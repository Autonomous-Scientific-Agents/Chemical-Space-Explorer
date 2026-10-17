//! Element table and stability rules, mirroring elements.py and rules.py.
//! Element codes are ranks of the symbols in Python string order, so code order == symbol order.

use crate::graph::{bits, edges, n_edges, simple_cycles};

pub const SYMBOLS: [&str; 11] = ["B", "Br", "C", "Cl", "F", "I", "N", "O", "P", "S", "Si"];
pub const B: u32 = 0;
pub const C: u32 = 2;
pub const N: u32 = 6;
pub const O: u32 = 7;
pub const P: u32 = 8;
pub const S: u32 = 9;
const MAX_RING: usize = 8;

pub fn code(symbol: &str) -> Option<u32> {
    SYMBOLS.iter().position(|&s| s == symbol).map(|i| i as u32)
}

pub fn is_halogen(e: u32) -> bool {
    matches!(e, 1 | 3 | 4 | 5)
}

#[derive(Clone, Copy, Debug)]
pub struct Cfg {
    pub hypervalent: bool,
    pub bond_order_limits: bool,
    pub strain: bool,
    pub heteroatom: bool,
    pub ban_n_o: bool,
    pub allow_rings: bool,
    pub allow_multiple_bonds: bool,
    // strict tier (ADR 0015)
    pub hetero_chains: bool,
    pub ban_s_o: bool,
    pub small_ring_heteroatoms: bool,
    pub tautomers: bool,
    pub thio_gem: bool,
    pub thiocarbonyl: bool,
}

pub fn is_hetero(e: u32) -> bool {
    e == N || e == O || e == S
}

impl Cfg {
    pub fn valences(&self, e: u32) -> &'static [u32] {
        if self.hypervalent && e == S {
            return &[2, 4, 6];
        }
        if self.hypervalent && e == P {
            return &[3, 5];
        }
        match e {
            B | N | P => &[3],
            C | 10 => &[4],
            O | S => &[2],
            _ => &[1],
        }
    }

    pub fn max_valence(&self, e: u32) -> u32 {
        *self.valences(e).last().unwrap()
    }

    pub fn max_bond_order(&self, e: u32) -> u32 {
        if !self.allow_multiple_bonds {
            return 1;
        }
        if !self.bond_order_limits {
            return 3;
        }
        match e {
            O | S => 2,
            e if is_halogen(e) => 1,
            _ => 3,
        }
    }

    pub fn pair_allowed(&self, a: u32, b: u32) -> bool {
        if self.ban_s_o && ((a == S && b == O) || (a == O && b == S)) {
            return false;
        }
        if !self.heteroatom {
            return true;
        }
        if a == O && b == O {
            return false;
        }
        if self.ban_n_o && ((a == N && b == O) || (a == O && b == N)) {
            return false;
        }
        for (x, y) in [(a, b), (b, a)] {
            if is_halogen(x) && (is_halogen(y) || y == N || y == O || y == S) {
                return false;
            }
        }
        true
    }

    /// Hereditary rules enforced while building skeletons.
    pub fn skeleton_ok(&self, labels: &[u32], adj: &[u64]) -> bool {
        if !self.allow_rings && n_edges(adj) >= labels.len() {
            return false;
        }
        if labels.iter().zip(adj).any(|(&l, a)| a.count_ones() > self.max_valence(l)) {
            return false;
        }
        if !edges(adj).iter().all(|&(u, v)| self.pair_allowed(labels[u], labels[v])) {
            return false;
        }
        (0..labels.len()).all(|v| self.hetero_ok(labels, adj, v))
    }

    /// Strict hereditary rules S1 and S3 around vertex v (see rules.hetero_ok).
    pub fn hetero_ok(&self, labels: &[u32], adj: &[u64], v: usize) -> bool {
        let het_v = is_hetero(labels[v]);
        if self.hetero_chains && het_v {
            let het_nbrs: Vec<usize> = bits(adj[v]).filter(|&u| is_hetero(labels[u])).collect();
            if het_nbrs.len() >= 2 {
                return false;
            }
            if het_nbrs.iter().any(|&u| bits(adj[u]).filter(|&w| is_hetero(labels[w])).count() >= 2) {
                return false;
            }
        }
        if self.small_ring_heteroatoms {
            let nbrs: Vec<usize> = bits(adj[v]).collect();
            for i in 0..nbrs.len() {
                for &b in &nbrs[i + 1..] {
                    let a = nbrs[i];
                    let base = het_v as u32 + is_hetero(labels[a]) as u32 + is_hetero(labels[b]) as u32;
                    if (adj[a] >> b) & 1 == 1 && base >= 2 {
                        return false; // triangle v-a-b
                    }
                    for c in bits(adj[a] & adj[b] & !(1u64 << v)) {
                        if base + is_hetero(labels[c]) as u32 >= 2 {
                            return false; // 4-ring v-a-c-b
                        }
                    }
                }
            }
        }
        true
    }
}

/// Ring facts needed by the strain rules.
pub struct RingInfo {
    n: usize,
    edge_ring: Vec<u8>, // n*n, smallest ring (<= MAX_RING) containing the edge, 0 = none
    pub atom_ring: Vec<u8>,
    pub bridgeheads: u64,
}

impl RingInfo {
    pub fn new(adj: &[u64]) -> RingInfo {
        let n = adj.len();
        let mut info = RingInfo { n, edge_ring: vec![0; n * n], atom_ring: vec![0; n], bridgeheads: 0 };
        if n_edges(adj) < n {
            return info;
        }
        let cycles = simple_cycles(adj, MAX_RING);
        for c in &cycles {
            let k = c.len() as u8;
            for i in 0..c.len() {
                let (u, v) = (c[i], c[(i + 1) % c.len()]);
                let e = u.min(v) * n + u.max(v);
                if info.edge_ring[e] == 0 || info.edge_ring[e] > k {
                    info.edge_ring[e] = k;
                }
                if info.atom_ring[u] == 0 || info.atom_ring[u] > k {
                    info.atom_ring[u] = k;
                }
            }
        }
        // Bredt: two rings <= 7 sharing a path of >= 3 atoms form a bridged system
        let small: Vec<(&Vec<usize>, u64)> = cycles
            .iter()
            .filter(|c| c.len() <= 7)
            .map(|c| (c, c.iter().fold(0u64, |m, &x| m | 1 << x)))
            .collect();
        for i in 0..small.len() {
            for j in i + 1..small.len() {
                let shared = small[i].1 & small[j].1;
                if shared.count_ones() < 3 {
                    continue;
                }
                for ring in [small[i].0, small[j].0] {
                    let k = ring.len();
                    for (q, &x) in ring.iter().enumerate() {
                        let prev = ring[(q + k - 1) % k];
                        let next = ring[(q + 1) % k];
                        if shared >> x & 1 == 1 && (shared >> prev & 1 == 0 || shared >> next & 1 == 0) {
                            info.bridgeheads |= 1 << x;
                        }
                    }
                }
            }
        }
        info
    }

    pub fn edge(&self, u: usize, v: usize) -> u8 {
        self.edge_ring[u.min(v) * self.n + u.max(v)]
    }
}

pub fn edge_cap(cfg: &Cfg, labels: &[u32], rings: &RingInfo, u: usize, v: usize) -> u32 {
    let mut cap = cfg.max_bond_order(labels[u]).min(cfg.max_bond_order(labels[v])).min(3);
    if cfg.strain {
        let size = rings.edge(u, v);
        if size != 0 {
            cap = cap.min(if size <= 4 { 1 } else if size < 8 { 2 } else { 3 });
        }
        if (rings.bridgeheads >> u | rings.bridgeheads >> v) & 1 == 1 {
            cap = 1;
        }
    }
    cap
}

/// Rules that need the full bond-order assignment; `order(u, v)` is the bond order of an edge.
pub fn molecule_ok(cfg: &Cfg, labels: &[u32], adj: &[u64], rings: &RingInfo, order: &dyn Fn(usize, usize) -> u32,
                   edge_list: &[(usize, usize)], orders: &[u32], h: &[u32]) -> bool {
    let n = labels.len();
    if cfg.strain {
        let mut multi = vec![0u32; n];
        for (&(u, v), &o) in edge_list.iter().zip(orders) {
            if o > 1 {
                multi[u] += 1;
                multi[v] += 1;
            }
        }
        for x in 0..n {
            if multi[x] >= 2 && rings.atom_ring[x] != 0 && rings.atom_ring[x] < 9 {
                return false; // cumulated bonds in a small ring
            }
        }
    }
    if cfg.tautomers || cfg.thio_gem || cfg.thiocarbonyl {
        for x in 0..n {
            if labels[x] == C && !strict_carbon_ok(cfg, labels, adj, x, order, h) {
                return false;
            }
        }
    }
    if cfg.heteroatom {
        for x in 0..n {
            if labels[x] != C {
                continue;
            }
            let nb: Vec<usize> = bits(adj[x]).collect();
            // enol: OH on a carbon that is double bonded to carbon
            if nb.iter().any(|&y| labels[y] == C && order(x, y) == 2)
                && nb.iter().any(|&y| labels[y] == O && h[y] > 0)
            {
                return false;
            }
            // hemiacetal / hemiaminal / gem-diol on an sp3 carbon
            if nb.iter().all(|&y| order(x, y) == 1) {
                let het: Vec<usize> = nb.iter().copied().filter(|&y| labels[y] == N || labels[y] == O).collect();
                if het.len() >= 2 && het.iter().any(|&y| h[y] > 0) {
                    return false;
                }
            }
        }
    }
    true
}

/// Strict molecule rules S4-S6 at carbon x (see rules.strict_carbon_ok).
fn strict_carbon_ok(cfg: &Cfg, labels: &[u32], adj: &[u64], x: usize, order: &dyn Fn(usize, usize) -> u32,
                    h: &[u32]) -> bool {
    let nbrs: Vec<usize> = bits(adj[x]).collect();
    let double: Vec<usize> = nbrs.iter().copied().filter(|&y| order(x, y) == 2).collect();
    let single_het: Vec<usize> = nbrs.iter().copied().filter(|&y| is_hetero(labels[y]) && order(x, y) == 1).collect();
    let single_het_h: Vec<usize> = single_het.iter().copied().filter(|&y| h[y] > 0).collect();
    if cfg.tautomers {
        let c_double = double.iter().any(|&y| labels[y] == C);
        if c_double && single_het_h.iter().any(|&y| labels[y] == S || labels[y] == N) {
            return false; // enethiol / NH-enamine
        }
        if !single_het.is_empty() && nbrs.iter().any(|&y| labels[y] == C && order(x, y) == 3) {
            return false; // ynol / ynamine / ynethiol
        }
        if double.iter().any(|&y| labels[y] == N) && single_het_h.iter().any(|&y| labels[y] == O || labels[y] == S) {
            return false; // imidic acid
        }
        if double.iter().any(|&y| labels[y] == S) && single_het_h.iter().any(|&y| labels[y] == O) {
            return false; // thiono acid
        }
        if c_double && double.iter().any(|&y| labels[y] == O) {
            return false; // ketene
        }
    }
    if cfg.thio_gem && double.is_empty() && nbrs.iter().all(|&y| order(x, y) == 1)
        && single_het.len() >= 2 && !single_het_h.is_empty() && single_het.iter().any(|&y| labels[y] == S)
    {
        return false; // hemithioacetal / gem-thiol
    }
    if cfg.thiocarbonyl && double.iter().any(|&y| labels[y] == S)
        && nbrs.iter().filter(|&&y| is_hetero(labels[y])).count() == 1
    {
        return false; // thioketone / thioaldehyde / thioketene
    }
    true
}
