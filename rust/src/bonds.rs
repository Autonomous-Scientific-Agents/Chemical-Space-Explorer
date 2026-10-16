//! Bond-order payload, mirroring bonds.py: distinct (element, H) labelings that pass all rules.

use crate::canon::{canonical_form, refined_colors};
use crate::graph::edges;
use crate::rules::{edge_cap, molecule_ok, Cfg, RingInfo};
use std::collections::HashSet;

/// One molecule: hydrogen count per atom and one passing bond order per edge (in `edges` order).
pub struct Molecule {
    pub h: Vec<u32>,
    pub orders: Vec<u32>,
}

pub fn molecules(labels: &[u32], adj: &[u64], cfg: &Cfg) -> (Vec<(usize, usize)>, Vec<Molecule>) {
    let edge_list = edges(adj);
    if !cfg.skeleton_ok(labels, adj) {
        return (edge_list, Vec::new());
    }
    let n = labels.len();
    let rings = RingInfo::new(adj);
    let deg: Vec<u32> = adj.iter().map(|a| a.count_ones()).collect();
    let mut spare: Vec<i64> = (0..n).map(|v| cfg.max_valence(labels[v]) as i64 - deg[v] as i64).collect();
    let caps: Vec<u32> = edge_list.iter().map(|&(u, v)| edge_cap(cfg, labels, &rings, u, v)).collect();
    let mut index = vec![usize::MAX; n * n];
    for (i, &(u, v)) in edge_list.iter().enumerate() {
        index[u * n + v] = i;
        index[v * n + u] = i;
    }
    let mut seen: HashSet<Vec<u32>> = HashSet::new();
    let mut accepted: Vec<Molecule> = Vec::new();
    let mut extra = vec![0u32; edge_list.len()];

    struct Ctx<'a> {
        labels: &'a [u32],
        adj: &'a [u64],
        cfg: &'a Cfg,
        rings: &'a RingInfo,
        edge_list: &'a [(usize, usize)],
        caps: &'a [u32],
        deg: &'a [u32],
        index: &'a [usize],
    }

    fn finish(c: &Ctx, extra: &[u32], seen: &mut HashSet<Vec<u32>>, accepted: &mut Vec<Molecule>) {
        let n = c.labels.len();
        let mut bond_sum: Vec<u32> = c.deg.to_vec();
        for (&(u, v), &x) in c.edge_list.iter().zip(extra) {
            bond_sum[u] += x;
            bond_sum[v] += x;
        }
        let h: Vec<u32> = (0..n)
            .map(|v| c.cfg.valences(c.labels[v]).iter().find(|&&val| val >= bond_sum[v]).unwrap() - bond_sum[v])
            .collect();
        if seen.contains(&h) {
            return;
        }
        let orders: Vec<u32> = extra.iter().map(|x| x + 1).collect();
        let order = |u: usize, v: usize| orders[c.index[u * n + v]];
        if molecule_ok(c.cfg, c.labels, c.adj, c.rings, &order, c.edge_list, &orders, &h) {
            seen.insert(h.clone());
            accepted.push(Molecule { h, orders });
        }
    }

    fn assign(c: &Ctx, i: usize, extra: &mut Vec<u32>, spare: &mut Vec<i64>, seen: &mut HashSet<Vec<u32>>,
              accepted: &mut Vec<Molecule>) {
        if i == c.edge_list.len() {
            finish(c, extra, seen, accepted);
            return;
        }
        let (u, v) = c.edge_list[i];
        let top = (c.caps[i] as i64 - 1).min(spare[u]).min(spare[v]);
        for x in 0..top + 1 {
            extra[i] = x as u32;
            spare[u] -= x;
            spare[v] -= x;
            assign(c, i + 1, extra, spare, seen, accepted);
            spare[u] += x;
            spare[v] += x;
        }
        extra[i] = 0;
    }

    let ctx = Ctx { labels, adj, cfg, rings: &rings, edge_list: &edge_list, caps: &caps, deg: &deg, index: &index };
    assign(&ctx, 0, &mut extra, &mut spare, &mut seen, &mut accepted);
    let colors = refined_colors(labels, adj);
    let mut distinct = colors.clone();
    distinct.sort_unstable();
    distinct.dedup();
    if distinct.len() == n {
        return (edge_list, accepted); // no symmetry: every hydrogen pattern is a distinct molecule
    }
    let mut forms = HashSet::new();
    let out = accepted
        .into_iter()
        .filter(|m| {
            let lab: Vec<u32> = labels.iter().zip(&m.h).map(|(&l, &h)| l * 16 + h).collect();
            forms.insert(canonical_form(&lab, adj))
        })
        .collect();
    (edge_list, out)
}

pub fn count_molecules(labels: &[u32], adj: &[u64], cfg: &Cfg) -> usize {
    molecules(labels, adj, cfg).1.len()
}
