//! Canonical labeling, a line-by-line port of canon.py so forms and permutations are identical.
//! Labels are u32 codes whose order matches the Python objects' sort order.

use crate::graph::nbr_lists;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

pub type Form = (Vec<u32>, Vec<u64>);

pub fn initial_colors(labels: &[u32]) -> Vec<u32> {
    let mut order: Vec<u32> = labels.to_vec();
    order.sort_unstable();
    order.dedup();
    labels.iter().map(|l| order.binary_search(l).unwrap() as u32).collect()
}

fn n_distinct(colors: &[u32]) -> usize {
    let mut c = colors.to_vec();
    c.sort_unstable();
    c.dedup();
    c.len()
}

/// Coarsest equitable refinement; colors are ranks of (color, sorted neighbor colors).
pub fn refine(nbrs: &[Vec<usize>], mut colors: Vec<u32>) -> Vec<u32> {
    let mut k = n_distinct(&colors);
    loop {
        let sigs: Vec<(u32, Vec<u32>)> = nbrs
            .iter()
            .enumerate()
            .map(|(v, nb)| {
                let mut s: Vec<u32> = nb.iter().map(|&u| colors[u]).collect();
                s.sort_unstable();
                (colors[v], s)
            })
            .collect();
        let mut uniq = sigs.clone();
        uniq.sort_unstable();
        uniq.dedup();
        colors = sigs.iter().map(|s| uniq.binary_search(s).unwrap() as u32).collect();
        if uniq.len() == k {
            return colors;
        }
        k = uniq.len();
    }
}

pub fn refined_colors(labels: &[u32], adj: &[u64]) -> Vec<u32> {
    refine(&nbr_lists(adj), initial_colors(labels))
}

/// Weisfeiler-Lehman hash; only equality matters (siblings with equal hashes are canonicalized).
pub fn invariant_hash(labels: &[u32], adj: &[u64]) -> u64 {
    fn h<T: Hash>(x: &T) -> u64 {
        let mut s = DefaultHasher::new();
        x.hash(&mut s);
        s.finish()
    }
    let nbrs = nbr_lists(adj);
    let mut hs: Vec<u64> = labels.iter().map(h).collect();
    for _ in 0..3 {
        hs = nbrs
            .iter()
            .enumerate()
            .map(|(v, nb)| {
                let mut s: Vec<u64> = nb.iter().map(|&u| hs[u]).collect();
                s.sort_unstable();
                h(&(hs[v], s))
            })
            .collect();
    }
    hs.sort_unstable();
    h(&hs)
}

fn orbit_of(v: usize, gens: &[&Vec<usize>], n: usize) -> Vec<bool> {
    let mut orbit = vec![false; n];
    orbit[v] = true;
    let mut stack = vec![v];
    while let Some(x) = stack.pop() {
        for g in gens {
            let y = g[x];
            if !orbit[y] {
                orbit[y] = true;
                stack.push(y);
            }
        }
    }
    orbit
}

struct Search<'a> {
    n: usize,
    nbrs: &'a [Vec<usize>],
    best: Option<(Vec<u64>, Vec<usize>)>,
    first: Option<(Vec<u64>, Vec<usize>)>,
    autos: Vec<Vec<usize>>,
}

impl<'a> Search<'a> {
    fn leaf(&mut self, colors: &[u32]) {
        let n = self.n;
        let perm: Vec<usize> = colors.iter().map(|&c| c as usize).collect();
        let mut inv = vec![0usize; n];
        for (v, &p) in perm.iter().enumerate() {
            inv[p] = v;
        }
        let cert: Vec<u64> = (0..n)
            .map(|p| self.nbrs[inv[p]].iter().map(|&u| 1u64 << perm[u]).sum())
            .collect();
        for r in [&self.first, &self.best] {
            if let Some((rc, rp)) = r {
                if *rc == cert {
                    let mut rinv = vec![0usize; n];
                    for (v, &p) in rp.iter().enumerate() {
                        rinv[p] = v;
                    }
                    self.autos.push((0..n).map(|v| rinv[perm[v]]).collect());
                    break;
                }
            }
        }
        if self.first.is_none() {
            self.first = Some((cert.clone(), perm.clone()));
        }
        if self.best.as_ref().map_or(true, |(bc, _)| cert < *bc) {
            self.best = Some((cert, perm));
        }
    }

    fn search(&mut self, colors: Vec<u32>, path: &mut Vec<usize>) {
        let n = self.n;
        let mut counts = vec![0usize; 2 * n + 2];
        for &c in &colors {
            counts[c as usize] += 1;
        }
        let target = (0..counts.len()).find(|&c| counts[c] > 1);
        let Some(target) = target else {
            self.leaf(&colors);
            return;
        };
        let target = target as u32;
        let mut tried = vec![false; n];
        for v in 0..n {
            if colors[v] != target || tried[v] {
                continue;
            }
            let child: Vec<u32> = colors
                .iter()
                .enumerate()
                .map(|(u, &c)| 2 * c + (c == target && u != v) as u32)
                .collect();
            let refined = refine(self.nbrs, child);
            path.push(v);
            self.search(refined, path);
            path.pop();
            let gens: Vec<&Vec<usize>> = self.autos.iter().filter(|g| path.iter().all(|&p| g[p] == p)).collect();
            for (w, inorb) in orbit_of(v, &gens, n).into_iter().enumerate() {
                if inorb {
                    tried[w] = true;
                }
            }
        }
    }
}

/// (form, perm): canonical form and perm[v] = canonical position of v.
pub fn canonicalize(labels: &[u32], adj: &[u64]) -> (Form, Vec<usize>) {
    let nbrs = nbr_lists(adj);
    let mut s = Search { n: labels.len(), nbrs: &nbrs, best: None, first: None, autos: Vec::new() };
    let start = refine(&nbrs, initial_colors(labels));
    s.search(start, &mut Vec::new());
    let (cert, perm) = s.best.unwrap();
    let mut sorted = labels.to_vec();
    sorted.sort_unstable();
    ((sorted, cert), perm)
}

pub fn canonical_form(labels: &[u32], adj: &[u64]) -> Form {
    canonicalize(labels, adj).0
}

/// True if some automorphism maps u to v.
pub fn same_orbit(labels: &[u32], adj: &[u64], u: usize, v: usize) -> bool {
    if u == v {
        return true;
    }
    let mu: Vec<u32> = labels.iter().enumerate().map(|(w, &l)| 2 * l + (w == u) as u32).collect();
    let mv: Vec<u32> = labels.iter().enumerate().map(|(w, &l)| 2 * l + (w == v) as u32).collect();
    canonical_form(&mu, adj) == canonical_form(&mv, adj)
}
