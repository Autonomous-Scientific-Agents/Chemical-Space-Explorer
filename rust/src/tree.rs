//! Generation tree (McKay canonical augmentation), a port of tree.py so children come out identical
//! and in the same order.

use crate::bonds::count_molecules;
use crate::canon::{canonical_form, canonicalize, initial_colors, invariant_hash, refine, same_orbit, Form};
use crate::graph::{nbr_lists, non_cut_vertices};
use crate::rules::Cfg;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

pub type Node = (Vec<u32>, Vec<u64>);

fn add_vertex(labels: &[u32], adj: &[u64], label: u32, mask: u64) -> Node {
    let n = labels.len();
    let mut l = labels.to_vec();
    l.push(label);
    let mut a: Vec<u64> = adj.iter().enumerate().map(|(v, &x)| if mask >> v & 1 == 1 { x | 1 << n } else { x }).collect();
    a.push(mask);
    (l, a)
}

/// Canonical-deletion test: is `new` in the orbit of the canonical non-cut vertex?
fn accept(labels: &[u32], adj: &[u64], new: usize) -> (bool, Option<Form>) {
    let deg: Vec<u32> = adj.iter().map(|a| a.count_ones()).collect();
    let key_new = (labels[new], deg[new]);
    if (0..labels.len()).any(|v| deg[v] == 1 && (labels[v], 1) < key_new) {
        return (false, None);
    }
    let nc = non_cut_vertices(adj);
    if nc.iter().map(|&v| (labels[v], deg[v])).min().unwrap() < key_new {
        return (false, None);
    }
    let cands: Vec<usize> = nc.iter().copied().filter(|&v| (labels[v], deg[v]) == key_new).collect();
    if cands.len() == 1 {
        return (true, None);
    }
    let colors = refine(&nbr_lists(adj), initial_colors(labels));
    let cmin = cands.iter().map(|&v| colors[v]).min().unwrap();
    if colors[new] != cmin {
        return (false, None);
    }
    let cands: Vec<usize> = cands.into_iter().filter(|&v| colors[v] == cmin).collect();
    if cands.len() == 1 {
        return (true, None);
    }
    let (form, perm) = canonicalize(labels, adj);
    let m = *cands.iter().min_by_key(|&&v| perm[v]).unwrap();
    (same_orbit(labels, adj, m, new), Some(form))
}

fn combinations(items: &[usize], k: usize, start: usize, cur: &mut Vec<usize>, f: &mut dyn FnMut(&[usize])) {
    if cur.len() == k {
        f(cur);
        return;
    }
    for i in start..items.len() {
        if items.len() - i < k - cur.len() {
            break;
        }
        cur.push(items[i]);
        combinations(items, k, i + 1, cur, f);
        cur.pop();
    }
}

pub fn children(labels: &[u32], adj: &[u64], elements: &[u32], cfg: &Cfg) -> Vec<Node> {
    let n = labels.len();
    let spare: Vec<i64> = (0..n).map(|v| cfg.max_valence(labels[v]) as i64 - adj[v].count_ones() as i64).collect();
    let mut index: HashMap<u64, usize> = HashMap::new();
    let mut buckets: Vec<Vec<(Node, Option<Form>)>> = Vec::new();
    for &e in elements {
        let eligible: Vec<usize> = (0..n).filter(|&v| spare[v] > 0 && cfg.pair_allowed(e, labels[v])).collect();
        let max_size = (cfg.max_valence(e) as usize).min(eligible.len()).min(if cfg.allow_rings { n } else { 1 });
        for size in 1..=max_size {
            combinations(&eligible, size, 0, &mut Vec::new(), &mut |subset| {
                let mask = subset.iter().fold(0u64, |m, &v| m | 1 << v);
                let child = add_vertex(labels, adj, e, mask);
                let (ok, mut form) = accept(&child.0, &child.1, n);
                if !ok {
                    return;
                }
                // siblings can be isomorphic; canonicalize only on hash collisions
                let hsh = invariant_hash(&child.0, &child.1);
                match index.get(&hsh) {
                    Some(&bi) => {
                        let f = form.take().unwrap_or_else(|| canonical_form(&child.0, &child.1));
                        let bucket = &mut buckets[bi];
                        for entry in bucket.iter_mut() {
                            if entry.1.is_none() {
                                entry.1 = Some(canonical_form(&entry.0 .0, &entry.0 .1));
                            }
                        }
                        if bucket.iter().any(|entry| entry.1.as_ref() == Some(&f)) {
                            return;
                        }
                        bucket.push((child, Some(f)));
                    }
                    None => {
                        index.insert(hsh, buckets.len());
                        buckets.push(vec![(child, form)]);
                    }
                }
            });
        }
    }
    buckets.into_iter().flatten().map(|(c, _)| c).collect()
}

/// Breadth-first expansion of one node by `depth` levels: per level (molecules, skeletons), last level.
fn expand_one(node: &Node, depth: usize, elements: &[u32], cfg: &Cfg, keep_last: bool) -> (Vec<(u64, u64)>, Vec<Node>) {
    let mut level = vec![node.clone()];
    let mut rows = Vec::with_capacity(depth);
    for _ in 0..depth {
        level = level.iter().flat_map(|(l, a)| children(l, a, elements, cfg)).collect();
        let mols: u64 = level.iter().map(|(l, a)| count_molecules(l, a, cfg) as u64).sum();
        rows.push((mols, level.len() as u64));
    }
    if !keep_last {
        level.clear();
    }
    (rows, level)
}

/// Expand every node by `depth` levels on `threads` scoped threads (all joined before returning, so
/// Python may fork afterwards); rows are summed, last levels concatenated in node order.
pub fn expand(nodes: &[Node], depth: usize, elements: &[u32], cfg: &Cfg, keep_last: bool,
              threads: usize) -> (Vec<(u64, u64)>, Vec<Node>) {
    if threads <= 1 {
        // shard workers call this once per process; no extra thread needed
        let parts = nodes.iter().map(|nd| expand_one(nd, depth, elements, cfg, keep_last));
        return combine(parts, depth);
    }
    let next = AtomicUsize::new(0);
    let slots: Vec<Mutex<Option<(Vec<(u64, u64)>, Vec<Node>)>>> = nodes.iter().map(|_| Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1).min(nodes.len().max(1)) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= nodes.len() {
                    break;
                }
                *slots[i].lock().unwrap() = Some(expand_one(&nodes[i], depth, elements, cfg, keep_last));
            });
        }
    });
    combine(slots.into_iter().map(|m| m.into_inner().unwrap().unwrap()), depth)
}

fn combine(parts: impl Iterator<Item = (Vec<(u64, u64)>, Vec<Node>)>, depth: usize) -> (Vec<(u64, u64)>, Vec<Node>) {
    let mut rows = vec![(0u64, 0u64); depth];
    let mut last = Vec::new();
    for (r, l) in parts {
        for (acc, x) in rows.iter_mut().zip(r) {
            acc.0 += x.0;
            acc.1 += x.1;
        }
        last.extend(l);
    }
    (rows, last)
}
