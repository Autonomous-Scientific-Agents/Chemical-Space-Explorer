//! Skeleton graphs as adjacency bitmasks (`adj[v]` = neighbors of v), mirroring graph.py.

#[inline]
pub fn bits(mut m: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        if m == 0 {
            None
        } else {
            let v = m.trailing_zeros() as usize;
            m &= m - 1;
            Some(v)
        }
    })
}

pub fn nbr_lists(adj: &[u64]) -> Vec<Vec<usize>> {
    adj.iter().map(|&a| bits(a).collect()).collect()
}

pub fn edges(adj: &[u64]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for u in 0..adj.len() {
        for v in bits(adj[u]) {
            if u < v {
                out.push((u, v));
            }
        }
    }
    out
}

pub fn n_edges(adj: &[u64]) -> usize {
    adj.iter().map(|a| a.count_ones() as usize).sum::<usize>() / 2
}

/// Vertices whose removal leaves the graph connected (complement of articulation points).
pub fn non_cut_vertices(adj: &[u64]) -> Vec<usize> {
    let n = adj.len();
    if n <= 2 {
        return (0..n).collect();
    }
    let nbrs = nbr_lists(adj);
    let mut disc = vec![-1i64; n];
    let mut low = vec![0i64; n];
    let mut cut = vec![false; n];
    disc[0] = 0;
    let mut t = 1i64;
    let mut root_children = 0;
    let mut stack: Vec<(usize, isize, usize)> = vec![(0, -1, 0)];
    while let Some(&(v, parent, idx)) = stack.last() {
        if idx < nbrs[v].len() {
            let w = nbrs[v][idx];
            stack.last_mut().unwrap().2 += 1;
            if disc[w] < 0 {
                disc[w] = t;
                low[w] = t;
                t += 1;
                stack.push((w, v as isize, 0));
                continue;
            }
            if w as isize != parent && disc[w] < low[v] {
                low[v] = disc[w];
            }
        } else {
            stack.pop();
            if parent >= 0 {
                let p = parent as usize;
                if low[v] < low[p] {
                    low[p] = low[v];
                }
                if p == 0 {
                    root_children += 1;
                } else if low[v] >= disc[p] {
                    cut[p] = true;
                }
            }
        }
    }
    cut[0] = root_children > 1;
    (0..n).filter(|&v| !cut[v]).collect()
}

/// All simple cycles with at most `max_len` vertices, each listed once (rooted at the smallest vertex).
pub fn simple_cycles(adj: &[u64], max_len: usize) -> Vec<Vec<usize>> {
    fn extend(adj: &[u64], s: usize, higher: u64, max_len: usize, path: &mut Vec<usize>, used: u64,
              out: &mut Vec<Vec<usize>>) {
        let v = *path.last().unwrap();
        if path.len() >= 3 && (adj[v] >> s) & 1 == 1 && path[1] < v {
            out.push(path.clone());
        }
        if path.len() == max_len {
            return;
        }
        for w in bits(adj[v] & higher & !used) {
            path.push(w);
            extend(adj, s, higher, max_len, path, used | 1 << w, out);
            path.pop();
        }
    }
    let mut out = Vec::new();
    for s in 0..adj.len() {
        let higher = if s + 1 >= 64 { 0 } else { u64::MAX << (s + 1) };
        let mut path = vec![s];
        extend(adj, s, higher, max_len, &mut path, 1 << s, &mut out);
    }
    out
}
