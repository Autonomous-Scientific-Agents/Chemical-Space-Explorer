//! Kekule SMILES writer, a port of export.to_smiles (output is identical).

use crate::graph::nbr_lists;
use crate::rules::SYMBOLS;

const BOND: [&str; 4] = ["", "", "=", "#"];

pub fn to_smiles(labels: &[u32], adj: &[u64], edge_list: &[(usize, usize)], orders: &[u32], h: &[u32]) -> String {
    let n = labels.len();
    let nbrs = nbr_lists(adj);
    let mut order = vec![0u32; n * n];
    for (&(u, v), &o) in edge_list.iter().zip(orders) {
        order[u * n + v] = o;
        order[v * n + u] = o;
    }
    let mut state = vec![0u8; n];
    let mut kids: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut opens: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    let mut closes: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];

    #[allow(clippy::too_many_arguments)]
    fn dfs(v: usize, parent: isize, nbrs: &[Vec<usize>], state: &mut [u8], kids: &mut [Vec<usize>],
           opens: &mut [Vec<(usize, usize)>], closes: &mut [Vec<(usize, usize)>]) {
        state[v] = 1;
        for &w in &nbrs[v] {
            if w as isize == parent {
                continue;
            }
            if state[w] == 0 {
                kids[v].push(w);
                dfs(w, v as isize, nbrs, state, kids, opens, closes);
            } else if state[w] == 1 {
                opens[w].push((w, v));
                closes[v].push((w, v));
            }
        }
        state[v] = 2;
    }
    dfs(0, -1, &nbrs, &mut state, &mut kids, &mut opens, &mut closes);

    struct Emit<'a> {
        labels: &'a [u32],
        h: &'a [u32],
        n: usize,
        order: Vec<u32>,
        kids: Vec<Vec<usize>>,
        opens: Vec<Vec<(usize, usize)>>,
        closes: Vec<Vec<(usize, usize)>>,
        digit: std::collections::HashMap<(usize, usize), u32>,
        free: Vec<u32>,
        next_digit: u32,
    }
    impl Emit<'_> {
        fn label(d: u32) -> String {
            if d < 10 { d.to_string() } else { format!("%{d}") }
        }
        fn take(&mut self) -> u32 {
            if !self.free.is_empty() {
                self.free.sort_unstable();
                return self.free.remove(0);
            }
            self.next_digit += 1;
            self.next_digit - 1
        }
        fn emit(&mut self, v: usize) -> String {
            let sym = SYMBOLS[self.labels[v] as usize];
            let mut out = if sym == "Si" {
                match self.h[v] {
                    0 => "[Si]".to_string(),
                    1 => "[SiH]".to_string(),
                    k => format!("[SiH{k}]"),
                }
            } else {
                sym.to_string()
            };
            let mut closed = Vec::new();
            for ring in self.closes[v].clone() {
                let d = self.digit.remove(&ring).unwrap();
                out += BOND[self.order[ring.0 * self.n + ring.1] as usize];
                out += &Self::label(d);
                closed.push(d);
            }
            for ring in self.opens[v].clone() {
                let d = self.take();
                self.digit.insert(ring, d);
                out += &Self::label(d);
            }
            self.free.extend(closed);
            let kids = self.kids[v].clone();
            for (i, &w) in kids.iter().enumerate() {
                let branch = format!("{}{}", BOND[self.order[v * self.n + w] as usize], self.emit(w));
                if i == kids.len() - 1 {
                    out += &branch;
                } else {
                    out += &format!("({branch})");
                }
            }
            out
        }
    }
    let mut e = Emit {
        labels, h, n, order, kids, opens, closes,
        digit: std::collections::HashMap::new(), free: Vec::new(), next_digit: 1,
    };
    e.emit(0)
}
