//! Python bindings. Skeletons cross the boundary as (labels: tuple[str], adj: tuple[int]) and the
//! rule configuration as the 7 RuleConfig booleans in field order.

pub mod bonds;
pub mod canon;
pub mod graph;
pub mod rules;
pub mod smiles;
pub mod tree;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rules::{Cfg, SYMBOLS};

// The hot paths allocate many small vectors; the macOS system allocator serializes threads on them.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// RuleConfig fields in declaration order (a list: PyO3 tuples stop at 12 elements).
type CfgTuple = Vec<bool>;
type PyNode = (Vec<&'static str>, Vec<u64>);

fn cfg_from(t: CfgTuple) -> Cfg {
    let f = |i: usize| t.get(i).copied().unwrap_or(false);
    Cfg {
        hypervalent: f(0),
        bond_order_limits: f(1),
        strain: f(2),
        heteroatom: f(3),
        ban_n_o: f(4),
        allow_rings: f(5),
        allow_multiple_bonds: f(6),
        hetero_chains: f(7),
        ban_s_o: f(8),
        small_ring_heteroatoms: f(9),
        tautomers: f(10),
        thio_gem: f(11),
        thiocarbonyl: f(12),
    }
}

fn codes(labels: &[String]) -> PyResult<Vec<u32>> {
    labels
        .iter()
        .map(|s| rules::code(s).ok_or_else(|| PyValueError::new_err(format!("unknown element {s:?}"))))
        .collect()
}

fn check(adj: &[u64]) -> PyResult<()> {
    if adj.len() > 63 {
        return Err(PyValueError::new_err("the Rust core supports at most 63 atoms"));
    }
    Ok(())
}

fn to_py(node: tree::Node) -> PyNode {
    (node.0.iter().map(|&c| SYMBOLS[c as usize]).collect(), node.1)
}

#[pyfunction]
fn children(labels: Vec<String>, adj: Vec<u64>, elements: Vec<String>, cfg: CfgTuple) -> PyResult<Vec<PyNode>> {
    check(&adj)?;
    let out = tree::children(&codes(&labels)?, &adj, &codes(&elements)?, &cfg_from(cfg));
    Ok(out.into_iter().map(to_py).collect())
}

#[pyfunction]
fn count_molecules(labels: Vec<String>, adj: Vec<u64>, cfg: CfgTuple) -> PyResult<usize> {
    check(&adj)?;
    Ok(bonds::count_molecules(&codes(&labels)?, &adj, &cfg_from(cfg)))
}

#[pyfunction]
fn molecule_smiles(labels: Vec<String>, adj: Vec<u64>, cfg: CfgTuple) -> PyResult<Vec<String>> {
    check(&adj)?;
    let labels = codes(&labels)?;
    let (edge_list, mols) = bonds::molecules(&labels, &adj, &cfg_from(cfg));
    Ok(mols.iter().map(|m| smiles::to_smiles(&labels, &adj, &edge_list, &m.orders, &m.h)).collect())
}

/// Expand nodes by `depth` levels on `threads` threads, releasing the GIL.
#[pyfunction]
#[pyo3(signature = (nodes, depth, elements, cfg, keep_last, threads=1))]
fn expand(py: Python<'_>, nodes: Vec<(Vec<String>, Vec<u64>)>, depth: usize, elements: Vec<String>, cfg: CfgTuple,
          keep_last: bool, threads: usize) -> PyResult<(Vec<(u64, u64)>, Vec<PyNode>)> {
    let nodes: Vec<tree::Node> = nodes.iter().map(|(l, a)| Ok((codes(l)?, a.clone()))).collect::<PyResult<_>>()?;
    let elements = codes(&elements)?;
    let cfg = cfg_from(cfg);
    let (rows, last) = py.allow_threads(|| tree::expand(&nodes, depth, &elements, &cfg, keep_last, threads));
    Ok((rows, last.into_iter().map(to_py).collect()))
}

#[pymodule]
fn chemuniverse_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(children, m)?)?;
    m.add_function(wrap_pyfunction!(count_molecules, m)?)?;
    m.add_function(wrap_pyfunction!(molecule_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(expand, m)?)?;
    Ok(())
}
