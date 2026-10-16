# ADR 0010: Rust core design: an exact port, sampling loops kept in Python, scoped threads, mimalloc

- Status: accepted. It refines ADR 0009, which named rayon.
- Date: 2026-10-16

## Context

ADR 0009 chose an optional Rust core. Building it raised four questions:
- How strict must parity with Python be?
- Where should the boundary between Rust and Python sit?
- How should threading coexist with the `fork`-based worker processes?
- Why did threads scale badly at first?

## Decision

1. **Exact port, not a reimplementation.**
   - `rust/src/{graph,canon,rules,bonds,tree,smiles}.rs` follow the Python code step by step: the same colour refinement, the same individualization–refinement search order, the same canonical deletion vertex, the same sibling bucketing and the same SMILES walk.
   - Element codes are the ranks of the symbols in Python string order, so every comparison orders the same way.
   - As a result the tree is identical, children come in the same order and SMILES strings match, not just the counts. `tests/test_backend_parity.py` checks this across 5 rule configurations and 4 element sets.
2. **Boundary.**
   - Rust provides `children`, `count_molecules`, `molecule_smiles` and a bulk `expand(nodes, depth, threads)` for exact mode.
   - The stratified and Knuth sampling loops stay in Python. That keeps the RNG (`random.Random`) and the merge order identical, so estimates are bit-identical across backends.
   - `CHEMUNIVERSE_BACKEND=auto|python|rust` selects the backend (`chemuniverse/backend.py`).
3. **Scoped `std::thread` workers with an atomic work index instead of rayon.**
   - All threads are joined before `expand` returns.
   - Python later `fork`s worker processes for the samplers, and forking while a rayon pool's idle threads are still alive risks deadlocks in the child.
4. **mimalloc as the global allocator.**
   - The hot paths allocate many small `Vec`s. With the macOS system allocator, 10 threads reached only 3.9×; with mimalloc, 8 threads reach 7.2×.
   - mimalloc also made one thread 1.5× faster.
5. **Optional, separately built package.**
   - Install with `pip install ./rust`. It is not a pip extra, because it isn't published on PyPI.
   - The pure-Python package still installs and runs without Rust.

## Alternatives

- Moving the sampling loops into Rust. This would be faster, but it needs a Rust port of CPython's Mersenne Twister and `random()` to keep estimates bit-identical. It may come later.
- A rayon global pool with fork protection. More fragile.
- Profile-guided reduction of allocations instead of mimalloc. Still worth doing later for more single-core speed.

## Consequences

- Measured single-core speedups of 15–23× and 10-core speedups of 75–107× against pure Python on one core (PRD 0002). That falls short of the 30–100× single-core estimate in ADR 0009, because the port keeps Python's data layout: sorted `Vec` signatures and `Vec`-keyed hash sets.
- Any change to the Python algorithm must be mirrored in Rust, or the parity test fails. The Python code stays the specification.
