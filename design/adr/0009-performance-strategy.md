# ADR 0009: Performance strategy of process parallelism plus an optional Rust core; no GPU, no Python threads

- Status: accepted
- Date: 2026-10-16

## Context

Measured on a 10-core Mac with Python 3.11:
- Exact CNOF runs at about 130 µs per molecule. The molecule count grows about 7.7× per added atom and the time about 10×, so 7 atoms take 8.7 s.
- One stratified pass to 12 atoms takes about 25 s. Strata grow only about 1.35× per level, but the time per level doubles, because the cost per node (payload backtracking plus ring enumeration) grows.
- The profile is dominated by interpreter overhead in tiny bit and integer functions. The payload takes about two thirds of the time.

## Decision

1. **Processes, not threads.** Stratified passes, Knuth probe batches and exact subtrees run in a `ProcessPoolExecutor` (`jobs`). Python threads are blocked by the GIL.
2. **Rust core** (PyO3 + maturin, rayon) as the optional package `chemuniverse-core`, chosen by the user on 2026-10-16. The pure-Python code remains the reference and the test oracle.
   - Rust was preferred over C for the safe Python bindings, data-race-free threading and memory safety in recursive search code.
   - C would have allowed direct use of nauty; we judged that unnecessary for molecules this small.
3. **No GPU.** The work is branchy and irregular (variable graphs, individualization–refinement, backtracking), and comparable tools such as nauty and surge are CPU-only.
4. **Algorithmic fix:** compute the smallest ring sizes without enumerating every simple cycle.
5. **Process start method `fork`** where available (macOS, Linux). `spawn` re-imports `__main__`, so it fails from notebooks and from scripts read on stdin, and it starts more slowly. The code is pure Python with no threads, so `fork` is safe. Windows falls back to the default.
6. **Results must not depend on `jobs`.** Each stratified pass and Knuth probe gets its own seed drawn from the caller's RNG. Exact subtrees are concatenated in their serial order, and export uses a fixed split into 128 subtree tasks. For a given seed, this changed the estimates once compared with the earlier single-RNG code; they still agree within error.

## Alternatives

- Cython or numba. Partial gains; numba handles recursive, object-heavy code poorly.
- Calling nauty only for canonical labeling. Gains at most about 1.5×, because the payload dominates.
- Free-threaded Python 3.13t. Immature, and the per-thread interpreter overhead remains.

## Consequences

- Stage 1 (processes) measured 5–6.5× on 10 cores (PRD 0002).

- Expected total speedup is about 400×: roughly 50× from Rust times roughly 8–10× from cores. This is to be measured with `scripts/bench.py` and recorded in PRD 0002.
- No implementation changes the exponential growth of exact enumeration; a 10× speedup buys only about one more atom.
- Rust needs a toolchain (rustup) to build, which pure-Python users don't need.
