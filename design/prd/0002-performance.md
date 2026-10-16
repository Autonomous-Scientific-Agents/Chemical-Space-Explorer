# PRD 0002: Performance

- Status: stages 1 and 2 done 2026-10-16; further single-core optimization of the Rust core is optional
- Date: 2026-10-16
- Decisions: ADR 0009

## Problem

Pure Python runs at about 130 µs per molecule. Exact enumeration stops at about 7 atoms; the site data stops at 12 atoms, which took about 15 minutes on 10 cores; dataset export to 11 atoms would take about a day.

## Goals

- Stage 1: use every core in `estimate()`, `enumerate_exact()` and `export` (`--jobs`), for about 8–10× on 10 cores.
- Stage 2: an optional Rust core (`chemuniverse-core`) with identical results, aiming for 30–100× on one core.
- Make the payload's ring analysis cheaper.

## Non-goals

- GPU support.
- Changing results. The Python backend stays the oracle.

## Requirements

- For a fixed seed, results are bit-identical across `jobs` values and across backends.
- `CHEMUNIVERSE_BACKEND=auto|python|rust`. The Rust core installs separately with `pip install ./rust`. It is not a pip extra, because it isn't on PyPI.
- `scripts/bench.py` reports molecules per second and seconds per pass, by backend and by jobs.

## Success metrics

The baseline (Python, one core, CNOF) is:
- exact enumeration to 7 atoms: 8.7 s;
- one stratified pass 6→12: about 25 s.

The table below is filled in from `scripts/bench.py` as each stage lands.

All timings are from `scripts/bench.py` on a 10-core Mac (8 performance + 2 efficiency cores), Python 3.11:

| Stage | exact CNOF ≤ 7 | exact CNOF ≤ 8 | 10 stratified passes 6→12 | export CNOF ≤ 7 | notes |
|---|---|---|---|---|---|
| baseline (jobs=1) | 9.8 s (7.9k mol/s) | 114.9 s | 193 s | 12.0 s (6.5k mol/s) | 1 core |
| stage 1 (jobs=10) | 1.9 s (40k mol/s) | 17.8 s | 32 s | 1.9 s (41k mol/s) | 5–6.5× |

Stage 1 gives 5–6.5×, short of the hoped-for 8–10×:
- The efficiency cores are about half speed, so 10 cores act like roughly 9.
- Serial passes share `cached_children`/`cached_count` across passes; separate processes don't, so a pass costs 32 s per core in a worker against 19 s in a serial run.
- Process start-up and the serial top of the tree add a fixed cost.
- Sending last-level nodes back to the parent costs almost nothing, and finer chunks gained about 10%, within noise.

### Stage 2: Rust core (ADR 0010)

Measured with `scripts/bench.py --jobs 1 10 --exact 8 --export 8`, with caches cleared before each task. Speedups are relative to Python on 1 core:

| Task (CNOF) | Python, 1 core | Python, 10 jobs | Rust, 1 core | Rust, 10 jobs |
|---|---|---|---|---|
| exact ≤ 8 (6.6×10⁵ molecules) | 112.7 s | 19.5 s | 4.8 s (23×) | 1.5 s (75×) |
| 10 stratified passes 6→12 | 192.6 s | 35.3 s | 13.0 s (15×) | 2.1 s (92×) |
| export ≤ 8 | 128.3 s | 18.8 s | 8.1 s (16×) | 1.2 s (107×) |

- The single-core gain (15–23×) is below the 30–100× estimate. The port keeps Python's data layout so that results stay identical: sorted `Vec` signatures in `refine` and `Vec`-keyed hash sets.
- In stratified sampling and export, the loops and the per-node calls across the language boundary still run in Python.
- Raw Rust `expand` (exact, no conversion back to Python) takes 4.1 s on 1 thread and 0.55 s on 10 threads for levels 6–8.
- **Next steps, if needed:**
  - Remove allocations from `refine` and `canonicalize`.
  - Move the export depth-first walk into Rust.
  - Port the sampling loops together with a bit-compatible Mersenne Twister.
