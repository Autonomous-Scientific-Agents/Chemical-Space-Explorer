# PRD 0003: Multi-node runs on Crux

- Status: in progress
- Date: 2026-10-16
- Decisions: ADR 0012 (ADR 0009 and ADR 0010 for performance)

## Problem

On a single 10-core Mac, the Rust core reaches exact C·N·O·F counts to about 10 atoms and sampling to about 17. A Crux allocation provides up to 184 × 128 = 23,552 cores, about 2,000× more.

## Goals

- Run any chemuniverse job type (exact counts, SMILES export, stratified passes) across many Crux nodes, with results identical to a single-machine run.
- Submit, monitor and collect results through the IRI API from a laptop agent, using the Genesis `iri-api` skill.
- Make runs resumable after a timeout or preemption.

## Non-goals

- GPUs; Crux has none.
- Moving TB-scale SMILES datasets off ALCF. That would use Globus later.
- Publishing a package to PyPI.

## Requirements

1. `chemuniverse shard plan|run|reduce` (done; tested against the single-machine functions).
2. `scripts/hpc/run_ensemble.py` (ensemble_launcher driver; done; tested locally).
3. `scripts/hpc/crux/{setup_env.sh,run_shards.sh,iri_job.py}` (written; not yet run on Crux).
4. A probe job (launch behavior), a setup job (environment, tests, benchmark) and then a 1→4-node scaling run in `debug`.

## Success metrics

- The setup job passes the test suite on Crux, and the parity tests agree with the Mac results.
- Scaling efficiency of at least 80% from 1 to 4 nodes in `debug`, and then up to tens of nodes in `workq`.
- First science run (target to be chosen by the user):
  - exact C·N·O·F counts to 13–14 atoms;
  - or a SMILES dataset;
  - or thousands of passes to 20+ atoms.

## Open questions

- Which ALCF project or allocation, and how many node-hours.
- The first science target.
- Whether IRI's generated batch script needs `ncpus=256` in `select`.
