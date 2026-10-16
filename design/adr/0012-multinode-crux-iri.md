# ADR 0012: Multi-node runs on ALCF Crux with sharding, ensemble_launcher and IRI submission

- Status: accepted (implementation in progress; not yet run on Crux)
- Date: 2026-10-16

## Context

The user has supercomputer access and chose ALCF **Crux**, a CPU-only HPE Cray EX system.
- Each node has 2× AMD EPYC 7742 (Rome/Zen2) CPUs, 128 cores, 256 GB of memory and Slingshot-11.
- The scheduler is PBS Pro: `workq` allows up to 184 nodes for 24 h, `debug` 1–4 nodes for 2 h, and `preemptable` up to 10 nodes for 72 h.
- The login nodes are Milan/Zen3.

The user asked to use:
- [ensemble_launcher](https://github.com/argonne-lcf/ensemble_launcher), ALCF's hierarchical task-ensemble launcher;
- the **Genesis skills** `iri-api` skill, to submit jobs through the **DOE IRI Facility API** at ALCF (`api.alcf.anl.gov`, v1). Crux resource ID: `8b9b42f7-572a-4909-8472-a0453436304c`.

Every unit of work is independent: exact and export subtrees (thanks to canonical augmentation), stratified passes, and element sets. No communication is needed apart from a final sum.

## Decision

1. **Launcher-independent shards** (`chemuniverse/shards.py`, `python -m chemuniverse shard plan|run|reduce`):
   - `plan` grows the tree serially to a split level and deals the split-level nodes to shards round-robin, ordered by a cost proxy (number of children).
   - For `passes`, `plan` stores one seed per pass, drawn exactly as `stratified_estimate` draws them.
   - Each shard writes `results/<i//1000>/<i>.json` atomically, so reruns skip finished shards. This makes the `preemptable` queue and resubmission after a timeout safe.
   - `reduce` gives results identical to `enumerate_exact`, `estimate` and `export`, which `tests/test_shards.py` checks.
2. **ensemble_launcher** (`scripts/hpc/run_ensemble.py`) runs one single-core `Task` per pending shard, using `ncpus=128` per node and nodes from `PBS_NODEFILE`. chemuniverse itself does not depend on it.
3. **Submission through IRI** with the Genesis `iri-api` scripts:
   - `scripts/hpc/crux/iri_job.py` writes the JobSpecs: `probe`, `setup` and `run`. Each runs `/bin/bash -c` (ALCF allows at most two arguments), uses absolute paths (IRI passes `$HOME` literally) and sets `custom_attributes.filesystems = home:grand`.
   - Small outputs go to `/home/<user>/chemuniverse-runs/<name>`, which the IRI filesystem API can read (it exposes Home and Eagle, not Grand). Large outputs go to `/lus/grand/projects/<project>/chemuniverse/<name>`.
4. **Environment built on a compute node** (`scripts/hpc/crux/setup_env.sh`):
   - It runs through the ALCF proxy, with Python 3.11 from `uv` (ensemble_launcher needs ≥ 3.10).
   - rustup and the Rust core are built with `target-cpu=znver2`, because the login nodes are Zen3 while the compute nodes are Zen2.
   - It then runs the tests and a one-node benchmark.
5. `run_shards.sh` lets only rank 0 drive the ensemble, in case the IRI launcher starts the command once per node. A `probe` job checks this behavior before the first multi-node run.

## Alternatives

- PBS job arrays or mpi4py master-worker. Workable, but the user asked for ensemble_launcher, which already provides hierarchical scheduling, work stealing and fault handling at 2,000+ nodes.
- SSH and `qsub` instead of IRI. That needs an interactive MFA session for every action, whereas IRI works with a Globus token.
- One task per split-level node. Too many tiny tasks: launcher overhead is about 5 ms per task on one master, so shards should average at least 10–30 s.

## Consequences

- Measured locally: 400 shards on 10 workers took 8.9 s of wall time for 64 core-seconds of work, with results identical to `enumerate_exact`.
- At 5,000 tasks of about 1 ms each, launcher overhead dominated (27 s), so shard sizes must be chosen for each run.
- Open items, to be settled by the probe and setup jobs:
  - how IRI launches multi-node jobs;
  - whether `ncpus=256` must be requested explicitly;
  - the Python and proxy behavior on compute nodes;
  - the real per-core speed on Zen2.
