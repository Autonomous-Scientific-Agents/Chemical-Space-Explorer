# ADR 0012: Multi-node runs on ALCF Crux with sharding, ensemble_launcher and IRI submission

- Status: accepted. In production on Crux since 2026-10-16: exact CNOF ≤ 12 and the CNOS ≤ 12 basic and strict datasets.
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

## Findings from the probe job (362402, 2 nodes, 2026-10-16)

- IRI runs the command **once**, on the first node (no PALS or PMI rank variables), and `PBS_NODEFILE` lists all nodes. So ensemble_launcher must reach the other nodes itself.
- The environment is bare: `PATH=/bin:/usr/bin`, no modules, **no `mpiexec`**, and system Python 3.6.
  - The JobSpecs therefore use a login shell (`bash -lc`).
  - `run_shards.sh` loads `cray-pals` if it's needed, and fails clearly on more than one node without `mpiexec`.
  - The driver passes `MPIConfig(flavor="cray-pals")` (`mpiexec -n … --ppn …`), because ensemble_launcher's default flavor is `mpich`.
- Nodes are AMD EPYC 7742 with `nproc=256` (SMT2). `/lus/grand/projects/IQC` exists.

## Findings from the first runs (2026-10-16)

| Job | Problem | Fix |
|---|---|---|
| 362409 | `pthread_create` EAGAIN. Crux defaults to `OMP_NUM_THREADS=256`, and numpy's OpenBLAS (imported through ensemble_launcher) starts 64 threads in each of 128 workers | `OMP/OPENBLAS/MKL_NUM_THREADS=1`; Rust `expand` runs inline for one thread |
| 362410, 362412 | All shards ran on one host. `PolicyConfig` defaults to `nchildren=1`, so a single worker owned every node | One worker per node (`nchildren = nnodes`); `fixed_leafs_children_policy` beyond 64 nodes |
| 362412 | Exit 141 from `find \| head` under `pipefail` | No `head` in pipelines |
| 362413 | Hung until walltime, with the logs on Grand (not IRI-readable) | Logs on `/home`; optional ensemble `timeout` |
| 362420, 362421 | Workers started by `mpiexec` (a bare `python`) could not import ensemble_launcher, even with `PATH` forwarded via PALS `--env` | Export `PYTHONPATH` = venv site-packages + checkout, forwarded with `--env` |
| 362422 | 2 nodes, 500 shards per host, counts exact | (works) |
| 362423 | CNOF ≤ 11 on 4 nodes: 62.9 s against 241.5 s on 1 node, so **3.84× (96%)**; 500 shards per host | — |
| 362424 | **CNOF ≤ 12 exact on 4 nodes: 7,605,193,456 at n = 12, 8.303×10⁹ total**. 310,119 core-seconds (86 core-hours) in 761 s. The slowest shard took 735 s, so efficiency was 80%, limited by one straggler | Split deeper by default (≥ 16 frontier nodes per shard; `--split-level`) |

`362410` doubles as a 1-node baseline: exact C·N·O·F up to 11 atoms took 27,180 core-seconds in 241.5 s on 128 cores, about 88% efficiency. It produced the first exact counts beyond the old 7-atom limit: **56,997,371** molecules at 10 atoms and **634,889,614** at 11 (6.981×10⁸ in total up to 11).

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
