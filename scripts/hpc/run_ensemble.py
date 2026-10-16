"""Run the pending shards of a chemuniverse run with ensemble_launcher, then reduce.

    python -m chemuniverse shard plan RUN_DIR --kind exact --max-atoms 13 --elements C N O F --shards 50000
    python scripts/hpc/run_ensemble.py RUN_DIR --ncpus 128      # inside a PBS job (nodes from PBS_NODEFILE)

Every task is one single-core shard (chemuniverse.shards.run_shard). Shards write
their own result files, so re-running after a timeout or preemption resumes.
Without PBS_NODEFILE the launcher runs on this machine only (for testing).
"""
import argparse
import json
import sys
import time
from pathlib import Path

from ensemble_launcher import EnsembleLauncher
from ensemble_launcher.config import LauncherConfig, PolicyConfig, SystemConfig
from ensemble_launcher.config.mpi_config import MPIConfig
from ensemble_launcher.ensemble import Task
from ensemble_launcher.helper_functions import get_nodes

from chemuniverse import shards


def launcher_config(nnodes, flavor):
    """None (auto-configure) unless an MPI flavor is given; then mirror the auto-configuration for
    serial tasks with that flavor's mpiexec flags."""
    if flavor is None or nnodes == 1:
        return None
    return LauncherConfig(
        child_executor_name="async_mpi",
        task_executor_name="async_processpool",
        comm_name="async_zmq",
        policy_config=PolicyConfig(nlevels=1 if nnodes <= 64 else 2),
        mpi_config=MPIConfig(flavor=flavor),
        master_logs=True,
        worker_logs=True,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run_dir")
    p.add_argument("--ncpus", type=int, default=128, help="worker processes per node (Crux: 128 physical cores)")
    p.add_argument("--limit", type=int, default=None, help="run at most this many pending shards")
    p.add_argument("--mpi-flavor", default=None,
                   help="mpiexec flavor for reaching other nodes, e.g. cray-pals on Crux (default: launcher's)")
    a = p.parse_args()
    run_dir = str(Path(a.run_dir).resolve())
    todo = shards.pending(run_dir)[: a.limit]
    print(f"{len(todo)} pending shards", flush=True)
    t0 = time.time()
    if todo:
        tasks = {f"shard-{i}": Task(task_id=f"shard-{i}", nnodes=1, ppn=1, executable=shards.run_shard,
                                    args=(run_dir, i)) for i in todo}
        EnsembleLauncher(ensemble_file=tasks, system_config=SystemConfig(name="chemuniverse", ncpus=a.ncpus),
                         launcher_config=launcher_config(len(get_nodes()), a.mpi_flavor)).run()
    left = shards.pending(run_dir)
    print(f"ensemble finished in {time.time() - t0:.1f} s; {len(left)} shards still pending", flush=True)
    if left:
        sys.exit(1)
    results, stats = shards.reduce(run_dir)
    print(json.dumps({"rows": [[r.n, r.count, r.se] for r in results], "stats": stats}), flush=True)


if __name__ == "__main__":
    main()
