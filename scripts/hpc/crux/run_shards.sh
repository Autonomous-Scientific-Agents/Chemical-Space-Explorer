#!/bin/bash
# Plan (if needed) and run a sharded chemuniverse job inside a PBS allocation on Crux.
#
#   run_shards.sh RUN_DIR [SUMMARY_DIR] -- [shard plan arguments]
#
# RUN_DIR      large outputs, e.g. /lus/grand/projects/<project>/chemuniverse/<name>
# SUMMARY_DIR  small copies (result.json, logs) on /home, readable through the IRI filesystem API
# Re-submitting the same command resumes: finished shards are skipped.
set -euo pipefail
# If a launcher (e.g. the IRI API's mpiexec) starts this once per node, only rank 0 drives the
# ensemble; ensemble_launcher itself reaches the other nodes listed in PBS_NODEFILE.
RANK=${PALS_RANKID:-${PMI_RANK:-0}}
if [ "$RANK" != 0 ]; then echo "rank $RANK on $(hostname): rank 0 drives the ensemble"; exit 0; fi
PREFIX=${CHEMUNIVERSE_PREFIX:-$HOME/chemuniverse-crux}
RUN_DIR=$1; shift
SUMMARY_DIR=$RUN_DIR
if [ "${1:-}" != "--" ] && [ $# -gt 0 ]; then SUMMARY_DIR=$1; shift; fi
[ "${1:-}" = "--" ] && shift
source "$PREFIX/venv/bin/activate"
mkdir -p "$RUN_DIR" "$SUMMARY_DIR"
LOG=$SUMMARY_DIR/run.log
{
  echo "== $(date -Is) job ${PBS_JOBID:-local} on $(wc -l < "${PBS_NODEFILE:-/dev/null}") nodes, chemuniverse $(git -C "$PREFIX/src" rev-parse --short HEAD)"
  if [ ! -f "$RUN_DIR/plan.json" ]; then python -m chemuniverse shard plan "$RUN_DIR" "$@"; fi
  cd "$RUN_DIR"   # ensemble_launcher writes logs/ into the working directory
  status=0
  python "$PREFIX/src/scripts/hpc/run_ensemble.py" "$RUN_DIR" --ncpus 128 || status=$?
  [ -f "$RUN_DIR/result.json" ] && [ "$SUMMARY_DIR" != "$RUN_DIR" ] && cp "$RUN_DIR/result.json" "$RUN_DIR/plan.json" "$SUMMARY_DIR/"
  [ $status -eq 0 ] && python -m chemuniverse shard reduce "$RUN_DIR"
  echo "== $(date -Is) exit $status"
  exit $status
} 2>&1 | tee -a "$LOG"
