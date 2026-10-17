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
# ensemble_launcher starts workers on the other nodes with mpiexec (cray-pals)
if ! command -v mpiexec >/dev/null; then
  set +u  # profile and module scripts reference unset variables
  source /etc/profile >/dev/null 2>&1 || true
  module load cray-pals >/dev/null 2>&1 || true
  set -u
fi
NNODES=$(wc -l < "${PBS_NODEFILE:-/dev/null}")
if [ "$NNODES" -gt 1 ] && ! command -v mpiexec >/dev/null; then
  echo "error: $NNODES nodes but no mpiexec; load cray-pals (module load cray-pals)" >&2; exit 2
fi
source "$PREFIX/venv/bin/activate"
# one thread per worker: Crux defaults to OMP_NUM_THREADS=256, and numpy's OpenBLAS (imported by
# ensemble_launcher) starts a thread pool in each of the 128 workers -> pthread_create EAGAIN (job 362409)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1
# Workers that mpiexec (PALS) starts as a bare `python` resolve to the venv's symlinked interpreter but
# do not see the venv (job 362421). Point them at its site-packages and at the editable checkout
# (.pth files are ignored for PYTHONPATH entries); run_ensemble.py forwards PYTHONPATH with --env.
export PYTHONPATH="$(python -c 'import site; print(":".join(site.getsitepackages()))'):$PREFIX/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$RUN_DIR" "$SUMMARY_DIR"
LOG=$SUMMARY_DIR/run.log
{
  echo "== $(date -Is) job ${PBS_JOBID:-local} on $(wc -l < "${PBS_NODEFILE:-/dev/null}") nodes, chemuniverse $(git -C "$PREFIX/src" rev-parse --short HEAD)"
  echo "== mpiexec: $(command -v mpiexec || echo none); nodes: $(sort -u "${PBS_NODEFILE:-/dev/null}" | cut -d. -f1 | tr '\n' ' ')"
  if [ ! -f "$RUN_DIR/plan.json" ]; then python -m chemuniverse shard plan "$RUN_DIR" "$@"; fi
  # ensemble_launcher writes logs/ into the working directory: keep them on /home (IRI-readable, live)
  cd "$SUMMARY_DIR"
  status=0
  # optional cap (seconds) so a hung ensemble still leaves time to collect logs
  ${CHEMUNIVERSE_TIMEOUT:+timeout $CHEMUNIVERSE_TIMEOUT} python -u "$PREFIX/src/scripts/hpc/run_ensemble.py" "$RUN_DIR" \
    --ncpus 128 --mpi-flavor cray-pals || status=$?
  [ -f "$RUN_DIR/result.json" ] && [ "$SUMMARY_DIR" != "$RUN_DIR" ] && cp "$RUN_DIR/result.json" "$RUN_DIR/plan.json" "$SUMMARY_DIR/"
  find "$SUMMARY_DIR/logs" -type f > "$SUMMARY_DIR/logs_index.txt" 2>/dev/null || true
  python -c "import json,collections,glob; c=collections.Counter(json.load(open(f))['host'] for f in glob.glob('$RUN_DIR/results/*/*.json')); print('== shards per host:', dict(c))" || true
  [ $status -eq 0 ] && python -m chemuniverse shard reduce "$RUN_DIR"
  # export runs: merge SMILES into ~2 GB checksummed chunks + dataset/manifest.json (ADR 0014)
  if [ $status -eq 0 ] && [ "$(python -c "import json; print(json.load(open('$RUN_DIR/plan.json'))['kind'])")" = export ]; then
    mkdir -p "$RUN_DIR/dataset"
    command -v lfs >/dev/null && lfs setstripe -c 8 "$RUN_DIR/dataset" || true
    python -m chemuniverse shard finalize "$RUN_DIR" --delete-parts --jobs 64 || status=$?
    [ -f "$RUN_DIR/dataset/manifest.json" ] && cp "$RUN_DIR/dataset/manifest.json" "$SUMMARY_DIR/"
  fi
  echo "== $(date -Is) exit $status"
  exit $status
} 2>&1 | tee -a "$LOG"
