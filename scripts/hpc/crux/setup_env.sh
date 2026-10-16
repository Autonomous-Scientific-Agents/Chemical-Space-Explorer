#!/bin/bash
# Build (or update) the chemuniverse environment for ALCF Crux. Run it on a COMPUTE node:
# the login nodes are Milan/Zen3 and the compute nodes Rome/Zen2, so the Rust core is
# compiled there with target-cpu=znver2. Idempotent; reruns update the checkout.
#
#   CHEMUNIVERSE_PREFIX  install prefix (default $HOME/chemuniverse-crux)
#   CHEMUNIVERSE_REF     git ref of Chemical-Universe to install (default main)
set -euo pipefail
PREFIX=${CHEMUNIVERSE_PREFIX:-$HOME/chemuniverse-crux}
REF=${CHEMUNIVERSE_REF:-main}
REPO=https://github.com/Autonomous-Scientific-Agents/Chemical-Universe.git

# compute nodes reach the internet only through the ALCF proxy
export http_proxy=http://proxy.alcf.anl.gov:3128 https_proxy=http://proxy.alcf.anl.gov:3128
export ftp_proxy=$http_proxy no_proxy=localhost,127.0.0.1,.alcf.anl.gov
export HTTP_PROXY=$http_proxy HTTPS_PROXY=$https_proxy

mkdir -p "$PREFIX/bin"
cd "$PREFIX"
echo "== host $(hostname)  $(grep -m1 'model name' /proc/cpuinfo)  nproc=$(nproc)"

# Python >= 3.10 via uv (independent of the cray-python version)
export PATH=$PREFIX/bin:$PATH UV_CACHE_DIR=$PREFIX/uv-cache UV_PYTHON_INSTALL_DIR=$PREFIX/python
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="$PREFIX/bin" INSTALLER_NO_MODIFY_PATH=1 sh
[ -d venv ] || uv venv --python 3.11 venv
source venv/bin/activate

# source checkout
if [ -d src/.git ]; then
  git -C src fetch -q origin && git -C src checkout -q "$REF" && git -C src pull -q --ff-only || true
else
  git clone -q "$REPO" src && git -C src checkout -q "$REF"
fi
echo "== chemuniverse $(git -C src rev-parse --short HEAD)"
uv pip install -q -e "src[test]" "ensemble_launcher @ git+https://github.com/argonne-lcf/ensemble_launcher.git@main"

# Rust toolchain inside the prefix, then the core built for Zen2
export RUSTUP_HOME=$PREFIX/rustup CARGO_HOME=$PREFIX/cargo PATH=$PREFIX/cargo/bin:$PATH
[ -x "$CARGO_HOME/bin/cargo" ] || curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --no-modify-path
RUSTFLAGS="-C target-cpu=znver2" uv pip install -q --reinstall ./src/rust

# check: Rust core present, test suite, a short single-node benchmark (one BLAS/OpenMP thread per process)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -c "import chemuniverse.backend as b; assert b.core, 'Rust core not importable'; print('== rust core ok')"
(cd src && CHEMUNIVERSE_BACKEND=rust python -m pytest -q -x)
(cd src && python scripts/bench.py --jobs 1 64 128 --exact 8 --export 8)
echo "== setup done: source $PREFIX/venv/bin/activate"
