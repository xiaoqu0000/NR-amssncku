#!/usr/bin/env bash
# AMSS-NCKU portable launcher.
# Reads config.env (if present) to pick the toolchain and MPI launcher, then
# runs the Python driver. Safe to commit; per-machine overrides go in config.env.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ -f config.env ]]; then
    # shellcheck disable=SC1091
    source config.env
fi

# Optional spack environment activation. Skipped unless SPACK_ENV_NAME is set
# AND spack is on PATH — open-source users never need this.
if [[ -n "${SPACK_ENV_NAME:-}" ]] && command -v spack >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(spack location -r)/share/spack/setup-env.sh"
    spack env activate "$SPACK_ENV_NAME"
fi

ulimit -s unlimited

export COMPILER="${COMPILER:-gnu}"
export MPI_LAUNCHER="${MPI_LAUNCHER:-auto}"
export ENABLE_CUDA="${ENABLE_CUDA:-no}"

python3 AMSS_NCKU_Program.py "$@"
