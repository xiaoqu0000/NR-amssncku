#!/usr/bin/env bash
# Preflight checker for AMSS-NCKU. Verifies the build/run prerequisites and
# prints the exact install commands for whatever is missing.
# Read-only: this script never installs anything itself.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ -f config.env ]]; then
    # shellcheck disable=SC1091
    source config.env
fi
COMPILER="${COMPILER:-gnu}"

missing_apt=()
missing_pip=()
fatal=0

check_cmd() {
    local cmd="$1" apt_pkg="$2" hint="${3:-}"
    if command -v "$cmd" >/dev/null 2>&1; then
        printf "  [ok]   %-16s -> %s\n" "$cmd" "$(command -v "$cmd")"
    else
        printf "  [MISS] %-16s (install: %s)%s\n" \
            "$cmd" "$apt_pkg" "${hint:+  $hint}"
        missing_apt+=("$apt_pkg")
        fatal=1
    fi
}

check_python_module() {
    local mod="$1" pip_pkg="${2:-$1}"
    if python3 -c "import $mod" 2>/dev/null; then
        printf "  [ok]   python:%-9s\n" "$mod"
    else
        printf "  [MISS] python:%-9s (install: pip install %s)\n" "$mod" "$pip_pkg"
        missing_pip+=("$pip_pkg")
        fatal=1
    fi
}

echo "AMSS-NCKU environment preflight"
echo "==============================="
echo "Configured COMPILER=$COMPILER"
echo
echo "Build tools:"
check_cmd make           build-essential
check_cmd python3        python3

if [[ "$COMPILER" == "gnu" ]]; then
    check_cmd gcc        build-essential
    check_cmd g++        build-essential
    check_cmd gfortran   gfortran
    check_cmd mpicxx     libopenmpi-dev   "or MPICH: libmpich-dev"
    check_cmd mpifort    libopenmpi-dev   "or MPICH: libmpich-dev"
elif [[ "$COMPILER" == "intel" ]]; then
    check_cmd icx        "intel-oneapi-compiler-dpcpp-cpp" \
        "see https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html"
    check_cmd ifx        "intel-oneapi-compiler-fortran"
    check_cmd mpiicpx    "intel-oneapi-mpi-devel"
else
    echo "  [WARN] Unknown COMPILER='$COMPILER' (expected 'gnu' or 'intel')."
    fatal=1
fi

echo
echo "MPI launcher:"
if command -v mpirun >/dev/null 2>&1; then
    printf "  [ok]   %-16s -> %s\n" "mpirun" "$(command -v mpirun)"
    mpirun --version 2>&1 | head -1 | sed 's/^/         /'
else
    printf "  [MISS] %-16s (install: libopenmpi-dev or libmpich-dev)\n" "mpirun"
    missing_apt+=("libopenmpi-dev")
    fatal=1
fi

echo
echo "Python packages:"
check_python_module numpy
check_python_module scipy
check_python_module matplotlib
check_python_module sympy

echo
if (( fatal == 0 )); then
    echo "All required tools and packages are present."
    echo "Next: ./run.sh"
    exit 0
fi

echo "Some prerequisites are missing. Suggested install commands:"
if (( ${#missing_apt[@]} )); then
    # Unique-ify apt list while preserving order.
    apt_uniq=$(printf '%s\n' "${missing_apt[@]}" | awk '!seen[$0]++' | tr '\n' ' ')
    echo "  sudo apt-get install ${apt_uniq}"
fi
if (( ${#missing_pip[@]} )); then
    pip_uniq=$(printf '%s\n' "${missing_pip[@]}" | awk '!seen[$0]++' | tr '\n' ' ')
    echo "  pip install ${pip_uniq}"
fi
exit 1
