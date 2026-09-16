#!/bin/bash
#SBATCH -p mt_module
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --job-name=MPMC_SCW_test
#SBATCH --err=%j_test.err
#SBATCH --output=%j_test.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${SCRIPT_DIR}"
mkdir -p logs results

PROFILE="${PROJECT_ROOT}/config/hpc.mk"
LOCAL_PROFILE="${PROJECT_ROOT}/config/hpc.local.mk"
read_make_var() {
    local var="$1"
    make -s -f - <<MAKEEOF
-include ${LOCAL_PROFILE}
include ${PROFILE}
all:
	@printf '%s' '\$(${var})'
MAKEEOF
}

PETSC_DIR="${PETSC_DIR:-$(read_make_var PETSC_DIR)}"
PETSC_ARCH="${PETSC_ARCH:-$(read_make_var PETSC_ARCH)}"
VTK_INSTALL="${VTK_INSTALL:-$(read_make_var VTK_INSTALL)}"
METIS_DIR="${METIS_DIR:-$(read_make_var METIS_DIR)}"
LOCAL_LIB_DIR="${LOCAL_LIB_DIR:-$(read_make_var LOCAL_LIB_DIR)}"
X11_LIB_DIR="${X11_LIB_DIR:-$(read_make_var X11_LIB_DIR)}"
MPI_MODULE="${MPI_MODULE:-$(read_make_var MPI_MODULE)}"
BLAS_MODULE="${BLAS_MODULE:-$(read_make_var BLAS_MODULE)}"
MODULES_ENABLED="${MODULES_ENABLED:-$(read_make_var MODULES_ENABLED)}"
MESH_DIR="${MESH_DIR:-$(read_make_var DQ_MESH_DIR)}"

if [[ "${MODULES_ENABLED}" == "1" ]]; then
    if ! command -v module >/dev/null 2>&1; then
        echo "[ERROR] module command is required by this profile but is unavailable" >&2
        exit 1
    fi
    [[ -n "${MPI_MODULE}" ]] && module load "${MPI_MODULE}"
    [[ -n "${BLAS_MODULE}" ]] && module load "${BLAS_MODULE}"
fi

export PETSC_DIR PETSC_ARCH
PETSC_LIB_DIR="${PETSC_DIR}/${PETSC_ARCH}/lib"

if [[ ! -e "${X11_LIB_DIR}/libX11.so.6" ]]; then
    echo "[ERROR] libX11.so.6 was not found in: ${X11_LIB_DIR}" >&2
    echo "        Override X11_LIB_DIR or config/hpc.local.mk." >&2
    exit 1
fi

export LD_LIBRARY_PATH="${LOCAL_LIB_DIR}:${PETSC_LIB_DIR}:${VTK_INSTALL}/lib:${METIS_DIR}/install/lib:${X11_LIB_DIR}:${LD_LIBRARY_PATH:-}"

NP="${NP:-${SLURM_NTASKS:-1}}"
PARTITION="${SLURM_JOB_PARTITION:-mt_module}"

echo "============================================================"
echo "MPMC_SCW unified test suite"
echo "PETSC_DIR      = ${PETSC_DIR}"
echo "PETSC_ARCH     = ${PETSC_ARCH}"
echo "MESH_DIR       = ${MESH_DIR}"
echo "Partition      = ${PARTITION}"
echo "MPI ranks      = ${NP}"
echo "Start          = $(date)"
echo "============================================================"

make clean
make run-full NP="${NP}" MESH_DIR="${MESH_DIR}" PARTITION="${PARTITION}"

echo "============================================================"
echo "All tests finished successfully."
echo "End = $(date)"
echo "============================================================"
