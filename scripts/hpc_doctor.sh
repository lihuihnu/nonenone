#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
profile="${root}/config/hpc.mk"
local_profile="${root}/config/hpc.local.mk"
mode="${1:-auto}"

if [[ ! -f "${profile}" ]]; then
    echo "[FAIL] missing ${profile}"
    exit 1
fi

# Local/private values are loaded first; config/hpc.mk contains ?= defaults.
read_make_var() {
    local var="$1"
    make -s -f - <<MAKEEOF
-include ${local_profile}
include ${profile}
all:
	@printf '%s' '\$(${var})'
MAKEEOF
}

MPMC_PLATFORM="$(read_make_var MPMC_PLATFORM)"
PETSC_DIR="$(read_make_var PETSC_DIR)"
PETSC_ARCH="$(read_make_var PETSC_ARCH)"
FMT_DIR="$(read_make_var FMT_DIR)"
METIS_DIR="$(read_make_var METIS_DIR)"
VTK_INSTALL="$(read_make_var VTK_INSTALL)"
X11_LIB_DIR="$(read_make_var X11_LIB_DIR)"
LOCAL_LIB_DIR="$(read_make_var LOCAL_LIB_DIR)"
MPI_LAUNCHER="$(read_make_var MPI_LAUNCHER)"
MPI_LAUNCHER_ARGS="$(read_make_var MPI_LAUNCHER_ARGS)"
MODULES_ENABLED="$(read_make_var MODULES_ENABLED)"
PETSC_CUSTOM_HOOK_LINKAGE="$(read_make_var PETSC_CUSTOM_HOOK_LINKAGE)"
SBATCH="$(read_make_var SBATCH)"
DQ_MESH_DIR="$(read_make_var DQ_MESH_DIR)"

if [[ "${mode}" != "auto" && "${mode}" != "local" && "${mode}" != "hpc" && "${mode}" != "full" ]]; then
    echo "Usage: $0 [auto|local|hpc|full]" >&2
    exit 2
fi
check_platform="${MPMC_PLATFORM}"
[[ "${mode}" == "local" || "${mode}" == "hpc" ]] && check_platform="${mode}"

fail=0
warn=0
ok()   { echo "[ OK ] $*"; }
warnf(){ echo "[WARN] $*"; warn=$((warn + 1)); }
failf(){ echo "[FAIL] $*"; fail=$((fail + 1)); }

check_command_or_path() {
    local label="$1" value="$2"
    if [[ -z "${value}" ]]; then
        failf "${label}: not configured"
    elif [[ "${value}" == */* ]]; then
        [[ -x "${value}" ]] && ok "${label}: ${value}" || failf "${label}: not executable: ${value}"
    elif command -v "${value}" >/dev/null 2>&1; then
        ok "${label}: $(command -v "${value}")"
    else
        failf "${label}: command not found: ${value}"
    fi
}

check_required_path() {
    local label="$1" path="$2"
    [[ -e "${path}" ]] && ok "${label}: ${path}" || failf "${label}: missing ${path}"
}

check_optional_path() {
    local label="$1" path="$2"
    if [[ -z "${path}" ]]; then
        warnf "${label}: not configured (optional for normal case builds)"
    elif [[ -e "${path}" ]]; then
        ok "${label}: ${path}"
    else
        warnf "${label}: missing ${path}"
    fi
}

echo "============================================================"
echo "MPMC_SCW environment doctor"
echo "Configured platform : ${MPMC_PLATFORM}"
echo "Check mode          : ${mode}"
echo "PETSC_DIR / ARCH    : ${PETSC_DIR} / ${PETSC_ARCH}"
echo "MPI launcher        : ${MPI_LAUNCHER} ${MPI_LAUNCHER_ARGS}"
echo "PETSc hook linkage  : ${PETSC_CUSTOM_HOOK_LINKAGE}"
echo "Local override      : $([[ -f "${local_profile}" ]] && echo "${local_profile}" || echo none)"
echo "============================================================"

command -v g++ >/dev/null 2>&1 && ok "g++: $(command -v g++)" || failf "g++: command not found"
command -v python3 >/dev/null 2>&1 && ok "python3: $(command -v python3)" || failf "python3: command not found"

check_required_path "PETSc variables" "${PETSC_DIR}/lib/petsc/conf/variables"
check_required_path "PETSc rules" "${PETSC_DIR}/lib/petsc/conf/rules"
check_required_path "PETSc arch lib" "${PETSC_DIR}/${PETSC_ARCH}/lib"
check_required_path "fmt include" "${FMT_DIR}/include"
check_required_path "METIS include" "${METIS_DIR}/include"
check_required_path "METIS library" "${METIS_DIR}/install/lib"
check_command_or_path "MPI launcher" "${MPI_LAUNCHER}"

if [[ "${PETSC_CUSTOM_HOOK_LINKAGE}" != "cxx" && "${PETSC_CUSTOM_HOOK_LINKAGE}" != "both" ]]; then
    failf "PETSC_CUSTOM_HOOK_LINKAGE must be cxx or both, got '${PETSC_CUSTOM_HOOK_LINKAGE}'"
fi

# The project's historical custom PETSc line-search hooks may be referenced with
# either C or C++ linkage depending on PETSc's --with-clanguage setting.  Detect
# the actual unresolved ABI before a long case link and verify the selected bridge.
petsc_shared="${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so"
if [[ -f "${petsc_shared}" ]] && command -v nm >/dev/null 2>&1; then
    petsc_undefined="$(nm -D --undefined-only "${petsc_shared}" 2>/dev/null || true)"
    needs_c=0
    needs_cxx=0
    if grep -Eq '[[:space:]]U[[:space:]]+(updateState|updateSol)$' <<<"${petsc_undefined}"; then
        needs_c=1
    fi
    if grep -Eq '[[:space:]]U[[:space:]]+_Z[0-9]+update(State|Sol)' <<<"${petsc_undefined}"; then
        needs_cxx=1
    fi

    if (( needs_c == 1 )); then
        if [[ "${PETSC_CUSTOM_HOOK_LINKAGE}" == "both" ]]; then
            ok "PETSc custom hooks: bare C updateState/updateSol detected; C bridge enabled"
        else
            failf "libpetsc.so requires bare C updateState/updateSol; set PETSC_CUSTOM_HOOK_LINKAGE := both in config/hpc.local.mk"
        fi
    fi
    if (( needs_cxx == 1 )); then
        ok "PETSc custom hooks: C++ updateState/updateSol ABI detected"
    fi
    if (( needs_c == 0 && needs_cxx == 0 )); then
        ok "PETSc custom hooks: libpetsc.so has no unresolved updateState/updateSol ABI requirement"
    fi
else
    warnf "PETSc custom hook ABI could not be inspected (missing libpetsc.so or nm)"
fi

launcher_base="$(basename "${MPI_LAUNCHER}")"
if [[ "${launcher_base}" == "mpiexec" || "${launcher_base}" == "mpirun" ]]; then
    if [[ " ${MPI_LAUNCHER_ARGS} " == *" --mpi="* || " ${MPI_LAUNCHER_ARGS} " == *" --mpi "* ]]; then
        failf "MPI_LAUNCHER_ARGS='${MPI_LAUNCHER_ARGS}' is a yhrun-style option and is invalid for ${launcher_base}; set MPI_LAUNCHER_ARGS :="
    else
        ok "MPI launcher arguments are compatible with ${launcher_base}"
    fi
fi

if [[ "${check_platform}" == "local" ]]; then
    if [[ "${MODULES_ENABLED}" == "1" ]]; then
        warnf "MODULES_ENABLED=1 on local platform; set MODULES_ENABLED := 0 unless Environment Modules is intentionally installed"
    else
        ok "Environment Modules disabled for local execution"
    fi
    for value in "${PETSC_DIR}" "${FMT_DIR}" "${METIS_DIR}" "${DQ_MESH_DIR}"; do
        [[ "${value}" == /vol8/* ]] && failf "local profile still resolves to cluster path: ${value}"
    done
else
    if [[ "${MODULES_ENABLED}" == "1" ]]; then
        if command -v module >/dev/null 2>&1; then ok "module command available"; else failf "module command not available"; fi
    fi
    check_command_or_path "sbatch" "${SBATCH}"
fi

# VTK and X11 are only required by the full VTK integration gate, not normal
# reservoir-case compilation/runs.
if [[ "${mode}" == "full" ]]; then
    check_required_path "VTK install" "${VTK_INSTALL}"
    check_required_path "X11 library directory" "${X11_LIB_DIR}"
else
    check_optional_path "VTK install" "${VTK_INSTALL}"
    check_optional_path "X11 library directory" "${X11_LIB_DIR}"
fi
check_optional_path "local library directory" "${LOCAL_LIB_DIR}"
check_optional_path "DQ mesh directory" "${DQ_MESH_DIR}"

if (( fail > 0 )); then
    echo "[FAIL] environment check failed (${fail} error(s), ${warn} warning(s))."
    exit 1
fi
echo "[PASS] environment check passed (${warn} warning(s))."
