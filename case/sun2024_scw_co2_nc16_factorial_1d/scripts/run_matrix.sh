#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
case_name="sun2024_scw_co2_nc16_factorial_1d"
mode="${1:-formal}"
np="${NP:-1}"
target_pv="${TARGET_PV:-4}"
formal_nx="${FORMAL_NX:-96}"
formal_dt_pv="${FORMAL_DT_PV:-0.01}"
rate_basis="${RATE_BASIS:-}"

if [[ -z "${rate_basis}" ]]; then
    echo "RATE_BASIS is required: in_situ_volume or reference_density" >&2
    exit 2
fi

rate_args="-rate_basis ${rate_basis}"
case "${rate_basis}" in
    in_situ_volume)
        ;;
    reference_density)
        : "${WATER_REFERENCE_DENSITY_KG_M3:?required for reference_density}"
        : "${CO2_REFERENCE_DENSITY_KG_M3:?required for reference_density}"
        rate_args+=" -water_reference_density_kg_m3 ${WATER_REFERENCE_DENSITY_KG_M3}"
        rate_args+=" -co2_reference_density_kg_m3 ${CO2_REFERENCE_DENSITY_KG_M3}"
        ;;
    *)
        echo "RATE_BASIS must be in_situ_volume or reference_density" >&2
        exit 2
        ;;
esac

base="${RESULT_BASE:-${repo_root}/results/sun2024-factorial/${rate_basis}}"

make -C "${repo_root}" case CASE="${case_name}" -j

run_case() {
    local id="$1"
    local nx="$2"
    local dt_pv="$3"
    local destination="$4"
    mkdir -p "${destination}"
    git -C "${repo_root}" rev-parse HEAD > "${destination}/input_git_commit.txt"
    printf '%q ' "$0" "$@" > "${destination}/driver_invocation.txt"
    printf '\n' >> "${destination}/driver_invocation.txt"
    {
        printf 'RATE_BASIS=%s\n' "${rate_basis}"
        printf 'WATER_REFERENCE_DENSITY_KG_M3=%s\n' "${WATER_REFERENCE_DENSITY_KG_M3:-not_applicable}"
        printf 'CO2_REFERENCE_DENSITY_KG_M3=%s\n' "${CO2_REFERENCE_DENSITY_KG_M3:-not_applicable}"
        printf 'TARGET_PV=%s\n' "${target_pv}"
        printf 'FORMAL_NX=%s\n' "${formal_nx}"
        printf 'FORMAL_DT_PV=%s\n' "${formal_dt_pv}"
        printf 'NP=%s\n' "${np}"
    } > "${destination}/driver_environment.txt"
    make -C "${repo_root}" run CASE="${case_name}" NP="${np}" \
        RESULT_DIR="${destination}" \
        RUN_ARGS="-experiment_id ${id} ${rate_args} -nx ${nx} -target_pv ${target_pv} -dt_pv ${dt_pv}"
}

case "${mode}" in
    formal)
        for id in R04 R06 R10 R12 C23 C24; do
            run_case "${id}" "${formal_nx}" "${formal_dt_pv}" \
                "${base}/formal/${id}"
        done
        ;;
    grid)
        run_case R12 48 0.005 "${base}/grid/R12_nx48"
        run_case R12 96 0.005 "${base}/grid/R12_nx96"
        run_case R12 192 0.005 "${base}/grid/R12_nx192"
        ;;
    time)
        run_case R12 192 0.02 "${base}/time/R12_dtpv0p02"
        run_case R12 192 0.01 "${base}/time/R12_dtpv0p01"
        run_case R12 192 0.005 "${base}/time/R12_dtpv0p005"
        ;;
    smoke)
        target_pv=0.0001
        run_case R12 48 0.0001 "${base}/smoke/R12"
        ;;
    *)
        echo "usage: $0 [formal|grid|time|smoke]" >&2
        exit 2
        ;;
esac
