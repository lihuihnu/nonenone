#!/usr/bin/env bash
set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${case_dir}/bin/h2o_co2_bsb_lumped_3d_lab"
manifest="${case_dir}/experiment_manifest.csv"
result_root="${1:-${case_dir}/results/pilot_matrix}"
pilot_pvi="${2:-0.0001}"
pilot_dt_pvi="${3:-${pilot_pvi}}"

mkdir -p "${result_root}"
while IFS=, read -r run_order run_id eos x_scw x_co2 y_h2o y_co2 rate target_pvi dt_pvi nx ny nz role; do
    [[ "${run_order}" == "run_order" ]] && continue
    run_dir="${result_root}/$(printf '%02d' "${run_order}")_${run_id}"
    mkdir -p "${run_dir}"
    echo "[PILOT ${run_order}/5] ${run_id}"
    "${binary}" \
        -eos "${eos}" -x_h2o_feed "${x_scw}" \
        -target_pv "${pilot_pvi}" -dt_pv "${pilot_dt_pvi}" \
        -result_dir "${run_dir}" \
        -snes_max_it 50 -snes_atol 1e-8 -snes_rtol 1e-8 -snes_stol 0 \
        -ksp_atol 1e-12 -ksp_rtol 1e-8 \
        >"${run_dir}/run.log" 2>&1
done < "${manifest}"

echo "Completed five-run pilot at ${pilot_pvi} PVI: ${result_root}"
