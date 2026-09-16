#!/usr/bin/env bash
set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${case_dir}/bin/h2o_co2_bsb_lumped_3d_lab"
manifest="${case_dir}/experiment_manifest.csv"
result_root="${1:-${case_dir}/results/formal_matrix}"

if [[ ! -x "${binary}" ]]; then
    echo "Missing executable: ${binary}" >&2
    echo "Build from the project root with: make -C case h2o_co2_bsb_lumped_3d_lab -j2" >&2
    exit 2
fi

mkdir -p "${result_root}"

while IFS=, read -r run_order run_id eos x_scw x_co2 y_h2o y_co2 rate target_pvi dt_pvi nx ny nz role; do
    if [[ "${run_order}" == "run_order" ]]; then
        continue
    fi
    run_dir="${result_root}/$(printf '%02d' "${run_order}")_${run_id}"
    mkdir -p "${run_dir}"
    echo "[RUN ${run_order}/5] ${run_id}: EOS=${eos}, x_SCW=${x_scw}, x_CO2=${x_co2}, grid=${nx}x${ny}x${nz}"
    "${binary}" \
        -eos "${eos}" \
        -x_h2o_feed "${x_scw}" \
        -nx "${nx}" -ny "${ny}" -nz "${nz}" \
        -target_pv "${target_pvi}" -dt_pv "${dt_pvi}" \
        -result_dir "${run_dir}" \
        -snes_max_it 50 -snes_atol 1e-8 -snes_rtol 1e-8 -snes_stol 0 \
        -ksp_atol 1e-12 -ksp_rtol 1e-8 \
        >"${run_dir}/run.log" 2>&1
done < "${manifest}"

python3 "${case_dir}/scripts/analyze_results.py" "${result_root}"
echo "Completed five-point comparison: ${result_root}"
