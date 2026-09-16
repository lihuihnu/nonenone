#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_ROOT="$(cd "${CASE_DIR}/.." && pwd)"
CASE_NAME="three_eos_3d_compare"

if [[ ! -x "${CASE_DIR}/bin/${CASE_NAME}" ]]; then
    echo "[ERROR] binary missing. Compile first:"
    echo "  cd ${CASE_ROOT} && make ${CASE_NAME} -j"
    exit 2
fi
if [[ ! -x "${CASE_DIR}/run.sh" ]]; then
    echo "[ERROR] generated run.sh missing. Prepare first:"
    echo "  cd ${CASE_ROOT} && make prepare CASE=${CASE_NAME} TASKS=8"
    exit 2
fi

for eos in pr sw cpa; do
    echo "[SBATCH] ${CASE_NAME} EOS=${eos} -> results/${eos}"
    (
        cd "${CASE_DIR}"
        sbatch --export="ALL,EOS=${eos},RESULT_DIR=./results/${eos}" run.sh
    )
done
