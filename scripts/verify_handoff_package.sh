#!/usr/bin/env bash
# 独立校验交接包校验和、源码快照与 Git 历史的一致性。
set -euo pipefail

if [[ "$#" -ne 1 || "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: scripts/verify_handoff_package.sh PACKAGE.zip"
    exit "$([[ "$#" -eq 1 ]] && echo 0 || echo 2)"
fi

for command_name in git tar python3 sha256sum diff mktemp; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "[FAIL] required command not found: ${command_name}" >&2
        exit 1
    fi
done

archive="$1"
if [[ ! -f "${archive}" ]]; then
    echo "[FAIL] package does not exist: ${archive}" >&2
    exit 1
fi

temporary_root="$(mktemp -d)"
cleanup() {
    rm -rf -- "${temporary_root}"
}
trap cleanup EXIT

python3 - "${archive}" "${temporary_root}/unpacked" <<'PY'
from pathlib import Path
import os
import sys
import zipfile

archive_path = Path(sys.argv[1])
destination = Path(sys.argv[2]).resolve()
destination.mkdir(parents=True, exist_ok=True)

with zipfile.ZipFile(archive_path) as archive:
    bad_member = archive.testzip()
    if bad_member is not None:
        raise RuntimeError(f"ZIP integrity check failed at {bad_member}")
    for member in archive.infolist():
        target = (destination / member.filename).resolve()
        if os.path.commonpath((str(destination), str(target))) != str(destination):
            raise RuntimeError(f"unsafe ZIP path: {member.filename}")
    archive.extractall(destination)
PY
mapfile -t roots < <(find "${temporary_root}/unpacked" -mindepth 1 -maxdepth 1 -type d)
if [[ "${#roots[@]}" -ne 1 ]]; then
    echo "[FAIL] package must contain exactly one top-level directory." >&2
    exit 1
fi
payload="${roots[0]}"

for required in source repository.bundle START_HERE.md HANDOFF_MANIFEST.txt SHA256SUMS; do
    if [[ ! -e "${payload}/${required}" ]]; then
        echo "[FAIL] package entry is missing: ${required}" >&2
        exit 1
    fi
done

(
    cd "${payload}"
    sha256sum -c SHA256SUMS >/dev/null
)

commit="$(awk -F': ' '$1 == "Commit" {print $2}' "${payload}/HANDOFF_MANIFEST.txt")"
if [[ ! "${commit}" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "[FAIL] manifest commit is invalid: ${commit}" >&2
    exit 1
fi

git clone -q "${payload}/repository.bundle" "${temporary_root}/clone"
git -C "${temporary_root}/clone" cat-file -e "${commit}^{commit}"
mkdir -p "${temporary_root}/expected"
git -C "${temporary_root}/clone" archive --format=tar "${commit}" |
    tar -xf - -C "${temporary_root}/expected"
if ! diff -qr "${payload}/source" "${temporary_root}/expected" >/dev/null; then
    echo "[FAIL] source snapshot does not match the manifest commit." >&2
    exit 1
fi

echo "[PASS] handoff package verified"
echo "       archive: ${archive}"
echo "       commit:  ${commit}"
