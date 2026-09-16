#!/usr/bin/env bash
# 生成不夹带机器私有文件的可继续开发交接包。
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/package_internal.sh [--force] [OUTPUT.zip]

The package contains:
  source/             committed HEAD source snapshot
  repository.bundle   all Git branches, tags and history
  START_HERE.md       continuation guide
  HANDOFF_MANIFEST.txt
  SHA256SUMS

Uncommitted working-tree changes are intentionally excluded.
EOF
}

force=0
output=""
while (($# > 0)); do
    case "$1" in
        --force) force=1 ;;
        -h|--help) usage; exit 0 ;;
        -*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
        *)
            if [[ -n "${output}" ]]; then
                echo "Only one output archive may be specified." >&2
                exit 2
            fi
            output="$1"
            ;;
    esac
    shift
done

for command_name in git tar python3 sha256sum mktemp; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "[FAIL] required command not found: ${command_name}" >&2
        exit 1
    fi
done

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
repo_name="$(basename "${repo_root}")"
commit="$(git -C "${repo_root}" rev-parse --verify HEAD)"
short_commit="$(git -C "${repo_root}" rev-parse --short=12 HEAD)"
branch="$(git -C "${repo_root}" symbolic-ref --quiet --short HEAD || echo DETACHED)"
describe="$(git -C "${repo_root}" describe --tags --always --dirty)"
package_name="${repo_name}-handoff-${short_commit}"

if [[ -z "${output}" ]]; then
    output="${repo_root}/../${package_name}.zip"
elif [[ "${output}" != /* ]]; then
    output="$(pwd)/${output}"
fi

if [[ -e "${output}" && "${force}" -ne 1 ]]; then
    echo "[FAIL] output already exists: ${output}" >&2
    echo "       pass --force to replace this exact archive" >&2
    exit 1
fi

dirty_count="$(git -C "${repo_root}" status --porcelain | wc -l | tr -d ' ')"
if [[ "${dirty_count}" -gt 0 ]]; then
    echo "[INFO] working tree has ${dirty_count} uncommitted path(s); they are excluded." >&2
fi

temporary_root="$(mktemp -d)"
cleanup() {
    rm -rf -- "${temporary_root}"
}
trap cleanup EXIT

payload="${temporary_root}/${package_name}"
source_dir="${payload}/source"
mkdir -p "${source_dir}"

git -C "${repo_root}" archive --format=tar "${commit}" |
    tar -xf - -C "${source_dir}"
git -C "${repo_root}" bundle create "${payload}/repository.bundle" --all

if [[ ! -f "${source_dir}/docs/CONTINUE_DEVELOPMENT.md" ]]; then
    echo "[FAIL] committed continuation guide is missing." >&2
    exit 1
fi
cp "${source_dir}/docs/CONTINUE_DEVELOPMENT.md" "${payload}/START_HERE.md"

created_utc="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
{
    echo "Format-Version: 1"
    echo "Project: ${repo_name}"
    echo "Created-UTC: ${created_utc}"
    echo "Commit: ${commit}"
    echo "Branch: ${branch}"
    echo "Describe: ${describe}"
    echo "Source-Snapshot: committed HEAD"
    echo "Git-History: repository.bundle (--all)"
    echo "Uncommitted-Paths-Excluded: ${dirty_count}"
    echo "Third-Party-Dependencies-Bundled: no"
    echo "External-Grid-Data-Bundled: no"
} >"${payload}/HANDOFF_MANIFEST.txt"

(
    cd "${payload}"
    sha256sum repository.bundle HANDOFF_MANIFEST.txt START_HERE.md >SHA256SUMS
    find source -type f -print0 | LC_ALL=C sort -z |
        xargs -0 sha256sum >>SHA256SUMS
)

mkdir -p "$(dirname "${output}")"
if [[ -e "${output}" ]]; then
    rm -f -- "${output}"
fi
python3 - "${temporary_root}" "${package_name}" "${output}" <<'PY'
from pathlib import Path
import sys
import zipfile

temporary_root = Path(sys.argv[1])
package_name = sys.argv[2]
output = Path(sys.argv[3])
package_root = temporary_root / package_name

with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(package_root.rglob("*")):
        archive.write(path, path.relative_to(temporary_root))

with zipfile.ZipFile(output) as archive:
    bad_member = archive.testzip()
    if bad_member is not None:
        raise RuntimeError(f"ZIP integrity check failed at {bad_member}")
PY

archive_hash="$(sha256sum "${output}" | awk '{print $1}')"
archive_size="$(wc -c <"${output}" | tr -d ' ')"
echo "[PASS] handoff package: ${output}"
echo "       commit: ${commit}"
echo "       bytes:  ${archive_size}"
echo "       sha256: ${archive_hash}"
