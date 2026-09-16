#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${root}"

python3 scripts/check_markdown_links.py
python3 scripts/check_source_comments.py
python3 scripts/check_vscode_config.py

bad_generated="$({ git ls-files | grep -E '(^|/)(results?|logs?|bin|build|\.build|reference_output)(/|$)|\.out$|\.log$|\.zip$|\.tar$|\.pyc$' || true; })"
if [[ -n "${bad_generated}" ]]; then
    echo "[FAIL] generated/runtime artifacts are tracked:" >&2
    echo "${bad_generated}" >&2
    exit 1
fi
echo "[PASS] no generated/runtime artifacts tracked"

# MATLAB 源码只允许作为明确归属某个正式算例的可复现辅助脚本；生成的
# 绘图脚本和散落在其它位置的 .m 文件仍视为运行产物。
bad_matlab="$({ git ls-files '*.m' | grep -Ev '^case/[^/]+/scripts/[^/]+\.m$' || true; })"
if [[ -n "${bad_matlab}" ]]; then
    echo "[FAIL] MATLAB files must live under case/<name>/scripts/:" >&2
    echo "${bad_matlab}" >&2
    exit 1
fi
echo "[PASS] tracked MATLAB helpers are confined to case-specific scripts/"

if [[ ! -f config/local.example.mk ]]; then
    echo "[FAIL] missing tracked local machine profile example: config/local.example.mk" >&2
    exit 1
fi
if ! grep -q '^# MPMC_SCW_GENERATED_RUN_SCRIPT=1$' case/run.sh.in; then
    echo "[FAIL] case/run.sh.in is missing the managed-script marker" >&2
    exit 1
fi
if grep -q '@MPI_LAUNCHER@ --mpi=' case/run.sh.in; then
    echo "[FAIL] run template hard-codes a yhrun-only --mpi option" >&2
    exit 1
fi
if ! grep -q 'MPI_LAUNCHER_ARGS_RAW="@MPI_LAUNCHER_ARGS@"' case/run.sh.in; then
    echo "[FAIL] run template does not consume profile MPI_LAUNCHER_ARGS" >&2
    exit 1
fi
echo "[PASS] local/HPC managed-run profile contract"

mapfile -t root_md < <(find . -maxdepth 1 -type f -name '*.md' -printf '%f\n' | sort)
expected_root_md="AGENTS.md README.md START_HERE.md"
if [[ "${root_md[*]:-}" != "${expected_root_md}" ]]; then
    echo "[FAIL] root Markdown must be limited to agent guidance and entry points; found: ${root_md[*]:-none}" >&2
    exit 1
fi
echo "[PASS] root documentation is limited to agent guidance and entry points"

# 只审计 Git 已跟踪的文本文件。
# 使用 git grep 而不是 grep -R，避免把本地编译生成的 .o、run.sh、logs/results
# 等未跟踪运行产物误判为“仓库中散落的超算绝对路径”。
bad_paths="$({ git grep -n -I -F '/vol8/home/' -- \
    . \
    ':(exclude)config/hpc.mk' \
    ':(exclude).vscode/c_cpp_properties.json' \
    ':(exclude)scripts/repo_audit.sh' || true; })"
if [[ -n "${bad_paths}" ]]; then
    echo "[FAIL] unexpected cluster-specific path outside config/hpc.mk:" >&2
    echo "${bad_paths}" >&2
    exit 1
fi
echo "[PASS] cluster dependency paths are centralized (VS Code defaults are audited mirrors)"

if git ls-files | grep -Eq 'AGENT_WORKLOG|_V[0-9]+\.md$|_VALIDATION\.md$'; then
    echo "[FAIL] version-process/project validation docs escaped docs/ consolidation:" >&2
    git ls-files | grep -E 'AGENT_WORKLOG|_V[0-9]+\.md$|_VALIDATION\.md$' >&2
    exit 1
fi
echo "[PASS] no top-level version-process documentation clutter"

obsolete_docs="$({ find docs -maxdepth 1 -type f -name 'V[0-9]*.md' -print 2>/dev/null; find docs/archive -type f -print 2>/dev/null || true; })"
if [[ -n "${obsolete_docs}" ]]; then
    echo "[FAIL] obsolete version/archive documentation exists in the current production package:" >&2
    echo "${obsolete_docs}" >&2
    exit 1
fi
echo "[PASS] no obsolete version/archive documentation in production package"
