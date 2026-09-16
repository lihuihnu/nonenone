#!/usr/bin/env python3
"""检查 VS Code IntelliSense 配置与项目/HPC 默认配置的一致性。"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / ".vscode" / "c_cpp_properties.json"
HPC = ROOT / "config" / "hpc.mk"

REQUIRED_PROJECT_INCLUDES = {
    "${workspaceFolder}/ad/include",
    "${workspaceFolder}/common/include",
    "${workspaceFolder}/grid/include",
    "${workspaceFolder}/indices/include",
    "${workspaceFolder}/models/include",
    "${workspaceFolder}/output/include",
    "${workspaceFolder}/AdaptiveTimeStepper/include",
    "${workspaceFolder}/case/include",
    "${workspaceFolder}/case/**",
    "${workspaceFolder}/test/include",
    "${workspaceFolder}/tools/include",
}


def make_var(name: str) -> str:
    proc = subprocess.run(
        ["make", "-s", "-f", "-"],
        cwd=ROOT,
        input=f"include {HPC}\nall:\n\t@printf '%s' '$({name})'\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return proc.stdout


def fail(message: str) -> None:
    print(f"[FAIL] {message}", file=sys.stderr)
    raise SystemExit(1)


if not CFG.is_file():
    fail(f"missing {CFG.relative_to(ROOT)}")

try:
    data = json.loads(CFG.read_text(encoding="utf-8"))
except json.JSONDecodeError as exc:
    fail(f"invalid JSON in {CFG.relative_to(ROOT)}: {exc}")

if data.get("version") != 4:
    fail("c_cpp_properties.json version must be 4")
configs = data.get("configurations")
if not isinstance(configs, list) or len(configs) != 1:
    fail("exactly one canonical VS Code configuration is required")
conf = configs[0]
if conf.get("cppStandard") != "c++17":
    fail("VS Code cppStandard must match the project C++17 baseline")
if "FMT_HEADER_ONLY" not in conf.get("defines", []):
    fail("VS Code defines must contain FMT_HEADER_ONLY")

env = data.get("env", {})
project_includes = set(env.get("mpmcProjectIncludes", []))
missing = sorted(REQUIRED_PROJECT_INCLUDES - project_includes)
if missing:
    fail("missing project include paths: " + ", ".join(missing))

petsc = make_var("PETSC_DIR")
arch = make_var("PETSC_ARCH")
fmt = make_var("FMT_DIR")
metis = make_var("METIS_DIR")
vtk = make_var("VTK_INSTALL")
expected_defaults = {
    f"{petsc}/include",
    f"{petsc}/{arch}/include",
    f"{fmt}/include",
    f"{metis}/include",
    f"{metis}/install/include",
    f"{vtk}/include/vtk-9.3",
}
defaults = set(env.get("mpmcHpcDefaultIncludes", []))
if defaults != expected_defaults:
    fail("VS Code default external include paths do not match config/hpc.mk")

expected_dynamic = {
    "${env:PETSC_DIR}/include",
    "${env:PETSC_DIR}/${env:PETSC_ARCH}/include",
    "${env:FMT_DIR}/include",
    "${env:METIS_DIR}/include",
    "${env:METIS_DIR}/install/include",
    "${env:VTK_INSTALL}/include/vtk-9.3",
}
if set(env.get("mpmcHpcEnvironmentIncludes", [])) != expected_dynamic:
    fail("VS Code environment include fallback is incomplete")

# The project config must be tracked; personal VS Code files must stay ignored.
tracked = subprocess.run(
    ["git", "ls-files", "--error-unmatch", ".vscode/c_cpp_properties.json"],
    cwd=ROOT,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
).returncode == 0
if not tracked:
    fail(".vscode/c_cpp_properties.json must be tracked by Git")

ignored_local = subprocess.run(
    ["git", "check-ignore", "-q", ".vscode/settings.json"], cwd=ROOT
).returncode == 0
if not ignored_local:
    fail("personal .vscode/settings.json must remain ignored")

print("[PASS] VS Code IntelliSense configuration is consistent")
