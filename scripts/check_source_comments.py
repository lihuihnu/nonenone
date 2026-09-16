#!/usr/bin/env python3
"""检查 C/C++ 源码是否保留统一的中文文件级注释。"""
from __future__ import annotations

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOTS = (
    "AdaptiveTimeStepper",
    "ad",
    "case",
    "common",
    "grid",
    "indices",
    "models",
    "output",
    "test",
    "tools",
)
EXTENSIONS = {".hpp", ".h", ".inc", ".cpp", ".cc", ".cxx", ".C"}
CHINESE = re.compile(r"[\u4e00-\u9fff]")


def main() -> int:
    files = sorted(
        path
        for root in SOURCE_ROOTS
        for path in (ROOT / root).rglob("*")
        if path.is_file() and path.suffix in EXTENSIONS
    )
    failures: list[str] = []

    for path in files:
        head = "\n".join(path.read_text(encoding="utf-8").splitlines()[:24])
        rel = path.relative_to(ROOT)
        if "@file" not in head or "@brief" not in head:
            failures.append(f"{rel}: 缺少文件级 @file/@brief")
            continue
        brief_lines = [line for line in head.splitlines() if "@brief" in line]
        if not brief_lines or not CHINESE.search(brief_lines[0]):
            failures.append(f"{rel}: 文件级 @brief 不是中文说明")

        # 生产源码的公开 Doxygen 摘要也必须直接给出中文语义；标准术语、
        # 变量名和公式可保留英文。测试代码只强制文件级说明，避免给测试样板
        # 增加无意义的注释负担。
        if rel.parts[0] != "test":
            for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if "@brief" in line and not CHINESE.search(line):
                    failures.append(f"{rel}:{line_no}: @brief 缺少中文接口说明")

    if failures:
        print("[FAIL] 源码中文注释规范检查失败:", file=sys.stderr)
        for item in failures:
            print(f"  - {item}", file=sys.stderr)
        return 1

    print(f"[PASS] {len(files)} 个 C/C++ 源文件均包含统一中文文件级注释")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
