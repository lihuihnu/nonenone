#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LINK_RE = re.compile(r"\[[^\]]*\]\(([^)]+)\)")
problems: list[str] = []
checked = 0

for md in sorted(ROOT.rglob("*.md")):
    if ".git" in md.parts:
        continue
    text = md.read_text(encoding="utf-8")
    # Strip fenced and inline code before interpreting Markdown link syntax.
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    text = re.sub(r"`[^`]*`", "", text)
    for raw in LINK_RE.findall(text):
        target = raw.strip()
        if not target or target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        target = target.split("#", 1)[0]
        if not target:
            continue
        checked += 1
        resolved = (md.parent / target).resolve()
        try:
            resolved.relative_to(ROOT.resolve())
        except ValueError:
            problems.append(f"{md.relative_to(ROOT)} -> outside repo: {raw}")
            continue
        if not resolved.exists():
            problems.append(f"{md.relative_to(ROOT)} -> missing: {raw}")

if problems:
    for problem in problems:
        print(f"[FAIL] {problem}")
    print(f"[FAIL] broken local Markdown links: {len(problems)} / {checked}")
    sys.exit(1)

print(f"[PASS] local Markdown links: {checked} checked")
