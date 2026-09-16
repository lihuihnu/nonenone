#!/usr/bin/env python3
"""Compare recomputed EOS results against preserved reference result sets."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EOS_ORDER = ("pr", "sw", "cpa")
PHYSICAL_FILES = (
    "solution_final.csv",
    "phase_state_final.csv",
    "reservoir_diagnostics.csv",
    "well_history.csv",
    "component_mass_balance.csv",
    "time_step_history.csv",
)
NON_PHYSICAL_COLUMN_PARTS = ("wall", "elapsed", "runtime")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--current", type=Path, default=ROOT / "results" / "recomputed_f405a18"
    )
    parser.add_argument(
        "--baseline", type=Path, default=ROOT / "results" / "full"
    )
    parser.add_argument(
        "--optimized-cpa",
        type=Path,
        default=ROOT / "results" / "optimization" / "full_optimized",
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        data = list(reader)
        return list(reader.fieldnames or ()), data


def finite_number(value: str) -> float | None:
    if value == "":
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def compare_file(reference: Path, current: Path) -> dict[str, object]:
    reference_fields, reference_rows = rows(reference)
    current_fields, current_rows = rows(current)
    if reference_fields != current_fields:
        raise ValueError(f"column mismatch: {reference} vs {current}")
    if len(reference_rows) != len(current_rows):
        raise ValueError(f"row-count mismatch: {reference} vs {current}")

    maximum_absolute = 0.0
    maximum_scaled = 0.0
    maximum_location = ""
    compared_numbers = 0
    string_mismatches = 0
    skipped_columns: list[str] = []
    for field in reference_fields:
        if any(part in field.lower() for part in NON_PHYSICAL_COLUMN_PARTS):
            skipped_columns.append(field)
            continue
        for index, (old_row, new_row) in enumerate(zip(reference_rows, current_rows)):
            old_text, new_text = old_row[field], new_row[field]
            old_number, new_number = finite_number(old_text), finite_number(new_text)
            if old_number is not None and new_number is not None:
                absolute = abs(new_number - old_number)
                scaled = absolute / max(abs(old_number), abs(new_number), 1.0e-12)
                compared_numbers += 1
                if absolute > maximum_absolute:
                    maximum_absolute = absolute
                    maximum_location = f"row {index + 2}, column {field}"
                maximum_scaled = max(maximum_scaled, scaled)
            elif old_text != new_text:
                string_mismatches += 1

    return {
        "rows": len(reference_rows),
        "numeric_values": compared_numbers,
        "max_abs_difference": maximum_absolute,
        "max_scaled_difference": maximum_scaled,
        "max_abs_location": maximum_location,
        "string_mismatches": string_mismatches,
        "excluded_nonphysical_columns": sorted(set(skipped_columns)),
    }


def summary(path: Path) -> dict[str, float | int]:
    _, data = rows(path / "simulation_summary.csv")
    row = data[0]
    fields = (
        "accepted_internal_steps",
        "rejected_internal_steps",
        "attempted_snes",
        "attempted_ksp",
        "total_snes_wall",
        "simulation_loop_wall",
    )
    return {field: float(row[field]) for field in fields}


def comparison(reference: Path, current: Path) -> dict[str, object]:
    return {
        name: compare_file(reference / name, current / name)
        for name in PHYSICAL_FILES
    }


def main() -> None:
    args = parse_args()
    report: dict[str, object] = {
        "current_root": str(args.current.resolve()),
        "baseline_root": str(args.baseline.resolve()),
        "physical_files": list(PHYSICAL_FILES),
        "timing_columns_excluded_from_physical_comparison": True,
        "current_summaries": {
            eos: summary(args.current / eos) for eos in EOS_ORDER
        },
        "baseline_to_current": {
            eos: comparison(args.baseline / eos, args.current / eos)
            for eos in EOS_ORDER
        },
        "optimized_cpa_to_current_cpa": comparison(
            args.optimized_cpa, args.current / "cpa"
        ),
        "optimized_cpa_summary": summary(args.optimized_cpa),
        "baseline_cpa_summary": summary(args.baseline / "cpa"),
    }
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
