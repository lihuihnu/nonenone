#!/usr/bin/env python3
"""Audit the experiment-driven five-component 2-D flow screening."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

SCOPE = "CONDITIONAL_5C_FLOW_SCREENING_NOT_FORMAL_VALIDATION"
COMPONENTS = (
    "H2O",
    "OIL_GASOLINE",
    "OIL_DIESEL",
    "OIL_MIDDLE",
    "OIL_HEAVY",
)
TARGETS = (0.25, 0.50, 0.75, 1.00)
MASS_LIMIT = 1e-6


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"empty table: {path}")
    return rows


def num(row, key):
    x = float(row[key])
    if not math.isfinite(x):
        raise ValueError(f"nonfinite {key}={x}")
    return x


def bracket(rows, target):
    p = [num(r, "pvi") for r in rows]
    if any(b <= a for a, b in zip(p, p[1:])):
        raise ValueError("PVI must increase strictly")
    if target < p[0] or target > p[-1]:
        raise ValueError(f"target {target} outside [{p[0]}, {p[-1]}]")
    for i, x in enumerate(p):
        if x == target:
            return rows[i], rows[i], 0.0
        if x > target:
            f = (target - p[i - 1]) / (x - p[i - 1])
            return rows[i - 1], rows[i], f
    raise AssertionError("unreachable")


def interp(rows, target, key):
    a, b, f = bracket(rows, target)
    x, y = float(a[key]), float(b[key])
    if f == 0.0:
        return x
    if not (math.isfinite(x) and math.isfinite(y)):
        return math.nan
    return x + f * (y - x)


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def audit(root: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    prod = read_csv(root / "producer_composition.csv")
    mb = read_csv(root / "component_mass_balance.csv")

    required = {
        "RF_total_hydrocarbon",
        "RF_Light_group",
        "RF_Heavy_group",
        "E_L_over_H_instant_mass",
        "E_L_over_H_cumulative_mass",
    }
    for comp in COMPONENTS:
        required.update(
            {
                f"m_dot_{comp}_produced_kg_s",
                f"Y_{comp}_mass_fraction",
                f"cum_{comp}_produced_kg",
                f"RF_{comp}",
            }
        )
    missing = sorted(required - set(prod[0]))
    if missing:
        raise ValueError(f"producer output contract missing columns: {missing}")

    final_pvi = num(prod[-1], "pvi")
    if final_pvi < 1.0:
        raise ValueError(f"trajectory stops before 1 PVI: {final_pvi}")

    mass_rows = []
    mass_gate = True
    present_components = {r["component"] for r in mb}
    if present_components != set(COMPONENTS):
        raise ValueError(
            f"mass ledger components {sorted(present_components)} "
            f"!= expected {sorted(COMPONENTS)}"
        )

    for comp in COMPONENTS:
        selected = [r for r in mb if r["component"] == comp]
        worst_rel = -1.0
        worst_abs = -1.0
        worst_step = -1
        for r in selected:
            initial = num(r, "initial_inventory_kg")
            injected = num(r, "cumulative_injected_kg")
            produced = num(r, "cumulative_produced_kg")
            current = num(r, "current_inventory_kg")
            error = current - initial - injected + produced
            budget = initial + injected
            rel = (
                abs(error) / budget
                if budget > 1e-24
                else (0.0 if abs(error) < 1e-24 else math.inf)
            )
            if rel > worst_rel:
                worst_rel = rel
                worst_abs = abs(error)
                worst_step = int(r["step"])
        mass_rows.append(
            {
                "scope": SCOPE,
                "component": comp,
                "max_relative_error": worst_rel,
                "max_absolute_error_kg": worst_abs,
                "worst_step": worst_step,
                "limit": MASS_LIMIT,
                "status": "PASS" if worst_rel <= MASS_LIMIT else "FAIL",
            }
        )
        mass_gate &= worst_rel <= MASS_LIMIT

    matched = []
    for target in TARGETS:
        row = {
            "scope": SCOPE,
            "target_pvi": target,
            "method": "BRACKETED_LINEAR_NO_EXTRAPOLATION",
        }
        for comp in COMPONENTS:
            for field in (
                f"RF_{comp}",
                f"Y_{comp}_mass_fraction",
                f"cum_{comp}_produced_kg",
            ):
                row[field] = interp(prod, target, field)
        for field in (
            "RF_Light_group",
            "RF_Heavy_group",
            "RF_total_hydrocarbon",
            "E_L_over_H_instant_mass",
            "E_L_over_H_cumulative_mass",
        ):
            row[field] = interp(prod, target, field)
        row["RF_Gasoline_minus_Heavy"] = (
            row["RF_OIL_GASOLINE"] - row["RF_OIL_HEAVY"]
        )
        row["RF_Diesel_minus_Heavy"] = (
            row["RF_OIL_DIESEL"] - row["RF_OIL_HEAVY"]
        )
        row["RF_Middle_minus_Heavy"] = (
            row["RF_OIL_MIDDLE"] - row["RF_OIL_HEAVY"]
        )
        matched.append(row)

    # Full producer trajectory with directly auditable selectivity diagnostics.
    trajectory = []
    for r in prod:
        rec = {
            "scope": SCOPE,
            "step": int(r["step"]),
            "time_s": num(r, "time_s"),
            "pvi": num(r, "pvi"),
            "RF_total_hydrocarbon": num(r, "RF_total_hydrocarbon"),
            "RF_Light_group": num(r, "RF_Light_group"),
            "RF_Heavy_group": num(r, "RF_Heavy_group"),
            "E_L_over_H_instant_mass": float(r["E_L_over_H_instant_mass"]),
            "E_L_over_H_cumulative_mass": float(r["E_L_over_H_cumulative_mass"]),
        }
        for comp in COMPONENTS:
            rec[f"RF_{comp}"] = float(r[f"RF_{comp}"])
            rec[f"Y_{comp}"] = float(r[f"Y_{comp}_mass_fraction"])
        trajectory.append(rec)

    write_csv(out / "five_component_mass_gate.csv", mass_rows)
    write_csv(out / "five_component_at_matched_pvi.csv", matched)
    write_csv(out / "five_component_producer_trajectory.csv", trajectory)

    status = {
        "scope": SCOPE,
        "final_pvi": final_pvi,
        "mass_limit": MASS_LIMIT,
        "mass_gate": "PASS" if mass_gate else "FAIL",
        "producer_output_contract": "PASS",
        "matched_pvi_coverage": "PASS",
        "physical_validation": "NOT_VALIDATED",
        "formal_360_380_pair": "NOT_STARTED",
    }
    (out / "status.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    return status


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    try:
        status = audit(a.root, a.out)
    except (OSError, KeyError, ValueError) as exc:
        print(f"FIVE_COMPONENT_SCREENING_AUDIT=FAIL: {exc}")
        raise SystemExit(2)
    print(json.dumps(status, indent=2))
    raise SystemExit(0 if status["mass_gate"] == "PASS" else 2)


if __name__ == "__main__":
    main()
