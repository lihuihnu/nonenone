#!/usr/bin/env python3
"""Compare the final SCW/CO2 co-injection and CO2-only result directories."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

NX, NY, NZ = 24, 16, 6


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def mean(rows: list[dict[str, str]], column: str) -> float:
    return sum(float(row[column]) for row in rows) / len(rows)


def saturation_column(rows: list[dict[str, str]], phase: str) -> str:
    candidates = {
        "oil": ("liquid_saturation", "oil_saturation"),
        "gas": ("vapor_saturation", "gas_saturation"),
        "water": ("water_saturation",),
    }[phase]
    for candidate in candidates:
        if candidate in rows[0]:
            return candidate
    raise KeyError(f"Cannot find {phase} saturation column")


def component_recovery(result_dir: Path, component: str) -> float:
    path = result_dir / "component_mass_balance.csv"
    rows = [row for row in read_rows(path) if row["component"] == component]
    if not rows:
        raise KeyError(f"No {component} rows in {path}")
    last = max(rows, key=lambda row: (float(row["time_day"]), int(row["step"])))
    initial = float(last["initial_inventory_kg"])
    return float(last["cumulative_produced_kg"]) / initial if initial > 0.0 else math.nan


def final_metrics(result_dir: Path) -> dict[str, float]:
    rows = read_rows(result_dir / "solution_final.csv")
    if len(rows) != NX * NY * NZ:
        raise ValueError(
            f"{result_dir}: expected {NX * NY * NZ} final cells, got {len(rows)}"
        )

    oil = saturation_column(rows, "oil")
    gas = saturation_column(rows, "gas")
    water = saturation_column(rows, "water")
    pressures = [float(row["pressure_Pa"]) for row in rows]
    gas_values = [float(row[gas]) for row in rows]
    water_values = [float(row[water]) for row in rows]

    gas_total = sum(gas_values)
    water_total = sum(water_values)
    gas_x = gas_y = gas_z = 0.0
    top_water = top_gas = 0.0
    plume_cells = 0
    for row, sg, sw in zip(rows, gas_values, water_values):
        cell = int(row["input_index"])
        i = cell % NX
        j = (cell // NX) % NY
        k = cell // (NX * NY)
        gas_x += sg * (i + 0.5) / NX
        gas_y += sg * (j + 0.5) / NY
        gas_z += sg * (k + 0.5) / NZ
        if k >= NZ // 2:
            top_water += sw
            top_gas += sg
        if sg > 1.0e-4:
            plume_cells += 1

    return {
        "pressure_mean_bar": sum(pressures) / len(pressures) / 1.0e5,
        "pressure_span_bar": (max(pressures) - min(pressures)) / 1.0e5,
        "oil_saturation_mean": mean(rows, oil),
        "gas_saturation_mean": mean(rows, gas),
        "water_saturation_mean": mean(rows, water),
        "gas_plume_cell_fraction": plume_cells / len(rows),
        "gas_centroid_x_fraction": gas_x / gas_total if gas_total > 0.0 else math.nan,
        "gas_centroid_y_fraction": gas_y / gas_total if gas_total > 0.0 else math.nan,
        "gas_centroid_z_fraction": gas_z / gas_total if gas_total > 0.0 else math.nan,
        "top_half_gas_share": top_gas / gas_total if gas_total > 0.0 else math.nan,
        "top_half_water_share": top_water / water_total if water_total > 0.0 else math.nan,
        "nC4_recovery_fraction": component_recovery(result_dir, "nC4"),
        "nC16_recovery_fraction": component_recovery(result_dir, "nC16"),
    }


def format_value(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.8g}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scw", type=Path, required=True,
                        help="SCW/CO2 co-injection result directory")
    parser.add_argument("--control", type=Path, required=True,
                        help="CO2-only result directory")
    parser.add_argument("--output", type=Path, default=Path("scw_effect_summary.csv"))
    args = parser.parse_args()

    scw = final_metrics(args.scw)
    control = final_metrics(args.control)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "scw_co2", "co2_only", "difference"))
        for metric in scw:
            writer.writerow((metric, format_value(scw[metric]),
                             format_value(control[metric]),
                             format_value(scw[metric] - control[metric])))

    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
