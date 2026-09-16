#!/usr/bin/env python3
"""Check newly generated Sun-2024 factorial grid/time convergence runs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


MASS_LIMIT = 1.0e-6
RECOVERY_LIMIT = 0.005
BREAKTHROUGH_LIMIT_PV = 0.02
BREAKTHROUGH_RATE_FRACTION = 0.01


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def scalar_metadata(directory: Path) -> dict[str, str]:
    rows = read_rows(directory / "experiment_metadata.csv")
    if len(rows) != 1:
        raise ValueError(f"expected one metadata row in {directory}")
    return rows[0]


def metrics(directory: Path) -> dict[str, float | str]:
    metadata = scalar_metadata(directory)
    rows = read_rows(directory / "component_mass_balance.csv")
    if not rows:
        raise ValueError(f"empty component mass balance in {directory}")

    by_component: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_component.setdefault(row["component"], []).append(row)
    if "nC16" not in by_component or "CO2" not in by_component:
        raise ValueError(f"CO2/nC16 rows are required in {directory}")

    nc16 = by_component["nC16"]
    co2 = by_component["CO2"]
    final_nc16 = nc16[-1]
    initial_nc16 = float(final_nc16["initial_inventory_kg"])
    recovery = float(final_nc16["cumulative_produced_kg"]) / initial_nc16
    def normalized_closure(row: dict[str, str]) -> float:
        scale = max(
            abs(float(row["initial_inventory_kg"])),
            abs(float(row["current_inventory_kg"])),
            abs(float(row["cumulative_injected_kg"])),
            abs(float(row["cumulative_produced_kg"])),
            1.0e-30,
        )
        return abs(float(row["balance_error_kg"])) / scale

    max_mass_error = max(normalized_closure(row) for row in rows)
    max_absolute_mass_error = max(
        abs(float(row["balance_error_kg"])) for row in rows
    )

    pore_volume = float(metadata["pore_volume_m3"])
    total_rate = float(metadata["reservoir_total_rate_m3_s"])
    target_pv = float(metadata["target_pv"])
    final_time_day = max(float(row["time_day"]) for row in rows)
    final_pv = final_time_day * 86400.0 * total_rate / pore_volume

    co2_injection_mass_rate = float(metadata["co2_mass_rate_kg_s"])
    breakthrough_pv = math.nan
    maximum_rate_ratio = 0.0
    if co2_injection_mass_rate > 0.0:
        for previous, row in zip(co2, co2[1:]):
            dt_seconds = (
                float(row["time_day"]) - float(previous["time_day"])
            ) * 86400.0
            if dt_seconds <= 0.0:
                continue
            produced_rate = max(
                0.0,
                (
                    float(row["cumulative_produced_kg"])
                    - float(previous["cumulative_produced_kg"])
                )
                / dt_seconds,
            )
            rate_ratio = produced_rate / co2_injection_mass_rate
            maximum_rate_ratio = max(maximum_rate_ratio, rate_ratio)
            if rate_ratio >= BREAKTHROUGH_RATE_FRACTION:
                breakthrough_pv = (
                    float(row["time_day"]) * 86400.0 * total_rate / pore_volume
                )
                break

    return {
        "directory": str(directory),
        "experiment_id": metadata["experiment_id"],
        "nx": float(metadata["nx"]),
        "dt_pv": float(metadata["dt_pv"]),
        "target_pv": target_pv,
        "final_pv": final_pv,
        "nc16_recovery": recovery,
        "co2_breakthrough_pv": breakthrough_pv,
        "max_co2_outlet_inlet_rate_ratio": maximum_rate_ratio,
        "max_component_relative_closure": max_mass_error,
        "max_component_absolute_error_kg": max_absolute_mass_error,
    }


def absolute_difference(a: float, b: float) -> float:
    if math.isnan(a) and math.isnan(b):
        return 0.0
    if math.isnan(a) or math.isnan(b):
        return math.inf
    return abs(a - b)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--medium", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records = [metrics(args.coarse), metrics(args.medium), metrics(args.fine)]
    medium, fine = records[1], records[2]
    recovery_delta = absolute_difference(
        float(medium["nc16_recovery"]), float(fine["nc16_recovery"])
    )
    breakthrough_delta = absolute_difference(
        float(medium["co2_breakthrough_pv"]),
        float(fine["co2_breakthrough_pv"]),
    )
    endpoint_ok = all(
        abs(float(record["final_pv"]) - float(record["target_pv"])) <= 1.0e-8
        for record in records
    )
    mass_ok = all(
        float(record["max_component_relative_closure"]) <= MASS_LIMIT
        for record in records
    )
    passed = (
        endpoint_ok
        and mass_ok
        and recovery_delta <= RECOVERY_LIMIT
        and breakthrough_delta <= BREAKTHROUGH_LIMIT_PV
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(records[0].keys()) + [
        "medium_fine_recovery_delta",
        "medium_fine_breakthrough_pv_delta",
        "gate_passed",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in records:
            row = dict(record)
            row["medium_fine_recovery_delta"] = recovery_delta
            row["medium_fine_breakthrough_pv_delta"] = breakthrough_delta
            row["gate_passed"] = passed
            writer.writerow(row)

    print(f"convergence gate: {'PASS' if passed else 'FAIL'}")
    print(f"medium/fine nC16 recovery delta: {recovery_delta:.8g}")
    print(f"medium/fine CO2 breakthrough delta [PV]: {breakthrough_delta:.8g}")
    print(
        "CO2 breakthrough definition: outlet/inlet component mass-rate ratio "
        f">= {BREAKTHROUGH_RATE_FRACTION:g}"
    )
    print(f"maximum component mass error limit: {MASS_LIMIT:.3g}")
    print(f"summary: {args.output}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
