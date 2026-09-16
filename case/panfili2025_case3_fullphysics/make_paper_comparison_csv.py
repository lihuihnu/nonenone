#!/usr/bin/env python3
"""Build a compact Panfili-2025 Fig.25/26 comparison table from MPMC CSV output.

Reads the standard MPMC output files in ./results and writes
./results/panfili_fig25_26_comparison.csv with one row per fixed output step.
No third-party Python packages are required.
"""

from __future__ import annotations

import csv
from pathlib import Path

SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.25
BAR_TO_PSI = 14.503773773
KG_PER_MTON = 1.0e9


def rows_by_step(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            yield int(row["step"]), row


def main() -> None:
    results = Path(__file__).resolve().parent / "results"
    balance_path = results / "component_mass_balance.csv"
    reservoir_path = results / "reservoir_diagnostics.csv"
    wells_path = results / "well_history.csv"
    for path in (balance_path, reservoir_path, wells_path):
        if not path.exists():
            raise SystemExit(f"missing required result file: {path}")

    pressure = {}
    time_days = {}
    for step, row in rows_by_step(reservoir_path):
        time_days[step] = float(row["time"])
        pressure[step] = float(row["pressure_avg"])

    co2_injected = {}
    for step, row in rows_by_step(balance_path):
        if row["component"] == "CO2":
            co2_injected[step] = float(row["cumulative_injected_kg"])

    water_prod_m3_s = {}
    with wells_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            step = int(row["step"])
            q = float(row["q_water_surface"])
            # Simulator sign: injection +, production -.  Fig.26 reports a
            # positive field-production magnitude.
            if q < 0.0:
                water_prod_m3_s[step] = water_prod_m3_s.get(step, 0.0) - q

    steps = sorted(set(pressure) & set(co2_injected))
    output = results / "panfili_fig25_26_comparison.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "step",
            "time_year",
            "cumulative_CO2_injected_kg",
            "cumulative_CO2_injected_Mton",
            "average_reservoir_pressure_bar",
            "average_reservoir_pressure_psi",
            "field_water_production_m3_s",
            "field_water_production_m3_day",
        ])
        for step in steps:
            days = time_days[step]
            injected = co2_injected[step]
            pw = pressure[step]
            qw = water_prod_m3_s.get(step, 0.0)
            writer.writerow([
                step,
                f"{days / DAYS_PER_YEAR:.12g}",
                f"{injected:.16g}",
                f"{injected / KG_PER_MTON:.16g}",
                f"{pw:.16g}",
                f"{pw * BAR_TO_PSI:.16g}",
                f"{qw:.16g}",
                f"{qw * SECONDS_PER_DAY:.16g}",
            ])

    print(output)


if __name__ == "__main__":
    main()
