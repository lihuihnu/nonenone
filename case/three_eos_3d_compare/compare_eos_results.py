#!/usr/bin/env python3
"""Merge PR/SW/CPA result CSVs from the 3-D EOS comparison case.

No third-party Python packages are required. Run after all three jobs finish:
    python3 compare_eos_results.py
"""
from __future__ import annotations

import csv
from pathlib import Path

EOS_NAMES = ("pr", "sw", "cpa")
ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "comparison"

RESERVOIR_COLUMNS = [
    "step", "time", "pressure_avg", "s_o_avg", "s_g_avg", "s_w_avg",
    "rho_o_avg", "rho_g_avg", "rho_w_avg", "mu_o_avg", "mu_g_avg", "mu_w_avg",
    "oil_only_cells", "gas_only_cells", "water_only_cells", "oil_gas_cells",
    "oil_water_cells", "gas_water_cells", "oil_gas_water_cells",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing required result file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_long_reservoir() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "eos_reservoir_history.csv"
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["eos", *RESERVOIR_COLUMNS])
        writer.writeheader()
        for eos in EOS_NAMES:
            rows = read_csv(RESULTS / eos / "reservoir_diagnostics.csv")
            for row in rows:
                writer.writerow({"eos": eos, **{key: row.get(key, "") for key in RESERVOIR_COLUMNS}})
    print(target)


def write_long_wells() -> None:
    target = OUT / "eos_well_history.csv"
    first = read_csv(RESULTS / EOS_NAMES[0] / "well_history.csv")
    if not first:
        raise SystemExit("empty PR well_history.csv")
    columns = list(first[0].keys())
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["eos", *columns])
        writer.writeheader()
        for eos in EOS_NAMES:
            rows = read_csv(RESULTS / eos / "well_history.csv")
            for row in rows:
                writer.writerow({"eos": eos, **row})
    print(target)


def write_solver_summary() -> None:
    """Merge per-EOS solver summaries so nonlinear robustness is compared separately."""
    target = OUT / "eos_solver_summary.csv"
    first = read_csv(RESULTS / EOS_NAMES[0] / "simulation_summary.csv")
    if not first:
        raise SystemExit("empty PR simulation_summary.csv")
    columns = list(first[0].keys())
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["eos", *columns])
        writer.writeheader()
        for eos in EOS_NAMES:
            rows = read_csv(RESULTS / eos / "simulation_summary.csv")
            if len(rows) != 1:
                raise SystemExit(f"expected one simulation_summary row for {eos}")
            writer.writerow({"eos": eos, **rows[0]})
    print(target)


def write_final_summary() -> None:
    target = OUT / "eos_final_summary.csv"
    fields = [
        "eos", "time_day", "pressure_avg_bar",
        "s_o_avg", "s_g_avg", "s_w_avg",
        "rho_o_avg_kg_m3", "rho_g_avg_kg_m3", "rho_w_avg_kg_m3",
        "mu_o_avg_Pa_s", "mu_g_avg_Pa_s", "mu_w_avg_Pa_s",
        "oil_gas_water_cells",
        "producer_bhp_bar", "producer_q_total_surface_m3_s",
        "producer_q_oil_surface_m3_s", "producer_q_gas_surface_m3_s",
        "producer_q_water_surface_m3_s",
        "accepted_internal_steps", "rejected_internal_steps",
        "attempted_snes", "attempted_ksp",
        "min_dt_day", "max_dt_day", "simulation_loop_wall_s",
    ]
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for eos in EOS_NAMES:
            reservoir = read_csv(RESULTS / eos / "reservoir_diagnostics.csv")[-1]
            wells = read_csv(RESULTS / eos / "well_history.csv")
            solver_rows = read_csv(RESULTS / eos / "simulation_summary.csv")
            if len(solver_rows) != 1:
                raise SystemExit(f"expected one simulation_summary row for {eos}")
            solver = solver_rows[0]
            producer_rows = [row for row in wells if row.get("name") == "PROD"]
            if not producer_rows:
                raise SystemExit(f"PROD rows missing for {eos}")
            producer = producer_rows[-1]
            writer.writerow({
                "eos": eos,
                "time_day": reservoir["time"],
                # reservoir_diagnostics.csv stores pressure directly in bar.
                "pressure_avg_bar": reservoir["pressure_avg"],
                "s_o_avg": reservoir["s_o_avg"],
                "s_g_avg": reservoir["s_g_avg"],
                "s_w_avg": reservoir["s_w_avg"],
                "rho_o_avg_kg_m3": reservoir["rho_o_avg"],
                "rho_g_avg_kg_m3": reservoir["rho_g_avg"],
                "rho_w_avg_kg_m3": reservoir["rho_w_avg"],
                "mu_o_avg_Pa_s": reservoir["mu_o_avg"],
                "mu_g_avg_Pa_s": reservoir["mu_g_avg"],
                "mu_w_avg_Pa_s": reservoir["mu_w_avg"],
                "oil_gas_water_cells": reservoir["oil_gas_water_cells"],
                # well_history.csv stores BHP directly in bar.
                "producer_bhp_bar": producer["bhp"],
                "producer_q_total_surface_m3_s": producer["q_total_surface"],
                "producer_q_oil_surface_m3_s": producer["q_oil_surface"],
                "producer_q_gas_surface_m3_s": producer["q_gas_surface"],
                "producer_q_water_surface_m3_s": producer["q_water_surface"],
                "accepted_internal_steps": solver["accepted_internal_steps"],
                "rejected_internal_steps": solver["rejected_internal_steps"],
                "attempted_snes": solver["attempted_snes"],
                "attempted_ksp": solver["attempted_ksp"],
                "min_dt_day": solver["min_dt_day"],
                "max_dt_day": solver["max_dt_day"],
                "simulation_loop_wall_s": solver["simulation_loop_wall"],
            })
    print(target)


def main() -> None:
    write_long_reservoir()
    write_long_wells()
    write_solver_summary()
    write_final_summary()


if __name__ == "__main__":
    main()
