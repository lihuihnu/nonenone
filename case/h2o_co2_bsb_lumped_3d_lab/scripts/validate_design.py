#!/usr/bin/env python3
"""Validate the five-point CO2/SCW design without reading run results."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_csv(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    plan = json.loads((ROOT / "RUN_PLAN.json").read_text(encoding="utf-8"))
    require(plan["design_status"] == "pilot_complete_formal_blocked",
            "wrong design status")

    components = read_csv("component_lumping.csv")
    names = [row["component"] for row in components]
    require(names == plan["components"], "component order differs from RUN_PLAN")
    require(math.isclose(sum(float(row["ppt_mole_fraction"]) for row in components),
                         1.0, rel_tol=0.0, abs_tol=1e-12),
            "PPT grouped fractions do not sum to one")
    hydrocarbon_rows = components[2:]
    require(math.isclose(
        sum(float(row["ppt_hydrocarbon_normalized_fraction"])
            for row in hydrocarbon_rows),
        1.0, rel_tol=0.0, abs_tol=1e-12),
        "normalized hydrocarbon fractions do not sum to one")

    initial = read_csv("initial_composition.csv")
    require([row["component"] for row in initial] == names,
            "initial component order mismatch")
    require(math.isclose(sum(float(row["overall_mole_fraction"]) for row in initial),
                         1.0, rel_tol=0.0, abs_tol=2e-15),
            "initial composition does not sum to one")

    geometry = plan["geometry_m"]
    grid = plan["grid"]
    require([grid["nx"], grid["ny"], grid["nz"]] == [60, 20, 1],
            "grid is not the unchanged nC10 60x20x1 grid")
    cell_size = [geometry["lx"] / grid["nx"],
                 geometry["ly"] / grid["ny"],
                 geometry["lz"] / grid["nz"]]
    require(all(math.isclose(value, 5.0, rel_tol=0.0, abs_tol=1e-15)
                for value in cell_size), "cells are not 5 m cubes")
    pore_volume = (geometry["lx"] * geometry["ly"] * geometry["lz"] *
                   plan["rock"]["porosity"])
    require(math.isclose(pore_volume, 30000.0, rel_tol=0.0, abs_tol=1e-12),
            "pore volume mismatch")
    rate = plan["injection"]["total_reservoir_rate_m3_s"]
    target_days = (pore_volume * plan["time"]["target_pvi"] / rate / 86400.0)
    require(math.isclose(target_days, plan["time"]["target_time_day"],
                         rel_tol=0.0, abs_tol=1e-10),
            "target time is inconsistent with rate and PVI")

    runs = read_csv("experiment_manifest.csv")
    require(len(runs) == 5, "the five explicitly listed mixtures must be retained")
    require([int(row["run_order"]) for row in runs] == list(range(1, 6)),
            "run order must be consecutive")
    require([row["run_id"] for row in runs] == plan["run_order"],
            "run order differs from RUN_PLAN")
    require({row["eos"] for row in runs} == {plan["primary_thermodynamic_model"]},
            "EOS must stay fixed across the primary comparison")
    levels = {float(row["x_scw_feed"]) for row in runs}
    require(levels == {0.0, 0.25, 0.5, 0.75, 1.0},
            "SCW composition levels are incomplete or duplicated")

    molar_mass_h2o = float(components[0]["M_kg_mol"])
    molar_mass_co2 = float(components[1]["M_kg_mol"])
    for row in runs:
        x_scw = float(row["x_scw_feed"])
        x_co2 = float(row["x_co2_feed"])
        require(math.isclose(x_scw + x_co2, 1.0, rel_tol=0.0, abs_tol=1e-15),
                f"feed mole fractions do not sum in {row['run_id']}")
        denominator = x_scw * molar_mass_h2o + x_co2 * molar_mass_co2
        expected_y = x_scw * molar_mass_h2o / denominator
        require(math.isclose(float(row["y_h2o_feed"]), expected_y,
                             rel_tol=0.0, abs_tol=5e-16),
                f"mass fraction mismatch in {row['run_id']}")
        require(math.isclose(float(row["y_h2o_feed"]) +
                             float(row["y_co2_feed"]),
                             1.0, rel_tol=0.0, abs_tol=5e-16),
                f"feed mass fractions do not sum in {row['run_id']}")
        require(math.isclose(float(row["total_reservoir_rate_m3_s"]), rate,
                             rel_tol=0.0, abs_tol=1e-20),
                f"rate mismatch in {row['run_id']}")
        require([int(row["grid_nx"]), int(row["grid_ny"]), int(row["grid_nz"])] ==
                [60, 20, 1], f"grid drift in {row['run_id']}")

    print("PASS: five-point CO2/SCW displacement design is internally consistent")
    print(f"components={len(names)} runs={len(runs)} cells={grid['nx']*grid['ny']*grid['nz']}")
    print(f"pore_volume_m3={pore_volume:.12g} target_days={target_days:.12g}")


if __name__ == "__main__":
    main()
