#!/usr/bin/env python3
"""Build the frozen common PR/SW/CPA calibration table from NIST ThermoML."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


SEED = 20260825
COMPONENTS = ("h2o", "co2", "ch4", "c2", "nc4")
WATER_MOLAR_MASS = 0.018015268


def variables(point: dict[str, Any]) -> dict[int, float]:
    return {
        int(item["nVarNumber"]): float(item["nVarValue"])
        for item in point["VariableValue"]
    }


def property_value(point: dict[str, Any]) -> tuple[float, float]:
    item = point["PropertyValue"][0]
    uncertainty = item.get("CombinedUncertainty", {}).get(
        "nCombExpandUncertValue", math.nan
    )
    return float(item["nPropValue"]), float(uncertainty)


def row_template(identifier: str, source: str, pair: str, kind: str) -> dict[str, Any]:
    row: dict[str, Any] = {
        "id": identifier,
        "source": source,
        "pair": pair,
        "kind": kind,
        "split": "",
        "temperature_K": 0.0,
        "pressure_Pa": 0.0,
        "salinity_molal": 0.0,
        "density_kg_m3": 0.0,
        "uncertainty": 0.0,
    }
    for phase in ("x", "y"):
        for component in COMPONENTS:
            row[f"{phase}_{component}"] = 0.0
    return row


def pair_by_variables(
    first: dict[str, Any], second: dict[str, Any]
) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    lookup = {
        tuple(round(value, 10) for value in variables(point).values()): point
        for point in second["NumValues"]
    }
    result = []
    for point in first["NumValues"]:
        key = tuple(round(value, 10) for value in variables(point).values())
        if key in lookup:
            result.append((point, lookup[key]))
    return result


def load_rows(reference_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    may = json.loads(
        (reference_dir / "may2015_pr_vle.json").read_text(encoding="utf-8-sig")
    )
    for pair, pressure_index, vapor_index, second_component in (
        ("ch4_c2", 0, 1, "c2"),
        ("ch4_nc4", 6, 7, "nc4"),
    ):
        for index, (pressure_point, vapor_point) in enumerate(
            pair_by_variables(
                may["PureOrMixtureData"][pressure_index],
                may["PureOrMixtureData"][vapor_index],
            )
        ):
            values = variables(pressure_point)
            pressure_kpa, pressure_u_kpa = property_value(pressure_point)
            y_second, y_u = property_value(vapor_point)
            x_ch4 = values[2]
            row = row_template(
                f"may2015_{pair}_{index:03d}", "May et al. (2015)", pair, "vle"
            )
            row.update(
                temperature_K=values[1],
                pressure_Pa=pressure_kpa * 1000.0,
                uncertainty=max(pressure_u_kpa * 1000.0, y_u),
                x_ch4=x_ch4,
                y_ch4=1.0 - y_second,
            )
            row[f"x_{second_component}"] = 1.0 - x_ch4
            row[f"y_{second_component}"] = y_second
            rows.append(row)

    petropoulou = json.loads(
        (reference_dir / "petropoulou2018_co2_ch4_vle.json").read_text(
            encoding="utf-8-sig"
        )
    )
    liquid, vapor = petropoulou["PureOrMixtureData"]
    # Liquid and gas capillary samples are stored in acquisition order.  Their
    # measured T/P values differ by up to 0.004 K and 0.4 kPa, so exact-key
    # matching would incorrectly discard valid tie lines.
    for index, (liquid_point, vapor_point) in enumerate(
        zip(liquid["NumValues"], vapor["NumValues"], strict=True)
    ):
        liquid_values = variables(liquid_point)
        vapor_values = variables(vapor_point)
        x_co2, x_u = property_value(liquid_point)
        y_co2, y_u = property_value(vapor_point)
        row = row_template(
            f"petropoulou2018_co2_ch4_{index:03d}",
            "Petropoulou et al. (2018)",
            "co2_ch4",
            "vle",
        )
        row.update(
            temperature_K=0.5 * (liquid_values[1] + vapor_values[1]),
            pressure_Pa=500.0 * (liquid_values[2] + vapor_values[2]),
            uncertainty=max(x_u, y_u),
            x_co2=x_co2,
            x_ch4=1.0 - x_co2,
            y_co2=y_co2,
            y_ch4=1.0 - y_co2,
        )
        rows.append(row)

    frost = json.loads(
        (reference_dir / "frost2014_methane_water_vle.json").read_text(
            encoding="utf-8-sig"
        )
    )
    pressure_data, vapor_data = frost["PureOrMixtureData"][:2]
    for index, (pressure_point, vapor_point) in enumerate(
        pair_by_variables(pressure_data, vapor_data)
    ):
        values = variables(pressure_point)
        pressure_kpa, pressure_u_kpa = property_value(pressure_point)
        y_h2o, y_u = property_value(vapor_point)
        x_ch4 = values[1]
        row = row_template(
            f"frost2014_h2o_ch4_{index:03d}",
            "Frost et al. (2014)",
            "h2o_ch4",
            "vle",
        )
        row.update(
            temperature_K=values[2],
            pressure_Pa=pressure_kpa * 1000.0,
            uncertainty=max(pressure_u_kpa * 1000.0, y_u),
            x_h2o=1.0 - x_ch4,
            x_ch4=x_ch4,
            y_h2o=y_h2o,
            y_ch4=1.0 - y_h2o,
        )
        rows.append(row)

    density = json.loads(
        (reference_dir / "yang2018_co2_ethane_density.json").read_text(
            encoding="utf-8-sig"
        )
    )["PureOrMixtureData"][0]
    for index, point in enumerate(density["NumValues"]):
        values = variables(point)
        rho, rho_u = property_value(point)
        x_c2 = values[1]
        row = row_template(
            f"yang2018_co2_c2_{index:03d}",
            "Yang et al. (2018)",
            "co2_c2",
            "density",
        )
        row.update(
            temperature_K=values[2],
            pressure_Pa=values[3] * 1000.0,
            density_kg_m3=rho,
            uncertainty=rho_u,
            y_co2=1.0 - x_c2,
            y_c2=x_c2,
        )
        rows.append(row)

    messabeb = json.loads(
        (reference_dir / "messabeb2016_sw_co2_brine.json").read_text(
            encoding="utf-8-sig"
        )
    )["PureOrMixtureData"][0]
    selected = [
        point for point in messabeb["NumValues"]
        if math.isclose(variables(point)[1], 1.0, abs_tol=1.0e-12)
    ]
    for index, point in enumerate(selected):
        values = variables(point)
        molality, molality_u = property_value(point)
        water_moles = 1.0 / WATER_MOLAR_MASS
        x_co2 = molality / (water_moles + molality)
        row = row_template(
            f"messabeb2016_h2o_co2_{index:03d}",
            "Messabeb et al. (2016)",
            "h2o_co2",
            "solubility",
        )
        row.update(
            temperature_K=values[2],
            pressure_Pa=values[3] * 1000.0,
            salinity_molal=1.0,
            uncertainty=molality_u,
            x_h2o=1.0 - x_co2,
            x_co2=x_co2,
            y_co2=1.0,
        )
        rows.append(row)

    return rows


def assign_split(rows: list[dict[str, Any]]) -> None:
    groups: dict[tuple[str, str, float, float], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        temperature = round(float(row["temperature_K"]), 2)
        composition = 0.0
        if row["kind"] == "density":
            composition = round(float(row["y_c2"]), 5)
        groups[(row["source"], row["pair"], temperature, composition)].append(row)

    for key, group in groups.items():
        group.sort(key=lambda item: (float(item["pressure_Pa"]), item["id"]))
        digest = hashlib.sha256(f"{SEED}:{key}".encode()).digest()
        offset = digest[0] % 4
        for index, row in enumerate(group):
            row["split"] = "validation" if index % 4 == offset else "training"

    for pair in sorted({row["pair"] for row in rows}):
        pair_rows = [row for row in rows if row["pair"] == pair]
        if not any(row["split"] == "validation" for row in pair_rows):
            pair_rows[-1]["split"] = "validation"
        if not any(row["split"] == "training" for row in pair_rows):
            pair_rows[0]["split"] = "training"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    rows = load_rows(args.reference_dir)
    assign_split(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = args.output_dir / "common_calibration_data.csv"
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in rows:
        counts[row["pair"]][row["split"]] += 1
    manifest = {
        "protocol_seed": SEED,
        "raw_data": "NIST ThermoML JSON; files and SHA-256 are in reference_data/manifest.json",
        "split_rule": "deterministic within-isotherm 3:1 training/validation split",
        "rows": len(rows),
        "counts": {pair: dict(value) for pair, value in counts.items()},
        "transformations": [
            "paired pressure and vapor-composition records by identical ThermoML variables",
            "converted kPa to Pa",
            "converted CO2 molality to salt-free H2O/CO2 liquid mole fraction",
            "used only the 1 mol/kg NaCl subset for the reservoir salinity",
        ],
        "limitations": [
            "the CO2-brine source does not report vapor water; y_CO2=1 is used only for its CO2 fugacity residual",
            "no public single-paper five-component data set contains every requested observable",
            "unmeasured binary pairs are not fitted",
        ],
    }
    (args.output_dir / "dataset_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
