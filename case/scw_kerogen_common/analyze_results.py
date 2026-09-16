#!/usr/bin/env python3
"""Reduce SCW/kerogen displacement outputs to viscosity and producer-composition histories."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {key} in output row")
    return value


def write_rows(path: Path, rows: list[dict[str, float | str]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def reduce_producer(
    rows: list[dict[str, str]], pore_volumes_per_day: float
) -> tuple[list[dict[str, float | str]], list[str]]:
    # Step zero is evaluated before the rate-controlled BHP unknowns have been
    # solved.  Its composition ratios can look plausible, but its absolute
    # rates are only initialization diagnostics and would distort the history.
    producer = [
        row for row in rows
        if row.get("name") == "PROD" and int(row["step"]) > 0
    ]
    if not producer:
        raise ValueError("well_history.csv contains no solved PROD rows")
    component_columns = [
        key for key in producer[0] if key.startswith("m_component_")
    ]
    if not component_columns:
        raise ValueError("well_history.csv predates per-component well-rate output")

    output: list[dict[str, float | str]] = []
    for row in producer:
        rates = {key: max(0.0, -number(row, key)) for key in component_columns}
        total = sum(rates.values())
        time_day = number(row, "time")
        record: dict[str, float | str] = {
            "step": int(row["step"]),
            "time_day": time_day,
            "pvi": time_day * pore_volumes_per_day,
            "producer_total_component_mass_rate_kg_s": total,
        }
        hydrocarbon_total = sum(
            rate for key, rate in rates.items() if key != "m_component_H2O"
        )
        for key, rate in rates.items():
            token = key.removeprefix("m_component_")
            record[f"{token}_mass_rate_kg_s"] = rate
            record[f"{token}_total_mass_fraction"] = rate / total if total > 0.0 else 0.0
            if token != "H2O":
                record[f"{token}_hydrocarbon_mass_fraction"] = (
                    rate / hydrocarbon_total if hydrocarbon_total > 0.0 else 0.0
                )
        output.append(record)
    return output, component_columns


def reduce_viscosity(
    rows: list[dict[str, str]], pore_volumes_per_day: float
) -> list[dict[str, float | str]]:
    output: list[dict[str, float | str]] = []
    for row in rows:
        so = number(row, "s_o_avg")
        sg = number(row, "s_g_avg")
        mu_o = number(row, "mu_o_avg")
        mu_g = number(row, "mu_g_avg")
        time_day = number(row, "time")
        denominator = so + sg
        mu_hc = (so * mu_o + sg * mu_g) / denominator if denominator > 0.0 else 0.0
        output.append(
            {
                "step": int(row["step"]),
                "time_day": time_day,
                "pvi": time_day * pore_volumes_per_day,
                "oil_viscosity_Pa_s": mu_o,
                "gas_viscosity_Pa_s": mu_g,
                "water_rich_viscosity_Pa_s": number(row, "mu_w_avg"),
                "hydrocarbon_phase_volume_weighted_viscosity_Pa_s": mu_hc,
                "oil_saturation": so,
                "gas_saturation": sg,
                "water_saturation": number(row, "s_w_avg"),
            }
        )
    return output


def first_threshold_crossing(
    rows: list[dict[str, float | str]], key: str, threshold: float
) -> dict[str, float | str | None]:
    """Locate a producer-composition threshold without using phase-slot labels.

    The interpolated PVI is a reporting aid between saved well-history records;
    it is not interpreted as sub-timestep resolution of the transport solution.
    """
    previous: dict[str, float | str] | None = None
    for row in rows:
        value = float(row[key])
        if value >= threshold:
            recorded_pvi = float(row["pvi"])
            interpolated_pvi = recorded_pvi
            if previous is not None:
                previous_value = float(previous[key])
                previous_pvi = float(previous["pvi"])
                if previous_value < threshold and value > previous_value:
                    fraction = (threshold - previous_value) / (value - previous_value)
                    interpolated_pvi = previous_pvi + fraction * (
                        recorded_pvi - previous_pvi
                    )
            return {
                "threshold": threshold,
                "first_recorded_step": int(row["step"]),
                "first_recorded_pvi": recorded_pvi,
                "linearly_interpolated_pvi": interpolated_pvi,
            }
        previous = row
    return {
        "threshold": threshold,
        "first_recorded_step": None,
        "first_recorded_pvi": None,
        "linearly_interpolated_pvi": None,
    }


def plot(output: Path, composition: list[dict[str, float | str]], viscosity: list[dict[str, float | str]]) -> bool:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    figure, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), constrained_layout=True)
    pvi = [float(row["pvi"]) for row in viscosity]
    axes[0].semilogy(pvi, [float(row["oil_viscosity_Pa_s"]) for row in viscosity], label="oil slot")
    axes[0].semilogy(pvi, [float(row["gas_viscosity_Pa_s"]) for row in viscosity], label="gas slot")
    axes[0].semilogy(pvi, [float(row["water_rich_viscosity_Pa_s"]) for row in viscosity], label="water-rich slot")
    axes[0].set(xlabel="injected pore volume (PVI)", ylabel="phase-volume-weighted viscosity (Pa s)")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend()

    composition_keys = [
        key for key in composition[0]
        if key.endswith("_total_mass_fraction")
    ]
    composition_pvi = [float(row["pvi"]) for row in composition]
    for key in composition_keys:
        axes[1].plot(
            composition_pvi,
            [float(row[key]) for row in composition],
            label=key.removesuffix("_total_mass_fraction"),
        )
    axes[1].set(xlabel="injected pore volume (PVI)", ylabel="producer total mass fraction", ylim=(0.0, 1.0))
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    figure.savefig(output / "scw_kerogen_diagnostics.png", dpi=180)
    plt.close(figure)
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--pore-volumes-per-day",
        type=float,
        default=0.25,
        help="constant injection rate used to convert time in days to PVI",
    )
    args = parser.parse_args()
    if not math.isfinite(args.pore_volumes_per_day) or args.pore_volumes_per_day <= 0.0:
        parser.error("--pore-volumes-per-day must be finite and positive")
    output = args.output or args.result_dir / "analysis"
    output.mkdir(parents=True, exist_ok=True)

    composition, source_columns = reduce_producer(
        read_rows(args.result_dir / "well_history.csv"),
        args.pore_volumes_per_day,
    )
    viscosity = reduce_viscosity(
        read_rows(args.result_dir / "reservoir_diagnostics.csv"),
        args.pore_volumes_per_day,
    )
    write_rows(output / "producer_composition.csv", composition)
    write_rows(output / "viscosity_history.csv", viscosity)
    plotted = plot(output, composition, viscosity)

    h2o_key = "H2O_total_mass_fraction"
    if h2o_key not in composition[0]:
        raise ValueError("producer output does not contain the H2O component rate")
    h2o_breakthrough = {
        "metric": "produced H2O component mass fraction",
        "initial_solved_value": float(composition[0][h2o_key]),
        "final_value": float(composition[-1][h2o_key]),
        "crossings": [
            first_threshold_crossing(composition, h2o_key, threshold)
            for threshold in (0.10, 0.50)
        ],
        "note": (
            "Thresholds use component mass rates and are independent of simulator "
            "phase-slot relabeling. Interpolated PVI is linear between saved records."
        ),
    }

    summary = {
        "result_dir": str(args.result_dir.resolve()),
        "pore_volumes_per_day": args.pore_volumes_per_day,
        "component_rate_columns": source_columns,
        "producer_final": composition[-1],
        "producer_h2o_breakthrough": h2o_breakthrough,
        "viscosity_final": viscosity[-1],
        "plot_written": plotted,
        "interpretation_note": (
            "Oil/gas/water are simulator phase roles. At strict-SCW conditions "
            "the hydrocarbon-enriched fluid can occupy the gas slot; inspect "
            "phase-state and saturation columns before calling it conventional oil."
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"wrote SCW/kerogen analysis to {output}")


if __name__ == "__main__":
    main()
