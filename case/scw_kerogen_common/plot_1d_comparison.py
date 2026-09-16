#!/usr/bin/env python3
"""Create provenance-bound comparison figures for the two SCW/kerogen 1-D cases."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image


PVI_PER_DAY = 0.25
CORE_LENGTH_M = 1.20
GRID_CELLS = 60
SNAPSHOT_STEPS = (0, 10, 20, 30, 40)

COLORS = {
    "H2O": "#0072B2",
    "Light_nC4": "#D55E00",
    "Middle_nC10": "#009E73",
    "Heavy_squalane": "#CC79A7",
    "hydrocarbon": "#D55E00",
    "water": "#0072B2",
    "gas": "#000000",
    "binary": "#0072B2",
    "lmh": "#D55E00",
}

LINESTYLES = ("-", "--", "-.", ":", (0, (5, 2, 1, 2)))
MARKERS = ("o", "s", "^", "D", "v")
TIME_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty CSV: {path}")
    return rows


def finite(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite value for {key}")
    return value


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def component_label(token: str) -> str:
    return {
        "H2O": "H₂O",
        "Light_nC4": "nC₄",
        "Middle_nC10": "nC₁₀",
        "Heavy_squalane": "Squalane",
    }.get(token, token)


def producer_composition(result_dir: Path, case: str) -> list[dict[str, object]]:
    source = read_rows(result_dir / "well_history.csv")
    producer = [
        row for row in source
        if row.get("name") == "PROD" and int(row["step"]) > 0
    ]
    component_columns = [
        key for key in producer[0] if key.startswith("m_component_")
    ]
    output: list[dict[str, object]] = []
    for row in producer:
        rates = {key: max(0.0, -finite(row, key)) for key in component_columns}
        total = sum(rates.values())
        if total <= 0.0:
            raise ValueError("non-positive producer component mass rate")
        time_day = finite(row, "time")
        for key, rate in rates.items():
            output.append(
                {
                    "case": case,
                    "step": int(row["step"]),
                    "time_day": time_day,
                    "pvi": time_day * PVI_PER_DAY,
                    "component": key.removeprefix("m_component_"),
                    "mass_rate_kg_s": rate,
                    "total_mass_fraction": rate / total,
                }
            )
    return output


def viscosity_saturation(result_dir: Path, case: str) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for row in read_rows(result_dir / "reservoir_diagnostics.csv"):
        time_day = finite(row, "time")
        so = finite(row, "s_o_avg")
        sg = finite(row, "s_g_avg")
        denominator = so + sg
        mu_o = finite(row, "mu_o_avg")
        mu_g = finite(row, "mu_g_avg")
        output.append(
            {
                "case": case,
                "step": int(row["step"]),
                "time_day": time_day,
                "pvi": time_day * PVI_PER_DAY,
                "hydrocarbon_viscosity_mPa_s": (
                    1000.0 * (so * mu_o + sg * mu_g) / denominator
                    if denominator > 0.0 else 0.0
                ),
                "water_viscosity_mPa_s": 1000.0 * finite(row, "mu_w_avg"),
                "oil_saturation": so,
                "gas_saturation": sg,
                "water_saturation": finite(row, "s_w_avg"),
            }
        )
    return output


def spatial_water_saturation(result_dir: Path, case: str) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for step in SNAPSHOT_STEPS:
        rows = read_rows(result_dir / f"solution_step_{step}.csv")
        if len(rows) != GRID_CELLS:
            raise ValueError(f"expected {GRID_CELLS} cells at step {step}")
        pvi = step * 0.1 * PVI_PER_DAY
        for row in rows:
            cell = int(row["input_index"])
            output.append(
                {
                    "case": case,
                    "step": step,
                    "pvi": pvi,
                    "cell": cell,
                    "x_m": (cell + 0.5) * CORE_LENGTH_M / GRID_CELLS,
                    "water_saturation": finite(row, "water_saturation"),
                    "oil_saturation": finite(row, "oil_saturation"),
                    "gas_saturation": finite(row, "gas_saturation"),
                }
            )
    return output


def mass_balance(result_dir: Path, case: str) -> list[dict[str, object]]:
    grouped: dict[tuple[int, float], list[float]] = defaultdict(list)
    for row in read_rows(result_dir / "component_mass_balance.csv"):
        grouped[(int(row["step"]), finite(row, "time_day"))].append(
            abs(finite(row, "relative_error"))
        )
    return [
        {
            "case": case,
            "step": step,
            "time_day": time_day,
            "pvi": time_day * PVI_PER_DAY,
            "max_abs_component_relative_error": max(errors),
        }
        for (step, time_day), errors in sorted(grouped.items())
    ]


def series(rows: Iterable[dict[str, object]], key: str) -> list[float]:
    return [float(row[key]) for row in rows]


def export_figure(fig: mpl.figure.Figure, output: Path, stem: str) -> list[str]:
    written: list[str] = []
    for suffix in ("png", "pdf", "svg"):
        target = output / f"{stem}.{suffix}"
        fig.savefig(
            target,
            dpi=600 if suffix == "png" else 300,
            facecolor="white",
            transparent=False,
            metadata={"Title": stem, "Creator": "Matplotlib"},
        )
        if suffix == "png":
            # Matplotlib writes an RGBA container even for an opaque white
            # canvas.  Strip the unused alpha channel so downstream checks do
            # not mistake the figure for a transparent scientific image.
            with Image.open(target) as source:
                rgb = source.convert("RGB")
            rgb.save(target, dpi=(600, 600), optimize=True)
        written.append(target.name)
    plt.close(fig)
    return written


def add_panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.11,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
    )


def make_spatial_figure(
    output: Path,
    spatial: dict[str, list[dict[str, object]]],
) -> list[str]:
    fig, axes = plt.subplots(2, 1, figsize=(7.09, 5.6), layout="constrained", sharex=True)
    for panel, (case, title) in enumerate(
        (("binary", "H₂O–squalane"), ("lmh", "H₂O–nC₄–nC₁₀–squalane"))
    ):
        ax = axes[panel]
        for index, step in enumerate(SNAPSHOT_STEPS):
            rows = [row for row in spatial[case] if int(row["step"]) == step]
            ax.plot(
                series(rows, "x_m"),
                series(rows, "water_saturation"),
                color=TIME_COLORS[index],
                linestyle=LINESTYLES[index],
                linewidth=1.5,
                label=f"{float(rows[0]['pvi']):.2f} PVI",
            )
        ax.set(ylabel="Water saturation", ylim=(0.0, 1.0), title=title)
        ax.grid(True, color="#D9D9D9", linewidth=0.6)
        ax.legend(ncol=5, fontsize=7.5, loc="upper right")
        add_panel_label(ax, f"({chr(97 + panel)})")
    axes[-1].set_xlabel("Distance from injector, x (m)")
    axes[-1].set_xlim(0.0, CORE_LENGTH_M)
    return export_figure(fig, output, "fig01_spatial_water_saturation")


def make_composition_figure(
    output: Path,
    composition: dict[str, list[dict[str, object]]],
) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.25), layout="constrained", sharey=True)
    order = {
        "binary": ("H2O", "Heavy_squalane"),
        "lmh": ("H2O", "Light_nC4", "Middle_nC10", "Heavy_squalane"),
    }
    titles = {"binary": "H₂O–squalane", "lmh": "H₂O–nC₄–nC₁₀–squalane"}
    for panel, case in enumerate(("binary", "lmh")):
        ax = axes[panel]
        for index, component in enumerate(order[case]):
            rows = [row for row in composition[case] if row["component"] == component]
            ax.plot(
                series(rows, "pvi"),
                series(rows, "total_mass_fraction"),
                color=COLORS[component],
                linestyle=LINESTYLES[index],
                marker=MARKERS[index],
                markevery=5,
                markersize=3.2,
                linewidth=1.5,
                label=component_label(component),
            )
        ax.set(xlabel="Injected pore volume (PVI)", xlim=(0.0, 1.0), title=titles[case])
        ax.grid(True, color="#D9D9D9", linewidth=0.6)
        ax.legend(fontsize=8)
        add_panel_label(ax, f"({chr(97 + panel)})")
    axes[0].set(ylabel="Producer component mass fraction", ylim=(0.0, 1.0))
    return export_figure(fig, output, "fig02_producer_composition")


def make_viscosity_saturation_figure(
    output: Path,
    histories: dict[str, list[dict[str, object]]],
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.6), layout="constrained", sharex="col")
    titles = {"binary": "H₂O–squalane", "lmh": "H₂O–nC₄–nC₁₀–squalane"}
    for column, case in enumerate(("binary", "lmh")):
        rows = histories[case]
        pvi = series(rows, "pvi")
        top = axes[0, column]
        top.plot(
            pvi,
            series(rows, "hydrocarbon_viscosity_mPa_s"),
            color=COLORS["hydrocarbon"],
            linestyle="-",
            marker="o",
            markevery=5,
            markersize=3.0,
            label="Hydrocarbon-rich phase",
        )
        top.plot(
            pvi,
            series(rows, "water_viscosity_mPa_s"),
            color=COLORS["water"],
            linestyle="--",
            marker="s",
            markevery=5,
            markersize=3.0,
            label="Water-rich phase",
        )
        top.set(yscale="log", ylim=(0.05, 1.2), title=titles[case])
        top.grid(True, which="both", color="#D9D9D9", linewidth=0.6)
        top.legend(fontsize=7.5)

        bottom = axes[1, column]
        bottom.plot(pvi, series(rows, "oil_saturation"), color=COLORS["hydrocarbon"], label="Oil slot")
        bottom.plot(pvi, series(rows, "water_saturation"), color=COLORS["water"], linestyle="--", label="Water slot")
        bottom.plot(pvi, series(rows, "gas_saturation"), color=COLORS["gas"], linestyle=":", label="Gas slot")
        bottom.set(xlabel="Injected pore volume (PVI)", xlim=(0.0, 1.0), ylim=(0.0, 1.0))
        bottom.grid(True, color="#D9D9D9", linewidth=0.6)
        bottom.legend(fontsize=7.5)
    axes[0, 0].set_ylabel("Average viscosity (mPa·s; log scale)")
    axes[1, 0].set_ylabel("Volume-averaged saturation")
    for index, ax in enumerate(axes.flat):
        add_panel_label(ax, f"({chr(97 + index)})")
    return export_figure(fig, output, "fig03_viscosity_and_saturation")


def make_balance_figure(
    output: Path,
    balances: dict[str, list[dict[str, object]]],
) -> list[str]:
    fig, ax = plt.subplots(figsize=(4.8, 3.35), layout="constrained")
    for index, (case, label) in enumerate(
        (("binary", "H₂O–squalane"), ("lmh", "H₂O–nC₄–nC₁₀–squalane"))
    ):
        rows = [
            row for row in balances[case]
            if float(row["max_abs_component_relative_error"]) > 0.0
        ]
        ax.semilogy(
            series(rows, "pvi"),
            series(rows, "max_abs_component_relative_error"),
            color=COLORS[case],
            linestyle=LINESTYLES[index],
            marker=MARKERS[index],
            markevery=5,
            markersize=3.2,
            label=label,
        )
    ax.set(
        xlabel="Injected pore volume (PVI)",
        ylabel="Maximum absolute component relative error",
        xlim=(0.0, 1.0),
    )
    ax.grid(True, which="both", color="#D9D9D9", linewidth=0.6)
    ax.legend(fontsize=8)
    return export_figure(fig, output, "fig04_component_mass_balance")


def first_threshold(rows: list[dict[str, object]], threshold: float) -> float | None:
    h2o = sorted(
        (row for row in rows if row["component"] == "H2O"),
        key=lambda row: float(row["pvi"]),
    )
    match = next(
        (row for row in h2o if float(row["total_mass_fraction"]) >= threshold),
        None,
    )
    return float(match["pvi"]) if match else None


def front_positions(rows: list[dict[str, object]]) -> dict[str, float | None]:
    output: dict[str, float | None] = {}
    for step in SNAPSHOT_STEPS[1:]:
        selected = [
            row for row in rows
            if int(row["step"]) == step and float(row["water_saturation"]) >= 0.5
        ]
        output[f"{step * 0.1 * PVI_PER_DAY:.2f}"] = (
            max(float(row["x_m"]) for row in selected) if selected else None
        )
    return output


def build_summary(
    composition: dict[str, list[dict[str, object]]],
    histories: dict[str, list[dict[str, object]]],
    spatial: dict[str, list[dict[str, object]]],
    balances: dict[str, list[dict[str, object]]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for case in ("binary", "lmh"):
        final_pvi = max(float(row["pvi"]) for row in composition[case])
        final_composition = {
            str(row["component"]): float(row["total_mass_fraction"])
            for row in composition[case]
            if math.isclose(float(row["pvi"]), final_pvi)
        }
        final_history = max(histories[case], key=lambda row: float(row["pvi"]))
        result[case] = {
            "final_pvi": final_pvi,
            "h2o_mass_fraction_50pct_first_pvi": first_threshold(composition[case], 0.5),
            "producer_final_mass_fractions": final_composition,
            "final_hydrocarbon_viscosity_mPa_s": final_history["hydrocarbon_viscosity_mPa_s"],
            "final_water_viscosity_mPa_s": final_history["water_viscosity_mPa_s"],
            "final_average_oil_saturation": final_history["oil_saturation"],
            "final_average_water_saturation": final_history["water_saturation"],
            "max_abs_component_relative_error_over_run": max(
                float(row["max_abs_component_relative_error"])
                for row in balances[case]
            ),
            "farthest_cell_center_with_Sw_at_least_0p5_m": front_positions(spatial[case]),
        }
    return result


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary-dir", type=Path, required=True)
    parser.add_argument("--lmh-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    expected = [
        args.output_dir / f"fig{index:02d}_{name}.png"
        for index, name in (
            (1, "spatial_water_saturation"),
            (2, "producer_composition"),
            (3, "viscosity_and_saturation"),
            (4, "component_mass_balance"),
        )
    ]
    if not args.force and any(path.exists() for path in expected):
        raise FileExistsError("output exists; use --force to replace generated artifacts")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    result_dirs = {"binary": args.binary_dir, "lmh": args.lmh_dir}
    composition = {
        case: producer_composition(path, case) for case, path in result_dirs.items()
    }
    histories = {
        case: viscosity_saturation(path, case) for case, path in result_dirs.items()
    }
    spatial = {
        case: spatial_water_saturation(path, case) for case, path in result_dirs.items()
    }
    balances = {
        case: mass_balance(path, case) for case, path in result_dirs.items()
    }

    write_rows(
        args.output_dir / "producer_composition_long.csv",
        composition["binary"] + composition["lmh"],
    )
    write_rows(
        args.output_dir / "viscosity_saturation_long.csv",
        histories["binary"] + histories["lmh"],
    )
    write_rows(
        args.output_dir / "spatial_saturation_long.csv",
        spatial["binary"] + spatial["lmh"],
    )
    write_rows(
        args.output_dir / "mass_balance_max_long.csv",
        balances["binary"] + balances["lmh"],
    )

    configure_style()
    figures: list[str] = []
    figures += make_spatial_figure(args.output_dir, spatial)
    figures += make_composition_figure(args.output_dir, composition)
    figures += make_viscosity_saturation_figure(args.output_dir, histories)
    figures += make_balance_figure(args.output_dir, balances)

    summary = build_summary(composition, histories, spatial, balances)
    (args.output_dir / "summary_metrics.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "generated_on": date.today().isoformat(),
        "status": "provisional internal scientific figures; human review required",
        "python": platform.python_version(),
        "matplotlib": mpl.__version__,
        "source_directories": {
            case: str(path.resolve()) for case, path in result_dirs.items()
        },
        "transformations": [
            "PVI = time_day * 0.25 PV/day",
            "producer mass fractions = absolute produced component mass rate / sum of absolute produced component mass rates",
            "hydrocarbon viscosity = saturation-volume-weighted oil/gas-slot viscosity; gas saturation is zero in these completed runs",
            "spatial curves use cell-center values without interpolation or smoothing",
            "mass-balance curve is the maximum absolute component relative error at each reported time",
            "zero mass-balance error at initial time is omitted from the logarithmic display only and retained in the CSV",
        ],
        "uncertainty": "not estimated; deterministic single-run outputs",
        "missing_data": "none in plotted completed histories",
        "figures": figures,
        "alt_text": {
            "fig01": "Two stacked line charts show water saturation along the 1.2 m one-dimensional core at 0, 0.25, 0.50, 0.75, and 1.00 injected pore volumes. In both fluids the high-water region advances from the injector at x=0 toward the producer at x=1.2 m; the four-component front is sharper and farther advanced at equal PVI.",
            "fig02": "Two line charts show producer component mass fractions versus injected pore volume. Water rises and hydrocarbon components fall after breakthrough. In the four-component case squalane remains the dominant hydrocarbon component, while nC4 and nC10 decline to small fractions.",
            "fig03": "A four-panel figure compares average hydrocarbon-rich and water-rich viscosities on logarithmic axes and average phase saturations on linear axes. The binary hydrocarbon-rich phase is more viscous than the four-component phase, while water-rich viscosity changes little. Average water saturation rises throughout both runs.",
            "fig04": "A logarithmic line chart shows maximum absolute component relative mass-balance error versus injected pore volume for both one-dimensional cases; errors remain small over the completed simulations.",
        },
    }
    (args.output_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {len(figures)} figure files to {args.output_dir}")


if __name__ == "__main__":
    main()
