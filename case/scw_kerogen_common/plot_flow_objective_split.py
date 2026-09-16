#!/usr/bin/env python3
"""Create separate single-heavy-component and multicomponent figure suites."""

from __future__ import annotations

import argparse
import json
import math
import platform
from collections import defaultdict
from datetime import date
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import plot_flow_objective_figures as base


SYSTEM_NAMES = {
    "binary": "H₂O–squalane",
    "lmh": "H₂O–nC₄–nC₁₀–squalane",
}


def figure_saturation_maps(
    runs: dict[str, Path], output: Path, system: str
) -> list[str]:
    eos_values = tuple(runs)
    fig, axes = plt.subplots(
        len(eos_values),
        4,
        figsize=(7.09, 1.28 + 1.08 * len(eos_values)),
        layout="constrained",
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    image = None
    for row_index, eos in enumerate(eos_values):
        for col_index, (step, pvi) in enumerate(base.SNAPSHOTS):
            field = base.reshape_column(
                runs[eos] / f"solution_step_{step}.csv", "water_saturation"
            )
            image = axes[row_index, col_index].imshow(
                field,
                origin="lower",
                extent=(0, base.LX, 0, base.LY),
                aspect="auto",
                interpolation="nearest",
                cmap="cividis",
                vmin=0.0,
                vmax=1.0,
            )
            axes[row_index, col_index].scatter(
                [0.01, 1.19],
                [0.0475, 0.0475],
                marker="x",
                s=16,
                color="white",
                linewidth=0.8,
            )
            axes[row_index, col_index].set_title(f"{pvi:.2f} PVI")
            if col_index == 0:
                axes[row_index, col_index].set_ylabel(f"{eos}\ny (m)")
            if row_index == len(eos_values) - 1:
                axes[row_index, col_index].set_xlabel("x (m)")
    fig.colorbar(
        image,
        ax=axes,
        label="Water saturation, $S_w$",
        fraction=0.025,
        pad=0.02,
    )
    fig.suptitle(f"{SYSTEM_NAMES[system]}: water-front evolution", fontsize=10)
    return base.export(fig, output, "fig02_2d_water_saturation")


def binary_overall_fields(path: Path) -> dict[str, np.ndarray]:
    rows = base.read_rows(path)
    output: dict[str, np.ndarray] = {}
    for index, component in enumerate(("H2O", "Heavy_squalane")):
        values = np.empty(base.NX * base.NY)
        column = f"z_{index}_{component}"
        for row in rows:
            values[int(row["input_index"])] = base.number(row, column)
        output[component] = values.reshape(base.NY, base.NX)
    return output


def figure_binary_composition_maps(
    runs: dict[str, Path], output: Path
) -> list[str]:
    fields = {
        eos: binary_overall_fields(path / "phase_state_step_400.csv")
        for eos, path in runs.items()
    }
    definitions = (
        (
            "H2O",
            "Overall H₂O (linear, 0.90–1.00)",
            mpl.colors.Normalize(0.90, 1.00),
            "cividis",
        ),
        (
            "Heavy_squalane",
            "Overall squalane (linear, 0–0.10)",
            mpl.colors.Normalize(0.0, 0.10),
            "magma",
        ),
    )
    eos_values = tuple(runs)
    fig, axes = plt.subplots(
        len(eos_values),
        2,
        figsize=(7.09, 1.25 + 1.15 * len(eos_values)),
        layout="constrained",
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for col, (key, title, norm, cmap) in enumerate(definitions):
        image = None
        for row, eos in enumerate(eos_values):
            image = axes[row, col].imshow(
                fields[eos][key],
                origin="lower",
                extent=(0, base.LX, 0, base.LY),
                aspect="auto",
                interpolation="nearest",
                norm=norm,
                cmap=cmap,
            )
            axes[row, col].scatter(
                [0.01, 1.19],
                [0.0475, 0.0475],
                marker="x",
                s=16,
                color="white",
                linewidth=0.8,
            )
            if col == 0:
                axes[row, col].set_ylabel(f"{eos}\ny (m)")
            if row == len(eos_values) - 1:
                axes[row, col].set_xlabel("x (m)")
            if row == 0:
                axes[row, col].set_title(title)
        fig.colorbar(image, ax=axes[:, col], fraction=0.035, pad=0.015)
    fig.suptitle("H₂O–squalane: overall composition at 1.00 PVI", fontsize=10)
    return base.export(fig, output, "fig03_2d_composition_at_1pvi")


def figure_viscosity(
    runs: dict[str, Path], output: Path, system: str
) -> tuple[list[str], list[dict[str, object]]]:
    data = {eos: base.diagnostics(path, system, eos) for eos, path in runs.items()}
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.05), layout="constrained")
    for col_index, phase in enumerate(("o", "w")):
        ax = axes[col_index]
        for eos, rows in data.items():
            valid = [r for r in rows if float(r[f"mu_{phase}_avg_mPa_s"]) > 0]
            if not valid:
                continue
            x = np.array([float(r["pvi"]) for r in valid])
            avg = np.array([float(r[f"mu_{phase}_avg_mPa_s"]) for r in valid])
            lo = np.array([float(r[f"mu_{phase}_min_mPa_s"]) for r in valid])
            hi = np.array([float(r[f"mu_{phase}_max_mPa_s"]) for r in valid])
            color = "#1B4F72" if phase == "o" else "#0072B2"
            ax.fill_between(
                x,
                np.maximum(lo, 1e-8),
                np.maximum(hi, 1e-8),
                color=color,
                alpha=0.08,
            )
            ax.plot(
                x,
                avg,
                color=color,
                linestyle=base.EOS_STYLES[eos],
                label=eos,
            )
        phase_name = "Hydrocarbon-rich" if phase == "o" else "Water-rich"
        ax.set(
            yscale="log",
            xlim=(0, 1),
            xlabel="Injected pore volume (PVI)",
            ylabel=f"{phase_name} viscosity (mPa·s; log)",
            title=f"{phase_name} phase",
        )
        ax.grid(True, which="both", color="#D9D9D9", linewidth=0.55)
        ax.legend(title="EOS", ncol=len(runs))
        base.panel(ax, f"({chr(97 + col_index)})")
    fig.suptitle(f"{SYSTEM_NAMES[system]}: viscosity response", fontsize=10)
    flat = [row for rows in data.values() for row in rows]
    return base.export(fig, output, "fig04_oil_water_viscosity"), flat


def figure_producer(
    runs: dict[str, Path], output: Path, system: str
) -> tuple[list[str], list[dict[str, object]]]:
    flat = [
        row
        for eos, path in runs.items()
        for row in base.producer_long(path, system, eos)
    ]
    if system == "binary":
        fig, ax = plt.subplots(1, 1, figsize=(5.0, 3.35), layout="constrained")
        axes = [ax]
        base.plot_component_lines(ax, flat, tuple(runs), normalized=False)
        titles = ("Total producer composition",)
    else:
        fig, axes_raw = plt.subplots(1, 3, figsize=(7.09, 2.9), layout="constrained", sharex=True)
        axes = list(axes_raw)
        base.plot_component_lines(axes[0], flat, tuple(runs), normalized=False)
        base.plot_component_lines(axes[1], flat, ("PR",), normalized=True)
        base.plot_component_lines(axes[2], flat, ("CPA",), normalized=True)
        titles = (
            "Total producer composition",
            "PR: HC-normalized",
            "CPA: HC-normalized",
        )
    for i, ax in enumerate(axes):
        ax.set(
            xlabel="Injected pore volume (PVI)",
            ylabel="Mass fraction",
            xlim=(0, 1),
            ylim=(0, 1),
            title=f"({chr(97+i)})  {titles[i]}",
        )
        ax.grid(True, color="#D9D9D9", linewidth=0.55)
        ax.legend(fontsize=6.1, ncol=2)
    fig.suptitle(f"{SYSTEM_NAMES[system]}: producer composition", fontsize=10)
    return base.export(fig, output, "fig05_producer_composition"), flat


def figure_recovery_balance(
    runs: dict[str, Path], output: Path, system: str
) -> tuple[list[str], list[dict[str, object]]]:
    flat = [
        row
        for eos, path in runs.items()
        for row in base.recovery_long(path, system, eos)
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.0), layout="constrained")
    components = (
        ("Heavy_squalane",)
        if system == "binary"
        else ("Light_nC4", "Middle_nC10", "Heavy_squalane")
    )
    for component in components:
        for eos in runs:
            rows = sorted(
                (
                    r
                    for r in flat
                    if r["eos"] == eos and r["component"] == component
                ),
                key=lambda r: float(r["pvi"]),
            )
            if rows:
                axes[0].plot(
                    [float(r["pvi"]) for r in rows],
                    [float(r["recovery_fraction_of_initial"]) for r in rows],
                    color=base.COMPONENT_COLORS[component],
                    linestyle=base.EOS_STYLES[eos],
                    label=f"{base.label_component(component)}–{eos}",
                )
    for eos in runs:
        grouped: dict[float, list[float]] = defaultdict(list)
        for row in flat:
            if row["eos"] == eos:
                grouped[float(row["pvi"])].append(
                    abs(float(row["relative_balance_error"]))
                )
        points = sorted(
            (pvi, max(values))
            for pvi, values in grouped.items()
            if max(values) > 0
        )
        axes[1].semilogy(
            [point[0] for point in points],
            [point[1] for point in points],
            color=base.SYSTEM_COLORS[system],
            linestyle=base.EOS_STYLES[eos],
            label=eos,
        )
    axes[0].set(
        title="Component recovery",
        ylabel="Cumulative recovery / initial inventory",
        ylim=(0, 1),
    )
    axes[1].set(
        title="Maximum component balance error",
        ylabel="Absolute relative error",
    )
    for i, ax in enumerate(axes):
        ax.set(xlabel="Injected pore volume (PVI)", xlim=(0, 1))
        ax.grid(True, which="both", color="#D9D9D9", linewidth=0.55)
        ax.legend(fontsize=6.5, ncol=2)
        base.panel(ax, f"({chr(97+i)})")
    fig.suptitle(f"{SYSTEM_NAMES[system]}: recovery and conservation", fontsize=10)
    return base.export(fig, output, "fig06_recovery_and_mass_balance"), flat


def figure_h2o_bhp(
    producer: list[dict[str, object]],
    runs: dict[str, Path],
    output: Path,
    system: str,
) -> tuple[list[str], list[dict[str, object]]]:
    seen: set[tuple[object, object]] = set()
    h2o_history = []
    for row in producer:
        key = (row["eos"], row["step"])
        if key not in seen:
            seen.add(key)
            h2o_history.append(
                {
                    field: row[field]
                    for field in (
                        "system",
                        "eos",
                        "step",
                        "pvi",
                        "producer_h2o_mass_fraction",
                        "reservoir_watercut",
                    )
                }
            )
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.0), layout="constrained")
    for eos, path in runs.items():
        rows = sorted(
            (r for r in h2o_history if r["eos"] == eos),
            key=lambda r: float(r["pvi"]),
        )
        axes[0].plot(
            [float(r["pvi"]) for r in rows],
            [float(r["producer_h2o_mass_fraction"]) for r in rows],
            color=base.SYSTEM_COLORS[system],
            linestyle=base.EOS_STYLES[eos],
            label=eos,
        )
        raw = [
            r
            for r in base.read_rows(path / "well_history.csv")
            if r.get("name") == "PROD" and int(r["step"]) > 0
        ]
        axes[1].plot(
            [base.number(r, "time") * base.PVI_PER_DAY for r in raw],
            [base.number(r, "bhp") for r in raw],
            color=base.SYSTEM_COLORS[system],
            linestyle=base.EOS_STYLES[eos],
            label=eos,
        )
    axes[0].set(
        ylabel="Produced H₂O mass fraction",
        ylim=(0, 1),
        title="Component-based water breakthrough",
    )
    axes[1].set(ylabel="Producer BHP (bar)", title="Producer pressure response")
    for i, ax in enumerate(axes):
        ax.set(xlabel="Injected pore volume (PVI)", xlim=(0, 1))
        ax.grid(True, color="#D9D9D9", linewidth=0.55)
        ax.legend(title="EOS", ncol=len(runs), fontsize=6.5)
        base.panel(ax, f"({chr(97+i)})")
    fig.suptitle(f"{SYSTEM_NAMES[system]}: producer response", fontsize=10)
    return base.export(fig, output, "fig07_producer_h2o_and_bhp"), h2o_history


def create_suite(
    root: Path,
    output: Path,
    system: str,
    runs: dict[str, Path],
    calibration: Path,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    figures: list[str] = []
    if system == "binary":
        figures += base.figure_calibration(calibration, output)
    figures += figure_saturation_maps(runs, output, system)
    if system == "binary":
        figures += figure_binary_composition_maps(runs, output)
    else:
        figures += base.figure_composition_maps(runs, output)
    names, viscosity = figure_viscosity(runs, output, system)
    figures += names
    names, producer = figure_producer(runs, output, system)
    figures += names
    names, recovery = figure_recovery_balance(runs, output, system)
    figures += names
    names, h2o_history = figure_h2o_bhp(producer, runs, output, system)
    figures += names

    base.write_rows(output / "viscosity_envelope_long.csv", viscosity)
    base.write_rows(output / "producer_composition_long.csv", producer)
    base.write_rows(output / "recovery_balance_long.csv", recovery)
    base.write_rows(output / "producer_h2o_breakthrough_long.csv", h2o_history)
    manifest = {
        "generated_on": date.today().isoformat(),
        "suite": system,
        "system": SYSTEM_NAMES[system],
        "status": "provisional internal figures; current uncalibrated flow models",
        "python": platform.python_version(),
        "matplotlib": mpl.__version__,
        "source_runs": {eos: str(path) for eos, path in runs.items()},
        "calibration_source": str(calibration) if system == "binary" else None,
        "figures": figures,
        "transformations": [
            "PVI = time_day * 0.25 PV/day",
            "2-D input_index reshaped as index = i + 60*j",
            "saturation maps use a common linear 0–1 color scale without interpolation",
            "producer mass fractions use absolute produced component mass rates",
            "hydrocarbon-normalized fractions exclude H2O from the denominator",
            "viscosity bands are grid-cell min/max, not uncertainty intervals",
            "component recovery = cumulative produced mass / initial inventory",
            "strict zero balance errors remain in CSV and are omitted only from log display",
        ],
        "unavailable_requested_output": {
            "local_2d_viscosity_maps": "cellwise viscosity is absent from current solution snapshots",
            "lmh_sw_full_run": (
                "four-component SW failed during initial phase transition"
                if system == "lmh"
                else None
            ),
        },
        "uncertainty": "none estimated; EOS curves and spatial min/max bands are not confidence intervals",
    }
    (output / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    if output.exists() and not args.force and any(output.iterdir()):
        raise FileExistsError("output directory is not empty; use --force")

    binary_root = root / "case/scw_kerogen_squalane_2d/results/current_20260916"
    lmh_root = root / "case/scw_kerogen_lmh_2d/results/current_20260916"
    binary_runs = {
        "PR": binary_root / "pr_np1_lu",
        "SW": binary_root / "sw_np1_lu",
        "CPA": binary_root / "cpa_np1_lu",
    }
    lmh_runs = {
        "PR": lmh_root / "pr_np1_lu",
        "CPA": lmh_root / "cpa_np1_lu",
    }
    calibration = root / "tools/example/scw_binary_calibration/results"

    base.configure_style()
    create_suite(root, output / "single_component", "binary", binary_runs, calibration)
    create_suite(root, output / "multicomponent", "lmh", lmh_runs, calibration)
    print(f"wrote separated suites to {output}")


if __name__ == "__main__":
    main()
