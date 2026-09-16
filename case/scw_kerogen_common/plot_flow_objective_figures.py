#!/usr/bin/env python3
"""Plot the SCW/kerogen flow-objective figure set from existing CSV outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from collections import defaultdict
from datetime import date
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


NX, NY = 60, 20
LX, LY = 1.20, 0.10
PVI_PER_DAY = 0.25
SNAPSHOTS = ((100, 0.25), (200, 0.50), (300, 0.75), (400, 1.00))

COMPONENTS = ("H2O", "Light_nC4", "Middle_nC10", "Heavy_squalane")
COMPONENT_COLORS = {
    "H2O": "#0072B2",
    "Light_nC4": "#D55E00",
    "Middle_nC10": "#009E73",
    "Heavy_squalane": "#CC79A7",
}
EOS_STYLES = {"PR": "-", "SW": "--", "CPA": "-."}
EOS_MARKERS = {"PR": "o", "SW": "s", "CPA": "^"}
SYSTEM_COLORS = {"binary": "#0072B2", "lmh": "#D55E00"}


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty CSV: {path}")
    return rows


def number(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {key}")
    return value


def label_component(token: str) -> str:
    return {
        "H2O": "H₂O",
        "Light_nC4": "nC₄",
        "Middle_nC10": "nC₁₀",
        "Heavy_squalane": "Squalane",
    }.get(token, token)


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"empty output table: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
        }
    )


def export(fig: mpl.figure.Figure, output: Path, stem: str) -> list[str]:
    names: list[str] = []
    for suffix in ("png", "pdf", "svg"):
        target = output / f"{stem}.{suffix}"
        fig.savefig(
            target,
            dpi=400 if suffix == "png" else 300,
            facecolor="white",
            transparent=False,
            metadata={"Title": stem, "Creator": "Matplotlib"},
        )
        if suffix == "png":
            with Image.open(target) as source:
                rgb = source.convert("RGB")
            rgb.save(target, dpi=(400, 400), optimize=True)
        names.append(target.name)
    plt.close(fig)
    return names


def panel(ax: mpl.axes.Axes, text: str) -> None:
    ax.text(-0.12, 1.04, text, transform=ax.transAxes, fontweight="bold", fontsize=10)


def figure_calibration(calibration: Path, output: Path) -> list[str]:
    pvt = read_rows(calibration / "squalane_pvt_viscosity_predictions.csv")
    lle = read_rows(calibration / "squalane_lle_predictions.csv")
    water = read_rows(calibration / "water_viscosity_check.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.8), layout="constrained")

    pressure_colors = {0.1: "#0072B2", 100.0: "#D55E00", 200.0: "#009E73"}
    for pressure in (0.1, 100.0, 200.0):
        rows = sorted(
            (r for r in pvt if math.isclose(number(r, "pressure_MPa"), pressure)),
            key=lambda r: number(r, "temperature_K"),
        )
        x = [number(r, "temperature_K") for r in rows]
        axes[0, 0].plot(x, [number(r, "target_density_kg_m3") for r in rows],
                        color=pressure_colors[pressure], marker="o", linestyle="none",
                        fillstyle="none", label=f"Exp., {pressure:g} MPa")
        axes[0, 0].plot(x, [number(r, "translated_PR_density_kg_m3") for r in rows],
                        color=pressure_colors[pressure], linestyle="-", label=f"VT-PR, {pressure:g} MPa")
        axes[0, 1].semilogy(x, [1000 * number(r, "target_viscosity_Pa_s") for r in rows],
                           color=pressure_colors[pressure], marker="o", linestyle="none", fillstyle="none")
        axes[0, 1].semilogy(x, [1000 * number(r, "calibrated_LBC_viscosity_Pa_s") for r in rows],
                           color=pressure_colors[pressure], linestyle="-")
    axes[0, 0].set(xlabel="Temperature (K)", ylabel="Squalane density (kg m⁻³)", title="Pure-squalane density")
    axes[0, 0].legend(ncol=2, fontsize=6.5)
    axes[0, 1].set(xlabel="Temperature (K)", ylabel="Squalane viscosity (mPa·s; log)", title="Pure-squalane viscosity")

    conditions = []
    current = {}
    fitted = {}
    for row in lle:
        key = row["id"]
        if key not in conditions:
            conditions.append(key)
        (current if row["variant"] == "current_case" else fitted)[key] = row
    xloc = np.arange(len(conditions))
    targets_wr = [number(current[k], "target_x_water_rich") for k in conditions]
    targets_sr = [number(current[k], "target_x_squalane_rich") for k in conditions]
    axes[1, 0].plot(xloc, targets_wr, "o", color="#0072B2", fillstyle="none", label="Exp., water-rich")
    axes[1, 0].plot(xloc, targets_sr, "s", color="#D55E00", fillstyle="none", label="Exp., squalane-rich")
    for variant, mapping, style in (("Current BIP", current, "--"), ("Fitted BIP", fitted, "-")):
        axes[1, 0].plot(xloc, [number(mapping[k], "predicted_x_water_rich") for k in conditions],
                        color="#0072B2", linestyle=style, label=f"{variant}, water-rich")
        axes[1, 0].plot(xloc, [number(mapping[k], "predicted_x_squalane_rich") for k in conditions],
                        color="#D55E00", linestyle=style, label=f"{variant}, squalane-rich")
    tick_labels = [f"{number(current[k], 'temperature_K'):.1f} K\n{number(current[k], 'pressure_MPa'):.2f} MPa" for k in conditions]
    axes[1, 0].set(xticks=xloc, xticklabels=tick_labels, ylabel="H₂O mole fraction", title="H₂O–squalane coexistence compositions", ylim=(0.65, 1.01))
    axes[1, 0].legend(fontsize=6.2, ncol=2)

    axes[1, 1].plot(xloc, [1000 * number(current[k], "water_rich_viscosity_Pa_s") for k in conditions],
                        color="#0072B2", marker="o", linestyle="-", label="Water-rich, current")
    axes[1, 1].plot(xloc, [1000 * number(current[k], "squalane_rich_viscosity_Pa_s") for k in conditions],
                        color="#D55E00", marker="s", linestyle="-", label="Squalane-rich, current")
    axes[1, 1].plot(xloc, [1000 * number(fitted[k], "water_rich_viscosity_Pa_s") for k in conditions],
                        color="#0072B2", marker="o", linestyle="--", label="Water-rich, fitted BIP")
    axes[1, 1].plot(xloc, [1000 * number(fitted[k], "squalane_rich_viscosity_Pa_s") for k in conditions],
                        color="#D55E00", marker="s", linestyle="--", label="Squalane-rich, fitted BIP")
    # Pure-water reference points provide context but are not mixture viscosity measurements.
    axes[1, 1].plot(xloc, [1000 * number(r, "iapws2008_viscosity_Pa_s") for r in water],
                        color="#000000", marker="^", linestyle=":", label="Pure H₂O, IAPWS-2008")
    axes[1, 1].set(xticks=xloc, xticklabels=tick_labels, ylabel="Viscosity (mPa·s)", yscale="log", title="Coexisting-phase viscosity")
    axes[1, 1].legend(fontsize=6.3)
    for i, ax in enumerate(axes.flat):
        ax.grid(True, which="both", color="#D9D9D9", linewidth=0.55)
        panel(ax, f"({chr(97+i)})")
    return export(fig, output, "fig01_pvt_lle_viscosity_calibration")


def reshape_column(path: Path, column: str) -> np.ndarray:
    rows = read_rows(path)
    if len(rows) != NX * NY:
        raise ValueError(f"expected {NX*NY} cells: {path}")
    values = np.empty(NX * NY)
    for row in rows:
        values[int(row["input_index"])] = number(row, column)
    return values.reshape(NY, NX)


def figure_saturation_maps(runs: dict[str, Path], output: Path) -> list[str]:
    fig, axes = plt.subplots(2, 4, figsize=(7.09, 3.1), layout="constrained", sharex=True, sharey=True)
    image = None
    for row_index, eos in enumerate(("PR", "CPA")):
        for col_index, (step, pvi) in enumerate(SNAPSHOTS):
            field = reshape_column(runs[eos] / f"solution_step_{step}.csv", "water_saturation")
            image = axes[row_index, col_index].imshow(
                field,
                origin="lower",
                extent=(0, LX, 0, LY),
                aspect="auto",
                interpolation="nearest",
                cmap="cividis",
                vmin=0.0,
                vmax=1.0,
            )
            axes[row_index, col_index].scatter([0.01, 1.19], [0.0475, 0.0475], marker="x", s=16, color="white", linewidth=0.8)
            axes[row_index, col_index].set_title(f"{pvi:.2f} PVI")
            if col_index == 0:
                axes[row_index, col_index].set_ylabel(f"{eos}\ny (m)")
            if row_index == 1:
                axes[row_index, col_index].set_xlabel("x (m)")
    fig.colorbar(image, ax=axes, label="Water saturation, $S_w$", fraction=0.025, pad=0.02)
    return export(fig, output, "fig02_2d_water_saturation_pr_cpa")


def phase_overall_fields(path: Path) -> dict[str, np.ndarray]:
    rows = read_rows(path)
    output: dict[str, np.ndarray] = {}
    for component in COMPONENTS:
        index = COMPONENTS.index(component)
        column = f"z_{index}_{component}"
        values = np.empty(NX * NY)
        for row in rows:
            values[int(row["input_index"])] = number(row, column)
        output[component] = values.reshape(NY, NX)
    denominator = sum(output[c] for c in COMPONENTS[1:])
    output["hc_squalane_fraction"] = np.divide(
        output["Heavy_squalane"], denominator, out=np.zeros_like(denominator), where=denominator > 0
    )
    return output


def figure_composition_maps(runs: dict[str, Path], output: Path) -> list[str]:
    fields = {eos: phase_overall_fields(path / "phase_state_step_400.csv") for eos, path in runs.items()}
    definitions = (
        ("H2O", "Overall H₂O (linear, 0.90–1.00)", mpl.colors.Normalize(0.90, 1.00), "cividis"),
        ("Light_nC4", "Overall nC₄ (log)", mpl.colors.LogNorm(1e-10, 3e-2), "viridis"),
        ("Middle_nC10", "Overall nC₁₀ (log)", mpl.colors.LogNorm(1e-10, 1.5e-1), "viridis"),
        ("hc_squalane_fraction", "Squalane / total HC (linear)", mpl.colors.Normalize(0, 1), "magma"),
    )
    fig, axes = plt.subplots(2, 4, figsize=(7.09, 3.2), layout="constrained", sharex=True, sharey=True)
    for col, (key, title, norm, cmap) in enumerate(definitions):
        images = []
        for row, eos in enumerate(("PR", "CPA")):
            image = axes[row, col].imshow(fields[eos][key], origin="lower", extent=(0, LX, 0, LY),
                                          aspect="auto", interpolation="nearest", norm=norm, cmap=cmap)
            images.append(image)
            axes[row, col].scatter([0.01, 1.19], [0.0475, 0.0475], marker="x", s=14, color="white", linewidth=0.7)
            if col == 0:
                axes[row, col].set_ylabel(f"{eos}\ny (m)")
            if row == 1:
                axes[row, col].set_xlabel("x (m)")
            if row == 0:
                axes[row, col].set_title(title, fontsize=8)
        fig.colorbar(images[-1], ax=axes[:, col], fraction=0.04, pad=0.015)
    fig.suptitle("Overall-composition fields at 1.00 PVI", fontsize=10)
    return export(fig, output, "fig03_2d_composition_at_1pvi")


def diagnostics(path: Path, system: str, eos: str) -> list[dict[str, object]]:
    result = []
    for row in read_rows(path / "reservoir_diagnostics.csv"):
        record: dict[str, object] = {"system": system, "eos": eos, "pvi": number(row, "time") * PVI_PER_DAY}
        for phase in ("o", "w"):
            for stat in ("min", "avg", "max"):
                record[f"mu_{phase}_{stat}_mPa_s"] = 1000 * number(row, f"mu_{phase}_{stat}")
        result.append(record)
    return result


def figure_viscosity_envelopes(all_runs: dict[tuple[str, str], Path], output: Path) -> tuple[list[str], list[dict[str, object]]]:
    data = {key: diagnostics(path, *key) for key, path in all_runs.items()}
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.2), layout="constrained", sharex=True)
    for row_index, system in enumerate(("binary", "lmh")):
        for col_index, phase in enumerate(("o", "w")):
            ax = axes[row_index, col_index]
            for eos in ("PR", "SW", "CPA"):
                rows = data.get((system, eos))
                if not rows:
                    continue
                valid = [r for r in rows if float(r[f"mu_{phase}_avg_mPa_s"]) > 0]
                if not valid:
                    continue
                x = np.array([float(r["pvi"]) for r in valid])
                avg = np.array([float(r[f"mu_{phase}_avg_mPa_s"]) for r in valid])
                lo = np.array([float(r[f"mu_{phase}_min_mPa_s"]) for r in valid])
                hi = np.array([float(r[f"mu_{phase}_max_mPa_s"]) for r in valid])
                color = "#1B4F72" if phase == "o" else "#0072B2"
                ax.fill_between(x, np.maximum(lo, 1e-8), np.maximum(hi, 1e-8), color=color, alpha=0.08)
                ax.plot(x, avg, color=color, linestyle=EOS_STYLES[eos], label=eos)
            ax.set(yscale="log", xlim=(0, 1), xlabel="Injected pore volume (PVI)")
            ax.set_ylabel(("Hydrocarbon-rich" if phase == "o" else "Water-rich") + " viscosity (mPa·s; log)")
            ax.set_title(("H₂O–squalane" if system == "binary" else "H₂O–nC₄–nC₁₀–squalane") + ("; oil slot" if phase == "o" else "; water slot"))
            ax.grid(True, which="both", color="#D9D9D9", linewidth=0.55)
            ax.legend(title="EOS", ncol=3)
            panel(ax, f"({chr(97 + row_index*2 + col_index)})")
    flat = [row for rows in data.values() for row in rows]
    return export(fig, output, "fig04_oil_water_viscosity_model_spread"), flat


def producer_long(path: Path, system: str, eos: str) -> list[dict[str, object]]:
    rows = [r for r in read_rows(path / "well_history.csv") if r.get("name") == "PROD" and int(r["step"]) > 0]
    columns = [key for key in rows[0] if key.startswith("m_component_")]
    result = []
    for row in rows:
        rates = {col.removeprefix("m_component_"): max(0.0, -number(row, col)) for col in columns}
        total = sum(rates.values())
        hc_total = sum(value for key, value in rates.items() if key != "H2O")
        qo = abs(number(row, "q_oil_reservoir"))
        qg = abs(number(row, "q_gas_reservoir"))
        qw = abs(number(row, "q_water_reservoir"))
        qtotal = qo + qg + qw
        h2o_mass_fraction = rates.get("H2O", 0.0) / total if total > 0 else 0.0
        for component, rate in rates.items():
            result.append(
                {
                    "system": system,
                    "eos": eos,
                    "step": int(row["step"]),
                    "pvi": number(row, "time") * PVI_PER_DAY,
                    "component": component,
                    "mass_rate_kg_s": rate,
                    "total_mass_fraction": rate / total if total > 0 else 0.0,
                    "hydrocarbon_normalized_mass_fraction": (rate / hc_total if component != "H2O" and hc_total > 0 else 0.0),
                    "producer_h2o_mass_fraction": h2o_mass_fraction,
                    "reservoir_watercut": qw / qtotal if qtotal > 0 else 0.0,
                }
            )
    return result


def plot_component_lines(ax: mpl.axes.Axes, rows: list[dict[str, object]], eos_values: tuple[str, ...], normalized: bool) -> None:
    value_key = "hydrocarbon_normalized_mass_fraction" if normalized else "total_mass_fraction"
    components = COMPONENTS[1:] if normalized else COMPONENTS
    for component in components:
        for eos in eos_values:
            selected = sorted(
                (r for r in rows if r["component"] == component and r["eos"] == eos),
                key=lambda r: float(r["pvi"]),
            )
            if not selected:
                continue
            ax.plot(
                [float(r["pvi"]) for r in selected],
                [float(r[value_key]) for r in selected],
                color=COMPONENT_COLORS[component],
                linestyle=EOS_STYLES[eos],
                marker=EOS_MARKERS[eos],
                markevery=50,
                markersize=2.7,
                label=f"{label_component(component)}–{eos}",
            )


def figure_producer_composition(all_runs: dict[tuple[str, str], Path], output: Path) -> tuple[list[str], list[dict[str, object]]]:
    all_rows = [producer_long(path, *key) for key, path in all_runs.items()]
    flat = [r for rows in all_rows for r in rows]
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.4), layout="constrained", sharex=True)
    binary = [r for r in flat if r["system"] == "binary"]
    lmh = [r for r in flat if r["system"] == "lmh"]
    plot_component_lines(axes[0, 0], binary, ("PR", "SW", "CPA"), normalized=False)
    plot_component_lines(axes[0, 1], lmh, ("PR", "CPA"), normalized=False)
    plot_component_lines(axes[1, 0], lmh, ("PR",), normalized=True)
    plot_component_lines(axes[1, 1], lmh, ("CPA",), normalized=True)
    titles = (
        "Binary: total producer composition",
        "Four-component: total producer composition",
        "Four-component PR: HC-normalized",
        "Four-component CPA: HC-normalized",
    )
    for i, ax in enumerate(axes.flat):
        ax.set(xlabel="Injected pore volume (PVI)", ylabel="Mass fraction", xlim=(0, 1), ylim=(0, 1))
        ax.grid(True, color="#D9D9D9", linewidth=0.55)
        ax.legend(fontsize=6.1, ncol=2)
        ax.set_title(f"({chr(97+i)})  {titles[i]}", loc="left", fontsize=8.5)
    return export(fig, output, "fig05_producer_total_and_hc_composition"), flat


def recovery_long(path: Path, system: str, eos: str) -> list[dict[str, object]]:
    output = []
    for row in read_rows(path / "component_mass_balance.csv"):
        initial = number(row, "initial_inventory_kg")
        produced = number(row, "cumulative_produced_kg")
        output.append(
            {
                "system": system,
                "eos": eos,
                "step": int(row["step"]),
                "pvi": number(row, "time_day") * PVI_PER_DAY,
                "component": row["component"],
                "cumulative_produced_kg": produced,
                "recovery_fraction_of_initial": produced / initial if initial > 0 else 0.0,
                "relative_balance_error": number(row, "relative_error"),
            }
        )
    return output


def figure_recovery_balance(all_runs: dict[tuple[str, str], Path], output: Path) -> tuple[list[str], list[dict[str, object]]]:
    flat = [r for key, path in all_runs.items() for r in recovery_long(path, *key)]
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.3), layout="constrained", sharex=True)
    for col, system in enumerate(("binary", "lmh")):
        eos_values = ("PR", "SW", "CPA") if system == "binary" else ("PR", "CPA")
        ax = axes[0, col]
        for component in COMPONENTS[1:]:
            for eos in eos_values:
                rows = sorted((r for r in flat if r["system"] == system and r["eos"] == eos and r["component"] == component), key=lambda r: float(r["pvi"]))
                if rows:
                    ax.plot([float(r["pvi"]) for r in rows], [float(r["recovery_fraction_of_initial"]) for r in rows],
                            color=COMPONENT_COLORS[component], linestyle=EOS_STYLES[eos], label=f"{label_component(component)}–{eos}")
        ax.set(title=("H₂O–squalane" if system == "binary" else "Four-component") + ": component recovery",
               ylabel="Cumulative recovery / initial inventory", ylim=(0, 1))
        ax.legend(fontsize=6.2, ncol=2)

        balance = axes[1, col]
        for eos in eos_values:
            grouped: dict[float, list[float]] = defaultdict(list)
            for row in flat:
                if row["system"] == system and row["eos"] == eos:
                    grouped[float(row["pvi"])].append(abs(float(row["relative_balance_error"])))
            points = sorted((pvi, max(vals)) for pvi, vals in grouped.items() if max(vals) > 0)
            balance.semilogy([p[0] for p in points], [p[1] for p in points], color=SYSTEM_COLORS[system],
                             linestyle=EOS_STYLES[eos], label=eos)
        balance.set(title="Maximum component balance error", ylabel="Absolute relative error")
        balance.legend(title="EOS", ncol=3)
    for i, ax in enumerate(axes.flat):
        ax.set(xlabel="Injected pore volume (PVI)", xlim=(0, 1))
        ax.grid(True, which="both", color="#D9D9D9", linewidth=0.55)
        panel(ax, f"({chr(97+i)})")
    return export(fig, output, "fig06_cumulative_recovery_and_mass_balance"), flat


def figure_h2o_bhp(producer: list[dict[str, object]], all_runs: dict[tuple[str, str], Path], output: Path) -> tuple[list[str], list[dict[str, object]]]:
    # One record per time is sufficient; component-long rows repeat the metric.
    # The component mass fraction is the primary breakthrough diagnostic because
    # it is invariant to oil/gas/water phase-slot reassignment.
    seen = set()
    h2o_history = []
    for row in producer:
        key = (row["system"], row["eos"], row["step"])
        if key not in seen:
            seen.add(key)
            h2o_history.append({k: row[k] for k in ("system", "eos", "step", "pvi", "producer_h2o_mass_fraction", "reservoir_watercut")})
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.1), layout="constrained")
    for system, title in (("binary", "H₂O–squalane"), ("lmh", "Four-component")):
        for eos in (("PR", "SW", "CPA") if system == "binary" else ("PR", "CPA")):
            rows = sorted((r for r in h2o_history if r["system"] == system and r["eos"] == eos), key=lambda r: float(r["pvi"]))
            axes[0].plot([float(r["pvi"]) for r in rows], [float(r["producer_h2o_mass_fraction"]) for r in rows],
                         color=SYSTEM_COLORS[system], linestyle=EOS_STYLES[eos], label=f"{title}–{eos}")
            raw = [r for r in read_rows(all_runs[(system, eos)] / "well_history.csv") if r.get("name") == "PROD" and int(r["step"]) > 0]
            axes[1].plot([number(r, "time") * PVI_PER_DAY for r in raw], [number(r, "bhp") for r in raw],
                         color=SYSTEM_COLORS[system], linestyle=EOS_STYLES[eos], label=f"{title}–{eos}")
    axes[0].set(ylabel="Produced H₂O mass fraction", ylim=(0, 1), title="Component-based water breakthrough")
    axes[1].set(ylabel="Producer BHP (bar)", title="Producer pressure response")
    for i, ax in enumerate(axes):
        ax.set(xlabel="Injected pore volume (PVI)", xlim=(0, 1))
        ax.grid(True, color="#D9D9D9", linewidth=0.55)
        ax.legend(fontsize=6.0)
        panel(ax, f"({chr(97+i)})")
    return export(fig, output, "fig07_producer_h2o_and_bhp"), h2o_history


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    output = args.output
    if output.exists() and not args.force and any(output.iterdir()):
        raise FileExistsError("output directory is not empty; use --force")
    output.mkdir(parents=True, exist_ok=True)
    root = args.root.resolve()

    calibration = root / "tools/example/scw_binary_calibration/results"
    binary = root / "case/scw_kerogen_squalane_2d/results/current_20260916"
    lmh = root / "case/scw_kerogen_lmh_2d/results/current_20260916"
    all_runs = {
        ("binary", "PR"): binary / "pr_np1_lu",
        ("binary", "SW"): binary / "sw_np1_lu",
        ("binary", "CPA"): binary / "cpa_np1_lu",
        ("lmh", "PR"): lmh / "pr_np1_lu",
        ("lmh", "CPA"): lmh / "cpa_np1_lu",
    }
    lmh_maps = {"PR": lmh / "pr_np1_lu", "CPA": lmh / "cpa_np1_lu"}

    configure_style()
    figures = []
    figures += figure_calibration(calibration, output)
    figures += figure_saturation_maps(lmh_maps, output)
    figures += figure_composition_maps(lmh_maps, output)
    names, viscosity = figure_viscosity_envelopes(all_runs, output)
    figures += names
    names, producer = figure_producer_composition(all_runs, output)
    figures += names
    names, recovery = figure_recovery_balance(all_runs, output)
    figures += names
    names, h2o_history = figure_h2o_bhp(producer, all_runs, output)
    figures += names

    write_rows(output / "viscosity_envelope_long.csv", viscosity)
    write_rows(output / "producer_composition_long.csv", producer)
    write_rows(output / "recovery_balance_long.csv", recovery)
    write_rows(output / "producer_h2o_breakthrough_long.csv", h2o_history)
    manifest = {
        "generated_on": date.today().isoformat(),
        "status": "provisional internal figures; current uncalibrated flow models",
        "python": platform.python_version(),
        "matplotlib": mpl.__version__,
        "source_runs": {f"{system}_{eos}": str(path) for (system, eos), path in all_runs.items()},
        "calibration_source": str(calibration),
        "figures": figures,
        "transformations": [
            "PVI = time_day * 0.25 PV/day",
            "2-D input_index reshaped as index = i + 60*j",
            "all saturation maps use a common linear 0–1 color scale without interpolation",
            "overall nC4 and nC10 maps use explicitly labelled logarithmic color scales",
            "producer mass fractions use absolute produced component mass rates",
            "hydrocarbon-normalized fractions exclude H2O from the denominator",
            "viscosity bands are reported grid-cell min/max from reservoir_diagnostics.csv; they are not uncertainty intervals",
            "component recovery = cumulative produced mass / initial component inventory",
            "strict zero balance errors are retained in CSV and omitted only from logarithmic display",
        ],
        "unavailable_requested_output": {
            "local_2d_viscosity_maps": "cellwise viscosity is not present in solution snapshots; only min/avg/max histories are available",
            "lmh_sw_full_run": "four-component SW failed during initial phase transition and therefore is not drawn as a completed-flow result",
        },
        "uncertainty": "none estimated; EOS curves and min/max spatial bands are model/spatial spread, not confidence intervals",
    }
    (output / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {len(figures)} figure files and four source tables to {output}")


if __name__ == "__main__":
    main()
