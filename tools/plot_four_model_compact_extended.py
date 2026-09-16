#!/usr/bin/env python3
"""Create compact extended comparisons for the four H2O-CO2-nC10 models."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize, PowerNorm, TwoSlopeNorm
from matplotlib.lines import Line2D
from PIL import Image

from plot_four_model_common_scale import (
    DEFAULT_SOURCE_ROOT,
    LX_M,
    LY_M,
    NX,
    NY,
    RUN_DIRECTORIES,
    self_color_colormap,
)


REFERENCE_MODEL = "New-PR"
MODELS = tuple(RUN_DIRECTORIES)
MODEL_COLORS = {
    "Traditional": "#111111",
    "New-PR": "#0072B2",
    "New-SW": "#009E73",
    "New-CPA": "#D55E00",
}
MODEL_MARKERS = {
    "Traditional": "o",
    "New-PR": "s",
    "New-SW": "^",
    "New-CPA": "D",
}
MODEL_LINESTYLES = {
    "Traditional": "-",
    "New-PR": "--",
    "New-SW": "-.",
    "New-CPA": ":",
}
TEXT_COLOR = "#233746"
AXIS_COLOR = "#7890A0"
GRID_COLOR = "#DCE6EC"
PANEL_BACKGROUND = "#FAFCFD"
DX_M = LX_M / NX
DY_M = LY_M / NY
WELL_ROW = 9
WELL_Y_M = (WELL_ROW + 0.5) * DY_M
PHASE_PRESENT_THRESHOLD = 1.0e-8
RASTER_DPI = 600
PRESENTATION_MODE = False
CORE_PROFILE_COLORS = {
    "Traditional": "#0072B2",
    "New-PR": "#D55E00",
    "New-SW": "#009E73",
    "New-CPA": "#CC79A7",
}
CORE_PROFILE_LINESTYLES = {
    "Traditional": "-",
    "New-PR": (0.0, (5.5, 2.2)),
    "New-SW": (0.0, (4.0, 1.6, 1.0, 1.6)),
    "New-CPA": (0.0, (1.0, 1.5)),
}
CORE_PROFILE_ZOOM_WINDOWS = {
    "pressure_mpa": (260.0, 300.0),
    "oil_saturation": (45.0, 70.0),
    "gas_saturation": (45.0, 70.0),
    "water_saturation": (45.0, 70.0),
}


def srgb_relative_luminance(rgb: np.ndarray | tuple[float, float, float]) -> float:
    channels = np.asarray(rgb, dtype=float)
    linear = np.where(
        channels <= 0.04045,
        channels / 12.92,
        ((channels + 0.055) / 1.055) ** 2.4,
    )
    return float(np.dot(linear, np.array([0.2126, 0.7152, 0.0722])))


def contrast_ratio(
    foreground: np.ndarray | tuple[float, float, float],
    background: np.ndarray | tuple[float, float, float],
) -> float:
    light = max(srgb_relative_luminance(foreground), srgb_relative_luminance(background))
    dark = min(srgb_relative_luminance(foreground), srgb_relative_luminance(background))
    return (light + 0.05) / (dark + 0.05)


def heatmap_text_color(cmap: ListedColormap, norm: Normalize, value: float) -> str:
    background = np.asarray(cmap(norm(value)))[:3]
    white = np.array([1.0, 1.0, 1.0])
    dark_text = np.asarray(mpl.colors.to_rgb(TEXT_COLOR))
    return "white" if contrast_ratio(white, background) > contrast_ratio(dark_text, background) else TEXT_COLOR


@dataclass(frozen=True)
class FieldSpec:
    key: str
    label: str
    unit: str
    filename: str
    phase: str | None = None
    difference_scale: float = 1.0
    difference_unit: str | None = None
    ppt_label: str | None = None


FIELDS = (
    FieldSpec("pressure_mpa", "Pressure", "MPa", "01_pressure", ppt_label=r"$p$ (MPa)"),
    FieldSpec("pressure_gradient_kpa_m", "Pressure-gradient magnitude", "kPa/m", "02_pressure_gradient", ppt_label=r"$|\nabla p|$ (kPa/m)"),
    FieldSpec("oil_saturation", "Oil saturation", "fraction", "03_oil_saturation", difference_scale=100.0, difference_unit="percentage points", ppt_label=r"$S_\mathrm{o}$"),
    FieldSpec("gas_saturation", "Gas saturation", "fraction", "04_gas_saturation", difference_scale=100.0, difference_unit="percentage points", ppt_label=r"$S_\mathrm{g}$"),
    FieldSpec("water_saturation", "Water saturation", "fraction", "05_water_saturation", difference_scale=100.0, difference_unit="percentage points", ppt_label=r"$S_\mathrm{w}$"),
    FieldSpec("oil_relperm", "Oil relative permeability", "fraction", "06_oil_relperm", ppt_label=r"$k_{r\mathrm{o}}$"),
    FieldSpec("gas_relperm", "Gas relative permeability", "fraction", "07_gas_relperm", ppt_label=r"$k_{r\mathrm{g}}$"),
    FieldSpec("water_relperm", "Water relative permeability", "fraction", "08_water_relperm", ppt_label=r"$k_{r\mathrm{w}}$"),
    FieldSpec("oil_x_h2o", "Oil-phase H2O mole fraction", "fraction", "09_oil_x_h2o", phase="oil", ppt_label=r"$x^\mathrm{o}_{\mathrm{H_2O}}$"),
    FieldSpec("oil_x_co2", "Oil-phase CO2 mole fraction", "fraction", "10_oil_x_co2", phase="oil", difference_scale=100.0, difference_unit="percentage points", ppt_label=r"$x^\mathrm{o}_{\mathrm{CO_2}}$"),
    FieldSpec("oil_x_nc10", "Oil-phase nC10 mole fraction", "fraction", "11_oil_x_nc10", phase="oil", difference_scale=100.0, difference_unit="percentage points", ppt_label=r"$x^\mathrm{o}_{\mathrm{nC}_{10}}$"),
    FieldSpec("gas_y_h2o", "Gas-phase H2O mole fraction", "fraction", "12_gas_y_h2o", phase="gas", ppt_label=r"$y^\mathrm{g}_{\mathrm{H_2O}}$"),
    FieldSpec("gas_y_co2", "Gas-phase CO2 mole fraction", "fraction", "13_gas_y_co2", phase="gas", difference_scale=100.0, difference_unit="percentage points", ppt_label=r"$y^\mathrm{g}_{\mathrm{CO_2}}$"),
    FieldSpec("gas_y_nc10", "Gas-phase nC10 mole fraction", "fraction", "14_gas_y_nc10", phase="gas", difference_scale=100.0, difference_unit="percentage points", ppt_label=r"$y^\mathrm{g}_{\mathrm{nC}_{10}}$"),
    FieldSpec("water_x_h2o", "Aqueous H2O mole fraction", "fraction", "15_water_x_h2o", phase="water", ppt_label=r"$x^\mathrm{w}_{\mathrm{H_2O}}$"),
    FieldSpec("water_x_co2", "Aqueous CO2 mole fraction", "fraction", "16_water_x_co2", phase="water", difference_scale=1.0e6, difference_unit="10^-6 mole fraction", ppt_label=r"$x^\mathrm{w}_{\mathrm{CO_2}}$"),
)
FIELD_BY_KEY = {field.key: field for field in FIELDS}

DIFFERENCE_FIELDS = (
    "pressure_mpa",
    "pressure_gradient_kpa_m",
    "oil_saturation",
    "gas_saturation",
    "water_saturation",
    "gas_relperm",
    "oil_x_co2",
    "water_x_co2",
)
FRONT_FOCUSED_DIFFERENCE_FIELDS = {
    "oil_saturation",
    "gas_saturation",
    "gas_relperm",
    "water_x_co2",
}

FULL_DOMAIN_ABSOLUTE_FIELDS = {
    "pressure_mpa",
    "pressure_gradient_kpa_m",
}
PRESENTATION_FOCUS_X_MAX_M = 100.0
ABSOLUTE_POWER_GAMMA = {
    "gas_saturation": 0.65,
    "oil_relperm": 0.75,
    "gas_relperm": 0.60,
    "oil_x_h2o": 0.60,
    "oil_x_co2": 0.65,
    "gas_y_h2o": 0.60,
    "gas_y_nc10": 0.60,
    "water_x_co2": 0.45,
}
ROBUST_DIFFERENCE_FIELDS = {
    "oil_saturation",
    "gas_saturation",
    "gas_relperm",
}
ROBUST_DIFFERENCE_QUANTILE = 0.99


PPT_DIFFERENCE_LABELS = {
    "pressure_mpa": r"$\Delta p$ vs PR (MPa)",
    "pressure_gradient_kpa_m": r"$\Delta|\nabla p|$ vs PR (kPa/m)",
    "oil_saturation": r"$\Delta S_\mathrm{o}$ vs PR (pp)",
    "gas_saturation": r"$\Delta S_\mathrm{g}$ vs PR (pp)",
    "water_saturation": r"$\Delta S_\mathrm{w}$ vs PR (pp)",
    "gas_relperm": r"$\Delta k_{r\mathrm{g}}$ vs PR",
    "oil_x_co2": r"$\Delta x^\mathrm{o}_{\mathrm{CO_2}}$ vs PR (pp)",
    "water_x_co2": r"$\Delta x^\mathrm{w}_{\mathrm{CO_2}}$ vs PR ($10^{-6}$)",
}


def field_label(field: FieldSpec) -> str:
    if PRESENTATION_MODE and field.ppt_label is not None:
        return field.ppt_label
    return field.label


def panel_title(panel: str, label: str) -> str:
    return label if PRESENTATION_MODE else f"({panel})  {label}"


def configure_matplotlib(presentation: bool = False) -> None:
    base_size = 12.5 if presentation else 9.5
    title_size = 13.5 if presentation else 10.5
    tick_size = 11.0 if presentation else 8.5
    legend_size = 11.0 if presentation else 8.5
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "mathtext.rm": "STIXGeneral",
            "mathtext.it": "STIXGeneral:italic",
            "mathtext.bf": "STIXGeneral:bold",
            "mathtext.default": "it",
            "font.size": base_size,
            "axes.titlesize": title_size,
            "axes.titleweight": "normal",
            "axes.titlecolor": TEXT_COLOR,
            "axes.labelsize": base_size,
            "axes.labelcolor": TEXT_COLOR,
            "axes.edgecolor": AXIS_COLOR,
            "axes.linewidth": 0.85,
            "xtick.labelsize": tick_size,
            "ytick.labelsize": tick_size,
            "xtick.color": TEXT_COLOR,
            "ytick.color": TEXT_COLOR,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 3.6,
            "ytick.major.size": 3.6,
            "xtick.major.width": 0.85,
            "ytick.major.width": 0.85,
            "legend.fontsize": legend_size,
            "legend.handlelength": 2.4,
            "legend.columnspacing": 1.2,
            "figure.dpi": 180,
            "savefig.dpi": RASTER_DPI,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.pad_inches": 0.04,
        }
    )


def corey_water(sw: np.ndarray) -> np.ndarray:
    result = np.zeros_like(sw)
    middle = (sw > 0.10) & (sw < 0.85)
    result[middle] = 0.30 * ((sw[middle] - 0.10) / 0.75) ** 2
    result[sw >= 0.85] = 0.30
    return result


def corey_gas(sg: np.ndarray) -> np.ndarray:
    result = np.zeros_like(sg)
    middle = (sg > 0.02) & (sg < 0.80)
    result[middle] = 0.80 * ((sg[middle] - 0.02) / 0.78) ** 2
    result[sg >= 0.80] = 0.80
    return result


def corey_oil(so: np.ndarray) -> np.ndarray:
    result = np.zeros_like(so)
    middle = (so > 0.10) & (so < 0.90)
    result[middle] = 0.80 * ((so[middle] - 0.10) / 0.80) ** 2
    result[so >= 0.90] = 0.80
    return result


def enrich_frame(model: str, frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.sort_values("input_index").reset_index(drop=True).copy()
    if len(frame) != NX * NY or not np.array_equal(frame["input_index"], np.arange(NX * NY)):
        raise ValueError(f"{model}: incomplete or unordered 60 x 20 field")

    frame["x_m"] = (frame["input_index"] % NX + 0.5) * DX_M
    frame["y_m"] = (frame["input_index"] // NX + 0.5) * DY_M
    frame["pressure_mpa"] = frame["pressure_Pa"] * 1.0e-6
    pressure_grid = frame["pressure_Pa"].to_numpy(dtype=float).reshape(NY, NX)
    gradient_y, gradient_x = np.gradient(pressure_grid, DY_M, DX_M)
    frame["pressure_gradient_kpa_m"] = np.hypot(gradient_x, gradient_y).ravel() * 1.0e-3

    is_traditional = "liquid_saturation" in frame
    if is_traditional:
        frame["oil_saturation"] = frame["liquid_saturation"]
        frame["gas_saturation"] = frame["vapor_saturation"]
        frame["oil_x_h2o"] = 0.0
        frame["oil_x_co2"] = frame["liquid_x_0_CO2"]
        frame["oil_x_nc10"] = 1.0 - frame["oil_x_co2"]
        frame["gas_y_h2o"] = 0.0
        frame["gas_y_co2"] = frame["vapor_y_0_CO2"]
        frame["gas_y_nc10"] = 1.0 - frame["gas_y_co2"]
        frame["water_x_h2o"] = 1.0
        frame["water_x_co2"] = 0.0
    else:
        frame["oil_x_h2o"] = frame["oil_x_0_H2O"]
        frame["oil_x_co2"] = frame["oil_x_1_CO2"]
        frame["oil_x_nc10"] = 1.0 - frame["oil_x_h2o"] - frame["oil_x_co2"]
        frame["gas_y_h2o"] = frame["gas_y_0_H2O"]
        frame["gas_y_co2"] = frame["gas_y_1_CO2"]
        frame["gas_y_nc10"] = 1.0 - frame["gas_y_h2o"] - frame["gas_y_co2"]
        frame["water_x_h2o"] = frame["water_x_0_H2O"]
        frame["water_x_co2"] = frame["water_x_1_CO2"]

    frame["oil_relperm"] = corey_oil(frame["oil_saturation"].to_numpy(dtype=float))
    frame["gas_relperm"] = corey_gas(frame["gas_saturation"].to_numpy(dtype=float))
    frame["water_relperm"] = corey_water(frame["water_saturation"].to_numpy(dtype=float))

    saturation_sum = frame[["oil_saturation", "gas_saturation", "water_saturation"]].sum(axis=1)
    if float(np.max(np.abs(saturation_sum - 1.0))) > 5.0e-8:
        raise ValueError(f"{model}: saturation closure exceeds tolerance")
    if not np.isfinite(frame[[field.key for field in FIELDS]].to_numpy(dtype=float)).all():
        raise ValueError(f"{model}: non-finite comparable quantity")
    return frame


def read_all(source_root: Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], float, dict[str, Path]]:
    data: dict[str, pd.DataFrame] = {}
    wells: dict[str, pd.DataFrame] = {}
    final_days: dict[str, float] = {}
    source_paths: dict[str, Path] = {}
    for model, directory in RUN_DIRECTORIES.items():
        run_directory = source_root / directory
        solution_path = run_directory / "solution_final.csv"
        data[model] = enrich_frame(model, pd.read_csv(solution_path))
        wells[model] = pd.read_csv(run_directory / "well_history.csv")
        summary = pd.read_csv(run_directory / "simulation_summary.csv")
        final_days[model] = float(summary["final_time_day"].iloc[-1])
        source_paths[model] = solution_path
    if max(final_days.values()) - min(final_days.values()) > 1.0e-10:
        raise ValueError(f"model final times differ: {final_days}")
    return data, wells, next(iter(final_days.values())), source_paths


def phase_mask(frame: pd.DataFrame, phase: str | None) -> np.ndarray:
    if phase is None:
        return np.ones(len(frame), dtype=bool)
    return frame[f"{phase}_saturation"].to_numpy(dtype=float) > PHASE_PRESENT_THRESHOLD


def masked_grid(frame: pd.DataFrame, field: FieldSpec) -> np.ma.MaskedArray:
    values = frame[field.key].to_numpy(dtype=float)
    valid = phase_mask(frame, field.phase)
    return np.ma.array(values.reshape(NY, NX), mask=(~valid).reshape(NY, NX))


def valid_values(frame: pd.DataFrame, field: FieldSpec) -> np.ndarray:
    values = frame[field.key].to_numpy(dtype=float)
    return values[phase_mask(frame, field.phase)]


def common_range(data: dict[str, pd.DataFrame], field: FieldSpec) -> tuple[float, float]:
    pooled = np.concatenate([valid_values(data[model], field) for model in MODELS])
    vmin, vmax = float(np.min(pooled)), float(np.max(pooled))
    if vmin == vmax:
        vmax = float(np.nextafter(vmax, np.inf))
    return vmin, vmax


def format_value(value: float) -> str:
    magnitude = abs(value)
    if magnitude == 0.0:
        return "0"
    if magnitude < 1.0e-3 or magnitude >= 1.0e4:
        return f"{value:.2e}"
    return f"{value:.3f}"


def format_showcase_value(value: float, unit: str) -> str:
    """Format values for display figures without debug-like precision."""
    if unit in {"kPa", "pp"}:
        return f"{value:.3g}"
    if abs(value) >= 100.0:
        return f"{value:.0f}"
    if abs(value) >= 10.0:
        return f"{value:.1f}"
    return f"{value:.2g}"


def add_wells(ax: plt.Axes) -> None:
    ax.scatter(
        [0.5 * DX_M], [WELL_Y_M], marker="^", s=28,
        facecolors="white", edgecolors="#2B8CBE", linewidths=1.0, zorder=4,
    )
    ax.scatter(
        [LX_M - 0.5 * DX_M], [WELL_Y_M], marker="v", s=28,
        facecolors="#173F5F", edgecolors="white", linewidths=0.7, zorder=4,
    )


def style_profile_axis(ax: plt.Axes, grid_axis: str = "both") -> None:
    """Apply a restrained manuscript style to line and dot panels."""
    ax.set_facecolor(PANEL_BACKGROUND)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(AXIS_COLOR)
    ax.spines["bottom"].set_color(AXIS_COLOR)
    ax.set_axisbelow(True)
    ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.55, alpha=0.9)


def style_showcase_axis(ax: plt.Axes, grid_axis: str = "y") -> None:
    """Apply a clean white presentation style for final report curves."""
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#9AA8B3")
    ax.spines["bottom"].set_color("#9AA8B3")
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)
    ax.set_axisbelow(True)
    ax.grid(True, axis=grid_axis, color="#E8EEF2", linewidth=0.75, alpha=0.95)


def style_box_axis(
    ax: plt.Axes,
    *,
    linewidth: float = 0.9,
    tick_labelsize: float | None = None,
) -> None:
    """Apply a grid-free four-sided frame to profile panels and insets."""
    ax.set_facecolor("white")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#6F7F8A")
        spine.set_linewidth(linewidth)
    ax.tick_params(
        direction="out",
        top=False,
        right=False,
        width=0.8,
        length=3.5,
        labelsize=tick_labelsize,
    )


def plot_model_profile(
    ax: plt.Axes,
    x: np.ndarray | pd.Series,
    y: np.ndarray | pd.Series,
    model: str,
    *,
    label: str | None = None,
) -> None:
    """Draw a model curve with color, dash and marker redundancy."""
    ax.plot(
        x,
        y,
        color=MODEL_COLORS[model],
        linestyle=MODEL_LINESTYLES[model],
        marker=MODEL_MARKERS[model],
        markerfacecolor="white",
        markeredgecolor=MODEL_COLORS[model],
        markeredgewidth=0.65,
        markersize=3.3,
        markevery=6,
        linewidth=1.55,
        label=model if label is None else label,
    )


def plot_direct_profile(
    ax: plt.Axes,
    x: np.ndarray | pd.Series,
    y: np.ndarray | pd.Series,
    model: str,
    *,
    linewidth: float = 2.55,
    alpha: float = 1.0,
    color: str | None = None,
    linestyle: str | tuple[float, tuple[float, ...]] | None = None,
) -> None:
    """Draw a cleaner direct-label curve for presentation figures."""
    ax.plot(
        x,
        y,
        color=MODEL_COLORS[model] if color is None else color,
        linestyle=MODEL_LINESTYLES[model] if linestyle is None else linestyle,
        linewidth=linewidth,
        solid_capstyle="round",
        dash_capstyle="round",
        alpha=alpha,
    )


def annotate_curve_at(
    ax: plt.Axes,
    x: np.ndarray | pd.Series,
    y: np.ndarray | pd.Series,
    model: str,
    *,
    x_target: float,
    dy_points: float = 0.0,
    dx_points: float = 6.0,
    label: str | None = None,
) -> None:
    """Place a compact model label on the nearest finite curve point."""
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    if not np.any(finite):
        return
    candidates = np.where(finite)[0]
    index = candidates[np.argmin(np.abs(x_values[candidates] - x_target))]
    ax.annotate(
        model if label is None else label,
        xy=(x_values[index], y_values[index]),
        xytext=(dx_points, dy_points),
        textcoords="offset points",
        color=MODEL_COLORS[model],
        fontsize=10.0 if PRESENTATION_MODE else 7.6,
        va="center",
        ha="left",
        clip_on=False,
    )


def save_figure(fig: plt.Figure, output: Path, stem: str) -> None:
    png = output / f"{stem}.png"
    pdf = output / f"{stem}.pdf"
    svg = output / f"{stem}.svg"
    fig.savefig(png, dpi=RASTER_DPI, facecolor="white")
    fig.savefig(pdf, dpi=RASTER_DPI, facecolor="white")
    fig.savefig(svg, dpi=RASTER_DPI, facecolor="white")
    plt.close(fig)
    with Image.open(png) as raster:
        raster.convert("RGB").save(
            png,
            dpi=(RASTER_DPI, RASTER_DPI),
            optimize=True,
        )


def plot_absolute_maps(data: dict[str, pd.DataFrame], output: Path, cmap: ListedColormap) -> list[dict[str, float | str]]:
    metrics: list[dict[str, float | str]] = []
    map_output = output / "field_maps"
    map_output.mkdir(parents=True, exist_ok=True)
    for field in FIELDS:
        vmin, vmax = common_range(data, field)
        gamma = ABSOLUTE_POWER_GAMMA.get(field.key)
        norm = (
            PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)
            if PRESENTATION_MODE and gamma is not None
            else Normalize(vmin=vmin, vmax=vmax)
        )
        front_focused = (
            PRESENTATION_MODE and field.key not in FULL_DOMAIN_ABSOLUTE_FIELDS
        )
        figure_size = (7.5, 6.15) if front_focused else (7.5, 3.05)
        fig, axes = plt.subplots(
            2,
            2,
            figsize=figure_size,
            sharex=True,
            sharey=True,
            layout="constrained",
        )
        image = None
        for column, model in enumerate(MODELS):
            ax = axes.flat[column]
            grid = masked_grid(data[model], field)
            image = ax.imshow(
                grid,
                origin="lower",
                extent=(0.0, LX_M, 0.0, LY_M),
                aspect="equal",
                interpolation="nearest",
                cmap=cmap,
                norm=norm,
            )
            values = valid_values(data[model], field)
            panel = chr(ord("a") + column)
            ax.set_title(panel_title(panel, model), pad=4, loc="left")
            if front_focused:
                ax.set_xlim(0.0, PRESENTATION_FOCUS_X_MAX_M)
                ax.set_xticks([0, 25, 50, 75, 100])
            else:
                ax.set_xticks([0, 100, 200, 300])
            ax.set_yticks([0, 50, 100])
            add_wells(ax)
            if column >= 2:
                ax.set_xlabel("x (m)", labelpad=1)
            if column % 2 == 0:
                ax.set_ylabel("y (m)", labelpad=1)
            metrics.append(
                {
                    "field": field.key,
                    "model": model,
                    "minimum": float(values.min()),
                    "maximum": float(values.max()),
                    "mean": float(values.mean()),
                    "valid_cells": int(values.size),
                    "display_gamma": gamma if PRESENTATION_MODE and gamma is not None else 1.0,
                    "display_x_max_m": (
                        PRESENTATION_FOCUS_X_MAX_M if front_focused else LX_M
                    ),
                    "display_colormap": "self_color",
                }
            )
        if field.phase is not None and not PRESENTATION_MODE:
            axes.flat[-1].text(
                0.98,
                0.05,
                "gray = phase absent",
                transform=axes.flat[-1].transAxes,
                ha="right",
                va="bottom",
                fontsize=7,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.5},
            )
        assert image is not None
        colorbar = fig.colorbar(
            image,
            ax=axes.ravel().tolist(),
            location="right",
            shrink=0.92,
            pad=0.015,
        )
        colorbar.set_label(
            (
                f"{field_label(field)}  [$\\gamma$={gamma:g}]"
                if PRESENTATION_MODE and gamma is not None
                else field_label(field)
            )
            if PRESENTATION_MODE
            else f"{field.label} ({field.unit})"
        )
        colorbar.outline.set_visible(False)
        save_figure(fig, map_output, field.filename)
    return metrics


def difference_arrays(data: dict[str, pd.DataFrame], field: FieldSpec, model: str) -> tuple[np.ma.MaskedArray, np.ndarray]:
    model_values = data[model][field.key].to_numpy(dtype=float)
    reference_values = data[REFERENCE_MODEL][field.key].to_numpy(dtype=float)
    valid = phase_mask(data[model], field.phase) & phase_mask(data[REFERENCE_MODEL], field.phase)
    difference = (model_values - reference_values) * field.difference_scale
    return np.ma.array(difference.reshape(NY, NX), mask=(~valid).reshape(NY, NX)), difference[valid]


def plot_difference_maps(data: dict[str, pd.DataFrame], output: Path, cmap: ListedColormap) -> list[dict[str, float | str]]:
    difference_output = output / "difference_maps"
    difference_output.mkdir(parents=True, exist_ok=True)
    comparisons = MODELS if PRESENTATION_MODE else ("Traditional", "New-SW", "New-CPA")
    metrics: list[dict[str, float | str]] = []
    for key in DIFFERENCE_FIELDS:
        field = FIELD_BY_KEY[key]
        front_focused = key in FRONT_FOCUSED_DIFFERENCE_FIELDS
        arrays = {model: difference_arrays(data, field, model) for model in comparisons}
        pooled_absolute = np.concatenate(
            [np.abs(values) for _, values in arrays.values()]
        )
        raw_limit = max(float(pooled_absolute.max()), np.finfo(float).eps)
        robust_scale = PRESENTATION_MODE and key in ROBUST_DIFFERENCE_FIELDS
        display_limit = (
            float(np.quantile(pooled_absolute, ROBUST_DIFFERENCE_QUANTILE))
            if robust_scale
            else raw_limit
        )
        display_limit = max(display_limit, np.finfo(float).eps)
        clipped_count = int(np.count_nonzero(pooled_absolute > display_limit))
        norm = TwoSlopeNorm(
            vmin=-display_limit,
            vcenter=0.0,
            vmax=display_limit,
        )
        if PRESENTATION_MODE:
            fig, axes = plt.subplots(
                2,
                2,
                figsize=(7.5, 6.15) if front_focused else (7.5, 3.05),
                sharex=True,
                sharey=True,
                layout="constrained",
            )
        else:
            fig, axes = plt.subplots(
                1, 3, figsize=(10.8, 2.55), sharex=True, sharey=True,
                layout="constrained",
            )
        flat_axes = np.ravel(axes)
        image = None
        for column, model in enumerate(comparisons):
            grid, values = arrays[model]
            ax = flat_axes[column]
            image = ax.imshow(
                grid,
                origin="lower",
                extent=(0.0, LX_M, 0.0, LY_M),
                aspect="equal",
                interpolation="nearest",
                cmap=cmap,
                norm=norm,
            )
            rmse = float(np.sqrt(np.mean(values * values)))
            panel = chr(ord("a") + column)
            ax.set_title(
                panel_title(panel, model if PRESENTATION_MODE else f"{model} - New-PR"),
                pad=4,
                loc="left",
            )
            if not PRESENTATION_MODE or column >= 2:
                ax.set_xlabel("x (m)", labelpad=1)
            if front_focused:
                ax.set_xlim(0.0, 100.0)
                ax.set_xticks([0, 25, 50, 75, 100])
            else:
                ax.set_xticks([0, 100, 200, 300])
            ax.set_yticks([0, 50, 100])
            add_wells(ax)
            if column == 0 or (PRESENTATION_MODE and column == 2):
                ax.set_ylabel("y (m)", labelpad=1)
            metrics.append(
                {
                    "field": field.key,
                    "model": model,
                    "difference_unit": field.difference_unit or field.unit,
                    "difference_minimum": float(values.min()),
                    "difference_maximum": float(values.max()),
                    "difference_max_abs": float(np.max(np.abs(values))),
                    "difference_rmse": rmse,
                    "difference_display_max_abs": display_limit,
                    "cells_above_display_maximum": int(
                        np.count_nonzero(np.abs(values) > display_limit)
                    ),
                    "compared_cells": int(values.size),
                }
            )
        assert image is not None
        unit = field.difference_unit or field.unit
        colorbar = fig.colorbar(
            image,
            ax=flat_axes.tolist(),
            location="right",
            shrink=0.92 if PRESENTATION_MODE else 0.88,
            pad=0.012,
            extend="both" if clipped_count else "neither",
        )
        colorbar.set_label(
            PPT_DIFFERENCE_LABELS[field.key]
            if PRESENTATION_MODE
            else f"Delta {field.label} ({unit})"
        )
        colorbar.outline.set_visible(False)
        save_figure(fig, difference_output, f"delta_{field.filename}")
    return metrics


def plot_core_centerlines(data: dict[str, pd.DataFrame], output: Path) -> None:
    fields = (
        FIELD_BY_KEY["pressure_mpa"],
        FIELD_BY_KEY["oil_saturation"],
        FIELD_BY_KEY["gas_saturation"],
        FIELD_BY_KEY["water_saturation"],
    )
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.1))
    fig.subplots_adjust(
        left=0.09,
        right=0.985,
        bottom=0.16,
        top=0.965,
        hspace=0.18,
        wspace=0.14,
    )
    for index, (ax, field) in enumerate(zip(axes.flat, fields)):
        rows = {}
        for model in MODELS:
            row = data[model].iloc[WELL_ROW * NX:(WELL_ROW + 1) * NX]
            rows[model] = row
            plot_direct_profile(
                ax,
                row["x_m"],
                row[field.key],
                model,
                linewidth=2.0,
                color=CORE_PROFILE_COLORS[model],
                linestyle=CORE_PROFILE_LINESTYLES[model],
            )
        panel = chr(ord("a") + index)
        ax.set_title(panel_title(panel, field_label(field)), loc="left")
        if PRESENTATION_MODE:
            ax.set_xlabel("x (m)" if index >= 2 else "")
        else:
            ax.set_xlabel(f"x (m), y={WELL_Y_M:g} m")
        ax.set_ylabel("" if PRESENTATION_MODE else field.unit)
        ax.set_xlim(0.0, LX_M)
        ax.set_xticks([0.0, 100.0, 200.0, 300.0])
        ax.set_box_aspect(0.62)
        style_box_axis(ax)
        if PRESENTATION_MODE and index < 2:
            ax.tick_params(labelbottom=False)

        x_min, x_max = CORE_PROFILE_ZOOM_WINDOWS[field.key]
        inset_bounds = {
            "pressure_mpa": [0.10, 0.08, 0.46, 0.40],
            "oil_saturation": [0.49, 0.08, 0.47, 0.44],
            "gas_saturation": [0.49, 0.52, 0.47, 0.44],
            "water_saturation": [0.49, 0.08, 0.47, 0.44],
        }[field.key]
        inset = ax.inset_axes(inset_bounds)
        zoom_values = []
        for model in MODELS:
            row = rows[model]
            x_values = row["x_m"].to_numpy(dtype=float)
            y_values = row[field.key].to_numpy(dtype=float)
            plot_direct_profile(
                inset,
                x_values,
                y_values,
                model,
                linewidth=1.25,
                color=CORE_PROFILE_COLORS[model],
                linestyle=CORE_PROFILE_LINESTYLES[model],
            )
            mask = (x_values >= x_min) & (x_values <= x_max)
            zoom_values.extend(y_values[mask & np.isfinite(y_values)])
        y_min = float(np.min(zoom_values))
        y_max = float(np.max(zoom_values))
        y_span = y_max - y_min
        y_pad = 0.08 * y_span if y_span > 0.0 else 0.02 * max(abs(y_min), 1.0)
        inset.set_xlim(x_min, x_max)
        inset.set_ylim(y_min - y_pad, y_max + y_pad)
        inset.set_xticks([x_min, 0.5 * (x_min + x_max), x_max])
        inset.yaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
        style_box_axis(inset, linewidth=0.75, tick_labelsize=6.3)
    handles = [
        Line2D(
            [0],
            [0],
            color=CORE_PROFILE_COLORS[model],
            linestyle=CORE_PROFILE_LINESTYLES[model],
            linewidth=2.2,
            solid_capstyle="round",
            dash_capstyle="round",
            label=model,
        )
        for model in MODELS
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.53, 0.015),
        handlelength=2.7,
    )
    save_figure(fig, output, "17_centerline_core_profiles")


def plot_core_centerline_differences(
    data: dict[str, pd.DataFrame],
    output: Path,
) -> None:
    """Plot centerline model-minus-New-PR differences with matched styling."""
    fields = (
        FIELD_BY_KEY["pressure_mpa"],
        FIELD_BY_KEY["oil_saturation"],
        FIELD_BY_KEY["gas_saturation"],
        FIELD_BY_KEY["water_saturation"],
    )
    difference_labels = {
        "pressure_mpa": r"$\Delta p$ (MPa)",
        "oil_saturation": r"$\Delta S_\mathrm{o}$",
        "gas_saturation": r"$\Delta S_\mathrm{g}$",
        "water_saturation": r"$\Delta S_\mathrm{w}$",
    }
    display_models = tuple(
        model for model in MODELS if model != REFERENCE_MODEL
    ) + (REFERENCE_MODEL,)
    reference_row = data[REFERENCE_MODEL].iloc[
        WELL_ROW * NX:(WELL_ROW + 1) * NX
    ]

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.1))
    fig.subplots_adjust(
        left=0.09,
        right=0.985,
        bottom=0.16,
        top=0.965,
        hspace=0.18,
        wspace=0.14,
    )
    for index, (ax, field) in enumerate(zip(axes.flat, fields)):
        x_values = reference_row["x_m"].to_numpy(dtype=float)
        reference_values = reference_row[field.key].to_numpy(dtype=float)
        differences = {}
        for model in display_models:
            row = data[model].iloc[WELL_ROW * NX:(WELL_ROW + 1) * NX]
            values = row[field.key].to_numpy(dtype=float)
            delta = values - reference_values
            differences[model] = delta
            plot_direct_profile(
                ax,
                x_values,
                delta,
                model,
                linewidth=2.0,
                color=CORE_PROFILE_COLORS[model],
                linestyle=CORE_PROFILE_LINESTYLES[model],
            )

        finite_values = np.concatenate([
            values[np.isfinite(values)] for values in differences.values()
        ])
        max_abs = float(np.max(np.abs(finite_values)))
        y_limit = 1.08 * max_abs if max_abs > 0.0 else 1.0
        panel = chr(ord("a") + index)
        ax.set_title(
            panel_title(panel, difference_labels[field.key]),
            loc="left",
        )
        if PRESENTATION_MODE:
            ax.set_xlabel("x (m)" if index >= 2 else "")
        else:
            ax.set_xlabel(f"x (m), y={WELL_Y_M:g} m")
        ax.set_ylabel("")
        ax.set_xlim(0.0, LX_M)
        ax.set_ylim(-y_limit, y_limit)
        ax.set_xticks([0.0, 100.0, 200.0, 300.0])
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
        ax.set_box_aspect(0.62)
        style_box_axis(ax)
        if PRESENTATION_MODE and index < 2:
            ax.tick_params(labelbottom=False)

    handles = [
        Line2D(
            [0],
            [0],
            color=CORE_PROFILE_COLORS[model],
            linestyle=CORE_PROFILE_LINESTYLES[model],
            linewidth=2.2,
            solid_capstyle="round",
            dash_capstyle="round",
            label=model,
        )
        for model in MODELS
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.53, 0.015),
        handlelength=2.7,
    )
    save_figure(fig, output, "17b_centerline_core_differences")


def plot_aqueous_co2_centerline(
    data: dict[str, pd.DataFrame],
    output: Path,
) -> None:
    """Plot the middle-line aqueous CO2 mole fraction on a portrait canvas."""
    field = FIELD_BY_KEY["water_x_co2"]
    fig, ax = plt.subplots(figsize=(3.8, 6.38))
    fig.subplots_adjust(
        left=0.20,
        right=0.965,
        bottom=0.19,
        top=0.96,
    )
    for model in MODELS:
        row = data[model].iloc[WELL_ROW * NX:(WELL_ROW + 1) * NX]
        values = row[field.key].to_numpy(dtype=float)
        values = np.where(phase_mask(row, field.phase), values, np.nan)
        plot_direct_profile(
            ax,
            row["x_m"],
            values * 1.0e3,
            model,
            linewidth=2.0,
            color=CORE_PROFILE_COLORS[model],
            linestyle=CORE_PROFILE_LINESTYLES[model],
        )

    ax.set_title(r"$x_{\mathrm{w,\,CO_2}}\;(10^{-3})$", loc="left")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("")
    ax.set_xlim(0.0, 90.0)
    ax.set_xticks([0.0, 30.0, 60.0, 90.0])
    ax.set_ylim(0.0, 13.0)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    style_box_axis(ax)
    for spine in ax.spines.values():
        spine.set_zorder(0)

    handles = [
        Line2D(
            [0],
            [0],
            color=CORE_PROFILE_COLORS[model],
            linestyle=CORE_PROFILE_LINESTYLES[model],
            linewidth=2.2,
            solid_capstyle="round",
            dash_capstyle="round",
            label=model,
        )
        for model in MODELS
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.53, 0.035),
        handlelength=2.7,
        columnspacing=1.5,
    )
    save_figure(fig, output, "17c_centerline_aqueous_co2_profile")


def plot_oil_h2o_centerline(
    data: dict[str, pd.DataFrame],
    output: Path,
) -> None:
    """Plot the middle-line oil-phase H2O mole fraction on a portrait canvas."""
    field = FIELD_BY_KEY["oil_x_h2o"]
    rows: dict[str, pd.DataFrame] = {}
    profiles: dict[str, np.ndarray] = {}
    pooled_values: list[np.ndarray] = []
    for model in MODELS:
        row = data[model].iloc[WELL_ROW * NX:(WELL_ROW + 1) * NX]
        values = row[field.key].to_numpy(dtype=float)
        values = np.where(phase_mask(row, field.phase), values, np.nan)
        rows[model] = row
        profiles[model] = values
        pooled_values.append(values[np.isfinite(values)])

    pooled_max = max(float(np.max(np.abs(values))) for values in pooled_values if values.size)
    exponent = int(np.floor(np.log10(pooled_max))) if pooled_max > 0.0 else 0
    scale = 10.0 ** (-exponent)

    fig, ax = plt.subplots(figsize=(3.8, 6.38))
    fig.subplots_adjust(
        left=0.20,
        right=0.965,
        bottom=0.19,
        top=0.96,
    )
    for model in MODELS:
        plot_direct_profile(
            ax,
            rows[model]["x_m"],
            profiles[model] * scale,
            model,
            linewidth=2.0,
            color=CORE_PROFILE_COLORS[model],
            linestyle=CORE_PROFILE_LINESTYLES[model],
        )

    ax.set_title(rf"$x_{{\mathrm{{o,\,H_2O}}}}\;(10^{{{exponent}}})$", loc="left")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("")
    ax.set_xlim(0.0, 300.0)
    ax.set_xticks([0.0, 100.0, 200.0, 300.0])
    ax.set_ylim(0.0, pooled_max * scale * 1.08)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    style_box_axis(ax)
    for spine in ax.spines.values():
        spine.set_zorder(0)

    handles = [
        Line2D(
            [0],
            [0],
            color=CORE_PROFILE_COLORS[model],
            linestyle=CORE_PROFILE_LINESTYLES[model],
            linewidth=2.2,
            solid_capstyle="round",
            dash_capstyle="round",
            label=model,
        )
        for model in MODELS
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.53, 0.035),
        handlelength=2.7,
        columnspacing=1.5,
    )
    save_figure(fig, output, "17d_centerline_oil_h2o_profile")


def plot_co2_centerlines(data: dict[str, pd.DataFrame], output: Path) -> None:
    fields = (
        (FIELD_BY_KEY["oil_x_co2"], 1.0, "mole fraction", 90.0),
        (FIELD_BY_KEY["gas_y_co2"], 1.0, "mole fraction", 60.0),
        (FIELD_BY_KEY["water_x_co2"], 1.0e6, r"mole fraction ($10^{-6}$)", 90.0),
    )
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.55))
    fig.subplots_adjust(left=0.065, right=0.988, bottom=0.18, top=0.72, wspace=0.30)
    for index, (ax, (field, scale, unit, xmax)) in enumerate(zip(axes, fields)):
        for model in MODELS:
            row = data[model].iloc[WELL_ROW * NX:(WELL_ROW + 1) * NX]
            values = row[field.key].to_numpy(dtype=float)
            valid = phase_mask(row, field.phase)
            values = np.where(valid, values, np.nan)
            plot_direct_profile(ax, row["x_m"], values * scale, model)
        panel = chr(ord("a") + index)
        ax.set_title(panel_title(panel, field_label(field)), loc="left")
        ax.set_xlabel("x (m)" if PRESENTATION_MODE else f"x (m), y={WELL_Y_M:g} m")
        ax.set_ylabel(unit)
        ax.set_xlim(0.0, xmax)
        ax.set_xticks(np.linspace(0.0, xmax, 4))
        ax.axvspan(55.0, min(70.0, xmax), color="#DCE6EC", alpha=0.26, linewidth=0.0, zorder=0)
        style_showcase_axis(ax, grid_axis="y")
    handles = [
        Line2D([0], [0], color=MODEL_COLORS[model], linestyle=MODEL_LINESTYLES[model], linewidth=2.5, label=model)
        for model in MODELS
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.52, 0.965), handlelength=2.7)
    save_figure(fig, output, "18_centerline_co2_profiles")


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(values)
    y = np.arange(1, len(x) + 1, dtype=float) / len(x)
    return x, y


def plot_ecdf(data: dict[str, pd.DataFrame], output: Path) -> None:
    fields = (
        FIELD_BY_KEY["pressure_mpa"],
        FIELD_BY_KEY["pressure_gradient_kpa_m"],
        FIELD_BY_KEY["gas_saturation"],
        FIELD_BY_KEY["water_saturation"],
    )
    fig, axes = plt.subplots(2, 2, figsize=(8.3, 5.6), layout="constrained")
    for index, (ax, field) in enumerate(zip(axes.flat, fields)):
        for model in MODELS:
            x, y = ecdf(valid_values(data[model], field))
            ax.plot(
                x, y, color=MODEL_COLORS[model],
                linestyle=MODEL_LINESTYLES[model], linewidth=1.55, label=model,
            )
        panel = chr(ord("a") + index)
        ax.set_title(panel_title(panel, field_label(field)), loc="left")
        ax.set_xlabel(field.unit)
        ax.set_ylabel("Cell fraction" if PRESENTATION_MODE else "Cumulative cell fraction")
        ax.set_ylim(0.0, 1.0)
        style_profile_axis(ax)
    axes[0, 0].legend(frameon=False, ncol=2, facecolor="white")
    save_figure(fig, output, "19_field_distributions_ecdf")


def plot_normalized_rmse_heatmap(difference_metrics: pd.DataFrame, data: dict[str, pd.DataFrame], output: Path, cmap: ListedColormap) -> None:
    comparisons = ["Traditional", "New-SW", "New-CPA"]
    rows = []
    labels = []
    for key in DIFFERENCE_FIELDS:
        field = FIELD_BY_KEY[key]
        vmin, vmax = common_range(data, field)
        scale = (vmax - vmin) * field.difference_scale
        row = []
        for model in comparisons:
            metric = difference_metrics[(difference_metrics.field == key) & (difference_metrics.model == model)].iloc[0]
            row.append(100.0 * float(metric.difference_rmse) / max(scale, np.finfo(float).eps))
        rows.append(row)
        labels.append(field_label(field))
    matrix = np.asarray(rows)
    fig, ax = plt.subplots(figsize=(6.6, 4.8), layout="constrained")
    norm = Normalize(vmin=0.0, vmax=max(float(matrix.max()), np.finfo(float).eps))
    image = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(len(comparisons)), comparisons)
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.set_xticks(np.arange(-0.5, len(comparisons), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            ax.text(
                column,
                row,
                f"{matrix[row, column]:.2f}%",
                ha="center",
                va="center",
                fontsize=7.5,
                color=heatmap_text_color(cmap, norm, float(matrix[row, column])),
            )
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("NRMSE (%)" if PRESENTATION_MODE else "RMSE / pooled field range (%)")
    colorbar.outline.set_visible(False)
    save_figure(fig, output, "20_normalized_rmse_heatmap")


def final_well_row(wells: pd.DataFrame, name: str) -> pd.Series:
    rows = wells[wells["name"] == name].sort_values("time")
    return rows.iloc[-1]


def plot_well_dots(wells: dict[str, pd.DataFrame], output: Path) -> pd.DataFrame:
    records = []
    for model in MODELS:
        injector = final_well_row(wells[model], "CO2_INJ")
        producer = final_well_row(wells[model], "PROD")
        records.append(
            {
                "model": model,
                "injector_bhp_bar": float(injector["bhp"]),
                "producer_total_reservoir_m3_s": -float(producer["q_total_reservoir"]),
                "producer_oil_reservoir_m3_s": -float(producer["q_oil_reservoir"]),
                "producer_water_reservoir_m3_s": -float(producer["q_water_reservoir"]),
            }
        )
    frame = pd.DataFrame(records)
    frame.to_csv(output / "well_final_metrics.csv", index=False)
    panels = (
        ("injector_bhp_bar", "Injector BHP", "bar", 1.0),
        ("producer_total_reservoir_m3_s", "Total rate" if PRESENTATION_MODE else "Producer total reservoir rate", r"$10^{-5}$ m³/s", 1.0e5),
        ("producer_oil_reservoir_m3_s", "Oil rate" if PRESENTATION_MODE else "Producer oil reservoir rate", r"$10^{-5}$ m³/s", 1.0e5),
        ("producer_water_reservoir_m3_s", "Water rate" if PRESENTATION_MODE else "Producer water reservoir rate", r"$10^{-7}$ m³/s", 1.0e7),
    )
    fig, axes = plt.subplots(2, 2, figsize=(8.3, 5.2), layout="constrained")
    y = np.arange(len(MODELS))
    for panel_index, (ax, (column, label, unit, scale)) in enumerate(zip(axes.flat, panels)):
        values = frame[column].to_numpy(dtype=float) * scale
        span = max(values.max() - values.min(), abs(values.mean()) * 1.0e-3, 1.0e-12)
        for model_index, model in enumerate(MODELS):
            ax.scatter(
                values[model_index], model_index,
                color=MODEL_COLORS[model], marker=MODEL_MARKERS[model],
                s=48 if PRESENTATION_MODE else 38, zorder=3,
            )
            ax.text(
                values[model_index] + 0.035 * span,
                model_index,
                f"{values[model_index]:.4g}",
                color=MODEL_COLORS[model],
                va="center",
                ha="left",
                fontsize=8.4 if PRESENTATION_MODE else 7.0,
            )
        ax.set_xlim(values.min() - 0.10 * span, values.max() + 0.26 * span)
        ax.set_yticks(y, MODELS)
        panel = chr(ord("a") + panel_index)
        ax.set_title(panel_title(panel, label), loc="left")
        ax.set_xlabel(unit)
        style_showcase_axis(ax, grid_axis="x")
    save_figure(fig, output, "21_final_well_response")
    return frame


def plot_plume_geometry(data: dict[str, pd.DataFrame], output: Path) -> pd.DataFrame:
    records = []
    cell_area = DX_M * DY_M
    for model in MODELS:
        sg = data[model]["gas_saturation"].to_numpy(dtype=float)
        x = data[model]["x_m"].to_numpy(dtype=float)
        record = {"model": model}
        for threshold in (0.01, 0.10):
            mask = sg > threshold
            record[f"area_sg_gt_{threshold:.2f}_m2"] = float(mask.sum() * cell_area)
            record[f"front_sg_gt_{threshold:.2f}_m"] = float(x[mask].max()) if np.any(mask) else 0.0
        records.append(record)
    frame = pd.DataFrame(records)
    frame.to_csv(output / "plume_geometry.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.2), layout="constrained")
    x_index = np.arange(len(MODELS))
    width = 0.34
    threshold_colors = ("#9ECAE1", "#3182BD")
    for offset, threshold in enumerate((0.01, 0.10)):
        values = frame[f"area_sg_gt_{threshold:.2f}_m2"]
        axes[0].bar(
            x_index + (offset - 0.5) * width,
            values,
            width,
            label=rf"$S_\mathrm{{g}}>{threshold:.2f}$",
            color=threshold_colors[offset],
            edgecolor="white",
            linewidth=0.5,
        )
        axes[1].plot(
            x_index,
            frame[f"front_sg_gt_{threshold:.2f}_m"],
            marker=["o", "s"][offset],
            markerfacecolor="white",
            markeredgecolor=threshold_colors[offset],
            color=threshold_colors[offset],
            linewidth=1.6,
            linestyle=["-", "--"][offset],
            label=rf"$S_\mathrm{{g}}>{threshold:.2f}$",
        )
    axes[0].set_xticks(x_index, MODELS, rotation=15)
    axes[0].set_ylabel("Area (m²)" if PRESENTATION_MODE else "Gas-plume area (m2)")
    axes[0].set_ylim(bottom=0.0)
    axes[0].set_title(panel_title("a", "Area" if PRESENTATION_MODE else "Plume area"), loc="left")
    axes[0].legend(frameon=False, facecolor="white")
    style_profile_axis(axes[0], grid_axis="y")
    axes[1].set_xticks(x_index, MODELS, rotation=15)
    axes[1].set_ylabel("Front x (m)" if PRESENTATION_MODE else "Farthest cell-center x (m)")
    front_columns = [column for column in frame if column.startswith("front_")]
    front_values = frame[front_columns].to_numpy(dtype=float)
    front_padding = max(2.5, 0.15 * float(front_values.max() - front_values.min()))
    axes[1].set_ylim(float(front_values.min()) - front_padding, float(front_values.max()) + front_padding)
    axes[1].set_title(panel_title("b", "Front" if PRESENTATION_MODE else "Front position"), loc="left")
    axes[1].legend(frameon=False, facecolor="white")
    style_profile_axis(axes[1], grid_axis="y")
    save_figure(fig, output, "22_gas_plume_geometry")
    return frame


def plot_model_spread_maps(
    data: dict[str, pd.DataFrame],
    output: Path,
    cmap: ListedColormap,
) -> pd.DataFrame:
    """Map the four-model max-minus-min spread without selecting a reference."""
    fields = (
        (FIELD_BY_KEY["pressure_mpa"], 1.0e3, "kPa"),
        (FIELD_BY_KEY["gas_saturation"], 100.0, "percentage points"),
        (FIELD_BY_KEY["oil_x_h2o"], 1.0e6, "10^-6 mole fraction"),
        (FIELD_BY_KEY["water_x_co2"], 1.0e6, "10^-6 mole fraction"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 4.35), layout="constrained")
    records: list[dict[str, float | int | str]] = []
    for index, (ax, (field, scale, unit)) in enumerate(zip(axes.flat, fields)):
        stack = np.vstack([
            data[model][field.key].to_numpy(dtype=float) for model in MODELS
        ])
        valid = np.logical_and.reduce([
            phase_mask(data[model], field.phase) for model in MODELS
        ])
        spread = (np.max(stack[:, valid], axis=0) - np.min(stack[:, valid], axis=0)) * scale
        full_spread = np.full(NX * NY, np.nan, dtype=float)
        full_spread[valid] = spread
        grid = np.ma.masked_invalid(full_spread.reshape(NY, NX))
        # A few isolated CPA phase-switch cells can be almost one saturation
        # unit apart from the other models.  Preserve those values in the
        # metrics, but keep the body of the front visible by using a clearly
        # extended, percentile-based display scale for this panel only.
        robust_scale = field.key == "gas_saturation"
        display_max = (
            float(np.quantile(spread, 0.99))
            if robust_scale
            else float(spread.max())
        )
        display_max = max(display_max, np.finfo(float).eps)
        clipped_count = int(np.count_nonzero(spread > display_max))
        spread_cmap = cmap.with_extremes(
            bad="#D9D9D9",
            over=tuple(np.asarray(cmap.colors)[-1]),
        )
        image = ax.imshow(
            grid,
            origin="lower",
            extent=(0.0, LX_M, 0.0, LY_M),
            aspect="equal",
            interpolation="nearest",
            cmap=spread_cmap,
            norm=Normalize(vmin=0.0, vmax=display_max),
        )
        panel = chr(ord("a") + index)
        ax.set_title(panel_title(panel, field_label(field)), loc="left")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_xticks([0, 100, 200, 300])
        ax.set_yticks([0, 50, 100])
        add_wells(ax)
        colorbar = fig.colorbar(
            image,
            ax=ax,
            shrink=0.86,
            pad=0.02,
            extend="max" if clipped_count else "neither",
        )
        colorbar.set_label(
            f"range ({'pp' if unit == 'percentage points' else '$10^{-6}$' if unit == '10^-6 mole fraction' else unit})"
            if PRESENTATION_MODE
            else f"max - min ({unit})"
        )
        colorbar.outline.set_visible(False)
        if robust_scale and not PRESENTATION_MODE:
            ax.text(
                0.02,
                0.03,
                f"99th-percentile scale; {clipped_count} cells above limit",
                transform=ax.transAxes,
                fontsize=6.5,
                color="#111111",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
            )
        records.append(
            {
                "field": field.key,
                "unit": unit,
                "spread_maximum": float(spread.max()),
                "spread_mean": float(spread.mean()),
                "spread_p95": float(np.quantile(spread, 0.95)),
                "display_maximum": display_max,
                "cells_above_display_maximum": clipped_count,
                "compared_cells": int(valid.sum()),
            }
        )
    save_figure(fig, output, "23_four_model_spatial_spread")
    frame = pd.DataFrame(records)
    frame.to_csv(output / "spatial_spread_metrics.csv", index=False)
    return frame


def plot_pairwise_rmse_heatmap(
    data: dict[str, pd.DataFrame],
    output: Path,
    cmap: ListedColormap,
) -> pd.DataFrame:
    """Compare all six model pairs on pooled-range-normalized field RMSE."""
    fields = (
        FIELD_BY_KEY["pressure_mpa"],
        FIELD_BY_KEY["oil_saturation"],
        FIELD_BY_KEY["gas_saturation"],
        FIELD_BY_KEY["water_saturation"],
        FIELD_BY_KEY["oil_x_h2o"],
        FIELD_BY_KEY["oil_x_co2"],
        FIELD_BY_KEY["water_x_co2"],
    )
    pairs = tuple(
        (MODELS[left], MODELS[right])
        for left in range(len(MODELS))
        for right in range(left + 1, len(MODELS))
    )
    matrix = np.zeros((len(fields), len(pairs)), dtype=float)
    records: list[dict[str, float | str | int]] = []
    for row, field in enumerate(fields):
        vmin, vmax = common_range(data, field)
        pooled_range = max(vmax - vmin, np.finfo(float).eps)
        for column, (left, right) in enumerate(pairs):
            valid = phase_mask(data[left], field.phase) & phase_mask(data[right], field.phase)
            delta = (
                data[left][field.key].to_numpy(dtype=float)
                - data[right][field.key].to_numpy(dtype=float)
            )[valid]
            rmse = float(np.sqrt(np.mean(delta * delta)))
            normalized = 100.0 * rmse / pooled_range
            matrix[row, column] = normalized
            records.append(
                {
                    "field": field.key,
                    "model_left": left,
                    "model_right": right,
                    "rmse": rmse,
                    "rmse_over_pooled_range_percent": normalized,
                    "compared_cells": int(valid.sum()),
                }
            )

    fig, ax = plt.subplots(figsize=(7.5, 5.1), layout="constrained")
    norm = Normalize(vmin=0.0, vmax=max(float(matrix.max()), np.finfo(float).eps))
    image = ax.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )
    pair_labels = [
        f"{left.replace('Traditional', 'Trad.').replace('New-', '')}–{right.replace('Traditional', 'Trad.').replace('New-', '')}"
        if PRESENTATION_MODE
        else f"{left}\nvs {right}"
        for left, right in pairs
    ]
    ax.set_xticks(np.arange(len(pairs)), pair_labels)
    ax.set_yticks(np.arange(len(fields)), [field_label(field) for field in fields])
    ax.set_xticks(np.arange(-0.5, len(pairs), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(fields), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.15)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            ax.text(
                column,
                row,
                f"{matrix[row, column]:.2f}%",
                ha="center",
                va="center",
                fontsize=7.0,
                color=heatmap_text_color(cmap, norm, float(matrix[row, column])),
            )
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("NRMSE (%)" if PRESENTATION_MODE else "Pairwise RMSE / pooled four-model range (%)")
    colorbar.outline.set_visible(False)
    save_figure(fig, output, "24_pairwise_normalized_rmse")
    frame = pd.DataFrame(records)
    frame.to_csv(output / "pairwise_rmse_metrics.csv", index=False)
    return frame


def plot_phase_exchange_profiles(
    data: dict[str, pd.DataFrame],
    output: Path,
) -> None:
    """Highlight H2O in oil and CO2 in water along the injector-producer line."""
    panels = (
        ("oil_x_h2o", 1.0e6, r"$x^\mathrm{o}_{\mathrm{H_2O}}$"),
        ("water_x_co2", 1.0e6, r"$x^\mathrm{w}_{\mathrm{CO_2}}$"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.35), sharex=True)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.16, top=0.78, wspace=0.23)
    for index, (ax, (key, scale, title)) in enumerate(zip(axes.flat, panels)):
        field = FIELD_BY_KEY[key]
        for model in MODELS:
            row = data[model].iloc[WELL_ROW * NX:(WELL_ROW + 1) * NX]
            values = row[key].to_numpy(dtype=float)
            valid = phase_mask(row, field.phase)
            values = np.where(valid, values, np.nan)
            ax.plot(
                row["x_m"],
                values * scale,
                color=MODEL_COLORS[model],
                linestyle=MODEL_LINESTYLES[model],
                linewidth=2.35 if model == "Traditional" else 2.75,
                alpha=0.72 if model == "Traditional" else 1.0,
                solid_capstyle="round",
                dash_capstyle="round",
                label=model,
            )
        ax.set_title(title, loc="left", pad=8)
        ax.set_xlabel("x (m)")
        ax.set_ylabel(r"mole fraction ($10^{-6}$)")
        ax.set_xlim(0.0, 90.0)
        ax.margins(y=0.10)
        ax.axvspan(55.0, 70.0, color="#DCE6EC", alpha=0.28, linewidth=0.0, zorder=0)
        style_showcase_axis(ax, grid_axis="y")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.965),
        handlelength=2.7,
        borderaxespad=0.0,
    )
    save_figure(fig, output, "25_phase_exchange_centerlines")


def plot_x_resolved_differences(
    data: dict[str, pd.DataFrame],
    output: Path,
) -> None:
    """Summarize how strongly each model departs from New-PR."""
    panels = (
        ("pressure_mpa", 1.0e3, r"$|\Delta p|$", "kPa"),
        ("gas_saturation", 100.0, r"$|\Delta S_\mathrm{g}|$", "pp"),
        ("oil_x_h2o", 1.0e6, r"$|\Delta x^\mathrm{o}_{\mathrm{H_2O}}|$", r"$10^{-6}$"),
        ("water_x_co2", 1.0e6, r"$|\Delta x^\mathrm{w}_{\mathrm{CO_2}}|$", r"$10^{-6}$"),
    )
    compared_models = ("Traditional", "New-SW", "New-CPA")
    y_positions = np.arange(len(compared_models))[::-1]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 5.6), layout="constrained")
    for index, (ax, (key, scale, title, unit)) in enumerate(zip(axes.flat, panels)):
        field = FIELD_BY_KEY[key]
        reference = data[REFERENCE_MODEL][key].to_numpy(dtype=float)
        reference_valid = phase_mask(data[REFERENCE_MODEL], field.phase)
        means = []
        for model in compared_models:
            values = data[model][key].to_numpy(dtype=float)
            valid = phase_mask(data[model], field.phase) & reference_valid
            delta = np.abs(values - reference)[valid]
            means.append(float(np.mean(delta) * scale))
        maximum = max(max(means), np.finfo(float).eps)
        for y, model, value in zip(y_positions, compared_models, means):
            ax.hlines(
                y,
                0.0,
                value,
                color=MODEL_COLORS[model],
                linewidth=3.0,
                alpha=0.25,
            )
            ax.scatter(
                [value],
                [y],
                s=66,
                color=MODEL_COLORS[model],
                edgecolor="white",
                linewidth=0.9,
                zorder=3,
            )
            ax.text(
                value + 0.035 * maximum,
                y,
                format_showcase_value(value, unit),
                color=MODEL_COLORS[model],
                va="center",
                ha="left",
                fontsize=9.0 if PRESENTATION_MODE else 7.2,
            )
        ax.set_title(title, loc="left")
        ax.set_yticks(y_positions, [m.replace("New-", "") for m in compared_models])
        ax.set_xlabel(unit)
        ax.set_xlim(0.0, maximum * 1.28)
        ax.set_ylim(-0.6, len(compared_models) - 0.4)
        style_showcase_axis(ax, grid_axis="x")
    save_figure(fig, output, "26_x_resolved_model_differences")


def plot_gas_front_contours(
    data: dict[str, pd.DataFrame],
    output: Path,
) -> None:
    """Overlay gas-front contours with color and line-style redundancy."""
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.12, top=0.77)
    x_edges = np.linspace(0.0, LX_M, NX + 1)
    y_edges = np.linspace(0.0, LY_M, NY + 1)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    for model in MODELS:
        grid = data[model]["gas_saturation"].to_numpy(dtype=float).reshape(NY, NX)
        for threshold, linewidth in ((0.01, 1.2), (0.10, 2.2)):
            ax.contour(
                x_centers,
                y_centers,
                grid,
                levels=[threshold],
                colors=[MODEL_COLORS[model]],
                linestyles=[MODEL_LINESTYLES[model]],
                linewidths=[linewidth],
            )
    add_wells(ax)
    ax.set_aspect("equal")
    ax.set_xlim(0.0, 90.0)
    ax.set_ylim(0.0, LY_M)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_xticks([0, 30, 60, 90])
    ax.set_yticks([0, 50, 100])
    ax.axvspan(55.0, 70.0, color="#DCE6EC", alpha=0.26, linewidth=0.0, zorder=0)
    style_showcase_axis(ax)
    model_handles = [
        Line2D(
            [0],
            [0],
            color=MODEL_COLORS[model],
            linestyle=MODEL_LINESTYLES[model],
            linewidth=1.7,
            label=model,
        )
        for model in MODELS
    ]
    threshold_handles = [
        Line2D([0], [0], color="#555555", linewidth=1.2, label=r"$S_\mathrm{g}=0.01$"),
        Line2D([0], [0], color="#555555", linewidth=2.2, label=r"$S_\mathrm{g}=0.10$"),
    ]
    fig.legend(handles=model_handles, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.55, 0.965), handlelength=2.7)
    ax.legend(handles=threshold_handles, frameon=False, loc="lower right")
    save_figure(fig, output, "27_gas_front_contours")


def write_provenance(
    output: Path,
    source_root: Path,
    source_paths: dict[str, Path],
    final_day: float,
    absolute_metrics: pd.DataFrame,
    difference_metrics: pd.DataFrame,
) -> None:
    absolute_metrics.to_csv(output / "field_metrics.csv", index=False)
    difference_metrics.to_csv(output / "difference_metrics.csv", index=False)
    manifest = {
        "source_root": str(source_root.resolve()),
        "source_files": {model: str(path.resolve()) for model, path in source_paths.items()},
        "final_time_day": final_day,
        "final_pvi": final_day / 3652.5,
        "grid": {"nx": NX, "ny": NY, "lx_m": LX_M, "ly_m": LY_M},
        "phase_present_threshold": PHASE_PRESENT_THRESHOLD,
        "transformations": [
            "Traditional liquid/vapor mapped to oil/gas",
            "Traditional oil/gas H2O mole fraction set to zero by its immiscible-water model definition",
            "Traditional aqueous phase set to pure H2O by its model definition",
            "nC10 phase mole fractions derived by closure",
            "phase compositions masked where the corresponding saturation is <= 1e-8",
            "pressure gradient from numpy.gradient on 5 m cell-center spacing: centered interior and one-sided boundary differences",
            "relative permeability derived from the exact shared Corey functions in benchmark_common.hpp",
            "field maps use nearest-neighbor pixels; no spatial smoothing or cell removal",
            "difference maps compare only cells where both model and New-PR phases are present",
            "presentation absolute maps except pressure and pressure gradient crop x to 0-100 m; the full-domain values still determine pooled color limits",
            "oil saturation, gas saturation, gas relative-permeability, and aqueous CO2 difference maps crop x to 0-100 m; all nonzero front differences lie inside this interval",
            "figure 17 centerline panels retain the full 0-300 m domain and use unsmoothed insets at x=260-300 m for pressure and x=45-70 m for saturations",
            "figure 17b centerline differences are model minus New-PR and use symmetric zero-centered y-limits over the full 0-300 m domain",
            "figure 17c shows the y=47.5 m aqueous CO2 mole-fraction profile scaled by 1e3 on a portrait canvas and crops to x=0-90 m because the remaining domain is at the numerical floor",
            "figure 17d shows the y=47.5 m oil-phase H2O mole-fraction profile on a portrait canvas; its scientific-notation exponent is selected from the pooled centerline magnitude and the full x=0-300 m domain is retained",
            "figure 18 centerline CO2-composition profiles crop to the 0-90 m front region",
            "gas-front contours crop to x=0-90 m because every 0.01 and 0.10 contour lies within that interval; contour coordinates are not smoothed",
            "presentation oil-saturation, gas-saturation, and gas-relative-permeability difference maps use a symmetric 99th-percentile display limit with explicit colorbar extensions; raw extrema are retained in difference_metrics.csv",
        ],
        "relative_permeability": {
            "water": "0 below Sw=0.10; 0.30*((Sw-0.10)/0.75)^2; capped 0.30 at Sw>=0.85",
            "gas": "0 below Sg=0.02; 0.80*((Sg-0.02)/0.78)^2; capped 0.80 at Sg>=0.80",
            "oil": "0 below So=0.10; 0.80*((So-0.10)/0.80)^2; capped 0.80 at So>=0.90",
        },
        "colormap": {
            "name": "self_color",
            "source": "user-supplied MRST self_color.m anchors embedded in plot_four_model_common_scale.py",
            "user_source_sha256": "E9111C8E7228339E9AC446690BC298A7A4A8C32076220FD712FA782B8760D51B",
            "interpolation": "PCHIP, 256 colors",
            "anchors_rgb_0_1": [
                [0.300, 0.550, 0.750],
                [0.550, 0.750, 0.850],
                [0.700, 0.800, 0.800],
                [0.950, 0.850, 0.700],
                [0.950, 0.650, 0.400],
                [0.950, 0.550, 0.300],
            ],
            "absolute_maps": "exact blue-to-orange self_color palette",
            "difference_maps": "the same exact self_color palette with zero mapped by TwoSlopeNorm",
            "spread_and_rmse_maps": "the same exact self_color palette",
            "missing_phase_color": "#D9D9D9",
        },
        "normalization": {
            "absolute_maps": "pooled minimum and maximum across all four models for each quantity; selected presentation fields use PowerNorm without clipping",
            "absolute_power_gamma": ABSOLUTE_POWER_GAMMA if PRESENTATION_MODE else {},
            "difference_maps": "symmetric about zero; selected presentation saturation/relative-permeability fields use an extended 99th-percentile display limit",
            "difference_robust_quantile": ROBUST_DIFFERENCE_QUANTILE,
            "spread_maps": "cellwise maximum minus minimum across all four models, starting at zero; gas saturation uses a colorbar extended above its 99th percentile and reports clipped-cell count",
            "pairwise_rmse": "RMSE divided by the pooled four-model field range",
        },
        "layout": {
            "presentation_mode": PRESENTATION_MODE,
            "absolute_field_maps": "2 x 2 model panels; physical x:y aspect preserved; presentation transport fields focus on x=0-100 m",
            "large_overall_titles": "none",
            "panel_titles": "panel letter plus model or quantity only; extrema and RMSE remain in CSV tables",
            "raster_dpi": RASTER_DPI,
            "formats": ["PNG", "PDF", "SVG"],
        },
        "styling": {
            "font": "Arial with DejaVu Sans fallback",
            "line_model_encoding": "color plus dash pattern plus marker shape",
            "heatmap_annotation_text": "per-cell white or navy selected by sRGB contrast ratio against the mapped self_color background",
            "profile_panel_background": PANEL_BACKGROUND,
            "grid_color": GRID_COLOR,
            "injector_marker": "open cyan upward triangle",
            "producer_marker": "filled navy downward triangle",
            "label_policy": (
                "compact symbols and units only; explanatory text moved to README/CSV"
                if PRESENTATION_MODE
                else "compact manuscript labels"
            ),
        },
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
        },
    }
    (output / "figure_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    readme_heading = (
        "# PPT-ready four-model comparison"
        if PRESENTATION_MODE
        else "# Compact extended four-model comparison"
    )
    presentation_note = (
        "\nThis export uses concise symbolic labels, larger projected-display text, "
        "and no in-panel explanatory sentences. Full definitions and numerical "
        "values remain in the CSV files and figure manifest.\n"
        if PRESENTATION_MODE
        else ""
    )
    (output / "README.md").write_text(
        readme_heading
        + """

No figure uses a large overall title. The 16 physical-field maps place the four
models in a compact 2 x 2 layout with one pooled color scale. Gray cells in phase
composition maps mean that the corresponding phase is absent, not zero.
Presentation transport maps focus on x = 0-100 m, where all active fronts occur;
pressure and pressure-gradient maps retain the full 300 m domain. Absolute maps
use the exact user-supplied MRST self_color palette. Selected wide-dynamic-range fields
use a labeled PowerNorm while retaining the same pooled four-model extrema.
"""
        + presentation_note
        + """

Additional figure types include New-PR difference maps, centerline profiles,
empirical cumulative distributions, reference-free spatial spread maps,
all-pairs normalized RMSE, H2O/CO2 phase-exchange profiles, x-resolved model
differences, gas-front contours, final-well dot plots, and gas-plume geometry.
Panel titles remain deliberately short; numerical extrema, RMSE values and
other summaries are retained in the accompanying CSV files.
The gas-saturation spread panel uses an extended 99th-percentile display scale
so the front remains visible while isolated phase-switch cells stay explicit.
The presentation oil/gas saturation and gas-relative-permeability difference
maps use the same explicit extension convention for isolated extrema.
Raw simulator CSV files are unchanged. See
`figure_manifest.json` for exact transformations and `field_metrics.csv` /
`difference_metrics.csv` for the underlying summaries.

The illustrated Chinese-language interpretation is available in
[`RESULT_ANALYSIS.md`](RESULT_ANALYSIS.md).
""",
        encoding="utf-8",
    )


def refresh_report_figures(output: Path) -> None:
    """Copy report figures to versioned names so Markdown previews refresh."""
    report_dir = output / "report_figures_10m6"
    report_dir.mkdir(parents=True, exist_ok=True)
    figure_pairs = (
        ("field_maps/01_pressure.png", "01_pressure_10m6.png"),
        ("21_final_well_response.png", "21_final_well_response_10m6.png"),
        ("field_maps/04_gas_saturation.png", "04_gas_saturation_10m6.png"),
        ("27_gas_front_contours.png", "27_gas_front_contours_10m6.png"),
        ("23_four_model_spatial_spread.png", "23_four_model_spatial_spread_10m6.png"),
        ("field_maps/09_oil_x_h2o.png", "09_oil_x_h2o_10m6.png"),
        ("field_maps/16_water_x_co2.png", "16_water_x_co2_10m6.png"),
        ("difference_maps/delta_16_water_x_co2.png", "delta_16_water_x_co2_10m6.png"),
        ("25_phase_exchange_centerlines.png", "25_phase_exchange_centerlines_10m6.png"),
        ("24_pairwise_normalized_rmse.png", "24_pairwise_normalized_rmse_10m6.png"),
        ("26_x_resolved_model_differences.png", "26_x_resolved_model_differences_10m6.png"),
    )
    for source, destination in figure_pairs:
        shutil.copy2(output / source, report_dir / destination)


def write_result_analysis(
    output: Path,
    source_root: Path,
    data: dict[str, pd.DataFrame],
    final_day: float,
) -> None:
    """Write a figure-linked, evidence-traceable analysis beside the exports."""
    fields = pd.read_csv(output / "field_metrics.csv")
    differences = pd.read_csv(output / "difference_metrics.csv")
    wells = pd.read_csv(output / "well_final_metrics.csv")
    plumes = pd.read_csv(output / "plume_geometry.csv")
    spreads = pd.read_csv(output / "spatial_spread_metrics.csv")
    pairs = pd.read_csv(output / "pairwise_rmse_metrics.csv")

    def field_value(field: str, model: str, column: str) -> float:
        row = fields[(fields["field"] == field) & (fields["model"] == model)]
        return float(row.iloc[0][column])

    def difference_value(field: str, model: str, column: str) -> float:
        row = differences[
            (differences["field"] == field) & (differences["model"] == model)
        ]
        return float(row.iloc[0][column])

    def spread_value(field: str, column: str) -> float:
        row = spreads[spreads["field"] == field]
        return float(row.iloc[0][column])

    def pair_value(field: str, left: str, right: str) -> float:
        row = pairs[
            (pairs["field"] == field)
            & (pairs["model_left"] == left)
            & (pairs["model_right"] == right)
        ]
        return float(row.iloc[0]["rmse_over_pooled_range_percent"])

    summaries: dict[str, pd.Series] = {}
    mass_errors: dict[str, float] = {}
    for model, run_directory in RUN_DIRECTORIES.items():
        run = source_root / run_directory
        summaries[model] = pd.read_csv(run / "simulation_summary.csv").iloc[0]
        mass = pd.read_csv(run / "component_mass_balance.csv")
        final_rows = mass[np.isclose(mass["time_day"], mass["time_day"].max())]
        mass_errors[model] = float(final_rows["relative_error"].abs().max())

    well_rows = {row["model"]: row for _, row in wells.iterrows()}
    bhp_values = wells["injector_bhp_bar"].to_numpy(dtype=float)
    rate_values = wells["producer_total_reservoir_m3_s"].to_numpy(dtype=float)
    bhp_span = float(np.ptp(bhp_values))
    rate_span = float(np.ptp(rate_values))
    bhp_span_percent = bhp_span / float(well_rows["New-PR"]["injector_bhp_bar"]) * 100.0
    rate_span_percent = (
        rate_span
        / float(well_rows["New-PR"]["producer_total_reservoir_m3_s"])
        * 100.0
    )

    cpa_phase_switch = (
        (data["New-CPA"]["oil_saturation"].to_numpy(dtype=float) <= PHASE_PRESENT_THRESHOLD)
        & (data["New-PR"]["oil_saturation"].to_numpy(dtype=float) > PHASE_PRESENT_THRESHOLD)
    )
    switch_coordinates = [
        f"({row.x_m:.1f}, {row.y_m:.1f})"
        for row in data["New-CPA"].loc[cpa_phase_switch, ["x_m", "y_m"]].itertuples()
    ]
    switch_coordinate_text = "、".join(switch_coordinates)

    pvi = final_day / 3652.5
    cpa_pr_cost_ratio = (
        float(summaries["New-CPA"]["simulation_loop_wall"])
        / float(summaries["New-PR"]["simulation_loop_wall"])
    )

    report = rf"""# H2O–CO2–nC10 二维算例结果分析

## 分析范围

本文比较 Traditional、New-PR、New-SW 和 New-CPA 在同一 60×20×1 规则二维算例中的结果。分析时刻为 {final_day:.2f} d，对应 {pvi:.3f} PVI。所有结论均由模拟 CSV 直接计算，图像只承担空间模式展示，不用于目测取数。

四模型绝对场采用每个物理量的共同最小值和最大值；除压力及压力梯度外，展示范围聚焦于全部驱替前缘所在的 x=0–100 m，统计量仍基于完整 300 m 计算域。全部空间色图使用用户提供的 6 锚点 MRST `self_color`，并按原 MATLAB 文件进行 256 色 PCHIP 插值；跨数量级场采用色条上明确标注的幂次归一化，不截断数据。差值场以 New-PR 为参照，并采用以零为中心的共同对称色标；饱和度及气相相对渗透率中的孤立极值使用带延伸端的第 99 百分位显示范围，原始极值保留在 CSV。相组成仅在对应相存在的单元中显示，灰色表示该相不存在，不表示组分摩尔分数为零。本分析属于对既有计算结果的描述性比较，不包含统计显著性检验。

## 主要发现

- 四组计算均到达 0.1 PVI，未发生井控切换；最终最大组分相对质量平衡误差不超过 {max(mass_errors.values()):.3e}。
- 压力场和主体气驱前缘具有较强一致性，四模型全域压力极差最大为 {spread_value('pressure_mpa', 'spread_maximum'):.2f} kPa。
- 井底压力差异很小，但生产井总地层流量仍出现 {rate_span_percent:.2f}% 的全幅差，说明热力学差异会通过相态、黏度和相对渗透率耦合传递到井响应。
- 模型差异主要集中在水–烃相间互溶；水相 CO2 的差异显著大于压力与主体饱和度差异。
- New-CPA 给出最高的油中水和水中 CO2，同时具有最高求解成本；修正后未再出现错误的 Gas+Water 特殊点。

## 1. 数值有效性

| 模型 | 接受/拒绝时间步 | 非线性求解 | 接受步 SNES | 接受步 KSP | 运行时间 (s) | 最终最大组分相对误差 |
|---|---:|---:|---:|---:|---:|---:|
| Traditional | {int(summaries['Traditional']['accepted_internal_steps'])}/{int(summaries['Traditional']['rejected_internal_steps'])} | {int(summaries['Traditional']['nonlinear_solves'])} | {int(summaries['Traditional']['accepted_snes'])} | {int(summaries['Traditional']['accepted_ksp'])} | {float(summaries['Traditional']['simulation_loop_wall']):.2f} | {mass_errors['Traditional']:.3e} |
| New-PR | {int(summaries['New-PR']['accepted_internal_steps'])}/{int(summaries['New-PR']['rejected_internal_steps'])} | {int(summaries['New-PR']['nonlinear_solves'])} | {int(summaries['New-PR']['accepted_snes'])} | {int(summaries['New-PR']['accepted_ksp'])} | {float(summaries['New-PR']['simulation_loop_wall']):.2f} | {mass_errors['New-PR']:.3e} |
| New-SW | {int(summaries['New-SW']['accepted_internal_steps'])}/{int(summaries['New-SW']['rejected_internal_steps'])} | {int(summaries['New-SW']['nonlinear_solves'])} | {int(summaries['New-SW']['accepted_snes'])} | {int(summaries['New-SW']['accepted_ksp'])} | {float(summaries['New-SW']['simulation_loop_wall']):.2f} | {mass_errors['New-SW']:.3e} |
| New-CPA | {int(summaries['New-CPA']['accepted_internal_steps'])}/{int(summaries['New-CPA']['rejected_internal_steps'])} | {int(summaries['New-CPA']['nonlinear_solves'])} | {int(summaries['New-CPA']['accepted_snes'])} | {int(summaries['New-CPA']['accepted_ksp'])} | {float(summaries['New-CPA']['simulation_loop_wall']):.2f} | {mass_errors['New-CPA']:.3e} |

四种模型的 `individual_well_control_switches` 均为 0，说明本次对比没有被不同井控路径混杂。质量守恒误差比本文讨论的场量差异低多个数量级，因此结果差异不能归因于可见的组分质量亏损。

## 2. 压力传播与井响应

四模型均形成由注入端向生产端平滑递减的压力场。全域四模型压力极差的平均值为 {spread_value('pressure_mpa', 'spread_mean'):.2f} kPa，最大值为 {spread_value('pressure_mpa', 'spread_maximum'):.2f} kPa。以 New-PR 为参照，Traditional、New-SW 和 New-CPA 的压力 RMSE 分别为 {difference_value('pressure_mpa', 'Traditional', 'difference_rmse') * 1e3:.2f}、{difference_value('pressure_mpa', 'New-SW', 'difference_rmse') * 1e3:.2f} 和 {difference_value('pressure_mpa', 'New-CPA', 'difference_rmse') * 1e3:.2f} kPa。New-CPA 在压力层面反而最接近 New-PR，表明 CPA 的主要影响没有表现为整体压力形态重构。

![四模型压力场，共同色标](report_figures_10m6/01_pressure_10m6.png)

最终注入井 BHP 为 {bhp_values.min():.3f}–{bhp_values.max():.3f} bar，全幅差为 {bhp_span:.3f} bar，相当于 New-PR 值的 {bhp_span_percent:.3f}%。生产井总地层流量为 {rate_values.min():.4e}–{rate_values.max():.4e} m³/s，全幅差相当于 New-PR 值的 {rate_span_percent:.2f}%。New-SW 的产量最高，Traditional 最低，New-CPA 位于 New-PR 与 New-SW 之间。

| 模型 | 注入井 BHP (bar) | 生产井总地层流量 (m³/s) | 油相流量 (m³/s) | 水相流量 (m³/s) |
|---|---:|---:|---:|---:|
| Traditional | {float(well_rows['Traditional']['injector_bhp_bar']):.4f} | {float(well_rows['Traditional']['producer_total_reservoir_m3_s']):.4e} | {float(well_rows['Traditional']['producer_oil_reservoir_m3_s']):.4e} | {float(well_rows['Traditional']['producer_water_reservoir_m3_s']):.4e} |
| New-PR | {float(well_rows['New-PR']['injector_bhp_bar']):.4f} | {float(well_rows['New-PR']['producer_total_reservoir_m3_s']):.4e} | {float(well_rows['New-PR']['producer_oil_reservoir_m3_s']):.4e} | {float(well_rows['New-PR']['producer_water_reservoir_m3_s']):.4e} |
| New-SW | {float(well_rows['New-SW']['injector_bhp_bar']):.4f} | {float(well_rows['New-SW']['producer_total_reservoir_m3_s']):.4e} | {float(well_rows['New-SW']['producer_oil_reservoir_m3_s']):.4e} | {float(well_rows['New-SW']['producer_water_reservoir_m3_s']):.4e} |
| New-CPA | {float(well_rows['New-CPA']['injector_bhp_bar']):.4f} | {float(well_rows['New-CPA']['producer_total_reservoir_m3_s']):.4e} | {float(well_rows['New-CPA']['producer_oil_reservoir_m3_s']):.4e} | {float(well_rows['New-CPA']['producer_water_reservoir_m3_s']):.4e} |

![最终井响应](report_figures_10m6/21_final_well_response_10m6.png)

## 3. 饱和度与气相前缘

四种模型的主体气相羽流形态相近。按气相饱和度大于 0.01 定义，影响面积为 {plumes['area_sg_gt_0.01_m2'].min():.0f}–{plumes['area_sg_gt_0.01_m2'].max():.0f} m²；按气相饱和度大于 0.10 定义，影响面积为 {plumes['area_sg_gt_0.10_m2'].min():.0f}–{plumes['area_sg_gt_0.10_m2'].max():.0f} m²。Traditional 的 0.10 前缘到达 57.5 m，其余三种模型到达 52.5 m。二者相差一个 5 m 网格，现阶段只能解释为单元尺度差异，不能据此断言连续尺度传播速度显著不同。

![四模型气相饱和度，共同色标](report_figures_10m6/04_gas_saturation_10m6.png)

![两级气相前缘轮廓](report_figures_10m6/27_gas_front_contours_10m6.png)

绝大多数单元的气相饱和度差异较小：四模型单元极差的 95% 分位数为 {spread_value('gas_saturation', 'spread_p95'):.2f} 个百分点。CPA 单一非水相修正后，最终场中 Gas+Water 单元数为 0；四模型气相饱和度最大单元极差为 {spread_value('gas_saturation', 'spread_maximum'):.2f} 个百分点。空间极差图将气相饱和度色标显示上限设为第 99 百分位，并保留超限单元，以避免局部极值掩盖主体前缘差异。

![四模型无参照空间极差](report_figures_10m6/23_four_model_spatial_spread_10m6.png)

该局部极差不宜直接解释为 CPA 预测了更快的整体前缘。更稳妥的判断是：不同热力学模型在相态边界附近产生了单元尺度响应差异。后续仍应通过相邻时间步追踪、网格加密和相态稳定性残差检查判断其网格依赖性。

## 4. 油中水

Traditional 按模型定义令油相 H2O 为零，因此其与三种新模型的差值同时包含“模型是否允许互溶”的结构性差异。New-PR、New-SW 和 New-CPA 的油相 H2O 范围分别为 {field_value('oil_x_h2o', 'New-PR', 'minimum') * 1e6:.0f}–{field_value('oil_x_h2o', 'New-PR', 'maximum') * 1e6:.0f}、{field_value('oil_x_h2o', 'New-SW', 'minimum') * 1e6:.0f}–{field_value('oil_x_h2o', 'New-SW', 'maximum') * 1e6:.0f} 和 {field_value('oil_x_h2o', 'New-CPA', 'minimum') * 1e6:.0f}–{field_value('oil_x_h2o', 'New-CPA', 'maximum') * 1e6:.0f} × 10^-6。New-CPA 的整体水平最高，并在注入影响区表现出更强的空间起伏。

![油相 H2O 摩尔分数，共同色标](report_figures_10m6/09_oil_x_h2o_10m6.png)

New-PR–New-SW 的成对归一化 RMSE 为 {pair_value('oil_x_h2o', 'New-PR', 'New-SW'):.2f}%，New-PR–New-CPA 为 {pair_value('oil_x_h2o', 'New-PR', 'New-CPA'):.2f}%。这说明普通 PR 与 SW 在油中水上的差异相对有限，而 CPA 缔合描述对该量的影响更明显。

## 5. 水中 CO2

水相 CO2 是四模型差异最强的物理量。Traditional 按模型定义为零；New-PR 的最大值为 {field_value('water_x_co2', 'New-PR', 'maximum') * 1e6:.0f} × 10^-6，New-SW 为 {field_value('water_x_co2', 'New-SW', 'maximum') * 1e6:.0f} × 10^-6，New-CPA 为 {field_value('water_x_co2', 'New-CPA', 'maximum') * 1e6:.0f} × 10^-6。富集主要位于注入端和气相前缘附近，说明局部相平衡处理直接控制 CO2 向水相的分配。

![水相 CO2 摩尔分数，共同色标](report_figures_10m6/16_water_x_co2_10m6.png)

![水相 CO2 相对 New-PR 的差值](report_figures_10m6/delta_16_water_x_co2_10m6.png)

New-PR–New-SW、New-PR–New-CPA 和 New-SW–New-CPA 的成对归一化 RMSE 分别为 {pair_value('water_x_co2', 'New-PR', 'New-SW'):.2f}%、{pair_value('water_x_co2', 'New-PR', 'New-CPA'):.2f}% 和 {pair_value('water_x_co2', 'New-SW', 'New-CPA'):.2f}%。因此，New-SW 与 New-CPA 在空间分布上更接近，但 New-CPA 仍给出更高的峰值与平均含量。对于溶解封存、含水相反应或水化学耦合问题，模型选择引起的不确定性远大于压力场差异。

![井连线上的油中水与水中 CO2](report_figures_10m6/25_phase_exchange_centerlines_10m6.png)

## 6. 综合差异

成对归一化 RMSE 定义为

$$
E_{{ab}}=\frac{{\sqrt{{N^{{-1}}\sum_{{i=1}}^N\left(u_{{a,i}}-u_{{b,i}}\right)^2}}}}
{{\max_{{m,i}}u_{{m,i}}-\min_{{m,i}}u_{{m,i}}}}\times100\%.
$$

该量使用四模型合并极差归一化，只用于比较空间分布差异，不是相对于某个“真实值”的预测误差。压力六组模型对的归一化 RMSE 为 {pairs.loc[pairs['field'] == 'pressure_mpa', 'rmse_over_pooled_range_percent'].min():.2f}%–{pairs.loc[pairs['field'] == 'pressure_mpa', 'rmse_over_pooled_range_percent'].max():.2f}%，气相饱和度为 {pairs.loc[pairs['field'] == 'gas_saturation', 'rmse_over_pooled_range_percent'].min():.2f}%–{pairs.loc[pairs['field'] == 'gas_saturation', 'rmse_over_pooled_range_percent'].max():.2f}%；相比之下，油中水最大达到 {pairs.loc[pairs['field'] == 'oil_x_h2o', 'rmse_over_pooled_range_percent'].max():.2f}%，水中 CO2 最大达到 {pairs.loc[pairs['field'] == 'water_x_co2', 'rmse_over_pooled_range_percent'].max():.2f}%。模型排序因此随观察量改变，不能用单一“整体最接近”概括所有物理量。

![全部模型对的归一化 RMSE](report_figures_10m6/24_pairwise_normalized_rmse_10m6.png)

![沿流动方向的平均绝对模型差](report_figures_10m6/26_x_resolved_model_differences_10m6.png)

## 7. 计算成本

New-CPA 的运行时间为 {float(summaries['New-CPA']['simulation_loop_wall']):.2f} s，是 New-PR 的 {cpa_pr_cost_ratio:.2f} 倍；其接受步 KSP 迭代为 {int(summaries['New-CPA']['accepted_ksp'])}，显著高于其余模型。该差异表明当前 CPA 雅可比与预条件组合具有更高代数求解成本。这里的时间仅对应本次单核运行、当前输出频率和求解器配置，不能直接外推为超算并行性能结论。

## 8. 结论与后续验证

本次 0.1 PVI 结果支持以下结论：压力传播、井控和主体气驱前缘对热力学模型较稳健；相间互溶对模型选择高度敏感，其中 New-CPA 给出最高的油中水和水中 CO2，New-SW 与 New-CPA 的水中 CO2 分布最接近。New-CPA 的计算成本，以及相态边界附近单元尺度差异，仍是进入 0.5 PVI 前最需要验证的两项内容。

当前结果只有 0.1 PVI 最终场。Traditional、New-PR 与 New-SW 保存两个固定输出区间，New-CPA 保存十个，因此现有文件不适合直接比较高分辨率瞬态演化。下一阶段应统一输出时刻并运行至 0.5 PVI，同时开展网格加密、时间步敏感性、孤立相态单元追踪和多进程性能测试。水中 CO2 与油中水的模型优劣最终还需要独立实验数据校核；本文只比较模型间差异，不把任何模型预设为真值。

## 数据与图件来源

| 内容 | 文件 |
|---|---|
| 绝对场范围与均值 | [field_metrics.csv](field_metrics.csv) |
| 相对 New-PR 的差值统计 | [difference_metrics.csv](difference_metrics.csv) |
| 最终井响应 | [well_final_metrics.csv](well_final_metrics.csv) |
| 气相羽流面积与前缘 | [plume_geometry.csv](plume_geometry.csv) |
| 四模型无参照空间极差 | [spatial_spread_metrics.csv](spatial_spread_metrics.csv) |
| 全部模型对归一化 RMSE | [pairwise_rmse_metrics.csv](pairwise_rmse_metrics.csv) |
| 源文件、掩膜、色标与软件版本 | [figure_manifest.json](figure_manifest.json) |

所有嵌入图均为同目录绘图流程生成的 600 dpi RGB PNG；PDF 和 SVG 版本与 PNG 同名，可用于论文排版或 PPT 后续编辑。原始模拟 CSV 未被修改。
"""
    (output / "RESULT_ANALYSIS.md").write_text(report, encoding="utf-8")


def main() -> None:
    global PRESENTATION_MODE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--presentation",
        action="store_true",
        help="use concise labels and larger type for PowerPoint/slides",
    )
    args = parser.parse_args()
    PRESENTATION_MODE = args.presentation
    configure_matplotlib(PRESENTATION_MODE)
    args.output.mkdir(parents=True, exist_ok=True)

    data, wells, final_day, source_paths = read_all(args.source_root)
    cmap = self_color_colormap().with_extremes(bad="#D9D9D9")
    absolute_metrics = pd.DataFrame(plot_absolute_maps(data, args.output, cmap))
    difference_metrics = pd.DataFrame(plot_difference_maps(data, args.output, cmap))
    plot_core_centerlines(data, args.output)
    plot_core_centerline_differences(data, args.output)
    plot_aqueous_co2_centerline(data, args.output)
    plot_oil_h2o_centerline(data, args.output)
    plot_co2_centerlines(data, args.output)
    plot_ecdf(data, args.output)
    plot_normalized_rmse_heatmap(difference_metrics, data, args.output, cmap)
    plot_well_dots(wells, args.output)
    plot_plume_geometry(data, args.output)
    plot_model_spread_maps(data, args.output, cmap)
    plot_pairwise_rmse_heatmap(data, args.output, cmap)
    plot_phase_exchange_profiles(data, args.output)
    plot_x_resolved_differences(data, args.output)
    plot_gas_front_contours(data, args.output)
    write_provenance(args.output, args.source_root, source_paths, final_day, absolute_metrics, difference_metrics)
    refresh_report_figures(args.output)
    write_result_analysis(args.output, args.source_root, data, final_day)
    print(f"wrote compact extended figures to {args.output}")


if __name__ == "__main__":
    main()
