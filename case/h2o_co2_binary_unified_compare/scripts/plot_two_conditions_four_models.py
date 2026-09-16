#!/usr/bin/env python3
"""Plot two-condition/four-model H2O-CO2 results in the established style."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize, PowerNorm
from matplotlib.lines import Line2D
from PIL import Image
from scipy.interpolate import PchipInterpolator


NX, NY = 60, 20
LX_M, LY_M = 300.0, 100.0
DX_M, DY_M = LX_M / NX, LY_M / NY
WELL_ROW = 9
WELL_Y_M = (WELL_ROW + 0.5) * DY_M
PHASE_EPS = 1.0e-8
RASTER_DPI = 600

MODELS = ("Traditional", "New-PR", "New-SW", "New-CPA")
MODEL_DIRS = {
    "Traditional": "traditional",
    "New-PR": "new_pr",
    "New-SW": "new_sw",
    "New-CPA": "new_cpa",
}
CONDITIONS = {
    "conventional": {
        "label": "常规工况  333.15 K, 5.16 MPa",
        "short": "333.15 K, 5.16 MPa",
        "pressure_mpa": 5.16,
    },
    "scw": {
        "label": "超临界水工况  653.15 K, 28 MPa",
        "short": "653.15 K, 28 MPa",
        "pressure_mpa": 28.0,
    },
}
MODEL_COLORS = {
    "Traditional": "#0072B2",
    "New-PR": "#D55E00",
    "New-SW": "#009E73",
    "New-CPA": "#CC79A7",
}
MODEL_LINESTYLES = {
    "Traditional": "-",
    "New-PR": (0.0, (5.5, 2.2)),
    "New-SW": (0.0, (4.0, 1.6, 1.0, 1.6)),
    "New-CPA": (0.0, (1.0, 1.5)),
}
SELF_COLOR_ANCHORS = np.array(
    [
        [0.300, 0.550, 0.750],
        [0.550, 0.750, 0.850],
        [0.700, 0.800, 0.800],
        [0.950, 0.850, 0.700],
        [0.950, 0.650, 0.400],
        [0.950, 0.550, 0.300],
    ],
    dtype=float,
)


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Noto Serif SC", "STSong", "DejaVu Serif"],
            "font.sans-serif": ["Microsoft YaHei", "Noto Sans SC", "SimHei"],
            "mathtext.fontset": "stix",
            "font.size": 9.5,
            "axes.titlesize": 10.0,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
            "axes.edgecolor": "#7890A0",
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.fontsize": 8.2,
            "figure.dpi": 180,
            "savefig.dpi": RASTER_DPI,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
        }
    )


def self_color_colormap() -> ListedColormap:
    anchor_index = np.arange(len(SELF_COLOR_ANCHORS), dtype=float)
    target_index = np.linspace(0.0, len(SELF_COLOR_ANCHORS) - 1.0, 256)
    colors = np.column_stack(
        [PchipInterpolator(anchor_index, SELF_COLOR_ANCHORS[:, i])(target_index) for i in range(3)]
    )
    return ListedColormap(np.clip(colors, 0.0, 1.0), name="self_color")


def read_frame(path: Path, model: str, initial_pressure_mpa: float) -> pd.DataFrame:
    frame = pd.read_csv(path).sort_values("input_index").reset_index(drop=True)
    if len(frame) != NX * NY or not np.array_equal(frame["input_index"], np.arange(NX * NY)):
        raise ValueError(f"{path}: incomplete 60 x 20 field")
    frame["x_m"] = (frame["input_index"] % NX + 0.5) * DX_M
    frame["y_m"] = (frame["input_index"] // NX + 0.5) * DY_M
    frame["pressure_mpa"] = frame["pressure_Pa"] * 1.0e-6
    frame["pressure_change_mpa"] = frame["pressure_mpa"] - initial_pressure_mpa
    if model == "Traditional":
        frame["gas_saturation"] = frame["vapor_saturation"]
        frame["water_x_co2"] = 0.0
        frame["gas_y_h2o"] = 0.0
        oil_saturation = frame["liquid_saturation"]
    else:
        frame["water_x_co2"] = frame["water_x_1_CO2"]
        frame["gas_y_h2o"] = frame["gas_y_0_H2O"]
        oil_saturation = frame["oil_saturation"]
    closure = oil_saturation + frame["gas_saturation"] + frame["water_saturation"]
    closure_error = float(np.max(np.abs(closure - 1.0)))
    if closure_error > 5.0e-8:
        raise ValueError(f"{path}: saturation closure {closure_error:.3e}")
    required = [
        "pressure_change_mpa", "gas_saturation", "water_saturation", "water_x_co2"
    ]
    if not np.isfinite(frame[required].to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: non-finite primary field")
    frame["gas_y_h2o"] = frame["gas_y_h2o"].where(frame["gas_saturation"] > PHASE_EPS)
    frame.attrs["saturation_closure_error"] = closure_error
    return frame


def read_all(root: Path):
    data, run_dirs = {}, {}
    for condition, meta in CONDITIONS.items():
        data[condition] = {}
        run_dirs[condition] = {}
        for model in MODELS:
            run_dir = root / "formal" / condition / MODEL_DIRS[model]
            run_dirs[condition][model] = run_dir
            data[condition][model] = read_frame(
                run_dir / "solution_final.csv", model, meta["pressure_mpa"]
            )
    return data, run_dirs


def add_wells(axis: plt.Axes, x_max: float) -> None:
    axis.scatter(
        [0.5 * DX_M], [WELL_Y_M], marker="^", s=25,
        facecolors="white", edgecolors="#2B8CBE", linewidths=0.95, zorder=4,
    )
    if x_max >= LX_M:
        axis.scatter(
            [LX_M - 0.5 * DX_M], [WELL_Y_M], marker="v", s=25,
            facecolors="#173F5F", edgecolors="white", linewidths=0.6, zorder=4,
        )


def style_map_axis(axis: plt.Axes, x_max: float, show_xlabel: bool, show_ylabel: bool) -> None:
    axis.set_xlim(0.0, x_max)
    axis.set_ylim(0.0, LY_M)
    axis.set_aspect("equal")
    axis.set_xticks([0, 25, 50, 75, 100] if x_max == 100.0 else [0, 100, 200, 300])
    axis.set_yticks([0, 50, 100])
    axis.set_xlabel("x (m)" if show_xlabel else "", labelpad=1)
    axis.set_ylabel("y (m)" if show_ylabel else "", labelpad=1)
    add_wells(axis, x_max)


def style_profile_axis(axis: plt.Axes) -> None:
    axis.set_facecolor("white")
    axis.grid(False)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("#6F7F8A")
        spine.set_linewidth(0.85)


def save_figure(figure: plt.Figure, output: Path, stem: str) -> None:
    for suffix in (".png", ".pdf", ".svg"):
        figure.savefig(output / f"{stem}{suffix}", dpi=RASTER_DPI, bbox_inches="tight")
    plt.close(figure)
    with Image.open(output / f"{stem}.png") as image:
        image.convert("RGB").save(output / f"{stem}.png", dpi=(RASTER_DPI, RASTER_DPI), optimize=True)


def field_limits(frames: dict[str, pd.DataFrame], field: str, scale: float):
    values = np.concatenate([frame[field].to_numpy(dtype=float) * scale for frame in frames.values()])
    values = values[np.isfinite(values)]
    vmin, vmax = float(np.min(values)), float(np.max(values))
    if vmin == vmax:
        vmax = float(np.nextafter(vmax, np.inf))
    return vmin, vmax


def plot_four_model_field(
    data, condition: str, output: Path, *, field: str, label: str, stem: str,
    scale: float = 1.0, x_max: float = 100.0, gamma: float | None = None,
) -> dict:
    frames = data[condition]
    vmin, vmax = field_limits(frames, field, scale)
    cmap = self_color_colormap()
    norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax) if gamma else Normalize(vmin=vmin, vmax=vmax)
    figure = plt.figure(figsize=(13.2, 4.45), facecolor="white")
    grid = figure.add_gridspec(
        2, 4, width_ratios=(1.0, 1.0, 0.055, 1.10), wspace=0.24, hspace=0.14
    )
    map_axes = [figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1]),
                figure.add_subplot(grid[1, 0]), figure.add_subplot(grid[1, 1])]
    image = None
    for index, (axis, model) in enumerate(zip(map_axes, MODELS)):
        values = frames[model][field].to_numpy(dtype=float).reshape(NY, NX) * scale
        image = axis.imshow(
            values, origin="lower", extent=(0, LX_M, 0, LY_M), interpolation="nearest",
            aspect="equal", cmap=cmap, norm=norm,
        )
        axis.set_title(model, loc="left", pad=2.5)
        style_map_axis(axis, x_max, show_xlabel=index >= 2, show_ylabel=index % 2 == 0)
        if index < 2:
            axis.tick_params(labelbottom=False)
        if index % 2 == 1:
            axis.tick_params(labelleft=False)
    assert image is not None
    colorbar_axis = figure.add_subplot(grid[:, 2])
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label(label)
    colorbar.outline.set_visible(False)

    profile = figure.add_subplot(grid[:, 3])
    center = slice(WELL_ROW * NX, (WELL_ROW + 1) * NX)
    for model in MODELS:
        row = frames[model].iloc[center]
        profile.plot(
            row["x_m"], row[field] * scale, color=MODEL_COLORS[model],
            linestyle=MODEL_LINESTYLES[model], linewidth=1.8, label=model,
            solid_capstyle="round", dash_capstyle="round",
        )
    profile.set_xlim(0.0, LX_M)
    profile.set_xticks([0, 100, 200, 300])
    profile.set_xlabel("x (m)")
    profile.set_title(f"Centerline  ·  {label}", loc="left", pad=4)
    style_profile_axis(profile)
    profile.legend(loc="best", frameon=False, ncol=2)

    figure.text(
        0.045, 0.985, f"{CONDITIONS[condition]['label']}  ·  {label}",
        ha="left", va="top", color="white", fontsize=14.0, fontweight="bold",
        fontfamily="Microsoft YaHei",
        bbox={"boxstyle": "square,pad=0.42", "facecolor": "#0072B2", "edgecolor": "none"},
    )
    figure.subplots_adjust(left=0.055, right=0.98, bottom=0.13, top=0.875)
    save_figure(figure, output, stem)
    return {
        "condition": condition, "field": field, "display_label": label,
        "common_vmin": vmin, "common_vmax": vmax, "display_scale": scale,
        "power_gamma": 1.0 if gamma is None else gamma, "display_x_max_m": x_max,
    }


def plot_centerline_summary(data, output: Path) -> None:
    specs = (
        ("pressure_change_mpa", 1.0, r"$\Delta p$ (MPa)"),
        ("gas_saturation", 1.0, r"$S_\mathrm{g}$"),
        ("water_x_co2", 1.0e3, r"$x^\mathrm{w}_{\mathrm{CO_2}}$ ($10^{-3}$)"),
    )
    center = slice(WELL_ROW * NX, (WELL_ROW + 1) * NX)
    figure, axes = plt.subplots(2, 3, figsize=(12.3, 6.1), sharex=True)
    for row_index, condition in enumerate(CONDITIONS):
        for column_index, (field, scale, label) in enumerate(specs):
            axis = axes[row_index, column_index]
            for model in MODELS:
                row = data[condition][model].iloc[center]
                axis.plot(
                    row["x_m"], row[field] * scale, color=MODEL_COLORS[model],
                    linestyle=MODEL_LINESTYLES[model], linewidth=1.75,
                    solid_capstyle="round", dash_capstyle="round",
                )
            axis.set_xlim(0, LX_M)
            axis.set_xticks([0, 100, 200, 300])
            axis.set_ylabel(label)
            axis.set_xlabel("x (m)" if row_index == 1 else "")
            axis.set_title(CONDITIONS[condition]["short"] if column_index == 0 else "", loc="left", pad=4)
            style_profile_axis(axis)
    handles = [
        Line2D([0], [0], color=MODEL_COLORS[m], linestyle=MODEL_LINESTYLES[m], linewidth=2.0, label=m)
        for m in MODELS
    ]
    figure.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.53, 0.005))
    figure.text(
        0.04, 0.985, "统一初始组成：两种工况、四种模型中心线对比",
        ha="left", va="top", color="white", fontsize=14.0, fontweight="bold",
        fontfamily="Microsoft YaHei",
        bbox={"boxstyle": "square,pad=0.42", "facecolor": "#0072B2", "edgecolor": "none"},
    )
    figure.subplots_adjust(left=0.085, right=0.985, bottom=0.13, top=0.88, hspace=0.22, wspace=0.27)
    save_figure(figure, output, "07_two_condition_centerline_summary")


def write_supporting_files(data, run_dirs, source_root: Path, output: Path, plot_metrics: list[dict]) -> None:
    source = pd.DataFrame({
        "input_index": np.arange(NX * NY),
        "x_m": data["conventional"]["Traditional"]["x_m"],
        "y_m": data["conventional"]["Traditional"]["y_m"],
    })
    field_rows, validation_rows = [], []
    for condition in CONDITIONS:
        for model in MODELS:
            frame = data[condition][model]
            prefix = f"{condition}_{MODEL_DIRS[model]}"
            for field in ("pressure_change_mpa", "gas_saturation", "water_saturation", "water_x_co2", "gas_y_h2o"):
                source[f"{prefix}_{field}"] = frame[field].to_numpy(dtype=float)
                values = frame[field].to_numpy(dtype=float)
                finite = values[np.isfinite(values)]
                field_rows.append({
                    "condition": condition, "model": model, "field": field,
                    "minimum": float(np.min(finite)) if len(finite) else np.nan,
                    "mean": float(np.mean(finite)) if len(finite) else np.nan,
                    "maximum": float(np.max(finite)) if len(finite) else np.nan,
                    "finite_cells": int(len(finite)),
                })
            run_dir = run_dirs[condition][model]
            summary = pd.read_csv(run_dir / "simulation_summary.csv").iloc[-1]
            balance = pd.read_csv(run_dir / "component_mass_balance.csv")
            balance = balance[balance["step"] == balance["step"].max()].copy()
            active = balance[balance["component"].isin(["H2O", "CO2"])]
            nc10 = balance[balance["component"] == "nC10"].iloc[0]
            total_rel = float(balance["balance_error_kg"].sum() / balance["initial_inventory_kg"].sum())
            gas_present = frame["gas_saturation"].to_numpy(dtype=float) > PHASE_EPS
            water_present = frame["water_saturation"].to_numpy(dtype=float) > PHASE_EPS
            validation_rows.append({
                "condition": condition, "model": model,
                "temperature_K": 333.15 if condition == "conventional" else 653.15,
                "initial_pressure_MPa": CONDITIONS[condition]["pressure_mpa"],
                "initial_z_H2O": 0.8, "initial_z_CO2": 0.2, "initial_z_nC10": 0.0,
                "final_time_day": float(summary["final_time_day"]),
                "accepted_steps": int(summary["accepted_internal_steps"]),
                "rejected_steps": int(summary["rejected_internal_steps"]),
                "well_control_switches": int(summary["individual_well_control_switches"]),
                "max_active_component_relative_error": float(active["relative_error"].abs().max()),
                "total_relative_error": total_rel,
                "inactive_nC10_absolute_trace_kg": float(abs(nc10["balance_error_kg"])),
                "gas_only_cells": int(np.sum(gas_present & ~water_present)),
                "water_only_cells": int(np.sum(~gas_present & water_present)),
                "gas_water_cells": int(np.sum(gas_present & water_present)),
                "max_saturation_closure_error": frame.attrs["saturation_closure_error"],
            })
    source.to_csv(output / "figure_source_data.csv", index=False)
    pd.DataFrame(field_rows).to_csv(output / "field_metrics.csv", index=False)
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(output / "validation_summary.csv", index=False)
    png_hashes = {}
    for path in sorted(output.glob("*.png")):
        png_hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "experiment": "Unified H2O-CO2 initial composition, two conditions, four models",
        "source_root": str(source_root.resolve()),
        "source_solution_files": {
            c: {m: str((run_dirs[c][m] / "solution_final.csv").resolve()) for m in MODELS}
            for c in CONDITIONS
        },
        "conditions": CONDITIONS,
        "initial_overall_mole_fraction": {"H2O": 0.8, "CO2": 0.2, "nC10": 0.0},
        "grid": {"nx": NX, "ny": NY, "nz": 1, "lx_m": LX_M, "ly_m": LY_M, "lz_m": 5.0},
        "final_time_day": 365.25,
        "final_pvi": 0.1,
        "transformations": [
            "cells sorted by input_index and reshaped to 20 x 60",
            "pressure converted Pa to MPa and referenced to each condition's own initial pressure",
            "Traditional vapor saturation mapped to gas saturation",
            "Traditional aqueous CO2 mole fraction set to structural zero because dissolution is disabled",
            "phase-specific gas composition masked where gas saturation <= 1e-8",
            "nearest-neighbor cell maps; no spatial interpolation",
            "common color limits across all four models within each condition and field",
        ],
        "plot_metrics": plot_metrics,
        "style": {
            "basis": "D:/DOCUMENTS/ChatGPT/MPMC_SCW_V2/RESULTS",
            "colormap": "six-anchor MRST self_color with PCHIP interpolation",
            "model_line_styles": {m: str(MODEL_LINESTYLES[m]) for m in MODELS},
            "raster_dpi": RASTER_DPI,
        },
        "software": {
            "python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "matplotlib": mpl.__version__,
        },
        "png_sha256": png_hashes,
    }
    (output / "figure_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    passed = bool(
        (validation["final_time_day"] == 365.25).all()
        and (validation["rejected_steps"] == 0).all()
        and (validation["max_active_component_relative_error"] <= 1.0e-6).all()
        and (validation["total_relative_error"].abs() <= 1.0e-6).all()
        and (validation["max_saturation_closure_error"] <= 5.0e-8).all()
    )
    (output / "README.md").write_text(
        "# 统一初始组成的两工况四模型结果\n\n"
        "初始总体摩尔组成为 H2O/CO2/nC10 = 0.80/0.20/0。两种工况分别为 "
        "333.15 K、5.16 MPa 和 653.15 K、28 MPa；终止时刻为 365.25 d（0.1 PVI）。\n\n"
        f"八组计算验收：{'通过' if passed else '未通过'}。图像继承旧 RESULTS 目录的配色、字体、井标记、"
        "四模型布局和线型，并同时导出 600 dpi PNG、PDF、SVG。所有地图使用原始网格值，不做空间插值。"
        "每张图内四模型共用色标。超临界工况的新模型为富水单相，因此其气相组成不被伪造。\n\n"
        "详细数值见 `validation_summary.csv`、`field_metrics.csv` 和 `figure_source_data.csv`。\n",
        encoding="utf-8",
    )
    (output / "ALT_TEXT.md").write_text(
        """# 图像说明\n\n01–06 每张图左侧为 Traditional、New-PR、New-SW、New-CPA 四个网格图，右侧为井排中心线剖面。常规工况显示气水两相置换；超临界工况下三种新模型保持富水单相，而 Traditional 因独立水相假设仍保留气相。07 将两种工况的压降、气饱和度和水相 CO2 摩尔分数中心线放在同一张六面板图中。\n""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)
    data, run_dirs = read_all(args.source_root)
    plot_metrics = []
    field_specs = (
        ("pressure_change_mpa", r"$\Delta p$ (MPa)", "pressure_change", 1.0, 300.0, None),
        ("gas_saturation", r"$S_\mathrm{g}$", "gas_saturation", 1.0, 100.0, None),
        ("water_x_co2", r"$x^\mathrm{w}_{\mathrm{CO_2}}$ ($10^{-3}$)", "aqueous_co2", 1.0e3, 100.0, 0.65),
    )
    index = 1
    for condition in CONDITIONS:
        for field, label, name, scale, x_max, gamma in field_specs:
            plot_metrics.append(plot_four_model_field(
                data, condition, args.output, field=field, label=label,
                stem=f"{index:02d}_{condition}_{name}_four_models",
                scale=scale, x_max=x_max, gamma=gamma,
            ))
            index += 1
    plot_centerline_summary(data, args.output)
    write_supporting_files(data, run_dirs, args.source_root, args.output, plot_metrics)


if __name__ == "__main__":
    main()
