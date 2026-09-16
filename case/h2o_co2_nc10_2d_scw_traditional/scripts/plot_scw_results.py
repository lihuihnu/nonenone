#!/usr/bin/env python3
"""Plot and validate the 60x20x1 supercritical-water displacement result."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image


NX, NY = 60, 20
LX_M, LY_M = 300.0, 100.0
BASE_PRESSURE_PA = 5.16e6
SCW_PRESSURE_PA = 28.0e6
BASE_TEMPERATURE_K = 333.15
SCW_TEMPERATURE_K = 653.15
WELL_ROW = 9


def configure_style() -> None:
    chinese_font = Path("C:/Windows/Fonts/msyh.ttc")
    font_family = "DejaVu Sans"
    if chinese_font.exists():
        mpl.font_manager.fontManager.addfont(chinese_font)
        font_family = mpl.font_manager.FontProperties(fname=chinese_font).get_name()
    mpl.rcParams.update(
        {
            "font.family": font_family,
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
        }
    )


def read_case(directory: Path, initial_pressure_pa: float) -> tuple[pd.DataFrame, float]:
    field_path = directory / "solution_final.csv"
    summary_path = directory / "simulation_summary.csv"
    frame = pd.read_csv(field_path).sort_values("input_index").reset_index(drop=True)
    summary = pd.read_csv(summary_path)
    if len(frame) != NX * NY:
        raise ValueError(f"{field_path}: expected {NX * NY} cells, found {len(frame)}")
    if not np.array_equal(frame["input_index"].to_numpy(), np.arange(NX * NY)):
        raise ValueError(f"{field_path}: incomplete or unordered input_index")
    required = {
        "pressure_Pa",
        "liquid_saturation",
        "vapor_saturation",
        "water_saturation",
        "liquid_x_0_CO2",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{field_path}: missing columns {sorted(missing)}")
    numeric = frame[list(required)].to_numpy(float)
    if not np.isfinite(numeric).all():
        raise ValueError(f"{field_path}: non-finite values")
    closure = (
        frame["liquid_saturation"]
        + frame["vapor_saturation"]
        + frame["water_saturation"]
    )
    if float(np.max(np.abs(closure - 1.0))) > 5.0e-8:
        raise ValueError(f"{field_path}: saturation closure exceeds 5e-8")
    frame["pressure_MPa"] = frame["pressure_Pa"] * 1.0e-6
    frame["pressure_change_MPa"] = (frame["pressure_Pa"] - initial_pressure_pa) * 1.0e-6
    frame["water_saturation_change_pp"] = (frame["water_saturation"] - 0.20) * 100.0
    frame["oil_saturation_change_pp"] = (frame["liquid_saturation"] - 0.80) * 100.0
    final_day = float(summary["final_time_day"].iloc[-1])
    return frame, final_day


def grid(frame: pd.DataFrame, field: str) -> np.ndarray:
    return frame[field].to_numpy(float).reshape(NY, NX)


def add_wells(ax: plt.Axes) -> None:
    y = (WELL_ROW + 0.5) * LY_M / NY
    ax.scatter(
        [(0.5) * LX_M / NX], [y], marker="o", s=34,
        facecolors="none", edgecolors="black", linewidths=1.1, zorder=5,
    )
    ax.scatter(
        [(59.5) * LX_M / NX], [y], marker="x", s=34,
        color="black", linewidths=1.1, zorder=5,
    )


def draw_map(
    ax: plt.Axes,
    values: np.ndarray,
    title: str,
    colorbar_label: str,
    cmap: str,
    norm: mpl.colors.Normalize,
) -> None:
    image = ax.imshow(
        values,
        origin="lower",
        extent=(0.0, LX_M, 0.0, LY_M),
        aspect="equal",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )
    add_wells(ax)
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_xticks([0, 100, 200, 300])
    ax.set_yticks([0, 50, 100])
    colorbar = ax.figure.colorbar(image, ax=ax, fraction=0.047, pad=0.025)
    colorbar.set_label(colorbar_label)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=400)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    plt.close(fig)
    with Image.open(stem.with_suffix(".png")) as image:
        image.convert("RGB").save(stem.with_suffix(".png"), dpi=(400, 400), optimize=True)


def final_field_figure(scw: pd.DataFrame, final_day: float, output: Path) -> None:
    specs = (
        ("pressure_MPa", "Pressure", "Pressure (MPa)", "cividis"),
        ("liquid_saturation", "Oil saturation", "Oil saturation", "viridis"),
        ("water_saturation", "SCW saturation", "SCW saturation", "viridis"),
        ("liquid_x_0_CO2", "Oil-phase CO2 mole fraction", "Mole fraction", "magma"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 5.9), layout="constrained")
    for label, (field, title, cbar, cmap) in zip("abcd", specs):
        values = grid(scw, field)
        norm = Normalize(float(np.min(values)), float(np.max(values)))
        draw_map(axes.flat["abcd".index(label)], values, f"({label}) {title}", cbar, cmap, norm)
    fig.suptitle(
        f"60x20x1 SCW displacement at {final_day:g} d (0.1 PVI)\n"
        f"T = {SCW_TEMPERATURE_K:.2f} K, initial p = {SCW_PRESSURE_PA * 1e-6:.1f} MPa",
        fontsize=13,
    )
    save_figure(fig, output / "01_scw_final_fields")


def comparison_figure(
    baseline: pd.DataFrame,
    scw: pd.DataFrame,
    final_day: float,
    output: Path,
) -> None:
    columns = (
        ("pressure_change_MPa", "Pressure change (MPa)"),
        ("water_saturation_change_pp", "Water/SCW saturation change (pp)"),
        ("liquid_x_0_CO2", "Oil-phase CO2 mole fraction"),
    )
    fig, axes = plt.subplots(3, 3, figsize=(13.2, 7.7), layout="constrained")
    row_data = (baseline, scw)
    row_names = (
        f"Baseline\n{BASE_TEMPERATURE_K:.2f} K, {BASE_PRESSURE_PA * 1e-6:.2f} MPa",
        f"SCW\n{SCW_TEMPERATURE_K:.2f} K, {SCW_PRESSURE_PA * 1e-6:.1f} MPa",
    )
    for column, (field, label) in enumerate(columns):
        pooled = np.concatenate([frame[field].to_numpy(float) for frame in row_data])
        if field == "liquid_x_0_CO2":
            common_norm: mpl.colors.Normalize = Normalize(float(np.min(pooled)), float(np.max(pooled)))
            common_cmap = "magma"
        else:
            limit = max(float(np.max(np.abs(pooled))), np.finfo(float).eps)
            common_norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
            common_cmap = "RdBu_r"
        for row, frame in enumerate(row_data):
            values = grid(frame, field)
            image = axes[row, column].imshow(
                values, origin="lower", extent=(0.0, LX_M, 0.0, LY_M),
                aspect="equal", interpolation="nearest", cmap=common_cmap,
                norm=common_norm,
            )
            add_wells(axes[row, column])
            axes[row, column].set_title(label if row == 0 else "")
            if column == 0:
                axes[row, column].set_ylabel(f"{row_names[row]}\ny (m)")
            axes[row, column].figure.colorbar(
                image, ax=axes[row, column], fraction=0.047, pad=0.02,
            )

        difference = grid(scw, field) - grid(baseline, field)
        limit = max(float(np.max(np.abs(difference))), np.finfo(float).eps)
        difference_image = axes[2, column].imshow(
            difference, origin="lower", extent=(0.0, LX_M, 0.0, LY_M),
            aspect="equal", interpolation="nearest", cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
        )
        add_wells(axes[2, column])
        if column == 0:
            axes[2, column].set_ylabel("SCW - baseline\ny (m)")
        axes[2, column].figure.colorbar(
            difference_image, ax=axes[2, column], fraction=0.047, pad=0.02,
        )

    for ax in axes.flat:
        ax.set_xticks([0, 100, 200, 300])
        ax.set_yticks([0, 50, 100])
    for ax in axes[-1, :]:
        ax.set_xlabel("x (m)")
    fig.suptitle(
        f"Thermodynamic-condition sensitivity at {final_day:g} d (0.1 PVI)",
        fontsize=13,
    )
    save_figure(fig, output / "02_baseline_vs_scw")


def centerline_figure(
    baseline: pd.DataFrame,
    scw: pd.DataFrame,
    final_day: float,
    output: Path,
) -> None:
    cell_slice = slice(WELL_ROW * NX, (WELL_ROW + 1) * NX)
    x = (np.arange(NX) + 0.5) * LX_M / NX
    fig, axes = plt.subplots(3, 1, figsize=(7.4, 7.5), sharex=True, layout="constrained")
    series = (
        ("pressure_change_MPa", "Pressure change (MPa)"),
        ("water_saturation_change_pp", "Water/SCW saturation change (pp)"),
        ("liquid_x_0_CO2", "Oil-phase CO2 mole fraction"),
    )
    styles = (
        (baseline, "Baseline: 333.15 K, 5.16 MPa", "#0072B2", "-", "o"),
        (scw, "SCW: 653.15 K, 28 MPa", "#D55E00", "--", "s"),
    )
    for ax, (field, ylabel) in zip(axes, series):
        for frame, label, color, linestyle, marker in styles:
            ax.plot(
                x, frame[field].to_numpy(float)[cell_slice], label=label,
                color=color, linestyle=linestyle, linewidth=1.7,
                marker=marker, markevery=6, markersize=3.5,
            )
        ax.axhline(0.0, color="#777777", linewidth=0.7, linestyle=":")
        ax.set_ylabel(ylabel)
        ax.grid(True, color="#D9E1E5", linewidth=0.6)
    axes[0].legend(frameon=False, loc="best")
    axes[-1].set_xlabel("x along injector-producer row (m)")
    fig.suptitle(f"Injector-producer centerline at {final_day:g} d", fontsize=13)
    save_figure(fig, output / "03_centerline_profiles")


def style_matched_overview(
    baseline: pd.DataFrame,
    scw: pd.DataFrame,
    final_day: float,
    output: Path,
) -> None:
    """Create a dashboard matching the earlier blue-header rounded-card style."""
    figure = plt.figure(figsize=(16.0, 5.6), facecolor="#FBF8F6")

    left_bounds = (0.018, 0.075, 0.585, 0.845)
    right_bounds = (0.620, 0.075, 0.362, 0.845)
    for bounds in (left_bounds, right_bounds):
        figure.add_artist(
            FancyBboxPatch(
                (bounds[0], bounds[1]), bounds[2], bounds[3],
                boxstyle="round,pad=0.008,rounding_size=0.018",
                transform=figure.transFigure,
                facecolor="white", edgecolor="#CBCBCB", linewidth=1.25,
                zorder=-5,
            )
        )

    def blue_header(x: float, y: float, width: float, text: str) -> None:
        figure.add_artist(
            Rectangle(
                (x, y), width, 0.055, transform=figure.transFigure,
                facecolor="#0072BC", edgecolor="none", zorder=8,
            )
        )
        figure.text(
            x + 0.010, y + 0.028, text, color="white", fontsize=14.0,
            fontweight="bold", ha="left", va="center", zorder=9,
        )

    blue_header(0.039, 0.895, 0.285, "653.15 K、28 MPa 超临界水最终场")
    blue_header(0.643, 0.895, 0.215, "注采井中心线对比")

    map_grid = figure.add_gridspec(
        2, 2, left=0.050, right=0.575, bottom=0.185, top=0.845,
        wspace=0.38, hspace=0.48,
    )
    map_axes = np.asarray(
        [[figure.add_subplot(map_grid[row, column]) for column in range(2)] for row in range(2)]
    )
    front_x_m = 120.0
    front_columns = int(round(front_x_m / LX_M * NX))
    field_specs = (
        ("pressure_MPa", "(a) 压力", "MPa", "cividis"),
        ("liquid_saturation", "(b) 油相饱和度", "饱和度", "YlOrBr"),
        ("water_saturation", "(c) 超临界水饱和度", "饱和度", "Blues"),
        ("liquid_x_0_CO2", r"(d) 油相 CO$_2$ 摩尔分数", "摩尔分数", "magma"),
    )
    for axis, (field, title, colorbar_label, colormap) in zip(map_axes.flat, field_specs):
        values = grid(scw, field)[:, :front_columns]
        image = axis.imshow(
            values, origin="lower", extent=(0.0, front_x_m, 0.0, LY_M),
            aspect="equal", interpolation="nearest", cmap=colormap,
            norm=Normalize(float(np.min(values)), float(np.max(values))),
        )
        axis.scatter(
            [0.5 * LX_M / NX], [(WELL_ROW + 0.5) * LY_M / NY],
            marker="^", s=22, facecolors="white", edgecolors="#0072BC",
            linewidths=0.9, zorder=5,
        )
        axis.set_title(title, loc="left", fontsize=9.6, pad=3)
        axis.set_xlabel("x (m)")
        axis.set_ylabel("y (m)")
        axis.set_xticks([0, 40, 80, 120])
        axis.set_yticks([0, 50, 100])
        axis.tick_params(labelsize=7.4, length=2.5)
        colorbar = figure.colorbar(image, ax=axis, fraction=0.045, pad=0.025)
        colorbar.set_label(colorbar_label, fontsize=7.8)
        colorbar.ax.tick_params(labelsize=7.0, length=2.0)

    profile_grid = figure.add_gridspec(
        3, 1, left=0.655, right=0.955, bottom=0.165, top=0.845,
        hspace=0.20,
    )
    profile_axes = [figure.add_subplot(profile_grid[index]) for index in range(3)]
    cell_slice = slice(WELL_ROW * NX, (WELL_ROW + 1) * NX)
    x = (np.arange(NX) + 0.5) * LX_M / NX
    profile_specs = (
        ("pressure_change_MPa", "压力变化 (MPa)"),
        ("water_saturation_change_pp", "水相饱和度变化 (百分点)"),
        ("liquid_x_0_CO2", r"油相 CO$_2$ 摩尔分数"),
    )
    profile_styles = (
        (baseline, "333.15 K、5.16 MPa", "#0072B2", "-", "o"),
        (scw, "653.15 K、28 MPa", "#D55E00", "--", "s"),
    )
    for axis, (field, ylabel) in zip(profile_axes, profile_specs):
        for frame, label, color, linestyle, marker in profile_styles:
            axis.plot(
                x, frame[field].to_numpy(float)[cell_slice],
                color=color, linestyle=linestyle, linewidth=1.55,
                marker=marker, markevery=7, markersize=3.0, label=label,
            )
        axis.axhline(0.0, color="#888888", linewidth=0.65, linestyle=":")
        axis.set_ylabel(ylabel, fontsize=8.4)
        axis.tick_params(labelsize=7.5, length=2.5)
        axis.grid(True, color="#DCE4E8", linewidth=0.55)
        axis.set_xlim(0.0, LX_M)
    profile_axes[0].legend(
        loc="best", frameon=False, fontsize=8.0, ncol=2,
        handlelength=2.8, columnspacing=1.1,
    )
    profile_axes[0].tick_params(labelbottom=False)
    profile_axes[1].tick_params(labelbottom=False)
    profile_axes[-1].set_xlabel("注入井—生产井方向 x (m)", fontsize=8.8)

    figure.text(
        0.050, 0.103,
        "二维场显示注入前缘 x = 0–120 m（无插值）；△为注入井位置，中心线覆盖完整 300 m。",
        fontsize=8.2, color="#555555", ha="left", va="center",
    )
    figure.text(
        0.973, 0.028, f"t = {final_day:g} d（0.1 PVI）",
        fontsize=8.3, color="#666666", ha="right", va="center",
    )
    save_figure(figure, output / "04_scw_style_matched_overview")


def write_provenance(
    output: Path,
    baseline_dir: Path,
    scw_dir: Path,
    baseline: pd.DataFrame,
    scw: pd.DataFrame,
    final_day: float,
) -> None:
    source = pd.DataFrame(
        {
            "input_index": scw["input_index"],
            "i": scw["input_index"] % NX,
            "j": scw["input_index"] // NX,
            "x_m": ((scw["input_index"] % NX) + 0.5) * LX_M / NX,
            "y_m": ((scw["input_index"] // NX) + 0.5) * LY_M / NY,
        }
    )
    for prefix, frame in (("baseline", baseline), ("scw", scw)):
        for field in (
            "pressure_MPa",
            "pressure_change_MPa",
            "liquid_saturation",
            "oil_saturation_change_pp",
            "water_saturation",
            "water_saturation_change_pp",
            "liquid_x_0_CO2",
        ):
            source[f"{prefix}_{field}"] = frame[field].to_numpy(float)
    source.to_csv(output / "figure_source_data.csv", index=False)

    manifest = {
        "audience_medium": "provisional general scientific figure for local review",
        "publisher_requirements": "not specified; no journal-compliance claim",
        "grid": {"nx": NX, "ny": NY, "lx_m": LX_M, "ly_m": LY_M},
        "final_time_day": final_day,
        "final_pvi": 0.1,
        "sources": {
            "baseline": str((baseline_dir / "solution_final.csv").resolve()),
            "scw": str((scw_dir / "solution_final.csv").resolve()),
        },
        "conditions": {
            "baseline": {"temperature_K": BASE_TEMPERATURE_K, "initial_pressure_Pa": BASE_PRESSURE_PA},
            "scw": {"temperature_K": SCW_TEMPERATURE_K, "initial_pressure_Pa": SCW_PRESSURE_PA},
        },
        "transformations": [
            "sort by input_index and reshape without interpolation to 20 rows x 60 columns",
            "pressure Pa converted to MPa",
            "pressure change computed relative to each case's own uniform initial pressure",
            "saturation change computed relative to initial saturation 0.20 or 0.80 and converted to percentage points",
            "SCW-minus-baseline panels are direct cellwise subtraction",
            "no smoothing, filtering, clipping, exclusions, binning, or normalization of raw values",
        ],
        "missing_values": "none in required final fields",
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
        },
        "outputs": [
            "01_scw_final_fields.png/pdf/svg",
            "02_baseline_vs_scw.png/pdf/svg",
            "03_centerline_profiles.png/pdf/svg",
            "04_scw_style_matched_overview.png/pdf/svg",
            "figure_source_data.csv",
        ],
    }
    (output / "figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output / "ALT_TEXT.md").write_text(
        """# Figure descriptions

## 01_scw_final_fields
Four cell-resolved maps at 365.25 days for the 653.15 K, 28 MPa case. Pressure
is nearly uniform around 27.56 MPa. The strongest changes occur near the left
CO2 injector: oil-phase CO2 increases, oil saturation increases locally, and
supercritical-water saturation decreases. The producer is at the far right.

## 02_baseline_vs_scw
Three-by-three comparison of the 333.15 K, 5.16 MPa baseline and the 653.15 K,
28 MPa SCW case. Columns show pressure change relative to each case's own
initial pressure, water/SCW saturation change, and oil-phase CO2 mole fraction.
The bottom row is the direct SCW-minus-baseline difference on a zero-centered
scale. Spatial patterns remain injector-centered, but magnitudes differ.

## 03_centerline_profiles
Line profiles along the injector-producer grid row. Solid circles show the
baseline and dashed squares show the SCW case, so model identity does not rely
on color alone. The largest differences are concentrated near the injector.

## 04_scw_style_matched_overview
Blue-header, rounded-card dashboard matching the earlier report style. The left
card shows four cell-resolved fields for the 653.15 K, 28 MPa case over the
injector-front window x=0-120 m without interpolation. The right card compares
the low-temperature baseline and SCW conditions along the complete 300 m
injector-producer row using both color and line style.
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--scw", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    configure_style()
    args.output.mkdir(parents=True, exist_ok=True)
    baseline, baseline_day = read_case(args.baseline, BASE_PRESSURE_PA)
    scw, scw_day = read_case(args.scw, SCW_PRESSURE_PA)
    if abs(baseline_day - scw_day) > 1.0e-8:
        raise ValueError(f"final times differ: baseline={baseline_day}, scw={scw_day}")
    final_field_figure(scw, scw_day, args.output)
    comparison_figure(baseline, scw, scw_day, args.output)
    centerline_figure(baseline, scw, scw_day, args.output)
    style_matched_overview(baseline, scw, scw_day, args.output)
    write_provenance(args.output, args.baseline, args.scw, baseline, scw, scw_day)
    print(f"wrote validated SCW figures to {args.output}")


if __name__ == "__main__":
    main()
