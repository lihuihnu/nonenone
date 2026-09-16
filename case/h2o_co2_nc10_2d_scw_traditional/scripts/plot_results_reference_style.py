#!/usr/bin/env python3
"""Plot the 653.15 K, 28 MPa case in the established RESULTS figure style."""

from __future__ import annotations

import argparse
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
BASE_PRESSURE_MPA = 5.16
SCW_PRESSURE_MPA = 28.0
RASTER_DPI = 600

TEXT_COLOR = "#233746"
AXIS_COLOR = "#7890A0"
BASE_COLOR = "#0072B2"
SCW_COLOR = "#D55E00"

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
            "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "mathtext.rm": "STIXGeneral",
            "mathtext.it": "STIXGeneral:italic",
            "mathtext.bf": "STIXGeneral:bold",
            "mathtext.default": "it",
            "font.size": 12.5,
            "axes.titlesize": 13.5,
            "axes.titleweight": "normal",
            "axes.titlecolor": TEXT_COLOR,
            "axes.labelsize": 12.5,
            "axes.labelcolor": TEXT_COLOR,
            "axes.edgecolor": AXIS_COLOR,
            "axes.linewidth": 0.85,
            "xtick.labelsize": 11.0,
            "ytick.labelsize": 11.0,
            "xtick.color": TEXT_COLOR,
            "ytick.color": TEXT_COLOR,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 3.6,
            "ytick.major.size": 3.6,
            "xtick.major.width": 0.85,
            "ytick.major.width": 0.85,
            "legend.fontsize": 11.0,
            "legend.handlelength": 2.7,
            "figure.dpi": 180,
            "savefig.dpi": RASTER_DPI,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
            "savefig.pad_inches": 0.04,
        }
    )


def self_color_colormap() -> ListedColormap:
    anchor_index = np.arange(len(SELF_COLOR_ANCHORS), dtype=float)
    target_index = np.linspace(anchor_index[0], anchor_index[-1], 256)
    colors = np.column_stack(
        [
            PchipInterpolator(anchor_index, SELF_COLOR_ANCHORS[:, channel])(target_index)
            for channel in range(3)
        ]
    )
    return ListedColormap(np.clip(colors, 0.0, 1.0), name="self_color")


def read_case(directory: Path, initial_pressure_mpa: float) -> tuple[pd.DataFrame, float]:
    solution_path = directory / "solution_final.csv"
    summary_path = directory / "simulation_summary.csv"
    frame = pd.read_csv(solution_path).sort_values("input_index").reset_index(drop=True)
    if len(frame) != NX * NY:
        raise ValueError(f"{solution_path}: expected {NX * NY} cells, found {len(frame)}")
    if not np.array_equal(frame["input_index"].to_numpy(), np.arange(NX * NY)):
        raise ValueError(f"{solution_path}: incomplete input_index")
    required = {
        "pressure_Pa",
        "liquid_saturation",
        "vapor_saturation",
        "water_saturation",
        "liquid_x_0_CO2",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{solution_path}: missing {sorted(missing)}")
    values = frame[list(required)].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"{solution_path}: non-finite values")
    closure = (
        frame["liquid_saturation"]
        + frame["vapor_saturation"]
        + frame["water_saturation"]
    )
    if float(np.max(np.abs(closure - 1.0))) > 5.0e-8:
        raise ValueError(f"{solution_path}: saturation closure exceeds 5e-8")

    frame["x_m"] = (frame["input_index"] % NX + 0.5) * DX_M
    frame["y_m"] = (frame["input_index"] // NX + 0.5) * DY_M
    frame["pressure_mpa"] = frame["pressure_Pa"] * 1.0e-6
    frame["pressure_change_mpa"] = frame["pressure_mpa"] - initial_pressure_mpa
    frame["oil_saturation"] = frame["liquid_saturation"]
    frame["gas_saturation"] = frame["vapor_saturation"]
    frame["scw_saturation"] = frame["water_saturation"]
    frame["oil_x_co2"] = frame["liquid_x_0_CO2"]
    summary = pd.read_csv(summary_path)
    return frame, float(summary["final_time_day"].iloc[-1])


def add_wells(axis: plt.Axes) -> None:
    axis.scatter(
        [0.5 * DX_M], [WELL_Y_M], marker="^", s=28,
        facecolors="white", edgecolors="#2B8CBE", linewidths=1.0, zorder=4,
    )
    axis.scatter(
        [LX_M - 0.5 * DX_M], [WELL_Y_M], marker="v", s=28,
        facecolors="#173F5F", edgecolors="white", linewidths=0.7, zorder=4,
    )


def style_box_axis(axis: plt.Axes, *, linewidth: float = 0.9, tick_size: float | None = None) -> None:
    axis.set_facecolor("white")
    axis.grid(False)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("#6F7F8A")
        spine.set_linewidth(linewidth)
    axis.tick_params(
        direction="out", top=False, right=False, width=0.8, length=3.5,
        labelsize=tick_size,
    )


def save_figure(figure: plt.Figure, stem: Path) -> None:
    figure.savefig(stem.with_suffix(".png"), dpi=RASTER_DPI, facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), dpi=RASTER_DPI, facecolor="white")
    figure.savefig(stem.with_suffix(".svg"), dpi=RASTER_DPI, facecolor="white")
    plt.close(figure)
    with Image.open(stem.with_suffix(".png")) as raster:
        raster.convert("RGB").save(
            stem.with_suffix(".png"), dpi=(RASTER_DPI, RASTER_DPI), optimize=True
        )


def plot_field_map(
    frame: pd.DataFrame,
    output: Path,
    cmap: ListedColormap,
    *,
    field: str,
    filename: str,
    label: str,
    front_focused: bool,
    gamma: float | None = None,
) -> dict[str, float | str]:
    field_values = frame[field].to_numpy(dtype=float)
    vmin = float(np.min(field_values))
    vmax = float(np.max(field_values))
    if vmin == vmax:
        vmax = float(np.nextafter(vmax, np.inf))
    norm: Normalize = (
        PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)
        if gamma is not None
        else Normalize(vmin=vmin, vmax=vmax)
    )
    figure_size = (7.5, 6.15) if front_focused else (7.5, 3.05)
    figure, axis = plt.subplots(figsize=figure_size, layout="constrained")
    image = axis.imshow(
        field_values.reshape(NY, NX), origin="lower",
        extent=(0.0, LX_M, 0.0, LY_M), aspect="equal",
        interpolation="nearest", cmap=cmap, norm=norm,
    )
    axis.set_title("653.15 K, 28 MPa", loc="left", pad=4)
    axis.set_xlabel("x (m)", labelpad=1)
    axis.set_ylabel("y (m)", labelpad=1)
    axis.set_yticks([0, 50, 100])
    if front_focused:
        axis.set_xlim(0.0, 100.0)
        axis.set_xticks([0, 25, 50, 75, 100])
    else:
        axis.set_xticks([0, 100, 200, 300])
    add_wells(axis)
    colorbar = figure.colorbar(image, ax=axis, location="right", shrink=0.92, pad=0.015)
    colorbar.set_label(f"{label}  [$\\gamma$={gamma:g}]" if gamma is not None else label)
    colorbar.outline.set_visible(False)
    save_figure(figure, output / filename)
    return {
        "field": field,
        "minimum": vmin,
        "maximum": float(np.max(field_values)),
        "mean": float(np.mean(field_values)),
        "display_gamma": 1.0 if gamma is None else gamma,
        "display_x_max_m": 100.0 if front_focused else 300.0,
    }


def plot_centerline_profiles(
    baseline: pd.DataFrame,
    scw: pd.DataFrame,
    output: Path,
) -> None:
    specs = (
        ("pressure_change_mpa", r"$\Delta p$ (MPa)", (260.0, 300.0), [0.10, 0.08, 0.46, 0.40]),
        ("oil_saturation", r"$S_\mathrm{o}$", (35.0, 115.0), [0.49, 0.08, 0.47, 0.44]),
        ("scw_saturation", r"$S_\mathrm{w}$", (35.0, 115.0), [0.49, 0.08, 0.47, 0.44]),
        ("oil_x_co2", r"$x^\mathrm{o}_{\mathrm{CO_2}}$", (35.0, 115.0), [0.49, 0.08, 0.47, 0.44]),
    )
    styles = (
        (baseline, "333.15 K, 5.16 MPa", BASE_COLOR, "-"),
        (scw, "653.15 K, 28 MPa", SCW_COLOR, (0.0, (5.5, 2.2))),
    )
    center = slice(WELL_ROW * NX, (WELL_ROW + 1) * NX)
    figure, axes = plt.subplots(2, 2, figsize=(8.6, 6.1))
    figure.subplots_adjust(
        left=0.09, right=0.985, bottom=0.16, top=0.965,
        hspace=0.18, wspace=0.14,
    )
    for index, (axis, (field, label, zoom_window, inset_bounds)) in enumerate(
        zip(axes.flat, specs)
    ):
        pooled_zoom: list[float] = []
        for frame, _, color, linestyle in styles:
            row = frame.iloc[center]
            x = row["x_m"].to_numpy(dtype=float)
            y = row[field].to_numpy(dtype=float)
            axis.plot(
                x, y, color=color, linestyle=linestyle, linewidth=2.0,
                solid_capstyle="round", dash_capstyle="round",
            )
        axis.set_title(label, loc="left")
        axis.set_xlim(0.0, LX_M)
        axis.set_xticks([0.0, 100.0, 200.0, 300.0])
        axis.set_xlabel("x (m)" if index >= 2 else "")
        axis.set_box_aspect(0.62)
        style_box_axis(axis)
        if index < 2:
            axis.tick_params(labelbottom=False)

        inset = axis.inset_axes(inset_bounds)
        for frame, _, color, linestyle in styles:
            row = frame.iloc[center]
            x = row["x_m"].to_numpy(dtype=float)
            y = row[field].to_numpy(dtype=float)
            inset.plot(
                x, y, color=color, linestyle=linestyle, linewidth=1.25,
                solid_capstyle="round", dash_capstyle="round",
            )
            mask = (x >= zoom_window[0]) & (x <= zoom_window[1]) & np.isfinite(y)
            pooled_zoom.extend(y[mask])
        zoom_values = np.asarray(pooled_zoom, dtype=float)
        y_min, y_max = float(np.min(zoom_values)), float(np.max(zoom_values))
        span = y_max - y_min
        pad = 0.08 * span if span > 0.0 else 0.02 * max(abs(y_min), 1.0)
        inset.set_xlim(*zoom_window)
        inset.set_ylim(y_min - pad, y_max + pad)
        inset.set_xticks([zoom_window[0], 0.5 * sum(zoom_window), zoom_window[1]])
        inset.yaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
        style_box_axis(inset, linewidth=0.75, tick_size=6.3)

    handles = [
        Line2D(
            [0], [0], color=color, linestyle=linestyle, linewidth=2.2,
            solid_capstyle="round", dash_capstyle="round", label=label,
        )
        for _, label, color, linestyle in styles
    ]
    figure.legend(
        handles=handles, loc="lower center", ncol=2, frameon=False,
        bbox_to_anchor=(0.53, 0.015), handlelength=2.9,
    )
    save_figure(figure, output / "17_centerline_core_profiles")


def write_supporting_files(
    output: Path,
    reference: Path,
    baseline_dir: Path,
    scw_dir: Path,
    baseline: pd.DataFrame,
    scw: pd.DataFrame,
    final_day: float,
    field_metrics: list[dict[str, float | str]],
) -> None:
    source = pd.DataFrame(
        {
            "input_index": scw["input_index"],
            "x_m": scw["x_m"],
            "y_m": scw["y_m"],
        }
    )
    for prefix, frame in (("baseline", baseline), ("scw", scw)):
        for field in (
            "pressure_mpa", "pressure_change_mpa", "oil_saturation",
            "gas_saturation", "scw_saturation", "oil_x_co2",
        ):
            source[f"{prefix}_{field}"] = frame[field].to_numpy(dtype=float)
    source.to_csv(output / "figure_source_data.csv", index=False)
    pd.DataFrame(field_metrics).to_csv(output / "field_metrics.csv", index=False)

    manifest = {
        "reference_directory": str(reference.resolve()),
        "source_files": {
            "baseline": str((baseline_dir / "solution_final.csv").resolve()),
            "scw": str((scw_dir / "solution_final.csv").resolve()),
        },
        "conditions": {
            "baseline": {"temperature_K": 333.15, "initial_pressure_MPa": 5.16},
            "scw": {"temperature_K": 653.15, "initial_pressure_MPa": 28.0},
        },
        "final_time_day": final_day,
        "final_pvi": 0.1,
        "grid": {"nx": NX, "ny": NY, "lx_m": LX_M, "ly_m": LY_M},
        "transformations": [
            "sort cells by input_index and reshape to 20 rows x 60 columns",
            "pressure converted from Pa to MPa",
            "pressure centerline expressed relative to each case's own initial pressure",
            "Traditional liquid/vapor mapped to oil/gas",
            "field maps use nearest-neighbor pixels with no spatial interpolation",
            "pressure map retains x=0-300 m; transport maps display x=0-100 m",
            "all color limits are calculated from the complete 1200-cell SCW field",
            "oil-phase CO2 map uses labeled PowerNorm gamma=0.65 without clipping",
            "centerline panels retain x=0-300 m and use unsmoothed local insets",
        ],
        "style": {
            "basis": "D:/DOCUMENTS/ChatGPT/MPMC_SCW_V2/RESULTS",
            "colormap": "exact six-anchor MRST self_color palette with 256-color PCHIP interpolation",
            "font": "Times New Roman/STIX serif fallback",
            "injector_marker": "open cyan upward triangle",
            "producer_marker": "filled navy downward triangle",
            "raster_dpi": RASTER_DPI,
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
        },
    }
    (output / "figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output / "README.md").write_text(
        """# 653.15 K, 28 MPa SCW figures in the established RESULTS style

The figures reproduce the visual rules used by the sibling `RESULTS` directory:
the exact MRST `self_color` palette, Times/STIX typography, compact symbolic
labels, cyan/navy well markers, nearest-neighbor cell display, 600 dpi PNG, and
parallel PDF/SVG exports. Pressure retains the full 300 m domain; transport maps
focus on x=0-100 m while computing color limits from all 1200 cells. The
centerline comparison retains the complete 300 m domain.

Only the Traditional independent-water SCW formulation is available for this
supercritical run. Fields that this formulation defines structurally as zero,
including aqueous CO2 dissolution and oil-phase H2O, are not presented as
informative maps.
""",
        encoding="utf-8",
    )
    (output / "ALT_TEXT.md").write_text(
        """# Figure descriptions

## Field maps
Cell-resolved pressure, oil saturation, supercritical-water saturation, and
oil-phase CO2 mole fraction for the 653.15 K, 28 MPa run at 365.25 days. The
pressure map spans the complete injector-producer domain. Transport maps focus
on the first 100 m, where the displacement front occurs. Upward and downward
triangles identify injector and producer locations.

## 17_centerline_core_profiles
Four centerline panels compare the 333.15 K, 5.16 MPa baseline with the 653.15 K,
28 MPa SCW run. Panels show pressure change relative to each case's initial
pressure, oil saturation, SCW/water saturation, and oil-phase CO2 mole fraction.
Solid blue and dashed orange lines provide redundant condition encoding, and
insets show unsmoothed local detail.
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--scw", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)
    map_output = args.output / "field_maps"
    map_output.mkdir(parents=True, exist_ok=True)
    baseline, baseline_day = read_case(args.baseline, BASE_PRESSURE_MPA)
    scw, scw_day = read_case(args.scw, SCW_PRESSURE_MPA)
    if abs(baseline_day - scw_day) > 1.0e-10:
        raise ValueError(f"final times differ: {baseline_day} versus {scw_day}")

    cmap = self_color_colormap()
    metrics = [
        plot_field_map(
            scw, map_output, cmap, field="pressure_mpa", filename="01_pressure",
            label=r"$p$ (MPa)", front_focused=False,
        ),
        plot_field_map(
            scw, map_output, cmap, field="oil_saturation", filename="03_oil_saturation",
            label=r"$S_\mathrm{o}$", front_focused=True,
        ),
        plot_field_map(
            scw, map_output, cmap, field="scw_saturation", filename="05_scw_saturation",
            label=r"$S_\mathrm{w}$", front_focused=True,
        ),
        plot_field_map(
            scw, map_output, cmap, field="oil_x_co2", filename="10_oil_x_co2",
            label=r"$x^\mathrm{o}_{\mathrm{CO_2}}$", front_focused=True, gamma=0.65,
        ),
    ]
    plot_centerline_profiles(baseline, scw, args.output)
    write_supporting_files(
        args.output, args.reference, args.baseline, args.scw,
        baseline, scw, scw_day, metrics,
    )
    print(f"wrote RESULTS-style figures to {args.output}")


if __name__ == "__main__":
    main()
