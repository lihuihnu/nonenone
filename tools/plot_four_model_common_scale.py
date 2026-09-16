#!/usr/bin/env python3
"""Plot one common-scale comparison figure per 2D physical field.

The top row contains the four absolute fields on one shared color scale.  The
bottom row contains differences from New-PR on one zero-centered scale.  Raw
simulation CSV files are read only and are never modified.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize, TwoSlopeNorm
from PIL import Image
from scipy.interpolate import PchipInterpolator


NX = 60
NY = 20
LX_M = 300.0
LY_M = 100.0
REFERENCE_MODEL = "New-PR"

DEFAULT_SOURCE_ROOT = Path(
    r"C:\Users\ncu55\Documents\Codex\2026-08-29"
    r"\referenced-chatgpt-conversation-this-is-an\outputs"
    r"\benchmark_runs\stage_0p1"
)

RUN_DIRECTORIES = {
    "Traditional": "traditional_phase_role",
    "New-PR": "pr_0p1",
    "New-SW": "sw_0p1_phase_order_fix",
    "New-CPA": "cpa_0p1_optimized4",
}

# Exact 6 anchors from the user-supplied MRST self_color.m.
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

FIELDS = (
    {
        "key": "pressure",
        "label": "Pressure",
        "absolute_unit": "MPa",
        "difference_unit": "kPa",
        "absolute_scale": 1.0e-6,
        "difference_scale": 1.0e-3,
        "filename": "01_pressure_four_models",
    },
    {
        "key": "oil_saturation",
        "label": "Oil saturation",
        "absolute_unit": "fraction",
        "difference_unit": "percentage points",
        "absolute_scale": 1.0,
        "difference_scale": 100.0,
        "filename": "02_oil_saturation_four_models",
    },
    {
        "key": "gas_saturation",
        "label": "Gas saturation",
        "absolute_unit": "fraction",
        "difference_unit": "percentage points",
        "absolute_scale": 1.0,
        "difference_scale": 100.0,
        "filename": "03_gas_saturation_four_models",
    },
    {
        "key": "water_saturation",
        "label": "Water saturation",
        "absolute_unit": "fraction",
        "difference_unit": "percentage points",
        "absolute_scale": 1.0,
        "difference_scale": 100.0,
        "filename": "04_water_saturation_four_models",
    },
)


def self_color_colormap() -> ListedColormap:
    """Reproduce the user-supplied MATLAB self_color.m palette with PCHIP."""
    x = np.arange(len(SELF_COLOR_ANCHORS), dtype=float)
    xi = np.linspace(x[0], x[-1], 256)
    colors = np.column_stack(
        [PchipInterpolator(x, SELF_COLOR_ANCHORS[:, channel])(xi) for channel in range(3)]
    )
    return ListedColormap(np.clip(colors, 0.0, 1.0), name="self_color")


def read_final_day(run_directory: Path) -> float:
    summary_path = run_directory / "simulation_summary.csv"
    summary = pd.read_csv(summary_path)
    if summary.empty or "final_time_day" not in summary:
        raise ValueError(f"missing final_time_day in {summary_path}")
    return float(summary["final_time_day"].iloc[-1])


def read_run(source_root: Path, model: str) -> tuple[pd.DataFrame, float, Path]:
    run_directory = source_root / RUN_DIRECTORIES[model]
    source_path = run_directory / "solution_final.csv"
    frame = pd.read_csv(source_path).sort_values("input_index").reset_index(drop=True)
    if len(frame) != NX * NY:
        raise ValueError(f"{model}: expected {NX * NY} cells, found {len(frame)}")
    expected_indices = np.arange(NX * NY)
    if not np.array_equal(frame["input_index"].to_numpy(), expected_indices):
        raise ValueError(f"{model}: input_index is not a complete 0..{NX * NY - 1} sequence")

    frame["pressure"] = frame["pressure_Pa"]
    if "oil_saturation" in frame:
        frame["oil_saturation"] = frame["oil_saturation"]
        frame["gas_saturation"] = frame["gas_saturation"]
    else:
        frame["oil_saturation"] = frame["liquid_saturation"]
        frame["gas_saturation"] = frame["vapor_saturation"]
    frame["water_saturation"] = frame["water_saturation"]

    comparable = frame[[field["key"] for field in FIELDS]].to_numpy(dtype=float)
    if not np.isfinite(comparable).all():
        raise ValueError(f"{model}: non-finite values found in comparable fields")
    saturation_sum = (
        frame["oil_saturation"]
        + frame["gas_saturation"]
        + frame["water_saturation"]
    )
    if np.max(np.abs(saturation_sum - 1.0)) > 5.0e-8:
        raise ValueError(f"{model}: saturation closure exceeds tolerance")
    return frame, read_final_day(run_directory), source_path


def field_grid(frame: pd.DataFrame, key: str, scale: float) -> np.ndarray:
    return (frame[key].to_numpy(dtype=float) * scale).reshape(NY, NX)


def add_well_markers(ax: plt.Axes) -> None:
    y = (9 + 0.5) * LY_M / NY
    ax.scatter(
        [(0 + 0.5) * LX_M / NX],
        [y],
        marker="o",
        s=30,
        facecolors="none",
        edgecolors="black",
        linewidths=1.0,
        zorder=5,
    )
    ax.scatter(
        [(59 + 0.5) * LX_M / NX],
        [y],
        marker="x",
        s=30,
        color="black",
        linewidths=1.0,
        zorder=5,
    )


def format_number(value: float) -> str:
    magnitude = abs(value)
    if magnitude == 0.0:
        return "0"
    if magnitude < 1.0e-3 or magnitude >= 1.0e4:
        return f"{value:.2e}"
    if magnitude < 0.1:
        return f"{value:.4f}"
    return f"{value:.3f}"


def plot_field(
    data: dict[str, pd.DataFrame],
    field: dict[str, object],
    final_day: float,
    output_directory: Path,
    cmap: ListedColormap,
) -> list[dict[str, float | str]]:
    models = list(RUN_DIRECTORIES)
    key = str(field["key"])
    absolute_scale = float(field["absolute_scale"])
    difference_scale = float(field["difference_scale"])

    absolute = {
        model: field_grid(data[model], key, absolute_scale) for model in models
    }
    all_absolute = np.concatenate([absolute[model].ravel() for model in models])
    absolute_min = float(np.min(all_absolute))
    absolute_max = float(np.max(all_absolute))
    if absolute_min == absolute_max:
        absolute_max = np.nextafter(absolute_max, np.inf)
    absolute_norm = Normalize(vmin=absolute_min, vmax=absolute_max)

    reference = data[REFERENCE_MODEL][key].to_numpy(dtype=float)
    differences = {
        model: (
            (data[model][key].to_numpy(dtype=float) - reference) * difference_scale
        ).reshape(NY, NX)
        for model in models
    }
    difference_limit = max(
        float(np.max(np.abs(differences[model]))) for model in models
    )
    difference_limit = max(difference_limit, np.finfo(float).eps)
    difference_norm = TwoSlopeNorm(
        vmin=-difference_limit, vcenter=0.0, vmax=difference_limit
    )

    fig, axes = plt.subplots(
        2,
        4,
        figsize=(16.0, 5.6),
        sharex=True,
        sharey=True,
        layout="constrained",
    )

    absolute_image = None
    difference_image = None
    contour_levels = np.linspace(absolute_min, absolute_max, 8)[1:-1]
    metrics: list[dict[str, float | str]] = []
    for column, model in enumerate(models):
        ax = axes[0, column]
        absolute_image = ax.imshow(
            absolute[model],
            origin="lower",
            extent=(0.0, LX_M, 0.0, LY_M),
            aspect="equal",
            interpolation="nearest",
            cmap=cmap,
            norm=absolute_norm,
        )
        if len(contour_levels):
            ax.contour(
                np.linspace(LX_M / (2 * NX), LX_M - LX_M / (2 * NX), NX),
                np.linspace(LY_M / (2 * NY), LY_M - LY_M / (2 * NY), NY),
                absolute[model],
                levels=contour_levels,
                colors="black",
                linewidths=0.35,
                alpha=0.35,
            )
        add_well_markers(ax)
        ax.set_title(model, fontweight="bold")
        ax.text(
            0.02,
            0.04,
            f"min {format_number(float(np.min(absolute[model])))}\n"
            f"max {format_number(float(np.max(absolute[model])))}",
            transform=ax.transAxes,
            fontsize=7.5,
            va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 2},
        )

        delta = differences[model]
        ax = axes[1, column]
        difference_image = ax.imshow(
            delta,
            origin="lower",
            extent=(0.0, LX_M, 0.0, LY_M),
            aspect="equal",
            interpolation="nearest",
            cmap=cmap,
            norm=difference_norm,
        )
        if np.min(delta) < 0.0 < np.max(delta):
            ax.contour(
                np.linspace(LX_M / (2 * NX), LX_M - LX_M / (2 * NX), NX),
                np.linspace(LY_M / (2 * NY), LY_M - LY_M / (2 * NY), NY),
                delta,
                levels=[0.0],
                colors="black",
                linewidths=0.55,
                alpha=0.65,
            )
        add_well_markers(ax)
        rmse = float(np.sqrt(np.mean(delta * delta)))
        max_abs = float(np.max(np.abs(delta)))
        ax.text(
            0.02,
            0.04,
            f"max |delta| {format_number(max_abs)}\nRMSE {format_number(rmse)}",
            transform=ax.transAxes,
            fontsize=7.5,
            va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 2},
        )
        metrics.append(
            {
                "field": key,
                "model": model,
                "absolute_min": float(np.min(absolute[model])),
                "absolute_max": float(np.max(absolute[model])),
                "difference_min_to_new_pr": float(np.min(delta)),
                "difference_max_to_new_pr": float(np.max(delta)),
                "difference_max_abs_to_new_pr": max_abs,
                "difference_rmse_to_new_pr": rmse,
            }
        )

    for ax in axes[1, :]:
        ax.set_xlabel("x (m)")
    axes[0, 0].set_ylabel("Absolute field\ny (m)")
    axes[1, 0].set_ylabel("Difference from New-PR\ny (m)")
    for ax in axes.flat:
        ax.set_xticks([0, 100, 200, 300])
        ax.set_yticks([0, 50, 100])

    assert absolute_image is not None and difference_image is not None
    absolute_colorbar = fig.colorbar(
        absolute_image, ax=axes[0, :], location="right", shrink=0.93, pad=0.015
    )
    absolute_colorbar.set_label(
        f"{field['label']} ({field['absolute_unit']})"
    )
    difference_colorbar = fig.colorbar(
        difference_image, ax=axes[1, :], location="right", shrink=0.93, pad=0.015
    )
    difference_colorbar.set_label(
        f"Delta vs New-PR ({field['difference_unit']})"
    )

    pvi = final_day / 3652.5
    fig.suptitle(
        f"{field['label']} at {final_day:g} d ({pvi:.3f} PVI): four-model comparison",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.004,
        "self_color: blue = low/negative, pale neutral = midpoint/zero, orange-red = high/positive; "
        "open circle = injector, x = producer",
        ha="center",
        fontsize=8,
        color="#333333",
    )

    stem = str(field["filename"])
    png_path = output_directory / f"{stem}.png"
    fig.savefig(png_path, dpi=300, facecolor="white")
    fig.savefig(output_directory / f"{stem}.pdf", facecolor="white")
    plt.close(fig)
    # Matplotlib writes an RGBA PNG even with a white figure face.  Flatten the
    # fully opaque raster to explicit RGB so downstream viewers cannot substitute
    # a different page background.
    with Image.open(png_path) as image:
        image.convert("RGB").save(png_path, dpi=(300, 300), optimize=True)
    return metrics


def write_provenance(
    output_directory: Path,
    source_root: Path,
    source_paths: dict[str, Path],
    final_day: float,
    metrics: list[dict[str, float | str]],
) -> None:
    pd.DataFrame(metrics).to_csv(output_directory / "comparison_metrics.csv", index=False)
    manifest = {
        "source_root": str(source_root.resolve()),
        "source_files": {model: str(path.resolve()) for model, path in source_paths.items()},
        "models": list(RUN_DIRECTORIES),
        "reference_model_for_differences": REFERENCE_MODEL,
        "grid": {"nx": NX, "ny": NY, "lx_m": LX_M, "ly_m": LY_M},
        "final_time_day": final_day,
        "final_pvi": final_day / 3652.5,
        "transformations": [
            "sort cells by input_index and reshape to 20 rows x 60 columns",
            "pressure: Pa to MPa for absolute maps and Pa to kPa for differences",
            "saturation differences: fraction to percentage points",
            "Traditional liquid/vapor saturation mapped to oil/gas saturation",
            "field pixels use nearest-neighbor rendering with no smoothing or cell removal",
            "thin contour guides are interpolated from the same cell-center values",
        ],
        "colormap": {
            "name": "self_color",
            "source": "user-supplied MRST self_color.m anchors embedded in this script",
            "user_source_sha256": "E9111C8E7228339E9AC446690BC298A7A4A8C32076220FD712FA782B8760D51B",
            "anchors_rgb_0_1": SELF_COLOR_ANCHORS.tolist(),
            "anchors_rgb_0_255": (SELF_COLOR_ANCHORS * 255.0).round().astype(int).tolist(),
            "interpolation": "PCHIP, 256 colors",
        },
        "normalization": {
            "absolute": "field-wise min and max pooled over all four models",
            "difference": "field-wise symmetric +/- maximum absolute difference over all four models",
        },
        "software": {
            "python": __import__("sys").version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
        },
    }
    (output_directory / "figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output_directory / "README.md").write_text(
        """# H2O-CO2-nC10 four-model field comparison

Each physical quantity has one figure. The upper row compares Traditional,
New-PR, New-SW, and New-CPA using the pooled minimum and maximum of all four
models. The lower row shows each model minus New-PR with one zero-centered,
symmetric difference range so small spatial discrepancies remain visible.

The four common fields are pressure, oil saturation, gas saturation, and water
saturation. PNG files are nominally 300 dpi; PDF files preserve vector text and contours.
The raw CSV files are unchanged. Exact sources, transformations, ranges, and the
MRST `self_color` anchors are recorded in `figure_manifest.json`; numeric
difference summaries are in `comparison_metrics.csv`.
""",
        encoding="utf-8",
    )


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.pad_inches": 0.08,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    configure_matplotlib()
    args.output.mkdir(parents=True, exist_ok=True)
    data: dict[str, pd.DataFrame] = {}
    final_days: dict[str, float] = {}
    source_paths: dict[str, Path] = {}
    for model in RUN_DIRECTORIES:
        data[model], final_days[model], source_paths[model] = read_run(
            args.source_root, model
        )
    if max(final_days.values()) - min(final_days.values()) > 1.0e-10:
        raise ValueError(f"model final times differ: {final_days}")
    final_day = next(iter(final_days.values()))

    cmap = self_color_colormap()
    metrics: list[dict[str, float | str]] = []
    for field in FIELDS:
        metrics.extend(plot_field(data, field, final_day, args.output, cmap))
    write_provenance(
        args.output, args.source_root, source_paths, final_day, metrics
    )
    print(f"wrote four common-scale field figures to {args.output}")


if __name__ == "__main__":
    main()
