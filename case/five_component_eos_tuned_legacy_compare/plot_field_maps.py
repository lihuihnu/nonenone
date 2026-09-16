#!/usr/bin/env python3
"""Plot vertically averaged day-30 fields for four flow-model packages."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_FULL = (
    ROOT.parent / "five_component_eos_tuned_compare" / "results"
    / "translated_density_30day_20260826"
)
DEFAULT_TRADITIONAL = ROOT / "results" / "translated_density_30day_20260826" / "pr"
DEFAULT_OUTPUT = ROOT / "figures" / "translated_density_30day_20260826" / "fields"

NX, NY, NZ = 20, 20, 5
LX, LY, LZ = 1000.0, 600.0, 50.0
DPI = 600
WIDTH_MM, HEIGHT_MM = 150.0, 92.0
FINAL_DAY = 30.0

MODELS = (
    ("our_pr", "Our PR", DEFAULT_FULL / "pr", False),
    ("our_sw", "Our SW", DEFAULT_FULL / "sw", False),
    ("our_cpa", "Our CPA", DEFAULT_FULL / "cpa", False),
    ("traditional", "Traditional", DEFAULT_TRADITIONAL, True),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-root", type=Path, default=DEFAULT_FULL)
    parser.add_argument("--traditional-root", type=Path, default=DEFAULT_TRADITIONAL)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def resolve(path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def relative_label(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.parent.parent).as_posix()
    except ValueError:
        return str(path.resolve())


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"missing source file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != NX * NY * NZ:
        raise ValueError(f"{path}: expected {NX * NY * NZ} rows, got {len(rows)}")
    indices = np.asarray([int(row["input_index"]) for row in rows])
    if not np.array_equal(indices, np.arange(NX * NY * NZ)):
        raise ValueError(f"{path}: input_index is not contiguous in k-j-i order")
    return rows


def volume(rows: list[dict[str, str]], column: str) -> np.ndarray:
    values = np.asarray([float(row[column]) for row in rows], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"non-finite values in {column}")
    return values.reshape(NZ, NY, NX)


def load_fields(full_root: Path, traditional_root: Path) -> dict[str, dict[str, object]]:
    roots = {
        "our_pr": (full_root / "pr", False),
        "our_sw": (full_root / "sw", False),
        "our_cpa": (full_root / "cpa", False),
        "traditional": (traditional_root, True),
    }
    result: dict[str, dict[str, object]] = {}
    for key, label, _, _ in MODELS:
        root, traditional = roots[key]
        initial = read_rows(root / "solution_step_0.csv")
        final = read_rows(root / "solution_final.csv")
        gas_saturation_column = "vapor_saturation" if traditional else "gas_saturation"
        gas_co2_column = "vapor_y_0_CO2" if traditional else "gas_y_1_CO2"
        pressure_change = (
            volume(final, "pressure_Pa") - volume(initial, "pressure_Pa")
        ).mean(axis=0) / 1.0e5
        water_saturation_change = 1.0e4 * (
            volume(final, "water_saturation")
            - volume(initial, "water_saturation")
        ).mean(axis=0)
        result[key] = {
            "label": label,
            "root": root,
            "pressure_change": pressure_change,
            "gas_saturation": volume(final, gas_saturation_column).mean(axis=0),
            "aqueous_co2": (
                np.zeros((NY, NX), dtype=float)
                if traditional
                else 1.0e3 * volume(final, "water_x_1_CO2").mean(axis=0)
            ),
            "gas_co2": volume(final, gas_co2_column).mean(axis=0),
            "water_saturation_change": water_saturation_change,
        }
    return result


def configure_fonts(language: str) -> None:
    family = "DejaVu Sans"
    if language == "zh":
        candidates = (
            Path("C:/Windows/Fonts/msyh.ttc"),
            Path("/mnt/c/Windows/Fonts/msyh.ttc"),
            Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"),
        )
        for candidate in candidates:
            if candidate.is_file():
                fm.fontManager.addfont(str(candidate))
                family = fm.FontProperties(fname=str(candidate)).get_name()
                break
        else:
            raise FileNotFoundError("Chinese font not found")
    mpl.rcParams.update({
        "font.family": family,
        "font.size": 9.4,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 8.8,
        "ytick.labelsize": 8.8,
        "axes.linewidth": 1.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.unicode_minus": False,
    })


def labels(language: str) -> dict[str, str]:
    if language == "zh":
        return {
            "x": "东西方向距离（m）",
            "y": "南北方向距离（m）",
            "pressure_change": "压力变化（bar）",
            "gas_saturation": "气相饱和度",
            "aqueous_co2": r"水相 CO$_2$ 摩尔分数（×10$^{-3}$）",
            "gas_co2": r"气相 CO$_2$ 摩尔分数",
            "water_saturation_change": r"水相饱和度变化（×10$^4$）",
        }
    return {
        "x": "West–east distance (m)",
        "y": "South–north distance (m)",
        "pressure_change": "Pressure change (bar)",
        "gas_saturation": "Gas saturation",
        "aqueous_co2": r"Aqueous CO$_2$ mole fraction (×10$^{-3}$)",
        "gas_co2": r"Gas-phase CO$_2$ mole fraction",
        "water_saturation_change": r"Water-saturation change (×10$^4$)",
    }


def pooled_limits(fields: dict[str, dict[str, object]], key: str) -> tuple[float, float]:
    arrays = [np.asarray(fields[model_key][key]) for model_key, _, _, _ in MODELS]
    return (
        min(float(np.min(array)) for array in arrays),
        max(float(np.max(array)) for array in arrays),
    )


def sequential_norm(fields: dict[str, dict[str, object]], key: str,
                    include_zero: bool = False) -> mpl.colors.Normalize:
    lower, upper = pooled_limits(fields, key)
    span = max(upper - lower, max(abs(lower), abs(upper), 1.0e-12) * 0.01)
    vmin = 0.0 if include_zero else lower - 0.04 * span
    vmax = upper + 0.04 * span
    return mpl.colors.Normalize(vmin=vmin, vmax=vmax)


def diverging_norm(fields: dict[str, dict[str, object]], key: str) -> mpl.colors.TwoSlopeNorm:
    lower, upper = pooled_limits(fields, key)
    if not lower < 0.0 < upper:
        raise ValueError(f"{key}: expected both negative and positive field values")
    return mpl.colors.TwoSlopeNorm(
        vmin=1.05 * lower, vcenter=0.0, vmax=1.05 * upper,
    )


def field_specs(fields: dict[str, dict[str, object]]) -> tuple[dict[str, object], ...]:
    return (
        {
            "key": "pressure_change",
            "number": "01",
            "cmap": "RdBu_r",
            "norm": diverging_norm(fields, "pressure_change"),
        },
        {
            "key": "gas_saturation",
            "number": "02",
            "cmap": "viridis",
            "norm": sequential_norm(fields, "gas_saturation"),
        },
        {
            "key": "aqueous_co2",
            "number": "03",
            "cmap": "cividis",
            "norm": sequential_norm(fields, "aqueous_co2", include_zero=True),
        },
        {
            "key": "gas_co2",
            "number": "04",
            "cmap": "magma",
            "norm": sequential_norm(fields, "gas_co2"),
        },
        {
            "key": "water_saturation_change",
            "number": "05",
            "cmap": "PuOr_r",
            "norm": diverging_norm(fields, "water_saturation_change"),
        },
    )


def add_wells(ax: plt.Axes) -> None:
    dx, dy = LX / NX, LY / NY
    injector = ((1 + 0.5) * dx, (4 + 0.5) * dy)
    producer = ((NX - 2 + 0.5) * dx, (NY - 5 + 0.5) * dy)
    ax.scatter(*injector, marker="*", s=105, facecolor="white", edgecolor="#111111",
               linewidth=1.0, zorder=4)
    ax.scatter(*producer, marker="X", s=64, facecolor="#111111", edgecolor="white",
               linewidth=0.8, zorder=4)
    annotation_box = {
        "facecolor": "white", "edgecolor": "none",
        "boxstyle": "round,pad=0.16", "alpha": 0.82,
    }
    ax.annotate("CO$_2$ INJ", injector, xytext=(7, 6), textcoords="offset points",
                ha="left", va="bottom", fontsize=7.8, color="#111111",
                bbox=annotation_box)
    ax.annotate("PROD", producer, xytext=(-7, -6), textcoords="offset points",
                ha="right", va="top", fontsize=7.8, color="#111111",
                bbox=annotation_box)


def export_field(data: np.ndarray, model_label: str, model_key: str,
                 spec: dict[str, object], language: str, output_root: Path,
                 force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = plt.subplots(
        figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), layout="constrained",
    )
    x_edges = np.linspace(0.0, LX, NX + 1)
    y_edges = np.linspace(0.0, LY, NY + 1)
    mesh = ax.pcolormesh(
        x_edges, y_edges, data, shading="flat", cmap=str(spec["cmap"]),
        norm=spec["norm"], rasterized=True,
    )
    add_wells(ax)
    ax.text(0.985, 0.965, model_label, transform=ax.transAxes, ha="right", va="top",
            fontsize=9.0, bbox={"facecolor": "white", "edgecolor": "#808080",
                                "boxstyle": "round,pad=0.25", "alpha": 0.92})
    ax.set(xlim=(0, LX), ylim=(0, LY), xlabel=text["x"], ylabel=text["y"],
           aspect="equal")
    ax.grid(False)
    ax.tick_params(direction="in", top=False, right=False, length=3.5, width=0.9)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)
        spine.set_color("#202020")
    colorbar = fig.colorbar(mesh, ax=ax, pad=0.025, fraction=0.046)
    colorbar.set_label(text[str(spec["key"])])

    target = output_root / language
    target.mkdir(parents=True, exist_ok=True)
    stem = f'{spec["number"]}_{spec["key"]}_{model_key}'
    paths = [target / f"{stem}.png", target / f"{stem}.pdf"]
    if not force and any(path.exists() for path in paths):
        raise FileExistsError(f"field figure exists; use --force: {paths[0]}")
    fig.savefig(paths[0], dpi=DPI, metadata={"Software": "Our field-map plotter"})
    fig.savefig(paths[1], dpi=DPI, metadata={
        "Creator": "Our field-map plotter",
        "Subject": "Vertically averaged day-30 reservoir field",
    })
    plt.close(fig)
    return paths


def write_source_table(fields: dict[str, dict[str, object]], output_root: Path,
                       force: bool) -> Path:
    path = output_root / "field_source_data.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        raise FileExistsError(f"field source table exists; use --force: {path}")
    columns = (
        "model", "i", "j", "x_center_m", "y_center_m",
        "pressure_change_bar", "gas_saturation",
        "aqueous_co2_mole_fraction_x1e3", "gas_co2_mole_fraction",
        "water_saturation_change_x1e4",
    )
    with path.open("w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for model_key, _, _, _ in MODELS:
            model = fields[model_key]
            for j in range(NY):
                for i in range(NX):
                    writer.writerow({
                        "model": model["label"],
                        "i": i,
                        "j": j,
                        "x_center_m": (i + 0.5) * LX / NX,
                        "y_center_m": (j + 0.5) * LY / NY,
                        "pressure_change_bar": model["pressure_change"][j, i],
                        "gas_saturation": model["gas_saturation"][j, i],
                        "aqueous_co2_mole_fraction_x1e3": model["aqueous_co2"][j, i],
                        "gas_co2_mole_fraction": model["gas_co2"][j, i],
                        "water_saturation_change_x1e4":
                            model["water_saturation_change"][j, i],
                    })
    return path


def main() -> None:
    args = parse_args()
    full_root = resolve(args.full_root)
    traditional_root = resolve(args.traditional_root)
    output_root = resolve(args.output_root)
    fields = load_fields(full_root, traditional_root)
    specs = field_specs(fields)
    generated: list[Path] = [write_source_table(fields, output_root, args.force)]
    for spec in specs:
        norm = spec["norm"]
        for model_key, _, _, _ in MODELS:
            data = fields[model_key][str(spec["key"])]
            data_min, data_max = float(np.min(data)), float(np.max(data))
            if data_min < norm.vmin or data_max > norm.vmax:
                raise ValueError(
                    f'{spec["key"]}/{model_key} range [{data_min}, {data_max}] '
                    f'exceeds shared color limits [{norm.vmin}, {norm.vmax}]'
                )
    for language in ("en", "zh"):
        configure_fonts(language)
        for spec in specs:
            for model_key, model_label, _, _ in MODELS:
                generated.extend(export_field(
                    fields[model_key][str(spec["key"])], model_label, model_key,
                    spec, language, output_root, args.force,
                ))

    manifest = {
        "purpose": "general provisional reservoir field figures; no target-journal compliance claimed",
        "source_roots": {
            "fully_compositional": relative_label(full_root),
            "traditional": relative_label(traditional_root),
        },
        "grid": {"cells": [NX, NY, NZ], "dimensions_m": [LX, LY, LZ]},
        "time_day": FINAL_DAY,
        "transformations": [
            "all displayed fields are unweighted arithmetic means over the five vertical layers",
            "pressure change is final pressure minus initial pressure, converted from Pa to bar",
            "water-saturation change is final minus initial and multiplied by 1e4",
            "aqueous CO2 mole fraction is multiplied by 1e3; Traditional is exactly zero by model definition",
            "no interpolation, smoothing, filtering, clipping, data normalization, or excluded cells",
        ],
        "shared_color_limits": {
            str(spec["key"]): (
                [spec["norm"].vmin, spec["norm"].vcenter, spec["norm"].vmax]
                if isinstance(spec["norm"], mpl.colors.TwoSlopeNorm)
                else [spec["norm"].vmin, spec["norm"].vmax]
            )
            for spec in specs
        },
        "figure_size_mm": [WIDTH_MM, HEIGHT_MM],
        "png_dpi": DPI,
        "generated": [relative_label(path) for path in generated],
    }
    manifest_path = output_root / "figure_manifest.json"
    if manifest_path.exists() and not args.force:
        raise FileExistsError(f"field manifest exists; use --force: {manifest_path}")
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(f"generated {len(generated) - 1} field figure files and one source table")


if __name__ == "__main__":
    main()
