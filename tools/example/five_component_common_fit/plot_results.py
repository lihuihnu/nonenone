#!/usr/bin/env python3
"""Draw bilingual single-axis figures for the common-data EOS comparison."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


MODELS = ("pr", "sw", "cpa")
LABEL = {"pr": "PR", "sw": "SW", "cpa": "CPA"}
COLOR = {"pr": "#0072B2", "sw": "#009E73", "cpa": "#D55E00"}
MARKER = {"pr": "o", "sw": "s", "cpa": "^"}
PAIR_LABEL = {
    "h2o_co2": r"H$_2$O–CO$_2$",
    "h2o_ch4": r"H$_2$O–CH$_4$",
    "co2_ch4": r"CO$_2$–CH$_4$",
    "co2_c2": r"CO$_2$–C$_2$H$_6$",
    "ch4_c2": r"CH$_4$–C$_2$H$_6$",
    "ch4_nc4": r"CH$_4$–n-C$_4$H$_{10}$",
}
PAIR_ORDER = tuple(PAIR_LABEL)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write an empty source-data table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style_axes(ax: Any) -> None:
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
        spine.set_color("#222222")
    ax.tick_params(direction="in", length=4, width=0.8, top=True, right=True)


def save(fig: Any, directory: Path, stem: str) -> list[str]:
    directory.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix in ("png", "pdf"):
        path = directory / f"{stem}.{suffix}"
        fig.savefig(path, dpi=400 if suffix == "png" else None, facecolor="white")
        outputs.append(str(path))
    plt.close(fig)
    return outputs


def validation_error(summary: list[dict[str, str]], directory: Path, zh: bool) -> list[str]:
    fig, ax = plt.subplots(figsize=(7.0, 4.3), layout="constrained")
    x = list(range(len(PAIR_ORDER)))
    offsets = {"pr": -0.16, "sw": 0.0, "cpa": 0.16}
    for model in MODELS:
        lookup = {
            row["pair"]: float(row["mean_abs_percent_error"])
            for row in summary
            if row["model"] == model and row["split"] == "validation"
        }
        ax.scatter(
            [value + offsets[model] for value in x],
            [lookup[pair] for pair in PAIR_ORDER],
            color=COLOR[model], marker=MARKER[model], s=45,
            edgecolors="#222222", linewidths=0.45, label=LABEL[model], zorder=3,
        )
    ax.set_yscale("log")
    ax.set_xticks(x, [PAIR_LABEL[pair] for pair in PAIR_ORDER], rotation=18, ha="right")
    ax.set_ylabel(
        "平均绝对归一化残差 [%]"
        if zh else "Mean absolute normalized residual [%]"
    )
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.legend(loc="upper right", frameon=True, edgecolor="#555555")
    style_axes(ax)
    return save(fig, directory, "01_validation_error")


def density_parity(predictions: list[dict[str, str]], directory: Path, zh: bool) -> list[str]:
    rows = [
        row for row in predictions
        if row["pair"] == "co2_c2" and row["split"] == "validation"
        and row["quantity"] == "mass_density_kg_m3"
    ]
    fig, ax = plt.subplots(figsize=(5.3, 4.4), layout="constrained")
    errors = [
        [100.0 * float(row["relative_error"]) for row in rows if row["model"] == model]
        for model in MODELS
    ]
    boxes = ax.boxplot(
        errors, patch_artist=True, widths=0.52, showfliers=True,
        medianprops={"color": "#222222", "linewidth": 1.5},
        whiskerprops={"color": "#444444", "linewidth": 1.0},
        capprops={"color": "#444444", "linewidth": 1.0},
        flierprops={"marker": "o", "markersize": 3.5, "markerfacecolor": "white",
                    "markeredgecolor": "#444444"},
    )
    for patch, model in zip(boxes["boxes"], MODELS, strict=True):
        patch.set_facecolor(COLOR[model])
        patch.set_alpha(0.78)
        patch.set_edgecolor("#222222")
    ax.axhline(0.0, color="#777777", linewidth=0.8)
    ax.set_xticks(range(1, len(MODELS) + 1), [LABEL[model] for model in MODELS])
    ax.set_ylabel("密度相对误差 [%]" if zh else "Density relative error [%]")
    style_axes(ax)
    return save(fig, directory, "02_density_parity")


def fitted_dry_bips(parameters: list[dict[str, str]], directory: Path, zh: bool) -> list[str]:
    dry = ("co2_ch4", "co2_c2", "ch4_c2", "ch4_nc4")
    fig, ax = plt.subplots(figsize=(6.5, 4.1), layout="constrained")
    x = list(range(len(dry)))
    offsets = {"pr": -0.14, "sw": 0.0, "cpa": 0.14}
    for model in MODELS:
        lookup = {
            row["pair"]: float(row["fitted"])
            for row in parameters if row["model"] == model
        }
        ax.scatter(
            [value + offsets[model] for value in x],
            [lookup[pair] for pair in dry], color=COLOR[model],
            marker=MARKER[model], s=48, edgecolors="#222222", linewidths=0.45,
            label=LABEL[model], zorder=3,
        )
    ax.axhline(0.0, color="#777777", linewidth=0.8)
    ax.set_xticks(x, [PAIR_LABEL[pair] for pair in dry], rotation=12, ha="right")
    ax.set_ylabel(r"拟合二元相互作用参数 $k_{ij}$" if zh else r"Fitted binary interaction $k_{ij}$")
    ax.legend(loc="upper right", frameon=True, edgecolor="#555555")
    style_axes(ax)
    return save(fig, directory, "03_fitted_dry_bips")


def initial_saturation(flow_root: Path, directory: Path, zh: bool) -> list[str]:
    phases = (("s_o_avg", "Oil", "#4477AA"), ("s_g_avg", "Gas", "#EE7733"),
              ("s_w_avg", "Water", "#66CCEE"))
    data: dict[str, dict[str, float]] = {}
    for model in MODELS:
        first = read_csv(flow_root / model / "reservoir_diagnostics.csv")[0]
        data[model] = {field: float(first[field]) for field, _, _ in phases}
    fig, ax = plt.subplots(figsize=(5.8, 4.2), layout="constrained")
    width = 0.23
    centers = list(range(len(MODELS)))
    for index, (field, label, color) in enumerate(phases):
        ax.bar(
            [center + (index - 1) * width for center in centers],
            [data[model][field] for model in MODELS],
            width=width, color=color, edgecolor="#222222", linewidth=0.6, label=label,
        )
    ax.set_xticks(centers, [LABEL[model] for model in MODELS])
    ax.set_ylabel("初始相饱和度" if zh else "Initial phase saturation")
    ax.set_ylim(0.0, 0.9)
    ax.legend(loc="upper right", frameon=True, edgecolor="#555555")
    style_axes(ax)
    return save(fig, directory, "04_initial_saturation")


def producer_mass_rate(flow_root: Path, directory: Path, zh: bool) -> list[str]:
    fig, ax = plt.subplots(figsize=(6.2, 4.2), layout="constrained")
    line_style = {"pr": "-", "sw": "--", "cpa": "-."}
    marker_start = {"pr": 1, "sw": 2, "cpa": 0}
    for model in MODELS:
        rows = [
            row for row in read_csv(flow_root / model / "well_history.csv")
            if row["name"] == "PROD"
        ]
        ax.plot(
            [float(row["time"]) for row in rows],
            [-float(row["m_total"]) for row in rows],
            color=COLOR[model], linestyle=line_style[model], marker=MARKER[model],
            markevery=(marker_start[model], 3),
            markersize=5.6, linewidth=1.8, label=LABEL[model],
        )
    ax.set_xlabel("时间 [天]" if zh else "Time [day]")
    ax.set_ylabel("生产井质量流量 [kg/s]" if zh else "Producer mass rate [kg/s]")
    ax.legend(loc="upper right", frameon=True, edgecolor="#555555")
    style_axes(ax)
    return save(fig, directory, "05_producer_mass_rate")


def write_source_data(
    summary: list[dict[str, str]],
    parameters: list[dict[str, str]],
    predictions: list[dict[str, str]],
    flow_root: Path,
    output_dir: Path,
) -> list[str]:
    source_dir = output_dir / "source_data"
    density = [
        row for row in predictions
        if row["pair"] == "co2_c2" and row["split"] == "validation"
        and row["quantity"] == "mass_density_kg_m3"
    ]
    saturation: list[dict[str, Any]] = []
    producer: list[dict[str, Any]] = []
    performance: list[dict[str, Any]] = []
    for model in MODELS:
        diagnostics = read_csv(flow_root / model / "reservoir_diagnostics.csv")
        for location, row in (("initial", diagnostics[0]), ("final", diagnostics[-1])):
            saturation.append({
                "model": LABEL[model],
                "state": location,
                "time_day": row["time"],
                "oil_saturation": row["s_o_avg"],
                "gas_saturation": row["s_g_avg"],
                "water_saturation": row["s_w_avg"],
            })
        for row in read_csv(flow_root / model / "well_history.csv"):
            if row["name"] == "PROD":
                producer.append({
                    "model": LABEL[model],
                    "time_day": row["time"],
                    "producer_mass_rate_kg_s": -float(row["m_total"]),
                })
        simulation = read_csv(flow_root / model / "simulation_summary.csv")[-1]
        mass_balance = read_csv(flow_root / model / "component_mass_balance.csv")
        performance.append({
            "model": LABEL[model],
            "accepted_steps": simulation["accepted_internal_steps"],
            "rejected_steps": simulation["rejected_internal_steps"],
            "snes_iterations": simulation["attempted_snes"],
            "ksp_iterations": simulation["attempted_ksp"],
            "simulation_wall_s": simulation["simulation_loop_wall"],
            "max_abs_mass_balance_relative_error": max(
                abs(float(row["relative_error"])) for row in mass_balance
            ),
        })

    tables = {
        "fit_summary.csv": summary,
        "fitted_parameters.csv": parameters,
        "density_validation.csv": density,
        "flow_saturation.csv": saturation,
        "producer_mass_rate.csv": producer,
        "flow_performance.csv": performance,
    }
    outputs = []
    for name, rows in tables.items():
        path = source_dir / name
        write_csv(path, rows)
        outputs.append(str(path))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit-results", type=Path, required=True)
    parser.add_argument("--flow-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.linewidth": 0.9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.transparent": False,
    })
    summary = read_csv(args.fit_results / "fit_summary.csv")
    parameters = read_csv(args.fit_results / "fitted_parameters.csv")
    predictions = read_csv(args.fit_results / "point_predictions.csv")
    outputs = []
    source_tables = write_source_data(
        summary, parameters, predictions, args.flow_root, args.output_dir
    )
    for language in ("en", "zh"):
        directory = args.output_dir / language
        zh = language == "zh"
        plt.rcParams["font.family"] = "Microsoft YaHei" if zh else "DejaVu Sans"
        outputs += validation_error(summary, directory, zh)
        outputs += density_parity(predictions, directory, zh)
        outputs += fitted_dry_bips(parameters, directory, zh)
        outputs += initial_saturation(args.flow_root, directory, zh)
        outputs += producer_mass_rate(args.flow_root, directory, zh)
    manifest = {
        "figures": outputs,
        "source_data": source_tables,
        "fit_results": str(args.fit_results),
        "flow_results": str(args.flow_root),
        "style": "single axes, no title, no grid, full frame, in-frame English legend",
        "formats": ["PNG 400 dpi", "vector PDF"],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
