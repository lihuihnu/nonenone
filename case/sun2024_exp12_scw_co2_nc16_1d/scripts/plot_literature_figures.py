#!/usr/bin/env python3
"""Redraw literature-based figures for the Sun et al. (2024) case report."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
from PIL import Image


CASE_DIR = Path(__file__).resolve().parents[1]
REPORT_DIR = CASE_DIR / "report"
DATA_DIR = REPORT_DIR / "data"
FIGURE_DIR = REPORT_DIR / "figures"

DOI = "10.3390/app14093588"
ACCESS_DATE = "2026-09-02"
MODEL_RECOVERY_PERCENT = 87.05100518460273
CRITICAL_TEMPERATURE_C = 373.946
CRITICAL_PRESSURE_MPA = 22.064

FAMILIES = ["Steam", "SCW", "SCW+N2", "SCW+CO2", "SCW+N2+CO2"]
COLORS = {
    "Steam": "#0072B2",
    "SCW": "#D55E00",
    "SCW+N2": "#009E73",
    "SCW+CO2": "#CC79A7",
    "SCW+N2+CO2": "#000000",
}
MARKERS = {
    "Steam": "o",
    "SCW": "s",
    "SCW+N2": "^",
    "SCW+CO2": "D",
    "SCW+N2+CO2": "X",
}
DISPLAY_NAMES = {
    "Steam": "蒸汽",
    "SCW": "超临界水",
    "SCW+N2": "超临界水 + N$_2$",
    "SCW+CO2": "超临界水 + CO$_2$",
    "SCW+N2+CO2": "超临界水 + N$_2$ + CO$_2$",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "mathtext.fontset": "dejavuserif",
            "axes.unicode_minus": False,
            "axes.linewidth": 0.9,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        0.015,
        0.985,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
        zorder=20,
    )


def export(fig: mpl.figure.Figure, stem: str) -> None:
    png_path = FIGURE_DIR / f"{stem}.png"
    fig.savefig(png_path, dpi=450)
    fig.savefig(FIGURE_DIR / f"{stem}.pdf")
    plt.close(fig)
    with Image.open(png_path) as image:
        image.convert("RGB").save(png_path, dpi=(450, 450))


def add_box(
    ax: mpl.axes.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    label: str,
    facecolor: str,
    edgecolor: str = "#233746",
    fontsize: float = 8.0,
) -> FancyBboxPatch:
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.25,rounding_size=0.7",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=1.1,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        label,
        ha="center",
        va="center",
        fontsize=fontsize,
    )
    return box


def add_arrow(
    ax: mpl.axes.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = "#233746",
    linestyle: str = "-",
    linewidth: float = 1.3,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=linewidth,
            linestyle=linestyle,
            color=color,
            shrinkA=0,
            shrinkB=0,
        )
    )


def plot_experimental_system() -> None:
    fig, ax = plt.subplots(figsize=(180 / 25.4, 86 / 25.4), layout="constrained")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 42)
    ax.axis("off")

    blue = "#DCEFF7"
    orange = "#FCE6D8"
    green = "#DDF1E9"
    gray = "#ECEFF1"
    beige = "#F4D3A7"

    add_box(ax, (2, 27), 11, 6, "供水与\n高压泵", blue)
    add_box(ax, (17, 27), 14, 6, "超临界水\n发生器", blue)
    add_arrow(ax, (13, 30), (17, 30), COLORS["SCW"])

    add_box(ax, (2, 14), 11, 6, "CO$_2$/N$_2$\n气源", orange)
    add_box(ax, (17, 14), 14, 6, "气体计量与\n增压单元", orange)
    add_arrow(ax, (13, 17), (17, 17), COLORS["SCW+CO2"])

    add_box(ax, (35, 20), 10, 7, "混合/切换\n阀组", gray)
    add_arrow(ax, (31, 30), (35, 24.5), COLORS["SCW"])
    add_arrow(ax, (31, 17), (35, 22.5), COLORS["SCW+CO2"])

    model_x, model_y, model_w, model_h = 49, 18, 25, 10
    ax.add_patch(
        FancyBboxPatch(
            (model_x - 1.0, model_y - 1.2),
            model_w + 2.0,
            model_h + 2.4,
            boxstyle="round,pad=0.2,rounding_size=0.8",
            facecolor="#FFF4E6",
            edgecolor=COLORS["SCW"],
            linewidth=1.6,
            linestyle="--",
        )
    )
    ax.add_patch(
        Rectangle(
            (model_x, model_y),
            model_w,
            model_h,
            facecolor=beige,
            edgecolor="#233746",
            linewidth=1.3,
        )
    )
    ax.text(
        model_x + model_w / 2,
        model_y + model_h / 2,
        "一维填砂管模型\n0.48 m × 0.039 m",
        ha="center",
        va="center",
        fontsize=8.5,
    )
    ax.text(
        model_x + model_w / 2,
        model_y - 2.6,
        "分段加热与保温套",
        color=COLORS["SCW"],
        ha="center",
        va="top",
        fontsize=7.4,
    )
    add_arrow(ax, (45, 23.5), (49, 23.5), "#233746")

    for index, probe_x in enumerate(np.linspace(model_x + 2.0, model_x + model_w - 2.0, 7), 1):
        ax.plot([probe_x, probe_x], [model_y + model_h, model_y + model_h + 4.2], color=COLORS["SCW+N2"], linewidth=0.9)
        ax.add_patch(Circle((probe_x, model_y + model_h), 0.38, facecolor=COLORS["SCW+N2"], edgecolor="white", linewidth=0.4))
        ax.text(probe_x, model_y + model_h + 4.8, f"T{index}", ha="center", fontsize=6.7)

    for sensor_x, sensor_label in [(model_x, "$P_{in}$"), (model_x + model_w, "$P_{out}$")]:
        ax.add_patch(Circle((sensor_x, model_y + 1.4), 0.75, facecolor="white", edgecolor=COLORS["Steam"], linewidth=1.2))
        ax.text(sensor_x, model_y + 1.4, sensor_label, ha="center", va="center", fontsize=6.5)

    add_box(ax, (51, 35), 21, 5, "温度/压力采集与控制", green, fontsize=8.1)
    for probe_x in np.linspace(model_x + 2.0, model_x + model_w - 2.0, 7):
        ax.plot([probe_x, probe_x], [model_y + model_h + 4.2, 35], color=COLORS["SCW+N2"], linewidth=0.55, linestyle=":")

    add_box(ax, (78, 20), 9, 7, "冷凝器与\n恒温水浴", blue)
    add_arrow(ax, (74, 23.5), (78, 23.5))
    add_box(ax, (90, 20), 8, 7, "背压阀\nBPR", gray)
    add_arrow(ax, (87, 23.5), (90, 23.5))

    add_box(ax, (78, 7), 10, 7, "气液\n分离器", green)
    add_arrow(ax, (94, 20), (88, 14))
    add_box(ax, (91, 9.5), 7, 4.5, "气体计量", orange, fontsize=7.2)
    add_arrow(ax, (88, 11.8), (91, 11.8), COLORS["SCW+CO2"])
    add_box(ax, (78, 0.5), 10, 4.5, "液体收集", blue, fontsize=7.2)
    add_arrow(ax, (83, 7), (83, 5), COLORS["Steam"])

    ax.text(40, 29.5, "伴热管线", ha="center", fontsize=7.2, color=COLORS["SCW"])
    ax.text(1.0, 39.5, "文献实验流程简化重绘", fontsize=11, fontweight="bold", ha="left")
    ax.text(
        1.0,
        -0.2,
        "依据 Sun et al. (2024) Fig. 1 简化重绘，非按比例；DOI: 10.3390/app14093588（CC BY 4.0）",
        fontsize=6.8,
        color="#455A64",
        ha="left",
        va="bottom",
    )
    export(fig, "04_literature_experimental_system")


def scatter_family_points(ax: mpl.axes.Axes, data: pd.DataFrame, y_column: str) -> None:
    for family in FAMILIES:
        group = data[data["fluid_family"] == family]
        ax.scatter(
            group["experiment"],
            group[y_column],
            s=35,
            color=COLORS[family],
            marker=MARKERS[family],
            edgecolors="white" if MARKERS[family] != "X" else COLORS[family],
            linewidths=0.55,
            label=DISPLAY_NAMES[family],
            zorder=3,
        )


def highlight_experiment_12(ax: mpl.axes.Axes, y: float) -> None:
    ax.scatter(
        [12],
        [y],
        s=115,
        facecolors="none",
        edgecolors="#D55E00",
        linewidths=1.6,
        zorder=5,
    )


def plot_experiment_scheme(data: pd.DataFrame) -> None:
    fig, axes = plt.subplots(
        1, 3, figsize=(180 / 25.4, 88 / 25.4), sharex=True, layout="constrained"
    )
    ax_t, ax_p, ax_gas = axes

    scatter_family_points(ax_t, data, "temperature_C")
    ax_t.axhline(
        CRITICAL_TEMPERATURE_C,
        color="#6D6D6D",
        linestyle="--",
        linewidth=1.0,
        label="水临界温度",
    )
    highlight_experiment_12(ax_t, 400.0)
    ax_t.set(ylabel="注入温度 (°C)", ylim=(300, 430))

    scatter_family_points(ax_p, data, "pressure_MPa")
    ax_p.axhline(
        CRITICAL_PRESSURE_MPA,
        color="#6D6D6D",
        linestyle="--",
        linewidth=1.0,
    )
    highlight_experiment_12(ax_p, 24.0)
    ax_p.set(ylabel="注入压力 (MPa)", ylim=(8, 27))

    experiments = data["experiment"].to_numpy()
    ax_gas.bar(
        experiments - 0.18,
        data["n2_rate_mL_min"],
        width=0.36,
        color=COLORS["SCW+N2"],
        edgecolor="#233746",
        linewidth=0.4,
        label="N$_2$",
    )
    ax_gas.bar(
        experiments + 0.18,
        data["co2_rate_mL_min"],
        width=0.36,
        color=COLORS["SCW+CO2"],
        edgecolor="#233746",
        linewidth=0.4,
        hatch="//",
        label="CO$_2$",
    )
    ax_gas.axvspan(11.55, 12.45, facecolor="#D55E00", alpha=0.08, edgecolor="#D55E00")
    ax_gas.set(ylabel="气体注入流量 (mL·min$^{-1}$)", ylim=(0, 2.35))
    ax_gas.legend(loc="upper left", fontsize=7.2)
    ax_gas.text(
        0.98,
        0.97,
        "各组水流量均为\n10 mL·min$^{-1}$",
        transform=ax_gas.transAxes,
        ha="right",
        va="top",
        fontsize=7.2,
    )

    for label, ax in zip(["(a)", "(b)", "(c)"], axes):
        panel_label(ax, label)
        ax.set_xlabel("文献实验编号")
        ax.set_xlim(0.3, 16.7)
        ax.set_xticks([1, 4, 8, 12, 16])
        for spine in ax.spines.values():
            spine.set_visible(True)

    handles, labels = ax_t.get_legend_handles_labels()
    family_handles = handles[: len(FAMILIES)]
    family_labels = labels[: len(FAMILIES)]
    fig.legend(
        family_handles,
        family_labels,
        loc="outside lower center",
        ncol=5,
        columnspacing=1.0,
        handletextpad=0.35,
        fontsize=7.2,
    )
    fig.text(
        0.5,
        0.005,
        "数据来源：Sun et al. (2024), Table 2；橙色外圈/浅色带标示本文采用的 Exp. 12 工况。",
        ha="center",
        fontsize=6.8,
        color="#455A64",
    )
    export(fig, "05_literature_experiment_scheme")


def plot_literature_results(scheme: pd.DataFrame, results: pd.DataFrame) -> None:
    data = results.merge(scheme[["experiment", "fluid_family"]], on="experiment", how="left")
    fig, axes = plt.subplots(
        1, 2, figsize=(180 / 25.4, 82 / 25.4), layout="constrained"
    )
    ax_eff, ax_heat = axes

    for family in FAMILIES:
        group = data[data["fluid_family"] == family]
        ax_eff.scatter(
            group["experiment"],
            group["displacement_efficiency_percent"],
            s=39,
            color=COLORS[family],
            marker=MARKERS[family],
            edgecolors="white" if MARKERS[family] != "X" else COLORS[family],
            linewidths=0.55,
            label=DISPLAY_NAMES[family],
            zorder=3,
        )
        ax_heat.scatter(
            group["heat_utilization_rate"],
            group["displacement_efficiency_percent"],
            s=39,
            color=COLORS[family],
            marker=MARKERS[family],
            edgecolors="white" if MARKERS[family] != "X" else COLORS[family],
            linewidths=0.55,
            zorder=3,
        )

    exp12 = data[data["experiment"] == 12].iloc[0]
    ax_eff.scatter(
        [12],
        [exp12["displacement_efficiency_percent"]],
        s=125,
        facecolors="none",
        edgecolors="#D55E00",
        linewidths=1.6,
        zorder=5,
    )
    ax_eff.scatter(
        [12],
        [MODEL_RECOVERY_PERCENT],
        s=70,
        marker="*",
        facecolor="white",
        edgecolor="#000000",
        linewidth=1.0,
        zorder=6,
        label="当前 Exp.12 派生模型（4 PV）",
    )
    ax_eff.annotate(
        "+4.09 个百分点",
        xy=(12, MODEL_RECOVERY_PERCENT),
        xytext=(9.2, 90.5),
        fontsize=7.0,
        arrowprops={"arrowstyle": "->", "color": "#000000", "linewidth": 0.8},
    )
    ax_eff.set(
        xlabel="文献实验编号",
        ylabel="终点驱替效率 (%)",
        xlim=(0.3, 16.7),
        ylim=(55, 95),
    )
    ax_eff.set_xticks([1, 4, 8, 12, 16])

    for _, row in data.iterrows():
        dx = 0.00012 if row["experiment"] % 2 else -0.00012
        ha = "left" if dx > 0 else "right"
        ax_heat.text(
            row["heat_utilization_rate"] + dx,
            row["displacement_efficiency_percent"] + 0.15,
            str(int(row["experiment"])),
            fontsize=6.2,
            ha=ha,
            va="bottom",
        )
    ax_heat.set(
        xlabel="热利用率（文献定义）",
        ylabel="终点驱替效率 (%)",
        xlim=(0.0162, 0.0283),
        ylim=(55, 95),
    )

    for label, ax in zip(["(a)", "(b)"], axes):
        panel_label(ax, label)
        for spine in ax.spines.values():
            spine.set_visible(True)

    handles, labels = ax_eff.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="outside lower center",
        ncol=3,
        columnspacing=1.0,
        handletextpad=0.35,
        fontsize=7.0,
    )
    fig.text(
        0.5,
        0.005,
        "文献数据：Sun et al. (2024), Table 5；星号为当前等温 nC16 派生模型，不能视为文献实验历史拟合。",
        ha="center",
        fontsize=6.8,
        color="#455A64",
    )
    export(fig, "06_literature_experimental_results")


def write_manifest() -> None:
    manifest = {
        "citation": {
            "authors": "Sun et al.",
            "title": "Analysis of Adaptability and Application Potential of Supercritical Multi-Source Multi-Component Thermal Fluid Technology for Offshore Heavy Oil in China",
            "journal": "Applied Sciences",
            "year": 2024,
            "volume": 14,
            "article": 3588,
            "doi": DOI,
            "license": "CC BY 4.0",
            "access_date": ACCESS_DATE,
        },
        "source_mapping": {
            "04_literature_experimental_system": "simplified redraw of Figure 1 and associated apparatus description, PDF page 5",
            "05_literature_experiment_scheme": "exact values transcribed from Table 2, PDF page 4",
            "06_literature_experimental_results": "exact values transcribed from Table 5, PDF page 9; current model endpoint from summary_metrics.json",
        },
        "transformations": [
            "no digitization of plotted literature curves",
            "no smoothing or interpolation",
            "experiment-family colors joined from Table 2 to Table 5 by experiment number",
            "current model recovery displayed only as a separately labeled comparison marker",
        ],
        "source_tables": [
            "data/sun2024_table2_experiment_scheme.csv",
            "data/sun2024_table5_experimental_results.csv",
        ],
        "generated_by": str(Path(__file__).relative_to(CASE_DIR)),
        "python": platform.python_version(),
        "matplotlib": mpl.__version__,
        "pandas": pd.__version__,
        "exports": ["PNG 450 dpi RGB opaque", "PDF vector"],
        "palette": [COLORS[family] for family in FAMILIES],
        "background": "#FFFFFF",
    }
    (REPORT_DIR / "literature_figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    configure_style()
    scheme = pd.read_csv(DATA_DIR / "sun2024_table2_experiment_scheme.csv")
    results = pd.read_csv(DATA_DIR / "sun2024_table5_experimental_results.csv")
    plot_experimental_system()
    plot_experiment_scheme(scheme)
    plot_literature_results(scheme, results)
    write_manifest()
    print("Generated three literature-based figures and a provenance manifest.")


if __name__ == "__main__":
    main()
