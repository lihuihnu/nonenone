#!/usr/bin/env python3
"""Compare fully compositional and traditional independent-water results.

The script reads simulator CSV outputs without smoothing or interpolation.
Production rates are multiplied by -1 so production is shown as a positive
magnitude; injector BHP at t=0 is omitted because it is only the initial guess.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_FULL = (
    ROOT.parent / "five_component_eos_tuned_compare" / "results"
    / "translated_density_30day_20260826"
)
DEFAULT_TRADITIONAL = ROOT / "results" / "translated_density_30day_20260826"
DEFAULT_FIGURES = ROOT / "figures" / "translated_density_30day_20260826"

EOS_ORDER = ("pr", "sw", "cpa")
EOS_LABEL = {"pr": "PR", "sw": "SW", "cpa": "CPA"}
COLORS = {"pr": "#0072B2", "sw": "#D55E00", "cpa": "#009E73"}
MARKERS = {"pr": "o", "sw": "s", "cpa": "^"}
SERIES_ORDER = (
    ("full", "pr"), ("full", "sw"), ("full", "cpa"),
    ("traditional", "pr"),
)
SERIES_LABEL = {
    ("full", "pr"): "Our PR",
    ("full", "sw"): "Our SW",
    ("full", "cpa"): "Our CPA",
    ("traditional", "pr"): "Traditional",
}
SERIES_COLOR = {
    ("full", "pr"): COLORS["pr"],
    ("full", "sw"): COLORS["sw"],
    ("full", "cpa"): COLORS["cpa"],
    ("traditional", "pr"): "#5F5F5F",
}
SERIES_MARKER = {
    ("full", "pr"): "o",
    ("full", "sw"): "s",
    ("full", "cpa"): "^",
    ("traditional", "pr"): "D",
}
SERIES_LINE = {
    ("full", "pr"): "-",
    ("full", "sw"): "--",
    ("full", "cpa"): "-.",
    ("traditional", "pr"): ":",
}

WIDTH_MM = 150.0
HEIGHT_MM = 92.0
DPI = 600
FINAL_DAY = 30.0
EXPECTED_STEPS = 120
TIME_TICKS = np.arange(0.0, FINAL_DAY + 0.1, 5.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-root", type=Path, default=DEFAULT_FULL)
    parser.add_argument("--traditional-root", type=Path, default=DEFAULT_TRADITIONAL)
    parser.add_argument("--figures-root", type=Path, default=DEFAULT_FIGURES)
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
    if not rows:
        raise ValueError(f"empty source file: {path}")
    return rows


def values(rows: list[dict[str, str]], column: str) -> np.ndarray:
    result = np.asarray([float(row[column]) for row in rows], dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"non-finite {column} values")
    return result


def load_model(
    root: Path, eos_order: tuple[str, ...],
) -> dict[str, dict[str, list[dict[str, str]]]]:
    model: dict[str, dict[str, list[dict[str, str]]]] = {}
    reference_time: np.ndarray | None = None
    for eos in eos_order:
        reservoir = read_rows(root / eos / "reservoir_diagnostics.csv")
        wells = read_rows(root / eos / "well_history.csv")
        solves = read_rows(root / eos / "nonlinear_solve_history.csv")
        summary = read_rows(root / eos / "simulation_summary.csv")
        solution_initial = read_rows(root / eos / "solution_step_0.csv")
        solution_final = read_rows(root / eos / "solution_final.csv")
        producer = [row for row in wells if row["name"] == "PROD"]
        injector = [row for row in wells if row["name"] == "CO2_INJ"]
        time = values(reservoir, "time")
        expected_records = EXPECTED_STEPS + 1
        if (len(reservoir) != expected_records or len(producer) != expected_records
                or len(injector) != expected_records):
            raise ValueError(
                f"{root}/{eos}: expected {expected_records} reservoir and well records"
            )
        summary_row = summary[0]
        accepted_solves = sum(row["converged"] == "1" for row in solves)
        rejected_solves = len(solves) - accepted_solves
        if int(summary_row["fixed_output_intervals"]) != EXPECTED_STEPS:
            raise ValueError(
                f"{root}/{eos}: expected {EXPECTED_STEPS} fixed output intervals"
            )
        if len(solves) != int(summary_row["nonlinear_solves"]):
            raise ValueError(f"{root}/{eos}: nonlinear solve history is incomplete")
        if accepted_solves != int(summary_row["accepted_internal_steps"]):
            raise ValueError(f"{root}/{eos}: accepted solve count is inconsistent")
        if rejected_solves != int(summary_row["rejected_internal_steps"]):
            raise ValueError(f"{root}/{eos}: rejected solve count is inconsistent")
        if len(solution_initial) != 2000 or len(solution_final) != 2000:
            raise ValueError(f"{root}/{eos}: expected 2000 solution records")
        if not np.isclose(float(summary_row["final_time_day"]), FINAL_DAY):
            raise ValueError(f"{root}/{eos}: final time is not {FINAL_DAY:g} day")
        if reference_time is None:
            reference_time = time
        elif not np.array_equal(time, reference_time):
            raise ValueError("output times differ between data sets")
        model[eos] = {
            "reservoir": reservoir,
            "producer": producer,
            "injector": injector,
            "solves": solves,
            "summary": summary,
            "solution_initial": solution_initial,
            "solution_final": solution_final,
        }
    return model


def configure_fonts(language: str) -> None:
    family = "DejaVu Sans"
    if language == "zh":
        candidates = (
            Path("/mnt/c/Windows/Fonts/msyh.ttc"),
            Path("C:/Windows/Fonts/msyh.ttc"),
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
        "font.size": 9.8,
        "axes.labelsize": 10.8,
        "xtick.labelsize": 9.3,
        "ytick.labelsize": 9.3,
        "legend.fontsize": 8.6,
        "axes.linewidth": 1.0,
        "lines.linewidth": 1.9,
        "lines.markersize": 4.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.unicode_minus": False,
    })


def new_axes() -> tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(
        figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
        layout="constrained",
    )
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#202020")
        spine.set_linewidth(1.0)
    ax.tick_params(direction="in", length=3.8, width=0.9, top=False, right=False)
    return fig, ax


def text_labels(language: str) -> dict[str, str]:
    if language == "zh":
        return {
            "time": "时间（天）",
            "pressure": "储层平均压力（bar）",
            "production": "生产井地面产量（m³/s）",
            "bhp": "注入井井底压力（bar）",
            "water_density": "水相平均密度（kg/m³）",
            "saturation": "平均饱和度",
            "runtime": "完整模拟耗时（s）",
            "water_co2": r"水相 CO$_2$ 平均摩尔分数（×10$^{-3}$）",
            "water_co2_change": r"水相 CO$_2$ 摩尔分数增量（×10$^6$）",
            "water_oil_partition": r"CO$_2$ 水/油相分配比",
            "gas_oil_partition": r"CO$_2$ 气/油相平衡比",
            "water_viscosity": "水相平均黏度（mPa·s）",
            "full_tick": "全组分",
            "traditional_tick": "独立水",
        }
    return {
        "time": "Time (day)",
        "pressure": "Mean reservoir pressure (bar)",
        "production": "Producer surface rate (m³/s)",
        "bhp": "Injector BHP (bar)",
        "water_density": "Mean water-phase density (kg/m³)",
        "saturation": "Mean saturation",
        "runtime": "Full simulation time (s)",
        "water_co2": r"Mean aqueous CO$_2$ mole fraction (×10$^{-3}$)",
        "water_co2_change": r"Aqueous CO$_2$ mole-fraction increase (×10$^6$)",
        "water_oil_partition": r"CO$_2$ water/oil partition ratio",
        "gas_oil_partition": r"CO$_2$ gas/oil equilibrium ratio",
        "water_viscosity": "Mean water-phase viscosity (mPa·s)",
        "full_tick": "Full",
        "traditional_tick": "Independent H$_2$O",
    }


def padded_limits(series: list[np.ndarray], fraction: float = 0.08) -> tuple[float, float]:
    lower = min(float(np.min(item)) for item in series)
    upper = max(float(np.max(item)) for item in series)
    span = max(upper - lower, max(abs(lower), abs(upper), 1.0) * 0.01)
    return lower - fraction * span, upper + fraction * span


def model_legend(ax: plt.Axes, location: str, columns: int = 2) -> None:
    handles = [
        Line2D([], [], color=SERIES_COLOR[key], marker=SERIES_MARKER[key],
               linestyle=SERIES_LINE[key], markerfacecolor=SERIES_COLOR[key],
               label=SERIES_LABEL[key])
        for key in SERIES_ORDER
    ]
    ax.legend(handles=handles, loc=location, ncol=columns, frameon=True,
              framealpha=0.97, facecolor="white", edgecolor="#A0A0A0")


def plot_time_series(ax: plt.Axes, data, source: str, column: str,
                     transform=None, skip_first: bool = False) -> list[np.ndarray]:
    plotted: list[np.ndarray] = []
    for model, eos in SERIES_ORDER:
        rows = data[model][eos][source]
        x = values(rows, "time")
        y = values(rows, column)
        if transform is not None:
            y = transform(y)
        if skip_first:
            x, y = x[1:], y[1:]
        plotted.append(y)
        key = (model, eos)
        ax.plot(x, y, color=SERIES_COLOR[key], linestyle=SERIES_LINE[key],
                marker=SERIES_MARKER[key], markevery=max(1, len(x) // 6),
                markerfacecolor=SERIES_COLOR[key],
                markeredgecolor=SERIES_COLOR[key], markeredgewidth=1.0, zorder=3)
    return plotted


def co2_composition(data, model: str, eos: str, phase: str,
                    step: str = "final") -> np.ndarray:
    rows = data[model][eos][f"solution_{step}"]
    if model == "traditional":
        if phase == "water":
            return np.zeros(len(rows), dtype=float)
        column = {"oil": "liquid_x_0_CO2", "gas": "vapor_y_0_CO2"}[phase]
    else:
        column = {
            "water": "water_x_1_CO2",
            "oil": "oil_x_1_CO2",
            "gas": "gas_y_1_CO2",
        }[phase]
    return values(rows, column)


def series_means(data, quantity) -> np.ndarray:
    return np.asarray([quantity(model, eos) for model, eos in SERIES_ORDER], dtype=float)


def model_handles() -> list[Line2D]:
    return [
        Line2D([], [], color=SERIES_COLOR[key], marker=SERIES_MARKER[key],
               linestyle="none", markerfacecolor=SERIES_COLOR[key],
               label=SERIES_LABEL[key])
        for key in SERIES_ORDER
    ]


def export(fig: plt.Figure, output_root: Path, language: str,
           stem: str, force: bool) -> list[Path]:
    target = output_root / language
    target.mkdir(parents=True, exist_ok=True)
    paths = [target / f"{stem}.png", target / f"{stem}.pdf"]
    if not force and any(path.exists() for path in paths):
        raise FileExistsError(f"figure exists; use --force: {paths[0]}")
    fig.savefig(paths[0], dpi=DPI, metadata={"Software": "Our flow-comparison plotter"})
    fig.savefig(paths[1], dpi=DPI, metadata={
        "Creator": "Our flow-comparison plotter",
        "Subject": "Fully compositional versus independent-water comparison",
    })
    plt.close(fig)
    return paths


def pressure_figure(data, output_root: Path, language: str, force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    plotted = plot_time_series(ax, data, "reservoir", "pressure_avg")
    ax.set(xlim=(0, FINAL_DAY), xticks=TIME_TICKS,
           ylim=padded_limits(plotted), xlabel=label["time"], ylabel=label["pressure"])
    model_legend(ax, "lower left", 2)
    return export(fig, output_root, language, "01_mean_pressure", force)


def production_figure(data, output_root: Path, language: str, force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    plotted = plot_time_series(
        ax, data, "producer", "q_total_surface", transform=lambda item: -item,
    )
    lower, upper = padded_limits(plotted)
    ax.set(xlim=(0, FINAL_DAY), xticks=TIME_TICKS,
           ylim=(lower, upper + 0.30 * (upper - lower)),
           xlabel=label["time"], ylabel=label["production"])
    model_legend(ax, "upper right", 2)
    return export(fig, output_root, language, "02_producer_rate", force)


def bhp_figure(data, output_root: Path, language: str, force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    plotted = plot_time_series(ax, data, "injector", "bhp", skip_first=True)
    lower, upper = padded_limits(plotted)
    span = upper - lower
    # Reserve in-frame space above every curve so the legend never masks data.
    ax.set(xlim=(0, FINAL_DAY), xticks=TIME_TICKS,
           ylim=(lower, upper + 0.35 * span),
           xlabel=label["time"], ylabel=label["bhp"])
    model_legend(ax, "upper center", 2)
    return export(fig, output_root, language, "03_injector_bhp", force)


def water_density_figure(data, output_root: Path, language: str, force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    records = list(SERIES_ORDER)
    initial = np.asarray([
        float(data[model][eos]["reservoir"][0]["rho_w_avg"])
        for model, eos in records
    ])
    final = np.asarray([
        float(data[model][eos]["reservoir"][-1]["rho_w_avg"])
        for model, eos in records
    ])
    y = np.arange(len(records))[::-1]
    for index, ((model, eos), y_value) in enumerate(zip(records, y)):
        key = (model, eos)
        ax.plot([initial[index], final[index]], [y_value, y_value],
                color="#A6A6A6", linewidth=1.4, zorder=1)
        ax.scatter(initial[index], y_value, s=34, facecolor="white",
                   edgecolor="#5F5F5F", linewidth=1.0, marker="o", zorder=2)
        ax.scatter(final[index], y_value, s=42,
                   facecolor=SERIES_COLOR[key], edgecolor=SERIES_COLOR[key],
                   linewidth=1.25, marker=SERIES_MARKER[key], zorder=3)
        ax.annotate(f"{final[index]:.2f}", (final[index], y_value),
                    xytext=(6, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=8.3)
    names = [SERIES_LABEL[key] for key in records]
    lower = min(float(initial.min()), float(final.min())) - 0.45
    upper = max(float(initial.max()), float(final.max())) + 0.95
    ax.set(xlim=(lower, upper), ylim=(-0.65, len(records) - 0.35),
           yticks=y, yticklabels=names, xlabel=label["water_density"])
    ax.legend(handles=(
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
               markeredgecolor="#5F5F5F", label="Initial"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="#303030",
               markeredgecolor="#303030", label="Day 30"),
    ), loc="lower right", ncol=2, frameon=True, framealpha=0.97,
       facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "04_water_density", force)


def saturation_figure(data, output_root: Path, language: str, force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    phase_specs = (
        ("Oil", "s_o_avg", "#CC6677", ""),
        ("Gas", "s_g_avg", "#4477AA", "//"),
        ("Water", "s_w_avg", "#228833", ".."),
    )
    x = np.arange(4)
    bottom = np.zeros(4)
    records = list(SERIES_ORDER)
    for phase, column, color, hatch in phase_specs:
        height = np.asarray([
            float(data[model][eos]["reservoir"][-1][column])
            for model, eos in records
        ])
        bars = ax.bar(x, height, bottom=bottom, width=0.68, color=color,
                      edgecolor="white", linewidth=0.8, hatch=hatch, label=phase)
        for bar, value, base in zip(bars, height, bottom):
            if value >= 0.06:
                ax.text(bar.get_x() + bar.get_width() / 2, base + value / 2,
                        f"{value:.3f}", ha="center", va="center", fontsize=7.8,
                        color="white" if value >= 0.15 else "#101010")
        bottom += height
    ticks = [
        SERIES_LABEL[(model, eos)] for model, eos in records
    ]
    # The empty right-hand part of the same axes is reserved for the legend.
    ax.set(xlim=(-0.6, 5.4), ylim=(0, 1), xticks=x, xticklabels=ticks,
           ylabel=label["saturation"])
    ax.legend(loc="center right", ncol=1, frameon=True, framealpha=0.97,
              facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "05_final_saturation", force)


def runtime_figure(data, output_root: Path, language: str, force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    records = list(SERIES_ORDER)
    times = np.asarray([
        float(data[model][eos]["summary"][0]["simulation_loop_wall"])
        for model, eos in records
    ])
    x = np.arange(4)
    bars = []
    for index, ((model, eos), wall) in enumerate(zip(records, times)):
        key = (model, eos)
        bars.append(ax.bar(
            index, wall, width=0.68, color=SERIES_COLOR[key],
            edgecolor=SERIES_COLOR[key], linewidth=1.4, zorder=3,
        )[0])
    for bar, wall in zip(bars, times):
        ax.text(bar.get_x() + bar.get_width() / 2, wall + 1.2, f"{wall:.1f}",
                ha="center", va="bottom", fontsize=8.4)
    ticks = [SERIES_LABEL[key] for key in records]
    ax.set(xlim=(-0.6, 3.6), ylim=(0, max(times) * 1.45), xticks=x,
           xticklabels=ticks, ylabel=label["runtime"])
    ax.legend(handles=model_handles(), loc="upper left", ncol=2, frameon=True,
              framealpha=0.97, facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "06_runtime", force)


def water_co2_figure(data, output_root: Path, language: str,
                     force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    means = 1.0e3 * series_means(
        data, lambda model, eos: np.mean(co2_composition(data, model, eos, "water")),
    )
    y = np.arange(len(SERIES_ORDER))[::-1]
    for key, position, value in zip(SERIES_ORDER, y, means):
        ax.barh(position, value, height=0.54, color=SERIES_COLOR[key],
                edgecolor=SERIES_COLOR[key], zorder=2)
        ax.text(value + 0.055, position, f"{value:.3f}", ha="left", va="center",
                fontsize=8.5)
    ax.set(xlim=(0, max(means) * 1.55), ylim=(-0.7, 4.25), yticks=[],
           xlabel=label["water_co2"])
    ax.legend(handles=model_handles(), loc="upper right", ncol=2, frameon=True,
              framealpha=0.97, facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "07_aqueous_co2_mole_fraction", force)


def water_co2_change_figure(data, output_root: Path, language: str,
                            force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    changes = 1.0e6 * series_means(
        data,
        lambda model, eos: (
            np.mean(co2_composition(data, model, eos, "water", "final"))
            - np.mean(co2_composition(data, model, eos, "water", "initial"))
        ),
    )
    x = np.arange(len(SERIES_ORDER))
    for key, position, value in zip(SERIES_ORDER, x, changes):
        ax.vlines(position, 0, value, color=SERIES_COLOR[key], linewidth=2.5, zorder=2)
        ax.scatter(position, value, s=52, color=SERIES_COLOR[key],
                   marker=SERIES_MARKER[key], edgecolor="white", linewidth=0.7, zorder=3)
        ax.text(position, value + 0.16, f"{value:.2f}", ha="center", va="bottom",
                fontsize=8.5)
    ax.axhline(0, color="#606060", linewidth=0.9, zorder=1)
    ax.set(xlim=(-0.55, 4.75), ylim=(-0.25, max(changes) * 1.60), xticks=[],
           ylabel=label["water_co2_change"])
    ax.legend(handles=model_handles(), loc="upper right", ncol=2, frameon=True,
              framealpha=0.97, facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "08_aqueous_co2_change", force)


def water_oil_partition_figure(data, output_root: Path, language: str,
                               force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    ratios = series_means(
        data,
        lambda model, eos: np.mean(
            co2_composition(data, model, eos, "water")
            / co2_composition(data, model, eos, "oil")
        ),
    )
    x = np.arange(len(SERIES_ORDER))
    for key, position, value in zip(SERIES_ORDER, x, ratios):
        ax.plot([position, position], [0, value], color=SERIES_COLOR[key],
                linewidth=2.2, zorder=2)
        ax.scatter(position, value, s=58, color=SERIES_COLOR[key],
                   marker=SERIES_MARKER[key], edgecolor="white", linewidth=0.7, zorder=3)
        ax.text(position, value + 0.0008, f"{value:.4f}", ha="center", va="bottom",
                fontsize=8.4)
    ax.set(xlim=(-0.55, 4.75), ylim=(0, max(ratios) * 1.62), xticks=[],
           ylabel=label["water_oil_partition"])
    ax.legend(handles=model_handles(), loc="upper right", ncol=2, frameon=True,
              framealpha=0.97, facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "09_co2_water_oil_partition", force)


def gas_oil_partition_figure(data, output_root: Path, language: str,
                             force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    ratios = series_means(
        data,
        lambda model, eos: np.mean(
            co2_composition(data, model, eos, "gas")
            / co2_composition(data, model, eos, "oil")
        ),
    )
    x = np.arange(len(SERIES_ORDER))
    ax.plot(x, ratios, color="#B5B5B5", linewidth=1.2, zorder=1)
    for key, position, value in zip(SERIES_ORDER, x, ratios):
        ax.scatter(position, value, s=70, color=SERIES_COLOR[key],
                   marker=SERIES_MARKER[key], edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(position, value + 0.008, f"{value:.3f}", ha="center", va="bottom",
                fontsize=8.5)
    lower, upper = padded_limits([ratios], fraction=0.16)
    ax.set(xlim=(-0.55, 4.75), ylim=(lower, upper + 0.06), xticks=[],
           ylabel=label["gas_oil_partition"])
    ax.legend(handles=model_handles(), loc="upper right", ncol=2, frameon=True,
              framealpha=0.97, facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "10_co2_gas_oil_partition", force)


def water_viscosity_figure(data, output_root: Path, language: str,
                           force: bool) -> list[Path]:
    label = text_labels(language)
    fig, ax = new_axes()
    viscosities = series_means(
        data,
        lambda model, eos: 1.0e3 * float(
            data[model][eos]["reservoir"][-1]["mu_w_avg"]
        ),
    )
    x = np.arange(len(SERIES_ORDER))
    for key, position, value in zip(SERIES_ORDER, x, viscosities):
        ax.bar(position, value, width=0.62, color=SERIES_COLOR[key],
               edgecolor=SERIES_COLOR[key], zorder=2)
        ax.text(position, value + 0.025, f"{value:.3f}", ha="center", va="bottom",
                fontsize=8.5)
    ax.set(xlim=(-0.55, 4.75), ylim=(0, max(viscosities) * 1.52), xticks=[],
           ylabel=label["water_viscosity"])
    ax.legend(handles=model_handles(), loc="upper right", ncol=2, frameon=True,
              framealpha=0.97, facecolor="white", edgecolor="#A0A0A0")
    return export(fig, output_root, language, "11_water_viscosity", force)


def write_source_table(data, output_root: Path, force: bool) -> Path:
    path = output_root / "source_data.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        raise FileExistsError(f"source table exists; use --force: {path}")
    fields = (
        "model", "eos", "time_day", "pressure_avg_bar", "s_o_avg", "s_g_avg",
        "s_w_avg", "rho_w_avg_kg_m3", "producer_rate_magnitude_m3_s",
        "injector_bhp_bar", "cumulative_nonlinear_solve_wall_s",
    )
    with path.open("w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for model, eos in SERIES_ORDER:
            item = data[model][eos]
            solve_end = values(item["solves"], "end_time")
            solve_wall = values(item["solves"], "wall_time")
            output_time = values(item["reservoir"], "time")
            cumulative = np.asarray([
                np.sum(solve_wall[solve_end <= time_day + 1.0e-12])
                for time_day in output_time
            ])
            for index, row in enumerate(item["reservoir"]):
                writer.writerow({
                    "model": model,
                    "eos": EOS_LABEL[eos],
                    "time_day": row["time"],
                    "pressure_avg_bar": row["pressure_avg"],
                    "s_o_avg": row["s_o_avg"],
                    "s_g_avg": row["s_g_avg"],
                    "s_w_avg": row["s_w_avg"],
                    "rho_w_avg_kg_m3": row["rho_w_avg"],
                    "producer_rate_magnitude_m3_s":
                        -float(item["producer"][index]["q_total_surface"]),
                    "injector_bhp_bar": item["injector"][index]["bhp"] if index else "",
                    "cumulative_nonlinear_solve_wall_s": cumulative[index],
                })
    return path


def write_composition_table(data, output_root: Path, force: bool) -> Path:
    path = output_root / "composition_statistics.csv"
    if path.exists() and not force:
        raise FileExistsError(f"composition table exists; use --force: {path}")
    fields = (
        "model", "time_day", "aqueous_co2_mole_fraction_mean",
        "oil_co2_mole_fraction_mean", "gas_co2_mole_fraction_mean",
        "co2_water_oil_partition_mean", "co2_gas_oil_partition_mean",
        "water_viscosity_mPa_s",
    )
    with path.open("w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for model, eos in SERIES_ORDER:
            for step, time_day, reservoir_index in (
                ("initial", 0.0, 0), ("final", FINAL_DAY, -1),
            ):
                water = co2_composition(data, model, eos, "water", step)
                oil = co2_composition(data, model, eos, "oil", step)
                gas = co2_composition(data, model, eos, "gas", step)
                writer.writerow({
                    "model": SERIES_LABEL[(model, eos)],
                    "time_day": time_day,
                    "aqueous_co2_mole_fraction_mean": np.mean(water),
                    "oil_co2_mole_fraction_mean": np.mean(oil),
                    "gas_co2_mole_fraction_mean": np.mean(gas),
                    "co2_water_oil_partition_mean": np.mean(water / oil),
                    "co2_gas_oil_partition_mean": np.mean(gas / oil),
                    "water_viscosity_mPa_s": 1.0e3 * float(
                        data[model][eos]["reservoir"][reservoir_index]["mu_w_avg"]
                    ),
                })
    return path


def main() -> None:
    args = parse_args()
    full_root = resolve(args.full_root)
    traditional_root = resolve(args.traditional_root)
    output_root = resolve(args.figures_root)
    data = {
        "full": load_model(full_root, EOS_ORDER),
        "traditional": load_model(traditional_root, ("pr",)),
    }
    generated: list[Path] = [
        write_source_table(data, output_root, args.force),
        write_composition_table(data, output_root, args.force),
    ]
    builders = (
        pressure_figure, production_figure, bhp_figure,
        water_density_figure, saturation_figure, runtime_figure,
        water_co2_figure, water_co2_change_figure,
        water_oil_partition_figure, gas_oil_partition_figure,
        water_viscosity_figure,
    )
    for language in ("en", "zh"):
        configure_fonts(language)
        for builder in builders:
            generated.extend(builder(data, output_root, language, args.force))
    manifest = {
        "purpose": "general provisional scientific figures; no target-journal compliance claimed",
        "source_roots": {
            "fully_compositional": relative_label(full_root),
            "traditional_independent_water": relative_label(traditional_root),
        },
        "transformations": [
            "producer q_total_surface multiplied by -1 to display production magnitude",
            "injector BHP at t=0 omitted because it is an unconstrained initial guess",
            "all attempted nonlinear-solve wall times cumulatively grouped by physical end time",
            "final saturation bars use the unmodified day-30 arithmetic cell means",
            "runtime bars use the unmodified simulation_loop_wall values",
            "water-density dumbbells use only the unmodified initial and day-30 arithmetic cell means",
            "aqueous, oil, and gas CO2 compositions use unweighted arithmetic means over all 2000 cells",
            "aqueous CO2 change is the day-30 cell mean minus the initial cell mean",
            "CO2 partition ratios are unweighted arithmetic means of the per-cell composition ratios",
            "water viscosity uses the unmodified day-30 arithmetic cell mean reported by the simulator",
            "no filtering, smoothing, interpolation, normalization, or exclusions",
        ],
        "missing_data": "none; injector BHP t=0 is intentionally blank in source_data.csv",
        "figure_size_mm": [WIDTH_MM, HEIGHT_MM],
        "png_dpi": DPI,
        "palette": COLORS,
        "redundant_encoding": "EOS uses color and marker; phase model uses line style and marker fill",
        "legend_series": [SERIES_LABEL[key] for key in SERIES_ORDER],
        "generated": [relative_label(path) for path in generated],
    }
    manifest_path = output_root / "figure_manifest.json"
    if manifest_path.exists() and not args.force:
        raise FileExistsError(f"manifest exists; use --force: {manifest_path}")
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    figure_count = len(generated) - 2
    print(f"generated {figure_count} figure files and two source-data tables")


if __name__ == "__main__":
    main()
