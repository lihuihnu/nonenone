#!/usr/bin/env python3
"""Create bilingual, single-axis scientific figures for the tuned EOS case.

Raw source files are read directly from the selected results root.  The only
transformations are: production rates are displayed as positive production
magnitudes, injector BHP excludes the unconstrained t=0 initial guess, and
nonlinear-solve wall times are cumulatively summed in chronological order.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Callable

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "full"
FIGURES = ROOT / "figures"
BASELINE_CPA = ROOT / "results" / "full" / "cpa"
OPTIMIZED_CPA = ROOT / "results" / "optimization" / "full_optimized"
EOS_ORDER = ("pr", "sw", "cpa")
EOS_LABEL = {"pr": "PR", "sw": "SW", "cpa": "CPA"}
COLORS = {"pr": "#0072B2", "sw": "#D55E00", "cpa": "#009E73"}
LINESTYLES = {"pr": "-", "sw": "--", "cpa": "-."}
MARKERS = {"pr": "o", "sw": "s", "cpa": "^"}
WIDTH_MM = 160.0
HEIGHT_MM = 100.0
DPI = 600


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="overwrite existing figures")
    parser.add_argument(
        "--results-root", type=Path, default=RESULTS,
        help="directory containing pr, sw and cpa result directories",
    )
    parser.add_argument(
        "--figures-root", type=Path, default=FIGURES,
        help="output directory for figures and plotting tables",
    )
    parser.add_argument(
        "--baseline-cpa", type=Path, default=BASELINE_CPA,
        help="preserved original CPA result directory",
    )
    parser.add_argument(
        "--optimized-cpa", type=Path, default=OPTIMIZED_CPA,
        help="preserved B2-reduced CPA result directory",
    )
    return parser.parse_args()


def resolve_case_path(path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def source_label(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"missing result file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty result file: {path}")
    return rows


def numeric(rows: list[dict[str, str]], column: str) -> np.ndarray:
    values = np.asarray([float(row[column]) for row in rows], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"non-finite values in column {column}")
    return values


def load_data() -> dict[str, dict[str, list[dict[str, str]]]]:
    data: dict[str, dict[str, list[dict[str, str]]]] = {}
    reference_times: np.ndarray | None = None
    for eos in EOS_ORDER:
        reservoir = read_rows(RESULTS / eos / "reservoir_diagnostics.csv")
        wells = read_rows(RESULTS / eos / "well_history.csv")
        solves = read_rows(RESULTS / eos / "nonlinear_solve_history.csv")
        producer = [row for row in wells if row["name"] == "PROD"]
        injector = [row for row in wells if row["name"] == "CO2_INJ"]
        if len(reservoir) != 21 or len(producer) != 21 or len(injector) != 21:
            raise ValueError(f"{eos}: expected 21 reservoir and per-well records")
        if len(solves) != 20 or any(row["converged"] != "1" for row in solves):
            raise ValueError(f"{eos}: expected 20 converged nonlinear solves")
        times = numeric(reservoir, "time")
        if reference_times is None:
            reference_times = times
        elif not np.array_equal(times, reference_times):
            raise ValueError("EOS output times are not identical")
        data[eos] = {
            "reservoir": reservoir,
            "producer": producer,
            "injector": injector,
            "solves": solves,
        }
    return data


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
            raise FileNotFoundError("Chinese font not found (Microsoft YaHei or Noto Sans CJK)")
    mpl.rcParams.update({
        "font.family": family,
        "font.size": 10.5,
        "axes.labelsize": 11.5,
        "xtick.labelsize": 10.0,
        "ytick.labelsize": 10.0,
        "legend.fontsize": 9.5,
        "axes.linewidth": 1.0,
        "lines.linewidth": 2.0,
        "lines.markersize": 5.0,
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
        spine.set_linewidth(1.0)
        spine.set_color("#202020")
    ax.tick_params(direction="in", length=4.0, width=0.9, top=False, right=False)
    ax.set_xlim(0.0, 5.0)
    ax.set_xticks(np.arange(0.0, 5.1, 1.0))
    return fig, ax


def plot_eos_lines(
    ax: plt.Axes,
    data: dict[str, dict[str, list[dict[str, str]]]],
    source: str,
    x_column: str,
    y_column: str,
    transform: Callable[[np.ndarray], np.ndarray] | None = None,
    skip_first: bool = False,
) -> None:
    for eos in EOS_ORDER:
        rows = data[eos][source]
        x = numeric(rows, x_column)
        y = numeric(rows, y_column)
        if transform is not None:
            y = transform(y)
        if skip_first:
            x, y = x[1:], y[1:]
        ax.plot(
            x,
            y,
            color=COLORS[eos],
            linestyle=LINESTYLES[eos],
            marker=MARKERS[eos],
            markevery=max(1, len(x) // 5),
            markerfacecolor="white",
            markeredgewidth=1.2,
            label=EOS_LABEL[eos],
            zorder=3,
        )


def padded_limits(values: list[np.ndarray], fraction: float = 0.08) -> tuple[float, float]:
    lower = min(float(np.min(value)) for value in values)
    upper = max(float(np.max(value)) for value in values)
    span = max(upper - lower, max(abs(lower), abs(upper), 1.0) * 0.01)
    return lower - fraction * span, upper + fraction * span


def labels(language: str) -> dict[str, str]:
    if language == "zh":
        return {
            "time": "时间（天）",
            "pressure": "储层压力（bar）",
            "production": "生产井地面产量（m³/s）",
            "water_density": "含水相平均密度（kg/m³）",
            "injector_bhp": "注入井井底压力（bar）",
            "solver_time": "累计非线性求解时间（s）",
            "runtime": "完整模拟耗时（s，对数坐标）",
            "saturation": "平均饱和度",
        }
    return {
        "time": "Time (day)",
        "pressure": "Reservoir pressure (bar)",
        "production": "Producer surface rate (m³/s)",
        "water_density": "Mean water-rich-phase density (kg/m³)",
        "injector_bhp": "Injector BHP (bar)",
        "solver_time": "Cumulative nonlinear-solve time (s)",
        "runtime": "Full simulation time (s, log scale)",
        "saturation": "Mean saturation",
    }


def export(fig: plt.Figure, stem: str, language: str, force: bool) -> list[Path]:
    out_dir = FIGURES / language
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = [out_dir / f"{stem}.png", out_dir / f"{stem}.pdf"]
    existing = [path for path in outputs if path.exists()]
    if existing and not force:
        raise FileExistsError(f"refusing to overwrite: {existing[0]} (use --force)")
    metadata = {
        "Creator": "Our compositional-flow result plotting script",
        "Subject": "PR/SW/CPA five-component flow comparison",
    }
    fig.savefig(outputs[0], dpi=DPI, metadata={"Software": metadata["Creator"]})
    fig.savefig(outputs[1], dpi=DPI, metadata=metadata)
    plt.close(fig)
    return outputs


def pressure_figure(data, language: str, force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = new_axes()
    all_values: list[np.ndarray] = []
    for eos in EOS_ORDER:
        rows = data[eos]["reservoir"]
        time = numeric(rows, "time")
        p_min = numeric(rows, "pressure_min")
        p_avg = numeric(rows, "pressure_avg")
        p_max = numeric(rows, "pressure_max")
        all_values.extend((p_min, p_max))
        ax.fill_between(time, p_min, p_max, color=COLORS[eos], alpha=0.09, linewidth=0)
        ax.plot(
            time, p_avg, color=COLORS[eos], linestyle=LINESTYLES[eos],
            marker=MARKERS[eos], markevery=4, markerfacecolor="white",
            markeredgewidth=1.2, label=EOS_LABEL[eos], zorder=3,
        )
    ax.set_ylim(*padded_limits(all_values, 0.06))
    ax.set_xlabel(text["time"])
    ax.set_ylabel(text["pressure"])
    ax.legend(loc="lower left", frameon=True, framealpha=0.96, edgecolor="#B0B0B0")
    return export(fig, "01_pressure_response", language, force)


def final_saturation_figure(data, language: str, force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = new_axes()
    # Keep a dedicated in-frame whitespace column for the phase legend.
    ax.set_xlim(-0.55, 3.55)
    ax.set_xticks(np.arange(3), [EOS_LABEL[eos] for eos in EOS_ORDER])
    phases = (
        ("Oil", "s_o_avg", "#CC6677", ""),
        ("Gas", "s_g_avg", "#4477AA", "//"),
        ("Water", "s_w_avg", "#228833", ".."),
    )
    bottom = np.zeros(3)
    for phase, column, color, hatch in phases:
        values = np.asarray([
            float(data[eos]["reservoir"][-1][column]) for eos in EOS_ORDER
        ])
        bars = ax.bar(
            np.arange(3), values, width=0.58, bottom=bottom, label=phase,
            color=color, edgecolor="white", linewidth=0.8, hatch=hatch,
        )
        for bar, value, base in zip(bars, values, bottom):
            label_color = "white" if value >= 0.12 else "#101010"
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                base + value / 2.0,
                f"{value:.3f}",
                ha="center", va="center", fontsize=9.0, color=label_color,
            )
        bottom += values
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(text["saturation"])
    ax.legend(loc="upper right", frameon=True, framealpha=0.96, edgecolor="#B0B0B0")
    return export(fig, "02_final_phase_saturation", language, force)


def production_figure(data, language: str, force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = new_axes()
    plot_eos_lines(
        ax, data, "producer", "time", "q_total_surface",
        transform=lambda values: -values,
    )
    values = [-numeric(data[eos]["producer"], "q_total_surface") for eos in EOS_ORDER]
    ax.set_ylim(*padded_limits(values, 0.08))
    ax.set_xlabel(text["time"])
    ax.set_ylabel(text["production"])
    ax.legend(loc="upper right", frameon=True, framealpha=0.96, edgecolor="#B0B0B0")
    return export(fig, "03_producer_surface_rate", language, force)


def water_density_figure(data, language: str, force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = new_axes()
    plot_eos_lines(ax, data, "reservoir", "time", "rho_w_avg")
    values = [numeric(data[eos]["reservoir"], "rho_w_avg") for eos in EOS_ORDER]
    ax.set_ylim(*padded_limits(values, 0.08))
    ax.set_xlabel(text["time"])
    ax.set_ylabel(text["water_density"])
    ax.legend(loc="center right", frameon=True, framealpha=0.96, edgecolor="#B0B0B0")
    return export(fig, "04_water_phase_density", language, force)


def injector_bhp_figure(data, language: str, force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = new_axes()
    plot_eos_lines(ax, data, "injector", "time", "bhp", skip_first=True)
    values = [numeric(data[eos]["injector"], "bhp")[1:] for eos in EOS_ORDER]
    ax.set_ylim(*padded_limits(values, 0.10))
    # The data end at day 5; reserve the right-hand in-frame whitespace for
    # the legend so it cannot obscure any BHP curve.
    ax.set_xlim(0.0, 6.5)
    ax.set_xticks(np.arange(0.0, 5.1, 1.0))
    ax.set_xlabel(text["time"])
    ax.set_ylabel(text["injector_bhp"])
    ax.legend(loc="lower right", frameon=True, framealpha=0.96, edgecolor="#B0B0B0")
    return export(fig, "05_injector_bhp", language, force)


def solver_time_figure(data, language: str, force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = new_axes()
    for eos in EOS_ORDER:
        rows = data[eos]["solves"]
        time = np.concatenate(([0.0], numeric(rows, "end_time")))
        cumulative = np.concatenate(([0.0], np.cumsum(numeric(rows, "wall_time"))))
        ax.plot(
            time, cumulative, color=COLORS[eos], linestyle=LINESTYLES[eos],
            marker=MARKERS[eos], markevery=4, markerfacecolor="white",
            markeredgewidth=1.2, label=EOS_LABEL[eos], zorder=3,
        )
    ax.set_ylim(0.0, None)
    ax.set_xlabel(text["time"])
    ax.set_ylabel(text["solver_time"])
    ax.legend(loc="upper left", frameon=True, framealpha=0.96, edgecolor="#B0B0B0")
    return export(fig, "06_cumulative_solver_time", language, force)


def simulation_loop_wall(path: Path) -> float:
    value = float(read_rows(path / "simulation_summary.csv")[0]["simulation_loop_wall"])
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"invalid simulation loop wall time in {path}")
    return value


def cpa_algorithm_figure(language: str, force: bool) -> list[Path]:
    text = labels(language)
    fig, ax = new_axes()
    values = np.asarray([
        simulation_loop_wall(BASELINE_CPA),
        simulation_loop_wall(OPTIMIZED_CPA),
        simulation_loop_wall(RESULTS / "cpa"),
    ])
    categories = (
        ("Original fixed-point", "B2-reduced", "Current backend")
        if language == "en"
        else ("原始定点法", "B2 约化", "当前后端")
    )
    colors = ("#666666", "#56B4E9", COLORS["cpa"])
    markers = ("s", "o", "D")
    y = np.arange(3)
    for index, value in enumerate(values):
        ax.plot(
            value, y[index], marker=markers[index], markersize=8.0,
            markerfacecolor="white", markeredgecolor=colors[index],
            markeredgewidth=2.0, linestyle="none", zorder=3,
        )
        label_on_left = value == float(np.max(values))
        ax.annotate(
            f"{value:.2f} s", xy=(value, y[index]),
            xytext=(-8 if label_on_left else 8, 0), textcoords="offset points",
            va="center", ha="right" if label_on_left else "left", color="#202020",
        )
    ax.set_xscale("log")
    ax.set_xlim(40.0, 1600.0)
    ax.set_ylim(-0.65, 2.65)
    ax.set_yticks(y, categories)
    ax.set_xlabel(text["runtime"])
    ax.set_ylabel("")
    ax.tick_params(axis="y", length=0)
    return export(fig, "07_cpa_algorithm_runtime", language, force)


def write_source_table(data, force: bool) -> Path:
    path = FIGURES / "source_data.csv"
    if path.exists() and not force:
        raise FileExistsError(f"refusing to overwrite: {path} (use --force)")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "eos", "time_day", "pressure_min_bar", "pressure_avg_bar", "pressure_max_bar",
        "s_o_avg", "s_g_avg", "s_w_avg", "rho_w_avg_kg_m3",
        "producer_surface_rate_magnitude_m3_s", "injector_bhp_bar",
        "cumulative_nonlinear_solve_wall_s",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for eos in EOS_ORDER:
            reservoir = data[eos]["reservoir"]
            producer = data[eos]["producer"]
            injector = data[eos]["injector"]
            cumulative = np.concatenate((
                [0.0], np.cumsum(numeric(data[eos]["solves"], "wall_time"))
            ))
            for index, row in enumerate(reservoir):
                writer.writerow({
                    "eos": EOS_LABEL[eos],
                    "time_day": row["time"],
                    "pressure_min_bar": row["pressure_min"],
                    "pressure_avg_bar": row["pressure_avg"],
                    "pressure_max_bar": row["pressure_max"],
                    "s_o_avg": row["s_o_avg"],
                    "s_g_avg": row["s_g_avg"],
                    "s_w_avg": row["s_w_avg"],
                    "rho_w_avg_kg_m3": row["rho_w_avg"],
                    "producer_surface_rate_magnitude_m3_s":
                        -float(producer[index]["q_total_surface"]),
                    "injector_bhp_bar": injector[index]["bhp"] if index > 0 else "",
                    "cumulative_nonlinear_solve_wall_s": cumulative[index],
                })
    return path


def write_algorithm_table(force: bool) -> Path:
    path = FIGURES / "cpa_algorithm_runtime.csv"
    if path.exists() and not force:
        raise FileExistsError(f"refusing to overwrite: {path} (use --force)")
    records = (
        ("original_fixed_point", BASELINE_CPA),
        ("b2_reduced", OPTIMIZED_CPA),
        ("current_backend", RESULTS / "cpa"),
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("algorithm_stage", "simulation_loop_wall_s", "source"))
        for stage, source in records:
            writer.writerow((stage, simulation_loop_wall(source), source_label(source)))
    return path


def main() -> None:
    global RESULTS, FIGURES, BASELINE_CPA, OPTIMIZED_CPA
    args = parse_args()
    RESULTS = resolve_case_path(args.results_root)
    FIGURES = resolve_case_path(args.figures_root)
    BASELINE_CPA = resolve_case_path(args.baseline_cpa)
    OPTIMIZED_CPA = resolve_case_path(args.optimized_cpa)
    data = load_data()
    generated: list[Path] = [
        write_source_table(data, args.force),
        write_algorithm_table(args.force),
    ]
    figure_builders = (
        pressure_figure,
        final_saturation_figure,
        production_figure,
        water_density_figure,
        injector_bhp_figure,
        solver_time_figure,
    )
    for language in ("en", "zh"):
        configure_fonts(language)
        for builder in figure_builders:
            generated.extend(builder(data, language, args.force))
        generated.extend(cpa_algorithm_figure(language, args.force))
    manifest = {
        "source_root": source_label(RESULTS),
        "source_files": [
            source_label(RESULTS / eos / name)
            for eos in EOS_ORDER
            for name in (
                "reservoir_diagnostics.csv", "well_history.csv",
                "nonlinear_solve_history.csv",
            )
        ] + [
            source_label(BASELINE_CPA / "simulation_summary.csv"),
            source_label(OPTIMIZED_CPA / "simulation_summary.csv"),
            source_label(RESULTS / "cpa" / "simulation_summary.csv"),
        ],
        "transformations": [
            "producer q_total_surface multiplied by -1 to show positive production magnitude",
            "injector BHP t=0 excluded because it is an unconstrained initial guess",
            "nonlinear solve wall_time cumulatively summed by accepted end_time",
            "pressure min/max shown as translucent spatial range around arithmetic cell mean",
            "CPA algorithm runtime shown as unaltered simulation-loop wall time on a logarithmic point plot",
        ],
        "missing_data": "none in source CSV; injector t=0 intentionally blank in source_data.csv",
        "figure_size_mm": [WIDTH_MM, HEIGHT_MM],
        "png_dpi": DPI,
        "palette": COLORS,
        "generated": [str(path.relative_to(ROOT)) for path in generated],
    }
    manifest_path = FIGURES / "figure_manifest.json"
    if manifest_path.exists() and not args.force:
        raise FileExistsError(f"refusing to overwrite: {manifest_path} (use --force)")
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"generated {len(generated) - 2} figure files and 2 source tables")


if __name__ == "__main__":
    main()
