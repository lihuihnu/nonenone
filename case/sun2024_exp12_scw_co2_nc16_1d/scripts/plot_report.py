#!/usr/bin/env python3
"""Generate report figures and derived tables for the Sun-2024 SCW case."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


CASE_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = CASE_DIR / "results" / "report_4pv"
REPORT_DIR = CASE_DIR / "report"
FIGURE_DIR = REPORT_DIR / "figures"
DATA_DIR = REPORT_DIR / "data"

PORE_VOLUME_M3 = 0.0002236273615419518
TOTAL_RATE_M3_S = 2.0e-7
SECONDS_PER_DAY = 86400.0
EXPERIMENTAL_RECOVERY_PERCENT = 82.96

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000"]
LINESTYLES = ["-", "--", "-.", ":", "-"]
MARKERS = ["o", "s", "^", "D", "x"]


def pv_from_days(time_day: pd.Series | np.ndarray | float) -> np.ndarray:
    return np.asarray(time_day, dtype=float) * SECONDS_PER_DAY * TOTAL_RATE_M3_S / PORE_VOLUME_M3


def first_crossing(x: np.ndarray, y: np.ndarray, threshold: float) -> float | None:
    """Linearly interpolate the first recorded threshold crossing."""
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    for index in range(1, len(x)):
        if y[index - 1] < threshold <= y[index]:
            if y[index] == y[index - 1]:
                return float(x[index])
            fraction = (threshold - y[index - 1]) / (y[index] - y[index - 1])
            return float(x[index - 1] + fraction * (x[index] - x[index - 1]))
    return None


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
            "xtick.major.size": 4.0,
            "ytick.major.size": 4.0,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
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
        fontsize=11,
        fontweight="bold",
    )


def export(fig: mpl.figure.Figure, stem: str) -> None:
    png_path = FIGURE_DIR / f"{stem}.png"
    fig.savefig(png_path, dpi=450)
    fig.savefig(FIGURE_DIR / f"{stem}.pdf")
    plt.close(fig)
    with Image.open(png_path) as image:
        image.convert("RGB").save(png_path, dpi=(450, 450))


def load_performance() -> tuple[pd.DataFrame, dict[str, float | int | None]]:
    mass = pd.read_csv(RAW_DIR / "component_mass_balance.csv")
    wells = pd.read_csv(RAW_DIR / "well_history.csv")
    diagnostics = pd.read_csv(RAW_DIR / "reservoir_diagnostics.csv")
    summary = pd.read_csv(RAW_DIR / "simulation_summary.csv").iloc[0]

    nc16 = mass[mass["component"] == "nC16"].copy().sort_values("step")
    co2 = mass[mass["component"] == "CO2"].copy().sort_values("step")
    grouped_balance = (
        mass.assign(abs_relative_error=mass["relative_error"].abs())
        .groupby(["step", "time_day"], as_index=False)["abs_relative_error"]
        .max()
    )
    injector = wells[wells["name"] == "SCW_CO2_INJ"].copy().sort_values("step")
    producer = wells[wells["name"] == "PROD"].copy().sort_values("step")

    history = pd.DataFrame(
        {
            "step": nc16["step"].to_numpy(),
            "time_day": nc16["time_day"].to_numpy(),
        }
    )
    history["injected_pv"] = pv_from_days(history["time_day"])
    history["nc16_recovery_percent"] = (
        nc16["cumulative_produced_kg"].to_numpy()
        / nc16["initial_inventory_kg"].to_numpy()
        * 100.0
    )
    oil_rate = np.abs(producer["q_oil_reservoir"].to_numpy())
    gas_rate = np.abs(producer["q_gas_reservoir"].to_numpy())
    water_rate = np.abs(producer["q_water_reservoir"].to_numpy())
    total_rate = oil_rate + gas_rate + water_rate
    history["producer_water_cut_percent"] = np.divide(
        water_rate * 100.0,
        total_rate,
        out=np.full_like(total_rate, np.nan),
        where=total_rate > 0.0,
    )
    history["pressure_drop_kpa"] = (
        injector["bhp"].to_numpy() - producer["bhp"].to_numpy()
    ) * 100.0
    injected_co2 = co2["cumulative_injected_kg"].to_numpy()
    history["co2_retained_percent_of_injected"] = np.divide(
        co2["current_inventory_kg"].to_numpy() * 100.0,
        injected_co2,
        out=np.full_like(injected_co2, np.nan),
        where=injected_co2 > 0.0,
    )
    history["co2_produced_percent_of_injected"] = np.divide(
        co2["cumulative_produced_kg"].to_numpy() * 100.0,
        injected_co2,
        out=np.full_like(injected_co2, np.nan),
        where=injected_co2 > 0.0,
    )
    history["max_component_relative_balance_error"] = grouped_balance[
        "abs_relative_error"
    ].to_numpy()
    history["average_oil_saturation"] = diagnostics["s_o_avg"].to_numpy()
    history["average_water_saturation"] = diagnostics["s_w_avg"].to_numpy()
    history.to_csv(DATA_DIR / "performance_history.csv", index=False)

    final_nc16 = nc16.iloc[-1]
    final_co2 = co2.iloc[-1]
    final_diag = diagnostics.iloc[-1]
    final_history = history.iloc[-1]
    metrics: dict[str, float | int | None] = {
        "final_time_day": float(summary["final_time_day"]),
        "final_injected_pv": float(final_history["injected_pv"]),
        "final_nc16_recovery_percent": float(final_history["nc16_recovery_percent"]),
        "experimental_recovery_percent": EXPERIMENTAL_RECOVERY_PERCENT,
        "recovery_difference_percentage_points": float(
            final_history["nc16_recovery_percent"] - EXPERIMENTAL_RECOVERY_PERCENT
        ),
        "final_nc16_remaining_kg": float(final_nc16["current_inventory_kg"]),
        "final_water_cut_percent": float(final_history["producer_water_cut_percent"]),
        "water_cut_50_percent_pv": first_crossing(
            history["injected_pv"].to_numpy(),
            history["producer_water_cut_percent"].to_numpy(),
            50.0,
        ),
        "water_cut_90_percent_pv": first_crossing(
            history["injected_pv"].to_numpy(),
            history["producer_water_cut_percent"].to_numpy(),
            90.0,
        ),
        "recovery_80_percent_pv": first_crossing(
            history["injected_pv"].to_numpy(),
            history["nc16_recovery_percent"].to_numpy(),
            80.0,
        ),
        "recovery_experimental_endpoint_pv": first_crossing(
            history["injected_pv"].to_numpy(),
            history["nc16_recovery_percent"].to_numpy(),
            EXPERIMENTAL_RECOVERY_PERCENT,
        ),
        "final_pressure_drop_kpa": float(final_history["pressure_drop_kpa"]),
        "maximum_pressure_drop_kpa": float(history["pressure_drop_kpa"].max()),
        "final_co2_injected_kg": float(final_co2["cumulative_injected_kg"]),
        "final_co2_produced_kg": float(final_co2["cumulative_produced_kg"]),
        "final_co2_retained_percent_of_injected": float(
            final_history["co2_retained_percent_of_injected"]
        ),
        "final_co2_produced_percent_of_injected": float(
            final_history["co2_produced_percent_of_injected"]
        ),
        "final_average_oil_saturation": float(final_diag["s_o_avg"]),
        "final_average_water_saturation": float(final_diag["s_w_avg"]),
        "final_max_component_relative_balance_error": float(
            final_history["max_component_relative_balance_error"]
        ),
        "maximum_component_relative_balance_error": float(
            history["max_component_relative_balance_error"].max()
        ),
        "accepted_internal_steps": int(summary["accepted_internal_steps"]),
        "rejected_internal_steps": int(summary["rejected_internal_steps"]),
        "maximum_snes_iterations": int(summary["max_snes_per_solve"]),
        "simulation_loop_wall_s": float(summary["simulation_loop_wall"]),
    }
    (DATA_DIR / "summary_metrics.json").write_text(
        json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return history, metrics


def plot_spatial_profiles() -> None:
    selected_steps = [0, 20, 40, 60, 80]
    profiles: list[tuple[float, pd.DataFrame]] = []
    for step in selected_steps:
        frame = pd.read_csv(RAW_DIR / f"solution_step_{step}.csv")
        frame["x_m"] = (frame["input_index"] + 0.5) * 0.48 / 48.0
        profiles.append((step * 4.0 / 80.0, frame))

    fig, axes = plt.subplots(
        2, 2, figsize=(180 / 25.4, 125 / 25.4), layout="constrained"
    )
    ax_sw, ax_so, ax_p, ax_co2 = axes.flat
    for index, (pv, frame) in enumerate(profiles):
        style = {
            "color": COLORS[index],
            "linestyle": LINESTYLES[index],
            "marker": MARKERS[index],
            "markevery": 8,
            "markersize": 3.6,
            "linewidth": 1.6,
            "label": f"{pv:g} PV",
        }
        ax_sw.plot(frame["x_m"], frame["water_saturation"], **style)
        ax_so.plot(frame["x_m"], frame["liquid_saturation"], **style)
        outlet_pressure = frame["pressure_Pa"].iloc[-1]
        ax_p.plot(
            frame["x_m"],
            (frame["pressure_Pa"] - outlet_pressure) / 1.0e3,
            **style,
        )
        ax_co2.plot(frame["x_m"], frame["liquid_x_0_CO2"], **style)

    ax_sw.set(ylabel=r"超临界水饱和度 $S_w$", ylim=(0.0, 1.0))
    ax_so.set(ylabel=r"油相饱和度 $S_o$", ylim=(0.0, 1.0))
    ax_p.set(ylabel=r"相对生产端压力 $p-p_{\mathrm{prod}}$ (kPa)")
    ax_co2.set(ylabel=r"油相 CO$_2$ 摩尔分数 $x_{o,\mathrm{CO_2}}$", ylim=(0.0, 0.35))
    for label, ax in zip(["(a)", "(b)", "(c)", "(d)"], axes.flat):
        panel_label(ax, label)
        ax.set_xlabel(r"轴向位置 $x$ (m)")
        ax.set_xlim(0.0, 0.48)
        for spine in ax.spines.values():
            spine.set_visible(True)

    handles, labels = ax_sw.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="outside lower center",
        ncol=5,
        handlelength=2.6,
        columnspacing=1.4,
    )
    export(fig, "01_spatial_profiles")


def plot_performance(history: pd.DataFrame, metrics: dict[str, float | int | None]) -> None:
    fig, axes = plt.subplots(
        2, 2, figsize=(180 / 25.4, 125 / 25.4), layout="constrained"
    )
    ax_rf, ax_wc, ax_dp, ax_co2 = axes.flat
    pv = history["injected_pv"]

    ax_rf.plot(
        pv,
        history["nc16_recovery_percent"],
        color=COLORS[0],
        linewidth=1.8,
        marker="o",
        markevery=8,
        markersize=3.8,
        label="数值结果",
    )
    ax_rf.axhline(
        EXPERIMENTAL_RECOVERY_PERCENT,
        color=COLORS[1],
        linestyle="--",
        linewidth=1.4,
        label="实验终点 82.96%",
    )
    ax_rf.scatter(
        [metrics["final_injected_pv"]],
        [metrics["final_nc16_recovery_percent"]],
        color=COLORS[0],
        edgecolor="white",
        linewidth=0.6,
        zorder=4,
    )
    ax_rf.set(ylabel="nC16 采收率 (%)", ylim=(0.0, 100.0))
    ax_rf.legend(loc="lower right")

    ax_wc.plot(
        pv,
        history["producer_water_cut_percent"],
        color=COLORS[2],
        linestyle="-",
        linewidth=1.8,
        marker="s",
        markevery=8,
        markersize=3.6,
    )
    ax_wc.set(ylabel="生产井含水率 (%)", ylim=(0.0, 102.0))

    ax_dp.plot(
        pv,
        history["pressure_drop_kpa"],
        color=COLORS[3],
        linestyle="-.",
        linewidth=1.8,
        marker="^",
        markevery=8,
        markersize=3.8,
    )
    ax_dp.set(ylabel=r"注采压差 $\Delta p$ (kPa)")

    ax_co2.plot(
        pv,
        history["co2_retained_percent_of_injected"],
        color=COLORS[0],
        linewidth=1.7,
        label="储层内滞留",
    )
    ax_co2.plot(
        pv,
        history["co2_produced_percent_of_injected"],
        color=COLORS[1],
        linestyle="--",
        linewidth=1.7,
        label="累计采出",
    )
    ax_co2.set(ylabel="占累计注入 CO$_2$ 比例 (%)", ylim=(-2.0, 102.0))
    ax_co2.legend(loc="center right")

    for label, ax in zip(["(a)", "(b)", "(c)", "(d)"], axes.flat):
        panel_label(ax, label)
        ax.set_xlabel("累计注入量 (PV)")
        ax.set_xlim(0.0, 4.0)
        for spine in ax.spines.values():
            spine.set_visible(True)
    export(fig, "02_performance_and_co2")


def plot_numerical_quality(history: pd.DataFrame) -> None:
    steps = pd.read_csv(RAW_DIR / "time_step_history.csv")
    accepted = steps[steps["status"] == "ACCEPT"].copy()
    rejected = steps[steps["status"] == "REJECT"].copy()
    accepted["pv"] = pv_from_days(accepted["time"])
    rejected["pv"] = pv_from_days(rejected["time"])

    fig, axes = plt.subplots(
        1, 2, figsize=(180 / 25.4, 64 / 25.4), layout="constrained"
    )
    ax_dt, ax_mb = axes
    ax_dt.step(
        accepted["pv"],
        accepted["dt"] * SECONDS_PER_DAY,
        where="post",
        color=COLORS[0],
        linewidth=1.5,
        label="接受时间步",
    )
    if not rejected.empty:
        ax_dt.scatter(
            rejected["pv"],
            rejected["dt"] * SECONDS_PER_DAY,
            color=COLORS[1],
            marker="x",
            s=34,
            linewidth=1.6,
            label="拒绝时间步",
            zorder=4,
        )
    ax_dt.set(xlabel="累计注入量 (PV)", ylabel="内部时间步长 (s)", xlim=(0.0, 4.0))
    ax_dt.legend(loc="lower right")

    ax_mb.plot(
        history["injected_pv"],
        history["max_component_relative_balance_error"] / 1.0e-13,
        color=COLORS[4],
        linewidth=1.6,
        marker="D",
        markevery=8,
        markersize=3.2,
    )
    ax_mb.set(
        xlabel="累计注入量 (PV)",
        ylabel=r"最大组分相对守恒误差 ($10^{-13}$)",
        xlim=(0.0, 4.0),
    )
    for label, ax in zip(["(a)", "(b)"], axes):
        panel_label(ax, label)
        for spine in ax.spines.values():
            spine.set_visible(True)
    export(fig, "03_numerical_quality")


def write_manifest() -> None:
    source_files = [
        "component_mass_balance.csv",
        "well_history.csv",
        "reservoir_diagnostics.csv",
        "simulation_summary.csv",
        "time_step_history.csv",
        *(f"solution_step_{step}.csv" for step in [0, 20, 40, 60, 80]),
    ]
    manifest = {
        "generated_by": str(Path(__file__).relative_to(CASE_DIR)),
        "python": platform.python_version(),
        "matplotlib": mpl.__version__,
        "pandas": pd.__version__,
        "raw_data_directory": str(RAW_DIR.relative_to(CASE_DIR)),
        "source_files": source_files,
        "transformations": [
            "injected PV = time_day * 86400 * 2e-7 / 0.0002236273615419518",
            "nC16 recovery = cumulative produced nC16 / initial nC16 inventory",
            "water cut uses absolute reservoir oil+gas+water production volumes",
            "threshold PV values use linear interpolation between adjacent output records",
            "spatial coordinate uses 48 uniform 1 cm cell centers",
            "spatial pressure = cell pressure - final-cell pressure at each snapshot",
            "no smoothing, filtering, or exclusion applied",
        ],
        "figure_size_mm": {
            "01_spatial_profiles": [180, 125],
            "02_performance_and_co2": [180, 125],
            "03_numerical_quality": [180, 64],
        },
        "exports": ["PNG 450 dpi", "PDF vector"],
        "palette": COLORS,
        "background": "#FFFFFF",
    }
    (REPORT_DIR / "figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def main() -> None:
    if not RAW_DIR.exists():
        raise FileNotFoundError(f"Missing raw result directory: {RAW_DIR}")
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    configure_style()
    history, metrics = load_performance()
    plot_spatial_profiles()
    plot_performance(history, metrics)
    plot_numerical_quality(history)
    write_manifest()
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
