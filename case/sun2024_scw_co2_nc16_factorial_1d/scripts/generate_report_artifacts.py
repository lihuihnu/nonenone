#!/usr/bin/env python3
"""Build auditable tables and publication-style figures for the Sun et al. case.

The script reads only the explicitly supplied newly generated run directories plus
the case-local experimental manifest and physical validation targets. It never
searches for or consumes pre-existing numerical experiment data.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties, fontManager


RUNS = ("R04", "R06", "R10", "R12", "C23", "C24")
PHYSICAL_RUNS = ("R04", "R06", "R10", "R12")
RUN_LABELS = {
    "R04": "R04  23 MPa, SCW",
    "R06": "R06  24 MPa, SCW",
    "R10": "R10  23 MPa, SCW+CO2",
    "R12": "R12  24 MPa, SCW+CO2",
    "C23": "C23  23 MPa, rate-matched SCW",
    "C24": "C24  24 MPa, rate-matched SCW",
}
COLORS = {
    "R04": "#6B6B6B",
    "R06": "#6B6B6B",
    "C23": "#0072B2",
    "C24": "#0072B2",
    "R10": "#D55E00",
    "R12": "#D55E00",
}
LINESTYLES = {"R04": "--", "R06": "--", "C23": "-.", "C24": "-.", "R10": "-", "R12": "-"}
MARKERS = {"R04": "o", "R06": "o", "C23": "s", "C24": "s", "R10": "^", "R12": "^"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: str) -> float:
    return float(value)


def configure_style() -> None:
    font_candidates = [
        Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\simhei.ttf"),
    ]
    for font_path in font_candidates:
        if font_path.exists():
            fontManager.addfont(str(font_path))
            family = FontProperties(fname=str(font_path)).get_name()
            break
    else:
        family = "DejaVu Sans"
    mpl.rcParams.update(
        {
            "font.family": family,
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.6,
            "lines.markersize": 4,
            "figure.dpi": 120,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )


def save_figure(fig: mpl.figure.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=600)
    fig.savefig(base.with_suffix(".pdf"))
    plt.close(fig)


def load_run(run_dir: Path) -> dict:
    meta = read_csv(run_dir / "experiment_metadata.csv")[0]
    dt_pv = as_float(meta["dt_pv"])
    balance = read_csv(run_dir / "component_mass_balance.csv")
    oil_rows = [r for r in balance if r["component"] == "nC16"]
    recovery = []
    for row in oil_rows:
        recovery.append(
            {
                "step": int(row["step"]),
                "pvi": int(row["step"]) * dt_pv,
                "recovery_fraction": as_float(row["cumulative_produced_kg"])
                / as_float(row["initial_inventory_kg"]),
            }
        )

    well = read_csv(run_dir / "well_history.csv")
    by_step: dict[int, dict[str, dict[str, str]]] = {}
    for row in well:
        by_step.setdefault(int(row["step"]), {})[row["name"]] = row
    pressure = []
    for step, rows in sorted(by_step.items()):
        if "SCW_CO2_INJ" not in rows or "PROD" not in rows:
            continue
        inj = rows["SCW_CO2_INJ"]
        prod = rows["PROD"]
        pressure.append(
            {
                "step": step,
                "pvi": step * dt_pv,
                "injector_bhp_MPa": as_float(inj["bhp"]) / 10.0,
                "producer_bhp_MPa": as_float(prod["bhp"]) / 10.0,
                "delta_p_MPa": (as_float(inj["bhp"]) - as_float(prod["bhp"])) / 10.0,
            }
        )

    max_closure = 0.0
    for row in balance:
        scale = max(
            abs(as_float(row["initial_inventory_kg"])),
            abs(as_float(row["current_inventory_kg"])),
            abs(as_float(row["cumulative_injected_kg"])),
            abs(as_float(row["cumulative_produced_kg"])),
            1.0e-30,
        )
        max_closure = max(max_closure, abs(as_float(row["balance_error_kg"])) / scale)

    # Component breakthrough is based on the finite-difference rate of cumulative
    # produced CO2, consistent with check_convergence.py. Phase well rates cannot
    # be used here because dissolved CO2 is not equivalent to the gas-phase rate.
    inlet_co2 = as_float(meta["co2_mass_rate_kg_s"])
    ratios: list[tuple[float, float]] = []
    if inlet_co2 > 0.0:
        co2_rows = [r for r in balance if r["component"] == "CO2"]
        total_rate = as_float(meta["reservoir_total_rate_m3_s"])
        pore_volume = as_float(meta["pore_volume_m3"])
        for previous, row in zip(co2_rows, co2_rows[1:]):
            dt_seconds = (as_float(row["time_day"]) - as_float(previous["time_day"])) * 86400.0
            if dt_seconds <= 0.0:
                continue
            produced_rate = max(
                0.0,
                (as_float(row["cumulative_produced_kg"]) - as_float(previous["cumulative_produced_kg"]))
                / dt_seconds,
            )
            pvi = as_float(row["time_day"]) * 86400.0 * total_rate / pore_volume
            ratios.append((pvi, produced_rate / inlet_co2))
    first_breakthrough = next((pvi for pvi, ratio in ratios if ratio >= 0.01), None)
    max_ratio = max((ratio for _, ratio in ratios), default=0.0)

    return {
        "meta": meta,
        "recovery": recovery,
        "pressure": pressure,
        "endpoint_recovery": recovery[-1]["recovery_fraction"],
        "max_delta_p_MPa": max(r["delta_p_MPa"] for r in pressure),
        "final_delta_p_MPa": pressure[-1]["delta_p_MPa"],
        "max_mass_closure_relative": max_closure,
        "breakthrough_pvi_1pct": first_breakthrough,
        "max_outlet_inlet_co2_mass_ratio": max_ratio,
    }


def plot_curves(data: dict[str, dict], figures: Path, key: str) -> None:
    ylabel = "n-C16 累计采收率（%）" if key == "recovery" else "注采压差（MPa）"
    filename = "fig01_recovery_curves" if key == "recovery" else "fig02_pressure_drop_curves"
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.05), sharey=False)
    groups = [("23 MPa", ("R04", "C23", "R10")), ("24 MPa", ("R06", "C24", "R12"))]
    for ax, (title, run_ids) in zip(axes, groups):
        for run_id in run_ids:
            if key == "recovery":
                x = [r["pvi"] for r in data[run_id]["recovery"]]
                y = [100.0 * r["recovery_fraction"] for r in data[run_id]["recovery"]]
            else:
                x = [r["pvi"] for r in data[run_id]["pressure"]]
                y = [r["delta_p_MPa"] for r in data[run_id]["pressure"]]
            markevery = max(1, len(x) // 10)
            ax.plot(
                x,
                y,
                color=COLORS[run_id],
                linestyle=LINESTYLES[run_id],
                marker=MARKERS[run_id],
                markevery=markevery,
                markerfacecolor="white",
                markeredgewidth=0.9,
                label=RUN_LABELS[run_id].split("  ", 1)[1],
            )
        ax.set_title(title)
        ax.set_xlabel("注入孔隙体积 PVI（–）")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 4)
        if key == "recovery":
            ax.set_ylim(bottom=0)
        ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.75)
        ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    save_figure(fig, figures / filename)


def plot_endpoint_and_effects(
    data: dict[str, dict], physical: dict[str, float], effects: list[dict], figures: Path
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.35))
    ax = axes[0]
    x = np.arange(len(PHYSICAL_RUNS))
    width = 0.36
    modeled = [100.0 * data[r]["endpoint_recovery"] for r in PHYSICAL_RUNS]
    observed = [physical[r] for r in PHYSICAL_RUNS]
    ax.bar(x - width / 2, modeled, width, color="#0072B2", label="当前等温代理模型")
    ax.bar(x + width / 2, observed, width, color="#E69F00", hatch="///", label="原论文物理实验")
    ax.set_xticks(x, PHYSICAL_RUNS)
    ax.set_ylabel("4 PV 位移效率/采收率（%）")
    ax.set_ylim(0, 100)
    ax.set_title("(a) 外部端点面对比（非标定）", pad=27)
    ax.grid(True, axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.75)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.07), ncol=2)

    ax = axes[1]
    labels = ["增流量\nSCW", "等流量\nCO2组分", "合并变化\nCO2+增流量"]
    c23 = [100.0 * e["effect_23MPa_fraction"] for e in effects]
    c24 = [100.0 * e["effect_24MPa_fraction"] for e in effects]
    x = np.arange(3)
    ax.bar(x - width / 2, c23, width, color="#56B4E9", label="23 MPa")
    ax.bar(x + width / 2, c24, width, color="#009E73", hatch="\\\\", label="24 MPa")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("4 PV 采收率差（百分点）")
    ax.set_title("(b) 数值效应分解")
    ax.grid(True, axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.75)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    save_figure(fig, figures / "fig03_endpoint_validation_and_effects")


def load_profile(
    run_dir: Path, step: int, nx: int, outlet_pressure_MPa: float, length_m: float = 0.48
) -> list[dict]:
    rows = read_csv(run_dir / f"solution_step_{step}.csv")
    profile = []
    for row in rows:
        index = int(row["input_index"])
        water = as_float(row["water_saturation"])
        liquid = as_float(row["liquid_saturation"])
        vapor = as_float(row["vapor_saturation"])
        profile.append(
            {
                "step": step,
                "pvi": step * 0.01,
                "x_m": (index + 0.5) / nx * length_m,
                "water_saturation": water,
                "oil_saturation": liquid,
                "vapor_saturation": vapor,
                "liquid_co2_mole_fraction": as_float(row["liquid_x_0_CO2"]),
                "pressure_MPa": as_float(row["pressure_Pa"]) / 1.0e6,
                "pressure_excess_kPa":
                    (as_float(row["pressure_Pa"]) / 1.0e6 - outlet_pressure_MPa) * 1000.0,
            }
        )
    return profile


def plot_profiles(profile_rows: list[dict], figures: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 5.45), sharex=True)
    panels = [
        ("water_saturation", "水相饱和度（–）", "(a) SCW 饱和度"),
        ("oil_saturation", "油相饱和度（–）", "(b) n-C16 相饱和度"),
        ("liquid_co2_mole_fraction", "液相 CO2 摩尔分数（–）", "(c) 液相 CO2 组分"),
        ("pressure_excess_kPa", "高于出口的压力（kPa）", "(d) 压力分布"),
    ]
    pvis = sorted({r["pvi"] for r in profile_rows})
    colors = ["#0072B2", "#009E73", "#E69F00", "#D55E00"]
    linestyles = [":", "--", "-.", "-"]
    markers = ["o", "s", "D", "^"]
    for ax, (key, ylabel, title) in zip(axes.flat, panels):
        for pvi, color, ls, marker in zip(pvis, colors, linestyles, markers):
            rows = [r for r in profile_rows if math.isclose(r["pvi"], pvi)]
            ax.plot(
                [r["x_m"] for r in rows],
                [r[key] for r in rows],
                color=color,
                linestyle=ls,
                marker=marker,
                markevery=max(1, len(rows) // 12),
                markerfacecolor="white",
                label=f"{pvi:g} PVI",
            )
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.75)
    for ax in axes[-1, :]:
        ax.set_xlabel("距注入口距离（m）")
    axes[0, 0].legend(frameon=False, ncol=2, loc="best")
    fig.tight_layout()
    save_figure(fig, figures / "fig04_R12_spatial_profiles")


def plot_convergence(grid_rows: list[dict], time_rows: list[dict], figures: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.25))
    ax = axes[0]
    gx = [int(float(r["nx"])) for r in grid_rows]
    gy = [100.0 * as_float(r["nc16_recovery"]) for r in grid_rows]
    ax.plot(gx, gy, color="#0072B2", marker="o", markerfacecolor="white")
    ax.set_xlabel("网格数 nx（–）")
    ax.set_ylabel("4 PV 采收率（%）")
    ax.set_title("(a) 网格收敛性")
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.75)
    ax2 = ax.twinx()
    ax2.semilogy(
        gx,
        [as_float(r["max_component_relative_closure"]) for r in grid_rows],
        color="#D55E00",
        marker="s",
        linestyle="--",
        markerfacecolor="white",
    )
    ax2.set_ylabel("最大相对质量闭合误差（–）", color="#D55E00")

    ax = axes[1]
    tx = [as_float(r["dt_pv"]) for r in time_rows]
    ty = [100.0 * as_float(r["nc16_recovery"]) for r in time_rows]
    order = np.argsort(tx)
    tx = [tx[i] for i in order]
    ty = [ty[i] for i in order]
    ax.plot(tx, ty, color="#009E73", marker="^", markerfacecolor="white")
    ax.set_xlabel("时间步长 ΔPVI（–）")
    ax.set_ylabel("4 PV 采收率（%）")
    ax.set_title("(b) 时间步收敛性")
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.75)
    ax.invert_xaxis()
    ax.text(
        0.04,
        0.06,
        "所有收敛算例均未达到\n1% CO2 出口/入口质量率阈值",
        transform=ax.transAxes,
        fontsize=7.5,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#BDBDBD"},
    )
    fig.tight_layout()
    save_figure(fig, figures / "fig05_numerical_convergence")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--formal-root", type=Path, required=True)
    parser.add_argument("--case-root", type=Path, required=True)
    parser.add_argument("--grid-summary", type=Path, required=True)
    parser.add_argument("--time-summary", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()

    figures = args.output_root / "figures"
    tables = args.output_root / "data"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    configure_style()

    data = {run: load_run(args.formal_root / run) for run in RUNS}
    recovery_rows: list[dict] = []
    pressure_rows: list[dict] = []
    endpoint_rows: list[dict] = []
    for run in RUNS:
        for row in data[run]["recovery"]:
            recovery_rows.append({"experiment_id": run, **row})
        for row in data[run]["pressure"]:
            pressure_rows.append({"experiment_id": run, **row})
        endpoint_rows.append(
            {
                "experiment_id": run,
                "role": data[run]["meta"]["role"],
                "pressure_MPa": as_float(data[run]["meta"]["pressure_Pa"]) / 1.0e6,
                "water_rate_reported_mL_min": data[run]["meta"]["water_rate_reported_mL_min"],
                "co2_rate_reported_mL_min": data[run]["meta"]["co2_rate_reported_mL_min"],
                "reservoir_total_rate_m3_s": data[run]["meta"]["reservoir_total_rate_m3_s"],
                "recovery_4PV_fraction": data[run]["endpoint_recovery"],
                "recovery_4PV_percent": 100.0 * data[run]["endpoint_recovery"],
                "max_delta_p_MPa": data[run]["max_delta_p_MPa"],
                "final_delta_p_MPa": data[run]["final_delta_p_MPa"],
                "max_relative_mass_closure": data[run]["max_mass_closure_relative"],
                "breakthrough_pvi_1pct": "not_reached" if data[run]["breakthrough_pvi_1pct"] is None else data[run]["breakthrough_pvi_1pct"],
                "max_outlet_inlet_co2_mass_ratio": data[run]["max_outlet_inlet_co2_mass_ratio"],
            }
        )
    write_csv(tables / "recovery_curves.csv", recovery_rows)
    write_csv(tables / "pressure_drop_curves.csv", pressure_rows)
    write_csv(tables / "endpoint_summary.csv", endpoint_rows)

    effects = [
        {
            "effect": "pure_SCW_flow_effect",
            "definition_23MPa": "C23-R04",
            "definition_24MPa": "C24-R06",
            "effect_23MPa_fraction": data["C23"]["endpoint_recovery"] - data["R04"]["endpoint_recovery"],
            "effect_24MPa_fraction": data["C24"]["endpoint_recovery"] - data["R06"]["endpoint_recovery"],
        },
        {
            "effect": "equal_reservoir_rate_CO2_composition_effect",
            "definition_23MPa": "R10-C23",
            "definition_24MPa": "R12-C24",
            "effect_23MPa_fraction": data["R10"]["endpoint_recovery"] - data["C23"]["endpoint_recovery"],
            "effect_24MPa_fraction": data["R12"]["endpoint_recovery"] - data["C24"]["endpoint_recovery"],
        },
        {
            "effect": "combined_CO2_plus_flow_effect",
            "definition_23MPa": "R10-R04",
            "definition_24MPa": "R12-R06",
            "effect_23MPa_fraction": data["R10"]["endpoint_recovery"] - data["R04"]["endpoint_recovery"],
            "effect_24MPa_fraction": data["R12"]["endpoint_recovery"] - data["R06"]["endpoint_recovery"],
        },
    ]
    write_csv(tables / "effect_decomposition.csv", effects)

    physical_rows = read_csv(args.case_root / "physical_validation_targets.csv")
    physical = {r["experiment_id"]: as_float(r["displacement_efficiency_percent"]) for r in physical_rows}
    comparison = []
    for run in PHYSICAL_RUNS:
        model = 100.0 * data[run]["endpoint_recovery"]
        comparison.append(
            {
                "experiment_id": run,
                "modeled_4PV_percent": model,
                "physical_experiment_percent": physical[run],
                "signed_difference_percentage_points": model - physical[run],
                "interpretation": "external_face_validity_not_calibration",
            }
        )
    write_csv(tables / "physical_endpoint_comparison.csv", comparison)

    profile_rows: list[dict] = []
    nx = int(data["R12"]["meta"]["nx"])
    outlet_pressure_MPa = as_float(data["R12"]["meta"]["pressure_Pa"]) / 1.0e6
    for step in (50, 100, 200, 400):
        profile_rows.extend(
            load_profile(args.formal_root / "R12", step, nx, outlet_pressure_MPa)
        )
    write_csv(tables / "R12_spatial_profiles.csv", profile_rows)

    grid_rows = read_csv(args.grid_summary)
    time_rows = read_csv(args.time_summary)
    write_csv(tables / "grid_convergence_summary.csv", grid_rows)
    write_csv(tables / "time_convergence_summary.csv", time_rows)

    plot_curves(data, figures, "recovery")
    plot_curves(data, figures, "pressure")
    plot_endpoint_and_effects(data, physical, effects, figures)
    plot_profiles(profile_rows, figures)
    plot_convergence(grid_rows, time_rows, figures)

    provenance = {
        "scope": "newly_generated_runs_only",
        "formal_run_root": str(args.formal_root.resolve()),
        "formal_runs": list(RUNS),
        "case_manifest": str((args.case_root / "experiment_manifest.csv").resolve()),
        "physical_targets": str((args.case_root / "physical_validation_targets.csv").resolve()),
        "grid_convergence_source": str(args.grid_summary.resolve()),
        "time_convergence_source": str(args.time_summary.resolve()),
        "transformations": {
            "pvi": "step * dt_pv",
            "recovery": "cumulative_produced_nC16 / initial_nC16_inventory",
            "pressure_drop": "injector_bhp - producer_bhp",
            "breakthrough": "first sampled PVI with abs(producer_CO2_mass_rate)/injector_CO2_mass_rate >= 0.01",
            "effects": ["C23-R04 / C24-R06", "R10-C23 / R12-C24", "R10-R04 / R12-R06"],
        },
        "rendering": {
            "software": f"Matplotlib {mpl.__version__}",
            "raster": "PNG 600 dpi",
            "vector": "PDF",
            "smoothing": "none",
            "palette": "Okabe-Ito-derived, with redundant line styles and markers",
        },
        "limitations": [
            "The model is isothermal and does not solve an energy equation or a 50–400 °C thermal front.",
            "n-C16 is a surrogate oil component; upgrading and reactions are not represented.",
            "The source did not resolve the injection-rate reference state; these formal runs use in_situ_volume.",
        ],
    }
    (tables / "figure_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(json.dumps({"endpoint_summary": endpoint_rows, "effects": effects}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
