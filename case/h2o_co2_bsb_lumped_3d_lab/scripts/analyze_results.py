#!/usr/bin/env python3
"""Summarize the five formal displacement runs using only standard-library Python."""

from __future__ import annotations

import csv
import html
import json
import math
import sys
from collections import defaultdict
from pathlib import Path


CASE_ROOT = Path(__file__).resolve().parents[1]
SECONDS_PER_DAY = 86400.0
HYDROCARBONS = {"L_C1_C6", "M_C7_C15", "H_C16_C27", "XH_C28plus"}
HEAVY = {"H_C16_C27", "XH_C28plus"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], key: str) -> float:
    return float(row[key])


def integrate_positive_production(
    rows: list[dict[str, str]], run_id: str
) -> tuple[float, list[dict[str, object]]]:
    producer = sorted(
        (row for row in rows if row["type"] == "PRODUCER"),
        key=lambda row: number(row, "time"),
    )
    cumulative = 0.0
    history: list[dict[str, object]] = []
    previous: dict[str, str] | None = None
    for row in producer:
        rate = max(0.0, -number(row, "q_oil_reservoir"))
        if previous is not None:
            previous_rate = max(0.0, -number(previous, "q_oil_reservoir"))
            dt_seconds = (number(row, "time") - number(previous, "time")) * SECONDS_PER_DAY
            cumulative += 0.5 * (previous_rate + rate) * dt_seconds
        history.append({
            "run_id": run_id,
            "time_day": number(row, "time"),
            "oil_rate_reservoir_m3_s": rate,
            "cumulative_oil_reservoir_m3": cumulative,
        })
        previous = row
    return cumulative, history


def breakthrough_pvi(
    balance: list[dict[str, str]], component: str, rate: float, pore_volume: float
) -> float | None:
    rows = sorted(
        (row for row in balance if row["component"] == component),
        key=lambda row: int(row["step"]),
    )
    for before, after in zip(rows, rows[1:]):
        injected = number(after, "cumulative_injected_kg") - number(
            before, "cumulative_injected_kg"
        )
        produced = number(after, "cumulative_produced_kg") - number(
            before, "cumulative_produced_kg"
        )
        if injected > 0.0 and produced >= 0.01 * injected:
            return number(after, "time_day") * SECONDS_PER_DAY * rate / pore_volume
    return None


def summarize_run(run_dir: Path, run: dict[str, str]) -> tuple[dict[str, object], list[dict[str, object]]]:
    required = [
        "design_metadata.csv",
        "component_mass_balance.csv",
        "well_history.csv",
        "reservoir_diagnostics.csv",
        "simulation_summary.csv",
        "well_control_switches.csv",
    ]
    missing = [name for name in required if not (run_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"{run['run_id']}: missing {', '.join(missing)}")

    metadata = read_rows(run_dir / "design_metadata.csv")[0]
    balance = read_rows(run_dir / "component_mass_balance.csv")
    wells = read_rows(run_dir / "well_history.csv")
    diagnostics = read_rows(run_dir / "reservoir_diagnostics.csv")
    simulation = read_rows(run_dir / "simulation_summary.csv")[0]
    switches = read_rows(run_dir / "well_control_switches.csv")

    final_step = max(int(row["step"]) for row in balance)
    final_balance = [row for row in balance if int(row["step"]) == final_step]
    initial_hc = sum(number(row, "initial_inventory_kg") for row in final_balance
                     if row["component"] in HYDROCARBONS)
    produced_hc = sum(number(row, "cumulative_produced_kg") for row in final_balance
                      if row["component"] in HYDROCARBONS)
    initial_heavy = sum(number(row, "initial_inventory_kg") for row in final_balance
                        if row["component"] in HEAVY)
    produced_heavy = sum(number(row, "cumulative_produced_kg") for row in final_balance
                         if row["component"] in HEAVY)
    hydrocarbon_recovery = produced_hc / initial_hc if initial_hc > 0.0 else math.nan
    heavy_recovery = produced_heavy / initial_heavy if initial_heavy > 0.0 else math.nan
    maximum_balance_error = max(abs(number(row, "relative_error")) for row in balance)

    cumulative_oil, history = integrate_positive_production(wells, run["run_id"])
    producer = max((row for row in wells if row["type"] == "PRODUCER"),
                   key=lambda row: number(row, "time"))
    injectors = [row for row in wells if row["type"] == "INJECTOR"]
    maximum_injector_bhp_pa = max(number(row, "bhp") for row in injectors) * 1.0e5
    final_diagnostics = max(diagnostics, key=lambda row: int(row["step"]))

    rate = number(metadata, "total_reservoir_rate_m3_s")
    pore_volume = number(metadata, "pore_volume_m3")
    target_pvi = number(metadata, "target_pvi")
    target_day = target_pvi * pore_volume / rate / SECONDS_PER_DAY
    final_day = number(simulation, "final_time_day")
    reached_target = final_day >= target_day * (1.0 - 1.0e-10)
    valid = (reached_target and maximum_balance_error <= 1.0e-6 and
             maximum_injector_bhp_pa <= 30.0e6 * (1.0 + 1.0e-10))

    summary: dict[str, object] = {
        "run_order": int(run["run_order"]),
        "run_id": run["run_id"],
        "co2_feed_mole_fraction": float(run["x_co2_feed"]),
        "scw_feed_mole_fraction": float(run["x_scw_feed"]),
        "final_pvi": final_day * SECONDS_PER_DAY * rate / pore_volume,
        "cumulative_hydrocarbon_produced_kg": produced_hc,
        "hydrocarbon_recovery_fraction": hydrocarbon_recovery,
        "cumulative_heavy_produced_kg": produced_heavy,
        "heavy_recovery_fraction": heavy_recovery,
        "cumulative_oil_reservoir_m3": cumulative_oil,
        "final_oil_rate_reservoir_m3_day": max(
            0.0, -number(producer, "q_oil_reservoir")
        ) * SECONDS_PER_DAY,
        "final_total_rate_reservoir_m3_day": max(
            0.0, -number(producer, "q_total_reservoir")
        ) * SECONDS_PER_DAY,
        "h2o_breakthrough_pvi": breakthrough_pvi(balance, "H2O", rate, pore_volume),
        "co2_breakthrough_pvi": breakthrough_pvi(balance, "CO2", rate, pore_volume),
        "final_mean_oil_saturation": number(final_diagnostics, "s_o_avg"),
        "final_mean_gas_saturation": number(final_diagnostics, "s_g_avg"),
        "final_mean_water_saturation": number(final_diagnostics, "s_w_avg"),
        "maximum_injector_bhp_Pa": maximum_injector_bhp_pa,
        "well_control_switch_count": len(switches),
        "maximum_component_relative_mass_balance_error": maximum_balance_error,
        "accepted_internal_steps": int(simulation["accepted_internal_steps"]),
        "rejected_internal_steps": int(simulation["rejected_internal_steps"]),
        "interpretation_gate": "PASS" if valid else "FAIL",
    }
    return summary, history


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: object, digits: int = 4) -> str:
    if value is None:
        return "—"
    number_value = float(value)
    if not math.isfinite(number_value):
        return "—"
    return f"{number_value:.{digits}g}"


def write_report(root: Path, summaries: list[dict[str, object]]) -> None:
    by_scw = sorted(summaries, key=lambda row: float(row["scw_feed_mole_fraction"]))
    baseline = by_scw[0]
    lines = [
        "# CO2/SCW 五点驱替结果比较",
        "",
        "下表只在每组 `interpretation_gate` 均为 PASS 时用于物理解读。",
        "",
        "| CO2/SCW (mol%) | 累计产烃 (kg) | 总烃采收率 (%) | 累计油相体积 (m3) | 最终油率 (m3/d) | 重质采收率 (%) | H2O突破 (PVI) | CO2突破 (PVI) | 门槛 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in by_scw:
        lines.append(
            f"| {100*float(row['co2_feed_mole_fraction']):.0f}/{100*float(row['scw_feed_mole_fraction']):.0f} "
            f"| {fmt(row['cumulative_hydrocarbon_produced_kg'], 6)} "
            f"| {100*float(row['hydrocarbon_recovery_fraction']):.4f} "
            f"| {fmt(row['cumulative_oil_reservoir_m3'], 6)} "
            f"| {fmt(row['final_oil_rate_reservoir_m3_day'], 6)} "
            f"| {100*float(row['heavy_recovery_fraction']):.4f} "
            f"| {fmt(row['h2o_breakthrough_pvi'])} "
            f"| {fmt(row['co2_breakthrough_pvi'])} "
            f"| {row['interpretation_gate']} |"
        )
    lines.extend(["", "## 相对纯 CO2 的效应", "",
                  "| SCW (mol%) | 产烃质量差 (kg) | 采收率差 (百分点) |",
                  "|---:|---:|---:|"])
    for row in by_scw[1:]:
        lines.append(
            f"| {100*float(row['scw_feed_mole_fraction']):.0f} "
            f"| {float(row['cumulative_hydrocarbon_produced_kg']) - float(baseline['cumulative_hydrocarbon_produced_kg']):.6g} "
            f"| {100*(float(row['hydrocarbon_recovery_fraction']) - float(baseline['hydrocarbon_recovery_fraction'])):.4f} |"
        )
    all_pass = all(row["interpretation_gate"] == "PASS" for row in summaries)
    lines.extend(["", f"总体解释门槛：**{'PASS' if all_pass else 'FAIL'}**。",
                  "" if all_pass else "至少一组未达到统一终点、质量闭合或井压门槛，因此暂不作配比优劣的物理解读。"])
    (root / "RESULTS_COMPARISON.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_svg(root: Path, summaries: list[dict[str, object]]) -> None:
    rows = sorted(summaries, key=lambda row: float(row["scw_feed_mole_fraction"]))
    width, height = 980, 430
    panel_width = 390
    x_starts = [70, 560]
    metrics = [
        ("cumulative_oil_reservoir_m3", "Cumulative oil, reservoir m3"),
        ("hydrocarbon_recovery_fraction", "Hydrocarbon recovery, %"),
    ]
    colors = ["#2F5597", "#4472C4", "#70AD47", "#ED7D31", "#C00000"]
    content = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
               '<rect width="100%" height="100%" fill="white"/>',
               '<style>text{font-family:Arial,sans-serif;fill:#222}.axis{stroke:#333;stroke-width:1}.grid{stroke:#ddd;stroke-width:1}</style>']
    for panel, ((key, title), x0) in enumerate(zip(metrics, x_starts)):
        values = [float(row[key]) * (100.0 if panel == 1 else 1.0) for row in rows]
        maximum = max(values) if max(values) > 0.0 else 1.0
        content.append(f'<text x="{x0 + panel_width/2}" y="28" text-anchor="middle" font-size="15" font-weight="bold">{html.escape(title)}</text>')
        content.append(f'<line class="axis" x1="{x0}" y1="350" x2="{x0+panel_width}" y2="350"/>')
        content.append(f'<line class="axis" x1="{x0}" y1="55" x2="{x0}" y2="350"/>')
        for tick in range(5):
            y = 350 - tick * 295 / 4
            tick_value = maximum * tick / 4
            content.append(f'<line class="grid" x1="{x0}" y1="{y:.1f}" x2="{x0+panel_width}" y2="{y:.1f}"/>')
            content.append(f'<text x="{x0-8}" y="{y+4:.1f}" text-anchor="end" font-size="10">{tick_value:.3g}</text>')
        bar_width = 52
        gap = (panel_width - bar_width * len(rows)) / (len(rows) + 1)
        for index, (row, value) in enumerate(zip(rows, values)):
            x = x0 + gap + index * (bar_width + gap)
            bar_height = 295 * value / maximum
            y = 350 - bar_height
            content.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width}" height="{bar_height:.1f}" fill="{colors[index]}"/>')
            label = f"{100*float(row['co2_feed_mole_fraction']):.0f}/{100*float(row['scw_feed_mole_fraction']):.0f}"
            content.append(f'<text x="{x+bar_width/2:.1f}" y="369" text-anchor="middle" font-size="10">{label}</text>')
            content.append(f'<text x="{x+bar_width/2:.1f}" y="{max(48,y-5):.1f}" text-anchor="middle" font-size="9">{value:.3g}</text>')
        content.append(f'<text x="{x0+panel_width/2}" y="397" text-anchor="middle" font-size="11">CO2/SCW feed (mol%)</text>')
    content.append('</svg>')
    (root / "final_oil_comparison.svg").write_text("\n".join(content), encoding="utf-8")


def main() -> int:
    result_root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else (CASE_ROOT / "results" / "formal_matrix")
    manifest = read_rows(CASE_ROOT / "experiment_manifest.csv")
    plan = json.loads((CASE_ROOT / "RUN_PLAN.json").read_text(encoding="utf-8"))
    summaries: list[dict[str, object]] = []
    histories: list[dict[str, object]] = []
    errors: list[str] = []
    for run in manifest:
        run_dir = result_root / f"{int(run['run_order']):02d}_{run['run_id']}"
        try:
            summary, history = summarize_run(run_dir, run)
            summaries.append(summary)
            histories.extend(history)
        except (FileNotFoundError, KeyError, ValueError) as error:
            errors.append(str(error))

    result_root.mkdir(parents=True, exist_ok=True)
    if errors:
        message = "Incomplete formal matrix:\n" + "\n".join(f"- {item}" for item in errors) + "\n"
        (result_root / "RESULTS_COMPARISON.md").write_text(message, encoding="utf-8")
        print(message, file=sys.stderr)
        return 2
    if len(summaries) != len(plan["run_order"]):
        raise RuntimeError("summary count differs from the registered matrix")

    write_csv(result_root / "comparison_summary.csv", summaries)
    write_csv(result_root / "production_history.csv", histories)
    write_report(result_root, summaries)
    write_svg(result_root, summaries)
    print(f"Wrote five-run comparison to {result_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
