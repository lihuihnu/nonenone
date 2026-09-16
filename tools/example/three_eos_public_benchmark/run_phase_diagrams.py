#!/usr/bin/env python3
"""Reconstruct public phase diagrams with the production PR, SW, and CPA paths."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import pathlib
import subprocess
import sys
from dataclasses import dataclass
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager
from matplotlib.lines import Line2D

from run_benchmark import load_reference_data, read_csv_by_id, run_mpmc_executable


PR_TEMPERATURE_K = 243.60
PR_ISOTHERM_TOLERANCE_K = 0.011
PR_COMPOSITION_MIN = 0.0
PR_COMPOSITION_MAX = 1.0
PR_CURVE_POINTS = 241
PR_FIGURE_X_LIMITS = (0.0, 0.76)
PR_FIGURE_PRESSURE_LIMITS_MPA = (1.0, 7.10)
CPA_TEMPERATURE_K = 333.15
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
LINESTYLES = ["-", "--", "-.", ":"]
MARKERS = ["o", "s", "^", "D"]
NEUTRAL = "#4B5563"


def register_chinese_font() -> str:
    candidates = [
        pathlib.Path("/mnt/c/Windows/Fonts/simhei.ttf"),
        pathlib.Path("C:/Windows/Fonts/simhei.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            font_manager.fontManager.addfont(candidate)
            return font_manager.FontProperties(fname=candidate).get_name()
    raise RuntimeError("A Chinese font is required for bilingual figure export")


CHINESE_FONT_FAMILY = register_chinese_font()

matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.9,
        "axes.labelsize": 9.5,
        "axes.titlesize": 8.8,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.0,
        "lines.solid_capstyle": "round",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.unicode_minus": False,
    }
)


@dataclass(frozen=True)
class RequestMeta:
    family: str
    kind: str
    coordinate: float
    observed_pressure_Pa: float | None = None
    observed_secondary: float | None = None
    observed_uncertainty: float | None = None
    observed_coordinate_uncertainty: float | None = None
    observed_secondary_uncertainty: float | None = None
    salinity_molal: float = 0.0
    temperature_K: float = 0.0


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def unique_sorted(values: Iterable[float]) -> list[float]:
    return sorted({round(float(value), 12) for value in values})


def interpolate_reference(
    coordinate: float,
    coordinates: list[float],
    pressures: list[float],
) -> float:
    order = np.argsort(np.asarray(coordinates))
    x = np.asarray(coordinates)[order]
    p = np.asarray(pressures)[order]
    return float(np.interp(coordinate, x, p))


def load_kurihara(path: pathlib.Path) -> tuple[dict[str, Any], list[dict[str, float]]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    points = [
        {
            "temperature_K": float(document["temperature_K"]),
            "pressure_Pa": float(point["pressure_kPa"]) * 1000.0,
            "x_methanol": float(point["x_methanol"]),
            "y_methanol": float(point["y_methanol"]),
        }
        for point in document["points"]
    ]
    return document, points


def add_request(
    rows: list[dict[str, Any]],
    metadata: dict[str, RequestMeta],
    request_id: str,
    model: str,
    operation: str,
    temperature_K: float,
    pressure_Pa: float,
    salinity_molal: float,
    z0: float,
    z1: float,
    meta: RequestMeta,
) -> None:
    rows.append(
        {
            "id": request_id,
            "model": model,
            "operation": operation,
            "temperature_K": temperature_K,
            "pressure_Pa": pressure_Pa,
            "salinity_molal": salinity_molal,
            "z0": z0,
            "z1": z1,
        }
    )
    metadata[request_id] = meta


def build_requests(
    pr_points: list[dict[str, Any]],
    sw_points: list[dict[str, Any]],
    cpa_points: list[dict[str, float]],
) -> tuple[list[dict[str, Any]], dict[str, RequestMeta]]:
    rows: list[dict[str, Any]] = []
    metadata: dict[str, RequestMeta] = {}

    pr_isotherm = [
        point
        for point in pr_points
        if abs(point["temperature_K"] - PR_TEMPERATURE_K)
        <= PR_ISOTHERM_TOLERANCE_K
    ]
    pr_x = [point["x_methane"] for point in pr_isotherm]
    pr_p = [point["pressure_Pa"] for point in pr_isotherm]
    # Scan the complete binary composition axis.  The benchmark executable
    # handles subcritical pure endpoints by direct liquid/vapour fugacity
    # equality.  Supercritical endpoints remain explicit missing boundaries.
    pr_grid = unique_sorted(
        list(np.linspace(PR_COMPOSITION_MIN, PR_COMPOSITION_MAX, PR_CURVE_POINTS))
        + pr_x
    )
    pr_observed = {round(point["x_methane"], 12): point for point in pr_isotherm}
    for index, x_methane in enumerate(pr_grid):
        observed = pr_observed.get(round(x_methane, 12))
        reference_pressure = (
            observed["pressure_Pa"]
            if observed is not None
            else interpolate_reference(x_methane, pr_x, pr_p)
        )
        add_request(
            rows,
            metadata,
            f"phase_pr_{index:03d}",
            "pr",
            "bubble_pressure",
            PR_TEMPERATURE_K,
            reference_pressure,
            0.0,
            x_methane,
            1.0 - x_methane,
            RequestMeta(
                family="PR",
                kind="experiment" if observed is not None else "curve",
                coordinate=x_methane,
                observed_pressure_Pa=(
                    None if observed is None else observed["pressure_Pa"]
                ),
                observed_secondary=(
                    None if observed is None else 1.0 - observed["y_ethane"]
                ),
                observed_uncertainty=(
                    None
                    if observed is None
                    else observed["pressure_measurement_u_Pa"]
                ),
                observed_coordinate_uncertainty=(
                    None if observed is None else observed["x_methane_uc"]
                ),
                observed_secondary_uncertainty=(
                    None if observed is None else observed["y_ethane_uc"]
                ),
                temperature_K=PR_TEMPERATURE_K,
            ),
        )

    sw_groups: dict[tuple[float, float], list[dict[str, Any]]] = {}
    for point in sw_points:
        key = (point["temperature_K"], point["salinity_molal"])
        sw_groups.setdefault(key, []).append(point)
    for group_index, ((temperature, salinity), points) in enumerate(
        sorted(sw_groups.items())
    ):
        pressures = [point["pressure_Pa"] for point in points]
        pressure_grid = unique_sorted(
            list(np.linspace(min(pressures), max(pressures), 81)) + pressures
        )
        observed = {round(point["pressure_Pa"], 6): point for point in points}
        for point_index, pressure in enumerate(pressure_grid):
            reference = observed.get(round(pressure, 6))
            add_request(
                rows,
                metadata,
                f"phase_sw_{group_index:02d}_{point_index:03d}",
                "sw",
                "sw_solubility",
                temperature,
                pressure,
                salinity,
                0.5,
                0.5,
                RequestMeta(
                    family="SW",
                    kind="experiment" if reference is not None else "curve",
                    coordinate=pressure,
                    observed_secondary=(
                        None if reference is None else reference["co2_molality"]
                    ),
                    observed_uncertainty=(
                        None if reference is None else reference["co2_molality_u"]
                    ),
                    salinity_molal=salinity,
                    temperature_K=temperature,
                ),
            )

    cpa_interior = [
        point for point in cpa_points if 0.0 < point["x_methanol"] < 1.0
    ]
    cpa_x = [point["x_methanol"] for point in cpa_interior]
    cpa_p = [point["pressure_Pa"] for point in cpa_interior]
    cpa_grid = unique_sorted(list(np.linspace(0.02, 0.965, 121)) + cpa_x)
    cpa_observed = {round(point["x_methanol"], 12): point for point in cpa_interior}
    for index, x_methanol in enumerate(cpa_grid):
        observed = cpa_observed.get(round(x_methanol, 12))
        reference_pressure = (
            observed["pressure_Pa"]
            if observed is not None
            else interpolate_reference(x_methanol, cpa_x, cpa_p)
        )
        add_request(
            rows,
            metadata,
            f"phase_cpa_{index:03d}",
            "cpa",
            "bubble_pressure",
            CPA_TEMPERATURE_K,
            reference_pressure,
            0.0,
            1.0 - x_methanol,
            x_methanol,
            RequestMeta(
                family="CPA",
                kind="experiment" if observed is not None else "curve",
                coordinate=x_methanol,
                observed_pressure_Pa=(
                    None if observed is None else observed["pressure_Pa"]
                ),
                observed_secondary=(
                    None if observed is None else observed["y_methanol"]
                ),
                temperature_K=CPA_TEMPERATURE_K,
            ),
        )

    return rows, metadata


def write_request_csv(path: pathlib.Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "id",
        "model",
        "operation",
        "temperature_K",
        "pressure_Pa",
        "salinity_molal",
        "z0",
        "z1",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def collect_rows(
    metadata: dict[str, RequestMeta],
    outputs: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for request_id, meta in metadata.items():
        output = outputs[request_id]
        converged = bool(output.get("converged", 0.0)) and finite(output.get("value"))
        predicted_secondary = math.nan
        predicted_value = math.nan
        if converged:
            predicted_value = float(output["value"])
            if meta.family == "PR":
                predicted_secondary = float(output["yg0"])
            elif meta.family == "CPA":
                predicted_secondary = float(output["yg1"])
            elif meta.family == "SW":
                predicted_secondary = predicted_value
        rows.append(
            {
                "id": request_id,
                "family": meta.family,
                "kind": meta.kind,
                "temperature_K": meta.temperature_K,
                "salinity_molal": meta.salinity_molal,
                "coordinate": meta.coordinate,
                "converged": int(converged),
                "predicted_value": predicted_value,
                "predicted_secondary": predicted_secondary,
                "observed_pressure_Pa": meta.observed_pressure_Pa,
                "observed_secondary": meta.observed_secondary,
                "observed_uncertainty": meta.observed_uncertainty,
                "observed_coordinate_uncertainty": meta.observed_coordinate_uncertainty,
                "observed_secondary_uncertainty": meta.observed_secondary_uncertainty,
            }
        )
    return rows


def write_rows(path: pathlib.Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_detail_rows(path: pathlib.Path) -> list[dict[str, Any]]:
    numeric_fields = {
        "temperature_K",
        "salinity_molal",
        "coordinate",
        "predicted_value",
        "predicted_secondary",
        "observed_pressure_Pa",
        "observed_secondary",
        "observed_uncertainty",
        "observed_coordinate_uncertainty",
        "observed_secondary_uncertainty",
    }
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row: dict[str, Any] = dict(raw)
            row["converged"] = int(raw["converged"])
            for field in numeric_fields:
                row[field] = None if raw[field] == "" else float(raw[field])
            rows.append(row)
    if not rows:
        raise RuntimeError(f"No phase-diagram rows found in {path}")
    return rows


def aard(predicted: list[float], observed: list[float]) -> float:
    pairs = [
        (p, o)
        for p, o in zip(predicted, observed)
        if finite(p) and finite(o) and abs(o) > 0.0
    ]
    return 100.0 * float(np.mean([abs((p - o) / o) for p, o in pairs]))


def mae(predicted: list[float], observed: list[float]) -> float:
    pairs = [
        (p, o)
        for p, o in zip(predicted, observed)
        if finite(p) and finite(o)
    ]
    return float(np.mean([abs(p - o) for p, o in pairs]))


def mean_relative_deviation(
    predicted: list[float], observed: list[float],
) -> float:
    pairs = [
        (p, o)
        for p, o in zip(predicted, observed)
        if finite(p) and finite(o) and abs(o) > 0.0
    ]
    return 100.0 * float(np.mean([(p - o) / o for p, o in pairs]))


def mean_deviation(predicted: list[float], observed: list[float]) -> float:
    pairs = [
        (p, o)
        for p, o in zip(predicted, observed)
        if finite(p) and finite(o)
    ]
    return float(np.mean([p - o for p, o in pairs]))


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for family in ("PR", "CPA"):
        family_rows = [row for row in rows if row["family"] == family]
        experimental = [
            row
            for row in rows
            if row["family"] == family and row["kind"] == "experiment"
        ]
        pressure_predicted = [row["predicted_value"] for row in experimental]
        pressure_observed = [row["observed_pressure_Pa"] for row in experimental]
        secondary_predicted = [
            row["predicted_secondary"] for row in experimental
        ]
        secondary_observed = [
            row["observed_secondary"] for row in experimental
        ]
        family_summary = {
            "requested_points": len(experimental),
            "converged_points": sum(row["converged"] for row in experimental),
            "calculation_points": len(family_rows),
            "calculation_converged_points": sum(
                row["converged"] for row in family_rows
            ),
            "coordinate_range": [
                min(row["coordinate"] for row in experimental),
                max(row["coordinate"] for row in experimental),
            ],
            "observed_pressure_range_Pa": [
                min(pressure_observed), max(pressure_observed),
            ],
            "observed_vapor_composition_range": [
                min(secondary_observed), max(secondary_observed),
            ],
            "pressure_AARD_percent": aard(pressure_predicted, pressure_observed),
            "pressure_mean_relative_deviation_percent":
                mean_relative_deviation(pressure_predicted, pressure_observed),
            "vapor_composition_MAE": mae(
                secondary_predicted, secondary_observed,
            ),
            "vapor_composition_mean_deviation": mean_deviation(
                secondary_predicted, secondary_observed,
            ),
        }
        if family == "PR":
            converged_curve = [row for row in family_rows if row["converged"]]
            family_summary["calculated_liquid_composition_range"] = [
                min(row["coordinate"] for row in converged_curve),
                max(row["coordinate"] for row in converged_curve),
            ]
            family_summary["calculated_vapor_composition_range"] = [
                min(row["predicted_secondary"] for row in converged_curve),
                max(row["predicted_secondary"] for row in converged_curve),
            ]
            family_summary["calculated_pressure_range_Pa"] = [
                min(row["predicted_value"] for row in converged_curve),
                max(row["predicted_value"] for row in converged_curve),
            ]
            uncertainties = [
                row["observed_uncertainty"] for row in experimental
                if finite(row["observed_uncertainty"])
            ]
            liquid_uncertainties = [
                row["observed_coordinate_uncertainty"] for row in experimental
                if finite(row["observed_coordinate_uncertainty"])
            ]
            vapor_uncertainties = [
                row["observed_secondary_uncertainty"] for row in experimental
                if finite(row["observed_secondary_uncertainty"])
            ]
            family_summary["pressure_standard_uncertainty_Pa"] = uncertainties[0]
            family_summary["liquid_composition_uc_range"] = [
                min(liquid_uncertainties), max(liquid_uncertainties),
            ]
            family_summary["vapor_composition_uc_range"] = [
                min(vapor_uncertainties), max(vapor_uncertainties),
            ]
        summary[family] = family_summary
    sw_rows = [row for row in rows if row["family"] == "SW"]
    sw_experimental = [
        row
        for row in rows
        if row["family"] == "SW" and row["kind"] == "experiment"
    ]
    sw_predicted = [row["predicted_secondary"] for row in sw_experimental]
    sw_observed = [row["observed_secondary"] for row in sw_experimental]
    sw_groups: list[dict[str, Any]] = []
    for temperature, salinity in sorted({
        (row["temperature_K"], row["salinity_molal"])
        for row in sw_experimental
    }):
        group = [
            row for row in sw_experimental
            if row["temperature_K"] == temperature
            and row["salinity_molal"] == salinity
        ]
        predicted = [row["predicted_secondary"] for row in group]
        observed = [row["observed_secondary"] for row in group]
        sw_groups.append({
            "temperature_K": temperature,
            "salinity_molal": salinity,
            "points": len(group),
            "pressure_range_Pa": [
                min(row["coordinate"] for row in group),
                max(row["coordinate"] for row in group),
            ],
            "observed_solubility_range_mol_kg": [min(observed), max(observed)],
            "AARD_percent": aard(predicted, observed),
            "mean_relative_deviation_percent":
                mean_relative_deviation(predicted, observed),
        })
    summary["SW"] = {
        "requested_points": len(sw_experimental),
        "converged_points": sum(row["converged"] for row in sw_experimental),
        "calculation_points": len(sw_rows),
        "calculation_converged_points": sum(
            row["converged"] for row in sw_rows
        ),
        "temperature_range_K": [
            min(row["temperature_K"] for row in sw_experimental),
            max(row["temperature_K"] for row in sw_experimental),
        ],
        "pressure_range_Pa": [
            min(row["coordinate"] for row in sw_experimental),
            max(row["coordinate"] for row in sw_experimental),
        ],
        "salinities_molal": sorted({
            row["salinity_molal"] for row in sw_experimental
        }),
        "solubility_AARD_percent": aard(sw_predicted, sw_observed),
        "solubility_MAE_mol_kg": mae(sw_predicted, sw_observed),
        "groups": sw_groups,
    }
    summary["curve_failures"] = sum(not row["converged"] for row in rows)
    summary["curve_requests"] = len(rows)
    return summary


def summarize_envelopes(
    envelopes: dict[str, list[dict[str, float]]],
) -> dict[str, float]:
    def maximum_difference(model_a: str, model_b: str, boundary: str) -> float:
        flag = f"has_{boundary}"
        value = f"{boundary}_pressure_bar"
        a = {
            row["temperature_K"]: row[value]
            for row in envelopes[model_a]
            if row.get(flag, 0.0) > 0.5
        }
        b = {
            row["temperature_K"]: row[value]
            for row in envelopes[model_b]
            if row.get(flag, 0.0) > 0.5
        }
        common = sorted(set(a) & set(b))
        return max(abs(a[temperature] - b[temperature]) for temperature in common)

    return {
        "PR_SW_max_bubble_difference_bar": maximum_difference("PR", "SW", "bubble"),
        "PR_SW_max_dew_difference_bar": maximum_difference("PR", "SW", "dew"),
        "PR_CPA_max_bubble_difference_bar": maximum_difference("PR", "CPA", "bubble"),
        "PR_CPA_max_dew_difference_bar": maximum_difference("PR", "CPA", "dew"),
    }


def style_axes(ax: plt.Axes) -> None:
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#111827")
        spine.set_linewidth(0.85)
    ax.tick_params(
        direction="in",
        top=True,
        right=True,
        length=4,
        width=0.8,
    )


def language_style(language: str) -> dict[str, Any]:
    if language == "zh":
        return {"font.family": CHINESE_FONT_FAMILY}
    return {"font.family": "DejaVu Sans"}


def localized(language: str, english: str, chinese: str) -> str:
    return chinese if language == "zh" else english


def localized_stem(stem: pathlib.Path, language: str) -> pathlib.Path:
    if language == "zh":
        return stem.with_name(f"{stem.name}_zh")
    return stem


def use_english_legend_font(legend: Any) -> None:
    """Keep bilingual figures' legends in a consistent Latin typeface."""
    legend.get_title().set_fontfamily("DejaVu Sans")
    for label in legend.get_texts():
        label.set_fontfamily("DejaVu Sans")


def save_figure(fig: plt.Figure, stem: pathlib.Path) -> list[pathlib.Path]:
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=300, facecolor="white", transparent=False)
    try:
        fig.savefig(pdf, facecolor="white", transparent=False)
    except PermissionError:
        pdf = stem.with_name(f"{stem.name}_latest").with_suffix(".pdf")
        fig.savefig(pdf, facecolor="white", transparent=False)
        print(f"[phase-diagram] locked PDF fallback: {pdf}")
    plt.close(fig)
    return [png, pdf]


def make_pr_figure(
    rows: list[dict[str, Any]],
    output_dir: pathlib.Path,
    language: str,
) -> pathlib.Path:
    selected = sorted(
        (row for row in rows if row["family"] == "PR"),
        key=lambda row: row["coordinate"],
    )
    experimental = [row for row in selected if row["kind"] == "experiment"]
    x = np.asarray([row["coordinate"] for row in selected])
    p = np.asarray(
        [row["predicted_value"] / 1.0e6 if row["converged"] else np.nan for row in selected]
    )
    y = np.asarray(
        [row["predicted_secondary"] if row["converged"] else np.nan for row in selected]
    )
    with plt.rc_context({**language_style(language), "font.size": 9.0}):
        fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
        observed_x = [row["coordinate"] for row in experimental]
        observed_y = [row["observed_secondary"] for row in experimental]
        observed_pressure = [
            row["observed_pressure_Pa"] / 1.0e6 for row in experimental
        ]
        finite_mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(p)
        finite_indices = np.flatnonzero(finite_mask)
        finite_segments = np.split(
            finite_indices,
            np.flatnonzero(np.diff(finite_indices) > 1) + 1,
        )
        two_phase = None
        for segment in finite_segments:
            if len(segment) < 2:
                continue
            order = segment[np.argsort(p[segment])]
            patch = ax.fill_betweenx(
                p[order],
                x[order],
                y[order],
                color="#DCEAF4",
                alpha=0.62,
                linewidth=0.0,
                zorder=1,
            )
            if two_phase is None:
                two_phase = patch
        liquid = ax.plot(
            observed_x,
            observed_pressure,
            linestyle="none",
            marker="o",
            markersize=6.2,
            markerfacecolor="white",
            markeredgecolor="#111111",
            markeredgewidth=1.1,
            zorder=4,
        )[0]
        bubble = ax.plot(
            x,
            p,
            color=COLORS[0],
            linewidth=2.1,
            zorder=2.5,
        )[0]
        dew = ax.plot(
            y,
            p,
            color=COLORS[1],
            linewidth=2.1,
            zorder=2.5,
        )[0]
        vapor = ax.plot(
            observed_y,
            observed_pressure,
            linestyle="none",
            marker="^",
            markersize=6.8,
            markerfacecolor="white",
            markeredgecolor="#111111",
            markeredgewidth=1.1,
            zorder=4,
        )[0]
        ax.set(
            xlabel=localized(
                language,
                "Methane mole fraction",
                "甲烷摩尔分数",
            ),
            ylabel=localized(
                language,
                "Equilibrium pressure, $P$ [MPa]",
                "平衡压力 $P$ [MPa]",
            ),
            xlim=PR_FIGURE_X_LIMITS,
            ylim=PR_FIGURE_PRESSURE_LIMITS_MPA,
        )
        ax.set_xticks(np.arange(0.0, 0.71, 0.1))
        ax.set_yticks(np.arange(1.0, 7.01, 1.0))
        style_axes(ax)
        labels = [
            "Our PR: liquid boundary",
            "Our PR: vapor boundary",
            "Our PR: two-phase region",
            "May et al. (2015): liquid data",
            "May et al. (2015): vapor data",
        ]
        legend = ax.legend(
            [bubble, dew, two_phase, liquid, vapor],
            labels,
            frameon=True,
            facecolor="white",
            edgecolor="#94A3B8",
            framealpha=0.96,
            ncol=1,
            loc="lower right",
            columnspacing=1.35,
            handlelength=2.0,
            fontsize=7.4,
        )
        legend.get_frame().set_linewidth(0.7)
        use_english_legend_font(legend)
        stem = localized_stem(
            output_dir / "01_pr_methane_ethane_pxy", language
        )
        save_figure(fig, stem)
    return stem.with_suffix(".png")


def make_sw_figures(
    rows: list[dict[str, Any]],
    output_dir: pathlib.Path,
    language: str,
) -> list[pathlib.Path]:
    salt_styles = {
        0.0: (COLORS[0], "-", "o"),
        1.0: (COLORS[1], "--", "s"),
        3.0: (COLORS[2], "-.", "^"),
        6.0: (COLORS[3], ":", "D"),
    }

    def series_for(
        selected: list[dict[str, Any]], salinity: float
    ) -> list[dict[str, Any]]:
        return sorted(
            (row for row in selected if row["salinity_molal"] == salinity),
            key=lambda row: row["coordinate"],
        )

    def experimental_for(series: list[dict[str, Any]]) -> list[dict[str, Any]]:
        return [row for row in series if row["kind"] == "experiment"]

    def solubility_aard(series: list[dict[str, Any]]) -> float:
        experimental = experimental_for(series)
        return aard(
            [row["predicted_secondary"] for row in experimental],
            [row["observed_secondary"] for row in experimental],
        )

    def plot_series(
        ax: plt.Axes,
        series: list[dict[str, Any]],
        salinity: float,
    ) -> tuple[Any, Any]:
        color, linestyle, marker = salt_styles[salinity]
        experimental = experimental_for(series)
        model_line = ax.plot(
            [row["coordinate"] / 1.0e6 for row in series],
            [
                row["predicted_secondary"] if row["converged"] else np.nan
                for row in series
            ],
            color=color,
            linestyle=linestyle,
            linewidth=2.0,
        )[0]
        experiment_points = ax.errorbar(
            [row["coordinate"] / 1.0e6 for row in experimental],
            [row["observed_secondary"] for row in experimental],
            yerr=[row["observed_uncertainty"] or 0.0 for row in experimental],
            fmt=marker,
            markersize=5.0,
            markerfacecolor="white",
            markeredgecolor=color,
            markeredgewidth=1.3,
            ecolor=color,
            elinewidth=0.8,
            capsize=2.0,
            linestyle="none",
            zorder=4,
        )
        return model_line, experiment_points

    paths: list[pathlib.Path] = []
    temperatures = sorted({row["temperature_K"] for row in rows if row["family"] == "SW"})
    for figure_index, temperature in enumerate(temperatures, start=2):
        selected = [
            row
            for row in rows
            if row["family"] == "SW" and row["temperature_K"] == temperature
        ]
        salinities = sorted({row["salinity_molal"] for row in selected})
        with plt.rc_context({**language_style(language), "font.size": 9.0}):
            if abs(temperature - 323.15) < 0.01:
                fig, ax = plt.subplots(
                    figsize=(6.4, 4.2), layout="constrained"
                )
                legend_handles: list[Line2D] = []
                legend_labels: list[str] = []
                for salinity in salinities:
                    series = series_for(selected, salinity)
                    plot_series(ax, series, salinity)
                    color, linestyle, marker = salt_styles[salinity]
                    legend_handles.append(
                        Line2D(
                            [0],
                            [0],
                            color=color,
                            linestyle=linestyle,
                            linewidth=2.0,
                            marker=marker,
                            markerfacecolor="white",
                            markeredgecolor=color,
                            markersize=5.0,
                        )
                    )
                    legend_labels.append(f"{salinity:g} mol/kg NaCl")
                ax.set(
                    xlabel=localized(language, "Pressure [MPa]", "压力 [MPa]"),
                    ylabel=localized(
                        language,
                        "Dissolved CO$_2$ [mol/kg H$_2$O]",
                        "溶解 CO$_2$ [mol/kg H$_2$O]",
                    ),
                    xlim=(4.3, 20.8),
                )
                style_axes(ax)
                legend = ax.legend(
                    legend_handles,
                    legend_labels,
                    ncol=2,
                    frameon=True,
                    facecolor="white",
                    edgecolor="#94A3B8",
                    framealpha=0.92,
                    loc="upper left",
                    columnspacing=1.2,
                    handlelength=2.2,
                )
                legend.get_frame().set_linewidth(0.7)
                use_english_legend_font(legend)

            elif abs(temperature - 373.15) < 0.01:
                fig, parity = plt.subplots(
                    figsize=(5.2, 4.6), layout="constrained"
                )
                all_experimental: list[dict[str, Any]] = []
                legend_handles = []
                for salinity in salinities:
                    experimental = experimental_for(series_for(selected, salinity))
                    all_experimental.extend(experimental)
                    color, _, marker = salt_styles[salinity]
                    parity.errorbar(
                        [row["observed_secondary"] for row in experimental],
                        [row["predicted_secondary"] for row in experimental],
                        xerr=[row["observed_uncertainty"] or 0.0 for row in experimental],
                        fmt=marker,
                        markersize=5.2,
                        markerfacecolor="white",
                        markeredgecolor=color,
                        markeredgewidth=1.3,
                        ecolor=color,
                        elinewidth=0.8,
                        capsize=2.0,
                    )
                    legend_handles.append(
                        Line2D(
                            [0],
                            [0],
                            color=color,
                            marker=marker,
                            linestyle="none",
                            markerfacecolor="white",
                            markeredgecolor=color,
                            markersize=5.5,
                        )
                    )
                observed_values = [row["observed_secondary"] for row in all_experimental]
                predicted_values = [row["predicted_secondary"] for row in all_experimental]
                low = 0.92 * min(observed_values + predicted_values)
                high = 1.06 * max(observed_values + predicted_values)
                parity.plot([low, high], [low, high], color=NEUTRAL, linestyle="--", linewidth=1.2)
                parity.set(
                    xlabel=localized(
                        language,
                        "Experiment [mol/kg H$_2$O]",
                        "实验值 [mol/kg H$_2$O]",
                    ),
                    ylabel=localized(
                        language,
                        "Our SW [mol/kg H$_2$O]",
                        "自主程序 SW [mol/kg H$_2$O]",
                    ),
                    xlim=(low, high),
                    ylim=(low, high),
                    aspect="equal",
                )
                reference_handle = Line2D(
                    [0], [0], color=NEUTRAL, linestyle="--", linewidth=1.2
                )
                legend = parity.legend(
                    [*legend_handles, reference_handle],
                    [
                        *[f"NaCl {salinity:g} mol/kg" for salinity in salinities],
                        "1:1 reference",
                    ],
                    frameon=True,
                    facecolor="white",
                    edgecolor="#94A3B8",
                    framealpha=0.92,
                    loc="upper left",
                )
                legend.get_frame().set_linewidth(0.7)
                use_english_legend_font(legend)
                style_axes(parity)

            else:
                fig, ax = plt.subplots(
                    figsize=(5.2, 4.0), layout="constrained"
                )
                metrics = [
                    solubility_aard(series_for(selected, salinity))
                    for salinity in salinities
                ]
                positions = np.arange(len(salinities))
                colors = [salt_styles[salinity][0] for salinity in salinities]
                ax.barh(positions, metrics, color=colors, height=0.58)
                ax.set_yticks(
                    positions,
                    labels=[f"NaCl {salinity:g} mol/kg" for salinity in salinities],
                )
                ax.invert_yaxis()
                ax.set(
                    xlabel=localized(
                        language,
                        "Absolute average relative deviation, AARD [%]",
                        "平均绝对相对偏差 AARD [%]",
                    ),
                    xlim=(0.0, max(metrics) * 1.30),
                )
                for position, value in zip(positions, metrics):
                    ax.text(
                        value + 0.03 * max(metrics),
                        position,
                        f"{value:.2f}%",
                        va="center",
                        fontsize=8.2,
                    )
                style_axes(ax)
            stem = localized_stem(
                output_dir
                / f"{figure_index:02d}_sw_co2_brine_{int(round(temperature))}",
                language,
            )
            save_figure(fig, stem)
        paths.append(stem.with_suffix(".png"))
    return paths


def make_cpa_figure(
    rows: list[dict[str, Any]],
    cpa_points: list[dict[str, float]],
    output_dir: pathlib.Path,
    language: str,
) -> pathlib.Path:
    selected = sorted(
        (row for row in rows if row["family"] == "CPA"),
        key=lambda row: row["coordinate"],
    )
    x = np.asarray([row["coordinate"] for row in selected])
    p = np.asarray(
        [row["predicted_value"] / 1000.0 if row["converged"] else np.nan for row in selected]
    )
    y = np.asarray(
        [row["predicted_secondary"] if row["converged"] else np.nan for row in selected]
    )
    with plt.rc_context({**language_style(language), "font.size": 9.0}):
        fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(p)
        fill_order = np.argsort(p[finite])
        two_phase = ax.fill_betweenx(
            p[finite][fill_order],
            x[finite][fill_order],
            y[finite][fill_order],
            color="#DCEAF4",
            alpha=0.72,
            linewidth=0.0,
            zorder=1,
        )
        bubble = ax.plot(
            x,
            p,
            color=COLORS[0],
            linewidth=2.1,
            zorder=2.5,
        )[0]
        dew = ax.plot(
            y,
            p,
            color=COLORS[1],
            linewidth=2.1,
            linestyle="--",
            zorder=2.5,
        )[0]
        liquid = ax.scatter(
            [point["x_methanol"] for point in cpa_points],
            [point["pressure_Pa"] / 1000.0 for point in cpa_points],
            marker="o",
            s=36,
            facecolors="white",
            edgecolors="#111111",
            linewidths=1.1,
            zorder=4,
        )
        vapor = ax.scatter(
            [point["y_methanol"] for point in cpa_points],
            [point["pressure_Pa"] / 1000.0 for point in cpa_points],
            marker="^",
            s=38,
            facecolors="white",
            edgecolors="#111111",
            linewidths=1.1,
            zorder=4,
        )
        ax.set(
            xlabel=localized(
                language,
                "Methanol mole fraction",
                "甲醇摩尔分数",
            ),
            ylabel=localized(
                language,
                "Equilibrium pressure, $P$ [kPa]",
                "平衡压力 $P$ [kPa]",
            ),
            xlim=(-0.02, 1.02),
        )
        style_axes(ax)
        legend = ax.legend(
            [bubble, dew, two_phase, liquid, vapor],
            [
                "Our CPA: liquid boundary",
                "Our CPA: vapor boundary",
                "Our CPA: two-phase region",
                "Kurihara et al.: liquid",
                "Kurihara et al.: vapor",
            ],
            frameon=True,
            facecolor="white",
            edgecolor="#94A3B8",
            framealpha=0.92,
            ncol=2,
            loc="upper left",
            columnspacing=1.35,
            handlelength=2.2,
        )
        legend.get_frame().set_linewidth(0.7)
        use_english_legend_font(legend)
        stem = localized_stem(
            output_dir / "05_cpa_methanol_water_pxy", language
        )
        save_figure(fig, stem)
    return stem.with_suffix(".png")


def read_same_fluid_envelopes(root: pathlib.Path) -> dict[str, list[dict[str, float]]]:
    result: dict[str, list[dict[str, float]]] = {}
    for model in ("pr", "sw", "cpa"):
        path = root / model / "pt_vle_envelope.csv"
        rows: list[dict[str, float]] = []
        with path.open(newline="", encoding="utf-8") as handle:
            for raw in csv.DictReader(handle):
                rows.append(
                    {
                        key: float(value)
                        for key, value in raw.items()
                        if value is not None and value != ""
                    }
                )
        result[model.upper()] = rows
    return result


def make_same_fluid_figure(
    envelopes: dict[str, list[dict[str, float]]],
    output_dir: pathlib.Path,
    language: str,
) -> pathlib.Path:
    model_styles = {
        "PR": (COLORS[0], "o"),
        "SW": (COLORS[1], "s"),
        "CPA": (COLORS[2], "^"),
    }
    with plt.rc_context({**language_style(language), "font.size": 9.0}):
        fig, ax = plt.subplots(figsize=(6.4, 4.4), layout="constrained")
        plotted: dict[tuple[str, str], Any] = {}
        for model in ("SW", "PR", "CPA"):
            rows = envelopes[model]
            color, marker = model_styles[model]
            temperatures = [row["temperature_K"] for row in rows]
            marker_size = 6.2 if model == "PR" else 4.4
            marker_face = "white" if model == "PR" else color
            plotted[(model, "bubble")] = ax.plot(
                temperatures,
                [
                    row["bubble_pressure_bar"]
                    if row.get("has_bubble", 0.0) > 0.5
                    else np.nan
                    for row in rows
                ],
                color=color,
                linewidth=1.9,
                linestyle="-",
                marker=marker,
                markersize=marker_size,
                markerfacecolor=marker_face,
                markeredgecolor=color,
                markeredgewidth=1.0,
            )[0]
            plotted[(model, "dew")] = ax.plot(
                temperatures,
                [
                    row["dew_pressure_bar"]
                    if row.get("has_dew", 0.0) > 0.5
                    else np.nan
                    for row in rows
                ],
                color=color,
                linewidth=1.9,
                linestyle="--",
                marker=marker,
                markersize=marker_size,
                markerfacecolor=marker_face,
                markeredgecolor=color,
                markeredgewidth=1.0,
            )[0]
        ax.set(
            xlabel=localized(language, "Temperature [K]", "温度 [K]"),
            ylabel=localized(
                language,
                "Equilibrium pressure [bar, log scale]",
                "平衡压力 [bar，对数坐标]",
            ),
            yscale="log",
        )
        handles = [
            plotted[(model, boundary)]
            for model in ("PR", "SW", "CPA")
            for boundary in ("bubble", "dew")
        ]
        labels = [
            f"{model} {'bubble' if boundary == 'bubble' else 'dew'}"
            for model in ("PR", "SW", "CPA")
            for boundary in ("bubble", "dew")
        ]
        legend = ax.legend(
            handles,
            labels,
            ncol=2,
            frameon=True,
            facecolor="white",
            edgecolor="#94A3B8",
            framealpha=0.92,
            loc="lower left",
            columnspacing=1.2,
            handlelength=2.8,
        )
        legend.get_frame().set_linewidth(0.7)
        use_english_legend_font(legend)
        style_axes(ax)
        stem = localized_stem(
            output_dir / "06_three_eos_same_fluid_pt_envelope", language
        )
        save_figure(fig, stem)
    return stem.with_suffix(".png")


def write_book(
    output_dir: pathlib.Path,
    image_paths: list[pathlib.Path],
    filename: str = "phase_diagram_figures.pdf",
) -> pathlib.Path:
    path = output_dir / filename
    try:
        with PdfPages(path) as pdf:
            for page_number, image_path in enumerate(image_paths, start=1):
                image = plt.imread(image_path)
                fig, ax = plt.subplots(figsize=(11.69, 8.27), layout="constrained")
                ax.imshow(image)
                ax.axis("off")
                fig.text(
                    0.975,
                    0.022,
                    f"{page_number} / {len(image_paths)}",
                    ha="right",
                    va="bottom",
                    fontsize=8,
                    color="#455A64",
                )
                pdf.savefig(fig, facecolor="white")
                plt.close(fig)
    except PermissionError:
        fallback = output_dir / "phase_diagram_figures_optimized.pdf"
        if fallback == path:
            raise
        return write_book(output_dir, image_paths, fallback.name)
    return path


def write_report(
    path: pathlib.Path,
    summary: dict[str, Any],
    image_paths_by_language: dict[str, list[pathlib.Path]],
    kurihara: dict[str, Any],
    reference_dir: pathlib.Path,
) -> None:
    pr = summary["PR"]
    sw = summary["SW"]
    cpa = summary["CPA"]
    same = summary["same_fluid"]
    sw_groups = {
        (group["temperature_K"], group["salinity_molal"]): group
        for group in sw["groups"]
    }
    image_paths = image_paths_by_language["en"]
    chinese_paths = image_paths_by_language["zh"]
    lines = [
        "# PR、SW 与 CPA 公开相图重建及模型比较",
        "",
        "> 本报告使用当前生产 EOS/Flash 逐点重建公开相平衡图。PR、SW、CPA 分别在其代表性体系中与公开数据比较；最后一张同流体图仅展示模型敏感性。",
        "",
        "## 实验设计与读图方法",
        "",
        "本报告分析的是静态热力学平衡，不包含储层网格、井、流速或时间推进。公开点直接来自三篇文献；自主程序使用预先冻结的纯组分参数、二元参数和生产 EOS/Flash 求解，不使用本报告中的点重新拟合。前三组实验分别检验各 EOS 的代表性适用体系，最后一组才把三种 EOS 放到同一工程流体上比较。",
        "",
        "| 图组 | 固定条件 | 扫描变量 | 公开点 | 模型设置 |",
        "|---|---|---|---:|---|",
        f"| PR 甲烷–乙烷 | T = {PR_TEMPERATURE_K:.2f} K | 液相甲烷摩尔分数与泡点压力 | {pr['requested_points']} | PR78，k_CH4,C2H6 = 0 |",
        "| SW CO2–盐水 | T = 323.15/373.15/423.15 K | P = 约 5–20 MPa；NaCl = 0/1/3/6 mol·kg⁻¹ H2O | 40 | SW 水相模型 + Chabab-2019 CO2–H2O 参数化 |",
        f"| SRK-CPA 甲醇–水 | T = {CPA_TEMPERATURE_K:.2f} K | 液相甲醇摩尔分数与泡点压力 | 12（10 个内部点计入误差） | 4C 水、2B 甲醇，k12 = -0.09 |",
        "| 三 EOS 同流体 | z(H2O/CO2/CH4/nC16) = 0.10/0.25/0.50/0.15 | T = 280–800 K；P = 1–600 bar | 无 | restricted O/G，同网格、同工程参数 |",
        "",
        "P–x–y 图中的液相边界表示给定液相组成 x 时开始出现第一泡气体的压力；气相边界表示与之平衡的气相组成 y。两条边界之间的色块是该温度下的气液两相区。圆点/三角形是公开液相/气相数据，线是自主程序计算。SW 图报告的是水相 CO2 平衡质量摩尔浓度，并不是二元烃类 P–x–y 包络。",
        "",
        "AARD 是逐点绝对相对误差的平均值；MAE 是绝对误差平均值。AARD 不显示高估或低估方向，因此各节同时报告平均有符号偏差。误差棒只表示相应数据源明确给出的测量不确定度，不代表模型置信区间；其方向和置信层级按各原文定义读取。",
        "",
        "## 状态方程参数设置",
        "",
        "以下参数表记录本次数值实验实际传入生产 EOS 的数值，而不是只列模型名称。`Tc`、`Pc`、`Vc`、`ω` 和 `M` 分别表示临界温度、临界压力、临界摩尔体积、偏心因子和摩尔质量。除特别说明外，二元作用系数采用经典对称混合规则 `a_ij = sqrt(a_i a_j)(1-k_ij)`。本组相图计算未启用体积平移，因此所有体积平移常数均为 0。",
        "",
        "### PR 甲烷–乙烷公开相图参数",
        "",
        "| 组分 | Tc [K] | Pc [MPa] | Vc [cm³/mol] | ω | M [g/mol] | PR alpha 系数 m_i |",
        "|---|---:|---:|---:|---:|---:|---:|",
        "| CH4 | 190.555 | 4.598837 | 97.7530 | 0.01131 | 16.043 | 0.392048 |",
        "| C2H6 | 305.400 | 4.883900 | 148.1772 | 0.09800 | 30.070 | 0.523189 |",
        "",
        "| PR 参数 | 本次取值 |",
        "|---|---:|",
        "| Ωa | 0.4572355289213822 |",
        "| Ωb | 0.07779607390388846 |",
        "| 立方根参数 δ1 / δ2 | 2.414213562373095 / -0.414213562373095 |",
        "| k_CH4,C2H6 | 0.0000 |",
        "| 最小组分摩尔分数 | 1.0×10⁻³⁰ |",
        "",
        "这里 `alpha_i(T) = [1 + m_i(1-sqrt(T/Tc_i))]²`。甲烷和乙烷的偏心因子均小于 0.49，因此 PR76 与 PR78 的高偏心因子分段在本体系中给出相同结果。",
        "",
        "### SW CO2–盐水公开实验参数",
        "",
        "| 组分 | Tc [K] | Pc [MPa] | Vc [cm³/mol] | ω | M [g/mol] |",
        "|---|---:|---:|---:|---:|---:|",
        "| CO2 | 304.1282 | 7.3773 | 94.0 | 0.22394 | 44.0095 |",
        "| H2O | 647.0960 | 22.0640 | 56.0 | 0.34430 | 18.01528 |",
        "",
        "SW 保留 PR 立方常数 `Ωa = 0.4572355289213822`、`Ωb = 0.07779607390388846` 和 `δ1/δ2 = 1±sqrt(2)`。非水相使用常数 `k_CO2,H2O = 0.1896`；水相中的该系数由 Chabab-2019 温度—盐度关联式覆盖。盐度 `m_s` 的单位为 mol NaCl/kg H2O，`Tr = T/304.13 K`：",
        "",
        "`k_CO2,H2O^aq = Tr(a + b Tr + c Tr m_s) + m_s²(d + e Tr) + f`",
        "",
        "| Chabab-2019 系数 | a | b | c | d | e | f |",
        "|---|---:|---:|---:|---:|---:|---:|",
        "| 数值 | 0.43575155 | -0.05766906744 | 0.00826464849 | 0.00129539193 | -0.0016698848 | -0.47866096 |",
        "",
        "水组分采用 SW alpha 修正：`alpha_H2O = {1 + 0.4530[1-Tr,w(1-0.0103 m_s^1.1)] + 0.0034(Tr,w^-3-1)}²`，其中 `Tr,w=T/647.096 K`。公开实验逐组传入 `m_s = 0、1、3、6 mol/kg H2O`；没有另外设置水相体积平移。",
        "",
        "### SRK-CPA 甲醇–水公开相图参数",
        "",
        "| 组分 | Tc [K] | Pc [MPa] | Vc [cm³/mol] | ω | M [g/mol] | a0 [Pa·m⁶/mol²] | b [m³/mol] | c1 | ε [J/mol] | β | 供体/受体位点 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        "| H2O（4C） | 647.30 | 22.0483 | 56.0 | 0.3443 | 18.01528 | 0.12277 | 1.4515×10⁻⁵ | 0.67359 | 16655 | 0.06920 | 2 / 2 |",
        "| CH3OH（2B） | 512.60 | 8.0959 | 118.0 | 0.5650 | 32.04186 | 0.40531 | 3.0978×10⁻⁵ | 0.43102 | 24591 | 0.01610 | 1 / 1 |",
        "",
        "SRK 立方项采用 `Ωa=0.42748`、`Ωb=0.08664`、`δ1=1`、`δ2=0`；`a_i(T)=a0_i[1+c1_i(1-sqrt(T/Tc_i))]²`。甲醇–水立方项 `k12=-0.09`。径向分布函数使用 simplified 形式；未显式设置交叉缔合矩阵时，CR-1 规则取 `ε_ij=(ε_i+ε_j)/2`、`β_ij=sqrt(β_iβ_j)`。",
        "",
        "### 三 EOS 同流体 P–T 比较参数",
        "",
        "同流体图使用另一套四组分工程参数。下表是 PR 和 SW 共用的纯组分元数据，也是 CPA 非缔合组分构造 SRK 立方参数的输入。",
        "",
        "| 组分 | z | Tc [K] | Pc [MPa] | Vc [cm³/mol] | ω | M [g/mol] | PR m_i |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        "| H2O | 0.10 | 647.096 | 22.0640 | 56.0 | 0.34430 | 18.01528 | 0.873643 |",
        "| CO2 | 0.25 | 304.1282 | 7.3773 | 94.0 | 0.22394 | 44.0095 | 0.706477 |",
        "| CH4 | 0.50 | 190.564 | 4.5992 | 99.0 | 0.01142 | 16.043 | 0.392217 |",
        "| nC16 | 0.15 | 723.0 | 1.4100 | 900.0 | 0.7420 | 226.441 | 1.397817 |",
        "",
        "| 对称 BIP `k_ij` | H2O | CO2 | CH4 | nC16 |",
        "|---|---:|---:|---:|---:|",
        "| H2O | 0 | 0.1896 | 0.4850 | 0.5000 |",
        "| CO2 | 0.1896 | 0 | 0.1200 | 0.0900 |",
        "| CH4 | 0.4850 | 0.1200 | 0 | 0.0000 |",
        "| nC16 | 0.5000 | 0.0900 | 0.0000 | 0 |",
        "",
        "| EOS | 同流体图中的附加设置 |",
        "|---|---|",
        "| PR | PR78；Ωa=0.45724，Ωb=0.07780，δ1/δ2=1±sqrt(2)；直接使用上表 BIP |",
        "| SW | 盐度为 0；水相 CO2–H2O 使用原始 SW 关联式，CH4–H2O 使用 SW 烃类关联式，nC16–H2O 固定为 0.5000；非水相仍使用上表 BIP |",
        "| CPA | simplified 径向分布；仅 H2O 采用标准 4C 缔合，CO2、CH4、nC16 不设缔合位点；所有立方项仍使用上表 BIP |",
        "",
        "CPA 同流体图中的实际 SRK 立方参数如下。水的数值由标准 4C 参数覆盖，其余三种组分由 `a0=0.42748 R²Tc²/Pc`、`b=0.08664 RTc/Pc` 和 `c1=0.480+1.574ω-0.176ω²` 计算。",
        "",
        "| 组分 | a0 [Pa·m⁶/mol²] | b [m³/mol] | c1 | ε [J/mol] | β | 供体/受体位点 |",
        "|---|---:|---:|---:|---:|---:|---:|",
        "| H2O（4C） | 0.12277 | 1.4515×10⁻⁵ | 0.67359 | 16655 | 0.06920 | 2 / 2 |",
        "| CO2 | 0.37051015 | 2.9696952×10⁻⁵ | 0.82365531 | 0 | 0 | 0 / 0 |",
        "| CH4 | 0.23333699 | 2.9847722×10⁻⁵ | 0.49795213 | 0 | 0 | 0 / 0 |",
        "| nC16 | 10.95573838 | 3.6937867×10⁻⁴ | 1.55100874 | 0 | 0 | 0 / 0 |",
        "",
        "参数表以 `tools/example/three_eos_public_benchmark/main.cpp`、`tools/example/h2o_co2_ch4_nc16_three_eos/main.cpp` 及生产关联式头文件为依据。它们是本次计算的冻结参数快照；若以后修改这些入口，应同步更新本报告生成器并重新计算结果。",
        "",
        "## 结论摘要",
        "",
        "| EOS | 公开体系 | 相图 | 公开点收敛 | 主要误差 |",
        "|---|---|---|---:|---:|",
        f"| PR | CH4–C2H6，{PR_TEMPERATURE_K:.2f} K | P–x–y | {pr['converged_points']}/{pr['requested_points']} | 压力 AARD {pr['pressure_AARD_percent']:.3f}%；气相组成 MAE {pr['vapor_composition_MAE']:.4f} |",
        f"| SW-2019 | CO2–NaCl–H2O，323.15–423.15 K | P–mCO2 | {sw['converged_points']}/{sw['requested_points']} | 溶解度 AARD {sw['solubility_AARD_percent']:.3f}%；MAE {sw['solubility_MAE_mol_kg']:.4f} mol/kg H2O |",
        f"| SRK-CPA | 甲醇–水，{CPA_TEMPERATURE_K:.2f} K | P–x–y | {cpa['converged_points']}/{cpa['requested_points']} | 压力 AARD {cpa['pressure_AARD_percent']:.3f}%；气相组成 MAE {cpa['vapor_composition_MAE']:.4f} |",
        "",
        f"三组公开重建曲线共请求 {summary['curve_requests']} 个生产计算点，其中 PR、SW、CPA 分别为 {pr['calculation_points']}、{sw['calculation_points']}、{cpa['calculation_points']} 点。共有 {summary['curve_failures']} 个状态未检测到可用边界或未收敛；这些状态保留在逐点 CSV 中，图线与两相色块均在相应位置断开，不跨越缺失区插值。同流体 P–T 包络由独立的共同网格扫描生成，不包含在上述点数中。",
        "",
        "六类图统一采用无网格完整四边框；图例放在坐标框内，并只保留曲线、点或分组等数据标签。温度、来源说明和误差指标集中写在报告中，避免与数据重叠。",
        "",
        "## PR：甲烷–乙烷等温 P–x–y 相图",
        "",
        f"**实验条件。** May 等在论文第 3612 页表 5 的 Isothermal Path 中给出 {pr['requested_points']} 组甲烷–乙烷 VLE 数据，实测温度为 243.60–243.61 K；本图不再混入同表 243.58 K 的 Isochoric Path 数据。液相甲烷摩尔分数 x_CH4 = {pr['coordinate_range'][0]:.4f}–{pr['coordinate_range'][1]:.4f}，平衡压力 P = {pr['observed_pressure_range_Pa'][0] / 1.0e6:.3f}–{pr['observed_pressure_range_Pa'][1] / 1.0e6:.3f} MPa，对应气相甲烷摩尔分数 y_CH4 = {pr['observed_vapor_composition_range'][0]:.4f}–{pr['observed_vapor_composition_range'][1]:.4f}。计算采用 PR78、甲烷/乙烷二元作用参数 k12 = 0。自主程序在 x_CH4 = {PR_COMPOSITION_MIN:.4f}–{PR_COMPOSITION_MAX:.4f} 的完整二元组成轴上设置 {PR_CURVE_POINTS} 个均匀点，并额外并入 {pr['requested_points']} 个文献液相组成点逐点求泡点；不使用文献曲线插值代替 EOS 结果。",
        "",
        f"**图像含义。** 图窗显示 x_CH4/y_CH4 = {PR_FIGURE_X_LIMITS[0]:.2f}–{PR_FIGURE_X_LIMITS[1]:.2f}、P = {PR_FIGURE_PRESSURE_LIMITS_MPA[0]:.2f}–{PR_FIGURE_PRESSURE_LIMITS_MPA[1]:.2f} MPa，完整保留从纯乙烷饱和端点到富甲烷高压终止区的两相包络，仅去除没有气液边界的右侧空白组成区；底层计算仍覆盖完整 0–1 组成轴。蓝色实线以液相组成 x_CH4 为横坐标，是自主程序泡点/液相边界；橙色实线以同一平衡状态的气相组成 y_CH4 为横坐标，是自主程序露点/气相边界。两条边界使用相同线型，仅以颜色区分。浅蓝区域表示两条已收敛边界之间的气液两相区。May 等的液相实验点用稍大的空心圆、气相实验点用稍大的空心三角形标明，图中仅显示文献数据的中心值，不绘制误差棒。自主计算实际得到的液相边界范围为 x_CH4 = {pr['calculated_liquid_composition_range'][0]:.4f}–{pr['calculated_liquid_composition_range'][1]:.4f}，对应气相范围 y_CH4 = {pr['calculated_vapor_composition_range'][0]:.4f}–{pr['calculated_vapor_composition_range'][1]:.4f}；边界之外不外推。",
        "",
        f"![PR methane-ethane P-x-y, English]({image_paths[0].name})",
        "",
        f"![PR 甲烷乙烷 P-x-y，中文]({chinese_paths[0].name})",
        "",
        f"**结果分析。** 完整计算边界的压力范围为 {pr['calculated_pressure_range_Pa'][0] / 1.0e6:.3f}–{pr['calculated_pressure_range_Pa'][1] / 1.0e6:.3f} MPa，其中 x_CH4=0 的纯乙烷饱和端点由同一 PR EOS 的液、气逸度相等条件直接求得。边界随后连续延伸至富甲烷的终止区。由于 {PR_TEMPERATURE_K:.2f} K 高于纯甲烷临界温度 190.555 K，x_CH4=1 不存在纯甲烷气液饱和点；程序在 x_CH4>{pr['calculated_liquid_composition_range'][1]:.4f} 未检测到两相边界，因此图线在这里终止而不是外推到 1。最后一个收敛网格点只界定终止区，不作为严格临界组成。文献覆盖区内，自主程序压力 AARD 为 {pr['pressure_AARD_percent']:.3f}%，平均有符号偏差为 {pr['pressure_mean_relative_deviation_percent']:.3f}%（负值表示整体略低估）；气相甲烷组成 MAE 为 {pr['vapor_composition_MAE']:.4f}，平均有符号偏差为 {pr['vapor_composition_mean_deviation']:.4f}。{pr['converged_points']} 个等温文献组成点均成功计算。",
        "",
        f"原论文第 3612 页表 5 同时提供了组成和压力不确定度，但本版图按简洁对比需求仅绘制实验中心值；原始不确定度数据仍保留在详细 CSV 和汇总数据中。论文报告的是乙烷气相摩尔分数 y_C2H6，本图转换为 y_CH4=1-y_C2H6。更稳妥的结论是 PR 在本温度和组成窗口内同时再现了压力趋势与气相富甲烷趋势，尚不能外推到其他温度或临界邻域。",
        "",
        "## SW：CO2–NaCl–H2O 平衡溶解度图",
        "",
        f"**实验条件。** Messabeb 等的数据覆盖 T = {sw['temperature_range_K'][0]:.2f}–{sw['temperature_range_K'][1]:.2f} K、P = {sw['pressure_range_Pa'][0] / 1.0e6:.2f}–{sw['pressure_range_Pa'][1] / 1.0e6:.2f} MPa。每个温度–盐度组合包含 4 个压力点。323.15 K 含 0、1、3、6 mol NaCl/kg H2O，共 16 点；373.15 K 和 423.15 K 各含 1、3、6 mol/kg，共 12 点。纵坐标是水相 CO2 质量摩尔浓度 m_CO2，单位 mol/kg H2O。",
        "",
        "计算采用 SW 水相 alpha 和 Chabab-2019 CO2–H2O 二元参数化；盐度直接作为 mol NaCl/kg H2O 输入。每组在实验压力范围内加入 81 点压力扫描，再把公开压力点并入计算。三种画法分别突出完整趋势、逐点一致性和高温下的盐度误差，不代表三组使用了不同模型。",
        "",
    ]
    sw_analysis = {
        323.15: (
            "**结果分析。** 四条曲线均随压力升高而上升并逐渐变缓，表示增压促进 CO2 溶解，但边际增幅减小。同一压力下曲线按 0、1、3、6 mol/kg 依次降低，直接显示盐析效应。分组 AARD 分别为 "
            f"{sw_groups[(323.15, 0.0)]['AARD_percent']:.2f}%、"
            f"{sw_groups[(323.15, 1.0)]['AARD_percent']:.2f}%、"
            f"{sw_groups[(323.15, 3.0)]['AARD_percent']:.2f}% 和 "
            f"{sw_groups[(323.15, 6.0)]['AARD_percent']:.2f}%。模型对纯水和 1 mol/kg 整体略低估，对 3 mol/kg 基本无系统偏差，对 6 mol/kg 整体高估约 "
            f"{sw_groups[(323.15, 6.0)]['mean_relative_deviation_percent']:.2f}%。"
        ),
        373.15: (
            "**结果分析。** 横坐标为实验溶解度，纵坐标为自主程序预测；点落在 1:1 线表示完全一致。三组数据都靠近对角线，1、3、6 mol/kg 的 AARD 分别为 "
            f"{sw_groups[(373.15, 1.0)]['AARD_percent']:.2f}%、"
            f"{sw_groups[(373.15, 3.0)]['AARD_percent']:.2f}% 和 "
            f"{sw_groups[(373.15, 6.0)]['AARD_percent']:.2f}%。1 和 3 mol/kg 平均略低于实验，6 mol/kg 平均高于实验 "
            f"{sw_groups[(373.15, 6.0)]['mean_relative_deviation_percent']:.2f}%；偏差方向随盐度改变，说明剩余误差主要与高盐相互作用有关，而不是统一的比例偏差。"
        ),
        423.15: (
            "**结果分析。** 柱长是该盐度下 4 个压力点的 AARD，不是溶解度本身。1、3、6 mol/kg 的 AARD 为 "
            f"{sw_groups[(423.15, 1.0)]['AARD_percent']:.2f}%、"
            f"{sw_groups[(423.15, 3.0)]['AARD_percent']:.2f}% 和 "
            f"{sw_groups[(423.15, 6.0)]['AARD_percent']:.2f}%。误差随盐度总体增大；6 mol/kg 平均高估约 "
            f"{sw_groups[(423.15, 6.0)]['mean_relative_deviation_percent']:.2f}%，是该温度最需要改进的区域。"
        ),
    }
    for image_path, chinese_path, temperature in zip(
        image_paths[1:4], chinese_paths[1:4], (323.15, 373.15, 423.15)
    ):
        lines.extend(
            [
                f"### {temperature:.2f} K",
                "",
                f"![SW CO2 brine validation {temperature:.2f} K, English]({image_path.name})",
                "",
                f"![SW CO2 盐水验证 {temperature:.2f} K，中文]({chinese_path.name})",
                "",
                sw_analysis[temperature],
                "",
            ]
        )
    lines.extend(
        [
            f"三温度合并后 40 个公开点全部收敛，溶解度 AARD 为 {sw['solubility_AARD_percent']:.3f}%，MAE 为 {sw['solubility_MAE_mol_kg']:.4f} mol/kg H2O。综合结论是模型正确重现增压溶解与盐析两个主趋势；误差在 6 mol/kg 高盐端较大。323.15 K 包含纯水基线，373.15 K 和 423.15 K 的公开数据只覆盖 1、3、6 mol/kg，因此没有补造 0 mol/kg 数据。",
            "",
            "## CPA：甲醇–水等温 P–x–y 相图",
            "",
            f"**实验条件。** Kurihara 等表 2 给出 T = {CPA_TEMPERATURE_K:.2f} K 的甲醇(1)–水(2) P–x–y 数据。图中保留 12 个公开点，包括纯水和纯甲醇端点；误差指标只使用 10 个内部混合物点，其 x_MeOH = {cpa['coordinate_range'][0]:.4f}–{cpa['coordinate_range'][1]:.4f}、P = {cpa['observed_pressure_range_Pa'][0] / 1000.0:.3f}–{cpa['observed_pressure_range_Pa'][1] / 1000.0:.3f} kPa、y_MeOH = {cpa['observed_vapor_composition_range'][0]:.4f}–{cpa['observed_vapor_composition_range'][1]:.4f}。完整公开压力含端点时约为 20.0–84.5 kPa。",
            "",
            "计算采用 SRK-CPA：水为 4C、甲醇为 2B 缔合方案，k12 = -0.09。生产曲线扫描 x_MeOH = 0.02–0.965；纯组分端点只显示公开数据，不强迫泡点求解器在 x = 0 或 1 的退化组成差上定义二元两相边界。",
            "",
            f"![CPA methanol-water P-x-y, English]({image_paths[4].name})",
            "",
            f"![CPA 甲醇水 P-x-y，中文]({chinese_paths[4].name})",
            "",
            "**图像含义。** 液相数据圆点位于较低甲醇摩尔分数一侧，气相数据三角形位于较高甲醇摩尔分数一侧；同一压力下 y_MeOH > x_MeOH，说明蒸气相显著富集甲醇。浅蓝区域较宽，体现甲醇–水强非理想性导致的较大气液组成差。",
            "",
            f"**结果分析。** 10 个内部点全部收敛。压力 AARD 为 {cpa['pressure_AARD_percent']:.3f}%，且平均有符号偏差也是 {cpa['pressure_mean_relative_deviation_percent']:.3f}%，说明内部泡点压力在所有点上均被系统低估，而不是正负误差相互抵消。气相甲醇组成 MAE 为 {cpa['vapor_composition_MAE']:.4f}，平均偏差为 {cpa['vapor_composition_mean_deviation']:.4f}。曲线抓住了压力随甲醇含量升高和蒸气富甲醇的主要形态，但液相边界仍存在约 4.6% 的系统性压力偏低，不能仅凭气相边界视觉接近就判为无偏。",
            "",
            "## 同一湿流体的三 EOS P–T 包络比较",
            "",
            "**计算条件。** 该图不是公开实验。四组分顺序为 H2O/CO2/CH4/nC16，总体摩尔组成为 0.10/0.25/0.50/0.15。三个模型使用相同 15 点线性温度网格（280–800 K）、17 点对数压力网格（1–600 bar）和相同边界细化容差 1e-5。相集合被限制为 O+G，水仍是混合物组分，但不允许独立水富集相出现。因此该图回答“同一工程参数下 EOS 结构会把 O/G 边界移到哪里”，不回答哪个模型最准确。",
            "",
            "同流体组使用独立工程参数集：SW 的 nC16–水参数取显式工程值 0.5，CPA 使用标准 4C 水和非缔合烃类 SRK fallback；三套参数没有针对同一 H2O–CO2–CH4–nC16 实验数据联合回归。",
            "",
            f"![Three-EOS same-fluid P-T envelope, English]({image_paths[5].name})",
            "",
            f"![三 EOS 同流体 P-T 包络，中文]({chinese_paths[5].name})",
            "",
            "**图像含义。** 纵轴是对数压力。实线为泡点，虚线为露点；在同一温度且两条边界都存在时，两线之间是 O/G 两相压力区。图中只连接扫描与细化实际找到的边界点：PR/SW 泡点出现在约 317–614 K，露点出现在约 503–614 K；CPA 泡点仅在约 577–651 K 被找到，露点约为 503–651 K。这些离散段不应被解释为严格闭合的完整包络或严格临界轨迹。",
            "",
            "**结果分析。** restricted O/G 投影中 PR 与 SW 几乎重合：共享温度点的最大泡点差为 "
            f"{same['PR_SW_max_bubble_difference_bar']:.3f} bar，最大露点差为 {same['PR_SW_max_dew_difference_bar']:.4f} bar；"
            f"CPA 与 PR 的最大泡点差达到 {same['PR_CPA_max_bubble_difference_bar']:.1f} bar，最大露点差为 {same['PR_CPA_max_dew_difference_bar']:.2f} bar。"
            "PR/SW 重合符合 SW 修正主要作用于水富集相、而本图禁止独立水相的模型结构。CPA 的露点变化较小而泡点明显移动，说明本参数包主要改变液相侧稳定性与液相组成，而不是简单把整个包络按同一比例平移。由于无公开实验锚点、有效边界温度范围也不完全重合，这一图只能支持模型敏感性结论；正式工程判断必须用目标流体的泡点、露点和相组成数据重新定参。",
            "",
            "## 参考文献（GB/T 7714—2015）",
            "",
            "[1] MAY E F, GUO J Y, OAKLEY J H, et al. Reference Quality Vapor–Liquid Equilibrium Data for the Binary Systems Methane + Ethane, + Propane, + Butane, and + 2-Methylpropane, at Temperatures from (203 to 273) K and Pressures to 9 MPa[J]. Journal of Chemical & Engineering Data, 2015, 60(12): 3606–3620. DOI: [10.1021/acs.jced.5b00610](https://doi.org/10.1021/acs.jced.5b00610).",
            "",
            "[2] MESSABEB H, CONTAMINE F, CÉZAC P, et al. Experimental Measurement of CO2 Solubility in Aqueous NaCl Solution at Temperature from 323.15 to 423.15 K and Pressure of up to 20 MPa[J]. Journal of Chemical & Engineering Data, 2016, 61(10): 3573–3584. DOI: [10.1021/acs.jced.6b00505](https://doi.org/10.1021/acs.jced.6b00505).",
            "",
            f"[3] KURIHARA K, MINOURA T, TAKEDA K, et al. Isothermal Vapor-Liquid Equilibria for Methanol + Ethanol + Water, Methanol + Water, and Ethanol + Water[J]. Journal of Chemical & Engineering Data, 1995, 40(3): 679–684. DOI: [{kurihara['doi']}](https://doi.org/{kurihara['doi']}).",
            "",
            "## 数据文件与溯源",
            "",
            "- PR 数值来自文献 [1] 的 [NIST ThermoML 机器可读记录](https://trc.nist.gov/ThermoML/10.1021/acs.jced.5b00610.html)。",
            "- SW 数值来自文献 [2] 的 [NIST ThermoML 机器可读记录](https://trc.nist.gov/ThermoML/10.1021/acs.jced.6b00505.html)。",
            f"- CPA 数值转录自文献 [3] 的 {kurihara['source_location']}。",
            "- NIST 说明：ThermoML 数据由 NIST/TRC 人员从原论文提取；原论文仍是最终权威来源。",
            "",
            "输入文件 SHA-256：",
            "",
            f"- `may2015_pr_vle.json`: `{sha256(reference_dir / 'may2015_pr_vle.json')}`",
            f"- `messabeb2016_sw_co2_brine.json`: `{sha256(reference_dir / 'messabeb2016_sw_co2_brine.json')}`",
            f"- `kurihara1995_cpa_meoh_water_333K.json`: `{sha256(reference_dir / 'kurihara1995_cpa_meoh_water_333K.json')}`",
            "",
            "## 比较边界",
            "",
            "- 三组公开验证使用不同体系，回答各 EOS 在代表性适用域内能否再现公开相平衡；不能据此把数值大小直接横向排序为‘最佳 EOS’。",
            "- SW 图是水相 CO2 平衡溶解度边界，不等同于常规烃类二元 P–x–y 包络。",
            "- Kurihara 表中的气相组成由作者根据实验结果计算；本报告按论文给出的公开表值比较，并未把它表述为独立直接测量。",
            "- 同流体 P–T 图使用工程参数集，目的是展示模型结构敏感性；PR/SW/CPA 的参数并未针对该四组分流体做统一回归。",
            "- 本报告不对平滑曲线做额外拟合；曲线只连接相邻成功的生产计算点。",
            "",
            "## 复现",
            "",
            "运行会生成请求 CSV、生产计算逐点 CSV、六张英文 PNG/PDF、六张中文 PNG/PDF、中英文汇总 PDF、指标 JSON 和本报告。生成物位于 `results/phase_diagrams/`，按仓库规则不进入 Git。",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpmc", type=pathlib.Path)
    parser.add_argument("--same-fluid-executable", type=pathlib.Path)
    parser.add_argument(
        "--redraw-only",
        action="store_true",
        help="reuse saved production CSV files and regenerate only presentation outputs",
    )
    parser.add_argument(
        "--reference-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent / "reference_data",
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent / "results" / "phase_diagrams",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    kurihara_path = args.reference_dir / "kurihara1995_cpa_meoh_water_333K.json"
    kurihara, cpa_points = load_kurihara(kurihara_path)
    same_fluid_root = args.output_dir / "same_fluid_raw"
    detail_path = args.output_dir / "phase_diagram_detail.csv"
    if args.redraw_only:
        detail_rows = read_detail_rows(detail_path)
    else:
        if args.mpmc is None or args.same_fluid_executable is None:
            parser.error(
                "--mpmc and --same-fluid-executable are required unless --redraw-only is used"
            )
        pr_points, sw_points, _ = load_reference_data(args.reference_dir)
        request_rows, metadata = build_requests(pr_points, sw_points, cpa_points)
        request_csv = args.output_dir / "phase_diagram_requests.csv"
        result_csv = args.output_dir / "phase_diagram_mpmc_results.csv"
        write_request_csv(request_csv, request_rows)
        run_mpmc_executable(args.mpmc, request_csv, result_csv)
        outputs = read_csv_by_id(result_csv)
        detail_rows = collect_rows(metadata, outputs)
        write_rows(detail_path, detail_rows)
        subprocess.run(
            [str(args.same_fluid_executable), str(same_fluid_root)], check=True
        )
    summary = summarize(detail_rows)
    envelopes = read_same_fluid_envelopes(same_fluid_root)
    summary["same_fluid"] = summarize_envelopes(envelopes)

    image_paths_by_language: dict[str, list[pathlib.Path]] = {}
    figure_books: dict[str, pathlib.Path] = {}
    for language in ("en", "zh"):
        image_paths = [make_pr_figure(detail_rows, args.output_dir, language)]
        image_paths.extend(
            make_sw_figures(detail_rows, args.output_dir, language)
        )
        image_paths.append(
            make_cpa_figure(detail_rows, cpa_points, args.output_dir, language)
        )
        image_paths.append(
            make_same_fluid_figure(envelopes, args.output_dir, language)
        )
        image_paths_by_language[language] = image_paths
        figure_books[language] = write_book(
            args.output_dir,
            image_paths,
            f"phase_diagram_figures_{language}.pdf",
        )

    (args.output_dir / "phase_diagram_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_report(
        args.output_dir / "phase_diagram_report.md",
        summary,
        image_paths_by_language,
        kurihara,
        args.reference_dir,
    )
    provenance = {
        "python": sys.version,
        "numpy": importlib.metadata.version("numpy"),
        "matplotlib": importlib.metadata.version("matplotlib"),
        "redraw_only": args.redraw_only,
        "detail_csv_sha256": sha256(detail_path),
        "figure_files": {
            language: [path.name for path in paths]
            for language, paths in image_paths_by_language.items()
        },
        "figure_books": {
            language: book.name for language, book in figure_books.items()
        },
    }
    if args.mpmc is not None:
        provenance["mpmc_executable"] = str(args.mpmc.resolve())
        provenance["mpmc_executable_sha256"] = sha256(args.mpmc)
    if args.same_fluid_executable is not None:
        provenance["same_fluid_executable"] = str(
            args.same_fluid_executable.resolve()
        )
        provenance["same_fluid_executable_sha256"] = sha256(
            args.same_fluid_executable
        )
    (args.output_dir / "phase_diagram_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"[phase-diagram] report: {args.output_dir / 'phase_diagram_report.md'}")
    for language, book in figure_books.items():
        print(f"[phase-diagram] {language} figure book: {book}")


if __name__ == "__main__":
    main()
