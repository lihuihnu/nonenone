#!/usr/bin/env python3
"""Draw four pressure-resolved parity panels for the CPA VLLE benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from plot_results import (
    AXES_BOX_ASPECT,
    COMPONENT_LABELS,
    COMPONENT_MARKERS,
    PHASE_LABELS,
    PHASE_STYLES,
    configure_style,
    read_rows,
    register_chinese_font,
)


FIGURE_SIZE_INCHES = (7.1, 7.0)
EXPECTED_PRESSURES = (19.0, 31.7, 46.0, 58.8)


def draw_panel(ax, rows: list[dict[str, str]], pressure: float) -> None:
    """Draw the nine phase-component comparisons at one pressure."""
    ax.plot(
        [0.0, 1.0],
        [0.0, 1.0],
        color="#4B5563",
        linewidth=1.0,
        linestyle="--",
        dashes=(4.0, 2.5),
        zorder=1,
    )
    for row in rows:
        phase = row["phase"]
        component = row["component"]
        phase_style = PHASE_STYLES[phase]
        experimental = float(row["experimental"])
        calculated = float(row["calculated"])
        ax.plot(
            [experimental, experimental],
            [experimental, calculated],
            color=phase_style["edgecolor"],
            linewidth=0.55,
            alpha=0.40,
            solid_capstyle="round",
            zorder=2,
        )
        ax.scatter(
            experimental,
            calculated,
            s=38,
            marker=COMPONENT_MARKERS[component],
            facecolor=phase_style["facecolor"],
            edgecolor=phase_style["edgecolor"],
            linewidth=1.05,
            hatch=phase_style["hatch"],
            alpha=0.97,
            zorder=3,
        )

    ax.set_xlim(-0.018, 1.018)
    ax.set_ylim(-0.018, 1.018)
    ticks = [index / 5.0 for index in range(6)]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_box_aspect(AXES_BOX_ASPECT)
    ax.set_title(f"{pressure:.1f} bar", fontsize=9.0, pad=5.0)
    ax.grid(False)
    ax.tick_params(
        direction="in",
        length=3.0,
        width=0.7,
        top=True,
        right=True,
        pad=2.5,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.85)


def add_phase_legend(ax) -> None:
    handles = [
        Patch(
            facecolor=PHASE_STYLES[phase]["facecolor"],
            edgecolor=PHASE_STYLES[phase]["edgecolor"],
            linewidth=1.05,
            hatch=PHASE_STYLES[phase]["hatch"],
            label=PHASE_LABELS[phase],
        )
        for phase in ("water-rich", "DME-rich", "vapor")
    ]
    legend = ax.legend(
        handles=handles,
        title="Phase",
        loc="upper left",
        fontsize=6.5,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="#9CA3AF",
        borderpad=0.45,
        handlelength=1.65,
        handletextpad=0.5,
        labelspacing=0.36,
    )
    legend.get_title().set_fontsize(6.6)


def add_component_legend(ax) -> None:
    handles = [
        Line2D(
            [],
            [],
            linestyle="none",
            marker=COMPONENT_MARKERS[component],
            markersize=5.2,
            markerfacecolor="#4B5563",
            markeredgecolor="white",
            label=COMPONENT_LABELS[component],
        )
        for component in ("H2O", "DME", "CO2")
    ]
    handles.append(
        Line2D(
            [],
            [],
            color="#4B5563",
            linewidth=1.0,
            linestyle="--",
            dashes=(4.0, 2.5),
            label="1:1",
        )
    )
    legend = ax.legend(
        handles=handles,
        title="Component",
        loc="lower right",
        fontsize=6.5,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="#9CA3AF",
        borderpad=0.45,
        handletextpad=0.5,
        labelspacing=0.36,
    )
    legend.get_title().set_fontsize(6.6)


def draw(
    rows: list[dict[str, str]],
    output_base: pathlib.Path,
    language: str,
    chinese_font: str,
) -> None:
    font_family = chinese_font if language == "zh" else "DejaVu Sans"
    configure_style(font_family)
    fig, axes = plt.subplots(
        2,
        2,
        figsize=FIGURE_SIZE_INCHES,
        sharex=True,
        sharey=True,
        layout="constrained",
    )

    pressures = tuple(sorted({float(row["pressure_bar"]) for row in rows}))
    if pressures != EXPECTED_PRESSURES:
        raise RuntimeError(f"Expected pressures {EXPECTED_PRESSURES}, found {pressures}.")

    for ax, pressure in zip(axes.flat, pressures):
        pressure_rows = [
            row for row in rows if abs(float(row["pressure_bar"]) - pressure) < 1.0e-9
        ]
        if len(pressure_rows) != 9:
            raise RuntimeError(
                f"Expected 9 phase-component values at {pressure:.1f} bar, "
                f"found {len(pressure_rows)}."
            )
        draw_panel(ax, pressure_rows, pressure)

    add_phase_legend(axes[0, 0])
    add_component_legend(axes[1, 1])
    fig.supxlabel(
        "实验摩尔分数" if language == "zh" else "Experimental mole fraction",
        fontsize=10.2,
    )
    fig.supylabel(
        "Our CPA 摩尔分数" if language == "zh" else "Our CPA mole fraction",
        fontsize=10.2,
    )

    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".png"), dpi=450, transparent=False)
    fig.savefig(output_base.with_suffix(".pdf"), transparent=False)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()

    rows = read_rows(args.input)
    chinese_font = register_chinese_font()
    stem = "02_cpa_vlle_composition_by_pressure"
    for language in ("en", "zh"):
        draw(rows, args.output_dir / language / stem, language, chinese_font)

    manifest = {
        "source_csv": str(args.input),
        "source_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "rows": len(rows),
        "panels": [f"{pressure:.1f} bar" for pressure in EXPECTED_PRESSURES],
        "transformation": "none; rows separated by pressure only",
        "encoding": {
            "phase": "color plus solid/open/hatched fill",
            "component": "marker shape",
            "reference": "1:1 dashed line in every panel",
            "residual": "vertical segment from y=x to the calculated value",
        },
        "figure_inches": list(FIGURE_SIZE_INCHES),
        "axes_box_aspect": AXES_BOX_ASPECT,
        "png_dpi": 450,
        "formats": ["png", "pdf"],
        "languages": ["en", "zh"],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / f"{stem}_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
