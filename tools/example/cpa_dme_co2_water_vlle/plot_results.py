#!/usr/bin/env python3
"""Draw a publication-style parity plot for the public CPA VLLE benchmark."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


PHASE_STYLES = {
    "water-rich": {"facecolor": "#0072B2", "edgecolor": "#005A8D", "hatch": None},
    "DME-rich": {"facecolor": "#FFF3E8", "edgecolor": "#D55E00", "hatch": None},
    "vapor": {"facecolor": "#E3F3ED", "edgecolor": "#009E73", "hatch": "///"},
}
COMPONENT_MARKERS = {"H2O": "o", "DME": "s", "CO2": "^"}
COMPONENT_LABELS = {"H2O": r"H$_2$O", "DME": "DME", "CO2": r"CO$_2$"}
FIGURE_SIZE_INCHES = (5.15, 4.45)
AXES_BOX_ASPECT = 1.0
PHASE_LABELS = {
    "water-rich": "Water-rich liquid",
    "DME-rich": "DME-rich liquid",
    "vapor": "Vapor",
}


def register_chinese_font() -> str:
    candidates = [
        pathlib.Path("/mnt/c/Windows/Fonts/msyh.ttc"),
        pathlib.Path("/mnt/c/Windows/Fonts/simhei.ttf"),
        pathlib.Path("C:/Windows/Fonts/msyh.ttc"),
        pathlib.Path("C:/Windows/Fonts/simhei.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            font_manager.fontManager.addfont(candidate)
            return font_manager.FontProperties(fname=candidate).get_name()
    raise RuntimeError("A Chinese font is required for bilingual figure export.")


def read_rows(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    required = {
        "pressure_bar",
        "phase",
        "component",
        "experimental",
        "calculated",
    }
    if not rows or not required.issubset(rows[0]):
        raise RuntimeError(f"Invalid phase-comparison CSV: {path}")
    if len(rows) != 36:
        raise RuntimeError(f"Expected 36 phase-composition values, found {len(rows)}.")
    return rows


def configure_style(font_family: str) -> None:
    matplotlib.rcParams.update(
        {
            "font.family": font_family,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.85,
            "axes.labelsize": 10.2,
            "xtick.labelsize": 8.6,
            "ytick.labelsize": 8.6,
            "legend.fontsize": 7.4,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.unicode_minus": False,
        }
    )


def draw(
    rows: list[dict[str, str]],
    output_base: pathlib.Path,
    language: str,
    chinese_font: str,
) -> None:
    font_family = chinese_font if language == "zh" else "DejaVu Sans"
    configure_style(font_family)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_INCHES, layout="constrained")

    ax.plot(
        [0.0, 1.0],
        [0.0, 1.0],
        color="#4B5563",
        linewidth=1.15,
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
        # The short residual segment makes the sign and magnitude of each
        # deviation legible without changing or displacing any observation.
        ax.plot(
            [experimental, experimental],
            [experimental, calculated],
            color=phase_style["edgecolor"],
            linewidth=0.65,
            alpha=0.38,
            solid_capstyle="round",
            zorder=2,
        )
        ax.scatter(
            experimental,
            calculated,
            s=46,
            marker=COMPONENT_MARKERS[component],
            facecolor=phase_style["facecolor"],
            edgecolor=phase_style["edgecolor"],
            linewidth=1.15,
            hatch=phase_style["hatch"],
            alpha=0.97,
            zorder=3,
        )

    # A small visual margin keeps near-zero trace-component markers inside the
    # frame while the displayed physical ticks remain exactly 0--1.
    ax.set_xlim(-0.018, 1.018)
    ax.set_ylim(-0.018, 1.018)
    ticks = [index / 5.0 for index in range(6)]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    # Keep the four-spine plotting frame itself square. The surrounding canvas
    # remains rectangular so axis labels have balanced publication margins.
    ax.set_box_aspect(AXES_BOX_ASPECT)
    ax.set_xlabel("实验摩尔分数" if language == "zh" else "Experimental mole fraction")
    ax.set_ylabel("Our CPA 摩尔分数" if language == "zh" else "Our CPA mole fraction")
    ax.grid(False)
    ax.tick_params(
        direction="in",
        length=3.2,
        width=0.75,
        top=True,
        right=True,
        pad=3.0,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)

    phase_handles = [
        Patch(
            facecolor=PHASE_STYLES[phase]["facecolor"],
            edgecolor=PHASE_STYLES[phase]["edgecolor"],
            linewidth=1.15,
            hatch=PHASE_STYLES[phase]["hatch"],
            label=PHASE_LABELS[phase],
        )
        for phase in ("water-rich", "DME-rich", "vapor")
    ]
    phase_legend = ax.legend(
        handles=phase_handles,
        title="Phase",
        loc="upper left",
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="#9CA3AF",
        borderpad=0.55,
        handlelength=1.75,
        handletextpad=0.55,
        labelspacing=0.42,
    )
    phase_legend.get_title().set_fontsize(7.4)
    ax.add_artist(phase_legend)

    component_handles = [
        Line2D(
            [],
            [],
            linestyle="none",
            marker=COMPONENT_MARKERS[component],
            markersize=6.0,
            markerfacecolor="#4B5563",
            markeredgecolor="white",
            label=COMPONENT_LABELS[component],
        )
        for component in ("H2O", "DME", "CO2")
    ]
    component_handles.append(
        Line2D(
            [],
            [],
            color="#4B5563",
            linewidth=1.25,
            linestyle="--",
            dashes=(4.0, 2.5),
            label="1:1",
        )
    )
    component_legend = ax.legend(
        handles=component_handles,
        title="Component",
        loc="lower right",
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="#9CA3AF",
        borderpad=0.55,
        handletextpad=0.55,
        labelspacing=0.42,
    )
    component_legend.get_title().set_fontsize(7.4)

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
    stem = "01_cpa_vlle_composition_parity"
    for language in ("en", "zh"):
        draw(rows, args.output_dir / language / stem, language, chinese_font)

    manifest = {
        "source_csv": str(args.input),
        "source_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "rows": len(rows),
        "transformation": "none; experimental values on x and calculated values on y",
        "encoding": {
            "phase": "color plus solid/open/hatched fill",
            "component": "marker shape",
            "reference": "1:1 dashed line",
            "residual": "vertical segment from y=x to the calculated value",
        },
        "figure_inches": list(FIGURE_SIZE_INCHES),
        "axes_box_aspect": AXES_BOX_ASPECT,
        "png_dpi": 450,
        "formats": ["png", "pdf"],
        "languages": ["en", "zh"],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
