#!/usr/bin/env python3
"""绘制 Heringer 2025 文献 PR 与程序 PR/SW/CPA 的三相闪蒸对比图。"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import LogLocator, PercentFormatter
from PIL import Image


COMPONENTS = [
    r"H$_2$O",
    r"CO$_2$",
    r"C$_1$",
    r"C$_{2-3}$",
    r"C$_{4-6}$",
    r"C$_{7-15}$",
    r"C$_{16-27}$",
    r"C$_{28+}$",
]
PHASES = ["oil", "gas", "water"]
PHASE_LABELS = {"oil": "油相", "gas": "气相", "water": "水相"}

MODEL_ROWS = {
    "程序 PR78": ("PR", "PR78"),
    "标准 SW": ("SW", "SW_water_alpha_paper_BIPs"),
    "标准 SRK-CPA": ("CPA", "SRK-CPA_4C-water"),
}

MODEL_STYLE = {
    "文献 PR": {"color": "#111111", "marker": "o", "linestyle": "-", "linewidth": 2.5},
    "程序 PR78": {"color": "#0072B2", "marker": "s", "linestyle": "--", "linewidth": 2.1},
    "标准 SW": {"color": "#D55E00", "marker": "^", "linestyle": "-.", "linewidth": 2.1},
    "标准 SRK-CPA": {"color": "#009E73", "marker": "D", "linestyle": ":", "linewidth": 2.4},
}

PHASE_STYLE = {
    "oil": {"color": "#D55E00", "hatch": "///"},
    "gas": {"color": "#0072B2", "hatch": "..."},
    "water": {"color": "#009E73", "hatch": "\\\\"},
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def find_summary_row(
    rows: list[dict[str, str]], backend: str, parameterization: str
) -> dict[str, str]:
    selected = [
        row
        for row in rows
        if row["backend"] == backend
        and row["parameterization"] == parameterization
        and row["path"] == "unrestricted"
    ]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one unrestricted {backend}/{parameterization} row, got {len(selected)}"
        )
    return selected[0]


def load_data(input_dir: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    summary_rows = read_csv(input_dir / "summary.csv")
    phase_rows = read_csv(input_dir / "phase_comparison.csv")

    pr_row = find_summary_row(summary_rows, "PR", "PR78")
    data: dict[str, Any] = {
        "文献 PR": {
            "beta": [
                float(pr_row["reference_beta_oil"]),
                float(pr_row["reference_beta_gas"]),
                float(pr_row["reference_beta_water"]),
            ],
            "composition_mae": None,
            "beta_mae": None,
            "composition": {phase: {} for phase in PHASES},
        }
    }

    for label, (backend, parameterization) in MODEL_ROWS.items():
        row = find_summary_row(summary_rows, backend, parameterization)
        data[label] = {
            "beta": [
                float(row["beta_oil"]),
                float(row["beta_gas"]),
                float(row["beta_water"]),
            ],
            "composition_mae": float(row["composition_mae"]),
            "beta_mae": float(row["beta_mae"]),
            "composition": {phase: {} for phase in PHASES},
        }

    reference_loaded = False
    for row in phase_rows:
        if row["backend"] != "PR" or row["parameterization"] != "PR78" or row["path"] != "unrestricted":
            continue
        data["文献 PR"]["composition"][row["phase"]][row["component"]] = float(
            row["reference_normalized"]
        )
        reference_loaded = True
    if not reference_loaded:
        raise ValueError("Could not locate normalized literature compositions")

    for label, (backend, parameterization) in MODEL_ROWS.items():
        for row in phase_rows:
            if (
                row["backend"] == backend
                and row["parameterization"] == parameterization
                and row["path"] == "unrestricted"
            ):
                data[label]["composition"][row["phase"]][row["component"]] = float(
                    row["calculated"]
                )

    for label, model in data.items():
        for phase in PHASES:
            if len(model["composition"][phase]) != len(COMPONENTS):
                raise ValueError(f"Incomplete composition for {label}/{phase}")
    return data, summary_rows


def configure_fonts() -> tuple[FontProperties, str]:
    candidates = [
        Path("C:/Windows/Fonts/msyh.ttc"),
        Path("/mnt/c/Windows/Fonts/msyh.ttc"),
    ]
    font_path = next((path for path in candidates if path.exists()), None)
    if font_path is None:
        return FontProperties(), "Matplotlib default"
    font = FontProperties(fname=str(font_path))
    mpl.font_manager.fontManager.addfont(str(font_path))
    family = font.get_name()
    mpl.rcParams["font.family"] = family
    mpl.rcParams["font.sans-serif"] = [family, "DejaVu Sans"]
    return font, f"{family} ({font_path})"


def plot_phase_fractions(ax: plt.Axes, data: dict[str, Any]) -> None:
    labels = list(data)
    y_positions = list(range(len(labels)))
    left = [0.0] * len(labels)
    for phase_index, phase in enumerate(PHASES):
        values = [data[label]["beta"][phase_index] for label in labels]
        style = PHASE_STYLE[phase]
        bars = ax.barh(
            y_positions,
            values,
            left=left,
            height=0.64,
            color=style["color"],
            edgecolor="white",
            linewidth=1.1,
            hatch=style["hatch"],
            label=PHASE_LABELS[phase],
            zorder=3,
        )
        for bar, value in zip(bars, values, strict=True):
            if value >= 0.075:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{100 * value:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=9.2,
                    color="white",
                    fontweight="bold",
                    zorder=4,
                )
        left = [a + b for a, b in zip(left, values, strict=True)]

    ax.set_xlim(0.0, 1.0)
    ax.set_yticks(y_positions, labels)
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xlabel("总进料的相摩尔分数")
    ax.set_title("a  三相分配：文献与三个状态方程", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#D9DEE7", linewidth=0.8, zorder=0)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)


def plot_errors(ax: plt.Axes, data: dict[str, Any]) -> None:
    labels = ["程序 PR78", "标准 SW", "标准 SRK-CPA"]
    y_positions = list(range(len(labels)))
    composition = [data[label]["composition_mae"] for label in labels]
    beta = [data[label]["beta_mae"] for label in labels]
    ax.scatter(
        composition,
        y_positions,
        s=76,
        marker="o",
        color="#0072B2",
        edgecolor="white",
        linewidth=0.8,
        label="三相组成 MAE",
        zorder=3,
    )
    ax.scatter(
        beta,
        y_positions,
        s=76,
        marker="s",
        color="#D55E00",
        edgecolor="white",
        linewidth=0.8,
        label="相分数 MAE",
        zorder=3,
    )
    for y, value in enumerate(beta):
        ax.annotate(
            f"{value:.1e}",
            (value, y),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            fontsize=8.3,
            color="#7A3200",
        )
    ax.set_xscale("log")
    ax.set_xlim(1.0e-5, 1.0e-1)
    ax.set_yticks(y_positions, labels)
    ax.invert_yaxis()
    ax.set_xlabel("相对文献 PR 的平均绝对误差（对数轴）")
    ax.set_title("b  误差量级", loc="left", fontweight="bold")
    ax.grid(axis="x", which="both", color="#D9DEE7", linewidth=0.75)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False, ncol=1)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)


def plot_composition(ax: plt.Axes, data: dict[str, Any], phase: str, panel: str) -> None:
    source_component_names = ["H2O", "CO2", "C1", "C2-3", "C4-6", "C7-15", "C16-27", "C28+"]
    x = list(range(len(source_component_names)))
    for label, model in data.items():
        values = [model["composition"][phase][component] for component in source_component_names]
        style = MODEL_STYLE[label]
        ax.plot(
            x,
            values,
            label=label,
            color=style["color"],
            marker=style["marker"],
            linestyle=style["linestyle"],
            linewidth=style["linewidth"],
            markersize=5.5,
            markerfacecolor="white" if label != "文献 PR" else style["color"],
            markeredgewidth=1.3,
            zorder=4 if label == "文献 PR" else 3,
        )

    ax.set_xticks(x, COMPONENTS, rotation=35, ha="right")
    ax.set_title(f"{panel}  {PHASE_LABELS[phase]}组成指纹", loc="left", fontweight="bold")
    ax.set_xlabel("组分")
    ax.grid(axis="y", color="#D9DEE7", linewidth=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if phase == "water":
        ax.set_yscale("log")
        ax.set_ylim(1.0e-10, 1.2)
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
        ax.set_ylabel(r"水相摩尔分数（log$_{10}$）")
        ax.text(
            0.02,
            0.05,
            "采用对数轴以保留痕量重组分",
            transform=ax.transAxes,
            fontsize=8.4,
            color="#555B66",
        )
    else:
        ax.set_ylim(0.0, 0.78)
        ax.yaxis.set_major_formatter(PercentFormatter(1.0))
        ax.set_ylabel(f"{PHASE_LABELS[phase]}摩尔分数")


def write_plot_data(path: Path, data: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["source", "phase", "component", "phase_beta", "phase_composition"])
        source_names = ["H2O", "CO2", "C1", "C2-3", "C4-6", "C7-15", "C16-27", "C28+"]
        for label, model in data.items():
            for phase_index, phase in enumerate(PHASES):
                for component in source_names:
                    writer.writerow(
                        [
                            label,
                            phase,
                            component,
                            f"{model['beta'][phase_index]:.17g}",
                            f"{model['composition'][phase][component]:.17g}",
                        ]
                    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    base = output_dir / "heringer2025_standard_eos_comparison"
    outputs = [base.with_suffix(suffix) for suffix in [".png", ".pdf", ".svg"]]
    outputs += [output_dir / "figure_data.csv", output_dir / "figure_manifest.json"]
    existing = [path for path in outputs if path.exists()]
    if existing and not args.force:
        raise FileExistsError("Refusing to overwrite: " + ", ".join(map(str, existing)))

    data, _ = load_data(args.input_dir.resolve())
    _, font_description = configure_fonts()
    mpl.rcParams.update(
        {
            "axes.unicode_minus": False,
            "axes.labelcolor": "#20242B",
            "axes.edgecolor": "#8B929E",
            "axes.linewidth": 0.8,
            "xtick.color": "#3F4650",
            "ytick.color": "#3F4650",
            "text.color": "#20242B",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    with mpl.rc_context():
        fig = plt.figure(figsize=(15.0, 9.1), facecolor="white")
        grid = fig.add_gridspec(
            2,
            3,
            height_ratios=[0.82, 1.28],
            width_ratios=[1, 1, 1],
            hspace=0.34,
            wspace=0.22,
        )
        fig.subplots_adjust(left=0.065, right=0.985, top=0.875, bottom=0.145)
        phase_ax = fig.add_subplot(grid[0, :2])
        error_ax = fig.add_subplot(grid[0, 2])
        oil_ax = fig.add_subplot(grid[1, 0])
        gas_ax = fig.add_subplot(grid[1, 1])
        water_ax = fig.add_subplot(grid[1, 2])

        plot_phase_fractions(phase_ax, data)
        plot_errors(error_ax, data)
        plot_composition(oil_ax, data, "oil", "c")
        plot_composition(gas_ax, data, "gas", "d")
        plot_composition(water_ax, data, "water", "e")

        handles, labels = oil_ax.get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.045),
            ncol=4,
            frameon=False,
            fontsize=10.5,
            handlelength=3.0,
        )
        fig.suptitle(
            "Heringer 2025 油–气–水三相闪蒸：PR、SW、CPA 与文献 PR 对比",
            fontsize=17,
            fontweight="bold",
            y=0.972,
        )
        fig.text(
            0.5,
            0.012,
            "650 K，390 bar｜Heringer et al. (2025), Table 6｜标准参数、无目标拟合｜各相组成按相内摩尔分数归一化",
            ha="center",
            va="bottom",
            fontsize=9.5,
            color="#555B66",
        )

        fig.savefig(outputs[0], dpi=300, facecolor="white", metadata={"Software": f"Matplotlib {mpl.__version__}"})
        fig.savefig(outputs[1], facecolor="white", metadata={"Creator": f"Matplotlib {mpl.__version__}"})
        fig.savefig(outputs[2], facecolor="white", metadata={"Creator": f"Matplotlib {mpl.__version__}"})
        plt.close(fig)

    # Matplotlib writes PNG as RGBA even with an opaque white canvas.  Convert
    # the generated presentation copy to true RGB so downstream page
    # backgrounds cannot alter its appearance.
    temporary_png = outputs[0].with_name(outputs[0].stem + ".rgb.png")
    with Image.open(outputs[0]) as png:
        png.convert("RGB").save(temporary_png, dpi=(300, 300), optimize=True)
    temporary_png.replace(outputs[0])

    write_plot_data(output_dir / "figure_data.csv", data)
    manifest = {
        "title": "Heringer 2025 standard-parameter PR/SW/CPA three-phase comparison",
        "source_files": [
            str((args.input_dir / "summary.csv").resolve()),
            str((args.input_dir / "phase_comparison.csv").resolve()),
        ],
        "transformations": [
            "Selected unrestricted PR78, SW_water_alpha_paper_BIPs and SRK-CPA_4C-water rows",
            "Used normalized Table 6 phase compositions already recorded in phase_comparison.csv",
            "Displayed water-phase composition on a base-10 logarithmic axis; no values were removed or replaced",
        ],
        "state": {"temperature_K": 650.0, "pressure_bar": 390.0},
        "outputs": [str(path) for path in outputs[:3]],
        "figure_inches": [15.0, 9.1],
        "png_dpi": 300,
        "matplotlib_version": mpl.__version__,
        "font": font_description,
        "background": "opaque white",
        "publisher_profile": None,
    }
    with (output_dir / "figure_manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
