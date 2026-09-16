#!/usr/bin/env python3
"""Create publication-style PR/SW/CPA comparison figures.

The script reads the unmodified CSV outputs written by this case.  It does not
smooth, interpolate, resample, or omit simulation output times.  Chinese and
English figures share identical data, dimensions, limits, and styles; only the
axis text changes.  Legend text remains English in both versions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import AutoMinorLocator, FuncFormatter, MaxNLocator


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
COMPARISON = RESULTS / "comparison"

EOS_ORDER = ("pr", "sw", "cpa")
EOS_LABELS = {"pr": "PR", "sw": "SW", "cpa": "CPA"}
EOS_COLORS = {"pr": "#0072B2", "sw": "#D55E00", "cpa": "#009E73"}
EOS_LINESTYLES = {"pr": "-", "sw": "--", "cpa": "-."}
EOS_MARKERS = {"pr": "o", "sw": "s", "cpa": "^"}

FIGURE_SIZE_MM = (120.0, 80.0)
DEFAULT_DPI = 300
MARK_EVERY = 12

matplotlib.rcParams.update(
    {
        "axes.unicode_minus": False,
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=RESULTS / "figures",
        help="output root (default: results/figures)",
    )
    parser.add_argument(
        "--language",
        choices=("all", "en", "zh"),
        default="all",
        help="figure language (default: all)",
    )
    parser.add_argument(
        "--formats",
        default="png,pdf",
        help="comma-separated output formats (default: png,pdf)",
    )
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace existing generated figures",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing required result file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"empty required result file: {path}")
    return rows


def numeric(row: dict[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, ValueError) as error:
        raise SystemExit(f"invalid numeric field {field!r}: {row.get(field)!r}") from error
    if not math.isfinite(value):
        raise SystemExit(f"non-finite numeric field {field!r}: {value}")
    return value


def grouped_by_eos(rows: Iterable[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {eos: [] for eos in EOS_ORDER}
    for row in rows:
        eos = row.get("eos", "").lower()
        if eos in grouped:
            grouped[eos].append(row)
    for eos in EOS_ORDER:
        grouped[eos].sort(key=lambda row: numeric(row, "time"))
        if not grouped[eos]:
            raise SystemExit(f"missing {eos.upper()} rows")
    return grouped


def register_chinese_font() -> str:
    candidates = (
        Path("C:/Windows/Fonts/simhei.ttf"),
        Path("/mnt/c/Windows/Fonts/simhei.ttf"),
    )
    for candidate in candidates:
        if candidate.exists():
            font_manager.fontManager.addfont(str(candidate))
            return font_manager.FontProperties(fname=str(candidate)).get_name()
    raise SystemExit("Chinese font not found (expected Windows SimHei)")


def validate_common_time_grid(
    reservoir: dict[str, list[dict[str, str]]],
    wells: dict[str, list[dict[str, str]]],
) -> None:
    reference = [numeric(row, "time") for row in reservoir["pr"]]
    if len(reference) != 121 or reference[0] != 0.0 or reference[-1] != 60.0:
        raise SystemExit("expected 121 common reservoir states from 0 to 60 days")
    for eos in EOS_ORDER:
        times = [numeric(row, "time") for row in reservoir[eos]]
        if times != reference:
            raise SystemExit(f"reservoir time grid differs for {eos.upper()}")
        for row in reservoir[eos]:
            saturation_sum = sum(numeric(row, field) for field in ("s_o_avg", "s_g_avg", "s_w_avg"))
            if abs(saturation_sum - 1.0) > 1.0e-10:
                raise SystemExit(
                    f"average saturation closure failed for {eos.upper()} at day {row['time']}"
                )
        for well_name in ("CO2_INJ", "PROD"):
            well_times = [
                numeric(row, "time") for row in wells[eos] if row.get("name") == well_name
            ]
            if well_times != reference:
                raise SystemExit(f"{well_name} time grid differs for {eos.upper()}")


def style_axes(ax: plt.Axes) -> None:
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#111827")
        spine.set_linewidth(0.9)
    ax.tick_params(
        which="major",
        direction="in",
        top=True,
        right=True,
        length=4.0,
        width=0.85,
        color="#111827",
    )
    ax.tick_params(
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=2.3,
        width=0.7,
        color="#111827",
    )
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))


def style_legend(legend: object) -> None:
    latin = font_manager.FontProperties(family="DejaVu Sans", size=8.2)
    for text in legend.get_texts():
        text.set_fontproperties(latin)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("#CBD5E1")
    frame.set_linewidth(0.7)
    frame.set_alpha(0.94)


def new_figure(language: str, chinese_font: str) -> tuple[plt.Figure, plt.Axes]:
    family = chinese_font if language == "zh" else "DejaVu Sans"
    width, height = (value / 25.4 for value in FIGURE_SIZE_MM)
    with matplotlib.rc_context(
        {
            "font.family": family,
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.4,
            "ytick.labelsize": 8.4,
            "legend.fontsize": 8.2,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    ):
        fig, ax = plt.subplots(figsize=(width, height), layout="constrained")
    style_axes(ax)
    return fig, ax


def localize(language: str, english: str, chinese: str) -> str:
    return chinese if language == "zh" else english


def plot_eos_lines(
    ax: plt.Axes,
    grouped: dict[str, list[dict[str, str]]],
    value: Callable[[dict[str, str]], float],
) -> None:
    for eos in EOS_ORDER:
        rows = grouped[eos]
        ax.plot(
            [numeric(row, "time") for row in rows],
            [value(row) for row in rows],
            color=EOS_COLORS[eos],
            linestyle=EOS_LINESTYLES[eos],
            linewidth=1.8,
            marker=EOS_MARKERS[eos],
            markersize=4.0,
            markerfacecolor="white",
            markeredgecolor=EOS_COLORS[eos],
            markeredgewidth=0.9,
            markevery=MARK_EVERY,
            label=EOS_LABELS[eos],
            zorder=3,
        )
    ax.set_xlim(0.0, 60.0)
    ax.set_xticks(range(0, 61, 10))
    legend = ax.legend(loc="best", ncols=1, borderpad=0.55, handlelength=2.8)
    style_legend(legend)


def atomic_save(fig: plt.Figure, target: Path, dpi: int, force: bool) -> None:
    if target.exists() and not force:
        raise SystemExit(f"refusing to overwrite existing figure: {target} (use --force)")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.tmp")
    try:
        fig.savefig(
            temporary,
            format=target.suffix.lstrip("."),
            dpi=dpi,
            facecolor="white",
            transparent=False,
        )
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def save_formats(
    fig: plt.Figure,
    stem: Path,
    formats: tuple[str, ...],
    dpi: int,
    force: bool,
) -> list[Path]:
    outputs = []
    for file_format in formats:
        target = stem.with_suffix(f".{file_format}")
        atomic_save(fig, target, dpi, force)
        outputs.append(target)
    plt.close(fig)
    return outputs


def reservoir_figure(
    rows: dict[str, list[dict[str, str]]],
    field: str,
    ylabel_en: str,
    ylabel_zh: str,
    language: str,
    chinese_font: str,
) -> plt.Figure:
    fig, ax = new_figure(language, chinese_font)
    plot_eos_lines(ax, rows, lambda row: numeric(row, field))
    ax.set_xlabel(localize(language, "Time (day)", "时间（天）"))
    ax.set_ylabel(localize(language, ylabel_en, ylabel_zh))
    return fig


def well_figure(
    well_rows: dict[str, list[dict[str, str]]],
    well_name: str,
    field: str,
    transform: Callable[[float], float],
    ylabel_en: str,
    ylabel_zh: str,
    language: str,
    chinese_font: str,
) -> plt.Figure:
    selected = {
        eos: [row for row in well_rows[eos] if row.get("name") == well_name]
        for eos in EOS_ORDER
    }
    fig, ax = new_figure(language, chinese_font)
    plot_eos_lines(ax, selected, lambda row: transform(numeric(row, field)))
    ax.set_xlabel(localize(language, "Time (day)", "时间（天）"))
    ax.set_ylabel(localize(language, ylabel_en, ylabel_zh))
    return fig


def load_mass_balance_history() -> dict[str, list[dict[str, float]]]:
    grouped: dict[str, list[dict[str, float]]] = {}
    for eos in EOS_ORDER:
        by_step: dict[int, dict[str, object]] = defaultdict(
            lambda: {"time": 0.0, "maximum": 0.0}
        )
        for row in read_csv(RESULTS / eos / "component_mass_balance.csv"):
            step = int(row["step"])
            relative_error = abs(numeric(row, "relative_error"))
            by_step[step]["time"] = numeric(row, "time_day")
            by_step[step]["maximum"] = max(float(by_step[step]["maximum"]), relative_error)
        grouped[eos] = [
            {"time": float(by_step[step]["time"]), "maximum": float(by_step[step]["maximum"])}
            for step in sorted(by_step)
        ]
    return grouped


def mass_balance_figure(
    rows: dict[str, list[dict[str, float]]],
    language: str,
    chinese_font: str,
) -> plt.Figure:
    fig, ax = new_figure(language, chinese_font)
    for eos in EOS_ORDER:
        ax.plot(
            [row["time"] for row in rows[eos]],
            [row["maximum"] for row in rows[eos]],
            color=EOS_COLORS[eos],
            linestyle=EOS_LINESTYLES[eos],
            linewidth=1.65,
            marker=EOS_MARKERS[eos],
            markersize=3.8,
            markerfacecolor="white",
            markeredgecolor=EOS_COLORS[eos],
            markeredgewidth=0.85,
            markevery=MARK_EVERY,
            label=EOS_LABELS[eos],
            zorder=3,
        )
    ax.set_yscale("symlog", base=10, linthresh=1.0e-16, linscale=0.6)
    ax.yaxis.set_major_formatter(
        FuncFormatter(
            lambda value, _position: (
                "0"
                if value == 0.0
                else rf"$10^{{{int(round(math.log10(abs(value))))}}}$"
            )
        )
    )
    ax.set_xlim(0.0, 60.0)
    ax.set_xticks(range(0, 61, 10))
    ax.set_xlabel(localize(language, "Time (day)", "时间（天）"))
    ax.set_ylabel(
        localize(
            language,
            "Maximum component mass-balance error",
            "最大组分质量守恒相对误差",
        )
    )
    legend = ax.legend(loc="lower right", borderpad=0.55, handlelength=2.8)
    style_legend(legend)
    return fig


def solver_steps_figure(
    rows: list[dict[str, str]], language: str, chinese_font: str
) -> plt.Figure:
    by_eos = {row["eos"].lower(): row for row in rows}
    accepted = [numeric(by_eos[eos], "accepted_internal_steps") for eos in EOS_ORDER]
    rejected = [numeric(by_eos[eos], "rejected_internal_steps") for eos in EOS_ORDER]
    x = list(range(len(EOS_ORDER)))

    fig, ax = new_figure(language, chinese_font)
    bars_accepted = ax.bar(
        x,
        accepted,
        width=0.56,
        color="#4477AA",
        edgecolor="#1F2937",
        linewidth=0.75,
        label="Accepted",
        zorder=3,
    )
    bars_rejected = ax.bar(
        x,
        rejected,
        width=0.56,
        bottom=accepted,
        color="#EE6677",
        edgecolor="#1F2937",
        linewidth=0.75,
        hatch="///",
        label="Rejected",
        zorder=3,
    )
    ax.set_xticks(x, [EOS_LABELS[eos] for eos in EOS_ORDER])
    ax.set_ylim(0.0, max(a + r for a, r in zip(accepted, rejected)) * 1.30)
    ax.set_xlabel(localize(language, "Equation of state", "状态方程"))
    ax.set_ylabel(localize(language, "Internal time-step attempts", "内部时间步尝试次数"))
    for accepted_bar, rejected_bar, accepted_value, rejected_value in zip(
        bars_accepted, bars_rejected, accepted, rejected
    ):
        ax.text(
            accepted_bar.get_x() + accepted_bar.get_width() / 2,
            accepted_value / 2,
            f"{int(accepted_value)}",
            ha="center",
            va="center",
            color="white",
            fontsize=8.0,
            fontweight="bold",
        )
        ax.text(
            rejected_bar.get_x() + rejected_bar.get_width() / 2,
            accepted_value + rejected_value / 2,
            f"{int(rejected_value)}",
            ha="center",
            va="center",
            color="#111827",
            fontsize=7.8,
            fontweight="bold",
        )
    legend = ax.legend(loc="upper center", ncols=2, borderpad=0.55)
    style_legend(legend)
    return fig


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    formats = tuple(part.strip().lower() for part in args.formats.split(",") if part.strip())
    unsupported = set(formats) - {"png", "pdf", "svg"}
    if not formats or unsupported:
        raise SystemExit(f"unsupported formats: {sorted(unsupported)}")
    if args.dpi < 150:
        raise SystemExit("DPI must be at least 150 for raster export")

    reservoir_path = COMPARISON / "eos_reservoir_history.csv"
    wells_path = COMPARISON / "eos_well_history.csv"
    solver_path = COMPARISON / "eos_solver_summary.csv"
    reservoir = grouped_by_eos(read_csv(reservoir_path))
    wells = grouped_by_eos(read_csv(wells_path))
    solver = read_csv(solver_path)
    validate_common_time_grid(reservoir, wells)
    mass_balance = load_mass_balance_history()
    chinese_font = register_chinese_font()

    languages = ("en", "zh") if args.language == "all" else (args.language,)
    outputs: list[Path] = []
    for language in languages:
        figures: tuple[tuple[str, plt.Figure], ...] = (
            (
                "01_average_pressure",
                reservoir_figure(
                    reservoir,
                    "pressure_avg",
                    "Average reservoir pressure (bar)",
                    "储层平均压力（bar）",
                    language,
                    chinese_font,
                ),
            ),
            (
                "02_average_oil_saturation",
                reservoir_figure(
                    reservoir,
                    "s_o_avg",
                    r"Average oil saturation, $S_o$",
                    r"平均油相饱和度，$S_o$",
                    language,
                    chinese_font,
                ),
            ),
            (
                "03_average_gas_saturation",
                reservoir_figure(
                    reservoir,
                    "s_g_avg",
                    r"Average gas saturation, $S_g$",
                    r"平均气相饱和度，$S_g$",
                    language,
                    chinese_font,
                ),
            ),
            (
                "04_average_water_saturation",
                reservoir_figure(
                    reservoir,
                    "s_w_avg",
                    r"Average water saturation, $S_w$",
                    r"平均水相饱和度，$S_w$",
                    language,
                    chinese_font,
                ),
            ),
            (
                "05_producer_total_rate",
                well_figure(
                    wells,
                    "PROD",
                    "q_total_surface",
                    lambda value: -value,
                    r"Production rate (m$^3$ s$^{-1}$)",
                    r"生产井总产量（m$^3$ s$^{-1}$）",
                    language,
                    chinese_font,
                ),
            ),
            (
                "06_injector_bhp",
                well_figure(
                    wells,
                    "CO2_INJ",
                    "bhp",
                    lambda value: value,
                    "Injector bottom-hole pressure (bar)",
                    "注入井井底压力（bar）",
                    language,
                    chinese_font,
                ),
            ),
            (
                "07_component_mass_balance",
                mass_balance_figure(mass_balance, language, chinese_font),
            ),
            (
                "08_time_step_attempts",
                solver_steps_figure(solver, language, chinese_font),
            ),
        )
        for name, figure in figures:
            outputs.extend(
                save_formats(
                    figure,
                    args.output_dir / language / name,
                    formats,
                    args.dpi,
                    args.force,
                )
            )

    source_paths = [
        reservoir_path,
        wells_path,
        solver_path,
        *(RESULTS / eos / "component_mass_balance.csv" for eos in EOS_ORDER),
    ]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "matplotlib_version": matplotlib.__version__,
        "figure_size_mm": list(FIGURE_SIZE_MM),
        "dpi": args.dpi,
        "formats": list(formats),
        "languages": list(languages),
        "source_files": {
            str(path.relative_to(ROOT)): sha256(path) for path in source_paths
        },
        "transformations": [
            "No smoothing, interpolation, resampling, or output-time omission.",
            "Production rate is plotted as positive withdrawal magnitude: -q_total_surface.",
            "Mass-balance trace is the maximum absolute relative error across five components at each output time.",
            "Mass-balance y-axis uses symlog with linthresh=1e-16 so exact zero remains visible.",
            "Time-step bars stack accepted and rejected nonlinear attempts from simulation_summary.csv.",
        ],
        "style": {
            "titles": "none",
            "grid": "none",
            "frame": "four-sided",
            "legend_language": "English",
            "palette": EOS_COLORS,
            "line_styles": EOS_LINESTYLES,
            "markers": EOS_MARKERS,
        },
        "outputs": {
            str(path.relative_to(args.output_dir)): sha256(path) for path in outputs
        },
    }
    manifest_path = args.output_dir / "figure_manifest.json"
    if manifest_path.exists() and not args.force:
        raise SystemExit(f"refusing to overwrite existing manifest: {manifest_path} (use --force)")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"validated common time grid: 121 states, 0--60 day")
    print(f"generated {len(outputs)} figure files in {args.output_dir}")
    print(manifest_path)


if __name__ == "__main__":
    main()
