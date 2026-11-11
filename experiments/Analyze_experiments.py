#!/usr/bin/env python3
"""
stacked_by_feed_2row_metrics.py

Expects files named like:
  phalloidin_starve.csv
  phalloidin_serum.csv
  pmlc2_starve.csv
  pmlc2_serum.csv

Each CSV is "wide":
  columns = conditions (e.g., NLS_bgal, WT, R5A/K6A/K7A, G14V, G14V-R5A/K6A/K7A)
  rows    = replicate values (blanks allowed)

Produces ONE figure per feed state (starve / serum) with 2 rows:
  Row 1 = phalloidin
  Row 2 = pmlc2

Each subplot:
- scatter (small points) + mean line + SEM error bar (GraphPad-ish)
- one-way ANOVA + Tukey HSD (all pairs)
- significance brackets ordered by leftmost condition
- Saves as PDF (with bbox_inches tight so brackets aren’t cut off)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})
try:
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
except ImportError as e:
    raise SystemExit("Missing statsmodels. Install: pip install statsmodels") from e


# =============================================================================
# CONFIG (edit these; no CLI args needed for sizing/fonts)
# =============================================================================
CONFIG = {
    "metrics": ["phalloidin", "pmlc2"],          # row order: top -> bottom
    "feed_states": ["starve", "serum"],          # makes one figure per feed state

    # Output naming
    "out_suffix": "stacked",                     # produces e.g. starve_stacked.pdf

    # Figure sizing/layout (tune to fit your manuscript panel)
    "figsize": (7, 4.0),                       # width, height inches for the *2-row* figure
    "hspace": 0.75,                              # space between rows

    # Typography
    "fontsize": {
        "axis_label": 8,
        "tick": 10,
        "panel_label": 10,                       # metric label ("phalloidin/RhoA", "pMLC2/RhoA")
        "sig": 12,                               # stars
    },
    "x_label_rewrites": {
        "WT": r"RhoA",                      # or just "RhoA"
        "R5A/K6A/K7A": r"RhoA$\mathrm{^{AAA}}$",
        "G14V": r"RhoA$\mathrm{^{G14V}}$",
        "G14V-R5A/K6A/K7A": r"RhoA$\mathrm{^{G14V/AAA}}$",
    },
    # Scatter
    "point_size": 7.0,
    "point_alpha": 0.95,
    "jitter": 0.14,
    "seed": 123,

    # Mean + SEM
    "mean_line_halfwidth": 0.18,
    "mean_line_lw": 1.2,
    "sem_lw": 1.0,
    "sem_capsize": 5,

    # Axes cosmetics
    "spine_w": 1.25,
    "tick_w": 1.25,
    "tick_l": 3.0,

    # Optional fixed y-zoom per metric (set ymax to e.g. 30.0, or None for auto)
    "ylims_by_metric": {
        "phalloidin": {"ymin": 0.0, "ymax": 30.0},
        "pmlc2": {"ymin": 0.0, "ymax": 30.0},
    },

    # Brackets: drawn ABOVE the axis; can extend outside axes
    "brackets": {
        "only_significant": True,                # set False to draw ns too (usually clutter)
        "base_pad": 0.06,                        # relative to visible y-range
        "step": 0.12,
        "h": 0.025,
        "lw": 0.5,
        "reserve_inside_axes": True,
        "n_reserved_levels": 6,     # max number of bracket rows you expect
        "top_pad_levels": 1,        # extra empty levels
    },

    # Savefig: prevent cropping of brackets above axes
    "savefig": {
        "bbox_inches": "tight",
        "pad_inches": 0.03,
    },

    # Metric -> y label
    "ylabel_by_metric": {
        "phalloidin": "phalloidin/RhoA",
        "pmlc2": "pMLC2/RhoA",
    },

    # Condition colors (mean line + SEM)
    "condition_colors": {
        "NLS_bgal": "#E69F00",              # orange
        "WT": "#0072B2",                    # blue
        "R5A/K6A/K7A": "#CC79A7",           # rose
        "G14V": "#7B61FF",                  # purple
        "G14V-R5A/K6A/K7A": "#0072B2",      # blue (same as WT)
    },

    # X-label formatting (display only; data keys unchanged)

    # Optional: force a consistent order across all plots (set None to use CSV order)
    # "condition_order": ["NLS_bgal","WT","R5A/K6A/K7A","G14V","G14V-R5A/K6A/K7A"],
    "condition_order": None,
}
# =============================================================================


# Matplotlib defaults for vector text
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["svg.fonttype"] = "none"


def style_axes(ax: plt.Axes):
    ax.grid(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines["left"].set_linewidth(CONFIG["spine_w"])
    ax.spines["bottom"].set_linewidth(CONFIG["spine_w"])
    ax.tick_params(
        axis="both",
        labelsize=CONFIG["fontsize"]["tick"],
        width=CONFIG["tick_w"],
        length=CONFIG["tick_l"],
        direction="out",
    )


def sem(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return 0.0
    return float(np.nanstd(x, ddof=1) / np.sqrt(len(x)))


def p_to_stars(p: float) -> str:
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def add_sig_bracket(ax: plt.Axes, x1: float, x2: float, y: float, h: float, text: str):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y],
            linewidth=CONFIG["brackets"]["lw"], color="black", clip_on=False)
    ax.text((x1 + x2) / 2, y -1, text,
            ha="center", va="bottom", fontsize=CONFIG["fontsize"]["sig"], color="black", clip_on=False)


def format_condition_label(label: str) -> str:
    return CONFIG["x_label_rewrites"].get(label, label)


def parse_metric_and_feed(filename: str) -> tuple[str | None, str | None]:
    stem = Path(filename).stem.lower()
    m = re.match(r"^(phalloidin|pmlc2)_(starve|serum)$", stem)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def load_wide_csv_to_long(csv_path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(csv_path)
    cols = list(df.columns)

    data = df.copy()
    for c in cols:
        data[c] = pd.to_numeric(data[c], errors="coerce")

    long_df = (data.melt(var_name="Condition", value_name="Value")
                  .dropna(subset=["Value"])
                  .reset_index(drop=True))
    return long_df, cols


def order_conditions(cols_from_csv: list[str]) -> list[str]:
    if CONFIG["condition_order"] is not None:
        base = [c for c in CONFIG["condition_order"] if c in cols_from_csv]
        extras = [c for c in cols_from_csv if c not in set(base)]
        return base + extras
    return cols_from_csv


def run_oneway_anova(long_df: pd.DataFrame, condition_order: list[str]) -> tuple[float, float]:
    groups = []
    for cond in condition_order:
        vals = long_df.loc[long_df["Condition"] == cond, "Value"].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            groups.append(vals)
    if len(groups) < 2:
        raise ValueError("Need at least 2 non-empty groups for ANOVA.")
    f_stat, p_val = stats.f_oneway(*groups)
    return float(f_stat), float(p_val)


def run_tukey(long_df: pd.DataFrame) -> pd.DataFrame:
    res = pairwise_tukeyhsd(endog=long_df["Value"].values,
                            groups=long_df["Condition"].values,
                            alpha=0.05)
    tbl = pd.DataFrame(res.summary().data[1:], columns=res.summary().data[0])
    tbl = tbl.rename(columns={"p-adj": "p_adj"})
    tbl["p_adj"] = tbl["p_adj"].astype(float)
    tbl["reject"] = tbl["reject"].astype(str).str.lower().isin(["true", "1", "yes"])
    tbl["stars"] = tbl["p_adj"].map(p_to_stars)
    return tbl


def plot_panel(ax: plt.Axes, *, metric: str, feed: str,
               long_df: pd.DataFrame, condition_order: list[str],
               show_xlabel: bool):
    rng = np.random.default_rng(CONFIG["seed"])
    ylabel = CONFIG["ylabel_by_metric"].get(metric, metric)

    x = np.arange(len(condition_order), dtype=float)

    # Scatter
    for i, cond in enumerate(condition_order):
        vals = long_df.loc[long_df["Condition"] == cond, "Value"].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        dx = rng.uniform(-CONFIG["jitter"], CONFIG["jitter"], size=len(vals))
        ax.scatter(np.full_like(vals, x[i], dtype=float) + dx, vals,
                   s=CONFIG["point_size"], color="black", alpha=CONFIG["point_alpha"], zorder=2)

    # Mean + SEM
    means = long_df.groupby("Condition")["Value"].mean()
    sems = long_df.groupby("Condition")["Value"].apply(lambda s: sem(s.to_numpy(dtype=float)))
    for i, cond in enumerate(condition_order):
        if cond not in means.index:
            continue
        m = float(means.loc[cond])
        se = float(sems.loc[cond]) if cond in sems.index else 0.0
        color = CONFIG["condition_colors"].get(cond, "black")

        ax.errorbar(x[i], m, yerr=se, fmt="none",
                    ecolor=color, elinewidth=CONFIG["sem_lw"],
                    capsize=CONFIG["sem_capsize"], capthick=CONFIG["sem_lw"], zorder=4)
        ax.hlines(m, x[i] - CONFIG["mean_line_halfwidth"], x[i] + CONFIG["mean_line_halfwidth"],
                  colors=color, linewidth=CONFIG["mean_line_lw"], zorder=5)

    # Labels
    ax.set_xticks(x)
    ax.set_xticklabels([format_condition_label(c) for c in condition_order], rotation=-45)
    if not show_xlabel:
        ax.set_xticklabels([])

    ax.set_ylabel(ylabel, fontsize=CONFIG["fontsize"]["axis_label"])
    style_axes(ax)

    # Panel label (metric) inside axis
    ax.text(0.01, 0.97, ylabel,
            transform=ax.transAxes, ha="left", va="top",
            fontsize=CONFIG["fontsize"]["panel_label"])

    # Y-lims (zoom) per metric
    yconf = CONFIG["ylims_by_metric"].get(metric, {})
    ymin_cfg = yconf.get("ymin", None)
    ymax_cfg = yconf.get("ymax", None)
    if ymin_cfg is not None or ymax_cfg is not None:
        ax.set_ylim(bottom=float(ymin_cfg) if ymin_cfg is not None else None,
                    top=float(ymax_cfg) if ymax_cfg is not None else None)

    # Stats + brackets
    f_stat, p_val = run_oneway_anova(long_df, condition_order)
    tukey = run_tukey(long_df)

    # Visible axis range (after ylim) controls bracket spacing
    ymin_plot, ymax_plot = ax.get_ylim()
    yrng = max(1e-12, ymax_plot - ymin_plot)

    base = ymax_plot + CONFIG["brackets"]["base_pad"] * yrng
    step = CONFIG["brackets"]["step"] * yrng
    h = CONFIG["brackets"]["h"] * yrng

    to_plot = tukey.copy()

    # Remove any comparison involving NLS_bgal
    # bgal_name = "NLS_bgal"
    # to_plot = to_plot[
    #     (to_plot["group1"] != bgal_name) &
    #     (to_plot["group2"] != bgal_name)
    # ].copy()

    if CONFIG["brackets"]["only_significant"]:
        to_plot = to_plot[to_plot["reject"]].copy()

    # Group brackets by leftmost condition
    cond_to_x = {c: i for i, c in enumerate(condition_order)}
    to_plot = to_plot[to_plot["group1"].isin(cond_to_x) & to_plot["group2"].isin(cond_to_x)].copy()
    if not to_plot.empty:
        def lr(row):
            a = cond_to_x[row["group1"]]
            b = cond_to_x[row["group2"]]
            return (min(a, b), max(a, b))
        lr_df = to_plot.apply(lr, axis=1, result_type="expand")
        to_plot["x_left"] = lr_df[0]
        to_plot["x_right"] = lr_df[1]
        to_plot["span"] = to_plot["x_right"] - to_plot["x_left"]
        to_plot = to_plot.sort_values(["x_left", "span", "x_right", "p_adj"],
                                      ascending=[True, True, True, True])

    k = 0
    for _, r in to_plot.iterrows():
        add_sig_bracket(ax, float(r["x_left"]), float(r["x_right"]),
                        y=base + k * step, h=h, text=str(r["stars"]))
        k += 1

    return {"anova_F": f_stat, "anova_p": p_val, "tukey": tukey}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", type=Path, default=Path("."), help="Directory containing the CSV files.")
    ap.add_argument("--outdir", type=Path, default=Path("figs_out"), help="Output directory.")
    args = ap.parse_args()

    indir = args.indir
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # Collect files
    files = list(indir.glob("*.csv"))
    table: dict[tuple[str, str], Path] = {}
    for f in files:
        metric, feed = parse_metric_and_feed(f.name)
        if metric and feed:
            table[(metric, feed)] = f

    # Build one figure per feed state: rows = metrics
    for feed in CONFIG["feed_states"]:
        needed = [(m, feed) for m in CONFIG["metrics"]]
        missing = [k for k in needed if k not in table]
        if missing:
            print(f"[SKIP] {feed}: missing files: {[f'{m}_{feed}.csv' for (m, _) in missing]}")
            continue

        # Decide a consistent condition order (use first metric’s CSV order unless forced)
        first_metric = CONFIG["metrics"][0]
        _, cols0 = load_wide_csv_to_long(table[(first_metric, feed)])
        cond_order = order_conditions(cols0)

        fig, axes = plt.subplots(
            nrows=len(CONFIG["metrics"]), ncols=1,
            figsize=(7,3.5),
            sharex=True,
        )

        if len(CONFIG["metrics"]) == 1:
            axes = [axes]

        stats_out = {}

        for i, metric in enumerate(CONFIG["metrics"]):
            csv_path = table[(metric, feed)]
            long_df, _ = load_wide_csv_to_long(csv_path)
            show_xlabel = (i == len(CONFIG["metrics"]) - 1)
            stats_dict = plot_panel(
                axes[i],
                metric=metric,
                feed=feed,
                long_df=long_df,
                condition_order=cond_order,
                show_xlabel=show_xlabel,
            )
            stats_out[(metric, feed)] = stats_dict

        # fig.tight_layout(pad=0.25)

        out_pdf = outdir / f"{feed}_{CONFIG['out_suffix']}.pdf"
        fig.savefig(out_pdf, **CONFIG["savefig"])
        plt.close(fig)
        print(f"[SAVED] {out_pdf}")

        # Save per-panel stats
        for metric in CONFIG["metrics"]:
            csv_path = table[(metric, feed)]
            stem = f"{metric}_{feed}"
            anova_txt = outdir / f"{stem}_anova.txt"
            posthoc_csv = outdir / f"{stem}_posthoc.csv"

            with open(anova_txt, "w", encoding="utf-8") as fh:
                fh.write(f"File: {csv_path.name}\n")
                fh.write("One-way ANOVA (scipy.stats.f_oneway)\n")
                fh.write(f"F = {stats_out[(metric, feed)]['anova_F']:.6g}\n")
                fh.write(f"p = {stats_out[(metric, feed)]['anova_p']:.6g}\n")

            stats_out[(metric, feed)]["tukey"].to_csv(posthoc_csv, index=False)

        print(f"[SAVED] stats for {feed}: (anova txt + tukey csv) per metric")


if __name__ == "__main__":
    main()