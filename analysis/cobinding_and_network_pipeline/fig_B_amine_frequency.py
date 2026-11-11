#!/usr/bin/env python3
"""Panel B (main figure): lipid co-binding occurs in ALL simulations, and the complex
membrane (GTP-2/2A/3A) adds an extra amine-lipid partner, PLA18 (plasmalogen-PE) — which
the simple-membrane refs (GTP-1/GDP) lack. Frame-matched full-res data.

Co-binding frequency = fraction of (post) frames with an acidic-residue <-> amine-lipid
co-binding event, per lipid type. From transition_states/gains/cobinding_frequency_amine.csv.

  python fig_B_amine_frequency.py -> ../paper_figures/panels/fig_B_amine_frequency.{png,svg}
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"    # editable text in Illustrator (SVG)
plt.rcParams["pdf.fonttype"] = 42        # editable text in PDF
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.linewidth"] = 2

def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure6"))
SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))
LIPCOL = {"POPS": "#E45756", "POPE": "#4C78A8", "POPC": "#54A24B", "PLA18": "#B279A2"}  # co-binding palette (match cobinding_occupancy.py)
MAIN_SIMS = ["GTP-2", "3A"]          # membrane-matched WT vs mutant (complex membrane)
SUPP_SIMS = ["GDP", "GTP-1"]         # simple-membrane references


def make(sims, outdir, stem, memb_label, cats):
    os.makedirs(outdir, exist_ok=True)
    d = pd.read_csv(f"{HERE}/transition_states/gains/cobinding_frequency_amine.csv").set_index("sim")
    fig, ax = plt.subplots(figsize=(1.2 * len(sims) + 1.4, 5.4), constrained_layout=True)  # narrow/tall (mockup ~0.68)
    x = np.arange(len(sims)); w = 0.8 / len(cats)
    catcol = {"total": "#333333", **LIPCOL}
    for k, cat in enumerate(cats):
        vals = [d.loc[s, "frac" if cat == "total" else f"frac_{cat}"] * 100 for s in sims]
        lab = "Total (Any Lipid)" if cat == "total" else cat
        ax.bar(x + (k - (len(cats) - 1) / 2) * w, vals, w, color=catcol[cat],
               edgecolor="k", linewidth=0.4, label=lab)
    if memb_label:
        ax.text(0.5, 0.97, memb_label, transform=ax.transAxes, ha="center", va="top",
                fontsize=10.5, color="0.35")
    ax.set_xticks(x); ax.set_xticklabels(
        [f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A") for s in sims],
        fontsize=12, rotation=15)
    ax.tick_params(axis="y", labelsize=12)
    ax.set_ylabel("Acidic↔Amine-Lipid Co-Binding\nFrequency (% of Frames)", fontsize=13)
    ax.set_ylim(0, 105)
    ax.legend(frameon=False, fontsize=9.5, title="Co-Bound Lipid", title_fontsize=10,
              ncol=2 if len(cats) >= 4 else len(cats),   # narrow panel -> 2x2 legend instead of a wide row
              loc="lower center", bbox_to_anchor=(0.5, 1.005), columnspacing=1.0, handletextpad=0.5)
    style_axes(ax)
    for e in ("png", "svg"):
        fig.savefig(f"{outdir}/{stem}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir}/{stem}.png")


def main():
    make(MAIN_SIMS, OUT, "fig_B_amine_frequency", "inner leaflet: POPC / POPE / POPS / PLA18",
         ["total", "POPS", "POPE", "PLA18"])
    make(SUPP_SIMS, SUPP, "supp_cobinding_frequency_nontransition",
         "", ["total", "POPS", "POPE"])   # no in-plot membrane label on the supp panel


if __name__ == "__main__":
    main()
