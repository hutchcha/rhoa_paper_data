#!/usr/bin/env python3
"""Quantitative dynetan allostery panel: (A) N-patch -> Switch I optimal-path length per sim
(lower = stronger MI communication; the channel strengthens with the mutation), and
(B) betweenness centrality of the salt-bridge residues D76/D78 (they become communication hubs
in the mutant). From the dynetan KSG-MI network outputs.

  python fig_dynetan_pathlength.py -> ../paper_figures/panels/fig_dynetan_pathlength.{png,svg}
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.linewidth"] = 2


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "transition_states", "dynetan")
OUT = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure7"))
ORDER = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
COL = {"GDP": "#B279A2", "GTP-1": "#72B7B2", "GTP-2": "#E45756", "2A": "#F58518", "3A": "#4C78A8"}


def path_length(sim):
    t = open(f"{D}/{sim}_path.txt").read()
    return float(t.split("length=")[1].split("\n")[0])


def path_str(sim):
    t = open(f"{D}/{sim}_path.txt").read()
    return t.split("path=")[1].strip()


SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))


def make(sims, outdir, path_stem, betw_stem, bracket=False):
    os.makedirs(outdir, exist_ok=True)
    x = np.arange(len(sims)); xt = [f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A") for s in sims]
    figw = 0.95 * len(sims) + 1.5          # narrower panel -> normal-thickness bars
    xpad = 0.7                              # whitespace on either side of the bar group

    figC, ax = plt.subplots(figsize=(figw, 4.8), constrained_layout=True)
    lens = [path_length(s) for s in sims]
    ax.bar(x, lens, width=0.55, color=[COL[s] for s in sims], edgecolor="k", linewidth=0.5)
    ax.set_xlim(-xpad, len(sims) - 1 + xpad)
    ax.set_xticks(x); ax.set_xticklabels(xt, fontsize=22, rotation=35, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="y", labelsize=24)
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax.set_ylabel("N-Patch→Switch I Path Length\nΣ −ln r$_{MI}$  (Lower = Stronger)", fontsize=17)
    ax.set_ylim(0, max(lens) * 1.2)
    style_axes(ax)
    for e in ("png", "svg"):
        figC.savefig(f"{outdir}/{path_stem}.{e}", dpi=200, bbox_inches="tight")
    plt.close(figC)

    figD, ax = plt.subplots(figsize=(figw, 4.8), constrained_layout=True)
    w = 0.3
    for k, res in enumerate([76, 78]):
        vals = [pd.read_csv(f"{D}/{s}_nodes.csv").set_index("resid").loc[res, "betweenness"] for s in sims]
        ax.bar(x + (k - 0.5) * w, vals, w, edgecolor="k", linewidth=0.5,
               color="#c0504d" if res == 76 else "#8064a2", label=f"D{res}")
    ax.set_xlim(-xpad, len(sims) - 1 + xpad)
    ax.set_xticks(x); ax.set_xticklabels(xt, fontsize=22, rotation=35, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="y", labelsize=24)
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax.set_ylabel("Betweenness Centrality", fontsize=24)
    ax.legend(frameon=False, fontsize=15, loc="upper left")
    style_axes(ax)
    for e in ("png", "svg"):
        figD.savefig(f"{outdir}/{betw_stem}.{e}", dpi=200, bbox_inches="tight")
    plt.close(figD)
    print(f"wrote {outdir}/{path_stem}.png + {betw_stem}.png")


def main():
    make(["GTP-2", "3A"], OUT, "fig_dynetan_pathlength", "fig_dynetan_betweenness", bracket=True)
    make(["GDP", "GTP-1"], SUPP, "supp_dynetan_pathlength_ref", "supp_dynetan_betweenness_ref")


if __name__ == "__main__":
    main()
