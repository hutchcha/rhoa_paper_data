#!/usr/bin/env python3
"""Panel C (main figure): lipid co-binding network across all 5 sims, using only CO-BINDABLE
lipids (POPS/POPE/PLA18 amine co-binding; POPC excluded). FIXED residue set + fixed node positions
in every panel so the networks are directly comparable — only the edges change. Residues with no
co-binding still show as nodes (no edge); weak edges are drawn. N-patch residues appear as Ala in
the mutants. Edge width = co-binding occupancy; nodes/edges colored by structural region.

Data: transition_states/gains/lipid_cobinding_enhancement_amine_min25.csv (full-res, frame-matched).

  python fig_cobinding_network.py -> ../paper_figures/panels/fig_C_cobinding_network.{png,svg}
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
NFR = 50000
SIMS = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
FIXED_BASIC = [5, 6, 7, 70, 183, 187]      # N-patch (R5/K6/K7) + SwII R70 + HVR R183/K187
FIXED_ACID = [54, 65, 67, 76, 78]          # interSw E54 + SwII D65/D67/D76/D78
REGION_COL = {"N-patch": "#8C564B", "SwI": "#C9A227", "interSw": "#57A79B",
              "SwII": "#8172B3", "HVR": "#4C72B0", "other": "#9AA0A6"}
AA1 = {"ARG": "R", "LYS": "K", "ASP": "D", "GLU": "E", "ALA": "A"}


def region_of(r):
    if 5 <= r <= 7: return "N-patch"
    if 28 <= r <= 42: return "SwI"
    if 43 <= r <= 55: return "interSw"
    if 56 <= r <= 78: return "SwII"
    if 180 <= r <= 190: return "HVR"
    return "other"


def resnum(x):
    return int("".join(c for c in str(x) if c.isdigit()))


def sim_resnames(sim):
    """per-sim resname for the fixed residues (Ala where the N-patch is mutated)."""
    u = mda.Universe(f"{HERE}/topologies/{sim}.psf")
    out = {}
    for r in FIXED_BASIC + FIXED_ACID:
        sel = u.select_atoms(f"protein and resid {r} and name CA")
        out[r] = AA1.get(sel.resnames[0], "?") if len(sel) else "?"
    return out


def make(sims, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    d = pd.read_csv(f"{HERE}/transition_states/gains/lipid_cobinding_enhancement_amine_min25.csv")
    d["b"] = d.basic.map(resnum); d["a"] = d.acidic.map(resnum)
    # fixed y-positions (same in every panel); compressed span -> tighter nodes, room for bigger labels
    yb = {r: (len(FIXED_BASIC) - 1 - i) / (len(FIXED_BASIC) - 1) * 7.0 + 1.3 for i, r in enumerate(FIXED_BASIC)}
    ya = {r: (len(FIXED_ACID) - 1 - i) / (len(FIXED_ACID) - 1) * 7.0 + 1.3 for i, r in enumerate(FIXED_ACID)}

    fig, axes = plt.subplots(1, len(sims), figsize=(3.6 * len(sims), 4.7), constrained_layout=True,
                             squeeze=False)   # wider + shorter (mockup ~1.29): "narrower height wise"
    axes = axes[0]
    for ax, sim in zip(axes, sims):
        rn = sim_resnames(sim)
        s = d[(d.sim == sim) & (d.b.isin(FIXED_BASIC)) & (d.a.isin(FIXED_ACID))]
        # edges (only between fixed residues; draw even weak ones)
        for _, r in s.iterrows():
            occ = r.co_frames / NFR
            reg = region_of(r.b)
            ax.plot([0.15, 0.85], [yb[r.b], ya[r.a]], color=REGION_COL[reg],
                    lw=0.5 + 8 * occ, alpha=0.8, solid_capstyle="round", zorder=1)
        # nodes (always all of them)
        for r in FIXED_BASIC:
            ax.scatter([0.15], [yb[r]], s=360, color=REGION_COL[region_of(r)], edgecolor="k",
                       linewidth=0.5, zorder=3)
            ax.text(0.075, yb[r], f"{rn[r]}{r}", ha="right", va="center", fontsize=14)
        for r in FIXED_ACID:
            ax.scatter([0.85], [ya[r]], s=360, color=REGION_COL[region_of(r)], edgecolor="k",
                       linewidth=0.5, zorder=3)
            ax.text(0.925, ya[r], f"{rn[r]}{r}", ha="left", va="center", fontsize=14)
        ax.set_xlim(-0.15, 1.15); ax.set_ylim(0.2, 9.9); ax.axis("off")
        ax.set_title(f"RhoA-{sim}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                     fontsize=15, fontweight="bold")
        ax.text(0.15, 9.2, "Basic", ha="center", fontsize=12, color="0.45")
        ax.text(0.85, 9.2, "Acidic", ha="center", fontsize=12, color="0.45")

    reg_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=REGION_COL[r], markeredgecolor="k",
                    markersize=13, label=r) for r in ["N-patch", "interSw", "SwII", "HVR"]]
    width_h = [Line2D([0], [0], color="0.5", lw=0.5 + 8 * f, label=f"{int(f*100)}% occ") for f in (0.1, 0.5, 1.0)]
    # stacked, centered legends (works for narrow 2-sim panels without overlap)
    leg1 = fig.legend(handles=reg_h, frameon=False, fontsize=14, ncol=4, loc="lower center",
                      bbox_to_anchor=(0.5, -0.10), title="Residue Region", title_fontsize=13.5,
                      handletextpad=0.4, columnspacing=1.4, markerscale=1.2)
    fig.add_artist(leg1)
    fig.legend(handles=width_h, frameon=False, fontsize=14, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.28), title="Edge = Co-Binding Occupancy", title_fontsize=13.5,
               handletextpad=0.7, columnspacing=1.8)
    for e in ("png", "svg"):
        fig.savefig(f"{outdir}/{stem}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir}/{stem}.png")


def main():
    make(["GTP-2", "3A"], OUT, "fig_C_cobinding_network")                       # main: WT vs mutant
    make(["GDP", "GTP-1"], SUPP, "supp_cobinding_network_nontransition")        # supplement


if __name__ == "__main__":
    main()
