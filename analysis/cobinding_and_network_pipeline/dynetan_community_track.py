#!/usr/bin/env python3
"""Cross-sim dynamical-community comparison as a SEQUENCE MEMBERSHIP TRACK (the clean way to show
community differences when the raw community indices are arbitrary per run).

Community indices from Louvain are not comparable across sims, so we colour by REFERENCE-MATCHED
identity: take one sim as reference (WT = GTP-2) and colour every residue by the reference
community that its own community overlaps most (max residue overlap). Same physical module -> same
colour in every sim, so splits/merges (e.g. the N-patch joining the switch module in the mutants)
are obvious, and each colour is labelled once from the reference.

One horizontal track per sim (residue number on x). Structural landmarks marked on top.

  python dynetan_community_track.py -> paper_figures/figure7/dynetan_community_track.{png,svg}
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
ORDER = ["GDP", "GTP-1", "GTP-2", "3A"]
REF = "GTP-2"
REGIONS = [("N-terminus", 4, 10, "#8C564B"), ("SwI", 28, 42, "#C9A227"), ("interSw", 43, 55, "#57A79B"),
           ("SwII", 56, 78, "#8172B3"), ("Rho ins.", 120, 137, "#9AA0A6"), ("HVR", 180, 189, "#4C72B0")]


def load():
    return {s: pd.read_csv(f"{D}/{s}_nodes.csv")[["resid", "community"]] for s in ORDER}


SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))


def make(plot_sims, outdir, stem):
    """Community track for `plot_sims`; colors are WT(GTP-2)-matched, computed over ALL sims in
    ORDER (which includes the GTP-2 reference) so subsets stay color-consistent."""
    os.makedirs(outdir, exist_ok=True)
    nd = load()
    lo = min(int(d.resid.min()) for d in nd.values())
    hi = max(int(d.resid.max()) for d in nd.values())
    resids = np.arange(lo, hi + 1)
    from dynetan_community_colors import community_colors
    _, _, per_sim, _ = community_colors(D, ORDER)         # colors over all (incl. REF)

    n = len(plot_sims)
    fig, ax = plt.subplots(figsize=(9.0, 0.75 * n + 1.1), constrained_layout=True)
    rowh = 0.8
    for yi, s in enumerate(plot_sims):
        y = n - 1 - yi
        rc_map = per_sim[s]
        for r in resids:
            if int(r) in rc_map:
                ax.add_patch(plt.Rectangle((r - 0.5, y - rowh / 2), 1, rowh,
                                           color=rc_map[int(r)], linewidth=0))
        ax.text(lo - 2, y, f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                ha="right", va="center", fontsize=13, fontweight="bold")

    ytop = n - 0.4
    for name, a, b, col in REGIONS:
        ax.plot([a, b], [ytop + 0.15, ytop + 0.15], color=col, lw=4, solid_capstyle="butt", clip_on=False)
        ax.text((a + b) / 2, ytop + 0.42, name, ha="center", va="bottom", fontsize=11,
                color=col, fontweight="bold", clip_on=False)
    ax.set_xlim(lo - 0.5, hi + 0.5); ax.set_ylim(-0.7, ytop + 0.9)
    ax.set_yticks([]); ax.set_xlabel("Residue Number", fontsize=13)
    ax.tick_params(axis="x", labelsize=11)
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    for e in ("png", "svg"):
        fig.savefig(f"{outdir}/{stem}.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir}/{stem}.png")


def main():
    make(["GTP-2", "3A"], OUT, "dynetan_community_track")            # main: WT vs mutant
    make(["GDP", "GTP-1"], SUPP, "supp_dynetan_community_track")     # supplement refs


if __name__ == "__main__":
    main()
