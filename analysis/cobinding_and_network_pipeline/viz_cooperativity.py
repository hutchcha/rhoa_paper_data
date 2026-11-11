#!/usr/bin/env python3
"""Bipartite network of basic->acidic COOPERATIVITY on shared POPS lipids.

Reads lipid_cobinding_enhancement.csv. Left nodes = basic residues (phosphate
binders), right nodes = acidic residues (amine binders). An edge means the pair
co-binds the same POPS more than chance; edge colour = enhancement
(P(acidic|basic)/P(acidic|no basic), log scale, red = cooperative >1, grey ~1,
blue = anti), edge width = co-binding frame count. One panel per sim.

  python viz_cooperativity.py [--min-enh 1.5] [--min-frames 100]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse
import numpy as np
import pandas as pd
import stable_state_config as C

OUT = "transition_states/gains"


def resid_of(lab):
    import re
    return int(re.match(r"[A-Za-z](\d+)", lab).group(1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-enh", type=float, default=1.5)
    ap.add_argument("--min-frames", type=int, default=100)
    ap.add_argument("--cutoff", type=int, default=25, help="read _min{cutoff}.csv (25 or 50)")
    ap.add_argument("--lipidset", default="", help="'' for POPS, or '_amine' for all amine lipids")
    args = ap.parse_args()
    df = pd.read_csv(f"{OUT}/lipid_cobinding_enhancement{args.lipidset}_min{args.cutoff}.csv")
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    
    C.apply_style()

    sims = [s for s in ["GDP", "GTP-1", "GTP-2", "2A", "3A"] if len(df) and s in set(df["sim"])]
    if not sims:
        print(f"no cooperative pairs to plot for lipidset '{args.lipidset}' cutoff {args.cutoff}"); return
    fig, axes = plt.subplots(1, len(sims), figsize=(4.6 * len(sims) + 2.0, 7.6), squeeze=False)
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)     # log10 enhancement
    import matplotlib.pyplot as _plt2; cmap = _plt2.get_cmap("RdBu_r")

    for ax, sim in zip(axes[0], sims):
        d = df[df.sim == sim].copy()
        d = d[d.co_frames >= args.min_frames]
        if not len(d):
            ax.axis("off"); ax.set_title(f"RhoA-{sim}\n(no pairs)"); continue
        basics = sorted(set(d["basic"]), key=resid_of)
        acids = sorted(set(d["acidic"]), key=resid_of)
        yb = {b: i for i, b in enumerate(np.linspace(0, max(len(basics), len(acids)) - 1, len(basics)))}
        by = {b: np.linspace(0, 10, len(basics))[i] for i, b in enumerate(basics)}
        ay = {a: np.linspace(0, 10, len(acids))[i] for i, a in enumerate(acids)}
        # edges
        for _, r in d.iterrows():
            enh = max(r.enhancement, 1e-3)
            col = cmap(norm(np.log10(enh)))
            lw = 0.7 + 3.5 * (r.co_frames / d.co_frames.max())
            ax.plot([0, 1], [by[r.basic], ay[r.acidic]], color=col, lw=lw,
                    alpha=0.85, zorder=1, solid_capstyle="round")
            # annotate strong ones
            if enh >= 3:
                ax.text(0.5, (by[r.basic] + ay[r.acidic]) / 2,
                        f"{enh:.0f}×", fontsize=6.5, ha="center", va="center",
                        color="0.15", zorder=3,
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))
        # nodes (labels padded well clear of the circles)
        for b in basics:
            ax.scatter(0, by[b], s=320, color="#3B6FB0", edgecolor="k", lw=0.5, zorder=2)
            ax.text(-0.20, by[b], b, ha="right", va="center", fontsize=10.5,
                    color="#1f4e79", fontweight="bold")
        for a in acids:
            ax.scatter(1, ay[a], s=320, color="#C0392B", edgecolor="k", lw=0.5, zorder=2)
            ax.text(1.20, ay[a], a, ha="left", va="center", fontsize=10.5,
                    color="#7b241c", fontweight="bold")
        ax.set_xlim(-0.7, 1.7); ax.set_ylim(-1, 11)
        ax.text(0, 10.6, "basic\n(→ phosphate /\nSer carboxylate)", ha="center", fontsize=8.5, color="#1f4e79", fontweight="bold")
        ax.text(1, 10.6, "acidic\n(→ amine)", ha="center", fontsize=9, color="#7b241c", fontweight="bold")
        ax.set_title(f"RhoA-{sim}", fontsize=13, fontweight="bold")
        ax.axis("off")

    # dedicated colorbar axis at the far right, clear of the last panel
    fig.subplots_adjust(left=0.02, right=0.90, top=0.86, bottom=0.04, wspace=0.05)
    import matplotlib.cm as _cm; sm = _cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cax = fig.add_axes([0.925, 0.28, 0.010, 0.44])
    cb = fig.colorbar(sm, cax=cax, ticks=[-1, 0, 1, 2, 3])
    cb.ax.set_yticklabels(["0.1×", "1×", "10×", "100×", "1000×"])
    cb.set_label("cooperativity enhancement  P(acidic|basic)/P(acidic|¬basic)", fontsize=9)
    lset = "amine lipids (POPS / POPE / PLA18)" if args.lipidset == "_amine" else "POPS lipids"
    fig.suptitle(f"Cooperative co-binding of single {lset}: basic residue recruits acidic partner\n"
                 "edge width = co-binding frames; red = cooperative, grey ≈ independent",
                 fontsize=14, fontweight="bold", y=0.97)
    cut = args.cutoff
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/cooperativity_network{args.lipidset}_min{cut}.{e}", dpi=160)
    plt.close(fig)
    print(f"wrote {OUT}/cooperativity_network{args.lipidset}_min{cut}.png")


if __name__ == "__main__":
    main()
