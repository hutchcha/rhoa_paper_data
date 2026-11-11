#!/usr/bin/env python3
"""Panel A (main figure): pre/cross/post sidechain-lipid H-bond occupancy for N-patch + effector
residues, ORGANIZED BY STRUCTURAL REGION with shaded backgrounds. Shows (i) both basic and acidic
effector-switch residues develop lipid H-bonds through the transition (unexpected acidic-lipid
electrostatics -> motivates co-binding), and (ii) relocation: N-patch basic-lipid interaction is
lost in the mutants (Ala) but HVR and Switch II R68/R70 persist/gain.

Sidechain lipid occupancy from region_hbonds/region_hbond_occupancy.csv (block mean±SEM).
Transitioning sims (GTP-2/2A/3A): pre/cross/post; refs (GTP-1/GDP): settled/post only.

  python fig_A_transition_states.py -> paper_figures/panels/fig_A_transition_states.{png,svg}
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
plt.rcParams["svg.fonttype"] = "none"    # editable text in Illustrator (SVG)
plt.rcParams["pdf.fonttype"] = 42        # editable text in PDF
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.linewidth"] = 2

def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")

OUT = "paper_figures/figure6"
AA1 = {"ARG": "R", "LYS": "K", "ASP": "D", "GLU": "E", "ALA": "A"}
STATE_COL = {"pre": "#4C78A8", "cross": "#F58518", "post": "#E45756"}
# residues organized by structural region (N->C); each region shaded in its color
REGION_GROUPS = [("N-patch", [5, 6, 7]),
                 ("SwI", [40]),
                 ("inter\nswitch", [51, 54]),
                 ("SwII", [65, 67, 68, 70, 76, 78]),
                 ("HVR", [183, 185, 187])]
REGION_COL = {"N-patch": "#8C564B", "SwI": "#C9A227", "inter\nswitch": "#57A79B",
              "SwII": "#8172B3", "HVR": "#4C72B0"}
TRANS_SIMS = ["GTP-2", "2A", "3A"]
ORDER = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
# main panel A = membrane-matched WT vs mutant only (2A dropped; GDP/GTP-1 -> supplement)
MAIN_SIMS = ["GTP-2", "3A"]
SUPP_SIMS = ["GDP", "GTP-1"]
SUPP_OUT = "paper_figures/supplementary"


def make(sims, outdir, stem, title=None, figw=None, figh=4.7, exact=False, csv=None):
    os.makedirs(outdir, exist_ok=True)
    d = pd.read_csv(csv or "transition_states/region_hbonds/region_hbond_occupancy.csv")
    d = d[d.part == "sc"]                    # sidechain = electrostatic interactions
    resids = [r for _, rs in REGION_GROUPS for r in rs]
    spans = []; i = 0                        # region span (start_idx, end_idx) for shading/labels
    for name, rs in REGION_GROUPS:
        spans.append((name, i, i + len(rs) - 1)); i += len(rs)

    fig, axes = plt.subplots(1, len(sims), figsize=(figw or 3.65 * len(sims), figh),
                             sharey=True, constrained_layout=True, squeeze=False)
    axes = axes[0]
    for ax, sim in zip(axes, sims):
        states = ["pre", "cross", "post"] if sim in TRANS_SIMS else ["post"]
        sd = d[d.sim == sim]
        x = np.arange(len(resids)); nb = len(states); w = 0.8 / nb
        for name, a, b in spans:
            ax.axvspan(a - 0.5, b + 0.5, color=REGION_COL[name], alpha=0.12, zorder=0)
            ax.text((a + b) / 2, 1.0, name, ha="center", va="top", fontsize=12.5, fontweight="bold",
                    color=REGION_COL[name], transform=ax.get_xaxis_transform())
        for k, st in enumerate(states):
            occ, sem = [], []
            for r in resids:
                row = sd[(sd.resid == r) & (sd.state == st)]
                occ.append(row["lipid"].iloc[0] if len(row) else np.nan)
                sem.append(row["lipid_sem"].iloc[0] if len(row) else np.nan)
            ax.bar(x + (k - (nb - 1) / 2) * w, occ, w, yerr=sem,
                   error_kw=dict(elinewidth=0.7, capsize=2, ecolor="0.35"),
                   color=STATE_COL[st], edgecolor="k", linewidth=0.4, zorder=2)
        labs = []                            # per-sim resname (R5/K6/K7 in WT, A5/A6/A7 in mutants)
        for r in resids:
            row = sd[sd.resid == r]
            rn = row["resname"].iloc[0] if len(row) else "?"
            labs.append(f"{AA1.get(rn, '?')}{r}")
        ax.set_xticks(x); ax.set_xticklabels(labs, rotation=45, fontsize=12.5)
        ax.tick_params(axis="y", labelsize=13)
        ax.set_title(f"RhoA-{sim}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                     fontsize=16.5, fontweight="bold", pad=6)
        style_axes(ax)
    axes[0].set_ylabel("Sidechain–Lipid H-Bond Occupancy", fontsize=15)
    axes[0].set_ylim(0, 1.05)

    # state legend up top, where the title was (only where >1 state is shown)
    if any(s in TRANS_SIMS for s in sims):
        state_h = [Patch(facecolor=STATE_COL[s], edgecolor="k", label=s.capitalize()) for s in ["pre", "cross", "post"]]
        # "outside" loc -> constrained_layout reserves space inside the canvas (no tight-bbox needed)
        fig.legend(handles=state_h, frameon=False, fontsize=13, ncol=3, loc="outside upper center",
                   title="Transition State", title_fontsize=12)
    if title:
        fig.suptitle(title, fontsize=13)
    # exact-size main panel (figw/figh already the target box) -> no tight-bbox trimming;
    # supplement (defaults) keeps tight bbox for a snug crop
    tight = None if exact else "tight"
    for e in ("png", "svg"):
        fig.savefig(f"{outdir}/{stem}.{e}", dpi=300, bbox_inches=tight)
    plt.close(fig)
    print(f"wrote {outdir}/{stem}.png")


def main():
    # main GTP-2/3A: pre/cross/post from the transition-containing analysis (the last-5 µs excerpt is
    # post-transition only, so pre/cross must come from the transition window — parent CSV)
    make(MAIN_SIMS, OUT, "fig_A_transition_states",
         figw=130 / 25.4, figh=85 / 25.4, exact=True)   # exact 130 x 85 mm box (matches Panel A width)
    # GDP/GTP-1 references are non-transitioning (post only) -> use the final-5 µs occupancy tally
    # no suptitle (described in caption); widen horizontally to fill the panel-A row
    make(SUPP_SIMS, SUPP_OUT, "supp_hbond_occupancy_nontransition", figw=11.0, figh=4.7,
         csv=os.path.join(os.path.dirname(os.path.abspath(__file__)), "rhoa_full_sampled_trajs",
                          "transition_states", "region_hbonds", "region_hbond_occupancy.csv"))


if __name__ == "__main__":
    main()
