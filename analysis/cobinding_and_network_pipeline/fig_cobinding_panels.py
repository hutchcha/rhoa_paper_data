#!/usr/bin/env python3
"""Co-binding main-figure panels from the frame-matched full-res data (19,500 frames/sim).

A: cooperativity enhancement matrix (WT) — basic x acidic residues, log enhancement, showing
   the N-patch (R5/K6/K7) -> Switch II D78 cooperative block.
C: the tether dose-response — Switch II (D76/D78) anionic-lipid co-binding occupancy, decomposed
   by the residue CLASS that shares the lipid (N-patch / HVR / switch-internal), WT -> 2A -> 3A.
   Shows the co-binding SOURCE switching from N-patch (WT) to HVR (3A) while the N-patch tether
   collapses 100 -> ~5 -> 0 %.

  python fig_cobinding_panels.py
Outputs -> ../paper_figures/panels/fig_cobind_A_enhancement.{png,svg}, fig_cobind_C_doseresponse.{png,svg}
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

HERE = os.path.dirname(os.path.abspath(__file__))
GAINS = os.path.join(HERE, "transition_states", "gains")
OUT = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "panels"))
NFRAMES = 19500
MINCO = 100          # co-frame reliability floor for the matrix

# residue classes (basic)
NPATCH = ["R5", "K6", "K7"]
HVR = ["R182", "R183", "K185", "K186", "K187"]
SWB = ["R68", "R70"]           # switch-internal basics
CLASS_OF = {**{r: "N-patch" for r in NPATCH}, **{r: "HVR" for r in HVR}, **{r: "switch" for r in SWB}}
CLASS_COL = {"N-patch": "#C44E52", "HVR": "#4C72B0", "switch": "#8172B3"}  # CVD-safe trio
SW2_ACID = ["D76", "D78"]


def load():
    return pd.read_csv(f"{GAINS}/lipid_cobinding_enhancement.csv")


def panel_A(d):
    w = d[(d.sim == "GTP-2") & (d.co_frames >= MINCO)].copy()
    basics = [b for b in NPATCH + SWB + HVR if b in set(w.basic)]
    acids = [a for a in ["E40", "D67", "D76", "D78", "E54"] if a in set(w.acidic)]
    M = np.full((len(basics), len(acids)), np.nan)
    for _, r in w.iterrows():
        if r.basic in basics and r.acidic in acids:
            M[basics.index(r.basic), acids.index(r.acidic)] = r.enhancement
    fig, ax = plt.subplots(figsize=(4.3, 4.6), constrained_layout=True)
    Mm = np.ma.masked_invalid(M)
    cmap = plt.cm.cividis.copy(); cmap.set_bad("#eeeeee")
    im = ax.imshow(Mm, cmap=cmap, norm=LogNorm(vmin=1, vmax=2000), aspect="auto")
    ax.set_xticks(range(len(acids))); ax.set_xticklabels(acids, fontsize=11)
    ax.set_yticks(range(len(basics))); ax.set_yticklabels(basics, fontsize=11)
    # annotate enhancement
    for i in range(len(basics)):
        for j in range(len(acids)):
            if not np.ma.is_masked(Mm[i, j]):
                v = M[i, j]
                ax.text(j, i, f"{v:.0f}" if v >= 10 else f"{v:.1f}", ha="center", va="center",
                        fontsize=8.5, color="white" if v > 8 else "#222")
    # bracket the N-patch rows
    npi = [basics.index(b) for b in NPATCH if b in basics]
    if npi:
        ax.add_patch(plt.Rectangle((-0.5, min(npi) - 0.5), len(acids), len(npi),
                                   fill=False, edgecolor="#C44E52", lw=2.2))
        ax.text(len(acids) - 0.4, np.mean(npi), "N-patch", color="#C44E52", fontsize=9,
                rotation=90, va="center", ha="left", fontweight="bold")
    for lab, col in [("Switch I", "#DFAF37"), ("Switch II", "#8172B3")]:
        pass
    ax.set_xlabel("acidic residue (effector switches)", fontsize=11)
    ax.set_ylabel("basic residue", fontsize=11)
    ax.set_title("Cooperative co-binding enhancement (WT)", fontsize=12, loc="left", fontweight="bold")
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label("enhancement  P(acidic|basic) / P(acidic|¬basic)", fontsize=9)
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/fig_cobind_A_enhancement.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_cobind_A_enhancement.png")


def panel_C(d):
    """Switch II D76/D78 co-binding occupancy by basic-residue class, per sim."""
    sims = ["GTP-2", "2A", "3A"]
    classes = ["N-patch", "HVR", "switch"]
    occ = {s: {c: 0.0 for c in classes} for s in sims}
    for s in sims:
        sub = d[(d.sim == s) & (d.acidic.isin(SW2_ACID))]
        for c in classes:
            # occupancy = max co-binding fraction among that class's residues (any of D76/D78)
            cl = sub[sub.basic.isin([b for b, cc in CLASS_OF.items() if cc == c])]
            occ[s][c] = (cl["co_frames"].max() / NFRAMES * 100) if len(cl) else 0.0
    fig, ax = plt.subplots(figsize=(5.0, 4.6), constrained_layout=True)
    x = np.arange(len(sims)); w = 0.26
    for k, c in enumerate(classes):
        vals = [occ[s][c] for s in sims]
        bars = ax.bar(x + (k - 1) * w, vals, w, color=CLASS_COL[c], edgecolor="k",
                      linewidth=0.4, label=f"{c} basics")
        for xi, v in zip(x + (k - 1) * w, vals):
            if v > 1:
                ax.text(xi, v + 1.5, f"{v:.0f}", ha="center", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels([f"RhoA-{s.replace('GTP-2','WT')}" for s in sims], fontsize=11)
    ax.set_ylabel("Switch II (D76/D78) anionic-lipid\nco-binding occupancy (%)", fontsize=11)
    ax.set_title("The lipid tether collapses and switches source", fontsize=12, loc="left", fontweight="bold")
    ax.set_ylim(0, 108)
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    ax.spines[["right", "top"]].set_visible(False)
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/fig_cobind_C_doseresponse.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_cobind_C_doseresponse.png")
    # print the numbers
    print("\nSwitch II co-binding occupancy (%) by class:")
    print(pd.DataFrame(occ).T.to_string(float_format=lambda v: f"{v:.1f}"))


def main():
    os.makedirs(OUT, exist_ok=True)
    d = load()
    panel_A(d)
    panel_C(d)


if __name__ == "__main__":
    main()
