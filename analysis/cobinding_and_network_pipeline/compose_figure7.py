#!/usr/bin/env python3
"""Effector-competence block for Figure 7: one ROW per effector crystal (ROCK1/1S1C, PKN/1CXZ,
anillin/4XOI); each row = [switch RMSD KDE | core-superposed representative snapshot]. Built for a
single switch at a time:
  - Switch I  -> the MAIN Figure 7 Panel A (fig7_panelA_swI_competence)
  - Switch II -> a SEPARATE SUPPLEMENTARY figure (supp_effector_competence_swII)

Reads the per-frame RMSD npz and the PyMOL renders (render_<eff>_<SwI|SwII>.png).

  python compose_figure7.py
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.linewidth"] = 2


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, "transition_states", "switches_fullres", "effector_rmsd_fullres.npz")
SNAP = os.path.join(HERE, "transition_states", "effector_snapshots")
F7 = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure7"))
SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))
ORDER = ["GDP", "GTP-1", "GTP-2", "3A"]   # 2A dropped (KDE panels)
COL = {"GDP": "#B279A2", "GTP-1": "#72B7B2", "GTP-2": "#E45756", "2A": "#F58518", "3A": "#4C78A8"}


def build_block(sw, xr, swname, outdir, stem, title, sims=None, render_tag="main", big=False):
    sims = sims or ORDER
    os.makedirs(outdir, exist_ok=True)
    z = np.load(NPZ, allow_pickle=True)
    eff_labels = list(z["eff_labels"])
    short = [e.split(" ")[0] for e in eff_labels]

    fig = plt.figure(figsize=(9.5, 11))
    gs = GridSpec(3, 2, figure=fig, width_ratios=([0.68, 1.45] if big else [1.0, 1.25]),
                  hspace=0.62 if big else 0.26, wspace=-0.02)
    for r, (lab, sh) in enumerate(zip(eff_labels, short)):
        ax = fig.add_subplot(gs[r, 0])
        grid = np.linspace(*xr, 400)
        for n in sims:
            ax.plot(grid, gaussian_kde(z[f"{n}|{sw}"][r])(grid), color=COL[n], lw=2,
                    label=f"RhoA-{n}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"))
        # effector identity carried by the x-axis label (no separate row title)
        ax.set_xlabel((f"RMSD to {sh} (Å)" if big else f"{swname} Cα RMSD to {sh} (Å)"),
                      fontsize=20 if big else 14)
        ax.set_ylabel("Density", fontsize=26 if big else 14)
        ax.tick_params(labelsize=28 if big else 12)
        style_axes(ax)
        if big:
            # range-aware sparse x-ticks (SwII spans ~2-10, SwI spans ~2-18)
            if xr[1] <= 11:
                ax.set_xticks([2.5, 5, 7.5, 10])
            else:
                ax.set_xticks([5, 10, 15])
            ax.yaxis.set_major_locator(plt.MaxNLocator(3))
        if r == 0:
            ax.legend(frameon=False, fontsize=16 if big else 12, loc="upper right")
        axs = fig.add_subplot(gs[r, 1])
        png = f"{SNAP}/render_{sh}_{sw}_{render_tag}.png"
        if os.path.exists(png):
            axs.imshow(plt.imread(png))
        axs.axis("off"); axs.set_anchor("W")           # pull render toward the plot
    # crystal-reference legend intentionally omitted (author labels the gray trace manually)
    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold")
    for e in ("png", "svg"):
        fig.savefig(f"{outdir}/{stem}.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir}/{stem}.png")


def _img_aspect(p):
    im = plt.imread(p); return im.shape[1] / im.shape[0]


def _place(ax, png, letter, fs=12):
    ax.imshow(plt.imread(png)); ax.axis("off"); ax.set_anchor("N")
    ax.annotate(letter, xy=(0, 1), xytext=(0, 3), xycoords="axes fraction",
                textcoords="offset points", fontsize=fs, fontweight="bold", va="bottom", ha="left")


def build_full():
    """Figure 7 main (GTP-2 vs 3A). LEFT = Panel A (SwI effector competence, KDE incl. GDP/GTP-1
    refs + no-2A snapshots). RIGHT = Panel B (community track, WT vs mutant) over C (community
    3-D networks, WT vs mutant)."""
    A = f"{F7}/fig7_panelA_swI_competence.png"
    B = f"{F7}/dynetan_community_track.png"
    C = f"{F7}/dynetan_community_network.png"
    aA, aB, aC = (_img_aspect(p) for p in (A, B, C))
    W = 7.0
    K = 1.0 / aB + 1.0 / aC
    wR = W / (1 + K * aA); wA = W - wR
    H = wA / aA
    fig = plt.figure(figsize=(W, H))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[wA, wR],
                  height_ratios=[wR / aB, wR / aC], wspace=0.04, hspace=0.16)
    _place(fig.add_subplot(gs[:, 0]), A, "A", fs=11)
    _place(fig.add_subplot(gs[0, 1]), B, "B", fs=11)
    _place(fig.add_subplot(gs[1, 1]), C, "C", fs=11)
    for e in ("png", "svg"):
        fig.savefig(f"{F7}/figure7_composite.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {F7}/figure7_composite.png  ({W}in x {H:.1f}in)")


def build_S8():
    """Fig S8: GTP-2 vs 3A network/competence quantification. (A) Switch II effector competence
    (tall, left); (B) N-patch→Switch I optimal-path length and (C) D76/D78 betweenness (communication
    hubs), stacked at right."""
    A = f"{SUPP}/supp_effector_competence_swII.png"          # SwII, GTP-2/3A
    B = f"{F7}/fig_dynetan_pathlength.png"
    C = f"{F7}/fig_dynetan_betweenness.png"
    aA, aB, aC = (_img_aspect(p) for p in (A, B, C))
    W = 7.0
    K = 1.0 / aB + 1.0 / aC
    wR = W / (1 + K * aA); wA = W - wR
    H = wA / aA
    fig = plt.figure(figsize=(W, H))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[wA, wR],
                  height_ratios=[wR / aB, wR / aC], wspace=0.05, hspace=0.20)
    _place(fig.add_subplot(gs[:, 0]), A, "A")
    _place(fig.add_subplot(gs[0, 1]), B, "B")
    _place(fig.add_subplot(gs[1, 1]), C, "C")
    for e in ("png", "svg"):
        fig.savefig(f"{SUPP}/supp_S8_swII_network_quant.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP}/supp_S8_swII_network_quant.png")


def build_S9():
    """Fig S9: both switches' effector competence for the non-transitioning references (GDP, GTP-1),
    each row = [switch RMSD KDE | core-superposed snapshot], as in main Fig. 7A. (A) Switch I,
    (B) Switch II."""
    A = f"{SUPP}/supp_swI_competence_refs.png"               # SwI, GDP/GTP-1
    B = f"{SUPP}/supp_swII_competence_refs.png"              # SwII, GDP/GTP-1
    aA, aB = (_img_aspect(p) for p in (A, B))
    W = 9.5
    wcol = W / 2.0
    H = max(wcol / aA, wcol / aB) + 0.3
    fig = plt.figure(figsize=(W, H))
    gs = GridSpec(1, 2, figure=fig, wspace=0.05)
    _place(fig.add_subplot(gs[0]), A, "A")
    _place(fig.add_subplot(gs[1]), B, "B")
    for e in ("png", "svg"):
        fig.savefig(f"{SUPP}/supp_S9_refs_switch_competence.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP}/supp_S9_refs_switch_competence.png")


def build_S10():
    """Fig S10: full network breakdown for the non-transitioning references (GDP, GTP-1) —
    community track (A), 3-D community networks (B), path length (C), betweenness (D)."""
    A = f"{SUPP}/supp_dynetan_community_track.png"
    B = f"{SUPP}/supp_dynetan_community_network.png"
    C = f"{SUPP}/supp_dynetan_pathlength_ref.png"
    Dimg = f"{SUPP}/supp_dynetan_betweenness_ref.png"
    aA, aB, aC, aD = (_img_aspect(p) for p in (A, B, C, Dimg))
    W = 8.0
    hA = W / aA                                  # track full width
    hB = W / aB                                  # renders full width
    h3 = W / (aC + aD)                           # path | betweenness share this height
    H = hA + hB + h3 + 0.6
    fig = plt.figure(figsize=(W, H))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[hA, hB, h3], width_ratios=[aC, aD],
                  hspace=0.16, wspace=0.06)
    _place(fig.add_subplot(gs[0, :]), A, "A")
    _place(fig.add_subplot(gs[1, :]), B, "B")
    _place(fig.add_subplot(gs[2, 0]), C, "C")
    _place(fig.add_subplot(gs[2, 1]), Dimg, "D")
    for e in ("png", "svg"):
        fig.savefig(f"{SUPP}/supp_S10_refs_network.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP}/supp_S10_refs_network.png")


def build_closest_overlay():
    """Standalone aligned image (NOT part of any composite): the single closest (min-Cα-RMSD)
    Switch I frame of RhoA-GTP-2 and RhoA-3A overlaid on each effector crystal — 'how close does
    each get at its best'."""
    from matplotlib.patches import Patch
    effs = ["ROCK1", "PKN", "anillin"]
    pngs = [(e, f"{SNAP}/render_{e}_SwI_closest.png") for e in effs]
    pngs = [(e, p) for e, p in pngs if os.path.exists(p)]
    if not pngs:
        print("no closest-overlay renders found"); return
    fig, axes = plt.subplots(1, len(pngs), figsize=(3.2 * len(pngs), 3.5))
    axes = np.atleast_1d(axes)
    for ax, (e, p) in zip(axes, pngs):
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(e, fontsize=14)
    handles = [Patch(color=COL["GTP-2"], label="RhoA-GTP-2, closest"),
               Patch(color=COL["3A"], label="RhoA-GTP-3A, closest"),
               Patch(color=(0.3, 0.3, 0.3), label="crystal Switch I")]
    fig.legend(handles=handles, frameon=False, ncol=3, loc="lower center", fontsize=11.5)
    fig.suptitle("Closest sampled Switch I conformation (single minimum-RMSD frame)", fontsize=13)
    fig.tight_layout(rect=[0, 0.07, 1, 0.94])
    for ext in ("png", "svg"):
        fig.savefig(f"{F7}/fig7_closest_structures.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {F7}/fig7_closest_structures.png  (standalone; not in the composite)")


def main():
    # main Fig. 7A: Switch I competence, clean GTP-2 vs 3A comparison (references moved to supp)
    build_block("SwI", (2, 18), "Switch I", F7, "fig7_panelA_swI_competence", "",
                sims=["GTP-2", "3A"], render_tag="main")
    # component competence blocks for the supplements
    build_block("SwII", (2, 10), "Switch II", SUPP, "supp_effector_competence_swII", "",
                sims=["GTP-2", "3A"], render_tag="main", big=True)       # -> S8
    build_block("SwI", (2, 18), "Switch I", SUPP, "supp_swI_competence_refs", "",
                sims=["GDP", "GTP-1"], render_tag="refs", big=True)      # -> S9
    build_block("SwII", (2, 10), "Switch II", SUPP, "supp_swII_competence_refs", "",
                sims=["GDP", "GTP-1"], render_tag="refs", big=True)      # -> S9
    build_full()          # main Fig. 7
    build_S8()            # GTP-2/3A: SwII competence + path length + betweenness
    build_S9()            # GDP/GTP-1: both switches' competence
    build_S10()           # GDP/GTP-1: full network breakdown
    build_closest_overlay()   # standalone: closest GTP-2/3A Switch I on each crystal


if __name__ == "__main__":
    main()
