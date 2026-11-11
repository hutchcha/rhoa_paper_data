#!/usr/bin/env python3
"""Assemble Figure 6 (lipid co-binding) from the finished panels:
  A  transition-state sidechain-lipid H-bond occupancy by region (5 sims)
  B  structural snapshot of the WT co-binding cluster (R5/K6/K7 + D78 bridging one POPS)
  C  acidic<->amine-lipid co-binding frequency (membrane-composition effect)
  D  co-binding network across all 5 sims (fixed residue set)

Reads the per-panel PNGs in paper_figures/figure6/ and lays them out; the structural render
sits beside the frequency bars. Panel letters added here.

  python compose_figure6.py -> paper_figures/figure6/figure6_composite.{png,svg}
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"   # keep text as editable text (Illustrator), not paths
plt.rcParams["pdf.fonttype"] = 42       # TrueType -> editable text in PDF too
plt.rcParams["ps.fonttype"] = 42
HERE = os.path.dirname(os.path.abspath(__file__))
F6 = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure6"))
SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))


def aspect(png, base=None):
    im = plt.imread(f"{base or F6}/{png}")
    return im.shape[1] / im.shape[0]        # width / height


def put(ax, png, letter, anchor=None, fill=False):
    # fill=True stretches the image to fill the whole cell (aspect="auto") -> guaranteed full-width
    # + left alignment; used for the wide bar panel A where a small horizontal stretch is invisible.
    ax.imshow(plt.imread(f"{F6}/{png}"), aspect="auto" if fill else None)
    ax.axis("off")
    if anchor:
        ax.set_anchor(anchor)               # 'S' -> image sits at bottom of its cell
    ax.annotate(letter, xy=(0, 1), xytext=(-2, 6), xycoords="axes fraction",
                textcoords="offset points", fontsize=17, fontweight="bold", va="bottom", ha="left")


def main():
    # membrane-matched WT vs mutant (GTP-2, 3A). Layout:
    #   row1: A transition (3A, narrow) | B H-bond occupancy (wider)
    #   row2: C structure | D frequency
    #   row3: E network (full width)
    W = 13.0
    aT = aspect("fig_transition_tiltrot.png")
    aOcc = aspect("fig_A_transition_states.png")
    aS = aspect("fig_struct_cobinding.png")
    aF = aspect("fig_B_amine_frequency.png")
    aN = aspect("fig_C_cobinding_network.png")
    h1 = W / (aT + aOcc)                        # A|B share this height
    h2 = W / (aS + aF)                          # C|D share this height
    hN = W / aN
    H = h1 + h2 + hN + 1.0
    fig = plt.figure(figsize=(W, H))
    # two independent 2-col rows + one full-width row, via nested gridspecs
    outer = GridSpec(3, 1, figure=fig, height_ratios=[h1, h2, hN], hspace=0.14)
    r1 = outer[0].subgridspec(1, 2, width_ratios=[aT, aOcc], wspace=0.04)
    r2 = outer[1].subgridspec(1, 2, width_ratios=[aS, aF], wspace=0.04)
    put(fig.add_subplot(r1[0]), "fig_transition_tiltrot.png", "A", anchor="N")
    put(fig.add_subplot(r1[1]), "fig_A_transition_states.png", "B", anchor="N")
    put(fig.add_subplot(r2[0]), "fig_struct_cobinding.png", "C", anchor="S")
    put(fig.add_subplot(r2[1]), "fig_B_amine_frequency.png", "D", anchor="S")
    put(fig.add_subplot(outer[2]), "fig_C_cobinding_network.png", "E")
    for e in ("png", "svg"):
        fig.savefig(f"{F6}/figure6_composite.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {F6}/figure6_composite.png")


def build_supp():
    """Figure 6 supplement (non-transitioning refs GDP, GTP-1): H-bond occupancy (A),
    co-binding frequency (B) + co-binding network (C)."""
    def put_s(ax, png, letter, anchor=None):
        ax.imshow(plt.imread(f"{SUPP}/{png}")); ax.axis("off")
        if anchor: ax.set_anchor(anchor)
        ax.annotate(letter, xy=(0, 1), xytext=(-2, 6), xycoords="axes fraction",
                    textcoords="offset points", fontsize=16, fontweight="bold", va="bottom", ha="left")
    W = 12.0
    aA = aspect("supp_hbond_occupancy_nontransition.png", SUPP)
    aB = aspect("supp_cobinding_frequency_nontransition.png", SUPP)
    aC = aspect("supp_cobinding_network_nontransition.png", SUPP)
    hA = W / aA
    h2 = W / (aB + aC)
    H = hA + h2 + 0.5
    fig = plt.figure(figsize=(W, H))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[hA, h2], width_ratios=[aB, aC],
                  hspace=0.12, wspace=0.03)
    put_s(fig.add_subplot(gs[0, :]), "supp_hbond_occupancy_nontransition.png", "A")
    put_s(fig.add_subplot(gs[1, 0]), "supp_cobinding_frequency_nontransition.png", "B", anchor="S")
    put_s(fig.add_subplot(gs[1, 1]), "supp_cobinding_network_nontransition.png", "C", anchor="S")
    for e in ("png", "svg"):
        fig.savefig(f"{SUPP}/supp_S7_refs_cobinding.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP}/supp_S7_refs_cobinding.png")


if __name__ == "__main__":
    main()
    build_supp()
