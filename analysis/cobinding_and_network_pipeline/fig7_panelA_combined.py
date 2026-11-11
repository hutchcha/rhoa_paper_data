#!/usr/bin/env python3
"""Figure 7 Panel A -- COMBINED replot. Instead of 3 stacked rows (one effector each, with a crystal
render), all three effector distributions are overlaid on a SINGLE plot: Switch I Cα RMSD-to-effector
KDE for RhoA-GTP-2 (solid) and RhoA-GTP-3A (dashed). Color = effector (legend); line style = system.
Crystal-structure renders removed. Does NOT overwrite the original fig7_panelA_swI_competence.*.

  conda activate memb_analysis ; python fig7_panelA_combined.py
"""
import os
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.linewidth"] = 1.2


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(1.2); ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(width=1.0, direction="out")


HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, "transition_states", "switches_fullres", "effector_rmsd_fullres.npz")
F7 = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure7"))

XR = (0, 7.5)                                  # Switch I RMSD range (x from 0; tails beyond ~7.5 are ~0)
SIMS = [("GTP-2", "-", "RhoA-GTP-2"),          # (key, linestyle, display label)
        ("3A",    "--", "RhoA-GTP-3A")]
EFF_COL = ["#1F77B4", "#F58518", "#54A24B"]    # ROCK1 / PKN / anillin (no lipids in panel -> palette reuse OK)
EFF_SHORT = ["ROCK1", "PKN", "Anillin"]
AXL, TK, LG = 8, 9, 6               # ticks bumped up a little; sized for a 3.2 x 1.7 in panel

z = np.load(NPZ, allow_pickle=True)
eff_labels = list(z["eff_labels"])             # ['ROCK1 (1S1C)', 'PKN (1CXZ)', 'anillin (4XOI)']
grid = np.linspace(*XR, 400)

fig, ax = plt.subplots(figsize=(3.2, 1.7), constrained_layout=True)
for r in range(len(eff_labels)):
    for key, ls, _ in SIMS:
        kde = gaussian_kde(z[f"{key}|SwI"][r])(grid)
        ax.plot(grid, kde, color=EFF_COL[r], lw=1.3, ls=ls)

ax.set_xlim(*XR); ax.set_xticks([0, 2.5, 5, 7.5])
ax.set_xlabel("Effector RMSD", fontsize=AXL)
ax.set_ylabel("Density", fontsize=AXL)
ax.tick_params(labelsize=TK)
ax.yaxis.set_major_locator(plt.MaxNLocator(3))
style_axes(ax)

# two compact legends: effector (color) and system (line style)
eff_handles = [Line2D([0], [0], color=EFF_COL[r], lw=1.8, label=EFF_SHORT[r]) for r in range(len(EFF_SHORT))]
sys_handles = [Line2D([0], [0], color="0.25", lw=1.5, ls=ls, label=lab) for _, ls, lab in SIMS]
leg1 = ax.legend(handles=eff_handles, title="Effector", frameon=False, fontsize=LG,
                 title_fontsize=LG, loc="upper right", handlelength=1.3,
                 labelspacing=0.25, handletextpad=0.4, borderpad=0.15)
ax.add_artist(leg1)
ax.legend(handles=sys_handles, frameon=False, fontsize=LG, loc="upper right",
          bbox_to_anchor=(1.0, 0.52), handlelength=1.8, labelspacing=0.25,
          handletextpad=0.4, borderpad=0.15)

STEM = "fig7_panelA_swI_competence_combined"
os.makedirs(F7, exist_ok=True)
for e in ("png", "svg", "pdf"):           # no bbox_inches="tight": keep canvas exactly 3.2 x 1.7 in
    fig.savefig(os.path.join(F7, f"{STEM}.{e}"), dpi=300)
plt.close(fig)
print(f"effectors: {eff_labels}")
for key, _, lab in SIMS:
    med = [float(np.median(z[f'{key}|SwI'][r])) for r in range(len(eff_labels))]
    print(f"  {lab:12s} median SwI RMSD -> " + "  ".join(f"{EFF_SHORT[r]} {med[r]:.1f}" for r in range(len(med))))
print(f"wrote {F7}/{STEM}.png / .svg / .pdf")
