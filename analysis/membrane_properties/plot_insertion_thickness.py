#!/usr/bin/env python3
"""Show that the anchor-insertion-depth difference between membranes is a thickness effect.
(A) bilayer thickness; (B) insertion depth measured from the bilayer CENTRE (suite metric) -- rises
with thickness; (C) insertion depth measured from the membrane SURFACE (phosphate plane) -- the
physically meaningful burial, essentially constant across membranes.

  conda activate memb_analysis ; python plot_insertion_thickness.py
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(f"{HERE}/insertion_thickness.csv").set_index("sim")

ORDER = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
def lab(s): return f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A")
LABELS = [lab(s) for s in ORDER]
C_SIMPLE, C_COMPLEX = "#4C78A8", "#E45756"
cols = [C_SIMPLE if d.loc[s, "membrane"] == "simple" else C_COMPLEX for s in ORDER]
x = np.arange(len(ORDER))

panels = [
    ("thick_PP",     "bilayer thickness, P–P (Å)",                 "A   Membrane thickness",              None),
    ("depth_center", "insertion depth from\nbilayer centre (Å)",   "B   Depth from centre (apparent)",    None),
    ("pen_surface",  "insertion depth from\nphosphate plane (Å)",  "C   Depth from surface (actual)",     None),
]
fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.4), constrained_layout=True)
for ax, (key, ylab, title, _) in zip(axes, panels):
    v = np.array([d.loc[s, f"{key}_mean"] for s in ORDER])
    e = np.array([d.loc[s, f"{key}_sd"] for s in ORDER])
    ax.bar(x, v, 0.68, yerr=e, color=cols, edgecolor="k", linewidth=0.6,
           error_kw=dict(elinewidth=0.9, capsize=3, ecolor="0.3"))
    for xi, val in zip(x, v):
        ax.text(xi, val + e[list(x).index(xi)] + (v.max() * 0.012), f"{val:.1f}",
                ha="center", fontsize=9.5, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(LABELS, fontsize=9, rotation=25, ha="right")
    ax.set_ylabel(ylab, fontsize=11.5)
    ax.set_title(title, fontsize=12, loc="left", fontweight="bold")
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.set_ylim(0, (v + e).max() * 1.15)

# shared membrane-type legend ABOVE the panels (clear of the bars)
fig.legend(handles=[Patch(facecolor=C_SIMPLE, edgecolor="k", label="simple membrane (POPC/POPE/POPS)"),
                    Patch(facecolor=C_COMPLEX, edgecolor="k", label="complex membrane (+PSM/PLA18)")],
           frameon=False, fontsize=10.5, ncol=2, loc="outside upper center",
           handlelength=1.4, columnspacing=2.0)
# make the "constant" point visually explicit on panel C: mean line across all sims
pc = np.array([d.loc[s, "pen_surface_mean"] for s in ORDER])
axes[2].axhline(pc.mean(), ls="--", lw=1.1, color="0.4", zorder=0)
axes[2].text(len(ORDER) - 1, pc.mean() + 0.3, f"mean {pc.mean():.1f} Å", ha="right", va="bottom",
             fontsize=8.5, color="0.4")

for e in ("png", "svg"):
    fig.savefig(f"{HERE}/fig_insertion_vs_thickness.{e}", dpi=250, bbox_inches="tight")
plt.close(fig)
print("wrote fig_insertion_vs_thickness.png / .svg")
