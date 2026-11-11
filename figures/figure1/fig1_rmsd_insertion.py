#!/usr/bin/env python3
"""Figure 1 rework (panels B, C): RhoA-GDP and RhoA-GTP-1 equilibration, 2 rows x 2 cols.
Row B = core Cα RMSD, row C = geranylgeranyl-anchor insertion depth; columns = RhoA-GDP / RhoA-GTP-1.
Raw trace (faint) + moving average (window 100). Sized to 69.96 x 77.29 mm for the main figure.

  conda activate memb_analysis ; python fig1_rmsd_insertion.py
"""
import pickle
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 1.3
for k, v in {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
             "mathtext.fontset": "custom", "mathtext.rm": "Arial",
             "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"}.items():
    matplotlib.rcParams[k] = v

OD = Path(__file__).resolve().parent / "data"
HERE = Path(__file__).resolve().parent
MM = 1 / 25.4
FIG_W, FIG_H = 69.96 * MM, 77.29 * MM
MA_WINDOW = 100
TOTAL_US = 20.0
TICK_FS, AXIS_FS, HEAD_FS, PANEL_FS = 10, 11, 11, 12

SIMS = ["RhoA-GDP", "RhoA-GTP"]                      # GTP = GTP-1
HEAD = {"RhoA-GDP": "RhoA-GDP", "RhoA-GTP": "RhoA-GTP-1"}
COL = {"RhoA-GDP": "#1a1a1a", "RhoA-GTP": "#1a1a1a"}  # both RMSD/I_D traces black

load = lambda p: pickle.load(open(OD / p, "rb"))
rmsd = load("rmsdrhoa.pkl")["rmsd"]
# insertion depth = membrane_COM_z - anchor_COM_z, anchor = geranylgeranyl ACYL carbons (CYSG C1-C20),
# not the whole CYSG residue (measure_anchor_insertion_acyl.py). center-referenced; ~7 A for all sims.
anchor = {k: v for k, v in pickle.load(open(OD / "anchor_insert_depth_acyl.pkl", "rb")).items()
          if k.startswith("RhoA")}


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(1.3); ax.spines["bottom"].set_linewidth(1.3)
    ax.tick_params(width=1.0, direction="out", labelsize=TICK_FS)


def line_plot(ax, x, y, color):
    y = np.asarray(y, float)
    ma = pd.Series(y).rolling(window=MA_WINDOW, center=True).mean()
    ax.plot(x, y, alpha=0.22, color=color, lw=0.5, zorder=1)
    ax.plot(x, ma, alpha=1.0, color=color, lw=1.2, zorder=2)


fig, axes = plt.subplots(2, 2, figsize=(FIG_W, FIG_H), sharex=True, constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.02, wspace=0.04, hspace=0.06)

# metric rows: (data dict, value-getter, y-label, panel-letter)
rmsd_get = lambda s: np.asarray(rmsd[s]["rmsd_A"], float)
anch_get = lambda s: np.asarray(anchor[s], float)
ROWS = [(rmsd_get, "RMSD (Å)", "B"), (anch_get, r"$\mathrm{I_D}$ (Å)", "C")]

for ri, (get, ylab, letter) in enumerate(ROWS):
    yvals = [get(s) for s in SIMS]
    ymax = max(np.nanmax(v) for v in yvals)
    ylo = -0.06 * ymax
    for ci, s in enumerate(SIMS):
        ax = axes[ri, ci]
        y = yvals[ci]
        line_plot(ax, np.linspace(0, TOTAL_US, len(y)), y, COL[s])
        ax.set_xlim(-0.4, TOTAL_US + 0.4); ax.set_xticks([0, 5, 10, 15, 20])
        ax.set_ylim(ylo, ymax * 1.06)
        style_axes(ax)
        if ri == 0:
            ax.set_title(HEAD[s], fontsize=HEAD_FS, fontweight="bold", pad=3)
        if ri == 1:
            ax.set_xlabel("Time (µs)", fontsize=AXIS_FS)
        if ci == 0:
            ax.set_ylabel(ylab, fontsize=AXIS_FS)
        else:
            ax.tick_params(labelleft=False)
    axes[ri, 0].text(-0.30, 1.17, letter, transform=axes[ri, 0].transAxes,
                     fontsize=PANEL_FS, fontweight="bold", va="top", ha="left")

for e in ("png", "svg"):
    fig.savefig(HERE / f"fig1_rmsd_insertion.{e}", dpi=600)
print(f"wrote fig1_rmsd_insertion.png  ({FIG_W/MM:.1f} x {FIG_H/MM:.1f} mm)")
