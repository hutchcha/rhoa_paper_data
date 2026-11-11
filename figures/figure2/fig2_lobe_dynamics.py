#!/usr/bin/env python3
"""Figure 2 rework: RhoA G-domain lobe dynamics, RhoA-GDP and RhoA-GTP-1 only.
Two metric columns: (A) Z_lb time evolution (GDP over GTP-1); (B) Z_lb1-vs-Z_lb2 density (GDP over
GTP-1). Rows = simulation. Sized 3.9 x 4.2 in; horizontal 0-1 density colorbar ABOVE column B.
(Protein snapshots added manually; panel C removed; GTP-2 -> Fig 3.)

  conda activate memb_analysis ; python fig2_lobe_dynamics.py
"""
import pickle
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.stats import gaussian_kde

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 1.1
for k, v in {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
             "mathtext.fontset": "custom", "mathtext.rm": "Arial",
             "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"}.items():
    matplotlib.rcParams[k] = v

OD = Path(__file__).resolve().parent / "data"
HERE = Path(__file__).resolve().parent
lobe = pickle.load(open(OD / "lobe_data.pkl", "rb"))

SIMS = ["RhoA-GDP", "RhoA-GTP"]                     # rows; GTP = GTP-1
HEAD = {"RhoA-GDP": "RhoA-GDP", "RhoA-GTP": "RhoA-GTP-1"}
C_LB1, C_LB2 = "#F58518", "#4C78A8"                 # Z_lb1 orange, Z_lb2 blue (lab lobe convention)
MA_WINDOW, TOTAL_US = 100, 20.0
TICK, AXIS, TITLE_FS, PANEL_FS, LEG_FS, CB_FS, CBT_FS = 11, 12, 11.5, 14, 10, 9, 8
ZLIM = (25, 65)


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(1.1); ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(width=0.9, direction="out", labelsize=TICK)


def line_plot(ax, x, y, color):
    y = np.asarray(y, float)
    ma = pd.Series(y).rolling(window=MA_WINDOW, center=True).mean()
    ax.plot(x, y, alpha=0.18, color=color, lw=0.4, zorder=1)
    ax.plot(x, ma, alpha=1.0, color=color, lw=1.2, zorder=2)


def draw_density(ax, df, nbins=60, levels=10, min_level=0.03):
    x, y = df["Lobe1"].to_numpy(), df["Lobe2"].to_numpy()
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    X, Y = np.meshgrid(np.linspace(x.min() - 5, x.max() + 5, nbins),
                       np.linspace(y.min() - 5, y.max() + 5, nbins))
    Z = gaussian_kde(np.vstack([x, y]))(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    Z = Z / np.nanmax(Z)
    Z = np.where(Z >= min_level, Z, np.nan)
    ax.contourf(X, Y, Z, levels=np.linspace(min_level, 1.0, levels), vmin=0, vmax=1, cmap="cividis")
    ax.set_xlim(*ZLIM); ax.set_ylim(*ZLIM)
    ax.set_xticks([30, 40, 50, 60]); ax.set_yticks([30, 40, 50, 60])
    style_axes(ax)


fig = plt.figure(figsize=(3.9, 4.2))
gs = fig.add_gridspec(2, 2, left=0.14, right=0.965, top=0.80, bottom=0.115, wspace=0.50, hspace=0.34)
axes = np.array([[fig.add_subplot(gs[r, c]) for c in range(2)] for r in range(2)])

for ri, s in enumerate(SIMS):
    # --- column A: Z_lb time evolution ---
    axT = axes[ri, 0]
    df = lobe[s]; t = np.linspace(0, TOTAL_US, len(df))
    line_plot(axT, t, df["Lobe1"], C_LB1)
    line_plot(axT, t, df["Lobe2"], C_LB2)
    axT.set_xlim(-0.4, TOTAL_US + 0.4); axT.set_xticks([0, 5, 10, 15, 20]); axT.set_ylim(*ZLIM)
    axT.set_ylabel("Z (Å)", fontsize=AXIS)
    axT.set_title(HEAD[s], fontsize=TITLE_FS, fontweight="bold", pad=2)
    style_axes(axT)
    if ri == 1:
        axT.set_xlabel("Time (µs)", fontsize=AXIS)
    else:
        axT.tick_params(labelbottom=False)
        axT.legend([Line2D([0], [0], color=C_LB1, lw=1.4), Line2D([0], [0], color=C_LB2, lw=1.4)],
                   [r"Z$_\mathrm{lb1}$", r"Z$_\mathrm{lb2}$"], frameon=False, fontsize=LEG_FS,
                   loc="upper center", ncol=2, handlelength=1.1, columnspacing=0.9,
                   borderaxespad=0.1, handletextpad=0.4)

    # --- column B: Z_lb1 vs Z_lb2 density ---
    axD = axes[ri, 1]
    draw_density(axD, df)
    axD.set_ylabel(r"Z$_\mathrm{lb2}$ (Å)", fontsize=AXIS)
    axD.set_title(HEAD[s], fontsize=TITLE_FS, fontweight="bold", pad=2)
    if ri == 1:
        axD.set_xlabel(r"Z$_\mathrm{lb1}$ (Å)", fontsize=AXIS)
    else:
        axD.tick_params(labelbottom=False)

# horizontal 0-1 colorbar ABOVE the top density panel (column B), clear of the titles
posB = axes[0, 1].get_position()
cax = fig.add_axes([posB.x0, 0.885, posB.width, 0.020])
cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap="cividis"), cax=cax, orientation="horizontal")
cb.set_ticks([0, 0.5, 1.0]); cb.ax.tick_params(labelsize=CBT_FS, pad=1, length=2)
cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
cb.set_label("Normalized Density", fontsize=CB_FS, labelpad=2)

# panel letters (left of each column, at the colorbar height)
fig.text(0.015, 0.915, "A", fontsize=PANEL_FS, fontweight="bold", va="top")
fig.text(posB.x0 - 0.075, 0.915, "B", fontsize=PANEL_FS, fontweight="bold", va="top")

for e in ("png", "svg"):
    fig.savefig(HERE / f"fig2_lobe_dynamics.{e}", dpi=600)
print("wrote fig2_lobe_dynamics.png (3.9 x 4.2 in)")
