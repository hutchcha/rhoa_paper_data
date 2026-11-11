#!/usr/bin/env python3
"""Figure 3 rework: RhoA-GTP-2 dynamics and G-domain reorientation (4 data panels; render added later).
(A) core Cα RMSD time series; (B) geranylgeranyl-anchor insertion depth; (C) Z_lb1/Z_lb2 time series
(the lobes swap membrane proximity at the ~5 µs reorientation); (D) Z_lb1-vs-Z_lb2 density.

  conda activate memb_analysis ; python fig3_gtp2_dynamics.py
"""
import pickle
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 1.1
for k, v in {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
             "mathtext.fontset": "custom", "mathtext.rm": "Arial",
             "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"}.items():
    matplotlib.rcParams[k] = v

OD = Path(__file__).resolve().parent / "data"
HERE = Path(__file__).resolve().parent

# core Cα RMSD (excl. switch I/II) referenced to frame 0 -> starts at 0 (see measure_gtp2_core_rmsd.py)
rmsd = pickle.load(open(OD / "gtp2_core_rmsd_ref0.pkl", "rb"))["rmsd_A"]
# insertion depth: acyl-carbon anchor (CYSG C1-C20), center-referenced (measure_anchor_insertion_acyl.py)
anchor = pickle.load(open(OD / "anchor_insert_depth_acyl.pkl", "rb"))["RhoA-GTP-2"]
lobe = pickle.load(open(OD / "lobe_data.pkl", "rb"))["RhoA-Resim"]      # GTP-2 = Resim

C_LINE = "#1a1a1a"                      # RMSD / insertion -> black
C_LB1, C_LB2 = "#F58518", "#4C78A8"     # Z_lb1 orange, Z_lb2 blue (lab lobe convention)
MA_WINDOW, TOTAL_US = 100, 20.0          # 100 frames @ ~1 ns/frame -> ~100 ns (as original)
RMSD_MA_WINDOW = 200                      # RMSD here is 1 ns/frame (vs 2 ns in rmsdrhoa.pkl); 200 -> ~200 ns, matching the original RMSD smoothing
DENS_START_US = 5.0                      # density over the last 15 µs only (post-reorientation)
TICK, AXIS, PANEL_FS, LEG_FS = 13, 14, 14, 11
ZLIM = (25, 65)


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(1.1); ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(width=0.9, direction="out", labelsize=TICK)


def line_plot(ax, x, y, color, lw=1.2, win=MA_WINDOW):
    y = np.asarray(y, float)
    ma = pd.Series(y).rolling(window=win, center=True).mean()
    ax.plot(x, y, alpha=0.18, color=color, lw=0.4, zorder=1)
    ax.plot(x, ma, alpha=1.0, color=color, lw=lw, zorder=2)


def draw_density(ax, df, nbins=60, levels=10, min_level=0.03):
    x, y = df["Lobe1"].to_numpy(), df["Lobe2"].to_numpy()
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    X, Y = np.meshgrid(np.linspace(x.min() - 5, x.max() + 5, nbins),
                       np.linspace(y.min() - 5, y.max() + 5, nbins))
    Z = gaussian_kde(np.vstack([x, y]))(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    Z = Z / np.nanmax(Z); Z = np.where(Z >= min_level, Z, np.nan)
    ax.contourf(X, Y, Z, levels=np.linspace(min_level, 1.0, levels), vmin=0, vmax=1, cmap="cividis")
    ax.set_xlim(*ZLIM); ax.set_ylim(*ZLIM)
    ax.set_xticks([30, 40, 50, 60]); ax.set_yticks([30, 40, 50, 60]); style_axes(ax)


fig = plt.figure(figsize=(5.65, 4.59), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.03, h_pad=0.03, wspace=0.05, hspace=0.06)
gs = fig.add_gridspec(2, 2)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])

# (A) RMSD  (win=200 -> ~200 ns smoothing, matching the 2 ns/frame original RMSD)
line_plot(axA, np.linspace(0, TOTAL_US, len(rmsd)), rmsd, C_LINE, win=RMSD_MA_WINDOW)
axA.set_xlim(-0.4, TOTAL_US + 0.4); axA.set_xticks([0, 5, 10, 15, 20]); axA.set_ylim(bottom=0)
axA.set_xlabel("Time (µs)", fontsize=AXIS); axA.set_ylabel("RMSD (Å)", fontsize=AXIS); style_axes(axA)

# (B) insertion depth
line_plot(axB, np.linspace(0, TOTAL_US, len(anchor)), anchor, C_LINE)
axB.set_xlim(-0.4, TOTAL_US + 0.4); axB.set_xticks([0, 5, 10, 15, 20]); axB.set_ylim(bottom=0)
axB.set_xlabel("Time (µs)", fontsize=AXIS); axB.set_ylabel(r"$\mathrm{I_D}$ (Å)", fontsize=AXIS); style_axes(axB)

# (C) Z_lb time series
t = np.linspace(0, TOTAL_US, len(lobe))
line_plot(axC, t, lobe["Lobe1"], C_LB1); line_plot(axC, t, lobe["Lobe2"], C_LB2)
axC.set_xlim(-0.4, TOTAL_US + 0.4); axC.set_xticks([0, 5, 10, 15, 20]); axC.set_ylim(*ZLIM)
axC.set_xlabel("Time (µs)", fontsize=AXIS); axC.set_ylabel(r"Z$_\mathrm{lb}$ (Å)", fontsize=AXIS)
axC.legend([Line2D([0], [0], color=C_LB1, lw=1.4), Line2D([0], [0], color=C_LB2, lw=1.4)],
           [r"Z$_\mathrm{lb1}$", r"Z$_\mathrm{lb2}$"], frameon=False, fontsize=LEG_FS, loc="lower left",
           ncol=2, handlelength=1.1, columnspacing=0.9, handletextpad=0.4, borderaxespad=0.2)
style_axes(axC)

# (D) Z_lb1 vs Z_lb2 density -- last 15 µs only (post-reorientation)
draw_density(axD, lobe[t >= DENS_START_US])
axD.set_xlabel(r"Z$_\mathrm{lb1}$ (Å)", fontsize=AXIS); axD.set_ylabel(r"Z$_\mathrm{lb2}$ (Å)", fontsize=AXIS)

for ax, L in [(axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")]:
    ax.text(-0.32, 1.10, L, transform=ax.transAxes, fontsize=PANEL_FS, fontweight="bold", va="top")

for e in ("png", "svg"):
    fig.savefig(HERE / f"fig3_gtp2_dynamics.{e}", dpi=600)
print("wrote fig3_gtp2_dynamics.png")
