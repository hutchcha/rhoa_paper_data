#!/usr/bin/env python3
"""Former Figure 3, panel B: time evolution of RhoA rotation (φ) and tilt (θ) angles, stacked
(φ over θ), three WT systems overlaid (raw + 100-frame moving average). Replot of the original
time_ev_tiltrot.py, same data (OneDrive 'Tilt and Rotation' angle txt files) and colors
(GDP orange, GTP-1 blue, GTP-2 green), with:
  - x-ticks [0, 5, 10, 15, 20]
  - y-ticks every 45°, φ in [-180, 180], θ in [0, 180]
  - sized to 94 x 60 mm (exact canvas; legend added separately in Illustrator)

  conda activate memb_analysis ; python tiltrot_timeseries.py
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 1.2
for k, v in {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
             "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic"}.items():
    matplotlib.rcParams[k] = v

OD = Path(__file__).resolve().parent / "data"
HERE = Path(__file__).resolve().parent
TOTAL_US = 20.0
MA_NS = 100.0                 # moving-average smoothing window in PHYSICAL time (~100 ns)
TK, AXL = 8, 9

# (label, color, rotation_txt, tilt_txt)  -- same colors as the original time_ev_tiltrot.py
SYS = [("GDP",   "orange", "RhoA_GDProtationangles.txt",     "RhoA_GDPtiltangles.txt"),
       ("GTP-1", "blue",   "RhoA_GTProtationangles.txt",     "RhoA_GTPtiltangles.txt"),
       ("GTP-2", "green",  "RhoA_GTPredorotationangles.txt", "RhoA_GTPredotiltangles.txt")]


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(1.2); ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(width=1.0, length=3, pad=2, direction="out", labelsize=TK)
    ax.yaxis.set_label_coords(-0.1, 0.5)


def ma_window(n):
    """frames spanning ~MA_NS of real time for a trace of n frames over TOTAL_US."""
    return max(1, round(MA_NS * n / (TOTAL_US * 1000.0)))


def line_plot(ax, x, y, color):
    y = np.asarray(y, float)
    ma = pd.Series(y).rolling(window=ma_window(len(y)), center=True).mean()
    ax.plot(x, y, alpha=0.25, color=color, lw=0.6, zorder=1)
    ax.plot(x, ma, alpha=1.0, color=color, lw=1.2, zorder=2)


fig, (axR, axT) = plt.subplots(2, 1, figsize=(94 / 25.4, 60 / 25.4), sharex=True)
fig.subplots_adjust(left=0.135, right=0.955, top=0.965, bottom=0.185, hspace=0.2)

for label, color, rotf, tiltf in SYS:
    rot = np.loadtxt(OD / rotf)
    tilt = np.loadtxt(OD / tiltf)
    line_plot(axR, np.linspace(0, TOTAL_US, len(rot)), rot, color)
    line_plot(axT, np.linspace(0, TOTAL_US, len(tilt)), tilt, color)

axR.set_ylabel("φ (°)", fontsize=AXL)
axR.set_ylim(-180, 180); axR.set_yticks(np.arange(-180, 181, 90))   # rotation: ticks every 90\u00b0
axR.tick_params(labelbottom=False)                       # x labels only on bottom panel
axT.set_ylabel("θ (°)", fontsize=AXL)
axT.set_ylim(0, 180); axT.set_yticks(np.arange(0, 181, 45))
for ax in (axR, axT):
    ax.set_xlim(0, TOTAL_US); ax.set_xticks([0, 5, 10, 15, 20])
    style_axes(ax)
axR.tick_params(labelbottom=False)
axT.set_xlabel("Time (µs)", fontsize=AXL, labelpad=2)

STEM = "tiltrot_timeseries"
for e in ("png", "svg", "pdf"):                          # no bbox_inches="tight": keep canvas 94x60 mm
    fig.savefig(HERE / f"{STEM}.{e}", dpi=600)
plt.close(fig)
print(f"wrote {STEM}.png / .svg / .pdf  (94 x 60 mm)")
