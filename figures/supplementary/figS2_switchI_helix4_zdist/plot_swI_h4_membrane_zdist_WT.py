#!/usr/bin/env python3
"""SI Figure S2 -- Time series of Switch I and Helix 4 distance from the membrane center
for the three wild-type systems (RhoA-GDP, RhoA-GTP-1, RhoA-GTP-2).

Distance to the membrane center = |z_region - z_mid|, where z_mid is the lipid (bilayer-COM)
midplane and z_region is the backbone COM of each region (Switch I resid 30-43, Helix 4
resid 141-151). The per-frame z-offsets were precomputed over the full 20-us WT trajectories
(see analysis/orientation_and_lobe/switchii_contacts.py) and the columns used here are bundled
in data/region_midplane_zmetrics_swI_helix4.csv. Stored d_a/d_b are signed z_region - z_mid
(negative, since the protein sits below the bilayer; distance = |d|). Stored time is 10x too
large, so real time = stored x 0.1 -> ~19.84 us.

Colors follow the lab lobe-COM convention (Switch I = Lobe 1 orange #F58518;
Helix 4 = Lobe 2 blue #4C78A8), matching main Figure 2.

  python plot_swI_h4_membrane_zdist_WT.py
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 1.2
for k, v in {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}.items():
    matplotlib.rcParams[k] = v

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "data", "region_midplane_zmetrics_swI_helix4.csv")

C_SWI, C_H4 = "#F58518", "#4C78A8"     # Switch I = Lobe 1 (orange); Helix 4 = Lobe 2 (blue)
TIME_SCALE = 0.1                        # stored time_us is 10x too large; real = stored x 0.1
MA_NS = 100.0                           # moving-average window in physical time (~100 ns)
SIMS = [("GDP", "RhoA-GDP"), ("GTP-1", "RhoA-GTP-1"), ("GTP-2", "RhoA-GTP-2")]
TK, AXL, LB = 11, 12, 12                   # tick / axis-label / in-panel-label font sizes


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(1.2); ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(width=1.0, length=4, pad=3, direction="out", labelsize=TK)


def ma(y, win):
    return pd.Series(y).rolling(window=max(1, win), center=True,
                                min_periods=max(1, win // 5)).mean().to_numpy()


def trace(ax, t, y, color):
    dt_ns = (t[-1] - t[0]) / (len(t) - 1) * 1000.0
    win = max(1, int(round(MA_NS / dt_ns)))
    ax.plot(t, y, alpha=0.20, color=color, lw=0.5, zorder=1)
    ax.plot(t, ma(y, win), alpha=1.0, color=color, lw=1.3, zorder=2)


df = pd.read_csv(CSV)
fig, axes = plt.subplots(len(SIMS), 1, figsize=(5.0, 6.0), sharex=True,
                         constrained_layout=True)

for ax, (key, label) in zip(axes, SIMS):
    t  = df[f"{key}_time_us_full"].to_numpy(float)
    da = df[f"{key}_d_a_full"].to_numpy(float)     # Switch I: z_region - z_mid (signed)
    db = df[f"{key}_d_b_full"].to_numpy(float)     # Helix 4
    m = np.isfinite(t) & np.isfinite(da) & np.isfinite(db)
    t = t[m] * TIME_SCALE
    trace(ax, t, np.abs(da[m]), C_SWI)             # distance to membrane center = |z_region - z_mid|
    trace(ax, t, np.abs(db[m]), C_H4)
    ax.set_xlim(0, 20); ax.set_xticks([0, 5, 10, 15, 20])
    ax.set_ylim(15, 80); ax.set_yticks([20, 40, 60])
    ax.text(0.015, 0.93, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=LB, fontweight="bold")
    style_axes(ax)

axes[-1].set_xlabel("Time (µs)", fontsize=AXL, labelpad=2)
fig.supylabel("Distance to Membrane Center (Å)", fontsize=AXL)

axes[0].legend(handles=[Line2D([0], [0], color=C_SWI, lw=1.6, label="Switch I"),
                        Line2D([0], [0], color=C_H4, lw=1.6, label="Helix 4")],
               frameon=False, fontsize=11, loc="upper right", ncol=2,
               handlelength=1.3, columnspacing=1.0, handletextpad=0.4, borderpad=0.15)

STEM = "swI_h4_membrane_zdist_WT"
for e in ("png", "svg", "pdf"):          # no bbox_inches="tight": keep canvas exactly 5.0 x 6.0 in
    fig.savefig(os.path.join(HERE, f"{STEM}.{e}"), dpi=600)
plt.close(fig)
for key, label in SIMS:
    da = np.abs(df[f"{key}_d_a_full"].dropna()); db = np.abs(df[f"{key}_d_b_full"].dropna())
    print(f"{label:12s} SwI {da.mean():.1f}  Helix4 {db.mean():.1f}  (A from membrane center)")
print(f"wrote {STEM}.png / .svg / .pdf  (5.0 x 6.0 in)")
