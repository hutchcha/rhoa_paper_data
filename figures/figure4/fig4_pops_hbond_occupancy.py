#!/usr/bin/env python3
"""Figure 4 rework (H-bond panel): REPLACE the GDP / GTP-1 / GTP-2 H-bond time series with the
mean N-patch (residues 5/6/7) hydrogen-bond occupancy to POPS only, one bar per simulation,
averaged over the LAST 15 US of each 20-us trajectory. (sum_POPS is a per-frame integer count
0-3, so a histogram looks odd; the mean occupancy is the clean summary.) Error bars = block-
averaged SEM (5 x 3-us blocks), which accounts for time correlation. GTP-2 engages POPS about
twice as much as GDP / GTP-1.

Data (per-frame, per-lipid, N-patch donors 5/6/7): WT_results_data/{GDP,GTP,Resim}_hbonds.csv
(mirrors .../Orientation Dynamics/RhoA/strided_trajs/hbonds_out/).

  conda activate memb_analysis ; python fig4_pops_hbond_occupancy.py
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
for k, v in {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
             "mathtext.fontset": "custom", "mathtext.rm": "Arial",
             "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"}.items():
    matplotlib.rcParams[k] = v


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FILES = {"RhoA-GDP":   f"{HERE}/data/GDP_hbonds.csv",
         "RhoA-GTP-1": f"{HERE}/data/GTP_hbonds.csv",
         "RhoA-GTP-2": f"{HERE}/data/Resim_hbonds.csv"}
COL = {"RhoA-GDP": "#B279A2", "RhoA-GTP-1": "#72B7B2", "RhoA-GTP-2": "#E45756"}   # per-system (FIGURE_STYLE)
TOTAL_US, LAST_US = 20.0, 15.0   # 20 us trajectories; analyse only the final 15 us
FIG_MM = (73, 55)                # exact placed size (mm)
AXL, YTK, XTK, VL = 11, 10, 8.5, 8    # y-label, y-ticks (larger); x-ticks, value labels

# NB: the CSV "time_us" column is mis-scaled; the Fig-4 pipeline maps frames to time with
# np.linspace(0, 20, n), so we reproduce that and keep t >= 5 us (the last 15 us).
NBLOCKS = 5                  # last-15-us split into 5 x 3-us blocks for the error bar


def block_sem(s, nblocks=NBLOCKS):
    """SEM of the mean from contiguous-block averaging (accounts for time correlation)."""
    bm = [b.mean() for b in np.array_split(s, nblocks)]
    return float(np.std(bm, ddof=1) / np.sqrt(nblocks))


means, sems = {}, {}
for sim, f in FILES.items():
    s_all = pd.read_csv(f)["sum_POPS"].to_numpy()
    t = np.linspace(0, TOTAL_US, len(s_all))
    s = s_all[t >= TOTAL_US - LAST_US]           # last 15 us only
    means[sim] = s.mean()
    sems[sim] = block_sem(s)

sims = list(FILES)
x = np.arange(len(sims))
fig, ax = plt.subplots(figsize=(FIG_MM[0] / 25.4, FIG_MM[1] / 25.4), constrained_layout=True)
ax.bar(x, [means[s] for s in sims], 0.52, yerr=[sems[s] for s in sims],
       color=[COL[s] for s in sims], edgecolor="k", linewidth=0.6,
       error_kw=dict(elinewidth=1.1, capsize=3, ecolor="0.2"))
for xi, s in zip(x, sims):
    ax.text(xi, means[s] + sems[s] + 0.04, f"{means[s]:.2f}", ha="center", va="bottom", fontsize=VL)
ax.set_xticks(x); ax.set_xticklabels(sims, fontsize=XTK, rotation=15)
ax.set_ylabel("POPS H-Bond Occupancy", fontsize=AXL)
ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax.set_ylim(0, 2.05)
ax.tick_params(axis="y", labelsize=YTK)
style_axes(ax)

STEM = "fig4_pops_hbond_occupancy"
for e in ("png", "svg", "pdf"):           # no bbox_inches="tight": keep the canvas exactly 73x55 mm
    fig.savefig(os.path.join(HERE, f"{STEM}.{e}"), dpi=300)
plt.close(fig)
for sim in FILES:
    print(f"{sim:12s} mean POPS H-bond occupancy {means[sim]:.3f} +/- {sems[sim]:.3f} (block SEM, last {LAST_US:.0f} us)")
print(f"wrote {STEM}.png / .svg / .pdf")
