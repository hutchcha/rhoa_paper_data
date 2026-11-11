#!/usr/bin/env python3
"""CHARMM-GUI *starting* (build-time) aggregate leaflet areas for the SIMPLE membrane.

Membrane Builder balances the two leaflets by cross-sectional AREA, not lipid count, then
fills each with as many lipids as that area holds. The simple-membrane upper leaflet is pure
POPC; the lower leaflet is a POPC/POPE/POPS mix. Because CHARMM-GUI gives POPC the largest
per-lipid area in the set (68.3 Å², confirmed from step3_size.inp of two independent builds),
the all-POPC upper leaflet needs ~10 FEWER lipids (190 vs 200) to tile the same area.

This figure shows the two starting leaflet areas side by side (lower leaflet stacked by lipid
contribution) to confirm the build was area-balanced.

  conda activate memb_analysis ; python charmmgui_simple_starting_area.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")


# ---- CHARMM-GUI build-time per-lipid areas (Å²) ----------------------------------------------
#   POPC confirmed exactly from step3_size.inp of the complex RhoA and a separate KRAS build
#   (102*68.3+100*55.4 = 12506.6 ; 83*68.3+80*55.4 = 10100.9).  POPE = CHARMM36 value (Klauda
#   2010, 59.2); POPS = 62.7 (CHARMM36). These reproduce the balanced leaflet areas below.
A_POPC, A_POPE, A_POPS = 68.3, 59.2, 62.7

# ---- SIMPLE-membrane leaflet composition, as built (initial-system PDB) -----------------------
UPPER = [("POPC", 190)]
LOWER = [("POPC", 100), ("POPE", 60), ("POPS", 40)]
APL = {"POPC": A_POPC, "POPE": A_POPE, "POPS": A_POPS}
COL = {"POPC": "#1F77B4", "POPE": "#FF7F0E", "POPS": "#2CA02C"}   # membrane OI_PALETTE

def area(comp): return sum(APL[l] * n for l, n in comp)
A_up, A_lo = area(UPPER), area(LOWER)
n_up, n_lo = sum(n for _, n in UPPER), sum(n for _, n in LOWER)

TK, AXL, VL, LG = 14, 15, 13, 12

fig, ax = plt.subplots(figsize=(4.6, 4.6), constrained_layout=True)

# upper leaflet: single pure-POPC block
ax.bar(0, A_up, 0.62, color=COL["POPC"], edgecolor="k", linewidth=0.6)
# lower leaflet: stacked by lipid contribution
bottom = 0.0
for lip, n in LOWER:
    h = APL[lip] * n
    ax.bar(1, h, 0.62, bottom=bottom, color=COL[lip], edgecolor="k", linewidth=0.6)
    bottom += h

# total-area labels above each bar
for x, A in [(0, A_up), (1, A_lo)]:
    ax.text(x, A + 230, f"{A:,.0f} Å²", ha="center", va="bottom", fontsize=VL, fontweight="bold")

# emphasise that the two starting areas match
ax.axhline(A_up, ls="--", lw=1.1, color="0.35", zorder=0)
ax.text(0.5, A_up - 350, f"Δ = {(A_lo/A_up - 1) * 100:+.1f}%", ha="center", va="top",
        fontsize=12, color="0.35")

ax.set_xlim(-0.65, 1.95); ax.set_ylim(0, max(A_up, A_lo) * 1.13)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"Upper\npure POPC\n({n_up} lipids)",
                    f"Lower\nPOPC/POPE/POPS\n({n_lo} lipids)"], fontsize=TK)
ax.set_ylabel("CHARMM-GUI Starting Leaflet Area (Å²)", fontsize=AXL)
ax.tick_params(labelsize=TK)
ax.legend(handles=[Patch(facecolor=COL[l], edgecolor="k", label=l) for l in ("POPC", "POPE", "POPS")],
          frameon=False, fontsize=LG, loc="center right", bbox_to_anchor=(1.0, 0.5),
          handlelength=1.2, labelspacing=0.4, borderpad=0.2)
style_axes(ax)

HERE = os.path.dirname(os.path.abspath(__file__))
STEM = "charmmgui_simple_starting_area"
for e in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(HERE, f"{STEM}.{e}"), dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"per-lipid APL (Å²): POPC {A_POPC}  POPE {A_POPE}  POPS {A_POPS}")
print(f"upper  {n_up} POPC                -> {A_up:,.1f} Å²  ({A_up/n_up:.1f} /lipid)")
print(f"lower  100 POPC / 60 POPE / 40 POPS -> {A_lo:,.1f} Å²  ({A_lo/n_lo:.1f} /lipid)")
print(f"starting-area imbalance: {(A_lo/A_up - 1) * 100:+.2f}%  ({A_lo - A_up:+.0f} Å²)")
print(f"wrote {STEM}.png / .svg / .pdf")
