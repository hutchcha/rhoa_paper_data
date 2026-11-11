#!/usr/bin/env python3
"""Two-panel cut of the combined SI membrane figure: just (A) POPC area-per-lipid by leaflet vs
experiment, and (C) bilayer thickness (P-P) -- Panel C restricted to the WT systems
RhoA-GDP / RhoA-GTP-1 / RhoA-GTP-2 (2A/3A dropped). Sized to 8.2 x 3.7 in.
Same data/colors as si_membrane_lipid_figure.py.

  conda activate memb_analysis ; python si_membrane_AC.py
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["hatch.linewidth"] = 1.0


def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")


def plabel(ax, L):
    ax.annotate(L, xy=(-0.2, 1.05), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=TS, fontweight="bold", annotation_clip=False)


ROOT = os.path.dirname(os.path.abspath(__file__))
APL = pd.read_csv(f"{ROOT}/data/apl_per_type.csv")
INS = pd.read_csv(f"{ROOT}/data/insertion_thickness.csv").set_index("sim")

def lab(s): return f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A")

# ---- data: APL (Panel A) ---------------------------------------------------------------------
def apl(sim, lip, leaf):
    r = APL[(APL.sim == sim) & (APL.leaflet == leaf) & (APL.lipid == lip)]
    return float(r.apl_mean.iloc[0]), float(r.apl_sem.iloc[0]), int(r.n.iloc[0])
pc_s, pc_s_e, _      = apl("GTP-1", "POPC", "upper")
pc_s_lo, pc_s_lo_e,_ = apl("GTP-1", "POPC", "lower")
pc_c, pc_c_e, _      = apl("GTP-2", "POPC", "upper")
pc_c_lo, pc_c_lo_e,_ = apl("GTP-2", "POPC", "lower")
sm_c, sm_c_e, _      = apl("GTP-2", "PSM", "upper")
APL_POPC_EXP, APL_PSM_EXP = 65.6, 61.9
C_POPC, C_PSM = "#1F77B4", "#9467BD"
H_LOWER, H_EXP = "///", "xx"

# ---- data: bilayer thickness (Panel C) -- WT systems only ------------------------------------
ORDER_C = ["GDP", "GTP-1", "GTP-2"]
LAB_C = [lab(s) for s in ORDER_C]
tv = np.array([INS.loc[s, "thick_PP_mean"] for s in ORDER_C])
te = np.array([INS.loc[s, "thick_PP_sd"] for s in ORDER_C])

# ===============================================================================================
fig, (axA, axC) = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True,
                               gridspec_kw=dict(width_ratios=[1.25, 1.0]))
TS, AXL, TK, VL = 16, 14, 12, 12

# ---- (A) POPC packing by leaflet vs experiment -----------------------------------------------
w = 0.17
barsA = [(0.00, pc_s,         pc_s_e,    C_POPC, ""),
         (0.18, pc_s_lo,      pc_s_lo_e, C_POPC, H_LOWER),
         (0.36, APL_POPC_EXP, 0,         C_POPC, H_EXP),
         (0.88, pc_c,         pc_c_e,    C_POPC, ""),
         (1.06, pc_c_lo,      pc_c_lo_e, C_POPC, H_LOWER),
         (1.24, APL_POPC_EXP, 0,         C_POPC, H_EXP),
         (1.42, sm_c,         sm_c_e,    C_PSM,  ""),
         (1.60, APL_PSM_EXP,  0,         C_PSM,  H_EXP)]
for xi, v, e, c, h in barsA:
    axA.bar(xi, v, w, yerr=e, color=c, edgecolor="k", linewidth=0.5, hatch=h,
            error_kw=dict(elinewidth=0.8, capsize=2, ecolor="0.25"))
axA.set_ylim(55, 68); axA.set_xlim(-0.22, 1.78)
axA.set_xticks([0.18, 1.24]); axA.set_xticklabels(["RhoA-GTP-1\nsimple", "RhoA-GTP-2\ncomplex"], fontsize=11)
axA.set_ylabel("Area Per Lipid (Å²)", fontsize=AXL)
plabel(axA, "A")
axA.legend(handles=[Patch(facecolor=C_POPC, edgecolor="k", label="POPC"),
                    Patch(facecolor=C_PSM, edgecolor="k", label="PSM"),
                    Patch(facecolor="0.85", edgecolor="k", label="upper leaflet"),
                    Patch(facecolor="0.85", edgecolor="k", hatch=H_LOWER, label="lower leaflet"),
                    Patch(facecolor="0.85", edgecolor="k", hatch=H_EXP, label="experiment")],
           frameon=False, fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3,
           handlelength=1.3, handleheight=1.0, labelspacing=0.3, columnspacing=1.0, borderpad=0.2)
style_axes(axA); axA.tick_params(labelsize=TK)

# ---- (C) bilayer thickness, P-P (WT only) ----------------------------------------------------
x = np.arange(len(ORDER_C))
axC.bar(x, tv, 0.55, yerr=te, color="#1F77B4", edgecolor="k", linewidth=0.5,
        error_kw=dict(elinewidth=0.8, capsize=2, ecolor="0.3"))
for xi, v in zip(x, tv):
    axC.text(xi, v + te[xi] + 0.4, f"{v:.1f}", ha="center", fontsize=VL - 0.5)
axC.set_xticks(x); axC.set_xticklabels(LAB_C, fontsize=TK, rotation=20, ha="right")
axC.set_ylabel("Bilayer Thickness, P–P (Å)", fontsize=AXL)
plabel(axC, "C")
axC.set_ylim(20, (tv + te).max() * 1.15)
style_axes(axC); axC.tick_params(labelsize=TK)

OUTDIR = os.path.dirname(os.path.abspath(__file__))
STEM = "si_membrane_AC"
for e in ("png", "svg", "pdf"):           # no bbox_inches="tight": keep canvas exactly 8.2 x 3.7 in
    fig.savefig(os.path.join(OUTDIR, f"{STEM}.{e}"), dpi=300)
plt.close(fig)
print(f"(C) thickness P-P: " + "  ".join(f"{s} {m:.2f}±{e:.2f}" for s, m, e in zip(ORDER_C, tv, te)))
print(f"wrote {STEM}.png / .svg / .pdf  (8.2 x 3.7 in)")
