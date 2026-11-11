#!/usr/bin/env python3
"""Combined SI figure (2 columns x 2 rows), sized for a Word page (6.5 in wide).
LEFT column  = lipid area-per-lipid (APL): (A) POPC packing by leaflet vs experiment,
               (B) complex-leaflet area, simulation vs experimental POPC/PSM mixture.
RIGHT column = membrane / anchor geometry: (C) bilayer thickness, (D) anchor insertion depth
               measured from the bilayer centre.

  conda activate memb_analysis ; python si_membrane_lipid_figure.py
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
matplotlib.rcParams["svg.fonttype"] = "none"   # editable text in Illustrator
matplotlib.rcParams["hatch.linewidth"] = 1.0   # crisper hatches (esp. for the PDF/vector export)

def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")

def plabel(ax, L):   # panel letter to the LEFT of the plot, freeing the top for a legend
    ax.annotate(L, xy=(-0.2, 1.05), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=TS, fontweight="bold", annotation_clip=False)
ROOT = os.path.dirname(os.path.abspath(__file__))
APL = pd.read_csv(f"{ROOT}/data/apl_per_type.csv")
AGG = pd.read_csv(f"{ROOT}/data/apl_aggregate.csv")
INS = pd.read_csv(f"{ROOT}/data/insertion_thickness.csv").set_index("sim")

def lab(s): return f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A")

# ---- data: APL -------------------------------------------------------------------------------
def apl(sim, lip, leaf):
    r = APL[(APL.sim == sim) & (APL.leaflet == leaf) & (APL.lipid == lip)]
    return float(r.apl_mean.iloc[0]), float(r.apl_sem.iloc[0]), int(r.n.iloc[0])
pc_s, pc_s_e, _      = apl("GTP-1", "POPC", "upper")
pc_s_lo, pc_s_lo_e,_ = apl("GTP-1", "POPC", "lower")
pc_c, pc_c_e, n_pc_c = apl("GTP-2", "POPC", "upper")
pc_c_lo, pc_c_lo_e,_ = apl("GTP-2", "POPC", "lower")
sm_c, sm_c_e, n_sm_c = apl("GTP-2", "PSM", "upper")
A_c = float(AGG[(AGG.sim == "GTP-2") & (AGG.leaflet == "upper")].agg_area_mean.iloc[0])
A_c_e = float(AGG[(AGG.sim == "GTP-2") & (AGG.leaflet == "upper")].agg_area_sem.iloc[0])
# experimental APLs (literature): POPC Kucerka 2011 (310 K); PSM Doktorova 2020 (fluid); mix Cardoso 2024
APL_POPC_EXP, APL_PSM_EXP = 65.6, 61.9      # experimental pure-component bilayer APLs (Kučerka; Doktorova)
ntot = n_pc_c + n_sm_c
A_ideal = APL_POPC_EXP * n_pc_c + APL_PSM_EXP * n_sm_c   # ideal, non-condensing mixing (no condensation)

C_POPC, C_PSM, C_EXP = "#1F77B4", "#9467BD", "#7F7F7F"
H_LOWER = "///"

# ---- data: insertion / thickness -------------------------------------------------------------
ORDER = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
C_SIMPLE, C_COMPLEX = "#1F77B4", "#D62728"
mcol = [C_SIMPLE if INS.loc[s, "membrane"] == "simple" else C_COMPLEX for s in ORDER]
LAB = [lab(s) for s in ORDER]

# ===============================================================================================
fig, ax = plt.subplots(2, 2, figsize=(7.0, 6.2), constrained_layout=True)
axA, axB = ax[0, 0], ax[0, 1]      # top row = lipid area-per-lipid
axC, axD = ax[1, 0], ax[1, 1]      # bottom row = membrane / anchor geometry
TS, AXL, TK, VL, LG = 16, 14, 12, 12, 12                # letter, axis, tick, value, legend sizes

# ---- (A) POPC packing by leaflet: measured (upper / lower) vs experimental bar per lipid -----
H_EXP = "xx"        # experimental bars = cross-hatch (vs // lower leaflet, solid upper)
w = 0.17
# (x, value, err, color, hatch)
barsA = [(0.00, pc_s,         pc_s_e,    C_POPC, ""),       # GTP-1 POPC upper
         (0.18, pc_s_lo,      pc_s_lo_e, C_POPC, H_LOWER),  # GTP-1 POPC lower
         (0.36, APL_POPC_EXP, 0,         C_POPC, H_EXP),    # experimental POPC
         (0.88, pc_c,         pc_c_e,    C_POPC, ""),       # GTP-2 POPC upper
         (1.06, pc_c_lo,      pc_c_lo_e, C_POPC, H_LOWER),  # GTP-2 POPC lower
         (1.24, APL_POPC_EXP, 0,         C_POPC, H_EXP),    # experimental POPC
         (1.42, sm_c,         sm_c_e,    C_PSM,  ""),       # GTP-2 PSM upper
         (1.60, APL_PSM_EXP,  0,         C_PSM,  H_EXP)]     # experimental PSM
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

# ---- (B) leaflet area: simulation vs experiment ----------------------------------------------
bvals = [A_ideal, A_c]; bcols = ["#1F77B4", "#1F77B4"]   # blue; ideal-mix vs sim told apart by hatch
axB.bar([0, 1], bvals, 0.5, yerr=[0, A_c_e], color=bcols, edgecolor="k", linewidth=0.5,
        error_kw=dict(elinewidth=0.8, capsize=2, ecolor="0.25"), hatch=["", H_LOWER])
for xi, v in zip([0, 1], bvals):
    axB.text(xi, v + 16, f"{v:,.0f}", ha="center", fontsize=VL)
axB.set_xlim(-0.6, 1.6); axB.set_ylim(11600, max(bvals) + 320)
axB.set_xticks([0, 1]); axB.set_xticklabels(["ideal\nmixing", "simulation"], fontsize=11)
axB.set_ylabel("Upper-Leaflet Area (Å²)", fontsize=AXL)
plabel(axB, "B")
style_axes(axB); axB.tick_params(labelsize=TK)

# ---- (C) bilayer thickness -------------------------------------------------------------------
x = np.arange(len(ORDER))
tv = np.array([INS.loc[s, "thick_PP_mean"] for s in ORDER]); te = np.array([INS.loc[s, "thick_PP_sd"] for s in ORDER])
axC.bar(x, tv, 0.68, yerr=te, color="#1F77B4", edgecolor="k", linewidth=0.5,
        error_kw=dict(elinewidth=0.8, capsize=2, ecolor="0.3"))
for xi, v in zip(x, tv):
    axC.text(xi, v + te[xi] + 0.4, f"{v:.1f}", ha="center", fontsize=VL - 0.5)
axC.set_xticks(x); axC.set_xticklabels(LAB, fontsize=TK, rotation=28, ha="right")
axC.set_ylabel("Bilayer Thickness, P–P (Å)", fontsize=AXL)
plabel(axC, "C")
axC.set_ylim(20, (tv + te).max() * 1.15)
style_axes(axC); axC.tick_params(labelsize=TK)

# ---- (D) insertion depth from the bilayer centre --------------------------------------------
C_CEN = "#1F77B4"                         # blue = from bilayer centre (center-referenced)
cen = np.array([INS.loc[s, "depth_center_mean"] for s in ORDER]); cen_e = np.array([INS.loc[s, "depth_center_sem"] for s in ORDER])
axD.bar(x, cen, 0.6, yerr=cen_e, color=C_CEN, edgecolor="k", linewidth=0.5,
        error_kw=dict(elinewidth=0.7, capsize=2, ecolor="0.3"))
axD.set_xticks(x); axD.set_xticklabels(LAB, fontsize=TK, rotation=28, ha="right")
axD.set_ylabel("Anchor Insertion Depth (Å)", fontsize=AXL)
plabel(axD, "D")
axD.set_ylim(0, (cen + cen_e).max() * 1.28)
style_axes(axD); axD.tick_params(labelsize=TK)

OUTDIR = os.path.dirname(os.path.abspath(__file__))
STEM = "si_membrane_lipid_figure"
# PDF keeps the bar hatches (SVG drops them); all three use Arial
for e in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(OUTDIR, f"{STEM}.{e}"), dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"POPC exp {APL_POPC_EXP}  PSM exp {APL_PSM_EXP}")
print(f"area sim {A_c:.0f} ({A_c/ntot:.1f})  ideal mix {A_ideal:.0f} ({A_ideal/ntot:.1f})  "
      f"condensation {(A_c/A_ideal-1)*100:+.1f}%")
print(f"centre depth mean {cen.mean():.2f} Å; thickness simple~{tv[:2].mean():.1f} complex~{tv[2:].mean():.1f}")
print("wrote si_membrane_lipid_figure.png / .svg")
