#!/usr/bin/env python3
"""Figure for the APL / PSM question (simple vs complex upper monolayer, protein-free).
(A) per-lipid-type area-per-lipid: POPC condenses when PSM is present, PSM ~= POPC.
(B) aggregate upper-monolayer area: PSM's size ALONE would not shrink the leaflet
    (dashed prediction ~= simple); only POPC condensation lowers the measured area.

  conda activate memb_analysis ; python fig_apl_psm_check.py
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
HERE = os.path.dirname(os.path.abspath(__file__))
dt = pd.read_csv(f"{HERE}/apl_per_type.csv")
da = pd.read_csv(f"{HERE}/apl_aggregate.csv")

SIMPLE, COMPLEX = "GTP-1", "GTP-2"          # representative WT pair (both GTP, 20 us)
def apl(sim, lip, leaf="upper"):
    r = dt[(dt.sim == sim) & (dt.leaflet == leaf) & (dt.lipid == lip)]
    return float(r.apl_mean.iloc[0]), float(r.apl_sem.iloc[0]), int(r.n.iloc[0])
def agg(sim):
    r = da[(da.sim == sim) & (da.leaflet == "upper")]
    return float(r.agg_area_mean.iloc[0]), float(r.agg_area_sem.iloc[0])

pc_s, pc_s_e, n_pc_s = apl(SIMPLE, "POPC", "upper")
pc_c, pc_c_e, n_pc_c = apl(COMPLEX, "POPC", "upper")
pc_s_lo, pc_s_lo_e, _ = apl(SIMPLE, "POPC", "lower")
pc_c_lo, pc_c_lo_e, _ = apl(COMPLEX, "POPC", "lower")
sm_c, sm_c_e, n_sm_c = apl(COMPLEX, "PSM", "upper")
A_s, A_s_e = agg(SIMPLE)
A_c, A_c_e = agg(COMPLEX)

# ---- EXPERIMENTAL area-per-lipid values (literature) -----------------------------------------
# POPC : 64.3 A^2 at 303 K (Kucerka, Nieh & Katsaras 2011, BBA 1808:2761); thermal expansion
#        ~0.19 A^2/K (alpha_A ~1.9e-3 nm^2/K) -> ~65.6 A^2 at 310 K (our sim temperature).
# PSM  : 61.9 A^2 fluid at 55 C (Doktorova et al. 2020, J Phys Chem B 124:5186). PSM Tm ~41 C, so
#        pure PSM is GEL at 310 K; the fluid value is the right reference for the mixed fluid leaflet.
# 50:50 POPC/PSM mixture : ~60.0 A^2 (vs pure POPC 64.7) -- Langmuir monolayer, Cardoso et al. 2024,
#        Molecules 29(24) -- experimental evidence of POPC/PSM condensation (monolayer, so indicative).
APL_POPC_EXP = 64.3 + 0.19 * (310 - 303)     # ~65.6 A^2 at 310 K
APL_PSM_EXP  = 61.9                           # fluid PSM, 55 C
APL_MIX_EXP  = 60.0                           # 50:50 POPC/PSM monolayer

A_ideal_exp = APL_POPC_EXP * n_pc_c + APL_PSM_EXP * n_sm_c   # experimental IDEAL mixing (no condensation)
A_exp_mix   = APL_MIX_EXP * (n_pc_c + n_sm_c)                # experimental mixture (condensed)

C_POPC, C_PSM, C_IDEAL, C_EXP = "#4C78A8", "#E45756", "0.72", "#54A24B"
H_LOWER = "///"                              # hatch = lower leaflet (upper = solid)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 4.5), constrained_layout=True)

# ---- Panel A: POPC packing by leaflet (upper has PSM in the complex), + PSM itself ----
# grouped by simulation; bar color = lipid, hatch = leaflet; PSM only exists in the complex upper
w = 0.26
# (x position, value, sem, color, hatch, value-label)
barsA = [
    (0.00, pc_s,    pc_s_e,    C_POPC, "",      f"{pc_s:.1f}"),     # GTP-1 POPC upper
    (0.28, pc_s_lo, pc_s_lo_e, C_POPC, H_LOWER, f"{pc_s_lo:.1f}"),  # GTP-1 POPC lower
    (0.92, pc_c,    pc_c_e,    C_POPC, "",      f"{pc_c:.1f}"),     # GTP-2 POPC upper
    (1.20, pc_c_lo, pc_c_lo_e, C_POPC, H_LOWER, f"{pc_c_lo:.1f}"),  # GTP-2 POPC lower
    (1.48, sm_c,    sm_c_e,    C_PSM,  "",      f"{sm_c:.1f}"),     # GTP-2 PSM upper
]
for xi, v, e, c, h, lab in barsA:
    axA.bar(xi, v, w, yerr=e, color=c, edgecolor="k", linewidth=0.6, hatch=h,
            error_kw=dict(elinewidth=0.9, capsize=3, ecolor="0.25"))
    axA.text(xi, v + 0.18, lab, ha="center", fontsize=10, fontweight="bold")
# experimental pure-component APL reference lines
axA.axhline(APL_POPC_EXP, ls="--", lw=1.2, color=C_POPC, alpha=0.9, zorder=0)
axA.axhline(APL_PSM_EXP, ls="--", lw=1.2, color=C_PSM, alpha=0.9, zorder=0)
axA.set_ylim(55, 69)
axA.set_xticks([0.14, 1.20]); axA.set_xticklabels(["RhoA-GTP-1\n(simple)", "RhoA-GTP-2\n(complex)"], fontsize=10.5)
axA.set_xlim(-0.3, 1.78)
axA.set_ylabel("area per lipid (Å²)", fontsize=12)
axA.set_title("A   POPC packing by leaflet", fontsize=12.5, loc="left", fontweight="bold")
axA.spines[["right", "top"]].set_visible(False)
axA.tick_params(labelsize=10.5)
legA = [Patch(facecolor=C_POPC, edgecolor="k", label="POPC"),
        Patch(facecolor=C_PSM, edgecolor="k", label="PSM"),
        Patch(facecolor="0.85", edgecolor="k", label="upper leaflet (+PSM)"),
        Patch(facecolor="0.85", edgecolor="k", hatch=H_LOWER, label="lower leaflet (no PSM)")]
axA.legend(handles=legA, frameon=False, fontsize=8.6, loc="upper right",
           handlelength=1.5, handleheight=1.1, labelspacing=0.35, borderpad=0.3)
# experimental reference-line labels, under the panel title (the dashed lines)
axA.text(0.015, 0.97, "experimental APL (dashed):", transform=axA.transAxes,
         fontsize=7.8, color="0.35", va="top", ha="left")
axA.text(0.015, 0.915, f"POPC  {APL_POPC_EXP:.1f} Å² (310 K)", transform=axA.transAxes,
         fontsize=8.4, color=C_POPC, va="top", ha="left", fontweight="bold")
axA.text(0.015, 0.86, f"PSM  {APL_PSM_EXP:.1f} Å² (fluid)", transform=axA.transAxes,
         fontsize=8.4, color=C_PSM, va="top", ha="left", fontweight="bold")

# ---- Panel B: simulation vs experimental POPC/PSM mixture area ----
ntot = n_pc_c + n_sm_c
blabs = ["experimental\nPOPC/PSM mixture", "simulation\n(RhoA-GTP-2,\nmeasured)"]
bvals = [A_exp_mix, A_c]; berrs = [0, A_c_e]
bcols = [C_EXP, C_POPC]
xb = np.arange(2)
bars = axB.bar(xb, bvals, 0.52, yerr=berrs, color=bcols, edgecolor="k", linewidth=0.6,
               error_kw=dict(elinewidth=0.9, capsize=3, ecolor="0.25"),
               hatch=["", "///"])
for xi, v in zip(xb, bvals):
    axB.text(xi, v + 8, f"{v:,.0f}", ha="center", fontsize=10.5, fontweight="bold")
    axB.text(xi, v - 42, f"{v/ntot:.1f} Å²/lipid", ha="center", fontsize=8.6, color="0.15",
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.8))
axB.set_xlim(-0.6, 1.6)
axB.set_ylim(11600, max(bvals) + 220)
axB.set_xticks(xb); axB.set_xticklabels(blabs, fontsize=10)
axB.set_ylabel("aggregate upper-monolayer area (Å²)", fontsize=12)
axB.set_title(f"B   Leaflet area: simulation vs experiment ({n_pc_c} POPC + {n_sm_c} PSM)",
              fontsize=11.5, loc="left", fontweight="bold")
axB.spines[["right", "top"]].set_visible(False)
axB.tick_params(labelsize=10.5)

for e in ("png", "svg"):
    fig.savefig(f"{HERE}/fig_apl_psm_check.{e}", dpi=250, bbox_inches="tight")
plt.close(fig)
print(f"POPC upper: simple={pc_s:.2f} complex={pc_c:.2f} (Δ{pc_c-pc_s:+.2f}); PSM upper={sm_c:.2f}")
print(f"POPC lower: simple={pc_s_lo:.2f} complex={pc_c_lo:.2f} (Δ{pc_c_lo-pc_s_lo:+.2f})")
print(f"complex intra-membrane POPC: lower {pc_c_lo:.2f} (no PSM) - upper {pc_c:.2f} (+PSM) = {pc_c_lo-pc_c:+.2f}")
print(f"EXP APLs: POPC(310K)={APL_POPC_EXP:.1f}  PSM(fluid)={APL_PSM_EXP:.1f}  mix(50:50)={APL_MIX_EXP:.1f}")
print(f"leaflet ({n_pc_c} POPC + {n_sm_c} PSM = {n_pc_c+n_sm_c}):")
print(f"  ideal mixing (exp)   = {A_ideal_exp:,.0f} A^2  ({A_ideal_exp/(n_pc_c+n_sm_c):.1f} A^2/lipid)")
print(f"  experimental mixture = {A_exp_mix:,.0f} A^2  ({APL_MIX_EXP:.1f} A^2/lipid)")
print(f"  simulation measured  = {A_c:,.0f} A^2  ({A_c/(n_pc_c+n_sm_c):.1f} A^2/lipid)")
print(f"  condensation (ideal->sim) = {A_c-A_ideal_exp:+,.0f} A^2 ({(A_c/A_ideal_exp-1)*100:+.1f}%)")
print("wrote fig_apl_psm_check.png / .svg")
