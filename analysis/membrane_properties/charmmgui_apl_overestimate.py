#!/usr/bin/env python3
"""Why the two leaflets hold different lipid counts: CHARMM-GUI's build-time per-lipid areas
(especially POPC) are overestimates.

(A) CHARMM-GUI *starting* leaflet areas, SIMPLE membrane. Membrane Builder balances the leaflets
    by AREA: the all-POPC upper leaflet and the POPC/POPE/POPS lower leaflet are built to nearly
    the same area (Δ = -0.7%), so the all-POPC leaflet (POPC = largest-area lipid) needs ~10
    FEWER lipids (190 vs 200) to tile it.
(B) Same, COMPLEX membrane. Here the upper leaflet is POPC+PSM: the small PSM (55.4 Å²) dilutes
    the big POPC, so the upper leaflet's mean area matches the lower and the counts are nearly
    equal (202 vs 200) -- i.e. the big simple-membrane asymmetry comes specifically from the
    upper leaflet being PURE POPC.
(C) Per-lipid area: CHARMM-GUI build estimate vs experiment vs equilibrated simulation. CHARMM-GUI
    gives POPC 68.3 Å² -- above experiment (65.6) and well above the relaxed simulation (64.1) --
    so it over-weights POPC at build time; the lipid COUNT is locked in before the membrane
    relaxes, making the leaflet asymmetry a build-time artifact of that high POPC estimate.

  conda activate memb_analysis ; python charmmgui_apl_overestimate.py
"""
import os
import numpy as np
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


def plabel(ax, L):
    ax.annotate(L, xy=(-0.18, 1.03), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=16, fontweight="bold", annotation_clip=False)


# =============================================================================================
# CHARMM-GUI build-time per-lipid areas (Å²).  POPC 68.3 and PSM 55.4 confirmed EXACTLY from
# step3_size.inp of the complex RhoA and a separate KRAS build; POPE 59.2 / POPS 62.7 are the
# standard CHARMM-GUI/CHARMM36 table values; PLA18 56.3 is then fixed by the complex lower-leaflet
# area (50 each: 50*(68.3+59.2+62.7+56.3)=12325).  The simple leaflets balance to -0.7% and the
# complex to -1.5% with these values -- a consistency check.
CG_A = {"POPC": 68.3, "POPE": 59.2, "POPS": 62.7, "PSM": 55.4, "PLA18": 56.3}
LIPCOL = {"POPC": "#1F77B4", "POPE": "#FF7F0E", "POPS": "#2CA02C",
          "PLA18": "#D62728", "PSM": "#9467BD"}                       # membrane OI_PALETTE
LEG_ORDER = ["POPC", "PSM", "POPE", "POPS", "PLA18"]
LIPLAB = {"PLA18": "PLA18 (pl-PE)"}   # pretty legend labels (default = name)

# per-lipid comparison (Panel C):
#   Experiment: POPC 65.6 (Kucerka 2011, 310 K), POPE 56.6 (303 K), POPS 62.7 (Pan 2014 Soft
#     Matter, 25 C), PSM 61.9 (Doktorova 2020, fluid SDP).
#   Simulation (equilibrated): POPC = simple pure-POPC upper leaflet; POPE/POPS = simple lower
#     leaflet; PSM = complex upper leaflet.  From apl_per_type.csv.
LIPIDS_C = ["POPC", "POPE", "POPS", "PSM"]
CG  = {l: CG_A[l] for l in LIPIDS_C}
EXP = {"POPC": 65.6, "POPE": 56.6, "POPS": 62.7, "PSM": 61.9}
SIM = {"POPC": 64.1, "POPE": 60.8, "POPS": 62.8, "PSM": 58.8}
C_CG, C_EXP, C_SIM = "#D4A017", "#7F7F7F", "#2E2E2E"   # build estimate / experiment / simulation

# ---- leaflet compositions, as built (initial-system PDBs) ------------------------------------
SIMPLE  = dict(upper=[("POPC", 190)],
               lower=[("POPC", 100), ("POPE", 60), ("POPS", 40)],
               xlabels=["Upper\npure POPC", "Lower\nPC/PE/PS"])
COMPLEX = dict(upper=[("POPC", 102), ("PSM", 100)],
               lower=[("POPC", 50), ("POPE", 50), ("POPS", 50), ("PLA18", 50)],
               xlabels=["Upper\nPOPC/PSM", "Lower\nPC/PE/PS/plPE"])

def leaf(comp): return sum(CG_A[l] * n for l, n in comp), sum(n for _, n in comp)

TK, AXL, VL, LG = 12, 13, 10, 11


def draw_starting(ax, memb, letter, show_ylabel, ymax):
    (A_up, n_up), (A_lo, n_lo) = leaf(memb["upper"]), leaf(memb["lower"])
    for xi, comp in [(0, memb["upper"]), (1, memb["lower"])]:
        bottom = 0.0
        for lip, n in comp:
            h = CG_A[lip] * n
            ax.bar(xi, h, 0.62, bottom=bottom, color=LIPCOL[lip], edgecolor="k", linewidth=0.6)
            bottom += h
    for xi, A in [(0, A_up), (1, A_lo)]:
        ax.text(xi, A + ymax * 0.015, f"{A:,.0f} Å²", ha="center", va="bottom",
                fontsize=10.5, fontweight="bold")
    ax.axhline(A_up, ls="--", lw=1.1, color="0.35", zorder=0)
    ax.text(0.5, ymax * 0.955, f"Δ = {(A_lo/A_up - 1) * 100:+.1f}%", ha="center", va="top",
            fontsize=11, color="k")
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, ymax)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{memb['xlabels'][0]}\n({n_up})", f"{memb['xlabels'][1]}\n({n_lo})"],
                       fontsize=TK)
    if show_ylabel:
        ax.set_ylabel("CHARMM-GUI Starting\nLeaflet Area (Å²)", fontsize=AXL)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(labelsize=TK)
    plabel(ax, letter); style_axes(ax)
    return (A_up, n_up, A_lo, n_lo)


# ===============================================================================================
fig = plt.figure(figsize=(9.2, 7.4), constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92])
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1], sharey=axA)
axC = fig.add_subplot(gs[1, :])

YMAX = 1.22 * max(leaf(SIMPLE["upper"])[0], leaf(SIMPLE["lower"])[0],
                  leaf(COMPLEX["upper"])[0], leaf(COMPLEX["lower"])[0])
sA = draw_starting(axA, SIMPLE,  "A", True,  YMAX)
sB = draw_starting(axB, COMPLEX, "B", False, YMAX)

# shared lipid legend across the top row (stacked leaflet composition, panels A & B)
fig.legend(handles=[Patch(facecolor=LIPCOL[l], edgecolor="k", label=LIPLAB.get(l, l)) for l in LEG_ORDER],
           frameon=False, fontsize=10.5, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 1.005),
           columnspacing=1.3, handlelength=1.1, handletextpad=0.5)

# ---- (C) per-lipid APL: CHARMM-GUI vs experiment vs simulation -------------------------------
x = np.arange(len(LIPIDS_C)); w = 0.26
for k, (lab, dvals, col) in enumerate([("CHARMM-GUI (build)", CG, C_CG),
                                       ("Experiment", EXP, C_EXP), ("Simulation", SIM, C_SIM)]):
    vals = [dvals[l] for l in LIPIDS_C]
    axC.bar(x + (k - 1) * w, vals, w, color=col, edgecolor="k", linewidth=0.5, label=lab)
    for xi, v in zip(x + (k - 1) * w, vals):
        axC.text(xi, v + 0.15, f"{v:.1f}", ha="center", va="bottom", fontsize=VL, color="0.15")
axC.set_xticks(x); axC.set_xticklabels(LIPIDS_C, fontsize=TK + 1)
axC.set_ylim(54, 70.5)
axC.set_ylabel("Area Per Lipid (Å²)", fontsize=AXL)
axC.tick_params(labelsize=TK)
axC.legend(frameon=False, fontsize=LG, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.09),
           columnspacing=1.6, handlelength=1.2, borderpad=0.2)
plabel(axC, "C"); style_axes(axC)

HERE = os.path.dirname(os.path.abspath(__file__))
STEM = "charmmgui_apl_overestimate"
for e in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(HERE, f"{STEM}.{e}"), dpi=300, bbox_inches="tight")
plt.close(fig)
for nm, s in [("simple ", sA), ("complex", sB)]:
    au, nu, al, nl = s
    print(f"(A/B) {nm}  upper {au:,.0f} ({nu})  lower {al:,.0f} ({nl})  Δ {(al/au-1)*100:+.2f}%")
for l in LIPIDS_C:
    print(f"(C) {l:5s}  CG {CG[l]:.1f}  exp {EXP[l]:.1f}  sim {SIM[l]:.1f}  "
          f"(CG-exp {CG[l]-EXP[l]:+.1f}, CG-sim {CG[l]-SIM[l]:+.1f})")
print(f"wrote {STEM}.png / .svg / .pdf")
