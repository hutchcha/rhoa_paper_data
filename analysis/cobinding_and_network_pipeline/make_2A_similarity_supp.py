#!/usr/bin/env python3
"""2A control supplement: RhoA-GTP-2A is globally similar in structure to the other simulations.
Core G-domain Cα RMSD of every sim (frame-matched last-1.95 µs excerpts, stride 10) to a common
reference (RhoA-GTP-2, first frame), after optimal core superposition. 2A's distribution overlaps
the others -> the double mutant is globally folded like the rest; it is used only as a control and
is otherwise omitted from the main figures.

  conda activate memb_analysis
  python make_2A_similarity_supp.py -> paper_figures/supplementary/supp_2A_global_similarity.{png,svg}
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))
# Fig S6: the three complex-membrane systems only (WT, and the two mutants) — are they all folded/stable?
ORDER = ["GTP-2", "2A", "3A"]
COL = {"GDP": "#B279A2", "GTP-1": "#72B7B2", "GTP-2": "#E45756", "2A": "#F58518", "3A": "#4C78A8"}
# G-domain core (exclude the mobile switches and Rho insert) -> reports fold, not loop motion
CORE = "protein and name CA and resid 4-179 and not (resid 28-42 or resid 56-78 or resid 120-137)"
STRIDE = 10


def core_coords(u, resids):
    ag = u.select_atoms("protein and name CA and resid " + " ".join(map(str, resids)))
    o = np.argsort(ag.resids)
    return ag.resids[o], ag.positions[o].astype(np.float64)


def main():
    os.makedirs(SUPP, exist_ok=True)
    uref = mda.Universe(f"{HERE}/topologies/GTP-2.psf", f"{HERE}/GTP-2_last5us.dcd")
    ref_resids = np.sort(uref.select_atoms(CORE).resids)
    uref.trajectory[0]
    _, ref = core_coords(uref, ref_resids); refc = ref.mean(0); ref0 = ref - refc

    dist = {}
    for name in ORDER:
        u = mda.Universe(f"{HERE}/topologies/{name}.psf", f"{HERE}/{name}_last5us.dcd")
        ag = u.select_atoms("protein and name CA and resid " + " ".join(map(str, ref_resids)))
        rr = np.argsort(ag.resids)
        vals = []
        for ts in u.trajectory[::STRIDE]:
            p = ag.positions[rr].astype(np.float64); com = p.mean(0)
            R, _ = rotation_matrix(p - com, ref0)
            d = ((p - com) @ R.T) - ref0
            vals.append(np.sqrt((d ** 2).sum(1).mean()))
        dist[name] = np.array(vals)
        print(f"  {name}: median core RMSD {np.median(dist[name]):.2f} Å ({len(vals)} frames)", flush=True)

    fig, ax = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
    grid = np.linspace(0, max(np.percentile(v, 99) for v in dist.values()) + 0.5, 400)
    for n in ORDER:
        lw = 3.0 if n == "2A" else 1.8
        z = 5 if n == "2A" else 2
        ax.plot(grid, gaussian_kde(dist[n])(grid), color=COL[n], lw=lw, zorder=z,
                label=f"RhoA-{n}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A")
                      + ("  (control)" if n == "2A" else ""))
    ax.set_xlabel("core G-domain Cα RMSD to RhoA-GTP-2 reference (Å)", fontsize=12)
    ax.set_ylabel("density", fontsize=12)
    ax.tick_params(labelsize=11)
    ax.legend(frameon=False, fontsize=10.5)
    ax.spines[["right", "top"]].set_visible(False)
    for e in ("png", "svg"):
        fig.savefig(f"{SUPP}/supp_2A_global_similarity.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP}/supp_2A_global_similarity.png")


if __name__ == "__main__":
    main()
