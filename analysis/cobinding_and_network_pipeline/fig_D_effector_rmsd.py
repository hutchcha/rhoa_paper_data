#!/usr/bin/env python3
"""Panel D (main figure): Switch I / Switch II Cα RMSD to the effector-bound RhoA conformation,
per sim, on the FRAME-MATCHED FULL-RES trajectories — the structural pivot showing the mutation
re-poises the effector switches. RMSD is to the NEAREST of three membrane-free effector crystals
(ROCK1 1S1C, PKN 1CXZ, anillin 4XOI), core-superposed (invariant G-domain core). Distributions
(KDE) per sim; lower = closer to the effector-competent conformation.

Effector crystals from ../effector_refs/. Saves per-frame arrays (npz) + the KDE panel.

  python fig_D_effector_rmsd.py [--stride 2]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"; plt.rcParams["svg.fonttype"] = "none"
HERE = os.path.dirname(os.path.abspath(__file__))
REFS_DIR = os.path.abspath(os.path.join(HERE, "..", "effector_refs"))
# NB the MAIN effector figure is now the per-effector composite (compose_figure7.py); this script
# still generates the per-frame RMSD npz + the supplementary per-effector grid. Its nearest-effector
# KDE (plot_main) now writes to figure7/ as a secondary/legacy view.
OUT = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure7"))
SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))
NPZ = os.path.join(HERE, "transition_states", "switches_fullres")
SWITCH1 = (28, 42); SWITCH2 = (56, 78); RHO_INSERT = (120, 137)
CORE_SEL = (f"protein and name CA and resid 4-179 and not "
            f"(resid {SWITCH1[0]}-{SWITCH1[1]} or resid {SWITCH2[0]}-{SWITCH2[1]} "
            f"or resid {RHO_INSERT[0]}-{RHO_INSERT[1]})")
S1 = f"protein and name CA and resid {SWITCH1[0]}-{SWITCH1[1]}"
S2 = f"protein and name CA and resid {SWITCH2[0]}-{SWITCH2[1]}"
EFFECTORS = {"ROCK1 (1S1C)": ("1S1C.pdb", "A"), "PKN (1CXZ)": ("1CXZ.pdb", "A"),
             "anillin (4XOI)": ("4XOI.pdb", "A")}
ORDER = ["GDP", "GTP-1", "GTP-2", "3A"]   # 2A dropped
COL = {"GDP": "#B279A2", "GTP-1": "#72B7B2", "GTP-2": "#E45756", "2A": "#F58518", "3A": "#4C78A8"}


def sorted_pos(ag):
    o = np.argsort(ag.resids)
    return ag.resids[o], ag.positions[o].astype(np.float64)


def load_refs(core_r, s1_r, s2_r):
    refs = {}
    for lab, (pdb, ch) in EFFECTORS.items():
        u = mda.Universe(f"{REFS_DIR}/{pdb}")
        chsel = f"(segid {ch} or chainID {ch}) and protein and name CA"

        def grab(want):
            ag = u.select_atoms(f"{chsel} and resid {' '.join(map(str, want))}")
            r, p = sorted_pos(ag); idx = {int(rr): pp for rr, pp in zip(r, p)}
            matched = [int(rr) for rr in want if rr in idx]
            return matched, np.array([idx[m] for m in matched])
        refs[lab] = {"core": grab(core_r), "s1": grab(s1_r), "s2": grab(s2_r)}
    return refs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stride", type=int, default=2)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True); os.makedirs(NPZ, exist_ok=True)

    # resid templates from the first sim
    u0 = mda.Universe(f"{HERE}/topologies/GTP-2.psf", f"{HERE}/GTP-2_last5us.dcd")
    core_r = np.sort(u0.select_atoms(CORE_SEL).resids)
    s1_r = np.sort(u0.select_atoms(S1).resids); s2_r = np.sort(u0.select_atoms(S2).resids)
    refs = load_refs(core_r, s1_r, s2_r)
    print("refs loaded:", {k: (len(v['core'][0]), len(v['s1'][0]), len(v['s2'][0])) for k, v in refs.items()})

    dist = {}
    for name in ORDER:
        u = mda.Universe(f"{HERE}/topologies/{name}.psf", f"{HERE}/{name}_last5us.dcd")
        # per-effector matched sim AtomGroups
        per_eff = {}
        for lab, R in refs.items():
            per_eff[lab] = (u.select_atoms("protein and name CA and resid " + " ".join(map(str, R["core"][0]))),
                            u.select_atoms("protein and name CA and resid " + " ".join(map(str, R["s1"][0]))),
                            u.select_atoms("protein and name CA and resid " + " ".join(map(str, R["s2"][0]))))
        frames = np.arange(0, len(u.trajectory), args.stride)
        d1 = np.empty((len(refs), len(frames))); d2 = np.empty((len(refs), len(frames)))
        for k, f in enumerate(frames):
            u.trajectory[int(f)]
            for e, (lab, R) in enumerate(refs.items()):
                core_ag, s1_ag, s2_ag = per_eff[lab]
                _, c = sorted_pos(core_ag); com = c.mean(0)
                refc = R["core"][1].mean(0)
                Rm, _ = rotation_matrix(c - com, R["core"][1] - refc)
                _, p1 = sorted_pos(s1_ag); _, p2 = sorted_pos(s2_ag)
                dd1 = ((p1 - com) @ Rm.T) - (R["s1"][1] - refc)
                dd2 = ((p2 - com) @ Rm.T) - (R["s2"][1] - refc)
                d1[e, k] = np.sqrt((dd1 ** 2).sum(1).mean()); d2[e, k] = np.sqrt((dd2 ** 2).sum(1).mean())
        dist[name] = {"SwI": d1, "SwII": d2}                   # full per-effector arrays (3, nframes)
        print(f"  [{name}] {len(frames)} frames | SwI med {np.median(d1.min(0)):.2f} "
              f"SwII med {np.median(d2.min(0)):.2f}", flush=True)
    eff_labels = list(refs.keys())
    np.savez_compressed(f"{NPZ}/effector_rmsd_fullres.npz",
                        eff_labels=np.array(eff_labels),
                        **{f"{n}|{sw}": dist[n][sw] for n in dist for sw in ("SwI", "SwII")})
    plot_main(dist, eff_labels)
    plot_supp(dist, eff_labels)


def _nearest_effector(arr, eff_labels):
    """dominant nearest effector (mode of per-frame argmin) + its fraction."""
    am = arr.argmin(0)
    vals, cts = np.unique(am, return_counts=True)
    i = vals[cts.argmax()]
    return eff_labels[i].split(" ")[0], cts.max() / len(am)


def plot_main(dist, eff_labels):
    os.makedirs(OUT, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    for ax, sw, xr in [(axes[0], "SwI", (2, 18)), (axes[1], "SwII", (2, 10))]:
        grid = np.linspace(*xr, 400)
        for n in ORDER:
            mn = dist[n][sw].min(0)
            near, frac = _nearest_effector(dist[n][sw], eff_labels)
            ax.plot(grid, gaussian_kde(mn)(grid), color=COL[n], lw=2,
                    label=f"RhoA-{n}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A")
                          + f"  (nearest: {near}, {frac*100:.0f}%)")
        ax.set_xlabel(f"{'Switch I' if sw=='SwI' else 'Switch II'} Cα RMSD to nearest effector (Å)", fontsize=11)
        ax.set_ylabel("density", fontsize=11)
        ax.set_title(f"{'A  Switch I' if sw=='SwI' else 'B  Switch II'}", fontsize=12, loc="left", fontweight="bold")
        ax.spines[["right", "top"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8.5)
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/fig_D_effector_rmsd.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_D_effector_rmsd.png")


def plot_supp(dist, eff_labels):
    """Supplementary: RMSD distribution to EACH effector separately (rows=switch, cols=effector)."""
    os.makedirs(SUPP, exist_ok=True)
    fig, axes = plt.subplots(2, len(eff_labels), figsize=(4.6 * len(eff_labels), 8),
                             constrained_layout=True)
    for row, (sw, xr) in enumerate([("SwI", (2, 18)), ("SwII", (2, 10))]):
        grid = np.linspace(*xr, 400)
        for col, lab in enumerate(eff_labels):
            ax = axes[row, col]
            for n in ORDER:
                ax.plot(grid, gaussian_kde(dist[n][sw][col])(grid), color=COL[n], lw=2,
                        label=f"RhoA-{n}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"))
            ax.set_title(f"{'Switch I' if sw=='SwI' else 'Switch II'} → {lab}", fontsize=11, loc="left",
                         fontweight="bold")
            ax.set_xlabel("Cα RMSD (Å)", fontsize=10); ax.set_ylabel("density", fontsize=10)
            ax.spines[["right", "top"]].set_visible(False)
            if row == 0 and col == 0:
                ax.legend(frameon=False, fontsize=9)
    fig.suptitle("Switch I / Switch II Cα RMSD to each effector-bound crystal (full-res)",
                 fontsize=13, fontweight="bold")
    for e in ("png", "svg"):
        fig.savefig(f"{SUPP}/supp_effector_rmsd_per_effector.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP}/supp_effector_rmsd_per_effector.png")


if __name__ == "__main__":
    main()
