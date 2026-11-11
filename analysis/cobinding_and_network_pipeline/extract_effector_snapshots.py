#!/usr/bin/env python3
"""For Figure 7: extract each sim's REPRESENTATIVE (effector-proximal Switch-I mode) frame relative to each
effector crystal, superpose it on that crystal's invariant G-domain core (same core used for the
RMSD), and write the transformed protein as a PDB. PyMOL then just loads crystal + the 5
core-aligned sim snapshots (no further alignment) and shows Switch I / Switch II.

Uses the per-frame RMSD arrays in transition_states/switches_fullres/effector_rmsd_fullres.npz to
pick the representative frame, so the snapshot matches the KDE it sits beside.

  conda activate memb_analysis
  python extract_effector_snapshots.py
"""
import warnings; warnings.filterwarnings("ignore")
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks


def closest_peak_rep(d1):
    """Frame index nearest the MEAN of the effector-proximal (lowest-RMSD) mode of the RMSD
    distribution, for visualization. For a bimodal distribution (e.g. RhoA-3A) this returns a
    frame from the effector-competent mode rather than an unrepresentative valley frame between
    modes; for a unimodal distribution it reduces to the ordinary mean of the single peak."""
    d1 = np.asarray(d1, float)
    xs = np.linspace(d1.min(), d1.max(), 400)
    ys = gaussian_kde(d1)(xs)
    peaks, _ = find_peaks(ys, prominence=ys.max() * 0.05)   # ignore shoulders/noise
    if len(peaks) == 0:
        peak_i = int(np.argmax(ys))
    else:
        peak_i = int(peaks[np.argmin(xs[peaks])])           # lowest-RMSD (effector-proximal) peak
    hi = d1.max()                                           # right edge of this mode
    later = [p for p in peaks if xs[p] > xs[peak_i]]
    if later:
        nxt = int(min(later, key=lambda p: xs[p]))
        valley = peak_i + int(np.argmin(ys[peak_i:nxt + 1]))
        hi = xs[valley]
    mode_mean = float(d1[d1 <= hi].mean())
    return int(np.argmin(np.abs(d1 - mode_mean))), mode_mean

HERE = os.path.dirname(os.path.abspath(__file__))
REFS_DIR = os.path.abspath(os.path.join(HERE, "..", "effector_refs"))
NPZ = os.path.join(HERE, "transition_states", "switches_fullres", "effector_rmsd_fullres.npz")
OUTDIR = os.path.join(HERE, "transition_states", "effector_snapshots")
SWITCH1 = (28, 42); SWITCH2 = (56, 78); RHO_INSERT = (120, 137)
CORE_SEL = (f"protein and name CA and resid 4-179 and not "
            f"(resid {SWITCH1[0]}-{SWITCH1[1]} or resid {SWITCH2[0]}-{SWITCH2[1]} "
            f"or resid {RHO_INSERT[0]}-{RHO_INSERT[1]})")
S1 = f"protein and name CA and resid {SWITCH1[0]}-{SWITCH1[1]}"
S2 = f"protein and name CA and resid {SWITCH2[0]}-{SWITCH2[1]}"
EFFECTORS = {"ROCK1 (1S1C)": ("1S1C.pdb", "A"), "PKN (1CXZ)": ("1CXZ.pdb", "A"),
             "anillin (4XOI)": ("4XOI.pdb", "A")}
ORDER = ["GDP", "GTP-1", "GTP-2", "3A"]   # 2A dropped


def sorted_pos(ag):
    o = np.argsort(ag.resids)
    return ag.resids[o], ag.positions[o].astype(np.float64)


def safe(lab):
    return lab.split(" ")[0]                     # "ROCK1 (1S1C)" -> "ROCK1"


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    z = np.load(NPZ, allow_pickle=True)
    eff_labels = list(z["eff_labels"])

    # crystal reference core/s1/s2 CA (matched resids), from the GTP-2 resid template
    u0 = mda.Universe(f"{HERE}/topologies/GTP-2.psf", f"{HERE}/GTP-2_last5us.dcd")
    core_r = np.sort(u0.select_atoms(CORE_SEL).resids)
    refs = {}
    for lab, (pdb, ch) in EFFECTORS.items():
        u = mda.Universe(f"{REFS_DIR}/{pdb}")
        chsel = f"(segid {ch} or chainID {ch}) and protein and name CA"
        ag = u.select_atoms(f"{chsel} and resid {' '.join(map(str, core_r))}")
        r, p = sorted_pos(ag); idx = {int(rr): pp for rr, pp in zip(r, p)}
        matched = [int(rr) for rr in core_r if rr in idx]
        refs[lab] = {"resids": matched, "pos": np.array([idx[m] for m in matched])}

    for lab in eff_labels:
        e = eff_labels.index(lab)
        R = refs[lab]; refc = R["pos"].mean(0)
        core_resid_sel = "protein and name CA and resid " + " ".join(map(str, R["resids"]))
        for name in ORDER:
            d1 = z[f"{name}|SwI"][e]                       # per-frame SwI RMSD to this effector
            stride = int(round(50000 / len(d1)))           # arrays were strided
            u = mda.Universe(f"{HERE}/topologies/{name}.psf", f"{HERE}/{name}_last5us.dcd")

            def super_write(fk, out):                      # superpose that frame on the crystal core -> PDB
                u.trajectory[fk * stride]
                _, c = sorted_pos(u.select_atoms(core_resid_sel)); com = c.mean(0)
                Rm, rc = rotation_matrix(c - com, R["pos"] - refc)
                prot = u.select_atoms("protein")
                prot.positions = (prot.positions.astype(np.float64) - com) @ Rm.T + refc
                prot.write(out)
                return rc

            rep_k, mode_mean = closest_peak_rep(d1)        # representative: mean of the effector-proximal mode
            rc = super_write(rep_k, f"{OUTDIR}/{safe(lab)}_{name}.pdb")
            print(f"  {safe(lab):6s} {name:5s} rep frame {rep_k*stride:5d} "
                  f"(SwI RMSD {d1[rep_k]:.2f}, mode-mean {mode_mean:.2f}, core-fit {rc:.2f})")
            if name in ("GTP-2", "3A"):                    # also the single closest (min-RMSD) frame, for a
                ck = int(np.argmin(d1))                    # separate 'how close does each get' overlay
                super_write(ck, f"{OUTDIR}/{safe(lab)}_{name}_closest.pdb")
                print(f"         {name:5s} closest frame {ck*stride:5d} (SwI RMSD {d1[ck]:.2f})")
    # record which crystals to load as reference
    with open(f"{OUTDIR}/effectors.txt", "w") as fh:
        for lab, (pdb, ch) in EFFECTORS.items():
            fh.write(f"{safe(lab)}\t{REFS_DIR}/{pdb}\t{ch}\n")
    print(f"\nwrote snapshots to {OUTDIR}")


if __name__ == "__main__":
    main()
