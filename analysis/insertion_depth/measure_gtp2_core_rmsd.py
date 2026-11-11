#!/usr/bin/env python3
"""Remeasure RhoA-GTP-2 catalytic-domain Cα RMSD referenced to the trajectory's FIRST frame,
so Fig 3A starts at 0 and rises to its plateau (the stored rmsdrhoa.pkl used a non-frame-0
reference, leaving the curve already elevated at t=0).

Selection matches the Fig 3 caption: G-domain Cα (resid 4-179) EXCLUDING switch I (28-42) and
switch II (56-78). Aligned on that same core, ref_frame = 0.

  conda activate memb_analysis ; python measure_gtp2_core_rmsd.py [--stride N]
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, pickle
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD

ROOT = Path("/home/hutchcha/rhoa_mut_analysis")
TOP = ROOT / "WT_trajs/RhoA-GTP-2.psf"
TRJ = ROOT / "WT_trajs/RhoA_GTP-2.dcd"
HERE = Path(__file__).resolve().parent
DT_NS = 1.008

CORE = ("protein and name CA and resid 4-179 and not "
        "(resid 28-42 or resid 56-78)")

ap = argparse.ArgumentParser()
ap.add_argument("--stride", type=int, default=1)
args = ap.parse_args()

u = mda.Universe(str(TOP), str(TRJ))
n = len(u.trajectory)
print(f"frames={n}  core atoms={u.select_atoms(CORE).n_atoms}  stride={args.stride}", flush=True)

r = RMSD(u, u, select=CORE, ref_frame=0)
r.run(step=args.stride, verbose=False)
arr = r.results.rmsd                      # cols: frame, time, rmsd
rmsd_A = arr[:, 2]
time_us = (arr[:, 0] * DT_NS) / 1000.0    # frame index -> µs

print(f"rmsd first3={np.round(rmsd_A[:3],3)}  min={rmsd_A.min():.3f}  "
      f"median={np.median(rmsd_A):.3f}  max={rmsd_A.max():.3f}", flush=True)

out = HERE / f"gtp2_core_rmsd_ref0{'_s'+str(args.stride) if args.stride>1 else ''}.pkl"
pickle.dump({"time_us": time_us, "rmsd_A": rmsd_A}, open(out, "wb"))
print(f"wrote {out.name} ({len(rmsd_A)} pts)", flush=True)
