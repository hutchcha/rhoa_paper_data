#!/usr/bin/env python3
"""Re-measure geranylgeranyl-anchor insertion depth using ONLY the acyl-chain carbons (CYSG
C1-C20) as the anchor, instead of the whole CYSG residue. Convention matches the suite /
insertion_vs_thickness: depth = membrane_heavy_COM_z - anchor_COM_z  (distance below the bilayer
centre; center-referenced). Full trajectories, every frame (same length as the original pkl so the
Fig-1/Fig-3 moving-average + time axis are unchanged).

  conda activate memb_analysis ; python measure_anchor_insertion_acyl.py
Writes anchor_insert_depth_acyl.pkl -> {RhoA-GDP, RhoA-GTP, RhoA-GTP-2}  (keys match fig1/fig3).
"""
import warnings; warnings.filterwarnings("ignore")
import os, pickle, numpy as np
import MDAnalysis as mda

HERE = os.path.dirname(os.path.abspath(__file__))
TRJ = "/home/hutchcha/rhoa_mut_analysis/WT_trajs"
PREN = "resname CYSG and name " + " ".join(f"C{i}" for i in range(1, 21))   # acyl chain C1-C20
MEMH = "resname POPC POPE POPS PLA18 PSM and not name H*"

SYS = [("RhoA-GDP",   "RhoA-GDP.psf",   "RhoA_GDP.dcd"),
       ("RhoA-GTP",   "RhoA-GTP.psf",   "RhoA_GTP-1.dcd"),      # GTP = GTP-1
       ("RhoA-GTP-2", "RhoA-GTP-2.psf", "RhoA_GTP-2.dcd")]

out = {}
for name, psf, dcd in SYS:
    u = mda.Universe(f"{TRJ}/{psf}", f"{TRJ}/{dcd}")
    anc = u.select_atoms(PREN)
    mem = u.select_atoms(MEMH)
    n = len(u.trajectory)
    d = np.empty(n, float)
    for i, ts in enumerate(u.trajectory):
        d[i] = mem.center_of_mass()[2] - anc.center_of_mass()[2]
        if i % 2000 == 0:
            print(f"  [{name}] frame {i}/{n}", flush=True)
    out[name] = d
    print(f"[{name}] done: n={n} mean={d.mean():.2f} last5us={d[int(n*0.75):].mean():.2f}", flush=True)

pickle.dump(out, open(f"{HERE}/anchor_insert_depth_acyl.pkl", "wb"))
print("wrote anchor_insert_depth_acyl.pkl")
for k, v in out.items():
    print(f"  {k}: full {v.mean():.2f}  last5us {v[int(len(v)*0.75):].mean():.2f} Å")
