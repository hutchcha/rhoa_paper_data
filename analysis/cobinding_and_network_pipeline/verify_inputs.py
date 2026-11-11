#!/usr/bin/env python3
"""Pre-flight check: confirm every sim's PSF matches its excerpt trajectory.

Loads each (topology, trajectory) pair from stable_state_config.SIMS and prints
the atom count + frame count. A PSF whose atom count doesn't match the DCD raises
immediately here (loudly) instead of failing hours into the H-bond search. Also
warns if the frame count != round(POST_US/DT_NS) (i.e. DT_NS/POST_US are off, which
would corrupt occupancy denominators).

  python verify_inputs.py
"""
import warnings; warnings.filterwarnings("ignore")
import os
import MDAnalysis as mda
import stable_state_config as C

expected = round(C.POST_US * 1000.0 / C.DT_NS)
print(f"expected frames per excerpt: round({C.POST_US} us / {C.DT_NS} ns) = {expected}\n")

ok = True
for s in C.SIMS:
    name, top, traj = s["name"], s["top"], s["traj"]
    if not os.path.exists(top):
        print(f"[MISSING PSF ] {name:6}: {top}"); ok = False; continue
    if not os.path.exists(traj):
        print(f"[MISSING TRAJ] {name:6}: {traj}"); ok = False; continue
    try:
        u = mda.Universe(top, traj)
        na, nf = u.atoms.n_atoms, len(u.trajectory)
        dur = nf * s["dt_ns"] / 1000.0
        note = "" if nf == expected else f"   [!] frames {nf} != expected {expected} -- check DT_NS/POST_US"
        if note:
            ok = False
        print(f"[OK] {name:6}: {na:>6} atoms, {nf:>6} frames ({dur:.3f} us){note}")
    except Exception as e:
        print(f"[FAIL] {name:6}: {type(e).__name__}: {e}")
        print(f"        -> PSF and trajectory atom counts likely differ; point "
              f"'{name}' top at the PSF used to write this excerpt")
        ok = False

import sys
print("\n" + ("ALL INPUTS OK -- safe to launch run_all.sh"
              if ok else "FIX THE ISSUES ABOVE before running run_all.sh"))
sys.exit(0 if ok else 1)
