#!/usr/bin/env python3
"""Write the WT (GTP-2) co-binding frame (chosen by scan_cobinding_frame.py) INCLUDING the
membrane, for a PyMOL session: protein + all lipids (membrane-context surface) + the bridging
POPS + GTP/MG. Frame/lipid come from transition_states/cobinding_struct/chosen.txt.

  conda activate memb_analysis
  python write_cobinding_session_pdb.py
"""
import warnings; warnings.filterwarnings("ignore")
import os
import MDAnalysis as mda
import stable_state_config as C

HERE = os.path.dirname(os.path.abspath(__file__))
D = f"{HERE}/transition_states/cobinding_struct"
SIM = "GTP-2"


def main():
    info = {}
    for ln in open(f"{D}/chosen.txt"):
        k, _, v = ln.strip().partition("=")
        info[k] = v
    fr = int(info["frame"])
    sim = dict(next(s for s in C.SIMS if s["name"] == SIM))
    if not os.path.exists(sim["traj"]):
        sim["traj"] = f"{HERE}/{os.path.basename(sim['traj'])}"
    if not os.path.exists(sim["top"]):
        sim["top"] = f"{HERE}/topologies/{SIM}.psf"
    u = mda.Universe(sim["top"], sim["traj"])
    u.trajectory[fr]
    sel = u.select_atoms("protein or resname POPC POPE POPS PLA18 PSM GTP MG")
    out = f"{D}/{SIM}_cobind_f{fr}_membrane.pdb"
    sel.write(out)
    print(f"wrote {out}  ({sel.n_atoms} atoms, {sel.n_residues} residues)")


if __name__ == "__main__":
    main()
