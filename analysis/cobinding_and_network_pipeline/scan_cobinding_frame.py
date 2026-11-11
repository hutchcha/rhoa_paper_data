#!/usr/bin/env python3
"""Pick a representative WT (GTP-2) lipid co-binding frame FROM THE EXISTING H-bond data
(the same precomputed bonds behind Fig 6A-C), not a fresh trajectory scan. Using load_bonds()
guarantees the render matches the panels' co-binding definition exactly.

We want the structural picture behind the network: N-patch basics (R5/K6/K7) and Switch II
acidic D78 (D76 shows ~no co-binding in GTP-2) H-bonding the SAME inner-leaflet amine lipid.
The occupancy table says this is a POPS in ~100% of frames. We choose the frame where the most
N-patch residues co-grip one POPS together with D78, then write just that frame as a PDB.

  conda activate memb_analysis
  python scan_cobinding_frame.py
"""
import warnings; warnings.filterwarnings("ignore")
import os
from collections import defaultdict
import MDAnalysis as mda
import stable_state_config as C
from lipid_hbond_analysis import load_bonds

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = f"{HERE}/transition_states/cobinding_struct"
SIM = "GTP-2"
NPATCH = {5, 6, 7}
ACID = {76, 78}
AMINE = {"POPS", "POPE", "PLA18"}


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sim = dict(next(s for s in C.SIMS if s["name"] == SIM))
    if not os.path.exists(sim["traj"]):        # config builds traj for repo-root cwd; fix here
        sim["traj"] = f"{HERE}/{os.path.basename(sim['traj'])}"
    if not os.path.exists(sim["top"]):
        sim["top"] = f"{HERE}/topologies/{SIM}.psf"
    df = load_bonds(sim)
    df = df[(df.lrn.isin(AMINE)) & (df.pclass != "other")]

    # per (frame, lipid molecule): which N-patch basics and which acidics grip it
    grip_np = defaultdict(set); grip_ac = defaultdict(set)
    for r in df.itertuples():
        key = (r.frame, r.lres, r.lrn)
        if r.pclass == "basic" and r.pres in NPATCH:
            grip_np[key].add(r.pres)
        elif r.pclass == "acidic" and r.pres in ACID:
            grip_ac[key].add(r.pres)

    cands = []  # (n_np, n_ac, has_pops, frame, lres, lrn, npatch, acids)
    for key in set(grip_np) & set(grip_ac):
        fr, lres, lrn = key
        nps, acs = grip_np[key], grip_ac[key]
        cands.append((len(nps), len(acs), lrn == "POPS", fr, lres, lrn,
                      tuple(sorted(nps)), tuple(sorted(acs))))
    if not cands:
        print("no frame with N-patch + D76/78 on the same amine lipid"); return
    # best tier = most N-patch residues, then most acidics, then POPS
    best_key = max((t[0], t[1], t[2]) for t in cands)
    tier = [t for t in cands if (t[0], t[1], t[2]) == best_key]
    npost = int(df["npost"].iloc[0])
    mid = npost // 2                          # representative mid-window frame (avoid frame 0)
    tier.sort(key=lambda t: abs(t[3] - mid))

    print(f"{len(cands)} co-binding frames total; best tier "
          f"({best_key[0]} N-patch + {best_key[1]} acid, POPS={best_key[2]}) has {len(tier)} frames.")
    print("nearest-to-midpoint candidates:")
    for n_np, n_ac, pops, fr, lres, lrn, nps, acs in tier[:5]:
        print(f"  frame {fr:5d}  lipid {lres} {lrn:5s}  N-patch {nps}  acid {acs}")

    n_np, n_ac, pops, fr, lres, lrn, nps, acs = tier[0]
    u = mda.Universe(sim["top"], sim["traj"])
    u.trajectory[int(fr)]
    sel = u.select_atoms(f"protein or (resname {lrn} and resid {lres}) or resname GTP MG")
    out = f"{OUTDIR}/{SIM}_cobind_f{fr}_lip{lres}.pdb"
    sel.write(out)
    with open(f"{OUTDIR}/chosen.txt", "w") as fh:
        fh.write(f"frame={fr}\nlipid_resid={lres}\nlipid_type={lrn}\n"
                 f"npatch={nps}\nacid={acs}\npdb={out}\n")
    print(f"\nCHOSEN frame {fr}, lipid {lres} ({lrn}), N-patch {nps}, acid {acs}\nwrote {out}")


if __name__ == "__main__":
    main()
