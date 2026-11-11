#!/usr/bin/env python3
"""Full-resolution protein<->lipid H-bond searches (ALL lipids), stride 1.

For each excerpt, searches G-domain<->lipid (gdlip) and anchor<->lipid (anclip)
hydrogen bonds over EVERY frame (0.1 ns) using MDAnalysis HydrogenBondAnalysis
(D-A <= 3.5 A, D-H-A angle >= 120 deg, polar donors/acceptors only). Writes
hbonds_stable/{name}_{anclip,gdlip}_raw_hbonds.npy (+ _meta.json). Resumable: a
search whose meta matches the current stride+frame count is skipped.

All lipids = POPS/POPE/POPC/PLA18 (config LIPID_SEL; PSM excluded = opposite
leaflet). The downstream cooperativity map/enhancement read these arrays and
filter to the amine-bearing subset themselves.

  python run_hbonds.py [--sims 2A,3A,GTP-2,GDP,GTP-1] [--stride 1]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse, json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import stable_state_config as C

HB = "hbonds_stable"
LIP_POLAR = f"({C.LIPID_SEL}) and not name C*"   # all lipids, H-bond-capable atoms


def done(meta_path, stride, n):
    if not os.path.exists(meta_path):
        return False
    try:
        m = json.load(open(meta_path))
        return m.get("stride") == stride and m.get("stop") == n
    except Exception:
        return False


def run_search(sim, tag, sel_a, sel_b, stride, kind):
    name = sim["name"]
    u = mda.Universe(sim["top"], sim["traj"])
    n = len(u.trajectory)
    mp = f"{HB}/{name}_{tag}_meta.json"
    if done(mp, stride, n):
        print(f"[skip] {name} {tag}: already done ({n} frames)", flush=True); return
    sel_all = f"({sel_a}) or ({sel_b})"
    hb = HBA(universe=u, between=[[sel_a, sel_b]],
             d_a_cutoff=C.HB_DA_CUTOFF, d_h_a_angle_cutoff=C.HB_ANGLE_CUTOFF)
    hb.hydrogens_sel = hb.guess_hydrogens(sel_all)
    hb.acceptors_sel = hb.guess_acceptors(sel_all)
    nf = len(range(0, n, stride))
    print(f"\n##### {name} {tag}: 0:{n}:{stride} = {nf} frames (all lipids) #####", flush=True)
    hb.run(start=0, stop=n, step=stride, verbose=True)
    np.save(f"{HB}/{name}_{tag}_raw_hbonds.npy", hb.results.hbonds)
    json.dump({"sim": name, "kind": kind, "start": 0, "stop": n, "stride": stride,
               "n_frames": nf, "traj": sim["traj"], "lipids": "POPS,POPE,POPC,PLA18"},
              open(mp, "w"), indent=2)
    print(f"[done] {name} {tag}: {hb.results.hbonds.shape[0]} bonds", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="2A,3A,GTP-2,GDP,GTP-1")
    ap.add_argument("--stride", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(HB, exist_ok=True)
    want = args.sims.split(",")
    for sim in C.SIMS:
        if sim["name"] not in want:
            continue
        if not os.path.exists(sim["traj"]):
            print(f"[SKIP] {sim['name']}: excerpt not found at {sim['traj']}", flush=True); continue
        # anchor is cheaper -> run first so partial results are useful early
        run_search(sim, "anclip", C.ANCHOR_SEL, LIP_POLAR, args.stride, "anchor_lipid")
        run_search(sim, "gdlip", C.GDOMAIN_SEL, LIP_POLAR, args.stride, "gdomain_lipid")
    print("\n=== ALL H-BOND SEARCHES DONE ===", flush=True)


if __name__ == "__main__":
    main()
