#!/usr/bin/env python3
"""Re-tally the cooperativity analysis over a SHIFTED post window, WITHOUT
re-searching and WITHOUT loading the giant trajectories.

The excerpt is the last 1.95 us; 2A's starts right at the transition tail, so its
last-1.95 window isn't fully settled. This re-tallies any sub-window (default the
last 1.5 us = drop the first 0.45 us) using only the PSF (atom names) + the
existing hbonds_stable/*.npy arrays, into a sub-directory of transition_states/gains.

  python replot_window.py --start-us 0.45 --end-us 1.95 --subdir last1.5us
"""
import warnings; warnings.filterwarnings("ignore")
import os, sys, argparse
import numpy as np, pandas as pd
import MDAnalysis as mda
import stable_state_config as C
import lipid_hbond_analysis as L


def load_bonds_topo_only(sim):
    """lipid_hbond_analysis.load_bonds but Universe(TOPOLOGY ONLY) -- re-tally
    needs atom names + the raw arrays, not the coordinates."""
    name = sim["name"]
    u = mda.Universe(sim["top"])                     # no trajectory
    resid, resname, names = u.atoms.resids, u.atoms.resnames, u.atoms.names
    lipset = set(C.LIPIDS)
    is_lip = np.array([rn in lipset for rn in resname])
    w = C.state_frames(name, "post", sim["dt_ns"])
    frames = []
    for tag in ("gdlip", "anclip"):
        p = f"{L.HB}/{name}_{tag}_raw_hbonds.npy"
        if not os.path.exists(p):
            continue
        hb = np.load(p)
        d = hb[:, 1].astype(int); a = hb[:, 3].astype(int); fr = hb[:, 0].astype(int)
        dl, al = is_lip[d], is_lip[a]
        keep = (dl ^ al) & (fr >= w[0]) & (fr < w[1])
        prot = np.where(dl, a, d)[keep]; lip = np.where(dl, d, a)[keep]
        role = np.where(dl[keep], "acceptor", "donor")
        frames.append(pd.DataFrame({
            "frame": fr[keep], "pres": resid[prot], "prn": resname[prot],
            "patom": names[prot], "role": role,
            "lres": resid[lip], "lrn": resname[lip], "latom": names[lip]}))
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    df["pclass"] = np.where(df["prn"].isin(L.BASIC), "basic",
                            np.where(df["prn"].isin(L.ACIDIC), "acidic", "other"))
    df["lgroup"] = df["latom"].map(L.lipgroup)
    df["npost"] = w[1] - w[0]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start-us", type=float, default=0.45)   # last 1.5 of the 1.95 excerpt
    ap.add_argument("--end-us", type=float, default=1.95)
    ap.add_argument("--subdir", default="last1.5us")
    a = ap.parse_args()

    names = [s["name"] for s in C.SIMS]
    C.STATES = {n: {"pre": None, "cross": None, "post": (a.start_us, a.end_us)} for n in names}
    nfr = round((a.end_us - a.start_us) * 1000.0 / C.DT_NS)
    print(f"window {a.start_us}-{a.end_us} us  ({nfr} frames)  -> transition_states/gains/{a.subdir}/\n")

    # topology-only re-tally + redirect all outputs to the sub-directory
    L.load_bonds = load_bonds_topo_only
    out = f"transition_states/gains/{a.subdir}"
    os.makedirs(out, exist_ok=True)
    L.OUT = out
    import cobinding_occupancy as O; O.OUT = out
    import viz_cooperativity as V; V.OUT = out

    for lipids in (["POPS"], ["POPS", "POPE", "PLA18"]):
        sys.argv = ["x", "--sims", "2A,3A,GTP-2,GDP,GTP-1", "--lipids", ",".join(lipids)]
        L.main()
    for ls in ("", "_amine"):
        for cut in (25, 50):
            sys.argv = ["x", "--lipidset", ls, "--cutoff", str(cut), "--min-frames", str(cut)]
            try:
                V.main()
            except SystemExit:
                pass
    sys.argv = ["x", "--sims", "GDP,GTP-1,GTP-2,2A,3A", "--lipids", "POPS", "--tag", "POPS"]; O.main()
    sys.argv = ["x", "--sims", "GDP,GTP-1,GTP-2,2A,3A", "--lipids", "POPS,POPE,PLA18"]; O.main()
    print(f"\nDONE -> transition_states/gains/{a.subdir}/")


if __name__ == "__main__":
    main()
