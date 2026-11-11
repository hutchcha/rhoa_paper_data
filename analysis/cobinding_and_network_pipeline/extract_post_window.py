#!/usr/bin/env python3
"""Write out the LAST `WINDOW_US` microseconds of each full-resolution trajectory,
with a matched (equal) number of frames, for post-state resampling.

Run this ON THE GPU MACHINE where the full trajectories live. It writes the
excerpts into the SAME directory as this script (no wrapping -- trajectories are
assumed already wrapped). rsync the excerpts back afterward.

Matched frames: N_FRAMES is computed once from WINDOW_US / DT_NS and applied to
every trajectory, so all outputs have the same frame count (assuming they share
the native DT_NS -- which they will if they came off the same run settings). If
a trajectory is shorter than the window, all of it is written and a warning is
printed. Set DRY_RUN = True first to just print the frame math without writing.

  python extract_post_window.py
"""
import os
import MDAnalysis as mda
import MDAnalysis.transformations as trn

# =========================== EDIT THIS BLOCK ===========================
WINDOW_US = 1.95        # keep the last 1.8 us of each trajectory
DT_NS     = 0.1        # native time PER FRAME (ns) of the FULL trajectories
OUT_EXT   = "dcd"      # output trajectory format (matches input)
DRY_RUN   = False      # True = only print frame counts, write nothing
N_THREADS = 10         # threads for the (parallelizable) wrap/center transforms

# On-the-fly wrap+center (only for trajectories with "wrap": True below):
# unwrap lipids -> center the CYSG Ca in the box by mass -> re-wrap lipid fragments.
CENTER_SEL = "resname CYSG and name CA"   # centering anchor (prenyl Cys Ca)
WRAP_SEL   = "not protein"                # atoms to unwrap / wrap (lipids, ions, water)

# GTP-2's full sim is split into 4 parts (nowater1..4.dcd) -- give them as a LIST
# and MDAnalysis reads them as one continuous trajectory (ChainReader).
_GTP2_DIR = "/mnt/nas2/chase_projects/rhoa/gtpsimredo/charmm36/charmm-gui-8253913405/namd/chasesim"
_GTP2_TRAJS = [f"{_GTP2_DIR}/nowater{i}.dcd" for i in (1, 2, 3, 4)]

# One entry per trajectory. "traj" may be a single path OR a list of paths.
# "wrap": True  -> apply the wrap/center workflow  (2A, 3A, GTP-2)
# "wrap": False -> already wrapped/centered, write as-is  (GDP, GTP-1)
SIMS = [
    {"name": "2A",    "top": "/mnt/data/rhoa/rhoa_2A_mutant_sim/chasesim/namd/namd_stages/nowat.psf",
                      "traj": "/mnt/data/rhoa/rhoa_2A_mutant_sim/anton_trajs/2A_nowater_10us.dcd", "wrap": True},
    {"name": "3A",    "top": "/mnt/data/rhoa/rhoa_3A_sim/anton_trajectories/nowat.psf",
                      "traj": "/mnt/data/rhoa/rhoa_3A_sim/anton_trajectories/3A_nowater_10us.dcd", "wrap": True},
    {"name": "GTP-2", "top": "/mnt/data_SSD/rhoa_redotraj/nowat.psf",
                      "traj": _GTP2_TRAJS, "wrap": True},
    {"name": "GDP",   "top": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/sys.psf",
                      "traj": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/cenwrapprod.dcd", "wrap": False},
    {"name": "GTP-1", "top": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/sys.psf",
                      "traj": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/cenwrapprod.dcd", "wrap": False},
]

# Optional: force an EXACT frame count (overrides the WINDOW_US/DT_NS math).
# e.g. set to 2A's post-state frame count. Leave as None to compute from the window.
N_FRAMES_OVERRIDE = None
# =======================================================================

OUTDIR = os.path.dirname(os.path.abspath(__file__))
N_FRAMES = int(N_FRAMES_OVERRIDE) if N_FRAMES_OVERRIDE else int(round(WINDOW_US * 1000.0 / DT_NS))


def main():
    print(f"output dir : {OUTDIR}")
    print(f"target     : last {WINDOW_US} us = {N_FRAMES} frames per trajectory "
          f"(dt = {DT_NS} ns/frame)")
    print(f"mode       : {'DRY RUN (nothing written)' if DRY_RUN else 'WRITING'}\n")

    for s in SIMS:
        name, top, traj = s["name"], s["top"], s["traj"]
        trajs = list(traj) if isinstance(traj, (list, tuple)) else [traj]   # allow multi-DCD
        missing = [p for p in [top] + trajs if not os.path.exists(p)]
        if missing:
            print(f"[SKIP] {name}: path(s) not found: {missing}\n"); continue
        # a list of trajectories is read as one continuous trajectory (ChainReader)
        u = mda.Universe(top, trajs)
        n = len(u.trajectory)
        if len(trajs) > 1:
            print(f"[{name}] {len(trajs)} trajectory parts chained")
        start = max(0, n - N_FRAMES)
        kept = n - start
        out = os.path.join(OUTDIR, f"{name}_last{WINDOW_US}us.{OUT_EXT}")

        do_wrap = s.get("wrap", False)
        if do_wrap:
            center = u.select_atoms(CENTER_SEL)
            ag = u.select_atoms(WRAP_SEL)
            if len(center) == 0:
                print(f"   [WARN] {name}: CENTER_SEL '{CENTER_SEL}' matched 0 atoms; "
                      f"writing WITHOUT wrap/center"); do_wrap = False
            else:
                u.trajectory.add_transformations(
                    trn.unwrap(ag, parallelizable=True, max_threads=N_THREADS),
                    trn.center_in_box(center, center="mass", parallelizable=True, max_threads=N_THREADS),
                    trn.wrap(ag, compound="fragments", parallelizable=True, max_threads=N_THREADS),
                )

        print(f"[{name}] {n} frames total ({n * DT_NS / 1000:.3f} us), {u.atoms.n_atoms} atoms"
              f"  | wrap/center: {'YES (on CYSG CA)' if do_wrap else 'no (already wrapped)'}")
        if kept < N_FRAMES:
            print(f"   [WARN] shorter than the window -> writing all {kept} frames "
                  f"({kept * DT_NS / 1000:.3f} us), fewer than the {N_FRAMES}-frame target")
        print(f"   keeping frames {start}:{n}  ({kept} frames, {kept * DT_NS / 1000:.3f} us)")
        print(f"   -> {out}")

        if DRY_RUN:
            print(); continue

        with mda.Writer(out, n_atoms=u.atoms.n_atoms) as W:
            for i, ts in enumerate(u.trajectory[start:n]):
                W.write(u.atoms)
                if i % 2000 == 0:
                    print(f"      {i}/{kept} frames written", flush=True)
        print(f"   done ({kept} frames)\n", flush=True)

    print("ALL DONE" + (" (dry run)" if DRY_RUN else ""))


if __name__ == "__main__":
    main()
