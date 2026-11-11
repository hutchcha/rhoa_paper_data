#!/usr/bin/env python3
"""
RhoA-only: per-frame SASA (Shrake-Rupley) for one residue using MDTraj,
then plot with your line_plot() using a 50 ns moving average.

NEW:
- Adds a "plot-only" option that loads previously saved CSVs and skips MDTraj.

Selection for the target residue is by:
(resSeq == TARGET_RESID) AND (residue.name == TARGET_RESNAME)

Given:
- dt = 0.1 ns per frame in the original trajectories
- moving average is set in ns and converted to frames based on dt_effective = dt_ns * stride_used

IMPORTANT:
- For plot-only mode to use the correct MA window, the script either:
  (A) infers dt from the CSV time axis, or
  (B) uses your DT_NS and STRIDE settings below.
  It will try (A) first.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import mdtraj as md
from tqdm import tqdm

# -----------------------
# Your style + line plot
# -----------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["svg.fonttype"] = "none"

def line_plot(xvalues, yvalues, ax, title, color, z, label, ma_window=5, ma_color="red", ma_z=2):
    """
    Plots the original data with transparency and overlays a moving average.
    """
    xvalues = np.array(xvalues)
    yvalues = np.array(yvalues)

    moving_avg = pd.Series(yvalues).rolling(window=ma_window, center=True).mean()

    ax.plot(xvalues, yvalues, alpha=0.25, color=color, zorder=z, label=None)
    ax.plot(xvalues, moving_avg, alpha=1.0, color=ma_color, zorder=ma_z)

    ax.set_title(title, fontsize=20)
    ax.tick_params(axis="both", labelsize=16)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")


# -----------------------
# RhoA sim paths (direct)
# -----------------------
systems = {
    "RhoA-GDP": {
        "top": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/sys.psf",
        "traj": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/cenwrapprod.dcd",
        "color": "#2ca02c",
    },
    "RhoA-GTP": {
        "top": "/mnt/data/rhoa/rhoagtpanalysis/20us/sys.psf",
        "traj": "/mnt/data/rhoa/rhoagtpanalysis/20us/cenwrapprod.dcd",
        "color": "#d62728",
    },
    "RhoA-Resim": {
        "top": "/mnt/nas/chase_projects/rhoa/rhoaredo/chasesim/gromacs/antonfiles/analysis/1us/nowat.psf",
        "traj": "/mnt/nas/chase_projects/rhoa/rhoaredo/chasesim/gromacs/antonfiles/analysis/20us/cenwrapprod.dcd",
        "color": "#9467bd",
    },
}

def mdtraj_res_index_by_resid_resname(topology: md.Topology, resid: int, resname: str) -> int:
    resname = resname.upper()
    hits = [r.index for r in topology.residues if (r.resSeq == resid and r.name.upper() == resname)]
    if len(hits) == 0:
        raise ValueError(f"Residue not found: {resname}{resid} (no matches for resSeq+name).")
    if len(hits) > 1:
        raise ValueError(f"Ambiguous: {resname}{resid} matched multiple residues: {hits}")
    return hits[0]


def sasa_timeseries_one_residue(
    top_path: str,
    traj_path: str,
    resid: int,
    resname: str,
    dt_ns: float = 0.1,
    stride: int = 1,
    chunk: int = 5000,
    probe_radius_nm: float = 0.14,
    n_sphere_points: int = 960,
    atom_selection: str = "protein",
):
    """
    Stream trajectory -> return time (ns) and SASA (Å^2) for ONE residue.
    SASA is computed using ONLY atoms in `atom_selection` (e.g., "protein" ignores lipids/water shielding).
    """
    full_top = md.load_topology(top_path)
    atom_indices = full_top.select(atom_selection)
    if atom_indices.size == 0:
        raise RuntimeError(f"Selection '{atom_selection}' returned 0 atoms.")
    sub_top = full_top.subset(atom_indices)

    resname_u = resname.upper()
    hits = [r.index for r in sub_top.residues if (r.resSeq == resid and r.name.upper() == resname_u)]
    if len(hits) == 0:
        raise ValueError(f"{resname}{resid} not found in sub_top for selection '{atom_selection}'.")
    if len(hits) > 1:
        raise ValueError(f"Ambiguous: {resname}{resid} matched multiple residues in sub_top: {hits}")
    res_i = hits[0]

    t_ns = []
    sasa_A2 = []
    frame0 = 0

    for trj in tqdm(
        md.iterload(traj_path, top=top_path, stride=stride, chunk=chunk, atom_indices=atom_indices),
        desc=f"SASA {resname}{resid} ({atom_selection})",
        unit="chunk",
    ):
        sr = md.shrake_rupley(trj, mode="residue", probe_radius=probe_radius_nm, n_sphere_points=n_sphere_points)
        sasa_nm2 = sr[:, res_i].astype(np.float64)

        n = trj.n_frames
        t0 = frame0 * dt_ns * stride
        t_ns.extend(t0 + np.arange(n) * dt_ns * stride)
        sasa_A2.extend(sasa_nm2 * 100.0)  # nm^2 -> Å^2
        frame0 += n

    return np.asarray(t_ns), np.asarray(sasa_A2)


def infer_dt_from_time_axis(time_ns: np.ndarray) -> float | None:
    """
    Infer dt (ns) from a CSV time axis. Robust to tiny float noise.
    Returns None if it can't infer reliably.
    """
    if time_ns is None or len(time_ns) < 5:
        return None
    diffs = np.diff(time_ns)
    diffs = diffs[np.isfinite(diffs)]
    if diffs.size == 0:
        return None
    # take median and round to 3 decimals (handles 0.1000000003)
    dt = float(np.median(diffs))
    if dt <= 0:
        return None
    return dt


def load_csvs(outdir: str, resid: int, resname: str, required_systems=None):
    resname = resname.upper()
    results = {}
    missing = []
    for sysname in systems.keys():
        if required_systems is not None and sysname not in required_systems:
            continue
        fn = f"{sysname.replace(' ', '_')}_{resname}{resid}_sasa.csv"
        path = os.path.join(outdir, fn)
        if not os.path.exists(path):
            missing.append(path)
            continue
        df = pd.read_csv(path)
        # basic column sanity
        if "time_ns" not in df.columns or "sasa_A2" not in df.columns:
            raise ValueError(f"CSV missing required columns time_ns/sasa_A2: {path}")
        results[sysname] = df
    return results, missing


def plot_results(results: dict, outdir: str, resid: int, resname: str, ma_ns: float, fallback_dt_ns: float, fallback_stride: int):
    if not results:
        raise RuntimeError("No results to plot.")

    os.makedirs(outdir, exist_ok=True)

    # Compute MA window per-system using inferred dt from CSV if possible
    ma_windows = {}
    for name, df in results.items():
        dt_eff = infer_dt_from_time_axis(df["time_ns"].values)
        if dt_eff is None:
            dt_eff = fallback_dt_ns * fallback_stride
        w = int(round(ma_ns / dt_eff))
        ma_windows[name] = max(w, 3)

    fig, axes = plt.subplots(len(results), 1, figsize=(12, 10), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)

    for ax, (name, df) in zip(axes, results.items()):
        line_plot(
            xvalues=np.linspace(0,20,len(df["sasa_A2"].values)),
            yvalues=df["sasa_A2"].values,
            ax=ax,
            title=name,
            color=systems[name]["color"],
            z=1,
            label=name,
            ma_window=ma_windows[name],
            ma_color="red",
            ma_z=3,
        )
        ax.set_ylabel(f"{resname.upper()}{resid} SASA (Å$^2$)", fontsize=18)

    axes[-1].set_xlabel("Time (μs)", fontsize=18)
    plt.tight_layout()

    out_pdf = os.path.join(outdir, f"rhoa_{resname.upper()}{resid}_sasa_ma{int(ma_ns)}ns.pdf")
    plt.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved plot: {out_pdf}")


def main():
    ap = argparse.ArgumentParser(description="Compute and/or plot residue SASA timeseries for RhoA sims.")
    ap.add_argument("--plot-only", action="store_true", help="Skip MDTraj SASA calculation and plot from saved CSVs.")
    ap.add_argument("--outdir", default="sasa_out_rhoa", help="Directory for CSVs and plot output.")
    ap.add_argument("--resid", type=int, default=7, help="Residue number (resSeq) to analyze.")
    ap.add_argument("--resname", default="LYS", help="Residue name (e.g., LYS).")
    ap.add_argument("--ma-ns", type=float, default=50.0, help="Moving average window in ns.")
    ap.add_argument("--dt-ns", type=float, default=0.1, help="Frame spacing in ns (fallback for MA window).")
    ap.add_argument("--stride", type=int, default=10, help="Stride used during calculation (fallback for MA window).")
    ap.add_argument("--chunk", type=int, default=5000, help="MDTraj iterload chunk size.")
    ap.add_argument("--atom-selection", default="protein", help='MDTraj selection used for SASA calculation (default "protein").')
    ap.add_argument("--resim-stride-div10", action="store_true",
                    help="If set, uses stride/10 for RhoA-Resim during calculation (keeps your old behavior).")
    args = ap.parse_args()

    outdir = args.outdir
    resid = args.resid
    resname = args.resname.upper()

    if args.plot_only:
        results, missing = load_csvs(outdir, resid, resname)
        if missing:
            print("Warning: missing CSVs (skipping these systems):")
            for m in missing:
                print(f"  {m}")
        plot_results(
            results=results,
            outdir=outdir,
            resid=resid,
            resname=resname,
            ma_ns=args.ma_ns,
            fallback_dt_ns=args.dt_ns,
            fallback_stride=args.stride,
        )
        return

    # ---- compute + save CSVs ----
    os.makedirs(outdir, exist_ok=True)

    results = {}
    for name, meta in systems.items():
        stride = args.stride
        if args.resim_stride_div10 and name == "RhoA-Resim":
            stride = max(1, int(stride / 10))

        t_ns, sasa_A2 = sasa_timeseries_one_residue(
            top_path=meta["top"],
            traj_path=meta["traj"],
            resid=resid,
            resname=resname,
            dt_ns=args.dt_ns,
            stride=stride,
            chunk=args.chunk,
            atom_selection=args.atom_selection,
        )

        df = pd.DataFrame({"time_ns": t_ns, "sasa_A2": sasa_A2, "resid": resid, "resname": resname})
        csv_path = os.path.join(outdir, f"{name.replace(' ', '_')}_{resname}{resid}_sasa.csv")
        df.to_csv(csv_path, index=False)
        print(f"[{name}] wrote {csv_path} (n_frames={len(df)})")
        results[name] = df

    # ---- plot ----
    plot_results(
        results=results,
        outdir=outdir,
        resid=resid,
        resname=resname,
        ma_ns=args.ma_ns,
        fallback_dt_ns=args.dt_ns,
        fallback_stride=args.stride,
    )


if __name__ == "__main__":
    main()