#!/usr/bin/env python3
import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis


# -------------------- Style --------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = "truetype"
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["svg.fonttype"] = "none"


# -------------------- Helpers --------------------
def style_axes(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")


def ensure_dir(d):
    os.makedirs(d, exist_ok=True)


def save_meta(outdir, meta):
    path = os.path.join(outdir, "meta.json")
    with open(path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[cache] wrote {path}")


def load_meta(outdir):
    path = os.path.join(outdir, "meta.json")
    if not os.path.exists(path):
        return None
    with open(path, "r") as f:
        return json.load(f)


def save_df(outdir, df, name):
    path = os.path.join(outdir, name)
    df.to_pickle(path)
    print(f"[cache] wrote {path}")


def load_df(outdir, name):
    path = os.path.join(outdir, name)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return pd.read_pickle(path)


# -------------------- Selections --------------------
DEFAULT_LIPID_ACCEPTOR_NAMES = {
    "POPC": ["P", "O11", "O12", "O13", "O14"],
    "POPE": ["P", "O11", "O12", "O13", "O14"],
    "POPS": ["P", "O11", "O12", "O13", "O14", "O13A", "O13B"],
    # You can add PLA1 here if you know exact acceptor atom names.
    # If not present, we will fall back to POPE names with a warning.
    # "PLA1": ["P", "O11", "O12", "O13", "O14"],
}


def build_name_clause(names):
    return "name " + " ".join(names)


def donor_sel_for_resids(resids):
    res_clause = " ".join(str(int(r)) for r in resids)
    return (
        f"protein and resid {res_clause} and "
        f"((resname ARG and name NH1 NH2) or (resname LYS and name NZ))"
    )


def acceptor_sel_for_lipids(lipid_resnames, lipid_acceptor_names):
    """
    Build an OR selection for acceptors across all lipid types.
    """
    parts = []
    for lip in lipid_resnames:
        names = lipid_acceptor_names.get(lip, None)
        if not names:
            # fallback: use POPE acceptors if present
            if "POPE" in lipid_acceptor_names:
                print(f"[WARNING] No acceptor names configured for '{lip}'. Falling back to POPE names.")
                names = lipid_acceptor_names["POPE"]
            else:
                raise ValueError(f"No acceptor atom names configured for lipid '{lip}' and no POPE fallback.")
        parts.append(f"(resname {lip} and ({build_name_clause(names)}))")
    return " or ".join(parts)


# -------------------- Angle math --------------------
def wrap_delta_deg(a, b):
    return ((a - b + 180.0) % 360.0) - 180.0


def assign_states(theta, phi, time_us,
                  A_theta, A_phi, B_theta, B_phi,
                  tol_deg=15.0,
                  transition_us_start=None, transition_us_end=None):
    """
    Returns labels: A, B, TRANS, UNASSIGNED
    """
    theta = np.asarray(theta, dtype=float)
    phi   = np.asarray(phi, dtype=float)
    t     = np.asarray(time_us, dtype=float)
    n = len(theta)
    if len(phi) != n or len(t) != n:
        raise ValueError("theta, phi, and time_us must have same length.")

    state = np.full(n, "UNASSIGNED", dtype=object)

    inA = (np.abs(theta - A_theta) <= tol_deg) & (np.abs(wrap_delta_deg(phi, A_phi)) <= tol_deg)
    inB = (np.abs(theta - B_theta) <= tol_deg) & (np.abs(wrap_delta_deg(phi, B_phi)) <= tol_deg)

    state[inA] = "A"
    state[inB] = "B"

    if transition_us_start is not None and transition_us_end is not None:
        inT = (t >= float(transition_us_start)) & (t <= float(transition_us_end))
        state[inT] = "TRANS"

    return state


# -------------------- H-bond compute (multi-lipid, cached) --------------------
def compute_hbonds_full(
    top, traj, sim_name, outdir,
    resids=(5, 6, 7),
    lipid_resnames=("POPC", "POPE", "POPS", "PLA1"),
    lipid_acceptor_names=None,
    d_a_cutoff=3.1,
    d_h_a_angle_cutoff=130.0,
    start=0, stop=None, step=1,
    force_recompute=False
):
    ensure_dir(outdir)
    if lipid_acceptor_names is None:
        lipid_acceptor_names = dict(DEFAULT_LIPID_ACCEPTOR_NAMES)

    hb_pkl = "hbonds_full.pkl"
    meta_wanted = {
        "sim_name": sim_name,
        "top": str(top),
        "traj": str(traj),
        "resids": [int(r) for r in resids],
        "lipid_resnames": list(lipid_resnames),
        "d_a_cutoff": float(d_a_cutoff),
        "d_h_a_angle_cutoff": float(d_h_a_angle_cutoff),
        "start": int(start) if start is not None else None,
        "stop": int(stop) if stop is not None else None,
        "step": int(step),
    }

    if (not force_recompute) and os.path.exists(os.path.join(outdir, hb_pkl)):
        meta_old = load_meta(outdir)
        if meta_old is not None:
            same = True
            for k, v in meta_wanted.items():
                if meta_old.get(k, None) != v:
                    same = False
                    break
            if same:
                print("[cache] loading cached hbonds_full.pkl")
                return load_df(outdir, hb_pkl), meta_old

    u = mda.Universe(top, traj)
    n_frames = u.trajectory.n_frames
    if stop is None:
        stop = n_frames

    donors_sel = donor_sel_for_resids(resids)
    acceptors_sel = acceptor_sel_for_lipids(lipid_resnames, lipid_acceptor_names)

    if u.select_atoms(donors_sel).n_atoms == 0:
        raise ValueError(f"[{sim_name}] Donor selection matched 0 atoms:\n{donors_sel}")

    h = HydrogenBondAnalysis(
        universe=u,
        donors_sel=donors_sel,
        acceptors_sel=acceptors_sel,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=True
    )

    print(f"[{sim_name}] Running full-trajectory H-bonds: frames {start}:{stop} step={step}")
    try:
        h.run(start=start, stop=stop, step=step, verbose=True)
    except TypeError:
        h.run(start=start, stop=stop, step=step)

    hb = h.results.hbonds

    try:
        sampled_frames = np.array(h.frames, dtype=int)
    except Exception:
        sampled_frames = np.arange(start, stop, step, dtype=int)
    sampled_frames = np.sort(sampled_frames)

    df = pd.DataFrame({"frame": sampled_frames})

    # per-res + per-lipid, plus per-lipid sums and overall sum_all
    lipid_resnames = tuple(lipid_resnames)
    for rid in resids:
        for lip in lipid_resnames:
            df[f"res{int(rid)}_{lip}"] = 0
    for lip in lipid_resnames:
        df[f"sum_{lip}"] = 0
    df["sum_all"] = 0

    if hb is not None and len(hb) > 0:
        donor_atom_indices = hb[:, 1].astype(int)
        acceptor_atom_indices = hb[:, 3].astype(int)
        frames_hb = hb[:, 0].astype(int)

        donors_atoms = u.atoms[donor_atom_indices]
        acceptor_atoms = u.atoms[acceptor_atom_indices]

        donor_resids = donors_atoms.resids.astype(int)
        acceptor_lipids = acceptor_atoms.resnames

        frame_to_row = {int(f): i for i, f in enumerate(sampled_frames)}
        resid_set = set(int(r) for r in resids)
        lipid_set = set(lipid_resnames)

        for f, rid, lip in zip(frames_hb, donor_resids, acceptor_lipids):
            rid = int(rid)
            lip = str(lip)
            if rid not in resid_set:
                continue
            if lip not in lipid_set:
                continue
            i = frame_to_row.get(int(f), None)
            if i is None:
                continue
            df.iat[i, df.columns.get_loc(f"res{rid}_{lip}")] += 1

    # lipid sums + total
    for lip in lipid_resnames:
        df[f"sum_{lip}"] = df[[f"res{int(r)}_{lip}" for r in resids]].sum(axis=1)
    df["sum_all"] = df[[f"sum_{lip}" for lip in lipid_resnames]].sum(axis=1)

    meta_used = dict(meta_wanted)
    meta_used.update({"n_frames": int(n_frames)})

    save_meta(outdir, meta_used)
    save_df(outdir, df, hb_pkl)
    df.to_csv(os.path.join(outdir, "hbonds_full.csv"), index=False)
    print(f"[out] wrote {os.path.join(outdir, 'hbonds_full.csv')}")

    return df, meta_used


# -------------------- Plot: mean hbonds per state heatmap (multi-lipid + TOTAL) --------------------
def plot_mean_hbonds_per_state_heatmap_multi(
    df: pd.DataFrame,
    lipid_resnames,
    outpath: str,
    state_col: str = "state",
    include_unassigned: bool = False,
    # order + label mapping you requested
    state_order=("A", "TRANS", "B"),
    state_label_map=None,
    figsize=(5.2, 2.1),
    cmap="viridis",
    annotate=True,
    min_frames_per_state=1,
):
    """
    Heatmap:
      rows = states (Start / Transition / OS_0), in order A, TRANS, B (top->bottom)
      cols = lipids + TOTAL
      color = mean hbonds/frame
      colorbar label: "Mean hbonds/frame" on right
      no title
    """
    if state_label_map is None:
        state_label_map = {"A": "Start", "TRANS": "Transition", "B": "OS_0", "UNASSIGNED": "UNASSIGNED"}

    lipid_resnames = list(lipid_resnames)
    cols = [f"sum_{lip}" for lip in lipid_resnames] + ["sum_all"]

    for c in cols:
        if c not in df.columns:
            raise ValueError(f"df missing required column '{c}'")

    if state_col not in df.columns:
        raise ValueError(f"df missing required column '{state_col}'")

    states = list(state_order)
    if include_unassigned and "UNASSIGNED" not in states:
        states.append("UNASSIGNED")

    # compute means
    M = np.full((len(states), len(cols)), np.nan, dtype=float)
    N = np.zeros((len(states),), dtype=int)

    for i, st in enumerate(states):
        sub = df.loc[df[state_col] == st]
        N[i] = int(len(sub))
        if N[i] < int(min_frames_per_state):
            continue
        for j, c in enumerate(cols):
            M[i, j] = float(np.nanmean(sub[c].to_numpy(dtype=float)))

    # write CSV summary (means + n)
    out_csv = os.path.splitext(outpath)[0] + "_summary.csv"
    recs = []
    for i, st in enumerate(states):
        row = {"state": st, "state_label": state_label_map.get(st, st), "n": int(N[i])}
        for lip in lipid_resnames:
            row[lip] = M[i, cols.index(f"sum_{lip}")]
        row["TOTAL"] = M[i, cols.index("sum_all")]
        recs.append(row)
    pd.DataFrame(recs).to_csv(out_csv, index=False)
    print(f"[out] wrote {out_csv}")

    # plot
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(M, aspect="auto", cmap=cmap)

    # ticks/labels
    ylabels = [state_label_map.get(st, st) for st in states]
    xlabels = lipid_resnames + ["TOTAL"]

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=10)

    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_xticklabels(xlabels, fontsize=10, rotation=0)

    # no title (as requested)
    # annotate values
    if annotate:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=9, color="black")

    # clean spines/ticks
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(width=1.0, length=0)

    # colorbar on right
    cbar = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cbar.set_label("Mean hbonds/frame", fontsize=10)
    cbar.ax.tick_params(labelsize=9, width=1.0, length=3)
    cbar.outline.set_linewidth(1.0)

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"[out] wrote {outpath}")


# -------------------- Main --------------------
def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--top", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--name", default="Resim")
    ap.add_argument("--outdir", default="transition_hbonds")

    # hbonds parameters
    ap.add_argument("--resids", nargs="+", type=int, default=[5, 6, 7])
    ap.add_argument("--lipids", nargs="+", default=["POPC", "POPE", "POPS", "PLA1"],
                    help="Lipid resnames to include (columns will be created for each + TOTAL).")
    ap.add_argument("--da", type=float, default=3.1)
    ap.add_argument("--angle", type=float, default=130.0)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=None)
    ap.add_argument("--step", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--plot-only", action="store_true")

    # orientation arrays
    ap.add_argument("--theta-txt", required=True, help="Tilt array (degrees), one per frame, txt.")
    ap.add_argument("--phi-txt", required=True, help="Rotation array (degrees), one per frame, txt.")
    ap.add_argument("--total-us", type=float, default=20.0, help="Total trajectory duration for linspace time axis.")

    # state centers + tolerance
    ap.add_argument("--A-theta", type=float, required=True)
    ap.add_argument("--A-phi", type=float, required=True)
    ap.add_argument("--B-theta", type=float, required=True)
    ap.add_argument("--B-phi", type=float, required=True)
    ap.add_argument("--tol-deg", type=float, default=15.0)

    # optional transition window (µs)
    ap.add_argument("--transition-us-start", type=float, default=None)
    ap.add_argument("--transition-us-end", type=float, default=None)

    args = ap.parse_args()
    ensure_dir(args.outdir)

    lipid_acceptor_names = dict(DEFAULT_LIPID_ACCEPTOR_NAMES)

    # 1) Hbonds (cached)
    if args.plot_only:
        print("[mode] plot-only: loading cached hbonds_full.pkl")
        hb_df = load_df(args.outdir, "hbonds_full.pkl")
    else:
        hb_df, _meta = compute_hbonds_full(
            top=args.top,
            traj=args.traj,
            sim_name=args.name,
            outdir=args.outdir,
            resids=tuple(args.resids),
            lipid_resnames=tuple(args.lipids),
            lipid_acceptor_names=lipid_acceptor_names,
            d_a_cutoff=args.da,
            d_h_a_angle_cutoff=args.angle,
            start=args.start,
            stop=args.stop,
            step=args.step,
            force_recompute=args.force,
        )

    # 2) Load angles
    theta = np.loadtxt(args.theta_txt, dtype=float)
    phi   = np.loadtxt(args.phi_txt, dtype=float)

    n_ang = len(theta)
    if len(phi) != n_ang:
        raise ValueError("theta and phi arrays have different lengths.")
    time_us = np.linspace(0.0, float(args.total_us), num=n_ang, endpoint=True)

    frames = hb_df["frame"].to_numpy(dtype=int)
    if np.max(frames) >= n_ang:
        raise ValueError(
            f"Angle arrays length ({n_ang}) is shorter than max frame in hbonds ({np.max(frames)}). "
            f"Make sure your angle arrays correspond to the same trajectory frames."
        )

    hb_df = hb_df.copy()
    hb_df["theta_deg"] = theta[frames]
    hb_df["phi_deg"]   = phi[frames]
    hb_df["time_us"]   = time_us[frames]

    # 3) Assign states
    hb_df["state"] = assign_states(
        theta=hb_df["theta_deg"].to_numpy(),
        phi=hb_df["phi_deg"].to_numpy(),
        time_us=hb_df["time_us"].to_numpy(),
        A_theta=args.A_theta, A_phi=args.A_phi,
        B_theta=args.B_theta, B_phi=args.B_phi,
        tol_deg=args.tol_deg,
        transition_us_start=args.transition_us_start,
        transition_us_end=args.transition_us_end,
    )

    # Save annotated table
    hb_df.to_csv(os.path.join(args.outdir, "hbonds_full_with_states.csv"), index=False)
    save_df(args.outdir, hb_df, "hbonds_full_with_states.pkl")

    # 4) Mean hbonds/state heatmap: lipids + TOTAL, rows Start/Transition/OS_0
    out_heat = os.path.join(args.outdir, f"{args.name}_mean_hbonds_by_state_heatmap.svg")
    plot_mean_hbonds_per_state_heatmap_multi(
        df=hb_df,
        lipid_resnames=args.lipids,
        outpath=out_heat,
        state_col="state",
        include_unassigned=False,
        state_order=("A", "TRANS", "B"),  # top->bottom
        state_label_map={"A": "Start", "TRANS": "Transition", "B": "OS_0", "UNASSIGNED": "UNASSIGNED"},
        figsize=(5.4, 2.15),
        cmap="viridis",
        annotate=True,
        min_frames_per_state=1,
    )


if __name__ == "__main__":
    main()
