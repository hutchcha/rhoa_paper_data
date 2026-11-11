#!/usr/bin/env python3
import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from MDAnalysis.analysis.rdf import InterRDF

# -------------------- Style --------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = "truetype"
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["svg.fonttype"] = "none"

# -------------------- OI palette --------------------
OI_PALETTE = [
    "#1F77B4",  # blue
    "#FF7F0E",  # orange
    "#2CA02C",  # green
    "#D62728",  # red
    "#9467BD",  # purple
    "#8C564B",  # brown
]

# -------------------- Analysis window controls --------------------
# For grouped stacked occupancy, RDFs, and heatmaps:
# use only the last 3/4 of the trajectory/data.
PLOT_KEEP_LAST_FRACTION = 0.75

def style_axes(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")

def moving_average(y, window):
    y = pd.Series(np.asarray(y, dtype=float))
    return y.rolling(window=window, center=True, min_periods=max(1, window//5)).mean().to_numpy()

def line_plot_raw_ma(ax, x, y, label, color, ma_window):
    ax.plot(x, y, alpha=0.25, color=color, linewidth=1.0, label=None)
    ax.plot(x, moving_average(y, ma_window), alpha=1.0, color=color, linewidth=2.0, label=label)

def frame_start_for_last_fraction(n_frames: int, fraction: float = PLOT_KEEP_LAST_FRACTION) -> int:
    """
    Return first frame index corresponding to keeping the last `fraction`
    of the trajectory.
    Example: fraction=0.75 -> start at n_frames//4.
    """
    if not (0 < float(fraction) <= 1.0):
        raise ValueError(f"fraction must be in (0, 1], got {fraction}")
    keep_count = int(np.ceil(float(n_frames) * float(fraction)))
    return max(0, int(n_frames) - keep_count)

def slice_last_fraction_df(df: pd.DataFrame, fraction: float = PLOT_KEEP_LAST_FRACTION, zero_time: bool = False) -> pd.DataFrame:
    """
    Keep only the last `fraction` of a dataframe.
    Optionally re-zero time_us.
    """
    start_idx = frame_start_for_last_fraction(len(df), fraction)
    out = df.iloc[start_idx:].copy().reset_index(drop=True)

    if zero_time and "time_us" in out.columns and len(out) > 0:
        out["time_us"] = out["time_us"] - float(out["time_us"].iloc[0])

    return out

def make_time_blocks_from_window(start_frame: int, stop_frame: int, n_blocks: int):
    """
    Split [start_frame, stop_frame) into n_blocks contiguous blocks.
    """
    if n_blocks < 1:
        raise ValueError("n_blocks must be >= 1")
    if stop_frame <= start_frame:
        raise ValueError(f"Invalid block window: start={start_frame}, stop={stop_frame}")

    edges = np.linspace(start_frame, stop_frame, n_blocks + 1)
    edges = np.round(edges).astype(int)
    edges[0] = int(start_frame)
    edges[-1] = int(stop_frame)

    blocks = []
    for i in range(n_blocks):
        s, e = int(edges[i]), int(edges[i + 1])
        if e <= s:
            e = min(int(stop_frame), s + 1)
        blocks.append((s, e))
    return blocks

# -------------------- Residue label helpers --------------------
def one_letter_from_name(resname: str) -> str:
    resname = (resname or "").upper()
    m = {
        "ARG": "R", "LYS": "K", "HIS": "H",
        "ASP": "D", "GLU": "E",
        "SER": "S", "THR": "T",
        "ASN": "N", "GLN": "Q",
        "ALA": "A", "VAL": "V", "LEU": "L", "ILE": "I",
        "GLY": "G", "PRO": "P",
        "TYR": "Y", "PHE": "F", "TRP": "W",
        "CYS": "C", "MET": "M",
    }
    return m.get(resname, "?")

def get_residue_label_map_from_universe(top, traj, resids):
    u = mda.Universe(top, traj)
    out = {}
    for rid in resids:
        ag = u.select_atoms(f"protein and resid {int(rid)}")
        if ag.n_atoms == 0:
            out[int(rid)] = f"res{int(rid)}"
            continue
        aa3 = ag.residues[0].resname
        out[int(rid)] = f"{one_letter_from_name(aa3)}{int(rid)}"
    return out

def parse_residue_labels_arg(s: str, resids):
    if not s:
        return None
    out = {}
    parts = [p.strip() for p in s.split(",") if p.strip()]
    for p in parts:
        rid_str, aa = [x.strip() for x in p.split(":")]
        rid = int(rid_str)
        out[rid] = f"{aa}{rid}"
    for rid in resids:
        out.setdefault(int(rid), f"res{int(rid)}")
    return out

# -------------------- Colors/styles --------------------
def make_lipid_color_map(lipid_resnames):
    colors = {}
    for i, lip in enumerate(lipid_resnames):
        colors[lip] = OI_PALETTE[i % len(OI_PALETTE)]
    return colors

def make_resid_styles(resids):
    """Make residue distinctions obvious in RDF via dash patterns + one marker series + TOTAL thicker."""
    rlist = [int(r) for r in resids]
    styles = {}

    dash_solid = "-"
    dash_long  = (0, (7, 2.5))
    dash_dot   = (0, (1.2, 2.2))
    dash_dashdot = (0, (6, 2, 1.5, 2))

    if len(rlist) >= 1:
        styles[rlist[0]] = dict(linestyle=dash_solid, lw=2.2, marker=None, markevery=None, ms=None)
    if len(rlist) >= 2:
        styles[rlist[1]] = dict(linestyle=dash_long,  lw=2.5, marker=None, markevery=None, ms=None)
    if len(rlist) >= 3:
        styles[rlist[2]] = dict(linestyle=dash_dot,   lw=2.2, marker="o", markevery=18, ms=3.2)

    cycle = [dash_solid, dash_long, dash_dot]
    for i, rid in enumerate(rlist):
        styles.setdefault(rid, dict(linestyle=cycle[i % len(cycle)], lw=2.2, marker=None, markevery=None, ms=None))

    styles["TOTAL"] = dict(linestyle=dash_dashdot, lw=3.2, marker=None, markevery=None, ms=None)
    return styles

# -------------------- Selections --------------------
def donor_sel_for_resids(resids):
    res_clause = " ".join(str(r) for r in resids)
    return f"protein and resid {res_clause} and ((resname ARG and name NH1 NH2) or (resname LYS and name NZ))"

DEFAULT_LIPID_ACCEPTOR_NAMES = {
    "POPC": ["O11","O12", "O13", "O14"],
    "POPE": ["O11","O12", "O13", "O14"],
    "POPS": ["O11","O12", "O13", "O14"],
}

DEFAULT_LIPID_P_NAMES = {
    "POPC": ["P"],
    "POPE": ["P"],
    "POPS": ["P"],
}

def build_name_clause(names):
    return "name " + " ".join(names)

def acceptor_sel_for_lipids(lipid_resnames, lipid_acceptor_names):
    parts = []
    for lip in lipid_resnames:
        names = lipid_acceptor_names.get(lip, None)
        if not names:
            raise ValueError(f"No acceptor atom names configured for lipid '{lip}'.")
        parts.append(f"(resname {lip} and ({build_name_clause(names)}))")
    return " or ".join(parts)

# -------------------- Cache helpers --------------------
def _ensure_dir(d):
    if d is None:
        return
    os.makedirs(d, exist_ok=True)

def save_meta(cache_dir, meta: dict):
    _ensure_dir(cache_dir)
    path = os.path.join(cache_dir, "meta.json")
    with open(path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[cache] wrote {path}")

def save_hbonds_pickle(cache_dir, sim_name, df: pd.DataFrame):
    _ensure_dir(cache_dir)
    path = os.path.join(cache_dir, f"hbonds_{sim_name}.pkl")
    df.to_pickle(path)
    print(f"[cache] wrote {path}")

def load_cached_hbonds(load_dir, sim_names):
    dfs = []
    for sim in sim_names:
        path = os.path.join(load_dir, f"hbonds_{sim}.pkl")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing cached hbonds pickle: {path}")
        dfs.append(pd.read_pickle(path))
    return dfs

def save_rdf_pickle(cache_dir, rdf_dict: dict, fname="rdf_results.pkl"):
    _ensure_dir(cache_dir)
    path = os.path.join(cache_dir, fname)
    pd.to_pickle(rdf_dict, path)
    print(f"[cache] wrote {path}")

def load_cached_rdfs(load_dir, fname="rdf_results.pkl"):
    path = os.path.join(load_dir, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing cached rdf pickle: {path}")
    return pd.read_pickle(path)

def save_hbonds_blocks_pickle(cache_dir, hb_blocks: dict, fname="hbonds_blocks.pkl"):
    _ensure_dir(cache_dir)
    path = os.path.join(cache_dir, fname)
    pd.to_pickle(hb_blocks, path)
    print(f"[cache] wrote {path}")

def load_hbonds_blocks_pickle(load_dir, fname="hbonds_blocks.pkl"):
    path = os.path.join(load_dir, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing cached hbonds blocks pickle: {path}")
    return pd.read_pickle(path)

# -------------------- H-bond analysis + aggregation --------------------
import inspect

def _build_resid_group_sel(resids):
    res_clause = " ".join(str(int(r)) for r in resids)
    return f"protein and resid {res_clause}"

def _build_lipid_group_sel(lipid_resnames):
    return " or ".join([f"resname {lip}" for lip in lipid_resnames])

def _make_hba_with_guessing(u, donors_group_sel, acceptors_group_sel,
                           d_a_cutoff, d_h_a_angle_cutoff):
    """
    Create HydrogenBondAnalysis in a version-tolerant way, preferring built-in guessing.
    """
    kwargs = dict(
        universe=u,
        donors_sel=donors_group_sel,
        acceptors_sel=acceptors_group_sel,
        d_a_cutoff=float(d_a_cutoff),
        d_h_a_angle_cutoff=float(d_h_a_angle_cutoff),
        update_selections=True,
    )

    sig = None
    try:
        sig = inspect.signature(HydrogenBondAnalysis.__init__)
    except Exception:
        sig = None

    if sig is not None:
        params = set(sig.parameters.keys())
        for k, v in [
            ("guess_hydrogens", True),
            ("guess_acceptors", True),
            ("guess_donors", True),
        ]:
            if k in params:
                kwargs[k] = v

        if ("between" in params) and ("donors_sel" not in params or "acceptors_sel" not in params):
            kwargs.pop("donors_sel", None)
            kwargs.pop("acceptors_sel", None)
            kwargs["between"] = [donors_group_sel, acceptors_group_sel]

    h = HydrogenBondAnalysis(**kwargs)

    for meth_name in ("guess_donors", "guess_acceptors", "guess_hydrogens"):
        if hasattr(h, meth_name):
            try:
                getattr(h, meth_name)()
            except TypeError:
                try:
                    getattr(h, meth_name)(donors_group_sel if "donor" in meth_name or "hydrogen" in meth_name else acceptors_group_sel)
                except Exception:
                    pass

    return h



def plot_rdf_single_residue_all_sims(
    rdf_dict,
    sim_names,
    resid,
    lipid_resnames,
    outpath,
    lipid_colors,
    resid_label_map=None,
    sim_linestyles=None,
    figsize=(5.2, 7.2),
    linewidth=2.0,
):
    """
    Plot RDFs for one residue across all simulations as stacked panels.

    New behavior:
      - one subplot per simulation
      - color = lipid
      - shared y-axis across all panels

    Signature is unchanged so the existing call in main() still works.
    """
    resid = int(resid)

    n_sims = len(sim_names)
    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=figsize,
        sharex=True,
        sharey=True,
    )
    if n_sims == 1:
        axes = [axes]

    found_any = False
    global_ymax = 0.0

    # First pass: determine global ymax so all panels truly share the same scale
    for sim in sim_names:
        for lip in lipid_resnames:
            key = (sim, resid, lip)
            if key not in rdf_dict:
                continue
            rdf = np.asarray(rdf_dict[key]["rdf"], dtype=float)
            if rdf.size > 0 and np.any(np.isfinite(rdf)):
                global_ymax = max(global_ymax, float(np.nanmax(rdf)))
                found_any = True

    if not found_any:
        raise ValueError(f"No RDF traces found for residue {resid} across sims={sim_names} and lipids={lipid_resnames}")

    # Plot
    for ax, sim in zip(axes, sim_names):
        sim_found = False

        for lip in lipid_resnames:
            key = (sim, resid, lip)
            if key not in rdf_dict:
                continue

            bins = np.asarray(rdf_dict[key]["bins"], dtype=float)
            rdf  = np.asarray(rdf_dict[key]["rdf"], dtype=float)

            ax.plot(
                bins, rdf,
                color=lipid_colors[lip],
                linestyle="-",
                linewidth=linewidth,
                label=lip if ax is axes[0] else None,
            )
            sim_found = True

        style_axes(ax)
        ax.tick_params(labelsize=16)
        ax.set_ylabel("g(r)", fontsize=18)

        ax.set_title(sim, fontsize=15)

        if global_ymax > 0:
            ax.set_ylim(0, global_ymax * 1.05)

        if not sim_found:
            ax.text(
                0.5, 0.5, "No data",
                transform=ax.transAxes,
                ha="center", va="center",
                fontsize=11
            )

    axes[-1].set_xlabel("r (Å)", fontsize=18)

    # Figure title from residue label
    if resid_label_map is not None and resid in resid_label_map:
        fig_title = resid_label_map[resid]
    else:
        fig_title = f"Residue {resid}"
    fig.suptitle(fig_title, fontsize=14, y=0.995)

    # Single legend on top panel only
    lipid_handles = [
        plt.Line2D([0], [0], color=lipid_colors[lip], lw=3, label=lip)
        for lip in lipid_resnames
    ]
    axes[0].legend(handles=lipid_handles, frameon=False, fontsize=10, loc="upper right", title="Lipid")

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")


def plot_hbond_timeseries_stacked_three(
    dfs, sim_names, lipid_resnames, outpath, lipid_colors, ma_window=500,
    figsize=(8, 6)
):
    """
    Stack three H-bond timeseries panels vertically in the order provided.
    Intended use:
        GDP, GTP, Resim
    with Resim last.

    Keeps the same plotting logic/style as plot_hbond_timeseries().
    """
    if len(dfs) != 3 or len(sim_names) != 3:
        raise ValueError("plot_hbond_timeseries_stacked_three requires exactly 3 dfs and 3 sim_names.")

    fig, axes = plt.subplots(
        nrows=3, ncols=1,
        figsize=figsize,
        sharex=True, sharey=False
    )

    for ax, df, sim in zip(axes, dfs, sim_names):
        t = np.linspace(0,20,len(df["time_us"]))
        for lip in lipid_resnames:
            col = f"sum_{lip}"
            if col not in df.columns:
                continue
            y = df[col].to_numpy(dtype=float)
            line_plot_raw_ma(ax, t, y, label=lip, color=lipid_colors[lip], ma_window=ma_window)

        ax.set_title(sim, fontsize=22)
        ax.set_ylabel("Hbond occupancy", fontsize=18)
        ax.tick_params(labelsize=16)
        style_axes(ax)
        ax.legend(frameon=False, fontsize=11, loc="upper right")

    axes[-1].set_xlabel("Time (µs)", fontsize=18)

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")

def run_hbonds_and_aggregate(
    top, traj, sim_name,
    resids=(5,6,7),
    lipid_resnames=("POPC","POPE","POPS"),
    dt_ns=1.0,
    d_a_cutoff=3.0,
    d_h_a_angle_cutoff=135.0,
    start=None, stop=None, step=None,
    outdir="hbonds_out",
    write_csv=True,
    resim_stride_factor=1,
):
    """
    Computes *frequency* (0/1) of H-bonding per donor atom per frame.
    """
    os.makedirs(outdir, exist_ok=True)

    u = mda.Universe(top, traj)

    dt_ns_eff = float(dt_ns)
    try:
        if getattr(u.trajectory, "dt", None) is not None:
            dt_ns_eff = float(u.trajectory.dt) / 1000.0
    except Exception:
        pass

    def _stride_fix(x):
        if x is None:
            return None
        x = int(x)
        if resim_stride_factor is None or int(resim_stride_factor) <= 1:
            return x
        f = int(resim_stride_factor)
        return max(0, x // f)

    start_eff = _stride_fix(start)
    stop_eff  = _stride_fix(stop)
    step_eff  = _stride_fix(step) if step is not None else None
    if step_eff == 0:
        step_eff = 1

    donors_group_sel   = _build_resid_group_sel(resids)
    acceptors_group_sel = _build_lipid_group_sel(lipid_resnames)

    if u.select_atoms(donors_group_sel).n_atoms == 0:
        raise ValueError(f"[{sim_name}] Donor GROUP selection matched 0 atoms:\n{donors_group_sel}")
    if u.select_atoms(acceptors_group_sel).n_atoms == 0:
        raise ValueError(f"[{sim_name}] Acceptor GROUP selection matched 0 atoms:\n{acceptors_group_sel}")

    h = _make_hba_with_guessing(
        u=u,
        donors_group_sel=donors_group_sel,
        acceptors_group_sel=acceptors_group_sel,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=d_h_a_angle_cutoff
    )

    print(f"\n[{sim_name}] Running H-bond analysis (guess donors/acceptors)...")
    try:
        h.run(start=start_eff, stop=stop_eff, step=step_eff, verbose=True)
    except TypeError:
        h.run(start=start_eff, stop=stop_eff, step=step_eff)

    hb = h.results.hbonds
    n_frames = u.trajectory.n_frames

    try:
        sampled_frames = np.array(h.frames, dtype=int)
    except Exception:
        s0 = 0 if start_eff is None else int(start_eff)
        e0 = n_frames if stop_eff is None else int(stop_eff)
        st = 1 if step_eff is None else int(step_eff)
        sampled_frames = np.arange(s0, e0, st, dtype=int)

    sampled_frames = np.sort(sampled_frames)
    df = pd.DataFrame({"frame": sampled_frames})
    df["time_us"] = (df["frame"].to_numpy(dtype=float) * dt_ns_eff) / 1000.0

    for rid in resids:
        for lip in lipid_resnames:
            df[f"res{int(rid)}_{lip}"] = 0

    if hb is None or len(hb) == 0:
        print(f"[{sim_name}] No hbonds detected under given cutoffs.")
        for lip in lipid_resnames:
            df[f"sum_{lip}"] = 0
        df["sum_all"] = 0
        if write_csv:
            out_csv = os.path.join(outdir, f"{sim_name}_hbonds.csv")
            df.to_csv(out_csv, index=False)
            print(f"[{sim_name}] Wrote {out_csv}")
        n_donors = float(len(resids))
        for lip in lipid_resnames:
            df[f"sum_{lip}_occ"] = df[f"sum_{lip}"] / n_donors
        df["sum_all_occ"] = df["sum_all"] / n_donors
        return df

    frames = hb[:, 0].astype(int)
    donor_atom_indices = hb[:, 1].astype(int)
    acceptor_atom_indices = hb[:, 3].astype(int)
    dist = hb[:, 4].astype(float) if hb.shape[1] >= 5 else np.zeros(len(hb), dtype=float)

    tmp = pd.DataFrame({
        "frame": frames,
        "donor_idx": donor_atom_indices,
        "acceptor_idx": acceptor_atom_indices,
        "dist": dist
    }).sort_values("dist", kind="mergesort").drop_duplicates(subset=["frame", "donor_idx"], keep="first")

    donors_atoms = u.atoms[tmp["donor_idx"].to_numpy(dtype=int)]
    acceptor_atoms = u.atoms[tmp["acceptor_idx"].to_numpy(dtype=int)]

    donor_resids = donors_atoms.resids.astype(int)
    acceptor_lipids = acceptor_atoms.resnames

    lipid_set = set(lipid_resnames)
    resid_set = set(int(r) for r in resids)
    frame_to_row = {int(f): i for i, f in enumerate(sampled_frames)}

    for f, rid, lip in zip(tmp["frame"].to_numpy(dtype=int), donor_resids, acceptor_lipids):
        if rid not in resid_set:
            continue
        if lip not in lipid_set:
            continue
        i = frame_to_row.get(int(f), None)
        if i is None:
            continue
        df.iat[i, df.columns.get_loc(f"res{int(rid)}_{lip}")] = 1

    for lip in lipid_resnames:
        df[f"sum_{lip}"] = df[[f"res{int(rid)}_{lip}" for rid in resids]].sum(axis=1)

    per_donor_any = []
    for rid in resids:
        cols = [f"res{int(rid)}_{lip}" for lip in lipid_resnames]
        per_donor_any.append(df[cols].max(axis=1))
    df["sum_all"] = pd.concat(per_donor_any, axis=1).sum(axis=1) / float(len(resids))

    if write_csv:
        out_csv = os.path.join(outdir, f"{sim_name}_hbonds.csv")
        df.to_csv(out_csv, index=False)
        print(f"[{sim_name}] Wrote {out_csv}")

    return df

# -------------------- Time-blocking --------------------
def make_time_blocks_from_nframes(n_frames: int, n_blocks: int):
    """
    Split [0, n_frames) into n_blocks contiguous blocks.
    Returns list of (start, stop) with stop exclusive.
    """
    if n_blocks < 1:
        raise ValueError("n_blocks must be >= 1")
    edges = np.linspace(0, n_frames, n_blocks + 1)
    edges = np.round(edges).astype(int)
    edges[0] = 0
    edges[-1] = n_frames
    blocks = []
    for i in range(n_blocks):
        s, e = int(edges[i]), int(edges[i+1])
        if e <= s:
            e = min(n_frames, s + 1)
        blocks.append((s, e))
    return blocks

def compute_hbond_blocks_from_df(df: pd.DataFrame, sim_name: str, n_frames: int, n_blocks: int,
                                 keep_last_fraction: float = PLOT_KEEP_LAST_FRACTION):
    """
    Slice an existing hbonds df into time blocks using only the last `keep_last_fraction`
    of the trajectory.
    Returns dict: (sim_name, block_index) -> df_slice
    """
    start_frame = frame_start_for_last_fraction(n_frames, keep_last_fraction)
    stop_frame = int(n_frames)

    blocks = make_time_blocks_from_window(start_frame, stop_frame, n_blocks)
    out = {}
    frames = df["frame"].to_numpy(dtype=int)

    for bi, (s, e) in enumerate(blocks, start=1):
        mask = (frames >= s) & (frames < e)
        dfi = df.loc[mask].copy()
        out[(sim_name, bi)] = dfi

    return out, blocks

# -------------------- Plotting: timeseries --------------------
def plot_hbond_timeseries(dfs, sim_names, lipid_resnames, outpath, lipid_colors, ma_window=500):
    n = len(dfs)
    fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(8.5, 2.8 * n), sharex=True, sharey=False)
    if n == 1:
        axes = [axes]

    for ax, df, sim in zip(axes, dfs, sim_names):
        t = df["time_us"].to_numpy()
        for lip in lipid_resnames:
            col = f"sum_{lip}"
            if col not in df.columns:
                continue
            y = df[col].to_numpy(dtype=float)
            line_plot_raw_ma(ax, t, y, label=lip, color=lipid_colors[lip], ma_window=ma_window)
        ax.set_ylim(0,10)
        ax.set_title(sim, fontsize=18)
        ax.set_ylabel("H-bond occupancy", fontsize=14)
        ax.tick_params(labelsize=12)
        style_axes(ax)
        ax.legend(frameon=False, fontsize=11, loc="upper right")

    axes[-1].set_xlabel("Time (µs)", fontsize=14)
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)

def block_mean_sem_from_df(df: pd.DataFrame, cols, n_blocks: int = 5):
    if n_blocks is None or int(n_blocks) <= 1:
        arr = df[cols].to_numpy(dtype=float)
        mean = np.nanmean(arr, axis=0)
        sem = np.zeros_like(mean)
        return mean, sem, 1

    n_blocks = int(n_blocks)
    idx_blocks = np.array_split(np.arange(len(df)), n_blocks)

    block_means = []
    for idx in idx_blocks:
        if idx.size == 0:
            continue
        block_arr = df.iloc[idx][cols].to_numpy(dtype=float)
        block_means.append(np.nanmean(block_arr, axis=0))

    if not block_means:
        mean = np.zeros(len(cols), dtype=float)
        sem  = np.zeros(len(cols), dtype=float)
        return mean, sem, 0

    bm = np.vstack(block_means)
    mean = bm.mean(axis=0)
    if bm.shape[0] > 1:
        sem = bm.std(axis=0, ddof=1) / np.sqrt(bm.shape[0])
    else:
        sem = np.zeros_like(mean)

    return mean, sem, bm.shape[0]

# -------------------- Plotting: stacked bar summary --------------------
def plot_grouped_stacked_summary(
    dfs, sim_names, resids, lipid_resnames, outpath, lipid_colors,
    resid_label_map,
    label_threshold_pct=6.0,
    sem_blocks: int = 5,
):
    # TOTAL removed per request
    res_keys   = list(resids)
    res_labels = [resid_label_map[int(r)] for r in resids]

    n_sims = len(dfs)
    n_bars = len(res_keys)
    n_lip  = len(lipid_resnames)

    means = np.zeros((n_sims, n_bars, n_lip), dtype=float)
    sem_total_height = np.zeros((n_sims, n_bars), dtype=float)

    for si, df in enumerate(dfs):
        for bi, key in enumerate(res_keys):
            lipid_cols = []
            for lip in lipid_resnames:
                col = f"res{int(key)}_{lip}"
                lipid_cols.append(col)

            df_work = df.copy()
            for col in lipid_cols:
                if col not in df_work.columns:
                    df_work[col] = 0.0

            mean_seg, sem_seg, _ = block_mean_sem_from_df(df_work, lipid_cols, n_blocks=sem_blocks)

            for li, lip in enumerate(lipid_resnames):
                means[si, bi, li] = float(mean_seg[li])

            df_work["_stack_sum"] = df_work[lipid_cols].sum(axis=1)
            mean_sum, sem_sum, _ = block_mean_sem_from_df(df_work, ["_stack_sum"], n_blocks=sem_blocks)
            sem_total_height[si, bi] = float(sem_sum[0])

    fig, ax = plt.subplots(figsize=(5.4, 6.2))
    style_axes(ax)

    group_x = np.arange(n_sims)
    bar_w   = 0.155
    offsets = (np.arange(n_bars) - (n_bars - 1) / 2.0) * (bar_w * 1.22)
    global_max = 0.0

    for bi, key in enumerate(res_keys):
        x = group_x + offsets[bi]
        bottom = np.zeros(n_sims, dtype=float)

        for li, lip in enumerate(lipid_resnames):
            h = means[:, bi, li]
            ax.bar(
                x, h, width=bar_w, bottom=bottom,
                color=lipid_colors[lip], edgecolor="none",
                label=lip if bi == 0 else None
            )

            for si in range(n_sims):
                occ = h[si]
                if occ < 0.06:
                    continue
                ax.text(
                    x[si], bottom[si] + 0.5 * h[si], f"{occ:.2f}",
                    ha="center", va="center", fontsize=12, color="black"
                )

            bottom += h

        ax.errorbar(
            x, bottom,
            yerr=sem_total_height[:, bi],
            fmt="none",
            ecolor="gray",
            elinewidth=1.6,
            capsize=3.5,
            capthick=1.6,
            zorder=1,
        )

        global_max = max(global_max, float(np.max(bottom)))

        for si in range(n_sims):
            ax.text(
                x[si], bottom[si] + 0.02 * max(1.0, global_max),
                res_labels[bi],
                ha="center", va="bottom", fontsize=15
            )

    ax.set_xticks(group_x)
    ax.set_xticklabels(sim_names, fontsize=18)
    ax.tick_params(axis="y", labelsize=16)

    ax.set_ylabel("H-bond Occupancy", fontsize=18)
    ax.set_ylim(0, max(1.0, global_max * 1.22))

    ax.legend(frameon=False, fontsize=11, loc="upper left", bbox_to_anchor=(0.02, 0.98))

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)

# -------------------- Heatmap --------------------
def compute_heatmap_mats(dfs, sim_names, resids, lipid_resnames):
    n_res = len(resids)
    n_sims = len(sim_names)
    mats = {}
    for lip in lipid_resnames:
        M = np.zeros((n_res, n_sims), dtype=float)
        for si, df in enumerate(dfs):
            for ri, rid in enumerate(resids):
                col = f"res{rid}_{lip}"
                M[ri, si] = float(df[col].mean()) if col in df.columns else 0.0
        mats[lip] = M

    Mtot = np.zeros((n_res, n_sims), dtype=float)
    for si, df in enumerate(dfs):
        for ri, rid in enumerate(resids):
            cols = [f"res{rid}_{lip}" for lip in lipid_resnames if f"res{rid}_{lip}" in df.columns]
            Mtot[ri, si] = float(df[cols].sum(axis=1).mean()) if cols else 0.0
    mats["TOTAL"] = Mtot
    return mats

def plot_hbond_heatmaps(dfs, sim_names, resids, lipid_resnames, outpath,
                        resid_label_map,
                        vmin=None, vmax=None):
    mats = compute_heatmap_mats(dfs, sim_names, resids, lipid_resnames)
    panels = list(lipid_resnames) + ["TOTAL"]
    n_pan = len(panels)

    all_vals = np.concatenate([mats[p].ravel() for p in panels])
    if vmin is None:
        vmin = float(np.nanmin(all_vals))
    if vmax is None:
        vmax = float(np.nanmax(all_vals))
        if vmax <= vmin:
            vmax = vmin + 1.0

    fig, axes = plt.subplots(1, n_pan, figsize=(2.75 * n_pan + 1.0, 3.6), sharey=True)
    if n_pan == 1:
        axes = [axes]

    cmap = "viridis"
    last_im = None

    yticklabels = [resid_label_map[int(r)] for r in resids]

    for ax, p in zip(axes, panels):
        M = mats[p]
        ny, nx = M.shape

        x = np.arange(nx + 1)
        y = np.arange(ny + 1)

        im = ax.pcolormesh(x, y, M, vmin=vmin, vmax=vmax, cmap=cmap, shading="flat")
        ax.set_xlim(0, nx)
        ax.set_ylim(ny, 0)

        last_im = im

        ax.set_title(p, fontsize=14)
        ax.set_xticks(np.arange(len(sim_names)))
        ax.set_xticklabels(sim_names, fontsize=12, rotation=30, ha="right")
        ax.set_yticks(np.arange(len(resids)))
        ax.set_yticklabels(yticklabels, fontsize=12)

        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(width=1.5, direction="out")

        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=9, color="black")

    cbar = fig.colorbar(
        last_im, ax=axes, orientation="vertical",
        fraction=0.028, pad=0.02, shrink=0.92, aspect=25
    )
    cbar.set_label("H-bond Occupancy", fontsize=12)
    cbar.ax.tick_params(labelsize=10, width=1.0, length=3)
    cbar.outline.set_linewidth(1.0)

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)

# -------------------- RDF compute --------------------
def compute_rdf_results(
    tops, trajs, sim_names, resids, lipid_resnames, lipid_p_names,
    r_range, nbins,
    include_total=True,
    start=None, stop=None, step=1,
    keep_last_fraction: float = PLOT_KEEP_LAST_FRACTION,
):
    """
    Computes RDFs for (sim, rid, lip) where rid is int residue or "TOTAL".
    Uses only the last `keep_last_fraction` of the trajectory by default.
    """
    rdf_dict = {}
    for top, traj, sim in zip(tops, trajs, sim_names):
        u = mda.Universe(top, traj)

        base_start = frame_start_for_last_fraction(u.trajectory.n_frames, keep_last_fraction)
        start_eff = base_start if start is None else max(int(start), base_start)
        stop_eff = u.trajectory.n_frames if stop is None else int(stop)
        step_eff = 1 if step is None else int(step)

        resid_items = [int(r) for r in resids]
        if include_total:
            resid_items = resid_items + ["TOTAL"]

        for rid in resid_items:
            if rid == "TOTAL":
                sel_res = donor_sel_for_resids(resids)
            else:
                sel_res = f"protein and resid {int(rid)} and ((resname ARG and name NH2) or (resname LYS and name NZ))"

            g1 = u.select_atoms(sel_res)
            if g1.n_atoms == 0:
                raise ValueError(f"[{sim}] RDF: selection matched 0 atoms for resid {rid}:\n{sel_res}")

            for lip in lipid_resnames:
                pnames = lipid_p_names.get(lip, None)
                if not pnames:
                    raise ValueError(f"No P atom names configured for lipid '{lip}'.")
                sel_p = f"resname {lip} and name P"
                g2 = u.select_atoms(sel_p)
                if g2.n_atoms == 0:
                    print(f"[WARNING] [{sim}] RDF: 0 atoms for lipid P selection: {sel_p}")
                    continue

                rdf = InterRDF(g1, g2, range=r_range, nbins=nbins, verbose=True)
                rdf.run(start=start_eff, stop=stop_eff, step=step_eff)

                rdf_dict[(sim, rid, lip)] = {"bins": np.array(rdf.bins), "rdf": np.array(rdf.rdf)}
    return rdf_dict

# -------------------- RDF plot --------------------
def plot_rdfs_all_residues_vertical(
    rdf_dict, sim_names, resids, lipid_resnames, outpath,
    lipid_colors, resid_styles, resid_label_map
):
    n_sims = len(sim_names)
    fig, axes = plt.subplots(n_sims, 1, figsize=(6.0, 2.7 * n_sims), sharex=True, sharey=True)
    if n_sims == 1:
        axes = [axes]

    plot_resids = [int(r) for r in resids] + ["TOTAL"]

    for ax, sim in zip(axes, sim_names):
        for lip in lipid_resnames:
            for rid in plot_resids:
                key = (sim, rid, lip)
                if key not in rdf_dict:
                    continue
                st = resid_styles[rid]
                ax.plot(
                    rdf_dict[key]["bins"], rdf_dict[key]["rdf"],
                    color=lipid_colors[lip],
                    linestyle=st["linestyle"],
                    linewidth=st["lw"],
                    marker=st["marker"],
                    markevery=st["markevery"],
                    markersize=st["ms"] if st["ms"] is not None else None,
                )

        ax.set_title(sim, fontsize=14, loc="left")
        style_axes(ax)
        ax.tick_params(labelsize=12)

    axes[-1].set_xlabel("r (Å)", fontsize=12)
    axes[0].set_ylabel("g(r)", fontsize=12)

    lipid_handles = [plt.Line2D([0], [0], color=lipid_colors[lip], lw=3, label=lip) for lip in lipid_resnames]
    resid_handles = []
    for r in resids:
        st = resid_styles[int(r)]
        resid_handles.append(
            plt.Line2D([0], [0], color="black", lw=st["lw"], linestyle=st["linestyle"],
                       marker=st["marker"], markersize=st["ms"] if st["ms"] else 0,
                       label=resid_label_map[int(r)])
        )
    stT = resid_styles["TOTAL"]
    resid_handles.append(
        plt.Line2D([0], [0], color="black", lw=stT["lw"], linestyle=stT["linestyle"], label="TOTAL")
    )

    leg1 = axes[0].legend(handles=lipid_handles, frameon=False, fontsize=10, loc="upper right", title="Lipid")
    axes[0].add_artist(leg1)
    axes[0].legend(handles=resid_handles, frameon=False, fontsize=10, loc="upper left", title="Residue")

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")

def aggregate_block_rdfs(rdf_blocks, sim, n_blocks, resids, lipid_resnames):
    """
    rdf_blocks keys: (sim, bi, rid, lip) -> {"bins":..., "rdf":...}
    Returns: dict[(rid, lip)] = {"bins": bins, "mean": mean_rdf, "sem": sem_rdf, "n": n_used}
    """
    out = {}
    for rid in resids:
        for lip in lipid_resnames:
            series = []
            bins_ref = None
            for bi in range(1, n_blocks + 1):
                k = (sim, bi, int(rid), lip)
                if k not in rdf_blocks:
                    continue
                bins = np.asarray(rdf_blocks[k]["bins"], dtype=float)
                rdf  = np.asarray(rdf_blocks[k]["rdf"], dtype=float)

                if bins_ref is None:
                    bins_ref = bins
                else:
                    if bins.shape != bins_ref.shape or np.max(np.abs(bins - bins_ref)) > 1e-8:
                        raise ValueError(
                            f"[aggregate_block_rdfs] Bin mismatch for {sim} rid={rid} lip={lip} block={bi}. "
                            f"Make sure all blocks use same r_range and nbins."
                        )
                series.append(rdf)

            if not series:
                continue

            arr = np.vstack(series)
            mean = arr.mean(axis=0)
            if arr.shape[0] > 1:
                sem = arr.std(axis=0, ddof=1) / np.sqrt(arr.shape[0])
            else:
                sem = np.zeros_like(mean)

            out[(int(rid), lip)] = {"bins": bins_ref, "mean": mean, "sem": sem, "n": int(arr.shape[0])}
    return out

def plot_blockavg_rdfs_overlay_two_sims_per_residue(
    aggA, aggB, sim_a, sim_b,
    resids, lipid_resnames,
    outpath, lipid_colors, resid_label_map,
    ls_a="-", ls_b=":", band_alpha=0.18
):
    n_res = len(resids)
    fig, axes = plt.subplots(n_res, 1, figsize=(3, 4), sharex=True, sharey=True)
    if n_res == 1:
        axes = [axes]

    for ax, rid in zip(axes, resids):
        for lip in lipid_resnames:
            k = (int(rid), lip)

            if k in aggA:
                bins = aggA[k]["bins"]
                mean = aggA[k]["mean"]
                sem  = aggA[k]["sem"]
                c = lipid_colors[lip]
                ax.plot(bins, mean, color=c, lw=1.5, ls=ls_a)
                ax.fill_between(bins, mean - sem, mean + sem, color=c, alpha=band_alpha, linewidth=0)

            if k in aggB:
                bins = aggB[k]["bins"]
                mean = aggB[k]["mean"]
                sem  = aggB[k]["sem"]
                c = lipid_colors[lip]
                ax.plot(bins, mean, color=c, lw=1.5, ls=ls_b)
                ax.fill_between(bins, mean - sem, mean + sem, color=c, alpha=band_alpha, linewidth=0)

        ax.set_title(resid_label_map[int(rid)], fontsize=14, loc="right", y=.7)
        style_axes(ax)
        ax.tick_params(labelsize=11)

    fig.supxlabel("r (Å)", fontsize=14, y=-.02)
    fig.supylabel("g(r)", fontsize=14)

    lipid_handles = [plt.Line2D([0], [0], color=lipid_colors[lip], lw=3, label=lip) for lip in lipid_resnames]
    sim_handles = [
        plt.Line2D([0], [0], color="black", lw=2.5, ls=ls_a, label=sim_a),
        plt.Line2D([0], [0], color="black", lw=2.5, ls=ls_b, label=sim_b),
    ]

    leg1 = axes[0].legend(handles=lipid_handles, frameon=False, fontsize=6, loc="upper right")
    axes[0].add_artist(leg1)
    axes[0].legend(handles=sim_handles, frameon=False, fontsize=6, loc="upper left")

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")

def plot_blockavg_rdfs_single_sim_per_residue(
    agg, sim, resids, lipid_resnames,
    outpath, lipid_colors, resid_label_map,
    ls="-", band_alpha=0.18
):
    n_res = len(resids)
    fig, axes = plt.subplots(n_res, 1, figsize=(6.6, 2.5 * n_res), sharex=True, sharey=True)
    if n_res == 1:
        axes = [axes]

    for ax, rid in zip(axes, resids):
        for lip in lipid_resnames:
            k = (int(rid), lip)
            if k not in agg:
                continue
            bins = agg[k]["bins"]
            mean = agg[k]["mean"]
            sem  = agg[k]["sem"]
            c = lipid_colors[lip]
            ax.plot(bins, mean, color=c, lw=2.4, ls=ls)
            ax.fill_between(bins, mean - sem, mean + sem, color=c, alpha=band_alpha, linewidth=0)

        ax.set_title(resid_label_map[int(rid)], fontsize=13, loc="left")
        style_axes(ax)
        ax.tick_params(labelsize=11)

    axes[-1].set_xlabel("r (Å)", fontsize=12)
    axes[0].set_ylabel("g(r)", fontsize=12)

    lipid_handles = [plt.Line2D([0], [0], color=lipid_colors[lip], lw=3, label=lip) for lip in lipid_resnames]
    axes[0].legend(handles=lipid_handles, frameon=False, fontsize=10, loc="upper right", title="Lipid")

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")

# -------------------- Combo figure: Heatmap + RDF --------------------
def plot_combo_heatmap_plus_rdf(
    dfs, sim_names, resids, lipid_resnames,
    rdf_dict,
    outpath,
    lipid_colors, resid_styles, resid_label_map,
    heat_vmin=None, heat_vmax=None
):
    mats = compute_heatmap_mats(dfs, sim_names, resids, lipid_resnames)
    heat_panels = list(lipid_resnames) + ["TOTAL"]
    n_pan = len(heat_panels)

    all_vals = np.concatenate([mats[p].ravel() for p in heat_panels])
    if heat_vmin is None:
        heat_vmin = float(np.nanmin(all_vals))
    if heat_vmax is None:
        heat_vmax = float(np.nanmax(all_vals))
        if heat_vmax <= heat_vmin:
            heat_vmax = heat_vmin + 1.0

    fig = plt.figure(figsize=(7.6, 3.3 + 2.4 * len(sim_names)))
    gs = GridSpec(2, 1, height_ratios=[1.05, 2.1], hspace=0.28, figure=fig)

    gs_heat = GridSpecFromSubplotSpec(1, n_pan, subplot_spec=gs[0], wspace=0.12)
    heat_axes = [fig.add_subplot(gs_heat[0, j]) for j in range(n_pan)]

    cmap = "viridis"
    yticklabels = [resid_label_map[int(r)] for r in resids]
    last_im = None

    for ax, p in zip(heat_axes, heat_panels):
        M = mats[p]
        im = ax.imshow(M, aspect="auto", origin="upper", vmin=heat_vmin, vmax=heat_vmax, cmap=cmap)
        last_im = im

        ax.set_title(p, fontsize=13)
        ax.set_xticks(np.arange(len(sim_names)))
        ax.set_xticklabels(sim_names, fontsize=11, rotation=30, ha="right")
        ax.set_yticks(np.arange(len(resids)))
        ax.set_yticklabels(yticklabels, fontsize=11)

        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(width=1.2, direction="out")

        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=8.5, color="black")

    cbar = fig.colorbar(last_im, ax=heat_axes, orientation="vertical",
                        fraction=0.025, pad=0.02, shrink=0.90, aspect=25)
    cbar.set_label("H-bond Occupancy", fontsize=11)
    cbar.ax.tick_params(labelsize=9, width=1.0, length=3)
    cbar.outline.set_linewidth(1.0)

    gs_rdf = GridSpecFromSubplotSpec(len(sim_names), 1, subplot_spec=gs[1], hspace=0.18)
    rdf_axes = [fig.add_subplot(gs_rdf[i, 0]) for i in range(len(sim_names))]

    plot_resids = [int(r) for r in resids] + ["TOTAL"]

    for ax, sim in zip(rdf_axes, sim_names):
        for lip in lipid_resnames:
            for rid in plot_resids:
                key = (sim, rid, lip)
                if key not in rdf_dict:
                    continue
                st = resid_styles[rid]
                ax.plot(
                    rdf_dict[key]["bins"], rdf_dict[key]["rdf"],
                    color=lipid_colors[lip],
                    linestyle=st["linestyle"],
                    linewidth=st["lw"],
                    marker=st["marker"],
                    markevery=st["markevery"],
                    markersize=st["ms"] if st["ms"] is not None else None,
                )
        ax.set_title(sim, fontsize=13, loc="left")
        style_axes(ax)
        ax.tick_params(labelsize=11)

    rdf_axes[-1].set_xlabel("r (Å)", fontsize=12)
    rdf_axes[0].set_ylabel("g(r)", fontsize=12)

    lipid_handles = [plt.Line2D([0], [0], color=lipid_colors[lip], lw=3, label=lip) for lip in lipid_resnames]
    resid_handles = []
    for r in resids:
        st = resid_styles[int(r)]
        resid_handles.append(
            plt.Line2D([0], [0], color="black", lw=st["lw"], linestyle=st["linestyle"],
                       marker=st["marker"], markersize=st["ms"] if st["ms"] else 0,
                       label=resid_label_map[int(r)])
        )
    stT = resid_styles["TOTAL"]
    resid_handles.append(plt.Line2D([0], [0], color="black", lw=stT["lw"], linestyle=stT["linestyle"], label="TOTAL"))

    leg1 = rdf_axes[0].legend(handles=lipid_handles, frameon=False, fontsize=9.5, loc="upper right", title="Lipid")
    rdf_axes[0].add_artist(leg1)
    rdf_axes[0].legend(handles=resid_handles, frameon=False, fontsize=9.5, loc="upper left", title="Residue")

    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")

def plot_combo_heatmap_left_plus_rdf_overlay(
    dfs, sim_names, resids, lipid_resnames,
    rdf_blocks, n_blocks,
    outpath,
    lipid_colors, resid_label_map,
    heat_vmin=None, heat_vmax=None,
    rdf_err="sem",
    band_alpha=0.18,
    ls_gdp="-", ls_gtp=":",
):
    if len(sim_names) != 2:
        raise ValueError("This combo layout requires exactly two sims (GDP/GTP).")
    simGDP, simGTP = sim_names[0], sim_names[1]

    mats = compute_heatmap_mats(dfs, sim_names, resids, lipid_resnames)
    panels = list(lipid_resnames) + ["TOTAL"]
    n_pan = len(panels)

    all_vals = np.concatenate([mats[p].ravel() for p in panels])
    if heat_vmin is None:
        heat_vmin = float(np.nanmin(all_vals))
    if heat_vmax is None:
        heat_vmax = float(np.nanmax(all_vals))
        if heat_vmax <= heat_vmin:
            heat_vmax = heat_vmin + 1.0

    def _aggregate(sim):
        agg = aggregate_block_rdfs(rdf_blocks, sim, n_blocks, resids, lipid_resnames)
        if rdf_err == "std":
            for k, v in agg.items():
                n = max(1, int(v.get("n", 1)))
                v["std"] = v["sem"] * np.sqrt(n)
        return agg

    aggGDP = _aggregate(simGDP)
    aggGTP = _aggregate(simGTP)

    fig = plt.figure(figsize=(11.6, 6.2))
    outer = GridSpec(1, 2, width_ratios=[1.05, 2.55], wspace=0.28, figure=fig)

    gs_left = GridSpecFromSubplotSpec(
        n_pan + 1, 1, subplot_spec=outer[0], height_ratios=[0.18] + [1.0] * n_pan, hspace=0.22
    )
    cax = fig.add_subplot(gs_left[0, 0])
    cax.axis("off")

    heat_axes = [fig.add_subplot(gs_left[i + 1, 0]) for i in range(n_pan)]

    cmap = "viridis"
    last_im = None
    yticklabels = [resid_label_map[int(r)] for r in resids]

    for i, (ax, p) in enumerate(zip(heat_axes, panels)):
        M = mats[p]
        ny, nx = M.shape

        x = np.arange(nx + 1)
        y = np.arange(ny + 1)

        im = ax.pcolormesh(
            x, y, M,
            vmin=heat_vmin, vmax=heat_vmax,
            cmap=cmap, shading="flat"
        )
        last_im = im

        ax.set_xticks(np.arange(nx) + 0.5)
        ax.set_yticks(np.arange(ny) + 0.5)

        ax.set_yticklabels(yticklabels, fontsize=10.5)

        if i == n_pan - 1:
            ax.set_xticklabels(sim_names, fontsize=10.5, rotation=30, ha="right")
        else:
            ax.set_xticklabels([])

        ax.set_xlim(0, nx)
        ax.set_ylim(ny, 0)

        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(width=1.0, length=3, direction="out")

        for r in range(ny):
            for c in range(nx):
                ax.text(
                    c + 0.5, r + 0.5, f"{M[r, c]:.2f}",
                    ha="center", va="center",
                    fontsize=8.2, color="black"
                )

        ax.text(
            1.02, 0.5, p,
            transform=ax.transAxes,
            rotation=-90,
            ha="left", va="center",
            fontsize=11
        )

    cax.axis("on")
    cbar = fig.colorbar(last_im, cax=cax, orientation="horizontal")
    cbar.set_label("H-bond Occupancy", fontsize=11, labelpad=6)
    cbar.ax.xaxis.set_label_position("top")
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.tick_params(labelsize=9.5, width=1.0, length=3)
    cbar.outline.set_linewidth(1.0)

    gs_rdf = GridSpecFromSubplotSpec(len(resids), 1, subplot_spec=outer[1], hspace=0.26)

    rdf_axes = []
    ref_ax = None
    for i in range(len(resids)):
        if ref_ax is None:
            ax = fig.add_subplot(gs_rdf[i, 0])
            ref_ax = ax
        else:
            ax = fig.add_subplot(gs_rdf[i, 0], sharey=ref_ax)
        rdf_axes.append(ax)

    for ax, rid in zip(rdf_axes, resids):
        rid = int(rid)

        for lip in lipid_resnames:
            k = (rid, lip)
            c = lipid_colors[lip]

            if k in aggGDP:
                bins = aggGDP[k]["bins"]
                mean = aggGDP[k]["mean"]
                errv = aggGDP[k]["std"] if (rdf_err == "std") else aggGDP[k]["sem"]
                ax.plot(bins, mean, color=c, lw=2.25, ls=ls_gdp)
                ax.fill_between(bins, mean - errv, mean + errv, color=c, alpha=band_alpha, linewidth=0)

            if k in aggGTP:
                bins = aggGTP[k]["bins"]
                mean = aggGTP[k]["mean"]
                errv = aggGTP[k]["std"] if (rdf_err == "std") else aggGTP[k]["sem"]
                ax.plot(bins, mean, color=c, lw=2.25, ls=ls_gtp)
                ax.fill_between(bins, mean - errv, mean + errv, color=c, alpha=band_alpha, linewidth=0)

        style_axes(ax)
        ax.tick_params(labelsize=11)

        ax.text(
            0.02, 0.88, resid_label_map[rid],
            transform=ax.transAxes,
            ha="left", va="top",
            fontsize=12,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="black", linewidth=1.2)
        )

    rdf_axes[0].set_title("RDFs (block mean ± error) — GDP solid, GTP dotted", fontsize=12.5, loc="left")
    rdf_axes[-1].set_xlabel("r (Å)", fontsize=12)

    mid = len(rdf_axes) // 2
    rdf_axes[mid].set_ylabel("g(r)", fontsize=12)

    lipid_handles = [plt.Line2D([0], [0], color=lipid_colors[lip], lw=3, label=lip) for lip in lipid_resnames]
    sim_handles = [
        plt.Line2D([0], [0], color="black", lw=2.4, ls=ls_gdp, label=simGDP),
        plt.Line2D([0], [0], color="black", lw=2.4, ls=ls_gtp, label=simGTP),
    ]
    leg1 = rdf_axes[0].legend(handles=lipid_handles, frameon=False, fontsize=10, loc="upper right", title="Lipid")
    rdf_axes[0].add_artist(leg1)
    rdf_axes[0].legend(handles=sim_handles, frameon=False, fontsize=10, loc="upper left", title="Simulation")

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")

# -------------------- Bundle plotting --------------------
def plot_bundle(
    dfs_ts_sub,
    dfs_stats_sub,
    names_sub, lipid_list_sub, tag,
    outdir,
    args,
    resid_label_map,
    rdf_dict_main=None,
    make_combo=False,
):
    lipid_colors_sub = make_lipid_color_map(lipid_list_sub)
    resid_styles = make_resid_styles(args.resids)

    ts_out = os.path.join(outdir, f"hbonds_timeseries_{tag}.svg")
    plot_hbond_timeseries(
        dfs_ts_sub, names_sub, tuple(lipid_list_sub), ts_out,
        lipid_colors=lipid_colors_sub, ma_window=args.ma_window
    )
    print(f"Wrote {ts_out}")

    heat_out = os.path.join(outdir, f"hbonds_heatmaps_{tag}.svg")
    plot_hbond_heatmaps(
        dfs_stats_sub, names_sub, tuple(args.resids), tuple(lipid_list_sub),
        heat_out, resid_label_map=resid_label_map,
        vmin=args.heatmap_vmin, vmax=args.heatmap_vmax
    )
    print(f"Wrote {heat_out}")

    if rdf_dict_main is not None:
        rdf_out = os.path.join(outdir, f"rdfs_all_residues_{tag}.svg")
        plot_rdfs_all_residues_vertical(
            rdf_dict_main, names_sub, tuple(args.resids), tuple(lipid_list_sub),
            rdf_out, lipid_colors=lipid_colors_sub,
            resid_styles=resid_styles, resid_label_map=resid_label_map
        )

        if make_combo:
            combo_out = os.path.join(outdir, f"combo_heatmap_rdf_{tag}.svg")
            plot_combo_heatmap_plus_rdf(
                dfs_stats_sub, names_sub, tuple(args.resids), tuple(lipid_list_sub),
                rdf_dict_main,
                combo_out,
                lipid_colors=lipid_colors_sub,
                resid_styles=resid_styles,
                resid_label_map=resid_label_map,
                heat_vmin=args.heatmap_vmin,
                heat_vmax=args.heatmap_vmax,
            )

def plot_resim_timeseries_plus_grouped_stacked(
    df_resim: pd.DataFrame,
    sim_name: str,
    resids,
    lipid_resnames,
    outpath: str,
    lipid_colors: dict,
    resid_label_map: dict,
    ma_window: int = 500,
    dt_ns: float | None = None,
    vlines_us=(4.8, 5.5),
    figsize=(3.9, 2.65),
    label_threshold_pct: float = 6.0,
):
    """
    2-row compact figure:
      Row 1: H-bond time evolution (sum_<lip>), raw + moving average
      Row 2: Grouped stacked bars for this sim: bars = residues + TOTAL, stacks = lipids.
    """
    if df_resim is None or len(df_resim) == 0:
        raise ValueError("df_resim is empty.")

    total_us = 20.0
    n = len(df_resim)
    t = np.linspace(0.0, total_us, num=n, endpoint=True)

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 1, height_ratios=[1.25, 1.0], hspace=0.28, figure=fig)
    ax_ts = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[1, 0])

    for lip in lipid_resnames:
        col = f"sum_{lip}"
        if col not in df_resim.columns:
            continue
        y = df_resim[col].to_numpy(dtype=float)
        line_plot_raw_ma(ax_ts, t, y, label=lip, color=lipid_colors[lip], ma_window=ma_window)

    for xvl in vlines_us:
        ax_ts.axvline(
            float(xvl),
            color="#1F77B4",
            linestyle=(0, (3.5, 3.5)),
            linewidth=1.6,
            alpha=1.0,
            zorder=5,
        )

    style_axes(ax_ts)
    ax_ts.set_ylabel("H-bond Occupancy", fontsize=10)
    ax_ts.tick_params(labelsize=9)
    ax_ts.set_title(sim_name, fontsize=11, loc="left", pad=2)
    ax_ts.legend(frameon=False, fontsize=8, loc="upper right", handlelength=2.0)
    ax_ts.set_xlim(0.0, total_us)

    style_axes(ax_bar)

    res_keys = list(resids) + ["TOTAL"]
    res_labels = [resid_label_map[int(r)] for r in resids] + ["TOTAL"]
    n_bars = len(res_keys)
    n_lip = len(lipid_resnames)

    means = np.zeros((n_bars, n_lip), dtype=float)
    for bi, key in enumerate(res_keys):
        for li, lip in enumerate(lipid_resnames):
            col = (f"sum_{lip}" if key == "TOTAL" else f"res{int(key)}_{lip}")
            means[bi, li] = float(df_resim[col].mean()) if col in df_resim.columns else 0.0

    total_all = float(np.sum(means))

    group_x = 0.0
    bar_w = 0.12
    offsets = (np.arange(n_bars) - (n_bars - 1) / 2.0) * (bar_w * 1.25)

    global_max = 0.0
    for bi, key in enumerate(res_keys):
        x = group_x + offsets[bi]
        bottom = 0.0

        for li, lip in enumerate(lipid_resnames):
            h = means[bi, li]
            ax_bar.bar(
                x, h, width=bar_w, bottom=bottom,
                color=lipid_colors[lip], edgecolor="none",
                label=None
            )

            if total_all > 1e-12 and h > 1e-12:
                pct = 100.0 * h / total_all
                if pct >= float(label_threshold_pct):
                    ax_bar.text(
                        x, bottom + 0.5 * h, f"{pct:.0f}%",
                        ha="center", va="center",
                        fontsize=8.5, color="black"
                    )

            bottom += h

        global_max = max(global_max, bottom)

        ax_bar.text(
            x, bottom + 0.03 * max(1.0, global_max),
            res_labels[bi],
            ha="center", va="bottom", fontsize=9
        )

    ax_bar.tick_params(axis="y", labelsize=9)
    fig.supylabel("H-bond Occupancy", fontsize=12)

    ax_bar.set_xticks([])
    ax_bar.set_xlim(-0.6, 0.6)
    ax_bar.set_ylim(0, max(1.0, global_max * 1.22))

    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outpath}")

def slice_last_half_for_plot(df: pd.DataFrame, zero_time: bool = True) -> pd.DataFrame:
    """
    Keep only the last half of the dataframe for plotting timeseries as before.
    """
    half_idx = len(df) // 2
    out = df.iloc[half_idx:].copy().reset_index(drop=True)

    if zero_time and "time_us" in out.columns and len(out) > 0:
        out["time_us"] = out["time_us"] - float(out["time_us"].iloc[0])

    return out

# -------------------- Main --------------------
def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--tops", nargs=3, required=False, help="Topology files for 3 sims (compute mode).")
    ap.add_argument("--trajs", nargs=3, required=False, help="Trajectory files for 3 sims (compute mode).")
    ap.add_argument("--names", nargs=3, required=True, help="Names for 3 sims (e.g., GDP GTP Resim).")

    ap.add_argument("--outdir", default="hbonds_out", help="Output directory for plots/CSVs.")
    ap.add_argument("--dt-ns", type=float, default=1.0, help="Time per sampled frame (ns).")
    ap.add_argument("--resids", nargs="+", type=int, default=[5,6,7], help="Residues to analyze.")
    ap.add_argument("--lipids", nargs="+", default=["POPC","POPE","POPS"], help="Base lipid resnames to analyze.")
    ap.add_argument("--da", type=float, default=3.1, help="Donor-acceptor distance cutoff (Å).")
    ap.add_argument("--angle", type=float, default=130.0, help="Minimum D-H-A angle cutoff (deg).")
    ap.add_argument("--ma-window", type=int, default=500, help="Moving average window (frames).")
    ap.add_argument("--rdf-rmax", type=float, default=20.0, help="RDF max distance (Å).")
    ap.add_argument("--rdf-nbins", type=int, default=200, help="RDF nbins.")
    ap.add_argument("--start", type=int, default=None, help="Start frame index for HBA.")
    ap.add_argument("--stop", type=int, default=None, help="Stop frame index for HBA.")
    ap.add_argument("--step", type=int, default=None, help="Step for HBA/RDF (frames).")
    ap.add_argument("--summary_label_threshold_pct", type=float, default=6.0,
                    help="Skip percent labels smaller than this (reduce clutter).")

    ap.add_argument("--heatmap-vmin", type=float, default=None, help="Heatmap vmin (optional).")
    ap.add_argument("--heatmap-vmax", type=float, default=None, help="Heatmap vmax (optional).")

    ap.add_argument("--resim-name", default="Resim", help="Sim name that should include extra lipid (default Resim).")
    ap.add_argument("--resim-extra-lipid", default="PLA18", help="Extra lipid resname present only in Resim.")

    ap.add_argument("--residue-labels", default=None,
                    help='Optional override like "5:R,6:K,7:K" to label residues as R5/K6/etc.')

    ap.add_argument("--save-cache", default=None,
                    help="Directory to save pickles (hbonds_<sim>.pkl and rdf_results.pkl + meta.json).")
    ap.add_argument("--load-dir", default=None,
                    help="Directory containing cached pickles to load (plot-only mode).")
    ap.add_argument("--plot-only", action="store_true",
                    help="Skip MDAnalysis H-bond compute; load hbonds from --load-dir and plot.")

    ap.add_argument("--combo-heatmap-rdf", action="store_true",
                    help="If set, also output a single svg combining heatmap (top) + RDF (bottom) stacked vertically.")

    ap.add_argument("--time-blocks", type=int, default=1,
                    help="Split each trajectory into this many time blocks for block-wise hbonds and RDFs (default 1 = off).")
    ap.add_argument("--compute-block-rdfs", action="store_true",
                    help="Compute block-wise RDFs for all sims. Requires --tops/--trajs even in plot-only mode.")
    ap.add_argument("--load-block-rdfs", action="store_true",
                    help="Load block-wise RDFs from cache (rdf_blocks.pkl) and plot them.")
    ap.add_argument("--rdf-blocks-pkl", default="rdf_blocks.pkl",
                    help="Filename to use inside cache dir for block RDF dict (default rdf_blocks.pkl).")
    ap.add_argument("--hbonds-blocks-pkl", default="hbonds_blocks.pkl",
                    help="Filename to use inside cache dir for block hbonds dict (default hbonds_blocks.pkl).")
    ap.add_argument("--bar-sem-blocks", type=int, default=1,
                help="Number of contiguous time blocks for stacked-bar block mean ± SEM. Use 1 to disable.")
    ap.add_argument("--resim-stride", type=int, default=10,
                    help="Stride factor already applied to the Resim trajectory on disk (default 10). "
                         "We divide start/stop/step by this factor for Resim so windows match original frames.")
    ap.add_argument("--resim-equil-us", type=float, default=5.0,
                    help="Equilibration time (µs) for Resim. Only used for reference line in compact Resim plot.")

    args = ap.parse_args()

    names = list(args.names)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    base_lipids = list(args.lipids)
    resim_name = args.resim_name
    extra_lipid = args.resim_extra_lipid

    plot_lipids_all = list(base_lipids)
    if extra_lipid not in plot_lipids_all:
        plot_lipids_all.append(extra_lipid)

    lipid_acceptor_names = dict(DEFAULT_LIPID_ACCEPTOR_NAMES)
    lipid_p_names = dict(DEFAULT_LIPID_P_NAMES)
    if extra_lipid not in lipid_acceptor_names:
        lipid_acceptor_names[extra_lipid] = lipid_acceptor_names["POPE"]
    if extra_lipid not in lipid_p_names:
        lipid_p_names[extra_lipid] = lipid_p_names["POPE"]

    resid_label_map = parse_residue_labels_arg(args.residue_labels, args.resids)

    dfs = None
    tops = args.tops
    trajs = args.trajs

    if args.plot_only:
        if args.load_dir is None:
            raise ValueError("--plot-only requires --load-dir <cache_dir>")
        dfs = load_cached_hbonds(args.load_dir, names)
        if resid_label_map is None:
            if tops is not None and trajs is not None:
                resid_label_map = get_residue_label_map_from_universe(tops[0], trajs[0], args.resids)
            else:
                resid_label_map = {int(r): f"res{int(r)}" for r in args.resids}
    else:
        if tops is None or trajs is None:
            raise ValueError("Compute mode requires --tops and --trajs (or use --plot-only with --load-dir).")

        if resid_label_map is None:
            resid_label_map = get_residue_label_map_from_universe(tops[0], trajs[0], args.resids)

        dfs = []
        per_sim_lipids = {}

        for top, traj, sim in zip(tops, trajs, names):
            lip_list = list(base_lipids)
            if sim == resim_name and extra_lipid not in lip_list:
                lip_list.append(extra_lipid)
            per_sim_lipids[sim] = lip_list

            stride_factor = (args.resim_stride if sim == resim_name else 1)

            df_full = run_hbonds_and_aggregate(
                top=top,
                traj=traj,
                sim_name=sim,
                resids=tuple(args.resids),
                lipid_resnames=tuple(lip_list),
                dt_ns=args.dt_ns,
                d_a_cutoff=args.da,
                d_h_a_angle_cutoff=args.angle,
                start=args.start,
                stop=args.stop,
                step=args.step,
                outdir=outdir,
                write_csv=True,
                resim_stride_factor=stride_factor,
            )

            dfs.append(df_full)

        if args.save_cache is not None:
            cache_dir = args.save_cache
            for sim, df in zip(names, dfs):
                save_hbonds_pickle(cache_dir, sim, df)

    rdf_dict = None
    if not args.plot_only:
        rdf_dict = {}
        for top, traj, sim in zip(tops, trajs, names):
            lip_list = list(base_lipids)
            if sim == resim_name and extra_lipid not in lip_list:
                lip_list.append(extra_lipid)

            rdf_part = compute_rdf_results(
                tops=[top], trajs=[traj], sim_names=[sim],
                resids=tuple(args.resids),
                lipid_resnames=tuple(lip_list),
                lipid_p_names=lipid_p_names,
                r_range=(0.0, args.rdf_rmax),
                nbins=args.rdf_nbins,
                include_total=True,
                start=None, stop=None,
                step=(args.step if args.step is not None else 1),
                keep_last_fraction=PLOT_KEEP_LAST_FRACTION,
            )
            rdf_dict.update(rdf_part)

        if args.save_cache is not None:
            save_rdf_pickle(args.save_cache, rdf_dict, fname="rdf_results.pkl")
    else:
        try:
            rdf_dict = load_cached_rdfs(args.load_dir, fname="rdf_results.pkl")
        except Exception:
            rdf_dict = None
            print("[plot-only] No rdf_results.pkl found; RDF plots (and combo) will be skipped unless you compute RDFs.")
        # -------------------- Single-residue RDF figures --------------------
    # Main text: residue 5 only
    # SI: same style for residues 6 and 7
    if rdf_dict is not None:
        rdf_sim_order = names[:]   # keep current sim ordering; should be GDP, GTP, Resim
        rdf_lipids_all = tuple(plot_lipids_all)
        rdf_colors_all = make_lipid_color_map(rdf_lipids_all)

        # Main text figure: residue 5
        out_rdf_res5 = os.path.join(outdir, "rdf_res5_all_sims_maintext.svg")
        plot_rdf_single_residue_all_sims(
            rdf_dict=rdf_dict,
            sim_names=rdf_sim_order,
            resid=5,
            lipid_resnames=rdf_lipids_all,
            outpath=out_rdf_res5,
            lipid_colors=rdf_colors_all,
            resid_label_map=resid_label_map,
            sim_linestyles={
                rdf_sim_order[0]: "-",
                rdf_sim_order[1]: "--",
                rdf_sim_order[2]: ":",
            } if len(rdf_sim_order) >= 3 else None,
        )

        # SI figures: residues 6 and 7
        for rid in (6, 7):
            out_rdf_si = os.path.join(outdir, f"rdf_res{rid}_all_sims_SI.svg")
            plot_rdf_single_residue_all_sims(
                rdf_dict=rdf_dict,
                sim_names=rdf_sim_order,
                resid=rid,
                lipid_resnames=rdf_lipids_all,
                outpath=out_rdf_si,
                lipid_colors=rdf_colors_all,
                resid_label_map=resid_label_map,
                sim_linestyles={
                    rdf_sim_order[0]: "-",
                    rdf_sim_order[1]: "--",
                    rdf_sim_order[2]: ":",
                } if len(rdf_sim_order) >= 3 else None,
            )

    # Timeseries data left as before
    dfs_plot = [slice_last_half_for_plot(df, zero_time=True) for df in dfs]
        # Combined stacked timeseries: GDP, GTP, Resim
    if len(names) == 3:
        ts_all_out = os.path.join(outdir, "hbonds_timeseries_all_three_stacked.pdf")
        plot_hbond_timeseries_stacked_three(
            dfs=dfs,
            sim_names=names,
            lipid_resnames=tuple(plot_lipids_all),
            outpath=ts_all_out,
            lipid_colors=make_lipid_color_map(plot_lipids_all),
            ma_window=args.ma_window,
        )

    # Summary / heatmap / occupancy use last 3/4
    dfs_stats = [slice_last_fraction_df(df, fraction=PLOT_KEEP_LAST_FRACTION, zero_time=False) for df in dfs]

    idx_resim = names.index(resim_name) if resim_name in names else None

    def subset(df_list, names_keep):
        return [df for df, nm in zip(df_list, names) if nm in names_keep], names_keep

    if idx_resim is not None:
        names_A = [n for n in names if n != resim_name]
        names_B = [resim_name]

        dfs_A_ts, names_A = subset(dfs_plot, names_A)
        dfs_B_ts, names_B = subset(dfs_plot, names_B)

        dfs_A_stats, _ = subset(dfs_stats, names_A)
        dfs_B_stats, _ = subset(dfs_stats, names_B)

        lip_A = list(base_lipids)
        lip_B = list(base_lipids) + [extra_lipid]

        plot_bundle(
            dfs_A_ts, dfs_A_stats, names_A, lip_A, tag="GDPGTP",
            outdir=outdir, args=args, resid_label_map=resid_label_map,
            rdf_dict_main=rdf_dict, make_combo=args.combo_heatmap_rdf
        )

        plot_bundle(
            dfs_B_ts, dfs_B_stats, names_B, lip_B, tag="Resim",
            outdir=outdir, args=args, resid_label_map=resid_label_map,
            rdf_dict_main=rdf_dict, make_combo=False
        )

        # Combined 3-group occupancy summary: GDP, GTP, Resim
        summary_all_out = os.path.join(outdir, "hbonds_grouped_stacked_summary_All3_last3q.svg")
        plot_grouped_stacked_summary(
            dfs_stats,
            names,
            tuple(args.resids),
            tuple(plot_lipids_all),
            summary_all_out,
            lipid_colors=make_lipid_color_map(plot_lipids_all),
            resid_label_map=resid_label_map,
            label_threshold_pct=args.summary_label_threshold_pct,
            sem_blocks=args.bar_sem_blocks,
        )
        print(f"Wrote {summary_all_out}")

        lipid_colors_resim = make_lipid_color_map(lip_B)
        df_resim_stats = dfs_stats[idx_resim]

        out_compact = os.path.join(outdir, "hbonds_timeseries_plus_summary_Resim_last3q.svg")
        plot_resim_timeseries_plus_grouped_stacked(
            df_resim=df_resim_stats,
            sim_name=f"{resim_name} (last 3/4)",
            resids=tuple(args.resids),
            lipid_resnames=tuple(lip_B),
            outpath=out_compact,
            lipid_colors=lipid_colors_resim,
            resid_label_map=resid_label_map,
            ma_window=args.ma_window,
            dt_ns=None,
            vlines_us=(args.resim_equil_us,),
            figsize=(3.9, 2.65),
            label_threshold_pct=args.summary_label_threshold_pct,
        )

    else:
        plot_bundle(
            dfs_plot, dfs_stats, names, list(base_lipids), tag="All",
            outdir=outdir, args=args, resid_label_map=resid_label_map,
            rdf_dict_main=rdf_dict, make_combo=args.combo_heatmap_rdf
        )

        summary_all_out = os.path.join(outdir, "hbonds_grouped_stacked_summary_All_last3q.svg")
        plot_grouped_stacked_summary(
            dfs_stats,
            names,
            tuple(args.resids),
            tuple(base_lipids),
            summary_all_out,
            lipid_colors=make_lipid_color_map(base_lipids),
            resid_label_map=resid_label_map,
            label_threshold_pct=args.summary_label_threshold_pct,
            sem_blocks=args.bar_sem_blocks,
        )
        print(f"Wrote {summary_all_out}")

    # -------------------- Time-blocking: hbonds + RDFs --------------------
    if args.time_blocks and args.time_blocks > 1:
        n_blocks = args.time_blocks

        if tops is None or trajs is None:
            raise ValueError("--time-blocks > 1 requires --tops/--trajs so we can define block boundaries from n_frames.")

        hb_blocks_all = {}
        blocks_meta = {}
        for top, traj, sim, df in zip(tops, trajs, names, dfs):
            u = mda.Universe(top, traj)

            hb_blocks, blocks = compute_hbond_blocks_from_df(
                df, sim, u.trajectory.n_frames, n_blocks,
                keep_last_fraction=PLOT_KEEP_LAST_FRACTION
            )
            hb_blocks_all.update(hb_blocks)
            blocks_meta[sim] = blocks

        if args.save_cache is not None:
            save_hbonds_blocks_pickle(args.save_cache, hb_blocks_all, fname=args.hbonds_blocks_pkl)

        rdf_blocks = None
        if args.compute_block_rdfs:
            rdf_blocks = {}
            step_rdf = (args.step if args.step is not None else 1)
            for top, traj, sim in zip(tops, trajs, names):
                u = mda.Universe(top, traj)

                lip_list = list(base_lipids)
                if sim == resim_name and extra_lipid not in lip_list:
                    lip_list.append(extra_lipid)

                start_global = frame_start_for_last_fraction(u.trajectory.n_frames, PLOT_KEEP_LAST_FRACTION)
                stop_global = u.trajectory.n_frames
                blocks = make_time_blocks_from_window(start_global, stop_global, n_blocks)

                for bi, (s, e) in enumerate(blocks, start=1):
                    print(f"[{sim}] Block {bi}/{n_blocks}: frames {s}:{e}")
                    part = compute_rdf_results(
                        tops=[top], trajs=[traj], sim_names=[sim],
                        resids=tuple(args.resids),
                        lipid_resnames=tuple(lip_list),
                        lipid_p_names=lipid_p_names,
                        r_range=(0.0, args.rdf_rmax),
                        nbins=args.rdf_nbins,
                        include_total=True,
                        start=int(s), stop=int(e),
                        step=step_rdf,
                        keep_last_fraction=1.0,
                    )
                    for (sim_k, rid, lip), v in part.items():
                        rdf_blocks[(sim_k, bi, rid, lip)] = v

            if args.save_cache is not None:
                save_rdf_pickle(args.save_cache, rdf_blocks, fname=args.rdf_blocks_pkl)

        if args.load_block_rdfs:
            if args.load_dir is None:
                raise ValueError("--load-block-rdfs requires --load-dir")
            rdf_blocks = load_cached_rdfs(args.load_dir, fname=args.rdf_blocks_pkl)

        if rdf_blocks is not None:
            blocks_outdir = "hbonds_out/block_averages"

            non_resim = [n for n in names if n != resim_name]
            if len(non_resim) != 2:
                raise ValueError("Combo heatmap+RDF figure requires exactly two non-Resim sims (GDP/GTP).")

            sim_a, sim_b = non_resim

            iA = names.index(sim_a)
            iB = names.index(sim_b)

            dfA = dfs_stats[iA]
            dfB = dfs_stats[iB]

            combo_out = os.path.join(
                blocks_outdir,
                f"combo_heatmap_rdf_blockavg_{sim_a}_vs_{sim_b}_N{n_blocks}.svg"
            )

            plot_combo_heatmap_left_plus_rdf_overlay(
                dfs=[dfA, dfB],
                sim_names=[sim_a, sim_b],
                resids=tuple(args.resids),
                lipid_resnames=tuple(base_lipids),
                rdf_blocks=rdf_blocks,
                n_blocks=n_blocks,
                outpath=combo_out,
                lipid_colors=make_lipid_color_map(base_lipids),
                resid_label_map=resid_label_map,
                heat_vmin=args.heatmap_vmin,
                heat_vmax=args.heatmap_vmax,
                rdf_err='sem',
            )

            if (resim_name in names) and (len([n for n in names if n != resim_name]) == 2):
                sim_a, sim_b = [n for n in names if n != resim_name]
                lipid_colors = make_lipid_color_map(base_lipids)

                aggA = aggregate_block_rdfs(rdf_blocks, sim_a, n_blocks, args.resids, base_lipids)
                aggB = aggregate_block_rdfs(rdf_blocks, sim_b, n_blocks, args.resids, base_lipids)

                outp = os.path.join(blocks_outdir, f"rdfs_blockavg_{sim_a}_vs_{sim_b}_N{n_blocks}.svg")
                plot_blockavg_rdfs_overlay_two_sims_per_residue(
                    aggA, aggB, sim_a, sim_b,
                    tuple(args.resids), tuple(base_lipids),
                    outp, lipid_colors=lipid_colors, resid_label_map=resid_label_map,
                    ls_a="-", ls_b=":", band_alpha=0.18
                )

            if resim_name in names:
                lipid_colors_R = make_lipid_color_map(base_lipids + [extra_lipid])
                aggR = aggregate_block_rdfs(rdf_blocks, resim_name, n_blocks, args.resids, base_lipids + [extra_lipid])

                outp = os.path.join(blocks_outdir, f"rdfs_blockavg_{resim_name}_N{n_blocks}.svg")
                plot_blockavg_rdfs_single_sim_per_residue(
                    aggR, resim_name,
                    tuple(args.resids), tuple(base_lipids + [extra_lipid]),
                    outp, lipid_colors=lipid_colors_R, resid_label_map=resid_label_map,
                    ls="-", band_alpha=0.18
                )

if __name__ == "__main__":
    main()
