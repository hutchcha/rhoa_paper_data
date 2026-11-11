#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib.patches as mpatches

import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from MDAnalysis.lib.distances import capped_distance
import matplotlib

# -----------------------
# Matplotlib style
# -----------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2

# Keep text as text (not paths) in vector exports
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Arial"
matplotlib.rcParams["mathtext.it"] = "Arial:italic"
matplotlib.rcParams["mathtext.bf"] = "Arial:bold"

# -----------------------
# User settings
# -----------------------
STEP = 1
CUTOFF = 4.0  # Angstrom
N_JOBS = max(1, cpu_count() - 1)
CHUNK_SIZE = 200  # frames per worker task

# Load existing CSV instead of recomputing
LOAD_FROM_CSV = True

# Selections
PROT_HEAVY_SEL = "(protein and not name H*) or (resname CYSF CYSG and not name H*)"
LIPID_HEAVY_SEL = "(resname POPC POPE CHL1 SAPI POPS) and not name H*"

# Output
OUT_CSV = "./Contacts/rcf_all_sims.csv"
PLOT_DIR = "./Contacts/plots"
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

# ---- figure font/tick preset ----
FS_SUP = 14     # supxlabel/supylabel
FS_ROW = 11     # row label box text (sim name)
FS_TICK = 11    # tick labels
FS_BAND = 10    # band labels

BAND_ALPHA = 0.22
SPINE_W = 1.25
TICK_W = 1.0
TICK_L = 3.5

# Bar styling
BAR_WIDTH = 0.75
BAR_COLOR = "black"

# -----------------------
# Regions to shade
# -----------------------
REGIONS = [
    ("NT",   1,   7,   "cyan"),
    ("SwI",  30,  41,  "yellow"),
    ("SwII", 56,  76,  "teal"),
    ("HVR",  180, 193, "blue"),
]

# -----------------------
# Residue callouts for top panel
# -----------------------
# tex_label can use mathtext, e.g. r"$K_{6}$" or r"$\\mathrm{Arg}^{5}$"
TOP_PANEL_CALLOUTS = [
    {"resid": 5, "tex_label": r"R$5$", "dx":  0.0, "dy": 0.4},
    {"resid": 6, "tex_label": r"K$6$", "dx":  5.0, "dy": 0.35},
    {"resid": 7, "tex_label": r"K$7$", "dx":  10.0, "dy": 0.24},
    {"resid": 51, "tex_label": r"K$51$", "dx":  0.0, "dy": 0.35},
    {"resid": 76, "tex_label": r"D$76$", "dx":  0.0, "dy": 0.36}
]

CALLOUT_FONTSIZE = 11
CALLOUT_LINEWIDTH = 1.1

# -----------------------
# Resim equilibration settings
# -----------------------
RESIM_NAME = "RhoA-Resim"
RESIM_EQ_START_US = 5.0  # analyze only t >= 5 us for Resim


def add_structure_bands(ax, alpha=BAND_ALPHA):
    """Add filled residue bands + return legend handles."""
    handles = []
    for label, start, end, color in REGIONS:
        ax.axvspan(start, end, ymin=0, ymax=1, facecolor=color, alpha=alpha, edgecolor=None)
        handles.append(mpatches.Patch(facecolor=color, alpha=alpha, label=label))
    return handles


def add_band_labels(ax, y=1.12, fontsize=FS_BAND, rotation=0, x_offset=0.1):
    """Label each shaded region above the top panel."""
    for label, start, end, color in REGIONS:
        x_mid = 0.5 * (start + end) + x_offset
        ax.text(
            x_mid,
            y,
            label,
            transform=ax.get_xaxis_transform(),  # x in data coords, y in axes coords
            ha="center",
            va="top",
            fontsize=fontsize,
            color=color,
            rotation=rotation,
            clip_on=False,
        )


def add_residue_callouts(ax, df, callouts, y_top, fontsize=CALLOUT_FONTSIZE):
    """
    Add TeX labels and leader lines to selected residues on the top panel.
    """
    df_idx = df.set_index("resid")

    for item in callouts:
        resid = int(item["resid"])
        tex_label = item.get("tex_label", rf"${resid}$")
        dx = float(item.get("dx", 0.0))
        dy = float(item.get("dy", 0.15))

        if resid not in df_idx.index:
            continue

        y_val = float(df_idx.loc[resid, "rcf"])
        x_text = resid + dx
        y_text = min(y_top - 0.03, y_val + dy)

        ax.annotate(
            tex_label,
            xy=(resid, y_val),
            xytext=(x_text, y_text),
            textcoords="data",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            arrowprops=dict(
                arrowstyle="-",
                lw=CALLOUT_LINEWIDTH,
                color="black",
                shrinkA=0,
                shrinkB=0,
            ),
            clip_on=False,
        )


def _infer_traj_dt_ps(u):
    """
    Best-effort infer timestep (ps) between frames.
    Uses MDAnalysis trajectory dt when available, else falls back to time differences.
    """
    dt = getattr(u.trajectory, "dt", None)
    if dt is not None and np.isfinite(dt) and dt > 0:
        return float(dt)  # MDAnalysis typically stores dt in ps

    try:
        if u.trajectory.n_frames >= 2:
            u.trajectory[0]
            t0 = float(getattr(u.trajectory.ts, "time", 0.0))
            u.trajectory[1]
            t1 = float(getattr(u.trajectory.ts, "time", 0.0))
            if np.isfinite(t1 - t0) and (t1 - t0) > 0:
                return float(t1 - t0)
    except Exception:
        pass

    raise ValueError(
        "Could not infer trajectory dt (ps). "
        "MDAnalysis did not provide trajectory.dt and time difference fallback failed."
    )


def _frame_index_for_time_us(u, target_us):
    """
    Convert a target time in microseconds to the first frame index at/after that time.
    Assumes ts.time (and dt) are in ps.
    """
    target_ps = float(target_us) * 1e6  # 1 us = 1e6 ps

    u.trajectory[0]
    t0 = float(getattr(u.trajectory.ts, "time", 0.0))  # ps
    dt_ps = _infer_traj_dt_ps(u)

    if target_ps <= t0:
        return 0

    idx = int(np.ceil((target_ps - t0) / dt_ps))
    idx = max(0, min(idx, u.trajectory.n_frames - 1))
    return idx


# -----------------------
# Worker: process chunk of frames for ONE trajectory
# -----------------------
def _chunk_contact_counts(args):
    """
    Returns:
      contacted_counts (np.ndarray [n_res]): # sampled frames where residue had >=1 contact
      n_frames (int): sampled frames processed
      unique_resids (np.ndarray [n_res]): residue ids in worker ordering
    """
    top, traj, step, cutoff, frame_indices = args

    u = mda.Universe(top, traj)
    prot = u.select_atoms(PROT_HEAVY_SEL)
    lip = u.select_atoms(LIPID_HEAVY_SEL)

    prot_resids = prot.resids
    unique_resids = np.unique(prot_resids)
    resid_to_idx = {rid: i for i, rid in enumerate(unique_resids)}
    prot_atom_res_idx = np.array([resid_to_idx[rid] for rid in prot_resids], dtype=np.int32)

    n_res = unique_resids.size
    contacted_counts = np.zeros(n_res, dtype=np.int64)
    n_frames = 0

    for fi in frame_indices:
        u.trajectory[fi]

        pairs = capped_distance(
            prot.positions,
            lip.positions,
            max_cutoff=cutoff,
            box=u.dimensions,
            return_distances=False,
        )

        if pairs is not None and len(pairs) > 0:
            if isinstance(pairs, tuple):
                prot_atom_idx = pairs[0]
            else:
                pairs = np.asarray(pairs)
                prot_atom_idx = pairs[:, 0]

            res_idx = prot_atom_res_idx[prot_atom_idx]
            contacted_counts[np.unique(res_idx)] += 1

        n_frames += 1

    return contacted_counts, n_frames, unique_resids


# -----------------------
# Compute RCF for one trajectory
# -----------------------
def compute_rcf_for_traj(
    top,
    traj,
    step=STEP,
    cutoff=CUTOFF,
    n_jobs=N_JOBS,
    chunk_size=CHUNK_SIZE,
    start_frame=None,
):
    """
    Returns DataFrame columns:
      resid, resname, label, rcf, n_frames_sampled
    """
    u0 = mda.Universe(top, traj)
    prot0 = u0.select_atoms(PROT_HEAVY_SEL)

    residues = prot0.residues
    meta = {int(r.resid): (r.resname, f"{r.resname}-{int(r.resid)}") for r in residues}

    if start_frame is None:
        start_frame = 0
    start_frame = int(max(0, start_frame))

    all_frames = np.arange(start_frame, u0.trajectory.n_frames, step, dtype=np.int64)
    if all_frames.size == 0:
        raise ValueError(
            f"No frames selected (n_frames={u0.trajectory.n_frames}, "
            f"start_frame={start_frame}, step={step})."
        )

    chunks = [all_frames[i : i + chunk_size] for i in range(0, all_frames.size, chunk_size)]
    args_list = [(top, traj, step, cutoff, chunk) for chunk in chunks]

    if n_jobs == 1:
        results = [_chunk_contact_counts(a) for a in args_list]
    else:
        with Pool(processes=n_jobs) as pool:
            results = pool.map(_chunk_contact_counts, args_list)

    total_counts = None
    total_frames = 0
    worker_unique_resids = results[0][2]

    for counts, n_frames, _ur in results:
        total_frames += n_frames
        if total_counts is None:
            total_counts = counts.copy()
        else:
            total_counts += counts

    rcf = total_counts / float(total_frames)

    out_resids = worker_unique_resids.astype(int)
    out_resnames = []
    out_labels = []
    for rid in out_resids:
        rn, lab = meta.get(int(rid), ("UNK", f"UNK-{int(rid)}"))
        out_resnames.append(rn)
        out_labels.append(lab)

    df = pd.DataFrame(
        {
            "resid": out_resids,
            "resname": out_resnames,
            "label": out_labels,
            "rcf": rcf,
            "n_frames_sampled": int(total_frames),
        }
    )
    return df


# -----------------------
# Plotting
# -----------------------
def plot_rcf_stacked_sims(all_df, sim_order, outpath, bar_color=BAR_COLOR):
    nrows = len(sim_order)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=1,
        figsize=(10.0, 3.5),
        sharex=True,
        sharey=True,
    )
    if nrows == 1:
        axes = [axes]

    fig.supylabel("RCF", fontsize=FS_SUP)
    fig.supxlabel("Residue #", fontsize=FS_SUP)

    max_resid = int(all_df["resid"].max())
    xtick_step = 20 if max_resid >= 200 else 10
    xticks = np.arange(0, max_resid + xtick_step, xtick_step)

    global_max = float(np.nanmax(all_df["rcf"].to_numpy()))
    y_top = max(1.2, global_max + 0.35)

    for ax_i, (ax, sim_name) in enumerate(zip(axes, sim_order)):
        df = all_df[all_df["sim"] == sim_name].copy()
        if df.empty:
            raise ValueError(f"No rows found for sim='{sim_name}'.")

        df = df.sort_values("resid")
        x = df["resid"].to_numpy()
        y = df["rcf"].to_numpy()

        ax.set_ylim(0, y_top)

        add_structure_bands(ax, alpha=BAND_ALPHA)
        ax.bar(x, y, color=bar_color, width=BAR_WIDTH)

        if ax is axes[0]:
            add_band_labels(ax, y=1.12, fontsize=FS_BAND)
            add_residue_callouts(ax, df, TOP_PANEL_CALLOUTS, y_top=y_top)

        bbox_props = dict(boxstyle="round,pad=0.25", fc="white", ec="black")
        ax.text(
            1.012,
            0.50,
            sim_name,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=FS_ROW,
            bbox=bbox_props,
            clip_on=False,
        )

        ax.grid(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.yaxis.tick_left()
        ax.xaxis.tick_bottom()
        ax.minorticks_off()

        ax.grid(axis="y", linestyle=":")

        ax.spines["left"].set_linewidth(SPINE_W)
        ax.spines["bottom"].set_linewidth(SPINE_W)

        ax.tick_params(axis="both", labelsize=FS_TICK, width=TICK_W, length=TICK_L, direction="out")

        ax.set_yticks(np.arange(0, 1.5, 0.5))
        ax.set_xticks(xticks)
        ax.set_xlim(0, max_resid + 5)

    fig.tight_layout(pad=0.35)
    fig.subplots_adjust(hspace=0.23, right=0.88)

    fig.savefig(outpath, dpi=300)
    plt.close(fig)


# -----------------------
# Main
# -----------------------
def main():
    sims = {
        "RhoA-GDP": ("/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/sys.psf", "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/cenwrapprod.dcd"),
        "RhoA-GTP": ("/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/sys.psf", "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/cenwrapprod.dcd"),
        "RhoA-Resim": ("/mnt/data_SSD/rhoa_redotraj/nowat.psf", "/mnt/data_SSD/rhoa_redotraj/cenwrapprod.dcd"),
    }
    sim_order = ["RhoA-GDP", "RhoA-GTP", "RhoA-Resim"]

    if LOAD_FROM_CSV and os.path.exists(OUT_CSV):
        print(f"Loading existing CSV: {OUT_CSV}")
        all_df = pd.read_csv(OUT_CSV)

        required_cols = {"sim", "resid", "resname", "label", "rcf", "n_frames_sampled"}
        missing = required_cols.difference(all_df.columns)
        if missing:
            raise ValueError(f"CSV is missing required columns: {sorted(missing)}")
    else:
        rows = []
        for sim_name, (top, traj) in sims.items():
            if not os.path.exists(top):
                raise FileNotFoundError(f"Missing topology for {sim_name}: {top}")
            if not os.path.exists(traj):
                raise FileNotFoundError(f"Missing trajectory for {sim_name}: {traj}")

            start_frame = None
            if sim_name == RESIM_NAME:
                u_tmp = mda.Universe(top, traj)
                start_frame = 2481
                print(
                    f"[{sim_name}] Using equilibrium window: "
                    f"t >= {RESIM_EQ_START_US} us (start_frame={start_frame})"
                )
                del u_tmp

            df = compute_rcf_for_traj(
                top,
                traj,
                step=STEP,
                cutoff=CUTOFF,
                n_jobs=N_JOBS,
                chunk_size=CHUNK_SIZE,
                start_frame=start_frame,
            )
            df.insert(0, "sim", sim_name)
            rows.append(df)

        all_df = pd.concat(rows, ignore_index=True)
        all_df.to_csv(OUT_CSV, index=False)
        print(f"Wrote CSV: {OUT_CSV}")

    out_svg = os.path.join(PLOT_DIR, "rcf_rhoa_stacked_sims.svg")
    plot_rcf_stacked_sims(all_df, sim_order=sim_order, outpath=out_svg)
    print(f"Wrote plot: {out_svg}")


if __name__ == "__main__":
    main()
