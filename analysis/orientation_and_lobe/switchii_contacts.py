#!/usr/bin/env python3

import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from MDAnalysis.lib.distances import capped_distance


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
# -------------------- User inputs --------------------

SIMULATIONS = [
    {
        "name": "GDP",
        "topology": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/sys.psf",
        "trajectory": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/cenwrapprod.dcd",
        "dt_ns": 1.0,   # ns between frames on disk
    },
    {
        "name": "GTP-1",
        "topology": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/sys.psf",
        "trajectory": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/cenwrapprod.dcd",
        "dt_ns": 1.0,
    },
    {
        "name": "GTP-2",
        "topology": "/mnt/data_SSD/rhoa_redotraj/nowat.psf",
        "trajectory": "/mnt/data_SSD/rhoa_redotraj/cenwrapprod.dcd",
        "dt_ns": 10.0,   # already strided by 10 on disk
    },
]

# -------------------- Region selections --------------------

REGION_A_NAME = "Switch I"
REGION_A_SELECTION = "protein and backbone and resid 30:43 and not name H*"

REGION_B_NAME = "Helix 4"
REGION_B_SELECTION = "protein and backbone and resid 141:151 and not name H*"

# Contact selections for lobe engagement
LOBE1_NAME = "Lobe 1"
LOBE1_CONTACT_SELECTION = "protein and resid 1:86 and not backbone and not name H*"

LOBE2_NAME = "Lobe 2"
LOBE2_CONTACT_SELECTION = "protein and resid 87:180 and not backbone and not name H*"

CONTACT_LIPID_SELECTION = "resname POPC POPE POPS PLA18 and not name H*"
CONTACT_CUTOFF_ANGSTROM = 4.0

# Membrane selection used to define the midplane COM
MEMBRANE_SELECTION = "resname POPC POPE POPS PLA18 and not name H*"

# -------------------- Plot controls --------------------

# Plot only the last N microseconds of the analyzed trajectory.
# Set to None to plot the full analyzed trace.
PLOT_LAST_US = 15.0

# Displayed time labels will be 5..20 while underlying plotted data remain 0..15
DISPLAY_TIME_START_US = 5.0
DISPLAY_TIME_END_US = 20.0
DISPLAY_TIME_TICKS = np.array([5, 10, 15, 20], dtype=float)

# Plot-only mode
PLOT_ONLY = True
OUTCSV = "region_midplane_zmetrics.csv"
LOADCSV = OUTCSV

# Plot only second half (equilibrium portion)
PLOT_EQUILIBRIUM_ONLY = False

# Analysis stride in ORIGINAL 1-ns frame units
# Desired sampling interval for analysis/plotting, in ns
TARGET_DT_NS = 0.10

# Autocorrelation lag window in ns
ACF_MAX_LAG_NS = 2000.0

# Histogram bins
HIST_BINS = 50

# -------------------- Output files --------------------

OUTFIG_CONTACTS = "lobe_contacts_time_evolution.png"
OUTFIG_CONTACTS_DELTA = "lobe_contacts_delta_centered.png"

OUTFIG_DIST = "region_midplane_distances_equilibrium.png"
OUTFIG_DELTA = "region_midplane_delta_equilibrium.png"
OUTFIG_ACF = "region_midplane_delta_autocorr.png"
OUTFIG_HIST = "region_midplane_delta_histogram.png"
OUTFIG_DELTA_CENTERED = "region_midplane_delta_centered_equilibrium.png"
OUTFIG_ZSCORE = "region_midplane_zscore_equilibrium.png"
OUTFIG_DELTA_ZSCORE = "region_midplane_delta_zscore_equilibrium.png"

# New combined figure
OUTFIG_COMBINED = "region_midplane_and_lobe_contacts_combined.svg"

# -------------------- Figure/layout controls --------------------

FIG_WIDTH_MAIN = 7.2
ROW_HEIGHT_MAIN = 2.5
FIG_WIDTH_HIST = 5.8
ROW_HEIGHT_HIST = 2.3

# Combined figure size
FIG_WIDTH_COMBINED = 5.5
ROW_HEIGHT_COMBINED = 2
COMBINED_WSPACE = 0.12
COMBINED_HSPACE = 0.08

TITLE_FONTSIZE = 14
AXIS_LABEL_FONTSIZE = 14
TICK_LABEL_FONTSIZE = 12

TRACE_LINEWIDTH = .85
TRACE_ALPHA = 0.95
AXHLINE_WIDTH = 1.1

# -------------------- Colors --------------------
# Switch I = Lobe 1 = orange
# Helix 4  = Lobe 2 = blue

COLOR_ORANGE = "#F58518"
COLOR_BLUE = "#4C78A8"
COLOR_GRAY = "#7F7F7F"

REGION_A_COLOR = COLOR_ORANGE
REGION_B_COLOR = COLOR_BLUE
LOBE1_COLOR = COLOR_ORANGE
LOBE2_COLOR = COLOR_BLUE

# -------------------- Style --------------------

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = "truetype"
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

def style_axes(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out", labelsize=TICK_LABEL_FONTSIZE)

def apply_shifted_time_ticks(ax):
    """
    Underlying x data remain 0..15 (zero-started last-15-us window),
    but displayed tick labels are forced to read 5..20 with no left offset.
    """
    internal_ticks = DISPLAY_TIME_TICKS - DISPLAY_TIME_START_US
    internal_end = DISPLAY_TIME_END_US - DISPLAY_TIME_START_US

    ax.set_xlim(0.0, internal_end)
    ax.set_xticks(internal_ticks)
    ax.set_xticklabels([str(int(t)) for t in DISPLAY_TIME_TICKS])

# -------------------- Helpers --------------------

def get_true_stride(sim, target_dt_ns):
    """
    Choose the frame stride needed to achieve approximately target_dt_ns
    given the on-disk frame spacing sim["dt_ns"].
    """
    base_dt_ns = sim["dt_ns"]
    return max(1, int(round(target_dt_ns / base_dt_ns)))

def center_signal(y):
    y = np.asarray(y, dtype=float)
    return y - np.nanmean(y)

def zscore_signal(y):
    y = np.asarray(y, dtype=float)
    mu = np.nanmean(y)
    sd = np.nanstd(y)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(y)
    return (y - mu) / sd

def line_plot(xvalues, yvalues, ax, title, color):
    xvalues = np.asarray(xvalues, dtype=float)
    yvalues = np.asarray(yvalues, dtype=float)

    ax.plot(
        xvalues,
        yvalues,
        color=color,
        linewidth=TRACE_LINEWIDTH,
        alpha=TRACE_ALPHA,
    )
    # ax.set_title(title, fontsize=TITLE_FONTSIZE)
    style_axes(ax)

def get_plot_window(time_us_full, *arrays, last_us=None, zero_start=True):
    """
    Slice one or more arrays to the last `last_us` microseconds using the actual time axis.
    Returns sliced_time, sliced_arrays...
    """
    t = np.asarray(time_us_full, dtype=float)

    if last_us is None:
        mask = np.isfinite(t)
    else:
        tmax = np.nanmax(t)
        tmin = tmax - float(last_us)
        mask = np.isfinite(t) & (t >= tmin)

    t_sel = t[mask]
    out = [np.asarray(arr, dtype=float)[mask] for arr in arrays]

    if zero_start and len(t_sel) > 0:
        t_sel = t_sel - t_sel[0]

    return (t_sel, *out)

def autocorrelation_normalized(y, max_lag=None):
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(y)
    y = y[valid]

    if len(y) < 3:
        return np.array([0]), np.array([np.nan])

    y = y - np.mean(y)
    var = np.var(y)
    if var == 0:
        lags = np.arange(len(y))
        acf = np.zeros(len(y))
        acf[0] = 1.0
        if max_lag is not None:
            lags = lags[:max_lag + 1]
            acf = acf[:max_lag + 1]
        return lags, acf

    corr = np.correlate(y, y, mode="full")
    corr = corr[len(y)-1:]
    norm = var * np.arange(len(y), 0, -1)
    acf = corr / norm
    lags = np.arange(len(acf))

    if max_lag is not None:
        lags = lags[:max_lag + 1]
        acf = acf[:max_lag + 1]

    return lags, acf

def compute_contacts_trace(
    u,
    region_sel,
    lipid_sel,
    cutoff=4.0,
    stride=1,
    show_progress=True,
    progress_desc="Computing contacts",
):
    region = u.select_atoms(region_sel)
    lipids = u.select_atoms(lipid_sel)

    if region.n_atoms == 0:
        raise ValueError(f"Region contact selection returned 0 atoms:\n{region_sel}")
    if lipids.n_atoms == 0:
        raise ValueError(f"Lipid contact selection returned 0 atoms:\n{lipid_sel}")

    n_total = len(u.trajectory)
    frame_indices = np.arange(0, n_total, stride, dtype=int)
    contacts = np.zeros(len(frame_indices), dtype=float)

    iterator = tqdm(
        frame_indices,
        total=len(frame_indices),
        desc=progress_desc,
        unit="frame",
        disable=not show_progress
    )

    for i, frame_idx in enumerate(iterator):
        ts = u.trajectory[frame_idx]

        pairs = capped_distance(
            region.positions,
            lipids.positions,
            max_cutoff=cutoff,
            box=ts.dimensions,
            return_distances=False,
        )
        contacts[i] = len(pairs)

    return frame_indices, contacts

def compute_region_midplane_metrics(
    u,
    region_a_sel,
    region_b_sel,
    membrane_sel,
    stride=1,
    show_progress=True,
    progress_desc="Computing COM z metrics",
):
    region_a = u.select_atoms(region_a_sel)
    region_b = u.select_atoms(region_b_sel)
    membrane = u.select_atoms(membrane_sel)

    if region_a.n_atoms == 0:
        raise ValueError(f"Region A selection returned 0 atoms:\n{region_a_sel}")
    if region_b.n_atoms == 0:
        raise ValueError(f"Region B selection returned 0 atoms:\n{region_b_sel}")
    if membrane.n_atoms == 0:
        raise ValueError(f"Membrane selection returned 0 atoms:\n{membrane_sel}")

    n_total = len(u.trajectory)
    frame_indices = np.arange(0, n_total, stride, dtype=int)

    z_a = np.zeros(len(frame_indices), dtype=float)
    z_b = np.zeros(len(frame_indices), dtype=float)
    z_mid = np.zeros(len(frame_indices), dtype=float)

    iterator = tqdm(
        frame_indices,
        total=len(frame_indices),
        desc=progress_desc,
        unit="frame",
        disable=not show_progress
    )

    for i, frame_idx in enumerate(iterator):
        u.trajectory[frame_idx]

        z_a[i] = region_a.center_of_mass()[2]
        z_b[i] = region_b.center_of_mass()[2]
        z_mid[i] = membrane.center_of_mass()[2]

    # Absolute values for easier viewing
    d_a = np.abs(z_a - z_mid)
    d_b = np.abs(z_b - z_mid)
    delta_z = d_a - d_b

    return frame_indices, z_a, z_b, z_mid, d_a, d_b, delta_z

def save_combined_csv(all_results, outcsv):
    max_len = max(len(v["delta_z_full"]) for v in all_results.values())
    csv_data = {}

    for name, result in all_results.items():
        n = len(result["delta_z_full"])
        pad = max_len - n

        for key in [
            "frame_indices", "time_us_full", "z_a_full", "z_b_full", "z_mid_full",
            "d_a_full", "d_b_full", "delta_z_full",
            "lobe1_contacts_full", "lobe2_contacts_full"
        ]:
            col = np.asarray(result[key], dtype=float)
            csv_data[f"{name}_{key}"] = np.pad(col, (0, pad), constant_values=np.nan)

    pd.DataFrame(csv_data).to_csv(outcsv, index=False)

def load_combined_csv(incsv, sim_names=None):
    df = pd.read_csv(incsv)

    if sim_names is None:
        sim_names = [
            col.replace("_delta_z_full", "")
            for col in df.columns
            if col.endswith("_delta_z_full")
        ]

    all_results = {}

    keys = [
        "frame_indices", "time_us_full", "z_a_full", "z_b_full", "z_mid_full",
        "d_a_full", "d_b_full", "delta_z_full",
        "lobe1_contacts_full", "lobe2_contacts_full"
    ]

    for name in sim_names:
        expected_cols = [f"{name}_{k}" for k in keys]
        missing = [c for c in expected_cols if c not in df.columns]
        if missing:
            raise ValueError(
                f"Missing expected columns for {name}: {missing}\n"
                f"Available columns:\n{list(df.columns)}"
            )

        result = {}
        for k in keys:
            result[k] = df[f"{name}_{k}"].to_numpy(dtype=float)

        valid = ~np.isnan(result["time_us_full"])
        for k in keys:
            result[k] = result[k][valid]

        result["half_idx"] = len(result["delta_z_full"]) // 2

        if len(result["time_us_full"]) > 1:
            result["effective_dt_ns"] = (result["time_us_full"][1] - result["time_us_full"][0]) * 1000.0
        else:
            result["effective_dt_ns"] = np.nan

        all_results[name] = result

    return all_results, list(sim_names)

# -------------------- Plot data prep helpers --------------------

def get_midplane_series(result):
    half_idx = result["half_idx"]

    if PLOT_EQUILIBRIUM_ONLY:
        d_a = result["d_a_full"][half_idx:]
        d_b = result["d_b_full"][half_idx:]
        t = result["time_us_full"][half_idx:]
    else:
        d_a = result["d_a_full"]
        d_b = result["d_b_full"]
        t = result["time_us_full"]

    x, d_a, d_b = get_plot_window(
        t,
        d_a,
        d_b,
        last_us=PLOT_LAST_US,
        zero_start=True
    )
    return x, d_a, d_b

def get_contact_series(result):
    half_idx = result["half_idx"]

    if PLOT_EQUILIBRIUM_ONLY:
        c1 = result["lobe1_contacts_full"][half_idx:]
        c2 = result["lobe2_contacts_full"][half_idx:]
        t = result["time_us_full"][half_idx:]
    else:
        c1 = result["lobe1_contacts_full"]
        c2 = result["lobe2_contacts_full"]
        t = result["time_us_full"]

    x, c1, c2 = get_plot_window(
        t,
        c1,
        c2,
        last_us=PLOT_LAST_US,
        zero_start=True
    )
    return x, c1, c2

# -------------------- Main --------------------

def main():
    print("Region A selection:", REGION_A_SELECTION)
    print("Region B selection:", REGION_B_SELECTION)
    print("Membrane selection:", MEMBRANE_SELECTION)
    print("Target dt (ns):    ", TARGET_DT_NS)
    print("PLOT_ONLY:         ", PLOT_ONLY)

    if PLOT_ONLY:
        print(f"\nPLOT_ONLY=True, loading saved data from: {LOADCSV}")
        all_results, sim_order = load_combined_csv(LOADCSV)
        print("Loaded sim order from CSV:", sim_order)

    else:
        sim_order = [sim["name"] for sim in SIMULATIONS]
        all_results = {}

        for sim in SIMULATIONS:
            name = sim["name"]
            top = sim["topology"]
            traj = sim["trajectory"]
            base_dt_ns = sim["dt_ns"]

            true_stride = get_true_stride(sim, TARGET_DT_NS)
            effective_dt_ns = base_dt_ns * true_stride

            print(f"\nLoading {name}")
            u = mda.Universe(top, traj)

            frame_indices, z_a, z_b, z_mid, d_a, d_b, delta_z = compute_region_midplane_metrics(
                u=u,
                region_a_sel=REGION_A_SELECTION,
                region_b_sel=REGION_B_SELECTION,
                membrane_sel=MEMBRANE_SELECTION,
                stride=true_stride,
                show_progress=True,
                progress_desc=f"{name} COM z metrics",
            )

            frame_idx_l1, lobe1_contacts = compute_contacts_trace(
                u=u,
                region_sel=LOBE1_CONTACT_SELECTION,
                lipid_sel=CONTACT_LIPID_SELECTION,
                cutoff=CONTACT_CUTOFF_ANGSTROM,
                stride=true_stride,
                show_progress=True,
                progress_desc=f"{name} {LOBE1_NAME} contacts",
            )

            frame_idx_l2, lobe2_contacts = compute_contacts_trace(
                u=u,
                region_sel=LOBE2_CONTACT_SELECTION,
                lipid_sel=CONTACT_LIPID_SELECTION,
                cutoff=CONTACT_CUTOFF_ANGSTROM,
                stride=true_stride,
                show_progress=True,
                progress_desc=f"{name} {LOBE2_NAME} contacts",
            )

            if not (np.array_equal(frame_indices, frame_idx_l1) and np.array_equal(frame_indices, frame_idx_l2)):
                raise ValueError(f"{name}: frame indices do not match between COM and contact analyses.")

            half_idx = len(frame_indices) // 2
            time_us_full = np.arange(len(frame_indices)) * effective_dt_ns / 1000.0

            print(f"{name}: true stride used = {true_stride}")
            print(f"{name}: analyzed frames = {len(frame_indices)}")
            print(f"{name}: base dt = {base_dt_ns} ns")
            print(f"{name}: effective dt = {effective_dt_ns} ns")

            all_results[name] = {
                "frame_indices": frame_indices,
                "time_us_full": time_us_full,
                "z_a_full": z_a,
                "z_b_full": z_b,
                "z_mid_full": z_mid,
                "d_a_full": d_a,
                "d_b_full": d_b,
                "delta_z_full": delta_z,
                "half_idx": half_idx,
                "effective_dt_ns": effective_dt_ns,
                "lobe1_contacts_full": lobe1_contacts,
                "lobe2_contacts_full": lobe2_contacts,
            }

        save_combined_csv(all_results, OUTCSV)
        print(f"\nSaved CSV: {OUTCSV}")

    n_sims = len(sim_order)

    # -------------------- Plot 1: distances to midplane --------------------

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=(FIG_WIDTH_MAIN, ROW_HEIGHT_MAIN * n_sims),
        sharex=True,
        sharey=True
    )
    if n_sims == 1:
        axes = [axes]

    for ax, name in zip(axes, sim_order):
        result = all_results[name]
        x, d_a, d_b = get_midplane_series(result)

        line_plot(x, d_a, ax, name, REGION_A_COLOR)
        line_plot(x, d_b, ax, name, REGION_B_COLOR)
        ax.set_ylabel("Z (Å)", fontsize=AXIS_LABEL_FONTSIZE)

    axes[-1].set_xlabel("Time (µs)", fontsize=AXIS_LABEL_FONTSIZE)
    for ax in axes:
        apply_shifted_time_ticks(ax)
    plt.tight_layout()
    plt.savefig(OUTFIG_DIST, dpi=300)
    print(f"Saved figure: {OUTFIG_DIST}")

    # -------------------- Plot 2: centered delta_z --------------------

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=(FIG_WIDTH_MAIN, ROW_HEIGHT_MAIN * n_sims),
        sharex=True,
        sharey=True
    )
    if n_sims == 1:
        axes = [axes]

    for ax, name in zip(axes, sim_order):
        result = all_results[name]
        x, d_a, d_b = get_midplane_series(result)

        d_a_center = center_signal(d_a)
        d_b_center = center_signal(d_b)
        delta_center = d_a_center - d_b_center

        line_plot(
            x, delta_center, ax,
            f"{name} : centered {REGION_A_NAME} - {REGION_B_NAME}",
            COLOR_GRAY
        )

        ax.axhline(0.0, color="black", linewidth=AXHLINE_WIDTH, linestyle="--", alpha=0.7)
        ax.set_ylabel("Centered Δz (Å)", fontsize=AXIS_LABEL_FONTSIZE)

    axes[-1].set_xlabel("Time (µs)", fontsize=AXIS_LABEL_FONTSIZE)
    for ax in axes:
        apply_shifted_time_ticks(ax)
    plt.tight_layout()
    plt.savefig(OUTFIG_DELTA_CENTERED, dpi=300)
    print(f"Saved figure: {OUTFIG_DELTA_CENTERED}")

    # -------------------- Plot 3: z-scored region distances --------------------

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=(FIG_WIDTH_MAIN, ROW_HEIGHT_MAIN * n_sims),
        sharex=True,
        sharey=True
    )
    if n_sims == 1:
        axes = [axes]

    for ax, name in zip(axes, sim_order):
        result = all_results[name]
        x, d_a, d_b = get_midplane_series(result)

        d_a_z = zscore_signal(d_a)
        d_b_z = zscore_signal(d_b)

        line_plot(x, d_a_z, ax, name, REGION_A_COLOR)
        line_plot(x, d_b_z, ax, name, REGION_B_COLOR)

        ax.set_ylabel("Z-scored distance", fontsize=AXIS_LABEL_FONTSIZE)

    axes[-1].set_xlabel("Time (µs)", fontsize=AXIS_LABEL_FONTSIZE)
    for ax in axes:
        apply_shifted_time_ticks(ax)
    plt.tight_layout()
    plt.savefig(OUTFIG_ZSCORE, dpi=300)
    print(f"Saved figure: {OUTFIG_ZSCORE}")

    # -------------------- Plot 4: autocorrelation of delta_z --------------------

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=(FIG_WIDTH_MAIN, ROW_HEIGHT_MAIN * n_sims),
        sharex=True,
        sharey=True
    )
    if n_sims == 1:
        axes = [axes]

    for ax, name in zip(axes, sim_order):
        result = all_results[name]
        half_idx = result["half_idx"]
        effective_dt_ns = result["effective_dt_ns"]

        if PLOT_EQUILIBRIUM_ONLY:
            delta_z = result["delta_z_full"][half_idx:]
        else:
            delta_z = result["delta_z_full"]

        max_lag_frames = max(1, int(round(ACF_MAX_LAG_NS / effective_dt_ns)))
        lags, acf = autocorrelation_normalized(delta_z, max_lag=max_lag_frames)
        lag_us = (lags * effective_dt_ns) / 1000.0

        ax.plot(lag_us, acf, color="black", linewidth=TRACE_LINEWIDTH)
        ax.axhline(0.0, color="black", linewidth=AXHLINE_WIDTH, linestyle="--", alpha=0.7)
        ax.set_title(f"{name} : autocorrelation of Δz", fontsize=TITLE_FONTSIZE)
        ax.set_ylabel("ACF", fontsize=AXIS_LABEL_FONTSIZE)
        style_axes(ax)

    axes[-1].set_xlabel("Lag (µs)", fontsize=AXIS_LABEL_FONTSIZE)
    plt.tight_layout()
    plt.savefig(OUTFIG_ACF, dpi=300)
    print(f"Saved figure: {OUTFIG_ACF}")

    # -------------------- Plot 5: histogram of delta_z --------------------

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=(FIG_WIDTH_HIST, ROW_HEIGHT_HIST * n_sims),
        sharex=True,
        sharey=True
    )
    if n_sims == 1:
        axes = [axes]

    for ax, name in zip(axes, sim_order):
        result = all_results[name]
        half_idx = result["half_idx"]

        if PLOT_EQUILIBRIUM_ONLY:
            delta_z = result["delta_z_full"][half_idx:]
        else:
            delta_z = result["delta_z_full"]

        ax.hist(delta_z[np.isfinite(delta_z)], bins=HIST_BINS, density=True, color="gray", alpha=0.85)
        ax.axvline(0.0, color="black", linewidth=AXHLINE_WIDTH, linestyle="--", alpha=0.7)
        ax.set_title(f"{name} : distribution of Δz", fontsize=TITLE_FONTSIZE)
        ax.set_ylabel("Density", fontsize=AXIS_LABEL_FONTSIZE)
        style_axes(ax)

    axes[-1].set_xlabel("Δz (Å)", fontsize=AXIS_LABEL_FONTSIZE)
    plt.tight_layout()
    plt.savefig(OUTFIG_HIST, dpi=300)
    print(f"Saved figure: {OUTFIG_HIST}")

    # -------------------- Plot 6: lobe sidechain contacts --------------------

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=(FIG_WIDTH_MAIN, ROW_HEIGHT_MAIN * n_sims),
        sharex=True,
        sharey=True
    )
    if n_sims == 1:
        axes = [axes]

    for ax, name in zip(axes, sim_order):
        result = all_results[name]
        x, c1, c2 = get_contact_series(result)

        line_plot(x, c1, ax, name, LOBE1_COLOR)
        line_plot(x, c2, ax, name, LOBE2_COLOR)
        ax.set_ylabel("# contacts", fontsize=AXIS_LABEL_FONTSIZE)

    axes[-1].set_xlabel("Time (µs)", fontsize=AXIS_LABEL_FONTSIZE)
    for ax in axes:
        apply_shifted_time_ticks(ax)
    plt.tight_layout()
    plt.savefig(OUTFIG_CONTACTS, dpi=300)
    print(f"Saved figure: {OUTFIG_CONTACTS}")

    # -------------------- Plot 7: centered lobe contact delta --------------------

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=(FIG_WIDTH_MAIN, ROW_HEIGHT_MAIN * n_sims),
        sharex=True,
        sharey=True
    )
    if n_sims == 1:
        axes = [axes]

    for ax, name in zip(axes, sim_order):
        result = all_results[name]
        x, c1, c2 = get_contact_series(result)

        c1_center = center_signal(c1)
        c2_center = center_signal(c2)
        delta_contacts = c1_center - c2_center

        line_plot(
            x,
            delta_contacts,
            ax,
            f"{name} : centered {LOBE1_NAME} - {LOBE2_NAME}",
            COLOR_GRAY
        )

        ax.axhline(0.0, color="black", linewidth=AXHLINE_WIDTH, linestyle="--", alpha=0.7)
        ax.set_ylabel("Centered Δ contacts", fontsize=AXIS_LABEL_FONTSIZE)

    axes[-1].set_xlabel("Time (µs)", fontsize=AXIS_LABEL_FONTSIZE)
    for ax in axes:
        apply_shifted_time_ticks(ax)
    plt.tight_layout()
    plt.savefig(OUTFIG_CONTACTS_DELTA, dpi=300)
    print(f"Saved figure: {OUTFIG_CONTACTS_DELTA}")

    # -------------------- Plot 8: combined two-column figure --------------------
    # Left column: region midplane distances
    # Right column: lobe contact time evolution

    fig = plt.figure(figsize=(FIG_WIDTH_COMBINED*1.2, ROW_HEIGHT_COMBINED * n_sims), constrained_layout=True)
    gs = gridspec.GridSpec(
        nrows=n_sims,
        ncols=2,
        figure=fig,
        wspace=COMBINED_WSPACE,
        hspace=COMBINED_HSPACE
    )

    left_axes = []
    right_axes = []

    for i, name in enumerate(sim_order):
        ax_left = fig.add_subplot(gs[i, 0], sharex=left_axes[0] if left_axes else None)
        ax_right = fig.add_subplot(gs[i, 1], sharex=right_axes[0] if right_axes else None)

        left_axes.append(ax_left)
        right_axes.append(ax_right)

        result = all_results[name]

        # Left column: midplane distances
        x_mid, d_a, d_b = get_midplane_series(result)
        line_plot(x_mid, d_a, ax_left, name, REGION_A_COLOR)
        line_plot(x_mid, d_b, ax_left, name, REGION_B_COLOR)
        ax_left.set_ylabel("Z (Å)", fontsize=AXIS_LABEL_FONTSIZE)

        # Right column: lobe contacts
        x_con, c1, c2 = get_contact_series(result)
        line_plot(x_con, c1, ax_right, name, LOBE1_COLOR)
        line_plot(x_con, c2, ax_right, name, LOBE2_COLOR)
        ax_right.set_ylabel("Contacts", fontsize=AXIS_LABEL_FONTSIZE)

        # Only bottom row gets x labels
        if i < n_sims - 1:
            ax_left.tick_params(labelbottom=False)
            ax_right.tick_params(labelbottom=False)

    # Column headers
    # left_axes[0].text(
    #     0.5, 1.18, "Region midplane distances",
    #     transform=left_axes[0].transAxes,
    #     ha="center", va="bottom",
    #     fontsize=TITLE_FONTSIZE + 1
    # )
    # right_axes[0].text(
    #     0.5, 1.18, "Lobe contact time evolution",
    #     transform=right_axes[0].transAxes,
    #     ha="center", va="bottom",
    #     fontsize=TITLE_FONTSIZE + 1
    # )

    # Apply shifted time ticks to all panels
    for ax in left_axes + right_axes:
        apply_shifted_time_ticks(ax)

    left_axes[-1].set_xlabel("Time (µs)", fontsize=AXIS_LABEL_FONTSIZE)
    right_axes[-1].set_xlabel("Time (µs)", fontsize=AXIS_LABEL_FONTSIZE)

    
    plt.savefig(OUTFIG_COMBINED, dpi=300)
    print(f"Saved figure: {OUTFIG_COMBINED}")

if __name__ == "__main__":
    main()
