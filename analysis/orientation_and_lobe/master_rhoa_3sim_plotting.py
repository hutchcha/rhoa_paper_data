#!/usr/bin/env python3

import pickle
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde
from tqdm import tqdm

try:
    import lipyphilic
except ImportError:
    lipyphilic = None

# =========================
# Global plotting style
# =========================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Arial"
matplotlib.rcParams["mathtext.it"] = "Arial:italic"
matplotlib.rcParams["mathtext.bf"] = "Arial:bold"

# =========================
# USER CONFIG
# =========================
TOTAL_TIME_US = 20.0
MAIN_MA_WINDOW = 100
FIG_OUTDIR = Path("RhoA_validation_figures")
MAIN_FIG = Path("RhoA_validation_summary.svg")
INDIVIDUAL_VALIDATION_FMT = "svg"
EXTRA_FMT = "svg"

# Easy font-size controls
TITLE_FONTSIZE = 21
AXIS_LABEL_FONTSIZE = 15
TICK_LABEL_FONTSIZE = 15
ROW_LABEL_FONTSIZE = 16
LEGEND_FONTSIZE = 10
PANEL_LABEL_FONTSIZE = 20
POLAR_TITLE_FONTSIZE = 16
POLAR_TICK_FONTSIZE = 11
LOBE_AXIS_LABEL_FONTSIZE = 13

# Layout controls
MAIN_FIG_WIDTH = 16.8
MAIN_FIG_HEIGHT = 17.5
ROW_HSPACE = 0.48
COL_WSPACE = 0.22

RHOA_SYSTEMS = {
    "RhoA-GDP": {
        "psf": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/sys.psf",
        "traj": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/cenwrapprod.dcd",
        "anchor_sel": "resname CYSG and not name CA HA C CB SG CT OT1 OT2 O H*",
        "anchor_rmsd_sel": "protein and name CA and not resid 25:43 58:78 120:140 178:191",
        "lobe1_sel": "protein and resid 5:86 and name CA",
        "lobe2_sel": "protein and resid 95:177 and name CA",
        "rotation_txt": "Tilt and Rotation/RhoA_GDProtationangles.txt",
        "tilt_txt": "Tilt and Rotation/RhoA_GDPtiltangles.txt",
        "color": "black",
        "stride_apl": 50,
        "stride_rmsd": 20,
        "stride_pp": 10,
        "stride_anchor": 10,
        "stride_lobe": 10,
    },
    "RhoA-GTP": {
        "psf": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/sys.psf",
        "traj": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/cenwrapprod.dcd",
        "anchor_sel": "resname CYSG and not name CA HA C CB SG CT OT1 OT2 O H*",
        "anchor_rmsd_sel": "protein and name CA and not resid 25:43 58:78 120:140 178:191",
        "lobe1_sel": "protein and resid 5:86 and name CA",
        "lobe2_sel": "protein and resid 95:177 and name CA",
        "rotation_txt": "Tilt and Rotation/RhoA_GTProtationangles.txt",
        "tilt_txt": "Tilt and Rotation/RhoA_GTPtiltangles.txt",
        "color": "tab:blue",
        "stride_apl": 50,
        "stride_rmsd": 20,
        "stride_pp": 10,
        "stride_anchor": 10,
        "stride_lobe": 10,
    },
    "RhoA-GTP-2": {
        "psf": "/mnt/data_SSD/rhoa_redotraj/nowat.psf",
        "traj": "/mnt/data_SSD/rhoa_redotraj/cenwrapprod.dcd",
        "anchor_sel": "resname CYSG and not name CA HA C CB SG CT OT1 OT2 O H*",
        "anchor_rmsd_sel": "protein and name CA and not resid 25:43 58:78 120:140 178:191",
        "lobe1_sel": "protein and resid 5:86 and name CA",
        "lobe2_sel": "protein and resid 95:177 and name CA",
        "rotation_txt": "Tilt and Rotation/RhoA_GTPredorotationangles.txt",
        "tilt_txt": "Tilt and Rotation/RhoA_GTPredotiltangles.txt",
        "color": "tab:orange",
        "stride_apl": 5,
        "stride_rmsd": 2,
        "stride_pp": 1,
        "stride_anchor": 1,
        "stride_lobe": 1,
    },
}

NAME_ALIASES = {
    "RhoA-GDP": ["RhoA-GDP"],
    "RhoA-GTP": ["RhoA-GTP"],
    "RhoA-GTP-2": ["RhoA-GTP-2", "RhoA-Resim", "RhoA_GTPredo", "rhoaredo", "RhoA redo"],
}

APL_PICKLE = Path("APL/all_sims_apl.pkl")
RMSD_PICKLE = Path("rmsd/rmsdrhoa.pkl")
ANCHOR_PICKLE = Path("anchor_insert_depth.pkl")
PP_PICKLE = Path("pp-dis/rhoa_ppdistances_val.pkl")
LOBE_PICKLE = Path("Lobe Center of mass/rhoa_lobe_data.pkl")

MEMBRANE_P_SEL = "resname POPC POPE PSM PLA18 POPS"
MEMBRANE_NONH_SEL = "resname POPC POPE POPS CHL1 SAPI PSM PLA18 and not name H*"
MEMBRANE_COM_SEL = "resname POPC POPE SAPI CHL1 and name C*"
APL_LIPIDS = ["POPC", "POPE", "POPS", "PLA18", "PSM"]

VALIDATION_ROWS = [
    ("apl", "Area per lipid", "APL (Å$^2$)"),
    ("rmsd", "RMSD", "RMSD (Å)"),
    ("pp", "P-P distance", "Bilayer thickness (Å)"),
    ("anchor", "Anchor insertion depth", r"$\mathrm{I}_{\mathrm{D}}$"),
]

ROW_PANEL_LABELS = ["A", "B", "C", "D"]

lipid_columns = [
    "Total", "POPC", "POPE", "POPS", "PLA18",
    "PSM", "SAPI", "CHL1", "Upper POPC", "Lower POPC"
]

lipid_colors = {
    "Total": "black",
    "POPC": "blue",
    "POPE": "green",
    "POPS": "gold",
    "PLA18": "purple",
    "PSM": "pink",
    "SAPI": "orange",
    "CHL1": "red",
    "Upper POPC": "cyan",
    "Lower POPC": "magenta"
}

# =========================
# Plotting helpers
# =========================
def line_plot(
    xvalues,
    yvalues,
    ax,
    title,
    color,
    z,
    label,
    ma_window=5,
    ma_color='red',
    ma_z=2,
    show_title=True,
    show_x_ticklabels=True,
    title_fontsize=TITLE_FONTSIZE,
    tick_fontsize=TICK_LABEL_FONTSIZE,
):
    xvalues = np.array(xvalues)
    yvalues = np.array(yvalues)

    y_series = pd.Series(yvalues)
    moving_avg = y_series.rolling(window=ma_window, center=True).mean()

    ax.plot(xvalues, yvalues, alpha=0.25, color=color, zorder=z, label=None)
    ax.plot(xvalues, moving_avg, alpha=1.0, color=ma_color, zorder=ma_z, label=label)

    if show_title and title:
        ax.set_title(title, fontsize=title_fontsize, pad=8)

    ax.tick_params(axis='both', labelsize=tick_fontsize)
    if not show_x_ticklabels:
        ax.tick_params(axis='x', labelbottom=False)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(width=1.5, direction='out')


def style_axes(ax, title=None):
    if title is not None:
        ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=8)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(width=1.5, direction='out')


def hide_y_axis(ax):
    ax.set_ylabel("")
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.spines['left'].set_visible(False)


def set_time_axis(ax):
    ax.set_xlim(0, TOTAL_TIME_US)
    ax.set_xticks([0, 5, 10, 15, 20])


def add_row_panel_label(ax, label):
    ax.annotate(
        label,
        xy=(-0.24, 1.10),
        xycoords='axes fraction',
        fontsize=PANEL_LABEL_FONTSIZE,
        fontweight='bold',
        ha='left',
        va='top',
        annotation_clip=False
    )


def plot_apl_panel(
    ax,
    apl_df,
    name,
    color,
    ma_window,
    show_legend=False,
    show_title=True,
    show_x_ticklabels=True,
):
    time = time_axis_from_length(len(apl_df))

    if "Total" in apl_df.columns:
        line_plot(
            time,
            apl_df["Total"].to_numpy(),
            ax,
            title=name,
            color="black",
            z=1,
            label="Total",
            ma_window=ma_window,
            ma_color="black",
            ma_z=2,
            show_title=show_title,
            show_x_ticklabels=show_x_ticklabels,
        )
    else:
        if show_title:
            ax.set_title(name, fontsize=TITLE_FONTSIZE, pad=8)

    for lip in lipid_columns:
        if lip == "Total":
            continue
        if lip not in apl_df.columns:
            continue

        arr = apl_df[lip].to_numpy()

        line_plot(
            time,
            arr,
            ax,
            title="" if "Total" in apl_df.columns else name,
            color=lipid_colors.get(lip, "gray"),
            z=1,
            label=lip,
            ma_window=ma_window,
            ma_color=lipid_colors.get(lip, "gray"),
            ma_z=2,
            show_title=False,
            show_x_ticklabels=show_x_ticklabels,
        )

    set_time_axis(ax)
    ax.set_ylabel("APL (Å$^2$)", fontsize=AXIS_LABEL_FONTSIZE)

    if show_legend:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.02),
            borderaxespad=0.0,
            fontsize=LEGEND_FONTSIZE,
            frameon=False
        )


def plot_rmsd_panel(ax, rmsd_dict, name, color, ma_window, show_title=True, show_x_ticklabels=True):
    time = np.linspace(0,20,len(rmsd_dict[name]["rmsd_A"]))
    vals = rmsd_dict[name]["rmsd_A"]
    line_plot(
        time, vals, ax, title=name, color=color, z=1,
        label="", ma_window=ma_window, ma_color=color, ma_z=2,
        show_title=show_title,
        show_x_ticklabels=show_x_ticklabels,
    )
    set_time_axis(ax)
    ax.set_ylabel("RMSD (Å)", fontsize=AXIS_LABEL_FONTSIZE)


def plot_pp_panel(ax, pp_data, name, color, ma_window, show_title=True, show_x_ticklabels=True):
    vals = np.asarray(pp_data[name])
    time = time_axis_from_length(len(vals))
    line_plot(
        time, vals, ax, title=name, color=color, z=1,
        label="", ma_window=ma_window, ma_color=color, ma_z=2,
        show_title=show_title,
        show_x_ticklabels=show_x_ticklabels,
    )
    set_time_axis(ax)
    ax.set_ylabel("Bilayer thickness (Å)", fontsize=AXIS_LABEL_FONTSIZE)


def plot_anchor_panel(ax, anchor_data, name, color, ma_window, show_title=True, show_x_ticklabels=True):
    vals = np.asarray(anchor_data[name])
    time = time_axis_from_length(len(vals))
    line_plot(
        time, vals, ax, title=name, color=color, z=1,
        label="", ma_window=ma_window, ma_color=color, ma_z=2,
        show_title=show_title,
        show_x_ticklabels=show_x_ticklabels,
    )
    set_time_axis(ax)
    ax.set_ylabel(r"$\mathrm{I}_{\mathrm{D}}$", fontsize=AXIS_LABEL_FONTSIZE)


def draw_lobe_contour(ax, df, title=None, levels=10, nbins=50, cmap="cividis", min_level=0.03
                      , bw_method="scott"):
    x = df["Lobe1"].to_numpy()
    y = df["Lobe2"].to_numpy()
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    x_min, x_max = np.nanmin(x), np.nanmax(x)
    y_min, y_max = np.nanmin(y), np.nanmax(y)
    X, Y = np.meshgrid(
        np.linspace(x_min - 5, x_max + 5, nbins),
        np.linspace(y_min - 5, y_max + 5, nbins)
    )
    coords = np.vstack([X.ravel(), Y.ravel()])
    xy = np.vstack([x, y])
    kde = gaussian_kde(xy, bw_method=bw_method)
    Z = kde(coords).reshape(X.shape)
    Z = Z / np.nanmax(Z)
    Z = np.where(Z >= min_level, Z, np.nan)

    ax.set_ylim(25, 65)
    ax.set_xlim(25, 65)
    levs = np.linspace(min_level, 1.0, levels)
    ax.contourf(X, Y, Z, levels=levs, vmin=0, vmax=1, cmap=cmap)
    style_axes(ax, title=title)
    ax.set_xlabel(r"Z$_\mathrm{lb1}$ (Å)", fontsize=LOBE_AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(r"Z$_\mathrm{lb2}$ (Å)", fontsize=LOBE_AXIS_LABEL_FONTSIZE)


def polar_density_histogram_smooth(rotationset, tiltset, ax, title=None):
    rotation_angles_rad = np.deg2rad(rotationset)
    tilt_angles_rad = np.deg2rad(tiltset)

    num_theta_bins = 100
    num_radius_bins = 100
    theta_bins = np.linspace(-np.pi, np.pi, num_theta_bins + 1)
    radius_bins = np.linspace(0, np.deg2rad(180), num_radius_bins + 1)

    histogram, _, _ = np.histogram2d(rotation_angles_rad, tilt_angles_rad, bins=[theta_bins, radius_bins])
    pad_width = 20
    histogram_padded = np.pad(histogram, ((pad_width, pad_width), (0, 0)), mode='wrap')
    histogram_smoothed = gaussian_filter(histogram_padded, sigma=(3, 3), mode='wrap')
    histogram_smoothed = histogram_smoothed[pad_width:-pad_width, :]
    histogram_smoothed = histogram_smoothed.T

    theta_centers = (theta_bins[:-1] + theta_bins[1:]) / 2
    radius_centers = (radius_bins[:-1] + radius_bins[1:]) / 2
    Theta, Radius = np.meshgrid(theta_centers, radius_centers)

    histogram_smoothed /= np.nanmax(histogram_smoothed)
    histogram_smoothed = histogram_smoothed.reshape(Theta.shape)
    histogram_smoothed[:, 0] = (histogram_smoothed[:, 0] + histogram_smoothed[:, -1]) / 2
    histogram_smoothed[:, -1] = histogram_smoothed[:, 0]

    levels = np.linspace(.02, 1, 6)

    ax.set_theta_zero_location('N')
    ax.tick_params(axis='x', which='both', labelsize=POLAR_TICK_FONTSIZE, pad=4)
    ax.tick_params(axis='y', which='both', labelsize=POLAR_TICK_FONTSIZE)
    ax.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '-135°', '-90°', '-45°'])
    ax.set_yticks(np.deg2rad([0, 30, 60, 90, 120, 150, 180]))
    ax.set_yticklabels(['', '', '', '', '', '', ''])
    ax.set_theta_direction('clockwise')
    ax.grid(linestyle='-', color='black', alpha=0.5, zorder=1)
    ax.contourf(Theta, Radius, histogram_smoothed, levels=levels, cmap='cividis', zorder=0)
    if title is not None:
        ax.set_title(title, fontsize=POLAR_TITLE_FONTSIZE)


# =========================
# Utility helpers
# =========================
def needs_recompute_apl(apl_df, expected_len):
    if apl_df is None:
        return True
    if not isinstance(apl_df, pd.DataFrame):
        return True
    if "Total" not in apl_df.columns:
        return True
    if len(apl_df["Total"]) != expected_len:
        return True
    return False


def needs_recompute_series(series_like, expected_len):
    if series_like is None:
        return True

    if isinstance(series_like, pd.Series):
        return len(series_like) != expected_len

    arr = np.asarray(series_like).squeeze()

    if arr.ndim != 1:
        return True

    return len(arr) != expected_len


def build_universes(config):
    universes = {}
    for name, cfg in config.items():
        universes[name] = mda.Universe(cfg["psf"], cfg["traj"])
    return universes


def time_axis_from_length(n, total_us=TOTAL_TIME_US):
    return np.linspace(0, total_us, n)


def last_half(arr):
    arr = np.asarray(arr)
    return arr[len(arr) // 2:]


def ensure_parent(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)


def load_pickle(path: Path):
    if path.exists():
        with open(path, "rb") as fh:
            return pickle.load(fh)
    return {}


def save_pickle(obj, path: Path):
    ensure_parent(path)
    with open(path, "wb") as fh:
        pickle.dump(obj, fh)


def find_alias_key(data, canonical_name):
    for alias in NAME_ALIASES[canonical_name]:
        if alias in data:
            return alias
    return None


def remap_top_level_keys(data):
    if not isinstance(data, dict):
        return data, False

    out = dict(data)
    changed = False
    for canonical in RHOA_SYSTEMS:
        alias = find_alias_key(out, canonical)
        if alias is not None and alias != canonical:
            out[canonical] = out[alias]
            changed = True
    return out, changed


def expected_samples(n_frames, stride):
    if stride <= 0:
        raise ValueError("stride must be a positive integer")
    return (n_frames + stride - 1) // stride


def shared_ylim_series(data_dict, key=None, pad_frac=0.05):
    ymin = np.inf
    ymax = -np.inf

    for name, v in data_dict.items():
        if isinstance(v, dict):
            if key is None:
                raise ValueError(f"{name} contains a dict; provide key=...")
            arr = np.asarray(v[key]).squeeze()
        else:
            arr = np.asarray(v).squeeze()

        ymin = min(ymin, np.nanmin(arr))
        ymax = max(ymax, np.nanmax(arr))

    pad = pad_frac * (ymax - ymin if ymax > ymin else 1.0)
    return ymin - pad, ymax + pad


def shared_ylim_apl(apl_data, columns=None, pad_frac=0.05):
    ymin = np.inf
    ymax = -np.inf

    for name, df in apl_data.items():
        cols = columns if columns is not None else df.columns
        for col in cols:
            if col not in df.columns:
                continue
            arr = np.asarray(df[col]).squeeze()
            arr = arr[np.isfinite(arr)]
            if arr.size == 0:
                continue
            ymin = min(ymin, np.nanmin(arr))
            ymax = max(ymax, np.nanmax(arr))

    if not np.isfinite(ymin) or not np.isfinite(ymax):
        raise ValueError("Could not determine APL y-limits from data.")

    pad = pad_frac * (ymax - ymin if ymax > ymin else 1.0)
    return ymin - pad, ymax + pad


# =========================
# Calculations
# =========================
def compute_apl_one(u, name, cfg):
    if lipyphilic is None:
        raise ImportError("lipyphilic is required to compute missing APL data.")

    step = cfg["stride_apl"]
    selection = "(name P) or (resname CHL1 and name O3)"

    lipid_ag = u.select_atoms(selection)
    if len(lipid_ag) == 0:
        raise ValueError(f"No lipid headgroup atoms found for selection: {selection}")

    residue_ids = lipid_ag.residues.resids
    resid_to_row = {resid: idx for idx, resid in enumerate(residue_ids)}

    u.trajectory[0]

    memcom_sel = u.select_atoms("resname POPC POPS POPE PLA18 PSM and name P")
    if len(memcom_sel) == 0:
        raise ValueError("No membrane phosphate atoms found for membrane COM calculation.")

    memcom_z = memcom_sel.center_of_mass()[2]

    headgroup_z = lipid_ag.positions[:, 2]
    leaflet_static = np.where(headgroup_z > memcom_z, 1, -1).astype(int)

    frame_indices = range(0, len(u.trajectory), step)
    n_frames = len(frame_indices)

    leaflets = np.repeat(leaflet_static[:, None], n_frames, axis=1)

    areas = lipyphilic.analysis.area_per_lipid.AreaPerLipid(
        universe=u,
        lipid_sel=selection,
        leaflets=leaflets
    )
    areas.run(step=step, verbose=True)

    apl_df = pd.DataFrame()
    apl_df["Total"] = areas.areas.mean(axis=0)

    for lip in APL_LIPIDS:
        if lip == "CHL1":
            ag_lip = u.select_atoms("resname CHL1 and name O3")
        else:
            ag_lip = u.select_atoms(f"resname {lip} and name P")

        lip_resids = ag_lip.residues.resids
        if lip_resids.size == 0:
            continue

        row_idxs = [resid_to_row[r] for r in lip_resids if r in resid_to_row]
        if not row_idxs:
            continue

        apl_df[lip] = areas.areas[row_idxs, :].mean(axis=0)

    popc_ag = u.select_atoms("resname POPC and name P")
    popc_resids = popc_ag.residues.resids
    popc_z = popc_ag.positions[:, 2]

    popcup_ids = popc_resids[popc_z > memcom_z]
    popclow_ids = popc_resids[popc_z < memcom_z]

    upper_rows = [resid_to_row[r] for r in popcup_ids if r in resid_to_row]
    lower_rows = [resid_to_row[r] for r in popclow_ids if r in resid_to_row]

    if upper_rows:
        apl_df["Upper POPC"] = areas.areas[upper_rows, :].mean(axis=0)
    if lower_rows:
        apl_df["Lower POPC"] = areas.areas[lower_rows, :].mean(axis=0)

    return apl_df


def ensure_apl_data(universes):
    data = load_pickle(APL_PICKLE)
    data, changed = remap_top_level_keys(data)

    for name, u in universes.items():
        alias = find_alias_key(data, name)
        expected_len = expected_samples(len(u.trajectory), RHOA_SYSTEMS[name]["stride_apl"])
        if alias is None or needs_recompute_apl(data[alias], expected_len):
            reason = "missing entry" if alias is None else "missing/incorrect Total column or stale length"
            print(f"Computing APL for {name} ({reason})")
            data[name] = compute_apl_one(u, name, RHOA_SYSTEMS[name])
            changed = True

    if changed:
        save_pickle(data, APL_PICKLE)

    return {name: data[find_alias_key(data, name)] for name in RHOA_SYSTEMS}


def aligned_roi_rmsd(u, roi_sel, align_sel="protein and name CA and not resid 30:43 56:76 131:141 178:191", ref_frame=0, stride=1):
    r = RMSD(
        u, u,
        select=align_sel,
        groupselections=[roi_sel],
        ref_frame=ref_frame
    ).run(step=stride, verbose=True)

    arr = r.results.rmsd
    time_ps = arr[:, 1]
    roi_rmsd_A = arr[:, 2]
    return time_ps / 1e6, roi_rmsd_A


def ensure_rmsd_data(universes):
    raw = load_pickle(RMSD_PICKLE)
    if "rmsd" in raw:
        data = raw["rmsd"]
    else:
        data = raw

    data, changed = remap_top_level_keys(data)

    for name, u in universes.items():
        alias = find_alias_key(data, name)
        expected_len = expected_samples(len(u.trajectory), RHOA_SYSTEMS[name]["stride_rmsd"])
        bad = alias is None
        if not bad:
            entry = data[alias]
            bad = (
                not isinstance(entry, dict)
                or "rmsd_A" not in entry
                or "time_us" not in entry
                or len(entry["rmsd_A"]) != expected_len
                or len(entry["time_us"]) != expected_len
            )

        if bad:
            print(f"Computing RMSD for {name}")
            t_us, rmsd_A = aligned_roi_rmsd(
                u,
                RHOA_SYSTEMS[name]["anchor_rmsd_sel"],
                align_sel="protein and name CA and not resid 25:43 58:78 120:140 178:191",
                ref_frame=1,
                stride=RHOA_SYSTEMS[name]["stride_rmsd"],
            )
            data[name] = {"time_us": t_us, "rmsd_A": rmsd_A}
            changed = True

    if changed:
        if "rmsd" in raw:
            raw["rmsd"] = data
            save_pickle(raw, RMSD_PICKLE)
        else:
            save_pickle({"rmsd": data}, RMSD_PICKLE)

    return {name: data[find_alias_key(data, name)] for name in RHOA_SYSTEMS}


def compute_pp_distance_one(u, cfg):
    zdist = []
    for _ts in tqdm(u.trajectory[::cfg["stride_pp"]], desc="P-P distance"):
        ag = u.select_atoms(MEMBRANE_P_SEL)
        mid_z = ag.center_of_mass()[2]
        top_z = u.select_atoms(f"({MEMBRANE_P_SEL} and name P) and prop z > {mid_z}").center_of_mass()[2]
        bot_z = u.select_atoms(f"({MEMBRANE_P_SEL} and name P) and prop z < {mid_z}").center_of_mass()[2]
        zdist.append(top_z - bot_z)
    return np.array(zdist)


def ensure_pp_data(universes):
    data = load_pickle(PP_PICKLE)
    data, changed = remap_top_level_keys(data)

    for name, u in universes.items():
        alias = find_alias_key(data, name)
        expected_len = expected_samples(len(u.trajectory), RHOA_SYSTEMS[name]["stride_pp"])
        if alias is None or needs_recompute_series(data[alias], expected_len):
            print(f"Computing P-P distance for {name}")
            data[name] = compute_pp_distance_one(u, RHOA_SYSTEMS[name])
            changed = True

    if changed:
        save_pickle(data, PP_PICKLE)

    return {name: data[find_alias_key(data, name)] for name in RHOA_SYSTEMS}


def compute_anchor_depth_one(u, cfg):
    depth = []
    for _ts in tqdm(u.trajectory[::cfg["stride_anchor"]], desc="Anchor insertion depth"):
        anchor_com = u.select_atoms(cfg["anchor_sel"]).center_of_mass()[2]
        mem_com = u.select_atoms(MEMBRANE_NONH_SEL).center_of_mass()[2]
        depth.append(mem_com - anchor_com)
    return np.array(depth)


def ensure_anchor_data(universes):
    data = load_pickle(ANCHOR_PICKLE)
    data, changed = remap_top_level_keys(data)

    for name, u in universes.items():
        alias = find_alias_key(data, name)
        expected_len = expected_samples(len(u.trajectory), RHOA_SYSTEMS[name]["stride_anchor"])
        if alias is None or needs_recompute_series(data[alias], expected_len):
            print(f"Computing anchor insertion depth for {name}")
            data[name] = compute_anchor_depth_one(u, RHOA_SYSTEMS[name])
            changed = True

    if changed:
        save_pickle(data, ANCHOR_PICKLE)

    return {name: data[find_alias_key(data, name)] for name in RHOA_SYSTEMS}


def compute_lobe_df_one(u, cfg):
    d1, d2 = [], []
    for _ts in tqdm(u.trajectory[::cfg["stride_lobe"]], desc="Lobe z COM"):
        lobe1 = u.select_atoms(cfg["lobe1_sel"])
        lobe2 = u.select_atoms(cfg["lobe2_sel"])
        memcen = u.select_atoms(MEMBRANE_COM_SEL)
        mem_z = memcen.center_of_mass()[2]
        d1.append(-1 * (lobe1.center_of_mass()[2] - mem_z))
        d2.append(-1 * (lobe2.center_of_mass()[2] - mem_z))
    df = pd.DataFrame({"Lobe1": np.array(d1), "Lobe2": np.array(d2)})
    half = len(df) // 2
    return df.iloc[half:].reset_index(drop=True)


def ensure_lobe_data(universes):
    data = load_pickle(LOBE_PICKLE)
    data, changed = remap_top_level_keys(data)

    for name, u in universes.items():
        alias = find_alias_key(data, name)
        if alias is None:
            print(f"Computing missing lobe COM data for {name}")
            data[name] = compute_lobe_df_one(u, RHOA_SYSTEMS[name])
            changed = True
        else:
            if alias != name:
                data[name] = data[alias]
                changed = True

    if changed:
        save_pickle(data, LOBE_PICKLE)

    return {name: data[find_alias_key(data, name)] for name in RHOA_SYSTEMS}


def load_polar_data():
    data = {}
    for name, cfg in RHOA_SYSTEMS.items():
        rot = np.loadtxt(cfg["rotation_txt"])
        tilt = np.loadtxt(cfg["tilt_txt"])
        data[name] = (last_half(rot), last_half(tilt))
    return data


# =========================
# Figure writers
# =========================
def save_validation_row_figure(row_key, row_title, ylabel, data_map, ma_window):
    names = list(RHOA_SYSTEMS.keys())
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6), constrained_layout=True)

    for col_idx, (ax, name) in enumerate(zip(axes, names)):
        color = RHOA_SYSTEMS[name]["color"]

        if row_key == "apl":
            plot_apl_panel(
                ax,
                data_map[name],
                name,
                color,
                ma_window,
                show_legend=(col_idx == len(names) - 1),
                show_title=True,
                show_x_ticklabels=True,
            )
        elif row_key == "rmsd":
            plot_rmsd_panel(
                ax, data_map, name, color, ma_window,
                show_title=True, show_x_ticklabels=True
            )
        elif row_key == "pp":
            plot_pp_panel(
                ax, data_map, name, color, ma_window,
                show_title=True, show_x_ticklabels=True
            )
        elif row_key == "anchor":
            plot_anchor_panel(
                ax, data_map, name, color, ma_window,
                show_title=True, show_x_ticklabels=True
            )

        ax.set_xlabel("Time (μs)", fontsize=AXIS_LABEL_FONTSIZE)

        if col_idx != 0:
            hide_y_axis(ax)

    axes[0].set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    add_row_panel_label(axes[0], row_title[0] if row_title else "")

    out = FIG_OUTDIR / f"validation_{row_key}.{INDIVIDUAL_VALIDATION_FMT}"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def write_extra_figures(lobe_data, polar_data):
    FIG_OUTDIR.mkdir(parents=True, exist_ok=True)

    for name in RHOA_SYSTEMS:
        fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.5))
        draw_lobe_contour(ax, lobe_data[name], title=name)
        fig.tight_layout()
        fig.savefig(FIG_OUTDIR / f"{name}_lobe_contour.{EXTRA_FMT}", bbox_inches="tight")
        plt.close(fig)

    for name in RHOA_SYSTEMS:
        fig, ax = plt.subplots(1, 1, figsize=(4.0, 3.5), subplot_kw={'projection': 'polar'})
        rot, tilt = polar_data[name]
        polar_density_histogram_smooth(rot, tilt, ax, title=name)
        fig.tight_layout()
        fig.savefig(FIG_OUTDIR / f"{name}_tilt_rotation_polar.{EXTRA_FMT}", bbox_inches="tight")
        plt.close(fig)


def make_main_validation_figure(apl_data, rmsd_data, pp_data, anchor_data, ma_window=MAIN_MA_WINDOW):
    names = list(RHOA_SYSTEMS.keys())

    fig, axes = plt.subplots(
        len(VALIDATION_ROWS),
        len(names),
        figsize=(MAIN_FIG_WIDTH, MAIN_FIG_HEIGHT),
        constrained_layout=False
    )

    row_data_lookup = {
        "apl": apl_data,
        "rmsd": rmsd_data,
        "pp": pp_data,
        "anchor": anchor_data,
    }

    apl_ylim = shared_ylim_apl(apl_data)
    rmsd_ylim = shared_ylim_series(rmsd_data, key='rmsd_A')
    pp_ylim = shared_ylim_series(pp_data)
    anchor_ylim = shared_ylim_series(anchor_data)

    n_rows = len(VALIDATION_ROWS)
    n_cols = len(names)

    for row_idx, (row_key, row_title, ylabel) in enumerate(VALIDATION_ROWS):
        for col_idx, name in enumerate(names):
            color = RHOA_SYSTEMS[name]["color"]
            ax = axes[row_idx, col_idx]

            show_title = (row_idx == 0)
            show_x_ticklabels = True
            show_apl_legend = (row_key == "apl" and col_idx == n_cols - 1)

            if row_key == "apl":
                plot_apl_panel(
                    ax,
                    apl_data[name],
                    name,
                    color,
                    ma_window,
                    show_legend=show_apl_legend,
                    show_title=show_title,
                    show_x_ticklabels=show_x_ticklabels,
                )
                ax.set_ylim(apl_ylim)

            elif row_key == "rmsd":
                plot_rmsd_panel(
                    ax,
                    rmsd_data,
                    name,
                    color,
                    ma_window,
                    show_title=show_title,
                    show_x_ticklabels=show_x_ticklabels,
                )
                ax.set_ylim(rmsd_ylim)

            elif row_key == "pp":
                plot_pp_panel(
                    ax,
                    pp_data,
                    name,
                    color,
                    ma_window,
                    show_title=show_title,
                    show_x_ticklabels=show_x_ticklabels,
                )
                ax.set_ylim(pp_ylim)

            elif row_key == "anchor":
                plot_anchor_panel(
                    ax,
                    anchor_data,
                    name,
                    color,
                    ma_window,
                    show_title=show_title,
                    show_x_ticklabels=show_x_ticklabels,
                )
                ax.set_ylim(anchor_ylim)

            ax.set_xlabel("Time (μs)", fontsize=AXIS_LABEL_FONTSIZE)

            if col_idx != 0:
                hide_y_axis(ax)


        axes[row_idx, 0].set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
        add_row_panel_label(axes[row_idx, 0], ROW_PANEL_LABELS[row_idx])

    fig.subplots_adjust(
        left=0.11,
        right=0.86,
        top=0.96,
        bottom=0.07,
        hspace=ROW_HSPACE,
        wspace=COL_WSPACE
    )

    fig.savefig(MAIN_FIG, bbox_inches="tight")
    plt.close(fig)

    for row_key, row_title, ylabel in VALIDATION_ROWS:
        save_validation_row_figure(row_key, row_title, ylabel, row_data_lookup[row_key], ma_window)


# =========================
# Main
# =========================
def main():
    FIG_OUTDIR.mkdir(parents=True, exist_ok=True)

    universes = build_universes(RHOA_SYSTEMS)

    apl_data = ensure_apl_data(universes)
    rmsd_data = ensure_rmsd_data(universes)
    pp_data = ensure_pp_data(universes)
    anchor_data = ensure_anchor_data(universes)
    lobe_data = ensure_lobe_data(universes)
    polar_data = load_polar_data()

    make_main_validation_figure(apl_data, rmsd_data, pp_data, anchor_data, ma_window=MAIN_MA_WINDOW)
    write_extra_figures(lobe_data, polar_data)

    print(f"Saved main validation figure -> {MAIN_FIG}")
    print(f"Saved individual validation figures -> {FIG_OUTDIR}")
    print(f"Saved lobe/polar figures -> {FIG_OUTDIR}")


if __name__ == "__main__":
    main()