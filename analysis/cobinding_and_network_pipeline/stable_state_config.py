#!/usr/bin/env python3
"""Shared config for the 'why is the stable state stable?' analyses.

Both the 3A (R5A/K6A/K7A) and the independent 2A (R5A/K6A) sims settle into the
same stable membrane-bound configuration -> it's a deep free-energy attractor.
We test two hypotheses for what holds it there once reached:
  H1  Lobe (protein-lipid) contacts pin it in place.
  H2  HVR/anchor interactions with the G-domain lock a single conformation.

Analyses (each its own script, all import this config):
  H1a  h1a_contacts.py      per-residue per-lipid CONTACT occupancy (4.0 A heavy)
  H1b  hbonds_h1b_h2b.py    per-residue per-lipid H-BOND occupancy (both lobes)
  H2a  h2a_anchor_rmsd.py   anchor backbone RMSD (own-frame superpose) KDE+ts
  H2b  hbonds_h1b_h2b.py    anchor<->G-domain H-bond occupancy per domain residue

Window: LAST HALF of every trajectory (characterise the settled state, not the
approach). Membrane note: only the complex-membrane sims (GTP-2 WT ref, 3A, 2A)
share a composition (POPS 50 / PLA18 50 / PSM 100), so per-lipid comparisons are
clean only within that trio. PSM is the opposite leaflet -> excluded from lipids.
"""
import os, re

# ---- simulations -----------------------------------------------------------
# 3A / 2A are the full 1 ns/frame extended runs (9921 frames = 10 us).
# GTP-2 (WT, same complex membrane) is only available 10 ns-strided; occupancy
# (a frame fraction) is robust to that coarser spacing.
# GTP-2 uses the ORIGINAL 20 us trajectory (19843 frames @ ~1.008 ns), NOT
# traj_strided/GTP-2_10ns.dcd -- that 10x-decimated copy exists only for PyLipID's
# memory limits, and using it here gave GTP-2 5x worse sampling than 3A/2A.
# ============================================================================
# FULL-RESOLUTION POST-STATE EXCERPTS (matched frame count, last ~1.95 us).
# EDIT DT_NS and POST_US to match your extraction (extract_post_window.py):
#   POST_US = the WINDOW_US you extracted (1.95 to match 2A's post state)
#   DT_NS   = the FULL trajectories' native per-frame time (ns) -- NOT 1.008
# so that round(POST_US*1000/DT_NS) == the excerpt frame count. This matters:
# occupancy denominators use POST_US/DT_NS, so they must equal the real #frames.
DT_NS   = 0.1     # <-- full-resolution per-frame time (ns)
POST_US = 5      # <-- 5 us extended excerpt (50000 frames); 2A dropped

# ---- PATHS: one row per sim -- edit these to wherever the files live ----------
# "top"  = topology. Defaults to the PSF bundled in topologies/ (synced with this
#          folder). Repoint to a GPU path if you'd rather not sync it.
# "traj" = the 1.95-us excerpt DCD. Default assumes it's in rhoa_full_sampled_trajs/
#          inside this bundle; point it anywhere (e.g. a scratch dir on the GPU)
#          so you don't have to move the big files around.
# Run  python verify_inputs.py  first -- it loads each PSF+trajectory pair and
# flags any atom-count mismatch or wrong frame count BEFORE the long run.
_HERE = os.path.dirname(os.path.abspath(__file__))
PATHS = {
    "2A":    {"top": f"{_HERE}/topologies/2A.psf",    "traj": f"{_HERE}/2A_last{POST_US}us.dcd"},
    "3A":    {"top": f"{_HERE}/topologies/3A.psf",    "traj": f"{_HERE}/3A_last{POST_US}us.dcd"},
    "GTP-2": {"top": f"{_HERE}/topologies/GTP-2.psf", "traj": f"{_HERE}/GTP-2_last{POST_US}us.dcd"},
    "GDP":   {"top": f"{_HERE}/topologies/GDP.psf",   "traj": f"{_HERE}/GDP_last{POST_US}us.dcd"},
    "GTP-1": {"top": f"{_HERE}/topologies/GTP-1.psf", "traj": f"{_HERE}/GTP-1_last{POST_US}us.dcd"},
}
SIMS = [{"name": n, "top": PATHS[n]["top"], "traj": PATHS[n]["traj"], "dt_ns": DT_NS}
        for n in ("2A", "3A", "GTP-2", "GDP", "GTP-1")]

# ---- lipids (PSM excluded: opposite leaflet) -------------------------------
LIPIDS = ["POPS", "POPE", "POPC", "PLA18"]
LIPID_SEL = "resname " + " ".join(LIPIDS)
LIPID_COLORS = {"POPS": "#E45756", "POPE": "#4C78A8", "POPC": "#54A24B", "PLA18": "#B279A2"}

# ---- protein / anchor selections -------------------------------------------
# protein = resid 4-189; CYSG (prenyl cys, resid 190) is a separate resname.
# IMPORTANT: guard every resid range with `protein` -- lipid molecules share the
# same resid numbers, so a bare `resid 180-189` also grabs membrane headgroups.
PROTEIN_SEL = "protein or resname CYSG"
ANCHOR_SEL  = "(protein and resid 180-189) or resname CYSG"   # HVR + prenyl Cys tip
GDOMAIN_SEL = "protein and resid 4-179"                        # catalytic / G-domain
LOBE1 = (4, 86)
LOBE2 = (87, 179)

# ---- H-bond geometric criteria (user: D-A <= 3.5 A, angle >= 120 deg) -------
HB_DA_CUTOFF = 3.5      # donor-acceptor distance, Angstrom
HB_ANGLE_CUTOFF = 120   # D-H-A angle, degrees

# ---- contact criterion (H1a) -----------------------------------------------
CONTACT_CUTOFF = 4.0    # Angstrom, heavy-atom to heavy-atom

# ---- structural regions (mutually exclusive) -------------------------------
REGIONS = {
    "N-term patch": [(5, 7)],
    "Switch I":     [(30, 40)],
    "beta1-3":      [(8, 11), (41, 48), (51, 55)],
    # extended 56-76 -> 56-78: D78 (and T77) were falling into no region at all,
    # despite being the partner we have been calling a "Switch II acidic"
    # throughout (D76/D78). beta4-6 still starts at 79, so no overlap.
    "Switch II":    [(56, 78)],
    "beta4-6":      [(79, 85), (112, 117), (155, 160)],
    "Helix 3":      [(86, 103)],
    "Helix 4":      [(138, 148)],
    "Helix 5":      [(164, 177)],
    "HVR+anchor":   [(180, 190)],
}
# order for heatmap rows (N->C)
REGION_ORDER = ["N-term patch", "beta1-3", "Switch I", "Switch II", "beta4-6",
                "Helix 3", "Helix 4", "Helix 5", "HVR+anchor"]

# ---- transition states -------------------------------------------------------
# All three sims move between the SAME two orientational states (~100 deg tilt /
# +65 deg rotation  ->  ~40 deg / -25 deg). Windows below are read off the native
# 1 ns tilt traces (rotation flips sign at the crossing and confirms them).
# NB rotation is DEGENERATE at low tilt (the azimuthal pointer collapses as
# tilt->0), so tilt is the state variable and rotation is only a cross-check.
# GTP-2 transitions at ~5 us but its saved H-bond raw starts at 10 us, so only
# its post state has bond-level data.
# Each excerpt IS the settled post state (the last POST_US of the full run), so
# every sim is post-only over its whole excerpt: post = (0, POST_US). No pre/cross.
# NB post = (0, POST_US) makes the occupancy denominator (POST_US*1000/DT_NS) equal
# the excerpt frame count -- keep DT_NS/POST_US consistent with the extraction.
STATES = {n: {"pre": None, "cross": None, "post": (0.0, POST_US)}
          for n in ("2A", "3A", "GTP-2", "GDP", "GTP-1")}
STATE_ORDER = ["pre", "cross", "post"]
STATE_COLORS = {"pre": "#4C78A8", "cross": "#F58518", "post": "#E45756"}

def state_frames(sim_name, state, dt_ns):
    """Frame index range (start, stop) for a state window, or None."""
    w = STATES.get(sim_name, {}).get(state)
    if w is None:
        return None
    return int(round(w[0] * 1000.0 / dt_ns)), int(round(w[1] * 1000.0 / dt_ns))

AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSG": "C",
          "GLN": "Q", "GLU": "E", "GLY": "G", "HSD": "H", "HIS": "H", "ILE": "I",
          "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
          "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

def res_label(resname, resid):
    """Compact residue label for plot ticks, e.g. ARG182 -> R182."""
    return f"{AA3TO1.get(str(resname), '?')}{int(resid)}"

def resolve(relpath):
    """Find a results file whether or not organize_results.py has moved it.

    Scripts write to flat dirs (rmsd_stable/, transition_states/, ...) but the
    organizer relocates csv/png outputs into results_by_sim/<sim>/<measure>/.
    Checks the original path first, then searches the organized tree by basename."""
    import glob as _glob
    if os.path.exists(relpath):
        return relpath
    hits = _glob.glob(os.path.join("results_by_sim", "*", "*", os.path.basename(relpath)))
    return hits[0] if hits else relpath

def region_of(resid):
    for name, ranges in REGIONS.items():
        if any(a <= resid <= b for a, b in ranges):
            return name
    return None

def parse_resid(residue_str):
    m = re.match(r"(\d+)", str(residue_str))
    return int(m.group(1)) if m else -1

# Analyse the LAST 2 us of every sim: the only window in which all three are in
# the SAME orientational state. 2A sits near-perpendicular (~100 deg tilt) for its
# first ~7 us and only joins the shared ~40 deg state at ~7-8 us, so any wider
# window averages two different states together for 2A. 3A settles by ~1-2 us and
# GTP-2 is stable throughout, so both tolerate a wider window -- 2A sets the limit.
# 2 us = 1984 frames/sim at ~1.008 ns.
WINDOW_US = 2.0

def window_slice(n_frames, stride, dt_ns=1.008):
    """Last-WINDOW_US window: (start, stop, step) for trajectory slicing.

    Gives every sim the same real duration AND (at 1 ns sampling) the same frame
    count: 3A/2A are 10 us so this is their whole trajectory; GTP-2 is 20 us so
    it takes the second half. Clamped to 0 if a trajectory is shorter."""
    n_window = int(round(WINDOW_US * 1000.0 / dt_ns))
    start = max(0, n_frames - n_window)
    return start, n_frames, stride

# plotting defaults shared by all figures
def apply_style():
    import matplotlib, matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["axes.linewidth"] = 1.5
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["svg.fonttype"] = "none"
