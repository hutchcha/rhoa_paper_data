#!/usr/bin/env python3
"""
parallel_rdfs.py

Compute block-averaged RDFs in parallel for multiple simulations.

Default use case here:
- Lipid phosphate atoms:
    POPS, POPC, POPE, PLA18  -> inner leaflet only, name P
- Basic sidechain terminal atoms:
    ARG5 -> CZ
    LYS6 -> NZ
    LYS7 -> NZ

Outputs a long-format CSV with columns:
    sim_name
    lipid_name
    residue_name
    ref_selection
    lipid_selection
    r_A
    rdf_mean
    rdf_sem
    n_blocks
    n_ref_atoms
    n_lip_atoms

This version:
- skips jobs where a lipid is absent in a given simulation
- skips jobs where the residue atom selection is absent
- uses rdf.results.bins and rdf.results.rdf (non-deprecated)
- uses only inner leaflet lipid phosphate atoms, defined as P atoms below
  the membrane COM z at the first analyzed frame
"""

from __future__ import annotations

import multiprocessing as mp
import traceback

import numpy as np
import pandas as pd

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
# =============================================================================
# USER INPUTS
# =============================================================================

SIMULATIONS = [
    {
        "name": "GDP",
        "topology": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/sys.psf",
        "trajectory": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/cenwrapprod.dcd",
    },
    {
        "name": "GTP",
        "topology": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/sys.psf",
        "trajectory": "/home/chutchins@uthouston.edu/projects/rhoa/rhoagtpanalysis/20us/cenwrapprod.dcd",
    },
    {
        "name": "Resim",
        "topology": "/mnt/data_SSD/rhoa_redotraj/nowat.psf",
        "trajectory": "/mnt/data_SSD/rhoa_redotraj/cenwrapprod.dcd",
    },
]

# Residue/reference selections to measure against lipid phosphates
REFERENCE_GROUPS = [
    {
        "residue_name": "Arg5",
        "selection": "resid 5 and resname ARG and name NH2",
    },
    {
        "residue_name": "Lys6",
        "selection": "resid 6 and resname LYS and name NZ",
    },
    {
        "residue_name": "Lys7",
        "selection": "resid 7 and resname LYS and name NZ",
    },
]

# Lipids to measure
LIPID_GROUPS = [
    {
        "lipid_name": "POPS",
        "selection": "resname POPS and name P",
    },
    {
        "lipid_name": "POPC",
        "selection": "resname POPC and name P",
    },
    {
        "lipid_name": "POPE",
        "selection": "resname POPE and name P",
    },
    {
        "lipid_name": "PLA18",
        "selection": "resname PLA18 and name P",
    },
]

# Fraction of trajectory to skip before analysis
FRACTION_TO_SKIP = 0.25

# Define the membrane phosphate atoms used to compute membrane COM z
MEMBRANE_PHOSPHATE_SELECTION = "resname POPS POPC POPE PLA18 and name P"

# Inner leaflet definition:
# inner leaflet = phosphate atoms with z < membrane_com_z
USE_INNER_LEAFLET_ONLY = True

# RDF settings
RDF_RANGE = (0.0, 20.0)   # Angstrom
NBINS = 200

# Block averaging
# If None or <= 0: use whole analyzed region as one block
# Otherwise, use this many frames per block
BLOCK_SIZE = 5000

# Parallelism
NPROC = 20

# Output
OUT_CSV = "parallel_rdfs_long.csv"
SKIPPED_CSV = "parallel_rdfs_skipped.csv"


# =============================================================================
# CORE
# =============================================================================

def _make_blocks(n_frames: int, block_size: int | None) -> list[tuple[int, int]]:
    """Return list of (start, stop) frame ranges after skipping an initial fraction."""
    if n_frames <= 0:
        return []

    analysis_start = int(n_frames * FRACTION_TO_SKIP)

    if block_size is None or block_size <= 0 or block_size >= (n_frames - analysis_start):
        return [(analysis_start, n_frames)]

    blocks: list[tuple[int, int]] = []
    start = analysis_start
    while start < n_frames:
        stop = min(start + block_size, n_frames)
        if stop > start:
            blocks.append((start, stop))
        start = stop
    return blocks


def _run_single_rdf_job(job: dict) -> dict | None:
    """
    Worker for one (simulation, residue, lipid) RDF job.

    Returns:
        dict with RDF results, or
        dict with status='skipped', or
        dict with status='failed'
    """
    sim_name = job["sim_name"]
    topology = job["topology"]
    trajectory = job["trajectory"]
    residue_name = job["residue_name"]
    lipid_name = job["lipid_name"]
    ref_selection = job["ref_selection"]
    lipid_selection = job["lipid_selection"]
    rdf_range = job["rdf_range"]
    nbins = job["nbins"]
    block_size = job["block_size"]

    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis.rdf import InterRDF

        u = mda.Universe(topology, trajectory)

        n_frames = len(u.trajectory)
        if n_frames == 0:
            return {
                "status": "skipped",
                "sim_name": sim_name,
                "residue_name": residue_name,
                "lipid_name": lipid_name,
                "ref_selection": ref_selection,
                "lipid_selection": lipid_selection,
                "reason": "trajectory contains 0 frames",
            }

        blocks = _make_blocks(n_frames, block_size)
        if len(blocks) == 0:
            return {
                "status": "skipped",
                "sim_name": sim_name,
                "residue_name": residue_name,
                "lipid_name": lipid_name,
                "ref_selection": ref_selection,
                "lipid_selection": lipid_selection,
                "reason": "no valid blocks generated",
            }

        # Use the first analyzed frame to define the membrane center and inner leaflet
        analysis_start = blocks[0][0]
        u.trajectory[analysis_start]

        ref_ag = u.select_atoms(ref_selection)
        if ref_ag.n_atoms == 0:
            return {
                "status": "skipped",
                "sim_name": sim_name,
                "residue_name": residue_name,
                "lipid_name": lipid_name,
                "ref_selection": ref_selection,
                "lipid_selection": lipid_selection,
                "reason": f"reference selection matched 0 atoms: {ref_selection}",
            }

        raw_lip_ag = u.select_atoms(lipid_selection)
        if raw_lip_ag.n_atoms == 0:
            return {
                "status": "skipped",
                "sim_name": sim_name,
                "residue_name": residue_name,
                "lipid_name": lipid_name,
                "ref_selection": ref_selection,
                "lipid_selection": lipid_selection,
                "reason": f"lipid selection matched 0 atoms: {lipid_selection}",
            }

        if USE_INNER_LEAFLET_ONLY:
            membrane_p = u.select_atoms(MEMBRANE_PHOSPHATE_SELECTION)
            if membrane_p.n_atoms == 0:
                return {
                    "status": "failed",
                    "sim_name": sim_name,
                    "residue_name": residue_name,
                    "lipid_name": lipid_name,
                    "ref_selection": ref_selection,
                    "lipid_selection": lipid_selection,
                    "reason": f"membrane phosphate selection matched 0 atoms: {MEMBRANE_PHOSPHATE_SELECTION}",
                }

            membrane_com_z = float(membrane_p.center_of_mass()[2])

            # Select only this lipid species' phosphate atoms below membrane COM z
            lip_ag = raw_lip_ag[raw_lip_ag.positions[:, 2] < membrane_com_z]

            leaflet_selection_desc = (
                f"({lipid_selection}) AND z < membrane_com_z({membrane_com_z:.3f} Å)"
            )
        else:
            lip_ag = raw_lip_ag
            leaflet_selection_desc = lipid_selection

        if lip_ag.n_atoms == 0:
            return {
                "status": "skipped",
                "sim_name": sim_name,
                "residue_name": residue_name,
                "lipid_name": lipid_name,
                "ref_selection": ref_selection,
                "lipid_selection": leaflet_selection_desc,
                "reason": (
                    f"no inner leaflet atoms found for lipid selection at analysis start: "
                    f"{lipid_selection}"
                ),
            }

        block_rdfs = []
        bins = None

        for start, stop in blocks:
            rdf = InterRDF(
                ref_ag,
                lip_ag,
                range=rdf_range,
                nbins=nbins,
                verbose=True,
            )
            rdf.run(start=start, stop=stop)

            this_bins = np.asarray(rdf.results.bins, dtype=float)
            this_rdf = np.asarray(rdf.results.rdf, dtype=float)

            if bins is None:
                bins = this_bins
            else:
                if len(this_bins) != len(bins) or not np.allclose(this_bins, bins):
                    return {
                        "status": "failed",
                        "sim_name": sim_name,
                        "residue_name": residue_name,
                        "lipid_name": lipid_name,
                        "ref_selection": ref_selection,
                        "lipid_selection": leaflet_selection_desc,
                        "reason": "RDF bins inconsistent across blocks",
                    }

            block_rdfs.append(this_rdf)

        if len(block_rdfs) == 0:
            return {
                "status": "skipped",
                "sim_name": sim_name,
                "residue_name": residue_name,
                "lipid_name": lipid_name,
                "ref_selection": ref_selection,
                "lipid_selection": leaflet_selection_desc,
                "reason": "no RDF blocks produced",
            }

        block_rdfs = np.asarray(block_rdfs, dtype=float)
        rdf_mean = np.mean(block_rdfs, axis=0)

        if block_rdfs.shape[0] > 1:
            rdf_sem = np.std(block_rdfs, axis=0, ddof=1) / np.sqrt(block_rdfs.shape[0])
        else:
            rdf_sem = np.zeros_like(rdf_mean)

        return {
            "status": "ok",
            "sim_name": sim_name,
            "residue_name": residue_name,
            "lipid_name": lipid_name,
            "ref_selection": ref_selection,
            "lipid_selection": leaflet_selection_desc,
            "r_A": bins,
            "rdf_mean": rdf_mean,
            "rdf_sem": rdf_sem,
            "n_blocks": int(block_rdfs.shape[0]),
            "n_ref_atoms": int(ref_ag.n_atoms),
            "n_lip_atoms": int(lip_ag.n_atoms),
        }

    except Exception as e:
        return {
            "status": "failed",
            "sim_name": sim_name,
            "residue_name": residue_name,
            "lipid_name": lipid_name,
            "ref_selection": ref_selection,
            "lipid_selection": lipid_selection,
            "reason": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
        }


def build_jobs() -> list[dict]:
    jobs: list[dict] = []

    for sim in SIMULATIONS:
        for ref in REFERENCE_GROUPS:
            for lip in LIPID_GROUPS:
                jobs.append({
                    "sim_name": sim["name"],
                    "topology": sim["topology"],
                    "trajectory": sim["trajectory"],
                    "residue_name": ref["residue_name"],
                    "lipid_name": lip["lipid_name"],
                    "ref_selection": ref["selection"],
                    "lipid_selection": lip["selection"],
                    "rdf_range": RDF_RANGE,
                    "nbins": NBINS,
                    "block_size": BLOCK_SIZE,
                })

    return jobs


def flatten_ok_results(results: list[dict]) -> pd.DataFrame:
    rows = []

    for res in results:
        if res["status"] != "ok":
            continue

        for r, mean, sem in zip(res["r_A"], res["rdf_mean"], res["rdf_sem"]):
            rows.append({
                "sim_name": res["sim_name"],
                "lipid_name": res["lipid_name"],
                "residue_name": res["residue_name"],
                "ref_selection": res["ref_selection"],
                "lipid_selection": res["lipid_selection"],
                "r_A": float(r),
                "rdf_mean": float(mean),
                "rdf_sem": float(sem),
                "n_blocks": int(res["n_blocks"]),
                "n_ref_atoms": int(res["n_ref_atoms"]),
                "n_lip_atoms": int(res["n_lip_atoms"]),
            })

    return pd.DataFrame(rows)


def collect_nonok_results(results: list[dict]) -> pd.DataFrame:
    rows = []

    for res in results:
        if res["status"] == "ok":
            continue

        rows.append({
            "status": res.get("status", "unknown"),
            "sim_name": res.get("sim_name"),
            "lipid_name": res.get("lipid_name"),
            "residue_name": res.get("residue_name"),
            "ref_selection": res.get("ref_selection"),
            "lipid_selection": res.get("lipid_selection"),
            "reason": res.get("reason", ""),
        })

    return pd.DataFrame(rows)


def main() -> None:
    jobs = build_jobs()
    print(f"[INFO] built {len(jobs)} RDF jobs")

    if len(jobs) == 0:
        raise RuntimeError("No RDF jobs were built.")

    with mp.Pool(processes=NPROC) as pool:
        results = pool.map(_run_single_rdf_job, jobs)

    ok_count = sum(r["status"] == "ok" for r in results)
    skipped_count = sum(r["status"] == "skipped" for r in results)
    failed_count = sum(r["status"] == "failed" for r in results)

    print(f"[INFO] completed jobs: ok={ok_count}, skipped={skipped_count}, failed={failed_count}")

    for r in results:
        if r["status"] == "skipped":
            print(f"[SKIP] [{r['sim_name']}] {r['residue_name']} vs {r['lipid_name']} -> {r['reason']}")
        elif r["status"] == "failed":
            print(f"[FAIL] [{r['sim_name']}] {r['residue_name']} vs {r['lipid_name']} -> {r['reason']}")

    ok_df = flatten_ok_results(results)
    if ok_df.empty:
        raise RuntimeError("All RDF jobs were skipped or failed; no output CSV written.")

    ok_df.to_csv(OUT_CSV, index=False)
    print(f"[SAVED] {OUT_CSV}")

    nonok_df = collect_nonok_results(results)
    if not nonok_df.empty:
        nonok_df.to_csv(SKIPPED_CSV, index=False)
        print(f"[SAVED] {SKIPPED_CSV}")


if __name__ == "__main__":
    main()
