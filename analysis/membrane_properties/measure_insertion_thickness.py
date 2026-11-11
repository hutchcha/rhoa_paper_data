#!/usr/bin/env python3
"""Is the anchor-insertion-depth difference between membranes just a thickness effect?

For each sim (final-5 us trajectory, stride 25) measure, per frame:
  z_mcom  = membrane heavy-atom COM z (bilayer center; suite convention)
  z_Plo   = lower-leaflet phosphate plane z   (the anchored, cytoplasmic leaflet)
  z_Pup   = upper-leaflet phosphate plane z
  z_anc   = geranylgeranyl anchor COM z (CYSG 190, C1-C20)   -- same selection as insertion_depth.py

Derived:
  depth_center = z_mcom - z_anc        suite "insertion depth" (higher = further below centre)
  thick_PP     = z_Pup - z_Plo         total bilayer thickness (P-P)
  half_lower   = z_mcom - z_Plo        lower monolayer thickness (centre -> lower P plane)
  pen_surface  = z_anc - z_Plo         penetration from the lower phosphate plane (surface-referenced)
  frac_pen     = pen_surface / half_lower   thickness-NORMALISED insertion (0 = at P plane, 1 = centre)

If depth_center differs between membranes but frac_pen (or pen_surface) does not, the raw difference
is a thickness/reference effect, not a real change in how deep the anchor buries.

  conda activate memb_analysis ; python measure_insertion_thickness.py
"""
import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, pandas as pd
import MDAnalysis as mda

HERE = os.path.dirname(os.path.abspath(__file__))
TRJ = "/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs"
SIMS = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
MEMBRANE = {"GDP": "simple", "GTP-1": "simple", "GTP-2": "complex", "2A": "complex", "3A": "complex"}
PREN = "resname CYSG and name " + " ".join(f"C{i}" for i in range(1, 21))     # geranylgeranyl C1-C20
MEMHEAVY = "resname POPC POPE POPS PLA18 PSM and not name H*"
STRIDE = 25

rows = []
for sim in SIMS:
    u = mda.Universe(f"{TRJ}/topologies/{sim}.psf", f"{TRJ}/{sim}_last5us.dcd")
    anc = u.select_atoms(PREN)
    memh = u.select_atoms(MEMHEAVY)
    P = u.select_atoms("resname POPC POPE POPS PLA18 PSM and name P")
    rec = {k: [] for k in ("depth_center", "thick_PP", "half_lower", "pen_surface", "frac_pen")}
    for ts in u.trajectory[::STRIDE]:
        pz = P.positions[:, 2]
        mid = np.median(pz)
        z_Plo = pz[pz < mid].mean(); z_Pup = pz[pz > mid].mean()
        z_mcom = memh.center_of_mass()[2]
        z_anc = anc.center_of_mass()[2]
        half_lower = z_mcom - z_Plo
        rec["depth_center"].append(z_mcom - z_anc)
        rec["thick_PP"].append(z_Pup - z_Plo)
        rec["half_lower"].append(half_lower)
        rec["pen_surface"].append(z_anc - z_Plo)
        rec["frac_pen"].append((z_anc - z_Plo) / half_lower)
    r = {"sim": sim, "membrane": MEMBRANE[sim], "nframes": len(rec["thick_PP"])}
    NBLOCKS = 5
    for k, v in rec.items():
        v = np.array(v)
        r[f"{k}_mean"] = v.mean(); r[f"{k}_sd"] = v.std(ddof=1)
        bmeans = np.array([b.mean() for b in np.array_split(v, NBLOCKS)])   # block-averaged SEM
        r[f"{k}_sem"] = bmeans.std(ddof=1) / np.sqrt(NBLOCKS)
    rows.append(r)
    print(f"[{sim}] depth_center={r['depth_center_mean']:.2f}  thick_PP={r['thick_PP_mean']:.2f}  "
          f"pen_surface={r['pen_surface_mean']:.2f}  frac_pen={r['frac_pen_mean']:.3f}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{HERE}/insertion_thickness.csv", index=False)
pd.set_option("display.width", 200)
cols = ["sim", "membrane", "depth_center_mean", "thick_PP_mean", "half_lower_mean",
        "pen_surface_mean", "frac_pen_mean"]
print("\n" + df[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\nwrote insertion_thickness.csv")
