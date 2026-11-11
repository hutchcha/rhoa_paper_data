# Analysis scripts

Upstream trajectory-analysis code grouped by type. These scripts read the production MD
trajectories (see `../DATA_AVAILABILITY.md`) and write the small processed tables that the
figure scripts in `../figures/` consume; the key output tables are bundled alongside the
scripts here so the provenance of every figure input is traceable. Absolute trajectory
paths inside some scripts point to the original compute hosts — repoint them to the
deposited trajectories to re-run.

| Subfolder | Produces data for | Key scripts |
|---|---|---|
| `membrane_properties/` | SI Fig. S1 (APL, P–P thickness, APL-vs-experiment) | `APL_background.py` (lipyphilic Voronoi APL), `measure_insertion_thickness.py`, `apl_voronoi_compare.py`, `charmmgui_apl_overestimate.py` |
| `orientation_and_lobe/` | Figs. 2, 3, 4; SI Fig. S2 | `time_ev_tiltrot.py` (tilt/rotation), `switchii_contacts.py` (Switch I / Helix 4 z-dist + lobe contacts), `tiltrotationbins.ipynb` (φ/θ density), `master_rhoa_3sim_plotting.py` (original master figure generator) |
| `insertion_depth/` | Figs. 1, 3 | `insertdepth.py`, `measure_anchor_insertion_acyl.py`, `measure_gtp2_core_rmsd.py` |
| `contacts_hbonds/` | Fig. 4 | `contacts_parralel_RhoA.py` (relative contact frequency), `hbonds_gr_analysis.py` + `parallel_rdfs.py` (R5/K6/K7 H-bonds & RDFs), `RhoA_Lys7_sasa.py` |
| `cobinding_and_network_pipeline/` | Figs. 6, 7 | self-contained pipeline — see its own README |
| `effector_refs/` | Fig. 7 | effector crystal structures (ROCK1 1S1C, PKN 1CXZ, anillin 4XOI) |
