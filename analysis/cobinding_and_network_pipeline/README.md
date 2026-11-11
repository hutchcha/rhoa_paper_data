# Co-binding & dynamical-network pipeline (main Figures 6 and 7)

This directory is a self-contained copy of the analysis pipeline that produces
**Figure 6** (cooperative basic–acidic lipid co-binding) and **Figure 7**
(dynamical-coupling network + Switch I effector competence) for the membrane-matched
wild type (RhoA-GTP-2) and the N-terminal triple mutant (RhoA-GTP-3A), together with
the reference systems (RhoA-GDP, RhoA-GTP-1) and the R5A/K6A control (RhoA-GTP-2A).

All small inputs needed to regenerate the figures are included here (`transition_states/`,
`topologies/`, `../effector_refs/`). The **production trajectories are not** — they are
too large for GitHub (see `../../DATA_AVAILABILITY.md`). Scripts expect the per-system
5 µs trajectory excerpts `{sim}_last5us.dcd` and the H-bond-window excerpts
`{sim}_last1.95us.dcd` in this directory (`sim` ∈ `GDP, GTP-1, GTP-2, 2A, 3A`).

Run scripts with this folder as the working directory.

## Legacy system names
`GDP` = RhoA-GDP · `GTP-1` = RhoA-GTP-1 (simple membrane) · `GTP-2` = RhoA-GTP-2
(complex membrane) · `2A` = RhoA-GTP-2A · `3A` = RhoA-GTP-3A.

## Shared modules
- `stable_state_config.py` — central config: topology/trajectory paths, state frame
  windows, atom selections. Imported by most stages.
- `dynetan_community_colors.py` — fixed community color map for the network renders.

## Figure 6 — lipid co-binding
Upstream data generation (needs trajectories):
1. `run_hbonds.py` → protein↔lipid H-bond search → `hbonds_stable/{sim}_{gdlip,anclip}_raw_hbonds.npy` (raw tables; on Zenodo, see DATA_AVAILABILITY).
2. `lipid_hbond_analysis.py` → headgroup map + per-pair cooperativity → `transition_states/gains/lipid_cobinding_*.csv`.
3. `cobinding_occupancy.py` → `transition_states/gains/cobinding_{occupancy,frequency}_*.csv`.
4. `plot_region_hbonds.py` → `transition_states/region_hbonds/region_hbond_occupancy.csv`.
   (`run_all.sh` / `run_all_5us.sh` chain these.)

Figure panels (regenerate from the bundled `transition_states/` CSVs — no trajectories needed):
- `fig_transition_tiltrot.py` → Fig 6A, RhoA-GTP-3A tilt/rotation + post-transition polar density (reads `../../figures/figure6` orientation txt / `orientation_out*`; see note).
- `fig_A_transition_states.py` → Fig 6B, sidechain–lipid H-bond occupancy by region (pre/cross/post) from `transition_states/region_hbonds/region_hbond_occupancy.csv`.
- `render_cobinding_struct.py` / `render_cobinding_panelB.py` (PyMOL) → Fig 6C co-binding snapshot from `transition_states/cobinding_struct/`.
- `fig_B_amine_frequency.py` → Fig 6D, co-binding frequency by lipid type from `transition_states/gains/cobinding_frequency_amine.csv`.
- `fig_cobinding_network.py` → Fig 6E, co-binding network from `transition_states/gains/lipid_cobinding_enhancement_amine_min25.csv` + `topologies/{sim}.psf`.
- `compose_figure6.py` → assembles the composite.

## Figure 7 — dynamical network + effector competence
Upstream data generation (needs trajectories + `../effector_refs/`):
- `fig_D_effector_rmsd.py` → computes Switch I/II Cα RMSD to effector crystals (ROCK1 1S1C, PKN 1CXZ, anillin 4XOI) → **writes `transition_states/switches_fullres/effector_rmsd_fullres.npz`** (the Fig 7D / S8 / S9 input).
- `dynetan_network.py` → generalized-correlation (dynetan) network → `transition_states/dynetan/{sim}_{nodes,edges}.csv`, `{sim}_path.txt`, `{sim}_{comm,btw}.pdb`.

Figure panels (regenerate from bundled `transition_states/dynetan/` + `switches_fullres/`):
- `render_community_network.py` (PyMOL) → Fig 7A, communities on the structure.
- `fig_dynetan_pathlength.py` → Fig 7B (D76/D78 betweenness) and Fig 7C (N-terminus→Switch I optimal-path length).
- `fig7_panelA_combined.py` → Fig 7D, combined Switch I competence KDE from `switches_fullres/effector_rmsd_fullres.npz`.
- `compose_figure7.py` → assembles the composite.
- `dynetan_community_track.py`, `render_network_on_structure.py`, `render_path_bundle.py`, `render_dynetan*.py`, `extract_effector_snapshots.py` + `render_effector_snapshots.py` → supporting renders / snapshot selection.

## Notes
- `dccm_mi.py` + `network_allostery.py` (outputs in `transition_states/dccm_mi/`,
  `transition_states/network/`) are an **alternate** mutual-information network kept for
  completeness; the figures in the paper use the **dynetan** outputs above.
- `environment.yml` / `environment_dynetan.yml` capture the conda environments used
  (dynetan needs its own env).
- `README_original_pipeline.md` is the author's original working-notes README.
