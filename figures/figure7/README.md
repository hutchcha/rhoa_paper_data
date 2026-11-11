# Figure 7 — Dynamical-coupling network and Switch I effector competence

Compares **RhoA-GTP-2** (wild type) and **RhoA-GTP-3A**. All code and small inputs are in
**`../../analysis/cobinding_and_network_pipeline/`** (shared with Figure 6); final
rendered panels are in `rendered_panels/` here.

| Panel | Content | Script (in the pipeline) |
|---|---|---|
| A | Dynamical communities drawn on the structure (GTP-2 top, GTP-3A bottom) | `render_community_network.py` (PyMOL) |
| B | D76/D78 node betweenness centrality | `fig_dynetan_pathlength.py` |
| C | N-terminus→Switch I optimal-path length | `fig_dynetan_pathlength.py` |
| D | Switch I Cα RMSD to effector crystals (ROCK1 1S1C, PKN 1CXZ, anillin 4XOI), GTP-2 vs. GTP-3A | `fig7_panelA_combined.py` |

Panel D reads `transition_states/switches_fullres/effector_rmsd_fullres.npz`; panels B/C
read `transition_states/dynetan/`; panel A reads `transition_states/dynetan/{sim}_comm.pdb`
— all bundled in the pipeline folder. Effector reference crystals are in
`../../analysis/effector_refs/`. The dynetan network and effector RMSD are computed from
the deposited 5 µs trajectory excerpts (see `../../DATA_AVAILABILITY.md`).
