# Figure 6 — Cooperative lipid co-binding of basic and acidic residues

Compares the membrane-matched wild type **RhoA-GTP-2** and the N-terminal triple mutant
**RhoA-GTP-3A** (with RhoA-GDP / RhoA-GTP-1 references).

All code and the small processed-data inputs live in
**`../../analysis/cobinding_and_network_pipeline/`** (shared with Figure 7). Final
rendered panels are in `rendered_panels/` here.

| Panel | Content | Script (in the pipeline) |
|---|---|---|
| A | RhoA-GTP-3A tilt/rotation + post-transition polar density | `fig_transition_tiltrot.py` |
| B | Sidechain–lipid H-bond occupancy by region (pre/cross/post) | `fig_A_transition_states.py` |
| C | Co-binding snapshot (R5/K6/K7 + D78 share one POPS) | `render_cobinding_struct.py` (PyMOL) |
| D | Co-binding frequency by lipid type | `fig_B_amine_frequency.py` |
| E | Residue co-binding network (edge width = occupancy) | `fig_cobinding_network.py` |

See `../../analysis/cobinding_and_network_pipeline/README.md` for the full run order and
the trajectory inputs (deposited; see `../../DATA_AVAILABILITY.md`).
