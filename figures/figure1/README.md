# Figure 1 — Simulation setup and global dynamics of GDP-/GTP-RhoA

**Panel B (scripted here):** time evolution of catalytic-domain (residues 1–180) Cα RMSD
excluding the switch regions, and geranylgeranyl insertion depth, for RhoA-GDP and
RhoA-GTP-1.

```bash
python fig1_rmsd_insertion.py      # -> fig1_rmsd_insertion.{png,svg,pdf}
```
Inputs (bundled in `data/`): `rmsdrhoa.pkl` (Cα RMSD per system),
`anchor_insert_depth_acyl.pkl` (prenyl-anchor acyl-chain insertion depth; produced by
`analysis/insertion_depth/measure_anchor_insertion_acyl.py`).

**Panels A and C** are PyMOL structural renders (system overview; representative
snapshots) made from the initial system + trajectory — not script-generated here.
