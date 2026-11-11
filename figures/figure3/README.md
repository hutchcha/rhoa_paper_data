# Figure 3 — Dynamics and G-domain reorientation during RhoA-GTP-2

**Panels A–D (scripted here):** catalytic-domain Cα RMSD (A), geranylgeranyl insertion
depth (B), Z_lb1 / Z_lb2 lobe-COM z-distance time series (C), and the normalized 2D
density of Z_lb1 vs. Z_lb2 (D) for RhoA-GTP-2, which reorients into the stable
configuration at ~5 µs.

```bash
python fig3_gtp2_dynamics.py       # -> fig3_gtp2_dynamics.{png,svg,pdf}
```
Inputs (bundled in `data/`): `gtp2_core_rmsd_ref0.pkl` (core Cα RMSD),
`anchor_insert_depth_acyl.pkl` (prenyl insertion depth), `lobe_data.pkl` (lobe-COM
z-distances; `RhoA-Resim` key = RhoA-GTP-2).

**Panel E** (starting vs. final configuration) is a PyMOL render.
