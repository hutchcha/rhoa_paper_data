# Figure 2 — RhoA G-domain dynamics on a bilayer surface

**Panels A and B (scripted here):** time evolution of the z-distance between the bilayer
COM and lobe 1 (Z_lb1, residues 1–86) and lobe 2 (Z_lb2, residues 87–177), and the
normalized 2D density of Z_lb1 vs. Z_lb2, for RhoA-GDP and RhoA-GTP-1.

```bash
python fig2_lobe_dynamics.py       # -> fig2_lobe_dynamics.{png,svg,pdf}
```
Input (bundled in `data/`): `lobe_data.pkl` (per-system lobe-COM z-distances; produced
by `analysis/orientation_and_lobe/`). Colors: Lobe 1 = orange `#F58518`, Lobe 2 = blue
`#4C78A8` (the lab lobe convention used throughout, incl. SI Fig. S2).

Representative membrane-proximal snapshots (right side of the figure) are PyMOL renders.
