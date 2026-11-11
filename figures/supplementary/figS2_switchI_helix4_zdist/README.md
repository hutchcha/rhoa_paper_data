# SI Figure S2 — Switch I / Helix 4 distance from the bilayer center

Time evolution of the COM z-distance of **Switch I** (resid 30–43) and **Helix 4**
(resid 141–151) — on opposite faces of the G-domain — from the bilayer center, for
**RhoA-GDP**, **RhoA-GTP-1**, and **RhoA-GTP-2**.

## Reproduce
```bash
python plot_swI_h4_membrane_zdist_WT.py
```
Writes `swI_h4_membrane_zdist_WT.{png,svg,pdf}` (5.0 × 6.0 in).

## Input (bundled)
`data/region_midplane_zmetrics_swI_helix4.csv` — per-frame signed z-offsets
`d_a = z_SwitchI − z_mid` and `d_b = z_Helix4 − z_mid` (Å) plus `time_us` for each
system. Distance to the membrane center is plotted as |d|. Columns extracted from the
full per-frame table produced by
`analysis/orientation_and_lobe/switchii_contacts.py` over the 20 µs trajectories.

Notes: stored time is scaled ×0.1 (→ ~19.84 µs); z_mid = lipid bilayer COM;
colors match main Fig. 2 (Switch I = Lobe 1 orange, Helix 4 = Lobe 2 blue).
