# SI Figure S1 — Bilayer structural properties

Area per lipid (APL) and bilayer thickness (P–P) across RhoA-GDP, RhoA-GTP-1, RhoA-GTP-2,
including comparison of POPC/PSM APL to pure-component scattering experiments.

**Scripted here** (the APL-vs-experiment and thickness comparison panels, C/D):
```bash
python si_membrane_lipid_figure.py    # combined A–D style figure -> .png/.svg
python si_membrane_AC.py              # the POPC-APL + bilayer-thickness cut -> .png/.svg/.pdf
```
Inputs (bundled in `data/`): `apl_per_type.csv`, `apl_aggregate.csv` (per-lipid and
aggregate APL; from `analysis/membrane_properties/apl_voronoi_compare.py` /
`APL_background.py`), `insertion_thickness.csv` (P–P thickness + insertion;
from `analysis/membrane_properties/measure_insertion_thickness.py`).

The APL/thickness **time-series** panels (A/B) are produced from
`analysis/membrane_properties/all_sims_apl.pkl` and `rhoa_ppdistances_val.pkl`
(APL via `lipyphilic` Voronoi tessellation; thickness = inter-leaflet phosphate COM
z-distance). Experimental reference values (POPC 65.6 Å², PSM 61.9 Å²) are cited in the
paper SI. See `analysis/membrane_properties/` for the APL-vs-experiment (`charmmgui_apl_overestimate.py`)
and PSM-condensation (`fig_apl_psm_check.py`) analyses.
