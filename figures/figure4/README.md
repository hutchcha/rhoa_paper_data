# Figure 4 — RhoA adopts a single membrane orientation

Membrane-referenced orientation (rotation φ, tilt θ), residue–lipid contacts, and
N-terminal R5/K6/K7–POPS hydrogen bonds across RhoA-GDP, RhoA-GTP-1, RhoA-GTP-2.

**Scripted here:**
- **Panel B** — φ/θ time series:
  ```bash
  python tiltrot_timeseries.py          # -> tiltrot_timeseries.{png,svg,pdf}
  ```
  Inputs (`data/`): `RhoA_{GDP,GTP,GTPredo}{tilt,rotation}angles.txt` (θ/φ per frame;
  `GTPredo` = RhoA-GTP-2). Computed by `analysis/orientation_and_lobe/time_ev_tiltrot.py`.
- **Panel E** — combined R5/K6/K7–POPS H-bond occupancy (last 15 µs, block-SEM):
  ```bash
  python fig4_pops_hbond_occupancy.py   # -> fig4_pops_hbond_occupancy.{png,svg,pdf}
  ```
  Inputs (`data/`): `{GDP,GTP,Resim}_hbonds.csv`.

**Other panels:**
- **A** — schematic/PyMOL definition of the φ/θ coordinates (see SI methods: VH5 from
  P168→D161; see `analysis/orientation_and_lobe/switchii_contacts.py` for the geometry).
- **C** — φ vs. θ 2D probability densities: `analysis/orientation_and_lobe/tiltrotationbins.ipynb`.
- **D** — residue–lipid relative contact frequency: `analysis/contacts_hbonds/contacts_parralel_RhoA.py`
  → `rcf_all_sims.csv` (bundled under `analysis/contacts_hbonds/`).
- **F** — PyMOL snapshot of an R5/K6/K7–POPS hydrogen-bonding example.
