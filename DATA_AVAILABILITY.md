# Data availability

This repository contains **all analysis/plotting code and the small processed data
needed to regenerate every quantitative figure panel** in the paper (see each
`figures/.../data/` and `analysis/.../` directory). The **raw MD trajectories are too
large for GitHub** and are catalogued below; they are available from the authors on
request and are being deposited on a public archive (Zenodo/figshare) — add the DOI
here on acceptance.

## What IS in this repository
- `initial_systems/` — all five fully assembled starting structures (PSF + PDB, ~215 MB).
- `figures/` — final plotting scripts + the exact processed-data inputs for each panel.
- `analysis/` — the upstream trajectory-analysis scripts and their small output tables.
- `experiments/` — cell-assay quantification code + source CSVs (Fig. 5 / SI Fig. S3).

## What is NOT in this repository (raw MD data)

### Production trajectories (the primary archival dataset, ~51 GB)
Frames every 100 ps; CHARMM36; topologies in `initial_systems/*.psf`.

| System | Length | Trajectory file (local name) | Size |
|---|---|---|---|
| RhoA-GDP    | 20 µs | `RhoA_GDP.dcd`    | 13 GB |
| RhoA-GTP-1  | 20 µs | `RhoA_GTP-1.dcd`  | 13 GB |
| RhoA-GTP-2  | 20 µs | `RhoA_GTP-2.dcd`  | 13 GB |
| RhoA-GTP-2A | 10 µs | `2A_10us.dcd`     | 6.2 GB |
| RhoA-GTP-3A | 10 µs | `3A_10us.dcd`     | 6.2 GB |

### Derived trajectory excerpts (regenerable; ~217 GB — not deposited)
The Fig. 6/7 (5 µs network/co-binding) and H-bond analyses read per-system excerpts
of the production runs, carried in the local working tree but **not** needed if you
start from the full trajectories above:
- `{GDP,GTP-1,GTP-2,2A,3A}_last5us.dcd` — final 5 µs (~31 GB each; used by the dynetan
  network, effector RMSD, and 2A-similarity analyses).
- `{GDP,GTP-1,GTP-2,2A,3A}_last1.95us.dcd` — H-bond analysis window (~12–13 GB each).

### Intermediate data (regenerable; ~0.7 GB — not deposited)
- `hbonds_stable/{sim}_{gdlip,anclip}_raw_hbonds.npy` — raw per-frame protein–lipid
  H-bond tables produced by `analysis/cobinding_and_network_pipeline/run_hbonds.py`.
  The derived co-binding/region-H-bond CSVs that the figures actually use **are**
  included under `analysis/cobinding_and_network_pipeline/transition_states/`.

### Other large, non-essential assets (not deposited)
PyMOL session files (`.pse`), movie trajectories, and high-resolution `.tif` figure
renders. The figure-rendering scripts regenerate the structural panels from the
PDB/PSF + trajectory inputs.

## Legacy file naming
In some raw analysis scripts the three wild-type systems appear as `GDP`, `GTP`, and
`Resim` → **RhoA-GDP**, **RhoA-GTP-1**, **RhoA-GTP-2**. Mutants are `2A` → RhoA-GTP-2A
and `3A` → RhoA-GTP-3A.

## Reproducing from raw trajectories
Scripts in `analysis/` that read trajectories expect the production `.dcd` (or the
`_last5us.dcd` / `_last1.95us.dcd` excerpts) alongside the matching topology. Obtain the
trajectories from the archive/authors, place them where each script expects (see that
subfolder's README), then run the analysis stage followed by the figure script.
