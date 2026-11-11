# Catalytic domain–lipid interactions negatively regulate RhoA signaling by modulating conformational dynamics

Analysis code, figure-generation scripts, assembled starting structures, and processed
data accompanying:

> C. M. Hutchins, Y. Zuo, S. Planque, J. Jakubec, J. Frost, A. A. Gorfe.
> *Catalytic domain–lipid interactions negatively regulate RhoA signaling by modulating
> conformational dynamics.* (2026).

This repository lets you **regenerate every quantitative figure panel from the included
processed data**, and documents how to recompute those data from the raw molecular-dynamics
(MD) trajectories (which are too large for GitHub — see
[`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)).

## Simulation systems

Five atomistic systems of geranylgeranylated human RhoA (residues 1–190) on a plasma-
membrane-mimetic bilayer, run on Anton 2/3 (CHARMM36, 310 K, frames every 100 ps).

| System | Nucleotide | Mutation | Bilayer | Length |
|---|---|---|---|---|
| RhoA-GDP    | GDP | WT | simple (POPC / POPC:POPE:POPS) | 20 µs |
| RhoA-GTP-1  | GTP | WT | simple | 20 µs |
| RhoA-GTP-2  | GTP | WT (HVR re-oriented) | complex (PSM:POPC / POPC:POPE:POPS:PlaPE) | 20 µs |
| RhoA-GTP-2A | GTP | R5A/K6A | complex | 10 µs |
| RhoA-GTP-3A | GTP | R5A/K6A/K7A | complex | 10 µs |

Starting coordinates from PDB **1FTN** (GDP) and **1A2B** (GTP); bilayers from
CHARMM-GUI Membrane Builder. Assembled PSF+PDB inputs are in
[`initial_systems/`](initial_systems/).

## Repository layout

```
initial_systems/      Assembled starting structures (PSF + PDB) for all 5 systems
figures/              One folder per paper figure: plotting script + bundled data + README
  figure1 … figure7     main figures
  supplementary/        figS1 (bilayer), figS2 (SwI/Helix4 z-dist), figS3 (cell assays)
analysis/             Upstream trajectory-analysis scripts, by type (see analysis/README.md)
  cobinding_and_network_pipeline/   self-contained Fig 6 & 7 pipeline + small data
  effector_refs/                    effector crystal structures (Fig 7)
experiments/          Cell-assay quantification + source CSVs (Fig 5 / SI S3)
common/               FIGURE_STYLE.md (plotting conventions)
requirements.txt      Python dependencies
DATA_AVAILABILITY.md  Catalogue of the large raw trajectories (deposited separately)
```

## Figure → code map

| Figure | What | Where |
|---|---|---|
| Fig 1B | Cα RMSD + geranylgeranyl insertion depth | `figures/figure1/` |
| Fig 2  | Lobe-COM z-distance + 2D density | `figures/figure2/` |
| Fig 3  | RhoA-GTP-2 reorientation (RMSD, insertion, lobe, density) | `figures/figure3/` |
| Fig 4B,E | Tilt/rotation time series; R5/K6/K7–POPS H-bonds | `figures/figure4/` |
| Fig 5  | NIH3T3 cell assays (F-actin, pMLC2) | `experiments/` |
| Fig 6  | Cooperative lipid co-binding | `analysis/cobinding_and_network_pipeline/` |
| Fig 7  | Dynamical network + Switch I effector competence | `analysis/cobinding_and_network_pipeline/` |
| SI S1  | Bilayer APL + thickness | `figures/supplementary/figS1_bilayer_properties/` |
| SI S2  | Switch I / Helix 4 distance from bilayer center | `figures/supplementary/figS2_switchI_helix4_zdist/` |
| SI S3  | GFP-RhoA F-actin activity | `experiments/` |

(Structural panels — system overviews, snapshots, co-binding and network renders — are
PyMOL figures; the rendering scripts are noted in each figure's README.)

## Reproducing a figure

```bash
pip install -r requirements.txt            # MDAnalysis, numpy, scipy, pandas, matplotlib, lipyphilic, dynetan, networkx
cd figures/figure2 && python fig2_lobe_dynamics.py
```
Each `figures/<fig>/` folder is self-contained: the plotting script reads its bundled
`data/` and writes the panel as `.png`/`.svg`/`.pdf`. See each folder's `README.md`.

## Software
MD: NAMD (equilibration), Anton 2/3 (production), CHARMM36 force field, CHARMM-GUI.
Analysis: MDAnalysis, lipyphilic, SciPy, NumPy, pandas, NetworkX, the Dynetan library
(generalized-correlation network), GROMACS tools; rendering with PyMOL and VMD; cell-assay
statistics with GraphPad Prism / SciPy.

## Naming note
Legacy names in some raw scripts: `GDP`→RhoA-GDP, `GTP`→RhoA-GTP-1, `Resim`→RhoA-GTP-2,
`2A`→RhoA-GTP-2A, `3A`→RhoA-GTP-3A.

## License
© 2026 Chase M. Hutchins, Alemayehu A. Gorfe, and co-authors. Released under the
**Creative Commons Attribution 4.0 International (CC BY 4.0)** license — see
[`LICENSE`](LICENSE). You may share and adapt this material for any purpose with
appropriate credit; please cite the paper above.
