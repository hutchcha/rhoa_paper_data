# Initial simulation systems

Fully assembled, solvated, and neutralized starting structures for the five RhoA
membrane systems in the paper, as built with CHARMM-GUI and used as input to the
production MD runs. Each system is provided as a matched **PSF topology** + **PDB
coordinate** pair (CHARMM36 force field).

| File stem | Nucleotide | Mutation | Bilayer | Production length |
|---|---|---|---|---|
| `RhoA-GDP_initial_system`    | GDP | WT | simple: 190 POPC (upper) / 98 POPC + 60 POPE + 40 POPS (lower) | 20 µs |
| `RhoA-GTP-1_initial_system`  | GTP | WT | simple (as above) | 20 µs |
| `RhoA-GTP-2_initial_system`  | GTP | WT (HVR rotated ~60° about membrane normal) | complex: 100 PSM + 102 POPC (upper) / 50 each POPC, POPE, POPS, PlaPE (lower) | 20 µs |
| `RhoA-GTP-2A_initial_system` | GTP | R5A/K6A | complex (as above) | 10 µs |
| `RhoA-GTP-3A_initial_system` | GTP | R5A/K6A/K7A | complex (as above) | 10 µs |

Each system contains one geranylgeranylated + carboxymethylated RhoA (residues
1–190, prenyl anchor on Cys190), Mg²⁺, GTP or GDP, TIP3P water, and Na⁺/Cl⁻ to
150 mM. The prenyl anchor is inserted in the intracellular (lower) leaflet.

## Provenance
- Catalytic-domain crystal structures: **PDB 1FTN** (GDP) and **PDB 1A2B** (GSP,
  converted to GTP). Missing N-/C-terminal residues modeled as extended chains;
  C-terminal Cys geranylgeranylated + carboxymethylated via CHARMM-GUI PDBReader.
- Bilayers built with the **CHARMM-GUI Membrane Builder**; leaflet lipid counts
  chosen to minimize monolayer area asymmetry (see paper SI, Fig. S1).
- Minimization/equilibration in **NAMD** (CHARMM36); production on **Anton 2/3**,
  NPT, 310 K, frames written every 100 ps.

> Naming note: in some of the raw analysis files the three wild-type systems appear
> under the legacy names `GDP`, `GTP`, and `Resim`, which correspond to
> **RhoA-GDP**, **RhoA-GTP-1**, and **RhoA-GTP-2** respectively.
