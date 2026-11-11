# RhoA full-resolution H-bond + cooperativity analysis

Re-runs the protein–lipid H-bond and basic/acidic **cooperativity** analysis on the
**full-sampling** (0.1 ns) matched post-state excerpts — equal frame count across
all five sims, all lipids included. Self-contained: rsync this folder to the GPU
machine, drop in the excerpts, and launch in the background.

## 1. One-time setup (on the GPU machine)
```bash
conda env create -f environment.yml
conda activate rhoa_analysis
```
(Or use any existing env with MDAnalysis ≥ 2.7, numpy, scipy, pandas, matplotlib.)

## 2. Topologies (PSFs) — already bundled
The five PSFs are in **`topologies/`** and are synced with this folder, so you
don't need to hunt for them on the GPU machine:
```
topologies/{2A,3A,GTP-2,GDP,GTP-1}.psf
```

## 3. Point the config at your files  (`stable_state_config.py`, `PATHS` block)
Each sim has an editable `top` (PSF) and `traj` (the 1.95-us excerpt). Defaults:
- `top`  → the bundled `topologies/<sim>.psf`
- `traj` → `rhoa_full_sampled_trajs/<sim>_last1.95us.dcd`

**You can point `traj` anywhere** (e.g. a scratch dir on the GPU where you wrote the
excerpts) so you don't move the big files around — just edit the path.

Also set at the top of the config:
- **`DT_NS`** — full trajectories' native per-frame time in ns (e.g. `0.1`). NOT 1.008.
- **`POST_US`** — the extracted window (`1.95`). Must satisfy
  `round(POST_US*1000/DT_NS) == excerpt frame count` (occupancy denominators use it;
  at DT_NS=0.1, POST_US=1.95 → 19500 frames).

## 4. Pre-flight check (do this before the long run)
```bash
python verify_inputs.py
```
Loads every PSF+trajectory pair and prints atom + frame counts. A PSF that doesn't
match its excerpt (wrong atom count) fails loudly here; it also warns if the frame
count != round(POST_US/DT_NS). `run_all.sh` runs this automatically and aborts if it fails.

## 5. Launch in the background
```bash
tmux new -d -s rhoa 'conda activate rhoa_analysis; bash run_all.sh 2>&1 | tee run_all.log'
tmux attach -t rhoa      # watch;  Ctrl-b d  to detach
# or:  nohup bash run_all.sh > run_all.log 2>&1 &
```
Resumable — relaunch after any interruption and it skips finished H-bond searches.

## Pipeline (`run_all.sh`)
1. `run_hbonds.py` — G-domain↔lipid + anchor↔lipid H-bonds, **all lipids, stride 1** → `hbonds_stable/`
2. `lipid_hbond_analysis.py` — headgroup map + per-pair cooperativity enhancement (POPS and all-amine)
3. `viz_cooperativity.py` — bipartite cooperativity networks (min25 + min50)
4. `cobinding_occupancy.py` — per-pair occupancy heatmaps + overall frequency

## Outputs
- `hbonds_stable/{sim}_{gdlip,anclip}_raw_hbonds.npy` — raw H-bond tables (re-tallyable)
- `transition_states/gains/` — maps, networks, occupancy/frequency figures + CSVs
  (`_amine` = all zwitterionic lipids; no suffix / `_POPS` = POPS only)

rsync `hbonds_stable/` and `transition_states/` back to compare against the
decimated (1 ns) results.

## Notes / expectations
- **Runtime:** all-lipid, stride-1 H-bond search over 5 × ~19.5k frames is the long
  part — expect **many hours** (order ~10–20 h depending on the machine). Background it.
- **RAM:** the cooperativity step loads all five sims' bonds into memory at once;
  the all-lipid full-res tables are large (a few GB). A workstation with ≥ 32 GB is comfortable.
- **Disk:** raw H-bond arrays can total a few GB.
- Only protein↔lipid H-bonds are computed (that's all the cooperativity needs) —
  no intra-protein or vdW-contact searches in this bundle.

---

## Allosteric analysis add-on (MI generalized correlation + network)

Separate from the H-bond pipeline above. Computes mutual-information generalized
correlation matrices (Lange–Grubmüller r_MI) of Cα motion for the network analysis.
Self-contained — uses only `stable_state_config.py` from this bundle, and resolves
trajectories at `<sim>_last<POST_US>us.dcd` (falls back automatically if the config
`traj` path is a dangling symlink).

```bash
conda activate rhoa_analysis          # needs MDAnalysis, numpy, scipy, matplotlib
# full 5-way run (uses all frames, frame-matched to the common minimum):
nohup python -u dccm_mi.py --sims 2A,3A,GTP-2,GTP-1,GDP > dccm_mi.log 2>&1 &
# quick check first (stride 5):
python dccm_mi.py --sims 2A,3A,GTP-2,GTP-1,GDP --stride 5
```

Requires all five `<sim>_last<POST_US>us.dcd` present (GTP-1/GTP-2 were NOT on the
analysis workstation — run this where all five excerpts live). `verify_inputs.py`
lists what is present.

Outputs → `transition_states/dccm_mi/`:
- `mi_<sim>.npy`, `mi_resids.npy` — matrices for the network step
- `mi_all.png` — all sims; `mi_diff.png` — difference vs WT (what the mutation rewires)
- `npatch_mi_profile.png` — N-patch(5-7) ↔ residue correlation, regions highlighted

**Runtime:** fast (minutes) — it's covariances + determinants, not an H-bond search.
The network step (betweenness / optimal N-patch→Switch I paths / communities) runs
off the saved `mi_*.npy` matrices and per-sim contact maps (4.5 Å / 75 % occupancy).

### Network step (`network_allostery.py`)
Runs off the saved `mi_*.npy` + per-sim contact maps (4.5 Å / 75 % occupancy).
```bash
python network_allostery.py --sims 2A,3A,GTP-2,GTP-1,GDP --contact-stride 10
```
Outputs → `transition_states/network/`: `betweenness.{png,csv}` (allosteric hubs),
`npatch_swI_path.{png}` + `npatch_swI_paths.csv` (optimal + suboptimal N-patch→Switch I
paths, per sim), `communities.{png,csv}` (dynamical modules). Edge distance = −ln r_MI.

---

## Porting the dynetan network analysis to a many-core machine

dynetan's generalized correlation (k-NN mutual information) is **CPU/numba-parallel** — it scales
with CPU cores, not the GPU — and is slow single-threaded, so run it where there are many cores.

**Files needed (all in this bundle):** `dynetan_network.py`, `render_dynetan.py`,
`environment_dynetan.yml`, `run_dynetan_gpu.sh`, `topologies/<sim>.psf`, and the five
`<sim>_last1.95us.dcd` excerpts (in this folder).

```bash
conda env create -f environment_dynetan.yml     # one-time
conda activate rhoa_dynetan
bash run_dynetan_gpu.sh                           # uses all cores; NCORES/SAMPLED/CSTRIDE overridable
```
`run_dynetan_gpu.sh` sets `NUMBA_NUM_THREADS=$(nproc)`, prints which excerpts it found, runs
`dynetan_network.py` over all 5 sims, then `render_dynetan.py`. Outputs land in
`transition_states/dynetan/` — copy that folder back to the analysis machine.

Outputs: `<sim>_nodes.csv` (resid, community, betweenness), `<sim>_edges.csv`, `<sim>_path.txt`
(N-patch→Switch I optimal path), `<sim>_{comm,btw}.pdb`, and the rendered comparison images
`dynetan_communities.png` / `dynetan_betweenness.png`.
