#!/usr/bin/env python3
"""Per-lipid-type area-per-lipid (APL) and aggregate monolayer areas for the simple (all-POPC-ish)
vs complex (POPC/PSM/...) bilayers, to answer whether the lower mean APL of the complex membrane is
driven by the small-headgroup sphingomyelin (PSM).

Method: PBC-aware 2-D Voronoi tessellation (freud) of the headgroup P atom of every lipid, split by
leaflet (z vs. the bilayer midplane). Each lipid's Voronoi cell area = its instantaneous APL; the
sum of a leaflet's cells = Lx*Ly (the aggregate monolayer area). Averaged over the final-5 us
sampled trajectory (stride 25). The UPPER leaflet is protein-free in both models, so its per-type
APL is a clean lipid-only comparison; the LOWER leaflet contains the protein, whose footprint is
distributed onto neighbouring lipids (flagged -> lower-leaflet APL is an upper bound).

  conda activate memb_analysis
  python apl_voronoi_compare.py
"""
import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, pandas as pd
import MDAnalysis as mda
import freud

HERE = os.path.dirname(os.path.abspath(__file__))
TRJ = "/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs"
SIMS = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]            # GDP/GTP-1 = simple ; GTP-2/2A/3A = complex
MEMBRANE = {"GDP": "simple", "GTP-1": "simple", "GTP-2": "complex", "2A": "complex", "3A": "complex"}
LIPIDS = ["POPC", "PSM", "POPE", "POPS", "PLA18"]
STRIDE = 50

rows = []          # per-(sim,leaflet,lipid) APL
agg = []           # per-(sim,leaflet) aggregate area + counts

for sim in SIMS:
    print(f"[{sim}] loading...", flush=True)
    u = mda.Universe(f"{TRJ}/topologies/{sim}.psf", f"{TRJ}/{sim}_last5us.dcd")
    P = u.select_atoms("resname POPC PSM POPE POPS PLA18 and name P")
    resn = np.array(P.resnames)
    voro = freud.locality.Voronoi()
    # accumulators: list of per-frame mean area for each (leaflet,lipid); and aggregate area
    acc_area = {}          # (leaf,lip) -> list of per-frame mean cell area
    acc_agg = {}           # (leaf) -> list of per-frame summed area
    acc_n = {}             # (leaf,lip) -> lipid count (constant)
    nfr = 0
    for ts in u.trajectory[::STRIDE]:
        Lx, Ly = ts.dimensions[0], ts.dimensions[1]
        pos = P.positions.copy()
        zmid = np.median(pos[:, 2])
        upper = pos[:, 2] > zmid
        box = freud.box.Box(Lx=Lx, Ly=Ly, is2D=True)
        for leaf, mask in (("upper", upper), ("lower", ~upper)):
            pts = np.zeros((mask.sum(), 3))
            pts[:, 0] = pos[mask, 0] - Lx / 2.0
            pts[:, 1] = pos[mask, 1] - Ly / 2.0
            pts = box.wrap(pts)
            voro.compute((box, pts))
            areas = np.asarray(voro.volumes)      # 2-D -> cell areas (A^2)
            rl = resn[mask]
            acc_agg.setdefault(leaf, []).append(areas.sum())
            for lip in LIPIDS:
                m = rl == lip
                if m.any():
                    acc_area.setdefault((leaf, lip), []).append(areas[m].mean())
                    acc_n[(leaf, lip)] = int(m.sum())
        nfr += 1
        if nfr % 200 == 0:
            print(f"  [{sim}] {nfr} frames...", flush=True)
    Lxy = np.mean([a for a in acc_agg["upper"]])   # ~ Lx*Ly (upper sum); lower sum equal within noise
    for (leaf, lip), vals in acc_area.items():
        v = np.array(vals)
        rows.append(dict(sim=sim, membrane=MEMBRANE[sim], leaflet=leaf, lipid=lip,
                         n=acc_n[(leaf, lip)], apl_mean=v.mean(), apl_sem=v.std(ddof=1) / np.sqrt(len(v))))
    for leaf, vals in acc_agg.items():
        v = np.array(vals)
        rows_n = sum(acc_n[(leaf, lip)] for lip in LIPIDS if (leaf, lip) in acc_n)
        agg.append(dict(sim=sim, membrane=MEMBRANE[sim], leaflet=leaf, n_lipids=rows_n,
                        agg_area_mean=v.mean(), agg_area_sem=v.std(ddof=1) / np.sqrt(len(v)),
                        gross_apl=v.mean() / rows_n))
    print(f"[{sim}] {nfr} frames  box~{Lx:.1f}x{Ly:.1f}  upper_agg={Lxy:.0f} A^2", flush=True)

df = pd.DataFrame(rows).sort_values(["membrane", "sim", "leaflet", "lipid"])
da = pd.DataFrame(agg).sort_values(["membrane", "sim", "leaflet"])
df.to_csv(f"{HERE}/apl_per_type.csv", index=False)
da.to_csv(f"{HERE}/apl_aggregate.csv", index=False)

pd.set_option("display.width", 160, "display.max_rows", 200)
print("\n================ PER-TYPE APL (A^2) ================")
print(df.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
print("\n================ AGGREGATE MONOLAYER AREA ================")
print(da.to_string(index=False, float_format=lambda x: f"{x:.1f}"))

# --- the boss's explicit test: upper leaflet, simple (POPC only) vs complex (POPC+PSM) ---
def upper(sim, lip):
    s = df[(df.sim == sim) & (df.leaflet == "upper") & (df.lipid == lip)]
    return (float(s.apl_mean.iloc[0]), int(s.n.iloc[0])) if len(s) else (np.nan, 0)

print("\n================ BOSS'S TEST: UPPER MONOLAYER ================")
for simple_sim in ["GTP-1", "GDP"]:
    for cx_sim in ["GTP-2"]:
        pc_s, n_pc_s = upper(simple_sim, "POPC")
        pc_c, n_pc_c = upper(cx_sim, "POPC")
        sm_c, n_sm_c = upper(cx_sim, "PSM")
        simple_agg = pc_s * n_pc_s
        complex_agg = pc_c * n_pc_c + sm_c * n_sm_c
        # predicted complex area if its POPC kept the simple-membrane APL (isolate PSM effect)
        pred_complex = pc_s * n_pc_c + sm_c * n_sm_c
        print(f"\n simple [{simple_sim}] upper : POPC APL={pc_s:.2f}  x {n_pc_s}  -> {simple_agg:.0f} A^2")
        print(f" complex[{cx_sim}] upper : POPC APL={pc_c:.2f} x {n_pc_c} + PSM APL={sm_c:.2f} x {n_sm_c}"
              f"  -> {complex_agg:.0f} A^2  ({n_pc_c + n_sm_c} lipids)")
        print(f"   PSM vs complex-POPC APL difference : {sm_c - pc_c:+.2f} A^2")
        print(f"   complex-POPC vs simple-POPC APL     : {pc_c - pc_s:+.2f} A^2 (condensation)")
        print(f"   aggregate area complex - simple     : {complex_agg - simple_agg:+.0f} A^2")
        print(f"   [APL(POPC_simple)x{n_pc_c} + APL(PSM)x{n_sm_c}] = {pred_complex:.0f} A^2"
              f"  vs simple {simple_agg:.0f} A^2  (PSM-only effect)")
print("\nwrote apl_per_type.csv + apl_aggregate.csv")
