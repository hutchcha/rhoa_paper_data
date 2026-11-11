#!/usr/bin/env python3
"""Whole-protein H-bond occupancy by region: pre / cross / post, lipid + intra,
backbone and sidechain, one figure per structural region.

Covers every residue 4-190. Regions come from stable_state_config, with HELIX 5
SPLIT using what the backbone-occupancy profile showed:
    Helix 5 (main)  164-175   fully helical, backbone occ ~1.0, state-independent
    Helix 5 hinge   176-183   a5 C-term through the disorder onset at R182/R183
    HVR+anchor      184-190
Residues in no defined region are grouped into contiguous "other N-M" blocks so
nothing is missed.

occupancy = fraction of that state's frames with >=1 H-bond.
  lipid : residue <-> any lipid   (gdlip array for 4-179, anclip for 180-190)
  intra : residue <-> other protein residue, |i-j| >= minsep

NB no atom-exclusivity dedup here: occupancy asks only whether >=1 bond exists,
and dedup can never remove every bond of a residue, so it cannot change the
answer -- skipping it makes the whole-protein pass fast.

  python plot_region_hbonds.py [--sims 3A,2A,GTP-2] [--minsep 2]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse, json
import numpy as np
import pandas as pd
import MDAnalysis as mda
import stable_state_config as C

OUT = "transition_states/region_hbonds"
HB = "hbonds_stable"
BACKBONE = {"N", "CA", "C", "O", "OT1", "OT2", "OXT"}
STATES3 = ("pre", "cross", "post")
COL = {("lipid", "pre"): "#F7C6C1", ("lipid", "cross"): "#E8756A", ("lipid", "post"): "#A93226",
       ("intra", "pre"): "#C3D9F0", ("intra", "cross"): "#5B8FC9", ("intra", "post"): "#17405F"}
LO, HI = 4, 190
BLOCK_NS = 50.0     # block-average window for error bars

def build_regions():
    """Ordered [(name, [resids])] covering LO..HI, with helix 5 split."""
    base = {k: list(v) for k, v in C.REGIONS.items()}
    base.pop("Helix 5", None); base.pop("HVR+anchor", None)
    base["Helix 5 (main)"] = [(164, 175)]
    base["Helix 5 hinge"] = [(176, 183)]
    base["HVR+anchor"] = [(184, 190)]
    assign = {}
    for name, ranges in base.items():
        for a, b in ranges:
            for r in range(a, b + 1):
                if LO <= r <= HI:
                    assign[r] = name
    # contiguous unassigned blocks -> "other a-b"
    r = LO
    while r <= HI:
        if r not in assign:
            start = r
            while r <= HI and r not in assign:
                r += 1
            for q in range(start, r):
                assign[q] = f"other {start}-{r-1}"
        else:
            r += 1
    order, seen = [], set()
    for r in range(LO, HI + 1):
        n = assign[r]
        if n not in seen:
            seen.add(n); order.append(n)
    return [(n, [r for r in range(LO, HI + 1) if assign[r] == n]) for n in order]

def residue_frames(sim, minsep, lipids=None):
    """-> dict[(resid, part, kind)] = unique frames; analysed grid; resnames.

    lipids: iterable of lipid resnames to count (default all C.LIPIDS). Pass
    ["POPS"] for a like-for-like comparison with the GDP/GTP-1 references,
    whose raw arrays were searched POPS-only."""
    name = sim["name"]
    lipset = set(lipids) if lipids else set(C.LIPIDS)
    u = mda.Universe(sim["top"], sim["traj"])
    resid, resname, names = u.atoms.resids, u.atoms.resnames, u.atoms.names
    is_lip = np.array([rn in lipset for rn in resname])
    out, analysed = {}, None
    # ---- lipid: gdlip covers 4-179, anclip covers 180-190 ----
    for suf in ("gdlip", "anclip"):
        p, m = f"{HB}/{name}_{suf}_raw_hbonds.npy", f"{HB}/{name}_{suf}_meta.json"
        if not (os.path.exists(p) and os.path.exists(m)):
            continue
        mm = json.load(open(m)); analysed = np.arange(mm["start"], mm["stop"], mm["stride"])
        hb = np.load(p)
        d = hb[:, 1].astype(int); a = hb[:, 3].astype(int); fr = hb[:, 0].astype(int)
        dl, al = is_lip[d], is_lip[a]; k = dl ^ al
        prot = np.where(dl, a, d)[k]; frk = fr[k]
        bb = np.isin(names[prot], list(BACKBONE))
        pr = resid[prot]
        df = pd.DataFrame({"r": pr, "bb": bb, "f": frk})
        for (r, isbb), g in df.groupby(["r", "bb"]):
            key = (int(r), "bb" if isbb else "sc", "lipid")
            prev = out.get(key, np.array([], int))
            out[key] = np.unique(np.concatenate([prev, g["f"].values]))
    # ---- intra-protein ---- (meta filename differs between pipelines:
    # the 3A/2A/GTP-2 runs wrote *_intra_hbonds_meta.json, run_gdp_gtp1.py
    # wrote *_intra_meta.json)
    p = f"{HB}/{name}_intra_raw_hbonds.npy"
    m = next((c for c in (f"{HB}/{name}_intra_hbonds_meta.json",
                          f"{HB}/{name}_intra_meta.json") if os.path.exists(c)), "")
    if os.path.exists(p) and os.path.exists(m):
        mm = json.load(open(m))
        if analysed is None:
            analysed = np.arange(mm["start"], mm["stop"], mm["stride"])
        hb = np.load(p)
        d = hb[:, 1].astype(int); a = hb[:, 3].astype(int); fr = hb[:, 0].astype(int)
        ri, rj = resid[d], resid[a]
        far = np.abs(ri - rj) >= minsep
        for side, me, other in (("d", d, rj), ("a", a, ri)):
            sel = far
            rr = resid[me][sel]; bb = np.isin(names[me][sel], list(BACKBONE)); ff = fr[sel]
            df = pd.DataFrame({"r": rr, "bb": bb, "f": ff})
            for (r, isbb), g in df.groupby(["r", "bb"]):
                key = (int(r), "bb" if isbb else "sc", "intra")
                prev = out.get(key, np.array([], int))
                out[key] = np.unique(np.concatenate([prev, g["f"].values]))
    resn = {}
    for r in range(LO, HI + 1):
        sel = u.select_atoms(f"(protein or resname CYSG) and resid {r} and name CA")
        resn[r] = sel.resnames[0] if len(sel) else None
    return out, analysed, resn

def occ_blocks(frames_hit, w, analysed, block_ns=50.0, dt_ns=1.008):
    """Occupancy computed in fixed-length blocks -> (mean, sem, n_blocks).

    Whole-window occupancy hides how variable the engagement is; 50 ns blocks
    give an error bar and are long enough to contain many H-bond make/break
    events but short enough to yield several blocks even in the pre/cross
    windows."""
    f0, f1 = w
    grid = analysed[(analysed >= f0) & (analysed < f1)]
    if grid.size == 0:
        return np.nan, np.nan, 0
    hit = np.asarray(frames_hit)
    inwin = np.isin(grid, hit) if hit.size else np.zeros(len(grid), bool)
    per_block = max(2, int(round(block_ns / dt_ns)))
    nb = len(grid) // per_block
    if nb < 1:
        return float(inwin.mean()), np.nan, 1
    vals = np.array([inwin[b * per_block:(b + 1) * per_block].mean() for b in range(nb)])
    sem = vals.std(ddof=1) / np.sqrt(nb) if nb > 1 else np.nan
    return float(vals.mean()), float(sem), nb

def occ(frames_hit, w, analysed):
    f0, f1 = w
    n = int(((analysed >= f0) & (analysed < f1)).sum())
    if n == 0:
        return np.nan
    h = np.asarray(frames_hit)
    return float(((h >= f0) & (h < f1)).sum()) / n if h.size else 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="3A,2A,GTP-2")
    ap.add_argument("--minsep", type=int, default=2)
    ap.add_argument("--block-ns", type=float, default=50.0)
    ap.add_argument("--lipid", default="all",
                    help="'all' (every membrane lipid) or a single resname like "
                         "POPS for a like-for-like comparison across all 5 sims")
    args = ap.parse_args()
    globals()["BLOCK_NS"] = args.block_ns
    os.makedirs(OUT, exist_ok=True)
    import matplotlib.pyplot as plt
    C.apply_style()
    sims = [s for s in C.SIMS if s["name"] in args.sims.split(",")]
    regions = build_regions()
    lipids = None if args.lipid == "all" else [args.lipid]
    suf = "" if args.lipid == "all" else f"_{args.lipid}"
    lipttl = "all lipids" if args.lipid == "all" else f"{args.lipid} only"

    rows, data, resn_all = [], {}, {}
    for sim in sims:
        got = residue_frames(sim, args.minsep, lipids)
        if got is None:
            continue
        fr, analysed, resn = got
        data[sim["name"]] = (fr, analysed)
        resn_all.update({k: v for k, v in resn.items() if v})
        for st in STATES3:
            w = C.state_frames(sim["name"], st, sim["dt_ns"])
            if not w:
                continue
            for r in range(LO, HI + 1):
                if resn.get(r) is None:
                    continue
                for part in ("bb", "sc"):
                    rec = {"sim": sim["name"], "resid": r, "resname": resn[r],
                           "part": part, "state": st}
                    for kind in ("lipid", "intra"):
                        m, e, nb = occ_blocks(fr.get((r, part, kind), np.array([], int)),
                                              w, analysed, BLOCK_NS, sim["dt_ns"])
                        rec[kind] = m; rec[f"{kind}_sem"] = e; rec[f"{kind}_nblocks"] = nb
                    rows.append(rec)
        print(f"[{sim['name']}] tallied", flush=True)
    if not rows:
        print("no data"); return
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/region_hbond_occupancy{suf}.csv", index=False)

    simnames = [s["name"] for s in sims if s["name"] in set(df["sim"])]
    cats = [(k, st) for k in ("lipid", "intra") for st in STATES3]

    # ONE large figure per sim: rows = regions, cols = backbone|sidechain,
    # x within a row = [ALL(region aggregate)] + each individual residue.
    maxgroups = max(len([r for r in res if r in resn_all]) + 1 for _, res in regions)
    for s in simnames:
        fr, analysed = data[s]
        nrow = len(regions)
        fig, axes = plt.subplots(nrow, 2,
                                 figsize=(max(10, 0.62 * maxgroups) * 2, 2.45 * nrow),
                                 squeeze=False)
        for i, (rname, resids) in enumerate(regions):
            present = [r for r in resids if r in resn_all]
            for j, part in enumerate(("bb", "sc")):
                ax = axes[i][j]
                if not present:
                    ax.axis("off"); continue
                groups = ["ALL"] + [f"{resn_all[r]}{r}" for r in present]
                x = np.arange(len(groups)); w = 0.13
                for k, (kind, st) in enumerate(cats):
                    wnd = C.state_frames(s, st, 1.008)
                    if not wnd:
                        continue
                    agg_frames = np.unique(np.concatenate(
                        [fr.get((r, part, kind), np.array([], int)) for r in present]
                        or [np.array([], int)]))
                    m0, e0, _ = occ_blocks(agg_frames, wnd, analysed, BLOCK_NS, 1.008)
                    vals, errs = [m0], [e0]
                    for r in present:
                        m1, e1, _ = occ_blocks(fr.get((r, part, kind), np.array([], int)),
                                               wnd, analysed, BLOCK_NS, 1.008)
                        vals.append(m1); errs.append(e1)
                    ax.bar(x + (k - (len(cats) - 1) / 2) * w, vals, width=w,
                           yerr=errs, error_kw=dict(elinewidth=0.5, capsize=1.2, ecolor="0.3"),
                           color=COL[(kind, st)], edgecolor="k", linewidth=0.3,
                           label=f"{kind} {st}")
                ax.set_xticks(x)
                ax.set_xticklabels(groups, fontsize=7.5, rotation=45, ha="right")
                ax.set_ylim(0, 1.05)
                ax.spines[["right", "top"]].set_visible(False)
                ax.set_title(f"{rname} — {'backbone' if part=='bb' else 'sidechain'}",
                             fontsize=10, loc="left")
                if j == 0:
                    ax.set_ylabel("H-bond occ.", fontsize=9)
                if i == 0 and j == 0:
                    ax.legend(frameon=False, fontsize=7, ncol=6)
        fig.suptitle(f"RhoA-{s}: H-bond occupancy by region — pre / cross / post "
                     f"(ALL = region aggregate, then residues; lipid = {lipttl}; "
                     f"{BLOCK_NS:.0f} ns block means +/- SEM)", fontsize=15, y=0.999)
        fig.tight_layout(rect=[0, 0, 1, 0.995])
        for e in ("png", "svg"):
            fig.savefig(f"{OUT}/ALL_regions_{s}{suf}.{e}", dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {OUT}/ALL_regions_{s}{suf}.png  ({nrow} regions x 2 parts)")

    # quick survey: biggest |post-pre| changes protein-wide
    piv = df.pivot_table(index=["sim", "resid", "resname", "part"], columns="state",
                         values=["lipid", "intra"])
    out = []
    for kind in ("lipid", "intra"):
        t = piv[kind].copy()
        t["delta"] = t["post"] - t["pre"]
        t["kind"] = kind
        out.append(t.reset_index())
    ch = pd.concat(out, ignore_index=True).dropna(subset=["delta"])
    ch.to_csv(f"{OUT}/region_hbond_changes{suf}.csv", index=False)
    print("\n=== largest |post - pre| occupancy changes (any residue, protein-wide) ===")
    top = ch.reindex(ch["delta"].abs().sort_values(ascending=False).index).head(25)
    for _, r in top.iterrows():
        print(f"  {r['sim']:>6} {r['resname']}{int(r['resid']):<4} {r['part']} {r['kind']:>5}: "
              f"{r['pre']:.2f} -> {r['post']:.2f}  ({r['delta']:+.2f})")

if __name__ == "__main__":
    main()
