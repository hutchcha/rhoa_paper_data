#!/usr/bin/env python3
"""How OFTEN do cooperative basic+acidic pairs co-grip a shared AMINE lipid, and
how frequent is the phenomenon overall?

Enhancement says how cooperative a pair is; this says how RELEVANT it is. The
acidic-to-amine mechanism works on any zwitterionic lipid that carries a
donatable amine: POPS (serine -NH3+), POPE (ethanolamine -NH3+), PLA18. POPC is
excluded (quaternary choline, no donatable H). So co-binding is counted on the
SAME lipid MOLECULE of ANY of these types.

Per pair: occupancy = fraction of POST frames co-gripping >=1 shared amine lipid.
Aggregate: fraction of POST frames with ANY basic+acidic co-binding (the avidity
signature -- specific pairs rotate, but the phenomenon is nearly always on), with
a POPS/POPE/PLA18 breakdown.

  python cobinding_occupancy.py [--lipids POPS,POPE,PLA18] [--min-occ 0.005]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import stable_state_config as C
from lipid_hbond_analysis import load_bonds, BASIC, ACIDIC

OUT = "transition_states/gains"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="GDP,GTP-1,GTP-2,2A,3A")
    ap.add_argument("--lipids", default="POPS,POPE,PLA18")
    ap.add_argument("--min-occ", type=float, default=0.005)
    ap.add_argument("--tag", default=None,
                    help="output filename suffix; default 'POPS' if only POPS else 'amine'")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    order = args.sims.split(",")
    lips = args.lipids.split(",")
    sims = sorted([s for s in C.SIMS if s["name"] in order], key=lambda s: order.index(s["name"]))
    ltag = "/".join(lips)
    suf = args.tag if args.tag else ("POPS" if lips == ["POPS"] else "amine")

    rows, occ_by_sim, lab, summary = [], {}, {}, []
    for sim in sims:
        name = sim["name"]
        df = load_bonds(sim)
        npost = int(df["npost"].iloc[0])
        pops = df[(df.lrn.isin(lips)) & (df.pclass != "other")]
        # molecule key includes resname so resids can't collide across lipid types
        bset, aset = defaultdict(set), defaultdict(set)
        for r in pops.itertuples():
            key = (r.frame, r.lres, r.lrn)
            (bset if r.pclass == "basic" else aset)[r.pres].add(key)
            lab[r.pres] = C.res_label(r.prn, r.pres)
        occ = {}
        for b, bs in bset.items():
            for a, as_ in aset.items():
                inter = bs & as_
                if inter:
                    frames = {fr for fr, lr, ln in inter}
                    occ[(b, a)] = len(frames) / npost
                    # per-lipid frame counts for this pair
                    perlip = {L: len({fr for fr, lr, ln in inter if ln == L}) for L in lips}
                    rows.append({"sim": name, "basic": lab[b], "acidic": lab[a],
                                 "basic_resid": b, "acidic_resid": a,
                                 "co_frames": len(frames), "n_post": npost,
                                 "occupancy": len(frames) / npost,
                                 **{f"frames_{L}": perlip[L] for L in lips}})
        occ_by_sim[name] = occ
        # aggregate: any co-binding on any amine lipid, per frame; + per-lipid breakdown
        fl_b = set().union(*bset.values()) if bset else set()
        fl_a = set().union(*aset.values()) if aset else set()
        co_fl = fl_b & fl_a
        union_frames = {fr for fr, lr, ln in co_fl}
        s = {"sim": name, "n_post": npost, "frames_with_cobinding": len(union_frames),
             "frac": len(union_frames) / npost, "total_events": len(co_fl)}
        for L in lips:
            uf = {fr for fr, lr, ln in co_fl if ln == L}
            s[f"frac_{L}"] = len(uf) / npost
        summary.append(s)
        print(f"[{name}] ANY co-binding in {len(union_frames)}/{npost} = "
              f"{len(union_frames)/npost*100:.1f}% of post frames | by lipid: "
              + ", ".join(f"{L} {s[f'frac_{L}']*100:.1f}%" for L in lips), flush=True)

    tab = pd.DataFrame(rows).sort_values(["sim", "occupancy"], ascending=[True, False])
    tab.to_csv(f"{OUT}/cobinding_occupancy_{suf}.csv", index=False)
    summ = pd.DataFrame(summary)
    summ.to_csv(f"{OUT}/cobinding_frequency_{suf}.csv", index=False)

    print("\n=== top pairs by occupancy (fraction of post, any amine lipid) ===")
    for name in order:
        d = tab[tab.sim == name].head(8)
        if not len(d):
            print(f"  {name}: none"); continue
        print(f"  {name}: " + " | ".join(
            f"{r.basic}+{r.acidic} {r.occupancy*100:.1f}%" for r in d.itertuples()))

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C.apply_style()
    simn = [s["name"] for s in sims]

    # ---------- Figure 1: per-pair occupancy heatmaps ----------
    vmax = max((max(o.values()) if o else 0) for o in occ_by_sim.values())
    vmax = max(0.05, vmax)
    fig, axes = plt.subplots(1, len(simn), figsize=(3.6 * len(simn) + 1.5, 6), squeeze=False)
    for ax, name in zip(axes[0], simn):
        occ = occ_by_sim[name]
        basics = sorted({b for b, a in occ}); acids = sorted({a for b, a in occ})
        if not basics:
            ax.axis("off"); ax.set_title(f"RhoA-{name}\n(no pairs)", fontsize=12); continue
        M = np.full((len(basics), len(acids)), np.nan)
        for i, b in enumerate(basics):
            for j, a in enumerate(acids):
                if (b, a) in occ:
                    M[i, j] = occ[(b, a)]
        im = ax.imshow(M, cmap="magma", vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(acids))); ax.set_xticklabels([lab[a] for a in acids],
                                                             rotation=90, fontsize=8, color="#7b241c")
        ax.set_yticks(range(len(basics))); ax.set_yticklabels([lab[b] for b in basics],
                                                              fontsize=8, color="#1f4e79")
        for i in range(len(basics)):
            for j in range(len(acids)):
                if np.isfinite(M[i, j]) and M[i, j] >= args.min_occ:
                    ax.text(j, i, f"{M[i,j]*100:.0f}", ha="center", va="center",
                            fontsize=7, color="white" if M[i, j] < vmax * 0.6 else "black")
        ax.set_title(f"RhoA-{name}", fontsize=12, fontweight="bold")
        if name == simn[0]:
            ax.set_ylabel("basic residue (→ phosphate)", fontsize=10)
    fig.suptitle(f"How often each basic+acidic pair co-grips a shared amine lipid "
                 f"({ltag}); numbers = % of post frames", fontsize=12, fontweight="bold")
    cax = fig.add_axes([0.93, 0.25, 0.010, 0.5])
    fig.colorbar(im, cax=cax, label="fraction of post frames co-binding")
    fig.subplots_adjust(left=0.04, right=0.91, top=0.88, bottom=0.12, wspace=0.35)
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/cobinding_occupancy_heatmap_{suf}.{e}", dpi=160)
    plt.close(fig)

    # ---------- Figure 2: overall frequency (grouped: ANY + per lipid; not stacked,
    # because a frame can co-bind two lipid types at once) ----------
    fig2, ax = plt.subplots(figsize=(1.7 * len(simn) + 3, 5))
    cats = ["ANY"] + lips
    COL = {"ANY": "#333333", "POPS": "#E45756", "POPE": "#4C78A8", "PLA18": "#B279A2"}
    x = np.arange(len(simn)); w = 0.8 / len(cats)
    for k, cat in enumerate(cats):
        field = "frac" if cat == "ANY" else f"frac_{cat}"
        vals = np.array([next(s[field] for s in summary if s["sim"] == n) for n in simn]) * 100
        bars = ax.bar(x + (k - (len(cats) - 1) / 2) * w, vals, w, color=COL.get(cat, "0.5"),
                      edgecolor="k", linewidth=0.4, label=cat)
        if cat == "ANY":
            for xi, v in zip(x, vals):
                ax.text(xi - (len(cats) - 1) / 2 * w, v + 1, f"{v:.0f}%", ha="center",
                        fontsize=9, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([f"RhoA-{n}" for n in simn], fontsize=11)
    ax.set_ylabel("% of post-state frames with ≥1\nbasic+acidic co-binding a shared amine lipid", fontsize=10)
    ax.set_title("Overall frequency of the cooperative co-binding phenomenon\n"
                 "(specific pairs interchange; the phenomenon persists = avidity)", fontsize=12, loc="left")
    ax.legend(frameon=False, fontsize=10, title="shared lipid (ANY = union)")
    ax.spines[["right", "top"]].set_visible(False)
    fig2.tight_layout()
    for e in ("png", "svg"):
        fig2.savefig(f"{OUT}/cobinding_frequency_{suf}.{e}", dpi=160, bbox_inches="tight")
    plt.close(fig2)
    print(f"\nwrote {OUT}/cobinding_*_{suf}.png and CSVs")


if __name__ == "__main__":
    main()
