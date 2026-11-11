#!/usr/bin/env python3
"""Map protein<->lipid H-bonds by headgroup chemistry, and test basic/acidic
COOPERATIVITY on single lipids.

Mechanism established: basic residues (Arg/Lys) DONATE to lipid PHOSPHATE oxygens;
acidic residues (Asp/Glu) ACCEPT from the lipid AMINE (-NH3+ of POPS serine /
POPE ethanolamine). Both salt-bridge to opposite ends of the zwitterionic head.

Part A -- MAP: per residue, post-state occupancy of H-bonding each lipid type's
  PHOSPHATE vs AMINE group -> two heatmaps (residues x lipid type).
Part B -- COOPERATIVITY: for each lipid MOLECULE per frame, does it have a basic
  partner (phosphate) AND an acidic partner (amine) at the same time? Test whether
  basic+acidic co-occur on the SAME lipid more than chance (odds ratio / lift /
  conditional probabilities), and rank the co-binding basic-acidic residue pairs.
  Hypothesis (user + PI): a basic residue damps the phosphate charge, reducing
  repulsion and enabling acidic residues to H-bond the amine -> POSITIVE association.

  python lipid_hbond_analysis.py [--sims 3A,2A,GTP-2] [--lipid POPS]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse, json
import numpy as np
import pandas as pd
import MDAnalysis as mda
import stable_state_config as C

HB = "hbonds_stable"
OUT = "transition_states/gains"
BASIC = {"ARG", "LYS", "HSP"}
ACIDIC = {"ASP", "GLU"}


def lipgroup(atom):
    # POPS headgroup has THREE charged groups: phosphate (-), serine amine (+),
    # serine CARBOXYLATE (-). O13A/O13B (q=-0.67 each) are the serine carboxylate,
    # NOT phosphate -- they must be split out (O13/O14 at -0.78 are phosphate).
    if atom == "N":
        return "amine"
    if atom in {"O13A", "O13B", "OC1", "OC2", "OT1", "OT2"}:
        return "carboxylate"
    if atom in {"O11", "O12", "O13", "O14", "P", "O1P", "O2P", "O3P", "O4P"}:
        return "phosphate"
    if atom in {"O21", "O22", "O31", "O32"}:
        return "ester"
    return "other"


def load_bonds(sim):
    """All protein<->lipid H-bonds in the POST window as a DataFrame."""
    name = sim["name"]
    u = mda.Universe(sim["top"], sim["traj"])
    resid, resname, names = u.atoms.resids, u.atoms.resnames, u.atoms.names
    lipset = set(C.LIPIDS)
    is_lip = np.array([rn in lipset for rn in resname])
    w = C.state_frames(name, "post", sim["dt_ns"])
    frames = []
    for tag in ("gdlip", "anclip"):
        p = f"{HB}/{name}_{tag}_raw_hbonds.npy"
        if not os.path.exists(p):
            continue
        hb = np.load(p)
        d = hb[:, 1].astype(int); a = hb[:, 3].astype(int); fr = hb[:, 0].astype(int)
        dl, al = is_lip[d], is_lip[a]
        keep = (dl ^ al) & (fr >= w[0]) & (fr < w[1])
        prot = np.where(dl, a, d)[keep]; lip = np.where(dl, d, a)[keep]
        role = np.where(dl[keep], "acceptor", "donor")
        frames.append(pd.DataFrame({
            "frame": fr[keep], "pres": resid[prot], "prn": resname[prot],
            "patom": names[prot], "role": role,
            "lres": resid[lip], "lrn": resname[lip], "latom": names[lip]}))
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    df["pclass"] = np.where(df["prn"].isin(BASIC), "basic",
                            np.where(df["prn"].isin(ACIDIC), "acidic", "other"))
    df["lgroup"] = df["latom"].map(lipgroup)
    df["npost"] = w[1] - w[0]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="3A,2A,GTP-2")
    ap.add_argument("--lipids", default="POPS",
                    help="one lipid (POPS) or a comma list (POPS,POPE,PLA18) for co-binding")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lips = args.lipids.split(",")
    suf = "" if lips == ["POPS"] else "_" + ("amine" if set(lips) >= {"POPS", "POPE"} else "".join(lips))
    ltag = "/".join(lips)
    sims = [s for s in C.SIMS if s["name"] in args.sims.split(",")]
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C.apply_style()

    dfs = {s["name"]: load_bonds(s) for s in sims}

    # ---------- Part A: MAP (residue x lipid type, phosphate vs amine) ----------
    maprows = []
    for name, df in dfs.items():
        npost = df["npost"].iloc[0]
        for (pres, prn, lrn, lg), g in df.groupby(["pres", "prn", "lrn", "lgroup"]):
            maprows.append({"sim": name, "resid": pres, "resname": prn,
                            "pclass": "basic" if prn in BASIC else ("acidic" if prn in ACIDIC else "other"),
                            "lipid": lrn, "group": lg, "occ": g["frame"].nunique() / npost})
    mp = pd.DataFrame(maprows)
    mp.to_csv(f"{OUT}/lipid_hbond_map{suf}.csv", index=False)

    # residues that reach >=0.1 occ to any charged group in any sim
    GROUPS = [("phosphate", "binds PHOSPHATE (basic → −)"),
              ("carboxylate", "binds Ser CARBOXYLATE (basic → −)\nPOPS only"),
              ("amine", "binds AMINE (acidic → +)")]
    key = mp[(mp.group.isin([g for g, _ in GROUPS])) & (mp.occ >= 0.10)]
    resids = sorted(set(key["resid"]), key=lambda r: (0 if key[key.resid==r]["pclass"].iloc[0]=="basic" else 1, r))
    reslab = {r: C.res_label(mp[mp.resid==r]["resname"].iloc[0], r) for r in resids}
    lipids = [l for l in C.LIPIDS if l in set(mp["lipid"])]
    simn = list(dfs.keys())
    # average occ across sims for the map
    fig, axes = plt.subplots(1, len(GROUPS), figsize=(2.4 * len(lipids) * len(GROUPS) * 0.5 + 4, 0.32 * len(resids) + 2),
                             squeeze=False, constrained_layout=True)
    for ax, (grp, ttl) in zip(axes[0], GROUPS):
        M = np.full((len(resids), len(lipids)), np.nan)
        for i, r in enumerate(resids):
            for j, L in enumerate(lipids):
                v = mp[(mp.resid == r) & (mp.lipid == L) & (mp.group == grp)]["occ"]
                if len(v):
                    M[i, j] = v.mean()
        im = ax.imshow(M, cmap="viridis", vmin=0, vmax=max(0.3, np.nanmax(M) if np.isfinite(M).any() else 0.3),
                       aspect="auto")
        ax.set_xticks(range(len(lipids))); ax.set_xticklabels(lipids, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(resids)))
        ax.set_yticklabels([reslab[r] for r in resids], fontsize=8)
        for i, r in enumerate(resids):
            ax.get_yticklabels()[i].set_color("#1f77b4" if mp[mp.resid==r]["pclass"].iloc[0]=="basic" else "#d62728")
        ax.set_title(ttl, fontsize=11)
        fig.colorbar(im, ax=ax, shrink=0.6, label="post occupancy")
    fig.suptitle("Protein–lipid H-bond map (blue=basic, red=acidic residues; mean over sims)",
                 fontsize=12, fontweight="bold")
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/lipid_hbond_map{suf}.{e}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/lipid_hbond_map{suf}.png")

    # ---------- Part B: COOPERATIVITY on single lipids ----------
    # molecule key = (frame, lres, lrn) so lipid resids can't collide across types
    print(f"\n=== basic/acidic COOPERATIVITY on single {ltag} lipids (POST) ===")
    coop_rows, pair_counter = [], {}
    for name, df in dfs.items():
        d = df[(df.lrn.isin(lips)) & (df.pclass != "other")]
        if not len(d):
            print(f"  {name}: no {ltag} basic/acidic bonds"); continue
        # per (frame, lipid molecule): has basic? has acidic?
        agg = d.groupby(["frame", "lres", "lrn"])["pclass"].agg(
            has_basic=lambda s: (s == "basic").any(),
            has_acidic=lambda s: (s == "acidic").any())
        # also engaged lipids with ONLY other residues (for n00)
        all_eng = df[df.lrn.isin(lips)].groupby(["frame", "lres", "lrn"]).size().index
        agg = agg.reindex(all_eng, fill_value=False)
        b = agg["has_basic"].values; a = agg["has_acidic"].values
        n11 = int((b & a).sum()); n10 = int((b & ~a).sum())
        n01 = int((~b & a).sum()); n00 = int((~b & ~a).sum())
        N = n11 + n10 + n01 + n00
        pB = (n11 + n10) / N; pA = (n11 + n01) / N; pAB = n11 / N
        lift = pAB / (pA * pB) if pA * pB > 0 else np.nan
        OR = (n11 * n00) / (n10 * n01) if n10 * n01 > 0 else np.nan
        pA_givenB = n11 / (n11 + n10) if (n11 + n10) else np.nan
        pA_givennotB = n01 / (n01 + n00) if (n01 + n00) else np.nan
        print(f"\n  -- RhoA-{name} ({N} engaged {ltag}-frames) --")
        print(f"     P(basic)={pB:.3f}  P(acidic)={pA:.3f}  P(both)={pAB:.3f}")
        print(f"     lift P(both)/P(b)P(a) = {lift:.2f}   odds ratio = {OR:.2f}")
        print(f"     P(acidic | basic)={pA_givenB:.3f}   vs   P(acidic | NO basic)={pA_givennotB:.3f}"
              f"   -> {'ENRICHED' if pA_givenB>pA_givennotB else 'not enriched'} "
              f"({pA_givenB/pA_givennotB:.2f}x)" if pA_givennotB else "")
        coop_rows.append({"sim": name, "N": N, "P_basic": pB, "P_acidic": pA, "P_both": pAB,
                          "lift": lift, "odds_ratio": OR,
                          "P_acidic_given_basic": pA_givenB, "P_acidic_given_nobasic": pA_givennotB})
        # co-binding residue pairs on the same lipid molecule
        both = agg[(agg.has_basic) & (agg.has_acidic)].index
        for fr, lr, ln in both:
            sub = d[(d.frame == fr) & (d.lres == lr) & (d.lrn == ln)]
            basics = sub[sub.pclass == "basic"]["pres"].unique()
            acids = sub[sub.pclass == "acidic"]["pres"].unique()
            for bb in basics:
                for aa in acids:
                    pair_counter[(name, int(bb), int(aa))] = pair_counter.get((name, int(bb), int(aa)), 0) + 1
    pd.DataFrame(coop_rows).to_csv(f"{OUT}/lipid_cobinding_stats{suf}.csv", index=False)

    # ---------- Part C: per-pair CONDITIONAL enhancement (the real test) ----------
    # For a specific (basic b, acidic a): over engaged lipid-frames, does b being
    # bound to a lipid make a MORE likely on that SAME lipid?
    #   enhancement = P(a | b) / P(a | not b)   ( >1 = cooperative for this pair )
    print(f"\n=== per-pair conditional enhancement on the same {ltag} lipid "
          f"(P(acidic|basic)/P(acidic|no basic); >1 = cooperative) ===")
    top_pairs = sorted(pair_counter.items(), key=lambda kv: -kv[1])
    seen, enh_rows = set(), []
    print(f"  {'sim':>6}  {'basic+acidic':>14}  {'P(a|b)':>7} {'P(a|¬b)':>8} {'enhance':>8}  {'co-frames':>9}")
    for (name, bb, aa), cnt in top_pairs:
        if cnt < 25 or (name, bb, aa) in seen:      # compute floor = 25; 50-subset written too
            continue
        seen.add((name, bb, aa))
        df = dfs[name]
        pops = df[df.lrn.isin(lips)]
        U = set(map(tuple, pops[["frame", "lres", "lrn"]].drop_duplicates().values))
        Sb = set(map(tuple, pops[pops.pres == bb][["frame", "lres", "lrn"]].drop_duplicates().values))
        Sa = set(map(tuple, pops[pops.pres == aa][["frame", "lres", "lrn"]].drop_duplicates().values))
        if not Sb or len(U) == len(Sb):
            continue
        pAgB = len(Sa & Sb) / len(Sb)
        pAgnotB = (len(Sa) - len(Sa & Sb)) / (len(U) - len(Sb))
        enh = pAgB / pAgnotB if pAgnotB > 0 else np.inf

        def lab(r):
            m = df[df.pres == r]
            return C.res_label(m["prn"].iloc[0], r) if len(m) else str(r)
        enh_rows.append({"sim": name, "basic": lab(bb), "acidic": lab(aa),
                         "P_a_given_b": pAgB, "P_a_given_notb": pAgnotB,
                         "enhancement": enh, "co_frames": cnt})
        if len(seen) <= 16:
            print(f"  {name:>6}  {lab(bb)+'+'+lab(aa):>14}  {pAgB:>7.3f} {pAgnotB:>8.3f} "
                  f"{enh:>8.2f}  {cnt:>9}")
    enh = pd.DataFrame(enh_rows, columns=["sim", "basic", "acidic",
                       "P_a_given_b", "P_a_given_notb", "enhancement", "co_frames"])
    enh.to_csv(f"{OUT}/lipid_cobinding_enhancement{suf}_min25.csv", index=False)
    enh[enh.co_frames >= 50].to_csv(f"{OUT}/lipid_cobinding_enhancement{suf}_min50.csv", index=False)
    if not suf:
        enh.to_csv(f"{OUT}/lipid_cobinding_enhancement.csv", index=False)  # back-compat = POPS min25

    if pair_counter:
        print(f"\n=== top co-binding basic+acidic residue pairs on the same {ltag} (frame-counts) ===")
        pr = sorted(pair_counter.items(), key=lambda kv: -kv[1])[:20]
        def lab(r):
            for df in dfs.values():
                m = df[df.pres == r]
                if len(m):
                    return C.res_label(m["prn"].iloc[0], r)
            return str(r)
        print(f"  {'sim':>6}  {'basic':>6} + {'acidic':>6}   frames")
        for (name, bb, aa), c in pr:
            print(f"  {name:>6}  {lab(bb):>6} + {lab(aa):>6}   {c}")
        pd.DataFrame([{"sim": n, "basic": lab(b), "acidic": lab(a), "frames": c}
                      for (n, b, a), c in pair_counter.items()]).to_csv(
            f"{OUT}/lipid_cobinding_pairs{suf}.csv", index=False)
    print(f"\nwrote {OUT}/lipid_cobinding_stats{suf}.csv and _pairs{suf}.csv")


if __name__ == "__main__":
    main()
