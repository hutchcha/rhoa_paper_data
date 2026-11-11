#!/usr/bin/env python3
"""Dynamical allosteric network analysis (Sethi/Luthey-Schulten) on the MI generalized-
correlation matrices — traces how a perturbation at the N-patch reaches Switch I.

Graph: node per residue (Cα). Edge between residues whose heavy atoms are within 4.5 Å for
>= 75 % of frames (per-sim contact map). Edge DISTANCE weight = -ln(r_MI) so strongly
correlated contacts are "short" (preferred communication). Then per sim:
  - betweenness centrality  -> communication bottleneck residues (allosteric hubs)
  - optimal + suboptimal N-patch(5-7) -> Switch I(28-42) paths -> the allosteric channel
  - greedy-modularity communities (weighted by r_MI) -> dynamically coupled modules
and compares them across sims (does the mutation lengthen/reroute the path, move the hubs,
or reassign the switches to different communities?).

Reads transition_states/dccm_mi/mi_<sim>.npy (+ mi_resids.npy) from dccm_mi.py, and the
trajectories (for contacts). Self-contained on the bundle config.

  python network_allostery.py [--sims 2A,3A,GTP-2,GTP-1,GDP] [--contact-stride 10] [--occ 0.75]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse
import numpy as np
import pandas as pd
import networkx as nx
import MDAnalysis as mda
from MDAnalysis.lib.distances import self_capped_distance
import stable_state_config as C
from dccm_mi import resolve_traj, OUT as MI_OUT, ALLCA, NPATCH, SWITCH1, SWITCH2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "transition_states", "network")
ORDER = ["GTP-2", "2A", "3A", "GTP-1", "GDP"]
CONTACT_CUT = 4.5
HL = [("N-patch", 5, 7), ("Switch I", 28, 42), ("interSw", 43, 55), ("Switch II", 56, 78),
      ("P-loop", 10, 17), ("Rho ins", 120, 137), ("α5", 164, 177), ("HVR", 180, 189)]


def contact_occupancy(sim, resids, stride):
    """Fraction of frames each residue pair has heavy atoms within CONTACT_CUT."""
    u = mda.Universe(sim["top"], resolve_traj(sim))
    heavy = u.select_atoms("protein and not name H* and resid 4-189")
    a2r = np.searchsorted(resids, heavy.resids)          # atom -> residue-index (resids sorted, contiguous)
    N = len(resids)
    occ = np.zeros((N, N))
    frames = np.arange(0, len(u.trajectory), stride)
    for f in frames:
        u.trajectory[int(f)]
        pairs = self_capped_distance(heavy.positions, CONTACT_CUT, box=u.dimensions,
                                     return_distances=False)
        ri, rj = a2r[pairs[:, 0]], a2r[pairs[:, 1]]
        m = ri != rj
        M = np.zeros((N, N), bool)
        M[ri[m], rj[m]] = True; M[rj[m], ri[m]] = True
        occ += M
    return occ / len(frames)


def build_graph(mi, occ, resids, occ_cut):
    N = len(resids)
    G = nx.Graph()
    for i in range(N):
        G.add_node(i, resid=int(resids[i]))
    iu = np.triu_indices(N, 1)
    for i, j in zip(*iu):
        if occ[i, j] >= occ_cut:
            r = max(mi[i, j], 1e-6)
            G.add_edge(i, j, dist=-np.log(r), corr=r)
    return G


def npatch_swI_paths(G, resids, ntol=0.5):
    r2i = {int(r): i for i, r in enumerate(resids)}
    src = [r2i[r] for r in range(NPATCH[0], NPATCH[1] + 1) if r in r2i]
    snk = [r2i[r] for r in range(SWITCH1[0], SWITCH1[1] + 1) if r in r2i]
    best = None
    for s in src:
        lengths = nx.single_source_dijkstra_path_length(G, s, weight="dist")
        for t in snk:
            if t in lengths and (best is None or lengths[t] < best[0]):
                best = (lengths[t], s, t)
    if best is None:
        return None
    L, s, t = best
    path = nx.dijkstra_path(G, s, t, weight="dist")
    # suboptimal degeneracy: count simple paths within ntol of optimal length
    nsub = 0
    for p in nx.shortest_simple_paths(G, s, t, weight="dist"):
        w = sum(G[p[k]][p[k+1]]["dist"] for k in range(len(p)-1))
        if w <= L * (1 + ntol):
            nsub += 1
            if nsub > 500:
                break
        else:
            break
    return dict(length=L, source=int(resids[s]), sink=int(resids[t]),
                path=[int(resids[n]) for n in path], n_suboptimal=nsub)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="2A,3A,GTP-2,GTP-1,GDP")
    ap.add_argument("--contact-stride", type=int, default=10)
    ap.add_argument("--occ", type=float, default=0.75)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    resids = np.load(f"{MI_OUT}/mi_resids.npy")

    want = [s for s in ORDER if s in args.sims.split(",")]
    betw = {}; paths = {}; comm = {}
    for name in want:
        mip = f"{MI_OUT}/mi_{name}.npy"
        if not os.path.exists(mip):
            print(f"!! no MI matrix for {name} (run dccm_mi.py) — skipping"); continue
        sim = [x for x in C.SIMS if x["name"] == name][0]
        if not os.path.exists(resolve_traj(sim)):
            print(f"!! no trajectory for {name} — skipping"); continue
        mi = np.load(mip)
        occ = contact_occupancy(sim, resids, args.contact_stride)
        G = build_graph(mi, occ, resids, args.occ)
        bc = nx.betweenness_centrality(G, weight="dist", normalized=True)
        betw[name] = np.array([bc.get(i, 0.0) for i in range(len(resids))])
        paths[name] = npatch_swI_paths(G, resids)
        cs = nx.community.greedy_modularity_communities(G, weight="corr")
        cassign = np.full(len(resids), -1)
        for ci, cset in enumerate(cs):
            for n in cset:
                cassign[n] = ci
        comm[name] = cassign
        p = paths[name]
        print(f"  [{name}] edges={G.number_of_edges()}  communities={len(cs)}  "
              f"N-patch→SwI: len={p['length']:.2f} via {p['source']}→{p['sink']} "
              f"({len(p['path'])} nodes, {p['n_suboptimal']} suboptimal)\n"
              f"          path: {p['path']}", flush=True)

    if not betw:
        print("nothing computed."); return
    # save
    pd.DataFrame({"resid": resids, **{n: betw[n] for n in betw}}).to_csv(f"{OUT}/betweenness.csv", index=False)
    pd.DataFrame([{"sim": n, **{k: (v if not isinstance(v, list) else "-".join(map(str, v)))
                                for k, v in paths[n].items()}} for n in paths]).to_csv(
        f"{OUT}/npatch_swI_paths.csv", index=False)
    pd.DataFrame({"resid": resids, **{n: comm[n] for n in comm}}).to_csv(f"{OUT}/communities.csv", index=False)

    plot(betw, paths, comm, resids)


def plot(betw, paths, comm, resids):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    COL = {"3A": "#4C78A8", "2A": "#F58518", "GTP-2": "#E45756", "GTP-1": "#72B7B2", "GDP": "#B279A2"}
    sims = [s for s in ORDER if s in betw]

    # betweenness profile
    fig, ax = plt.subplots(figsize=(15, 4.6))
    for name in sims:
        ax.plot(resids, betw[name], lw=1.4, color=COL.get(name), label=f"RhoA-{name}")
    for nm, lo, hi in HL:
        ax.axvspan(lo, hi, color="#ffe08a" if "Switch" in nm or nm == "interSw" else "#e8e8e8", alpha=0.4)
        ax.text((lo+hi)/2, ax.get_ylim()[1]*0.96, nm, ha="center", va="top", fontsize=8)
    ax.set_xlabel("residue", fontsize=12); ax.set_ylabel("betweenness centrality", fontsize=12)
    ax.set_title("Allosteric communication bottlenecks (edge-weighted betweenness)", fontsize=13, loc="left")
    ax.legend(frameon=False, ncol=5, fontsize=10); ax.spines[["right", "top"]].set_visible(False)
    fig.tight_layout()
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/betweenness.{e}", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # N-patch->SwI optimal path length + degeneracy
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    x = np.arange(len(sims))
    axes[0].bar(x, [paths[n]["length"] for n in sims], color=[COL[n] for n in sims], edgecolor="k", linewidth=0.4)
    axes[0].set_xticks(x); axes[0].set_xticklabels([f"RhoA-{n}" for n in sims], rotation=20, fontsize=10)
    axes[0].set_ylabel("optimal path length  Σ-ln r$_{MI}$", fontsize=10)
    axes[0].set_title("N-patch → Switch I communication distance", fontsize=11, loc="left")
    axes[1].bar(x, [paths[n]["n_suboptimal"] for n in sims], color=[COL[n] for n in sims], edgecolor="k", linewidth=0.4)
    axes[1].set_xticks(x); axes[1].set_xticklabels([f"RhoA-{n}" for n in sims], rotation=20, fontsize=10)
    axes[1].set_ylabel("# suboptimal paths (≤1.5× optimal)", fontsize=10)
    axes[1].set_title("channel degeneracy", fontsize=11, loc="left")
    for ax in axes:
        ax.spines[["right", "top"]].set_visible(False)
    fig.suptitle("N-patch → Switch I allosteric path", fontsize=13, fontweight="bold")
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/npatch_swI_path.{e}", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # communities as colored residue bands per sim
    fig, ax = plt.subplots(figsize=(15, 0.7 * len(sims) + 1.2))
    for row, name in enumerate(sims):
        ca = comm[name]
        ncom = ca.max() + 1
        colors = cm.tab20(np.linspace(0, 1, max(ncom, 1)))
        for i, r in enumerate(resids):
            ax.barh(row, 1, left=r, height=0.8, color=colors[ca[i] % 20], edgecolor="none")
    ax.set_yticks(range(len(sims))); ax.set_yticklabels([f"RhoA-{n}" for n in sims])
    for _, lo, hi in HL:
        ax.axvline(lo, color="k", lw=0.4, alpha=0.3); ax.axvline(hi, color="k", lw=0.4, alpha=0.3)
    ax.set_xlabel("residue", fontsize=12)
    ax.set_title("Dynamical communities (color = module; do the switches share the N-patch's module?)",
                 fontsize=12, loc="left")
    ax.set_xlim(resids.min(), resids.max())
    fig.tight_layout()
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/communities.{e}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/betweenness.png, npatch_swI_path.png, communities.png")


if __name__ == "__main__":
    main()
