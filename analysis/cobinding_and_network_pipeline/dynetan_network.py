#!/usr/bin/env python3
"""dynetan (Generalized-correlation Dynamical Network Analysis, Melo/Luthey-Schulten) on the
full-res post excerpts — MI network + communities + N-patch→Switch I optimal path, per sim,
with outputs ready for VMD/PyMOL structure rendering.

For each sim: generalized-correlation (KSG mutual information) network over residue nodes,
contacts 4.5 Å / 75 %, Louvain-style communities, node betweenness, and the optimal
N-patch(5-7)→Switch I(28-42) path. Writes, per sim, into transition_states/dynetan/:
  <sim>_nodes.csv     resid, community, betweenness   (community/betweenness -> B-factor for coloring)
  <sim>.pdb           aligned protein frame (cartoon target for rendering)
  <sim>_path.txt      N-patch→Switch I optimal path (resids)
  <sim>_edges.csv     network edges (resid_i, resid_j, gen-correlation)  for drawing the network

  python dynetan_network.py [--sims 2A,3A,GTP-2,GTP-1,GDP] [--sampled 4000] [--contact-stride 10] [--ncores 6]
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse
import numpy as np
import pandas as pd
import networkx as nx
import MDAnalysis as mda
from dynetan.proctraj import DNAproc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "transition_states", "dynetan")
NPATCH = (5, 7); SWITCH1 = (28, 42)
ORDER = ["GTP-2", "2A", "3A", "GTP-1", "GDP"]


def resolve_traj(name):
    """Locate <name>_last5us.dcd. Checks $TRAJ_DIR first (set it on the GPU box to
    point at wherever the excerpts live, no need to move the big files), then the bundle."""
    cands = []
    if os.environ.get("TRAJ_DIR"):
        cands.append(os.path.join(os.environ["TRAJ_DIR"], f"{name}_last5us.dcd"))
    cands.append(os.path.join(HERE, f"{name}_last5us.dcd"))
    for p in cands:
        if os.path.exists(p):
            return p
    return None


def process(name, sampled, cstride, ncores, debug=False, corr_frames=10000):
    psf = os.path.join(HERE, "topologies", f"{name}.psf")
    dcd = resolve_traj(name)
    dnap = DNAproc(notebookMode=False)
    dnap.setNumWinds(1)
    dnap.setNumSampledFrames(sampled)
    dnap.loadSystem(psf, dcd)
    # NOTE: this dynetan build's calcCor ignores numSampledFrames and runs the KSG MI on the
    # FULL in-memory trajectory (win_len = n_frames/numWinds). MI cost is ~O(N_frames^2), so the
    # 50k-frame 5 us excerpts would take ~weeks. Subsample the trajectory to ~corr_frames here so
    # the correlation runs on a converged-but-tractable set (~0.5 ns spacing over 5 us).
    nfr = dnap.workU.trajectory.n_frames
    fstride = max(1, round(nfr / corr_frames)) if corr_frames else 1
    if fstride > 1:
        dnap.workU.transfer_to_memory(step=fstride, verbose=False)
        print(f"  [{name}] subsampled {nfr} -> {dnap.workU.trajectory.n_frames} frames "
              f"for correlation (step {fstride})", flush=True)
    dnap.setSegIDs(["PROA"])
    dnap.selectSystem(withSolvent=True, inputSelStr="protein and resid 4-189")
    dnap.setCutoffDist(4.5)
    dnap.setContactPersistence(0.75)
    dnap.prepareNetwork(autocomp_groups=True)
    dnap.alignTraj(inMemory=True)
    dnap.findContacts(stride=cstride, n_cores=ncores)
    dnap.filterContacts(notSameRes=True, notConsecutiveRes=False, removeIsolatedNodes=True)
    dnap.calcCor(ncores=ncores, corr_type="gencor")
    dnap.calcGraphInfo()
    dnap.calcCommunities()
    try:
        dnap.calcOptPaths(ncores=ncores)
        dnap.calcBetween(ncores=ncores)
    except Exception as e:
        print(f"  [{name}] optPaths/between warn: {e}", flush=True)

    G = dnap.nxGraphs[0]
    resids = np.array(dnap.nodesAtmSel.resids)
    if debug:
        print("DNAproc attrs:", [a for a in dir(dnap) if not a.startswith('_')][:60], flush=True)
        n0 = list(G.nodes(data=True))[0]; e0 = list(G.edges(data=True))[0]
        print("node attr sample:", n0, flush=True)
        print("edge attr sample:", e0, flush=True)

    # community per node (dynetan stores 'modularity' on nodes)
    comm = {}
    for n, d in G.nodes(data=True):
        comm[n] = d.get("modularity", d.get("community", -1))
    # weighted betweenness from the graph (edge 'weight' = gen-correlation -> distance = -ln w)
    Gd = nx.Graph()
    for i, j, d in G.edges(data=True):
        w = d.get("weight", d.get("corr", None))
        if w is None:
            # fall back to correlation matrix
            w = dnap.corrMatAll[0][i, j] if hasattr(dnap, "corrMatAll") else 0.1
        w = float(max(w, 1e-6))
        Gd.add_edge(i, j, dist=-np.log(w), corr=w)
    btw = nx.betweenness_centrality(Gd, weight="dist", normalized=True)

    os.makedirs(OUT, exist_ok=True)
    nodes_df = pd.DataFrame({"node": list(range(len(resids))), "resid": resids,
                             "community": [comm.get(i, -1) for i in range(len(resids))],
                             "betweenness": [btw.get(i, 0.0) for i in range(len(resids))]})
    nodes_df.to_csv(f"{OUT}/{name}_nodes.csv", index=False)
    edges = [(int(resids[i]), int(resids[j]), float(Gd[i][j]["corr"])) for i, j in Gd.edges()]
    pd.DataFrame(edges, columns=["resid_i", "resid_j", "gencorr"]).to_csv(f"{OUT}/{name}_edges.csv", index=False)

    # optimal N-patch -> Switch I path (min -ln corr), from the same graph
    r2n = {int(r): i for i, r in enumerate(resids)}
    src = [r2n[r] for r in range(NPATCH[0], NPATCH[1] + 1) if r in r2n]
    snk = [r2n[r] for r in range(SWITCH1[0], SWITCH1[1] + 1) if r in r2n]
    best = None
    for s in src:
        try:
            L = nx.single_source_dijkstra_path_length(Gd, s, weight="dist")
        except Exception:
            continue
        for t in snk:
            if t in L and (best is None or L[t] < best[0]):
                best = (L[t], s, t)
    path_res = []
    if best:
        L, s, t = best
        path = nx.dijkstra_path(Gd, s, t, weight="dist")
        path_res = [int(resids[n]) for n in path]
        with open(f"{OUT}/{name}_path.txt", "w") as fh:
            fh.write(f"length={L:.3f}\npath={'-'.join(map(str, path_res))}\n")

    # aligned protein frame -> PDB (cartoon target); write community & betweenness into B-factors of two copies
    u = dnap.getU()
    if not hasattr(u.atoms, "tempfactors"):     # PSF-loaded Universe has no B-factor attr
        u.add_TopologyAttr("tempfactors")
    prot = u.select_atoms("protein and resid 4-189")
    u.trajectory[0]
    # community B-factor
    bfac_comm = np.zeros(len(prot))
    bfac_btw = np.zeros(len(prot))
    r2comm = dict(zip(resids, [comm.get(i, -1) for i in range(len(resids))]))
    r2btw = dict(zip(resids, [btw.get(i, 0.0) for i in range(len(resids))]))
    for k, a in enumerate(prot):
        bfac_comm[k] = r2comm.get(a.resid, -1)
        bfac_btw[k] = r2btw.get(a.resid, 0.0) * 1000  # scale for visibility
    try:
        prot.tempfactors = bfac_comm
        prot.write(f"{OUT}/{name}_comm.pdb")
        prot.tempfactors = bfac_btw
        prot.write(f"{OUT}/{name}_btw.pdb")
    except Exception as e:
        print(f"  [{name}] PDB write skipped ({e}) — CSVs still saved", flush=True)

    ncomm = len(set(comm.values()))
    print(f"  [{name}] nodes={len(resids)} edges={Gd.number_of_edges()} communities={ncomm} "
          f"path({len(path_res)}): {'-'.join(map(str, path_res))}", flush=True)
    return nodes_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="2A,3A,GTP-2,GTP-1,GDP")
    ap.add_argument("--sampled", type=int, default=4000)
    ap.add_argument("--contact-stride", type=int, default=10)
    ap.add_argument("--ncores", type=int, default=6)
    ap.add_argument("--corr-frames", type=int, default=10000,
                    help="subsample the trajectory to ~this many frames for the KSG MI (0 = use all)")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    for name in [s for s in ORDER if s in args.sims.split(",")]:
        if resolve_traj(name) is None:
            print(f"!! no trajectory for {name} — skipping"); continue
        process(name, args.sampled, args.contact_stride, args.ncores, debug=args.debug,
                corr_frames=args.corr_frames)
    print("DONE")


if __name__ == "__main__":
    main()
