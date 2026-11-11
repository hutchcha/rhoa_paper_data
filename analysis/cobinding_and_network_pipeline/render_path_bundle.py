#!/usr/bin/env python3
"""Suboptimal-path BUNDLE (WISP-style): the N-patch -> Switch I allosteric channel drawn on the
structure as the optimal path plus the top-k near-optimal suboptimal paths, with edge/node width =
how many paths use that connection/residue. Shows the communication CHANNEL (and its width /
degeneracy), not the whole network.

Step 1 (conda env): compute the bundle per sim from <sim>_edges.csv + <sim>_path.txt (networkx),
  write <sim>_bundle_edges.csv / <sim>_bundle_nodes.csv.
Step 2 (OUTSIDE conda): /usr/bin/pymol -cq transition_states/dynetan/_render_bundle.py
Step 3 (conda): python render_path_bundle.py --compose-only

  python render_path_bundle.py [--k 25]
"""
import os, argparse, subprocess
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "transition_states", "dynetan")
OUT = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "panels"))
ORDER = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
REF = "GTP-2"
REGION_RGB = {"N-patch": (0.55,0.34,0.29), "SwI": (0.79,0.64,0.15), "interSw": (0.34,0.65,0.60),
              "SwII": (0.51,0.45,0.70), "HVR": (0.30,0.45,0.69), "other": (0.60,0.63,0.65)}


def compute_bundle(sim, k):
    import networkx as nx
    import pandas as pd
    e = pd.read_csv(f"{D}/{sim}_edges.csv")
    G = nx.Graph()
    for _, r in e.iterrows():
        w = max(float(r["gencorr"]), 1e-6)
        G.add_edge(int(r["resid_i"]), int(r["resid_j"]), dist=-np.log(w))
    t = open(f"{D}/{sim}_path.txt").read()
    p0 = [int(x) for x in t.split("path=")[1].strip().split("-")]
    src, snk = p0[0], p0[-1]
    ecount = defaultdict(int); ncount = defaultdict(int); npath = 0
    for p in nx.shortest_simple_paths(G, src, snk, weight="dist"):
        for rr in p: ncount[rr] += 1
        for a, b in zip(p, p[1:]): ecount[tuple(sorted((a, b)))] += 1
        npath += 1
        if npath >= k: break
    with open(f"{D}/{sim}_bundle_edges.csv", "w") as fh:
        fh.write("i,j,count\n")
        for (a, b), c in ecount.items(): fh.write(f"{a},{b},{c}\n")
    with open(f"{D}/{sim}_bundle_nodes.csv", "w") as fh:
        fh.write("resid,count\n")
        for rr, c in ncount.items(): fh.write(f"{rr},{c}\n")
    return npath, src, snk


def region_of(r):
    if 5 <= r <= 7: return "N-patch"
    if 28 <= r <= 42: return "SwI"
    if 43 <= r <= 55: return "interSw"
    if 56 <= r <= 78: return "SwII"
    if 180 <= r <= 190: return "HVR"
    return "other"


def pml_script(sims, k):
    body = f'''
from pymol import cmd
from pymol.cgo import CYLINDER, COLOR, SPHERE
import csv
D = {D!r}; REF = {REF!r}; SIMS = {sims!r}; K = {k}
REGION_RGB = {REGION_RGB!r}
def region_of(r):
    if 5<=r<=7: return "N-patch"
    if 28<=r<=42: return "SwI"
    if 43<=r<=55: return "interSw"
    if 56<=r<=78: return "SwII"
    if 180<=r<=190: return "HVR"
    return "other"
def safe(s): return "o"+s.replace("-","").replace(" ","")
cmd.bg_color("white"); cmd.set("auto_zoom",0); cmd.set("ray_opaque_background",0); cmd.set("ray_shadows",0)
cmd.set("cartoon_transparency",0.8); cmd.set("label_size",16); cmd.set("label_color","black")
for s in SIMS: cmd.load(f"{{D}}/{{s}}_comm.pdb", safe(s))
for s in SIMS:
    if s!=REF: cmd.align(safe(s), safe(REF))
cmd.hide("everything"); cmd.show("cartoon"); cmd.color("gray85")
for s in SIMS:
    o=safe(s)
    coords={{}}
    cmd.iterate_state(1, f"{{o}} and name CA and polymer", "coords[int(resi)]=(x,y,z)", space={{'coords':coords}})
    ec=[]
    with open(f"{{D}}/{{s}}_bundle_edges.csv") as fh:
        for row in csv.DictReader(fh):
            i,j,c=int(row['i']),int(row['j']),int(row['count'])
            if i in coords and j in coords:
                frac=c/K; rad=0.15+1.4*frac
                x1,y1,z1=coords[i]; x2,y2,z2=coords[j]
                ec+=[CYLINDER, x1,y1,z1, x2,y2,z2, rad, 1.0,0.5,0.0, 1.0,0.5,0.0]
    cmd.load_cgo(ec, f"be_{{o}}")
    nc=[]
    with open(f"{{D}}/{{s}}_bundle_nodes.csv") as fh:
        for row in csv.DictReader(fh):
            r,c=int(row['resid']),int(row['count'])
            if r in coords:
                frac=c/K; rad=0.4+1.1*frac
                cr,cg,cb=REGION_RGB[region_of(r)]
                x,y,z=coords[r]; nc+=[COLOR,cr,cg,cb, SPHERE,x,y,z,rad]
    cmd.load_cgo(nc, f"bn_{{o}}")
    for rr,txt in [(6,"N-patch"),(37,"T37"),(76,"D76/78")]:
        if rr in coords:
            cmd.label(f"{{o}} and name CA and resi {{rr}}", repr(txt))
cmd.orient(safe(REF)); cmd.zoom(safe(REF), 4)
for s in SIMS:
    o=safe(s); cmd.disable("all")
    cmd.enable(o); cmd.enable("be_"+o); cmd.enable("bn_"+o)
    cmd.ray(1000,1000); cmd.png(f"{{D}}/_{{s}}_bundle.png", dpi=150)
'''
    open(f"{D}/_render_bundle.py", "w").write(body)


def compose(sims):
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, len(sims), figsize=(3.3*len(sims), 3.9), constrained_layout=True)
    for ax, s in zip(np.atleast_1d(axes), sims):
        png = f"{D}/_{s}_bundle.png"
        if os.path.exists(png): ax.imshow(plt.imread(png))
        t = open(f"{D}/{s}_path.txt").read(); pth = t.split("path=")[1].strip() if "path=" in t else ""
        ax.set_title(f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                     fontsize=13, fontweight="bold")
        ax.set_xlabel("optimal: " + "→".join(pth.split("-")), fontsize=8.5, color="0.3")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
    h = [Line2D([0],[0],marker="o",color="w",markerfacecolor=REGION_RGB[r],markeredgecolor="k",markersize=10,label=r)
         for r in ["N-patch","SwI","interSw","SwII"]]
    h += [Line2D([0],[0],color=(1.0,0.5,0.0),lw=5,label="path bundle (width = # paths)")]
    fig.legend(handles=h, frameon=False, fontsize=10, ncol=5, loc="lower center", bbox_to_anchor=(0.5,-0.03))
    fig.suptitle("N-patch→Switch I communication channel (optimal + top suboptimal paths; "
                 "orange width = number of paths using that link)", fontsize=12, fontweight="bold")
    for e in ("png","svg"):
        fig.savefig(f"{OUT}/dynetan_path_bundle.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/dynetan_path_bundle.png")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--compose-only", action="store_true")
    args = ap.parse_args()
    sims = [s for s in ORDER if os.path.exists(f"{D}/{s}_path.txt")]
    if not args.compose_only:
        for s in sims:
            n, src, snk = compute_bundle(s, args.k)
            print(f"  [{s}] {n} paths, {src}->{snk}", flush=True)
        pml_script(sims, args.k)
        print("run: /usr/bin/pymol -cq transition_states/dynetan/_render_bundle.py  (outside conda), then --compose-only")
    else:
        compose(sims)


if __name__ == "__main__":
    main()
