#!/usr/bin/env python3
"""Coarse-grained COMMUNITY network on the structure (the canonical dynamical-network-analysis
image): cartoon colored by dynamical community; one big node per community at its centroid;
inter-community edges weighted by cross-community MI coupling. The communities containing the
N-patch, Switch I and Switch II are labeled — so "N-patch is its own module in WT, merges with the
interswitch/switch module in the mutants" reads directly.

Reads transition_states/dynetan/<sim>_nodes.csv (resid,community), <sim>_edges.csv, <sim>_comm.pdb.
Renders each sim (aligned) with PyMOL, then composes a labeled grid.

Run the PyMOL step OUTSIDE the conda env (its python shadows pymol):
  python -c "import render_community_network as R; R.pml_script(['GDP','GTP-1','GTP-2','2A','3A'])"
  /usr/bin/pymol -cq transition_states/dynetan/_render_comm.py
  python render_community_network.py --compose-only
"""
import os, argparse, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "transition_states", "dynetan")
OUT = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure7"))
SUPP = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "supplementary"))
ORDER = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
REF = "GTP-2"
PALETTE = [(0.12,0.47,0.71),(1.0,0.5,0.05),(0.17,0.63,0.17),(0.84,0.15,0.16),(0.58,0.40,0.74),
           (0.55,0.34,0.29),(0.89,0.47,0.76),(0.6,0.6,0.6),(0.74,0.74,0.13),(0.09,0.75,0.81),
           (0.68,0.78,0.91),(1.0,0.73,0.47),(0.6,0.85,0.55),(0.9,0.6,0.6)]


def pml_script(sims):
    # WT-matched community -> RGB per sim (identical palette to the sequence track)
    import pandas as pd
    from dynetan_community_colors import community_colors
    _, _, per_sim, _ = community_colors(D, sims)
    COMMCOL = {}
    for s in sims:
        nodes = pd.read_csv(f"{D}/{s}_nodes.csv")
        cc = {}
        for r, c in zip(nodes.resid, nodes.community):
            if int(r) in per_sim[s]:
                cc[int(c)] = [float(x) for x in per_sim[s][int(r)][:3]]
        COMMCOL[s] = cc
    body = f'''
from pymol import cmd
from pymol.cgo import CYLINDER, COLOR, SPHERE
import csv
from collections import defaultdict
D = {D!r}; REF = {REF!r}; SIMS = {sims!r}
COMMCOL = {COMMCOL!r}
GRAY = [0.6, 0.6, 0.6]

def safe(s): return "o" + s.replace("-", "").replace(" ", "")

cmd.bg_color("white"); cmd.set("auto_zoom", 0)
cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0)
cmd.set("cartoon_transparency", 0.2)
cmd.set("cartoon_loop_radius", 0.35)
cmd.set("label_size", 18); cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
for s in SIMS: cmd.load(f"{{D}}/{{s}}_comm.pdb", safe(s))
for s in SIMS:
    if s != REF: cmd.align(safe(s), safe(REF))
cmd.hide("everything"); cmd.show("cartoon")

for s in SIMS:
    o = safe(s)
    coords = {{}}
    cmd.iterate_state(1, f"{{o}} and name CA and polymer", "coords[int(resi)]=(x,y,z)", space={{'coords':coords}})
    comm = {{}}
    with open(f"{{D}}/{{s}}_nodes.csv") as fh:
        for row in csv.DictReader(fh):
            comm[int(row['resid'])] = int(row['community'])
    members = defaultdict(list)
    for r,c in comm.items():
        if r in coords: members[c].append(r)
    cent = {{}}; size = {{}}
    for c,rs in members.items():
        xs = [coords[r] for r in rs]
        cent[c] = (sum(p[0] for p in xs)/len(xs), sum(p[1] for p in xs)/len(xs), sum(p[2] for p in xs)/len(xs))
        size[c] = len(rs)
    # color cartoon by WT-matched community color
    for c,rs in members.items():
        col = COMMCOL[s].get(c, GRAY)
        cmd.set_color(f"cc_{{o}}_{{c}}", list(col))
        cmd.color(f"cc_{{o}}_{{c}}", f"{{o}} and resi " + "+".join(str(r) for r in rs))
    # inter-community coupling (sum gencorr of crossing edges)
    coup = defaultdict(float)
    with open(f"{{D}}/{{s}}_edges.csv") as fh:
        for row in csv.DictReader(fh):
            i,j,w = int(row['resid_i']), int(row['resid_j']), float(row['gencorr'])
            if i in comm and j in comm and comm[i] != comm[j]:
                coup[tuple(sorted((comm[i],comm[j])))] += w
    maxc = max(coup.values()) if coup else 1.0
    ce = []
    for (a,b),w in coup.items():
        if w < 0.15*maxc: continue               # drop weak inter-community couplings (declutter)
        if a in cent and b in cent:
            x1,y1,z1 = cent[a]; x2,y2,z2 = cent[b]
            rad = 0.15 + 0.9*(w/maxc)             # thinner edges, scaled by coupling strength
            ce += [CYLINDER, x1,y1,z1, x2,y2,z2, rad, 0.35,0.35,0.38, 0.35,0.35,0.38]
    cmd.load_cgo(ce, f"cedges_{{o}}")
    cn = []
    maxsz = max(size.values()) if size else 1
    for c,(x,y,z) in cent.items():
        col = COMMCOL[s].get(c, GRAY)
        rad = 0.9 + 2.0*(size[c]/maxsz)          # node radius linear in community size (0.9..2.9)
        cn += [COLOR, col[0],col[1],col[2], SPHERE, x,y,z, rad]
    cmd.load_cgo(cn, f"cnodes_{{o}}")
    # (region text labels removed -- the sequence track now serves as the colour/region key)

cmd.orient(safe(REF)); cmd.zoom(safe(REF), 4)
for s in SIMS:
    o = safe(s)
    cmd.disable("all")
    for pre in ("", "cedges_", "cnodes_"):
        cmd.enable(pre + o)
    cmd.ray(1000, 1000)
    cmd.png(f"{{D}}/_{{s}}_commnet.png", dpi=150)
'''
    path = f"{D}/_render_comm.py"
    open(path, "w").write(body)
    return path


def compose(sims, outdir, stem, ncols=None):
    ncols = ncols or len(sims); nrows = -(-len(sims) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 3.4 * nrows), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.01, h_pad=0.01, wspace=0.0, hspace=0.02)
    axes = np.atleast_1d(axes).ravel()
    for ax in axes[len(sims):]:
        ax.axis("off")
    for ax, s in zip(axes, sims):
        png = f"{D}/_{s}_commnet.png"
        if os.path.exists(png):
            ax.imshow(plt.imread(png))
        ax.set_title(f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                     fontsize=14, fontweight="bold", pad=2)   # community count -> caption
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    os.makedirs(outdir, exist_ok=True)
    for e in ("png", "svg"):
        fig.savefig(f"{outdir}/{stem}.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outdir}/{stem}.png ({nrows}x{ncols})")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--compose-only", action="store_true")
    args = ap.parse_args()
    have = [s for s in ORDER if os.path.exists(f"{D}/{s}_comm.pdb")]
    if not args.compose_only:
        pml_script(have)
        subprocess.run(["/usr/bin/pymol", "-cq", f"{D}/_render_comm.py"],
                       stdout=open(f"{D}/_render_comm.log", "w"), stderr=subprocess.STDOUT)
    compose(["GTP-2", "3A"], OUT, "dynetan_community_network")                    # main
    compose(["GDP", "GTP-1"], SUPP, "supp_dynetan_community_network")             # supplement


if __name__ == "__main__":
    main()
