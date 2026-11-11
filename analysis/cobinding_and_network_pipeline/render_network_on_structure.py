#!/usr/bin/env python3
"""PSN-style node+edge network drawn ON the protein structure, highlighting the N-patch -> Switch I
optimal path. Renders each sim (aligned to a common view) with PyMOL, then composes a labeled grid.

For each sim: a semi-transparent cartoon for context; network edges (dynetan contact+MI edges) as
thin gray cylinders between residue Cα; nodes as spheres colored by structural region; and the
N-patch->Switch I optimal path drawn boldly (thick orange cylinders + large path-node spheres).

Reads transition_states/dynetan/<sim>_edges.csv, <sim>_path.txt, <sim>_comm.pdb.

  python render_network_on_structure.py [--edge-min 0.45]
"""
import os, argparse, subprocess, csv
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
# region colors as (r,g,b) 0-1 for PyMOL
REGION_RGB = {"N-patch": (0.55, 0.34, 0.29), "SwI": (0.79, 0.64, 0.15), "interSw": (0.34, 0.65, 0.60),
              "SwII": (0.51, 0.45, 0.70), "HVR": (0.30, 0.45, 0.69), "other": (0.60, 0.63, 0.65)}


def present():
    return [s for s in ORDER if os.path.exists(f"{D}/{s}_comm.pdb")]


def pml_script(sims, edge_min):
    """PyMOL python script: load all, align, draw network+path CGO per sim, render each separately."""
    body = f'''
from pymol import cmd
from pymol.cgo import CYLINDER, COLOR, SPHERE
import csv
D = {D!r}
REF = {REF!r}
SIMS = {sims!r}
EDGE_MIN = {edge_min}
REGION_RGB = {REGION_RGB!r}

def region_of(r):
    if 5 <= r <= 7: return "N-patch"
    if 28 <= r <= 42: return "SwI"
    if 43 <= r <= 55: return "interSw"
    if 56 <= r <= 78: return "SwII"
    if 180 <= r <= 190: return "HVR"
    return "other"

def safe(s):
    return "o" + s.replace("-", "").replace(" ", "")     # hyphen-free object name (o2A, oGTP2, ...)

ROI = set(range(5,8)) | set(range(28,79)) | set(range(180,191))   # N-patch+SwI+interSw+SwII+HVR (drop scaffold)
KEY_LABELS = [(37,"T37"),(77,"D76/78"),(6,"N-patch")]
cmd.bg_color("white")
cmd.set("auto_zoom", 0)
cmd.set("ray_opaque_background", 0)
cmd.set("ray_shadows", 0)
cmd.set("cartoon_transparency", 0.82)
cmd.set("label_size", 15); cmd.set("label_color", "black")
for s in SIMS:
    cmd.load(f"{{D}}/{{s}}_comm.pdb", safe(s))
for s in SIMS:
    if s != REF:
        cmd.align(safe(s), safe(REF))
cmd.hide("everything"); cmd.show("cartoon"); cmd.color("skyblue")

for s in SIMS:
    o = safe(s)
    # Ca coords for this (aligned) object
    coords = {{}}
    cmd.iterate_state(1, f"{{o}} and name CA and polymer", "coords[int(resi)]=(x,y,z)", space={{'coords':coords}})
    # edges
    net = []
    with open(f"{{D}}/{{s}}_edges.csv") as fh:
        for row in csv.DictReader(fh):
            i, j, w = int(row['resid_i']), int(row['resid_j']), float(row['gencorr'])
            if w >= EDGE_MIN and i in ROI and j in ROI and i in coords and j in coords:
                x1,y1,z1 = coords[i]; x2,y2,z2 = coords[j]
                net += [CYLINDER, x1,y1,z1, x2,y2,z2, 0.10, 0.55,0.55,0.58, 0.55,0.55,0.58]
    cmd.load_cgo(net, f"net_{{o}}")
    # nodes (ROI only)
    nodes = []
    for r,(x,y,z) in coords.items():
        if r not in ROI: continue
        cr,cg,cb = REGION_RGB[region_of(r)]
        nodes += [COLOR, cr,cg,cb, SPHERE, x,y,z, 0.6]
    cmd.load_cgo(nodes, f"nodes_{{o}}")
    # path (bold)
    path = []
    with open(f"{{D}}/{{s}}_path.txt") as fh:
        t = fh.read()
    pres = [int(p) for p in t.split("path=")[1].strip().split("-")] if "path=" in t else []
    pcgo = []
    for k in range(len(pres)-1):
        a,b = pres[k], pres[k+1]
        if a in coords and b in coords:
            x1,y1,z1 = coords[a]; x2,y2,z2 = coords[b]
            pcgo += [CYLINDER, x1,y1,z1, x2,y2,z2, 0.35, 1.0,0.45,0.0, 1.0,0.45,0.0]
    for r in pres:
        if r in coords:
            x,y,z = coords[r]
            cr,cg,cb = REGION_RGB[region_of(r)]
            pcgo += [COLOR, cr,cg,cb, SPHERE, x,y,z, 1.05]
    cmd.load_cgo(pcgo, f"path_{{o}}")
    for rr,txt in KEY_LABELS:
        cmd.label(f"{{o}} and name CA and resi {{rr}}", repr(txt))

cmd.set("ray_trace_mode", 0)
cmd.orient(safe(REF)); cmd.zoom(safe(REF), 4)      # full-protein view, set AFTER all CGO loaded
for s in SIMS:
    o = safe(s)
    cmd.disable("all")
    cmd.enable(o); cmd.enable("net_"+o); cmd.enable("nodes_"+o); cmd.enable("path_"+o)
    cmd.ray(1000, 1000)
    cmd.png(f"{{D}}/_{{s}}_netstruct.png", dpi=150)
'''
    path = f"{D}/_render_netstruct.py"
    open(path, "w").write(body)
    return path


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--edge-min", type=float, default=0.45)
    args = ap.parse_args()
    sims = present()
    if not sims:
        print("no dynetan outputs"); return
    print("rendering network-on-structure:", sims)
    script = pml_script(sims, args.edge_min)
    with open(f"{D}/_render_netstruct.log", "w") as lg:      # capture_output breaks pymol ray; redirect instead
        subprocess.run(["pymol", "-cq", script], stdout=lg, stderr=lg)

    fig, axes = plt.subplots(1, len(sims), figsize=(3.3 * len(sims), 3.7), constrained_layout=True)
    for ax, s in zip(np.atleast_1d(axes), sims):
        png = f"{D}/_{s}_netstruct.png"
        if os.path.exists(png):
            ax.imshow(plt.imread(png))
        t = open(f"{D}/{s}_path.txt").read()
        pth = t.split("path=")[1].strip() if "path=" in t else ""
        ax.set_title(f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                     fontsize=13, fontweight="bold")
        ax.set_xlabel("N-patch→SwI: " + "→".join(pth.split("-")), fontsize=8, color="0.3")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.suptitle("Co-binding/MI network on the structure — N-patch→Switch I path in orange "
                 "(edges = MI-correlated contacts; spheres = residue nodes by region)",
                 fontsize=12, fontweight="bold")
    out = f"{OUT}/dynetan_network_on_structure.png"
    for e in ("png", "svg"):
        fig.savefig(out.replace(".png", f".{e}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
