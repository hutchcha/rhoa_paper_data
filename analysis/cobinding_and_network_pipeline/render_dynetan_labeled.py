#!/usr/bin/env python3
"""Labeled dynetan network renders: render each sim's structure SEPARATELY (aligned to a common
view) with PyMOL, then compose a labeled grid (sim titles + shared colorbar) with matplotlib —
so you can actually tell which structure is which.

betweenness: cartoon putty, spectrum by node betweenness (same scale across sims), N-patch->Switch I
             optimal path drawn as black spheres.
communities: cartoon colored by dynamical community (per-sim modules; colors are per-sim).

Reads transition_states/dynetan/<sim>_{btw,comm}.pdb + <sim>_nodes.csv + <sim>_path.txt.

  python render_dynetan_labeled.py [--mode betweenness|communities]
"""
import os, argparse, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

plt.rcParams["font.family"] = "Arial"
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "transition_states", "dynetan")
ORDER = ["GDP", "GTP-1", "GTP-2", "2A", "3A"]
REF = "GTP-2"
COMM_PALETTE = ["marine", "orange", "forest", "purpleblue", "yellow", "salmon", "cyan",
                "deepteal", "wheat", "slate", "firebrick", "limegreen", "violet", "lightpink"]


def present():
    return [s for s in ORDER if os.path.exists(f"{D}/{s}_btw.pdb")]


def path_resis(sim):
    p = f"{D}/{sim}_path.txt"
    if os.path.exists(p):
        t = open(p).read()
        if "path=" in t:
            return t.split("path=")[1].strip().split("-")
    return []


def render_pngs(sims, mode):
    """One PyMOL session: load all, align to REF, set common view, render each sim to its own PNG."""
    vmax = max(pd.read_csv(f"{D}/{s}_nodes.csv")["betweenness"].max() for s in sims) * 1000
    L = ["bg_color white", "set ray_opaque_background, 0", "set ray_shadows, 0",
         "set cartoon_putty_scale, 2", "set cartoon_putty_radius, 0.4"]
    suff = "btw" if mode == "betweenness" else "comm"
    for s in sims:
        L.append(f"load {D}/{s}_{suff}.pdb, {s}")
    for s in sims:
        if s != REF:
            L.append(f"align {s}, {REF}")
    L += ["hide everything", "show cartoon"]
    for s in sims:
        if mode == "betweenness":
            L.append(f"cartoon putty, {s}")
            L.append(f"spectrum b, blue_white_red, {s}, minimum=0, maximum={vmax:.2f}")
            pr = path_resis(s)
            if pr:
                sel = "+".join(pr)
                L.append(f"show spheres, {s} and name CA and resi {sel}")
                L.append(f"set sphere_scale, 0.6, {s} and name CA and resi {sel}")
                L.append(f"color black, {s} and name CA and resi {sel}")
        else:
            nd = pd.read_csv(f"{D}/{s}_nodes.csv")
            for k, c in enumerate(sorted(nd["community"].unique())):
                col = COMM_PALETTE[k % len(COMM_PALETTE)]
                L.append(f"color {col}, {s} and b>{c-0.5} and b<{c+0.5}")
    L += [f"orient {REF}", "set ray_trace_mode, 1"]
    outs = {}
    for s in sims:
        png = f"{D}/_{s}_{suff}_render.png"
        L += ["disable all", f"enable {s}", "ray 1000, 1000", f"png {png}, dpi=150"]
        outs[s] = png
    pml = f"{D}/_render_labeled_{mode}.pml"
    open(pml, "w").write("\n".join(L))
    subprocess.run(["pymol", "-cq", pml], capture_output=True, text=True)
    return outs, vmax


def compose(sims, outs, mode, vmax):
    fig, axes = plt.subplots(1, len(sims), figsize=(3.2 * len(sims), 3.6), constrained_layout=True)
    for ax, s in zip(np.atleast_1d(axes), sims):
        if os.path.exists(outs[s]):
            ax.imshow(plt.imread(outs[s]))
        ax.set_title(f"RhoA-{s}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                     fontsize=13, fontweight="bold")
        p = path_resis(s)
        if mode == "betweenness" and p:
            ax.set_xlabel("N-patch→SwI path: " + "→".join(p), fontsize=8, color="0.3")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    if mode == "betweenness":
        sm = cm.ScalarMappable(norm=colors.Normalize(0, vmax / 1000), cmap="bwr")
        cb = fig.colorbar(sm, ax=list(np.atleast_1d(axes)), shrink=0.6, location="bottom", pad=0.02)
        cb.set_label("node betweenness centrality (communication bottleneck)", fontsize=10)
    title = ("Allosteric communication bottlenecks (dynetan MI network); black = N-patch→Switch I path"
             if mode == "betweenness" else "Dynamical communities (per-sim modules)")
    fig.suptitle(title, fontsize=13, fontweight="bold")
    out = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "panels", f"dynetan_{mode}_labeled.png"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for e in ("png", "svg"):
        fig.savefig(out.replace(".png", f".{e}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="betweenness", choices=["betweenness", "communities"])
    args = ap.parse_args()
    sims = present()
    if not sims:
        print("no dynetan PDBs found"); return
    print("rendering:", sims)
    outs, vmax = render_pngs(sims, args.mode)
    compose(sims, outs, args.mode, vmax)


if __name__ == "__main__":
    main()
