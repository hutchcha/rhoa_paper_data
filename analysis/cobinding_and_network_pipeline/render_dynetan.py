#!/usr/bin/env python3
"""Render the dynetan networks on the protein structure (PyMOL, headless) into side-by-side
comparison images across sims.

Reads dynetan_network.py outputs (transition_states/dynetan/<sim>_{comm,btw}.pdb,
<sim>_nodes.csv, <sim>_path.txt). Produces, into the same dir:
  dynetan_communities.png   grid: each sim's cartoon colored by dynamical community (modules)
  dynetan_betweenness.png   grid: each sim's cartoon, putty+spectrum by betweenness (same scale)
                            with the N-patch→Switch I optimal path drawn as spheres/sticks

All sims are superposed on the first present sim so orientation is comparable.

  python render_dynetan.py [--sims 2A,3A,GTP-2,GTP-1,GDP]
"""
import os, argparse, subprocess
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "transition_states", "dynetan")
ORDER = ["GTP-2", "2A", "3A", "GTP-1", "GDP"]
PALETTE = ["marine", "orange", "forest", "purple", "yellow", "salmon", "cyan", "deepteal",
           "wheat", "slate", "firebrick", "limegreen", "violet", "lightblue", "gray70"]


def present(sims):
    return [s for s in ORDER if s in sims and os.path.exists(f"{D}/{s}_comm.pdb")]


def communities_pml(sims):
    lines = ["bg_color white", "set ray_opaque_background, 0", "set grid_mode, 1",
             "set cartoon_transparency, 0", "set ray_shadows, 0"]
    ref = sims[0]
    for s in sims:
        lines.append(f"load {D}/{s}_comm.pdb, {s}")
    for s in sims[1:]:
        lines.append(f"align {s}, {ref}")
    lines += ["hide everything", "show cartoon", "set cartoon_transparency, 0"]
    for s in sims:
        nd = pd.read_csv(f"{D}/{s}_nodes.csv")
        for k, c in enumerate(sorted(nd["community"].unique())):
            col = PALETTE[k % len(PALETTE)]
            lines.append(f"color {col}, {s} and b>{c-0.5} and b<{c+0.5}")
        lines.append(f"set_name {s}, {s}")
    lines += [f"orient {ref}", "set ray_trace_mode, 1", "ray 2400, 1400",
              f"png {D}/dynetan_communities.png, dpi=150"]
    return "\n".join(lines)


def betweenness_pml(sims):
    # common betweenness scale
    vmax = max(pd.read_csv(f"{D}/{s}_nodes.csv")["betweenness"].max() for s in sims) * 1000
    lines = ["bg_color white", "set ray_opaque_background, 0", "set grid_mode, 1",
             "set ray_shadows, 0", "set cartoon_putty_scale, 2",
             "set cartoon_putty_radius, 0.4"]
    ref = sims[0]
    for s in sims:
        lines.append(f"load {D}/{s}_btw.pdb, {s}")
    for s in sims[1:]:
        lines.append(f"align {s}, {ref}")
    lines += ["hide everything", "show cartoon"]
    for s in sims:
        lines.append(f"cartoon putty, {s}")
        lines.append(f"spectrum b, blue_white_red, {s}, minimum=0, maximum={vmax:.2f}")
        # draw optimal path as spheres, if present
        pf = f"{D}/{s}_path.txt"
        if os.path.exists(pf):
            txt = open(pf).read()
            if "path=" in txt:
                pres = txt.split("path=")[1].strip().split("-")
                sel = "+".join(pres)
                lines.append(f"show spheres, {s} and name CA and resi {sel}")
                lines.append(f"set sphere_scale, 0.6, {s} and name CA and resi {sel}")
                lines.append(f"color black, {s} and name CA and resi {sel}")
    lines += [f"orient {ref}", "ray 2400, 1400", f"png {D}/dynetan_betweenness.png, dpi=150"]
    return "\n".join(lines)


def run_pml(script, tag):
    path = f"{D}/_render_{tag}.pml"
    with open(path, "w") as fh:
        fh.write(script)
    r = subprocess.run(["pymol", "-cq", path], capture_output=True, text=True)
    print(r.stdout[-500:] if r.stdout else "", r.stderr[-500:] if r.stderr else "")
    print(f"wrote {D}/dynetan_{tag}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="2A,3A,GTP-2,GTP-1,GDP")
    args = ap.parse_args()
    sims = present(args.sims.split(","))
    if not sims:
        print("no dynetan PDB outputs found — run dynetan_network.py first"); return
    print("rendering sims:", sims)
    run_pml(communities_pml(sims), "communities")
    run_pml(betweenness_pml(sims), "betweenness")


if __name__ == "__main__":
    main()
