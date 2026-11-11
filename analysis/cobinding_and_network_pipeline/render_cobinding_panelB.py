#!/usr/bin/env python3
"""Render Fig 6 Panel C (co-binding structural snapshot) from the session (.pse): zoomed-in on the
R5/K6/K7 + D78 / POPS co-binding cluster, membrane drawn as a SMOOTH sand-colored surface (large
probe radius so headgroup detail reads), interacting residues (util.cbac) and POPS headgroup as
sticks with yellow polar-contact dashes, and residue + lipid labels.

  python render_cobinding_panelB.py    # runs anaconda pymol
-> paper_figures/figure6/fig_struct_cobinding.png
"""
import os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
D = f"{HERE}/transition_states/cobinding_struct"
PSE = f"{D}/GTP-2_cobinding_session.pse"
OUT = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure6", "fig_struct_cobinding.png"))
TAILS = ("C22+C23+C24+C25+C26+C27+C28+C29+C210+C211+C212+C213+C214+C215+C216+C217+C218+"
         "C32+C33+C34+C35+C36+C37+C38+C39+C310+C311+C312+C313+C314+C315+C316")
LABELS = {5: "R5", 6: "K6", 7: "K7", 78: "D78"}


def main():
    labcmds = "\n    ".join(
        f'cmd.label("prot and resi {r} and name CA", {(chr(34)+t+chr(34))!r})'
        for r, t in LABELS.items())
    body = f'''from pymol import cmd

def run():
    cmd.load({PSE!r})
    cmd.bg_color("white")
    cmd.set_color("mem_sand", [0.91, 0.87, 0.78])
    cmd.set("transparency_mode", 3)
    cmd.set("ray_interior_color", "gray90")
    cmd.set("backface_cull", 0)
    cmd.set("two_sided_lighting", 1)

    # smooth, sand membrane: big probe + forced surface recompute merges lipids into a slab
    cmd.color("mem_sand", "membrane")
    cmd.hide("surface", "membrane")
    cmd.set("surface_quality", 2)             # finer triangulation -> smoother slab
    cmd.set("solvent_radius", 6.0)            # large probe -> smooth continuous membrane
    cmd.set("surface_mode", 0)
    cmd.set("cavity_cull", 8)                 # drop small disconnected islands
    cmd.set("transparency", 0.5, "membrane")
    cmd.show("surface", "membrane")
    cmd.rebuild()

    # residue + lipid labels (small; author can reposition in Illustrator)
    cmd.set("label_size", 15)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_bg_color", "white"); cmd.set("label_bg_transparency", 0.3)
    {labcmds}
    cmd.label("bridge and name P", '"POPS"')

    # zoom in tight on the cluster (interacting sidechains + POPS headgroup, tails excluded)
    cmd.zoom("inter or (bridge and not (name {TAILS}))", 3.2)
    cmd.ray(1500, 1300)
    cmd.png({OUT!r}, dpi=200)
    print("PMLOK wrote {OUT}")


try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
'''
    path = f"{D}/_render_panelB.py"
    open(path, "w").write(body)
    subprocess.run(["pymol", "-cq", path],
                   stdout=open(f"{D}/_render_panelB.log", "w"), stderr=subprocess.STDOUT)
    print(open(f"{D}/_render_panelB.log").read()[-400:])


if __name__ == "__main__":
    main()
