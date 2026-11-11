#!/usr/bin/env python3
"""Build an EDITABLE PyMOL session (.pse) of the WT co-binding interaction, using the author's
lighting script (paper_figures/figure6/'set lights2.pml'). Objects are split for easy editing:
  prot      -- protein cartoon, faint
  bridge    -- the shared POPS headgroup, sticks
  membrane  -- all other lipids, transparent surface (full membrane context)
  lig       -- GTP + Mg
Interacting residues R5/K6/K7 + D78 shown as sticks, coloured with util.cbac (cyan carbons).
Basic-amine -> POPS phosphate and D78-carboxylate -> POPS serine-amine polar contacts drawn as
dashes (yellow, per the lighting script). NO text labels (author adds them in Illustrator).

Anaconda pymol only (see render_cobinding_struct.py notes); emitted script self-reports errors.

  python build_cobinding_session.py           # writes .pml, runs anaconda pymol, saves .pse
"""
import os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
D = f"{HERE}/transition_states/cobinding_struct"
LIGHTS = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure6", "set lights2.pml"))
HEAD = ("resn POPS and name N+HN1+HN2+HN3+C11+C12+C13+O13A+O13B+P+O11+O12+O13+O14+"
        "C1+C2+C3+HA+HB+HS+HX+HY+O21+C21+O22+O31+C31+O32")


def main():
    info = {}
    for ln in open(f"{D}/chosen.txt"):
        k, _, v = ln.strip().partition("=")
        info[k] = v
    fr = int(info["frame"])
    pdb = f"{D}/GTP-2_cobind_f{fr}_membrane.pdb"
    pse = f"{D}/GTP-2_cobinding_session.pse"

    body = f'''from pymol import cmd, util

def run():
    cmd.set("auto_zoom", 0)
    # author lighting/palette (has a space in the filename -> feed line by line)
    for ln in open({LIGHTS!r}):
        s = ln.split("#", 1)[0].strip()
        if s:
            cmd.do(s)

    cmd.load({pdb!r}, "frame")
    cmd.create("prot", "frame and polymer")
    cmd.create("bridge", "frame and resn POPS and resi 323")
    cmd.create("membrane", "frame and (resn POPC+POPE+POPS+PLA18+PSM) and not (resn POPS and resi 323)")
    cmd.create("lig", "frame and resn GTP+MG")
    cmd.delete("frame")

    cmd.hide("everything")
    # protein: faint cartoon
    cmd.show("cartoon", "prot"); cmd.color("gray80", "prot")
    cmd.set("cartoon_transparency", 0.5, "prot")
    cmd.set("cartoon_side_chain_helper", 1)
    # interacting residues as sticks, cyan carbons (util.cbac)
    cmd.select("inter", "prot and resi 5+6+7+78")
    cmd.show("sticks", "inter and (sidechain or name CA)")
    util.cbac("inter")
    cmd.set("stick_radius", 0.22, "prot")
    # shared POPS: headgroup sticks (tails present but hidden -> easy to reveal)
    cmd.show("sticks", "bridge and ({HEAD})")
    util.cbaw("bridge"); cmd.color("wheat", "bridge and elem C")
    cmd.set("stick_radius", 0.22, "bridge")
    # ligand
    cmd.show("sticks", "lig and resn GTP"); cmd.color("gray70", "lig and elem C")
    cmd.show("spheres", "lig and resn MG"); cmd.color("green", "lig and resn MG")
    # membrane context: transparent surface
    cmd.show("surface", "membrane"); cmd.color("mem_bg", "membrane")
    cmd.set("transparency", 0.6, "membrane")
    cmd.set("surface_quality", 1)
    cmd.hide("everything", "hydro")

    # polar contacts (dashes yellow via lights script); labels hidden
    cmd.distance("d_bas", "(inter and resi 5+6+7 and name NZ+NH1+NH2+NE)",
                 "bridge and name O11+O12+O13+O14+P", 4.0, mode=0)
    cmd.distance("d_aci", "(inter and resi 78 and name OD1+OD2)",
                 "bridge and name N+HN1+HN2+HN3", 4.0, mode=0)
    cmd.hide("labels", "d_bas"); cmd.hide("labels", "d_aci")

    # alternate angle: tilt membrane out of plane so the bridge into the headgroup reads
    cmd.orient("inter or (bridge and ({HEAD}))")
    cmd.turn("x", -55); cmd.turn("y", 25)
    cmd.zoom("inter or (bridge and ({HEAD}))", 9.0)
    cmd.set("field_of_view", 25)

    cmd.save({pse!r})
    print("PMLOK wrote {pse}")


try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
'''
    path = f"{D}/_build_session.py"
    open(path, "w").write(body)
    subprocess.run(["pymol", "-cq", path],
                   stdout=open(f"{D}/_build_session.log", "w"), stderr=subprocess.STDOUT)
    print(open(f"{D}/_build_session.log").read()[-400:])
    print(f"session: {pse}")


if __name__ == "__main__":
    main()
