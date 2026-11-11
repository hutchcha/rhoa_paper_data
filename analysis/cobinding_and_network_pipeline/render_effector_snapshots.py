"""Render one structural snapshot per effector for Figure 7: the crystal's Switch I/II (gray
reference target) with each sim's core-aligned representative switches overlaid, colored by sim
(same palette as the RMSD KDEs). Switch I is drawn bold, Switch II faint, so the eye reads the
Switch-I competence story (3A's Switch I sits closest to the gray crystal).

Snapshots come from extract_effector_snapshots.py (already superposed on the crystal core).
Uses the lab lighting script. Anaconda pymol only; emitted script self-reports errors.

  conda activate memb_analysis   # extract step
  python render_effector_snapshots.py     # writes the .pml then runs anaconda pymol
"""
import os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
SNAP = os.path.join(HERE, "transition_states", "effector_snapshots")
LIGHTS = os.path.abspath(os.path.join(HERE, "..", "paper_figures", "figure6", "set lights2.pml"))
# one render per sim GROUP so a panel only shows the structures for the sims it actually plots:
#   main -> GTP-2/3A (main Fig. 7A, S8);   refs -> GDP/GTP-1 (S9)
GROUPS = {"main": ["GTP-2", "3A"], "refs": ["GDP", "GTP-1"]}
# match fig_D KDE palette
SIM_COL = {"GDP": (0.70, 0.475, 0.635), "GTP-1": (0.447, 0.718, 0.698),
           "GTP-2": (0.894, 0.341, 0.339), "2A": (0.961, 0.522, 0.094), "3A": (0.298, 0.471, 0.706)}


def main():
    effs = [ln.split("\t") for ln in open(f"{SNAP}/effectors.txt").read().splitlines() if ln.strip()]
    body = f'''from pymol import cmd

GROUPS = {GROUPS!r}
SIM_COL = {SIM_COL!r}
EFFS = {effs!r}
SNAP = {SNAP!r}

def run():
    cmd.set("auto_zoom", 0)
    for ln in open({LIGHTS!r}):
        s = ln.split("#", 1)[0].strip()
        if s:
            try: cmd.do(s)
            except Exception: pass
    cmd.set("cartoon_loop_radius", 0.35)
    cmd.set("cartoon_transparency", 0.0)
    for name, col in SIM_COL.items():
        cmd.set_color("sc_" + name, list(col))

    SW = {{"SwI": "resi 28-42", "SwII": "resi 56-78"}}   # main = SwI only; SwII -> supplement
    for eff, pdb, ch in EFFS:
        for swname, swsel in SW.items():
            for group, gsims in GROUPS.items():            # one render per sim group
                cmd.delete("all")
                cmd.load(pdb, "xtal")
                cmd.hide("everything")
                rho = f"xtal and chain {{ch}} and polymer"
                cmd.show("cartoon", rho)                   # faint core context
                cmd.color("gray90", rho)
                cmd.set("cartoon_transparency", 0.82, rho)
                xsw = f"xtal and chain {{ch}} and {{swsel}}"  # gray crystal reference switch
                cmd.set("cartoon_transparency", 0.0, xsw)
                cmd.color("gray30", xsw)
                for name in gsims:                         # only the sims plotted in this panel
                    o = "o" + name.replace("-", "")
                    cmd.load(f"{{SNAP}}/{{eff}}_{{name}}.pdb", o)
                    cmd.hide("everything", o)
                    cmd.show("cartoon", f"{{o}} and {{swsel}}")
                    cmd.color("sc_" + name, f"{{o}} and {{swsel}}")
                cmd.orient(xsw)
                cmd.zoom(xsw, 3.5)
                cmd.ray(1500, 1350)
                out = f"{{SNAP}}/render_{{eff}}_{{swname}}_{{group}}.png"
                cmd.png(out, dpi=200)
                print("PMLOK", out)

    # separate 'how close does each get' overlay: the single closest (min-RMSD) Switch I frame
    # for GTP-2 and 3A on the gray crystal (not used in the KDE composites)
    for eff, pdb, ch in EFFS:
        swsel = "resi 28-42"
        cmd.delete("all")
        cmd.load(pdb, "xtal")
        cmd.hide("everything")
        rho = f"xtal and chain {{ch}} and polymer"
        cmd.show("cartoon", rho); cmd.color("gray90", rho)
        cmd.set("cartoon_transparency", 0.82, rho)
        xsw = f"xtal and chain {{ch}} and {{swsel}}"
        cmd.set("cartoon_transparency", 0.0, xsw); cmd.color("gray30", xsw)
        for name in ("GTP-2", "3A"):
            o = "c" + name.replace("-", "")
            cmd.load(f"{{SNAP}}/{{eff}}_{{name}}_closest.pdb", o)
            cmd.hide("everything", o)
            cmd.show("cartoon", f"{{o}} and {{swsel}}")
            cmd.color("sc_" + name, f"{{o}} and {{swsel}}")
        cmd.orient(xsw); cmd.zoom(xsw, 3.5)
        cmd.ray(1500, 1350)
        out = f"{{SNAP}}/render_{{eff}}_SwI_closest.png"
        cmd.png(out, dpi=200)
        print("PMLOK", out)

try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
'''
    path = f"{SNAP}/_render_effsnap.py"
    open(path, "w").write(body)
    subprocess.run(["pymol", "-cq", path],
                   stdout=open(f"{SNAP}/_render_effsnap.log", "w"), stderr=subprocess.STDOUT)
    print(open(f"{SNAP}/_render_effsnap.log").read()[-500:])


if __name__ == "__main__":
    main()
