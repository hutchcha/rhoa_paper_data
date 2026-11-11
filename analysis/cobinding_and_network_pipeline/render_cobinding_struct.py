#!/usr/bin/env python3
"""Render the WT (GTP-2) lipid co-binding snapshot chosen by scan_cobinding_frame.py:
protein cartoon faint, the shared inner-leaflet POPS as sticks, the N-patch basics (R5/K6/K7)
and Switch II acidic (D78) as sticks, with polar-contact dashes to the lipid. Structural
companion to the co-binding network panel (Fig 6C): basic residues grip the POPS phosphate,
the acidic D78 grips the POPS serine amine -- a lipid-bridged basic+acidic cluster.

Anaconda's PyMOL (python 3.13) is the only working build here; its traceback formatter is buggy,
so the emitted script wraps everything in try/except and prints the real error itself.

  conda activate memb_analysis      # for scan step only
  python render_cobinding_struct.py --write-pml
  pymol -cq transition_states/cobinding_struct/_render_cobind.py   # anaconda pymol
"""
import os, argparse, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
D = f"{HERE}/transition_states/cobinding_struct"
AA1 = {5: "R5", 6: "K6", 7: "K7", 76: "D76", 78: "D78"}


def read_chosen():
    info = {}
    for ln in open(f"{D}/chosen.txt"):
        k, _, v = ln.strip().partition("=")
        info[k] = v
    return info


def parse_ids(s):                       # "(5, 6, 7)" -> [5, 6, 7]
    return [int(x) for x in s.strip("() ").replace(" ", "").split(",") if x]


def write_pml():
    info = read_chosen()
    pdb = info["pdb"]; lr = info["lipid_resid"]; lt = info["lipid_type"]
    bas_ids = parse_ids(info["npatch"]); aci_ids = parse_ids(info["acid"])
    basics = "+".join(map(str, bas_ids)); acids = "+".join(map(str, aci_ids))
    labels = "\n    ".join(                 # 2nd arg is a python expr -> needs embedded quotes
        f'cmd.label("resi {r} and name CA and polymer", {(chr(34) + AA1.get(r, str(r)) + chr(34))!r})'
        for r in bas_ids + aci_ids)
    cmds = f'''
    cmd.load({pdb!r}, "sys")
    cmd.hide("everything")
    cmd.show("cartoon", "polymer"); cmd.color("gray80", "polymer")
    cmd.set("cartoon_transparency", 0.6)
    cmd.set("cartoon_side_chain_helper", 1)

    # POPS headgroup only (serine + phosphate + glycerol/ester carbonyls); acyl tails hidden
    head = ("resn {lt} and resi {lr} and name "
            "N+HN1+HN2+HN3+C11+C12+C13+O13A+O13B+P+O11+O12+O13+O14+"
            "C1+C2+C3+HA+HB+HS+HX+HY+O21+C21+O22+O31+C31+O32")
    cmd.select("lip", head)
    cmd.show("sticks", "lip"); cmd.color("wheat", "lip and elem C")
    cmd.select("bas", "resi {basics} and sidechain")
    cmd.select("aci", "resi {acids} and sidechain")
    cmd.show("sticks", "bas or aci")
    cmd.color("marine", "bas and elem C"); cmd.color("firebrick", "aci and elem C")
    cmd.util.cnc("lip or bas or aci")
    cmd.set("stick_radius", 0.22)
    cmd.hide("everything", "hydro")          # heavy atoms only -> declutter

    {labels}

    cmd.distance("d_bas", "(bas and name NZ+NH1+NH2+NE)", "lip", 4.0, mode=0)
    cmd.distance("d_aci", "(aci and name OD1+OD2+OE1+OE2)", "lip", 4.0, mode=0)
    cmd.hide("labels", "d_bas"); cmd.hide("labels", "d_aci")
    cmd.color("marine", "d_bas"); cmd.color("firebrick", "d_aci")
    cmd.set("dash_radius", 0.06); cmd.set("dash_gap", 0.35)

    cmd.orient("lip or bas or aci"); cmd.zoom("lip or bas or aci", 3.0); cmd.turn("y", 15)
    cmd.ray(1400, 1200)
    cmd.png("{D}/GTP-2_cobinding.png", dpi=200)
    print("PMLOK wrote {D}/GTP-2_cobinding.png")'''
    body = f'''from pymol import cmd
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0)
cmd.set("label_size", 22); cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
cmd.set("dash_color", "gray40")


def run():{cmds}


try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
'''
    path = f"{D}/_render_cobind.py"
    open(path, "w").write(body)
    print(f"wrote {path}\nchosen: frame {info['frame']} lipid {lr} ({lt}) "
          f"N-patch {info['npatch']} acid {info['acid']}")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-pml", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.parse_args()
    path = write_pml()
    # run via anaconda pymol (python 3.13); do NOT capture_output (breaks ray)
    subprocess.run(["pymol", "-cq", path],
                   stdout=open(f"{D}/_render_cobind.log", "w"), stderr=subprocess.STDOUT)
    print(f"see {D}/_render_cobind.log")


if __name__ == "__main__":
    main()
