from pymol import cmd

GROUPS = {'main': ['GTP-2', '3A'], 'refs': ['GDP', 'GTP-1']}
SIM_COL = {'GDP': (0.7, 0.475, 0.635), 'GTP-1': (0.447, 0.718, 0.698), 'GTP-2': (0.894, 0.341, 0.339), '2A': (0.961, 0.522, 0.094), '3A': (0.298, 0.471, 0.706)}
EFFS = [['ROCK1', '/home/hutchcha/rhoa_mut_analysis/effector_refs/1S1C.pdb', 'A'], ['PKN', '/home/hutchcha/rhoa_mut_analysis/effector_refs/1CXZ.pdb', 'A'], ['anillin', '/home/hutchcha/rhoa_mut_analysis/effector_refs/4XOI.pdb', 'A']]
SNAP = '/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/effector_snapshots'

def run():
    cmd.set("auto_zoom", 0)
    for ln in open('/home/hutchcha/rhoa_mut_analysis/paper_figures/figure6/set lights2.pml'):
        s = ln.split("#", 1)[0].strip()
        if s:
            try: cmd.do(s)
            except Exception: pass
    cmd.set("cartoon_loop_radius", 0.35)
    cmd.set("cartoon_transparency", 0.0)
    for name, col in SIM_COL.items():
        cmd.set_color("sc_" + name, list(col))

    SW = {"SwI": "resi 28-42", "SwII": "resi 56-78"}   # main = SwI only; SwII -> supplement
    for eff, pdb, ch in EFFS:
        for swname, swsel in SW.items():
            for group, gsims in GROUPS.items():            # one render per sim group
                cmd.delete("all")
                cmd.load(pdb, "xtal")
                cmd.hide("everything")
                rho = f"xtal and chain {ch} and polymer"
                cmd.show("cartoon", rho)                   # faint core context
                cmd.color("gray90", rho)
                cmd.set("cartoon_transparency", 0.82, rho)
                xsw = f"xtal and chain {ch} and {swsel}"  # gray crystal reference switch
                cmd.set("cartoon_transparency", 0.0, xsw)
                cmd.color("gray30", xsw)
                for name in gsims:                         # only the sims plotted in this panel
                    o = "o" + name.replace("-", "")
                    cmd.load(f"{SNAP}/{eff}_{name}.pdb", o)
                    cmd.hide("everything", o)
                    cmd.show("cartoon", f"{o} and {swsel}")
                    cmd.color("sc_" + name, f"{o} and {swsel}")
                cmd.orient(xsw)
                cmd.zoom(xsw, 3.5)
                cmd.ray(1500, 1350)
                out = f"{SNAP}/render_{eff}_{swname}_{group}.png"
                cmd.png(out, dpi=200)
                print("PMLOK", out)

    # separate 'how close does each get' overlay: the single closest (min-RMSD) Switch I frame
    # for GTP-2 and 3A on the gray crystal (not used in the KDE composites)
    for eff, pdb, ch in EFFS:
        swsel = "resi 28-42"
        cmd.delete("all")
        cmd.load(pdb, "xtal")
        cmd.hide("everything")
        rho = f"xtal and chain {ch} and polymer"
        cmd.show("cartoon", rho); cmd.color("gray90", rho)
        cmd.set("cartoon_transparency", 0.82, rho)
        xsw = f"xtal and chain {ch} and {swsel}"
        cmd.set("cartoon_transparency", 0.0, xsw); cmd.color("gray30", xsw)
        for name in ("GTP-2", "3A"):
            o = "c" + name.replace("-", "")
            cmd.load(f"{SNAP}/{eff}_{name}_closest.pdb", o)
            cmd.hide("everything", o)
            cmd.show("cartoon", f"{o} and {swsel}")
            cmd.color("sc_" + name, f"{o} and {swsel}")
        cmd.orient(xsw); cmd.zoom(xsw, 3.5)
        cmd.ray(1500, 1350)
        out = f"{SNAP}/render_{eff}_SwI_closest.png"
        cmd.png(out, dpi=200)
        print("PMLOK", out)

try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
