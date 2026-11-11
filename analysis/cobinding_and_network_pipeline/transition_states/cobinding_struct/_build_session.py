from pymol import cmd, util

def run():
    cmd.set("auto_zoom", 0)
    # author lighting/palette (has a space in the filename -> feed line by line)
    for ln in open('/home/hutchcha/rhoa_mut_analysis/paper_figures/figure6/set lights2.pml'):
        s = ln.split("#", 1)[0].strip()
        if s:
            cmd.do(s)

    cmd.load('/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/cobinding_struct/GTP-2_cobind_f9750_membrane.pdb', "frame")
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
    cmd.show("sticks", "bridge and (resn POPS and name N+HN1+HN2+HN3+C11+C12+C13+O13A+O13B+P+O11+O12+O13+O14+C1+C2+C3+HA+HB+HS+HX+HY+O21+C21+O22+O31+C31+O32)")
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
    cmd.orient("inter or (bridge and (resn POPS and name N+HN1+HN2+HN3+C11+C12+C13+O13A+O13B+P+O11+O12+O13+O14+C1+C2+C3+HA+HB+HS+HX+HY+O21+C21+O22+O31+C31+O32))")
    cmd.turn("x", -55); cmd.turn("y", 25)
    cmd.zoom("inter or (bridge and (resn POPS and name N+HN1+HN2+HN3+C11+C12+C13+O13A+O13B+P+O11+O12+O13+O14+C1+C2+C3+HA+HB+HS+HX+HY+O21+C21+O22+O31+C31+O32))", 9.0)
    cmd.set("field_of_view", 25)

    cmd.save('/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/cobinding_struct/GTP-2_cobinding_session.pse')
    print("PMLOK wrote /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/cobinding_struct/GTP-2_cobinding_session.pse")


try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
