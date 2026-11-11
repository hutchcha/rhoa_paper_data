from pymol import cmd
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0)
cmd.set("label_size", 22); cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
cmd.set("dash_color", "gray40")


def run():
    cmd.load('/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/cobinding_struct/GTP-2_cobind_f9750_lip323.pdb', "sys")
    cmd.hide("everything")
    cmd.show("cartoon", "polymer"); cmd.color("gray80", "polymer")
    cmd.set("cartoon_transparency", 0.6)
    cmd.set("cartoon_side_chain_helper", 1)

    # POPS headgroup only (serine + phosphate + glycerol/ester carbonyls); acyl tails hidden
    head = ("resn POPS and resi 323 and name "
            "N+HN1+HN2+HN3+C11+C12+C13+O13A+O13B+P+O11+O12+O13+O14+"
            "C1+C2+C3+HA+HB+HS+HX+HY+O21+C21+O22+O31+C31+O32")
    cmd.select("lip", head)
    cmd.show("sticks", "lip"); cmd.color("wheat", "lip and elem C")
    cmd.select("bas", "resi 5+6+7 and sidechain")
    cmd.select("aci", "resi 78 and sidechain")
    cmd.show("sticks", "bas or aci")
    cmd.color("marine", "bas and elem C"); cmd.color("firebrick", "aci and elem C")
    cmd.util.cnc("lip or bas or aci")
    cmd.set("stick_radius", 0.22)
    cmd.hide("everything", "hydro")          # heavy atoms only -> declutter

    cmd.label("resi 5 and name CA and polymer", '"R5"')
    cmd.label("resi 6 and name CA and polymer", '"K6"')
    cmd.label("resi 7 and name CA and polymer", '"K7"')
    cmd.label("resi 78 and name CA and polymer", '"D78"')

    cmd.distance("d_bas", "(bas and name NZ+NH1+NH2+NE)", "lip", 4.0, mode=0)
    cmd.distance("d_aci", "(aci and name OD1+OD2+OE1+OE2)", "lip", 4.0, mode=0)
    cmd.hide("labels", "d_bas"); cmd.hide("labels", "d_aci")
    cmd.color("marine", "d_bas"); cmd.color("firebrick", "d_aci")
    cmd.set("dash_radius", 0.06); cmd.set("dash_gap", 0.35)

    cmd.orient("lip or bas or aci"); cmd.zoom("lip or bas or aci", 3.0); cmd.turn("y", 15)
    cmd.ray(1400, 1200)
    cmd.png("/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/cobinding_struct/GTP-2_cobinding.png", dpi=200)
    print("PMLOK wrote /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/cobinding_struct/GTP-2_cobinding.png")


try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
