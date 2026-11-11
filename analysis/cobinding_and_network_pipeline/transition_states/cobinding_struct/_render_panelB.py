from pymol import cmd

def run():
    cmd.load('/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/cobinding_struct/GTP-2_cobinding_session.pse')
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
    cmd.label("prot and resi 5 and name CA", '"R5"')
    cmd.label("prot and resi 6 and name CA", '"K6"')
    cmd.label("prot and resi 7 and name CA", '"K7"')
    cmd.label("prot and resi 78 and name CA", '"D78"')
    cmd.label("bridge and name P", '"POPS"')

    # zoom in tight on the cluster (interacting sidechains + POPS headgroup, tails excluded)
    cmd.zoom("inter or (bridge and not (name C22+C23+C24+C25+C26+C27+C28+C29+C210+C211+C212+C213+C214+C215+C216+C217+C218+C32+C33+C34+C35+C36+C37+C38+C39+C310+C311+C312+C313+C314+C315+C316))", 3.2)
    cmd.ray(1500, 1300)
    cmd.png('/home/hutchcha/rhoa_mut_analysis/paper_figures/figure6/fig_struct_cobinding.png', dpi=200)
    print("PMLOK wrote /home/hutchcha/rhoa_mut_analysis/paper_figures/figure6/fig_struct_cobinding.png")


try:
    run()
except Exception as e:
    print("PMLERR:", type(e).__name__, str(e))
