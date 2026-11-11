bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set cartoon_putty_scale, 2
set cartoon_putty_radius, 0.4
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/GDP_btw.pdb, GDP
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/GTP-1_btw.pdb, GTP-1
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/GTP-2_btw.pdb, GTP-2
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/2A_btw.pdb, 2A
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/3A_btw.pdb, 3A
align GDP, GTP-2
align GTP-1, GTP-2
align 2A, GTP-2
align 3A, GTP-2
hide everything
show cartoon
cartoon putty, GDP
spectrum b, blue_white_red, GDP, minimum=0, maximum=111.16
show spheres, GDP and name CA and resi 7+56+42
set sphere_scale, 0.6, GDP and name CA and resi 7+56+42
color black, GDP and name CA and resi 7+56+42
cartoon putty, GTP-1
spectrum b, blue_white_red, GTP-1, minimum=0, maximum=111.16
show spheres, GTP-1 and name CA and resi 5+43+42
set sphere_scale, 0.6, GTP-1 and name CA and resi 5+43+42
color black, GTP-1 and name CA and resi 5+43+42
cartoon putty, GTP-2
spectrum b, blue_white_red, GTP-2, minimum=0, maximum=111.16
show spheres, GTP-2 and name CA and resi 5+43+42
set sphere_scale, 0.6, GTP-2 and name CA and resi 5+43+42
color black, GTP-2 and name CA and resi 5+43+42
cartoon putty, 2A
spectrum b, blue_white_red, 2A, minimum=0, maximum=111.16
show spheres, 2A and name CA and resi 7+56+42
set sphere_scale, 0.6, 2A and name CA and resi 7+56+42
color black, 2A and name CA and resi 7+56+42
cartoon putty, 3A
spectrum b, blue_white_red, 3A, minimum=0, maximum=111.16
show spheres, 3A and name CA and resi 7+56+42
set sphere_scale, 0.6, 3A and name CA and resi 7+56+42
color black, 3A and name CA and resi 7+56+42
orient GTP-2
set ray_trace_mode, 1
disable all
enable GDP
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_GDP_btw_render.png, dpi=150
disable all
enable GTP-1
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_GTP-1_btw_render.png, dpi=150
disable all
enable GTP-2
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_GTP-2_btw_render.png, dpi=150
disable all
enable 2A
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_2A_btw_render.png, dpi=150
disable all
enable 3A
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_3A_btw_render.png, dpi=150