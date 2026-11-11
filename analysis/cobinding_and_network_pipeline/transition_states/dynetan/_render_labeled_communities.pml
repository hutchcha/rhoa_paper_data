bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set cartoon_putty_scale, 2
set cartoon_putty_radius, 0.4
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/GDP_comm.pdb, GDP
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/GTP-1_comm.pdb, GTP-1
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/GTP-2_comm.pdb, GTP-2
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/2A_comm.pdb, 2A
load /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/3A_comm.pdb, 3A
align GDP, GTP-2
align GTP-1, GTP-2
align 2A, GTP-2
align 3A, GTP-2
hide everything
show cartoon
color marine, GDP and b>-0.5 and b<0.5
color orange, GDP and b>0.5 and b<1.5
color forest, GDP and b>1.5 and b<2.5
color purpleblue, GDP and b>2.5 and b<3.5
color yellow, GDP and b>3.5 and b<4.5
color salmon, GDP and b>4.5 and b<5.5
color cyan, GDP and b>5.5 and b<6.5
color deepteal, GDP and b>6.5 and b<7.5
color wheat, GDP and b>7.5 and b<8.5
color slate, GDP and b>8.5 and b<9.5
color marine, GTP-1 and b>-0.5 and b<0.5
color orange, GTP-1 and b>0.5 and b<1.5
color forest, GTP-1 and b>1.5 and b<2.5
color purpleblue, GTP-1 and b>2.5 and b<3.5
color yellow, GTP-1 and b>3.5 and b<4.5
color salmon, GTP-1 and b>4.5 and b<5.5
color cyan, GTP-1 and b>5.5 and b<6.5
color deepteal, GTP-1 and b>6.5 and b<7.5
color wheat, GTP-1 and b>7.5 and b<8.5
color marine, GTP-2 and b>-0.5 and b<0.5
color orange, GTP-2 and b>0.5 and b<1.5
color forest, GTP-2 and b>1.5 and b<2.5
color purpleblue, GTP-2 and b>2.5 and b<3.5
color yellow, GTP-2 and b>3.5 and b<4.5
color salmon, GTP-2 and b>4.5 and b<5.5
color cyan, GTP-2 and b>5.5 and b<6.5
color deepteal, GTP-2 and b>6.5 and b<7.5
color wheat, GTP-2 and b>7.5 and b<8.5
color marine, 2A and b>-0.5 and b<0.5
color orange, 2A and b>0.5 and b<1.5
color forest, 2A and b>1.5 and b<2.5
color purpleblue, 2A and b>2.5 and b<3.5
color yellow, 2A and b>3.5 and b<4.5
color salmon, 2A and b>4.5 and b<5.5
color cyan, 2A and b>5.5 and b<6.5
color deepteal, 2A and b>6.5 and b<7.5
color marine, 3A and b>-0.5 and b<0.5
color orange, 3A and b>0.5 and b<1.5
color forest, 3A and b>1.5 and b<2.5
color purpleblue, 3A and b>2.5 and b<3.5
color yellow, 3A and b>3.5 and b<4.5
color salmon, 3A and b>4.5 and b<5.5
color cyan, 3A and b>5.5 and b<6.5
color deepteal, 3A and b>6.5 and b<7.5
color wheat, 3A and b>7.5 and b<8.5
orient GTP-2
set ray_trace_mode, 1
disable all
enable GDP
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_GDP_comm_render.png, dpi=150
disable all
enable GTP-1
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_GTP-1_comm_render.png, dpi=150
disable all
enable GTP-2
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_GTP-2_comm_render.png, dpi=150
disable all
enable 2A
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_2A_comm_render.png, dpi=150
disable all
enable 3A
ray 1000, 1000
png /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan/_3A_comm_render.png, dpi=150