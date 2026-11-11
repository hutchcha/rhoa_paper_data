bg_color white
set ray_opaque_background, 0
set grid_mode, 1
set cartoon_transparency, 0
set ray_shadows, 0
load /mnt/data_SSD/rhoa_1.8us_trajs/transition_states/dynetan/GTP-2_comm.pdb, GTP-2
load /mnt/data_SSD/rhoa_1.8us_trajs/transition_states/dynetan/2A_comm.pdb, 2A
load /mnt/data_SSD/rhoa_1.8us_trajs/transition_states/dynetan/3A_comm.pdb, 3A
load /mnt/data_SSD/rhoa_1.8us_trajs/transition_states/dynetan/GTP-1_comm.pdb, GTP-1
load /mnt/data_SSD/rhoa_1.8us_trajs/transition_states/dynetan/GDP_comm.pdb, GDP
align 2A, GTP-2
align 3A, GTP-2
align GTP-1, GTP-2
align GDP, GTP-2
hide everything
show cartoon
set cartoon_transparency, 0
color marine, GTP-2 and b>-0.5 and b<0.5
color orange, GTP-2 and b>0.5 and b<1.5
color forest, GTP-2 and b>1.5 and b<2.5
color purple, GTP-2 and b>2.5 and b<3.5
color yellow, GTP-2 and b>3.5 and b<4.5
color salmon, GTP-2 and b>4.5 and b<5.5
color cyan, GTP-2 and b>5.5 and b<6.5
color deepteal, GTP-2 and b>6.5 and b<7.5
color wheat, GTP-2 and b>7.5 and b<8.5
set_name GTP-2, GTP-2
color marine, 2A and b>-0.5 and b<0.5
color orange, 2A and b>0.5 and b<1.5
color forest, 2A and b>1.5 and b<2.5
color purple, 2A and b>2.5 and b<3.5
color yellow, 2A and b>3.5 and b<4.5
color salmon, 2A and b>4.5 and b<5.5
color cyan, 2A and b>5.5 and b<6.5
color deepteal, 2A and b>6.5 and b<7.5
set_name 2A, 2A
color marine, 3A and b>-0.5 and b<0.5
color orange, 3A and b>0.5 and b<1.5
color forest, 3A and b>1.5 and b<2.5
color purple, 3A and b>2.5 and b<3.5
color yellow, 3A and b>3.5 and b<4.5
color salmon, 3A and b>4.5 and b<5.5
color cyan, 3A and b>5.5 and b<6.5
color deepteal, 3A and b>6.5 and b<7.5
color wheat, 3A and b>7.5 and b<8.5
set_name 3A, 3A
color marine, GTP-1 and b>-0.5 and b<0.5
color orange, GTP-1 and b>0.5 and b<1.5
color forest, GTP-1 and b>1.5 and b<2.5
color purple, GTP-1 and b>2.5 and b<3.5
color yellow, GTP-1 and b>3.5 and b<4.5
color salmon, GTP-1 and b>4.5 and b<5.5
color cyan, GTP-1 and b>5.5 and b<6.5
color deepteal, GTP-1 and b>6.5 and b<7.5
color wheat, GTP-1 and b>7.5 and b<8.5
set_name GTP-1, GTP-1
color marine, GDP and b>-0.5 and b<0.5
color orange, GDP and b>0.5 and b<1.5
color forest, GDP and b>1.5 and b<2.5
color purple, GDP and b>2.5 and b<3.5
color yellow, GDP and b>3.5 and b<4.5
color salmon, GDP and b>4.5 and b<5.5
color cyan, GDP and b>5.5 and b<6.5
color deepteal, GDP and b>6.5 and b<7.5
color wheat, GDP and b>7.5 and b<8.5
color slate, GDP and b>8.5 and b<9.5
set_name GDP, GDP
orient GTP-2
set ray_trace_mode, 1
ray 2400, 1400
png /mnt/data_SSD/rhoa_1.8us_trajs/transition_states/dynetan/dynetan_communities.png, dpi=150