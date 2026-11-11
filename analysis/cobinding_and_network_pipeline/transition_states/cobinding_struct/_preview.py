from pymol import cmd
cmd.load("GTP-2_cobinding_session.pse")
cmd.set("ray_trace_mode", 1)
cmd.ray(1400, 1150)
cmd.png("_session_preview.png", dpi=150)
print("PREVIEW_DONE")
