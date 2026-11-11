
from pymol import cmd
from pymol.cgo import CYLINDER, COLOR, SPHERE
import csv
D = '/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan'
REF = 'GTP-2'
SIMS = ['GDP', 'GTP-1', 'GTP-2', '2A', '3A']
EDGE_MIN = 0.45
REGION_RGB = {'N-patch': (0.55, 0.34, 0.29), 'SwI': (0.79, 0.64, 0.15), 'interSw': (0.34, 0.65, 0.6), 'SwII': (0.51, 0.45, 0.7), 'HVR': (0.3, 0.45, 0.69), 'other': (0.6, 0.63, 0.65)}

def region_of(r):
    if 5 <= r <= 7: return "N-patch"
    if 28 <= r <= 42: return "SwI"
    if 43 <= r <= 55: return "interSw"
    if 56 <= r <= 78: return "SwII"
    if 180 <= r <= 190: return "HVR"
    return "other"

def safe(s):
    return "o" + s.replace("-", "").replace(" ", "")     # hyphen-free object name (o2A, oGTP2, ...)

ROI = set(range(5,8)) | set(range(28,79)) | set(range(180,191))   # N-patch+SwI+interSw+SwII+HVR (drop scaffold)
KEY_LABELS = [(37,"T37"),(77,"D76/78"),(6,"N-patch")]
cmd.bg_color("white")
cmd.set("auto_zoom", 0)
cmd.set("ray_opaque_background", 0)
cmd.set("ray_shadows", 0)
cmd.set("cartoon_transparency", 0.82)
cmd.set("label_size", 15); cmd.set("label_color", "black")
for s in SIMS:
    cmd.load(f"{D}/{s}_comm.pdb", safe(s))
for s in SIMS:
    if s != REF:
        cmd.align(safe(s), safe(REF))
cmd.hide("everything"); cmd.show("cartoon"); cmd.color("skyblue")

for s in SIMS:
    o = safe(s)
    # Ca coords for this (aligned) object
    coords = {}
    cmd.iterate_state(1, f"{o} and name CA and polymer", "coords[int(resi)]=(x,y,z)", space={'coords':coords})
    # edges
    net = []
    with open(f"{D}/{s}_edges.csv") as fh:
        for row in csv.DictReader(fh):
            i, j, w = int(row['resid_i']), int(row['resid_j']), float(row['gencorr'])
            if w >= EDGE_MIN and i in ROI and j in ROI and i in coords and j in coords:
                x1,y1,z1 = coords[i]; x2,y2,z2 = coords[j]
                net += [CYLINDER, x1,y1,z1, x2,y2,z2, 0.10, 0.55,0.55,0.58, 0.55,0.55,0.58]
    cmd.load_cgo(net, f"net_{o}")
    # nodes (ROI only)
    nodes = []
    for r,(x,y,z) in coords.items():
        if r not in ROI: continue
        cr,cg,cb = REGION_RGB[region_of(r)]
        nodes += [COLOR, cr,cg,cb, SPHERE, x,y,z, 0.6]
    cmd.load_cgo(nodes, f"nodes_{o}")
    # path (bold)
    path = []
    with open(f"{D}/{s}_path.txt") as fh:
        t = fh.read()
    pres = [int(p) for p in t.split("path=")[1].strip().split("-")] if "path=" in t else []
    pcgo = []
    for k in range(len(pres)-1):
        a,b = pres[k], pres[k+1]
        if a in coords and b in coords:
            x1,y1,z1 = coords[a]; x2,y2,z2 = coords[b]
            pcgo += [CYLINDER, x1,y1,z1, x2,y2,z2, 0.35, 1.0,0.45,0.0, 1.0,0.45,0.0]
    for r in pres:
        if r in coords:
            x,y,z = coords[r]
            cr,cg,cb = REGION_RGB[region_of(r)]
            pcgo += [COLOR, cr,cg,cb, SPHERE, x,y,z, 1.05]
    cmd.load_cgo(pcgo, f"path_{o}")
    for rr,txt in KEY_LABELS:
        cmd.label(f"{o} and name CA and resi {rr}", repr(txt))

cmd.set("ray_trace_mode", 0)
cmd.orient(safe(REF)); cmd.zoom(safe(REF), 4)      # full-protein view, set AFTER all CGO loaded
for s in SIMS:
    o = safe(s)
    cmd.disable("all")
    cmd.enable(o); cmd.enable("net_"+o); cmd.enable("nodes_"+o); cmd.enable("path_"+o)
    cmd.ray(1000, 1000)
    cmd.png(f"{D}/_{s}_netstruct.png", dpi=150)
