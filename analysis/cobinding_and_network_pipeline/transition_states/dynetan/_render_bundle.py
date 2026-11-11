
from pymol import cmd
from pymol.cgo import CYLINDER, COLOR, SPHERE
import csv
D = '/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan'; REF = 'GTP-2'; SIMS = ['GDP', 'GTP-1', 'GTP-2', '2A', '3A']; K = 25
REGION_RGB = {'N-patch': (0.55, 0.34, 0.29), 'SwI': (0.79, 0.64, 0.15), 'interSw': (0.34, 0.65, 0.6), 'SwII': (0.51, 0.45, 0.7), 'HVR': (0.3, 0.45, 0.69), 'other': (0.6, 0.63, 0.65)}
def region_of(r):
    if 5<=r<=7: return "N-patch"
    if 28<=r<=42: return "SwI"
    if 43<=r<=55: return "interSw"
    if 56<=r<=78: return "SwII"
    if 180<=r<=190: return "HVR"
    return "other"
def safe(s): return "o"+s.replace("-","").replace(" ","")
cmd.bg_color("white"); cmd.set("auto_zoom",0); cmd.set("ray_opaque_background",0); cmd.set("ray_shadows",0)
cmd.set("cartoon_transparency",0.8); cmd.set("label_size",16); cmd.set("label_color","black")
for s in SIMS: cmd.load(f"{D}/{s}_comm.pdb", safe(s))
for s in SIMS:
    if s!=REF: cmd.align(safe(s), safe(REF))
cmd.hide("everything"); cmd.show("cartoon"); cmd.color("gray85")
for s in SIMS:
    o=safe(s)
    coords={}
    cmd.iterate_state(1, f"{o} and name CA and polymer", "coords[int(resi)]=(x,y,z)", space={'coords':coords})
    ec=[]
    with open(f"{D}/{s}_bundle_edges.csv") as fh:
        for row in csv.DictReader(fh):
            i,j,c=int(row['i']),int(row['j']),int(row['count'])
            if i in coords and j in coords:
                frac=c/K; rad=0.15+1.4*frac
                x1,y1,z1=coords[i]; x2,y2,z2=coords[j]
                ec+=[CYLINDER, x1,y1,z1, x2,y2,z2, rad, 1.0,0.5,0.0, 1.0,0.5,0.0]
    cmd.load_cgo(ec, f"be_{o}")
    nc=[]
    with open(f"{D}/{s}_bundle_nodes.csv") as fh:
        for row in csv.DictReader(fh):
            r,c=int(row['resid']),int(row['count'])
            if r in coords:
                frac=c/K; rad=0.4+1.1*frac
                cr,cg,cb=REGION_RGB[region_of(r)]
                x,y,z=coords[r]; nc+=[COLOR,cr,cg,cb, SPHERE,x,y,z,rad]
    cmd.load_cgo(nc, f"bn_{o}")
    for rr,txt in [(6,"N-patch"),(37,"T37"),(76,"D76/78")]:
        if rr in coords:
            cmd.label(f"{o} and name CA and resi {rr}", repr(txt))
cmd.orient(safe(REF)); cmd.zoom(safe(REF), 4)
for s in SIMS:
    o=safe(s); cmd.disable("all")
    cmd.enable(o); cmd.enable("be_"+o); cmd.enable("bn_"+o)
    cmd.ray(1000,1000); cmd.png(f"{D}/_{s}_bundle.png", dpi=150)
