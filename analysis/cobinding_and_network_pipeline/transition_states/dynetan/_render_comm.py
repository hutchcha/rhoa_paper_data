
from pymol import cmd
from pymol.cgo import CYLINDER, COLOR, SPHERE
import csv
from collections import defaultdict
D = '/home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs/transition_states/dynetan'; REF = 'GTP-2'; SIMS = ['GDP', 'GTP-1', 'GTP-2', '2A', '3A']
COMMCOL = {'GDP': {0: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765], 2: [1.0, 0.4980392156862745, 0.054901960784313725], 4: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 6: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765], 8: [1.0, 0.4980392156862745, 0.054901960784313725], 10: [0.8901960784313725, 0.4666666666666667, 0.7607843137254902], 1: [0.5803921568627451, 0.403921568627451, 0.7411764705882353], 3: [0.5490196078431373, 0.33725490196078434, 0.29411764705882354], 7: [0.7372549019607844, 0.7411764705882353, 0.13333333333333333], 5: [0.4980392156862745, 0.4980392156862745, 0.4980392156862745], 9: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765]}, 'GTP-1': {0: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765], 1: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 4: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 6: [0.7372549019607844, 0.7411764705882353, 0.13333333333333333], 7: [1.0, 0.4980392156862745, 0.054901960784313725], 9: [0.5490196078431373, 0.33725490196078434, 0.29411764705882354], 5: [0.5803921568627451, 0.403921568627451, 0.7411764705882353], 2: [0.8901960784313725, 0.4666666666666667, 0.7607843137254902], 3: [0.4980392156862745, 0.4980392156862745, 0.4980392156862745], 8: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765]}, 'GTP-2': {7: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765], 1: [1.0, 0.4980392156862745, 0.054901960784313725], 3: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 8: [0.8392156862745098, 0.15294117647058825, 0.1568627450980392], 5: [0.5803921568627451, 0.403921568627451, 0.7411764705882353], 0: [0.5490196078431373, 0.33725490196078434, 0.29411764705882354], 4: [0.8901960784313725, 0.4666666666666667, 0.7607843137254902], 2: [0.4980392156862745, 0.4980392156862745, 0.4980392156862745], 6: [0.7372549019607844, 0.7411764705882353, 0.13333333333333333]}, '2A': {0: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765], 4: [1.0, 0.4980392156862745, 0.054901960784313725], 2: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 3: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 5: [0.4980392156862745, 0.4980392156862745, 0.4980392156862745], 6: [0.5490196078431373, 0.33725490196078434, 0.29411764705882354], 7: [0.8901960784313725, 0.4666666666666667, 0.7607843137254902], 1: [0.7372549019607844, 0.7411764705882353, 0.13333333333333333]}, '3A': {6: [0.12156862745098039, 0.4666666666666667, 0.7058823529411765], 3: [1.0, 0.4980392156862745, 0.054901960784313725], 0: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 4: [0.17254901960784313, 0.6274509803921569, 0.17254901960784313], 7: [1.0, 0.4980392156862745, 0.054901960784313725], 2: [0.4980392156862745, 0.4980392156862745, 0.4980392156862745], 1: [0.5490196078431373, 0.33725490196078434, 0.29411764705882354], 5: [0.8901960784313725, 0.4666666666666667, 0.7607843137254902], 8: [0.7372549019607844, 0.7411764705882353, 0.13333333333333333]}}
GRAY = [0.6, 0.6, 0.6]

def safe(s): return "o" + s.replace("-", "").replace(" ", "")

cmd.bg_color("white"); cmd.set("auto_zoom", 0)
cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0)
cmd.set("cartoon_transparency", 0.2)
cmd.set("cartoon_loop_radius", 0.35)
cmd.set("label_size", 18); cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
for s in SIMS: cmd.load(f"{D}/{s}_comm.pdb", safe(s))
for s in SIMS:
    if s != REF: cmd.align(safe(s), safe(REF))
cmd.hide("everything"); cmd.show("cartoon")

for s in SIMS:
    o = safe(s)
    coords = {}
    cmd.iterate_state(1, f"{o} and name CA and polymer", "coords[int(resi)]=(x,y,z)", space={'coords':coords})
    comm = {}
    with open(f"{D}/{s}_nodes.csv") as fh:
        for row in csv.DictReader(fh):
            comm[int(row['resid'])] = int(row['community'])
    members = defaultdict(list)
    for r,c in comm.items():
        if r in coords: members[c].append(r)
    cent = {}; size = {}
    for c,rs in members.items():
        xs = [coords[r] for r in rs]
        cent[c] = (sum(p[0] for p in xs)/len(xs), sum(p[1] for p in xs)/len(xs), sum(p[2] for p in xs)/len(xs))
        size[c] = len(rs)
    # color cartoon by WT-matched community color
    for c,rs in members.items():
        col = COMMCOL[s].get(c, GRAY)
        cmd.set_color(f"cc_{o}_{c}", list(col))
        cmd.color(f"cc_{o}_{c}", f"{o} and resi " + "+".join(str(r) for r in rs))
    # inter-community coupling (sum gencorr of crossing edges)
    coup = defaultdict(float)
    with open(f"{D}/{s}_edges.csv") as fh:
        for row in csv.DictReader(fh):
            i,j,w = int(row['resid_i']), int(row['resid_j']), float(row['gencorr'])
            if i in comm and j in comm and comm[i] != comm[j]:
                coup[tuple(sorted((comm[i],comm[j])))] += w
    maxc = max(coup.values()) if coup else 1.0
    ce = []
    for (a,b),w in coup.items():
        if w < 0.15*maxc: continue               # drop weak inter-community couplings (declutter)
        if a in cent and b in cent:
            x1,y1,z1 = cent[a]; x2,y2,z2 = cent[b]
            rad = 0.15 + 0.9*(w/maxc)             # thinner edges, scaled by coupling strength
            ce += [CYLINDER, x1,y1,z1, x2,y2,z2, rad, 0.35,0.35,0.38, 0.35,0.35,0.38]
    cmd.load_cgo(ce, f"cedges_{o}")
    cn = []
    maxsz = max(size.values()) if size else 1
    for c,(x,y,z) in cent.items():
        col = COMMCOL[s].get(c, GRAY)
        rad = 0.9 + 2.0*(size[c]/maxsz)          # node radius linear in community size (0.9..2.9)
        cn += [COLOR, col[0],col[1],col[2], SPHERE, x,y,z, rad]
    cmd.load_cgo(cn, f"cnodes_{o}")
    # (region text labels removed -- the sequence track now serves as the colour/region key)

cmd.orient(safe(REF)); cmd.zoom(safe(REF), 4)
for s in SIMS:
    o = safe(s)
    cmd.disable("all")
    for pre in ("", "cedges_", "cnodes_"):
        cmd.enable(pre + o)
    cmd.ray(1000, 1000)
    cmd.png(f"{D}/_{s}_commnet.png", dpi=150)
