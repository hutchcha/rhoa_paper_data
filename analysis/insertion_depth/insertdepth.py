import matplotlib.pyplot as plt 
from lipyphilic.lib.assign_leaflets import AssignLeaflets
from lipyphilic.lib.area_per_lipid import AreaPerLipid
import scipy.sparse as sparse
import seaborn as sns
import pandas as pd
import numpy as np
import MDAnalysis as mda 
import matplotlib.ticker as tkr
import matplotlib.pylab as plb
from tqdm import tqdm


rheb = mda.Universe("rhebanchor/sys.psf", "rhebanchor/cenwrapprod.dcd")
rhoa = mda.Universe("rhoanchor/2us/sys.psf", "rhoanchor/2us/cenwrapprod.dcd")
diras = mda.Universe("diras3nte/nowat.psf", "diras3nte/../rotate10/cenwrapprodcom.dcd")
names = {0:"rheb", 1:"rhoa", 2:"diras3nte"}
unis = [rheb, rhoa, diras]
insert_dic = {}
for i, u in enumerate(unis):
    
    insert = []
    print(names[i])
    for ts in tqdm(u.trajectory[::1]):
        memcom = u.select_atoms("resname POPC POPS POPE SAPI CHL1 PLA18 PSM and not name H*").center_of_mass()[2]
        if u is rheb:
            anchor = u.select_atoms("resname CYSF and not name CA HA C CB SG CT OT1 OT2 O H*").center_of_mass()[2]
        elif u is rhoa:
            anchor = u.select_atoms("resname CYSG and name CA HA C CB SG CT OT1 OT2 O H*").center_of_mass()[2]
        elif u is diras:
            anchor = u.select_atoms("resname GLYM and not name CA HA C CB SG CT OT1 OT2 O H*").center_of_mass()[2]

        insert_depth = memcom - anchor
        insert.append(insert_depth)
    insert_dic[f'{names[i]}'] = np.array(insert)
import pickle
with open('insertdepthdictionary.pkl', 'wb') as file:
    pickle.dump(insert_dic, file)