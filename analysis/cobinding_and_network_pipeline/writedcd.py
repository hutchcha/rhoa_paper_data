import scipy.sparse as sparse
import seaborn as sns
import pandas as pd
import numpy as np
import MDAnalysis as mda 
from MDAnalysis.analysis import contacts
import matplotlib.pyplot as plt
import MDAnalysis.transformations as trn
from tqdm import tqdm



u = mda.Universe("nowat.psf", "3A_nowater_10us.dcd")
protein = u.atoms
prot = u.select_atoms("resname CYSG and name CA")
ag = u.select_atoms("not protein")
# we will use mass as weights for the center calculation
workflow = (trn.unwrap(ag, parallelizable=True, max_threads=10),
                   trn.center_in_box(prot, center='mass', parallelizable=True, max_threads=10),
                   trn.wrap(ag, compound='fragments', parallelizable=True, max_threads=10))

u.trajectory.add_transformations(*workflow)
with mda.Writer("3A_cenwrapprod8us.dcd", protein.n_atoms) as W:
    for ts in tqdm(u.trajectory[::10]):
        W.write(protein)
