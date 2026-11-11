import warnings; warnings.filterwarnings("ignore")
import os, time, numba
print("NUMBA_NUM_THREADS =", numba.config.NUMBA_NUM_THREADS, flush=True)
from dynetan.proctraj import DNAproc
d = DNAproc(notebookMode=False); d.setNumWinds(1); d.setNumSampledFrames(1000)
d.loadSystem("topologies/3A.psf", "3A_last1.95us.dcd"); d.setSegIDs(["PROA"])
d.selectSystem(withSolvent=True, inputSelStr="protein and resid 4-189")
d.setCutoffDist(4.5); d.setContactPersistence(0.75)
d.prepareNetwork(autocomp_groups=True); d.alignTraj(inMemory=True)
d.findContacts(stride=30, n_cores=8); d.filterContacts()
t = time.time(); d.calcCor(ncores=8, corr_type="gencor")
print("CALCCOR_DONE %.0fs" % (time.time() - t), flush=True)
