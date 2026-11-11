#!/bin/bash
cd /home/hutchcha/rhoa_mut_analysis/rhoa_full_sampled_trajs
source ~/anaconda3/etc/profile.d/conda.sh
conda activate memb_analysis
export NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8
python -u dynetan_network.py --sims GTP-2,2A,3A,GTP-1,GDP --sampled 1200 --contact-stride 20 --ncores 8 --debug > dynetan_full.log 2>&1
echo "DYNETAN_BATCH_DONE exit=$?" >> dynetan_full.log
