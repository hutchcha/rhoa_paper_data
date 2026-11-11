#!/bin/bash
# Run the dynetan MI generalized-correlation network analysis on a many-core machine.
#
# dynetan's KSG mutual-information is CPU/numba-parallel (scales with cores, NOT the GPU),
# so it just wants lots of CPU threads. This drives dynetan_network.py over all 5 sims and
# then renders the structure images with render_dynetan.py.
#
# Setup (once):   conda env create -f environment_dynetan.yml
# Run:            bash run_dynetan_gpu.sh
#   env knobs:    TRAJ_DIR=/path/to/excerpts   (if the *_last1.95us.dcd files aren't in this dir)
#                 NCORES=32  SAMPLED=2000  CSTRIDE=10
#
# Needs, alongside this script: dynetan_network.py, render_dynetan.py, topologies/<sim>.psf,
# and the five <sim>_last1.95us.dcd excerpts (here, or under $TRAJ_DIR).
set -e
cd "$(dirname "$0")"

# activate env (edit if your conda lives elsewhere)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rhoa_dynetan

NC=${NCORES:-$(nproc)}
export NUMBA_NUM_THREADS=$NC OMP_NUM_THREADS=$NC MKL_NUM_THREADS=$NC
SAMPLED=${SAMPLED:-2000}
CSTRIDE=${CSTRIDE:-10}

echo "dynetan: $NC cores | sampled=$SAMPLED | contact-stride=$CSTRIDE | TRAJ_DIR=${TRAJ_DIR:-<bundle>}"

# pre-flight: list which excerpts are found
python - <<'PY'
import os
from dynetan_network import resolve_traj, ORDER
for n in ORDER:
    p = resolve_traj(n)
    print(f"  {n:6} {'OK  '+p if p else 'MISSING'}")
PY

python -u dynetan_network.py --sims GTP-2,2A,3A,GTP-1,GDP \
    --sampled "$SAMPLED" --contact-stride "$CSTRIDE" --ncores "$NC" --debug \
    > dynetan_full.log 2>&1
echo "DYNETAN_BATCH_DONE exit=$?" | tee -a dynetan_full.log

# structure rendering (needs pymol; skips cleanly if absent)
python -u render_dynetan.py --sims GTP-2,2A,3A,GTP-1,GDP >> dynetan_full.log 2>&1 \
    && echo "rendered dynetan_communities.png / dynetan_betweenness.png" \
    || echo "render step skipped (pymol not found?) — copy transition_states/dynetan/ back and render locally"

echo "DONE. Copy transition_states/dynetan/ back to the analysis machine."
