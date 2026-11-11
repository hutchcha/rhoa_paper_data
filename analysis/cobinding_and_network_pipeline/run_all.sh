#!/bin/bash
# Full-resolution H-bond + cooperativity pipeline. Runs everything in sequence.
# Each stage is resumable (H-bond searches skip completed sims), so if this is
# interrupted, just relaunch and it continues.
#
# Launch in the background (survives logout) with EITHER:
#   conda activate rhoa_analysis
#   nohup bash run_all.sh > run_all.log 2>&1 &
# OR with tmux (preferred -- lets you re-attach and watch progress):
#   tmux new -d -s rhoa 'conda activate rhoa_analysis; bash run_all.sh 2>&1 | tee run_all.log'
#   tmux attach -t rhoa      # watch;  Ctrl-b then d  to detach
set -u
cd "$(dirname "$0")"
echo "===== START $(date) ====="

echo "=== [0/4] verifying inputs (PSF <-> trajectory match, frame counts) ==="
python -u verify_inputs.py || { echo "INPUT VERIFICATION FAILED -- aborting (fix paths in stable_state_config.py)"; exit 1; }

echo "=== [1/4] H-bond searches (all lipids, stride 1) ==="
python -u run_hbonds.py

echo "=== [2/4] cooperativity map + enhancement (POPS, then all amine lipids) ==="
python -u lipid_hbond_analysis.py --sims 2A,3A,GTP-2,GDP,GTP-1 --lipids POPS
python -u lipid_hbond_analysis.py --sims 2A,3A,GTP-2,GDP,GTP-1 --lipids POPS,POPE,PLA18

echo "=== [3/4] cooperativity networks (POPS + amine, min25 + min50) ==="
python -u viz_cooperativity.py --lipidset ''     --cutoff 25 --min-frames 25
python -u viz_cooperativity.py --lipidset ''     --cutoff 50 --min-frames 50
python -u viz_cooperativity.py --lipidset _amine --cutoff 25 --min-frames 25
python -u viz_cooperativity.py --lipidset _amine --cutoff 50 --min-frames 50

echo "=== [4/4] occupancy / frequency (POPS + all amine) ==="
python -u cobinding_occupancy.py --sims GDP,GTP-1,GTP-2,2A,3A --lipids POPS --tag POPS
python -u cobinding_occupancy.py --sims GDP,GTP-1,GTP-2,2A,3A --lipids POPS,POPE,PLA18

echo "===== ALL DONE $(date) ====="
