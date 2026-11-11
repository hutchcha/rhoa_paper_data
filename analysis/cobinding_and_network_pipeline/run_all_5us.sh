#!/bin/bash
# Full 5 us re-analysis pipeline (GTP-2, 3A, GDP, GTP-1 -- 2A dropped) + figure regeneration.
# Resumable: each stage writes a marker in _pipeline_5us/ and is skipped if already done.
# Detached-safe: run with  setsid nohup bash run_all_5us.sh > _pipeline_5us/driver.log 2>&1 &
#
# Envs: memb_analysis for python analyses; base (anaconda pymol) for structure renders.
set -u
cd "$(dirname "$0")"
MARK=_pipeline_5us
mkdir -p "$MARK"
SIMS="GDP,GTP-1,GTP-2,3A"
source "$HOME/anaconda3/etc/profile.d/conda.sh"

log(){ echo "[$(date '+%F %T')] $*"; }

# run a stage once: stage <name> <conda-env> <command...>
stage(){
  local name="$1" env="$2"; shift 2
  if [ -f "$MARK/$name.done" ]; then log "SKIP $name (done)"; return 0; fi
  log "START $name"
  conda activate "$env"
  if "$@" > "$MARK/$name.log" 2>&1; then
    touch "$MARK/$name.done"; log "DONE  $name"
  else
    log "FAIL  $name (see $MARK/$name.log)"; return 1
  fi
}

# ---- co-binding ----
stage hbonds        memb_analysis python run_hbonds.py --sims "$SIMS" --stride 1 || exit 1
stage cobind_occ    memb_analysis python cobinding_occupancy.py --sims "$SIMS" --lipids POPS,POPE,PLA18 || exit 1
stage cobind_enh    memb_analysis python lipid_hbond_analysis.py --sims "$SIMS" --lipids POPS,POPE,PLA18 || exit 1
# regenerate Fig 6 co-binding panels + composite
stage fig6          memb_analysis bash -c "python fig_B_amine_frequency.py && python fig_cobinding_network.py && python compose_figure6.py"

# ---- effector competence ----
stage effector_rmsd memb_analysis python fig_D_effector_rmsd.py --stride 2 || exit 1
stage effector_snap_extract memb_analysis python extract_effector_snapshots.py || exit 1
stage effector_snap_render  base          python render_effector_snapshots.py

# ---- dynamical network (per-sim so a crash resumes at the right sim) ----
export NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
conda activate memb_analysis
for s in GDP GTP-1 GTP-2 3A; do
  if [ -f "$MARK/dynetan_$s.done" ]; then log "SKIP dynetan_$s"; continue; fi
  log "START dynetan_$s"
  if python dynetan_network.py --sims "$s" --sampled 2000 --contact-stride 10 --ncores 8 \
        --corr-frames 10000 \
        > "$MARK/dynetan_$s.log" 2>&1; then
    touch "$MARK/dynetan_$s.done"; log "DONE  dynetan_$s"
  else
    log "FAIL  dynetan_$s"; exit 1
  fi
done

# regenerate Fig 7 network pieces + composites (needs effector snapshots + dynetan)
stage fig7_track    memb_analysis python dynetan_community_track.py
stage fig7_render   base          python render_community_network.py
stage fig7_quant    memb_analysis python fig_dynetan_pathlength.py
stage fig7_compose  memb_analysis python compose_figure7.py

# ---- 2A global-similarity control supplement (5 us) ----
stage supp_2A       memb_analysis python make_2A_similarity_supp.py

log "ALL DONE"
touch "$MARK/ALL_DONE"
