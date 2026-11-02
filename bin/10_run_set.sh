#!/usr/bin/env bash
# 10_run_set.sh SPLIT MATRIX    e.g.  ALLOW_EVAL=1 bash bin/10_run_set.sh eval liquid
# Runs 06_trim_kma.sh (the registered Arm A pipeline) once for every run of one set of the sample sheet, one run at a time.
# Finished runs are skipped. A run that fails is reported at the end and NOT retried automatically; its partial outputs are
# kept for inspection (nothing is overwritten). Log: logs/run_<split>_<matrix>.log. One job at a time (flock).
set -uo pipefail
SPLIT=${1:?split (dev|eval)}; MATRIX=${2:?matrix (liquid|sediment)}
REPO=${REPO:-$HOME/projects/Wastewater_AMR_Metagenomics}; cd "$REPO" || exit 1
SHEET=phaseA_audit/out/sample_sheet_draft.tsv; mkdir -p logs
exec 8>logs/.run_set.lock; flock -n 8 || { echo "another run job is active"; exit 1; }
LOG=logs/run_${SPLIT}_${MATRIX}.log; exec > >(tee -a "$LOG") 2>&1
RUNS=$(awk -F'\t' -v s="$SPLIT" -v m="$MATRIX" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["split"]==s && $h["matrix"]==m {print $h["run_accession"]}' "$SHEET")
N=$(echo "$RUNS" | wc -l); echo "$(date +%F\ %T) set $SPLIT $MATRIX: $N runs"
failed=""; k=0
for RUN in $RUNS; do
  k=$((k+1)); echo "$(date +%T) ===== run $k of $N: $RUN"
  bash bin/06_trim_kma.sh "$RUN" || { failed="$failed $RUN"; echo "$(date +%T) FAILED $RUN"; }
done
done_n=0; for RUN in $RUNS; do [ -s results/kma/${RUN}_resfinder_2.1.1.res ] && [ -s results/kma/${RUN}_silva_138.1_NR99.res ] && done_n=$((done_n+1)); done
echo "$(date +%F\ %T) DONE: $done_n of $N runs complete for $SPLIT $MATRIX; failed:${failed:- none}"
[ "$done_n" = "$N" ]
