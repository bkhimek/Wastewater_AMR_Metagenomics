#!/usr/bin/env bash
# 05_qc_pilot.sh - per-sample QC for verified DEVELOPMENT (Day 1) samples only.
# For each sample with both mates verified: seqkit stats, FastQC, bbduk 38.90 dry run
# (authors' parameters + our assumed adapter options, D11), fastp 0.23.4 dry run (comparison only).
# No trimmed reads are written. Safe to re-run: steps with existing output are skipped.
# Writes results/qc_summary_pilot.tsv. Never touches runs outside ALLOW (evaluation data stay closed).
# Usage: bash bin/05_qc_pilot.sh            (ALLOW="SRR... SRR..." to change the list; THREADS=4)
set -uo pipefail
REPO=${REPO:-$HOME/projects/Wastewater_AMR_Metagenomics}
RAW=$REPO/data/fastq_raw; VER=$RAW/verified.tsv; RES=$REPO/results
ALLOW=${ALLOW:-"SRR22752025 SRR22752015 SRR22752048 SRR22752047 SRR22752038 SRR22752037"}
THREADS=${THREADS:-4}
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/miniconda3}
case "$REPO" in /mnt/c*|/mnt/d*) echo "ERROR: repo must be on the WSL disk"; exit 1;; esac
mkdir -p "$RES/qc_raw_dev" "$RES/qc_trim_dev" "$RES/qc_stats"
for run in $ALLOW; do
  r1=$RAW/${run}_1.fastq.gz; r2=$RAW/${run}_2.fastq.gz
  if [ -f "$r1" ] && [ -f "$r2" ] && grep -q -F "${run}_1.fastq.gz" "$VER" && grep -q -F "${run}_2.fastq.gz" "$VER"; then :; else
    echo "$(date +%T) [$run] skip (both mates not verified yet)"; continue; fi
  echo "$(date +%T) [$run] start"
  st=$RES/qc_stats/seqkit_$run.tsv
  [ -s "$st" ] || { seqkit stats -a -T -j 2 "$r1" "$r2" > "$st.tmp" && mv "$st.tmp" "$st"; }
  [ -s "$RES/qc_raw_dev/${run}_2_fastqc.zip" ] || fastqc -q -t 2 -o "$RES/qc_raw_dev" "$r1" "$r2"
  bl=$RES/qc_trim_dev/bbduk_dryrun_$run.log
  grep -q '^Result:' "$bl" 2>/dev/null || mamba run -n wwamr-bbduk bbduk.sh in="$r1" in2="$r2" out=null \
      ref=adapters ktrim=r k=23 mink=11 hdist=1 tbo tpe maxns=0 minlength=75 qtrim=r trimq=20 \
      threads=$THREADS -Xmx8g 2> "$bl"
  fj=$RES/qc_trim_dev/fastp_$run.json
  [ -s "$fj" ] || mamba run -n wwamr-qc fastp -i "$r1" -I "$r2" --detect_adapter_for_pe --length_required 50 \
      -w $THREADS -j "$fj" -h "$RES/qc_trim_dev/fastp_$run.html" 2> "$RES/qc_trim_dev/fastp_$run.log"
  echo "$(date +%T) [$run] done"
done
python3 -I "$(dirname "$0")/05_collect_qc.py" "$REPO" $ALLOW
