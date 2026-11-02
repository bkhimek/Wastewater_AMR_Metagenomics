#!/usr/bin/env bash
# 13_explore_raw.sh RUN [RUN ...]
# EXPLORATORY (post-freeze, not part of the R1 verdict; review step C, variant "raw"):
# KMA 1.4.9 against ResFinder 2.1.1 with the SAME flags as 06_trim_kma.sh, but on the RAW reads (no bbduk at all).
# Only one thing differs from the registered run: the trimming step is skipped. Output goes to results/explore/ with the
# suffix _raw; the registered outputs in results/kma/ are never touched or overwritten.
# Usage: bash bin/13_explore_raw.sh SRR22752046 SRR22752045      (E021 = discrepant, E026 = matched control)
set -uo pipefail
REPO=${REPO:-$HOME/projects/Wastewater_AMR_Metagenomics}; cd "$REPO" || exit 1
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/miniconda3}
IDX=refdb/kma_index; OUTD=results/explore; mkdir -p "$OUTD" logs
for RUN in "$@"; do
  R1=data/fastq_raw/${RUN}_1.fastq.gz; R2=data/fastq_raw/${RUN}_2.fastq.gz
  grep -q -F "${RUN}_1.fastq.gz" data/fastq_raw/verified.tsv && grep -q -F "${RUN}_2.fastq.gz" data/fastq_raw/verified.tsv || { echo "ERROR: $RUN not verified"; exit 1; }
  OUT=$OUTD/${RUN}_raw_resfinder_2.1.1
  if [ -s "$OUT.res" ]; then echo "$(date +%T) [$RUN] raw: already done (skip, nothing overwritten)"; else
    echo "$(date +%T) [$RUN] KMA ResFinder on RAW reads start"
    mamba run -n wwamr-kma kma -ipe "$R1" "$R2" -o "$OUT" -t_db "$IDX/resfinder_2.1.1" -mem_mode -ef -1t1 -cge -nf -vcf -t ${THREADS:-4} \
      > logs/kma_explore_${RUN}_raw.log 2>&1 || { echo "ERROR: KMA failed, see logs/kma_explore_${RUN}_raw.log"; exit 1; }
    echo "$(date +%T) [$RUN] done"
  fi
  echo "--- $RUN raw reads: erm(B) templates and neighbours (readCount, covered positions, readCountAln)   [registered, trimmed run in brackets]"
  for pre in "erm(B)_" "erm(G)_" "erm(F)_" "erm(Q)_"; do
    awk -F'\t' -v p="$pre" 'substr($1,1,length(p))==p && $2>=50 {printf "  %-22s readCount %6d  covPos %5d  aln %4d\n",$1,$2,$5,$14}' "$OUT.mapstat" | head -6
  done
  echo "  [trimmed, registered run:]"
  for pre in "erm(B)_" "erm(G)_" "erm(F)_" "erm(Q)_"; do
    awk -F'\t' -v p="$pre" 'substr($1,1,length(p))==p && $2>=50 {printf "  %-22s readCount %6d  covPos %5d  aln %4d\n",$1,$2,$5,$14}' "results/kma/${RUN}_resfinder_2.1.1.mapstat" | head -6
  done
done
