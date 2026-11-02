#!/usr/bin/env bash
# 14_explore_single_template.sh TEMPLATE RUN [RUN ...]
# EXPLORATORY (post-freeze, not part of the R1 verdict; review step C, variant "single-template index"):
# builds a KMA index that contains ONLY one ResFinder 2.1.1 template (so no other template can compete for reads) and
# maps the RAW reads of each run with the same KMA flags as the registered pipeline. The only difference to the registered
# run: the database (one template instead of all of ResFinder). Raw reads are used because trimming was shown not to
# change the E021 outcome. Output: results/explore/, registered outputs are never touched.
# Usage: bash bin/14_explore_single_template.sh erm\(B\)_2_K00551 SRR22752046 SRR22752045
set -uo pipefail
T=${1:?template name, e.g. 'erm(B)_2_K00551'}; shift
REPO=${REPO:-$HOME/projects/Wastewater_AMR_Metagenomics}; cd "$REPO" || exit 1
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/miniconda3}
FSA=refdb/resfinder_db_2.1.1/all.fsa; OUTD=results/explore; mkdir -p "$OUTD" logs
SAFE=$(echo "$T" | tr -c 'A-Za-z0-9_\n' '_'); IDX=$OUTD/idx_single_$SAFE
if [ ! -s "$IDX.name" ]; then
  awk -v t="$T" '/^>/{keep=(substr($0,2)==t)} keep' "$FSA" > "$IDX.fsa"
  [ -s "$IDX.fsa" ] || { echo "ERROR: template $T not found in $FSA"; exit 1; }
  echo "$(date +%T) index with one template: $(grep -c '>' "$IDX.fsa") sequence(s), $(grep -v '>' "$IDX.fsa" | tr -d '\n' | wc -c) bp"
  mamba run -n wwamr-kma kma index -i "$IDX.fsa" -o "$IDX" > logs/kma_explore_index_$SAFE.log 2>&1 || { echo "ERROR: index failed"; exit 1; }
fi
for RUN in "$@"; do
  R1=data/fastq_raw/${RUN}_1.fastq.gz; R2=data/fastq_raw/${RUN}_2.fastq.gz
  grep -q -F "${RUN}_1.fastq.gz" data/fastq_raw/verified.tsv && grep -q -F "${RUN}_2.fastq.gz" data/fastq_raw/verified.tsv || { echo "ERROR: $RUN not verified"; exit 1; }
  OUT=$OUTD/${RUN}_raw_single_$SAFE
  if [ ! -s "$OUT.res" ]; then
    echo "$(date +%T) [$RUN] KMA on raw reads, single-template index, start"
    mamba run -n wwamr-kma kma -ipe "$R1" "$R2" -o "$OUT" -t_db "$IDX" -mem_mode -ef -1t1 -cge -nf -vcf -t ${THREADS:-4} \
      > logs/kma_explore_${RUN}_single_$SAFE.log 2>&1 || { echo "ERROR: KMA failed"; exit 1; }
  fi
  echo "--- $RUN: $T alone in the database:"
  awk -F'\t' '!/^#/ {printf "  readCount %6d  fragmentCount %6d  covPos %5d  readCountAln %4d\n",$2,$3,$5,$14}' "$OUT.mapstat"
  echo "  (registered run, full ResFinder, trimmed: $(awk -F'\t' -v t="$T" '$1==t {printf "readCount %d covPos %d aln %d", $2,$5,$14}' results/kma/${RUN}_resfinder_2.1.1.mapstat))"
done
