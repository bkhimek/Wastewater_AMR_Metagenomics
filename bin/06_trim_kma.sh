#!/usr/bin/env bash
# 06_trim_kma.sh <run_accession> - reproduction arm for ONE development sample:
#   1) bbduk 38.90 with the authors' stated parameters, writing PAIRED trimmed reads (unpaired reads dropped)
#   2) KMA 1.4.9 against ResFinder db tag 2.1.1 and SILVA 138.1 NR99 (DNA) with the authors' flags
#      -ipe -mem_mode -ef -1t1 -cge -nf -vcf   (supplement MOESM1, "as in Munk et al.")
# Assumptions (to be stated in the registration): bbduk adapter options (ref=adapters ktrim=r k=23 mink=11 hdist=1 tbo tpe),
# SILVA variant NR99, KMA index built with defaults from all.fsa / the DNA-converted fasta.
# Idempotent: finished steps are skipped. Only development runs are accepted (evaluation data stay closed).
# Usage: bash bin/06_trim_kma.sh SRR22752025        (THREADS=4; KEEP_TRIMMED=1 keeps the trimmed reads)
set -uo pipefail
RUN=${1:?usage: 06_trim_kma.sh <run_accession>}
REPO=${REPO:-$HOME/projects/Wastewater_AMR_Metagenomics}
DEV="SRR22752025 SRR22752015 SRR22752048 SRR22752047 SRR22752038 SRR22752037 SRR22752014 SRR22752036 SRR22752035"
case " $DEV " in *" $RUN "*) ;; *) echo "ERROR: $RUN is not in the development list"; exit 1;; esac
THREADS=${THREADS:-4}
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/miniconda3}
RAW=$REPO/data/fastq_raw; TRIM=$REPO/data/trimmed; KM=$REPO/results/kma; IDX=$REPO/refdb/kma_index
mkdir -p "$TRIM" "$KM" "$REPO/logs"
R1=$RAW/${RUN}_1.fastq.gz; R2=$RAW/${RUN}_2.fastq.gz
for f in "$R1" "$R2"; do
  [ -f "$f" ] && grep -q -F "$(basename "$f")" "$RAW/verified.tsv" || { echo "ERROR: $f not verified"; exit 1; }
done
T1=$TRIM/${RUN}_1.trim.fq.gz; T2=$TRIM/${RUN}_2.trim.fq.gz; BL=$REPO/logs/bbduk_${RUN}.log
KEEP_TRIMMED=${KEEP_TRIMMED:-0}   # 0 (default): delete the trimmed reads after both KMA runs (about 11 GB per sample; they can be regenerated)
if [ -s "$KM/${RUN}_resfinder_2.1.1.res" ] && [ -s "$KM/${RUN}_silva_138.1_NR99.res" ]; then echo "$(date +%T) [$RUN] all done (skip)"; exit 0; fi
if [ -s "$T1" ] && [ -s "$T2" ] && grep -q '^Result:' "$BL" 2>/dev/null; then echo "$(date +%T) [$RUN] bbduk: skip (done)"; else
  echo "$(date +%T) [$RUN] bbduk start"
  rm -f "$T1.tmp" "$T2.tmp"
  mamba run -n wwamr-bbduk bbduk.sh in="$R1" in2="$R2" out="$T1.tmp" out2="$T2.tmp" \
    ref=adapters ktrim=r k=23 mink=11 hdist=1 tbo tpe maxns=0 minlength=75 qtrim=r trimq=20 \
    threads=$THREADS -Xmx8g 2> "$BL" || { echo "ERROR: bbduk failed, see $BL"; exit 1; }
  mv "$T1.tmp" "$T1" && mv "$T2.tmp" "$T2"
  grep -E '^(Input|Result|Total Removed):' "$BL"
fi
for DB in resfinder_2.1.1 silva_138.1_NR99; do
  OUT=$KM/${RUN}_${DB}
  if [ -s "$OUT.res" ]; then echo "$(date +%T) [$RUN] KMA $DB: skip (done)"; continue; fi
  echo "$(date +%T) [$RUN] KMA $DB start"
  mamba run -n wwamr-kma kma -ipe "$T1" "$T2" -o "$OUT" -t_db "$IDX/$DB" -mem_mode -ef -1t1 -cge -nf -vcf -t $THREADS \
    > "$REPO/logs/kma_${RUN}_${DB}.log" 2>&1 || { echo "ERROR: KMA $DB failed, see logs/kma_${RUN}_${DB}.log"; exit 1; }
  echo "$(date +%T) [$RUN] KMA $DB done: $(($(wc -l < "$OUT.res") - 1)) templates in .res"
done
if [ "$KEEP_TRIMMED" != 1 ]; then rm -f "$T1" "$T2"; echo "$(date +%T) [$RUN] trimmed reads deleted (KEEP_TRIMMED=1 keeps them)"; fi
