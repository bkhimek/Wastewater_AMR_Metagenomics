#!/usr/bin/env bash
# 03c_fetch_set.sh SPLIT MATRIX      e.g.  bash bin/03c_fetch_set.sh eval liquid
# Downloads every FASTQ file of one set of the sample sheet with 03b_fetch_chunked.sh (chunked, size and md5 checked
# against ENA), skips files already listed in data/fastq_raw/verified.tsv, appends each verified file to verified.tsv,
# and makes up to ROUNDS passes (default 3) over the files that failed. Transfer and checksum only: the reads are
# never opened, counted or QC'd. One job at a time (flock). Run from the project root, preferably inside tmux.
set -uo pipefail
SPLIT=${1:?split (dev|eval)}; MATRIX=${2:?matrix (liquid|sediment)}; ROUNDS=${ROUNDS:-3}
SHEET=phaseA_audit/out/sample_sheet_draft.tsv; RAW=data/fastq_raw; VER=$RAW/verified.tsv
mkdir -p logs "$RAW"; touch "$VER"
exec 9>"$RAW/.fetch.lock"; flock -n 9 || { echo "another fetch job is running"; exit 1; }
LOG=logs/fetch_${SPLIT}_${MATRIX}.log; exec > >(tee -a "$LOG") 2>&1

# list: the paired files only (RUN_1 and RUN_2); unpaired files are not used by the registered pipeline.
# md5 is checked against the live ENA record inside 03b_fetch_chunked.sh (the file is moved into place only if it matches).
LIST=$(awk -F'\t' -v s="$SPLIT" -v m="$MATRIX" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
  $h["split"]==s && $h["matrix"]==m {print $h["run_accession"]"_1.fastq.gz\t"$h["run_accession"]; print $h["run_accession"]"_2.fastq.gz\t"$h["run_accession"]}' "$SHEET")
TOTAL=$(echo "$LIST" | wc -l); echo "$(date +%F\ %T) set $SPLIT $MATRIX: $TOTAL files"

is_ok() { grep -q "^$1	" "$VER"; }
record() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$(md5sum "$RAW/$1" | cut -d' ' -f1)" "$(stat -c %s "$RAW/$1")" "$2" "$(date +%FT%T)" >> "$VER"; }
for ((r=1; r<=ROUNDS; r++)); do
  bad=0
  while IFS=$'\t' read -r F RUN; do
    [ -n "$F" ] || continue
    is_ok "$F" && continue
    # a file under its final name exists only after 03b verified it against ENA and moved it: record it, do not fetch again
    if [ -s "$RAW/$F" ]; then record "$F" "$RUN"; echo "$(date +%T) recorded existing verified file $F"; continue; fi
    echo "$(date +%T) round $r: $F"
    if bash bin/03b_fetch_chunked.sh "$F" && [ -s "$RAW/$F" ]; then record "$F" "$RUN"
    else bad=$((bad+1)); echo "$(date +%T) FAILED $F (will retry in the next round)"; fi
  done <<< "$LIST"
  [ "$bad" = 0 ] && break
done
ok=$(echo "$LIST" | while IFS=$'\t' read -r F R; do is_ok "$F" && [ -s "$RAW/$F" ] && echo 1; done | wc -l)
echo "$(date +%F\ %T) DONE: $ok of $TOTAL files verified for $SPLIT $MATRIX"
[ "$ok" = "$TOTAL" ]
