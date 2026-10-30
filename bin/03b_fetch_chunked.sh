#!/usr/bin/env bash
# 03b_fetch_chunked.sh FILE_NAME   e.g. SRR22752036_1.fastq.gz
# Looks up URL, size and md5 at ENA, downloads with chunk_fetch.sh (no resume splicing), verifies md5,
# and only then places the file in data/fastq_raw/. Run from the project root (~/projects/Wastewater_AMR_Metagenomics).
set -uo pipefail
F=${1:?file name}; RUN=${F%%_*}; OUT=data/fastq_raw
J=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp,fastq_bytes,fastq_md5&format=tsv") || { echo "ENA lookup failed"; exit 1; }
row=$(echo "$J" | tail -n1)
IFS=$'\t' read -r _ ftp bytes md5 <<< "$row"
for k in 1 2 3; do
  u=$(echo "$ftp" | cut -d';' -f$k)
  [ "${u##*/}" = "$F" ] && { URL="https://$u"; BYTES=$(echo "$bytes" | cut -d';' -f$k); MD5=$(echo "$md5" | cut -d';' -f$k); }
done
[ -n "${URL:-}" ] || { echo "file $F not found in ENA record"; exit 1; }
echo "$F  $BYTES bytes  md5 $MD5"
bash bin/chunk_fetch.sh "$URL" "$BYTES" "$OUT/$F.chunked" || exit 1
got=$(md5sum "$OUT/$F.chunked" | cut -d' ' -f1)
if [ "$got" = "$MD5" ]; then mv "$OUT/$F.chunked" "$OUT/$F"; echo "VERIFIED $F"; else mv "$OUT/$F.chunked" "$OUT/$F.bad_chunked_$(date +%H%M%S)"; echo "MD5 MISMATCH $F (got $got)"; exit 1; fi
