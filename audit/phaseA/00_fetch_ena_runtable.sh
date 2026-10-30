#!/usr/bin/env bash
# Phase A step 1 (2026-10-27): fetch the ENA run table for PRJNA904380 and keep it UNMODIFIED.
# Run from WSL:  bash /mnt/c/Users/krist/OneDrive/Documents/Projects/metagenomics/Wastewater_AMR_Metagenomics/phaseA_audit/00_fetch_ena_runtable.sh
# Output goes next to this script (raw/). Nothing is downloaded except metadata text (no reads).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pwd
RAW="$HERE/raw"; mkdir -p "$RAW"
ACC="PRJNA904380"
API="https://www.ebi.ac.uk/ena/portal/api"
STAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# 1. Field names ENA accepts for read_run (saved, so the field list below can be checked against it)
curl -fsS "$API/returnFields?dataPortal=ena&result=read_run&format=tsv" -o "$RAW/ena_returnFields_read_run.tsv"
echo "returnFields rows: $(wc -l < "$RAW/ena_returnFields_read_run.tsv")"

# 2. Run table
FIELDS="study_accession,secondary_study_accession,run_accession,experiment_accession,sample_accession,secondary_sample_accession,submission_accession,run_alias,experiment_alias,sample_alias,sample_title,sample_description,experiment_title,study_title,library_name,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,scientific_name,tax_id,base_count,read_count,fastq_bytes,fastq_md5,fastq_ftp,first_public,collection_date,country,isolation_source,host,strain"
URL="$API/filereport?accession=$ACC&result=read_run&fields=$FIELDS&format=tsv&download=true"
curl -fsS "$URL" -o "$RAW/ena_filereport_${ACC}.tsv" || { echo "FAILED: a field name may be invalid; compare FIELDS with raw/ena_returnFields_read_run.tsv"; exit 1; }

# 3. Provenance record
{
  echo "accession: $ACC"; echo "fetched_utc: $STAMP"; echo "url: $URL"
  echo "rows_excluding_header: $(($(wc -l < "$RAW/ena_filereport_${ACC}.tsv") - 1))"
  echo "sha256: $(sha256sum "$RAW/ena_filereport_${ACC}.tsv" | cut -d' ' -f1)"
} | tee "$RAW/ena_filereport_${ACC}.provenance.txt"
echo "Done. Raw table is kept unmodified; later scripts read it and write elsewhere."
