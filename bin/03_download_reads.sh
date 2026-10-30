#!/usr/bin/env bash
# 03_download_reads.sh v5.3 (2026-10-28; a partial file that does not start with a gzip header (for example an 'Index of' web page returned by the http route) is removed; use the default https; after an HTTP error such as 403 wait 3 min before the next try; v5 + a file failing the md5 check is kept as <name>.bad instead of deleted, so the damage can be located with cmp): download FASTQ files from ENA in chunks (SET=...); verify each against ENA's md5 and byte size.
# Copy to <repo>/bin/ and run from WSL.  Modes:
#   plan   (default) list files and sizes, check free space; no download
#   test   one stream, first 50 MB, estimate total time
#   ptest  JOBS parallel streams, first 100 MB of different files, aggregate speed   (JOBS=4 bash bin/03_download_pilot.sh ptest)
#   run    download + verify, JOBS files at a time (default 1)                       (JOBS=4 nohup bash bin/03_download_pilot.sh run > logs/download_pilot.log 2>&1 &)
# Chunks (SET, default pilot), all derived from the ENA table only:
#   pilot        Day 1 replicates A and B of influent, effluent, pond (6 samples, 30.3 GB)
#   day1_rest    the other 3 Day 1 samples (replicate C)               (3 samples, 16.2 GB)
#   eval_liquid  liquid samples of Days 2 and 3                        (18 samples, 87.0 GB)  - evaluation data: download only, do not open before the freeze
#   sediment     the 6 sediment samples (all evaluation)               (6 samples, 31.4 GB)
#   all          all 33 samples                                        (164.9 GB)
#   or give your own list: SAMPLES="R001 R006 ..."   Example:  SET=day1_rest bash bin/03_download_reads.sh plan
# Only one run per output directory (lock file). PROTO=http uses http:// instead of https:// (md5 still verifies every file).
# Stop any time (kill the job); re-running continues partial files and skips verified ones.
# Safety: refuses to write under /mnt/c or /mnt/d (OneDrive); resumable; verified files are skipped. Internal mode: one <file>.
set -euo pipefail
MODE="${1:-plan}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TABLE="${ENA_TABLE:-$REPO/audit/phaseA/raw/ena_filereport_PRJNA904380.tsv}"
OUT="${DATA_DIR:-$REPO/data/fastq_raw}"
PROTO="${PROTO:-https}"; SET="${SET:-pilot}"; SAMPLES="${SAMPLES:-${PILOT_SAMPLES:-}}"
JOBS="${JOBS:-1}"
[ "$MODE" = one ] || echo "pwd: $(pwd) | repo: $REPO | mode: $MODE | jobs: $JOBS"
[ -f "$TABLE" ] || { echo "ERROR: ENA table not found: $TABLE"; exit 1; }
mkdir -p "$OUT"; OUTR="$(readlink -f "$OUT")"
case "$OUTR" in /mnt/c/*|/mnt/d/*) echo "ERROR: output $OUTR is on a Windows drive (OneDrive). Use the WSL repo (copy this script to <repo>/bin/)."; exit 1;; esac
if [ "$MODE" = run ]; then
  exec 9>"$OUT/.run.lock"
  flock -n 9 || { echo "ERROR: another download is already running for $OUT (see: ps -eo pid,etime,args | grep 03_download). Stop it first."; exit 1; }
fi
PLAN="$OUT/download_plan.tsv"; VER="$OUT/verified.tsv"; touch "$VER"

if [ "$MODE" != one ]; then
python3 - "$TABLE" "$SAMPLES" "$SET" "$PROTO" > "$PLAN" <<'PY'
import csv, re, sys, collections
rows = {r["sample_alias"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t") if r["library_source"] == "METAGENOMIC"}
paper = sorted(a for a in rows if re.fullmatch(r"[REPS]\d{3}", a) and not a.endswith("000"))
assert len(paper) == 33, len(paper)
rank = {}
for L in "REP":
    ds = sorted({rows[a]["collection_date"] for a in paper if a[0] == L})
    rank.update({a: ds.index(rows[a]["collection_date"]) + 1 for a in paper if a[0] == L})
pilot = "E001 E006 R001 R006 P001 P006".split()
day1 = [a for a in paper if a[0] in "REP" and rank[a] == 1]
sets = {"pilot": pilot, "day1_rest": [a for a in day1 if a not in pilot],
        "eval_liquid": [a for a in paper if a[0] in "REP" and rank[a] > 1],
        "sediment": [a for a in paper if a[0] == "S"], "all": paper}
chosen = sys.argv[2].split() if sys.argv[2].strip() else sets[sys.argv[3]]
print("\t".join(["alias", "run", "idx", "url", "bytes", "md5", "fname"]))
for a in chosen:
    r = rows[a]  # KeyError = alias not found among metagenomic runs
    urls, sizes, md5s = r["fastq_ftp"].split(";"), r["fastq_bytes"].split(";"), r["fastq_md5"].split(";")
    assert len(urls) == len(sizes) == len(md5s) == 2, (a, urls)
    for i, (u, s_, m) in enumerate(zip(urls, sizes, md5s), 1):
        print("\t".join([a, r["run_accession"], str(i), sys.argv[4] + "://" + u, s_, m, u.rsplit("/", 1)[1]]))
PY
NF=$(( $(wc -l < "$PLAN") - 1 )); TOT=$(awk -F'\t' 'NR>1{s+=$5} END{print s+0}' "$PLAN")
AVAIL=$(df --output=avail -B1 "$OUT" | tail -1)
echo "files: $NF | total: $(awk -v t=$TOT 'BEGIN{printf "%.1f", t/1e9}') GB | free on target: $(awk -v t=$AVAIL 'BEGIN{printf "%.0f", t/1e9}') GB"
[ "${SKIP_SPACE_CHECK:-0}" = 1 ] || [ "$AVAIL" -gt $(( TOT * 12 / 10 )) ] || { echo "ERROR: less than 1.2 x total free"; exit 1; }
fi

download_one() {
  local f="$1" row alias run idx url bytes md5 fname SZ SUM
  row=$(awk -F'\t' -v f="$f" '$7==f' "$PLAN"); [ -n "$row" ] || { echo "ERROR: $f not in plan"; return 1; }
  IFS=$'\t' read -r alias run idx url bytes md5 fname <<< "$row"
  if [ -f "$OUT/$fname" ] && grep -q -F -- "$fname"$'\t'"$md5" "$VER"; then echo "$(date -u +%T) [$fname] skip (verified)"; return 0; fi
  echo "$(date -u +%T) [$fname] download $(awk -v b=$bytes 'BEGIN{printf "%.2f", b/1e9}') GB"
  local a=0 cur rc W=20
  while :; do
    cur=$(stat -c %s "$OUT/$fname.part" 2>/dev/null || echo 0)
    [ "$cur" = "$bytes" ] && break
    [ "$cur" -gt "$bytes" ] && { echo "[$fname] partial file larger than expected; removed"; rm -f "$OUT/$fname.part"; }
    [ "$a" -ge 40 ] && { echo "[$fname] FAILED after $a attempts (re-run to resume)"; return 1; }
    a=$((a+1)); [ "$a" -gt 1 ] && { echo "$(date -u +%T) [$fname] retry $a (resume at $cur bytes; waited ${W:-20} s)"; sleep "${W:-20}"; }
    rc=0; curl -fsSL --connect-timeout 30 --speed-time 120 --speed-limit 20000 -C - -o "$OUT/$fname.part" "$url" || rc=$?
    if [ "$rc" = 22 ]; then W=180; else W=20; fi   # curl 22 = HTTP error (e.g. 403): wait 3 min instead of hammering the server
    if [ -s "$OUT/$fname.part" ] && [ "$(head -c 2 "$OUT/$fname.part" | od -An -tx1 | tr -d ' \n')" != "1f8b" ]; then
      echo "[$fname] partial file does not start with a gzip header (a web page?); removed"; rm -f "$OUT/$fname.part"; W=60; fi
  done
  SZ=$(stat -c %s "$OUT/$fname.part")
  [ "$SZ" = "$bytes" ] || { echo "[$fname] SIZE MISMATCH: got $SZ expected $bytes (partial kept; re-run to resume)"; return 1; }
  SUM=$(md5sum "$OUT/$fname.part" | cut -d' ' -f1)
  [ "$SUM" = "$md5" ] || { echo "[$fname] MD5 MISMATCH: got $SUM expected $md5 (kept as $fname.bad for diagnosis; will be downloaded again)"; mv -f "$OUT/$fname.part" "$OUT/$fname.bad"; return 1; }
  mv "$OUT/$fname.part" "$OUT/$fname"; printf '%s\t%s\t%s\t%s\t%s\n' "$fname" "$md5" "$bytes" "$run" "$(date -u +%FT%TZ)" >> "$VER"
  echo "$(date -u +%T) [$fname] verified"
}

case "$MODE" in
  plan) cut -f1,2,3,5,7 "$PLAN" | column -t ;;
  test)
    URL=$(awk -F'\t' 'NR==2{print $4}' "$PLAN")
    curl -fsSL -r 0-52428799 -o /dev/null -w '%{size_download} %{speed_download}\n' "$URL" > "$OUT/.speedtest"
    read -r SZ SP < "$OUT/.speedtest"
    awk -v sz="$SZ" -v sp="$SP" -v t="$TOT" 'BEGIN{printf "downloaded %.0f MB at %.2f MB/s; estimated total %.1f hours for %.1f GB\n", sz/1e6, sp/1e6, t/sp/3600, t/1e9}' ;;
  ptest)
    t0=$(date +%s.%N); i=0
    while IFS=$'\t' read -r alias run idx url rest; do
      i=$((i+1)); if [ "$i" -gt "$JOBS" ]; then break; fi
      ( curl -fsSL -r 0-104857599 -o /dev/null -w '%{size_download}\n' "$url" > "$OUT/.pt$i" ) &
    done < <(tail -n +2 "$PLAN")
    wait; t1=$(date +%s.%N)
    B=$(cat "$OUT"/.pt* | awk '{s+=$1} END{print s+0}'); rm -f "$OUT"/.pt*
    awk -v b="$B" -v t0="$t0" -v t1="$t1" -v n="$JOBS" -v t="$TOT" 'BEGIN{r=b/(t1-t0); printf "%d parallel streams: %.0f MB in %.1f s = %.2f MB/s aggregate; estimated %.1f hours for %.1f GB\n", n, b/1e6, t1-t0, r/1e6, t/r/3600, t/1e9}' ;;
  one) download_one "$2"; exit $? ;;
  run)
    export DATA_DIR="$OUT" ENA_TABLE="$TABLE"; BAD=0
    tail -n +2 "$PLAN" | cut -f7 | xargs -P "$JOBS" -I{} bash "$0" one {} || BAD=1
    echo "done: $(wc -l < "$VER") verified entries of $NF; problems: $BAD"; [ "$BAD" -eq 0 ] ;;
  *) echo "ERROR: unknown mode $MODE"; exit 1 ;;
esac
