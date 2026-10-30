#!/usr/bin/env bash
# 04_network_diag.sh (2026-10-27): is ENA slow, or is the local connection / WSL network slow?  Run from WSL:  bash bin/04_network_diag.sh
# Each test downloads at most 30-50 MB and stops after 90 s. Compare the MB/s lines. Nothing is written except /dev/null.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TABLE="${ENA_TABLE:-$REPO/audit/phaseA/raw/ena_filereport_PRJNA904380.tsv}"
FTP=$(python3 - "$TABLE" <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    if r["sample_alias"] == "E001" and r["library_source"] == "METAGENOMIC": print(r["fastq_ftp"].split(";")[0])
PY
)
echo "curl: $(curl --version | head -1)"; echo "target: $FTP"
t() { # label, url, range-end
  local r; r=$(curl -fsSL -m 90 -r 0-"$3" -o /dev/null -w '%{size_download} %{speed_download} %{time_total}' "$2" 2>&1) || r="$r FAILED"
  echo "$1: $(echo "$r" | awk '{ if (NF>=3 && $1+0>0) printf "%.0f MB in %.1f s = %.2f MB/s", $1/1e6, $3, $2/1e6; else print $0 }')"
}
t "ENA https  " "https://$FTP" 31457279
t "ENA ftp    " "ftp://$FTP" 31457279
t "Cloudflare " "https://speed.cloudflare.com/__down?bytes=52428800" 52428799
t "ENA https again" "https://$FTP" 31457279
echo "--- WSL network"; ip -o link show | awk '{print $2, $3, $4, $5}' | grep -i -E "eth|mtu" | head -3
cat /etc/resolv.conf | grep -i nameserver | head -2
