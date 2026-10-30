#!/usr/bin/env python3
"""Phase A steps 2-4 (2026-10-27): classify PRJNA904380 runs, test mapping of the 33 paper sample IDs, size.
Usage: python3 01_classify_runs.py <ena_filereport.tsv> <MOESM3.xlsx> <outdir>
Reads the raw ENA table (never modified) and the sample IDs from the header of Table S6 (IDs only).
Writes: out/run_classification.tsv (one row per run with class and paper-ID match), out/classification_summary.txt
"""
import sys, csv, collections, openpyxl, os
tsv, s6, out = sys.argv[1:4]
rows = list(csv.DictReader(open(tsv), delimiter="\t"))
ws = openpyxl.load_workbook(s6, data_only=True).active
vals = [[c.value for c in r] for r in ws.iter_rows()]
hdr = [i for i, r in enumerate(vals) if r and r[0] == "gene"][0]
paper_ids = [x for x in vals[hdr][2:35]]
assert len(paper_ids) == 33
L = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); L.append(s)
p("ENA rows:", len(rows), "| paper IDs (S6 header):", len(paper_ids))
def cls(r):
    if r["library_source"] == "METAGENOMIC": return "metagenome"
    if r["instrument_model"] == "MinION": return "isolate_nanopore"
    return "isolate_illumina"
for r in rows: r["class"] = cls(r)
p("class counts:", dict(collections.Counter(r["class"] for r in rows)))
meta = [r for r in rows if r["class"] == "metagenome"]
p("metagenome platform:", dict(collections.Counter(r["instrument_model"] for r in meta)),
  "| layout:", dict(collections.Counter(r["library_layout"] for r in meta)))
# map by sample_alias
by_alias = collections.defaultdict(list)
for r in meta: by_alias[r["sample_alias"]].append(r["run_accession"])
matched = [i for i in paper_ids if i in by_alias]
p("paper IDs found as sample_alias among metagenomes:", len(matched), "of 33")
p("paper IDs NOT found:", [i for i in paper_ids if i not in by_alias])
extra = sorted(a for a in by_alias if a not in paper_ids)
p("metagenome aliases not in the paper's 33:", extra)
for a in extra:
    for r in meta:
        if r["sample_alias"] == a:
            p("  extra:", a, r["run_accession"], r["sample_title"], "| base_count", r["base_count"], "| collection_date", r["collection_date"], "| desc:", r["sample_description"][:80])
dups = {a: v for a, v in by_alias.items() if len(v) > 1}
p("aliases with more than one run:", dups)
# size / depth of the 33
sel = [r for r in meta if r["sample_alias"] in paper_ids]
gb = sum(int(r["base_count"]) for r in sel) / 1e9
fb = sum(sum(int(x) for x in r["fastq_bytes"].split(";") if x) for r in sel) / 1e9
p("33 mapped runs: total bases %.1f Gbp | FASTQ.gz bytes %.1f GB" % (gb, fb))
bc = [int(r["base_count"]) / 1e9 for r in sel]
p("per-sample Gbp: min %.1f max %.1f" % (min(bc), max(bc)))
rc = [int(r["read_count"]) / 1e6 for r in sel]
p("per-sample reads (M, ENA read_count): min %.1f max %.1f total %.0f" % (min(rc), max(rc), sum(rc)))
alle = sum(int(r["base_count"]) for r in meta) / 1e9
p("all 39 metagenome runs: %.1f Gbp" % alle)
with open(os.path.join(out, "run_classification.tsv"), "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["run_accession", "sample_accession", "sample_alias", "class", "in_paper_33", "instrument_model", "base_count", "read_count", "fastq_bytes"])
    for r in rows:
        w.writerow([r["run_accession"], r["sample_accession"], r["sample_alias"], r["class"], r["sample_alias"] in paper_ids and r["class"] == "metagenome", r["instrument_model"], r["base_count"], r["read_count"], r["fastq_bytes"]])
open(os.path.join(out, "classification_summary.txt"), "w").write("\n".join(L) + "\n")
