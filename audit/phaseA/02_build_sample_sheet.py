#!/usr/bin/env python3
"""Phase A step 3 (2026-10-27): draft sample sheet from the raw ENA table plus type/day from the paper's S4 (IDs, type, day only).
Usage: python3 02_build_sample_sheet.py <ena_filereport.tsv> <MOESM1.docx> <outdir>
Checks: replicate letter in ENA sample_title (A/B/C) versus replicate index by ID order (the S5 assumption).
Licence note: contains sample IDs/type/day from S4 (CC BY-NC-ND 4.0); do not publish this file without checking.
"""
import sys, csv, collections, zipfile, re, os
tsv, docx, out = sys.argv[1:4]
rows = {r["sample_alias"]: r for r in csv.DictReader(open(tsv), delimiter="\t") if r["library_source"] == "METAGENOMIC"}
x = zipfile.ZipFile(docx).read("word/document.xml").decode()
s4 = {}
for tr in re.findall(r"<w:tr[ >].*?</w:tr>", x, flags=re.S):
    f = [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<w:tc>.*?</w:tc>", tr, flags=re.S)]
    if len(f) == 4 and re.fullmatch(r"[REPS]\d{3}", f[0]): s4[f[0]] = (f[2], int(f[3]))
assert len(s4) == 33
grp = collections.defaultdict(list)
for i, (t, d) in s4.items(): grp[(t, d)].append(i)
sheet = []; bad = 0
for (t, d), ids in sorted(grp.items()):
    for n, i in enumerate(sorted(ids)):
        r = rows[i]; m = re.fullmatch(r"(\d{6})_([REPS])_([ABC])(\d)", r["sample_title"])
        letter = m.group(3) if m else "?"
        exp = "ABC"[n]
        if letter != exp: bad += 1
        sheet.append([r["run_accession"], r["sample_accession"], i, t, d, r["collection_date"], r["sample_title"], letter, n + 1,
                      "dev" if d == 1 else "eval", "sediment" if t == "SED" else "liquid", round(int(r["base_count"]) / 1e9, 2), r["fastq_bytes"], r["fastq_md5"]])
with open(os.path.join(out, "sample_sheet_draft.tsv"), "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["run_accession", "sample_accession", "sample_id", "type", "day", "collection_date", "ena_sample_title", "replicate_letter", "replicate_index_by_id_order", "split", "matrix", "gbp", "fastq_bytes", "fastq_md5"])
    w.writerows(sheet)
print("rows", len(sheet), "| split:", dict(collections.Counter((r[9], r[10]) for r in sheet)))
print("replicate letter disagrees with ID-order index:", bad, "of", len(sheet))
print("titles by type/day:", {k: [rows[i]["sample_title"] for i in sorted(v)] for k, v in sorted(grp.items())})
