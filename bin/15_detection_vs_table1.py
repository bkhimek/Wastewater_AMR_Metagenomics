#!/usr/bin/env python3
"""
15_detection_vs_table1.py [--out results/detection_vs_table1.tsv]

DESCRIPTIVE (not part of the R1 verdict). Metagenomic detection of the five genes per sample type, from our KMA runs,
compared with Table 1 of Taylor et al. 2025 (MGS rows, % of samples where the ARG was detected).
A gene counts as detected in a sample when at least one ResFinder 2.1.1 template of that gene family has KMA readCount >= 50
(the same >= 50 rule as the registered pipeline; family = base names of the published subtypes, as in 09/11).
Sample types and sizes: INF 9, EFF 9, POND 9 (Day 1 = dev, Days 2-3 = eval) and SED 6 (eval, sediment, descriptive).
Standard library only. Run from the project root.
"""
import csv, sys
from pathlib import Path

KMA = Path("results/kma")
GENES = ["ermB", "sul1", "tetA", "tetQ", "tetW"]
TYPES = ["INF", "EFF", "POND", "SED"]
# Table 1, MGS rows (% of samples where the ARG was detected), as printed in the paper
TABLE1_MGS = {
    "INF": [100, 100, 100, 100, 100], "EFF": [100, 100, 100, 100, 100],
    "POND": [0, 100, 77.77, 0, 0], "SED": [0, 100, 66.6, 0, 16.6],
}
OUT = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else Path("results/detection_vs_table1.tsv")


def counts(path):
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        out[p[0]] = float(p[1])
    return out


def main():
    rd = list(csv.reader(open("refdata/paper_S6_subtype_read_counts.csv")))
    fam = {g: set() for g in GENES}
    for r in rd[1:]:
        fam[r[0]].add(r[1].split("_")[0])
    sheet = list(csv.DictReader(open("phaseA_audit/out/sample_sheet_draft.tsv"), delimiter="\t"))
    det = {t: {g: [] for g in GENES} for t in TYPES}
    for s in sheet:
        f = KMA / f"{s['run_accession']}_resfinder_2.1.1.mapstat"
        if not f.exists() or f.stat().st_size == 0:
            continue
        c = counts(f)
        for g in GENES:
            tot = sum(v for n, v in c.items() if n.split("_")[0] in fam[g] and v >= 50)
            det[s["type"]][g].append((s["sample_id"], tot))
    rows = []
    print("\nMGS detection (gene detected = any template of the gene with readCount >= 50): ours vs Table 1 of the paper")
    print(f"{'type':6s}{'n':>3s}" + "".join(f"{g:>24s}" for g in GENES))
    print(f"{'':9s}" + "".join(f"{'ours n/N  %  | paper %':>24s}" for _ in GENES))
    for t in TYPES:
        n = len(det[t][GENES[0]])
        line = f"{t:6s}{n:>3d}"
        for i, g in enumerate(GENES):
            k = sum(1 for _, v in det[t][g] if v > 0)
            pct = 100.0 * k / n if n else float("nan")
            line += f"{f'{k}/{n}  {pct:5.1f}  | {TABLE1_MGS[t][i]:g}':>24s}"
            rows.append([t, g, n, k, round(pct, 2), TABLE1_MGS[t][i]])
        print(line)
    print("\nPer-sample gene-family totals (readCount sum of templates >= 50) for POND and SED:")
    for t in ("POND", "SED"):
        print(f"  {t}")
        ids = [sid for sid, _ in det[t][GENES[0]]]
        print("    " + f"{'sample':8s}" + "".join(f"{g:>8s}" for g in GENES))
        for j, sid in enumerate(ids):
            print("    " + f"{sid:8s}" + "".join(f"{det[t][g][j][1]:8.0f}" for g in GENES))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["type", "gene", "n_samples", "n_detected", "pct_ours", "pct_table1"]); w.writerows(rows)
    print(f"\nTable written to {OUT}")


if __name__ == "__main__":
    main()
