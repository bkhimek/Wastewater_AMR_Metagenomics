#!/usr/bin/env python3
"""
07_compare_paper_counts.py RUN SAMPLE_ID [PAPER_CSV]

Compares our KMA output for one sample with the authors' published subtype read counts
(Supplementary Table S6, saved as CSV) and computes the ARG / 16S ratio under several
possible definitions. Reads only files in results/kma/ (run from the project root).
Standard library only. Development use: run it on Day-1 development samples only.

  python3 bin/07_compare_paper_counts.py SRR22752025 E001
"""
import csv, sys, statistics as st
from pathlib import Path

RUN = sys.argv[1]
SAMPLE = sys.argv[2]
PAPER = Path(sys.argv[3] if len(sys.argv) > 3 else "refdata/paper_S6_subtype_read_counts.csv")
KMA = Path("results/kma")
GENES = ["ermB", "sul1", "tetA", "tetQ", "tetW"]
COUNTS = ["readCount", "fragmentCount", "readCountAln", "fragmentCountAln"]
MIN_COUNT = 50          # the supplement filters subtype counts at 50 (smallest published non-zero value is 50)
COV_CUTOFFS = [None, 90.0, 98.98]


def read_res(path):
    out = {}
    with open(path) as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        header = [h.strip() for h in header]
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            row = dict(zip(header, [x.strip() for x in p]))
            out[row["Template"]] = row
    return out


def read_mapstat(path):
    out, header = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                header = [h.strip() for h in header]
                continue
            p = line.rstrip("\n").split("\t")
            row = dict(zip(header, p))
            out[row["refSequence"]] = row
    return out


def fnum(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def load(db):
    res = read_res(KMA / f"{RUN}_{db}.res")
    ms = read_mapstat(KMA / f"{RUN}_{db}.mapstat")
    merged = {}
    for t, r in res.items():
        m = ms.get(t, {})
        merged[t] = dict(
            cov=fnum(r["Template_Coverage"]), ident=fnum(r["Template_Identity"]),
            qid=fnum(r["Query_Identity"]), depth=fnum(r["Depth"]),
            **{c: fnum(m.get(c)) for c in COUNTS})
    return merged


# ---- paper table ---------------------------------------------------------
with open(PAPER) as f:
    rd = list(csv.reader(f))
hdr, rows = rd[0], rd[1:]
if SAMPLE not in hdr:
    sys.exit(f"sample {SAMPLE} not found in {PAPER}; columns: {hdr[2:35]}")
ci = hdr.index(SAMPLE)
paper = {r[1]: (r[0], int(float(r[ci] or 0))) for r in rows}

arg = load("resfinder_2.1.1")
s16 = load("silva_138.1_NR99")

print(f"\nRun {RUN} = sample {SAMPLE}; paper file {PAPER}")
print(f"ResFinder templates in .res: {len(arg)}; SILVA templates: {len(s16)}\n")

# ---- Table A: subtype by subtype -----------------------------------------
print("TABLE A  authors' published count (Table S6) vs our KMA counts, per ARG subtype")
print(f"{'subtype':22s}{'paper':>7s}{'readCnt':>9s}{'fragCnt':>9s}{'readAln':>9s}{'fragAln':>9s}{'cov%':>7s}{'id%':>7s}{'depth':>8s}")
for sub, (gene, pc) in paper.items():
    a = arg.get(sub)
    if pc == 0 and (a is None or max(a[c] for c in COUNTS) < 1):
        continue
    if a is None:
        print(f"{sub:22s}{pc:7d}   (not in our .res)")
        continue
    print(f"{sub:22s}{pc:7d}" + "".join(f"{a[c]:9.0f}" for c in COUNTS) +
          f"{a['cov']:7.1f}{a['ident']:7.1f}{a['depth']:8.1f}")

# ---- Table B: which count column reproduces the published numbers? --------
print("\nTABLE B  which of our counts reproduces the published subtype counts?")
print("(over the 28 published subtypes; 'agree' = our count >= 50 exactly when the paper lists a count)")
print(f"{'our column':18s}{'agree':>8s}{'median our/paper':>20s}{'spearman-like (rank corr)':>28s}")
def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
    for k, i in enumerate(o): r[i] = k
    return r
for c in COUNTS:
    agree, ratios, xs, ys = 0, [], [], []
    for sub, (gene, pc) in paper.items():
        a = arg.get(sub); ours = a[c] if a else 0.0
        agree += (ours >= MIN_COUNT) == (pc > 0)
        xs.append(pc); ys.append(ours)
        if pc > 0 and ours > 0:
            ratios.append(ours / pc)
    rx, ry = rank(xs), rank(ys)
    try:
        corr = st.correlation(rx, ry)
    except Exception:
        corr = float("nan")
    med = st.median(ratios) if ratios else float("nan")
    print(f"{c:18s}{agree:>5d}/28{med:>20.2f}{corr:>28.2f}")

# ---- Table C: ARG / 16S under different definitions -----------------------
print("\nTABLE C  total of the 5 genes of interest divided by 16S, under different definitions")
paper_total = sum(pc for (g, pc) in paper.values())
print(f"paper: sum of published subtype counts for {SAMPLE} = {paper_total} (16S counts are not published)")
print(f"{'count':18s}{'ARG (subtype>=50)':>19s}" + "".join(
    f"{('16S all' if c is None else '16S cov>=%g' % c):>17s}" for c in COV_CUTOFFS))
for c in COUNTS:
    num = sum(arg[s][c] for s in paper if s in arg and arg[s][c] >= MIN_COUNT)
    line = f"{c:18s}{num:19.0f}"
    for cut in COV_CUTOFFS:
        den = sum(t[c] for t in s16.values() if cut is None or t["cov"] >= cut)
        line += f"{(num / den if den else float('nan')):17.4f}"
    print(line)
print("\nDenominator sums (16S):")
for c in COUNTS:
    print(f"  {c:18s}" + "".join(
        f"  {('all' if cut is None else 'cov>=%g' % cut)}: {sum(t[c] for t in s16.values() if cut is None or t['cov'] >= cut):,.0f}"
        for cut in COV_CUTOFFS))
print("\nNOTE: development use only; numbers are for deciding the counting definition, not results.")
