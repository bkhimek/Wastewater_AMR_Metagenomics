#!/usr/bin/env python3
"""
08_scan_16S_definitions.py [PAPER_CSV]

Tries many possible rules for counting "16S reads" (the denominator of ARG per 16S) and checks
which rule reproduces the ratios read from the paper's Fig. 1c (MGS), Day-1 bars, on all
development samples that have KMA results. Run from the project root; standard library only.
DEVELOPMENT USE ONLY (Day-1 samples). The figure values are approximate readings by eye (+-0.1),
so a good rule must match the PATTERN (INF ~2.3, EFF ~2.0, POND ~0.5) and not a single number.

Numerator = sum of KMA readCount over the 28 published subtypes of the 5 genes, counts >= 50
(this reproduced the published counts in 07_compare_paper_counts.py).
"""
import csv, math, sys
from pathlib import Path

PAPER = Path(sys.argv[1] if len(sys.argv) > 1 else "refdata/paper_S6_subtype_read_counts.csv")
KMA = Path("results/kma")
COUNTS = ["readCount", "fragmentCount", "readCountAln", "fragmentCountAln"]

# run -> (sample id, ratio read from Fig. 1c, MGS panel, Day 1; approximate)
DEV = {
    "SRR22752048": ("R001", 2.25), "SRR22752047": ("R006", 2.45), "SRR22752036": ("R011", 2.30),
    "SRR22752025": ("E001", 1.95), "SRR22752015": ("E006", 2.05), "SRR22752014": ("E011", 2.10),
    "SRR22752038": ("P001", 0.40), "SRR22752037": ("P006", 0.55), "SRR22752035": ("P011", 0.55),
}


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def read_res(path):
    out = {}
    with open(path) as f:
        header = [h.strip() for h in f.readline().lstrip("#").rstrip("\n").split("\t")]
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= len(header):
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
                header = [h.strip() for h in line.lstrip("#").rstrip("\n").split("\t")]
                continue
            p = line.rstrip("\n").split("\t")
            out[p[0]] = dict(zip(header, p))
    return out


def load(run, db):
    res = read_res(KMA / f"{run}_{db}.res")
    ms = read_mapstat(KMA / f"{run}_{db}.mapstat")
    t = {}
    for name, r in res.items():
        m = ms.get(name, {})
        t[name] = dict(cov=fnum(r["Template_Coverage"]), ident=fnum(r["Template_Identity"]),
                       qid=fnum(r["Query_Identity"]), depth=fnum(r["Depth"]),
                       **{c: fnum(m.get(c)) for c in COUNTS})
    return t


# ---- candidate 16S rules -------------------------------------------------
def taxon_ok(name, mode):
    if mode == "all":
        return True
    if mode == "bacteria":
        return " Bacteria;" in name
    if mode == "bact_noorg":
        return " Bacteria;" in name and "Chloroplast" not in name and "Mitochondria" not in name
    if mode == "bact_arch":
        return (" Bacteria;" in name or " Archaea;" in name) and "Chloroplast" not in name and "Mitochondria" not in name
    raise ValueError(mode)


RULES = []
for tax in ["all", "bacteria", "bact_noorg"]:
    RULES.append((f"{tax}", lambda n, t, tax=tax: taxon_ok(n, tax)))
    for cut in [50, 80, 90, 95, 98.98, 99.5]:
        RULES.append((f"{tax} + Template_Coverage>={cut}", lambda n, t, tax=tax, cut=cut: taxon_ok(n, tax) and t["cov"] >= cut))
        RULES.append((f"{tax} + Template_Identity>={cut}", lambda n, t, tax=tax, cut=cut: taxon_ok(n, tax) and t["ident"] >= cut))
    for d in [1, 5, 10, 50, 100]:
        RULES.append((f"{tax} + Depth>={d}", lambda n, t, tax=tax, d=d: taxon_ok(n, tax) and t["depth"] >= d))
    for q in [95, 99]:
        RULES.append((f"{tax} + Query_Identity>={q}", lambda n, t, tax=tax, q=q: taxon_ok(n, tax) and t["qid"] >= q))
    RULES.append((f"{tax} + Coverage>=98.98 & Identity>=80.65", lambda n, t, tax=tax: taxon_ok(n, tax) and t["cov"] >= 98.98 and t["ident"] >= 80.65))
TOPN = [1, 5, 10, 20, 50, 100, 500]


def denominators(s16, col):
    """dict rule-name -> denominator for one sample and one count column"""
    out = {}
    for name, fn in RULES:
        out[name] = sum(t[col] for n, t in s16.items() if fn(n, t))
    vals = sorted((t[col] for t in s16.values()), reverse=True)
    for k in TOPN:
        out[f"top {k} templates by {col}"] = sum(vals[:k])
    return out


# ---- data ----------------------------------------------------------------
with open(PAPER) as f:
    rd = list(csv.reader(f))
subtypes = [r[1] for r in rd[1:]]

samples = []
for run, (sid, target) in DEV.items():
    need = [KMA / f"{run}_{d}.{e}" for d in ("resfinder_2.1.1", "silva_138.1_NR99") for e in ("res", "mapstat")]
    if all(p.exists() and p.stat().st_size > 0 for p in need):
        samples.append((run, sid, target))
if len(samples) < 2:
    sys.exit("need KMA results for at least 2 development samples")
print(f"samples used ({len(samples)}): " + ", ".join(f"{sid}" for _, sid, _ in samples))

num, den = {}, {}
for run, sid, target in samples:
    arg = load(run, "resfinder_2.1.1")
    num[sid] = sum(arg[s]["readCount"] for s in subtypes if s in arg and arg[s]["readCount"] >= 50)
    s16 = load(run, "silva_138.1_NR99")
    for col in COUNTS:
        for rule, v in denominators(s16, col).items():
            den.setdefault((col, rule), {})[sid] = v

print("\nARG numerator (readCount, 28 subtypes, >=50): " + ", ".join(f"{k}={v:.0f}" for k, v in num.items()))
print("Fig. 1c target ratios (approximate):          " + ", ".join(f"{sid}={t}" for _, sid, t in samples))

# ---- score: mean absolute log2 error against the figure ---------------------
scores = []
for (col, rule), d in den.items():
    errs, ok = [], True
    for _, sid, target in samples:
        v = d.get(sid, 0)
        if v <= 0 or num[sid] <= 0:
            ok = False
            break
        errs.append(abs(math.log2((num[sid] / v) / target)))
    if ok:
        scores.append((sum(errs) / len(errs), max(errs), col, rule))
scores.sort()

print("\nBEST 15 DEFINITIONS (error = mean |log2(our ratio / figure ratio)|; 0 = perfect, 0.25 = 19% off, 1 = factor 2)")
print(f"{'mean err':>9s}{'worst':>7s}  {'count column':18s} rule")
for e, w, col, rule in scores[:15]:
    print(f"{e:9.2f}{w:7.2f}  {col:18s} {rule}")

if scores:
    e, w, col, rule = scores[0]
    print(f"\nBest rule in detail: {col}, {rule}")
    print(f"{'sample':8s}{'ARG':>8s}{'16S':>10s}{'ratio':>8s}{'figure':>8s}")
    for _, sid, target in samples:
        v = den[(col, rule)][sid]
        print(f"{sid:8s}{num[sid]:8.0f}{v:10.0f}{num[sid] / v:8.2f}{target:8.2f}")
print("\nFor reference, the simplest rule (all templates, readCount):")
for _, sid, target in samples:
    v = den[("readCount", "all")][sid]
    print(f"  {sid}: 16S={v:,.0f} ratio={num[sid] / v:.4f} (figure {target})")
print("\nNOTE: development use only. A rule is plausible only if it fits all sample types, and still needs a reason in the paper's text.")
