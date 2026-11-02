#!/usr/bin/env python3
"""
16_threshold_vs_qpcr.py [--out results/threshold_vs_qpcr.tsv]

EXPLORATORY / descriptive (RQ3, detection level). Not part of the R1 verdict.
How does the metagenomic presence/absence call depend on the detection rule, and how well does it agree with the qPCR
presence/absence calls of the paper (Supplementary Table S5, saved as refdata/paper_S5_presence_absence.csv)?
qPCR is the comparison, not "truth": agreement is reported, with positives among qPCR-positives (sensitivity-like) and
positives among qPCR-negatives (extra calls).

MGS call per sample and gene (gene family = base names of the published S6 subtypes, as in scripts 09/11/15):
  rule A(T):  any template of the gene with KMA readCount >= T                 T = 1, 5, 10, 20, 50 (registered), 100, 200
  rule B(T):  as A(T), and the template covers >= 90 % of its length (Template_Coverage from the ResFinder .res)
Mapping of S5 rows to our samples: type, day and replicate number = order of the run accessions within type and day
(draft sample sheet). Check printed: our A(50) calls against the MGS columns of S5 (must agree if the mapping is right).
Standard library only. Run from the project root.
"""
import csv, sys
from pathlib import Path

KMA = Path("results/kma")
GENES = ["ermB", "sul1", "tetA", "tetQ", "tetW"]
THR = [1, 5, 10, 20, 50, 100, 200]
OUT = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else Path("results/threshold_vs_qpcr.tsv")


def mapstat(path):
    return {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(path) if not l.startswith("#")}


def cov(path):
    out = {}
    with open(path) as f:
        h = [x.strip() for x in f.readline().lstrip("#").rstrip("\n").split("\t")]
        ic = h.index("Template_Coverage")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > ic:
                out[p[0].strip()] = float(p[ic])
    return out


def main():
    rd = list(csv.reader(open("refdata/paper_S6_subtype_read_counts.csv")))
    fam = {g: set(r[1].split("_")[0] for r in rd[1:] if r[0] == g) for g in GENES}
    s5 = {}
    for r in csv.DictReader(open("refdata/paper_S5_presence_absence.csv")):
        key = (r["type"], int(r["day"]), int(r["replicate"]))
        s5[key] = ([int(r["mgs_" + g]) for g in GENES], [int(r["qpcr_" + g]) for g in GENES])
    sheet = list(csv.DictReader(open("phaseA_audit/out/sample_sheet_draft.tsv"), delimiter="\t"))
    samples = []
    for s in sheet:
        run = s["run_accession"]
        m = KMA / f"{run}_resfinder_2.1.1.mapstat"
        r = KMA / f"{run}_resfinder_2.1.1.res"
        if not (m.exists() and r.exists()):
            continue
        key = (s["type"], int(s["day"]), int(s["replicate_index_by_id_order"]))
        if key not in s5:
            print("WARNING: no S5 row for", key, s["sample_id"]); continue
        c, cv = mapstat(m), cov(r)
        calls = {}
        for rule in ("A", "B"):
            for T in THR:
                calls[(rule, T)] = [int(any(v >= T and n.split("_")[0] in fam[g] and (rule == "A" or cv.get(n, 0.0) >= 90.0)
                                            for n, v in c.items())) for g in GENES]
        samples.append((s["sample_id"], key, calls, s5[key][0], s5[key][1]))
    if len(samples) != 33:
        print(f"WARNING: {len(samples)} samples matched (expected 33)")

    # mapping check
    chk = sum(1 for _, _, calls, mg, _ in samples for a, b in zip(calls[("A", 50)], mg) if a == b)
    print(f"\nMapping check: our A(50) calls equal the MGS columns of Table S5 in {chk} of {5 * len(samples)} sample-gene pairs")
    bad = [(sid, key, g) for sid, key, calls, mg, _ in samples for g, a, b in zip(GENES, calls[("A", 50)], mg) if a != b]
    if bad:
        print("  differences:", bad)

    def stats(sel, rule, T):
        tp = fp = fn = tn = 0
        for sid, key, calls, mg, qp in samples:
            if not sel(key):
                continue
            for a, q in zip(calls[(rule, T)], qp):
                tp += a and q; fp += a and not q; fn += (not a) and q; tn += (not a) and (not q)
        return tp, fp, fn, tn

    rows = []
    for label, sel in (("all 33 samples", lambda k: True), ("POND + SED only (15 samples)", lambda k: k[0] in ("POND", "SED"))):
        print(f"\nAgreement of the metagenomic call with the qPCR call, {label}")
        print(f"{'rule':8s}{'T':>5s}{'qPCR+ & MGS+':>14s}{'qPCR+ & MGS-':>14s}{'qPCR- & MGS+':>14s}{'qPCR- & MGS-':>14s}{'MGS+ among qPCR+':>18s}{'agree':>8s}")
        for rule in ("A", "B"):
            for T in THR:
                tp, fp, fn, tn = stats(sel, rule, T)
                n = tp + fp + fn + tn
                sens = tp / (tp + fn) if tp + fn else float("nan")
                print(f"{rule:8s}{T:>5d}{tp:>14d}{fn:>14d}{fp:>14d}{tn:>14d}{sens:>18.3f}{(tp + tn) / n:>8.3f}")
                rows.append([label, rule, T, tp, fn, fp, tn, round(sens, 4), round((tp + tn) / n, 4)])
    print("\nWhere the two methods disagree at the registered rule A(50): sample, gene, qPCR call, MGS call")
    for sid, key, calls, mg, qp in samples:
        for g, a, q in zip(GENES, calls[("A", 50)], qp):
            if a != q:
                print(f"  {sid} {key[0]} day {key[1]} rep {key[2]}  {g:5s} qPCR {q}  MGS {a}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["subset", "rule", "T", "qpcrpos_mgspos", "qpcrpos_mgsneg", "qpcrneg_mgspos", "qpcrneg_mgsneg", "mgs_pos_among_qpcr_pos", "agreement"])
        w.writerows(rows)
    print(f"\nTable written to {OUT}")


if __name__ == "__main__":
    main()
