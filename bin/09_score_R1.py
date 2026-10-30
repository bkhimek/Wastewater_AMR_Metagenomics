#!/usr/bin/env python3
"""
09_score_R1.py  --set dev|eval  [--detail] [--selftest]

Scoring script for registration R1 (reproduction arm). Reads only results/kma/*.mapstat, the
published Table S6 (CSV) and the sample sheet. Standard library only. Run from the project root.

REGISTRATION R1 (header copied from the registration draft; thresholds are those registered in the
frozen version of REGISTRATION_R1_reproduction_arm.md; keep both in sync).

Definitions (fixed):
  - Our subtype count = KMA readCount (column 2 of the .mapstat) of the ResFinder 2.1.1 template whose
    name equals the published subtype name; 0 if the template is absent from the .mapstat.
  - A subtype is "called" when our count >= 50 (inclusive; Table S6 contains an exact 50: R006,
    tet(W)_1_DQ060146). A published count of 0 means "not reported" (below the authors' filter).
  - Published-positive = published count > 0.
  - P1 ratio (per sample, and per gene) = sum of our counts over CALLED subtypes / sum of published counts.
  - P2 = P1 ratio < 1.00.
  - P3 (positive-call recovery) per sample: TP = published-positive and called; FN = published-positive and
    not called; FP = published 0 and called ("extra positive"). recall = TP / (TP + FN).
  - P4 per sample: over published-positive subtypes, ratio = our raw count / published count (our count 0
    gives ratio 0, so misses lower the median); report the median and the share of subtypes within +-10%.
  - P5 (secondary, outside the verdict): excluded-template share of the SILVA readCount =
    (Eukaryota + Chloroplast + Mitochondria) / all templates; ordering POND > EFF > INF per day (mean of
    the replicates).
  - Primary 16S denominator (D13) and the sensitivity rules are computed by a separate script and are not
    part of the verdict.
"""
import csv, statistics as st, sys, tempfile, os
from pathlib import Path

KMA = Path("results/kma")
PAPER = Path("refdata/paper_S6_subtype_read_counts.csv")
SHEET = Path("phaseA_audit/out/sample_sheet_draft.tsv")
MIN_CALL = 50
GENES = ["ermB", "sul1", "tetA", "tetQ", "tetW"]
TYPES = ["INF", "EFF", "POND"]

# ---- registered thresholds (to be fixed at the freeze; None = not yet registered) ---------------
REG = dict(
    P1_low=0.93, P1_high=1.07, P1_min_overall=16, P1_min_per_type=5,   # of 18 samples, of 6 per type
    P2_min=15,
    P3_gene_band=(0.90, 1.10), P3_gene_min_share=0.90,                # P3: gene-sample pairs (published > 0) in band
    P3b_false_gene_max=2,                                             # P3b: genes called where published 0
    P3c_recall_min=0.90, P3c_extra_max=18,                            # P3c: subtype level, reported, outside verdict
    P4_low=0.90, P4_high=1.10, P4_min=16,
)


def called(count):
    return count >= MIN_CALL


def load_paper(path):
    with open(path) as f:
        rd = list(csv.reader(f))
    hdr, rows = rd[0], rd[1:]
    sub = {}
    for r in rows:
        sub[r[1]] = (r[0], {s: int(float(r[hdr.index(s)] or 0)) for s in hdr[2:35]})
    return sub


def load_sheet(path):
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        return list(rd)


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
            row = dict(zip(header, p))
            out[row["refSequence"]] = float(row["readCount"])
    return out


def score_sample(paper, sid, ours):
    """ours: dict subtype -> readCount. Returns metrics dict."""
    m = {}
    called_sum, paper_sum = 0.0, 0.0
    per_gene = {g: [0.0, 0.0] for g in GENES}
    tp = fn = fp = 0
    ratios = []
    for name, (gene, counts) in paper.items():
        pc, oc = counts[sid], ours.get(name, 0.0)
        paper_sum += pc
        per_gene[gene][1] += pc
        if called(oc):
            called_sum += oc
            per_gene[gene][0] += oc
        if pc > 0:
            ratios.append(oc / pc)
            if called(oc): tp += 1
            else: fn += 1
        elif called(oc):
            fp += 1
    m["ratio"] = called_sum / paper_sum if paper_sum else float("nan")
    m["gene_ratio"] = {g: (a / b if b else float("nan")) for g, (a, b) in per_gene.items()}
    m["tp"], m["fn"], m["fp"] = tp, fn, fp
    m["recall"] = tp / (tp + fn) if (tp + fn) else float("nan")
    m["median_ratio"] = st.median(ratios) if ratios else float("nan")
    m["within10"] = sum(1 for r in ratios if 0.9 <= r <= 1.1) / len(ratios) if ratios else float("nan")
    m["n_pos"] = tp + fn
    m["paper_sum"], m["our_sum"] = paper_sum, called_sum
    # gene-family numerator (G1): all templates whose base name belongs to the gene (derived from the published
    # subtype names, e.g. ermB -> erm(B)); every template with count >= MIN_CALL counts, published or not.
    bases = {g: set() for g in GENES}
    for name, (gene, _) in paper.items():
        bases[gene].add(name.split("_")[0])
    fam = {g: 0.0 for g in GENES}
    for name, c in ours.items():
        for g in GENES:
            if name.split("_")[0] in bases[g] and called(c):
                fam[g] += c
    m["fam"] = fam
    m["fam_gene_ratio"] = {g: (fam[g] / per_gene[g][1] if per_gene[g][1] else float("nan")) for g in GENES}
    m["fam_ratio"] = sum(fam.values()) / paper_sum if paper_sum else float("nan")
    return m


def excl_share(path):
    """Share of the SILVA readCount (mapstat, column 2) on Eukaryota, Chloroplast or Mitochondria templates."""
    tot = ex = 0.0
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            n = p[0].split(" ", 1)[1] if " " in p[0] else p[0]
            c = float(p[1])
            tot += c
            if n.startswith("Eukaryota") or "Chloroplast" in n or "Mitochondria" in n:
                ex += c
    return ex / tot if tot else float("nan")


def selftest():
    assert not called(49) and called(50) and called(51), "boundary rule"
    paper = {"gA_1": ("ermB", {"S1": 100, "S2": 0}), "gB_1": ("sul1", {"S1": 50, "S2": 0}),
             "gC_1": ("tetA", {"S1": 200, "S2": 60})}
    m = score_sample(paper, "S1", {"gA_1": 100, "gB_1": 49, "gC_1": 220})
    assert (m["tp"], m["fn"], m["fp"]) == (2, 1, 0), m            # gB_1 at 49 is a miss
    assert abs(m["ratio"] - 320 / 350) < 1e-9, m
    assert abs(m["median_ratio"] - st.median([1.0, 0.98, 1.1])) < 1e-9, m
    m = score_sample(paper, "S2", {"gA_1": 51, "gC_1": 50})
    assert (m["tp"], m["fn"], m["fp"]) == (1, 0, 1), m             # gC_1 at exactly 50 is called (TP); gA_1 at 51 vs published 0 is an extra positive
    assert abs(m["ratio"] - (51 + 50) / 60) < 1e-9, m
    m = score_sample(paper, "S1", {})
    assert m["ratio"] == 0 and m["recall"] == 0 and m["fp"] == 0, m
    # a pipeline that detects nothing: agreement is high but recall is 0
    print("selftest passed (boundaries 49/50/51, extra positives, total miss)")


def main():
    if "--selftest" in sys.argv:
        selftest(); return
    which = sys.argv[sys.argv.index("--set") + 1]
    expect = {"dev": 9, "eval": 18}[which]
    paper = load_paper(PAPER)
    sheet = [r for r in load_sheet(SHEET) if r["split"] == which and r["matrix"] == "liquid"]
    ids = [r["sample_id"] for r in sheet]
    if len(sheet) != expect or len(set(ids)) != expect:
        sys.exit(f"expected {expect} unique {which} liquid samples, found {len(sheet)} ({len(set(ids))} unique): refusing to score")
    for r in sheet:
        for need in (KMA / f"{r['run_accession']}_resfinder_2.1.1.mapstat", KMA / f"{r['run_accession']}_silva_138.1_NR99.mapstat"):
            if not need.exists() or need.stat().st_size == 0:
                sys.exit(f"missing or empty result file: {need}")
    res = []
    for r in sorted(sheet, key=lambda x: (TYPES.index(x["type"]), int(x["day"]), x["sample_id"])):
        ours = read_mapstat(KMA / f"{r['run_accession']}_resfinder_2.1.1.mapstat")
        m = score_sample(paper, r["sample_id"], ours)
        m["excl"] = excl_share(KMA / f"{r['run_accession']}_silva_138.1_NR99.mapstat")
        m.update(sid=r["sample_id"], type=r["type"], day=int(r["day"]))
        m["ours"] = ours
        res.append(m)

    print(f"\nR1 scoring, set = {which}, n = {len(res)}; call threshold: count >= {MIN_CALL}\n")
    print(f"{'sample':7s}{'type':6s}{'day':>4s}{'ratio':>8s}{'paper':>8s}{'ours':>8s}{'posPub':>8s}{'TP':>4s}{'FN':>4s}{'FP':>4s}{'recall':>8s}{'medRat':>8s}{'in+-10%':>9s}{'excl16S':>9s}")
    for m in res:
        print(f"{m['sid']:7s}{m['type']:6s}{m['day']:4d}{m['ratio']:8.3f}{m['paper_sum']:8.0f}{m['our_sum']:8.0f}{m['n_pos']:8d}{m['tp']:4d}{m['fn']:4d}{m['fp']:4d}{m['recall']:8.2f}{m['median_ratio']:8.2f}{m['within10']:9.2f}{m['excl']:9.2f}")

    print("\nPer-gene ratios (ours called / published):")
    print(f"{'sample':7s}" + "".join(f"{g:>8s}" for g in GENES))
    for m in res:
        print(f"{m['sid']:7s}" + "".join(f"{m['gene_ratio'][g]:8.2f}" for g in GENES))

    print("\nGene-family numerator G1 (all templates of the gene with count >= 50) / published sum of the gene:")
    print(f"{'sample':7s}{'total':>8s}" + "".join(f"{g:>8s}" for g in GENES))
    for m in res:
        print(f"{m['sid']:7s}{m['fam_ratio']:8.3f}" + "".join(f"{m['fam_gene_ratio'][g]:8.2f}" for g in GENES))
    pairs = [(m, g) for m in res for g in GENES if m["fam_gene_ratio"][g] == m["fam_gene_ratio"][g]]   # published > 0
    inb = [(m["sid"], g, m["fam_gene_ratio"][g]) for m, g in pairs if 0.90 <= m["fam_gene_ratio"][g] <= 1.10]
    out = [(m["sid"], g, m["fam_gene_ratio"][g]) for m, g in pairs if not (0.90 <= m["fam_gene_ratio"][g] <= 1.10)]
    gfp = [(m["sid"], g, m["fam"][g]) for m in res for g in GENES
           if m["fam_gene_ratio"][g] != m["fam_gene_ratio"][g] and m["fam"][g] > 0]
    print(f"  G1: gene-sample pairs with a published value: {len(pairs)}; within 0.90-1.10: {len(inb)}; outside: {out}")
    print(f"  G2: gene-level false positives (published 0, ours called): {len(gfp)} {gfp}")
    share = len(inb) / len(pairs) if pairs else float("nan")
    print(f"  P3  gene-level pairs in band: {len(inb)} of {len(pairs)} = {share:.3f} (registered: >= {REG['P3_gene_min_share']})")
    print(f"  P3b false genes: {len(gfp)} (registered: <= {REG['P3b_false_gene_max']})")
    tpt, fnt, fpt = sum(m["tp"] for m in res), sum(m["fn"] for m in res), sum(m["fp"] for m in res)
    print(f"  P3c subtype level (reported, outside verdict): pooled recall {tpt / (tpt + fnt):.3f} (registered: >= {REG['P3c_recall_min']}); extra positives {fpt} (registered: <= {REG['P3c_extra_max']})")

    if "--detail" in sys.argv:
        print("\nMisses (published-positive, our count < 50) and extra positives (published 0, our count >= 50):")
        print(f"{'sample':7s}{'kind':8s}{'subtype':24s}{'published':>10s}{'ours':>8s}")
        for m in res:
            for name, (gene, counts) in paper.items():
                pc, oc = counts[m["sid"]], m["ours"].get(name, 0.0)
                if pc > 0 and not called(oc):
                    print(f"{m['sid']:7s}{'MISS':8s}{name:24s}{pc:10d}{oc:8.0f}")
                elif pc == 0 and called(oc):
                    print(f"{m['sid']:7s}{'EXTRA':8s}{name:24s}{pc:10d}{oc:8.0f}")

    R = REG
    p1 = [R["P1_low"] <= m["ratio"] <= R["P1_high"] for m in res]
    p2 = [m["ratio"] < 1.0 for m in res]
    p4 = [R["P4_low"] <= m["median_ratio"] <= R["P4_high"] for m in res]
    print("\nCounts of samples meeting each registered criterion (verdict only for set = eval, n = 18):")
    print(f"  P1 total ratio in [{R['P1_low']}, {R['P1_high']}]: {sum(p1)} of {len(res)} (registered: >= {R['P1_min_overall']} of 18)")
    for t in TYPES:
        sel = [ok for ok, m in zip(p1, res) if m["type"] == t]
        print(f"      {t}: {sum(sel)} of {len(sel)} (registered: >= {R['P1_min_per_type']} of 6 per type)")
    print(f"  P2 ratio < 1.00:                 {sum(p2)} of {len(res)} (registered: >= {R['P2_min']} of 18)")
    print("  P3c subtype-level detail per sample type (reported, outside the verdict):")
    for t in TYPES:
        sel = [m for m in res if m["type"] == t]
        print(f"      {t}: recall {min(m['recall'] for m in sel):.2f} to {max(m['recall'] for m in sel):.2f}; "
              f"misses (FN) {sum(m['fn'] for m in sel)}; extra positives (FP) {sum(m['fp'] for m in sel)}; published positives {sum(m['n_pos'] for m in sel)}")
    print(f"  P4 median ours/published in [{R['P4_low']}, {R['P4_high']}]: {sum(p4)} of {len(res)} (registered: >= {R['P4_min']} of 18)")
    print("  P5 (secondary, outside verdict) excluded-template share of 16S readCount, mean per type and day:")
    for d in sorted({m["day"] for m in res}):
        parts = []
        for t in TYPES:
            v = [m["excl"] for m in res if m["type"] == t and m["day"] == d]
            parts.append((t, sum(v) / len(v) if v else float("nan")))
        order = " > ".join(t for t, _ in sorted(parts, key=lambda x: -x[1]))
        print(f"      day {d}: " + ", ".join(f"{t} {v:.2f}" for t, v in parts) + f"  -> {order}")
    print("\nNOTE: set = dev is for testing the script and setting thresholds. The eval run happens once, after the freeze.")


if __name__ == "__main__":
    main()
