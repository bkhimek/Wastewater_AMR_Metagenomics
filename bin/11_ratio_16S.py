#!/usr/bin/env python3
"""
11_ratio_16S.py [--out results/ratio_16S.tsv]

ARG / 16S ratio for every sample of the sample sheet that has KMA results (dev + eval, liquid + sediment).
Descriptive analysis; NOT part of the R1 verdict. Standard library only. Run from the project root.

Numerator (registration section 4): per gene (ermB, sul1, tetA, tetQ, tetW) the sum of KMA counts over ALL
ResFinder 2.1.1 templates of that gene family whose count >= 50 (family = base names of the published subtypes,
e.g. ermB -> erm(B)); total = sum over the 5 genes. The >= 50 rule is always applied to readCount; for the
alternative count columns the same templates are summed.

Denominator, PRIMARY (decision D13): sum of readCount over SILVA 138.1 NR99 templates labelled Bacteria or Archaea,
excluding Chloroplast, Mitochondria (and Eukaryota, which are not Bacteria/Archaea).

Sensitivity rules (all reported, fixed before the evaluation):
  D13      primary
  all      all SILVA templates
  bact     Bacteria only (no Archaea), organelles excluded
  org      Bacteria + Archaea INCLUDING Chloroplast and Mitochondria templates
  frag     D13 taxa, fragmentCount instead of readCount (numerator also fragmentCount)
  aln      D13 taxa, readCountAln instead of readCount (numerator also readCountAln)
  cov90    D13 taxa, only templates with Template_Coverage >= 90 % (from the SILVA .res)

Column "fig1c" (Day-1 development samples only): approximate readings by eye of the paper's Fig. 1c (+-0.1).
No rule is chosen by agreement with Fig. 1c.
"""
import csv, sys, statistics as st
from pathlib import Path

KMA = Path("results/kma")
PAPER = Path("refdata/paper_S6_subtype_read_counts.csv")
SHEET = Path("phaseA_audit/out/sample_sheet_draft.tsv")
GENES = ["ermB", "sul1", "tetA", "tetQ", "tetW"]
MIN_CALL = 50
FIG1C = {"R001": 2.25, "R006": 2.45, "R011": 2.30, "E001": 1.95, "E006": 2.05, "E011": 2.10,
         "P001": 0.40, "P006": 0.55, "P011": 0.55}
OUT = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else Path("results/ratio_16S.tsv")


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


def read_cov(path):
    cov = {}
    with open(path) as f:
        header = [h.strip() for h in f.readline().lstrip("#").rstrip("\n").split("\t")]
        ic, it = header.index("Template_Coverage"), header.index("Template")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > ic:
                cov[p[it].strip()] = float(p[ic])
    return cov


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def taxon(name):
    """SILVA header: 'ACC.start.stop Bacteria;Phylum;...' -> lineage string"""
    return name.split(" ", 1)[1] if " " in name else name


def is_org(lin):
    return "Chloroplast" in lin or "Mitochondria" in lin


def d13(lin):
    return (lin.startswith("Bacteria") or lin.startswith("Archaea")) and not is_org(lin)


def families():
    with open(PAPER) as f:
        rd = list(csv.reader(f))
    fam = {g: set() for g in GENES}
    for r in rd[1:]:
        fam[r[0]].add(r[1].split("_")[0])
    return fam


def numerator(arg, fam, col):
    tot = {g: 0.0 for g in GENES}
    for name, row in arg.items():
        if num(row["readCount"]) < MIN_CALL:
            continue
        base = name.split("_")[0]
        for g in GENES:
            if base in fam[g]:
                tot[g] += num(row[col])
    return tot


def main():
    fam = families()
    sheet = list(csv.DictReader(open(SHEET), delimiter="\t"))
    rows, skipped = [], []
    for s in sheet:
        run = s["run_accession"]
        paths = [KMA / f"{run}_{d}.{e}" for d in ("resfinder_2.1.1", "silva_138.1_NR99") for e in ("res", "mapstat")]
        if not all(p.exists() and p.stat().st_size > 0 for p in paths):
            skipped.append(s["sample_id"]); continue
        arg = read_mapstat(KMA / f"{run}_resfinder_2.1.1.mapstat")
        sil = read_mapstat(KMA / f"{run}_silva_138.1_NR99.mapstat")
        cov = read_cov(KMA / f"{run}_silva_138.1_NR99.res")
        n_rc = numerator(arg, fam, "readCount"); n_fc = numerator(arg, fam, "fragmentCount"); n_al = numerator(arg, fam, "readCountAln")

        def den(col, keep, extra=None):
            return sum(num(r[col]) for n, r in sil.items() if keep(taxon(n)) and (extra is None or extra(n)))
        D = {
            "D13": (sum(n_rc.values()), den("readCount", d13)),
            "all": (sum(n_rc.values()), den("readCount", lambda l: True)),
            "bact": (sum(n_rc.values()), den("readCount", lambda l: l.startswith("Bacteria") and not is_org(l))),
            "org": (sum(n_rc.values()), den("readCount", lambda l: l.startswith("Bacteria") or l.startswith("Archaea"))),
            "frag": (sum(n_fc.values()), den("fragmentCount", d13)),
            "aln": (sum(n_al.values()), den("readCountAln", d13)),
            "cov90": (sum(n_rc.values()), den("readCount", d13, lambda n: cov.get(n, 0.0) >= 90.0)),
        }
        row = dict(sample=s["sample_id"], run=run, matrix=s["matrix"], type=s["type"], day=s["day"], split=s["split"],
                   arg_total=sum(n_rc.values()), d13_den=D["D13"][1])
        for g in GENES:
            row["ratio_" + g] = n_rc[g] / D["D13"][1] if D["D13"][1] else float("nan")
        for k, (a, b) in D.items():
            row[k] = a / b if b else float("nan")
        row["fig1c"] = FIG1C.get(s["sample_id"], "")
        rows.append(row)

    if not rows:
        sys.exit("no samples with complete KMA results")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

    RULES = ["D13", "all", "bact", "org", "frag", "aln", "cov90"]
    print(f"\nARG / 16S under the primary rule (D13) and sensitivity rules; n = {len(rows)} samples (skipped: {skipped or 'none'})")
    print(f"{'sample':7s}{'matrix':9s}{'type':5s}{'day':>4s}{'split':>6s}{'ARG':>8s}{'16S(D13)':>10s}" + "".join(f"{r:>8s}" for r in RULES) + f"{'fig1c':>7s}")
    key = lambda r: (r["matrix"], ["INF", "EFF", "POND", "SED"].index(r["type"]) if r["type"] in ("INF", "EFF", "POND", "SED") else 9, int(r["day"]), r["sample"])
    for r in sorted(rows, key=key):
        print(f"{r['sample']:7s}{r['matrix']:9s}{r['type']:5s}{int(r['day']):4d}{r['split']:>6s}{r['arg_total']:8.0f}{r['d13_den']:10.0f}" +
              "".join(f"{r[k]:8.3f}" for k in RULES) + f"{(str(r['fig1c']) if r['fig1c'] != '' else '-'):>7s}")
    print("\nMean ratio by group and rule (last column: range of the primary rule D13):")
    groups = sorted({(r["matrix"], r["type"], r["day"]) for r in rows})
    print(f"{'matrix':9s}{'type':5s}{'day':>4s}{'n':>3s}" + "".join(f"{k:>8s}" for k in RULES) + f"{'D13 range':>14s}")
    for g in groups:
        sel = [r for r in rows if (r["matrix"], r["type"], r["day"]) == g]
        print(f"{g[0]:9s}{g[1]:5s}{int(g[2]):4d}{len(sel):3d}" + "".join(f"{st.mean(r[k] for r in sel):8.3f}" for k in RULES) +
              f"  {min(r['D13'] for r in sel):.3f}-{max(r['D13'] for r in sel):.3f}")
    dev = [r for r in rows if r["fig1c"] != ""]
    if dev:
        print("\nDay-1 development samples: primary ratio vs Fig. 1c (approximate reading, +-0.1); descriptive, no rule is chosen by this:")
        for r in dev:
            print(f"  {r['sample']}: ours {r['D13']:.4f}   figure {r['fig1c']}   ratio ours/figure {r['D13'] / float(r['fig1c']):.4f}")
    print(f"\nTable written to {OUT}")


if __name__ == "__main__":
    main()
