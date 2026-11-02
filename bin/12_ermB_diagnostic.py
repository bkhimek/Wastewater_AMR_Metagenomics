#!/usr/bin/env python3
"""
12_ermB_diagnostic.py [--out results/ermB_diagnostic.tsv]

EXPLORATORY (post-freeze, not part of the R1 verdict). Step A of the independent review of 2026-11-02:
for the nine effluent samples (E001 E006 E011 E016 E021 E026 E031 E036 E041), list every erm(B) template with any reads,
and, in a separate block, the neighbouring erm templates that are not erm(B) (never added to ermB totals).
Columns: sample, published ermB total (Table S6), our ermB total (readCount >= 50), template, readCount, fragmentCount,
readCountAln, fragmentCountAln, covered positions (mapstat), Score, Template_length, Template_Identity, Template_Coverage,
Query_Identity, Depth (.res). A template missing from the .res but present in the mapstat is printed with '-' (absent row).
Standard library only; reads results/kma/*.mapstat and *.res, refdata S6 CSV and the sample sheet. Run from the project root.
"""
import csv, sys
from pathlib import Path

KMA = Path("results/kma")
SAMPLES = ["E001", "E006", "E011", "E016", "E021", "E026", "E031", "E036", "E041"]
OUT = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else Path("results/ermB_diagnostic.tsv")
NEIGH = ("erm(F)", "erm(G)", "erm(Q)", "erm(35)", "erm(38)", "erm(50)", "erm(T)", "erm(C)", "erm(X)")


def mapstat(path):
    out, header = {}, None
    for line in open(path):
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = [h.strip() for h in line.lstrip("#").rstrip("\n").split("\t")]
            continue
        p = line.rstrip("\n").split("\t")
        out[p[0]] = dict(zip(header, p))
    return out


def res(path):
    out = {}
    with open(path) as f:
        header = [h.strip() for h in f.readline().lstrip("#").rstrip("\n").split("\t")]
        for line in f:
            p = [x.strip() for x in line.rstrip("\n").split("\t")]
            if len(p) >= len(header):
                out[p[0]] = dict(zip(header, p))
    return out


def f(x, d=0):
    try:
        return f"{float(x):.{d}f}"
    except (TypeError, ValueError):
        return "-"


def main():
    rd = list(csv.reader(open("refdata/paper_S6_subtype_read_counts.csv")))
    hdr = rd[0]
    sh = {r["sample_id"]: r["run_accession"] for r in csv.DictReader(open("phaseA_audit/out/sample_sheet_draft.tsv"), delimiter="\t")}
    rows = []
    for s in SAMPLES:
        run = sh[s]
        ms = mapstat(KMA / f"{run}_resfinder_2.1.1.mapstat")
        rs = res(KMA / f"{run}_resfinder_2.1.1.res")
        i = hdr.index(s)
        pub = sum(float(r[i] or 0) for r in rd[1:] if r[0] == "ermB")
        ours = sum(float(r["readCount"]) for n, r in ms.items() if n.startswith("erm(B)_") and float(r["readCount"]) >= 50)
        for block, test in (("ermB", lambda n: n.startswith("erm(B)_")), ("neighbour", lambda n: n.startswith(NEIGH) and not n.startswith("erm(B)"))):
            for n in sorted(ms, key=lambda k: -float(ms[k]["readCount"])):
                m = ms[n]
                if not test(n) or float(m["readCount"]) <= 0:
                    continue
                if block == "neighbour" and float(m["readCount"]) < 50:
                    continue
                r = rs.get(n, {})
                rows.append([s, run, block, f(pub), f(ours), n, f(m["readCount"]), f(m["fragmentCount"]), f(m.get("readCountAln")),
                             f(m.get("fragmentCountAln")), f(m.get("refCoveredPositions")), f(r.get("Score")), f(r.get("Template_length")),
                             f(r.get("Template_Identity"), 1), f(r.get("Template_Coverage"), 1), f(r.get("Query_Identity"), 1), f(r.get("Depth"), 1)])
    cols = ["sample", "run", "block", "ermB_published", "ermB_ours_ge50", "template", "readCount", "fragmentCount", "readCountAln",
            "fragmentCountAln", "coveredPos", "Score", "Template_length", "Template_Identity", "Template_Coverage", "Query_Identity", "Depth"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(cols); w.writerows(rows)
    print(f"\n{'sample':7s}{'blk':5s}{'pub':>5s}{'ours':>5s}  {'template':22s}{'rdCnt':>6s}{'frag':>6s}{'rdAln':>6s}{'frAln':>6s}{'cov.pos':>8s}{'len':>6s}{'id%':>6s}{'cov%':>6s}{'depth':>7s}")
    last = None
    for r in rows:
        if r[0] != last:
            print("-" * 112); last = r[0]
        print(f"{r[0]:7s}{r[2][:4]:5s}{r[3]:>5s}{r[4]:>5s}  {r[5]:22s}{r[6]:>6s}{r[7]:>6s}{r[8]:>6s}{r[9]:>6s}{r[10]:>8s}{r[12]:>6s}{r[13]:>6s}{r[14]:>6s}{r[16]:>7s}")
    print(f"\nTable written to {OUT}   (exploratory; E011/E021/E036 are the samples with the ermB shortfall, E016/E026/E031/E041 the matched controls)")


if __name__ == "__main__":
    main()
