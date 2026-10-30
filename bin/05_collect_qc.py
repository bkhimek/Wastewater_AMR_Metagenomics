"""Collect per-sample QC numbers into results/qc_summary_pilot.tsv (counts only; no reads are opened)."""
import csv, json, os, re, sys
repo, runs = sys.argv[1], sys.argv[2:]
res = os.path.join(repo, "results")
cols = ["run","pairs","avg_len","gb_pair","q30_r1","q30_r2","gc","fastqc_fail_modules",
        "bbduk_reads_in","bbduk_qtrimmed_pct","bbduk_ktrimmed_pct","bbduk_removed_pct","bbduk_remaining_pct",
        "fastp_pass_pct","fastp_lowq_pct","fastp_short_pct","fastp_adapter_pct"]
def seqkit(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    a, b = rows[0], rows[1]
    return {"pairs": a["num_seqs"], "avg_len": a["avg_len"],
            "gb_pair": round((int(a["sum_len"]) + int(b["sum_len"])) / 1e9, 2),
            "q30_r1": a["Q30(%)"], "q30_r2": b["Q30(%)"], "gc": a["GC(%)"]}
def bbduk(path):
    t = open(path).read()
    g = lambda lab: re.search(lab + r":\s+(\d+) reads \(([\d.]+)%\)", t)
    i = re.search(r"Input:\s+(\d+) reads", t)
    if not (i and g("QTrimmed") and g("Total Removed") and g("Result")): return {}
    kt = g("KTrimmed")
    return {"bbduk_reads_in": i.group(1), "bbduk_qtrimmed_pct": g("QTrimmed").group(2),
            "bbduk_ktrimmed_pct": kt.group(2) if kt else "0",
            "bbduk_removed_pct": g("Total Removed").group(2), "bbduk_remaining_pct": g("Result").group(2)}
def fastp(path):
    j = json.load(open(path)); n = j["summary"]["before_filtering"]["total_reads"]
    f = j["filtering_result"]; p = lambda x: round(100.0 * x / n, 2)
    return {"fastp_pass_pct": p(f["passed_filter_reads"]), "fastp_lowq_pct": p(f["low_quality_reads"]),
            "fastp_short_pct": p(f["too_short_reads"]), "fastp_adapter_pct": p(j["adapter_cutting"]["adapter_trimmed_reads"])}
def fastqc_fails(run):
    import zipfile
    out = []
    for m in (1, 2):
        z = os.path.join(res, "qc_raw_dev", f"{run}_{m}_fastqc.zip")
        if not os.path.exists(z): return "NA"
        zf = zipfile.ZipFile(z); name = [x for x in zf.namelist() if x.endswith("summary.txt")][0]
        out += [f"R{m}:{l.split(chr(9))[1]}" for l in zf.read(name).decode().splitlines() if l.startswith("FAIL")]
    return ";".join(out) if out else "none"
dest = os.path.join(res, "qc_summary_pilot.tsv")
with open(dest, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="NA"); w.writeheader()
    for run in runs:
        row = {"run": run}; sk = os.path.join(res, "qc_stats", f"seqkit_{run}.tsv")
        if not os.path.exists(sk): continue
        row.update(seqkit(sk)); row["fastqc_fail_modules"] = fastqc_fails(run)
        b = os.path.join(res, "qc_trim_dev", f"bbduk_dryrun_{run}.log")
        if os.path.exists(b): row.update(bbduk(b))
        f = os.path.join(res, "qc_trim_dev", f"fastp_{run}.json")
        if os.path.exists(f): row.update(fastp(f))
        w.writerow(row)
print("wrote", dest)
