# Workflow and script inventory

Scripts in `bin/` are numbered in the order they were written. The registered verdict uses only the scripts marked "yes".
See the [README](../README.md) for the overview and the reproduction steps.

| Script | Purpose | Part of the registered verdict |
| --- | --- | --- |
| `03c_fetch_set.sh` (with `03b_fetch_chunked.sh`, `chunk_fetch.sh`) | Recommended download: chunked, checksum-verified, from ENA | transfer only |
| `03_download_reads.sh` | Superseded older download script (resume-based; could corrupt files) | transfer only |
| `04_network_diag.sh` | Test whether a slow download is ENA or the local link | no |
| `05_qc_pilot.sh`, `05_collect_qc.py` | Per-sample QC on development samples | no |
| `06_trim_kma.sh` | **Registered pipeline** for one sample: bbduk trimming, then KMA against ResFinder and SILVA | yes |
| `10_run_set.sh` | Runs `06_trim_kma.sh` once per run of a set | yes |
| `09_score_R1.py` | **Frozen scoring script** (P1 to P4, P3b, P3c, P5); `--set dev` or `--set eval`, `--selftest` | yes |
| `07_compare_paper_counts.py`, `08_scan_16S_definitions.py` | Development-only comparison with the paper; scan of 16S counting rules | no |
| `11_ratio_16S.py` | ARG per 16S for all samples under the primary rule and six alternatives | exploratory |
| `12_ermB_diagnostic.py`, `13_explore_raw.sh`, `14_explore_single_template.sh` | *ermB* investigation: per-template table, raw reads, single-template index | exploratory |
| `15_detection_vs_table1.py`, `16_threshold_vs_qpcr.py` | Detection per sample type; threshold sensitivity against published qPCR calls | exploratory |

Audit scripts in `audit/phaseA/`: `00_fetch_ena_runtable.sh` (ENA run table), `01_classify_runs.py` (needs `openpyxl`), `02_build_sample_sheet.py <ena_table.tsv> <MOESM1.docx> <outdir>` (sample sheet).

Known inconsistency: downstream scripts read `phaseA_audit/out/sample_sheet_draft.tsv`, while the audit scripts are stored in `audit/phaseA/`. Standardising the paths is on the to-do list.
