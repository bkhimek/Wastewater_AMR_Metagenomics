# Wastewater AMR metagenomics: a preregistered reproduction study

**Can published wastewater resistance-gene counts be reproduced from public sequencing reads and a written description of the original analysis?**

This project reproduced the metagenomic analysis of Taylor et al. (2025) from the public reads. The pipeline met every registered evaluation criterion on 18 held-out wastewater samples. Gene-level agreement was strong. The work also showed two limits: reads cannot always be assigned to the same closely related resistance-gene subtypes as in the paper, and the paper's ARG/16S normalisation could not be reproduced.

The repository holds the analysis scripts, the preregistration, the evaluation results and the exploratory investigations. It currently documents the **R1 reproduction phase** of a wider project on how sequencing depth, normalisation and alternative analytical methods affect wastewater ARG results.

Author: Krzysztof Gizynski. Status: the registered reproduction (R1) is complete and scored; everything after the scoring is exploratory and labelled as such.

> Taylor W, Bohm K, Dyet K, Weaver L, Pattis I. 2025. Comparative analysis of qPCR and metagenomics for detecting antimicrobial resistance in wastewater: a case study. *BMC Research Notes* 18:5. [doi:10.1186/s13104-024-07027-9](https://doi.org/10.1186/s13104-024-07027-9). Reads: BioProject [PRJNA904380](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA904380).

Preregistration: [REGISTRATION_R1_reproduction_arm.md](registration/REGISTRATION_R1_reproduction_arm.md).

## Headline results

A pipeline built from the tools and settings reported by the paper, with documented assumptions where details were not stated, was run on 18 held-out liquid samples (Days 2 and 3) and scored once with a frozen script.

| Criterion (registered before the data were opened) | Evaluation result |
| --- | --- |
| P1: total ratio, ours / published, within 0.93 to 1.07 in at least 16 of 18 samples **and** in at least 5 of 6 samples of each sample type | 17 of 18 (range 0.911 to 0.991); per-type condition met |
| P2: ratio below 1.00 in at least 15 of 18 | 18 of 18 |
| P3: gene-sample pairs within 0.90 to 1.10 in at least 90% | 68 of 70 (97.1%) |
| P3b: at most 2 additional gene detections | 0 |
| P4: median subtype ratio within 0.90 to 1.10 in at least 16 of 18 | 17 of 18 |

**Why 70 and not 90?** 18 samples × 5 genes = 90 pairs, but a ratio can only be calculated where the paper reports a positive count. That is the case for 70 pairs; the other 20 have no published value. Genes that we detect but the paper does not are counted separately under P3b.

**Practical conclusion.** Gene-level agreement was stronger than subtype-level agreement. Short-read data support detecting a resistance-gene family more reliably than naming the exact allele, because counts swap between near-identical subtypes.

### What did not reproduce

- ***ermB* in effluent.** Gene ratios 0.73, 0.22 and 0.71 in three of nine effluent samples (the other four genes agree). Exploratory tests show that competing, closely related reference sequences change how KMA assigns reads, a plausible explanation; why our assignments differ from the authors' remains unresolved because their reference files were not available. Trimming did not explain the missing block in the samples tested.
- **ARG per 16S (the paper's Fig. 1c).** Our ratio is about 100 times smaller. Scaling by 100 matches influent within about 5%, but effluent and pond deviate in opposite directions. The factor of 100 is a possible explanation, not a verified unit conversion; the authors' 16S counting rule is unknown.
- **Subtype labels.** See the practical conclusion above.

### Detection against qPCR (exploratory)

Metagenomic presence/absence reproduces the paper's Table 1 in 20 of 20 cells and its Table S5 in 165 of 165 sample-gene pairs. Against the published qPCR calls, metagenomics finds 26 of 36 positive calls in pond and sediment. For the other 10, the analysis assigned no reads to the gene family, and no read threshold from 1 to 100 changes this.

Two cautions. The paper's qPCR detection limit was about 50 copies/µL, which is not equivalent to the 50-read threshold used for metagenomic calls. And no assigned reads does not prove that the gene is biologically absent.

## Design in brief

- **Development set** (Day 1, 9 liquid samples): used to set open choices such as the 16S counting rule.
- **Evaluation set** (Days 2 and 3, 18 liquid samples): downloaded but not opened until the registration was frozen (tag `R1-registered`), then run and scored once.
- **Sediment** (6 samples, Days 2 and 3): analysed descriptively only, because no sediment was available for development.
- All predictions, including failures, are reported. Changes after the freeze are listed in [R1_post_freeze_deviations.md](registration/R1_post_freeze_deviations.md).

## Where to look

| What | Link |
| --- | --- |
| Preregistration (version 2): pipeline, assumptions, criteria, counting definitions | [REGISTRATION_R1_reproduction_arm.md](registration/REGISTRATION_R1_reproduction_arm.md) |
| Changes made after the freeze | [R1_post_freeze_deviations.md](registration/R1_post_freeze_deviations.md) |
| Evaluation scoring and raw output | [R1_results_eval.md](registration/R1_results_eval.md), [score_R1_eval_output.txt](registration/score_R1_eval_output.txt) |
| Exploratory *ermB* investigation | [R1_ermB_diagnostic_note.md](registration/R1_ermB_diagnostic_note.md), [ermB_diagnostic_effluent.tsv](registration/ermB_diagnostic_effluent.tsv) |
| Exploratory tables: ARG per 16S under seven counting rules; detection per sample type | [ratio_16S_all_samples.tsv](registration/ratio_16S_all_samples.tsv), [detection_vs_table1.tsv](registration/detection_vs_table1.tsv) |
| Scripts, numbered in the order they were written | [bin/](bin/) and the inventory in [docs/WORKFLOW.md](docs/WORKFLOW.md) |
| Dataset audit: ENA run table, run classification, sample-sheet checks | [audit/phaseA/](audit/phaseA/) |
| Scientific report (Word file): methods, results, limitations | [Reproducing_ARG_metagenomics_vs_qPCR_in_wastewater_report.docx](docs/Reproducing_ARG_metagenomics_vs_qPCR_in_wastewater_report.docx) |

The main scripts are `06_trim_kma.sh` (the registered pipeline for one sample), `09_score_R1.py` (the frozen scoring script) and `11` to `16` (exploratory analyses). Exploratory scripts write to their own folders and never overwrite registered outputs.

## Tags

| Tag | State |
| --- | --- |
| `phaseA-go` | Dataset audit complete, go decision |
| `R1-registered` | Registration frozen, before any evaluation or sediment file was opened |
| `R1-eval-scored` | Evaluation scored once, post-freeze changes logged |
| `R1-exploratory-ermB` | Exploratory ARG per 16S and *ermB* diagnostic added |

## Software

| Tool | Version |
| --- | --- |
| KMA | 1.4.9 (flags `-ipe -mem_mode -ef -1t1 -cge -nf -vcf`) |
| bbduk (BBTools) | 38.90 (`maxns=0 minlength=75 qtrim=r trimq=20`, adapter options assumed) |
| ResFinder database | tag 2.1.1 (`all.fsa`, 3,146 sequences) |
| SILVA | 138.1 SSURef NR99 (RNA converted to DNA) |
| Python | 3; the scripts in `bin/` use only the standard library (`audit/phaseA/01_classify_runs.py` also needs `openpyxl`) |

## Reproducing the analysis

The numerical results were reproduced on the author's machine. Rebuilding them from a fresh clone has **not yet been tested**; a clean-environment reproduction test is a planned milestone. The steps below are what the scripts expect.

**Before you start**

- Use a Linux or WSL disk, not a Windows-mounted path.
- Clone the repository to `$HOME/projects/Wastewater_AMR_Metagenomics`. `06_trim_kma.sh` assumes this location; to use another, set `REPO=/your/path` before running it.
- The scripts read the sample sheet from `phaseA_audit/out/sample_sheet_draft.tsv`, although the audit scripts live in `audit/phaseA/`. Step 1 creates that file in the expected place.
- Raw reads need roughly 200 GB of disk at peak. Building the SILVA index needs about 37.6 GB of RAM.

**1. Software environments.** The scripts call three conda environments by name. No environment files are exported yet; a plausible installation (versions from the table above, not tested from scratch) is:

```bash
mamba create -n wwamr-kma   -c conda-forge -c bioconda kma=1.4.9
mamba create -n wwamr-bbduk -c conda-forge -c bioconda bbmap=38.90
mamba create -n wwamr-qc    -c conda-forge -c bioconda fastp seqkit fastqc
```

**2. Sample sheet.** Download the paper's Supplementary File 1 (the `...MOESM1_ESM.docx` file) and fetch the ENA run table, then build the sheet where the scripts expect it:

```bash
bash audit/phaseA/00_fetch_ena_runtable.sh      # writes the ENA run table; see the script header for the output path
mkdir -p phaseA_audit/out
python3 audit/phaseA/02_build_sample_sheet.py \
  audit/phaseA/raw/ena_filereport_PRJNA904380.tsv \
  /path/to/13104_2024_7027_MOESM1_ESM.docx \
  phaseA_audit/out
```

**3. Reference data.** Obtain Supplementary Tables S6 (MOESM3) and S5 (MOESM2) from the paper and save them as `refdata/paper_S6_subtype_read_counts.csv` and `refdata/paper_S5_presence_absence.csv` (column layouts are described in the headers of scripts `09` and `16`). They are not redistributed here.

**4. Databases.** Build both KMA indexes under `refdb/`:

```bash
# ResFinder database, tag 2.1.1 (commit 1a53e55)
mkdir -p refdb/kma_index && cd refdb
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
git -C resfinder_db checkout 2.1.1
kma index -i resfinder_db/all.fsa -o kma_index/resfinder_2.1.1     # 3,146 sequences, seconds

# SILVA 138.1 SSURef NR99
mkdir -p silva_138.1 && cd silva_138.1
curl -fLO https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
md5sum SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz   # expect 0b47532bddb5cb2765f38d7420345cbd (197,133,995 bytes)
zcat SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz | sed '/^>/!y/U/T/' > SILVA_138.1_SSURef_NR99_DNA.fasta   # 510,508 sequences
kma index -i SILVA_138.1_SSURef_NR99_DNA.fasta -o ../kma_index/silva_138.1_NR99
```

The SILVA index failed at 28 GB of RAM locally; the one used here was built on a 64 GB cloud instance (about 7 hours). The SILVA release page serves HTML for wrong addresses, so check the file size and md5, not only the download's exit code. The paper does not state the SILVA variant or the ResFinder file used; NR99 and `all.fsa` are documented assumptions.

**5. Reads.** `bash bin/03c_fetch_set.sh` is the recommended download method (chunked, checked against the ENA md5). `bin/03_download_reads.sh` is **superseded**: its resume-based method could corrupt files. Usage is in the script headers.

**6. Pipeline.** `bash bin/06_trim_kma.sh <run_accession>` for development runs. Evaluation runs need `ALLOW_EVAL=1`, for example `ALLOW_EVAL=1 bash bin/10_run_set.sh eval liquid`.

**7. Scoring.** `python3 bin/09_score_R1.py --set dev` or `--set eval` (add `--detail` for per-sample tables; `--selftest` checks the scorer).

**8. Exploratory analyses.** Run scripts 11 to 16 as described in their headers.

## Data and licence notes

- Raw reads, indexes and result folders are not committed (size); they are rebuilt from public sources.
- The article is published under CC BY-NC-ND 4.0. The licence of the individual supplementary files has not been verified. As a precaution, their content (Tables S5 and S6 as CSV, and the sample sheet, which carries sample IDs, types and days from Table S4) is not republished here.
- The code and the report are released under the [MIT licence](LICENSE), copyright Krzysztof Gizynski. The licence covers the author's own work only; the paper and its supplementary tables remain under their own terms (see above).
- The paper's Fig. 1c values used for comparison were read off the published graph (about ±0.1) and are approximate.

## Limits

The hold-out is imperfect: the published values were known beforehand, and all samples come from one plant. Influent and effluent were sampled on 19 to 21 November 2019 and pond and sediment on 10 to 12 December 2019, so sample type and sampling date are confounded. The 18 evaluation samples are not independent experiments. The qPCR comparison is at presence/absence level, because per-sample quantitative qPCR values are not public.
