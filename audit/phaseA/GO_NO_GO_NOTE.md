# Phase A go/no-go note (2026-10-27)

Dataset: Taylor et al. (2025), BMC Research Notes 18:5, doi 10.1186/s13104-024-07027-9; BioProject PRJNA904380. Research use only.
Evidence: `raw/ena_filereport_PRJNA904380.tsv` (fetched 2026-10-27T15:02:09Z, sha256 in the provenance file), `out/classification_summary.txt`, `out/sample_sheet_checks.txt`, the paper PDF and supplementary files, and the pages named in step 6. Results are own calculations.

## Decision: GO, with RQ3 reduced to detection level
Conditions: (a) RQ3 quantitative comparison only if the authors supply per-sample values; (b) the points under "Consequences for the plan" are written into proposal v1.3 before registration; (c) the licence points in step 6 are respected.

## Acceptance criteria (proposal section 3)
| Criterion | Result |
|---|---|
| Raw FASTQs with stable accessions and read-quality information | Met for accessions (33 paired NovaSeq runs in ENA, 164.9 GB FASTQ.gz). Read quality not yet inspected; done in WP1 on Day 1. |
| Complete sample sheet (run, matrix, stage, date, platform, replicates) | Met from ENA: `sample_alias` equals the paper ID for 33 of 33; type, date and A/B/C replicate letter are in ENA. Whether replicates are independent composites or subsamples: unknown. |
| qPCR outputs per sample | Partly met: per-sample detection calls (S5) are public; per-sample quantitative values are not (requested 2026-10-27). |
| Enough comparable samples without conflating matrix, place and date | Met with caveat: influent/effluent (19 to 21 Nov 2019) and pond/sediment (10 to 12 Dec 2019) are separate campaigns, so stage and date are confounded. Sediment has no Day 1 samples. |
| Storage and compute budget; rights recorded | Met: WSL has 619 GB free (2026-10-27, `df -h ~`); 164.9 GB FASTQ.gz; MEGAHIT on a subset; metaSPAdes only with a costed cloud run. Rights: step 6. |

## Step 5: assay data (decision on RQ3)
- Public: S5 (presence/absence per sample, MGS and qPCR; 165 cells, 154 agree, 11 discordant), S6 (subtype read counts, 33 samples x 28 subtypes), S1 to S4 (primers, standards, sampling).
- qPCR: 16S, ermB, sul1, tetA, tetQ, tetW; 5 ng input DNA; 4 technical replicates; standards from 5 x 10^6 to 5 x 10^1 copies/ul; signals below about 50 copies/ul excluded. ARG per 16S is the ARG copy number divided by the 16S copy number; for MGS it is ARG read count divided by 16S read count (paper methods and S1).
- Not public: per-sample copy numbers or ARG-per-16S values; how the 4 technical replicates form one detection call; replicate type. Non-detects are treated as censored (registered rule), not zero.
- Decision: RQ3 is detection-level (S5 versus the registered MGS rule); quantitative RQ3 is conditional on the authors' reply (no reply is the planning assumption).

## Step 6: rights (read directly 2026-10-27)
- Article: CC BY-NC-ND 4.0 (paper PDF, licence statement). The statement says adapted material derived from the article may not be shared. Whether the supplementary files carry the same licence was not confirmed from the journal page (not reachable here); treated as the same, as the proposal does.
- Reads (ENA/INSDC): INSDC policy sets no licence requirement or fee on use or redistribution of sequence records; credit by citing the original submission. The ENA policy page notes possible benefit-sharing obligations for commercial users (Cali Fund). Source: ENA and INSDC Policies page.
- NCBI: no restrictions on use of data in its databases; submitters or countries may claim IP rights, which NCBI cannot assess. Source: NCBI Website and Data Usage Policies page.
- Citation: cite Taylor et al. 2025 and the BioProject PRJNA904380 in every output.
- Practical rules: (1) public files must not contain S4, S5 or S6 content; build the public sample sheet from ENA fields only (type from the ID letter, date and replicate letter from the ENA record); (2) report own aggregate statistics with citation, not copies of the tables; (3) the `sample_sheet_draft.tsv` in `out/` uses S4 for type and day and stays local; (4) whether a public portfolio on a consultancy site counts as non-commercial use of the NC-licensed supplement is unsettled; the email asked about reuse conditions; if the answer matters, ask the publisher. This is a record of the sources read, not legal advice.

## Step 7: Alicante contingency
Not adopted (6 SRA experiments; project intro). Paper not opened in this audit. Taylor passes the gate, so the contingency is closed unless the plan changes.

## Step 8: pilot choice and consequences for the plan
- Pilot (Day 1, 6 samples, replicates A and B of influent, effluent and pond): E001, E006, R001, R006, P001, P006; 30.3 GB FASTQ.gz, 67.9 Gbp. Day 1 holds 9 samples (46.5 GB); the remaining three (replicate C) follow after the pilot.
- Day 1 has no sediment, so sediment settings cannot be tuned on development data; sediment is evaluation-only and uses settings frozen on liquid samples.
- All 11 discordant cells are in Days 2 and 3; Day 1 has 0 of 45. Development data therefore carry no information about where the methods disagree, and thresholds are set without sight of the discordances.
- Between-stage comparisons are confounded with date (two campaigns three weeks apart); every such claim says so.
- Six extra metagenomic runs (3 blanks, 3 `*000`) are outside the 33; their use as descriptive controls is decided before registration.
- Proposal v1.3 corrections: download estimate (164.9 GB), "no controls" risk (blanks exist, very shallow), the mapping statement, the Knight et al. entry (confirmed below).

## Knight et al. reference (D4): confirmed
Knight ME, Webster G, Perry WB, Baldwin A, Rushton L, Pass DA, Cross G, Durance I, Muziasari W, Kille P, Farkas K, Weightman AJ, Jones DL. National-scale antimicrobial resistance surveillance in wastewater: a comparative analysis of HT qPCR and metagenomic approaches. Water Research 2024;262:121989. doi 10.1016/j.watres.2024.121989; PMID 39018584. Wales, 47 treatment plants (weekly composite influent) plus hospital effluent; HT qPCR (73 genes) and metagenomics (545 genes). Details taken from the Bangor University research portal record, not from the journal page; accession numbers were not shown there. Reserve source only; the placeholder in v1.0 is replaced by this reference.
