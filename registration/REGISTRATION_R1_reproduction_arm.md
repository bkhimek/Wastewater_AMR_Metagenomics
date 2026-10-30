# Registration R1: reproduction arm (version 2, candidate for freeze)

Status: version 2 of 2026-10-30, written after an independent review of the first draft and after the Day-1 development results. Becomes binding when committed and tagged `R1-registered` BEFORE any Day 2/3 or sediment file is opened. Any change after the tag goes into the post-freeze deviation table (section 6b).

## 1. Header (copied to the top of the evaluation scoring script)
```
# REGISTRATION R1 v2 (2026-10-30)
# Question:      Does a pipeline built from the authors' stated tools and settings reproduce the
#                published metagenomic ARG read counts (Taylor et al. 2025, Table S6) on held-out samples?
# Arms:          A = reproduction arm (section 2). Exploratory, outside the verdict: fastp instead of bbduk;
#                ResFinder head (2.6.0) instead of 2.1.1; 16S sensitivity rules (section 4).
# Data:          dev  = Day 1, 9 liquid samples (INF R001/R006/R011, EFF E001/E006/E011, POND P001/P006/P011).
#                eval = Days 2 and 3, 18 liquid samples (verdict). Sediment (6 samples): descriptive only.
# Setting rule:  open choices were set on dev only (16S rule D13, thresholds below); nothing changes after the freeze.
# Freeze:        tag R1-registered, created before eval data are opened.
# Predictions:   P1 to P4 decide the verdict; P3c and P5 are reported and do not decide it (section 3).
# Eval policy:   eval scored once per arm; outputs never overwritten; failures reported.
# Not in verdict: ARG/16S ratio versus Fig. 1c; sediment; subtype-level recall (P3c); P5; rank correlation;
#                fastp and ResFinder-head runs.
```

## 2. Arm A: pipeline and stated assumptions
Per sample: raw reads -> bbduk 38.90 (`maxns=0 minlength=75 qtrim=r trimq=20 ref=adapters`, paired outputs so unpaired reads are dropped) -> KMA 1.4.9 (`-ipe -mem_mode -ef -1t1 -cge -nf -vcf`, 4 threads) against ResFinder db 2.1.1 (tag, `all.fsa`, 3,146 sequences) and SILVA 138.1 SSURef NR99 (RNA converted to DNA). Script: `06_trim_kma.sh` (commit hash added at freeze).

The SILVA index was built on a temporary cloud machine (37.6 GB needed; local cap 28 GB) and verified by checksum after transfer. The cross-machine equivalence control (cloud build versus local build) was done with the ResFinder index only; for SILVA the evidence is the checksum and a working end-to-end run on all 9 dev samples.

Assumptions the authors do not state (recorded, not tested):
1. Adapter options for bbduk (`ktrim=r k=23 mink=11 hdist=1 tbo tpe`); the paper gives only `ref=adapters.fa`.
2. SILVA variant NR99 (not full Ref); ResFinder index built from `all.fsa` of tag 2.1.1; KMA index settings are the defaults.
3. "Read count" = KMA `readCount` (28/28 agreement on E001; HANDOVER addenda 17 and 20).
4. Subtype filter: count >= 50, inclusive. Evidence: Table S6 contains an exact 50 (R006, tet(W)_1_DQ060146) while the supplement text says "greater than 50"; the table is taken as authoritative.
5. No template coverage or identity filter on the ARG side (subtypes are published with large counts at 8% template coverage).

## 3. Predictions (numeric, written before any eval result)
Unit: eval liquid samples (n = 18). "Paper" = Supplementary Table S6: 28 subtypes of ermB, sul1, tetA, tetQ, tetW. Definitions are those in the docstring of `09_score_R1.py`.

Two numerators are used and kept apart:
- S6 numerator (reproduction of Table S6): KMA `readCount` of the template named as each published subtype, called when >= 50.
- G numerator (own rule, used for the gene level and for ARG/16S): per gene family, the sum of `readCount` over ALL ResFinder templates of that gene with `readCount` >= 50, whether or not the subtype is among the 28 published.

| ID | Prediction (verdict criteria) | Dev result (Day 1, n = 9) | Eval result | Met? |
|---|---|---|---|---|
| P1 | Per sample, (sum of our S6 counts >= 50) / (sum of published counts) lies in 0.93 to 1.07 in at least 16 of 18 samples, and in at least 5 of 6 samples of each type (INF, EFF, POND). Per-gene ratios reported. | 0.967 to 0.986 (9 of 9) | | |
| P2 | The same ratio is below 1.00 in at least 15 of 18 samples. Direction only; the attribution to the assumed adapter options stays a hypothesis (testing it needs a controlled trimming comparison). | 9 of 9 | | |
| P3 | Gene level: over gene-sample pairs with a published value, the ratio (G numerator / published sum of that gene) lies in 0.90 to 1.10 in at least 90% of the pairs. | 35 of 36 pairs (97%); one outside: E011 ermB 0.73 | | |
| P3b | Gene-level false positives (G numerator > 0 where the paper reports nothing for that gene in that sample): at most 2 in total. | 0 | | |
| P4 | Median of (our S6 count / published count) over published-positive subtypes (a miss counts as ratio 0) lies in 0.90 to 1.10 in at least 16 of 18 samples; the share of subtypes within +-10% is reported. | 0.98 to 1.05 (median) | | |

Reported, outside the verdict:
- P3c, subtype level: pooled recall of published-positive subtypes (registered level >= 0.90; dev 0.955, 5 misses in 9 samples) and extra positive subtypes (registered level <= 18 in total; dev 5). Reason it does not decide the verdict: the misses are large counts (about 660 to 1,700 reads) swapped between near-identical alleles (tetA(P)_1 vs tetA(P)_3 differ at 2 of 1,263 positions; tet(A)_4 vs tet(A)_6 are 95% similar; erm(B)_2 vs erm(B)_18 98% similar), so subtype labels inside these groups are unstable with short reads. A null pipeline would fail P3 and P1.
- P5, secondary: share of the SILVA readCount on Eukaryota, Chloroplast and Mitochondria templates ("excluded-template share"), expected order POND > EFF > INF on Day 2 and Day 3 (mean of replicates). Dev: INF 0.12, EFF 0.29, POND 0.57. Motivates D13. The INF/EFF and POND samples come from different campaign dates, so the order is not attributed to treatment stage.

Reporting rules: every prediction is reported whether met or not; sample-level values go in a table; P1 to P4 also per sample type.

Known dev observation, recorded before the freeze: E011 ermB is 21% below the published value even when all 601 `erm(B)` reads in the sample are counted; the cause is unexplained. It is the only pair outside the band on dev.

## 4. Counting definitions (decisions D12 and D13)
- ARG numerator for ARG/16S: the G numerator of section 3. The supplement describes summing at gene level and does not say whether the 50 filter was applied before summing; the rule here is explicit and sensitivity variants are reported.
- 16S denominator (primary): sum of `readCount` over SILVA templates labelled Bacteria or Archaea, excluding Chloroplast, Mitochondria and Eukaryota.
- Sensitivity rules, all reported (descriptive): all templates; Bacteria only; fragmentCount; readCountAln; 16S template coverage >= 90%; organelles included.
- Rules were fixed before any eval result. No denominator is chosen by agreement with Fig. 1c on eval.
- The published ratio (Fig. 1c, about 2 for EFF) could not be reproduced with any of about 100 rules tried on dev. This is an unresolved discrepancy and will be stated prominently in the final report as a limit of the reproduction.

## 5. Matrices
Sediment (6 samples) is excluded from the verdict because the development set contains no sediment, so nothing was tuned or checked on that matrix. It is analysed and described only.

## 6a. Development decisions and pre-freeze changes
| What | Why | Timing | Effect on eval | Where |
|---|---|---|---|---|
| Dev samples inspected and compared with Table S6 | needed to set thresholds on dev only | before registration (dev only) | not applicable | this file; HANDOVER addenda 17 to 22 |
| SRR22752036_1 obtained with a chunked download after md5 failures | server error text inside the data stream (root cause found) | before analysis | none expected: md5 matches ENA | HANDOVER addendum 11 |
| SILVA index built on a temporary cloud machine | needs 37.6 GB | before analysis | unknown; checksum verified, ResFinder control identical | HANDOVER addenda 12, 13 |
| Fig. 1c ratio not reproduced; own 16S rule (D13) | no tested rule fits; rule chosen by reasoning | on dev, before eval | unknown | HANDOVER addenda 20 to 22 |
| P3 redefined (positive-call recovery, then gene level) | the first definition could be passed by a pipeline that detects nothing; subtype swaps inside near-identical alleles | after dev results, before eval | unknown | HANDOVER addenda 24, 25 |

## 6b. Post-freeze deviations
(empty; every entry: what changed, why, when relative to seeing eval results, effect, where disclosed)

## 7. Freeze checklist
1. `09_score_R1.py` passes `--selftest` including the synthetic test of counts 0, 49, 50, 51; it refuses to score unless exactly 18 unique expected liquid eval samples with complete result files are present (no missing, duplicate or unexpected samples).
2. Record versions and checksums: scripts, KMA, bbduk, indexes, `all.fsa`, Table S6 CSV, sample sheet.
3. Commit scripts, this registration and the handover files; tag `phaseA-go` and `R1-registered`; push. The Table S6 CSV is not committed (licence unconfirmed).
4. Download and checksum of the eval data are allowed before the freeze (they are transfer, not analysis); do not open, count or QC them.
5. Run the pipeline on eval once; run `09_score_R1.py --set eval` once; fill the eval column; report every prediction.

## 8. Limits stated in advance
- The hold-out is limited: Table S6 and the paper's results were known before registration, and Day 1 and Days 2/3 come from the same plant and season.
- The 18 eval samples are not independent experiments (3 types x 3 replicates x 2 days, collected on different campaign dates).
- Subtype labels among near-identical alleles are not reliably identifiable from short reads.
- The authors' 16S denominator is not reproduced.
