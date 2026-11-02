# R1 results: evaluation set (Days 2 and 3, 18 liquid samples), scored once

Scored 2026-11-02 ~06:45 with `bin/09_score_R1.py --set eval --detail` (identical to the committed script, checked with git diff against ad97ff6; selftest passed). Raw output: `score_R1_eval_output.txt` (same content as `results/score_R1_eval.txt`). Pipeline: `10_run_set.sh` / `06_trim_kma.sh`, once per sample, 24 of 24 runs complete, no failures, no reruns.

| ID | Registered criterion | Dev (n=9) | Eval (n=18) | Met? |
|---|---|---|---|---|
| P1 | total ratio in 0.93-1.07 in >= 16 of 18 and >= 5 of 6 per type | 9 of 9 | 17 of 18 (INF 6/6, EFF 5/6, POND 6/6); range 0.911-0.991; the one outside: E021 0.911 | MET |
| P2 | ratio < 1.00 in >= 15 of 18 | 9 of 9 | 18 of 18 | MET |
| P3 | gene-sample pairs in 0.90-1.10 >= 90% | 35/36 (97%) | 68/70 (97.1%); outside: E021 ermB 0.22, E036 ermB 0.71 | MET |
| P3b | false genes <= 2 | 0 | 0 | MET |
| P4 | median ours/published in 0.90-1.10 in >= 16 of 18 | 9 of 9 | 17 of 18; the one outside: P036 (median 0.00; 3 published positives, 2 swapped between sibling alleles) | MET |
| P3c (reported, outside verdict) | recall >= 0.90, extras <= 18 | 0.955, 5 | 0.934, 9 | (level reached) |
| P5 (secondary, outside verdict) | POND > EFF > INF on both days | one day | Day 2: INF 0.10, EFF 0.29, POND 0.58; Day 3: 0.11, 0.27, 0.57; ordering as predicted on both days | MET |

## Observations to report (not tested predictions)
- All numerator-level predictions were met. Total ratios are 0.91-0.99 of the published counts, always below 1 (P2).
- Subtype labels swap between near-identical alleles, as seen on dev: tet(W)_1 vs tet(W)_5 (E016), sul1_9 vs sul1_2 and tet(A)_6 vs tet(A)_4 (P036), erm(B)_2 vs erm(B)_10/_20 (E036, E021). P036 has only 3 published positives, so its subtype-level median collapses to 0.
- Recurrent ermB shortfall in effluent: gene ratio 0.73 (E011, dev), 0.22 (E021), 0.71 (E036); other effluent samples 0.94-0.98; influent and pond not affected. Cause unexplained. E011 and E021 both have published erm(B)_2 = 659 in Table S6 (to be checked against the supplement file).
- Influent erm(B)_1 and erm(B)_10: published counts (e.g. 1575 in R031) appear as small counts here (22, 14), but the gene totals match (ermB 0.96-0.99), so the reads sit on other erm(B) templates.
- Share of subtypes within +-10% of the published count (reported, not a criterion): INF 0.47-0.68, EFF 0.62-0.93, POND 0-0.75; the median is close to 1, individual subtype counts vary more.
- Sediment (6 samples) was run but not scored (descriptive only).
