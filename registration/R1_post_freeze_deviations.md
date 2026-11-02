# R1 post-freeze changes and deviations (section 6b of the registration; kept separate so the tagged text stays unchanged)

Tag `R1-registered` = commit ad97ff6 (2026-10-30 09:33). All entries below were made before any evaluation result was seen.

| # | Date | What changed | Why | Relative to results | Effect on evaluation | Where |
|---|---|---|---|---|---|---|
| 1 | 2026-10-30 to 11-01 | Evaluation reads downloaded with `03c_fetch_set.sh` (chunk size 32 MB; wrapper bug fixes), each file checked against the live ENA md5 | download over mobile broadband failed with the earlier method; transfer is allowed by the freeze checklist | before any result | none: files identical to ENA (md5) | HANDOVER addenda 27, 28 |
| 2 | 2026-11-01 | `06_trim_kma.sh`: guard changed so evaluation runs are accepted when `ALLOW_EVAL=1` and the sample sheet marks them `eval`; processing steps unchanged (git diff against ad97ff6 shows only the guard lines) | the freeze script accepted development runs only | before any result | none on counts | this file; commit following the evaluation run |
| 3 | 2026-11-01 | `10_run_set.sh` added: runs `06_trim_kma.sh` once per run of a set, one at a time | automation of the run | before any result | none | this file |
| 4 | | (append here any later change: what, why, timing, effect) | | | | |
