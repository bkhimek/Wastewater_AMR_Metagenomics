# ermB diagnostic note (exploratory, post-freeze; not part of the R1 verdict)

Date: 2026-11-02. Follows the independent review of 2026-11-02 (steps A to C, stop rule). Registered outputs in `results/kma/` were never changed; all follow-up runs wrote to `results/explore/`. Versions: KMA 1.4.9 (flags `-ipe -mem_mode -ef -1t1 -cge -nf -vcf`), ResFinder db tag 2.1.1 (`all.fsa`, 3,146 sequences), bbduk 38.90 for the registered trimmed runs. Scripts: `12_ermB_diagnostic.py`, `13_explore_raw.sh`, `14_explore_single_template.sh` (commit following this note).

## Question
In 3 of 9 effluent samples (E011 dev, E021, E036) our ermB gene total is far below the published value (gene ratios 0.73, 0.22, 0.71); the other six agree within 2-4% and the other four genes agree in the same samples.

## Findings
1. **Step A (existing outputs, all nine effluent samples).** Two kinds of erm(B) counts. Whole-gene support: E001 erm(B)_18 570 reads (coverage 100%, depth 4.6), E011 erm(B)_18 557 (99.7%, 4.3), E036 erm(B)_10 548 (99.2%, 4.1), E031 erm(B)_10 97 (98.9%), E041 erm(B)_1 117 (98.6%). Thin blocks: erm(B)_2 with 672-780 reads but only 26-52% coverage, depth 0.4-1.0 and 2-5 aligned reads (E006, E016, E026, E031, E041). Non-ermB templates show the same pattern in every sample (erm(38) 2,856-3,534 reads at 14-26% coverage; erm(Q) 1,460-1,980; erm(F), erm(G) 600-850), so `readCount` behaves as a voting count for these templates, not as reads aligned along the gene.
2. **Step C1 (raw reads, no trimming), E021 vs control E026.** E021: erm(B)_2 still absent (erm(B)_18 112 vs 101 trimmed; erm(G)_1 815 vs 775). E026: erm(B)_2 816 raw vs 698 trimmed. Trimming changes counts by about 2-17% and does not create or remove the block.
3. **Step C2 (single-template index with erm(B)_2 only, raw reads).** E021: readCount 899, fragmentCount 893, 738 of 738 positions covered, 21 aligned reads. E026: 906, 897, 738 of 738, 25 aligned. With no competing template the two samples are alike; in the full database E021 has 0 reads on erm(B)_2 and E026 has 698.

## Status of each hypothesis
| Hypothesis | Status | Evidence |
|---|---|---|
| KMA reference competition / ambiguous allele assignment | **Supported** | C2: reads exist in E021 at the level of the control when erm(B)_2 stands alone; in the full database they are labelled erm(B)_18/_20, erm(G)_1, erm(F)_3 |
| Assumed bbduk adapter options | **Not supported for the missing block** | C1: raw reads give the same absence in E021 (counts +2-17%) |
| Reference content / index differences from the authors' run | **Unresolved** | cannot be tested without the authors' exact files; would change the competition outcome |
| Sample/run mix-up, hidden unpaired reads | **Not supported** | other four genes agree; all runs have exactly two files; md5 verified |
| Error in Table S6 | **Not supported** | published ermB/sul1 ratio of E021 (0.61) lies inside the normal range 0.44-0.78; the shared value 659 (E011, E021) is unexplained but not evidence of a duplicated column |

## Caution
The single-template count is an upper bound: with one template, any sufficiently similar read (including reads from related erm genes) can map to it. It shows that sequence similar to erm(B)_2 is present, not its true abundance. The thin erm(B)_2 blocks (about 700 reads at 26-52% coverage) are not whole-gene evidence; the paper text says a whole-gene coverage threshold was required, yet its published counts contain such low-support blocks.

## Consequence for the report
The R1 verdict stands unchanged. The ermB shortfall is reported as a template-competition effect among near-identical erm references, with the whole-gene-supported reads (about 550 in E011 and E036) present but labelled on a sibling template. Subtype-level and ermB-level counts below whole-gene support should not be interpreted as abundance. Stop rule applied: explanation recorded as supported (competition), not supported (trimming), unresolved (reference/index difference).
