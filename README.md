# HIV drug repositioning: signature reversal and metabolic modelling (GSE199911)

Dataset: **GSE199911**, acute HIV-1 infection vs healthy controls (PBMC RNA-seq)
Status: all workflow steps (1, 2a, 2b, 3, 4) and the Step 5 sensitivity analysis have been run. Results are computational hypotheses; validation steps are listed at the end.

---

## 1. Workflow

```
GSE199911 (acute HIV vs healthy control)
        |
Differential expression analysis                        [Step 1  - done]
        |
   +----+------------------------------+
   |                                   |
Transcriptomic signature          Enzyme / gene-set mapping
   |                                   |
CMap / iLINCS / L1000             Human1 GEMs (GIMME / iMAT)
(signature reversion)             (metabolic flux simulation)
   |                                   |
Reversing compounds               Target enzymes & pathways
   [Step 2a - done]                    [Step 2b - done]
   |                                   |
   +-----------------+-----------------+
                     |
         Candidate repurposed drugs                     [Steps 3 and 4 - done]
```

---

## 2. Data

- **GSE199911**: PBMC RNA-seq (50 bp single-end, aligned to hg19). The GEO count file holds 75 samples: **29 untreated HIV-infected (about 1 month after screening) and 46 uninfected controls**.
- **Differential expression** (Step 1) was provided by the client: `Supplementary_-_1.xlsx`, sheet "HIV DEGs". It lists 379 genes (|log2FC| ≥ 1.5, padj < 0.05): 344 up, 35 down in HIV.
- **Group labels** were set from the `studygroup` column in the count file (1 = HIV, 2 = control). This matches the GEO sample annotation ("disease state: HIV infected" vs "HIV non-infected") for the 73 samples whose annotation could be read; the other 2 samples were labelled by `studygroup` alone.
- **Client drug list**: `NO_CI_DRUGS.xlsx` (15 unique drugs with targets, mechanism, network proximity and contraindications).
- **Metabolic model**: Human-GEM (Human1) version 2.1.1, 12,748 reactions.
- **Drug-gene interactions**: DGIdb (87,667 interactions, 18,005 drugs; includes approval and anti-cancer flags).

---

## 3. Step 2a: signature reversal (L1000 / SigCom LINCS)

**Method.** Genes were ranked by the DESeq2 Wald statistic and only genes present in the L1000 gene space were kept (this removes LOC/LINC/antisense transcripts). The final signature was **150 up and 19 down genes**. It was searched against the LINCS L1000 chemical perturbation library through the SigCom LINCS API. A signature counted as a *reverser* when z_up ≤ -5 and FDR < 0.05 (the up-regulated HIV genes sit at the bottom of the compound's expression profile). CMap/CLUE (clue.io) was retired on 31 January 2026, so SigCom LINCS, which searches the same LINCS data, was used.

**Results.**

| Item | Value |
|---|---|
| Signature hits returned | 10,000 (5,000 reversing + 5,000 mimicking) |
| Distinct compounds / cell lines | 2,969 / 191 |
| Compounds with ≥1 reversing signature | 592 |
| Reversing in ≥3 / ≥10 cell lines | 309 / 27 |
| "Reverser" compounds (≥3 reversing signatures and more than twice as many reversing as mimicking) | 268 |

Top reversers (by number of cell lines): palbociclib (CDK4/6 inhibitor, 76 cell lines), PHA-793887 (CDK), R-547 (CDK), nutlin-3 and AMG-232 (MDM2), vorinostat (HDAC), osimertinib (EGFR), teniposide (topoisomerase II), JNJ-7706621 (CDK/Aurora), rucaparib (PARP).

**Why these results.** The up-regulated genes are dominated by cell-cycle and DNA-replication genes (CDC45, MCM2/4/10, CCNA2, TOP2A, PLK1, CDK1, FOXM1, E2F2/7/8) plus interferon-stimulated genes (IFI27, ISG15, USP18), CD38, CD8A/B and immunoglobulin variable segments. Compounds that stop cells dividing push the cell-cycle genes down, so anti-proliferative drugs rank at the top whatever the disease. The 19 down genes were not enriched in any L1000 signature (FDR ≥ 0.92), so the ranking rests on the up genes.

---

## 4. Step 2b: enzyme mapping and Human1 metabolic modelling

**Methods.**
1. **DEG mapping.** The 379 DEGs were mapped to Human1 genes, reactions and subsystems.
2. **Reaction activity.** Counts were normalised (median-of-ratios, log2), genes with fewer than 10 counts in 3 samples were removed (2,271 Human1 genes retained). For every Human1 reaction, an activity score was computed per sample from its gene rule (AND = lowest gene, OR = highest gene) and compared between HIV and control with a Welch t-test and Benjamini-Hochberg correction. Pathway enrichment counted reactions sharing the same gene rule once.
3. **Flux models.** Context-specific models for HIV and control were built with GIMME (group-median expression, 25th-percentile threshold, 10% of maximum growth required) in a defined medium, and the two flux solutions were compared.
4. **Candidate enzymes.** Each enzyme was scored by how many evidence sources agree on direction (DEG, reaction activity, GIMME flux).

**Results.**

| Item | Value |
|---|---|
| DEGs present in Human1 reactions | 17 of 379 (all up-regulated, none of the 35 down genes) |
| Reactions tested for activity | 6,752: 566 more active, 470 less active in HIV (FDR < 0.05, \|log2 diff\| ≥ 0.5) |
| GIMME flux changes | 106 of 12,748 reactions (46 higher, 60 lower in HIV) |
| Candidate enzymes with any evidence | 231 (103 up, 116 down, 12 conflicting) |

**Pathways more active in HIV** (FDR < 0.05): **pyrimidine metabolism** (FDR 0.0001), **nucleotide metabolism** (0.004) and the **mitochondrial carnitine shuttle** (0.004). No pathway was significantly less active (closest: starch and sucrose metabolism, FDR 0.079; xenobiotics metabolism, 0.097). The strongest single reaction is ribonucleotide reduction (RRM1/RRM2/RRM2B with thioredoxin), about 3.7-fold higher in HIV (log2 difference +1.89, FDR ~10⁻¹¹).

**Best-supported enzymes:** TK1 (thymidine kinase; DEG + reaction activity + flux, all up) and, with two agreeing sources, TYMS (thymidylate synthase), CD38, AICDA, SLC27A2, SCD, CMPK2 (all up) and AK5, SLC4A8, SLC4A10 (down). The DEG-mapped enzymes also include RRM2 and FABP5. Human1 contains some DNA and protein-modification reactions, so a few of the 17 mapped DEGs (for example POLQ, UBE2C, USP18) are not classic metabolic enzymes.

**Interpretation.**
- The metabolic signal is **nucleotide synthesis** (RRM2, TYMS, TK1, plus DNA replication genes). This matches Step 2a: expanding, proliferating immune cells need more nucleotides. The raised carnitine shuttle and fatty-acid oxidation reactions point to a change in lipid handling in activated cells. CD38 regulates NAD+ metabolism, and CMPK2 is an interferon-stimulated gene, in line with the interferon response in acute HIV.
- **The GIMME flux comparison is not reliable evidence and was not used as such.** Only 9 of the 106 flux changes agree in direction with the expression data and 2 oppose it. The largest changes are transport and isoenzyme swaps, and total nucleotide-metabolism flux comes out *lower* in HIV, the opposite of the expression evidence. This is typical of picking one of many equivalent optimal flux solutions. The reaction-activity analysis is the statistically sound part of Step 2b.
- Bulk PBMC data mix cell types. Part of the signal may reflect a shift in cell composition (for example expansion of proliferating plasmablasts and T cells) rather than metabolic reprogramming inside each cell.

---

## 5. Step 3: linking drugs to the metabolic findings

### 3a. Client drug list (15 drugs)
`results_3/5_drug_integration_client_list.csv` combines the client's signature score and network proximity, contraindications, SigCom reversal and Step 2b support for the targets.

| Drug | Main finding |
|---|---|
| **Fedratinib** (JAK2/FLT3 inhibitor) | Best supported: reverses the signature in 23 of 25 SigCom signatures across 12 cell lines, client NES -1.55 (FDR 0.027), proximity z -3.63, no contraindications. No Human1 metabolic target. |
| Dexamethasone | Mimics rather than reverses (13 mimicking vs 1 reversing), conflicting with its negative score in the client table. |
| Lenalidomide | Mimics (12 of 12 signatures) and has 10 contraindications. |
| Benzbromarone | Targets ABCG2 and AKR1C1, but their reaction activity is mixed or weak. |
| Other 10 drugs | Few or no SigCom signatures; 6 of them are not among the top 10,000 hits, so were not assessed. |

None of the 15 drugs has metabolic-target support: their targets are mostly kinases, receptors and transporters, while the HIV metabolic signal sits in nucleotide synthesis.

### 3b. Drugs that act on the HIV-changed enzymes (DGIdb)
The 36 Step 2b target genes (33 up, 3 down in HIV) were looked up in DGIdb, checking whether a drug acts in the opposing direction (an inhibitor of a gene that is up in HIV).

| Item | Value |
|---|---|
| Drugs in DGIdb hitting these genes | 197 (most targeted genes: TYMS 47 drugs, ADA 23, CTPS1 19, NT5E 19, RRM2 13, CD38 8, SCD 7, TK1 6) |
| With a documented inhibitory action on the gene | 56 (the other 139 have no action type recorded) |
| Also present in the SigCom hits | 28; 9 of them are signature reversers |
| **Reverse the signature and inhibit the target** | **1: fluorouracil (TYMS inhibitor): 8 reversing, 0 mimicking signatures, 3 cell lines** |

Other direct nucleotide-synthesis inhibitors did **not** reverse the signature consistently: gemcitabine (3 reversing vs 4 mimicking), clofarabine (2 vs 3), floxuridine (4 vs 7), pemetrexed (1 vs 11), raltitrexed (1 vs 3). One possible explanation, not tested, is that cells arrested in S phase by antimetabolites keep the replication genes switched on. Most of the remaining 48 direct inhibitors (for example hydroxyurea, capecitabine, CD38 antibodies) have no L1000 signature among the top hits, so their effect on this signature is **unknown**, not negative.

**Important caveats.** (i) TYMS, RRM2 and TK1 are themselves part of the up-regulated signature, so a drug that inhibits them and reverses the signature is not independent confirmation. (ii) Some DGIdb "inhibitor" claims come from binding screens and are doubtful (for example ROCK inhibitors listed against CTPS1). (iii) Reversers with a blank interaction type (etoposide, topotecan and others linked to TYMS) are not direct inhibitors of that gene and were not counted as consistent.

---

## 6. Step 4: final candidate ranking

`results_4/7_final_candidate_ranking.csv` ranks 337 candidates (SigCom reversers, drugs hitting Step 2b enzymes, and the client list). Each drug gets simple evidence points (not a statistical score): reversal (2 points, +1 if in ≥10 cell lines), consistent action on a Step 2b enzyme (2 points, +1 if the gene is best supported), and for client-list drugs one point each for a negative signature score, network proximity z ≤ -1.96 and no contraindications. Drug class comes from DGIdb flags (approved, anti-cancer, immunotherapy).

| Tier | Meaning | Drugs |
|---|---|---|
| 1 | Reverses the signature and has independent support | **2: fedratinib, fluorouracil** |
| 2 | Strong reversal only (≥10 cell lines) | 23 (palbociclib, PHA-793887, R-547, nutlin-3, vorinostat, AMG-232, JNJ-7706621, rucaparib and others) |
| 2b | Reversal in fewer cell lines | 243 |
| 3 | Acts on a Step 2b enzyme, reversal not shown | 55 |
| 4 | Given list, no reversal | 14 |

**Drug classes among the 268 signature reversers:** 62 anti-cancer drugs, 20 approved non-oncology drugs, 57 investigational or unannotated, and 129 not found in DGIdb (mostly research compounds or lab codes).

**Repurposing shortlist** (`results_4/7_repurposing_shortlist.csv`): the 20 approved, non-anti-cancer drugs that reverse the signature. The strongest are pyrvinium pamoate (20 reversing vs 0 mimicking, 12 cell lines), pitavastatin (16 vs 1, 6 cell lines), hexachlorophene, sertraline, niclosamide, etravirine, lomitapide, fluspirilene and ivermectin. Etravirine is an HIV antiretroviral, but the cell lines are not HIV-infected, so its reversal is probably an off-target effect, and it was not robust in the sensitivity analysis (Section 7). Of the 19 distinct shortlist drugs, 12 stay reversers in every sensitivity run.

**Interpretation.**
1. **Fedratinib is the single best-supported candidate.** It reverses the signature consistently, agrees with the client's independent score, has good network proximity and no recorded contraindications. It acts on JAK/STAT signalling, which fits the interferon and cytokine biology of acute infection. It is approved for myelofibrosis and carries a safety warning that must be checked against current labelling before any suggestion of repurposing.
2. **Fluorouracil is mechanistically coherent** (it inhibits TYMS, which is raised in HIV, and reverses the signature) but it is a cytotoxic cancer drug. The argument is partly circular (see the Step 3b caveats), although its reversal persists when TYMS and RRM2 are removed from the signature (Section 7).
3. **The non-oncology shortlist is exploratory.** Many of these drugs (pyrvinium, niclosamide, hexachlorophene, sertraline, fluspirilene, amiodarone) are known to affect many cell types through general growth inhibition or cell stress, so their reversal may reflect that rather than anything HIV-specific. The DGIdb anti-cancer flag is also imperfect (for example ixabepilone and colchicine appear in the shortlist), so each drug needs manual review.
4. **No drug has strong, HIV-specific, independent evidence.** The conclusions are hypotheses for follow-up, not recommendations.

---

## 7. Step 5: sensitivity analysis

**Why.** The signature is dominated by cell-cycle genes and contains 30 immunoglobulin segment genes that L1000 measures poorly. The question was whether the top candidates still reverse the signature without those genes.

**Method.** The Step 2a search was repeated with a modified signature. Each time the signature was rebuilt from the next-ranked genes, so it still had 150 up genes. A candidate counts as *robust* if it is a reverser in the baseline and in all three variant runs.
- `noIG`: the 30 immunoglobulin genes removed.
- `noProlif`: 43 core cell-cycle genes removed (published S-phase and G2/M sets), including TYMS and RRM2.
- `noIG_noProlif`: both (73 of the 150 original up genes removed).

**Results: overall reversers.**

| Run | Reversers | Shared with baseline | Baseline reversers retained | New reversers |
|---|---|---|---|---|
| baseline | 268 | 268 | 100% | 0 |
| noIG | 269 | 234 | 87.3% | 35 |
| noProlif | 288 | 214 | 79.9% | 74 |
| noIG_noProlif | 277 | 194 | 72.4% | 83 |

**Results: candidates.**

| Group | Robust in all 4 runs |
|---|---|
| Tier 1 (fedratinib, fluorouracil) | 2 of 2 |
| Tier 2 (strong reversal, ≥10 cell lines) | 23 of 23 |
| Repurposing shortlist | 12 of 19 |

- **Fedratinib** reverses in every run: 20–23 reversing signatures in 10–12 cell lines each time, with 0–2 mimicking signatures. The reversal does not depend on the immunoglobulin or core cell-cycle genes.
- **Fluorouracil** still reverses (6–9 reversing signatures, 2–3 cell lines) even with TYMS and RRM2 removed from the signature. This answers part of the circularity concern from Step 3b: its reversal is not simply due to its own target gene being in the signature.
- **Tier 2 drugs all stay reversers**, although the strength falls for some (palbociclib goes from 76 to 42 cell lines in the combined run).
- **Shortlist, robust (12):** pitavastatin, hexachlorophene, sertraline, niclosamide, ivermectin, ingenol mebutate, bazedoxifene, amiodarone, nifedipine, perphenazine, doconexent, raloxifene. **Not robust (7):** lomitapide (3 of 4 runs), fluspirilene, colchicine and ixabepilone (2 of 4), and etravirine, nortriptyline and carvedilol (baseline only). Etravirine losing its reversal supports the view that it was an off-target effect.

**Interpretation.**
- The leading candidates do not depend on the immunoglobulin genes or on the 43 core cell-cycle genes, so their reversal is not an artefact of those genes.
- This does **not** show that they act on anything other than proliferation. Many proliferation genes remain in the rebuilt signature (for example KIF15, PBK, MELK, MCM10), and the strongest reversers are broadly growth-inhibitory drugs. A stricter test removes the whole proliferation programme (see next steps).
- About a fifth to a quarter of the baseline reversers drop out and 35–83 new ones appear in each variant, so the full reverser list is moderately sensitive to the signature; the robust candidates above are the ones to rely on.

---

## 8. Limitations

- **Cell-line data** (L1000) and **bulk PBMC** expression: results may not transfer to specific immune cell types or to HIV-infected cells.
- **Proliferation signal dominates**, so top reversers are generic anti-proliferative agents.
- **SigCom returns only the top 5,000 hits per direction**, so per-compound counts are lower bounds and drugs absent from the hits were not actually tested.
- **Immunoglobulin genes** (30 of the 150 up genes) are poorly measured in L1000; the sensitivity analysis in Step 5 tests the effect of removing them.
- **Medium and growth requirement are modelling assumptions**; the flux comparison is one optimal solution per condition and was not used as evidence.
- **Expression is not length-normalised** and the reaction-activity test is not adjusted for age or sex (both are available in the count file).
- **DGIdb annotations are incomplete and sometimes noisy**; 129 reversers could not be matched by name.
- **Computational predictions only**; they do not show efficacy or safety in HIV.

---

## 9. Files

| File | Description |
|---|---|
| `step2a_signature_reversal.py` | Step 2a: signature search in SigCom LINCS |
| `step2b_pipeline.py` | Step 2b: GEO download, Human1 mapping, reaction activity, GIMME/iMAT, candidate enzymes |
| `step3_integration.py` | Step 3a: client drug list integrated with Steps 2a and 2b |
| `step3b_drug_target_link.py` | Step 3b: DGIdb drug-target linking |
| `step4_final_ranking.py` | Step 4: final ranking and repurposing shortlist |
| `step5_sensitivity_compare.py` | Step 5: compares baseline and sensitivity runs |
| `results_v2/` | Step 2a results (compound summary, raw hits, signature gene lists) |
| `results_2b/` | Step 2b results: `1_*` DEG mapping, `2_*` reaction activity, `3_*` flux (exploratory), `4_candidate_target_enzymes.csv`, `groups_used.csv` |
| `results_3/` | `5_drug_integration_client_list.csv`; `6_drugs_targeting_HIV_metabolic_genes.csv`; `6_sigcom_reversers_with_targets.csv` |
| `results_4/` | `7_final_candidate_ranking.csv`; `7_repurposing_shortlist.csv` |
| `results_5/` | `8_sensitivity_robustness.csv`; `8_sensitivity_run_summary.csv` |
| `make_exclude_file.py` | Helper: builds a broad proliferation gene list from an MSigDB Hallmark file for a stricter sensitivity run |

## 10. How to reproduce

```bash
pip install numpy pandas scipy pyyaml requests openpyxl
python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --outdir results_v2
python step2b_pipeline.py --degs Supplementary_-_1.xlsx --hiv-code 1
python step3_integration.py
python step3b_drug_target_link.py
python step4_final_ranking.py

# Step 5: sensitivity analysis
python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --exclude-ig --outdir results_noIG
python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --exclude-proliferation --outdir results_noProlif
python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --exclude-ig --exclude-proliferation --outdir results_noIG_noProlif
python step5_sensitivity_compare.py \
    --variant noIG=results_noIG/sigcom_compound_summary.csv \
    --variant noProlif=results_noProlif/sigcom_compound_summary.csv \
    --variant noIG_noProlif=results_noIG_noProlif/sigcom_compound_summary.csv
```

Options: `--method both` adds iMAT in Step 2b (slower); `--growth-frac 0.05` or `0.5` tests the sensitivity of the flux model.

## 11. Next steps

1. **Stricter proliferation test.** Remove the whole proliferation programme (MSigDB Hallmark E2F targets, G2M checkpoint and mitotic spindle, a few hundred genes) with `make_exclude_file.py` and `--exclude-file`, and see which candidates still reverse the remaining signature. This is the real test of whether anything beyond cell division is being reversed.
2. **Test the drugs that were never assessed.** Look up the L1000 signatures of the 48 enzyme-targeting drugs that are missing from the top hits (for example hydroxyurea, capecitabine) instead of relying on the top-10,000 list.
3. **Cross-check with a second tool** (L2S2 or iLINCS), and resolve the dexamethasone and lenalidomide discrepancy with the client's original signature.
4. **Replicate the HIV signature in an independent cohort** before drawing biological conclusions.
5. **Manual review** of the top candidates: literature on HIV, safety and contraindications, and drug exposure in immune cells.
