# HIV drug repositioning: signature reversal and metabolic modelling (GSE199911)

Dataset: **GSE199911**, acute HIV-1 infection vs healthy controls (PBMC RNA-seq)
Status: Steps 1, 2a and 2b complete. Step 3 (integration) started: client-list integration done, wider drug-target search pending.

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
         Candidate repurposed drugs                     [Step 3 - in progress]
```

---

## 2. Data

- **GSE199911**: PBMC RNA-seq (50 bp single-end, aligned to hg19). The GEO count file holds 75 samples: **29 untreated HIV-infected (about 1 month after screening) and 46 uninfected controls**.
- **Differential expression** (Step 1) was provided by the client: `Supplementary_-_1.xlsx`, sheet "HIV DEGs". It lists 379 genes (|log2FC| ≥ 1.5, padj < 0.05): 344 up, 35 down in HIV.
- **Group labels** were set from the `studygroup` column in the count file (1 = HIV, 2 = control). This matches the GEO sample annotation ("disease state: HIV infected" vs "HIV non-infected") for the 73 samples whose annotation could be read; the remaining 2 samples were labelled by `studygroup` alone.
- **Client drug list**: `NO_CI_DRUGS.xlsx` (15 unique drugs, with targets, mechanism, network proximity and contraindications).
- **Metabolic model**: Human-GEM (Human1) version 2.1.1, 12,748 reactions.

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

**Best-supported enzymes:** TK1 (thymidine kinase; DEG + reaction activity + flux, all up) and, with two agreeing sources, TYMS (thymidylate synthase), CD38, AICDA, SLC27A2, SCD, CMPK2 (all up) and AK5, SLC4A8, SLC4A10 (down). The DEG-mapped enzymes also include RRM2 and FABP5. Note that Human1 contains some DNA and protein-modification reactions, so a few of the 17 mapped DEGs (for example POLQ, UBE2C, USP18) are not classic metabolic enzymes.

**Interpretation.**
- The metabolic signal is **nucleotide synthesis** (RRM2, TYMS, TK1, plus DNA replication genes). This matches Step 2a: expanding, proliferating immune cells need more nucleotides. The raised carnitine shuttle and fatty-acid oxidation reactions point to a change in lipid handling in activated cells. CD38 is a regulator of NAD+ metabolism, and CMPK2 is an interferon-stimulated gene, in line with the interferon response in acute HIV.
- **The GIMME flux comparison is not reliable evidence and was not used as such.** Only 9 of the 106 flux changes agree in direction with the expression data and 2 oppose it. The largest changes are transport and adenylate kinase / pyruvate kinase isoenzyme swaps, and total nucleotide-metabolism flux comes out *lower* in HIV, the opposite of the expression evidence. This is typical of picking one of many equivalent optimal flux solutions. The reaction-activity analysis is the statistically sound part of Step 2b.
- Bulk PBMC data mix cell types. Part of the signal may reflect a shift in cell composition (for example expansion of proliferating plasmablasts and T cells) rather than metabolic reprogramming inside each cell.

---

## 5. Step 3 (in progress): integrating reversal and metabolic evidence for the client drug list

For the 15 drugs in `NO_CI_DRUGS.xlsx`, `results_3/5_drug_integration_client_list.csv` combines the client's signature score and network proximity, contraindications, SigCom reversal (Step 2a) and metabolic evidence for their targets (Step 2b). Five yes/no flags are used: SigCom reversal, client signature score, network proximity, no contraindications, and metabolic-target support.

| Drug | Main finding |
|---|---|
| **Fedratinib** (JAK2/FLT3 inhibitor) | Best supported: **4 of 5 flags**. Reverses the signature in 23 of 25 SigCom signatures across 12 cell lines, client NES -1.55 (FDR 0.027), proximity z -3.63, no contraindications. No Human1 metabolic target support. |
| Dexamethasone | Mimics rather than reverses (13 mimicking vs 1 reversing signature), which conflicts with its negative score in the client table. |
| Lenalidomide | Mimics (12 of 12 signatures) and has 10 contraindications. |
| Benzbromarone | Targets ABCG2 and AKR1C1, but their reaction activity is mixed or weak (one supporting source), so no metabolic support. |
| Others (oxybenzone, nicergoline, ciprofibrate, bromhexine, anastrozole) | 1–2 SigCom signatures, all mimicking. |
| Six drugs (naftidrofuryl, mephentermine, flupirtine, colforsin, zoxazolamine, ciclacillin) | Not among the top 10,000 SigCom hits, so not assessed in Step 2a. |

**Key finding:** none of the 15 drugs has metabolic-target support from Step 2b. Their targets are mostly kinases, receptors and transporters, while the HIV metabolic signal sits in nucleotide synthesis, and only 6 of their 324 targets fall in significantly changed reactions. Linking the metabolic findings to drugs needs a drug-target database that covers enzymes such as TYMS, RRM1/RRM2, TK1, CD38, CMPK2 and SCD (for example DrugCentral or DGIdb). This is the next sub-step.

---

## 6. Limitations

- **Cell-line data** (L1000) and **bulk PBMC** expression: results may not transfer to specific immune cell types or to infected cells.
- **Proliferation signal dominates.** Top reversers are generic anti-proliferative agents, not HIV-specific, and many are cytotoxic oncology drugs.
- **SigCom returns only the top 5,000 hits per direction**, so per-compound counts are lower bounds.
- **Immunoglobulin genes** (30 of the 150 up genes) are poorly measured in L1000; a sensitivity run without them is recommended.
- **Medium and growth requirement are modelling assumptions** (defined medium; 10% of maximum growth). The flux comparison is one optimal solution per condition.
- **Expression is not length-normalised** and the reaction-activity test is not adjusted for age or sex (both are available in the count file).
- **Computational predictions only**; they do not show efficacy or safety in HIV.

---

## 7. Files

| File | Description |
|---|---|
| `step2a_signature_reversal.py` | Step 2a: signature search in SigCom LINCS |
| `step2b_pipeline.py` | Step 2b: GEO download, Human1 mapping, reaction activity, GIMME/iMAT, candidate enzymes |
| `step3_integration.py` | Step 3 (part 1): integrates client drug list with Steps 2a and 2b |
| `results_v2/` | Step 2a results (`sigcom_compound_summary.csv`, `sigcom_raw_signature_hits.csv`, signature gene lists) |
| `results_2b/1_metabolic_DEGs.csv`, `1_subsystem_enrichment_DEGs.csv` | DEGs found in Human1 and their pathway enrichment |
| `results_2b/2_reaction_activity_HIV_vs_control.csv`, `2_subsystem_activity_enrichment.csv` | Reaction-level and pathway-level activity |
| `results_2b/3_flux_comparison_gimme.csv`, `3_subsystem_flux_gimme.csv` | GIMME flux comparison (exploratory) |
| `results_2b/4_candidate_target_enzymes.csv` | Candidate enzymes with evidence sources |
| `results_2b/groups_used.csv`, `sample_annotation_from_GEO.csv` | Sample labels used |
| `results_3/5_drug_integration_client_list.csv` | Per-drug evidence table for the client list |

## 8. How to reproduce

```bash
pip install numpy pandas scipy pyyaml requests openpyxl
python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --outdir results_v2
python step2b_pipeline.py --degs Supplementary_-_1.xlsx --hiv-code 1
python step3_integration.py
```

Options: `--method both` adds iMAT in Step 2b (slower); `--growth-frac 0.05` or `0.5` tests how sensitive the flux model is to the growth requirement.

## 9. Next steps

1. Add a drug-target database (DrugCentral or DGIdb) and link the 592 SigCom reversers and the candidate enzymes (TYMS, RRM1/RRM2, TK1, CD38, CMPK2, SCD) to drugs.
2. Sensitivity checks: exclude immunoglobulin genes from the signature, adjust reaction activity for age and sex, vary the growth fraction.
3. Cross-check reversers with a second tool (L2S2 or iLINCS) and resolve the dexamethasone and lenalidomide discrepancy with the client's original signature.
4. Rank final candidates, with cytotoxic or contraindicated drugs flagged separately.