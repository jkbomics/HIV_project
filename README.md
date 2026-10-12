# HIV_project# HIV drug repositioning: transcriptomic signature reversal (Steps 1–2a)

Dataset: **GSE199911** (acute HIV vs healthy control)
Status: Steps 1 and 2a complete. Step 2b (Human1 metabolic modelling) and Step 3 (integration) pending.

---

## 1. Workflow overview

```
GSE199911 (acute HIV vs healthy control)
        |
Differential expression analysis                 [Step 1  - done]
        |
   +----+-----------------------------+
   |                                  |
Transcriptomic signature         Enzyme / gene-set mapping
   |                                  |
CMap / iLINCS / L1000            Human1 GEMs (GIMME / iMAT)
(signature reversion)            (metabolic flux simulation)
   |                                  |
Reversing compounds              Target enzymes & pathways
   [Step 2a - done]                   [Step 2b - pending]
   |                                  |
   +----------------+-----------------+
                    |
        Candidate repurposed drugs               [Step 3 - pending]
```

This document covers the left branch (Steps 1 and 2a).

---

## 2. What was done

### Step 1: Differential expression
- Differential expression results for GSE199911 (acute HIV vs healthy control) were taken from `Supplementary_-_1.xlsx`, sheet "HIV DEGs" (DESeq2-style output with Entrez gene IDs).
- **379 differentially expressed genes**: 344 up-regulated and 35 down-regulated in acute HIV (all with |log2FC| ≥ 1.5 and padj < 0.05).

### Step 2a: Signature reversion against L1000
The goal is to find compounds whose effect on gene expression is the *opposite* of the HIV signature.

1. **Gene ID conversion.** Entrez IDs were converted to official gene symbols using the NCBI `Homo_sapiens.gene_info` file.
2. **Signature construction.** Genes were ranked by the DESeq2 Wald statistic. Only genes present in the L1000 gene space (as resolved by the SigCom LINCS gene dictionary) were kept, so that the signature is not filled with genes that cannot be measured (LOC*, LINC*, antisense transcripts, pseudogenes). The final signature is:
   - **150 up-regulated genes** (top 150 by Wald stat)
   - **19 down-regulated genes** (all down genes that exist in the L1000 gene space; 16 of the 35 down-regulated DEGs were lost because they are non-coding or unmeasured)
3. **Signature search.** The signature was submitted to the **SigCom LINCS** API (two-sided rank enrichment, `ranktwosided`) against the LINCS L1000 **chemical perturbation** library (`l1000_cp`). The API returns the 5,000 strongest reversing and 5,000 strongest mimicking signatures.
4. **Classification of hits.** A signature was called a **reverser** when `z_up ≤ -5` and `FDR_up < 0.05` (the up-regulated HIV genes sit at the *bottom* of that compound's expression signature, meaning the compound pushes them down). A **mimicker** is the opposite (`z_up ≥ 5`).
5. **Per-compound summary.** For each compound: number of reversing signatures, number of mimicking signatures, number of distinct cell lines in which it reverses, and median z_up of its reversing signatures. Compounds were ranked by the number of cell lines in which they reverse, then by median z_up.
6. **Cross-check with the client drug list.** Results were matched to `NO_CI_DRUGS.xlsx` by compound ID/name.

> CMap/CLUE (clue.io) was originally planned but the site and its query tools were retired on 31 January 2026. SigCom LINCS searches the same underlying LINCS L1000 data and was used instead.

---

## 3. Results

| Item | Value |
|---|---|
| DEGs used | 379 (344 up, 35 down) |
| Signature submitted | 150 up + 19 down genes (all found in L1000 space) |
| Signature hits returned | 10,000 (5,000 reversing + 5,000 mimicking) |
| Distinct compounds / cell lines in hits | 2,969 compounds / 191 cell lines |
| Compounds with ≥1 reversing signature | 592 |
| Compounds reversing in ≥3 cell lines | 309 |
| Compounds reversing in ≥10 cell lines | 27 |

**Top reversing compounds (by number of cell lines):**

| Compound | Class | Reversing / mimicking signatures | Cell lines reversing | Median z_up (reversers) |
|---|---|---|---|---|
| palbociclib | CDK4/6 inhibitor | 262 / 0 | 76 | -11.7 |
| PHA-793887 | CDK inhibitor | 64 / 0 | 39 | -10.6 |
| R-547 | CDK inhibitor | 82 / 0 | 38 | -11.9 |
| nutlin-3 | MDM2 inhibitor | 69 / 1 | 31 | -11.1 |
| AMG-232 | MDM2 inhibitor | 51 / 3 | 20 | -11.6 |
| vorinostat | HDAC inhibitor | 71 / 3 | 20 | -10.6 |
| AZD-9291 (osimertinib) | EGFR inhibitor | 24 / 4 | 16 | -11.5 |
| teniposide | Topoisomerase II inhibitor | 26 / 0 | 16 | -11.0 |
| JNJ-7706621 | CDK/Aurora inhibitor | 16 / 2 | 16 | -10.4 |
| rucaparib | PARP inhibitor | 16 / 1 | 16 | -9.6 |

The top 50 are dominated by CDK inhibitors, MDM2/p53 activators, HDAC inhibitors, topoisomerase poisons, PARP inhibitors, PI3K/AKT/mTOR inhibitors, MEK/RAF inhibitors, HSP90 and proteasome inhibitors. The full table is in `results_v2/sigcom_compound_summary.csv`.

**Overlap with the client drug list (`NO_CI_DRUGS.xlsx`):**

| Drug | Result in this analysis | Existing values in client table |
|---|---|---|
| **Fedratinib (TG-101348)**, JAK2/FLT3 inhibitor | **Reverses**: 23 reversing vs 2 mimicking signatures, 12 cell lines, median z_up -11.4 | `hiv_mean_nes_up` -1.55, FDR 0.027 |
| Dexamethasone | Mostly **mimics** (13 mimicking vs 1 reversing) | `hiv_mean_nes_up` -1.70, FDR 0.010 |
| Lenalidomide | **Mimics** only (12 of 12 signatures) | `hiv_mean_nes_up` -1.34, FDR 0.013 |
| Ruxolitinib, anastrozole | Mimic only (3 of 3, 2 of 2) | n/a |

Only two compounds from the client list appear among the reversers in the cleaned results, and only fedratinib is a consistent reverser.

---

## 4. Why we got these results

- **The signature is dominated by proliferation genes.** The most up-regulated genes include CDC45, CDCA7, NUSAP1, CCNA2, ASPM, NCAPG, MCM2/4/10, TOP2A, PLK1, AURKB, CDK1, CCNB1/2, FOXM1, MYBL2 and E2F2/7/8. These are DNA replication, mitosis and E2F-target genes.
- **Compounds that stop the cell cycle reverse such a signature, whatever the disease.** CDK inhibitors, topoisomerase poisons, MDM2/p53 activators, HDAC inhibitors and PARP inhibitors all switch these genes off, so they rank at the top. This explains why the leading hits are anti-proliferative agents.
- **The signature also carries immune biology.** CD38, CD8A/CD8B, LAG3, CCL5 and several interferon-stimulated genes (IFI27, ISG15, USP18, RTP4) are up, consistent with an activated, proliferating lymphocyte response and an interferon response in acute HIV infection. The many immunoglobulin variable-region genes (IGHV/IGKV/IGLV, 30 of the 150 up genes) are consistent with expansion of antibody-producing B cells.
- **The down-regulated side adds almost nothing.** Only 19 down genes could be used, and none of the L1000 signatures were enriched for them (FDR_down ≥ 0.92 throughout). The ranking therefore rests on the up-regulated genes. This is why scoring is based on the up set.

---

## 5. Interpretation

1. **The top reversers mainly show that the signature is a proliferation signature.** They are strong, reproducible reversers across many cell lines, but this pattern would be seen for almost any disease with a strong cell-cycle signal. They are not evidence of HIV-specific activity and are not recommended as repurposing candidates on this basis alone (many are cytotoxic oncology agents).
2. **Fedratinib is the most interesting hit for the client list.** It reverses consistently across 12 cell lines, it comes from a drug list already filtered for contraindications, and it agrees with the independent reversal value already in the client table. It also acts on JAK/STAT signalling, which is plausible in the context of interferon and cytokine signalling in acute infection. This is a hypothesis for follow-up, not a conclusion.
3. **Dexamethasone and lenalidomide do not reverse this signature here.** Their values in the client table were negative (suggesting reversal), but in this run they behave as mimickers. The difference may come from using a different signature, different gene-set size, or a different scoring method. It should be resolved before these two drugs are carried forward.
4. **Step 2a alone cannot nominate final drug candidates.** The metabolic branch (Step 2b) is needed to separate generic anti-proliferative effects from compounds acting on HIV-relevant metabolic enzymes and pathways.

---

## 6. Limitations

- **Cell-line data, not primary immune cells.** The L1000 signatures come mostly from cancer and other cell lines (the largest groups are MCF7, VCAP, A549, A375, MCF10A), not from PBMCs or T/B cells infected with HIV.
- **Truncated results.** The API returns only the top 5,000 hits in each direction, so per-compound signature counts are lower bounds. Use "reversing vs mimicking" as a consistency check rather than exact percentages.
- **Immunoglobulin genes.** 30 of the 150 up genes are IG variable segments. They resolved in the SigCom gene dictionary, but their measurement in L1000 is unreliable, and they may add noise. A sensitivity run without them is recommended.
- **Short down-regulated set** (19 genes), so the search is effectively one-sided.
- **Thresholds are set by the analyst** (z_up ≤ -5, FDR < 0.05, top 150 up genes). Results should be checked for stability to these choices.
- **Computational prediction only.** Reversal of a transcriptomic signature in vitro does not demonstrate clinical efficacy or safety in HIV.

---

## 7. Files

| File | Description |
|---|---|
| `step2a_signature_reversal.py` | Script for Step 2a (v2) |
| `results_v2/signature_up_symbols.txt` | 150 up genes used |
| `results_v2/signature_down_symbols.txt` | 19 down genes used |
| `results_v2/sigcom_raw_signature_hits.csv` | All 10,000 signature hits with z-scores, FDR, cell line, dose, time |
| `results_v2/sigcom_compound_summary.csv` | Per-compound summary (592 compounds), ranked |
| `Supplementary_-_1.xlsx` | Input DEG table (GSE199911) |
| `NO_CI_DRUGS.xlsx` | Client drug list used for the cross-check |

## 8. How to reproduce

```bash
pip install pandas openpyxl requests
python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --outdir results_v2
```

Optional arguments: `--zcut` (default 5), `--fdr` (default 0.05), `--n-up` (default 150), `--n-down`, `--database` (default `l1000_cp`), `--require-down`, `--list-databases`. The SigCom LINCS API is queried live, so results may change slightly if the database is updated.

## 9. Next steps

1. Sensitivity run excluding immunoglobulin genes; check that the top reversers and fedratinib stay stable.
2. Cross-check with an independent tool (L2S2 and/or iLINCS) using the same signature, and keep compounds supported by more than one tool.
3. Resolve the discordance for dexamethasone and lenalidomide by confirming which signature produced the existing values in `NO_CI_DRUGS.xlsx`.
4. Step 2b: map DEGs to Human1 enzymes and run GIMME/iMAT flux modelling (needs the full expression matrix, not only the 379 DEGs).
5. Step 3: integrate reversing compounds with target enzymes/pathways, network proximity and contraindication data to rank final candidates.