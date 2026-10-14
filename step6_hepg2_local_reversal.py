#!/usr/bin/env python3
"""
Step 6: local signature reversal with the client's LINCS gene-set files (500 up / 500 down genes per compound).
Priority rule from the client: use HepG2; for a drug with no HepG2 signature, use HA1E for that drug only.

For every compound signature, the HIV up-regulated genes are tested for over-representation among the compound's
DOWN genes (reversal) and among its UP genes (mimicry), with a hypergeometric test against the 12,328-gene L1000 space
(the union of genes in the two files). The down-regulated HIV genes are reported but not used for calling, because
only about 13 of them are in the L1000 gene space.

Works offline (no API), so it covers drugs that were missing from the SigCom top hits.

Run:
    python step6_hepg2_local_reversal.py
    python step6_hepg2_local_reversal.py --exclude-proliferation --outdir results_6_noProlif   (sensitivity)

Outputs (results_6/):
    9_local_reversal_all_signatures.csv     every signature (all cell lines) with overlaps, fold enrichment and FDR
    9_local_focus_candidates.csv            candidates from Steps 2-4 with the cell line used and SigCom status alongside
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

# Core proliferation genes (Tirosh et al. 2016 S-phase and G2/M sets), same lists as Step 2a
S_PHASE = ("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 CENPU HELLS RFC2 RPA2 NASP "
           "RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM "
           "CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8").split()
G2M = ("HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 PIMREG SMC4 CCNB2 "
       "CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 CDC20 TTK CDC25C KIF2C RANGAP1 "
       "NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 "
       "CENPA").split()
SALTS = r"(hydrochloride|hcl|sodium|potassium|mesylate|acetate|phosphate|sulfate|sulphate|citrate|maleate|" \
        r"tartrate|besylate|tosylate|bromide|chloride|fumarate|succinate|calcium|disodium|dihydrate|hydrate)"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", re.sub(SALTS, " ", str(s).lower()))


def read_gmt(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        d[f[0][:-3]] = (f[1], set(x for x in f[2:] if x))      # strip the _UP / _DN suffix
    return d


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = p[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(n); q[o] = np.minimum(r, 1)
    return q


def score_cell_line(name, up_path, dn_path, up_genes, dn_genes, min_fold, fdr):
    """Score every signature of one cell line against the HIV signature."""
    U, D = read_gmt(up_path), read_gmt(dn_path)
    sigs = [s for s in U if s in D]
    univ = set().union(*[U[s][1] for s in sigs], *[D[s][1] for s in sigs]); N = len(univ)
    up, dn = up_genes & univ, dn_genes & univ
    n_set = len(next(iter(U.values()))[1])
    print(f"[{name}] {len(sigs)} signatures; gene space {N}; HIV signature in gene space: {len(up)} of "
          f"{len(up_genes)} up, {len(dn)} of {len(dn_genes)} down (expected overlap by chance {len(up) * n_set / N:.1f})")
    rows = [{"cell_line": name, "signature": s, "drug": U[s][0], "up_in_drugDN": len(up & D[s][1]),
             "up_in_drugUP": len(up & U[s][1]), "down_in_drugUP": len(dn & U[s][1]),
             "down_in_drugDN": len(dn & D[s][1])} for s in sigs]
    R = pd.DataFrame(rows)
    exp = len(up) * n_set / N
    R["fold_reversal"] = R["up_in_drugDN"] / exp
    R["fold_mimic"] = R["up_in_drugUP"] / exp
    R["p_reversal"] = hypergeom.sf(R["up_in_drugDN"] - 1, N, len(up), n_set)
    R["p_mimic"] = hypergeom.sf(R["up_in_drugUP"] - 1, N, len(up), n_set)
    R["fdr_reversal"] = bh(R["p_reversal"]); R["fdr_mimic"] = bh(R["p_mimic"])
    rev = (R["fdr_reversal"] < fdr) & (R["fold_reversal"] >= min_fold)
    mim = (R["fdr_mimic"] < fdr) & (R["fold_mimic"] >= min_fold)
    R["call"] = np.where(rev & ~mim, "reverser", np.where(mim & ~rev, "mimicker", "neutral"))
    R["key"] = R["drug"].map(norm)
    print(f"[{name}] calls: {R['call'].value_counts().to_dict()}")
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cell", action="append",
                    help="NAME:up.gmt:dn.gmt, repeat for each cell line IN PRIORITY ORDER. Default: HEPG2 first, "
                         "then HA1E (used only for drugs with no HepG2 signature), as instructed")
    ap.add_argument("--sig-up", default="results_v2/signature_up_symbols.txt")
    ap.add_argument("--sig-down", default="results_v2/signature_down_symbols.txt")
    ap.add_argument("--ranking", default="results_4/7_final_candidate_ranking.csv")
    ap.add_argument("--shortlist", default="results_4/7_repurposing_shortlist.csv")
    ap.add_argument("--client", default="results_3/5_drug_integration_list.csv")
    ap.add_argument("--metab", default="results_3/6_drugs_targeting_HIV_metabolic_genes.csv")
    ap.add_argument("--sigcom", default="results_v2/sigcom_compound_summary.csv")
    ap.add_argument("--exclude-proliferation", action="store_true",
                    help="sensitivity: drop core cell-cycle genes from the HIV signature before testing")
    ap.add_argument("--min-fold", type=float, default=2.0)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--outdir", default="results_6")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    cells = []
    if a.cell:
        for c in a.cell:
            n, u, d = c.split(":")
            cells.append((n, u, d))
    else:
        cells = [("HEPG2", "HEPG2_gene_sets_up_500.gmt", "HEPG2_gene_sets_dn_500.gmt"),
                 ("HA1E", "HA1E_gene_sets_up_500.gmt", "HA1E_gene_sets_dn_500.gmt")]
    have = []
    for n, u, d in cells:
        if Path(u).exists() and Path(d).exists():
            have.append((n, u, d))
        else:
            print(f"NOTE: files for {n} not found ({u}, {d}); skipping this cell line")
    if not have:
        raise SystemExit("No signature files found")

    up_all = set(open(a.sig_up).read().split()); dn_all = set(open(a.sig_down).read().split())
    if a.exclude_proliferation:
        up_all -= set(S_PHASE) | set(G2M)
        print("Sensitivity: core proliferation genes removed from the HIV signature")

    allR = [score_cell_line(n, u, d, up_all, dn_all, a.min_fold, a.fdr) for n, u, d in have]
    R = pd.concat(allR, ignore_index=True)
    order = [n for n, _, _ in have]
    R.sort_values(["cell_line", "p_reversal"]).drop(columns="key").to_csv(out / "9_local_reversal_all_signatures.csv", index=False)

    # ---- focus candidates: HepG2 first, next cell line only if the drug has no signature in the earlier ones
    S = pd.read_csv(a.sigcom); S["key"] = S["compound"].map(norm); S = S.drop_duplicates("key").set_index("key")
    S["sigcom_reverser"] = (S["n_reversing"] >= 3) & (S["n_reversing"] > 2 * S["n_mimicking"])
    foc = []
    rk = pd.read_csv(a.ranking); rk["key"] = rk["drug"].map(norm)
    for _, r in rk.iterrows():
        foc.append((r["key"], r["drug"], "ranking: " + str(r["tier"])))
    for d in pd.read_csv(a.shortlist)["drug"]:
        foc.append((norm(d), d, "repurposing shortlist"))
    for d in pd.read_csv(a.client)["Name"]:
        foc.append((norm(d), d, "client list"))
    me = pd.read_csv(a.metab); me = me[me["n_consistent_genes"] > 0]
    for d in me["drug"]:
        foc.append((norm(d), d, "direct inhibitor of a Step 2b enzyme"))
    F = pd.DataFrame(foc, columns=["key", "name", "group"])
    F = F.groupby("key").agg(name=("name", "first"), groups=("group", lambda x: "; ".join(sorted(set(x))))).reset_index()

    def best_per_key(cell_df):
        m = F[["key", "name"]].merge(cell_df, on="key")
        m["_exact"] = (m["name"].str.lower() == m["drug"].str.lower()).astype(int)   # prefer the exact name
        return m.sort_values(["_exact", "p_reversal"], ascending=[False, True]).drop_duplicates("key") \
                .drop(columns=["name", "_exact", "drug", "signature"]).set_index("key")
    per_cell = {n: best_per_key(R[R["cell_line"] == n]) for n in order}
    chosen = []
    for k in F["key"]:
        for n in order:
            if k in per_cell[n].index:
                row = per_cell[n].loc[k].to_dict(); row["key"] = k; chosen.append(row); break
    C = pd.DataFrame(chosen)
    F = F.merge(C, on="key", how="left")
    for n in order:                                   # informational: call in every cell line that has the drug
        F[f"call_in_{n}"] = F["key"].map(per_cell[n]["call"]) if len(per_cell[n]) else None
    F = F.join(S[["n_reversing", "n_mimicking", "n_cell_lines_reversing", "sigcom_reverser"]], on="key")
    F = F.rename(columns={"n_reversing": "sigcom_n_reversing", "n_mimicking": "sigcom_n_mimicking",
                          "n_cell_lines_reversing": "sigcom_cell_lines_reversing", "cell_line": "cell_line_used"})
    F["covered"] = F["call"].notna()
    cov = F["covered"]
    F["agreement"] = np.where(~cov, "no signature in any cell line",
                     np.where(F["sigcom_reverser"].isna(), "not in SigCom summary",
                     np.where(F["sigcom_reverser"] & (F["call"] == "reverser"), "both reverse",
                     np.where(F["sigcom_reverser"] & (F["call"] == "mimicker"), "CONFLICT (SigCom reverses, local mimics)",
                     np.where(F["sigcom_reverser"], "SigCom reverses, local neutral",
                     np.where(F["call"] == "reverser", "local reverses, SigCom not flagged", "neither"))))))
    cols = ["name", "groups", "cell_line_used", "call", "up_in_drugDN", "up_in_drugUP", "fold_reversal", "fold_mimic",
            "fdr_reversal", "fdr_mimic"] + [f"call_in_{n}" for n in order] + \
           ["sigcom_n_reversing", "sigcom_n_mimicking", "sigcom_cell_lines_reversing", "agreement"]
    F = F.sort_values(["covered", "call", "fold_reversal"], ascending=[False, False, False])
    F[cols].to_csv(out / "9_local_focus_candidates.csv", index=False)

    pd.set_option("display.width", 230); pd.set_option("display.max_columns", 20)
    print(f"\nFocus candidates: {len(F)}; with a signature in at least one cell line: {int(cov.sum())}")
    print("Cell line used:", F.loc[cov, "cell_line_used"].value_counts().to_dict())
    print(F.loc[cov, "call"].value_counts().to_string())
    print("\nAgreement with SigCom:")
    print(F["agreement"].value_counts().to_string())
    cl = F[F["groups"].str.contains("client list")]
    print(f"\nDrug-list ({int(cl['covered'].sum())} with a signature):")
    print(cl[["name", "cell_line_used", "call", "up_in_drugDN", "up_in_drugUP", "fold_reversal"]].round(2).to_string(index=False))
    print("\nFiles written to", out.resolve())


if __name__ == "__main__":
    main()
