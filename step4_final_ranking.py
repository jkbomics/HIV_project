#!/usr/bin/env python3
"""
Step 4: final candidate ranking. Combines
  Step 2a  SigCom signature reversal              (results_v2/sigcom_compound_summary.csv)
  Step 2b/3b  direct drug action on HIV-changed metabolic enzymes   (results_3/6_drugs_targeting_...csv)
  Client list evidence                            (results_3/5_drug_integration_client_list.csv)
  DGIdb drug class flags: approved / anti-neoplastic / immunotherapy  (results_3/DGIdb_interactions.tsv)

Outputs (results_4/):
  7_final_candidate_ranking.csv      all candidates with evidence, class and tier
  7_repurposing_shortlist.csv        approved, non-anti-cancer candidates that reverse the signature

Run:  python step4_final_ranking.py
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

SALTS = r"(hydrochloride|hcl|sodium|potassium|mesylate|acetate|phosphate|sulfate|sulphate|citrate|maleate|" \
        r"tartrate|besylate|tosylate|bromide|chloride|fumarate|succinate|calcium|disodium|dihydrate|hydrate)"


def norm(s):
    s = re.sub(SALTS, " ", str(s).lower())
    return re.sub(r"[^a-z0-9]", "", s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigcom", default="results_v2/sigcom_compound_summary.csv")
    ap.add_argument("--metab", default="results_3/6_drugs_targeting_HIV_metabolic_genes.csv")
    ap.add_argument("--client", default="results_3/5_drug_integration_list.csv")
    ap.add_argument("--dgidb", default="results_3/DGIdb_interactions.tsv")
    ap.add_argument("--outdir", default="results_4")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(exist_ok=True, parents=True)

    # ---- SigCom
    S = pd.read_csv(a.sigcom); S["key"] = S["compound"].map(norm)
    S = S.drop_duplicates("key").set_index("key")
    S["reverser_flag"] = (S["n_reversing"] >= 3) & (S["n_reversing"] > 2 * S["n_mimicking"])

    # ---- DGIdb class flags (per drug)
    D = pd.read_csv(a.dgidb, sep="\t", dtype=str, usecols=["drug_name", "gene_name", "approved", "immunotherapy",
                                                           "anti_neoplastic"]).dropna(subset=["drug_name"])
    D["key"] = D["drug_name"].map(norm)
    cls = D.groupby("key").agg(approved=("approved", lambda x: (x == "TRUE").any()),
                               anti_neoplastic=("anti_neoplastic", lambda x: (x == "TRUE").any()),
                               immunotherapy=("immunotherapy", lambda x: (x == "TRUE").any()),
                               n_known_targets=("gene_name", "nunique"))

    # ---- metabolic target evidence (direction-consistent drugs)
    M = pd.read_csv(a.metab); M["key"] = M["drug"].map(norm)
    M = M[M["n_consistent_genes"] > 0].drop_duplicates("key").set_index("key")

    # ---- client list
    C = pd.read_csv(a.client); C["key"] = C["Name"].map(norm)
    C = C.drop_duplicates("key").set_index("key")

    # client-list drugs sometimes carry a lab code (e.g. TG-101348): also look them up by generic name
    for k in C.index:
        alt = norm(C.loc[k, "drug_name"]) if "drug_name" in C.columns else None
        if k not in cls.index and alt in cls.index:
            cls.loc[k] = cls.loc[alt]

    keys = set(S.index[S["reverser_flag"]]) | set(M.index) | set(C.index)
    rows = []
    for k in keys:
        r = {"key": k}
        if k in S.index:
            s = S.loc[k]
            r.update(drug=s["compound"], sigcom_n_reversing=s["n_reversing"], sigcom_n_mimicking=s["n_mimicking"],
                     sigcom_cell_lines_reversing=s["n_cell_lines_reversing"],
                     sigcom_median_z_up=s["median_z_up_reversers"], sigcom_reverser=bool(s["reverser_flag"]))
        if k in M.index:
            r.setdefault("drug", M.loc[k, "drug"])
            r.update(metabolic_target_genes=M.loc[k, "genes"], metabolic_target_tierA=M.loc[k, "tier_A_genes"],
                     metabolic_consistent=True)
        if k in C.index:
            c = C.loc[k]
            r.setdefault("drug", c["Name"])
            r.update(client_list=True, client_nes=c["hiv_mean_nes_up"], client_fdr=c["hiv_fdr_up"],
                     client_proximity_z=c["Proximity Z_score"], client_contraindications=c["n_contraindications"])
        if k in cls.index:
            r.update(approved=bool(cls.loc[k, "approved"]), anti_neoplastic=bool(cls.loc[k, "anti_neoplastic"]),
                     immunotherapy=bool(cls.loc[k, "immunotherapy"]), n_known_targets=int(cls.loc[k, "n_known_targets"]))
        rows.append(r)
    R = pd.DataFrame(rows)
    for c in ["sigcom_reverser", "metabolic_consistent", "client_list", "approved", "anti_neoplastic", "immunotherapy"]:
        R[c] = R[c].fillna(False).astype(bool) if c in R else False
    R["in_DGIdb"] = R["key"].isin(cls.index)

    # ---- transparent evidence points (not a statistical score)
    pts = pd.DataFrame(index=R.index)
    pts["reversal"] = R["sigcom_reverser"].astype(int) * 2 + (R["sigcom_cell_lines_reversing"].fillna(0) >= 10).astype(int)
    pts["metabolic_target"] = R["metabolic_consistent"].astype(int) * 2 + R["metabolic_target_tierA"].fillna("").ne("").astype(int)
    pts["client_signature"] = ((R["client_nes"] < 0) & (R["client_fdr"] < 0.05)).astype(int)
    pts["client_proximity"] = (R["client_proximity_z"] <= -1.96).astype(int)
    pts["client_safety"] = (R["client_contraindications"] == 0).astype(int)
    R["evidence_points"] = pts.sum(axis=1)
    R["evidence_detail"] = pts.apply(lambda r: ",".join(f"{k}:{v}" for k, v in r.items() if v), axis=1)

    def drug_class(r):
        if r["anti_neoplastic"]:
            return "anti-cancer (cytotoxic/targeted)"
        if r["immunotherapy"]:
            return "immunotherapy"
        if r["approved"]:
            return "approved, non-oncology"
        return "investigational / not annotated" if r["in_DGIdb"] else "not found in DGIdb"
    R["drug_class"] = R.apply(drug_class, axis=1)

    def tier(r):
        two = r["sigcom_reverser"] and (r["metabolic_consistent"] or (r["client_list"] and r["evidence_points"] >= 5))
        if two:
            return "1: reversal + independent support"
        if r["sigcom_reverser"] and r["sigcom_cell_lines_reversing"] >= 10:
            return "2: strong reversal only"
        if r["metabolic_consistent"] and not r["sigcom_reverser"]:
            return "3: enzyme-targeting, reversal not shown"
        if r["sigcom_reverser"]:
            return "2b: reversal (<10 cell lines)"
        return "4: client list, no reversal"
    R["tier"] = R.apply(tier, axis=1)
    R = R.sort_values(["tier", "evidence_points", "sigcom_cell_lines_reversing"], ascending=[True, False, False])
    cols = ["drug", "tier", "evidence_points", "evidence_detail", "drug_class", "approved", "anti_neoplastic",
            "sigcom_n_reversing", "sigcom_n_mimicking", "sigcom_cell_lines_reversing", "sigcom_median_z_up",
            "metabolic_target_genes", "metabolic_target_tierA", "client_list", "client_nes", "client_proximity_z",
            "client_contraindications", "n_known_targets"]
    R[[c for c in cols if c in R]].to_csv(out / "7_final_candidate_ranking.csv", index=False)

    short = R[R["sigcom_reverser"] & R["approved"] & ~R["anti_neoplastic"] & ~R["immunotherapy"]]
    short = short.sort_values(["sigcom_cell_lines_reversing", "sigcom_median_z_up"], ascending=[False, True])
    short[[c for c in cols if c in short]].to_csv(out / "7_repurposing_shortlist.csv", index=False)

    pd.set_option("display.width", 230); pd.set_option("display.max_columns", 20)
    print(f"Candidates: {len(R)}")
    print(R["tier"].value_counts().to_string())
    print("\nDrug class among SigCom reversers:")
    print(R[R["sigcom_reverser"]]["drug_class"].value_counts().to_string())
    show = ["drug", "tier", "evidence_points", "drug_class", "sigcom_n_reversing", "sigcom_n_mimicking",
            "sigcom_cell_lines_reversing", "metabolic_target_genes"]
    print("\nTier 1 and the 15 best others:")
    print(R[R["tier"].str.startswith("1")][show].to_string(index=False))
    print(f"\nApproved, non-anti-cancer reversers (shortlist): {len(short)}")
    print(short[["drug", "sigcom_n_reversing", "sigcom_n_mimicking", "sigcom_cell_lines_reversing",
                 "sigcom_median_z_up", "n_known_targets"]].head(25).to_string(index=False))
    print("\nFiles written to", out.resolve())


if __name__ == "__main__":
    main()
