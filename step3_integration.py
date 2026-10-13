#!/usr/bin/env python3
"""
Step 3 (part 1): integrate Step 2a (signature reversal) and Step 2b (Human1 metabolic evidence)
for the client's candidate drugs (NO_CI_DRUGS.xlsx). Local files only, no downloads.

Run:  python step3_integration.py
      (expects the files below in the current folder, or pass paths)
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs", default="NO_CI_DRUGS.xlsx")
    ap.add_argument("--sigcom-raw", default="results_v2/sigcom_raw_signature_hits.csv")
    ap.add_argument("--activity", default="results_2b/2_reaction_activity_HIV_vs_control.csv")
    ap.add_argument("--degs", default="results_2b/1_metabolic_DEGs.csv")
    ap.add_argument("--candidates", default="results_2b/4_candidate_target_enzymes.csv")
    ap.add_argument("--genes", default="results_2b/genes.tsv", help="Human-GEM genes.tsv")
    ap.add_argument("--outdir", default="results_3")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(exist_ok=True, parents=True)

    n = pd.read_excel(a.drugs)
    drugs = n.drop_duplicates("Name")[["Name", "drug_name", "MOA", "sig_count", "hiv_mean_nes_up", "hiv_fdr_up",
                                       "Proximity Z_score", "n_contraindications", "In HIV module"]].copy()
    drugs["MOA"] = drugs["MOA"].astype(str).str.split(",").str[0]
    targets = n.groupby("Name")["Target"].apply(list)

    # ---- Step 2a: SigCom LINCS evidence
    r = pd.read_csv(a.sigcom_raw)
    r["compound"] = r["compound"].astype(str).str.lower()
    agg = r.groupby("compound").apply(lambda g: pd.Series({
        "sigcom_n_signatures": len(g), "sigcom_n_reversing": int(g["reverser"].sum()),
        "sigcom_n_mimicking": int(g["mimicker"].sum()),
        "sigcom_cell_lines_reversing": g.loc[g["reverser"], "cell_line"].astype(str).nunique(),
        "sigcom_median_z_up_reversers": g.loc[g["reverser"], "z_up"].median()}), include_groups=False).reset_index()
    agg = agg.set_index("compound")
    rows = []
    for _, d in drugs.iterrows():
        hit = None
        for key in (str(d["Name"]).lower(), str(d["drug_name"]).lower()):
            if key in agg.index:
                hit = agg.loc[key]; break
        rows.append(hit if hit is not None else pd.Series(np.nan, index=agg.columns))
    drugs = pd.concat([drugs.reset_index(drop=True), pd.DataFrame(rows).reset_index(drop=True)], axis=1)

    # ---- Step 2b: metabolic evidence for each drug's Human1 targets
    h1 = set(pd.read_csv(a.genes, sep="\t")["geneSymbols"])
    act = pd.read_csv(a.activity)
    gdir = {}
    for _, x in act.iterrows():
        for g in str(x["genes"]).split(";"):
            gdir.setdefault(g, []).append(x["direction"])
    deg = set(pd.read_csv(a.degs)["symbol"]); cand = pd.read_csv(a.candidates).set_index("gene")

    def metab(name):
        ts = [t for t in targets[name] if t in h1]
        sig = [t for t in ts if t in gdir and any(d != "ns" for d in gdir[t])]
        return pd.Series({
            "n_targets_total": len(targets[name]), "n_targets_in_Human1": len(ts),
            "human1_targets_with_activity_change": ";".join(sig),
            "human1_targets_that_are_DEGs": ";".join(t for t in ts if t in deg),
            "human1_targets_in_candidate_list": ";".join(f"{t}({cand.loc[t,'direction_in_HIV']},{cand.loc[t,'n_agreeing_sources']})"
                                                         for t in ts if t in cand.index)})
    drugs = pd.concat([drugs, drugs["Name"].apply(metab)], axis=1)

    # ---- evidence flags (simple, transparent; not a score)
    drugs["flag_reversal_SigCom"] = (drugs.sigcom_n_reversing.fillna(0) >= 3) & \
        (drugs.sigcom_n_reversing.fillna(0) > 2 * drugs.sigcom_n_mimicking.fillna(0))
    drugs["flag_client_signature"] = (drugs.hiv_mean_nes_up < 0) & (drugs.hiv_fdr_up < 0.05)
    drugs["flag_network_proximity"] = drugs["Proximity Z_score"] <= -1.96
    drugs["flag_no_contraindications"] = drugs["n_contraindications"] == 0
    # metabolic support = a Human1 target backed by >=2 agreeing evidence sources in Step 2b
    drugs["flag_metabolic_target_2b"] = drugs["human1_targets_in_candidate_list"].astype(str).apply(
        lambda s: any(int(x.split(",")[1].rstrip(")")) >= 2 for x in s.split(";") if "," in x))
    fl = [c for c in drugs.columns if c.startswith("flag_")]
    drugs["n_flags"] = drugs[fl].sum(axis=1)
    drugs = drugs.sort_values(["n_flags", "sigcom_n_reversing"], ascending=False)
    drugs.to_csv(out / "5_drug_integration_client_list.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(drugs[["Name", "hiv_mean_nes_up", "Proximity Z_score", "n_contraindications", "sigcom_n_reversing",
                 "sigcom_n_mimicking", "sigcom_cell_lines_reversing", "flag_metabolic_target_2b", "n_flags"]].to_string(index=False))
    print("Written:", (out / "5_drug_integration_client_list.csv").resolve())


if __name__ == "__main__":
    main()
