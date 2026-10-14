#!/usr/bin/env python3
"""
Step 3 (part 2): link the HIV metabolic target genes (Step 2b) to drugs, and to the
signature reversers (Step 2a), using the DGIdb drug-gene interaction table.

Outputs (results_3/):
  6_drugs_targeting_HIV_metabolic_genes.csv   drugs hitting Step 2b genes, with SigCom evidence
  6_sigcom_reversers_with_targets.csv         SigCom reversers annotated with their DGIdb targets

Run (needs internet the first time, to download DGIdb):
    python step3b_drug_target_link.py
or with a file you already downloaded:
    python step3b_drug_target_link.py --dgidb interactions.tsv
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

DGIDB_URL = "https://www.dgidb.org/data/latest/interactions.tsv"
INHIBIT = {"inhibitor", "antagonist", "blocker", "antibody", "negative modulator", "suppressor",
           "inverse agonist", "antisense oligonucleotide", "cleavage", "inhibitory allosteric modulator",
           "antagonist/inhibitor", "inhibits"}
ACTIVATE = {"agonist", "activator", "positive modulator", "inducer", "stimulator", "potentiator",
            "partial agonist", "positive allosteric modulator", "enhancer", "activates"}
SALTS = r"(hydrochloride|hcl|sodium|potassium|mesylate|acetate|phosphate|sulfate|sulphate|citrate|maleate|" \
        r"tartrate|besylate|tosylate|bromide|chloride|fumarate|succinate|calcium|disodium|dihydrate|hydrate)"


def norm(s):
    s = str(s).lower()
    s = re.sub(SALTS, " ", s)
    return re.sub(r"[^a-z0-9]", "", s)


def get_dgidb(path):
    if path:
        return pd.read_csv(path, sep="\t", dtype=str)
    cache = Path("results_3/DGIdb_interactions.tsv")
    if not cache.exists():
        cache.parent.mkdir(exist_ok=True, parents=True)
        print("Downloading", DGIDB_URL)
        r = requests.get(DGIDB_URL, timeout=600)
        r.raise_for_status()
        cache.write_bytes(r.content)
    return pd.read_csv(cache, sep="\t", dtype=str)


def target_genes(cand, deg, act, sub):
    """Genes from Step 2b with the direction they move in HIV, and an evidence tier."""
    T = {}
    for _, r in cand.iterrows():
        if r["n_agreeing_sources"] >= 2 and r["direction_in_HIV"] in ("UP", "DOWN"):
            T[r["gene"]] = (r["direction_in_HIV"], "A: >=2 agreeing sources")
    for g in deg["symbol"].dropna():
        T.setdefault(g, ("UP", "B: differentially expressed enzyme"))
    for _, s in sub[sub["fdr"] < 0.05].iterrows():
        rx = act[(act["subsystem"] == s["subsystem"]) & (act["direction"] == s["direction"])]
        for gs in rx["genes"]:
            gl = [g for g in str(gs).split(";") if g]
            if 0 < len(gl) <= 3:
                for g in gl:
                    T.setdefault(g, (s["direction"], f"B: in significant pathway ({s['subsystem']})"))
    return T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dgidb", help="DGIdb interactions.tsv (downloaded automatically if omitted)")
    ap.add_argument("--sigcom", default="results_v2/sigcom_compound_summary.csv")
    ap.add_argument("--candidates", default="results_2b/4_candidate_target_enzymes.csv")
    ap.add_argument("--degs", default="results_2b/1_metabolic_DEGs.csv")
    ap.add_argument("--activity", default="results_2b/2_reaction_activity_HIV_vs_control.csv")
    ap.add_argument("--subsys", default="results_2b/2_subsystem_activity_enrichment.csv")
    ap.add_argument("--outdir", default="results_3")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(exist_ok=True, parents=True)

    cand = pd.read_csv(a.candidates); deg = pd.read_csv(a.degs)
    act = pd.read_csv(a.activity); sub = pd.read_csv(a.subsys)
    T = target_genes(cand, deg, act, sub)
    print(f"Step 2b target genes: {len(T)} "
          f"({sum(v[0]=='UP' for v in T.values())} up, {sum(v[0]=='DOWN' for v in T.values())} down in HIV)")

    D = get_dgidb(a.dgidb)
    D.columns = [c.strip().lower() for c in D.columns]
    need = {"gene_name", "drug_name"}
    if not need.issubset(D.columns):
        raise SystemExit(f"DGIdb file lacks columns {need}; found {list(D.columns)}")
    D["interaction_type"] = D.get("interaction_type", pd.Series("", index=D.index)).fillna("").str.lower()
    # the real DGIdb table has rows with a missing gene or drug name: drop them
    D = D.dropna(subset=["gene_name", "drug_name"]).copy()
    D["gene_name"] = D["gene_name"].astype(str).str.strip()
    D["drug_name"] = D["drug_name"].astype(str).str.strip()
    D = D[(D["drug_name"] != "") & (D["gene_name"] != "")].copy()
    D["drug_key"] = D["drug_name"].map(norm)
    print(f"DGIdb: {len(D)} interactions, {D['drug_key'].nunique()} drugs")

    S = pd.read_csv(a.sigcom).copy()
    S["drug_key"] = S["compound"].map(norm)
    S = S.drop_duplicates("drug_key").set_index("drug_key")
    S["reverser_flag"] = (S["n_reversing"] >= 3) & (S["n_reversing"] > 2 * S["n_mimicking"])

    # ---- A. drugs hitting Step 2b genes (direction-aware)
    H = D[D["gene_name"].isin(T)].copy()
    H["gene_direction_in_HIV"] = H["gene_name"].map(lambda g: T[g][0])
    H["gene_evidence"] = H["gene_name"].map(lambda g: T[g][1])

    def effect(r):
        t = r["interaction_type"]
        if t in INHIBIT or "inhibit" in t or "antagon" in t or "block" in t:
            return "opposes_gene"
        if t in ACTIVATE or "agonist" in t or "activat" in t:
            return "enhances_gene"
        return "unknown_effect"
    H["effect"] = H.apply(effect, axis=1)
    H["direction_consistent"] = ((H["gene_direction_in_HIV"] == "UP") & (H["effect"] == "opposes_gene")) | \
                                ((H["gene_direction_in_HIV"] == "DOWN") & (H["effect"] == "enhances_gene"))
    g = H.groupby("drug_key").agg(
        drug=("drug_name", lambda x: sorted(set(x))[0]),
        genes=("gene_name", lambda x: ";".join(sorted(set(x)))),
        n_genes=("gene_name", "nunique"),
        n_consistent_genes=("gene_name", lambda x: H.loc[x.index][H.loc[x.index, "direction_consistent"]]["gene_name"].nunique()),
        interaction_types=("interaction_type", lambda x: ";".join(sorted(set(t for t in x if t)))),
        tier_A_genes=("gene_name", lambda x: ";".join(sorted(set(
            n for n in x if T[n][1].startswith("A"))))),
        n_sources=("interaction_source_db_name", lambda x: x.nunique()) if "interaction_source_db_name" in H else ("gene_name", "size"),
    ).reset_index()
    g = g.join(S[["n_signatures", "n_reversing", "n_mimicking", "n_cell_lines_reversing",
                  "median_z_up_reversers", "reverser_flag"]], on="drug_key")
    g["in_SigCom_hits"] = g["n_signatures"].notna()
    g["reverser_flag"] = g["reverser_flag"].fillna(False).astype(bool)
    g["supported_both"] = g["reverser_flag"] & (g["n_consistent_genes"] > 0)
    g = g.sort_values(["supported_both", "n_consistent_genes", "n_reversing"], ascending=False)
    g.drop(columns="drug_key").to_csv(out / "6_drugs_targeting_HIV_metabolic_genes.csv", index=False)

    # ---- B. SigCom reversers and their (any) DGIdb targets
    R = S[S["reverser_flag"]].copy().reset_index()
    tg = D.groupby("drug_key")["gene_name"].apply(lambda x: sorted(set(x)))
    R["n_known_targets"] = R["drug_key"].map(tg).apply(lambda l: len(l) if isinstance(l, list) else 0)
    R["targets_in_step2b"] = R["drug_key"].map(tg).apply(
        lambda l: ";".join(x for x in l if x in T) if isinstance(l, list) else "")
    R["annotated_in_DGIdb"] = R["drug_key"].isin(tg.index)
    R = R.sort_values(["n_cell_lines_reversing", "median_z_up_reversers"], ascending=[False, True])
    R.drop(columns="drug_key").to_csv(out / "6_sigcom_reversers_with_targets.csv", index=False)

    print(f"\nDrugs in DGIdb hitting Step 2b genes: {len(g)}")
    print(f"  of which reverse the signature in SigCom (flag): {int(g['reverser_flag'].sum())}")
    print(f"  of which reverse AND act on the gene in the opposing direction: {int(g['supported_both'].sum())}")
    print(f"SigCom reversers (flag): {len(R)}; annotated in DGIdb: {int(R['annotated_in_DGIdb'].sum())}; "
          f"with a Step 2b target: {int((R['targets_in_step2b'] != '').sum())}")
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 20)
    cols = ["drug", "genes", "interaction_types", "n_consistent_genes", "n_reversing", "n_mimicking",
            "n_cell_lines_reversing", "supported_both"]
    print("\nTop drugs hitting Step 2b genes:")
    print(g[cols].head(25).to_string(index=False))
    print("\nFiles written to", out.resolve())


if __name__ == "__main__":
    main()