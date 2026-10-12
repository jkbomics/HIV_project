#!/usr/bin/env python3
"""
Step 2a: transcriptomic signature -> L1000 signature reversion (SigCom LINCS API)

Input : DEG table from GSE199911 (Supplementary_-_1.xlsx, sheet "HIV DEGs")
Output: results/
          signature_up_symbols.txt, signature_down_symbols.txt
          sigcom_raw_signature_hits.csv     one row per matching L1000 signature
          sigcom_compound_summary.csv       one row per compound (reversers ranked first)

Setup (Codespace or local):
    pip install pandas openpyxl requests
Run:
    python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx
If the database name is rejected, list the valid ones with:
    python step2a_signature_reversal.py --list-databases
"""
import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# If these URLs stop working, the same APIs are also served under https://ldp3.cloud/
DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"
META_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
GENE_INFO_URL = ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
                 "Homo_sapiens.gene_info.gz")


# ----------------------------------------------------------------- helpers
def post(url, payload, tries=3):
    for i in range(tries):
        try:
            r = requests.post(url, json=payload, timeout=180)
            r.raise_for_status()
            return r.json()
        except Exception as e:  # noqa: BLE001
            if i == tries - 1:
                sys.exit(f"Request failed: {url}\n{e}")
            time.sleep(3 * (i + 1))


def list_databases():
    r = requests.get(DATA_API + "/listdata", timeout=60)
    r.raise_for_status()
    print(json.dumps(r.json(), indent=2))


# ----------------------------------------------------------------- step 1: load DEGs
def load_degs(path, sheet):
    d = pd.read_excel(path, sheet_name=sheet) if str(path).endswith(".xlsx") \
        else pd.read_csv(path)
    need = {"GeneID", "log2FoldChange", "stat"}
    if not need.issubset(d.columns):
        sys.exit(f"DEG table needs columns {need}; found {list(d.columns)}")
    d = d.dropna(subset=["GeneID"]).copy()
    d["GeneID"] = d["GeneID"].astype(int)
    return d


def entrez_to_symbol(ids, outdir):
    cache = outdir / "Homo_sapiens.gene_info.gz"
    if not cache.exists():
        print("Downloading NCBI gene_info for Entrez -> symbol mapping ...")
        r = requests.get(GENE_INFO_URL, timeout=300)
        r.raise_for_status()
        cache.write_bytes(r.content)
    gi = pd.read_csv(cache, sep="\t", compression="gzip",
                     usecols=["GeneID", "Symbol"], dtype={"GeneID": int})
    m = gi.drop_duplicates("GeneID").set_index("GeneID")["Symbol"]
    return pd.Series(list(ids)).map(m).values


def build_signature(d, n_up, n_down, outdir):
    d = d.copy()
    d["symbol"] = entrez_to_symbol(d["GeneID"], outdir)
    d = d.dropna(subset=["symbol"])
    # keep only genes that exist in the L1000 gene space (drops LOC*/LINC*/IG* etc.),
    # THEN take the top N, so the signature is not filled with unmeasurable genes
    gmap = resolve_genes(list(d["symbol"]))
    d = d[d["symbol"].isin(gmap)]
    up = d[d["log2FoldChange"] > 0].sort_values("stat", ascending=False).head(n_up)
    dn = d[d["log2FoldChange"] < 0].sort_values("stat", ascending=True).head(n_down)
    (outdir / "signature_up_symbols.txt").write_text("\n".join(up["symbol"]) + "\n")
    (outdir / "signature_down_symbols.txt").write_text("\n".join(dn["symbol"]) + "\n")
    print(f"Signature (L1000-mappable genes only): {len(up)} up, {len(dn)} down genes")
    if len(dn) < 10:
        print("WARNING: <10 down genes usable -> down side is weak; use up-only scoring.")
    return list(up["symbol"]), list(dn["symbol"]), gmap


# ----------------------------------------------------------------- step 2: SigCom LINCS
def resolve_genes(symbols):
    payload = {"filter": {"where": {"meta.symbol": {"inq": symbols}},
                          "fields": ["id", "meta.symbol"]}}
    res = post(META_API + "/entities/find", payload)
    return {e["meta"]["symbol"]: e["id"] for e in res}


def search(up_ids, down_ids, database, limit):
    q = {"up_entities": up_ids, "down_entities": down_ids,
         "limit": limit, "database": database}
    return post(DATA_API + "/enrich/ranktwosided", q)["results"]


def fetch_signature_meta(uuids, chunk=200):
    out = {}
    for i in range(0, len(uuids), chunk):
        part = uuids[i:i + chunk]
        payload = {"filter": {"where": {"id": {"inq": part}}}}
        for s in post(META_API + "/signatures/find", payload):
            out[s["id"]] = s.get("meta", {})
    return out


def collect(up, down, gmap, database, limit):
    up_ids = [gmap[g] for g in up if g in gmap]
    dn_ids = [gmap[g] for g in down if g in gmap]
    missing = [g for g in up + down if g not in gmap]
    print(f"Genes found in L1000 gene space: up {len(up_ids)}/{len(up)}, "
          f"down {len(dn_ids)}/{len(down)}")
    if missing:
        print("  not found:", ", ".join(missing[:30]), "..." if len(missing) > 30 else "")

    hits = search(up_ids, dn_ids, database, limit)
    meta = fetch_signature_meta([h["uuid"] for h in hits])
    rows = []
    for h in hits:
        m = meta.get(h["uuid"], {})
        rows.append({
            "uuid": h["uuid"],
            "compound": m.get("pert_name") or m.get("pert_iname"),
            "cell_line": m.get("cell_line"),
            "dose": m.get("pert_dose") or m.get("pert_idose"),
            "time": m.get("pert_time") or m.get("pert_itime"),
            "z_up": h.get("z-up"),
            # flip z-down so reversers have negative z_up AND negative z_down
            "z_down": -h["z-down"] if h.get("z-down") is not None else None,
            "fdr_up": h.get("fdr-up"),
            "fdr_down": h.get("fdr-down"),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- step 3: summarise
def summarise(df, fdr, min_sigs, zcut, require_down=False):
    df = df.dropna(subset=["compound", "z_up"]).copy()
    # some compound / cell-line names come back as numbers from the API -> force text,
    # otherwise groupby fails with "'<' not supported between 'str' and 'int'"
    df["compound"] = df["compound"].astype(str)
    df["cell_line"] = df["cell_line"].astype(str)
    # Reverser = up-regulated disease genes sit at the BOTTOM of the L1000 signature.
    df["reverser"] = (df["z_up"] <= -zcut) & (df["fdr_up"] < fdr)
    df["mimicker"] = (df["z_up"] >= zcut) & (df["fdr_up"] < fdr)
    if require_down:
        df["reverser"] &= df["z_down"] < 0
    g = df.groupby("compound")
    s = pd.DataFrame({
        "n_signatures": g.size(),
        "n_reversing": g["reverser"].sum(),
        "n_mimicking": g["mimicker"].sum(),
        "median_z_up_reversers": df[df["reverser"]].groupby("compound")["z_up"].median(),
        "best_z_up": g["z_up"].min(),
    })
    s["n_cell_lines_reversing"] = (df[df["reverser"]].groupby("compound")["cell_line"]
                                   .nunique().reindex(s.index).fillna(0).astype(int))
    s["n_reversing"] = s["n_reversing"].astype(int)
    s["n_mimicking"] = s["n_mimicking"].astype(int)
    # NOTE: the API returns only the top hits per direction, so n_signatures is truncated;
    # use n_reversing vs n_mimicking as a consistency check, not as exact percentages.
    s = s[(s["n_signatures"] >= min_sigs) & (s["n_reversing"] > 0)]
    s = s.sort_values(["n_cell_lines_reversing", "median_z_up_reversers"],
                      ascending=[False, True])
    return df, s.reset_index()


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--degs", default="Supplementary_-_1.xlsx")
    ap.add_argument("--sheet", default="HIV DEGs")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--n-up", type=int, default=150)
    ap.add_argument("--n-down", type=int, default=150)   # you only have 35 down genes
    ap.add_argument("--database", default="l1000_cp",
                    help="chemical perturbation library; check with --list-databases")
    ap.add_argument("--limit", type=int, default=5000, help="top signatures to return")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-sigs", type=int, default=3,
                    help="min signatures per compound to keep in summary")
    ap.add_argument("--zcut", type=float, default=5.0,
                    help="|z_up| needed to call a signature a reverser/mimicker")
    ap.add_argument("--require-down", action="store_true",
                    help="also require z_down<0 (only sensible with >=10 usable down genes)")
    ap.add_argument("--list-databases", action="store_true")
    a = ap.parse_args()

    if a.list_databases:
        return list_databases()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    d = load_degs(a.degs, a.sheet)
    up, down, gmap = build_signature(d, a.n_up, a.n_down, out)
    raw = collect(up, down, gmap, a.database, a.limit)
    raw, summ = summarise(raw, a.fdr, a.min_sigs, a.zcut, a.require_down)

    raw.to_csv(out / "sigcom_raw_signature_hits.csv", index=False)
    summ.to_csv(out / "sigcom_compound_summary.csv", index=False)
    print(f"\n{len(raw)} signature hits, {len(summ)} compounds summarised")
    print(summ.head(20).to_string(index=False))
    print(f"\nFiles written to {out.resolve()}")


if __name__ == "__main__":
    main()
