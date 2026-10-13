#!/usr/bin/env python3
"""
Step 2b (complete): enzyme/gene-set mapping + Human1 metabolic flux simulation

Stages
  0  Download GSE199911 counts + sample annotation from GEO (or use local files)
  1  Map your DEGs to Human1 (Human-GEM) enzymes, reactions and subsystems
  2  Reaction-level activity per sample (from GPR rules) and HIV-vs-control test
  3  Context-specific flux models for HIV and control (GIMME, optional iMAT)
  4  Compare flux, rank target enzymes and pathways

Setup (Codespace or local; no COBRApy or solver install needed, only scipy/HiGHS):
    pip install numpy pandas scipy pyyaml requests openpyxl

Run (everything, GIMME):
    python step2b_pipeline.py --degs Supplementary_-_1.xlsx

Run with your own count matrix / group file instead of downloading:
    python step2b_pipeline.py --counts my_counts.csv --groups my_groups.csv

Add iMAT (slow MILP, has a time limit per condition):
    python step2b_pipeline.py --method both --imat-time-limit 900
"""
import argparse, gzip, io, pickle, re, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import scipy.sparse as sp
import yaml
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.stats import hypergeom, ttest_ind

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199911/"
COUNTS_URL = GEO + "suppl/GSE199911_processeddata_genecounts.csv.gz"
MATRIX_URL = GEO + "matrix/GSE199911_series_matrix.txt.gz"
GEM = "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/"

# Nutrients allowed to enter the model (everything else: no uptake, secretion allowed).
# The last three vitamins/cofactors are needed for the generic biomass reaction to be feasible.
MEDIUM = ['glucose', 'o2', 'h2o', 'h+', 'co2', 'pi', 'nh3', 'nh4+', 'sulfate', 'hco3-', 'na+', 'k+',
          'ca2+', 'mg2+', 'fe2+', 'fe3+', 'chloride', 'zinc', 'cu2+', 'iodide', 'alanine', 'arginine',
          'asparagine', 'aspartate', 'cysteine', 'glutamine', 'glutamate', 'glycine', 'histidine',
          'isoleucine', 'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine',
          'threonine', 'tryptophan', 'tyrosine', 'valine', 'linoleate', 'linolenate', 'choline',
          'inositol', 'thiamin', 'riboflavin', 'folate', 'pyridoxine', 'pantothenate', 'nicotinamide',
          'biotin', 'cholesterol', 'oleate', 'palmitate', 'retinol', 'alpha-tocopherol',
          'hypoxanthine', 'thymidine', 'pyruvate', 'l-lactate', 'acetate', 'taurine', 'cystine',
          'gamma-tocotrienol', 'lipoic acid', 'aquacob(iii)alamin']


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = p[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(n); q[o] = np.minimum(r, 1)
    return q


def download(url, path):
    if path.exists():
        return path
    log("Downloading", url)
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    path.write_bytes(r.content)
    return path


# =============================================================== Stage 0: data
def parse_series_matrix(path):
    """Return a DataFrame with one row per sample: gsm, title, text (all annotation)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    keys = {}
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                keys.setdefault(parts[0], []).append([p.strip('"') for p in parts[1:]])
    gsm = keys["!Sample_geo_accession"][0]
    title = keys["!Sample_title"][0]
    ann = [""] * len(gsm)
    for k, rows in keys.items():
        if k in ("!Sample_characteristics_ch1", "!Sample_source_name_ch1", "!Sample_title"):
            for row in rows:
                ann = [a + " | " + v for a, v in zip(ann, row)]
    return pd.DataFrame({"gsm": gsm, "title": title, "text": ann})


def infer_group(text):
    t = text.lower()
    ctrl = re.search(r"control|healthy|uninfected|negative|\bhc\b|\bnc\b|seroneg", t)
    hiv = re.search(r"hiv|infected|acute|positive|\bpos\b", t)
    if ctrl and not re.search(r"hiv[- ]?(positive|infected)", t):
        return "control"
    if hiv:
        return "HIV"
    return "unknown"


def load_counts_and_groups(args, out):
    cpath = Path(args.counts) if args.counts else download(COUNTS_URL, out / "GSE199911_genecounts.csv.gz")
    raw = pd.read_csv(cpath)
    log(f"Count file: {raw.shape[0]} rows x {raw.shape[1]} columns")
    log("  first columns:", [str(c) for c in raw.columns[:6]])
    samples_in_rows = (args.layout == "samples-rows") or (
        args.layout == "auto" and raw.shape[0] < raw.shape[1])
    sample_meta = None
    if samples_in_rows:
        sc = args.sample_col or ("sample" if "sample" in raw.columns else raw.columns[0])
        ens = [c for c in raw.columns if str(c).startswith("ENSG")]
        if len(ens) > 1000:
            gene_cols = ens
        else:
            skip = {sc, "unnamed: 0", "studygroup", "age", "sex"}
            gene_cols = [c for c in raw.select_dtypes(include=[np.number]).columns if str(c).lower() not in skip]
        meta_cols = [c for c in raw.columns if c not in gene_cols and c != sc]
        log(f"Layout: samples in rows (sample IDs from column '{sc}'); {len(gene_cols)} gene columns; "
            f"sample-info columns kept aside: {meta_cols}")
        raw[sc] = raw[sc].astype(str)
        sample_meta = raw.set_index(sc)[meta_cols]
        X = raw.set_index(sc)[gene_cols].T
        X.index = X.index.astype(str)
    else:
        log("Layout: genes in rows (gene IDs taken from the first column)")
        X = raw.set_index(raw.columns[0]).select_dtypes(include=[np.number])
        X.index = X.index.astype(str)
    log(f"Counts: {X.shape[0]} genes x {X.shape[1]} samples; first gene IDs {list(X.index[:3])}; "
        f"first sample IDs {list(X.columns[:3])}")

    meta = None
    if args.groups:
        g = pd.read_csv(args.groups)
        g.columns = [c.lower() for c in g.columns]
        grp = dict(zip(g["sample"].astype(str), g["group"]))
    else:
        mpath = download(MATRIX_URL, out / "GSE199911_series_matrix.txt.gz")
        meta = parse_series_matrix(mpath)
        meta["group"] = meta["text"].map(infer_group)
        meta.to_csv(out / "sample_annotation_from_GEO.csv", index=False)
        grp = {}
        for _, r in meta.iterrows():
            grp[r["gsm"]] = r["group"]; grp[r["title"]] = r["group"]
    G = pd.DataFrame({"sample": list(X.columns)})
    G["group"] = [grp.get(c) if grp.get(c) in ("HIV", "control") else "unknown" for c in G["sample"]]
    if sample_meta is not None and "studygroup" in sample_meta.columns:
        G["studygroup"] = sample_meta.loc[G["sample"], "studygroup"].values
        log("Cross-check, GEO-annotation group vs studygroup code in the count file:")
        log("\n" + pd.crosstab(G["group"], G["studygroup"]).to_string())
        if args.hiv_code is None and G["group"].nunique() < 3 and (G["group"] == "control").sum() == 0:
            log("GEO text did not separate HIV from control (everything looks like HIV). "
                "Re-run with --hiv-code <value of studygroup that means HIV>, e.g. --hiv-code 1")
        if args.hiv_code is not None:
            code = G["studygroup"].astype(str) == str(args.hiv_code)
            G["group"] = np.where(code, "HIV", "control")
            log(f"Groups set from studygroup: code {args.hiv_code} = HIV, all others = control")
    G.to_csv(out / "groups_used.csv", index=False)
    log("Group counts:", G["group"].value_counts().to_dict(),
        f"(written to {out/'groups_used.csv'} - CHECK IT)")
    if (G["group"] == "unknown").all() and meta is not None:
        log("No count-file sample ID matched a GEO sample title or GSM accession, or the GEO text had no")
        log("  recognisable HIV/control wording. First rows of the GEO annotation:")
        log("\n" + meta[["gsm", "title", "text"]].head(4).to_string())
        log("  Open results_2b/sample_annotation_from_GEO.csv, then either pass --hiv-code <studygroup value that")
        log("  means HIV> or build a groups CSV (columns sample,group) and pass --groups.")
    if (G["group"] == "HIV").sum() < 3 or (G["group"] == "control").sum() < 3:
        sys.exit("Could not assign enough samples to HIV/control (see messages above).")
    keep = G[G["group"].isin(["HIV", "control"])]
    return X[keep["sample"]], keep.set_index("sample")["group"]


def normalise(X):
    """DESeq2-style median-of-ratios size factors, then log2(x+1)."""
    X = X.loc[(X >= 10).sum(axis=1) >= 3]
    lg = np.log(X.where(X > 0))
    ref = lg.mean(axis=1)
    ok = ref.notna() & np.isfinite(ref) & (X > 0).all(axis=1)
    sf = np.exp((lg[ok].sub(ref[ok], axis=0)).median(axis=0))
    return np.log2(X.div(sf, axis=1) + 1)


# =============================================================== Human-GEM
def load_model(out):
    cache = out / "humangem_parsed.pkl"
    if cache.exists():
        return pickle.load(open(cache, "rb"))
    yml = download(GEM + "Human-GEM.yml", out / "Human-GEM.yml")
    genes = download(GEM + "genes.tsv", out / "genes.tsv")
    log("Parsing Human-GEM (about a minute)...")
    doc = yaml.load(open(yml), Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader))
    sec = {}
    for item in doc:
        k, v = item if isinstance(item, tuple) else list(item.items())[0]
        sec[k] = v
    meta = dict(sec["metaData"])
    mets = {dict(m)["id"]: str(dict(m).get("name")).lower() for m in sec["metabolites"]}
    mi = {m: i for i, m in enumerate(mets)}
    rows, cols, vals = [], [], []
    ids, gpr, ss, lb, ub, exch, bio = [], [], [], [], [], {}, None
    for j, r in enumerate(sec["reactions"]):
        r = dict(r)
        ids.append(r["id"]); gpr.append(r.get("gene_reaction_rule") or "")
        s = r.get("subsystem") or ""
        ss.append(";".join(s) if isinstance(s, list) else s)
        lb.append(r["lower_bound"]); ub.append(r["upper_bound"])
        m = dict(r["metabolites"])
        for k, c in m.items():
            rows.append(mi[k]); cols.append(j); vals.append(c)
        if len(m) == 1:
            exch[j] = mets[list(m)[0]]
        if r.get("objective_coefficient"):
            bio = j
    M = dict(S=sp.csr_matrix((vals, (rows, cols)), shape=(len(mets), len(ids))), ids=ids, gpr=gpr,
             subsystem=ss, lb=np.array(lb, float), ub=np.array(ub, float), exch=exch, bio=bio, meta=meta,
             genes=pd.read_csv(genes, sep="\t"))
    pickle.dump(M, open(cache, "wb"))
    return M


def medium_bounds(M):
    lb = M["lb"].copy()
    for j, nm in M["exch"].items():
        if nm not in MEDIUM:
            lb[j] = 0.0   # no uptake of anything outside the medium
    return lb, M["ub"].copy()


# =============================================================== GPR evaluation
def parse_gpr(s):
    toks = re.findall(r"\(|\)|\bor\b|\band\b|ENSG\d+", s)
    pos = [0]

    def factor():
        t = toks[pos[0]]; pos[0] += 1
        if t == "(":
            n = expr(); pos[0] += 1; return n
        return t

    def term():
        parts = [factor()]
        while pos[0] < len(toks) and toks[pos[0]] == "and":
            pos[0] += 1; parts.append(factor())
        return ("and", parts) if len(parts) > 1 else parts[0]

    def expr():
        parts = [term()]
        while pos[0] < len(toks) and toks[pos[0]] == "or":
            pos[0] += 1; parts.append(term())
        return ("or", parts) if len(parts) > 1 else parts[0]

    return expr() if toks else None


def eval_gpr(node, E):
    """AND -> min, OR -> max over the genes that have data; None if no gene has data."""
    if node is None:
        return None
    if isinstance(node, str):
        return E.get(node)
    vals = [v for v in (eval_gpr(c, E) for c in node[1]) if v is not None]
    if not vals:
        return None
    return np.minimum.reduce(vals) if node[0] == "and" else np.maximum.reduce(vals)


def reaction_matrix(M, E, nsamp):
    cache, R = {}, np.full((len(M["ids"]), nsamp), np.nan)
    for j, g in enumerate(M["gpr"]):
        if not g:
            continue
        if g not in cache:
            cache[g] = eval_gpr(parse_gpr(g), E)
        if cache[g] is not None:
            R[j] = cache[g]
    return R


def gene_expression_for_model(logx, M):
    """Map count-matrix row IDs (Entrez, Ensembl or symbol) to Human-GEM Ensembl IDs."""
    G = M["genes"].copy()
    G["geneEntrezID"] = pd.to_numeric(G["geneEntrezID"], errors="coerce")
    idx = pd.Series(logx.index.astype(str), index=logx.index)
    if idx.str.fullmatch(r"\d+").mean() > 0.8:
        m = G.dropna(subset=["geneEntrezID"]).assign(k=lambda d: d["geneEntrezID"].astype(int).astype(str))
        kind = "Entrez"
    elif idx.str.startswith("ENSG").mean() > 0.8:
        m = G.assign(k=G["genes"]); idx = idx.str.replace(r"\..*", "", regex=True); kind = "Ensembl"
    else:
        m = G.assign(k=G["geneSymbols"]); kind = "symbol"
    log(f"Count matrix gene IDs detected as {kind}")
    mp = dict(zip(m["k"], m["genes"]))
    tmp = logx.copy(); tmp["ensg"] = idx.map(mp).values
    tmp = tmp.dropna(subset=["ensg"]).groupby("ensg").max()   # one row per Human1 gene
    log(f"{len(tmp)} Human1 genes have expression data")
    return tmp


def sanity_check_groups(logx, groups):
    """Interferon-stimulated genes are strongly up in acute HIV; warn/stop if the labels look reversed."""
    isg = {"IFI27": "ENSG00000165949", "ISG15": "ENSG00000187608", "MX1": "ENSG00000157601",
           "USP18": "ENSG00000084110", "IFI44L": "ENSG00000137959"}
    ids = pd.Series(logx.index.astype(str).str.replace(r"\..*", "", regex=True), index=logx.index)
    isH = (groups.loc[logx.columns] == "HIV").values
    diffs = {}
    for name, e in isg.items():
        hit = ids[ids == e]
        if len(hit):
            row = logx.loc[hit.index[0]].values
            diffs[name] = row[isH].mean() - row[~isH].mean()
    if not diffs:
        log("Sanity check skipped (interferon-stimulated genes not found by Ensembl ID)")
        return
    log("Sanity check, log2 difference HIV minus control: " + ", ".join(f"{k} {v:+.2f}" for k, v in diffs.items()))
    if np.mean([v > 0 for v in diffs.values()]) < 0.5:
        sys.exit("STOP: interferon-stimulated genes are LOWER in the samples labelled HIV, so the labels look "
                 "reversed. Re-run with the other --hiv-code value.")


# =============================================================== Stage 1: DEG mapping
def stage1_deg_mapping(args, M, out):
    D = pd.read_excel(args.degs, sheet_name=args.sheet) if str(args.degs).endswith(".xlsx") else pd.read_csv(args.degs)
    D["GeneID"] = D["GeneID"].astype(int)
    G = M["genes"].copy()
    G["geneEntrezID"] = pd.to_numeric(G["geneEntrezID"], errors="coerce")
    e2g = G.dropna(subset=["geneEntrezID"]).drop_duplicates("geneEntrezID").set_index("geneEntrezID")
    D["ensg"] = D["GeneID"].map(e2g["genes"]); D["symbol"] = D["GeneID"].map(e2g["geneSymbols"])
    sym = dict(zip(G["genes"], G["geneSymbols"]))
    gl = [sorted(set(re.findall(r"ENSG\d+", g))) for g in M["gpr"]]
    allg = set(x for l in gl for x in l)
    met = D[D["ensg"].isin(allg)].copy()
    log(f"Stage 1: {len(D)} DEGs -> {len(met)} are Human1 enzymes/transporters "
        f"({(met.log2FoldChange > 0).sum()} up, {(met.log2FoldChange < 0).sum()} down)")
    gene_rx = {}
    for j, l in enumerate(gl):
        for x in l:
            gene_rx.setdefault(x, []).append(j)
    met["n_reactions"] = met["ensg"].map(lambda x: len(gene_rx[x]))
    met["subsystems"] = met["ensg"].map(lambda x: ";".join(sorted({M["subsystem"][j] for j in gene_rx[x]} - {""})))
    met.sort_values("padj").to_csv(out / "1_metabolic_DEGs.csv", index=False)
    N, K, rows = len(allg), len(set(met["ensg"])), []
    sub = {}
    for j, l in enumerate(gl):
        if l and M["subsystem"][j]:
            sub.setdefault(M["subsystem"][j], set()).update(l)
    for s, gs in sub.items():
        k = len(gs & set(met["ensg"]))
        if k:
            rows.append((s, len(gs), k, hypergeom.sf(k - 1, N, K, len(gs)),
                         ";".join(sorted(sym.get(x, x) for x in gs & set(met["ensg"])))))
    S = pd.DataFrame(rows, columns=["subsystem", "n_genes", "n_deg", "p", "deg_genes"])
    S["fdr"] = bh(S["p"]); S.sort_values("p").to_csv(out / "1_subsystem_enrichment_DEGs.csv", index=False)
    return met, gl, sym


# =============================================================== Stage 2: reaction activity
def stage2_reaction_activity(M, gx, groups, gl, sym, out, min_diff):
    E = {g: gx.loc[g].values for g in gx.index}
    R = reaction_matrix(M, E, gx.shape[1])
    isH = (groups.loc[gx.columns] == "HIV").values
    have = ~np.isnan(R).any(axis=1)
    idx = np.where(have)[0]
    a, b = R[idx][:, isH], R[idx][:, ~isH]
    t, p = ttest_ind(a, b, axis=1, equal_var=False)
    df = pd.DataFrame({"rxn": np.array(M["ids"])[idx], "subsystem": np.array(M["subsystem"])[idx],
                       "gpr": np.array(M["gpr"])[idx],
                       "genes": [";".join(sorted(sym.get(x, x) for x in gl[j])) for j in idx],
                       "mean_HIV": a.mean(1), "mean_control": b.mean(1), "log2_diff": a.mean(1) - b.mean(1),
                       "p": p})
    df["fdr"] = bh(np.nan_to_num(df["p"], nan=1.0))
    df["direction"] = np.where((df.fdr < 0.05) & (df.log2_diff >= min_diff), "UP",
                               np.where((df.fdr < 0.05) & (df.log2_diff <= -min_diff), "DOWN", "ns"))
    df.sort_values("p").to_csv(out / "2_reaction_activity_HIV_vs_control.csv", index=False)
    log(f"Stage 2: {len(df)} reactions with expression data; "
        f"{(df.direction == 'UP').sum()} more active, {(df.direction == 'DOWN').sum()} less active in HIV")
    # pathway enrichment; reactions sharing the same GPR within a subsystem count once
    u = df[df.subsystem != ""].drop_duplicates(["subsystem", "gpr"])
    n_all = len(u); rows = []
    for d in ("UP", "DOWN"):
        k_all = (u.direction == d).sum()
        for s, g in u.groupby("subsystem"):
            k = (g.direction == d).sum()
            if k:
                rows.append((s, d, len(g), k, hypergeom.sf(k - 1, n_all, k_all, len(g)), g.log2_diff.mean()))
    S = pd.DataFrame(rows, columns=["subsystem", "direction", "n_reactions", "n_significant", "p", "mean_log2_diff"])
    S["fdr"] = S.groupby("direction")["p"].transform(lambda x: bh(x))
    S.sort_values("p").to_csv(out / "2_subsystem_activity_enrichment.csv", index=False)
    return df, R


# =============================================================== Stage 3: flux models
def build_split(M, lb, ub):
    n = M["S"].shape[1]
    rev = np.where(lb < 0)[0]
    A = sp.hstack([M["S"], -M["S"][:, rev]]).tocsr()
    lo = np.concatenate([np.maximum(lb, 0), np.maximum(-ub[rev], 0)])
    hi = np.concatenate([np.maximum(ub, 0), -lb[rev]])
    return A, rev, lo, hi


def max_biomass(M, lb, ub):
    A, rev, lo, hi = build_split(M, lb, ub)
    c = np.zeros(A.shape[1]); c[M["bio"]] = -1
    res = linprog(c, A_eq=A, b_eq=np.zeros(A.shape[0]), bounds=np.c_[lo, hi], method="highs")
    if res.status != 0:
        sys.exit("Model cannot grow in the defined medium: " + res.message)
    return -res.fun


def pfba(M, lb, ub, growth, weights=None, w_budget=None):
    """Parsimonious flux with biomass fixed at its lower bound; optional GIMME penalty budget."""
    lb = lb.copy(); lb[M["bio"]] = max(lb[M["bio"]], growth)
    A, rev, lo, hi = build_split(M, lb, ub)
    n = len(lb)
    c = np.ones(A.shape[1])
    A_ub = b_ub = None
    if weights is not None:
        w = np.concatenate([weights, weights[rev]])
        A_ub, b_ub = sp.csr_matrix(w.reshape(1, -1)), [w_budget]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A, b_eq=np.zeros(A.shape[0]), bounds=np.c_[lo, hi],
                  method="highs")
    if res.status != 0:
        return None
    v = res.x[:n].copy(); v[rev] -= res.x[n:]
    return v


def gimme(M, expr, thr, frac, lb, ub):
    gmax = max_biomass(M, lb, ub); growth = frac * gmax
    w = np.where(np.isnan(expr), 0.0, np.maximum(thr - np.nan_to_num(expr, nan=thr), 0.0))
    lb2 = lb.copy(); lb2[M["bio"]] = max(lb2[M["bio"]], growth)
    A, rev, lo, hi = build_split(M, lb2, ub)
    cw = np.concatenate([w, w[rev]])
    r1 = linprog(cw, A_eq=A, b_eq=np.zeros(A.shape[0]), bounds=np.c_[lo, hi], method="highs")
    if r1.status != 0:
        sys.exit("GIMME LP infeasible: " + r1.message)
    incons = r1.fun
    v = pfba(M, lb, ub, growth, weights=w, w_budget=incons * (1 + 1e-4) + 1e-6)
    return v, incons, growth, gmax


def imat(M, expr, lo_t, hi_t, eps, growth, lb, ub, tlimit):
    n = len(lb)
    high = np.where(~np.isnan(expr) & (expr >= hi_t))[0]
    low = np.where(~np.isnan(expr) & (expr <= lo_t))[0]
    high = high[np.isin(high, np.where((ub >= eps) | (lb <= -eps))[0])]
    lb = lb.copy(); lb[M["bio"]] = max(lb[M["bio"]], growth)
    hp = [j for j in high if ub[j] >= eps]
    hn = [j for j in high if lb[j] <= -eps]
    nb = len(hp) + len(hn) + len(low); nv = n + nb
    rows, cols, vals, lo_c, hi_c = [], [], [], [], []
    r = 0; k = n

    def add(coefs, lo_, hi_):
        nonlocal r
        for c_, v_ in coefs:
            rows.append(r); cols.append(c_); vals.append(v_)
        lo_c.append(lo_); hi_c.append(hi_); r += 1

    yp, yn, zl = {}, {}, {}
    for j in hp:                       # v_j >= eps when y=1
        yp[j] = k; add([(j, 1.0), (k, -(eps - lb[j]))], lb[j], np.inf); k += 1
    for j in hn:                       # v_j <= -eps when y=1
        yn[j] = k; add([(j, 1.0), (k, ub[j] + eps)], -np.inf, ub[j]); k += 1
    for j in low:                      # v_j = 0 when z=1
        zl[j] = k
        add([(j, 1.0), (k, ub[j])], -np.inf, ub[j])
        add([(j, 1.0), (k, lb[j])], lb[j], np.inf); k += 1
    for j in set(yp) & set(yn):        # one direction only
        add([(yp[j], 1.0), (yn[j], 1.0)], -np.inf, 1.0)
    Ac = sp.csr_matrix((vals, (rows, cols)), shape=(r, nv))
    Aeq = sp.hstack([M["S"], sp.csr_matrix((M["S"].shape[0], nb))]).tocsr()
    cons = [LinearConstraint(Aeq, 0, 0), LinearConstraint(Ac, lo_c, hi_c)]
    c = np.zeros(nv); c[n:] = -1
    integ = np.r_[np.zeros(n), np.ones(nb)]
    bnds = Bounds(np.r_[lb, np.zeros(nb)], np.r_[ub, np.ones(nb)])
    res = milp(c, constraints=cons, integrality=integ, bounds=bnds,
               options={"time_limit": tlimit, "mip_rel_gap": 0.02, "disp": False})
    if res.x is None:
        log("  iMAT: no feasible solution found within the time limit")
        return None, None
    y = np.round(res.x[n:])
    sat_hi = sum(y[i - n] > 0.5 for i in list(yp.values()) + list(yn.values()))
    sat_lo = sum(y[i - n] > 0.5 for i in zl.values())
    v = None
    # Re-solve a parsimonious flux with the binary decisions fixed. Big-M solutions can violate a
    # strict re-solve by solver tolerance, so relax slightly step by step if needed.
    for hi_scale, lo_tol in [(1.0, 0.0), (0.5, 1e-4), (0.1, 1e-3)]:
        lb2, ub2 = lb.copy(), ub.copy()
        for j, i in yp.items():
            if y[i - n] > 0.5: lb2[j] = max(lb2[j], eps * hi_scale)
        for j, i in yn.items():
            if y[i - n] > 0.5: ub2[j] = min(ub2[j], -eps * hi_scale)
        for j, i in zl.items():
            if y[i - n] > 0.5: lb2[j] = -lo_tol; ub2[j] = lo_tol
        v = pfba(M, lb2, ub2, growth)
        if v is not None:
            break
    if v is None:
        v = res.x[:n].copy()
        log("  iMAT: parsimonious re-solve failed, using the raw MILP flux (less unique)")
    info = dict(high=len(high), low=len(low), sat_high=int(sat_hi), sat_low=int(sat_lo), status=res.message[:40])
    return v, info


def stage3_models(args, M, gx, groups, gl, sym, out):
    lb, ub = medium_bounds(M)
    E = {}
    res = {}
    for g in ("HIV", "control"):
        cols = groups[groups == g].index.intersection(gx.columns)
        E[g] = {x: np.array([gx.loc[x, cols].median()]) for x in gx.index}
    Rg = {g: reaction_matrix(M, E[g], 1)[:, 0] for g in E}
    pooled = np.concatenate([Rg[g][~np.isnan(Rg[g])] for g in Rg])
    lo_t, hi_t = np.quantile(pooled, args.low_q), np.quantile(pooled, args.high_q)
    log(f"Stage 3: thresholds (log2 expression) low={lo_t:.2f} (q{args.low_q}), high={hi_t:.2f} (q{args.high_q}); "
        f"{(~np.isnan(Rg['HIV'])).sum()} reactions with data")
    out_flux = {}
    for method in (["gimme", "imat"] if args.method == "both" else [args.method]):
        for g in ("HIV", "control"):
            t0 = time.time()
            if method == "gimme":
                v, inc, growth, gmax = gimme(M, Rg[g], lo_t, args.growth_frac, lb, ub)
                log(f"  GIMME {g}: inconsistency score {inc:.3f}, max biomass {gmax:.1f}, {time.time()-t0:.0f}s")
            else:
                gmax = max_biomass(M, lb, ub)
                v, info = imat(M, Rg[g], lo_t, hi_t, args.imat_eps, args.growth_frac * gmax, lb, ub,
                               args.imat_time_limit)
                log(f"  iMAT {g}: {info}, {time.time()-t0:.0f}s")
            if v is None:
                log(f"  {method} {g}: failed; skipping"); continue
            out_flux[(method, g)] = v
    return out_flux, Rg


def stage4_compare(M, out_flux, Rg, act, met, gl, sym, out, args):
    ids, ss = np.array(M["ids"]), np.array(M["subsystem"])
    genes = np.array([";".join(sorted(sym.get(x, x) for x in l)) for l in gl])
    deg_genes = set(met["symbol"]) if met is not None else set()
    summ = []
    for method in sorted({m for m, _ in out_flux}):
        if (method, "HIV") not in out_flux or (method, "control") not in out_flux:
            continue
        vh, vc = np.abs(out_flux[(method, "HIV")]), np.abs(out_flux[(method, "control")])
        d = pd.DataFrame({"rxn": ids, "subsystem": ss, "genes": genes, "abs_flux_HIV": vh,
                          "abs_flux_control": vc, "diff": vh - vc,
                          "log2_fold": np.log2((vh + 1e-3) / (vc + 1e-3))})
        d["flux_change"] = np.where((d["diff"].abs() >= args.flux_min_diff) & (d.log2_fold.abs() >= 1),
                                    np.where(d["diff"] > 0, "UP", "DOWN"), "ns")
        d = d.merge(act[["rxn", "log2_diff", "fdr", "direction"]].rename(
            columns={"log2_diff": "expr_log2_diff", "fdr": "expr_fdr", "direction": "expr_direction"}),
            on="rxn", how="left")
        d["has_DEG_enzyme"] = [any(x in deg_genes for x in str(g).split(";")) for g in d["genes"]]
        d.sort_values("diff", key=abs, ascending=False).to_csv(out / f"3_flux_comparison_{method}.csv", index=False)
        s = d[d.subsystem != ""].groupby("subsystem").agg(
            n_reactions=("rxn", "size"), flux_HIV=("abs_flux_HIV", "sum"), flux_control=("abs_flux_control", "sum"),
            n_up=("flux_change", lambda x: (x == "UP").sum()), n_down=("flux_change", lambda x: (x == "DOWN").sum()))
        s["flux_diff"] = s.flux_HIV - s.flux_control
        s.sort_values("flux_diff", key=abs, ascending=False).to_csv(out / f"3_subsystem_flux_{method}.csv")
        agree = d[(d.flux_change != "ns") & (d.expr_direction.isin(["UP", "DOWN"])) & (d.flux_change == d.expr_direction)]
        log(f"Stage 4 [{method}]: {(d.flux_change=='UP').sum()} reactions with higher flux in HIV, "
            f"{(d.flux_change=='DOWN').sum()} lower; {len(agree)} also supported by differential expression")
        summ.append((method, d))
    return summ


def candidate_targets(summ, act, met, out):
    """One row per enzyme gene. Evidence sources: DEG, reaction activity, GIMME flux, iMAT flux.
    A source votes UP/DOWN (in HIV vs control) only if all its signals for that gene agree."""
    votes = {}                                   # gene -> source -> set of directions

    def add(g, src, d):
        votes.setdefault(g, {}).setdefault(src, set()).add(d)

    if met is not None:
        for _, r in met.iterrows():
            add(r["symbol"], "DEG", "UP" if r["log2FoldChange"] > 0 else "DOWN")
    for _, r in act[act.direction != "ns"].iterrows():
        gs = [g for g in str(r["genes"]).split(";") if g]
        if 0 < len(gs) <= 3:                     # skip large promiscuous gene sets
            for g in gs:
                add(g, "reaction_activity", r["direction"])
    for method, d in summ:
        for _, r in d[d.flux_change != "ns"].iterrows():
            gs = [g for g in str(r["genes"]).split(";") if g]
            if 0 < len(gs) <= 3:
                for g in gs:
                    add(g, method + "_flux", r["flux_change"])
    rows = []
    for g, srcs in votes.items():
        v = {s: (list(d)[0] if len(d) == 1 else "mixed") for s, d in srcs.items()}
        up = sorted(s for s, d in v.items() if d == "UP")
        dn = sorted(s for s, d in v.items() if d == "DOWN")
        mixed = sorted(s for s, d in v.items() if d == "mixed")
        direction = "UP" if len(up) > len(dn) else "DOWN" if len(dn) > len(up) else "conflict"
        n = max(len(up), len(dn)) if direction != "conflict" else 0
        rows.append((g, direction, n, ";".join(up if direction != "DOWN" else dn),
                     ";".join(dn if direction != "DOWN" else up), ";".join(mixed)))
    T = pd.DataFrame(rows, columns=["gene", "direction_in_HIV", "n_agreeing_sources",
                                    "agreeing_sources", "opposing_sources", "mixed_sources"])
    T = T.sort_values(["n_agreeing_sources", "gene"], ascending=[False, True])
    T.to_csv(out / "4_candidate_target_enzymes.csv", index=False)
    log(f"Stage 4: {len(T)} enzymes with any evidence; {(T.n_agreeing_sources >= 3).sum()} supported by >=3 "
        f"agreeing sources, {(T.n_agreeing_sources >= 2).sum()} by >=2")
    return T


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--degs", default="Supplementary_-_1.xlsx")
    ap.add_argument("--sheet", default="HIV DEGs")
    ap.add_argument("--counts", help="local count matrix (genes x samples, first column = gene ID)")
    ap.add_argument("--groups", help="CSV with columns sample,group (HIV/control)")
    ap.add_argument("--layout", choices=["auto", "samples-rows", "genes-rows"], default="auto",
                    help="orientation of the count file (auto: fewer rows than columns = samples in rows)")
    ap.add_argument("--hiv-code", help="value of the 'studygroup' column that means HIV (others = control)")
    ap.add_argument("--sample-col", help="column holding sample IDs when samples are in rows (default: first)")
    ap.add_argument("--outdir", default="results_2b")
    ap.add_argument("--method", choices=["gimme", "imat", "both"], default="gimme")
    ap.add_argument("--growth-frac", type=float, default=0.1, help="required fraction of max biomass")
    ap.add_argument("--low-q", type=float, default=0.25, help="low-expression quantile (GIMME threshold, iMAT low)")
    ap.add_argument("--high-q", type=float, default=0.75, help="iMAT high-expression quantile")
    ap.add_argument("--imat-eps", type=float, default=0.1)
    ap.add_argument("--imat-time-limit", type=int, default=900, help="seconds per iMAT condition")
    ap.add_argument("--min-log2-diff", type=float, default=0.5, help="reaction activity effect size cutoff")
    ap.add_argument("--flux-min-diff", type=float, default=0.01)
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    M = load_model(out)
    log(f"Human-GEM {M['meta'].get('version')}: {len(M['ids'])} reactions, biomass reaction {M['ids'][M['bio']]}")

    met = gl = sym = None
    try:
        met, gl, sym = stage1_deg_mapping(a, M, out)
    except FileNotFoundError:
        log(f"DEG file {a.degs} not found; skipping Stage 1")
    if gl is None:
        gl = [sorted(set(re.findall(r"ENSG\d+", g))) for g in M["gpr"]]
        sym = dict(zip(M["genes"]["genes"], M["genes"]["geneSymbols"]))

    X, groups = load_counts_and_groups(a, out)
    logx = normalise(X)
    sanity_check_groups(logx, groups)
    gx = gene_expression_for_model(logx, M)
    gx.to_csv(out / "expression_log2norm_human1_genes.csv.gz")

    act, _ = stage2_reaction_activity(M, gx, groups, gl, sym, out, a.min_log2_diff)
    flux, Rg = stage3_models(a, M, gx, groups, gl, sym, out)
    summ = stage4_compare(M, flux, Rg, act, met, gl, sym, out, a)
    candidate_targets(summ, act, met, out)
    log("Done. Results in", out.resolve())


if __name__ == "__main__":
    main()