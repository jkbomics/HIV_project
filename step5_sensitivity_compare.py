#!/usr/bin/env python3
"""
Step 5: sensitivity analysis. Compare the signature-reversal results of the baseline run with runs that
use a modified signature (no immunoglobulin genes, no core proliferation genes, or both), and report which
candidates stay reversers.

First produce the variant runs with step2a_signature_reversal.py, for example:
    python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --exclude-ig --outdir results_noIG
    python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --exclude-proliferation --outdir results_noProlif
    python step2a_signature_reversal.py --degs Supplementary_-_1.xlsx --exclude-ig --exclude-proliferation --outdir results_noIG_noProlif

Then compare:
    python step5_sensitivity_compare.py \\
        --variant noIG=results_noIG/sigcom_compound_summary.csv \\
        --variant noProlif=results_noProlif/sigcom_compound_summary.csv \\
        --variant noIG_noProlif=results_noIG_noProlif/sigcom_compound_summary.csv

Output (results_5/):
    8_sensitivity_robustness.csv   one row per focus candidate, with its reversal status in every run
    8_sensitivity_run_summary.csv  one row per run: number of reversers and overlap with the baseline
"""
import argparse
import re
from pathlib import Path

import pandas as pd

SALTS = r"(hydrochloride|hcl|sodium|potassium|mesylate|acetate|phosphate|sulfate|sulphate|citrate|maleate|" \
        r"tartrate|besylate|tosylate|bromide|chloride|fumarate|succinate|calcium|disodium|dihydrate|hydrate)"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", re.sub(SALTS, " ", str(s).lower()))


def load(path):
    S = pd.read_csv(path)
    S["key"] = S["compound"].map(norm)
    S = S.drop_duplicates("key").set_index("key")
    S["reverser"] = (S["n_reversing"] >= 3) & (S["n_reversing"] > 2 * S["n_mimicking"])
    return S


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", default="results_v2/sigcom_compound_summary.csv")
    ap.add_argument("--variant", action="append", required=True, help="name=path/to/sigcom_compound_summary.csv")
    ap.add_argument("--ranking", default="results_4/7_final_candidate_ranking.csv")
    ap.add_argument("--shortlist", default="results_4/7_repurposing_shortlist.csv")
    ap.add_argument("--outdir", default="results_5")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    runs = {"baseline": load(a.baseline)}
    for v in a.variant:
        name, path = v.split("=", 1)
        runs[name] = load(path)

    # ---- run-level summary
    base_rev = set(runs["baseline"].index[runs["baseline"]["reverser"]])
    rows = []
    for name, S in runs.items():
        rev = set(S.index[S["reverser"]])
        rows.append({"run": name, "n_reversers": len(rev),
                     "n_shared_with_baseline": len(rev & base_rev),
                     "pct_of_baseline_reversers_retained": round(100 * len(rev & base_rev) / max(len(base_rev), 1), 1),
                     "n_new_reversers_not_in_baseline": len(rev - base_rev)})
    runsum = pd.DataFrame(rows)
    runsum.to_csv(out / "8_sensitivity_run_summary.csv", index=False)

    # ---- focus candidates
    R = pd.read_csv(a.ranking)
    R["key"] = R["drug"].map(norm)
    focus = R[R["tier"].str.startswith(("1", "2:"))].copy()
    focus["group"] = focus["tier"].str.split(":").str[0].map({"1": "tier 1", "2": "tier 2 (strong reversal)"})
    short = pd.read_csv(a.shortlist)
    short["key"] = short["drug"].map(norm)
    short["group"] = "repurposing shortlist"
    focus = pd.concat([focus, short]).drop_duplicates("key")

    recs = []
    for _, f in focus.iterrows():
        rec = {"drug": f["drug"], "group": f["group"], "drug_class": f.get("drug_class")}
        n_rev_runs = 0
        for name, S in runs.items():
            if f["key"] in S.index:
                s = S.loc[f["key"]]
                rec[f"{name}_n_reversing"] = s["n_reversing"]; rec[f"{name}_n_mimicking"] = s["n_mimicking"]
                rec[f"{name}_cell_lines"] = s["n_cell_lines_reversing"]; rec[f"{name}_reverser"] = bool(s["reverser"])
                n_rev_runs += bool(s["reverser"])
            else:
                rec[f"{name}_n_reversing"] = 0; rec[f"{name}_n_mimicking"] = None
                rec[f"{name}_cell_lines"] = 0; rec[f"{name}_reverser"] = False
        rec["n_runs_reverser"] = n_rev_runs
        rec["robust_all_runs"] = n_rev_runs == len(runs)
        recs.append(rec)
    T = pd.DataFrame(recs).sort_values(["robust_all_runs", "n_runs_reverser", "group"], ascending=[False, False, True])
    T.to_csv(out / "8_sensitivity_robustness.csv", index=False)

    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 20)
    print("Run summary:")
    print(runsum.to_string(index=False))
    print("\nRobustness of focus candidates (reverser in how many of the", len(runs), "runs):")
    print(T.groupby("group")["robust_all_runs"].agg(["sum", "count"]).rename(
        columns={"sum": "reverser_in_all_runs", "count": "n_candidates"}).to_string())
    print("\nTier 1 candidates:")
    print(T[T["group"] == "tier 1"][["drug", "n_runs_reverser"] + [f"{n}_reverser" for n in runs]].to_string(index=False))
    print("\nFiles written to", out.resolve())


if __name__ == "__main__":
    main()
