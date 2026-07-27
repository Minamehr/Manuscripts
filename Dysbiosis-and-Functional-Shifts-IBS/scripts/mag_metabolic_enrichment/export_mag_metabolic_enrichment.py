#!/usr/bin/env python3
"""Transparent export of anvi'o metabolic-enrichment results."""

from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd


def sep(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","


def choose(columns, candidates, explicit=None):
    if explicit:
        if explicit not in columns:
            raise ValueError(f"Column not found: {explicit}")
        return explicit
    lookup = {c.lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    raise ValueError(f"Could not find one of: {candidates}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--q-column")
    p.add_argument("--module-column")
    p.add_argument("--group-column")
    p.add_argument("--p-column")
    p.add_argument("--n-hc-column")
    p.add_argument("--n-ibs-column")
    p.add_argument("--expected-significant", type=int, default=24)
    p.add_argument("--expected-hc-units", type=int, default=94)
    p.add_argument("--expected-ibs-units", type=int, default=127)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, sep=sep(args.input))
    columns = list(df.columns)
    qcol = choose(columns, ["adjusted_q_value", "q_value", "qvalue", "q"], args.q_column)
    mcol = choose(columns, ["MODULE", "module", "module_name", "metabolic_module"], args.module_column)
    gcol = choose(columns, ["enriched_group", "group", "enrichment_group"], args.group_column)
    try:
        pcol = choose(columns, ["unadjusted_p_value", "p_value", "pvalue", "p"], args.p_column)
    except ValueError:
        pcol = None

    df[qcol] = pd.to_numeric(df[qcol], errors="coerce")
    df = df.sort_values([qcol, mcol], na_position="last")
    significant = df.loc[df[qcol] < .05].copy()
    df.to_csv(args.output_dir / "metabolic_enrichment_all.csv", index=False)
    significant.to_csv(args.output_dir / "metabolic_enrichment_q_lt_0.05.csv", index=False)

    validation = {
        "input_rows": int(len(df)),
        "significant_modules_q_lt_0.05": int(len(significant)),
        "expected_significant_modules": args.expected_significant,
        "significant_count_matches_lock": bool(len(significant) == args.expected_significant),
        "significant_by_enriched_group": {
            str(k): int(v) for k, v in significant[gcol].astype(str).value_counts().to_dict().items()
        },
        "q_column": qcol, "module_column": mcol, "group_column": gcol, "p_column": pcol,
        "significance_rule": "BH-adjusted q < 0.05 only",
        "unsupported_filters_applied": False,
        "expected_analysis_units": {
            "HC": args.expected_hc_units,
            "IBS": args.expected_ibs_units,
            "total": args.expected_hc_units + args.expected_ibs_units,
        },
    }
    if args.n_hc_column and args.n_ibs_column:
        h = pd.to_numeric(df[args.n_hc_column], errors="coerce").dropna().unique()
        i = pd.to_numeric(df[args.n_ibs_column], errors="coerce").dropna().unique()
        validation["observed_n_hc_values"] = [float(x) for x in h]
        validation["observed_n_ibs_values"] = [float(x) for x in i]
        validation["analysis_unit_counts_match_lock"] = bool(
            args.expected_hc_units in h and args.expected_ibs_units in i)
    (args.output_dir / "metabolic_enrichment_validation.json").write_text(
        json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
