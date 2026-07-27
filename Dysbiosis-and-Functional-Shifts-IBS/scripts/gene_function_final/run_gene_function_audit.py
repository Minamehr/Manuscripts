#!/usr/bin/env python3
"""Final sample-level functional audit for HC versus IBS."""

from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal


def sep(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    return out


def orient(df, orientation):
    return df if orientation == "features-by-samples" else df.T


def aggregate_genes(table, mapping, gene_column, function_column, method):
    table = table.copy()
    table.index = table.index.astype(str)
    mapping = mapping[[gene_column, function_column]].dropna().copy()
    mapping[gene_column] = mapping[gene_column].astype(str)
    mapping[function_column] = mapping[function_column].astype(str)
    merged = mapping.merge(table.reset_index(names=gene_column), on=gene_column, how="inner")
    numeric = [c for c in merged if c not in {gene_column, function_column}]
    grouped = merged.groupby(function_column, sort=False)[numeric]
    return grouped.sum() if method == "sum" else grouped.mean()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--abundance", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--analysis-name", required=True)
    p.add_argument("--orientation", choices=["features-by-samples", "samples-by-features"], default="features-by-samples")
    p.add_argument("--mapping", type=Path)
    p.add_argument("--gene-column", default="gene_id")
    p.add_argument("--function-column", default="function")
    p.add_argument("--aggregation", choices=["sum", "mean"], default="sum")
    p.add_argument("--sample-column", default="SampleID")
    p.add_argument("--group-column", default="GeneralGroup")
    p.add_argument("--subtype-column", default="Subtype")
    p.add_argument("--hc-label", default="HC")
    p.add_argument("--ibs-label", default="IBS")
    p.add_argument("--exclude-sample", action="append", default=[])
    p.add_argument("--minimum-overall-prevalence", type=float, default=.10)
    p.add_argument("--minimum-group-prevalence", type=float, default=0)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(args.abundance, sep=sep(args.abundance), index_col=0)
    table = orient(table, args.orientation)
    if args.mapping:
        mapping = pd.read_csv(args.mapping, sep=sep(args.mapping))
        table = aggregate_genes(table, mapping, args.gene_column, args.function_column, args.aggregation)

    meta = pd.read_csv(args.metadata, sep=sep(args.metadata))
    meta[args.sample_column] = meta[args.sample_column].astype(str)
    meta = meta.set_index(args.sample_column, drop=False)
    table.columns = table.columns.astype(str)
    samples = [s for s in meta.index if s in table.columns and s not in set(args.exclude_sample)]
    meta, table = meta.loc[samples], table.loc[:, samples]
    valid = meta[args.group_column].isin([args.hc_label, args.ibs_label])
    meta, table = meta.loc[valid], table.loc[:, meta.index[valid]]

    hc_samples = meta.index[meta[args.group_column] == args.hc_label]
    ibs_samples = meta.index[meta[args.group_column] == args.ibs_label]
    overall = table.gt(0).mean(axis=1)
    hprev = table[hc_samples].gt(0).mean(axis=1)
    iprev = table[ibs_samples].gt(0).mean(axis=1)
    keep = ((overall >= args.minimum_overall_prevalence) &
            (hprev >= args.minimum_group_prevalence) &
            (iprev >= args.minimum_group_prevalence))
    filtered = table.loc[keep]

    rows = []
    for feature, values in filtered.iterrows():
        h = values[hc_samples].to_numpy(float)
        i = values[ibs_samples].to_numpy(float)
        stat, pv = mannwhitneyu(h, i, alternative="two-sided")
        rows.append({
            "feature": feature,
            "n_hc": len(h), "n_ibs": len(i),
            "prevalence_hc": np.mean(h > 0), "prevalence_ibs": np.mean(i > 0),
            "mean_hc": np.mean(h), "mean_ibs": np.mean(i),
            "median_hc": np.median(h), "median_ibs": np.median(i),
            "mean_difference_ibs_minus_hc": np.mean(i) - np.mean(h),
            "median_difference_ibs_minus_hc": np.median(i) - np.median(h),
            "mannwhitney_u": stat, "p_value": pv,
        })
    results = pd.DataFrame(rows)
    if not results.empty:
        results["bh_q"] = bh(results["p_value"])
        results = results.sort_values(["bh_q", "p_value", "feature"])
    results.to_csv(args.output_dir / f"{args.analysis_name}_hc_vs_ibs_all.csv", index=False)
    results.loc[results["bh_q"] < .05].to_csv(
        args.output_dir / f"{args.analysis_name}_hc_vs_ibs_q_lt_0.05.csv", index=False)

    subtype = []
    if args.subtype_column in meta.columns:
        labels = list(meta[args.subtype_column].dropna().astype(str).unique())
        for feature, values in filtered.iterrows():
            arrays = [values[meta.index[meta[args.subtype_column].astype(str) == label]].to_numpy(float)
                      for label in labels]
            arrays = [x for x in arrays if len(x)]
            if len(arrays) >= 2:
                stat, pv = kruskal(*arrays)
                subtype.append({"feature": feature, "kruskal_h": stat, "p_value": pv})
    subtype_df = pd.DataFrame(subtype)
    if not subtype_df.empty:
        subtype_df["bh_q"] = bh(subtype_df["p_value"])
        subtype_df = subtype_df.sort_values(["bh_q", "p_value", "feature"])
    subtype_df.to_csv(args.output_dir / f"{args.analysis_name}_subtype_kruskal_all.csv", index=False)

    summary = {
        "analysis_name": args.analysis_name,
        "input_features": int(table.shape[0]),
        "retained_features": int(filtered.shape[0]),
        "samples": int(filtered.shape[1]),
        "hc_samples": int(len(hc_samples)),
        "ibs_samples": int(len(ibs_samples)),
        "minimum_overall_prevalence": args.minimum_overall_prevalence,
        "minimum_group_prevalence": args.minimum_group_prevalence,
        "aggregation": args.aggregation if args.mapping else "input_already_aggregated",
        "hc_vs_ibs_features_q_lt_0.05": int((results["bh_q"] < .05).sum()) if not results.empty else 0,
    }
    (args.output_dir / f"{args.analysis_name}_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
