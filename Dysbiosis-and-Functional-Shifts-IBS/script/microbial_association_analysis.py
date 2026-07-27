#!/usr/bin/env python3
"""HC-versus-IBS microbial association analysis with resampling and sensitivity checks."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, pearsonr


def sep(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    return out


def orient(df: pd.DataFrame, orientation: str) -> pd.DataFrame:
    return df.T if orientation == "features-by-samples" else df


def samplewise_multiplicative(matrix: np.ndarray, delta: float | None = None) -> np.ndarray:
    rows = []
    for row in np.asarray(matrix, float):
        if np.any(row < 0) or row.sum() <= 0:
            raise ValueError("Each sample must contain non-negative values and at least one positive value.")
        x = row / row.sum()
        zero = x == 0
        z = int(zero.sum())
        if z:
            d = len(x)
            replacement = 1.0 / (d * d) if delta is None else float(delta)
            if z * replacement >= 1:
                raise ValueError("Zero-replacement value is too large.")
            x[zero] = replacement
            x[~zero] *= (1 - z * replacement) / x[~zero].sum()
        rows.append(x)
    return np.vstack(rows)


def feature_half_minimum(matrix: np.ndarray) -> np.ndarray:
    x = np.asarray(matrix, float).copy()
    for j in range(x.shape[1]):
        positive = x[:, j][x[:, j] > 0]
        replacement = float(positive.min() / 2.0) if len(positive) else 1e-12
        x[x[:, j] == 0, j] = replacement
    x = x / x.sum(axis=1, keepdims=True)
    return x


def clr_transform(matrix: np.ndarray, zero_method: str, delta: float | None = None) -> np.ndarray:
    if zero_method == "multiplicative":
        closed = samplewise_multiplicative(matrix, delta)
    elif zero_method == "feature-half-minimum":
        closed = feature_half_minimum(matrix)
    else:
        raise ValueError(zero_method)
    logged = np.log(closed)
    return logged - logged.mean(axis=1, keepdims=True)


def correlation_matrix(x: np.ndarray, estimator: str) -> np.ndarray:
    values = np.apply_along_axis(rankdata, 0, x) if estimator == "spearman" else x
    corr = np.corrcoef(values, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def upper(corr: np.ndarray):
    tri = np.triu_indices_from(corr, k=1)
    return corr[tri], tri


def domain_mask(features, domain_map, tri):
    if not domain_map:
        return None
    domains = np.array([domain_map.get(f, "Unknown") for f in features], dtype=object)
    return domains[tri[0]] != domains[tri[1]]


def network_metrics(corr: np.ndarray, threshold: float, cross=None):
    vals, _ = upper(corr)
    out = {
        "mean_abs_correlation": float(np.mean(np.abs(vals))),
        "edge_density": float(np.mean(np.abs(vals) >= threshold)),
        "positive_edge_density": float(np.mean(vals >= threshold)),
        "negative_edge_density": float(np.mean(vals <= -threshold)),
    }
    if cross is not None and cross.any():
        cv = vals[cross]
        out["cross_domain_mean_abs_correlation"] = float(np.mean(np.abs(cv)))
        out["cross_domain_edge_density"] = float(np.mean(np.abs(cv) >= threshold))
    return out


def metric_differences(hc_corr, ibs_corr, threshold, cross=None):
    hc = network_metrics(hc_corr, threshold, cross)
    ibs = network_metrics(ibs_corr, threshold, cross)
    out = {name: hc[name] - ibs[name] for name in hc}
    hv, _ = upper(hc_corr)
    iv, _ = upper(ibs_corr)
    out["global_rms_matrix_difference"] = float(np.sqrt(np.mean((hv - iv) ** 2)))
    return out


def prepare_table(raw, metadata, args, prevalence):
    raw = raw.copy()
    raw.index = raw.index.astype(str)
    metadata = metadata.copy()
    metadata[args.sample_column] = metadata[args.sample_column].astype(str)
    metadata = metadata.set_index(args.sample_column, drop=False)
    samples = [s for s in metadata.index if s in raw.index and s not in set(args.exclude_sample)]
    x = raw.loc[samples]
    meta = metadata.loc[samples]
    valid = meta[args.group_column].isin([args.hc_label, args.ibs_label])
    x, meta = x.loc[valid], meta.loc[valid]
    groups = meta[args.group_column].astype(str)
    hc = x.loc[groups == args.hc_label]
    ibs = x.loc[groups == args.ibs_label]
    keep = (hc.gt(0).mean(axis=0) >= prevalence) & (ibs.gt(0).mean(axis=0) >= prevalence)
    x = x.loc[:, keep]
    if x.shape[1] < 2:
        raise ValueError("Fewer than two features remain after prevalence filtering.")
    return x, groups, list(x.columns)


def group_corrs(clr, groups, args, estimator):
    return (
        correlation_matrix(clr[groups == args.hc_label], estimator),
        correlation_matrix(clr[groups == args.ibs_label], estimator),
    )


def permutation_test(clr, groups, args, cross, rng):
    hc, ibs = group_corrs(clr, groups, args, "pearson")
    observed = metric_differences(hc, ibs, args.threshold, cross)
    null = {name: np.empty(args.permutations) for name in observed}
    for i in range(args.permutations):
        permuted = rng.permutation(groups)
        ph, pi = group_corrs(clr, permuted, args, "pearson")
        diff = metric_differences(ph, pi, args.threshold, cross)
        for name, value in diff.items():
            null[name][i] = value
    rows = []
    for name, value in observed.items():
        values = null[name]
        if name == "global_rms_matrix_difference":
            p = (1 + np.sum(values >= value)) / (args.permutations + 1)
        else:
            p = (1 + np.sum(np.abs(values) >= abs(value))) / (args.permutations + 1)
        rows.append({
            "metric": name,
            "observed_hc_minus_ibs": value,
            "permutation_p": p,
            "null_mean": float(values.mean()),
            "null_sd": float(values.std(ddof=1)),
        })
    result = pd.DataFrame(rows)
    result["bh_q"] = bh(result["permutation_p"])
    return result


def resampling(clr, groups, args, cross, rng, iterations, replace):
    hc_idx = np.flatnonzero(groups == args.hc_label)
    ibs_idx = np.flatnonzero(groups == args.ibs_label)
    n = min(len(hc_idx), len(ibs_idx))
    rows = []
    for i in range(iterations):
        h = rng.choice(hc_idx, n, replace=replace)
        b = rng.choice(ibs_idx, n, replace=replace)
        diff = metric_differences(
            correlation_matrix(clr[h], "pearson"),
            correlation_matrix(clr[b], "pearson"),
            args.threshold,
            cross,
        )
        diff["iteration"] = i + 1
        rows.append(diff)
    return pd.DataFrame(rows)


def summarise_resampling(table, name):
    rows = []
    for column in table.columns:
        if column == "iteration":
            continue
        values = table[column].to_numpy(float)
        rows.append({
            "analysis": name,
            "metric": column,
            "mean": float(values.mean()),
            "median": float(np.median(values)),
            "ci_low_2.5": float(np.quantile(values, 0.025)),
            "ci_high_97.5": float(np.quantile(values, 0.975)),
            "fraction_positive": float(np.mean(values > 0)),
        })
    return pd.DataFrame(rows)


def edge_stability(clr_group, features, args, rng):
    tri = np.triu_indices(len(features), k=1)
    boot = np.empty((args.edge_bootstraps, len(tri[0])), dtype=np.float32)
    for i in range(args.edge_bootstraps):
        idx = rng.choice(np.arange(len(clr_group)), len(clr_group), replace=True)
        boot[i] = correlation_matrix(clr_group[idx], "pearson")[tri]
    mean = boot.mean(axis=0)
    low = np.quantile(boot, 0.025, axis=0)
    high = np.quantile(boot, 0.975, axis=0)
    positive_fraction = (boot >= args.threshold).mean(axis=0)
    negative_fraction = (boot <= -args.threshold).mean(axis=0)
    same_sign_threshold_fraction = np.maximum(positive_fraction, negative_fraction)
    stable = (same_sign_threshold_fraction >= args.edge_stability_fraction) & ((low > 0) | (high < 0))
    return pd.DataFrame({
        "feature_1": np.asarray(features)[tri[0]],
        "feature_2": np.asarray(features)[tri[1]],
        "bootstrap_mean_r": mean,
        "ci_low_2.5": low,
        "ci_high_97.5": high,
        "same_sign_threshold_fraction": same_sign_threshold_fraction,
        "stable_edge": stable,
    })


def matrix_concordance(path, observed, features):
    ref = pd.read_csv(path, sep=sep(path), index_col=0)
    common = [f for f in features if f in ref.index and f in ref.columns]
    if len(common) < 3:
        raise ValueError("Fewer than three common features for matrix concordance.")
    idx = [features.index(f) for f in common]
    obs = observed[np.ix_(idx, idx)]
    reference = ref.loc[common, common].to_numpy(float)
    tri = np.triu_indices(len(common), k=1)
    value, p = pearsonr(obs[tri], reference[tri])
    return {"n_common_features": len(common), "matrix_pearson_r": float(value), "p_value": float(p)}


def sensitivity_analysis(raw, meta, args, rng):
    rows = []
    prevalences = [float(x) for x in args.prevalence_sensitivity.split(",")]
    thresholds = [float(x) for x in args.thresholds.split(",")]
    for prevalence in prevalences:
        try:
            x, groups_series, features = prepare_table(raw, meta, args, prevalence)
        except ValueError:
            continue
        groups = groups_series.to_numpy(str)
        mean_abundance = x.mean(axis=0)
        q25, q75 = mean_abundance.quantile([0.25, 0.75])
        strata = {
            "all": list(x.columns),
            "lowest_abundance_quartile": list(mean_abundance.index[mean_abundance <= q25]),
            "highest_abundance_quartile": list(mean_abundance.index[mean_abundance >= q75]),
        }
        for stratum, selected in strata.items():
            if len(selected) < 2:
                continue
            sx = x.loc[:, selected]
            for zero_method in ["multiplicative", "feature-half-minimum"]:
                transformed = clr_transform(sx.to_numpy(float), zero_method, args.zero_delta)
                for estimator in ["pearson", "spearman"]:
                    hc, ibs = group_corrs(transformed, groups, args, estimator)
                    for threshold in thresholds:
                        for metric, value in metric_differences(hc, ibs, threshold).items():
                            rows.append({
                                "prevalence": prevalence,
                                "abundance_stratum": stratum,
                                "zero_method": zero_method,
                                "estimator": estimator,
                                "threshold": threshold,
                                "n_features": len(selected),
                                "metric": metric,
                                "hc_minus_ibs": value,
                            })
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--abundance", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--layer-name", required=True)
    p.add_argument("--orientation", choices=["features-by-samples", "samples-by-features"], default="features-by-samples")
    p.add_argument("--sample-column", default="SampleID")
    p.add_argument("--group-column", default="GeneralGroup")
    p.add_argument("--hc-label", default="HC")
    p.add_argument("--ibs-label", default="IBS")
    p.add_argument("--exclude-sample", action="append", default=[])
    p.add_argument("--feature-metadata", type=Path)
    p.add_argument("--feature-column", default="feature")
    p.add_argument("--domain-column", default="domain")
    p.add_argument("--prevalence", type=float, default=0.50)
    p.add_argument("--prevalence-sensitivity", default="0.20,0.25,0.30,0.40,0.60,0.70,0.80")
    p.add_argument("--threshold", type=float, default=0.40)
    p.add_argument("--thresholds", default="0.35,0.40,0.50,0.60")
    p.add_argument("--zero-delta", type=float)
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--equal-size-iterations", type=int, default=500)
    p.add_argument("--balanced-bootstraps", type=int, default=500)
    p.add_argument("--edge-bootstraps", type=int, default=500)
    p.add_argument("--edge-stability-fraction", type=float, default=0.70)
    p.add_argument("--seed", type=int, default=20260722)
    p.add_argument("--sparcc-hc", type=Path)
    p.add_argument("--sparcc-ibs", type=Path)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    raw = pd.read_csv(args.abundance, sep=sep(args.abundance), index_col=0)
    raw = orient(raw, args.orientation)
    meta = pd.read_csv(args.metadata, sep=sep(args.metadata))

    x, group_series, features = prepare_table(raw, meta, args, args.prevalence)
    groups = group_series.to_numpy(str)
    transformed = clr_transform(x.to_numpy(float), "multiplicative", args.zero_delta)
    hc_corr, ibs_corr = group_corrs(transformed, groups, args, "pearson")
    _, tri = upper(hc_corr)

    domain_map = None
    if args.feature_metadata:
        fm = pd.read_csv(args.feature_metadata, sep=sep(args.feature_metadata))
        domain_map = dict(zip(fm[args.feature_column].astype(str), fm[args.domain_column].astype(str)))
    cross = domain_mask(features, domain_map, tri)

    descriptive = []
    for label, matrix in [(args.hc_label, hc_corr), (args.ibs_label, ibs_corr)]:
        for threshold in [float(x) for x in args.thresholds.split(",")]:
            row = network_metrics(matrix, threshold, cross)
            row.update({"group": label, "threshold": threshold, "n_features": len(features)})
            descriptive.append(row)
    pd.DataFrame(descriptive).to_csv(args.output_dir / "descriptive_metrics.csv", index=False)

    permutation_test(transformed, groups, args, cross, rng).to_csv(
        args.output_dir / "primary_permutation_tests.csv", index=False)

    equal = resampling(transformed, groups, args, cross, rng, args.equal_size_iterations, replace=False)
    balanced = resampling(transformed, groups, args, cross, rng, args.balanced_bootstraps, replace=True)
    equal.to_csv(args.output_dir / "equal_size_iterations.csv", index=False)
    balanced.to_csv(args.output_dir / "balanced_bootstrap_iterations.csv", index=False)
    pd.concat([
        summarise_resampling(equal, "equal_size_without_replacement"),
        summarise_resampling(balanced, "balanced_bootstrap"),
    ], ignore_index=True).to_csv(args.output_dir / "resampling_summary.csv", index=False)

    edge_stability(transformed[groups == args.hc_label], features, args, rng).to_csv(
        args.output_dir / "edge_stability_hc.csv.gz", index=False, compression="gzip")
    edge_stability(transformed[groups == args.ibs_label], features, args, rng).to_csv(
        args.output_dir / "edge_stability_ibs.csv.gz", index=False, compression="gzip")

    sensitivity_analysis(raw, meta, args, rng).to_csv(
        args.output_dir / "sensitivity_analysis.csv", index=False)

    pd.DataFrame(hc_corr, index=features, columns=features).to_csv(args.output_dir / "pearson_hc.csv")
    pd.DataFrame(ibs_corr, index=features, columns=features).to_csv(args.output_dir / "pearson_ibs.csv")

    concordance = {}
    if args.sparcc_hc and args.sparcc_ibs:
        concordance["HC"] = matrix_concordance(args.sparcc_hc, hc_corr, features)
        concordance["IBS"] = matrix_concordance(args.sparcc_ibs, ibs_corr, features)

    summary = {
        "layer": args.layer_name,
        "samples_total": int(len(groups)),
        "hc_samples": int(np.sum(groups == args.hc_label)),
        "ibs_samples": int(np.sum(groups == args.ibs_label)),
        "features_retained": int(len(features)),
        "primary_prevalence": args.prevalence,
        "primary_threshold": args.threshold,
        "permutations": args.permutations,
        "equal_size_iterations": args.equal_size_iterations,
        "balanced_bootstraps": args.balanced_bootstraps,
        "edge_bootstraps": args.edge_bootstraps,
        "sparcc_concordance": concordance,
    }
    (args.output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2))
    pd.Series(features, name="feature").to_csv(args.output_dir / "retained_features.csv", index=False)


if __name__ == "__main__":
    main()
