#!/usr/bin/env python3
"""Clean reproducibility implementation of the final HC-versus-IBS network analysis."""

from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, pearsonr


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
    return df.T if orientation == "features-by-samples" else df


def mult_replace(row, delta=None):
    x = np.asarray(row, float)
    if np.any(x < 0) or x.sum() <= 0:
        raise ValueError("Each sample must contain non-negative values and at least one positive value.")
    x = x / x.sum()
    zero = x == 0
    z = int(zero.sum())
    if z == 0:
        return x
    d = len(x)
    delta = 1 / (d * d) if delta is None else float(delta)
    if z * delta >= 1:
        raise ValueError("Zero replacement delta is too large.")
    x[zero] = delta
    x[~zero] *= (1 - z * delta) / x[~zero].sum()
    return x


def clr(matrix, delta=None):
    replaced = np.vstack([mult_replace(row, delta) for row in matrix])
    logged = np.log(replaced)
    return logged - logged.mean(axis=1, keepdims=True)


def corr(x, estimator):
    if estimator == "spearman":
        x = np.apply_along_axis(rankdata, 0, x)
    out = np.corrcoef(x, rowvar=False)
    out = np.nan_to_num(out, nan=0.0)
    np.fill_diagonal(out, 1.0)
    return out


def upper(c):
    tri = np.triu_indices_from(c, k=1)
    return c[tri], tri


def cross_mask(features, domain_map, tri):
    if not domain_map:
        return None
    domains = np.array([domain_map.get(f, "Unknown") for f in features], object)
    return domains[tri[0]] != domains[tri[1]]


def metrics(c, threshold, cmask=None):
    vals, _ = upper(c)
    out = {
        "mean_abs_correlation": float(np.mean(np.abs(vals))),
        "edge_density": float(np.mean(np.abs(vals) >= threshold)),
        "positive_edge_density": float(np.mean(vals >= threshold)),
        "negative_edge_density": float(np.mean(vals <= -threshold)),
    }
    if cmask is not None and cmask.any():
        cross = vals[cmask]
        out["cross_domain_mean_abs_correlation"] = float(np.mean(np.abs(cross)))
        out["cross_domain_edge_density"] = float(np.mean(np.abs(cross) >= threshold))
    return out


def differences(ch, ci, threshold, cmask=None):
    mh, mi = metrics(ch, threshold, cmask), metrics(ci, threshold, cmask)
    out = {k: mh[k] - mi[k] for k in mh}
    vh, _ = upper(ch)
    vi, _ = upper(ci)
    out["global_correlation_matrix_distance"] = float(np.linalg.norm(vh - vi))
    return out


def prepare(raw, metadata, args, prevalence, abundance_quantile):
    raw.index = raw.index.astype(str)
    metadata = metadata.copy()
    metadata[args.sample_column] = metadata[args.sample_column].astype(str)
    metadata = metadata.set_index(args.sample_column, drop=False)
    samples = [s for s in metadata.index if s in raw.index and s not in set(args.exclude_sample)]
    x = raw.loc[samples]
    m = metadata.loc[samples]
    valid = m[args.group_column].isin([args.hc_label, args.ibs_label])
    x, m = x.loc[valid], m.loc[valid]
    g = m[args.group_column].astype(str)
    hc = x.loc[g == args.hc_label]
    ibs = x.loc[g == args.ibs_label]
    keep = (hc.gt(0).mean() >= prevalence) & (ibs.gt(0).mean() >= prevalence)
    x = x.loc[:, keep]
    if abundance_quantile > 0:
        means = x.mean()
        x = x.loc[:, means >= means.quantile(abundance_quantile)]
    if x.shape[1] < 2:
        raise ValueError("Fewer than two features remain after filtering.")
    return x, g, list(x.columns)


def group_corrs(z, groups, args, estimator):
    return corr(z[groups == args.hc_label], estimator), corr(z[groups == args.ibs_label], estimator)


def permutation(z, groups, args, cmask, rng):
    ch, ci = group_corrs(z, groups, args, "pearson")
    obs = differences(ch, ci, args.threshold, cmask)
    null = {k: [] for k in obs}
    for _ in range(args.permutations):
        pg = rng.permutation(groups)
        ph, pi = group_corrs(z, pg, args, "pearson")
        d = differences(ph, pi, args.threshold, cmask)
        for k, v in d.items():
            null[k].append(v)
    rows = []
    for k, observed in obs.items():
        values = np.asarray(null[k])
        if k == "global_correlation_matrix_distance":
            p = (1 + np.sum(values >= observed)) / (len(values) + 1)
        else:
            p = (1 + np.sum(np.abs(values) >= abs(observed))) / (len(values) + 1)
        rows.append({"metric": k, "observed_hc_minus_ibs": observed, "permutation_p": p,
                     "null_mean": values.mean(), "null_sd": values.std(ddof=1)})
    out = pd.DataFrame(rows)
    out["bh_q"] = bh(out["permutation_p"])
    return out


def resample(z, groups, args, cmask, rng, iterations, replace):
    hi = np.flatnonzero(groups == args.hc_label)
    ii = np.flatnonzero(groups == args.ibs_label)
    n = min(len(hi), len(ii))
    rows = []
    for i in range(iterations):
        h = rng.choice(hi, n, replace=replace)
        b = rng.choice(ii, n, replace=replace)
        d = differences(corr(z[h], "pearson"), corr(z[b], "pearson"), args.threshold, cmask)
        d["iteration"] = i + 1
        rows.append(d)
    return pd.DataFrame(rows)


def resample_summary(df, name):
    rows = []
    for col in df.columns:
        if col == "iteration":
            continue
        v = df[col].to_numpy(float)
        rows.append({"analysis": name, "metric": col, "mean": v.mean(), "median": np.median(v),
                     "ci_low_2.5": np.quantile(v, .025), "ci_high_97.5": np.quantile(v, .975),
                     "fraction_positive": np.mean(v > 0)})
    return pd.DataFrame(rows)


def edge_stability(zgroup, features, args, rng):
    p = len(features)
    tri = np.triu_indices(p, k=1)
    nedge = len(tri[0])
    boot = np.empty((args.edge_stability_iterations, nedge), dtype=np.float32)
    for i in range(args.edge_stability_iterations):
        idx = rng.choice(np.arange(len(zgroup)), len(zgroup), replace=True)
        boot[i] = corr(zgroup[idx], "pearson")[tri]
    mean = boot.mean(0)
    low = np.quantile(boot, .025, axis=0)
    high = np.quantile(boot, .975, axis=0)
    sign = np.maximum((boot > 0).mean(0), (boot < 0).mean(0))
    thresh = (np.abs(boot) >= args.threshold).mean(0)
    return pd.DataFrame({
        "feature_1": np.asarray(features)[tri[0]],
        "feature_2": np.asarray(features)[tri[1]],
        "bootstrap_mean_r": mean,
        "ci_low_2.5": low,
        "ci_high_97.5": high,
        "sign_stability": sign,
        "threshold_stability": thresh,
    })


def concordance(path, observed, features):
    ref = pd.read_csv(path, sep=sep(path), index_col=0)
    common = [f for f in features if f in ref.index and f in ref.columns]
    idx = [features.index(f) for f in common]
    o = observed[np.ix_(idx, idx)]
    r = ref.loc[common, common].to_numpy(float)
    tri = np.triu_indices(len(common), k=1)
    value, p = pearsonr(o[tri], r[tri])
    return {"n_common_features": len(common), "matrix_pearson_r": value, "p_value": p}


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
    p.add_argument("--prevalence", type=float, default=.50)
    p.add_argument("--abundance-quantile", type=float, default=0)
    p.add_argument("--threshold", type=float, default=.40)
    p.add_argument("--thresholds", default=".35,.40,.50,.60")
    p.add_argument("--prevalence-sensitivity", default=".50,.60,.70,.80")
    p.add_argument("--abundance-quantile-sensitivity", default="0,.25,.50,.75")
    p.add_argument("--zero-delta", type=float)
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--equal-size-iterations", type=int, default=500)
    p.add_argument("--balanced-bootstrap-iterations", type=int, default=500)
    p.add_argument("--edge-stability-iterations", type=int, default=500)
    p.add_argument("--seed", type=int, default=20260722)
    p.add_argument("--sparcc-hc", type=Path)
    p.add_argument("--sparcc-ibs", type=Path)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    raw = pd.read_csv(args.abundance, sep=sep(args.abundance), index_col=0)
    raw = orient(raw, args.orientation)
    meta = pd.read_csv(args.metadata, sep=sep(args.metadata))
    x, g, features = prepare(raw.copy(), meta, args, args.prevalence, args.abundance_quantile)
    groups = g.to_numpy(str)
    z = clr(x.to_numpy(float), args.zero_delta)
    ch, ci = group_corrs(z, groups, args, "pearson")
    _, tri = upper(ch)

    dmap = None
    if args.feature_metadata:
        fm = pd.read_csv(args.feature_metadata, sep=sep(args.feature_metadata))
        dmap = dict(zip(fm[args.feature_column].astype(str), fm[args.domain_column].astype(str)))
    cmask = cross_mask(features, dmap, tri)

    desc = []
    for label, matrix in [(args.hc_label, ch), (args.ibs_label, ci)]:
        for threshold in map(float, args.thresholds.split(",")):
            row = metrics(matrix, threshold, cmask)
            row.update({"layer": args.layer_name, "group": label, "threshold": threshold, "n_features": len(features)})
            desc.append(row)
    pd.DataFrame(desc).to_csv(args.output_dir / "descriptive_metrics.csv", index=False)

    permutation(z, groups, args, cmask, rng).to_csv(args.output_dir / "permutation_tests.csv", index=False)
    eq = resample(z, groups, args, cmask, rng, args.equal_size_iterations, False)
    bs = resample(z, groups, args, cmask, rng, args.balanced_bootstrap_iterations, True)
    eq.to_csv(args.output_dir / "equal_size_iterations.csv", index=False)
    bs.to_csv(args.output_dir / "balanced_bootstrap_iterations.csv", index=False)
    pd.concat([resample_summary(eq, "equal_size"), resample_summary(bs, "balanced_bootstrap")]).to_csv(
        args.output_dir / "resampling_summary.csv", index=False)

    if args.edge_stability_iterations:
        edge_stability(z[groups == args.hc_label], features, args, rng).to_csv(
            args.output_dir / "edge_stability_hc.csv.gz", index=False, compression="gzip")
        edge_stability(z[groups == args.ibs_label], features, args, rng).to_csv(
            args.output_dir / "edge_stability_ibs.csv.gz", index=False, compression="gzip")

    sens = []
    for prev in map(float, args.prevalence_sensitivity.split(",")):
        for aq in map(float, args.abundance_quantile_sensitivity.split(",")):
            try:
                sx, sg, sf = prepare(raw.copy(), meta, args, prev, aq)
            except ValueError:
                continue
            sz = clr(sx.to_numpy(float), args.zero_delta)
            sgroups = sg.to_numpy(str)
            for est in ("pearson", "spearman"):
                sh, si = group_corrs(sz, sgroups, args, est)
                for threshold in map(float, args.thresholds.split(",")):
                    for metric, value in differences(sh, si, threshold).items():
                        sens.append({"prevalence": prev, "abundance_quantile": aq, "estimator": est,
                                     "threshold": threshold, "n_features": len(sf), "metric": metric,
                                     "hc_minus_ibs": value})
    pd.DataFrame(sens).to_csv(args.output_dir / "sensitivity_analysis.csv", index=False)

    sp_h, sp_i = group_corrs(z, groups, args, "spearman")
    pd.DataFrame(ch, index=features, columns=features).to_csv(args.output_dir / "pearson_hc.csv")
    pd.DataFrame(ci, index=features, columns=features).to_csv(args.output_dir / "pearson_ibs.csv")
    pd.DataFrame(sp_h, index=features, columns=features).to_csv(args.output_dir / "spearman_hc.csv")
    pd.DataFrame(sp_i, index=features, columns=features).to_csv(args.output_dir / "spearman_ibs.csv")

    sc = {}
    if args.sparcc_hc and args.sparcc_ibs:
        sc["HC"] = concordance(args.sparcc_hc, ch, features)
        sc["IBS"] = concordance(args.sparcc_ibs, ci, features)

    summary = {"layer": args.layer_name, "samples_total": len(groups),
               "hc_samples": int(np.sum(groups == args.hc_label)),
               "ibs_samples": int(np.sum(groups == args.ibs_label)),
               "features_retained": len(features), "shared_prevalence": args.prevalence,
               "primary_threshold": args.threshold, "permutations": args.permutations,
               "equal_size_iterations": args.equal_size_iterations,
               "balanced_bootstrap_iterations": args.balanced_bootstrap_iterations,
               "edge_stability_iterations": args.edge_stability_iterations,
               "seed": args.seed, "sparcc_concordance": sc}
    (args.output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2))
    pd.Series(features, name="feature").to_csv(args.output_dir / "retained_features.csv", index=False)


if __name__ == "__main__":
    main()
