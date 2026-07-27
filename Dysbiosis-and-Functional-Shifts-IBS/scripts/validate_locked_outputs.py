#!/usr/bin/env python3
"""Validate key results against the locked manuscript values."""

import argparse, json
from pathlib import Path


def load(path):
    return json.loads(Path(path).read_text())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--species-summary", type=Path)
    p.add_argument("--strain-summary", type=Path)
    p.add_argument("--kegg-summary", type=Path)
    p.add_argument("--functional-summary", action="append", default=[], type=Path)
    p.add_argument("--mag-summary", type=Path)
    args = p.parse_args()
    checks = []
    for name, path, expected in [
        ("species_features", args.species_summary, 676),
        ("strain_features", args.strain_summary, 132),
        ("kegg_taxa_features", args.kegg_summary, 64),
    ]:
        if path:
            observed = load(path).get("features_retained")
            checks.append({"check": name, "expected": expected, "observed": observed,
                           "pass": observed == expected})
    for path in args.functional_summary:
        observed = load(path).get("hc_vs_ibs_features_q_lt_0.05")
        checks.append({"check": f"functional_zero_fdr::{path.name}",
                       "expected": 0, "observed": observed, "pass": observed == 0})
    if args.mag_summary:
        observed = load(args.mag_summary).get("significant_modules_q_lt_0.05")
        checks.append({"check": "mag_significant_modules", "expected": 24,
                       "observed": observed, "pass": observed == 24})
    print(json.dumps({"all_pass": all(x["pass"] for x in checks), "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
