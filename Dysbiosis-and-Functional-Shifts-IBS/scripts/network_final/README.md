# Final microbial association analysis

`run_final_network_analysis.py` is a cleaned reproducibility implementation of the locked HC-versus-IBS analysis. It does not reproduce the superseded hub, centrality, fragmentation, or subtype-network claims.

## Primary locked settings

- final cohort: HC n=17, IBS n=25;
- exclude `IBSC02` when present in older input tables;
- shared prevalence >=50% in both groups;
- sample-wise multiplicative zero replacement;
- closure and CLR transformation;
- Pearson correlation as the primary estimator;
- primary edge threshold: |r| >=0.40;
- threshold sensitivity: 0.35, 0.40, 0.50, 0.60;
- 1,000 group-label permutations;
- 500 equal-size comparisons;
- 500 balanced bootstraps;
- 500 edge-stability bootstraps;
- prevalence, abundance-quantile, Spearman, and optional SparCC sensitivity checks.

Locked retained feature counts:

- species: 676;
- strains: 132;
- KEGG-annotated taxa: 64.

## Example

```bash
python run_final_network_analysis.py   --abundance Species_Level_Data.csv   --metadata metadata_for_analysis.csv   --output-dir results/species   --layer-name species   --orientation features-by-samples   --sample-column SampleID   --group-column GeneralGroup   --exclude-sample IBSC02   --feature-metadata species_domains.csv
```

Run the script separately for species, strain, and KEGG-annotated taxon layers.
