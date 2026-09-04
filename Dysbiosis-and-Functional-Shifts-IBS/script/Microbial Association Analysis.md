# Microbial association analysis

This analysis compares microbial association structure between healthy controls (HC) and participants with irritable bowel syndrome (IBS).

## Input format

### Abundance table

The first column contains feature identifiers. The remaining columns contain sample abundances.

Default orientation:

```text
feature,S01,S02,S03,...
taxon_1,0.10,0.00,0.03,...
taxon_2,0.00,0.02,0.01,...
```

Use `--orientation samples-by-features` when samples are rows instead.

### Metadata

The metadata table must contain:

```text
SampleID,GeneralGroup
S01,HC
S02,IBS
```

Default column names:

- sample identifier: `SampleID`
- comparison group: `GeneralGroup`
- group labels: `HC` and `IBS`

### Optional feature-domain metadata

For cross-domain summaries:

```text
feature,domain
taxon_1,Bacteria
taxon_2,Viruses
```

### Optional SparCC matrices

Previously generated HC and IBS SparCC correlation matrices may be supplied with:

```text
--sparcc-hc HC_sparcc_matrix.tsv
--sparcc-ibs IBS_sparcc_matrix.tsv
```

The script then reports matrix concordance with the CLR-Pearson results.

## Primary analysis

The primary settings are:

- shared prevalence of at least 50% in both HC and IBS;
- sample-wise multiplicative zero replacement;
- closure and centred log-ratio transformation;
- Pearson correlations;
- edge-density summaries at `|r| >= 0.40`;
- 1,000 group-label permutations;
- 500 equal-size HC-versus-IBS comparisons;
- 500 balanced bootstrap comparisons;
- 500 edge-stability bootstraps;
- Benjamini-Hochberg correction of permutation p-values.

Sensitivity analyses include:

- alternative prevalence thresholds;
- correlation thresholds of 0.35, 0.40, 0.50, and 0.60;
- Spearman correlations;
- feature-wise half-minimum zero replacement;
- low- and high-abundance feature strata;
- optional CLR-Pearson versus SparCC matrix concordance.

The primary retained feature sets in the study contained:

- 676 species-level features;
- 132 strain-level features;
- 64 KEGG-annotated taxon features.

## Species-level run

```bash
python script/microbial_association_analysis.py \
  --abundance /path/to/Species_Level_Data.csv \
  --metadata /path/to/metadata_for_analysis.csv \
  --output-dir network_results/species \
  --layer-name species \
  --orientation features-by-samples \
  --sample-column SampleID \
  --group-column GeneralGroup
```

When the abundance table still contains the previously excluded sample, add:

```text
--exclude-sample IBSC02
```

When a domain mapping is available, add:

```text
--feature-metadata /path/to/species_domains.tsv \
--feature-column feature \
--domain-column domain
```

## Strain-level run

```bash
python script/microbial_association_analysis.py \
  --abundance /path/to/Strain_Level_Data.csv \
  --metadata /path/to/metadata_for_analysis.csv \
  --output-dir network_results/strain \
  --layer-name strain \
  --orientation features-by-samples \
  --sample-column SampleID \
  --group-column GeneralGroup
```

## KEGG-annotated taxon run

```bash
python script/microbial_association_analysis.py \
  --abundance /path/to/KEGG_annotated_taxon_abundance.csv \
  --metadata /path/to/metadata_for_analysis.csv \
  --output-dir network_results/kegg_taxa \
  --layer-name kegg_taxa \
  --orientation features-by-samples \
  --sample-column SampleID \
  --group-column GeneralGroup
```

## Main outputs

For each layer, the script writes:

- `descriptive_metrics.csv`
- `primary_permutation_tests.csv`
- `equal_size_iterations.csv`
- `balanced_bootstrap_iterations.csv`
- `resampling_summary.csv`
- `edge_stability_hc.csv.gz`
- `edge_stability_ibs.csv.gz`
- `sensitivity_analysis.csv`
- `pearson_hc.csv`
- `pearson_ibs.csv`
- `run_summary.json`
- `retained_features.csv`

The analysis evaluates group-level differences in overall association strength, edge density, positive and negative edge density, matrix-wide correlation differences, and-when feature-domain information is supplied-cross-domain association structure.
