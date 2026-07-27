# Final gene-function audit

`run_gene_function_audit.py` implements the final sample-level HC-versus-IBS functional audit.

It can read an already aggregated function-by-sample table or aggregate a gene-by-sample table using a gene-to-function mapping. The default aggregation is a sum across all genes assigned to the same function.

The script applies prevalence filtering, two-sided Mann-Whitney U tests, and BH correction. It does not present raw p-value-only findings as significant.

Run independently for CAZyme, KEGG, KOfam, COG, Pfam, ARG, and virulence-factor tables.

## Example

```bash
python run_gene_function_audit.py   --abundance Merged_CAZyme_Abundance.tsv   --metadata metadata_for_analysis.csv   --output-dir results/cazyme   --analysis-name CAZyme   --orientation features-by-samples   --exclude-sample IBSC02
```

The locked sample-level audit found zero BH-significant features in CAZyme, KEGG, KOfam, COG, and Pfam.
