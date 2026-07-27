# Final MAG-level metabolic enrichment export

`export_mag_metabolic_enrichment.py` reads the direct output of `anvi-compute-metabolic-enrichment` and exports all results plus modules passing BH-adjusted q < 0.05.

It deliberately does not apply:

- `enrichment_score > 5`;
- `N_IBS >= 50`;
- any other undocumented filtering rule.

The locked analysis used 221 group-specific bin/MAG entries:

- HC: 94;
- IBS: 127.

The locked HC-versus-IBS analysis contained 24 modules passing q < 0.05.

## Example

```bash
python export_mag_metabolic_enrichment.py   --input metabolic_enrichment.tsv   --output-dir results/mag_metabolism   --n-hc-column N_HC   --n-ibs-column N_IBS
```
