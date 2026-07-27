# Galaxy record: short-read preprocessing and taxonomic profiling

## Execution environment

- Platform: Galaxy
- Public instance: UseGalaxy.eu
- Published taxonomy workflow:
  https://usegalaxy.eu/published/workflow?id=7491883694fff308

## Analysis steps

1. Raw paired-end reads were quality-filtered and adapters were removed using Trim Galore v0.6.6.
2. Human-derived reads were removed by alignment against GRCh38.p13 using BBMap.
3. Host-filtered reads were taxonomically classified using Kraken2.
4. Bracken was used for abundance estimation.
5. MetaPhlAn was used for complementary taxonomic profiling and cross-method comparison.

## Reproducibility information to add

For each Galaxy tool, add from the original Galaxy history:

- Galaxy wrapper version;
- underlying software version;
- complete parameter settings;
- input dataset identifiers;
- output dataset identifiers;
- database build or release;
- exported workflow or history link, when available.

## Scope note

This file documents analyses actually performed within Galaxy. It is not a reconstructed command-line substitute.
