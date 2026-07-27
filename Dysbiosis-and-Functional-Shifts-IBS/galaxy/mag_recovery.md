# Galaxy record: group-wise co-assembly and initial MAG recovery

## Galaxy-executed steps

Gut microbial genomes from healthy controls and patients with IBS were reconstructed using group-wise co-assembly and multiple binning tools implemented within Galaxy.

The Galaxy-based recovery steps were:

1. Separate co-assembly of clean reads from:
   - HC;
   - IBS-C;
   - IBS-D;
   - IBS-M.
2. Co-assembly using metaSPAdes.
3. Genome binning using:
   - MetaBAT2;
   - MaxBin2;
   - CONCOCT.
4. Bin refinement using the metaWRAP Bin Refinement module with:
   - minimum completeness threshold: 70%;
   - maximum contamination threshold: 5%.
5. Initial completeness and contamination assessment using CheckM.

## Downstream analyses performed separately

The resulting bins were exported for downstream processing, including:

- anvi'o-based genome assessment and manual inspection;
- assessment using 71 universal bacterial single-copy genes;
- pooling of bins from the four study groups;
- duplicate identification at >=99% ANI;
- retention of a nonredundant MAG set;
- taxonomy assignment using GTDB-Tk;
- functional annotation using COG, Pfam, KOfam, and CAZyme resources;
- metabolic reconstruction and enrichment analysis;
- phylogenomic reconstruction using anvi'o and MUSCLE;
- strain-level haplotyping using anvi'o and DESMAN.

The final analysis retained 154 nonredundant quality-filtered MAGs.

## Reproducibility information to add

From the original Galaxy histories, add:

- metaSPAdes wrapper and software versions;
- MetaBAT2 wrapper and software versions;
- MaxBin2 wrapper and software versions;
- CONCOCT wrapper and software versions;
- metaWRAP wrapper and software versions;
- CheckM wrapper and software versions;
- exact parameters and database releases;
- exported histories or workflows, when available.

## Scope note

This file documents only the Galaxy-executed recovery stage. Downstream anvi'o, GTDB-Tk, DESMAN, metabolic, statistical, and figure-generation code belongs under `scripts/mag_analysis_final/` after validation.
