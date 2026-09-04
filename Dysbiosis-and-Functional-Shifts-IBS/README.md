# Shotgun metagenomic analysis of the gut microbiome in irritable bowel syndrome

This repository documents the final analysis workflow used for 42 gut metagenomes: 17 healthy controls and 25 participants with irritable bowel syndrome, including IBS-C, IBS-D, and IBS-M.

Raw sequencing data are available from the European Nucleotide Archive under accession **PRJEB104707**.

## Repository organisation

### `galaxy/`

This directory documents analyses performed in Galaxy:

1. `preprocessing_and_taxonomy.md`
   - paired-end quality trimming;
   - removal of human reads;
   - Kraken2/Bracken taxonomic profiling;
   - MetaPhlAn taxonomic profiling.

2. `gene_catalogue.md`
   - individual metagenomic assembly;
   - gene prediction;
   - removal of genes shorter than 100 bp;
   - MMseqs2 construction of the nonredundant gene catalogue.

3. `mag_recovery.md`
   - group-wise co-assembly for HC, IBS-C, IBS-D, and IBS-M;
   - MetaBAT2, MaxBin2, and CONCOCT binning;
   - metaWRAP bin refinement;
   - CheckM quality assessment.

### `script/`

This directory contains analyses performed outside Galaxy:

1. `Gene Annotation and Abundance`
   - UniProt TrEMBL taxonomic annotation;
   - KOBAS, eggNOG-mapper, dbCAN/HMMER, Pfam, and PathoFact annotation;
   - BWA-MEM mapping to the nucleotide gene catalogue;
   - Samsum abundance estimation and FPKM output.

2. `MAG Downstream Analysis`
   - anvi'o database construction and genome assessment;
   - COG, KOfam, Pfam, and CAZyme annotation;
   - final cross-group dereplication at 99% ANI;
   - GTDB-Tk taxonomy;
   - abundance/recruitment summaries and phylogenomics.

3. `MAG Metabolic Enrichment`
   - anvi'o metabolic-module estimation;
   - HC-versus-IBS metabolic enrichment.

4. `Strain-Level Haplotyping`
   - anvi'o nucleotide-variation profiling;
   - DESMAN haplotype inference and model selection.

4. `microbial_association_analysis.py` and `Microbial Association Analysis.md`
   - final species-, strain-, and KEGG-annotated taxon association analyses;
   - gene-catalogue-derived input matrices;
   - shared-prevalence filtering, CLR-Pearson primary analysis, permutation/resampling/bootstrap tests, and SparCC sensitivity analysis.

Local database paths, Galaxy history links, and sample-specific paths must be adapted before applying the workflow to another dataset.

## Main outputs

- Nonredundant gene catalogue: **934,495 genes**.
- Nonredundant quality-filtered MAG set: **154 MAGs**.
- MAG metabolic enrichment input: **221 group-specific genome/bin entries** comprising 94 HC and 127 IBS entries.
