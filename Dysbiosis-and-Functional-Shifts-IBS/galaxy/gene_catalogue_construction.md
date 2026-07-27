# Galaxy record: metagenomic assembly and gene-catalogue construction

## Galaxy-executed steps

Metagenomic assembly, gene prediction, sequence-length filtering, and sequence clustering for gene-catalogue construction were performed using versioned tools implemented within Galaxy.

The Galaxy-based steps were:

1. Per-sample metagenomic assembly using metaSPAdes.
2. Gene prediction from assembled contigs using Prodigal in metagenomic mode (`-p meta`).
3. Removal of genes shorter than the threshold reported in the final manuscript Methods.
4. Clustering of retained genes using MMseqs2 at:
   - >=95% sequence identity;
   - >=90% alignment coverage.

The resulting nonredundant catalogue contained 934,495 genes.

## Downstream analyses performed separately

The following downstream steps should not be implied to have been performed in Galaxy unless confirmed by the original histories:

- taxonomic assignment against UniProt TrEMBL;
- KEGG annotation using KOBAS;
- COG annotation using eggNOG-mapper;
- CAZyme annotation using dbCAN and HMMER;
- protein-domain annotation using Pfam;
- ARG and virulence-factor identification using PathoFact;
- read alignment to the gene catalogue using BWA-MEM;
- read counting using Samsum;
- FPKM normalization;
- sample-level functional aggregation and statistical testing.

## Reproducibility information to add

From the original Galaxy history, add:

- metaSPAdes wrapper and software versions;
- Prodigal wrapper and software versions;
- MMseqs2 wrapper and software versions;
- the exact sequence-length filtering implementation;
- all nondefault parameters;
- input and output dataset identifiers;
- exported history or workflow link, when available.

## Scope note

This file documents the Galaxy-executed portion of gene-catalogue construction. Downstream annotation and abundance scripts belong under `scripts/gene_function_final/` only after validation against the final locked results.
