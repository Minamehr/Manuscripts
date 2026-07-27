# Individual assembly and gene-catalogue construction

All steps in this section were performed in Galaxy.

## Input

Host-filtered paired FASTQ files for each of the 42 samples.

## 1. Individual metagenomic assembly

**Tool:** metaSPAdes
**Design:** one assembly per sample
**Input:** host-filtered paired reads
**Output:** one contig assembly per sample

## 2. Gene prediction

**Tool:** Prodigal
**Mode:** metagenomic
**Input:** the individual sample assemblies
**Output:** predicted nucleotide genes, predicted proteins, and gene-coordinate files

Gene prediction was performed separately for each individual assembly.

## 3. Sequence-length filtering

Predicted nucleotide genes shorter than **100 bp** were removed before catalogue clustering. The corresponding predicted proteins were retained for the genes passing this filter.

## 4. Nonredundant catalogue construction

**Tool:** MMseqs2
**Input:** combined predicted proteins from all samples
**Minimum amino-acid sequence identity:** 0.95
**Minimum alignment coverage:** 0.90
**Coverage mode:** coverage relative to the shorter sequence
**Output:** representative proteins and membership assignments

The representative protein identifiers were used to recover the matching nucleotide gene sequences.

## Final catalogue

The final catalogue contained **934,495 nonredundant genes**.

## Downstream processing

Taxonomic annotation, functional annotation, read mapping, Samsum abundance estimation, and FPKM output were performed outside Galaxy and are documented in:

`../script/Gene Annotation and Abundance`
