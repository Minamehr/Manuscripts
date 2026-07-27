# Preprocessing and taxonomic profiling

All steps in this section were performed in Galaxy.

## Input

Paired-end raw shotgun metagenomic FASTQ files.

## 1. Quality trimming

**Tool:** Trim Galore
**Software version:** 0.6.6
**Mode:** paired-end
**Input:** raw R1 and R2 FASTQ files
**Output:** paired trimmed FASTQ files

Default paired-end trimming settings were used unless otherwise recorded in the Galaxy history.

## 2. Human-read removal

**Tool:** BBMap
**Software version:** 38.95
**Reference:** GRCh38.p13
**Input:** paired trimmed FASTQ files
**Output:** paired reads that did not map to the human reference

Only the paired unmapped reads were retained for downstream taxonomic, assembly, and genome-recovery analyses.

## 3. Kraken2 and Bracken profiling

**Input:** host-filtered paired FASTQ files

Kraken2 produced per-read classifications and per-sample taxonomic reports. Bracken used the Kraken2 reports to estimate species-level abundance.

The Kraken2 database selection, Bracken database, and read-length setting must match those recorded in the Galaxy workflow/history.

## 4. MetaPhlAn profiling

**Input:** host-filtered paired FASTQ files
**Output:** species- and strain-level taxonomic abundance tables

MetaPhlAn was used as a complementary marker-gene-based taxonomic profiler.

## Public Galaxy workflow

https://usegalaxy.eu/published/workflow?id=7491883694fff308

The shared workflow provides the connected Galaxy tools and workflow structure for taxonomic profiling.
