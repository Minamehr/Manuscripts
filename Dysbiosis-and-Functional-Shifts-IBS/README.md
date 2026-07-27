# Shotgun Metagenomic Analysis of the Gut Microbiome in Irritable Bowel Syndrome

This repository contains documentation, cleaned reproducibility scripts, supporting metadata, and archived exploratory code associated with a shotgun metagenomic study comparing patients with irritable bowel syndrome (IBS) and healthy controls.

The final study cohort comprised 42 metagenomes:

- 17 healthy controls
- 25 patients with IBS
- IBS-C, IBS-D, and IBS-M subtype groups

Raw sequencing data are available from the European Nucleotide Archive under accession **PRJEB104707**.

## Repository principles

This repository distinguishes three types of material:

1. **Galaxy-executed analyses**
   These are documented using the Galaxy instance, tool names, tool versions, parameters, inputs, outputs, and public workflow or history links where available. They are not represented by invented command-line scripts.

2. **Cleaned reproducibility implementations**
   These scripts implement the final analyses reported in the revised manuscript and are validated against the final locked results.

3. **Archived exploratory analyses**
   Earlier scripts are retained for provenance but are not used to support the final manuscript conclusions.

## Analysis overview

### Short-read preprocessing and taxonomic profiling

Read preprocessing and read-based taxonomic profiling were performed using versioned tools on the European Galaxy server, UseGalaxy.eu.

The main steps included:

- adapter and quality trimming with Trim Galore;
- removal of human-derived reads using BBMap against GRCh38.p13;
- taxonomic classification with Kraken2;
- abundance estimation with Bracken;
- complementary taxonomic profiling with MetaPhlAn.

The initial version of the Galaxy Metagenomic Taxonomy Analysis workflow is publicly available at:

https://usegalaxy.eu/published/workflow?id=7491883694fff308

See `galaxy/preprocessing_and_taxonomy.md`.

### Gene-catalogue construction

Metagenomic assembly, gene prediction, sequence-length filtering, and sequence clustering for gene-catalogue construction were performed using versioned tools implemented within Galaxy.

The Galaxy-based steps included:

- individual assembly using metaSPAdes;
- gene prediction using Prodigal in metagenomic mode;
- removal of genes shorter than the threshold reported in the final Methods;
- sequence clustering using MMseqs2 at >=95% sequence identity and >=90% alignment coverage.

Downstream taxonomic and functional annotation and gene-abundance estimation were performed using the software and databases described in the manuscript, including UniProt TrEMBL, KOBAS, eggNOG-mapper, dbCAN/HMMER, Pfam, PathoFact, BWA-MEM, Samsum, and FPKM normalization.

The final nonredundant gene catalogue contained **934,495 genes**.

See `galaxy/gene_catalogue_construction.md`.

### Genome recovery and downstream MAG analysis

Gut microbial genomes were reconstructed using group-wise co-assembly and multiple binning tools implemented within Galaxy.

The Galaxy-based recovery steps included:

- separate co-assemblies for HC, IBS-C, IBS-D, and IBS-M using metaSPAdes;
- binning with MetaBAT2, MaxBin2, and CONCOCT;
- bin refinement with the metaWRAP Bin Refinement module;
- initial completeness and contamination assessment with CheckM.

The resulting bins were exported for downstream analysis using anvi'o, GTDB-Tk, MUSCLE, DESMAN, and associated statistical procedures.

The final analysis retained **154 nonredundant quality-filtered MAGs**.

See `galaxy/mag_recovery.md`.

### Final network analysis

The final microbial association analysis used:

- shared-prevalence filtering;
- multiplicative zero replacement, closure, and CLR transformation;
- Pearson correlation as the primary estimator;
- label permutation testing;
- equal-size HC-versus-IBS comparisons;
- balanced bootstrap assessment;
- edge-stability bootstrapping;
- BH multiple-testing correction;
- prevalence, threshold, Spearman, and SparCC sensitivity analyses.

The final analysis did not provide robust evidence of widespread microbial-network fragmentation, altered cross-domain organization, or global functional-network restructuring in IBS.

Earlier SparCC hub, centrality, fragmentation, and subtype-network scripts are retained only under `archive/`.

## Directory structure

```text
galaxy/
    preprocessing_and_taxonomy.md
    gene_catalogue_construction.md
    mag_recovery.md

scripts/
    diversity_and_taxonomy/
    gene_function_final/
    mag_analysis_final/
    network_final/
    figures/

archive/
    legacy_command_examples/
    legacy_network_analysis/
```

## Galaxy citations

Analyses performed within Galaxy should cite:

- The Galaxy Community. Galaxy for accessible, reproducible, and collaborative data analyses: 2026 update. *Nucleic Acids Research*. 2026;54(W1):W105-W116. doi:10.1093/nar/gkag469.
- Blankenberg D, Von Kuster G, Bouvier E, Baker D, Afgan E, Stoler N, et al. Dissemination of scientific software with Galaxy ToolShed. *Genome Biology*. 2014;15:403. doi:10.1186/gb4161.

The original publication for every individual scientific tool should also be cited.

## Repository status

The Galaxy documentation is complete at the level currently supported by the manuscript record. Exact Galaxy wrapper versions and exported histories should be added where available.

Cleaned reproducibility implementations are included for the final network, gene-function, and MAG metabolic-enrichment analyses. These scripts encode the locked analytical logic reported in the revised manuscript but are not presented as byte-for-byte copies of the original historical scripts. Archived exploratory scripts should not be interpreted as final analytical implementations.
