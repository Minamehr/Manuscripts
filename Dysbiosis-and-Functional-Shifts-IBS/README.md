# Metagenomic Insights into Dysbiosis and Functional Shifts in Irritable Bowel Syndrome

This repository contains scripts, workflows, and metadata associated with the analysis presented in the manuscript:  
**“Metagenomic Insights into Dysbiosis and Functional Shifts in Irritable Bowel Syndrome”**

The analysis explores taxonomic, functional, and genomic shifts in the gut microbiome of individuals with IBS (including subtypes IBS-C, IBS-D, and IBS-M) compared to healthy controls, using whole metagenome shotgun sequencing.

---

## Repository Overview

This repository is organized into several parts reflecting key stages of the analysis:

### **Part 1: Preprocessing**
- **Quality control** of raw reads using TrimGalore
- **Host read removal** using BBMap against the GRCh38.p13 human reference genome
- **Initial taxonomic profiling** with:
  - Kraken2 + Bracken
  - MetaPhlAn (for validation)

### **Part 2: Assembly and Gene Catalog Construction**
- **Individual metagenomic assemblies** using MetaSPAdes
- **Gene prediction** with Prodigal (discarding genes <100 bp)
- **Clustering of predicted genes** into a non-redundant gene catalog using MMseqs2 (95% identity, 90% overlap)
- **Taxonomic annotation** via UniProt TrEMBL
- **Functional annotation** with:
  - KEGG (via KOBAS)
  - COG (via eggNOG-mapper)
  - CAZymes (via dbCAN + HMMER)
  - ARGs and virulence factors (via PathoFact)
- **Gene abundance estimation** using BWA-MEM and samsum, normalized as FPKM

### **Part 3: Genome Reconstruction**
- **Group-wise co-assemblies** (IBS-C, IBS-D, IBS-M, HC)
- **MAG binning** with MetaBAT2, MaxBin2, and CONCOCT
- **Bin refinement** via metaWRAP (Bin_refinement module)
- **Quality assessment** with CheckM and anvi’o (≥70% completeness, ≤5% redundancy)
- **Dereplication** using `anvi-dereplicate-genomes`
- **Taxonomy assignment** using GTDB-Tk
- **Functional annotation** of MAGs with anvi’o (COG, KEGG KO, Pfam, CAZy)

### **Part 4: Strain-Level Haplotyping**
- **SNV calling and haplotype inference** using anvi’o and DESMAN
- **Optimal haplotype number** based on posterior deviance and SNV uncertainty (<10%)

### **Part 5: Analyses**
- **Functional Enrichment** ausing `anvi-estimate-metabolism` and `anvi-compute-metabolic-enrichment`
- **Microbial network analysis** with:
  - **SparCC** for correlation inference (|r| ≥ 0.5)
  - **NetworkX** for centrality and connectivity metrics
  - **Cytoscape** for visualization and modularity detection


