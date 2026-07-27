# Group-wise co-assembly and MAG recovery

All assembly, binning, refinement, and initial quality-assessment steps in this section were performed in Galaxy.

## Input groups

Host-filtered paired reads were analysed separately for:

- HC;
- IBS-C;
- IBS-D;
- IBS-M.

## 1. Group-wise co-assembly

**Tool:** metaSPAdes
**Design:** one co-assembly for each of the four study groups
**Output:** HC, IBS-C, IBS-D, and IBS-M co-assemblies

## 2. Coverage profiles

Reads from each sample were mapped to the appropriate group co-assembly to generate the contig-coverage information required for binning.

## 3. Independent genome binning

Each group co-assembly was binned using:

- MetaBAT2;
- MaxBin2;
- CONCOCT.

## 4. Bin refinement

**Tool:** metaWRAP Bin Refinement
**Minimum completeness:** 70%
**Maximum contamination:** 5%

The three independent bin sets were supplied to metaWRAP for refinement.

## 5. Initial genome-quality assessment

**Tool:** CheckM
**Input:** refined group-specific bins
**Output:** completeness and contamination estimates

## Downstream processing

The refined bins were exported from Galaxy for anvi'o processing, 99% ANI dereplication, GTDB-Tk classification, genome annotation, abundance/recruitment analysis, phylogenomics, metabolic enrichment, and DESMAN haplotyping.

These steps are documented in:

- `../script/MAG Downstream Analysis`
- `../script/MAG Metabolic Enrichment`
- `../script/Strain-Level Haplotyping`
