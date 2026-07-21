# TGF-β and IL-2 Differentially Shape T Follicular Regulatory Cell Differentiation and Stability In Vitro

RNA-seq differential expression analysis pipeline for the study:

> Bach L, Chang Y, Arteaga Transito O, Ghasemi M, Steinheuer LM, Steffen T, De Domenico E, Ulas T, Wunderlich FT, Beyer MD, Thurley K, et al. Baumjohann D. **TGF-β and IL-2 differentially shape T follicular regulatory cell differentiation and stability in vitro.** *Cellular & Molecular Immunology* (2026). https://doi.org/10.1038/s41423-026-01440-9

Raw sequencing data are deposited in GEO under accession [GSE306188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306188).

## Overview

Naïve CD4⁺ T cells from GREAT.Smart-17A.*FoxP3*ʰᶜᴰ² reporter mice were differentiated in vitro into Treg, Tfr, and Tfh (TGF-β) subsets, sorted by flow cytometry, and profiled by bulk RNA-seq (QuantSeq 3′ mRNA-seq, single-end, Illumina NovaSeq 6000). This repository contains the R code used to go from the gene-level count matrix to the differential expression results and figures reported in the paper: protein-coding gene filtering, DESeq2-based differential expression between Treg/Tfr/Tfh, PCA/QC, volcano plots, and overlap analysis of differentially expressed genes (DEGs) across the three pairwise comparisons.

**Scope of this repository:** the pipeline here starts from the annotated gene count matrix. Upstream processing of the raw FASTQ files (FastQC, adapter/quality trimming with Cutadapt, alignment to the *Mus musculus* GRCm39 reference genome with HISAT2, and gene-level quantification with featureCounts) was performed as described in the Methods section of the paper and is not included as executable code here.

## Repository structure

```
TGFb-IL2-Tfr-RNAseq/
├── data/
│   ├── raw_counts/       # annotated_count_data_2.xlsx (Ensembl IDs x samples, last row = group labels)
│   └── metadata/         # Curated_List.xlsx (curated gene set used to label volcano plots)
├── scripts/
│   ├── 01_gene_filtering_deseq2/   # load data, protein-coding filter, expression filter, DESeq2
│   ├── 02_pca_qc/                  # PCA, scree plot, VST expression boxplot
│   ├── 03_volcano_plots/           # plain, curated-list-labeled, and threshold-labeled volcano plots
│   └── 04_deg_overlap/             # Venn diagram + overlapping/exclusive gene lists per comparison
├── results/
│   ├── tables/    # DEG_Significant_Results.xlsx, Euler_Overlapping_Genes.xlsx, and cached .rds objects
│   └── figures/   # PNG/SVG plots, organized by parameter tag
└── docs/
```

Each `results/tables/{tag}` and `results/figures/{tag}` subfolder is named after the filtering/significance parameters used for that run (e.g. `Thr5_Frac0.75_FDR0.1_LFC0.5`), so different parameter choices don't overwrite each other.

## Requirements

- R (analysis in the paper was run with v4.4.2)
- R packages: `dplyr`, `readxl`, `ggplot2`, `ggrepel`, `ggfortify`, `ggvenn`, `biomaRt`, `DESeq2`, `openxlsx`, `svglite`
- Internet access when running script 01 (queries the Ensembl BioMart database to identify protein-coding genes)

## Input data

Place the following files under `data/`:
- `data/raw_counts/annotated_count_data_2.xlsx` — gene count matrix with Ensembl gene IDs in the first column, one column per sample, gene symbols in the last column, and a final row giving each sample's condition (`Treg` / `Tfr` / `Tfh`).
- `data/metadata/Curated_List.xlsx` — single-column list of gene symbols to highlight in the curated-list volcano plots (script 03).

## How to run

Run the scripts in order from the project root (e.g. via the `.Rproj` file, or `setwd()` to the repo root):

1. **`scripts/01_gene_filtering_deseq2/`** — loads the count matrix, restricts to protein-coding genes (via Ensembl BioMart), filters low-expression genes, fits DESeq2 across all pairwise comparisons of Treg/Tfr/Tfh, and writes the significant DEGs to Excel. Also caches the fitted `dds` object, the DEG results, and the Ensembl ID ↔ gene symbol mapping as `.rds` files so the later scripts don't need to repeat this step.
2. **`scripts/02_pca_qc/`** — PCA (PC1 vs PC2, PC1 vs PC3), a scree plot, and a boxplot of VST-normalized expression per condition.
3. **`scripts/03_volcano_plots/`** — three volcano plot variants per comparison: unlabeled, labeled with the curated gene list, and labeled using a separate, stricter FDR/log2FC cutoff (used only to decide which points get a text label, not which are called significant).
4. **`scripts/04_deg_overlap/`** — a Venn diagram of DEGs across the three comparisons, plus an Excel workbook listing which genes fall into each overlap region (triple overlap, pairwise overlaps, and comparison-exclusive genes).

## Analysis parameters

As described in the paper's Methods:

| Parameter | Value | Meaning |
|---|---|---|
| `expr_thresh` | 5 | minimum raw count for a gene to be considered expressed in a sample |
| `min_frac` | 0.75 | fraction of samples within a group that must clear `expr_thresh` |
| `alpha_thresh` | 0.1 | FDR (Benjamini–Hochberg) significance threshold |
| `lfc_thresh` | 0.5 | minimum absolute log2 fold-change for significance |

Starting from 20,605 annotated genes, restricting to protein-coding genes (Ensembl BioMart) leaves 13,311 genes; the expression filter above further reduces this to 5,520 genes used for all downstream differential expression analyses. Differential expression was tested with DESeq2's Wald test and Benjamini–Hochberg correction for each pairwise comparison among Treg, Tfr, and Tfh.

## Citation

If you use this code, please cite the paper:

Bach L, Chang Y, Arteaga Transito O, Ghasemi M, Steinheuer LM, Steffen T, De Domenico E, Ulas T, Wunderlich FT, Beyer MD, Thurley K, et al. Baumjohann D. TGF-β and IL-2 differentially shape T follicular regulatory cell differentiation and stability in vitro. *Cell Mol Immunol* (2026). https://doi.org/10.1038/s41423-026-01440-9
