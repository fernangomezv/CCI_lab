# CCI_lab

# Phase 3 (GEO): Normal → Dysplasia → Cancer Progression

This repository contains the **R script used to analyze oral mucosa progression** (Normal / Dysplasia / Cancer) using public GEO datasets:

- **GSE30784**
- **GSE25099**

The pipeline performs:

- Download and preprocessing of GEO series
- Automatic recoding of lesions into **Normal / Dysplasia / Cancer**
- **Differential expression analysis** (pairwise: CvsN, DvsN, CvsD)
- **GSVA** enrichment for **IFN type I / II / III** signatures  
- **PCA (robust)**, volcano plots, heatmaps
- Export of **Excel summaries** per dataset

> 📝 *Nota (ES): Este script corresponde exactamente al flujo de análisis utilizado en el manuscrito. Cualquier actualización mayor será versionada en este repositorio.*

---

## 1. Repository contents

- `phase3_GEO_IFN_pipeline.R`  
  Main script implementing the Phase 3 pipeline described above.  
  - Generic functions to process any GEO series with Normal / Dysplasia / Cancer labels  
  - Dedicated section with a more detailed workflow for **GSE30784**

You can rename the script in the repo if you prefer a different filename; just keep this README in sync.

---

## 2. Analysis overview

For each GEO series (GSE30784 and GSE25099), the pipeline:

1. **Downloads and loads ExpressionSet objects** using `GEOquery`.
2. **Standardizes lesion labels** (Normal / Dysplasia / Cancer) from GEO metadata using robust regex-based parsing.
3. **Preprocesses expression data**:
   - Detects **RNA-seq (counts)** vs **microarray (intensities)**.
   - For RNA-seq:
     - Filters low-count features.
     - Runs **DESeq2** with `Lesion` as design factor.
     - Uses **VST**-transformed counts for PCA and heatmaps.
   - For microarrays:
     - Applies **log2 transform + quantile normalization** (`limma`).
     - Collapses probes to **gene level** (median) using GPL and Bioconductor annotation (Illumina v2/v3/v4).
4. **Differential expression (AED / DGEA)**:
   - Pairwise comparisons:
     - **Cancer vs Normal (CvsN)**
     - **Dysplasia vs Normal (DvsN)**
     - **Cancer vs Dysplasia (CvsD)**
   - Uses:
     - **DESeq2** for count-like data
     - **limma** for microarray data
5. **IFN signatures and GSVA**:
   - Defines three curated gene sets:
     - `IFN_Type_I`
     - `IFN_Type_II`
     - `IFN_Type_III`
   - Runs **GSVA** (with compatibility for old/new GSVA APIs).
   - Generates **grouped ggbetweenstats plots** by lesion and IFN set.
6. **Visualization**:
   - **Volcano plots** for each contrast (logFC vs –log10 FDR).
   - **Heatmaps** (top 50 genes by minimum FDR across contrasts).
   - **PCA plots** (robust, with variance-based filtering of problematic samples).
   - **GSVA IFN plots** (Normal vs Dysplasia vs Cancer).
7. **Export**:
   - Writes **Excel workbooks** per GEO series with:
     - GSVA scores per sample (IFN signatures)
     - DE tables for each contrast
     - Summary tables with counts of over-/under-expressed genes
   - Writes **CSV** files for:
     - DGE results per contrast
     - GSVA scores (IFN signatures)
     - Probe → gene symbol mappings (for GSE30784)
   - Saves all **figures as PNG** (and selected plots also as PDF).

---

## 3. Requirements

### 3.1. R version

- Recommended: **R ≥ 4.2**

### 3.2. Required packages

The script takes care of installing missing packages automatically via `BiocManager` and `install.packages`. It relies on:

**Bioconductor packages**

- `BiocManager`
- `GEOquery`
- `DESeq2`
- `limma`
- `GSVA`
- `apeglm`
- `Biobase`
- `AnnotationDbi`
- `illuminaHumanv3.db`
- `illuminaHumanv2.db`
- `illuminaHumanv4.db`

**CRAN packages**

- `tidyverse`
- `openxlsx`
- `pheatmap`
- `RColorBrewer`
- `ggplot2`
- `ggrepel`
- `ggstatsplot`
- `stringr`
- `tibble`
- `rlang`
- `dplyr`

If you want to pre-install everything manually:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_needed <- c(
  "GEOquery","DESeq2","limma","GSVA","apeglm","Biobase",
  "AnnotationDbi","illuminaHumanv3.db","illuminaHumanv2.db","illuminaHumanv4.db"
)
for (pkg in bioc_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cran_needed <- c(
  "tidyverse","openxlsx","pheatmap","RColorBrewer",
  "ggplot2","ggrepel","ggstatsplot","stringr","tibble","rlang","dplyr"
)
for (pkg in cran_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 4. how to run
# 1) Set your working directory to the repo
setwd("/path/to/this/repo")

# 2) Source the script
source("phase3_GEO_IFN_pipeline.R")

# 3) Run Phase 3 for both GEO series
process_gse_phase3("GSE30784")
process_gse_phase3("GSE25099")

# (Optional) The dedicated section for GSE30784 is executed automatically
# if you run the full script top-to-bottom.

# 5. Reproducibilitty
set.seed(123)

# 6. Citation

# If you use this code or derivative pipelines in your work, please cite: The GEO datasets:
# GSE30784 – Oral mucosa progression (Normal / Dysplasia / Cancer)
# GSE25099 – Complementary oral lesion dataset (Normal / Dysplasia / Cancer)
# The analytical tools:
# Love MI et al., Genome Biol, 2014 (DESeq2)
# Ritchie ME et al., Nucleic Acids Res, 2015 (limma)
# Hänzelmann S et al., BMC Bioinformatics, 2013 (GSVA)
# and any relevant Bioconductor annotation packages.

# 7. Contact
# For questions, suggestions, or bug reports, please open an Issue in this repository or contact:
# Fernán Gómez Valenzuela
# Cancer & Chronic Inflammation (CCI) Lab
# Center for Precision Oncology
# fernan.gomez@umayor.cl
# https://orcid.org/0000-0001-6889-976X
