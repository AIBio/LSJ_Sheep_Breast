# 🐑 Sheep Mammary Gland scRNA-seq Analysis

This repository contains the full R code and workflow for single-cell RNA sequencing (scRNA-seq) analysis of sheep mammary gland tissues. The analysis was performed using Seurat and multiple complementary downstream tools to uncover cellular composition, interactions, and biological functions at single-cell resolution.

## 📁 Project Structure

The analysis is organized into the following parts:

### 1️⃣ Seurat Pipeline  
Basic preprocessing and quality control, normalization, dimensionality reduction, clustering, and visualization using the Seurat framework.

### 2️⃣ Immune Cluster Analysis  
Subclustering and annotation of immune-related cells, identification of marker genes, and functional enrichment.

### 3️⃣ Epithelial Cluster Analysis  
Refined analysis of epithelial subpopulations, including luminal and basal cell types, with emphasis on lineage hierarchy and function.

### 4️⃣ Final Annotation  
Comprehensive cell type annotation based on canonical markers, literature references, and cross-validation with multiple databases.

### 5️⃣ Trajectory Analysis  
Pseudotime and lineage trajectory inference to explore cellular differentiation dynamics using tools such as Monocle or Slingshot.

### 6️⃣ CellChat Analysis  
Cell-cell communication inference using [CellChat](https://github.com/sqjin/CellChat), focusing on signaling pathway activities across cell types.

### 7️⃣ Integrative Analysis  
Integration of multiple samples or conditions using Seurat’s integration framework, identifying conserved and condition-specific features.

### 8️⃣ Metabolism Analysis  
Analysis of cell-type-specific metabolic signatures using curated gene sets and pathway enrichment tools such as GSVA or clusterProfiler.

### 9️⃣ Fibroblast Analysis  
Focused characterization of fibroblast populations, including subcluster identification and functional profiling.

## 📦 Requirements

This project was built using:
- R (≥ 4.1)
- [Seurat](https://satijalab.org/seurat/)
- [CellChat](https://github.com/sqjin/CellChat)
- [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) or [Slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html)
- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
- Other packages as noted in script headers

You can install all required R packages using:

```r
install.packages("Seurat")
# Other packages can be installed similarly or via Bioconductor
