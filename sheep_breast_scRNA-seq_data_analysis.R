# ========================
# - Title: Sheep breast scRNAseq data
# - Dataset: ***
# - Article: ***
# - Content:
# -         1st part: Seurat pipeline;
# -         2nd part: Immune cluster analysis;
# -         3rd part: Epithelial cluster analysis;
# -         4th part: Final annotation;
# -         5th part: Trajectory analysis;
# -         6th part: Cell Chat analysis;
# -         7th part: Integrative analysis;
# -         8th part: Metabolism analysis;
# -         9th part: Fibroblast analysis;
# ========================



# ========================
# 1st part: Global setting ----
# ========================

### >>> 1. Setting workding directory
setwd("/home/yhw/bioinfo/proj12th")
dir.create("R/CodeData", recursive = T)
dir.create("R/Graphs", recursive = T)
dir.create("R/Table", recursive = T)



# =================
# 2nd part: Library ----
# =================

### >>> 1. Load packages
pks <- c("forcats", "dplyr", "tidyr", "stringr", "circlize",
         "ggplot2", "ggalluvial", "ggrepel", "RColorBrewer", "cowplot", "Seurat",
         "harmony", "clustree", "moonBook", "webr")
for (i in pks) {
  library(i, character.only = T)
}


### >>> 2. Functions
theme_dp <- theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 14, angle = 0),
                  axis.title.y = element_text(face="plain", colour = "#000000", size = 14, angle = 90),
                  axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
                  axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0))
cols.gp <- c(`LL-rep1` = "#A6CEE3", `LL-rep2` = "#1F78B4", `PL-rep1` = "#FB9A99", `PL-rep2` = "#E31A1C")
pd.col <- c("#607d8b","#795548","#ff5722","#ffc107","#cddc39","#4caf50","#009688",
            "#00bcd4","#2196f3","#3f51b5","#673ab7","#9c27b0","#e91e63","#f44336",
            "#b0bec5","#bcaaa4","#ffab91","#ffe082","#e6ee9c","#a5d6a7","#80cbc4",
            "#80deea","#90caf9","#9fa8da","#b39ddb","#ce93d8","#f48fb1","#ef9a9a",
            "#37474f","#4e342e","#d84315","#ff8f00","#9e9d24","#2e7d32","#00695c",
            "#00838f","#1565c0","#283593","#4527a0","#6a1b9a","#ad1457","#c62828")



# =========================
# 3rd part: Seurat pipeline ----
# =========================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Seurat_output")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Load data
raw.file <- "/home/yhw/bioinfo/project-temp/jobs/proj12th/analysis"
bc.dir <- grep("rep", list.files(raw.file, "filtered_feature_bc_matrix", full.names = T, recursive = T), value = T)
sr.list <- list()
for (dir in bc.dir) {
  sr.list[[dir]] <- LoadSCdata(mode = "10x", dir = gsub(".h5", "", dir),
                               proj.name = str_split_fixed(dir, "/", 10)[,9])
  sr.list[[dir]]$Group <- str_split_fixed(dir, "/", 10)[,9]
}
names(sr.list) <- str_split_fixed(bc.dir, "/", 10)[,9]
sr.sheep <- Reduce(function(x, y) merge(x, y, add.cell.ids = c(x@project.name, y@project.name)), sr.list)


### >>> 3. Seurat pipeline
# quality control
mt.gene <- read.table("/home/yhw/bioinfo/project-temp/jobs/proj12th/analysis/metadata/mt_gene.txt")
sr.sheep[["percent.mt"]] <- round(PercentageFeatureSet(sr.sheep, features = mt.gene$V1), 1)
pdf(file.path(res.out, "Violin_plot_to_show_quality_control_before_filtering.pdf"), height = 4.5, width = 10)
VlnPlot(sr.sheep, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "Group", cols = cols.gp)
dev.off()
# filter cells
sr.sheep <- subset(sr.sheep, nFeature_RNA >= 500 & nFeature_RNA <= 5000 & nCount_RNA >= 1000 & percent.mt <= 30)
pdf(file.path(res.out, "Violin_plot_to_show_quality_control_after_filtering.pdf"), height = 4.5, width = 10)
VlnPlot(sr.sheep, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "Group", cols = cols.gp)
dev.off()
# pca
sr.sheep <- NormalizeData(sr.sheep) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.sheep), verbose = FALSE)
sr.sheep <- RunPCA(sr.sheep, features = VariableFeatures(sr.sheep), verbose = FALSE)
pdf(file.path(res.out, "Elbowplot_after_filtering.pdf"), height = 4, width = 4)
ElbowPlot(sr.sheep)
dev.off()
# cluster cells + UMAP/tSNE
dim.n <- 12
sr.sheep <- RunUMAP(sr.sheep, dims = 1:dim.n) %>% RunTSNE(dims = 1:dim.n)
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_before_correction.pdf"), height = 5, width = 6)
DimPlot(sr.sheep, reduction = "umap", group.by = c("Group"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.05, cols = cols.gp)
dev.off()
# remove batch effect
sr.sheep <- RunHarmony(sr.sheep, group.by.vars = "Group")
sr.sheep <- RunUMAP(sr.sheep, reduction = "harmony", dims = 1:dim.n)
sr.sheep <- FindNeighbors(sr.sheep, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = seq(0, 2, 0.2))
pdf(file.path(res.out, "Tree_plot_to_show_clustering_after_correction.pdf"), height = 10, width = 12.5)
clustree(sr.sheep, prefix = "RNA_snn_res.",
         node_label = "percent.mt", node_label_aggr = "median") +
  scale_color_brewer(palette = "Blues") +
  scale_fill_brewer(palette = "Blues")
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_after_correction.pdf"), height = 5, width = 13)
p1 <- DimPlot(sr.sheep, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.05, cols = cols.gp)
p2 <- DimPlot(sr.sheep, reduction = "umap", group.by = c("RNA_snn_res.2"),
              label = TRUE, repel = TRUE, pt.size = 0.05, cols = pd.col)
p1 + p2
dev.off()
sr.sheep$Final_clusters <- paste0("C", sr.sheep$RNA_snn_res.2)
# find markers
Idents(sr.sheep) <- sr.sheep$RNA_snn_res.2
markers.all <- FindAllMarkers(sr.sheep, assay = "RNA", only.pos = T, logfc.threshold = log2(1.25))
write.csv(markers.all, file.path(res.out, "Marker_genes_of_42_populations.csv"), row.names = F)
markers.go <- list()
for (type in unique(as.character(markers$cluster))) {
  pd.gene <- markers %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == type)
  markers.go[[type]] <- Pipe.GO(species = "human", genelist = pd.gene$gene, basename = type,
                                genetype = "SYMBOL", res.out = file.path(sr.out, "Marker_GO"))
}
#
pdf(file.path(outdir, "Hs_Ss_heart_cell_type_marker_gene_expression_after_imputation.pdf"), height = 4, width = 12)
marker.gene <- c("EPCAM",
                 # myoepithelial lineage
                 "KRT17", "KRT5", "ACTA2", "MYL9", "MYLK", "MYH11",
                 # - Luminal cells (Krt19, Krt18, and Krt8)
                 "KRT19", "KRT18", "KRT8",
                 # HS cells
                 "PRLR", "ESR1", "CITED1", "PROM1",
                 # AV cells
                 "ENSOARG00020009476", "CSN3", "WFDC8", "LTF", "ELF5",
                 # luminal progenitor
                 "KIT", "ALDH1A3", "CD14",
                 # mature alveolar
                 "WAPL",
                 # - fibroblasts
                 "COL1A1", "COL1A2", "COL3A1", "FN1",
                 # - vascular/lymphatic cells
                 "PECAM1", "CDH5", "ENG",
                 # vascular endothelial cells
                 "SOX17", "SELE",
                 # pericytes
                 "RGS5", "DES", "NOTCH3",
                 # lymphatic endothelial cells
                 "MMRN1", "PROX1", "FLT4",
                 # - immune cells
                 "PTPRC",
                 # myeloid cells
                 "ENSOARG00020013258",
                 # dendritic cells
                 "NAPSA", "TRAF1", "FLT3",
                 # Ma macrophages
                 "CSF1R", "FCGR3A", "ENSOARG00020017289", "CD163",
                 # Mb macrophages
                 "MMP12", "MMP13", "SPIC",
                 # lymphocytes
                 "CD3D", "CD3E", "CD3G",
                 # T cells
                 "CD8A", "CD8B1", "CD4",
                 # B cells
                 "CD79A", "CD79B",
                 # - adipocyte
                 "ADIPOQ", "ADRB1", "ACSL9", "CYP2E1", "ALDH1A1")
marker.gene <- intersect(marker.gene, rownames(sr.sheep))
DotPlot(sr.sheep, features = marker.gene, cols = c("#eaeaea", "#fc0330"),
        col.min = 0, col.max = 2, group.by = "RNA_snn_res.2",
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker.gene) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
markers.all %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster == 30) -> pd.gene
FeaturePlot(sr.sheep, features = pd.gene$gene[1:20], ncol = 5)
# Epithelial cells --> cluster 4,5,22,27,29,31,33,40
FeaturePlot(sr.sheep, features = c("EPCAM"), order = T)
FeaturePlot(sr.sheep, features = c("PEBP4", "FOLR1", "AQP4", "CLDN18", # mature alveolar
                                   "KRT19", "CD24", "ITGA6", "KIT", "SLPI" # luminal progenitor
))
FeaturePlot(sr.sheep, features = c("MARCO", "ALDH2", "APOE", "TREM2", # mature alveolar
                                   "FABP4" # luminal progenitor
))
# myoepithelial lineage --> cluster 11
FeaturePlot(sr.sheep, features = c("KRT17", "KRT5", "ACTA2", "MYL9", "MYLK", "MYH11",
                                   "TP63", "KRT17", "MYH11", "CALD1"), order = T)
# Vascular endothelial cells --> cluster 15
FeaturePlot(sr.sheep, features = c("SOX17", "SELE"), order = T)
# Lymphatic endothelial cells --> cluster 20
FeaturePlot(sr.sheep, features = c("MMRN1", "PROX1", "FLT4"), order = T)
# Pericytes --> cluster 1,2,12,24
FeaturePlot(sr.sheep, features = c("RGS5", "DES", "NOTCH3"), order = T)
# Smooth muscle cells --> cluster 7,8,9,19
FeaturePlot(sr.sheep, features = c("ACTA2", "MYL9", "MYLK", "MYH11", "KRT5", "TAGLN"))
# Fibroblasts --> cluster 10,16,21,26,36,41
FeaturePlot(sr.sheep, features = c("COL1A1", "COL1A2", "COL3A1", "FN1"))
# Immune cell --> cluster 0,3,6,13,14,17,18,23,25,28,32,34,35,37,38
FeaturePlot(sr.sheep, features = c("PTPRC"))
# Adipose-related genes --> cluster *
FeaturePlot(sr.sheep, features = c("FABP4", "PPARG", "TMEM26", "PLIN1", "ACSL1", # Fat cell
                                   "TCF21", "CEBPA", "HOX3", "LEP", # White fat cell
                                   "PAT2", "P2RX5", "UCP1", # Brown fat cell
                                   "KAZALD1", "LTBP4", "LHFP", "C1QTNF7", "CYGB"))
FeaturePlot(sr.sheep, features = c("ADIPOR1", "ADIPOR2", "PLIN2", "EBF1", "ACSL1", # Fat cell
                                   "FOXO1", "LPIN1", "GHR", "PLIN1", "PPARG", "ADIPOQ"))
FeaturePlot(sr.sheep, features = c("PDGFRA", "SCA1", "CD29", "CD24", "CD34", # mesenchymal cell
                                   "DPP4", "WNT2", "BMP7", "PIL6", "CDH11", "PDGFRB"), order = T)
FeaturePlot(sr.sheep, features = c("ICAM1", "PPARG", "FABP4", "CD36", "DLK1",
                                   "CLEC11A", "FMO2", "ADIPOQ", "PLIN1", "CAR3"))
FeaturePlot(sr.sheep, features = c("APOD", "GPX3"))
# new fibroblast cells --> cluster *
FeaturePlot(sr.sheep, features = c("CD140B", "EGFR", "LEPR", "TUBB2A", "TUBB2B",
                                   "TUBB6", "FBLN2", "FBLN5", "SPRY1"))

# rename cells
new.cluster.ids <- c("Immune cells","Pericyte","Pericyte","Immune cells","Epithelial cells",
                     "Epithelial cells","Immune cells","Smooth muscle cells","Smooth muscle cells",
                     "Smooth muscle cells","Fibroblast","Myoepithelial cells","Pericyte",
                     "Immune cells","Immune cells","Vascular endothelial cells","Fibroblast",
                     "Immune cells","Immune cells","Smooth muscle cells","Lymphatic endothelial cells",
                     "Fibroblast","Epithelial cells","Immune cells","Pericyte","Immune cells",
                     "Fibroblast","Epithelial cells","Immune cells","Epithelial cells","Doublets",
                     "Epithelial cells","Immune cells","Epithelial cells","Immune cells","Immune cells",
                     "Fibroblast","Immune cells","Immune cells","Doublets","Epithelial cells","Doublets")
names(new.cluster.ids) <- levels(sr.sheep)
sr.sheep <- RenameIdents(sr.sheep, new.cluster.ids)
sr.sheep$CellType <- as.character(Idents(sr.sheep))
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_after_correction_with_annotation.pdf"), height = 5, width = 14)
p1 <- DimPlot(sr.sheep, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.25, cols = cols.gp)
p2 <- DimPlot(sr.sheep, reduction = "umap", group.by = c("CellType"),
              label = TRUE, repel = TRUE, pt.size = 0.25, cols = pd.col)
p1 + p2
dev.off()



# =================================
# 4th part: Immune cluster analysis ----
# =================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Immune_cluster")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Load data
# filter cells
sr.sheep.immune <- subset(sr.sheep, RNA_snn_res.2 %in% c(0,3,6,13,14,17,18,23,25,28,32,34,35,37,38))
# pca
sr.sheep.immune <- NormalizeData(sr.sheep.immune) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.sheep.immune), verbose = FALSE)
sr.sheep.immune <- RunPCA(sr.sheep.immune, features = VariableFeatures(sr.sheep.immune), verbose = FALSE)
pdf(file.path(res.out, "Elbowplot_after_filtering.pdf"), height = 4, width = 4)
ElbowPlot(sr.sheep.immune)
dev.off()
# cluster cells + UMAP/tSNE
dim.n <- 6
sr.sheep.immune <- RunUMAP(sr.sheep.immune, dims = 1:dim.n) %>% RunTSNE(dims = 1:dim.n)
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_before_correction.pdf"), height = 5, width = 6)
DimPlot(sr.sheep.immune, reduction = "umap", group.by = c("Group"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.05, cols = cols.gp)
dev.off()
# remove batch effect
sr.sheep.immune <- RunHarmony(sr.sheep.immune, group.by.vars = "Group")
sr.sheep.immune <- RunUMAP(sr.sheep.immune, reduction = "harmony", dims = 1:dim.n)
sr.sheep.immune <- FindNeighbors(sr.sheep.immune, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = seq(0, 2, 0.2))
pdf(file.path(res.out, "Tree_plot_to_show_clustering_after_correction.pdf"), height = 10, width = 12.5)
clustree(sr.sheep.immune, prefix = "RNA_snn_res.",
         node_label = "percent.mt", node_label_aggr = "median") +
  scale_color_brewer(palette = "Blues") +
  scale_fill_brewer(palette = "Blues")
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_after_correction.pdf"), height = 5, width = 12)
p1 <- DimPlot(sr.sheep.immune, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = cols.gp)
p2 <- DimPlot(sr.sheep.immune, reduction = "umap", group.by = c("RNA_snn_res.2"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = pd.col)
p1 + p2
dev.off()
# find markers
Idents(sr.sheep.immune) <- sr.sheep.immune$RNA_snn_res.2
markers <- FindAllMarkers(sr.sheep.immune, assay = "RNA", only.pos = T, logfc.threshold = log2(1.25))
write.csv(markers, file.path(res.out, "Marker_genes_of_30_immune_subpopulations.csv"), row.names = F)
markers.go <- list()
for (type in unique(as.character(markers$cluster))) {
  pd.gene <- markers %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == type)
  markers.go[[type]] <- Pipe.GO(species = "human", genelist = pd.gene$gene, basename = type,
                                genetype = "SYMBOL", res.out = file.path(sr.out, "Marker_GO"))
}
#
sr.sheep.immune <- SeuratWrappers::RunALRA(sr.sheep.immune)
sr.sheep.immune <- ScaleData(sr.sheep.immune, features = rownames(sr.sheep.immune))

pdf(file.path(outdir, "Hs_Ss_heart_cell_type_marker_gene_expression_after_imputation.pdf"), height = 4, width = 12)
marker.gene <- c(# - immune cells
                 "PTPRC",
                 # myeloid cells
                 "ENSOARG00020013258",
                 # dendritic cells
                 "NAPSA", "TRAF1", "FLT3",
                 # Ma macrophages
                 "CSF1R", "FCGR3A", "ENSOARG00020017608", "ENSOARG00020017289", "CD163",
                 # Mb macrophages
                 "MMP12", "MMP13", "SPIC",
                 # lymphocytes
                 "CD3D", "CD3E", "CD3G",
                 # T cells
                 "CD8A", "CD8B1", "CD4", "CXCR6",
                 # B cells
                 "CD79A", "CD79B")
marker.gene <- intersect(marker.gene, rownames(sr.sheep.immune))
DotPlot(sr.sheep.immune, features = marker.gene, cols = c("#eaeaea", "#fc0330"),
        col.min = 0, col.max = 2, group.by = "RNA_snn_res.2",
        scale = T, scale.min = 0, scale.max = 100, assay = "alra") + RotatedAxis() +
  scale_x_discrete(labels = marker.gene) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster == 26) -> pd.gene
FeaturePlot(sr.sheep.immune, features = pd.gene$gene[1:20], ncol = 5)
# Neutrophil --> cluster 22
FeaturePlot(sr.sheep.immune, features = c("HPSE", "CSF3R", "CD93", "ALPL", "S100A12", "MGAM"), ncol = 4)
# B cell --> cluster 27,24,(a part of 24)
FeaturePlot(sr.sheep.immune, features = c("CD79A", "CD79B", "CD19", "MS4A1"), ncol = 4)
# Mast cell --> cluster 28
FeaturePlot(sr.sheep.immune, features = c("CPA3", "MS4A2", "KIT", "TPSB2"), ncol = 4)
# Ma Macrophage --> cluster 17,18
FeaturePlot(sr.sheep.immune, features = c("CSF1R", "FCGR3A", "ENSOARG00020017608", "ENSOARG00020017289", "CD163",
                                          "C1QA", "C1QB", "C1QC"), ncol = 4)
# Mb Macrophage --> cluster 15
FeaturePlot(sr.sheep.immune, features = c("CSF1R", "FCGR3A", "ENSOARG00020017608", "ENSOARG00020017289", "CD163",
                                          "C1QA", "C1QB", "C1QC"), ncol = 4)
# CD4+ T cell --> cluster 2,3,8,10,11,14
FeaturePlot(sr.sheep.immune, features = c("CD3D", "CD3E", "CD3G", "CD4"), ncol = 4)
# CD8+ T cell --> cluster 0,1,7,16,21,23
FeaturePlot(sr.sheep.immune, features = c("CD3D", "CD3E", "CD3G", "CD8A"), ncol = 4)
# NK cell --> cluster 4,5,6,9,12,13,19,20,25,26
FeaturePlot(sr.sheep.immune, features = c("ATP9A", "BATF", "ANP32B", "CENPF", "CXCR6"), ncol = 4)
# Dendritic cell --> cluster 29,30,(a part of 24)
FeaturePlot(sr.sheep.immune, features = c("NAPSA", "TRAF1", "CD209A", "FLT3"), ncol = 4)
FeaturePlot(sr.sheep.immune, features = c("EBI3", "FABP4", "CCL22", "CCL17"), ncol = 4)
# Myeloid cells --> cluster *
FeaturePlot(sr.sheep.immune, features = c("ENSOARG00020013258", "FCGR3A", "ENSOARG00020016111"), ncol = 3)
# plasma cells --> cluster *
FeaturePlot(sr.sheep.immune, features = c("ENSOARG00020001068", "ENSOARG00020025540",
                                          "ENSOARG00020025540", "ENSOARG00020025540"), ncol = 2)
# rename cells
new.cluster.ids <- c("CD8+ T cells","CD8+ T cells","CD4+ T cells","CD4+ T cells","NK cells","NK cells","NK cells",
                     "CD8+ T cells","CD4+ T cells","NK cells","CD4+ T cells","CD4+ T cells","NK cells","NK cells",
                     "CD4+ T cells","Mb Macrophage","CD8+ T cells","Ma Macrophage","Ma Macrophage","NK cells",
                     "NK cells","CD8+ T cells","Neutrophil","CD8+ T cells","B cells","NK cells","NK cells","B cells",
                     "Mast cells","Dendritic cells","Dendritic cells")
names(new.cluster.ids) <- levels(sr.sheep.immune)
sr.sheep.immune <- RenameIdents(sr.sheep.immune, new.cluster.ids)
sr.sheep.immune$CellType <- as.character(Idents(sr.sheep.immune))
# refined cell types
cell.pos <- data.frame(group = "highlighted",
                       xmin = -10, xmax = -5, ymin = 6.5, ymax = 15)
tmp <- ExtractCellByPos(object = sr.sheep.immune, object.type = "seurat",
                        dim.name = "umap", group.coor = cell.pos,
                        group.name = "tmp", pt.size = 1)
keep.cell <- intersect(rownames(subset(sr.sheep.immune@meta.data, CellType == "B cells")), tmp$id$highlighted)
sr.sheep.immune$CellType[colnames(sr.sheep.immune) %in% keep.cell] <- "Dendritic cells"
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_after_correction_with_annotation.pdf"), height = 5, width = 13)
p1 <- DimPlot(sr.sheep.immune, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = cols.gp)
p2 <- DimPlot(sr.sheep.immune, reduction = "umap", group.by = c("CellType"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = pd.col)
p1 + p2
dev.off()
#
pd.gene <- c(
  # Ma Macrophage
  "C1QA", "C1QB", "C1QC", "CSF1R", "CD163",
  # Mb Macrophage
  "MMP12", "SPIC",
  # Dendritic cells
  "TRAF1", "FLT3", "FABP4", "CCL22", "CCL17",
  # B cells
  "CD79A", "CD79B", "CD19", "MS4A1",
  # Neutrophil
  "CSF3R", "ALPL", "S100A12", "MGAM",
  # Mast cells
  "CPA3", "MS4A2", "KIT", "TPSB2",
  # lymphocytes
  "CD3D", "CD3E", "CD3G",
  # T cells
  "CD4", "CD8A",
  # NK cells
  "GZMA", "AFAP1L2", "SH2D1B"
)
pd.cell <- c("Ma Macrophage","Mb Macrophage","Dendritic cells","B cells","Neutrophil",
             "Mast cells","CD8+ T cells","CD4+ T cells","NK cells")
pdf(file.path(res.out, "Dot_plot_to_show_marker_gene_expression_with_annotation.pdf"), height = 5, width = 12)
DotPlot(sr.sheep.immune, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        col.min = -0.5, col.max = 1, group.by = "CellType",
        scale = T, scale.min = 0, scale.max = 75) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  scale_y_discrete(limit = pd.cell) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_marker_gene_expression_with_annotation.pdf"), height = 60, width = 30)
FeaturePlot(sr.sheep.immune, features = c("PTPRC", pd.gene), cols = c("#E7E7E7", "#1A5A9B"), pt.size = 0.25)
dev.off()
DefaultAssay(sr.sheep.immune) <- "alra"
pdf(file.path(res.out, "Scatter_plot_to_show_marker_gene_expression_with_annotation_update.pdf"), height = 10, width = 22)
FeaturePlot(sr.sheep.immune, features = c("ENSOARG00020016111", "CSF1R", "FCGR3A", "ENSOARG00020017608", 
                                          "ENSOARG00020017289", "MRC1", 'CD163', "MMP12",  
                                          "SPIC", "CD209"), cols = c("#E7E7E7", "#1A5A9B"), 
            pt.size = 0.25, order = T, ncol = 4)
dev.off()
# find markers (cluster)
Idents(sr.sheep.immune) <- sr.sheep.immune$CellType
markers.immune <- FindAllMarkers(sr.sheep.immune, assay = "RNA", only.pos = T, logfc.threshold = log2(1.25))
write.csv(markers.immune, file.path(res.out, "Marker_genes_of_final_clusters.csv"), row.names = F)
pd.gene <- markers.immune[-grep("ENSOARG", markers.immune$gene), ] %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC)
# plot heatmap (cluster)
tmp <- subset(sr.sheep.immune, downsample = 100)
tmp <- ScaleData(tmp, features = rownames(tmp))
library(SeuratWrappers)
tmp <- SeuratWrappers::RunALRA(tmp)
tmp <- ScaleData(tmp, features = rownames(tmp))
DefaultAssay(tmp) <- "RNA"
pdf(file.path(res.out, "Heatmap_to_show_top5_marker_gene_expression_cluster_update_RNA.pdf"), height = 10, width = 10)
DoHeatmap(tmp, features = unique(pd.gene$gene))
dev.off()
DefaultAssay(tmp) <- "alra"
pdf(file.path(res.out, "Heatmap_to_show_top5_marker_gene_expression_cluster_update_ALRA.pdf"), height = 10, width = 10)
DoHeatmap(tmp, features = unique(pd.gene$gene))
dev.off()


### >>> 3. T cell cluster analysis
# filter cells
sr.sheep.immune.tcell <- subset(sr.sheep.immune, RNA_snn_res.2 %in% c(21,1,11,23,0,2,8,3,10,7,16,14,19,9,20,4,12,13,5,6,25,26))
# pca
sr.sheep.immune.tcell <- NormalizeData(sr.sheep.immune.tcell) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.sheep.immune.tcell), verbose = FALSE)
sr.sheep.immune.tcell <- RunPCA(sr.sheep.immune.tcell, features = VariableFeatures(sr.sheep.immune.tcell), verbose = FALSE)
pdf(file.path(res.out, "Elbowplot_after_filtering_Tcell.pdf"), height = 4, width = 4)
ElbowPlot(sr.sheep.immune.tcell)
dev.off()
# cluster cells + UMAP/tSNE
dim.n <- 10
sr.sheep.immune.tcell <- RunUMAP(sr.sheep.immune.tcell, dims = 1:dim.n) %>% RunTSNE(dims = 1:dim.n)
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_before_correction.pdf"), height = 5, width = 6)
DimPlot(sr.sheep.immune.tcell, reduction = "umap", group.by = c("Group"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.05, cols = cols.gp)
dev.off()
# remove batch effect
sr.sheep.immune.tcell <- RunHarmony(sr.sheep.immune.tcell, group.by.vars = "Group")
sr.sheep.immune.tcell <- RunUMAP(sr.sheep.immune.tcell, reduction = "harmony", dims = 1:dim.n)
sr.sheep.immune.tcell <- FindNeighbors(sr.sheep.immune.tcell, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = seq(0, 2, 0.2))
pdf(file.path(res.out, "Tree_plot_to_show_clustering_after_correction.pdf"), height = 10, width = 12.5)
clustree(sr.sheep.immune.tcell, prefix = "RNA_snn_res.",
         node_label = "percent.mt", node_label_aggr = "median") +
  scale_color_brewer(palette = "Blues") +
  scale_fill_brewer(palette = "Blues")
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_after_correction_Tcell.pdf"), height = 5, width = 12)
p1 <- DimPlot(sr.sheep.immune.tcell, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = cols.gp)
p2 <- DimPlot(sr.sheep.immune.tcell, reduction = "umap", group.by = c("RNA_snn_res.2"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = pd.col)
p1 + p2
dev.off()



# =====================================
# 5th part: Epithelial cluster analysis ----
# =====================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Epithelial_cluster")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Load data
# filter cells
sr.sheep.epithelial <- subset(sr.sheep, RNA_snn_res.2 %in% c(4,5,11,22,27,29,31,33,40))
# pca
sr.sheep.epithelial <- NormalizeData(sr.sheep.epithelial) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.sheep.epithelial), verbose = FALSE)
sr.sheep.epithelial <- RunPCA(sr.sheep.epithelial, features = VariableFeatures(sr.sheep.epithelial), verbose = FALSE)
pdf(file.path(res.out, "Elbowplot_after_filtering.pdf"), height = 4, width = 4)
ElbowPlot(sr.sheep.epithelial)
dev.off()
# cluster cells + UMAP/tSNE
dim.n <- 10
sr.sheep.epithelial <- RunUMAP(sr.sheep.epithelial, dims = 1:dim.n) %>% RunTSNE(dims = 1:dim.n)
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_before_correction.pdf"), height = 5, width = 6)
DimPlot(sr.sheep.epithelial, reduction = "umap", group.by = c("Group"),
        ncol = 1, label = TRUE, repel = TRUE, pt.size = 0.05, cols = cols.gp)
dev.off()
# remove batch effect
sr.sheep.epithelial <- RunHarmony(sr.sheep.epithelial, group.by.vars = "Group")
sr.sheep.epithelial <- RunUMAP(sr.sheep.epithelial, reduction = "harmony", dims = 1:dim.n)
sr.sheep.epithelial <- FindNeighbors(sr.sheep.epithelial, reduction = "harmony", dims = 1:dim.n) %>%
  FindClusters(resolution = seq(0, 2, 0.2))
pdf(file.path(res.out, "Tree_plot_to_show_clustering_after_correction.pdf"), height = 10, width = 12.5)
clustree(sr.sheep.epithelial, prefix = "RNA_snn_res.",
         node_label = "percent.mt", node_label_aggr = "median") +
  scale_color_brewer(palette = "Blues") +
  scale_fill_brewer(palette = "Blues")
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_after_correction.pdf"), height = 5, width = 12)
p1 <- DimPlot(sr.sheep.epithelial, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = cols.gp)
p2 <- DimPlot(sr.sheep.epithelial, reduction = "umap", group.by = c("RNA_snn_res.2"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = pd.col)
p1 + p2
dev.off()
# find markers
Idents(sr.sheep.epithelial) <- sr.sheep.epithelial$RNA_snn_res.2
markers.epi <- FindAllMarkers(sr.sheep.epithelial, assay = "RNA", only.pos = T, logfc.threshold = log2(1.25))
write.csv(markers.epi, file.path(res.out, "Marker_genes_of_28_epithelial_subpopulations.csv"), row.names = F)
marker.gene <- c("EPCAM",
                 # myoepithelial lineage
                 "KRT17", "KRT5", "ACTA2", "MYL9", "MYLK", "MYH11",
                 # Luminal cells
                 "KRT19", "KRT18", "KRT8",
                 # HS cells
                 "PRLR", "ESR1", "CITED1", "PROM1",
                 # AV cells
                 "ENSOARG00020009476", "CSN3", "WFDC8", "LTF", "ELF5",
                 # luminal progenitor
                 "KIT", "ALDH1A3", "CD14",
                 # mature luminal cells
                 "FOXA1", "ESR1",
                 # mature alveolar
                 "WAPL")
marker.gene <- intersect(marker.gene, rownames(sr.sheep.epithelial))
DotPlot(sr.sheep.epithelial, features = marker.gene, cols = c("#eaeaea", "#fc0330"),
        col.min = 0, col.max = 2, group.by = "RNA_snn_res.2",
        scale = T, scale.min = 0, scale.max = 100) + RotatedAxis() +
  scale_x_discrete(labels = marker.gene) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
markers.epi %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster == 1) -> pd.gene
FeaturePlot(sr.sheep.epithelial, features = pd.gene$gene[1:20], ncol = 5)
# Epithelial cells
FeaturePlot(sr.sheep.epithelial, features = c("EPCAM"), order = T)
# Luminal progenitor --> cluster 2,5,13,17
FeaturePlot(sr.sheep.epithelial, features = c("KIT", "ALDH1A3", "CD14"))
# Luminal cells --> cluster 4,6,7,8,11,12,15,16,19,24
FeaturePlot(sr.sheep.epithelial, features = c("KRT19", "KRT18", "KRT8", "EPCAM"))
# HS cells --> cluster 3,10,23
FeaturePlot(sr.sheep.epithelial, features = c("PRLR", "ESR1", "CITED1", "PROM1"))
# AV cells --> cluster 1,22
FeaturePlot(sr.sheep.epithelial, features = c("ENSOARG00020009476", "CSN3", "WFDC8", "LTF", "ELF5"))
# Myoepithelial cells --> cluster 0,9,18,20,26
FeaturePlot(sr.sheep.epithelial, features = c("KRT17", "KRT5", "ACTA2", "MYL9", "MYLK", "MYH11"), ncol = 3)
# HS-AV --> cluster 14,21,25,27
FeaturePlot(sr.sheep.epithelial, features = c("PTPRC", "CD3D", "CD3E", "CD3G"), ncol = 2)

# hormone-responsive (HR) mature LCs --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("ESR1", "PRLR"))
# secretory mature LCs --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("ELF5", "ALDH1A3", "KIT", "SLPI", "KRT23"), ncol = 3)
# new AV cells --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("CD14", "FOLR1", "KRT5", "KRT14", "KRT17", "CAV1", "ANXA8",
                                              "CD95", "CD47", "RANK", "CD47", "TAPBP", "CTSS"))
# new HS cells --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("ESR1", "PDZK1", "SERPINA1", "SPDEF", "CITED1",
                                              "CNK1", "AKAP13", "WNT4"))
# BS cells --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("K14", "K17", "CD10", "IEGs", "FOS", "JUN", "EGR1", "IER2"))
# Mature alveolar --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("FOXA1", "ESR1", "WAPL", "ENSOARG00020004485"))
# milk biosynthesis-related genes --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("ENSOARG00020009476", "CSN3", "LTF"))
# luminal lineage 1 --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("CLDN4", "JUN", "KLF6"))
# luminal lineage 2 --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("XDH", "CSN3"))
# mature luminal cells --> cluster *
FeaturePlot(sr.sheep.epithelial, features = c("ESR1", "PGR", "FOXA1"))

# rename cells
new.cluster.ids <- c("Myoepithelial cells","AV cells","Luminal progenitor","HS cells","Luminal cells",
                     "Luminal progenitor","Luminal cells","Luminal cells","Luminal cells","Myoepithelial cells",
                     "HS cells","Luminal cells","Luminal cells","Luminal progenitor","HS-AV","Luminal cells",
                     "Luminal cells","Luminal progenitor","Myoepithelial cells","Luminal cells","Myoepithelial cells",
                     "HS-AV","AV cells","HS cells","Luminal cells","HS-AV","Myoepithelial cells","HS-AV")
names(new.cluster.ids) <- levels(sr.sheep.epithelial)
sr.sheep.epithelial <- RenameIdents(sr.sheep.epithelial, new.cluster.ids)
sr.sheep.epithelial$CellType <- as.character(Idents(sr.sheep.epithelial))
sr.sheep.epithelial$CellType <- gsub("HS-AV", "HS-AV cells", sr.sheep.epithelial$CellType)
sr.sheep.epithelial$Stage <- gsub("-rep.", "", sr.sheep.epithelial$Group)
# refined cell types
cell.pos <- data.frame(group = "highlighted",
                       xmin = 7.5, xmax = 15, ymin = -2.5, ymax = 5)
tmp <- ExtractCellByPos(object = sr.sheep.epithelial, object.type = "seurat",
                        dim.name = "umap", group.coor = cell.pos,
                        group.name = "tmp", pt.size = 1)
keep.cell <- intersect(rownames(subset(sr.sheep.epithelial@meta.data, CellType == "Luminal cells")), tmp$id$highlighted)
sr.sheep.epithelial$CellType[colnames(sr.sheep.epithelial) %in% keep.cell] <- "Myoepithelial cells"
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_after_correction_with_annotation.pdf"), height = 5, width = 13)
p1 <- DimPlot(sr.sheep.epithelial, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = cols.gp)
p2 <- DimPlot(sr.sheep.epithelial, reduction = "umap", group.by = c("CellType"),
              label = TRUE, repel = TRUE, pt.size = 0.5, cols = pd.col)
p1 + p2
dev.off()
#
pd.gene <- c(# luminal progenitor
  "KIT", "ALDH1A3", "CD14",
  # Luminal cells
  "KRT19", "KRT18", "KRT8",
  # HS cells
  "PRLR", "ESR1", "PROM1",
  # AV cells
  "ENSOARG00020009476", "LTF", "ELF5",
  # myoepithelial lineage
  "KRT17", "KRT5", "MYLK")
pd.cell <- c("Luminal progenitor","Luminal cells","HS-AV cells","AV cells","HS cells","Myoepithelial cells")
pdf(file.path(res.out, "Dot_plot_to_show_marker_gene_expression_with_annotation.pdf"), height = 5, width = 10)
DotPlot(sr.sheep.epithelial, features = pd.gene, cols = c("#eaeaea", "#fc0330"),
        col.min = -0.5, col.max = 1, group.by = "CellType",
        scale = T, scale.min = 0, scale.max = 50) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  scale_y_discrete(limit = pd.cell) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
pdf(file.path(res.out, "Scatter_plot_to_show_marker_gene_expression_with_annotation.pdf"), height = 26.5, width = 30)
FeaturePlot(sr.sheep.epithelial, features = c("EPCAM", pd.gene), cols = c("#E7E7E7", "#1A5A9B"), pt.size = 0.25)
dev.off()
Vis.Annotation.Ratio(sr.meta = sr.sheep.epithelial@meta.data, annotation = c("CellType", "Stage"),
                     pd.title = "CellType by Stage", pd.height = 5, pd.width = 10, pd.col = cols.stage,
                     res.out = file.path(res.out, "Stat_ratio/Cell_type_stat_by_cell_type"))
Vis.Annotation.Ratio(sr.meta = sr.sheep.epithelial@meta.data, annotation = c("Stage", "CellType"),
                     pd.title = "CellType by Stage", pd.height = 5, pd.width = 15, pd.col = pd.col.fix,
                     res.out = file.path(res.out, "Stat_ratio/Cell_type_stat_by_stage"))
Vis.Annotation.Ratio(sr.meta = sr.sheep.epithelial@meta.data, annotation = c("Group", "CellType"),
                     pd.title = "CellType by Sample", pd.height = 5, pd.width = 15, pd.col = pd.col.fix,
                     res.out = file.path(res.out, "Stat_ratio/Cell_type_stat_by_Sample"))



# ==========================
# 6th part: Final annotation ----
# ==========================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Final_annotation")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Merge data
# prepare data
sr.list.anno <- list()
sr.list.anno$data1 <- sr.sheep[, sr.sheep$CellType %in% c("Fibroblast", "Lymphatic endothelial cells",
                                                          "Pericyte", "Smooth muscle cells",
                                                          "Vascular endothelial cells")]
sr.list.anno$data2 <- sr.sheep.immune
sr.list.anno$data3 <- sr.sheep.epithelial
# merge
sr.sheep.final <- Reduce(function(x, y) merge(x, y), sr.list.anno)

# pca
sr.sheep.final <- NormalizeData(sr.sheep.final) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.sheep.final), verbose = FALSE)
sr.sheep.final <- RunPCA(sr.sheep.final, features = VariableFeatures(sr.sheep.final), verbose = FALSE)
ElbowPlot(sr.sheep.final)
dev.off()
dim.n <- 12
sr.sheep.final <- RunUMAP(sr.sheep.final, dims = 1:dim.n)
sr.sheep.final <- RunHarmony(sr.sheep.final, group.by.vars = "Group")
sr.sheep.final <- RunUMAP(sr.sheep.final, reduction = "harmony", dims = 1:dim.n)
pd.col.fix <- pd.col[sample(1:length(pd.col), length(pd.col), replace = F)]
pd.col.fix <- pd.col.fix[1:length(unique(sr.sheep.final$CellType))]
names(pd.col.fix) <- sort(unique(sr.sheep.final$CellType))
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_with_annotation_after_correction.pdf"), height = 5, width = 14)
p1 <- DimPlot(sr.sheep.final, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.05, cols = cols.gp)
p2 <- DimPlot(sr.sheep.final, reduction = "umap", group.by = c("CellType"),
              label = TRUE, repel = TRUE, pt.size = 0.05, cols = pd.col.fix)
p1 + p2
dev.off()
cols.stage <- cols.gp[c(2, 4)]
names(cols.stage) <- c("LL", "PL")
# modify cell types
table(sr.sheep.final$CellType)
sr.sheep.final@meta.data <- sr.sheep.final@meta.data %>% 
  mutate(Cluster = case_when(CellType %in% c("B cells", "CD4+ T cells", "CD8+ T cells", "Dendritic cells",
                                             "Ma Macrophage", "Mast cells", "Mb Macrophage",
                                             "Neutrophil", "NK cells") ~ "Immune cells",
                             CellType %in% c("AV cells", "HS cells", "HS-AV cells", "Luminal cells", 
                                             "Luminal progenitor", "Myoepithelial cells") ~ "Epithelial cells",
                             TRUE ~ CellType))
table(sr.sheep.final@meta.data$Cluster)
pdf(file.path(res.out, "Scatter_plot_to_show_clustering_with_annotation2_after_correction.pdf"), height = 5, width = 14)
p1 <- DimPlot(sr.sheep.final, reduction = "umap", group.by = c("Group"),
              label = TRUE, repel = TRUE, pt.size = 0.05, cols = cols.gp)
p2 <- DimPlot(sr.sheep.final, reduction = "umap", group.by = c("Cluster"),
              label = TRUE, repel = TRUE, pt.size = 0.05, cols = pd.col)
p1 + p2
dev.off()
# find markers (cell type)
Idents(sr.sheep.final) <- sr.sheep.final$CellType
sr.sheep.final$Stage <- gsub("-rep.", "", sr.sheep.final$Group)
table(sr.sheep.final$Stage)
markers.final <- FindAllMarkers(sr.sheep.final, assay = "RNA", only.pos = T, logfc.threshold = log2(1.25))
write.csv(markers.final, file.path(res.out, "Marker_genes_of_final_cell_types.csv"), row.names = F)
pd.gene <- markers.final[-grep("ENSOARG", markers.final$gene), ] %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC)
# plot heatmap (cell type)
tmp <- subset(sr.sheep.final, downsample = 100)
tmp <- ScaleData(tmp, features = rownames(tmp))
library(SeuratWrappers)
tmp <- SeuratWrappers::RunALRA(tmp)
tmp <- ScaleData(tmp, features = rownames(tmp))
DefaultAssay(tmp) <- "RNA"
pdf(file.path(res.out, "Heatmap_to_show_top5_marker_gene_expression.pdf"), height = 15, width = 20)
DoHeatmap(tmp, features = unique(pd.gene$gene))
dev.off()
DefaultAssay(tmp) <- "alra"
pdf(file.path(res.out, "Heatmap_to_show_top5_marker_gene_expression_ALRA.pdf"), height = 15, width = 20)
DoHeatmap(tmp, features = unique(pd.gene$gene))
dev.off()
# marker gene GO (cell type)
markers.go <- list()
for (type in unique(as.character(markers$cluster))) {
  pd.gene <- markers.final %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(cluster == type)
  markers.go[[type]] <- Pipe.GO(species = "human", genelist = pd.gene$gene, basename = type,
                                genetype = "SYMBOL", res.out = file.path(sr.out, "Marker_GO"))
}
# find markers (cluster)
Idents(sr.sheep.final) <- sr.sheep.final$Cluster
markers.final.2 <- FindAllMarkers(sr.sheep.final, assay = "RNA", only.pos = T, logfc.threshold = log2(1.25))
write.csv(markers.final.2, file.path(res.out, "Marker_genes_of_final_clusters.csv"), row.names = F)
pd.gene <- markers.final.2[-grep("ENSOARG", markers.final.2$gene), ] %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC)
# plot heatmap (cluster)
tmp <- subset(sr.sheep.final, downsample = 100)
tmp <- ScaleData(tmp, features = rownames(tmp))
library(SeuratWrappers)
tmp <- SeuratWrappers::RunALRA(tmp)
tmp <- ScaleData(tmp, features = rownames(tmp))
DefaultAssay(tmp) <- "RNA"
pdf(file.path(res.out, "Heatmap_to_show_top5_marker_gene_expression_cluster.pdf"), height = 10, width = 20)
DoHeatmap(tmp, features = unique(pd.gene$gene))
dev.off()
DefaultAssay(tmp) <- "alra"
pdf(file.path(res.out, "Heatmap_to_show_top5_marker_gene_expression_cluster_ALRA.pdf"), height = 10, width = 20)
DoHeatmap(tmp, features = unique(pd.gene$gene))
dev.off()
# plot cell composition (cell type)
sr.sheep.final@meta.data %>% 
  group_by(Group, CellType) %>% 
  summarise(Number = length(CellType)) %>%
  group_by(Group) %>%
  summarise(CellType = CellType,
            Number = Number,
            Percentage = Number/sum(Number)) -> pd
write.csv(pd, file.path(res.out, "Stat_ratio_cell_type_stat_by_sample.csv"), row.names = F)
sr.sheep.final@meta.data %>% 
  group_by(Stage, CellType) %>% 
  summarise(Number = length(CellType)) %>%
  group_by(Stage) %>%
  summarise(CellType = CellType,
            Number = Number,
            Percentage = Number/sum(Number)) -> pd
write.csv(pd, file.path(res.out, "Stat_ratio_cell_type_stat_by_stage.csv"), row.names = F)
Vis.Annotation.Ratio(sr.meta = sr.sheep.final@meta.data, annotation = c("CellType", "Stage"),
                     pd.title = "CellType by Stage", pd.height = 5, pd.width = 10, pd.col = cols.stage,
                     res.out = file.path(res.out, "Stat_ratio/Cell_type_stat_by_cell_type"))
Vis.Annotation.Ratio(sr.meta = sr.sheep.final@meta.data, annotation = c("Stage", "CellType"),
                     pd.title = "CellType by Stage", pd.height = 5, pd.width = 15, pd.col = pd.col.fix,
                     res.out = file.path(res.out, "Stat_ratio/Cell_type_stat_by_stage"))
Vis.Annotation.Ratio(sr.meta = sr.sheep.final@meta.data, annotation = c("Group", "CellType"),
                     pd.title = "CellType by Sample", pd.height = 5, pd.width = 15, pd.col = pd.col.fix,
                     res.out = file.path(res.out, "Stat_ratio/Cell_type_stat_by_Sample"))
# plot cell composition (cluster)
sr.sheep.final@meta.data %>% 
  group_by(Group, Cluster) %>% 
  summarise(Number = length(Cluster)) %>%
  group_by(Group) %>%
  summarise(Cluster = Cluster,
            Number = Number,
            Percentage = Number/sum(Number)) -> pd
write.csv(pd, file.path(res.out, "Stat_ratio_cluster_stat_by_sample.csv"), row.names = F)
sr.sheep.final@meta.data %>% 
  group_by(Stage, CellType) %>% 
  summarise(Number = length(CellType)) %>%
  group_by(Stage) %>%
  summarise(CellType = CellType,
            Number = Number,
            Percentage = Number/sum(Number)) -> pd
write.csv(pd, file.path(res.out, "Stat_ratio_cluster_stat_by_stage.csv"), row.names = F)
Vis.Annotation.Ratio(sr.meta = sr.sheep.final@meta.data, annotation = c("Cluster", "Stage"),
                     pd.title = "Cluster by Stage", pd.height = 5, pd.width = 10, pd.col = cols.stage,
                     res.out = file.path(res.out, "Stat_ratio/Stage_stat_by_Cluster"))
Vis.Annotation.Ratio(sr.meta = sr.sheep.final@meta.data, annotation = c("Stage", "Cluster"),
                     pd.title = "Cluster by Stage", pd.height = 5, pd.width = 15, pd.col = pd.col,
                     res.out = file.path(res.out, "Stat_ratio/Cluster_stat_by_stage"))
Vis.Annotation.Ratio(sr.meta = sr.sheep.final@meta.data, annotation = c("Group", "Cluster"),
                     pd.title = "Cluster by Sample", pd.height = 5, pd.width = 15, pd.col = pd.col,
                     res.out = file.path(res.out, "Stat_ratio/Cluster_stat_by_Sample"))
# plot cell relationship (cell type)
Vis.Annotation.Relationship(sr.meta = sr.sheep.final@meta.data,
                            annotation = c("Stage", "CellType"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = FALSE,
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(res.out, "Stat_relation/Cell_type_stat_by_Stage"))
Vis.Annotation.Relationship(sr.meta = sr.sheep.final@meta.data,
                            annotation = c("Group", "CellType"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = FALSE,
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(res.out, "Stat_relation/Cell_type_stat_by_Sample"))
# plot cell relationship (cluster)
Vis.Annotation.Relationship(sr.meta = sr.sheep.final@meta.data,
                            annotation = c("Stage", "Cluster"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = FALSE,
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(res.out, "Stat_relation/Cluster_stat_by_Stage"))
Vis.Annotation.Relationship(sr.meta = sr.sheep.final@meta.data,
                            annotation = c("Group", "Cluster"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = FALSE,
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(res.out, "Stat_relation/Cluster_stat_by_Sample"))
# 
sr.sheep.final$Replicate <- str_split_fixed(sr.sheep.final$Group, "-", 2)[,2]
pd <- sr.sheep.final@meta.data %>% 
  mutate(cell.type = CellType)
RatioPlot(meta.data = pd, group = "Stage", rep = "Replicate", cell = "CellType", pt.col = pd.col)
# plot marker genes
pd.gene <- c(
  # Lymphatic endothelial cells
  "MMRN1", "PROX1", "FLT4",

  # Vascular endothelial cells
  "SOX17", "SELE",

  # Pericytes
  "RGS5", "DES", "NOTCH3",

  # Smooth muscle cells
  "ACTA2", "MYL9", "MYLK", "MYH11", "KRT5", "TAGLN",

  # Fibroblast
  "COL1A1", "COL1A2", "COL3A1", "FN1",

  # Immune cell
  "PTPRC",
  # Ma Macrophage
  "C1QA", "C1QB", "C1QC", "CSF1R", "CD163",
  # Mb Macrophage
  "MMP12", "SPIC",
  # Dendritic cells
  "TRAF1", "FLT3", "FABP4", "CCL22", "CCL17",
  # B cells
  "CD79A", "CD79B", "CD19", "MS4A1",
  # Neutrophil
  "CSF3R", "ALPL", "S100A12", "MGAM",
  # Mast cells
  "CPA3", "MS4A2", "KIT", "TPSB2",
  # Lymphocytes
  "CD3D", "CD3E", "CD3G",
  # T cells
  "CD4", "CD8A",
  # NK cells
  "GZMA", "AFAP1L2", "SH2D1B",

  # Epithelial cells
  "EPCAM",
  # luminal progenitor
  "KIT", "ALDH1A3", "CD14",
  # Luminal cells
  "KRT19", "KRT18", "KRT8",
  # HS cells
  "PRLR", "ESR1", "PROM1",
  # AV cells
  "ENSOARG00020009476", "LTF", "ELF5",
  # myoepithelial lineage
  "KRT17", "KRT5", "MYLK"
) %>% unique()
pd.cell <- c("Lymphatic endothelial cells", "Vascular endothelial cells", "Pericyte",
             "Smooth muscle cells", "Fibroblast", "Ma Macrophage", "Mb Macrophage",
             "Dendritic cells", "B cells", "Neutrophil", "Mast cells", "CD4+ T cells", "CD8+ T cells",
             "NK cells", "Luminal progenitor", "Luminal cells", "HS cells", "AV cells", "HS-AV cells",
             "Myoepithelial cells")
pdf(file.path(res.out, "Dot_plot_to_show_marker_gene_expression_with_annotation.pdf"), height = 7.5, width = 20)
DotPlot(sr.sheep.final, features = unique(pd.gene), cols = c("#eaeaea", "#fc0330"),
        col.min = -0.5, col.max = 1, group.by = "CellType",
        scale = T, scale.min = 0, scale.max = 50) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  scale_y_discrete(limit = pd.cell) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
pd.cell <- c("Lymphatic endothelial cells", "Vascular endothelial cells", "Pericyte",
             "Smooth muscle cells", "Fibroblast", "Immune cells", "Epithelial cells")
pdf(file.path(res.out, "Dot_plot_to_show_marker_gene_expression_with_cluster_annotation.pdf"), height = 5, width = 20)
DotPlot(sr.sheep.final, features = unique(pd.gene), cols = c("#eaeaea", "#fc0330"),
        col.min = -0.5, col.max = 1, group.by = "Cluster",
        scale = T, scale.min = 0, scale.max = 50) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(pd.gene)) +
  scale_y_discrete(limit = pd.cell) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
pd.gene <- c(
  # AV cells
  "ENSOARG00020009476", "LTF", "ELF5",
  # Immune cell
  "PTPRC",
  # B cells
  "CD79A", "CD79B", "CD19", "MS4A1",
  # Lymphocytes
  "CD3D", "CD3E", "CD3G",
  # T cells
  "CD4", "CD8A",
  # Dendritic cells
  "TRAF1", "FLT3", "FABP4", "CCL22", "CCL17",
  # Fibroblast
  "COL1A1", "COL1A2", "COL3A1", "FN1",
  # Epithelial cells
  "EPCAM",
  # HS cells
  "PRLR", "ESR1", "PROM1",
  # Luminal cells
  "KRT19", "KRT18", "KRT8",
  # luminal progenitor
  "KIT", "ALDH1A3", "CD14",
  # Lymphatic endothelial cells
  "MMRN1", "PROX1", "FLT4",
  # Ma Macrophage
  "C1QA", "C1QB", "C1QC", "CSF1R", "CD163",
  # Mast cells
  "CPA3", "MS4A2", "KIT", "TPSB2",
  # Mb Macrophage
  "MMP12", "SPIC",
  # myoepithelial lineage
  "KRT17", "KRT5", "MYLK",
  # Neutrophil
  "CSF3R", "ALPL", "S100A12", "MGAM",
  # NK cells
  "GZMA", "AFAP1L2", "SH2D1B",
  # Pericytes
  "RGS5", "DES", "NOTCH3",
  # Smooth muscle cells
  "ACTA2", "MYL9", "MYLK", "MYH11", "KRT5", "TAGLN",
  # Vascular endothelial cells
  "SOX17", "SELE"
) %>% unique()
DefaultAssay(tmp) <- "RNA"
pdf(file.path(res.out, "Heatmap_to_show_cell_marker_gene_expression.pdf"), height = 15, width = 20)
DoHeatmap(tmp, features = unique(pd.gene))
dev.off()
DefaultAssay(tmp) <- "alra"
pdf(file.path(res.out, "Heatmap_to_show_cell_marker_gene_expression_ALRA.pdf"), height = 15, width = 20)
DoHeatmap(tmp, features = unique(pd.gene))
dev.off()
DefaultAssay(tmp) <- "RNA"
pdf(file.path(res.out, "Heatmap_to_show_cell_marker_gene_expression_cluster.pdf"), height = 10, width = 20)
DoHeatmap(tmp, features = unique(pd.gene))
dev.off()
DefaultAssay(tmp) <- "alra"
pdf(file.path(res.out, "Heatmap_to_show_cell_marker_gene_expression_cluster_ALRA.pdf"), height = 10, width = 20)
DoHeatmap(tmp, features = unique(pd.gene))
dev.off()
# plot expression
marker.gene <- list(cell.marker = unique(pd.gene))
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.sheep.final))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(res.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 5, width = 5.5)
    print(FeaturePlot(sr.sheep.final, features = marker.gene[[type]][g], order = T,
                      cols = c("#E7E7E7", "#1A5A9B"), pt.size = 0.1, slot = "data", min.cutoff = 0.1))
    dev.off()
  }; rm(g)
}


### >>> 3. Differential expression analysis
deg.final <- list()
for (i in unique(sr.sheep.final$CellType)) {
  tmp <- subset(sr.sheep.final, CellType == i)
  test <- SeuratToEdger(sr.obj = tmp, group.by = "Stage",
                        comparison = c("LL", "PL"), sample.n = 300)
  print(table(tmp$CellType, tmp$Stage))
  tmp <- gsub(" ", "_", i)
  deg.final[[tmp]] <- edgeR.scRNA(count = test$count, meta = test$meta, expr = test$expr,
                                  sample.n = NULL, g1 = "1_LL", g2 = "2_PL",
                                  lfc = log2(2), sig = 0.05,
                                  res.out = file.path(res.out, paste0("DEGs/PL_vs_LL/", tmp)))
}


### >>> 4. Plotting
sheep.tfs <- read.table("/home/yhw/document/AnimalTFDBs/Ovis_aries_TF.txt", header = T, sep = "\t")
sheep.tfs$Symbol[sheep.tfs$Symbol == "-"] <- sheep.tfs$Ensembl[sheep.tfs$Symbol == "-"]
for (i in names(deg.final)) {
  # all
  write.csv(subset(deg.final[[i]]$all, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t1_2_PL_Vs_1_LL_all_TFs_results.csv")),
            row.names = F)
  pdf(file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/Volcano_plot_to_show_all_genes.pdf")),
      height = 5, width = 7)
  print(VisDEG.volcano(deg.data = deg.final[[i]]$all, geneset = NULL,
                       p.col = "PValue", lfc.col = "logFC", sig = 0.05, lfc = 1,
                       title = paste0(i, ": all genes"), up.col = "#B71C1C", down.col = "#01579B"))
  dev.off()
  pdf(file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/Volcano_plot_to_show_all_TFs.pdf")),
      height = 5, width = 7)
  print(VisDEG.volcano(deg.data = deg.final[[i]]$all, geneset = sheep.tfs$Symbol,
                       p.col = "PValue", lfc.col = "logFC", sig = 0.05, lfc = 1,
                       title = paste0(i, ": all TFs"), up.col = "#B71C1C", down.col = "#01579B"))
  dev.off()
  # all sig
  write.csv(subset(deg.final[[i]]$sig, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t2_2_PL_Vs_1_LL_all_sig_TFs_results.csv")),
            row.names = F)
  # up sig
  write.csv(subset(deg.final[[i]]$up, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t3_2_PL_Vs_1_LL_up_sig_TFs_lfc1_sig0.05.csv")),
            row.names = F)
  # down sig
  write.csv(subset(deg.final[[i]]$down, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t4_2_PL_Vs_1_LL_down_sig_TFs_lfc1_sig0.05.csv")),
            row.names = F)

}


### >>> 5. GO and GSEA analysis
go.res <- list()
gsea.res <- list()
for (i in names(deg.final)) {
  gsea.res[[i]] <- Pipe.GSEA(deg.obj = deg.final[[i]]$sig[-grep("ENSOARG", deg.final[[i]]$sig$SYMBOL), ],
                             deg.type = "edger", lfc = log2(1.5), sig = 0.05,
                             species = "human", basename = i,
                             genetype = "SYMBOL", gene.col = "SYMBOL",
                             outdir = file.path(res.out, paste0("GSEA/PL_vs_LL/", i)))
  go.res[[paste0(i, ".up")]] <- Pipe.GO(species = "human",
                                        genelist = deg.final[[i]]$up$SYMBOL,
                                        basename = paste0(i, "_up"),
                                        genetype = "SYMBOL",
                                        res.out = file.path(res.out, paste0("GO/PL_vs_LL/", i)))
  go.res[[paste0(i, ".down")]] <- Pipe.GO(species = "human",
                                          genelist = deg.final[[i]]$down$SYMBOL,
                                          basename = paste0(i, "_down"),
                                          genetype = "SYMBOL",
                                          res.out = file.path(res.out, paste0("GO/PL_vs_LL/", i)))
}
library(aPEAR)
for (i in names(go.res)) {
  for (j in names(go.res$Fibroblast.up)) {
    pd <- na.omit(subset(go.res[[i]][[j]], pvalue <= 0.05)[1:20, ])
    if (nrow(pd) >= 5) {
      pdf(file.path(res.out, paste0("GO/Network_plot_of_", i, "_", j, ".pdf")), height = 10, width = 10)
      print(enrichmentNetwork(pd, colorBy = "pvalue", colorType = "pval", pCutoff = -5, drawEllipses = TRUE))
      dev.off()
    }
  }
}



# =============================
# 7th part: Trajectory analysis ----
# =============================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Trajectory_analysis")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. URD analysis
pd.gene <- c("EPCAM","KIT","ALDH1A3","CD14","KRT19","KRT18","KRT8","PRLR","ESR1",
             "PROM1","ENSOARG00020009476","LTF","ELF5","KRT17","KRT5","MYLK") %>% unique()
urd <- list()
tree <- list()
# all epithelial cells
epi.reduction <- Embeddings(sr.sheep.epithelial, reduction = "umap") %>% as.data.frame()
urd$all.epithelial <- URD.Btree(sr.obj = sr.sheep.epithelial, 
                                group.by = c("CellType", "Stage"), 
                                tip.gp = c("AV cells", "HS cells", "HS-AV cells"), 
                                start.cell = "Luminal progenitor", 
                                urd.knn = 86, urd.sigma = NULL, 
                                embedding = epi.reduction,
                                group.col = list(CellType = pd.col.fix[names(pd.col.fix) %in% unique(sr.sheep.epithelial$CellType)],
                                                 Stage = cols.stage),
                                pt.gene = pd.gene, 
                                expr.col = c("#ebebeb", "#7d3c98"), 
                                pt.size = 1, 
                                res.out = file.path(res.out, "URD/all_epithelial"))
tree$all.epithelial <- URD.Atree(urd = urd$all.epithelial, 
                                 tip.gp = c("AV cells", "HS cells", "HS-AV cells"), 
                                 tip.cluster = c(2, 6, 5), 
                                 group.by = c("CellType", "Stage"), 
                                 group.col = list(CellType = pd.col.fix[names(pd.col.fix) %in% unique(sr.sheep.epithelial$CellType)],
                                                  Stage = cols.stage),
                                 expr.col = c("#ebebeb", "#7d3c98"), 
                                 pt.gene = pd.gene, 
                                 pt.size = 1, 
                                 res.out = file.path(res.out, "URD/all_epithelial"))
# only luminal cells
epi.reduction <- Embeddings(subset(sr.sheep.epithelial, CellType != "Myoepithelial cells"), reduction = "umap") %>% as.data.frame()
urd$only.luminal <- URD.Btree(sr.obj = subset(sr.sheep.epithelial, CellType != "Myoepithelial cells"), 
                              group.by = c("CellType", "Stage"), 
                              tip.gp = c("AV cells", "HS cells", "HS-AV cells"), 
                              start.cell = "Luminal progenitor", 
                              urd.knn = 79, urd.sigma = NULL, 
                              embedding = epi.reduction,
                              group.col = list(CellType = pd.col.fix[names(pd.col.fix) %in% unique(sr.sheep.epithelial$CellType)],
                                               Stage = cols.stage),
                              pt.gene = pd.gene, 
                              expr.col = c("#ebebeb", "#7d3c98"), 
                              pt.size = 1, 
                              res.out = file.path(res.out, "URD/only_luminal"))
tree$only.luminal <- URD.Atree(urd = urd$only.luminal, 
                               tip.gp = c("AV cells", "HS cells", "HS-AV cells"), 
                               tip.cluster = c(2, 6, 5), 
                               group.by = c("CellType", "Stage"), 
                               group.col = list(CellType = pd.col.fix[names(pd.col.fix) %in% unique(sr.sheep.epithelial$CellType)],
                                                Stage = cols.stage),
                               expr.col = c("#ebebeb", "#7d3c98"), 
                               pt.gene = pd.gene, 
                               pt.size = 1, 
                               res.out = file.path(res.out, "URD/only_luminal"))


### >>> 3. Monocle2 analysis
mono2 <- list()
# all epithelial
mono2$all.epithelial <- TI.Monocle2(sr.obj = sr.sheep.epithelial, 
                                    outdir = file.path(res.out, "monocle2/all_epithelial"), 
                                    sp.name = "all_epithelial", 
                                    gp.name = "CellType", 
                                    pd.color = pd.col.fix[names(pd.col.fix) %in% unique(sr.sheep.epithelial$CellType)], 
                                    point.size = 0.75)
mono2$all_epithelial <- TI.Monocle2.Vis(mn.obj = mono2$all.epithelial$data, 
                                        sp.name = "all_epithelial", 
                                        root = 3, 
                                        branch = 1, 
                                        useless.gene = "ENSOARG", 
                                        point.size = 1,
                                        outdir = file.path(res.out, "monocle2/all_epithelial"))
mono2$all_epithelial_go <- TI.Monocle2.GO(mn.heat = mono2$all_epithelial$plot, 
                                          species = "human", sp.name = "all_epithelial", 
                                          genetype = "SYMBOL", 
                                          outdir = file.path(res.out, "monocle2/all_epithelial/GO"))
# only luminal
mono2$only_luminal <- TI.Monocle2(sr.obj = subset(sr.sheep.epithelial, CellType != "Myoepithelial cells"), 
                                  outdir = file.path(res.out, "monocle2/only_luminal"), 
                                  sp.name = "only_luminal", 
                                  gp.name = "CellType", 
                                  pd.color = pd.col.fix[names(pd.col.fix) %in% unique(sr.sheep.epithelial$CellType)], 
                                  point.size = 0.75)
mono2$only_luminal <- TI.Monocle2.Vis(mn.obj = mono2$only_luminal$data, 
                                      sp.name = "only_luminal", 
                                      root = 2, 
                                      branch = 1, 
                                      useless.gene = "ENSOARG", 
                                      point.size = 1,
                                      outdir = file.path(res.out, "monocle2/only_luminal"))
mono2$only_luminal_go <- TI.Monocle2.GO(mn.heat = mono2$only_luminal$plot, 
                                        species = "human", sp.name = "only_luminal", 
                                        genetype = "SYMBOL", 
                                        outdir = file.path(res.out, "monocle2/only_luminal/GO"))


### >>> 4. RNA velocity
# all epithelial
for (i in unique(sr.sheep.epithelial$Group)) {
  SeuratToVelocyto(sr.ob = subset(sr.sheep.epithelial, Group == i), 
                   reduction = "umap", 
                   out.dir = file.path(res.out, "Velocyte/all_epithelial"), 
                   sample.name = i, 
                   str.in.barcode = "Seurat.*rep._")
}
# only luminal
for (i in unique(sr.sheep.epithelial$Group)) {
  SeuratToVelocyto(sr.ob = subset(sr.sheep.epithelial, Group == i & CellType != "Myoepithelial cells"), 
                   reduction = "umap", 
                   out.dir = file.path(res.out, "Velocyte/only_luminal"), 
                   sample.name = i, 
                   str.in.barcode = "Seurat.*rep._")
}



# ============================
# 8th part: Cell Chat analysis ----
# ============================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Cell_communication")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Epithelial cells
library(AfterChat)
library(CellChat)
db.cellchat.hs <- CellChatDB.human
ct.epi <- list()
ct.epi <- Pipe.CellChat(sc.obj = sr.sheep.epithelial, cell.db = db.cellchat.hs,
                        group.by = "CellType", split.by = "Stage")
for (i in names(ct.epi)) {
  ct.epi[[i]] <- computeCommunProb(ct.epi[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
for (i in names(ct.epi)[2]) {
  AfterChat::PathCentrality(ct.obj = ct.epi[[i]], outdir = file.path(res.out, "Epithelial"), file.prefix = i)
  AfterChat::PathInteracion(ct.obj = ct.epi[[i]], outdir = file.path(res.out, "Epithelial"), file.prefix = i)
  AfterChat::LRsContribution(ct.obj = ct.epi[[i]], outdir = file.path(res.out, "Epithelial"), file.prefix = i)
  #LRsInteraction(ct.obj = ct.epi[[i]], outdir = file.path(res.out, "Epithelial"), file.prefix = i,
  #               cell.source = c(), cell.target = c())
  AfterChat::PathClustering(ct.obj = ct.epi[[i]], outdir = file.path(res.out, "Epithelial"), file.prefix = i)
}


### >>> 3. Immune cells
ct.immune <- list()
sr.sheep.immune$Stage <- gsub("-rep.", "", sr.sheep.immune$Group)
ct.immune <- Pipe.CellChat(sc.obj = sr.sheep.immune, cell.db = db.cellchat.hs,
                           group.by = "CellType", split.by = "Stage")
for (i in names(ct.immune)) {
  ct.immune[[i]] <- computeCommunProb(ct.immune[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
for (i in names(ct.immune)) {
  AfterChat::PathCentrality(ct.obj = ct.immune[[i]], outdir = file.path(res.out, "Immune"), file.prefix = i)
  AfterChat::PathInteracion(ct.obj = ct.immune[[i]], outdir = file.path(res.out, "Immune"), file.prefix = i)
  AfterChat::LRsContribution(ct.obj = ct.immune[[i]], outdir = file.path(res.out, "Immune"), file.prefix = i)
  #LRsInteraction(ct.obj = ct.immune[[i]], outdir = file.path(res.out, "Immune"), file.prefix = i,
  #               cell.source = c(), cell.target = c())
  AfterChat::PathClustering(ct.obj = ct.immune[[i]], outdir = file.path(res.out, "Immune"), file.prefix = i)
}


### >>> 4. Others
ct.others <- list()
sr.sheep.others <- subset(sr.sheep.final, Cluster %in% c("Fibroblast", "Lymphatic endothelial cells",
                                                         "Pericyte", "Smooth muscle cells", "Vascular endothelial cells"))

ct.others <- Pipe.CellChat(sc.obj = sr.sheep.others, cell.db = db.cellchat.hs,
                           group.by = "CellType", split.by = "Stage")
for (i in names(ct.others)) {
  ct.others[[i]] <- computeCommunProb(ct.others[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
for (i in names(ct.others)) {
  AfterChat::PathCentrality(ct.obj = ct.others[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
  AfterChat::PathInteracion(ct.obj = ct.others[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
  AfterChat::LRsContribution(ct.obj = ct.others[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
  #LRsInteraction(ct.obj = ct.others[[i]], outdir = file.path(res.out, "othersthelial"), file.prefix = i,
  #               cell.source = c(), cell.target = c())
  AfterChat::PathClustering(ct.obj = ct.others[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
}


### >>> 5. Immune to stromal cells
ct.cross1 <- list()
sr.sheep.others <- subset(sr.sheep.final, Cluster %in% c("Immune cells",
                                                         "Fibroblast", "Lymphatic endothelial cells",
                                                         "Pericyte", "Smooth muscle cells", "Vascular endothelial cells"))
ct.cross1 <- Pipe.CellChat(sc.obj = sr.sheep.others, cell.db = db.cellchat.hs,
                           group.by = "CellType", split.by = "Stage")
for (i in names(ct.cross1)) {
  ct.cross1[[i]] <- computeCommunProb(ct.cross1[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
table(sr.sheep.others$CellType)
for (i in names(ct.cross1)) {
  AfterChat::PathCentrality(ct.obj = ct.cross1[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
  AfterChat::PathInteracion(ct.obj = ct.cross1[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
  AfterChat::LRsContribution(ct.obj = ct.cross1[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
  LRsInteraction(ct.obj = ct.cross1[[i]], outdir = file.path(res.out, "Others"), file.prefix = i,
                 cell.source = c("B cells", "CD4+ T cells", "CD8+ T cells", "Dendritic cells",
                                 "Ma Macrophage", "Mast cells", "Mb Macrophage",
                                 "Neutrophil", "NK cells"), 
                 cell.target = c("Fibroblast", "Lymphatic endothelial cells",
                                 "Pericyte", "Smooth muscle cells", "Vascular endothelial cells"))
  AfterChat::PathClustering(ct.obj = ct.cross1[[i]], outdir = file.path(res.out, "Others"), file.prefix = i)
}


### >>> 6. Immune to epithelial cells
library(CellChat)
library(AfterChat)
res.out <- file.path(getwd(), "R/Graphs/Cell_communication")
ct.cross2 <- list()
sr.sheep.others <- subset(sr.sheep.final, Cluster %in% c("Immune cells", "Epithelial cells"))
table(sr.sheep.others$CellType)
sr.sheep.others@meta.data <- sr.sheep.others@meta.data %>% 
  mutate(CellType2 = case_when(CellType == "AV cells" ~ "AV",
                               CellType == "HS cells" ~ "HS",
                               CellType == "HS-AV cells" ~ "HS-AV",
                               CellType == "Luminal cells" ~ "Luminal",
                               CellType == "Luminal progenitor" ~ "Prog",
                               CellType == "Myoepithelial cells" ~ "Myoepi",
                               CellType == "Dendritic cells" ~ "Dend",
                               CellType == "Neutrophil" ~ "Neut",
                               TRUE ~ CellType))
sr.sheep.others$CellType2 <- gsub(" cells", "", sr.sheep.others$CellType2)
sr.sheep.others$CellType2 <- gsub(" Macrophage", "", sr.sheep.others$CellType2)
table(sr.sheep.others$CellType2)
ct.cross2 <- Pipe.CellChat(sc.obj = sr.sheep.others, cell.db = db.cellchat.hs,
                           group.by = "CellType2", split.by = "Stage")
for (i in names(ct.cross2)) {
  ct.cross2[[i]] <- computeCommunProb(ct.cross2[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
for (i in names(ct.cross2)) {
  AfterChat::PathCentrality(ct.obj = ct.cross2[[i]], outdir = file.path(res.out, "Immune_to_Epithelial"), file.prefix = i)
  AfterChat::PathInteracion(ct.obj = ct.cross2[[i]], outdir = file.path(res.out, "Immune_to_Epithelial"), file.prefix = i)
  AfterChat::LRsContribution(ct.obj = ct.cross2[[i]], outdir = file.path(res.out, "Immune_to_Epithelial"), file.prefix = i)
  LRsInteraction(ct.obj = ct.cross2[[i]], outdir = file.path(res.out, "Immune_to_Epithelial"), file.prefix = i,
                 cell.source = c("B", "CD4+ T", "CD8+ T", "Dendritic", "Ma Macrophage", 
                                 "Mast", "Mb Macrophage", "Neutrophil", "NK"), 
                 cell.target = c("AV", "HS", "HS-AV", "Luminal", 
                                 "Progenitor", "Myoepithelial"))
  AfterChat::PathClustering(ct.obj = ct.cross2[[i]], outdir = file.path(res.out, "Immune_to_Epithelial"), file.prefix = i)
}


### >>> 6. Stromal cells to epithelial cells
library(CellChat)
library(AfterChat)
sr.sheep.others <- subset(sr.sheep.final, Cluster %in% c("Fibroblast", "Lymphatic endothelial cells",
                                                         "Pericyte", "Smooth muscle cells", "Vascular endothelial cells", 
                                                         "Epithelial cells"))
table(sr.sheep.others$CellType)
sr.sheep.others@meta.data <- sr.sheep.others@meta.data %>% 
  mutate(CellType2 = case_when(CellType == "AV cells" ~ "AV",
                               CellType == "HS cells" ~ "HS",
                               CellType == "HS-AV cells" ~ "HS-AV",
                               CellType == "Luminal cells" ~ "Luminal",
                               CellType == "Luminal progenitor" ~ "Prog",
                               CellType == "Myoepithelial cells" ~ "Myoepi",
                               CellType == "Fibroblast" ~ "Fibro",
                               CellType == "Lymphatic endothelial cells" ~ "LymEndo",
                               CellType == "Smooth muscle cells" ~ "SmoMuscle",
                               CellType == "Vascular endothelial cells" ~ "VasEndo",
                               CellType == "Pericyte" ~ "Pericyte",
                               TRUE ~ CellType))
table(sr.sheep.others$CellType2)
ct.cross3 <- Pipe.CellChat(sc.obj = sr.sheep.others, cell.db = db.cellchat.hs,
                           group.by = "CellType2", split.by = "Stage")
for (i in names(ct.cross3)) {
  ct.cross3[[i]] <- computeCommunProb(ct.cross3[[i]], type = "truncatedMean", trim = 0.1) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}
LRsInteraction.new <- edit(LRsInteraction)
for (i in names(ct.cross3)) {
  AfterChat::PathCentrality(ct.obj = ct.cross3[[i]], outdir = file.path(res.out, "Stromal_to_Epithelial"), file.prefix = i)
  AfterChat::PathInteracion(ct.obj = ct.cross3[[i]], outdir = file.path(res.out, "Stromal_to_Epithelial"), file.prefix = i)
  AfterChat::LRsContribution(ct.obj = ct.cross3[[i]], outdir = file.path(res.out, "Stromal_to_Epithelial"), file.prefix = i)
  LRsInteraction.new(ct.obj = ct.cross3[[i]], outdir = file.path(res.out, "Stromal_to_Epithelial"), file.prefix = i,
                     cell.source = c("Fibro", "LymEndo", "SmoMuscle", "VasEndo", "Pericyte"), 
                     cell.target = c("AV", "HS", "HS-AV", "Luminal", "Prog", "Myoepi"))
  AfterChat::PathClustering(ct.obj = ct.cross3[[i]], outdir = file.path(res.out, "Stromal_to_Epithelial"), file.prefix = i)
}



# ==============================
# 9th part: Integrative analysis ----
# ==============================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Integrative_Analysis")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Load human data
sr.human <- readRDS("/home/yhw/document/public_data/GSE195665_human_breast/scRNA-seq.all.cell.rds")
Idents(sr.human) <- sr.human$author_cell_type
table(sr.human$author_cell_type) %>% as.data.frame()
sr.human <- subset(sr.human, downsample = 1000)
table(sr.human$author_cell_type) %>% as.data.frame()


### >>> 3. Process human data
# filter homologous genes
hs.oar <- read.table("/home/yhw/document/ensembl/Homologous_Genes/release108_GRCh38_Oar_v3.1_mart_export.txt", header = T, sep = "\t")
hs.oar <- subset(hs.oar, Sheep.homology.type == "ortholog_one2one")
keep.gene <- intersect(hs.oar$Gene.stable.ID, rownames(sr.human))
hs.oar <- subset(hs.oar, Gene.stable.ID %in% keep.gene)
keep.gene <- intersect(hs.oar$Sheep.gene.name, rownames(sr.sheep.final))
hs.oar <- subset(hs.oar, Sheep.gene.name %in% keep.gene)
hs.oar <- hs.oar[hs.oar$Gene.name == hs.oar$Sheep.gene.name, ]
rownames(hs.oar) <- hs.oar$Gene.stable.ID
# filter genes in human data
count <- GetAssayData(sr.human, slot = "count")
count <- count[rownames(hs.oar), ]
rownames(count) <- hs.oar$Gene.name
meta <- sr.human@meta.data
hs.ref <- CreateSeuratObject(counts = count, meta.data = meta)
rm(count, meta)
hs.ref <- NormalizeData(hs.ref) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(hs.ref), verbose = FALSE)
hs.ref <- RunPCA(hs.ref, features = VariableFeatures(hs.ref), verbose = FALSE)
ElbowPlot(hs.ref)
dev.off()
hs.ref <- RunUMAP(hs.ref, dims = 1:10)
DimPlot(hs.ref, group.by = "author_cell_type")
unique(hs.ref$author_cell_type)
hs.ref@meta.data <- hs.ref@meta.data %>% 
  mutate(Cluster = case_when(author_cell_type %in% c("LummHR-SCGB", "LummHR-major", "LummHR-active", "Lumsec-HLA", 
                                                     "Lumsec-myo", "Lumsec-prol", "Lumsec-basal", "Lumsec-KIT", 
                                                     "Lumsec-major", "Lumsec-lac", "basal") ~ "Epithelial",
                             author_cell_type %in% c("Fibro-prematrix", "Fibro-matrix", "Fibro-SFRP4", "Fibro-major", 
                                                     "Lymph-immune", "Lymph-major", "Lymph-valve2", "Lymph-valve1", 
                                                     "Vas-capillary", "Vas-venous", "Vas-arterial",
                                                     "vsmc", "pericytes") ~ "Stromal",
                             author_cell_type %in% c("Mast", "b_naive", "bmem_switched", "plasma_IgG",
                                                     "plasma_IgA", "bmem_unswitched", "CD8-Tem", "NKT", "GD", "CD4-naive", "NK",
                                                     "NK-ILCs", "CD4-activated", "CD4-Treg", "CD4-Th", "CD8-activated", "T_prol",
                                                     "CD4-Th-like", "CD4-Tem", "CD8-Trm", "Macro-lipo", "mDC", "Mono-non-classical", 
                                                     "Mono-classical", "cDC2", "cDC1", "Macro-m2", "Macro-IFN", "Macro-m1", 
                                                     "Macro-m1-CCL", "Macro-m2-CXCL", "mye-prol", "pDC", "Neutrophil") ~ "Immune"))
# split data
hs.ref.list <- SplitObject(hs.ref, split.by = "Cluster")
# pca
hs.ref.list <- lapply(X = hs.ref.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
})
p1 <- ElbowPlot(hs.ref.list$Epithelial)
p2 <- ElbowPlot(hs.ref.list$Immune)
p3 <- ElbowPlot(hs.ref.list$Stromal)
plot_grid(p1, p2, p3, ncol = 3)
# cluster cells + UMAP/tSNE
dim.n <- c(15, 15, 15)
names(dim.n) <- names(hs.ref.list)
for (i in names(hs.ref.list)) {
  hs.ref.list[[i]]@project.name <- i
}
hs.ref.list <- lapply(X = hs.ref.list, FUN = function(x) {
  dims <- unname(dim.n[x@project.name])
  x <- RunUMAP(x, dims = 1:dim.n) %>% RunTSNE(dims = 1:dim.n)
})
DimPlot(hs.ref.list$Epithelial, group.by = "author_cell_type", cols = pd.col)
DimPlot(hs.ref.list$Immune, group.by = "author_cell_type", cols = pd.col)
DimPlot(hs.ref.list$Stromal, group.by = "author_cell_type", cols = pd.col)


### >>> 4. Process data
# filter genes in sheep data
count <- GetAssayData(sr.sheep.final, slot = "count")
table(rownames(count) %in% hs.oar$Sheep.gene.name)
count <- count[hs.oar$Sheep.gene.name, ]
meta <- sr.sheep.final@meta.data
query.sheep <- CreateSeuratObject(counts = count, meta.data = meta)
rm(count, meta)
query.sheep <- NormalizeData(query.sheep) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(query.sheep), verbose = FALSE)
query.sheep <- RunPCA(query.sheep, features = VariableFeatures(query.sheep), verbose = FALSE)
ElbowPlot(query.sheep)
dev.off()
query.sheep <- RunUMAP(query.sheep, dims = 1:12)
DimPlot(query.sheep, group.by = "CellType")
query.sheep@meta.data <- query.sheep@meta.data %>% 
  mutate(Cluster2 = case_when(Cluster %in% c("Epithelial cells") ~ "Epithelial",
                              Cluster %in% c("Fibroblast", "Lymphatic endothelial cells", "Pericyte", 
                                             "Smooth muscle cells", "Vascular endothelial cells") ~ "Stromal",
                              Cluster %in% c("Immune cells") ~ "Immune",))
# split data
query.sheep.list <- SplitObject(query.sheep, split.by = "Cluster2")
# pca
query.sheep.list <- lapply(X = query.sheep.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
})
p1 <- ElbowPlot(query.sheep.list$Epithelial)
p2 <- ElbowPlot(query.sheep.list$Immune)
p3 <- ElbowPlot(query.sheep.list$Stromal)
plot_grid(p1, p2, p3, ncol = 3)
# cluster cells + UMAP/tSNE
dim.n <- c(15, 15, 15)
names(dim.n) <- names(query.sheep.list)
for (i in names(query.sheep.list)) {
  query.sheep.list[[i]]@project.name <- i
}
query.sheep.list <- lapply(X = query.sheep.list, FUN = function(x) {
  dims <- unname(dim.n[x@project.name])
  x <- RunUMAP(x, dims = 1:dim.n) %>% RunTSNE(dims = 1:dim.n)
})
DimPlot(query.sheep.list$Epithelial, group.by = "CellType", cols = pd.col)
DimPlot(query.sheep.list$Immune, group.by = "CellType", cols = pd.col)
DimPlot(query.sheep.list$Stromal, group.by = "CellType", cols = pd.col)


### >>> 5. Annotate cells
tmp <- AnnotateCell(ref = hs.ref, query = query.sheep, 
                    gene.set = NULL, group.by = "author_cell_type", dim.n = 30)
table(tmp$query@meta.data$Cluster)
# Epithelial
sheep.anno <- list()
sheep.anno$Epithelial <- AnnotateCell(ref = hs.ref.list$Epithelial, 
                                      query = query.sheep.list$Epithelial,
                                      gene.set = NULL, group.by = "author_cell_type", dim.n = 30)
pdf(file.path(res.out, "Sscatter_to_show_epithelial_annotation_results_by_three_methods.pdf"), height = 14, width = 18)
DimPlot(sheep.anno$Epithelial$query, group.by = c("CellType", "Seurat.id", "scPred.prediction", "singleR.labels"), 
        cols = pd.col, pt.size = 1)
dev.off()
Vis.Annotation.Relationship(sr.meta = sheep.anno$Epithelial$query@meta.data,
                            annotation = c("CellType", "scPred.prediction"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = F, 
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(res.out, "Sankey_plot_of_epithelial_annotation"))
# Immune
sheep.anno$Immune <- AnnotateCell(ref = hs.ref.list$Immune, 
                                  query = query.sheep.list$Immune,
                                  gene.set = NULL, group.by = "author_cell_type", dim.n = 30)
pdf(file.path(res.out, "Sscatter_to_show_immune_annotation_results_by_three_methods.pdf"), height = 14, width = 22)
DimPlot(sheep.anno$Immune$query, group.by = c("CellType", "Seurat.id", "scPred.prediction", "singleR.labels"), 
        cols = pd.col, pt.size = 0.75)
dev.off()
Vis.Annotation.Relationship(sr.meta = sheep.anno$Immune$query@meta.data,
                            annotation = c("CellType", "scPred.prediction", "Seurat.id", "singleR.labels"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = F, 
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(res.out, "Sankey_plot_of_immune_annotation"))
# Stromal
sheep.anno$Stromal <- AnnotateCell(ref = hs.ref.list$Stromal, 
                                   query = query.sheep.list$Stromal,
                                   gene.set = NULL, group.by = "author_cell_type", dim.n = 30)
pdf(file.path(res.out, "Sscatter_to_show_stromal_annotation_results_by_three_methods.pdf"), height = 14, width = 18)
DimPlot(sheep.anno$Stromal$query, group.by = c("CellType", "Seurat.id", "scPred.prediction", "singleR.labels"), 
        cols = pd.col, pt.size = 1)
dev.off()
Vis.Annotation.Relationship(sr.meta = sheep.anno$Stromal$query@meta.data,
                            annotation = c("CellType", "scPred.prediction"),
                            pd.title = "Relationship between various annotations",
                            pd.col = NULL, split.meta = F, 
                            pd.height = 15, pd.width = 15,
                            res.out = file.path(res.out, "Sankey_plot_of_stromal_annotation"))


### >>> 6. Integrate datasets across species
# Epithelial
inte.list <- list(hs.ref.list$Epithelial, sheep.anno$Epithelial$query)
inte.list[[1]]@meta.data <- inte.list[[1]]@meta.data %>% 
  mutate(Species = "Human",
         Group = "Human.Breast",
         Stage = "Mature",
         CellType.final = author_cell_type)
inte.list[[1]]@meta.data <- inte.list[[1]]@meta.data[, c("Group", "Stage", "Species", "Cluster", "CellType.final")]
inte.list[[2]]@meta.data <- inte.list[[2]]@meta.data %>% 
  mutate(Species = "Sheep",
         CellType.final = scPred.prediction)
inte.list[[2]] <- subset(inte.list[[2]], CellType.final != "unassigned")
inte.list[[2]]@meta.data <- inte.list[[2]]@meta.data[, c("Group", "Stage", "Species", "Cluster", "CellType.final")]
anchors <- FindIntegrationAnchors(object.list = inte.list, dims = 1:20)
inte.list$merged <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(inte.list$merged) <- "integrated"
inte.list$merged <- ScaleData(inte.list$merged, verbose = FALSE)
inte.list$merged <- RunPCA(inte.list$merged, npcs = 30, verbose = FALSE)
ElbowPlot(inte.list$merged)
inte.list$merged <- RunUMAP(inte.list$merged, reduction = "pca", dims = 1:15)
DimPlot(inte.list$merged, reduction = "umap", group.by = c("Species", "CellType.final"), 
        cols = pd.col, label = TRUE, repel = TRUE)
inte.list$merged <- RunHarmony(inte.list$merged, group.by.vars = "Species")
inte.list$merged <- RunUMAP(inte.list$merged, reduction = "harmony", dims = 1:15)
pdf(file.path(res.out, "Scatter_to_show_epithelial_annotation_results_by_final_integration.pdf"), height = 6, width = 15)
DimPlot(inte.list$merged, reduction = "umap", group.by = c("Species", "CellType.final"), 
        cols = pd.col, label = TRUE, repel = TRUE, pt.size = 0.75)
dev.off()
pdf(file.path(res.out, "Scatter_to_show_epithelial_annotation_results_by_final_integration_splited.pdf"), height = 6, width = 12)
DimPlot(inte.list$merged, reduction = "umap", group.by = "CellType.final", 
        cols = pd.col, label = TRUE, repel = TRUE, pt.size = 0.75, split.by = "Species")
dev.off()
# Immune
inte.list2 <- list(hs.ref.list$Immune, sheep.anno$Immune$query)
inte.list2[[1]]@meta.data <- inte.list2[[1]]@meta.data %>% 
  mutate(Species = "Human",
         Group = "Human.Breast",
         Stage = "Mature",
         CellType.final = author_cell_type)
inte.list2[[1]]@meta.data <- inte.list2[[1]]@meta.data[, c("Group", "Stage", "Species", "Cluster", "CellType.final")]
inte.list2[[2]]@meta.data <- inte.list2[[2]]@meta.data %>% 
  mutate(Species = "Sheep",
         CellType.final = scPred.prediction)
inte.list2[[2]] <- subset(inte.list2[[2]], CellType.final != "unassigned")
inte.list2[[2]]@meta.data <- inte.list2[[2]]@meta.data[, c("Group", "Stage", "Species", "Cluster", "CellType.final")]
anchors <- FindIntegrationAnchors(object.list = inte.list2, dims = 1:30)
inte.list2$merged <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(inte.list2$merged) <- "integrated"
inte.list2$merged <- ScaleData(inte.list2$merged, verbose = FALSE)
inte.list2$merged <- RunPCA(inte.list2$merged, npcs = 30, verbose = FALSE)
ElbowPlot(inte.list2$merged)
inte.list2$merged <- RunUMAP(inte.list2$merged, reduction = "pca", dims = 1:15)
DimPlot(inte.list2$merged, reduction = "umap", group.by = c("Species", "CellType.final"), 
        cols = pd.col, label = TRUE, repel = TRUE)
inte.list2$merged <- RunHarmony(inte.list2$merged, group.by.vars = "Species")
inte.list2$merged <- RunUMAP(inte.list2$merged, reduction = "harmony", dims = 1:15)
pdf(file.path(res.out, "Scatter_to_show_immune_annotation_results_by_final_integration.pdf"), height = 6, width = 17)
DimPlot(inte.list2$merged, reduction = "umap", group.by = c("Species", "CellType.final"), 
        cols = pd.col, label = TRUE, repel = TRUE, pt.size = 0.25)
dev.off()
pdf(file.path(res.out, "Scatter_to_show_immune_annotation_results_by_final_integration_splited.pdf"), height = 6, width = 15)
DimPlot(inte.list2$merged, reduction = "umap", group.by = "CellType.final", 
        cols = pd.col, label = TRUE, repel = TRUE, pt.size = 0.25, split.by = "Species")
dev.off()
# Stromal
inte.list3 <- list(hs.ref.list$Stromal, sheep.anno$Stromal$query)
inte.list3[[1]]@meta.data <- inte.list3[[1]]@meta.data %>% 
  mutate(Species = "Human",
         Group = "Human.Breast",
         Stage = "Mature",
         CellType.final = author_cell_type)
inte.list3[[1]]@meta.data <- inte.list3[[1]]@meta.data[, c("Group", "Stage", "Species", "Cluster", "CellType.final")]
inte.list3[[2]]@meta.data <- inte.list3[[2]]@meta.data %>% 
  mutate(Species = "Sheep",
         CellType.final = scPred.prediction)
inte.list3[[2]] <- subset(inte.list3[[2]], CellType.final != "unassigned")
inte.list3[[2]]@meta.data <- inte.list3[[2]]@meta.data[, c("Group", "Stage", "Species", "Cluster", "CellType.final")]
anchors <- FindIntegrationAnchors(object.list = inte.list3, dims = 1:20)
inte.list3$merged <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(inte.list3$merged) <- "integrated"
inte.list3$merged <- ScaleData(inte.list3$merged, verbose = FALSE)
inte.list3$merged <- RunPCA(inte.list3$merged, npcs = 30, verbose = FALSE)
ElbowPlot(inte.list3$merged)
inte.list3$merged <- RunUMAP(inte.list3$merged, reduction = "pca", dims = 1:10)
DimPlot(inte.list3$merged, reduction = "umap", group.by = c("Species", "CellType.final"), 
        cols = pd.col, label = TRUE, repel = TRUE)
inte.list3$merged <- RunHarmony(inte.list3$merged, group.by.vars = "Species")
inte.list3$merged <- RunUMAP(inte.list3$merged, reduction = "harmony", dims = 1:10)
pdf(file.path(res.out, "Scatter_to_show_stromal_annotation_results_by_final_integration.pdf"), height = 6, width = 15)
DimPlot(inte.list3$merged, reduction = "umap", group.by = c("Species", "CellType.final"), 
        cols = pd.col, label = TRUE, repel = TRUE, pt.size = 0.5)
dev.off()
pdf(file.path(res.out, "Scatter_to_show_stromal_annotation_results_by_final_integration_splited.pdf"), height = 6, width = 12.5)
DimPlot(inte.list3$merged, reduction = "umap", group.by = "CellType.final", 
        cols = pd.col, label = TRUE, repel = TRUE, pt.size = 0.5, split.by = "Species")
dev.off()


### >>> 7. Correlation analysis
# method 1: Epithelial
expr.hs = AverageExpression(subset(inte.list$merged, Species == "Human"), 
                            assays = "RNA", group.by = "CellType.final", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.sheep = AverageExpression(subset(inte.list$merged, Species == "Sheep"), 
                               assays = "RNA", group.by = "CellType.final", slot = "data") %>% as.data.frame()
colnames(expr.sheep) = paste("Sheep", gsub("RNA.", "", colnames(expr.sheep)), sep = "_")
deg.hs = VariableFeatures(subset(inte.list$merged, Species == "Human"))
deg.sheep = VariableFeatures(subset(inte.list$merged, Species == "Sheep"))
spe.score.seurat <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                    ExpressionTableSpecies2 = expr.sheep, DEgenesSpecies2 = deg.sheep,
                                    Species1 = "human", Species2 = "human",
                                    filename = file.path(res.out, "Epithelial_Specificity_score_correlation"))
# method 1: Immune
expr.hs = AverageExpression(subset(inte.list2$merged, Species == "Human"), 
                            assays = "RNA", group.by = "CellType.final", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.sheep = AverageExpression(subset(inte.list2$merged, Species == "Sheep"), 
                               assays = "RNA", group.by = "CellType.final", slot = "data") %>% as.data.frame()
colnames(expr.sheep) = paste("Sheep", gsub("RNA.", "", colnames(expr.sheep)), sep = "_")
deg.hs = VariableFeatures(subset(inte.list2$merged, Species == "Human"))
deg.sheep = VariableFeatures(subset(inte.list2$merged, Species == "Sheep"))
spe.score.seurat <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                    ExpressionTableSpecies2 = expr.sheep, DEgenesSpecies2 = deg.sheep,
                                    Species1 = "human", Species2 = "human",
                                    filename = file.path(res.out, "Immune_Specificity_score_correlation"))
# method 1: Stromal
expr.hs = AverageExpression(subset(inte.list3$merged, Species == "Human"), 
                            assays = "RNA", group.by = "CellType.final", slot = "data") %>% as.data.frame()
colnames(expr.hs) = paste("Hs", gsub("RNA.", "", colnames(expr.hs)), sep = "_")
expr.sheep = AverageExpression(subset(inte.list3$merged, Species == "Sheep"), 
                               assays = "RNA", group.by = "CellType.final", slot = "data") %>% as.data.frame()
colnames(expr.sheep) = paste("Sheep", gsub("RNA.", "", colnames(expr.sheep)), sep = "_")
deg.hs = VariableFeatures(subset(inte.list3$merged, Species == "Human"))
deg.sheep = VariableFeatures(subset(inte.list3$merged, Species == "Sheep"))
spe.score.seurat <- CorrComparePlot(ExpressionTableSpecies1 = expr.hs, DEgenesSpecies1 = deg.hs,
                                    ExpressionTableSpecies2 = expr.sheep, DEgenesSpecies2 = deg.sheep,
                                    Species1 = "human", Species2 = "human",
                                    filename = file.path(res.out, "Stromal_Specificity_score_correlation"))
# - method 2: Epithelial
ct.score <- list()
ct.score[["Epithelial"]] <- Pipe.MetaNeighbor(sc.obj = inte.list$merged,
                                              study.col = "Species",
                                              celltype.col = "CellType.final",
                                              var.genes = VariableFeatures(inte.list$merged),
                                              res.out = file.path(res.out, "MetaNeighbor_Epithelial"))
g1 <- str_split_fixed(rownames(ct.score[["Epithelial"]]$AUROC.scores), "\\|", 2)[,1]
g1.col <- c(Human = "#E41A1C", Sheep = "#377EB8")
g2 <- str_split_fixed(rownames(ct.score[["Epithelial"]]$AUROC.scores), "\\|", 2)[,2]
g2.col <- pd.col[1:length(unique(g2))]
names(g2.col) <- sort(unique(g2))
Plot.MetaNeighbor(score = ct.score[["Epithelial"]]$AUROC.scores, g1 = g1, g2 = g2,
                  g1.col = g1.col, g2.col = g2.col,
                  res.out = file.path(res.out, "MetaNeighbor_Epithelial"))
rm(tmp, g1, g2, g1.col, g2.col)
# - method 2: Immune
tmp <- inte.list2$merged[, sample(1:ncol(inte.list2$merged), 5000)]
ct.score[["Immune"]] <- Pipe.MetaNeighbor(sc.obj = tmp,
                                          study.col = "Species",
                                          celltype.col = "CellType.final",
                                          var.genes = VariableFeatures(inte.list2$merged),
                                          res.out = file.path(res.out, "MetaNeighbor_Immune"))
g1 <- str_split_fixed(rownames(ct.score[["Immune"]]$AUROC.scores), "\\|", 2)[,1]
g1.col <- c(Human = "#E41A1C", Sheep = "#377EB8")
g2 <- str_split_fixed(rownames(ct.score[["Immune"]]$AUROC.scores), "\\|", 2)[,2]
g2.col <- pd.col[1:length(unique(g2))]
names(g2.col) <- sort(unique(g2))
Plot.MetaNeighbor(score = ct.score[["Immune"]]$AUROC.scores, g1 = g1, g2 = g2,
                  g1.col = g1.col, g2.col = g2.col,
                  res.out = file.path(res.out, "MetaNeighbor_Immune"))
rm(tmp, g1, g2, g1.col, g2.col)
# - method 2: Stromal
tmp <- inte.list3$merged[, sample(1:ncol(inte.list3$merged), 5000)]
ct.score[["Stromal"]] <- Pipe.MetaNeighbor(sc.obj = tmp,
                                           study.col = "Species",
                                           celltype.col = "CellType.final",
                                           var.genes = VariableFeatures(inte.list3$merged),
                                           res.out = file.path(res.out, "MetaNeighbor_Stromal"))
g1 <- str_split_fixed(rownames(ct.score[["Stromal"]]$AUROC.scores), "\\|", 2)[,1]
g1.col <- c(Human = "#E41A1C", Sheep = "#377EB8")
g2 <- str_split_fixed(rownames(ct.score[["Stromal"]]$AUROC.scores), "\\|", 2)[,2]
g2.col <- pd.col[1:length(unique(g2))]
names(g2.col) <- sort(unique(g2))
Plot.MetaNeighbor(score = ct.score[["Stromal"]]$AUROC.scores, g1 = g1, g2 = g2,
                  g1.col = g1.col, g2.col = g2.col,
                  res.out = file.path(res.out, "MetaNeighbor_Stromal"))
rm(tmp, g1, g2, g1.col, g2.col)


### >>> 8. Cell Cycle analysis
exp.mat <- read.table(file = "/home/yhw/document/Seurat/nestorawa_forcellcycle_expressionMatrix.txt", 
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
inte.list$merged <- CellCycleScoring(inte.list$merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
inte.list2$merged <- CellCycleScoring(inte.list2$merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
inte.list3$merged <- CellCycleScoring(inte.list3$merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


### >>> 9. Plot gene expression
breast.markers <- read.csv("/home/yhw/document/public_data/GSE195665_human_breast/Cell_Markers.csv")
marker.gene <- split(breast.markers[, 5], f = breast.markers$CellType)
# Epithelial
library(SCP)
ht <- GroupHeatmap(
  srt = inte.list$merged,
  features = unique(c(marker.gene$`LumHR-SCGB`, marker.gene$`LumHR-major`, marker.gene$`LumHR-active`, 
                      marker.gene$`LumSec-HLA`, marker.gene$`LumSec-myo`, marker.gene$`LumSec-prol`, 
                      marker.gene$`LumSec-basal`, marker.gene$`LumSec-KIT`, marker.gene$`LumSec-major`, 
                      marker.gene$`LumSec-lac`, marker.gene$Basal)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_epithelial_cell_type_markers.pdf"), height = 15, width = 13)
print(ht$plot)
dev.off()
# Immune
inte.list2$merged@meta.data <- inte.list2$merged@meta.data %>% 
  mutate(Cluster2 = case_when(CellType.final %in% c("Macro-lipo", "mDC", "Mono-non-classical", 
                                                    "Mono-classical", "cDC2", "cDC1", "Macro-m2", 
                                                    "Macro-IFN", "Macro-m1", "Macro-m1-CCL", "Macro-m2-CXCL", 
                                                    "mye-prol", "pDC", "Neutrophil", "Mast") ~ "Myeloid cells",
                              CellType.final %in% c("CD8-Tem", "NKT", "GD", "CD4-naive", "NK",
                                                    "NK-ILCs", "CD4-activated", "CD4-Treg", "CD4-Th", 
                                                    "CD8-activated", "T_prol", "CD4-Th-like", 
                                                    "CD4-Tem", "CD8-Trm") ~ "T cells",
                              CellType.final %in% c("b_naive", "bmem_switched", "plasma_IgG",
                                                    "plasma_IgA", "bmem_unswitched") ~ "B cells"))

ht <- GroupHeatmap(
  srt = subset(inte.list2$merged, Cluster2 == "B cells"),
  features = unique(c(marker.gene$`B naive`, marker.gene$`B mem-switched`, 
                      marker.gene$`B plasma-IgG`, marker.gene$`B plasma-IgA`,
                      marker.gene$`B mem-unswitched`)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(10, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_immune_B_cells_cell_type_markers.pdf"), height = 7, width = 13)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = subset(inte.list2$merged, Cluster2 == "T cells"),
  features = unique(c(marker.gene$`CD8+ TEM`, marker.gene$NKT, 
                      marker.gene$`GD-T`, marker.gene$`CD4+ Naive`,
                      marker.gene$NK, marker.gene$`NK/ILCs`, marker.gene$`CD4+ Activated`,
                      marker.gene$`CD4+ Reg`, marker.gene$`CD4+ Th`, marker.gene$`CD8+ Activated`, 
                      marker.gene$`T-prol`, marker.gene$`CD4+ Th-like`, 
                      marker.gene$`CD4+ TEM`, marker.gene$`CD8+ TRM`)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(10, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_immune_T_cells_cell_type_markers.pdf"), height = 15, width = 15)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = subset(inte.list2$merged, Cluster2 == "Myeloid cells"),
  features = unique(c(marker.gene$`Macro-lipo`, marker.gene$mDC, 
                      marker.gene$`Mono-classical`, marker.gene$`Mono-nonclassical`,
                      marker.gene$cDC2, marker.gene$cDC1, marker.gene$`Macro-m2`,
                      marker.gene$`Macro-IFN`, marker.gene$`Macro-m1`, marker.gene$`Macro-m1-CLL`, marker.gene$`Macro-m2-CLL`,
                      marker.gene$`Mye-prol`, marker.gene$pDC, marker.gene$Neutrophil, marker.gene$Mast)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_immune_myeloid_cells_cell_type_markers.pdf"), height = 15, width = 15)
print(ht$plot)
dev.off()
# Stromal
ht <- GroupHeatmap(
  srt = inte.list3$merged,
  features = unique(c(marker.gene$Fibroblast, marker.gene$`Fibro-major`, marker.gene$`Fibro-matrix`, 
                      marker.gene$`Fibro-prematrix`, marker.gene$`Fibro-SFRP4`, 
                      marker.gene$Lymphatic, marker.gene$`Lym-immune`, marker.gene$`Lym-major`, 
                      marker.gene$`Lym-valve1`, marker.gene$`Lym-valve2`, 
                      marker.gene$Vascular, marker.gene$`Vas-arterial`, 
                      marker.gene$`Vas-capillary`, marker.gene$`Vas-venous`, 
                      marker.gene$`Vas-venous`, 
                      marker.gene$Perivascular, marker.gene$Pericyte, marker.gene$VSMC)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(6, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_stromal_cell_type_markers.pdf"), height = 15, width = 13)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = inte.list3$merged,
  features = unique(c(marker.gene$Fibroblast, marker.gene$`Fibro-major`, marker.gene$`Fibro-matrix`, 
                      marker.gene$`Fibro-prematrix`, marker.gene$`Fibro-SFRP4`)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_stromal_fibroblast_cell_type_markers.pdf"), height = 7, width = 13)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = inte.list3$merged,
  features = unique(c(marker.gene$Lymphatic, marker.gene$`Lym-immune`, marker.gene$`Lym-major`, 
                      marker.gene$`Lym-valve1`, marker.gene$`Lym-valve2`)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_stromal_lymphatic_cell_type_markers.pdf"), height = 7, width = 13)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = inte.list3$merged,
  features = unique(c(marker.gene$Vascular, marker.gene$`Vas-arterial`, 
                      marker.gene$`Vas-capillary`, marker.gene$`Vas-venous`, 
                      marker.gene$`Vas-venous`, 
                      marker.gene$Perivascular, marker.gene$Pericyte, marker.gene$VSMC)),
  group.by = c("Species", "CellType.final"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "S.Score", "G2M.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, assay = "RNA", dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "Dotplot_to_show_stromal_vascular_cell_type_markers.pdf"), height = 10, width = 13)
print(ht$plot)
dev.off()
# 
sr.out <- file.path(res.out, "Marker_genes_of_ORG_merged_ORG")
dir.create(sr.out, recursive = T)
DefaultAssay(sr.foxh1.org.merge2) <- "RNA"
for (type in names(marker.gene)) {
  marker.gene[[type]] <- intersect(marker.gene[[type]], rownames(sr.foxh1.org.merge2))
}
for (type in names(marker.gene)) {
  tmp.out <- file.path(sr.out, type)
  dir.create(tmp.out, recursive = T)
  for (g in seq(1, length(marker.gene[[type]]))) {
    pdf(file.path(tmp.out, paste0("Scatter_plot_to_show_gene_expression_of_", marker.gene[[type]][g], ".pdf")),
        height = 4.5, width = 5)
    print(FeatureDimPlot(
      srt = pancreas_sub, features = c("Sox9", "Neurog3", "Fev", "Rbp4"),
      reduction = "UMAP", theme_use = "theme_blank"
    ))
    dev.off()
  }; rm(g)
}


### >>> 9. 
# Epithelial
table(inte.list$merged$Species, inte.list$merged$CellType.final)
pd <- AverageExpression(inte.list$merged, group.by = c("Species", "CellType.final"), assays = "RNA", slot = "data")
write.csv(pd$RNA, file.path(res.out, "Epithelial_RNA_expression.csv"), row.names = T,quote = F)
pd <- AverageExpression(inte.list$merged, group.by = c("Species", "CellType.final"), assays = "integrated", slot = "data")
write.csv(pd$integrated, file.path(res.out, "Epithelial_Integrated_expression.csv"), row.names = T, quote = F)
# Immune
table(inte.list2$merged$Species, inte.list2$merged$CellType.final)
pd <- AverageExpression(inte.list2$merged, group.by = c("Species", "CellType.final"), assays = "RNA", slot = "data")
write.csv(pd$RNA, file.path(res.out, "Immune_RNA_expression.csv"), row.names = T,quote = F)
pd <- AverageExpression(inte.list2$merged, group.by = c("Species", "CellType.final"), assays = "integrated", slot = "data")
write.csv(pd$integrated, file.path(res.out, "Immune_Integrated_expression.csv"), row.names = T, quote = F)
# Stromal
table(inte.list3$merged$Species, inte.list3$merged$CellType.final)
pd <- AverageExpression(inte.list3$merged, group.by = c("Species", "CellType.final"), assays = "RNA", slot = "data")
write.csv(pd$RNA, file.path(res.out, "Stromal_RNA_expression.csv"), row.names = T,quote = F)
pd <- AverageExpression(inte.list3$merged, group.by = c("Species", "CellType.final"), assays = "integrated", slot = "data")
write.csv(pd$integrated, file.path(res.out, "Stromal_Integrated_expression.csv"), row.names = T, quote = F)



# ==============================
# 10th part: Metabolism analysis ----
# ==============================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Metabolism_Analysis")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Run metabolism analysis
table(sr.sheep.final$CellType)
table(sr.sheep.final$Cluster)
meta.sheep <- Pipe.Metabolism(sr.obj = sr.sheep.final, species = "human", group.by = "CellType", 
                              pt.size = 0.25, threads = 16, cluster.n = 5, prefix = "Sheep", 
                              outdir = file.path(res.out, "Sheep_all_cells"))


### >>> 3. Differential expression analysis
# cell type
meta.sheep <- RunDEtest(srt = meta.sheep, group_by = "CellType", fc.threshold = 1, only.pos = FALSE)
pdf(file.path(res.out, "Volcano_plot_to_show_differential_expression_analysis_in_cell_type_level.pdf"), height = 20, width = 30)
VolcanoPlot(srt = meta.sheep, group_by = "CellType", pt.size = 4)
dev.off()
# cluster
meta.sheep <- RunDEtest(srt = meta.sheep, group_by = "Cluster", fc.threshold = 1, only.pos = FALSE)
pdf(file.path(res.out, "Volcano_plot_to_show_differential_expression_analysis_in_cluster_level.pdf"), height = 18, width = 20)
VolcanoPlot(srt = meta.sheep, group_by = "Cluster", pt.size = 4)
dev.off()


### >>> 4. Visualize the metabolism
# scatter plot
dir.create(file.path(res.out, "Scatter_plots"))
for (i in rownames(meta.sheep)) {
  pdf(file.path(res.out, paste0("Scatter_plots/Pathway_", gsub("/", "", i), ".pdf")), height = 6, width = 6)
  print(FeatureDimPlot(srt = meta.sheep, features = i,
                       reduction = "UMAP", theme_use = "theme_blank", 
                       assay = "metabolism", slot = "data"))
  dev.off()
}
# violin plot: cell type
dir.create(file.path(res.out, "Violin_plots"))
for (i in rownames(meta.sheep)) {
  pdf(file.path(res.out, paste0("Violin_plots/Pathway_", gsub("/", "", i), "_in_cell_type_level.pdf")), height = 6, width = 20)
  print(FeatureStatPlot(srt = meta.sheep, stat.by = i, add_trend = T, add_box = T,
                        pairwise_method = "t.test", group.by = "CellType", 
                        split.by = "Stage", comparisons = TRUE))
  dev.off()
}
# violin plot: cluster
for (i in rownames(meta.sheep)) {
  pdf(file.path(res.out, paste0("Violin_plots/Pathway_", gsub("/", "", i), "_in_cluster_level.pdf")), height = 6, width = 10)
  print(FeatureStatPlot(srt = meta.sheep, stat.by = i, add_trend = T,  add_box = T,
                        pairwise_method = "t.test", group.by = "Cluster", 
                        split.by = "Stage", comparisons = TRUE))
  dev.off()
}
# heatmap: cell type
ht1 <- GroupHeatmap(meta.sheep, cluster_rows = T, cluster_columns = T, show_row_names = T,
                    features = rownames(meta.sheep),
                    group.by = c("CellType"), assay = "metabolism", slot = "data")
dev.off()
pdf(file.path(res.out, "Heatmap_plot_to_show_all_metabolism_in_cell_type_level.pdf"), height = 20, width = 20)
ht1
dev.off()
ht1 <- GroupHeatmap(meta.sheep, cluster_rows = T, cluster_columns = T, show_row_names = T,
                    features = rownames(meta.sheep), split.by = "Stage", 
                    group.by = c("CellType"), assay = "metabolism", slot = "data")
dev.off()
pdf(file.path(res.out, "Heatmap_plot_to_show_all_metabolism_in_cell_type_level_splited.pdf"), height = 20, width = 20)
ht1
dev.off()
# heatmap: cluster
ht1 <- GroupHeatmap(meta.sheep, cluster_rows = T, cluster_columns = T, show_row_names = T,
                    features = rownames(meta.sheep),
                    group.by = c("Cluster"), assay = "metabolism", slot = "data")
dev.off()
pdf(file.path(res.out, "Heatmap_plot_to_show_all_metabolism_in_cluster_level.pdf"), height = 20, width = 15)
ht1
dev.off()
ht1 <- GroupHeatmap(meta.sheep, cluster_rows = T, cluster_columns = T, show_row_names = T,
                    features = rownames(meta.sheep), split.by = "Stage", 
                    group.by = c("Cluster"), assay = "metabolism", slot = "data")
dev.off()
pdf(file.path(res.out, "Heatmap_plot_to_show_all_metabolism_in_cluster_level_splited.pdf"), height = 20, width = 15)
ht1
dev.off()
# dotplot: cluster
ht <- GroupHeatmap(srt = meta.sheep, features = rownames(meta.sheep),
                   group.by = "Cluster", heatmap_palette = "YlOrRd", 
                   cluster_rows = T, cluster_columns = T, exp_cutoff = 0.025,
                   add_dot = TRUE, add_reticle = TRUE, assay = "metabolism", slot = "data")
print(ht$plot)



# ==============================
# 11th part: Fibroblast analysis ----
# ==============================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Fibroblast_Analysis")
if (!dir.exists(res.out)) {
  dir.create(res.out, recursive = T)
}


### >>> 2. Clustering
sr.fibro <- subset(inte.list3$merged, CellType.final %in% c("Fibro-major", "Fibro-matrix", 
                                                            "Fibro-prematrix", "Fibro-SFRP4"))
DefaultAssay(sr.fibro) <- "RNA"
sr.list <- SplitObject(sr.fibro, split.by = "Stage")
for (i in 1:length(sr.list)) {
  sr.list[[i]] <- NormalizeData(sr.list[[i]], verbose = FALSE)
  sr.list[[i]] <- FindVariableFeatures(sr.list[[i]], selection.method = "vst", 
                                       nfeatures = 1000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = sr.list, dims = 1:15)
sr.fibro <- IntegrateData(anchorset = anchors, dims = 1:15)
DefaultAssay(sr.fibro) <- "integrated"
sr.fibro <- ScaleData(sr.fibro, verbose = FALSE)
sr.fibro <- RunPCA(sr.fibro, npcs = 15, verbose = FALSE)
ElbowPlot(sr.fibro)
sr.fibro <- RunUMAP(sr.fibro, reduction = "pca", dims = 1:6)
pdf(file.path(res.out, "Scatter_plot_after_correction_after_annotation.pdf"), height = 5, width = 10)
CellDimPlot(srt = sr.fibro, group.by = c("Stage", "CellType.final"), 
            label_repel = T, label = T, pt.size = 0.5, label_insitu = F,
            reduction = "UMAP", theme_use = "theme_blank", palette = "Paired")
dev.off()
# cell cycle
s.genes <- intersect(cc.genes$s.genes, rownames(sr.fibro))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(sr.fibro))
sr.fibro <- CellCycleScoring(sr.fibro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
pdf(file.path(res.out, "Stacked_barplot_to_show_cell_type_count_after_annotation.pdf"), 
    height = 8, width = 5)
CellStatPlot(sr.fibro, stat.by = "CellType.final", stat_type = "percent", 
             group.by = "Stage", plot_type = "trend", label = T)
dev.off()



### >>> Cluster-specific genes and GO
# COSG
Idents(sr.fibro) <- sr.fibro$CellType.final
library("COSG")
fibro.markers <- cosg(sr.fibro, groups = 'all', assay = 'integrated', 
                      slot = 'data', mu = 1, n_genes_user = 300)
write.csv(fibro.markers, file.path(res.out, "All_marker_genes_COSG.csv"), quote = F, row.names = F)
pd.gene <- c(fibro.markers$names$`Fibro-prematrix`[1:10], 
             fibro.markers$names$`Fibro-matrix`[1:10],
             fibro.markers$names$`Fibro-SFRP4`[1:10], 
             fibro.markers$names$`Fibro-major`[1:10])
ht <- GroupHeatmap(srt = sr.fibro, flip = F, 
                   limits = c(0, 2), exp_cutoff = 0, exp_method = "zscore", 
                   assay = "integrated", slot = "data",
                   cluster_columns = F, cluster_rows = F,
                   features = pd.gene,
                   group.by = c("CellType.final"), heatmap_palette = "Reds",
                   cell_split_palette = "Blues", cell_annotation = c("G2M.Score"), 
                   cell_annotation_palette = c("Paired"),
                   show_row_names = TRUE, row_names_side = "left",
                   add_dot = TRUE, add_reticle = TRUE, dot_size = unit(5, "mm"))
dev.off()
pdf(file.path(res.out, paste0("Dot_plot_show_cell_markers_COSG_after_annotation_", i, ".pdf")), 
    height = 8, width = 5)
print(ht$plot)
dev.off()
fibro.go <- list()
for (i in colnames(fibro.markers$names)) {
  fibro.go[[i]] <- Pipe.GO(species = "human", 
                           genelist = fibro.markers$names[, i],
                           basename = i, 
                           genetype = "SYMBOL",
                           res.out = file.path(res.out, "GO_COSG"))
}
# Seurat
Idents(sr.fibro) <- sr.fibro$CellType.final
fibro.markers2 <- FindAllMarkers(sr.fibro, assay = "integrated", only.pos = T, logfc.threshold = log2(1.25))
write.csv(fibro.markers2, file.path(res.out, "All_marker_genes_Seurat.csv"), quote = F, row.names = F)

pd.gene <- fibro.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ht <- GroupHeatmap(srt = sr.fibro, flip = F, 
                   limits = c(0, 2), exp_cutoff = 0, exp_method = "zscore", 
                   assay = "integrated", slot = "data",
                   cluster_columns = F, cluster_rows = F,
                   features = pd.gene$gene,
                   group.by = c("CellType.final"), heatmap_palette = "Reds",
                   cell_split_palette = "Blues", cell_annotation = c("G2M.Score"), 
                   cell_annotation_palette = c("Paired"),
                   show_row_names = TRUE, row_names_side = "left",
                   add_dot = TRUE, add_reticle = TRUE, dot_size = unit(5, "mm"))
dev.off()
pdf(file.path(res.out, paste0("Dot_plot_show_cell_markers_Seurat_after_annotation_", i, ".pdf")), 
    height = 8, width = 5)
print(ht$plot)
dev.off()
fibro.go2 <- list()
for (i in unique(fibro.markers2$cluster)) {
  fibro.go2[[i]] <- Pipe.GO(species = "human", 
                            genelist = subset(fibro.markers2, cluster == i)$gene,
                            basename = i, 
                            genetype = "SYMBOL",
                            res.out = file.path(res.out, "GO_Seurat"))
}


### >>> 8. Differential expression analysis
table(sr.fibro$CellType.final, sr.fibro$Stage)
deg.sc2 <- list()
for (i in unique(sr.fibro$CellType.final)) {
  tmp <- subset(sr.fibro, CellType.final == i)
  test <- SeuratToEdger(sr.obj = tmp, group.by = "Stage",
                        comparison = c("LL", "PL"), sample.n = 300)
  print(table(tmp$CellType.final, tmp$Stage))
  tmp <- gsub(" ", "_", i)
  deg.sc2[[tmp]] <- edgeR.scRNA(count = test$count, meta = test$meta, expr = test$expr,
                                sample.n = NULL, g1 = "1_LL", g2 = "2_PL",
                                lfc = log2(2), sig = 0.05,
                                res.out = file.path(res.out, paste0("DEGs/PL_vs_LL/", tmp)))
}
rm(tmp, tmp.sr, i, tmp.name)
# plotting
for (i in names(deg.sc2)) {
  # all
  write.csv(subset(deg.sc2[[i]]$all, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t1_2_PL_Vs_1_LL_all_TFs_results.csv")),
            row.names = F)
  pdf(file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/Volcano_plot_to_show_all_genes.pdf")),
      height = 5, width = 7)
  print(VisDEG.volcano(deg.data = deg.sc2[[i]]$all, geneset = NULL,
                       p.col = "PValue", lfc.col = "logFC", sig = 0.05, lfc = 1,
                       title = paste0(i, ": all genes"), up.col = "#B71C1C", down.col = "#01579B"))
  dev.off()
  pdf(file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/Volcano_plot_to_show_all_TFs.pdf")),
      height = 5, width = 7)
  print(VisDEG.volcano(deg.data = deg.sc2[[i]]$all, geneset = sheep.tfs$Symbol,
                       p.col = "PValue", lfc.col = "logFC", sig = 0.05, lfc = 1,
                       title = paste0(i, ": all TFs"), up.col = "#B71C1C", down.col = "#01579B"))
  dev.off()
  # all sig
  write.csv(subset(deg.sc2[[i]]$sig, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t2_2_PL_Vs_1_LL_all_sig_TFs_results.csv")),
            row.names = F)
  # up sig
  write.csv(subset(deg.sc2[[i]]$up, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t3_2_PL_Vs_1_LL_up_sig_TFs_lfc1_sig0.05.csv")),
            row.names = F)
  # down sig
  write.csv(subset(deg.sc2[[i]]$down, SYMBOL %in% sheep.tfs$Symbol),
            file.path(res.out, paste0("DEGs/PL_vs_LL/", i, "/edgeR_t4_2_PL_Vs_1_LL_down_sig_TFs_lfc1_sig0.05.csv")),
            row.names = F)
  
}
# plot multiple plots
FeatureDimPlot(sr.sema6a.hs, c("GAD1", "GAD2", "GSH2", "DARPP32", "BCL11B",
                               "NKX2-1", "CHAT", "OLIG2", "ISL1", "SLC18A3",
                               "NKX2-2", "NKX6-1", "DLX1", "DLX2", "DLX5"), 
               reduction = "UMAP", theme_use = "theme_blank", assay = "alra")
table(GetAssayData(sr.sema6a.hs, slot = "data", assay = "alra")["CHAT", ] > 0)
pd.gene <- c("FN1","FBN1","FBLN1","FBLN5","FBLN7","HAS2","HAS3","LUM","COL14A1","COL18A1",
             "COL5A3","COL5A1","COL5A2","COL12A1","COL27A1","COL4A5","COL9A3","COL4A1",
             "COL15A1","COL3A1","COL1A2","COL4A2","COL6A3","COL6A2","COL12A1","COL16A1",
             "COL9A1","COL4A4","COL19A1","COL21A1","COL23A1","COL6A1","COL11A2","COL5A2",
             "COL17A1","COL8A1","COL6A6","COL10A1","COL22A1","COL11A1")
pd.gene <- intersect(pd.gene, rownames(sr.sheep))
pd <- list()
for (i in names(deg.sc2)) {
  pd[[gsub(" ", "_", gsub("-", "_", gsub("\\/", "-", i)))]] <- deg.sc2[[i]]$all
}
pdf(file.path(res.out, "DEG_All_DEG_vs_Ctrl_volcano_custom.pdf"), height = 8, width = 8)
VisDEG.volcano.multi(deg.list = pd, geneset = pd.gene, 
                     lfc = log2(1.5), sig = 0.05,
                     top.n = 10, pt.size = 3, jitter.width = 0.35, label.gene = NULL,
                     high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#B2B2B2")
dev.off()
pdf(file.path(res.out, "DEG_All_DEG_vs_Ctrl_volcano_custom2.pdf"), height = 8, width = 8)
VisDEG.volcano.multi(deg.list = pd, geneset = NULL, 
                     lfc = log2(1.5), sig = 0.05,
                     top.n = 10, pt.size = 1, jitter.width = 0.35, label.gene = NULL,
                     high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#B2B2B2")
dev.off()


### >>> 5. GO and GSEA analysis
go.res2 <- list()
gsea.res2 <- list()
for (i in names(deg.sc2)) {
  gsea.res2[[i]] <- Pipe.GSEA(deg.obj = deg.sc2[[i]]$all,
                              deg.type = "edger", lfc = log2(1.25), sig = 0.05,
                              species = "human", basename = i,
                              genetype = "SYMBOL", gene.col = "SYMBOL",
                              outdir = file.path(res.out, paste0("GSEA/PL_vs_LL/", i)))
  go.res2[[paste0(i, ".up")]] <- Pipe.GO(species = "human",
                                         genelist = subset(deg.sc2[[i]]$all, PValue <= 0.05 & logFC >= log2(1.25))$SYMBOL,
                                         basename = paste0(i, "_up"),
                                         genetype = "SYMBOL",
                                         res.out = file.path(res.out, paste0("GO/PL_vs_LL/", i)))
  go.res2[[paste0(i, ".down")]] <- Pipe.GO(species = "human",
                                           genelist = subset(deg.sc2[[i]]$all, PValue <= 0.05 & logFC <= -log2(1.25))$SYMBOL,
                                           basename = paste0(i, "_down"),
                                           genetype = "SYMBOL",
                                           res.out = file.path(res.out, paste0("GO/PL_vs_LL/", i)))
}
library(aPEAR)
for (i in names(go.res2)) {
  for (j in names(go.res2$Fibroblast.up)) {
    pd <- na.omit(subset(go.res2[[i]][[j]], pvalue <= 0.05)[1:20, ])
    if (nrow(pd) >= 5) {
      pdf(file.path(res.out, paste0("GO/Network_plot_of_", i, "_", j, ".pdf")), height = 10, width = 10)
      print(enrichmentNetwork(pd, colorBy = "pvalue", colorType = "pval", pCutoff = -5, drawEllipses = TRUE))
      dev.off()
    }
  }
}



# ====================
# Last part: Save Data ----
# ====================
library(qs)
rm(tmp, test, p1, p2, p3, pd, pd.gene, sr.list.anno, sr.list, test, i, j)
qsave(sr.sheep, "R/CodeData/proj12.qs")
saveRDS(sr.sheep, "R/CodeData/proj12.rds")
save.image("R/CodeData/proj12.RData")
load("R/CodeData/proj12.RData")
