# -------------------------------
# Load Required Libraries
# -------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(SingleR)
  library(celldex)
})

set.seed(123)

# -------------------------------
# Load Cell Ranger Data
# -------------------------------
data_dir <- "D:/Breast_Cancer_Single_Cell/filtered_feature_bc_matrix"

raw_counts <- Read10X(data.dir = data_dir)

# Create Seurat object
sce <- CreateSeuratObject(
  counts = raw_counts,
  project = "IDC_scRNA",
  min.cells = 3,
  min.features = 200
)

# -------------------------------
# 3. Quality Control
# -------------------------------
sce[["percent.mt"]] <- PercentageFeatureSet(
  sce,
  pattern = "^MT-"
)

sce <- subset(
  sce,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

# -------------------------------
# SCTransform
# -------------------------------
sce <- SCTransform(
  sce,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

DefaultAssay(sce) <- "SCT"

# -------------------------------
# 5. PCA, Clustering, UMAP
# -------------------------------
sce <- RunPCA(sce, verbose = FALSE)

sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.4)

sce <- RunUMAP(sce, dims = 1:20)

DimPlot(sce, reduction = "umap", label = TRUE) + NoLegend()

# ============================================================
# AUTOMATIC CELL-TYPE ANNOTATION (SingleR)
# ============================================================

# Convert Seurat to SingleCellExperiment
sce_sce <- as.SingleCellExperiment(sce, assay = "SCT")

# Load reference atlas
ref <- HumanPrimaryCellAtlasData()

# Run SingleR
singleR_pred <- SingleR(
  test = sce_sce,
  ref = ref,
  labels = ref$label.main
)

# Add annotations to Seurat metadata
sce$SingleR_label <- singleR_pred$labels

# Visualize annotation
DimPlot(
  sce,
  reduction = "umap",
  group.by = "SingleR_label",
  label = TRUE,
  repel = TRUE
) + NoLegend()

# ============================================================
# EPITHELIAL CELL EXTRACTION
# ============================================================

# Identify epithelial cells using markers + annotation
epi_cells <- WhichCells(
  sce,
  expression = EPCAM > 1 | KRT8 > 1 | KRT18 > 1
)

epi <- subset(sce, cells = epi_cells)

DimPlot(epi, reduction = "umap", label = TRUE)

# ============================================================
# TUMOR vs NORMAL EPITHELIAL CLASSIFICATION
# ============================================================

# Heuristic rule:
# Normal-like: low nCount_RNA + low EMT markers
# Tumor-like: high nCount_RNA or VIM+

epi$epi_state <- ifelse(
  epi$nCount_RNA > median(epi$nCount_RNA) | FetchData(epi, "VIM") > 1,
  "Tumor_Epithelial",
  "Normal_Epithelial"
)

epi$epi_state <- factor(
  epi$epi_state,
  levels = c("Normal_Epithelial", "Tumor_Epithelial")
)

DimPlot(
  epi,
  reduction = "umap",
  group.by = "epi_state",
  cols = c("steelblue", "firebrick"),
  pt.size = 1
)

# ============================================================
# DIFFERENTIAL EXPRESSION: Tumor vs Normal Epithelium
# ============================================================

Idents(epi) <- "epi_state"

tumor_vs_normal <- FindMarkers(
  epi,
  ident.1 = "Tumor_Epithelial",
  ident.2 = "Normal_Epithelial",
  assay = "SCT",
  logfc.threshold = 0.25,
  min.pct = 0.25
)

write.csv(
  tumor_vs_normal,
  file = "Tumor_vs_Normal_Epithelial_DEGs.csv"
)

# Top DEGs visualization
top_genes <- rownames(
  tumor_vs_normal %>%
    arrange(desc(avg_log2FC)) %>%
    head(20)
)

DoHeatmap(
  epi,
  features = top_genes,
  group.by = "epi_state",
  assay = "SCT"
)

#============================
# Define Metabolic Gene Sets
#============================

metabolic_programs <- list(
  
  Glycolysis = c(
    "HK2","HK1","PFKP","ALDOA","GAPDH",
    "ENO1","PKM","LDHA"
  ),
  
  TCA_Cycle = c(
    "CS","ACO2","IDH3A","OGDH",
    "SUCLG1","SDHA","FH","MDH2"
  ),
  
  OxPhos = c(
    "NDUFA1","NDUFB8","SDHB",
    "UQCRC1","COX5A","ATP5F1A"
  ),
  
  Pentose_Phosphate = c(
    "G6PD","PGLS","PGD","TKT","TALDO1"
  ),
  
  Fatty_Acid_Synthesis = c(
    "ACLY","ACACA","FASN","SCD"
  ),
  
  Glutamine_Metabolism = c(
    "GLS","GLUD1","SLC1A5"
  )
)

# ===============================
# FDCA: Flux Score Computation
#================================

# Keep genes present in dataset
metabolic_programs <- lapply(
  metabolic_programs,
  function(g) intersect(g, rownames(epi))
)

# Remove empty pathways
metabolic_programs <- metabolic_programs[
  sapply(metabolic_programs, length) > 3
]

#===============================
# Add Flux Scores (FDCA Core)
#===============================

epi <- AddModuleScore(
  epi,
  features = metabolic_programs,
  assay = "SCT",
  name = "FDCA_Flux"
)

#==================================
# Visualize Flux Landscapes (FBCA)
#==================================

FeaturePlot(
  epi,
  features = grep("^FDCA_Flux", colnames(epi@meta.data), value = TRUE),
  reduction = "umap",
  ncol = 3
)

#====================================
# Flux Comparison: Tumor vs Normal
#====================================

flux_df <- epi@meta.data %>%
  select(epi_state, starts_with("FDCA_Flux")) %>%
  pivot_longer(
    cols = starts_with("FDCA_Flux"),
    names_to = "Pathway",
    values_to = "Flux"
  )

ggplot(flux_df, aes(x = epi_state, y = Flux, fill = epi_state)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~ Pathway, scales = "free") +
  theme_minimal() +
  labs(
    title = "FDCA Metabolic Flux Comparison",
    subtitle = "Tumor vs Normal Epithelial Cells"
  )

#=========================================
# FBCA: Cluster-Level Flux Dynamics
#=========================================

epi$cluster <- Idents(epi)

cluster_flux <- epi@meta.data %>%
  group_by(cluster) %>%
  summarise(across(starts_with("FDCA_Flux"), mean))

cluster_flux


#=================================
# Save FDCA / FBCA Results
#=================================
write.csv(
  epi@meta.data %>% select(epi_state, starts_with("FDCA_Flux")),
  file = "FDCA_CellLevel_FluxScores.csv"
)

write.csv(
  cluster_flux,
  file = "FBCA_ClusterLevel_FluxScores.csv"
)
