# LOAD LIBRARIES ####
# Restart Rstudio or R

# install.packages('ggplot2')
# install.packages('cowplot')
# install.packages('Matrix')
# install.packages('ggridges')
# install.packages('ggrepel')
# install.packages('dplyr')
# install.packages('Seurat')
# install.packages('monocle3')
# install.packages('plotly')
# install.packages('clustree')
# install.packages('patchwork')
# install.packages('future')

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("DoubletFinder")
# BiocManager::install("glmGamPoi")

# Run the following code once you have Seurat installed
suppressWarnings(
  {
    library(Rcpp)
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(Seurat)
    library(monocle3)
    library(harmony)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
    library(glmGamPoi)
  }
)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle3")
packageVersion("harmony")

# Set global environment parameter
#options(future.globals.maxSize = 8000 * 1024^2)

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
{
  HP2022801.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\1_220624_Fahd_GEX1_F1_HP-20228-01\filtered_feature_bc_matrix)")
  SAMN15877725.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\2_220624_Fahd_GEX2_F2_SAMN15877725\filtered_feature_bc_matrix)")
  HP2024001.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\3_220624_Fahd_GEX3_F3_HP-20240-01\filtered_feature_bc_matrix)")
  HP2031401.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\4_220624_Fahd_GEX4_F4_HP-20314-01\filtered_feature_bc_matrix)")
  HP2105501.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\5_210406 GEX_F5_HP-21055-01\filtered_feature_bc_matrix)")
  HP2106201.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\6_210406 GEX_F6_HP-21062-01\filtered_feature_bc_matrix)")
  HP2107001.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\7_210406 GEX_F7a_HP-21070-01\filtered_feature_bc_matrix)")
  HP2107901.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\9_210714 GEX_F9a_HP-21079-01\filtered_feature_bc_matrix)")
  HP2108601.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\10_211108 GEX_F10a_HP-21086-01\filtered_feature_bc_matrix)")
  HP2108901.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\11_211108 GEX_F11a_HP-21089-01\filtered_feature_bc_matrix)")
  HP2110001.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\12_211108 GEX_F12a_HP-21100-01\filtered_feature_bc_matrix)")
  HP2121601.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\13_211108 GEX_F13a_HP-21216-01\filtered_feature_bc_matrix)")
  HP2123201.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\14_220124 GEX_F14_HP21232-01\filtered_feature_bc_matrix)")
  HP2132801.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\15_220124 GEX_F15_HP-21328-01\filtered_feature_bc_matrix)")
  HP2202101.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\16_220324 GEX_F16_HP-22021-01\filtered_feature_bc_matrix)")


# STEP 2: Create Seurat objects ####

  HP2022801 <- CreateSeuratObject(counts = HP2022801.data, min.features = 500)
  SAMN15877725 <- CreateSeuratObject(counts = SAMN15877725.data, min.features = 500)
  HP2024001 <- CreateSeuratObject(counts = HP2024001.data, min.features = 500)
  HP2031401 <- CreateSeuratObject(counts = HP2031401.data, min.features = 500)
  HP2105501 <- CreateSeuratObject(counts = HP2105501.data, min.features = 500)
  HP2106201 <- CreateSeuratObject(counts = HP2106201.data, min.features = 500)
  HP2107001 <- CreateSeuratObject(counts = HP2107001.data, min.features = 500)
  HP2107901 <- CreateSeuratObject(counts = HP2107901.data, min.features = 500)
  HP2108601 <- CreateSeuratObject(counts = HP2108601.data, min.features = 500)
  HP2108901 <- CreateSeuratObject(counts = HP2108901.data, min.features = 500)
  HP2110001 <- CreateSeuratObject(counts = HP2110001.data, min.features = 500)
  HP2121601 <- CreateSeuratObject(counts = HP2121601.data, min.features = 500)
  HP2123201 <- CreateSeuratObject(counts = HP2123201.data, min.features = 500)
  HP2132801 <- CreateSeuratObject(counts = HP2132801.data, min.features = 500)
  HP2202101 <- CreateSeuratObject(counts = HP2202101.data, min.features = 500)
  }

# Sample specific Metadata addition
{
  HP2022801$sample <- "HP2022801"
  SAMN15877725$sample <- "SAMN15877725"
  HP2024001$sample <- "HP2024001"
  HP2031401$sample <- "HP2031401"
  HP2105501$sample <- "HP2105501"
  HP2106201$sample <- "HP2106201"
  HP2107001$sample <- "HP2107001"
  HP2107901$sample <- "HP2107901"
  HP2108601$sample <- "HP2108601"
  HP2108901$sample <- "HP2108901"
  HP2110001$sample <- "HP2110001"
  HP2121601$sample <- "HP2121601"
  HP2123201$sample <- "HP2123201"
  HP2132801$sample <- "HP2132801"
  HP2202101$sample <- "HP2202101"


# Sex specific Metadata addition
  HP2022801$sex <- "male"
  SAMN15877725$sex <- "male"
  HP2024001$sex <- "female"
  HP2031401$sex <- "male"
  HP2105501$sex <- "female"
  HP2106201$sex <- "female"
  HP2107001$sex <- "male"
  HP2107901$sex <- "male"
  HP2108601$sex <- "female"
  HP2108901$sex <- "female"
  HP2110001$sex <- "male"
  HP2121601$sex <- "female"
  HP2123201$sex <- "male"
  HP2132801$sex <- "female"
  HP2202101$sex <- "female"

# Ancestry specific Metadata addition
  HP2022801$ancestry <- "white"
  SAMN15877725$ancestry <- "white"
  HP2024001$ancestry <- "white"
  HP2031401$ancestry <- "black"
  HP2105501$ancestry <- "white"
  HP2106201$ancestry <- "black"
  HP2107001$ancestry <- "white"
  HP2107901$ancestry <- "white"
  HP2108601$ancestry <- "white"
  HP2108901$ancestry <- "white"
  HP2110001$ancestry <- "black"
  HP2121601$ancestry <- "black"
  HP2123201$ancestry <- "black"
  HP2132801$ancestry <- "black"
  HP2202101$ancestry <- "black"

# Ancestry and sex specific Metadata addition
  HP2022801$ancestry_sex <- "white_male"
  SAMN15877725$ancestry_sex <- "white_male"
  HP2024001$ancestry_sex <- "white_female"
  HP2031401$ancestry_sex <- "black_male"
  HP2105501$ancestry_sex <- "white_female"
  HP2106201$ancestry_sex <- "black_female"
  HP2107001$ancestry_sex <- "white_male"
  HP2107901$ancestry_sex <- "white_male"
  HP2108601$ancestry_sex <- "white_female"
  HP2108901$ancestry_sex <- "white_female"
  HP2110001$ancestry_sex <- "black_male"
  HP2121601$ancestry_sex <- "black_female"
  HP2123201$ancestry_sex <- "black_male"
  HP2132801$ancestry_sex <- "black_female"
  HP2202101$ancestry_sex <- "black_female"
  }

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
{
  HP2022801[["percent.mt"]] <- PercentageFeatureSet(object = HP2022801, pattern = "^MT-")
  SAMN15877725[["percent.mt"]] <- PercentageFeatureSet(object = SAMN15877725, pattern = "^MT-")
  HP2024001[["percent.mt"]] <- PercentageFeatureSet(object = HP2024001, pattern = "^MT-")
  HP2031401[["percent.mt"]] <- PercentageFeatureSet(object = HP2031401, pattern = "^MT-")
  HP2105501[["percent.mt"]] <- PercentageFeatureSet(object = HP2105501, pattern = "^MT-")
  HP2106201[["percent.mt"]] <- PercentageFeatureSet(object = HP2106201, pattern = "^MT-")
  HP2107001[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001, pattern = "^MT-")
  HP2107901[["percent.mt"]] <- PercentageFeatureSet(object = HP2107901, pattern = "^MT-")
  HP2108601[["percent.mt"]] <- PercentageFeatureSet(object = HP2108601, pattern = "^MT-")
  HP2108901[["percent.mt"]] <- PercentageFeatureSet(object = HP2108901, pattern = "^MT-")
  HP2110001[["percent.mt"]] <- PercentageFeatureSet(object = HP2110001, pattern = "^MT-")
  HP2121601[["percent.mt"]] <- PercentageFeatureSet(object = HP2121601, pattern = "^MT-")
  HP2123201[["percent.mt"]] <- PercentageFeatureSet(object = HP2123201, pattern = "^MT-")
  HP2132801[["percent.mt"]] <- PercentageFeatureSet(object = HP2132801, pattern = "^MT-")
  HP2202101[["percent.mt"]] <- PercentageFeatureSet(object = HP2202101, pattern = "^MT-")

# QC information before thresholding
  summary(head(HP2022801@meta.data))
  summary(head(SAMN15877725@meta.data))
  summary(head(HP2024001@meta.data))
  summary(head(HP2031401@meta.data))
  summary(head(HP2105501@meta.data))
  summary(head(HP2106201@meta.data))
  summary(head(HP2107001@meta.data))
  summary(head(HP2107901@meta.data))
  summary(head(HP2108601@meta.data))
  summary(head(HP2108901@meta.data))
  summary(head(HP2110001@meta.data))
  summary(head(HP2121601@meta.data))
  summary(head(HP2123201@meta.data))
  summary(head(HP2132801@meta.data))
  summary(head(HP2202101@meta.data))

# Visualize QC metrics as a violin plot
  VlnPlot(object = HP2022801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = SAMN15877725, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2024001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2031401, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2105501, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2106201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2110001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2121601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2123201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2132801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2202101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
  HP2022801 <- subset(x = HP2022801, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  SAMN15877725 <- subset(x = SAMN15877725, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2024001 <- subset(x = HP2024001, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2031401 <- subset(x = HP2031401, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2105501 <- subset(x = HP2105501, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2106201 <- subset(x = HP2106201, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2107001 <- subset(x = HP2107001, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2107901 <- subset(x = HP2107901, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2108601 <- subset(x = HP2108601, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2108901 <- subset(x = HP2108901, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2110001 <- subset(x = HP2110001, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2121601 <- subset(x = HP2121601, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2123201 <- subset(x = HP2123201, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2132801 <- subset(x = HP2132801, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2202101 <- subset(x = HP2202101, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# QC information after thresholding
  summary(head(HP2022801@meta.data))
  summary(head(SAMN15877725@meta.data))
  summary(head(HP2024001@meta.data))
  summary(head(HP2031401@meta.data))
  summary(head(HP2105501@meta.data))
  summary(head(HP2106201@meta.data))
  summary(head(HP2107001@meta.data))
  summary(head(HP2107901@meta.data))
  summary(head(HP2108601@meta.data))
  summary(head(HP2108901@meta.data))
  summary(head(HP2110001@meta.data))
  summary(head(HP2121601@meta.data))
  summary(head(HP2123201@meta.data))
  summary(head(HP2132801@meta.data))
  summary(head(HP2202101@meta.data))

# Visualize QC metrics post thresholding as a violin plot
  VlnPlot(object = HP2022801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = SAMN15877725, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2024001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2031401, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2105501, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2106201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2110001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2121601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2123201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2132801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2202101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Step 4: Add cell IDs ####
# Add cell IDs
  HP2022801 <- RenameCells(HP2022801, add.cell.id = "HP2022801")
  SAMN15877725 <- RenameCells(SAMN15877725, add.cell.id = "SAMN15877725")
  HP2024001 <- RenameCells(HP2024001, add.cell.id = "HP2024001")
  HP2031401 <- RenameCells(HP2031401, add.cell.id = "HP2031401")
  HP2105501 <- RenameCells(HP2105501, add.cell.id = "HP2105501")
  HP2106201 <- RenameCells(HP2106201, add.cell.id = "HP2106201")
  HP2107001 <- RenameCells(HP2107001, add.cell.id = "HP2107001")
  HP2107901 <- RenameCells(HP2107901, add.cell.id = "HP2107901")
  HP2108601 <- RenameCells(HP2108601, add.cell.id = "HP2108601")
  HP2108901 <- RenameCells(HP2108901, add.cell.id = "HP2108901")
  HP2110001 <- RenameCells(HP2110001, add.cell.id = "HP2110001")
  HP2121601 <- RenameCells(HP2121601, add.cell.id = "HP2121601")
  HP2123201 <- RenameCells(HP2123201, add.cell.id = "HP2123201")
  HP2132801 <- RenameCells(HP2132801, add.cell.id = "HP2132801")
  HP2202101 <- RenameCells(HP2202101, add.cell.id = "HP2202101")
  }

# Step 5: creating a list of all datasets
{
  pancreas.list <- list("HP2022801" = HP2022801, "SAMN15877725" = SAMN15877725, "HP2107001" = HP2107001, "HP2107901" = HP2107901,
                        "HP2024001" = HP2024001, "HP2105501" = HP2105501, "HP2108601" = HP2108601, "HP2108901" = HP2108901, 
                        "HP2031401" = HP2031401, "HP2110001" = HP2110001, "HP2123201" = HP2123201,
                        "HP2106201" = HP2106201, "HP2121601" = HP2121601, "HP2132801" = HP2132801, "HP2202101" = HP2202101
                         )
}

# Normalization for visualization
pancreas.combined <- lapply(X = pancreas.list, FUN = function(x) {
  x <- NormalizeData(x)
 })

# normalize and identify variable features for each dataset independently
pancreas.list <- lapply(X = pancreas.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = features)
pancreas.list <- lapply(X = pancreas.list, FUN = RunPCA, features = features)

# Perform integration note k.anchors = 20 increase in integration strength (5 by default)
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT",
                                           anchor.features = features, dims = 1:30, reduction = "rpca", 
                                           k.anchor = 5)
pancreas.combined <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", dims = 1:30,
                                   verbose = TRUE)

# Post integration UMAP
pancreas.combined <- RunPCA(pancreas.combined, npcs = 50, verbose = TRUE)
pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", dims = 1:50)

# Step 9a: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
DefaultAssay(pancreas.combined) <- "integrated"
pancreas.combined <- FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:50)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.1)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.2)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.3)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.4)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.5)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.combined, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.3
pancreas.combined <- FindClusters(pancreas.combined, resolution = 0.3)

# Alternatively build a cluster tree
DefaultAssay(object = pancreas.combined) <- "integrated"
pancreas.combined = BuildClusterTree(pancreas.combined, slot = "scale.data")
PlotClusterTree(pancreas.combined)

# Visualization
DimPlot(pancreas.combined, reduction = "umap", group.by = "ancestry_sex")
DimPlot(pancreas.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(pancreas.combined, reduction = "umap", group.by = "integrated_snn_res.0.3", label = TRUE, repel = TRUE)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.combined, prefix = "integrated_snn_res.")

# Discovery based Plotting
DefaultAssay(pancreas.combined) <- "RNA"
FeaturePlot(object = pancreas.combined,
            features = c("SOX10"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 100,
            slot = 'counts',
            order = TRUE)

#Rename Idents
pancreas.combined <- RenameIdents(pancreas.combined, 
                                    "0" = "Alpha-GCGhi", 
                                    "1" = "Beta-INShi",
                                    "2" = "Activated-Stellate", 
                                    "3" = "Ductal",
                                    "4" = "Alpha-GCGhi", 
                                    "5" = "Acinar",
                                    "6" = "Transdifferentiating-Endo", 
                                    "7" = "Beta-INSlow",
                                    "8" = "Alpha-GCGlow", 
                                    "9" = "Endothelial",
                                    "10" = "Delta", 
                                    "11" = "Quiescent-Stellate1",
                                    "12" = "Quiescent-Stellate2",
                                    "13" = "Acinar",
                                    "14" = "Alpha-GCGlow",
                                    "15" = "Beta-ERStress",
                                    "16" = "Macrophages",
                                    "17" = "Gamma",
                                    "18" = "Mast",
                                    "19" = "Ductal",
                                    "20" = "Delta",
                                    "21" = "Transdifferentiating-Exo",
                                    "22" = "T-cells",
                                    "23" = "Endothelial",
                                    "24" = "Schwann"
                                    )

# Manual allocation of GHRL epsilon cells
DefaultAssay(object = pancreas.combined) <- "RNA"
Idents(pancreas.combined, WhichCells(object = pancreas.combined, expression = GHRL > 100, slot = 'counts')) <- 'Epsilon'
DimPlot(pancreas.combined, reduction = "umap", label = TRUE)

# Saving this information in the metadata slot
table(Idents(pancreas.combined))
pancreas.combined$celltype <- Idents(pancreas.combined)
summary(pancreas.combined@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("Beta-INShi", "Beta-INSlow", "Beta-ERStress", "Transdifferentiating-Endo", "Alpha-GCGhi", "Alpha-GCGlow", "Delta", "Gamma", "Epsilon",
               "Ductal", "Transdifferentiating-Exo", "Acinar", 
               "Quiescent-Stellate1", "Quiescent-Stellate2", "Activated-Stellate",
               "Macrophages", "T-cells", "Mast",
               "Schwann", "Endothelial")
head(pancreas.combined@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.combined@meta.data$celltype <- factor(x = pancreas.combined@meta.data$celltype, levels = my_levels)
Idents(pancreas.combined) <- "celltype"

# Observing cells
DimPlot(pancreas.combined, 
        split.by = "ancestry_sex", group.by = "celltype", 
        label = FALSE, ncol = 2,  
        cols = c("red4", "red3", "grey40", "orange", "lightgoldenrod3", "yellow4", "indianred", "orangered", "black",
                 "royalblue2", "steelblue1", "darkcyan",
                 "springgreen4", "green3", "darkturquoise",
                 "purple4", "purple", "deeppink",
                 "violetred", "violet"
))

DimPlot(pancreas.combined, 
        group.by = "celltype", 
        label = FALSE, ncol = 1,  
        cols = c("red4", "red3", "grey40", "orange", "lightgoldenrod3", "yellow4", "indianred", "orangered", "black",
                 "royalblue2", "steelblue1", "darkcyan",
                 "springgreen4", "green3", "darkturquoise",
                 "purple4", "purple", "deeppink",
                 "violetred", "violet"
        ))

markers <- FindAllMarkers(pancreas.combined, assay = "RNA",
                          logfc.threshold = 1,
                          test.use = "wilcox",
                          slot = "data",
                          min.pct = 0.5)
# Saving and load
saveRDS(pancreas.combined, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.combined.rds)")
#pancreas.combined <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.combined.rds)")

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.combined) <- "celltype"

# Select only beta cells
pancreas.combined$celltype.sample <- paste(Idents(pancreas.combined),pancreas.combined$ancestry_sex, sep = "_")
table(pancreas.combined@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta-INShi_white_male", "Beta-INShi_white_female", "Beta-INShi_black_male", "Beta-INShi_black_female",
                "Beta-INSlow_white_male", "Beta-INSlow_white_female", "Beta-INSlow_black_male", "Beta-INSlow_black_female",
                "Beta-ERStress_white_male", "Beta-ERStress_white_female", "Beta-ERStress_black_male", "Beta-ERStress_black_female",
                "Transdifferentiating-Endo_white_male", "Transdifferentiating-Endo_white_female", "Transdifferentiating-Endo_black_male", "Transdifferentiating-Endo_black_female",
                "Alpha-GCGhi_white_male", "Alpha-GCGhi_white_female", "Alpha-GCGhi_black_male", "Alpha-GCGhi_black_female",
                "Alpha-GCGlow_white_male", "Alpha-GCGlow_white_female", "Alpha-GCGlow_black_male", "Alpha-GCGlow_black_female",
                "Delta_white_male", "Delta_white_female", "Delta_black_male", "Delta_black_female",
                "Gamma_white_male", "Gamma_white_female", "Gamma_black_male", "Gamma_black_female",
                "Epsilon_white_male", "Epsilon_white_female", "Epsilon_black_male", "Epsilon_black_female",
                "Ductal_white_male", "Ductal_white_female", "Ductal_black_male", "Ductal_black_female",
                "Transdifferentiating-Exo_white_male", "Transdifferentiating-Exo_white_female", "Transdifferentiating-Exo_black_male", "Transdifferentiating-Exo_black_female",
                "Acinar_white_male", "Acinar_white_female", "Acinar_black_male", "Acinar_black_female",
                "Quiescent-Stellate1_white_male", "Quiescent-Stellate1_white_female", "Quiescent-Stellate1_black_male", "Quiescent-Stellate1_black_female",
                "Quiescent-Stellate2_white_male", "Quiescent-Stellate2_white_female", "Quiescent-Stellate2_black_male", "Quiescent-Stellate2_black_female",
                "Activated-Stellate_white_male", "Activated-Stellate_white_female", "Activated-Stellate_black_male", "Activated-Stellate_black_female",
                "Macrophages_white_male", "Macrophages_white_female", "Macrophages_black_male", "Macrophages_black_female",
                "T-cells_white_male", "T-cells_white_female", "T-cells_black_male", "T-cells_black_female",
                "Mast_white_male", "Mast_white_female", "Mast_black_male", "Mast_black_female",
                "Schwann_white_male", "Schwann_white_female", "Schwann_black_male", "Schwann_black_female",
                "Endothelial_white_male", "Endothelial_white_female", "Endothelial_black_male", "Endothelial_black_female")

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.combined@meta.data$celltype.sample <- factor(x = pancreas.combined@meta.data$celltype.sample, levels = my_levels2)
table(pancreas.combined@meta.data$celltype.sample)

# Re select organized idents
Idents(pancreas.combined) <- "celltype.sample"
DefaultAssay(object = pancreas.combined) <- "SCT"

# Selected genes
markers.to.plot <- c("INS", "IAPP", "NKX6-1", "MAFA", "MAFB", "GCG", "DPP4", "GC", "LEPR", "SST", "FRZB", "PPY", "CALB1", "THSD7A", "GHRL", "PHGR1",
                     "CFTR", "KRT19", "MMP7", "CELA2A", "CELA2B", "CELA3A", "RGS5", "CSRP2", "FABP4", "COL3A1", "FMOD", "PDGFRB", "MKI67", "HIST1H4C", "STMN1", 
                     "CD86", "CSF1R", "SDS", "NKG7", "IL2RB", "CCL5", "RGS13", "TPSB2", "TPSAB1", "SOX10", "CDH19", "NGFR", "CD34", "ENG", "VWF", "UCN3")

# Dotplot
DotPlot(pancreas.combined,  
        dot.scale = 8,
        col.min = 0, #minimum level
        col.max = 1,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid = c("white"), high =c("red3")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# End #
