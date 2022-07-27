# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed

install.packages('clustree')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

suppressWarnings(
  {
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(Seurat)
    library(monocle3)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
  }
)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle3")

# Set global environment parameter
options(future.globals.maxSize = 8000 * 1024^2)

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
}

# STEP 2: Create Seurat objects ####
{
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
HP2022801$ancestry <- "white_male"
SAMN15877725$ancestry <- "white_male"
HP2024001$ancestry <- "white_female"
HP2031401$ancestry <- "black_male"
HP2105501$ancestry <- "white_female"
HP2106201$ancestry <- "black_female"
HP2107001$ancestry <- "white_male"
HP2107901$ancestry <- "white_male"
HP2108601$ancestry <- "white_female"
HP2108901$ancestry <- "white_female"
HP2110001$ancestry <- "black_male"
HP2121601$ancestry <- "black_female"
HP2123201$ancestry <- "black_male"
HP2132801$ancestry <- "black_female"
HP2202101$ancestry <- "black_female"
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
}

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

# Step 5: Merge Datasets
# Based on comment to Issue #4753 https://github.com/satijalab/seurat/issues/4753
# We use RPCA to yield conserved mapping and set Tx as control refrence samples
# Merge panc_sex datasets
pancreas.list <- list("HP2022801" = HP2022801, "SAMN15877725" = SAMN15877725, "HP2024001" = HP2024001, "HP2031401" = HP2031401,
                      "HP2105501" = HP2105501, "HP2106201" = HP2106201, "HP2107001" = HP2107001, "HP2107901" = HP2107901,
                      "HP2108601" = HP2108601, "HP2108901" = HP2108901, "HP2110001" = HP2110001, "HP2121601" = HP2121601,
                      "HP2123201" = HP2123201, "HP2132801" = HP2132801, "HP2202101" = HP2202101)

# Step 6: Data normalization
#Normalise data
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# Step 7: Feature selection
# Select features for downstream integration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list)
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- ScaleData(x, features = pancreas.features, verbose = FALSE)
  x <- RunPCA(x, features = pancreas.features, verbose = FALSE)
})

# Step 8: Anchor identification and data integration
# Identify anchors and integrate dataset
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = pancreas.features, reduction = "rpca", dims = 1:50)

# this command creates an 'integrated' data assay
pancreas.combined <- IntegrateData(anchorset = pancreas.anchors, dims = 1:50)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(pancreas.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.combined <- ScaleData(pancreas.combined, verbose = FALSE)
pancreas.combined <- RunPCA(pancreas.combined, npcs = 50, verbose = FALSE)
pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", dims = 1:50)
pancreas.combined <- FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:50)
pancreas.combined <- FindClusters(pancreas.combined, resolution = 0.5)

# Visualization
DimPlot(pancreas.combined, reduction = "umap", group.by = "seurat_clusters")

# Step 9: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Scaling this is weird, but as done in https://satijalab.org/seurat/articles/integration_large_datasets.html
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)

# Dimensionality assessment using PCA analysis
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))

# Examine data dimensionality
ElbowPlot(pancreas.integrated)

# Step 9a: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.1)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.2)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.3)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.4)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.5)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.integrated, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.5
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)

# Alternatively build a cluster tree
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated=BuildClusterTree(pancreas.integrated, slot = "scale.data")
PlotClusterTree(pancreas.integrated)

# Step 10: non-linear dimensionality assessment ####
# Run PCA and UMAP calculations
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)

# Change default assay to integrated, to view dimensionality
Idents(pancreas.integrated) <- "treatment"
DimPlot(pancreas.integrated, reduction = "umap", 
        cols = c('black', 'red'), 
        label = FALSE,
        order = FALSE)
Idents(pancreas.integrated) <- "sex"
Idents(pancreas.integrated) <- "sample"
Idents(pancreas.integrated) <- "seurat_clusters"
#Idents(pancreas.integrated) <- "integrated_snn_res.0.3"
DimPlot(pancreas.integrated, reduction = "umap", label = FALSE)

#Visualize gene expression
DefaultAssay(object = pancreas.integrated) <- "RNA"
DefaultAssay(object = pancreas.integrated)
FeaturePlot(object = pancreas.integrated,
            features = c("INS", "GCG", "SST", "PPY", "GHRL",
                         "KRT19", "CPA1",
                         "COL1A1", "VWF", "SOX10",
                         "TPSAB1", "SDS", "TRAC"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            #max.cutoff = 1,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.combined,
            features = c("SST", "INS", "GCG", "PPY", "GHRL"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            max.cutoff = 10000,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.combined,
            features = c("KRT19", "CPA1", "SPP1", "VWF", "SDS", "COL1A1"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            max.cutoff = 10000,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.combined,
            features = c("DDIT3"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            max.cutoff = 10000,
            slot = 'counts',
            order = TRUE)

VlnPlot(
  object = pancreas.integrated,
  features = c("AR"),
  assay = 'RNA',
  slot = 'counts',
  cols = c('red',
           'red4',
           'orange',
           'lightgoldenrod3',
           'sienna',
           'indianred',
           'orangered1',
           'black',
           'darkturquoise',
           'paleturquoise',
           'lightgreen',
           'springgreen4',
           'darkolivegreen',
           'purple4',
           'purple',
           'deeppink',
           'violetred',
           'violet'),
  #y.max = 3,
  pt.size = 1
)

#Rename Idents
pancreas.integrated <- RenameIdents(pancreas.integrated, 
                                    "0" = "Beta INS-hi", 
                                    "1" = "Alpha GCG-hi",
                                    "2" = "Ductal", 
                                    "3" = "Transdifferentiating Beta",
                                    "4" = "Beta INS-low", 
                                    "5" = "Acinar",
                                    "6" = "Activated Stellate", 
                                    "7" = "Ductal",
                                    "8" = "Endothelial", 
                                    "9" = "Delta",
                                    "10" = "Alpha GCG-hi", 
                                    "11" = "Quiescent Stellate",
                                    "12" = "Alpha GCG-low",
                                    "13" = "Ductal",
                                    "14" = "Quiescent Stellate",
                                    "15" = "Gamma",
                                    "16" = "Macrophage",
                                    "17" = "Proliferating Stellate",
                                    "18" = "Schwann",
                                    "19" = "Mast",
                                    "20" = "T-Lymphocyte"
)

#plot <- DimPlot(pancreas.integrated, reduction = "umap")
DefaultAssay(object = pancreas.integrated) <- "RNA"
Idents(pancreas.integrated, WhichCells(object = pancreas.integrated, expression = GHRL > 1, slot = 'counts')) <- 'Epsilon'
#pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Epsilon")

DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)

# Saving this information in the metadata slot
table(Idents(pancreas.integrated))
pancreas.integrated$celltype <- Idents(pancreas.integrated)
head(pancreas.integrated@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("Beta INS-hi", "Beta INS-low", "Transdifferentiating Beta", "Alpha GCG-hi", "Alpha GCG-low", "Delta", "Gamma", "Epsilon",
               "Ductal", "Acinar", 
               "Quiescent Stellate", "Activated Stellate", "Proliferating Stellate",
               "Macrophage", "T-Lymphocyte", "Mast",
               "Schwann", "Endothelial")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype <- factor(x = pancreas.integrated@meta.data$celltype, levels = my_levels)
Idents(pancreas.integrated) <- "celltype"

# Observing cells
DimPlot(pancreas.integrated, split.by = "sample", group.by = "celltype", label = FALSE, ncol = 2,  cols = c("red",
                                                                                                            "red4",
                                                                                                            "orange",
                                                                                                            "lightgoldenrod3",
                                                                                                            "sienna",
                                                                                                            "indianred",
                                                                                                            "orangered1",
                                                                                                            "black",
                                                                                                            "darkturquoise",
                                                                                                            "paleturquoise",
                                                                                                            "lightgreen",
                                                                                                            "springgreen4",
                                                                                                            "darkolivegreen",
                                                                                                            "purple4",
                                                                                                            "purple",
                                                                                                            "deeppink",
                                                                                                            "violetred",
                                                                                                            "violet"
))

Idents(pancreas.integrated) <- "treatment"
DHT <- subset(pancreas.integrated, idents = "DHT[10nM]")
DimPlot(DHT, group.by = "treatment", cols = "red")
EtOH <- subset(pancreas.integrated, idents = "EtOH")
DimPlot(EtOH, group.by = "treatment", cols = "blue")

DimPlot(pancreas.integrated, group.by = "treatment")
UMAPPlot(pancreas.integrated, reduction = "umap",
         pt.size = .75,
         cols = c("red",
                  "red4",
                  "orange",
                  "lightgoldenrod3",
                  "sienna",
                  "indianred",
                  "orangered1",
                  "black",
                  "darkturquoise",
                  "paleturquoise",
                  "lightgreen",
                  "springgreen4",
                  "darkolivegreen",
                  "purple4",
                  "purple",
                  "deeppink",
                  "violetred",
                  "violet"
         ),
         label = FALSE)
#beta.hi <- subset(pancreas.integrated, idents = "Beta INS-hi")

#table(beta.hi$treatment)

DotPlot(pancreas.integrated,
        group.by = "treatment",
        #split.by = "treatment",
        features = c("AR"), 
        cols = c("yellow", "red"), 
        col.min = -10, 
        col.max = 10)

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "celltype"

# Select only beta cells
beta.hi <- subset(pancreas.integrated, idents = "Beta INS-hi")

pancreas.integrated$celltype.sample <- paste(Idents(pancreas.integrated),pancreas.integrated$treatment, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta INS-hi_EtOH", "Beta INS-hi_DHT[10nM]", "Beta INS-low_EtOH", "Beta INS-low_DHT[10nM]", "Transdifferentiating Beta_EtOH", "Transdifferentiating Beta_DHT[10nM]", "Alpha GCG-hi_EtOH", "Alpha GCG-hi_DHT[10nM]", "Alpha GCG-low_EtOH", "Alpha GCG-low_DHT[10nM]", "Delta_EtOH", "Delta_DHT[10nM]", "Gamma_EtOH", "Gamma_DHT[10nM]", "Epsilon_EtOH", "Epsilon_DHT[10nM]",
                "Ductal_EtOH", "Ductal_DHT[10nM]", "Acinar_EtOH", "Acinar_DHT[10nM]", 
                "Quiescent Stellate_EtOH", "Quiescent Stellate_DHT[10nM]", "Activated Stellate_EtOH", "Activated Stellate_DHT[10nM]", "Proliferating Stellate_EtOH", "Proliferating Stellate_DHT[10nM]",
                "Macrophage_EtOH", "Macrophage_DHT[10nM]", "T-Lymphocyte_EtOH", "T-Lymphocyte_DHT[10nM]", "Mast_EtOH", "Mast_DHT[10nM]", "Schwann_EtOH", "Schwann_DHT[10nM]", "Endothelial_EtOH", "Endothelial_DHT[10nM]")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.sample <- factor(x = pancreas.integrated@meta.data$celltype.sample, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.sample)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.sample"
DefaultAssay(object = pancreas.integrated) <- "RNA"

# Selected genes
markers.to.plot <- c("MT-ATP8", "MT-ATP6",
                     "MT-CO1", "MT-CO2", "MT-CO3",
                     "MT-CYB",
                     "MT-ND6", "MT-ND5", "MT-ND4", "MT-ND4L", "MT-ND3", "MT-ND2", "MT-ND1")

# Dotplot
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        col.min = -1, #minimum level
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



# Select only beta cells
Idents(pancreas.integrated) <- "celltype.sample"
pathway <- subset(pancreas.integrated, idents = c("Beta INS-hi_EtOH", "Beta INS-hi_DHT[10nM]", "Beta INS-low_EtOH", "Beta INS-low_DHT[10nM]", 
                                                  "Alpha GCG-hi_EtOH", "Alpha GCG-hi_DHT[10nM]", "Alpha GCG-low_EtOH", "Alpha GCG-low_DHT[10nM]"))

# Selected genes
Idents(pancreas.integrated) <- "celltype.sample"
markers.to.plot <- c("MT-ATP8", "MT-ATP6",
                     "MT-CO1", "MT-CO2", "MT-CO3",
                     "MT-CYB",
                     "MT-ND6", "MT-ND5", "MT-ND4", "MT-ND4L", "MT-ND3", "MT-ND2", "MT-ND1")

markers.to.plot <- c("GCK", "ALDOA", "BPGM",
                     "ENO1", "ENO2", "GAPDH",
                     "GPI",
                     "PFKL", "PFKM", "PGAM1", "PGK1", "PKM", "TPI1")
markers.to.plot <- c("AR", "SLC25A4")
# Dotplot
DotPlot(pathway,
        dot.scale = 8,
        col.min = -1, #minimum level
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

VlnPlot(pathway, features = markers.to.plot, slot = "counts", #"data" 
        split.by = 'treatment', ncol = 1)

# Visualize co-expression of two features simultaneously
FeaturePlot(pancreas.integrated, features = c("GLP1R", "AR"), blend = TRUE, order = TRUE, blend.threshold = 0.01, pt.size = 2)
FeatureScatter(pancreas.integrated, feature1 = "GLP1R", feature2 = "AR", pt.size = 3,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))

# Calculate %of GLP1R+ Beta-hi cells expressing AR
Idents(pancreas.integrated) <- "celltype"
beta.hi <- subset(pancreas.integrated, idents = "Beta INS-hi")

Idents(beta.hi) <- "sample"
HP2107001_ctrl_beta.hi <- subset(beta.hi, idents = "HP2107001_ctrl")
Idents(HP2107001_ctrl_beta.hi) <- "celltype"
HP2107001_ctrl_beta.hi.GLP1R <- subset(HP2107001_ctrl_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_ctrl_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_ctrl_beta.hi@meta.data)*100

HP2107001_DHT_beta.hi <- subset(beta.hi, idents = "HP2107001_DHT")
Idents(HP2107001_DHT_beta.hi) <- "celltype"
HP2107001_DHT_beta.hi.GLP1R <- subset(HP2107001_DHT_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_DHT_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_DHT_beta.hi@meta.data)*100

HP2107701_ctrl_beta.hi <- subset(beta.hi, idents = "HP2107701_ctrl")
Idents(HP2107701_ctrl_beta.hi) <- "celltype"
HP2107701_ctrl_beta.hi.GLP1R <- subset(HP2107701_ctrl_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_ctrl_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_ctrl_beta.hi@meta.data)*100

HP2107701_DHT_beta.hi <- subset(beta.hi, idents = "HP2107701_DHT")
Idents(HP2107701_DHT_beta.hi) <- "celltype"
HP2107701_DHT_beta.hi.GLP1R <- subset(HP2107701_DHT_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_DHT_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_DHT_beta.hi@meta.data)*100

HP2107901_ctrl_beta.hi <- subset(beta.hi, idents = "HP2107901_ctrl")
Idents(HP2107901_ctrl_beta.hi) <- "celltype"
HP2107901_ctrl_beta.hi.GLP1R <- subset(HP2107901_ctrl_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_ctrl_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_ctrl_beta.hi@meta.data)*100

HP2107901_DHT_beta.hi <- subset(beta.hi, idents = "HP2107901_DHT")
Idents(HP2107901_DHT_beta.hi) <- "celltype"
HP2107901_DHT_beta.hi.GLP1R <- subset(HP2107901_DHT_beta.hi, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_DHT_beta.hi.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_DHT_beta.hi@meta.data)*100

# Calculate %of GLP1R+ Beta-low cells expressing AR
Idents(pancreas.integrated) <- "celltype"
beta.low <- subset(pancreas.integrated, idents = "Beta INS-low")

Idents(beta.low) <- "sample"
HP2107001_ctrl_beta.low <- subset(beta.low, idents = "HP2107001_ctrl")
Idents(HP2107001_ctrl_beta.low) <- "celltype"
HP2107001_ctrl_beta.low.GLP1R <- subset(HP2107001_ctrl_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_ctrl_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_ctrl_beta.low@meta.data)*100

HP2107001_DHT_beta.low <- subset(beta.low, idents = "HP2107001_DHT")
Idents(HP2107001_DHT_beta.low) <- "celltype"
HP2107001_DHT_beta.low.GLP1R <- subset(HP2107001_DHT_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_DHT_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_DHT_beta.low@meta.data)*100

HP2107701_ctrl_beta.low <- subset(beta.low, idents = "HP2107701_ctrl")
Idents(HP2107701_ctrl_beta.low) <- "celltype"
HP2107701_ctrl_beta.low.GLP1R <- subset(HP2107701_ctrl_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_ctrl_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_ctrl_beta.low@meta.data)*100

HP2107701_DHT_beta.low <- subset(beta.low, idents = "HP2107701_DHT")
Idents(HP2107701_DHT_beta.low) <- "celltype"
HP2107701_DHT_beta.low.GLP1R <- subset(HP2107701_DHT_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_DHT_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_DHT_beta.low@meta.data)*100

HP2107901_ctrl_beta.low <- subset(beta.low, idents = "HP2107901_ctrl")
Idents(HP2107901_ctrl_beta.low) <- "celltype"
HP2107901_ctrl_beta.low.GLP1R <- subset(HP2107901_ctrl_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_ctrl_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_ctrl_beta.low@meta.data)*100

HP2107901_DHT_beta.low <- subset(beta.low, idents = "HP2107901_DHT")
Idents(HP2107901_DHT_beta.low) <- "celltype"
HP2107901_DHT_beta.low.GLP1R <- subset(HP2107901_DHT_beta.low, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_DHT_beta.low.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_DHT_beta.low@meta.data)*100

# Calculate %of GLP1R+ Trans-beta cells expressing AR
Idents(pancreas.integrated) <- "celltype"
trans.beta <- subset(pancreas.integrated, idents = "Transdifferentiating Beta")

Idents(trans.beta) <- "sample"
HP2107001_ctrl_trans.beta <- subset(trans.beta, idents = "HP2107001_ctrl")
Idents(HP2107001_ctrl_trans.beta) <- "celltype"
HP2107001_ctrl_trans.beta.GLP1R <- subset(HP2107001_ctrl_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_ctrl_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_ctrl_trans.beta@meta.data)*100

HP2107001_DHT_trans.beta <- subset(trans.beta, idents = "HP2107001_DHT")
Idents(HP2107001_DHT_trans.beta) <- "celltype"
HP2107001_DHT_trans.beta.GLP1R <- subset(HP2107001_DHT_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107001_DHT_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107001_DHT_trans.beta@meta.data)*100

HP2107701_ctrl_trans.beta <- subset(trans.beta, idents = "HP2107701_ctrl")
Idents(HP2107701_ctrl_trans.beta) <- "celltype"
HP2107701_ctrl_trans.beta.GLP1R <- subset(HP2107701_ctrl_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_ctrl_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_ctrl_trans.beta@meta.data)*100

HP2107701_DHT_trans.beta <- subset(trans.beta, idents = "HP2107701_DHT")
Idents(HP2107701_DHT_trans.beta) <- "celltype"
HP2107701_DHT_trans.beta.GLP1R <- subset(HP2107701_DHT_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107701_DHT_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107701_DHT_trans.beta@meta.data)*100

HP2107901_ctrl_trans.beta <- subset(trans.beta, idents = "HP2107901_ctrl")
Idents(HP2107901_ctrl_trans.beta) <- "celltype"
HP2107901_ctrl_trans.beta.GLP1R <- subset(HP2107901_ctrl_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_ctrl_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_ctrl_trans.beta@meta.data)*100

HP2107901_DHT_trans.beta <- subset(trans.beta, idents = "HP2107901_DHT")
Idents(HP2107901_DHT_trans.beta) <- "celltype"
HP2107901_DHT_trans.beta.GLP1R <- subset(HP2107901_DHT_trans.beta, subset = GLP1R > 0)
length(WhichCells(object = HP2107901_DHT_trans.beta.GLP1R, expression = AR > 0 & GLP1R > 0))/nrow(HP2107901_DHT_trans.beta@meta.data)*100

FeatureScatter(beta.hi, feature1 = "GLP1R", feature2 = "AR", pt.size = 2,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))

FeatureScatter(beta.low, feature1 = "GLP1R", feature2 = "AR", pt.size = 2,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))

tranbeta <- subset(pancreas.integrated, idents = "Transdifferentiating Beta")
FeatureScatter(tranbeta, feature1 = "GLP1R", feature2 = "AR", pt.size = 2,
               cols = c("red",
                        "red4",
                        "orange",
                        "lightgoldenrod3",
                        "sienna",
                        "indianred",
                        "orangered1",
                        "black",
                        "darkturquoise",
                        "paleturquoise",
                        "lightgreen",
                        "springgreen4",
                        "darkolivegreen",
                        "purple4",
                        "purple",
                        "deeppink",
                        "violetred",
                        "violet"
               ))
# Save file this will change but for showing them on 07132021 its fine
#saveRDS(pancreas.integrated, file = r"(C:/Users/mqadir/Box/Lab 2301/RNAseq DHT data/wkdir/pancreas.integrated.rds)")
#pancreas.integrated <- readRDS(r"(C:/Users/mqadir/Box/Lab 2301/RNAseq DHT data/wkdir/pancreas.integrated.rds)")

# Identify conserved cell markers
DefaultAssay(pancreas.integrated) <- "RNA"
markers.beta.hi <- FindConservedMarkers(pancreas.integrated, ident.1 = "Beta INS-hi", grouping.var = "treatment", verbose = TRUE)
head(markers.beta.hi)

markers.beta.low <- FindConservedMarkers(pancreas.integrated, ident.1 = "Beta INS-low", grouping.var = "treatment", verbose = TRUE)
head(markers.beta.low)

markers.alpha.hi <- FindConservedMarkers(pancreas.integrated, ident.1 = "Alpha GCG-hi", grouping.var = "treatment", verbose = TRUE)
head(markers.alpha.hi)

markers.alpha.low <- FindConservedMarkers(pancreas.integrated, ident.1 = "Alpha GCG-low", grouping.var = "treatment", verbose = TRUE)
head(markers.alpha.low)

markers.transbeta <- FindConservedMarkers(pancreas.integrated, ident.1 = "Transdifferentiating Beta", grouping.var = "treatment", verbose = TRUE)
head(markers.transbeta)

markers.delta <- FindConservedMarkers(pancreas.integrated, ident.1 = "Delta", grouping.var = "treatment", verbose = TRUE)
head(markers.delta)

markers.gamma <- FindConservedMarkers(pancreas.integrated, ident.1 = "Gamma", grouping.var = "treatment", verbose = TRUE)
head(markers.gamma)

markers.epsilon <- FindConservedMarkers(pancreas.integrated, ident.1 = "Epsilon", grouping.var = "treatment", verbose = TRUE)
head(markers.epsilon)

markers.ductal <- FindConservedMarkers(pancreas.integrated, ident.1 = "Ductal", grouping.var = "treatment", verbose = TRUE)
head(markers.ductal)

markers.acinar <- FindConservedMarkers(pancreas.integrated, ident.1 = "Acinar", grouping.var = "treatment", verbose = TRUE)
head(markers.acinar)

markers.quiescentstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Quiescent Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.quiescentstellate)

markers.activatedstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Activated Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.activatedstellate)

markers.prolifstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Proliferating Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.prolifstellate)

markers.macrophage <- FindConservedMarkers(pancreas.integrated, ident.1 = "Macrophage", grouping.var = "treatment", verbose = TRUE)
head(markers.macrophage)

markers.tlympho <- FindConservedMarkers(pancreas.integrated, ident.1 = "T-Lymphocyte", grouping.var = "treatment", verbose = TRUE)
head(markers.tlympho)

markers.mast <- FindConservedMarkers(pancreas.integrated, ident.1 = "Mast", grouping.var = "treatment", verbose = TRUE)
head(markers.mast)

markers.schwann <- FindConservedMarkers(pancreas.integrated, ident.1 = "Schwann", grouping.var = "treatment", verbose = TRUE)
head(markers.schwann)

markers.endothelial <- FindConservedMarkers(pancreas.integrated, ident.1 = "Endothelial", grouping.var = "treatment", verbose = TRUE)
head(markers.endothelial)

# Now over-write the SCT assay with new-analyzed data from this subsetted data
pancreas.integrated <- SCTransform(pancreas.integrated, assay = "RNA", new.assay.name = "SCT", verbose = TRUE, return.only.var.genes = TRUE)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
# Here we define a DE gene as a gene which has:
# Fold Change of >1.1x
# Atleast 10% of cells express that gene
Idents(object = pancreas.integrated) <- "celltype"
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated.markers <- FindAllMarkers(object = pancreas.integrated, 
                                              features = VariableFeatures(pancreas.integrated, assay = 'integrated'), 
                                              only.pos = TRUE, 
                                              min.pct = 0.1, 
                                              logfc.threshold = 0.137504, 
                                              assay = 'RNA',
                                              slot = c('data'))

pancreas.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(pancreas.integrated.markers, r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\pancreas.integrated.markers.RNA.csv)")

# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to SCT, save information in the "SCT" assay
# You can toggle between integrated, SCT and RNA to see different expression profiles/different normalizations
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Create heatmap using doheatmap
top10.nomes <- pancreas.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = pancreas.integrated, 
          features = top10.nomes$gene, 
          disp.min = -1, 
          disp.max = 1,
          label = FALSE) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", 
                                                                            "#fbfcbd", 
                                                                            "#ff0000"))(256))
# Identify conserved cell markers
Idents(pancreas.integrated) <- factor(Idents(pancreas.integrated), levels = c("Beta INS-hi", "Beta INS-low", "Transdifferentiating Beta", "Alpha GCG-hi", "Alpha GCG-low", "Delta", "Gamma", "Epsilon",
                                                                              "Ductal", "Acinar", 
                                                                              "Quiescent Stellate", "Activated Stellate", "Proliferating Stellate", 
                                                                              "Macrophage", "T-Lymphocyte", "Mast", "Schwann", "Endothelial"))
markers.to.plot <- c("INS", "IAPP", "NKX6-1", "MAFA", "MAFB", "GCG", "DPP4", "GC", "LEPR", "SST", "FRZB", "PPY", "CALB1", "THSD7A", "GHRL", "PHGR1",
                     "CFTR", "KRT19", "MMP7", "CELA2A", "CELA2B", "CELA3A", "RGS5", "CSRP2", "FABP4", "COL3A1", "FMOD", "PDGFRB", "MKI67", "HIST1H4C", "STMN1", 
                     "CD86", "CSF1R", "SDS", "NKG7", "IL2RB", "CCL5", "RGS13", "TPSB2", "TPSAB1", "SOX10", "CDH19", "NGFR", "CD34", "ENG", "VWF", "UCN3")
markers.to.plot <- c("INS", "IAPP", "PDX1", "MAFA", "MAFB", "GCG", "DPP4", "GC", "LEPR", "SST", "PPY", "THSD7A", "GHRL", "FRZB",
                     "CFTR", "MMP7", "CELA2A", "CELA3A", "RGS5", "FABP4", "COL3A1", "FMOD", "MKI67", "STMN1", 
                     "CSF1R", "SDS", "NKG7", "CCL5", "TPSB2", "TPSAB1", "SOX10", "NGFR", "ENG", "VWF")

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "celltype"
pancreas.integrated$celltype.sample <- paste(Idents(pancreas.integrated),pancreas.integrated$treatment, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta INS-hi_EtOH", "Beta INS-hi_DHT[10nM]", "Beta INS-low_EtOH", "Beta INS-low_DHT[10nM]", "Transdifferentiating Beta_EtOH", "Transdifferentiating Beta_DHT[10nM]", "Alpha GCG-hi_EtOH", "Alpha GCG-hi_DHT[10nM]", "Alpha GCG-low_EtOH", "Alpha GCG-low_DHT[10nM]", "Delta_EtOH", "Delta_DHT[10nM]", "Gamma_EtOH", "Gamma_DHT[10nM]", "Epsilon_EtOH", "Epsilon_DHT[10nM]",
                "Ductal_EtOH", "Ductal_DHT[10nM]", "Acinar_EtOH", "Acinar_DHT[10nM]", 
                "Quiescent Stellate_EtOH", "Quiescent Stellate_DHT[10nM]", "Activated Stellate_EtOH", "Activated Stellate_DHT[10nM]", "Proliferating Stellate_EtOH", "Proliferating Stellate_DHT[10nM]",
                "Macrophage_EtOH", "Macrophage_DHT[10nM]", "T-Lymphocyte_EtOH", "T-Lymphocyte_DHT[10nM]", "Mast_EtOH", "Mast_DHT[10nM]", "Schwann_EtOH", "Schwann_DHT[10nM]", "Endothelial_EtOH", "Endothelial_DHT[10nM]")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.sample <- factor(x = pancreas.integrated@meta.data$celltype.sample, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.sample)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.sample"
Idents(pancreas.integrated) <- "celltype"
DefaultAssay(object = pancreas.integrated) <- "RNA"
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
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
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Older dotplot configuration, shows only percentage not expression
Idents(pancreas.integrated) <- "celltype"
DotPlot(pancreas.integrated, features = rev(markers.to.plot), 
        cols = c("blue", "red"), 
        dot.scale = 8, 
        split.by = "treatment") + 
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme_light() + 
  #coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =10, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =8, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10))

# Diff gene testing across conditions
# pancreas.integrated$treatment.dht <- paste(Idents(pancreas.integrated), pancreas.integrated$treatment, sep = "_")
# pancreas.integrated$celltype.split <- Idents(pancreas.integrated)
# choosing only those genes which are differentially expressed
# Optimise idents and assay
Idents(pancreas.integrated) <- "celltype.sample"
DefaultAssay(object = pancreas.integrated) <- "RNA"

# 1.Beta-cells (INS Hi)
beta.INSHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                       ident.1 = "Beta INS-hi_DHT[10nM]", ident.2 = "Beta INS-hi_EtOH", 
                                       test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                       min.pct = 0.1,
                                       logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                       pseudocount.use = 1,
                                       verbose = TRUE)
head(beta.INSHi.DHT.response, n = 15)
write.csv(beta.INSHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\beta.INShi.DHT.response.csv)")

# 2.Beta-cells (INS low)
beta.INSLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Beta INS-low_DHT[10nM]", ident.2 = "Beta INS-low_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0.1,
                                        logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(beta.INSLow.DHT.response, n = 15)
write.csv(beta.INSLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\beta.INSLow.DHT.response.csv)")

# 3.Alpha-cells (GCG hi)
alpha.GCGHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Alpha GCG-hi_DHT[10nM]", ident.2 = "Alpha GCG-hi_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0.1,
                                        logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(alpha.GCGHi.DHT.response, n = 15)
write.csv(alpha.GCGHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\alpha.GCGHi.DHT.response.csv)")

# 4.Alpha-cells (GCG low)
alpha.GCGLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                         ident.1 = "Alpha GCG-low_DHT[10nM]", ident.2 = "Alpha GCG-low_EtOH", 
                                         test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                         min.pct = 0.1,
                                         logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                         pseudocount.use = 1,
                                         verbose = TRUE)
head(alpha.GCGLow.DHT.response, n = 15)
write.csv(alpha.GCGLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\alpha.GCGLow.DHT.response.csv)")

# 5.Trandifferentiating Endocrine-Cells
tranbeta.DHT.response <- FindMarkers(pancreas.integrated, 
                                     ident.1 = "Transdifferentiating Beta_DHT[10nM]", ident.2 = "Transdifferentiating Beta_EtOH", 
                                     test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                     min.pct = 0.1,
                                     logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                     pseudocount.use = 1,
                                     verbose = TRUE)
head(tranbeta.DHT.response, n = 15)
write.csv(tranbeta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\tranbeta.DHT.response.csv)")

# 6.Delta-Cells
delta.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Delta_DHT[10nM]", ident.2 = "Delta_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0.1,
                                  logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(delta.DHT.response, n = 15)
write.csv(delta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\delta.DHT.response.csv)")

# 7.Gamma-Cells
gamma.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Gamma_DHT[10nM]", ident.2 = "Gamma_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0.1,
                                  logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(gamma.DHT.response, n = 15)
write.csv(gamma.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\gamma.DHT.response.csv)")

# 7.Epsilon-Cells
epsilon.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "Epsilon_DHT[10nM]", ident.2 = "Epsilon_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0.1,
                                    logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(epsilon.DHT.response, n = 15)
write.csv(epsilon.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\epsilon.DHT.response.csv)")

# 8.Ductal-Cells
ductal.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Ductal_DHT[10nM]", ident.2 = "Ductal_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0.1,
                                   logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(ductal.DHT.response, n = 15)
write.csv(ductal.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\ductal.DHT.response.csv)")

# 9.Acinar-Cells
acinar.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Acinar_DHT[10nM]", ident.2 = "Acinar_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0.1,
                                   logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(acinar.DHT.response, n = 15)
write.csv(acinar.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\acinar.DHT.response.csv)")

# 10.Quiescent Stellate-Cells
qstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Quiescent Stellate_DHT[10nM]", ident.2 = "Quiescent Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(qstellate.DHT.response, n = 15)
write.csv(qstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\qstellate.DHT.response.csv)")

# 11.Activated Stellate-Cells
astellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Activated Stellate_DHT[10nM]", ident.2 = "Activated Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(astellate.DHT.response, n = 15)
write.csv(astellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\astellate.DHT.response.csv)")

# 12.Proliferating Stellate-Cells
pstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Proliferating Stellate_DHT[10nM]", ident.2 = "Proliferating Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0.1,
                                      logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(pstellate.DHT.response, n = 15)
write.csv(pstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\pstellate.DHT.response.csv)")

# 13.Macrophage-Cells
macrophage.DHT.response <- FindMarkers(pancreas.integrated, 
                                       ident.1 = "Macrophage_DHT[10nM]", ident.2 = "Macrophage_EtOH", 
                                       test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                       min.pct = 0.1,
                                       logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                       pseudocount.use = 1,
                                       verbose = TRUE)
head(macrophage.DHT.response, n = 15)
write.csv(macrophage.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\macrophage.DHT.response.csv)")

# 14.T Lymphocyte-Cells
tlympho.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "T-Lymphocyte_DHT[10nM]", ident.2 = "T-Lymphocyte_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0.1,
                                    logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(tlympho.DHT.response, n = 15)
write.csv(tlympho.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\tlympho.DHT.response.csv)")

# 15.Mast-Cells
mast.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Mast_DHT[10nM]", ident.2 = "Mast_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0.1,
                                 logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = TRUE)
head(mast.DHT.response, n = 15)
write.csv(mast.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\mast.DHT.response.csv)")

# 16.Schwann-Cells
schwann.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "Schwann_DHT[10nM]", ident.2 = "Schwann_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0.1,
                                    logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(schwann.DHT.response, n = 15)
write.csv(schwann.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\schwann.DHT.response.csv)")

# 17.Endothelial-Cells
endothelial.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Endothelial_DHT[10nM]", ident.2 = "Endothelial_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0.1,
                                        logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(endothelial.DHT.response, n = 15)
write.csv(endothelial.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\endothelial.DHT.response.csv)")

#
#
#
#
#
# Running DE for all genes irrlevant of FC and PCT filtering
# Optimise idents
Idents(pancreas.integrated) <- "celltype.sample"
DefaultAssay(object = pancreas.integrated) <- "RNA"

# 1.Beta-cells (INS Hi)
beta.INSHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                       ident.1 = "Beta INS-hi_DHT[10nM]", ident.2 = "Beta INS-hi_EtOH", 
                                       test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                       min.pct = 0,
                                       logfc.threshold = 0, 
                                       pseudocount.use = 1,
                                       verbose = TRUE)
head(beta.INSHi.DHT.response, n = 15)
write.csv(beta.INSHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1beta.INSHi.DHT.response.csv)")

# 2.Beta-cells (INS low)
beta.INSLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Beta INS-low_DHT[10nM]", ident.2 = "Beta INS-low_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0,
                                        logfc.threshold = 0, 
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(beta.INSLow.DHT.response, n = 15)
write.csv(beta.INSLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1beta.INSLow.DHT.response.csv)")

# 3.Alpha-cells (GCG hi)
alpha.GCGHi.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Alpha GCG-hi_DHT[10nM]", ident.2 = "Alpha GCG-hi_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0,
                                        logfc.threshold = 0,
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(alpha.GCGHi.DHT.response, n = 15)
write.csv(alpha.GCGHi.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1alpha.GCGHi.DHT.response.csv)")

# 4.Alpha-cells (GCG low)
alpha.GCGLow.DHT.response <- FindMarkers(pancreas.integrated, 
                                         ident.1 = "Alpha GCG-low_DHT[10nM]", ident.2 = "Alpha GCG-low_EtOH", 
                                         test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                         min.pct = 0,
                                         logfc.threshold = 0, 
                                         pseudocount.use = 1,
                                         verbose = TRUE)
head(alpha.GCGLow.DHT.response, n = 15)
write.csv(alpha.GCGLow.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1alpha.GCGLow.DHT.response.csv)")

# 5.Trandifferentiating Endocrine-Cells
tranbeta.DHT.response <- FindMarkers(pancreas.integrated, 
                                     ident.1 = "Transdifferentiating Beta_DHT[10nM]", ident.2 = "Transdifferentiating Beta_EtOH", 
                                     test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                     min.pct = 0,
                                     logfc.threshold = 0,
                                     pseudocount.use = 1,
                                     verbose = TRUE)
head(tranbeta.DHT.response, n = 15)
write.csv(tranbeta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1tranbeta.DHT.response.csv)")

# 6.Delta-Cells
delta.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Delta_DHT[10nM]", ident.2 = "Delta_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0,
                                  logfc.threshold = 0, 
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(delta.DHT.response, n = 15)
write.csv(delta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1delta.DHT.response.csv)")

# 7.Gamma-Cells
gamma.DHT.response <- FindMarkers(pancreas.integrated, 
                                  ident.1 = "Gamma_DHT[10nM]", ident.2 = "Gamma_EtOH", 
                                  test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                  min.pct = 0,
                                  logfc.threshold = 0, 
                                  pseudocount.use = 1,
                                  verbose = TRUE)
head(gamma.DHT.response, n = 15)
write.csv(gamma.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1gamma.DHT.response.csv)")

# 7.Epsilon-Cells
epsilon.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "Epsilon_DHT[10nM]", ident.2 = "Epsilon_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0,
                                    logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(epsilon.DHT.response, n = 15)
write.csv(epsilon.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\epsilon.DHT.response.csv)")

# 8.Ductal-Cells
ductal.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Ductal_DHT[10nM]", ident.2 = "Ductal_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0,
                                   logfc.threshold = 0, 
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(ductal.DHT.response, n = 15)
write.csv(ductal.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1ductal.DHT.response.csv)")

# 9.Acinar-Cells
acinar.DHT.response <- FindMarkers(pancreas.integrated, 
                                   ident.1 = "Acinar_DHT[10nM]", ident.2 = "Acinar_EtOH", 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0,
                                   logfc.threshold = 0, 
                                   pseudocount.use = 1,
                                   verbose = TRUE)
head(acinar.DHT.response, n = 15)
write.csv(acinar.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1acinar.DHT.response.csv)")

# 10.Quiescent Stellate-Cells
qstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Quiescent Stellate_DHT[10nM]", ident.2 = "Quiescent Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0,
                                      logfc.threshold = 0, 
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(qstellate.DHT.response, n = 15)
write.csv(qstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1qstellate.DHT.response.csv)")

# 11.Activated Stellate-Cells
astellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Activated Stellate_DHT[10nM]", ident.2 = "Activated Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0,
                                      logfc.threshold = 0, 
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(astellate.DHT.response, n = 15)
write.csv(astellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1astellate.DHT.response.csv)")

# 12.Proliferating Stellate-Cells
pstellate.DHT.response <- FindMarkers(pancreas.integrated, 
                                      ident.1 = "Proliferating Stellate_DHT[10nM]", ident.2 = "Proliferating Stellate_EtOH", 
                                      test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                      min.pct = 0,
                                      logfc.threshold = 0, 
                                      pseudocount.use = 1,
                                      verbose = TRUE)
head(pstellate.DHT.response, n = 15)
write.csv(pstellate.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1pstellate.DHT.response.csv)")

# 13.Macrophage-Cells
macrophage.DHT.response <- FindMarkers(pancreas.integrated, 
                                       ident.1 = "Macrophage_DHT[10nM]", ident.2 = "Macrophage_EtOH", 
                                       test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                       min.pct = 0,
                                       logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                       pseudocount.use = 1,
                                       verbose = TRUE)
head(macrophage.DHT.response, n = 15)
write.csv(macrophage.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1macrophage.DHT.response.csv)")

# 14.T Lymphocyte-Cells
tlympho.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "T-Lymphocyte_DHT[10nM]", ident.2 = "T-Lymphocyte_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0,
                                    logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(tlympho.DHT.response, n = 15)
write.csv(tlympho.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1tlympho.DHT.response.csv)")

# 15.Mast-Cells
mast.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Mast_DHT[10nM]", ident.2 = "Mast_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0,
                                 logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = TRUE)
head(mast.DHT.response, n = 15)
write.csv(mast.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1mast.DHT.response.csv)")

# 16.Schwann-Cells
schwann.DHT.response <- FindMarkers(pancreas.integrated, 
                                    ident.1 = "Schwann_DHT[10nM]", ident.2 = "Schwann_EtOH", 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0,
                                    logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                    pseudocount.use = 1,
                                    verbose = TRUE)
head(schwann.DHT.response, n = 15)
write.csv(schwann.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1schwann.DHT.response.csv)")

# 17.Endothelial-Cells
endothelial.DHT.response <- FindMarkers(pancreas.integrated, 
                                        ident.1 = "Endothelial_DHT[10nM]", ident.2 = "Endothelial_EtOH", 
                                        test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                        min.pct = 0,
                                        logfc.threshold = 0, # based on output log2 so 0.137504 is ~1.1 FC
                                        pseudocount.use = 1,
                                        verbose = TRUE)
head(endothelial.DHT.response, n = 15)
write.csv(endothelial.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\1endothelial.DHT.response.csv)")


# Plotting DE genes ###
Idents(pancreas.integrated) <- "celltype"
beta.cells <- subset(pancreas.integrated, idents = "Beta INS-hi")
plots <- VlnPlot(beta.cells, features = c("INS", "DDIT3", "MIF", "DEPP1", "PLCG2", "IAPP"), group.by = "treatment", 
                 pt.size = 0, combine = TRUE)
plots <- VlnPlot(beta.cells, features = c("MT-CO3", "MT-ND1", "MT-ND4", "MT-ATP6", "MT-CO1", "MT-CYB"), group.by = "treatment", 
                 pt.size = 0, combine = TRUE)
wrap_plots(plots = plots, nrow = 1, ncol = 1)

# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1alpha.GCGHi.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('ATP5F1E', 'ATP1B1', 'COX17', 'MT-ND4L',
                              'ATP5MC1', 'ATP5MD', 'COX7A1', 'ATP5ME', 'UQCR10',
                              'COX7A2', 'NDUFC1', 'LDHA', 'COX8A', 'COX7B', 'COX6C', 'NDUFC2', 
                              'NDUFAB1', 'UQCRQ',
                              'MT-CO3', 'MT-ND1','MT-ND4', 'MT-Co1', 'MT-ATP6', 'MT-CYB', 'MT-ND2',
                              'MT-CO2', 'MT-ND5', 'MT-ND3', 'PCSK', 'SOD2', 'MT-ND6', 'ATP9A'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1,0.8),
                ylim = c(0,300),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))


# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1alpha.GCGLow.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('PCSK1N', 'COX4I1', 'ATP5F1B', 'MT-CO3',
                              'MT-ND4L', 'MT-ATP8', 'MT-ND5', 'COX20', 'MT-ND6',
                              'ABCC8', 'ATP1B1', 'ATP5MC1'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1,1.6),
                ylim = c(0,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))

# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1beta.INSHi.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('TPI1', 'IAPP', 'PGK1', 'NPY', 'CTNNB1', 'VEGFA', 'PKM',
                              'MT-CO3', 'MT-ND1', 'MT-ND4', 'MT-CO1', 'MT-ATP6', 'MT-CYB', 'MT-ND2', 'MT-CO2', 'MT-ND3', 'MT-ND6', 'ABCC8',
                              'MAFA', 'PDX1', 'PCSK1', 'MAFB', 'ACTB', 'GSN'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1.5,1.5),
                ylim = c(0,300),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))

# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1delta.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c(''), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-2.2,2.2),
                ylim = c(0,50),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                lengthConnectors = unit(0.02, 'npc'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))


# Load data
volcanodat <- read.csv(r"(C:\Users\mqadir\Box\!FAHD\1. AR-DHT Project\DHT_scRNAseq_Islets\1. DGE_analysis\1. All Genes\1alpha.GCGLow.DHT.response.csv)",
                       header = TRUE, sep = ",", row.names = 1)
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$avg_log2FC > 0.137504 & volcanodat$p_val < 0.05)] <- 'high'

# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$avg_log2FC < -0.137504 & volcanodat$p_val < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'avg_log2FC',
                y = 'p_val',
                #selectLab = FALSE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                selectLab = c('PCSK1N', 'COX4I1', 'ATP5F1B', 'MT-CO3',
                              'MT-ND4L', 'MT-ATP8', 'MT-ND5', 'COX20', 'MT-ND6',
                              'ABCC8', 'ATP1B1', 'ATP5MC1'), # use this for labelling genes on plot
                #boxedLabels = TRUE,
                xlim = c(-1,1.6),
                ylim = c(0,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 0.137504,
                #pointSize = c(ifelse(volcanodat$avg_log2FC < -1 | volcanodat$avg_log2FC > 1, 6, 4)),
                pointSize = 3,
                labSize = 5,
                labFace = 'bold',
                #boxedLabels = TRUE,
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))





# Calculating percentages
X1 <- NULL
table(x = FetchData(pancreas.integrated, vars = c('celltype', 'sample')))
x1 <- subset(pancreas.integrated, subset = (celltype == c("beta", "alpha")) & (sex == "Male"))
table(x = FetchData(x1, vars = c('celltype', 'sex')))
x2 <- subset(pancreas.integrated, subset = (celltype != c("beta")) & (sex != "Male")) # wont run because you cant subset a vector with no cells which is what is left once all male cells are removed :)
table(x = FetchData(x2, vars = c('celltype', 'sex')))






theme_set(theme_cowplot())
beta.cells <- subset(pancreas.integrated, idents = "beta")
Idents(beta.cells) <- "treatment"
avg.beta.cells <- log1p(AverageExpression(beta.cells, verbose = FALSE)$RNA)
avg.beta.cells$gene <- rownames(avg.beta.cells)

alpha.cells <- subset(pancreas.integrated, idents = "alpha")
Idents(alpha.cells) <- "treatment"
avg.alpha.cells <- log1p(AverageExpression(alpha.cells, verbose = FALSE)$RNA)
avg.alpha.cells$gene <- rownames(avg.alpha.cells)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.beta.cells, aes("ctrl", "DHT[10nM]")) + geom_point() + ggtitle("Beta Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)

# SAVING GEM ####
pancreas_integrated_GEM <- GetAssayData(object = pancreas.integrated, slot = "counts")
write.csv(pancreas_integrated_GEM, file = r"(C:\Users\mqadir\Box\FMJ lab\2. Ongoing papers\BARKO paper\GEO Upload\pancreas_integrated_GEM.csv)")

# Extra
# Installation of the latest released version
install.packages('GOplot')
library(GOplot)
packageVersion("GOplot")

gene.data <- read.csv("D:/R-Projects/DHT/Data output/beta.DHT.de.data.csv")
up.go <- read.csv("D:/R-Projects/DHT/Data output/up/DHTGOup.csv")
down.go <- read.csv("D:/R-Projects/DHT/Data output/Down/DHTGOdown.csv")

head(gene.data)
head(up.go)
circ <- circle_dat(up.go, gene.data)
circ
process <- List('cellular response to decreased oxygen levels', "cellular response to hypoxia",
                'response to unfolded protein', 'canonical glycolysis',
                'glucose catabolic process to pyruvate', 'glycolytic process through glucose-6-phosphate',
                'gluconeogenesis', 'cellular response to oxidative stress',
                'protein stabilization ', 'amino acid transport ')

process
chord <- chord_dat(data = circ, genes = gene.data)
chord <- chord_dat(data = circ, process = process)
chord <- chord_dat(data = circ, genes = gene.data, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)


GOBubble(circ, labels = 1)


