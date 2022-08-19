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
  }
)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle3")
packageVersion("harmony")

# Enable parallelization
plan("multiprocess", workers = 1)

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

# Step 5: Merge Datasets
# Based on comment to Issue #4753 https://github.com/satijalab/seurat/issues/4753
# we first integrate all of the similar samples
{
  whitemale.list <- list("HP2022801" = HP2022801, "SAMN15877725" = SAMN15877725, "HP2107001" = HP2107001, "HP2107901" = HP2107901
                         )
  
  whitefemale.list <- list("HP2024001" = HP2024001, "HP2105501" = HP2105501, "HP2108601" = HP2108601, "HP2108901" = HP2108901 
                           )
  
  blackmale.list <- list("HP2031401" = HP2031401, "HP2110001" = HP2110001, "HP2123201" = HP2123201
                         )
  
  blackfemale.list <- list("HP2106201" = HP2106201, "HP2121601" = HP2121601, "HP2132801" = HP2132801, "HP2202101" = HP2202101
                           )

# Normalization using SCT look at https://github.com/satijalab/seurat/issues/4623
  whitemale.list <- lapply(X = whitemale.list,
                           FUN = function(x) {
                             x <- SCTransform(x, vars.to.regress = "percent.mt")
                            })

  whitefemale.list <- lapply(X = whitefemale.list,
                             FUN = function(x) {
                               x <- SCTransform(x, vars.to.regress = "percent.mt")
                              })
  
  blackmale.list <- lapply(X = blackmale.list,
                           FUN = function(x) {
                             x <- SCTransform(x, vars.to.regress = "percent.mt")
                           })
  
  blackfemale.list <- lapply(X = blackfemale.list,
                             FUN = function(x) {
                               x <- SCTransform(x, vars.to.regress = "percent.mt")
                             })

}

# Saving
# saveRDS(whitemale.list, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.list.rds)")


# Staged loading
#pancreas.list <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.list.rds)")

# Find features
{
features.whitemale <- SelectIntegrationFeatures(object.list = whitemale.list, nfeatures = 3000)
features.whitefemale <- SelectIntegrationFeatures(object.list = whitefemale.list, nfeatures = 3000)
features.blackmale <- SelectIntegrationFeatures(object.list = blackmale.list, nfeatures = 3000)
features.blackfemale <- SelectIntegrationFeatures(object.list = blackfemale.list, nfeatures = 3000)

# prep SCT integration
whitemale.list <- PrepSCTIntegration(object.list = whitemale.list, anchor.features = features.whitemale)
whitefemale.list <- PrepSCTIntegration(object.list = whitefemale.list, anchor.features = features.whitefemale)
blackmale.list <- PrepSCTIntegration(object.list = blackmale.list, anchor.features = features.blackmale)
blackfemale.list <- PrepSCTIntegration(object.list = blackfemale.list, anchor.features = features.blackfemale)

# Find integration anchors
anchors.whitemale <- FindIntegrationAnchors(object.list = whitemale.list,
                                            anchor.features = features.whitemale,
                                            normalization.method = "SCT")

anchors.whitefemale <- FindIntegrationAnchors(object.list = whitefemale.list, 
                                              anchor.features = features.whitefemale,
                                              normalization.method = "SCT")

anchors.blackmale <- FindIntegrationAnchors(object.list = blackmale.list, 
                                            anchor.features = features.blackmale,
                                            normalization.method = "SCT")

anchors.blackfemale <- FindIntegrationAnchors(object.list = blackfemale.list, 
                                              anchor.features = features.blackfemale,
                                              normalization.method = "SCT")
}

# Integrate data
whitemale.combined <- IntegrateData(anchorset = anchors.whitemale, normalization.method = "SCT")
whitefemale.combined <- IntegrateData(anchorset = anchors.whitefemale, normalization.method = "SCT")
blackmale.combined <- IntegrateData(anchorset = anchors.blackmale, normalization.method = "SCT")
blackfemale.combined <- IntegrateData(anchorset = anchors.blackfemale, normalization.method = "SCT")

# Save individually integrated samples
{
  saveRDS(whitemale.combined, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\whitemale.combined.rds)")
  saveRDS(whitefemale.combined, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\whitefemale.combined.rds)")
  saveRDS(blackmale.combined, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\blackmale.combined.rds)")
  saveRDS(blackfemale.combined, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\blackfemale.combined.rds)")
  }

# Load
{
  whitemale.combined <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\whitemale.combined.rds)")
  whitefemale.combined <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\whitefemale.combined.rds)")
  blackmale.combined <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\blackmale.combined.rds)")
  blackfemale.combined <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\blackfemale.combined.rds)")
  }

# Analysis of new integrated list
pancreas.list <- list("whitemale.combined" = whitemale.combined, "whitefemale.combined"=whitefemale.combined, 
                      "blackmale.combined"=blackmale.combined, "blackfemale.combined"=blackfemale.combined)
# Find features
features.pancreas <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)

# prep SCT integration
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = features.pancreas)

# Find integration anchors
anchors.pancreas <- FindIntegrationAnchors(object.list = pancreas.list,
                                            anchor.features = features.pancreas,
                                            normalization.method = "SCT")

# Integrate data
pancreas.combined <- IntegrateData(anchorset = anchors.pancreas, normalization.method = "SCT")

# PCA
pancreas.combined <- RunPCA(pancreas.combined, verbose = TRUE)
pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", dims = 1:30)
pancreas.combined <- FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:30)
pancreas.combined <- FindClusters(pancreas.combined, resolution = 0.5)
p1 <- DimPlot(pancreas.combined, reduction = "umap", group.by = "ancestry_sex")
p2 <- DimPlot(pancreas.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
              repel = TRUE)
p1 + p2

# Save pancreas.combined.cca
saveRDS(pancreas.combined, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.combined.cca.rds)")


# Run RPCA
{
  whitemale.list <- lapply(X = whitemale.list, FUN = RunPCA, 
                         verbose = TRUE, 
                         features = anchors.wm)
  whitefemale.list <- lapply(X = whitfeemale.list, FUN = RunPCA, 
                           verbose = TRUE, 
                           features = anchors.wf)
  blackmale.list <- lapply(X = blackmale.list, FUN = RunPCA, 
                           verbose = TRUE, 
                           features = anchors.bm)
  whitefemale.list <- lapply(X = whitefemale.list, FUN = RunPCA, 
                             verbose = TRUE, 
                             features = anchors.bf)
  }
