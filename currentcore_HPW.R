# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 03/23/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

# LOAD LIBRARIES ####
# Restart Rstudio or R

install.packages('ggplot2')
install.packages('cowplot')
install.packages('Matrix')
install.packages('ggridges')
install.packages('ggrepel')
install.packages('dplyr')
#install.packages('Seurat')
install.packages('plotly')
install.packages('clustree')
install.packages('patchwork')
install.packages('future')
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
BiocManager::install("EnhancedVolcano")
BiocManager::install("DoubletFinder")
BiocManager::install("glmGamPoi")
BiocManager::install("GOSemSim")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationHub")
BiocManager::install("GenomeInfoDb")
BiocManager::install("MeSHDbi")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("dittoSeq")
BiocManager::install("escape")
BiocManager::install("ComplexHeatmap")
BiocManager::install(c("DropletUtils", "Nebulosa"))
BiocManager::install("hdf5r", force = TRUE)
BiocManager::install("DESeq2")

# install Seurat from Github (automatically updates sctransform)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("gaospecial/ggVennDiagram")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("harmony")
BiocManager::install("EnrichmentBrowser")
install.packages('SoupX')
install.packages('tidyverse')
install.packages("viridis")
install.packages("circlize")
install.packages("scCustomize")
install.packages("devtools")
install.packages("archive")
install.packages("R.utils")
install.packages("qs")


# Run the following code once you have Seurat installed
suppressWarnings(
  {
    library(leiden)
    library(stringr)
    library(hdf5r)
    library(SoupX)
    library(Rcpp)
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(tidyverse)
    library(data.table)
    library(reticulate)
    library(Seurat)
    library(monocle3)
    library(harmony)
    library(Signac)
    library(EnsDb.Hsapiens.v86)
    library(GenomeInfoDb)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
    library(glmGamPoi)
    library(GOSemSim)
    library(org.Hs.eg.db)
    library(AnnotationHub)
    library(MeSHDbi)
    library(clusterProfiler)
    library(DOSE)
    library(dittoSeq)
    library(escape)
    library(EnrichmentBrowser)
    library(viridisLite)
    library(viridis)
    library(ComplexHeatmap)
    library(circlize)
    library(scCustomize)
    library(Nebulosa)
    library(DropletUtils)
    library(ggvenn)
    library(ggVennDiagram)
    library(devtools)
    library(R.utils)
    library(qs)
    library(parallel)
    library(readr)
    library(DESeq2)
    library(beeswarm)
    library(limma)
    library(edgeR)
    library(GenomicFeatures)
    library(data.table)
    library(RColorBrewer)
    library(pheatmap)
  }
)

# Set global environment parameter par-proc
#options(future.globals.maxSize = 8000 * 1024^2)
set.seed(1234)

# Python env
if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/mqadir/AppData/Local/r-miniconda/envs/r-reticulate",Sys.getenv()["PATH"],sep=";"))
py_config()

# WD
setwd(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\WD)")
(WD <- getwd())


# Check package versions
packageVersion("clusterProfiler")
packageVersion("dittoSeq")
packageVersion("escape")
packageVersion("seurat")
packageVersion("signac")
packageVersion("EnrichmentBrowser")

############################ STAGE ############################
############################   1   ############################
# OBJECT SETUP AND NORMALIZATION ###
# STEP 1: Load 
system.time({
#user    system elapsed 
#8328.56 254.45 8669.53 
  {
    samples <- c("HP2022801", "SAMN15877725", "HP2107001", "HP2107901",
                 "HP2024001", "HP2105501", "HP2108601", "HP2108901",
                 "HP2031401", "HP2110001", "HP2123201",
                 "HP2106201", "HP2121601", "HP2132801", "HP2202101")
  }
  
  {
    for (sample in samples){
      wd <- sprintf('E:/1.SexbasedStudyrawdata/Cellranger_raw_data/scRNAseq/%s', samples)
    }
    }
  {
    # Automated SoupX created
    for (x in wd){
    sample_name <- str_split_fixed(x, "/", n=5)[5] #Adjust this to output your sample name
    rna_counts <- Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5'))
    data <- CreateSeuratObject(counts=rna_counts)
    data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
    data <- subset(x = data, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 20)
  
    #Running sctransform takes into account sequencing depth at each cell
    #data <- SCTransform(data, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE, vst.flavor = "v2")
    #data <- SCTransform(data, verbose = FALSE)
  
    #Log normalization alternative to sctransform
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
  
    data <- RunPCA(data, verbose = FALSE)
    data <- RunUMAP(data, dims = 1:30, verbose = FALSE)
    data <- FindNeighbors(data, dims = 1:30, verbose = FALSE)
    data <- FindClusters(data, algorithm=4, resolution = 1, verbose=FALSE)
  
    #Read in RNA assay counts from our filtered seurat object
    DefaultAssay(data) <- 'RNA'
    toc <- GetAssayData(object = data, slot = "counts") #with nFeature >500 filter
    tod <- Seurat::Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5')) #raw count matrix
  
    #Pull out the required metadata from the clustered filtered adata object
    #We need the UMAP coordinates (RD1 and RD2) and the cluster assignments at minimum
    metadata <- (cbind(as.data.frame(data[["umap"]]@cell.embeddings),
                       as.data.frame(Idents(data)),
                       as.data.frame(Idents(data))))
    colnames(metadata) <- c("RD1","RD2","Cluster","Annotation")
  
    #Create the SoupChannel Object
    sc <- SoupChannel(tod,toc)
  
    #Add in the metadata (dimensionality reduction UMAP coords and cluster assignments)
    sc <- setDR(sc,metadata[colnames(sc$toc),c("RD1","RD2")])
    sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))
    sc <- autoEstCont(sc)
    out <- adjustCounts(sc, roundToInt = TRUE) # very important step, use roundToInt = TRUE
  
    #Create post-SoupX Seurat Object
    data2 <- CreateSeuratObject(out)
    data2[['percent.mt']] <- PercentageFeatureSet(data2, pattern = '^MT-')
    #data2 <- NormalizeData(data2, normalization.method = "LogNormalize", scale.factor = 10000)  #Can be changed to sctransform
    #data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 2000)
    #all.genes <- rownames(data2)
    #data2 <- ScaleData(data2, features = all.genes)
    data2 <- SCTransform(data2, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE, vst.flavor = "v2")
    data2 <- RunPCA(data2, verbose = FALSE)
    data2 <- RunUMAP(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindNeighbors(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindClusters(data2, algorithm=4, resolution = 1, verbose=FALSE)
    qsave(data2, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/RDS files/current/1_soupx_rounttoint/%s.qs",sample_name))
    }
  }
})
  

############################ STAGE ############################
############################   2   ############################
# Load data
system.time({
#user     system elapsed 
#30790.28 562.56 31476.69 
{
  HP2022801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2022801.qs)")
  SAMN15877725 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\SAMN15877725.qs)")
  HP2024001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2024001.qs)")
  HP2031401 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2031401.qs)")
  HP2105501 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2105501.qs)")
  HP2106201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2106201.qs)")
  HP2107001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2107001.qs)")
  HP2107901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2107901.qs)")
  HP2108601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2108601.qs)")
  HP2108901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2108901.qs)")
  HP2110001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2110001.qs)")
  HP2121601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2121601.qs)")
  HP2123201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2123201.qs)")
  HP2132801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2132801.qs)")
  HP2202101 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2202101.qs)")
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
                                
  # Ancestry and sex specific UNION Metadata addition
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
  
# Doublet removal
# Optimization
{
  sweep.res <- paramSweep_v3(HP2022801, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2022801) * 0.05)  # expect 4% doublets
  HP2022801 <- doubletFinder_v3(HP2022801, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(SAMN15877725, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(SAMN15877725) * 0.05)  # expect 4% doublets
  SAMN15877725 <- doubletFinder_v3(SAMN15877725, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2024001, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2024001) * 0.05)  # expect 4% doublets
  HP2024001 <- doubletFinder_v3(HP2024001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2031401, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2031401) * 0.05)  # expect 4% doublets
  HP2031401 <- doubletFinder_v3(HP2031401, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2105501, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2105501) * 0.05)  # expect 4% doublets
  HP2105501 <- doubletFinder_v3(HP2105501, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2106201, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2106201) * 0.05)  # expect 4% doublets
  HP2106201 <- doubletFinder_v3(HP2106201, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2107001, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2107001) * 0.05)  # expect 4% doublets
  HP2107001 <- doubletFinder_v3(HP2107001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2107901, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2107901) * 0.05)  # expect 4% doublets
  HP2107901 <- doubletFinder_v3(HP2107901, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2108601, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2108601) * 0.05)  # expect 4% doublets
  HP2108601 <- doubletFinder_v3(HP2108601, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2108901, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2108901) * 0.05)  # expect 4% doublets
  HP2108901 <- doubletFinder_v3(HP2108901, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2110001, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2110001) * 0.05)  # expect 4% doublets
  HP2110001 <- doubletFinder_v3(HP2110001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2121601, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2121601) * 0.05)  # expect 4% doublets
  HP2121601 <- doubletFinder_v3(HP2121601, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2123201, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2123201) * 0.05)  # expect 4% doublets
  HP2123201 <- doubletFinder_v3(HP2123201, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2132801, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2132801) * 0.05)  # expect 4% doublets
  HP2132801 <- doubletFinder_v3(HP2132801, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2202101, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2202101) * 0.05)  # expect 4% doublets
  HP2202101 <- doubletFinder_v3(HP2202101, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  }

# Setup one metadata column
HP2022801$doublets <- HP2022801$DF.classifications_0.25_0.09_213
SAMN15877725$doublets <- SAMN15877725$DF.classifications_0.25_0.09_201
HP2107001$doublets <- HP2107001$DF.classifications_0.25_0.09_212
HP2107901$doublets <- HP2107901$DF.classifications_0.25_0.09_165
HP2024001$doublets <- HP2024001$DF.classifications_0.25_0.09_151
HP2105501$doublets <- HP2105501$DF.classifications_0.25_0.09_156
HP2108601$doublets <- HP2108601$DF.classifications_0.25_0.09_273
HP2108901$doublets <- HP2108901$DF.classifications_0.25_0.09_214
HP2031401$doublets <- HP2031401$DF.classifications_0.25_0.09_233
HP2110001$doublets <- HP2110001$DF.classifications_0.25_0.09_290
HP2123201$doublets <- HP2123201$DF.classifications_0.25_0.09_78
HP2106201$doublets <- HP2106201$DF.classifications_0.25_0.09_325
HP2121601$doublets <- HP2121601$DF.classifications_0.25_0.09_179
HP2132801$doublets <- HP2132801$DF.classifications_0.25_0.09_120
HP2202101$doublets <- HP2202101$DF.classifications_0.25_0.09_200

# Step 4: Add cell IDs ####
# Add cell IDs
{
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

# Save point
{
  qsave(HP2022801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2022801.qs)")
  qsave(SAMN15877725, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\SAMN15877725.qs)")
  qsave(HP2024001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2024001.qs)")
  qsave(HP2031401, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2031401.qs)")
  qsave(HP2105501, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2105501.qs)")
  qsave(HP2106201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2106201.qs)")
  qsave(HP2107001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107001.qs)")
  qsave(HP2107901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107901.qs)")
  qsave(HP2108601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108601.qs)")
  qsave(HP2108901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108901.qs)")
  qsave(HP2110001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2110001.qs)")
  qsave(HP2121601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2121601.qs)")
  qsave(HP2123201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2123201.qs)")
  qsave(HP2132801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2132801.qs)")
  qsave(HP2202101, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2202101.qs)")
}
})

############################ STAGE ############################
############################   3   ############################

# Analysis of Ruth's data
# Hpap dataset
system.time({
  #user    system elapsed 
  #9294.29 614.52 9908.50 (~2.75hrs)
hpap <- readRDS(file=r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\Raw_hpap_file\hpap_islet_scRNAseq.rds)")
table(hpap@meta.data[["Library"]])

Idents(hpap) <- "Library"
hpap$ancestry <- plyr::mapvalues(
  x= hpap$Library,
  from = c("HPAP-019",
           "HPAP-020",
           "HPAP-021",
           "HPAP-022",
           "HPAP-023",
           "HPAP-024",
           "HPAP-026",
           "HPAP-028",
           "HPAP-029",
           "HPAP-032",
           "HPAP-034",
           "HPAP-035",
           "HPAP-036",
           "HPAP-037",
           "HPAP-038",
           "HPAP-039",
           "HPAP-040",
           "HPAP-042",
           "HPAP-043",
           "HPAP-044",
           "HPAP-045",
           "HPAP-047",
           "HPAP-049",
           "HPAP-050",
           "HPAP-051",
           "HPAP-052",
           "HPAP-053",
           "HPAP-054",
           "HPAP-055",
           "HPAP-056",
           "HPAP-057",
           "HPAP-058",
           "HPAP-059",
           "HPAP-061",
           "HPAP-063",
           "HPAP-064",
           "HPAP-065",
           "HPAP-070",
           "HPAP-071",
           "HPAP-072",
           "HPAP-074",
           "HPAP-075",
           "HPAP-077",
           "HPAP-079",
           "HPAP-080",
           "HPAP-081",
           "HPAP-082",
           "HPAP-083",
           "HPAP-084",
           "HPAP-085",
           "HPAP-087",
           "HPAP-088",
           "HPAP-090",
           "HPAP-091",
           "HPAP-092",
           "HPAP-099",
           "HPAP-100",
           "HPAP-101",
           "HPAP-103",
           "HPAP-104",
           "HPAP-105",
           "HPAP-106",
           "HPAP-107",
           "HPAP-108",
           "HPAP-109"),
  to = c("white", #"HPAP-019",
         "white", #"HPAP-020",
         "white", #"HPAP-021",
         "white", #"HPAP-022",
         "white", #"HPAP-023",
         "white", #"HPAP-024",
         "white", #"HPAP-026",
         "white", #"HPAP-028",
         "white", #"HPAP-029",
         "white", #"HPAP-032",
         "white", #"HPAP-034",
         "white", #"HPAP-035",
         "white", #"HPAP-036",
         "white", #"HPAP-037",
         "white", #"HPAP-038",
         "white", #"HPAP-039",
         "white", #"HPAP-040",
         "white", #"HPAP-042",
         "hispanic", #"HPAP-043",
         "white", #"HPAP-044",
         "white", #"HPAP-045",
         "white", #"HPAP-047",
         "white",  #"HPAP-049",
         "hispanic", #"HPAP-050",
         "black", #"HPAP-051",
         "black", #"HPAP-052",
         "white", #"HPAP-053",
         "white", #"HPAP-054",
         "hispanic", #"HPAP-055",
         "white", #"HPAP-056",
         "white", #"HPAP-057",
         "black", #"HPAP-058",
         "white", #"HPAP-059",
         "black", #"HPAP-061",
         "white", #"HPAP-063",
         "black", #"HPAP-064",
         "white", #"HPAP-065",
         "black", #"HPAP-070",
         "white", #"HPAP-071",
         "hispanic", #"HPAP-072",
         "white", #"HPAP-074",
         "white", #"HPAP-075",
         "white", #"HPAP-077",
         "hispanic", #"HPAP-079",
         "hispanic", #"HPAP-080",
         "white", #"HPAP-081",
         "white", #"HPAP-082",
         "black", #"HPAP-083",
         "white", #"HPAP-084",
         "white", #"HPAP-085",
         "white", #"HPAP-087",
         "white", #"HPAP-088",
         "white", #"HPAP-090",
         "hispanic", #"HPAP-091",
         "hispanic", #"HPAP-092",
         "hispanic",  #"HPAP-099",
         "white", #"HPAP-100",
         "hispanic", #"HPAP-101",
         "white", #"HPAP-103",
         "hispanic", #"HPAP-104",
         "hispanic", #"HPAP-105",
         "white", #"HPAP-106",
         "white", #"HPAP-107",
         "black", #"HPAP-108",
         "hispanic" #"HPAP-109"
  )
)

table(hpap@meta.data[["Library"]])
table(hpap@meta.data[["Sex"]])
table(hpap@meta.data[["ancestry"]])

# Subsetting data
Idents(hpap) <- "Library"
{
  HPAP_019 <- subset(hpap, idents = c("HPAP-019"))
  HPAP_020 <- subset(hpap, idents = c("HPAP-020"))
  HPAP_021 <- subset(hpap, idents = c("HPAP-021"))
  HPAP_022 <- subset(hpap, idents = c("HPAP-022"))
  HPAP_023 <- subset(hpap, idents = c("HPAP-023"))
  HPAP_024 <- subset(hpap, idents = c("HPAP-024"))
  HPAP_026 <- subset(hpap, idents = c("HPAP-026"))
  HPAP_028 <- subset(hpap, idents = c("HPAP-028"))
  HPAP_029 <- subset(hpap, idents = c("HPAP-029"))
  HPAP_032 <- subset(hpap, idents = c("HPAP-032"))
  HPAP_034 <- subset(hpap, idents = c("HPAP-034"))
  HPAP_035 <- subset(hpap, idents = c("HPAP-035"))
  HPAP_036 <- subset(hpap, idents = c("HPAP-036"))
  HPAP_037 <- subset(hpap, idents = c("HPAP-037"))
  HPAP_038 <- subset(hpap, idents = c("HPAP-038"))
  HPAP_039 <- subset(hpap, idents = c("HPAP-039"))
  HPAP_040 <- subset(hpap, idents = c("HPAP-040"))
  HPAP_042 <- subset(hpap, idents = c("HPAP-042"))
  HPAP_043 <- subset(hpap, idents = c("HPAP-043"))
  HPAP_044 <- subset(hpap, idents = c("HPAP-044"))
  HPAP_045 <- subset(hpap, idents = c("HPAP-045"))
  HPAP_047 <- subset(hpap, idents = c("HPAP-047"))
  HPAP_049 <- subset(hpap, idents = c("HPAP-049"))
  HPAP_050 <- subset(hpap, idents = c("HPAP-050"))
  HPAP_051 <- subset(hpap, idents = c("HPAP-051"))
  HPAP_052 <- subset(hpap, idents = c("HPAP-052"))
  HPAP_053 <- subset(hpap, idents = c("HPAP-053"))
  HPAP_054 <- subset(hpap, idents = c("HPAP-054"))
  HPAP_055 <- subset(hpap, idents = c("HPAP-055"))
  HPAP_056 <- subset(hpap, idents = c("HPAP-056"))
  HPAP_057 <- subset(hpap, idents = c("HPAP-057"))
  HPAP_058 <- subset(hpap, idents = c("HPAP-058"))
  HPAP_059 <- subset(hpap, idents = c("HPAP-059"))
  HPAP_061 <- subset(hpap, idents = c("HPAP-061"))
  HPAP_063 <- subset(hpap, idents = c("HPAP-063"))
  HPAP_064 <- subset(hpap, idents = c("HPAP-064"))
  HPAP_065 <- subset(hpap, idents = c("HPAP-065"))
  HPAP_070 <- subset(hpap, idents = c("HPAP-070"))
  HPAP_071 <- subset(hpap, idents = c("HPAP-071"))
  HPAP_072 <- subset(hpap, idents = c("HPAP-072"))
  HPAP_074 <- subset(hpap, idents = c("HPAP-074"))
  HPAP_075 <- subset(hpap, idents = c("HPAP-075"))
  HPAP_077 <- subset(hpap, idents = c("HPAP-077"))
  HPAP_079 <- subset(hpap, idents = c("HPAP-079"))
  HPAP_080 <- subset(hpap, idents = c("HPAP-080"))
  HPAP_081 <- subset(hpap, idents = c("HPAP-081"))
  HPAP_082 <- subset(hpap, idents = c("HPAP-082"))
  HPAP_083 <- subset(hpap, idents = c("HPAP-083"))
  HPAP_084 <- subset(hpap, idents = c("HPAP-084"))
  HPAP_085 <- subset(hpap, idents = c("HPAP-085"))
  HPAP_087 <- subset(hpap, idents = c("HPAP-087"))
  HPAP_088 <- subset(hpap, idents = c("HPAP-088"))
  HPAP_090 <- subset(hpap, idents = c("HPAP-090"))
  HPAP_091 <- subset(hpap, idents = c("HPAP-091"))
  HPAP_092 <- subset(hpap, idents = c("HPAP-092"))
  HPAP_099 <- subset(hpap, idents = c("HPAP-099"))
  HPAP_100 <- subset(hpap, idents = c("HPAP-100"))
  HPAP_101 <- subset(hpap, idents = c("HPAP-101"))
  HPAP_103 <- subset(hpap, idents = c("HPAP-103"))
  HPAP_104 <- subset(hpap, idents = c("HPAP-104"))
  HPAP_105 <- subset(hpap, idents = c("HPAP-105"))
  HPAP_106 <- subset(hpap, idents = c("HPAP-106"))
  HPAP_107 <- subset(hpap, idents = c("HPAP-107"))
  HPAP_108 <- subset(hpap, idents = c("HPAP-108"))
  HPAP_109 <- subset(hpap, idents = c("HPAP-109"))
}

pancreas.list <- list("HPAP_019" = HPAP_019, "HPAP_020" = HPAP_020, "HPAP_021" = HPAP_021, "HPAP-022" = HPAP_022, "HPAP_023" = HPAP_023,
                      "HPAP_024" = HPAP_024, "HPAP_026" = HPAP_026, "HPAP_028" = HPAP_028, "HPAP_029" = HPAP_029, "HPAP_032" = HPAP_032,
                      "HPAP_034" = HPAP_034, "HPAP_035" = HPAP_035, "HPAP_036" = HPAP_036, "HPAP_037" = HPAP_037, "HPAP_038" = HPAP_038,
                      "HPAP_039" = HPAP_039, "HPAP_040" = HPAP_040, "HPAP_042" = HPAP_042, "HPAP_043" = HPAP_043, "HPAP_044" = HPAP_044,
                      "HPAP_045" = HPAP_045, "HPAP_047" = HPAP_047, "HPAP_049" = HPAP_049, "HPAP_050" = HPAP_050, "HPAP_051" = HPAP_051,
                      "HPAP_052" = HPAP_052, "HPAP_053" = HPAP_053, "HPAP_054" = HPAP_054, "HPAP_055" = HPAP_055, "HPAP_056" = HPAP_056,
                      "HPAP_057" = HPAP_057, "HPAP_058" = HPAP_058, "HPAP_059" = HPAP_059, "HPAP_061" = HPAP_061, "HPAP_063" = HPAP_063,
                      "HPAP_064" = HPAP_064, "HPAP_065" = HPAP_065, "HPAP_070" = HPAP_070, "HPAP_071" = HPAP_071, "HPAP_072" = HPAP_072,
                      "HPAP_074" = HPAP_074, "HPAP_075" = HPAP_075, "HPAP_077" = HPAP_077, "HPAP_079" = HPAP_079, "HPAP_080" = HPAP_080,
                      "HPAP_081" = HPAP_081, "HPAP_082" = HPAP_082, "HPAP_083" = HPAP_083, "HPAP_084" = HPAP_084, "HPAP_085" = HPAP_085,
                      "HPAP_087" = HPAP_087, "HPAP_088" = HPAP_088, "HPAP_090" = HPAP_090, "HPAP_091" = HPAP_091, "HPAP_092" = HPAP_092,
                      "HPAP_099" = HPAP_099, "HPAP_100" = HPAP_100, "HPAP_101" = HPAP_101, "HPAP_103" = HPAP_103, "HPAP_104" = HPAP_104,
                      "HPAP_105" = HPAP_105, "HPAP_106" = HPAP_106, "HPAP_107" = HPAP_107, "HPAP_108" = HPAP_108, "HPAP_109" = HPAP_109
)

# Pulling data for analysis
pancreas_subset <- pancreas.list[c("HPAP-022", "HPAP_026", "HPAP_035", "HPAP_036", "HPAP_037", "HPAP_040", "HPAP_044", "HPAP_051", 
                                     "HPAP_052", "HPAP_053", "HPAP_054", "HPAP_056", "HPAP_057", "HPAP_058", "HPAP_059", "HPAP_061", 
                                     "HPAP_063", "HPAP_065", "HPAP_070", "HPAP_074", "HPAP_075", "HPAP_077", "HPAP_079", "HPAP_080", 
                                     "HPAP_081", "HPAP_082", "HPAP_083", "HPAP_085", "HPAP_088", "HPAP_090", "HPAP_091", "HPAP_099",
                                     "HPAP_100", "HPAP_101", "HPAP_103", "HPAP_105", "HPAP_106", "HPAP_108", "HPAP_109")
]

pancreas_subset <- lapply(X = pancreas_subset, 
                            FUN = function(x){
                              x[['percent.mt']] <- PercentageFeatureSet(x, 
                                                                        pattern = '^MT-')
                              x <- SCTransform(x,
                                               vars.to.regress = "percent.mt", 
                                               verbose = TRUE, 
                                               return.only.var.genes = FALSE, 
                                               vst.flavor = "v2")
                            })

# Read in SoupX corrected/doublet classified data
{
  HP2022801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2022801.qs)")
  SAMN15877725 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\SAMN15877725.qs)")
  HP2024001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2024001.qs)")
  HP2031401 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2031401.qs)")
  HP2105501 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2105501.qs)")
  HP2106201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2106201.qs)")
  HP2107001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107001.qs)")
  HP2107901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107901.qs)")
  HP2108601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108601.qs)")
  HP2108901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108901.qs)")
  HP2110001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2110001.qs)")
  HP2121601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2121601.qs)")
  HP2123201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2123201.qs)")
  HP2132801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2132801.qs)")
  HP2202101 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2202101.qs)")
}


# Adjust and add metadata
# Sex specific Metadata addition
{
  HP2022801$Sex <- "M"
  SAMN15877725$Sex <- "M"
  HP2024001$Sex <- "F"
  HP2031401$Sex <- "M"
  HP2105501$Sex <- "F"
  HP2106201$Sex <- "F"
  HP2107001$Sex <- "M"
  HP2107901$Sex <- "M"
  HP2108601$Sex <- "F"
  HP2108901$Sex <- "F"
  HP2110001$Sex <- "M"
  HP2121601$Sex <- "F"
  HP2123201$Sex <- "M"
  HP2132801$Sex <- "F"
  HP2202101$Sex <- "F"
  
  
  HP2022801@meta.data[["Tissue Source"]] <- "Tulane"
  SAMN15877725@meta.data[["Tissue Source"]] <- "Tulane"
  HP2024001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2031401@meta.data[["Tissue Source"]] <- "Tulane"
  HP2105501@meta.data[["Tissue Source"]] <- "Tulane"
  HP2106201@meta.data[["Tissue Source"]] <- "Tulane"
  HP2107001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2107901@meta.data[["Tissue Source"]] <- "Tulane"
  HP2108601@meta.data[["Tissue Source"]] <- "Tulane"
  HP2108901@meta.data[["Tissue Source"]] <- "Tulane"
  HP2110001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2121601@meta.data[["Tissue Source"]] <- "Tulane"
  HP2123201@meta.data[["Tissue Source"]] <- "Tulane"
  HP2132801@meta.data[["Tissue Source"]] <- "Tulane"
  HP2202101@meta.data[["Tissue Source"]] <- "Tulane"
  
  
  HP2022801@meta.data[["Library"]] <- deparse(substitute(HP2022801))
  SAMN15877725@meta.data[["Library"]] <- deparse(substitute(SAMN15877725))
  HP2024001@meta.data[["Library"]] <- deparse(substitute(HP2024001))
  HP2031401@meta.data[["Library"]] <- deparse(substitute(HP2031401))
  HP2105501@meta.data[["Library"]] <- deparse(substitute(HP2105501))
  HP2106201@meta.data[["Library"]] <- deparse(substitute(HP2106201))
  HP2107001@meta.data[["Library"]] <- deparse(substitute(HP2107001))
  HP2107901@meta.data[["Library"]] <- deparse(substitute(HP2107901))
  HP2108601@meta.data[["Library"]] <- deparse(substitute(HP2108601))
  HP2108901@meta.data[["Library"]] <- deparse(substitute(HP2108901))
  HP2110001@meta.data[["Library"]] <- deparse(substitute(HP2110001))
  HP2121601@meta.data[["Library"]] <- deparse(substitute(HP2121601))
  HP2123201@meta.data[["Library"]] <- deparse(substitute(HP2123201))
  HP2132801@meta.data[["Library"]] <- deparse(substitute(HP2132801))
  HP2202101@meta.data[["Library"]] <- deparse(substitute(HP2202101))
  
  
  HP2022801@meta.data[["sample"]] <- NULL
  SAMN15877725@meta.data[["sample"]] <- NULL
  HP2024001@meta.data[["sample"]] <- NULL
  HP2031401@meta.data[["sample"]] <- NULL
  HP2105501@meta.data[["sample"]] <- NULL
  HP2106201@meta.data[["sample"]] <- NULL
  HP2107001@meta.data[["sample"]] <- NULL
  HP2107901@meta.data[["sample"]] <- NULL
  HP2108601@meta.data[["sample"]] <- NULL
  HP2108901@meta.data[["sample"]] <- NULL
  HP2110001@meta.data[["sample"]] <- NULL
  HP2121601@meta.data[["sample"]] <- NULL
  HP2123201@meta.data[["sample"]] <- NULL
  HP2132801@meta.data[["sample"]] <- NULL
  HP2202101@meta.data[["sample"]] <- NULL
  
  
  HP2022801@meta.data[["Chemistry"]] <- "10Xv3"
  SAMN15877725@meta.data[["Chemistry"]] <- "10Xv3"
  HP2024001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2031401@meta.data[["Chemistry"]] <- "10Xv3"
  HP2105501@meta.data[["Chemistry"]] <- "10Xv3"
  HP2106201@meta.data[["Chemistry"]] <- "10Xv3"
  HP2107001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2107901@meta.data[["Chemistry"]] <- "10Xv3"
  HP2108601@meta.data[["Chemistry"]] <- "10Xv3"
  HP2108901@meta.data[["Chemistry"]] <- "10Xv3"
  HP2110001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2121601@meta.data[["Chemistry"]] <- "10Xv3"
  HP2123201@meta.data[["Chemistry"]] <- "10Xv3"
  HP2132801@meta.data[["Chemistry"]] <- "10Xv3"
  HP2202101@meta.data[["Chemistry"]] <- "10Xv3"
  
  
  HP2022801@meta.data[["Diabetes Status"]] <- "ND"
  SAMN15877725@meta.data[["Diabetes Status"]] <- "ND"
  HP2024001@meta.data[["Diabetes Status"]] <- "ND"
  HP2031401@meta.data[["Diabetes Status"]] <- "ND"
  HP2105501@meta.data[["Diabetes Status"]] <- "ND"
  HP2106201@meta.data[["Diabetes Status"]] <- "ND"
  HP2107001@meta.data[["Diabetes Status"]] <- "ND"
  HP2107901@meta.data[["Diabetes Status"]] <- "ND"
  HP2108601@meta.data[["Diabetes Status"]] <- "ND"
  HP2108901@meta.data[["Diabetes Status"]] <- "ND"
  HP2110001@meta.data[["Diabetes Status"]] <- "ND"
  HP2121601@meta.data[["Diabetes Status"]] <- "ND"
  HP2123201@meta.data[["Diabetes Status"]] <- "ND"
  HP2132801@meta.data[["Diabetes Status"]] <- "ND"
  HP2202101@meta.data[["Diabetes Status"]] <- "ND"
}

# Make list
pancreas_qadir <- list("HP2022801" = HP2022801, "SAMN15877725" = SAMN15877725, "HP2107001" = HP2107001, "HP2107901" = HP2107901,
                       "HP2024001" = HP2024001, "HP2105501" = HP2105501, "HP2108601" = HP2108601, "HP2108901" = HP2108901, 
                       "HP2031401" = HP2031401, "HP2110001" = HP2110001, "HP2123201" = HP2123201,
                       "HP2106201" = HP2106201, "HP2121601" = HP2121601, "HP2132801" = HP2132801, "HP2202101" = HP2202101
)

# Subset out all single cells
pancreas_qadir <- lapply(X = pancreas_qadir, 
                         FUN = function(x){
                           
                           # Subset out all singlets
                           Idents(x) <- "doublets"
                           x <- subset(x = x, idents = c("Singlet"))
                         })

# Join datasets
pancreas_combined <- c(pancreas_subset, pancreas_qadir)

# merge data sets
sampleset <- c("HPAP-022", "HPAP_026", "HPAP_035", "HPAP_036", "HPAP_037", "HPAP_040", "HPAP_044", "HPAP_051", 
               "HPAP_052", "HPAP_053", "HPAP_054", "HPAP_056", "HPAP_057", "HPAP_058", "HPAP_059", "HPAP_061", 
               "HPAP_063", "HPAP_065", "HPAP_070", "HPAP_074", "HPAP_075", "HPAP_077", "HPAP_079", "HPAP_080", 
               "HPAP_081", "HPAP_082", "HPAP_083", "HPAP_085", "HPAP_088", "HPAP_090", "HPAP_091", "HPAP_099",
               "HPAP_100", "HPAP_101", "HPAP_103", "HPAP_105", "HPAP_106", "HPAP_108", "HPAP_109", "HP2022801", 
               "SAMN15877725", "HP2107001", "HP2107901", "HP2024001", "HP2105501", "HP2108601", "HP2108901", 
               "HP2031401", "HP2110001", "HP2123201", "HP2106201", "HP2121601", "HP2132801", "HP2202101")

pancreas_rna <- merge(pancreas_combined[[sampleset[[1]]]], 
                     y=pancreas_combined[sampleset[2:length(sampleset)]], 
                     project='pancreas', 
                     merge.data = TRUE)

# Inspect data
Idents(pancreas_rna) <- "Tissue Source"
VlnPlot(object = pancreas_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Subset
pancreas_rna <- subset(x = pancreas_rna, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 15)

# Inspect data
pancreas_rna
table(pancreas_rna@meta.data[["Chemistry"]])
table(pancreas_rna@meta.data[["Sex"]])
table(pancreas_rna@meta.data[["Library"]])
table(pancreas_rna@meta.data[["Tissue Source"]])
table(pancreas_rna@meta.data[["Diabetes Status"]])
table(pancreas_rna@meta.data[["ancestry"]])
table(pancreas_rna@meta.data[["doublets"]])
table(pancreas_rna@meta.data[["Cell Type"]])

# #Split data on basis of disease status and calculate integration features on object list
# DefaultAssay(pancreas_rna) <- "SCT"
# integrationfeatures <- SelectIntegrationFeatures(pancreas_combined, nfeatures = 3000, verbose = TRUE)
# 
# # replacing variable features with integrationfeatures
# VariableFeatures(pancreas_rna, assay = "SCT") <- integrationfeatures
# VariableFeatures(pancreas_rna, assay = "RNA") <- integrationfeatures


#Perform basic threshold filtering and log normalization
DefaultAssay(pancreas_rna) <- "RNA"
pancreas_rna <- NormalizeData(pancreas_rna, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas_rna <- FindVariableFeatures(pancreas_rna, selection.method = "vst", nfeatures = 2000)
pancreas_rna <- ScaleData(pancreas_rna, verbose = FALSE) %>% 
  RunPCA(pc.genes = pancreas_rna@assays$RNA@var.features, npcs = 20, verbose = FALSE)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
Idents(pancreas_rna) <- "ancestry"
pancreas_rna$ancestry_sex <- paste(Idents(pancreas_rna), pancreas_rna$Sex, sep = "_")
table(pancreas_rna@meta.data[["ancestry_sex"]])

#Run Harmony batch correction with library and tissue source covariates
Idents(pancreas_rna) <- "Library"
pancreas_rna <- RunHarmony(pancreas_rna, 
                           assay.use = "RNA",
                           reduction = "pca",
                           dims.use = 1:20,
                           group.by.vars = c('Library','Tissue Source', 'Chemistry'),
                           kmeans_init_nstart=20, kmeans_init_iter_max=100,
                           plot_convergence = TRUE)

unique(pancreas_rna@meta.data[["Library"]])
unique(pancreas_rna@meta.data[["Tissue Source"]])
unique(pancreas_rna@meta.data[["ancestry_sex"]])
unique(pancreas_rna@meta.data[["Cell Type"]])
unique(pancreas_rna@meta.data[["Chemistry"]])
table(pancreas_rna@meta.data[["Chemistry"]])

# Run UMAP
pancreas_rna <- RunUMAP(pancreas_rna, reduction = "harmony", dims = 1:20, return.model = TRUE)
DimPlot(pancreas_rna, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)

# Clustering
pancreas_rna <- pancreas_rna %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 4, 4.6, 5, 6, 7, 8, 9, 10), method = 'igraph') #25 res

# Save file
qsave(pancreas_rna, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\pancreas_rna.qs)")
qsave(pancreas_rna, file = r"(E:\2.SexbasedStudyCurrent\QS files\pancreas_rna.qs)")
})


############################ STAGE ############################
############################   4   ############################

system.time({
# Load data
pancreas_rna <- qread(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\pancreas_rna.qs)")

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
#clustree(pancreas_rna, prefix = "RNA_snn_res.")

# View clustering
DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
FeaturePlot(object = pancreas_rna,
            features = c("MKI67"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 1000,
            slot = 'counts',
            order = TRUE,
            raster=FALSE)

# Subclustering
Idents(pancreas_rna) <- "RNA_snn_res.6"
subset_clust <- subset(pancreas_rna, idents = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", 
                                                "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", 
                                                "41", "42", "43", "44", #"45", 
                                                "46", "47", "48", "49", "50", 
                                                #"51", 
                                                "52", "53", "54", "55", "56", "57", "58", "59", "60",
                                                "61", "62", "63", "64", #"65", 
                                                "66", #"67", 
                                                "68", "69", "70", 
                                                "71", "72", "73", "74", #"75", 
                                                "76", #"77", 
                                                "78", "79", #"80", 
                                                "81", "82", "83", "84", #"85", 
                                                "86", "87", #"88", 
                                                "89", "90" 
                                                ))

# Checking cluster loss
DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
DimPlot(subset_clust, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)

# As cells were subsetted re-run batch correction and cluster assignmnet
subset_clust <- RunHarmony(subset_clust, 
                           assay.use = "RNA",
                           reduction = "pca",
                           dims.use = 1:20,
                           group.by.vars = c('Library','Tissue Source', 'Chemistry'),
                           kmeans_init_nstart=20, kmeans_init_iter_max=100,
                           plot_convergence = TRUE)


# UMAP 
subset_clust <- RunUMAP(subset_clust, reduction = "harmony", dims = 1:20, return.model = TRUE)

#Neighbours + Clustering
subset_clust <- subset_clust %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(0.5, 6), method = 'igraph')

DefaultAssay(subset_clust) <- "RNA"
DimPlot(subset_clust, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
FeaturePlot(object = subset_clust,
            features = c("HSPB1", "DNAJB6", "HSPH1", "GADD45B"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 100,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=TRUE)

FeaturePlot(object = subset_clust,
            features = c("INS"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=TRUE)

DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
FeaturePlot(object = pancreas_rna,
            features = c("TRAC"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=FALSE)

# Cluster assignment
table(subset_clust@meta.data$RNA_snn_res.6)
Idents(subset_clust) <- "RNA_snn_res.6"
subset_clust <- RenameIdents(subset_clust, 
                             "1" = "alpha",
                             "2" = "alpha",
                             "3" = "alpha",
                             "4" = "alpha",
                             "5" = "alpha",
                             "6" = "alpha",
                             "7" = "activated_stellate",
                             "8" = "alpha",
                             "9" = "acinar",
                             "10" = "alpha",
                             "11" = "beta",
                             "12" = "beta",
                             "13" = "ductal",
                             "14" = "beta",
                             "15" = "beta",
                             "16" = "alpha",
                             "17" = "alpha",
                             "18" = "acinar",
                             "19" = "ductal",
                             "20" = "beta",
                             "21" = "beta+alpha",
                             "22" = "activated_stellate",
                             "23" = "alpha",
                             "24" = "acinar",
                             "25" = "ductal",
                             "26" = "beta",
                             "27" = "quiescent_stellate",
                             "28" = "acinar",
                             "29" = "endothelial",
                             "30" = "alpha",
                             "31" = "endothelial",
                             "32" = "beta+delta",
                             "33" = "alpha",
                             "34" = "acinar",
                             "35" = "acinar",
                             "36" = "ductal",
                             "37" = "beta+alpha",
                             "38" = "beta",
                             "39" = "beta+alpha",
                             "40" = "alpha",
                             "41" = "beta",
                             "42" = "acinar",
                             "43" = "beta",
                             "44" = "beta",
                             "45" = "delta",
                             "46" = "beta",
                             "47" = "delta",
                             "48" = "activated_stellate",
                             "49" = "acinar",
                             "50" = "alpha",
                             "51" = "beta",
                             "52" = "alpha",
                             "53" = "gamma",
                             "54" = "mast",
                             "55" = "macrophages",
                             "56" = "beta+alpha",
                             "57" = "beta",
                             "58" = "acinar",
                             "59" = "gamma",
                             "60" = "quiescent_stellate",
                             "61" = "acinar",
                             "62" = "ductal",
                             "63" = "activated_stellate",
                             "64" = "acinar",
                             "65" = "alpha",
                             "66" = "endothelial",
                             "67" = "ductal",
                             "68" = "ductal",
                             "69" = "acinar",
                             "70" = "cycling_endo",
                             "71" = "acinar",
                             "72" = "endothelial",
                             "73" = "acinar",
                             "74" = "acinar",
                             "75" = "endothelial",
                             "76" = "acinar",
                             "77" = "beta",
                             "78" = "alpha",
                             "79" = "alpha",
                             "80" = "alpha",
                             "81" = "schwann",
                             "82" = "activated_stellate",
                             "83" = "beta+alpha",
                             "84" = "beta+alpha",
                             "85" = "ductal",
                             "86" = "acinar",
                             "87" = "ductal"
                             )

# Check renaming
table(subset_clust@active.ident)
unique(subset_clust@active.ident)


############################## #
########## SUBSET1 ########### #
############################## #

# Subset epsilon cells
DefaultAssay(object = subset_clust) <- "RNA"
Idents(subset_clust, WhichCells(object = subset_clust, expression = GHRL > 50, slot = 'counts')) <- 'epsilon'
subset_clust$celltype_qadir <- Idents(subset_clust)
Idents(subset_clust) <- "celltype_qadir"
DimPlot(subset_clust, reduction = "umap", label = TRUE)

############################## #
########## SUBSET2 ########### #
############################## #

# Subsetting Lymphocytes from Mast cells
# First subset all cells, they are called "Mast" in the primary object
mast <- subset(subset_clust, idents = "mast")
table(mast@active.ident)
mast <- mast %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(0.5), method = 'igraph')

# Run UMAP
mast <- RunUMAP(mast, reduction = "harmony", dims = 1:20, return.model = TRUE)
DimPlot(mast, reduction = 'umap', label = FALSE, pt.size = 1, raster=FALSE)
# Cluster assignment
table(mast@meta.data$RNA_snn_res.0.5)
Idents(mast) <- "RNA_snn_res.0.5"
DimPlot(mast, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE, pt.size = 1, raster=FALSE)
FeaturePlot(object = mast,
            features = c("TRAC"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=FALSE)

table(mast@active.ident)
mast <- RenameIdents(mast, 
                     "1" = "mast",
                     "2" = "lymphocyte",
                     "3" = "mast",
                     "4" = "lymphocyte",
                     "5" = "lymphocyte"
)

table(mast@active.ident)

# Generate a new column called celltype_qadir in the metadata copying all Ident info there, this is active idents so check
table(subset_clust@active.ident)
subset_clust$celltype_qadir <- as.character(Idents(subset_clust)) #as.character imp
table(subset_clust$celltype_qadir)

# Change the information of cells containing sub-cluster information
subset_clust$celltype_qadir[Cells(mast)] <- paste(Idents(mast))
table(subset_clust$celltype_qadir)
DimPlot(subset_clust, 
        #split.by = "Tissue Source", 
        group.by = "celltype_qadir", 
        label = FALSE, 
        ncol = 1,  
        cols = c("dodgerblue3",      #beta
                 "turquoise2",       #beta+alpha
                 "lightseagreen",    #alpha
                 "darkseagreen2",    #cycling_endo
                 "khaki2",           #epsilon 
                 "springgreen4",     #gamma
                 "chartreuse3",      #delta
                 "burlywood3",       #beta+delta
                 "darkorange",       #ductal
                 "salmon3",          #acinar
                 "orangered",        #activated_setallate
                 "salmon",           #quiescent_stellate
                 "red",              #endothelial
                 "magenta3",         #macrophages
                 "orchid1",          #lymphocytes
                 "red4",             #mast
                 "grey30"            #schwann
        )
)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "beta+alpha", "alpha", "cycling_endo", "epsilon", "gamma", "delta", "beta+delta",
               "ductal", "acinar", 
               "activated_stellate", "quiescent_stellate", "endothelial",
               "macrophages", "lymphocyte", "mast",
               "schwann")
table(subset_clust@meta.data$celltype_qadir)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
subset_clust@meta.data$celltype_qadir <- factor(x = subset_clust@meta.data$celltype_qadir, levels = my_levels)
table(subset_clust@meta.data$celltype_qadir)
unique(subset_clust@meta.data$celltype_qadir)

# Set celltype_qadir as default
Idents(subset_clust) <- "celltype_qadir"

# Change the information of cells containing sub-cluster information
DimPlot(subset_clust, 
        #split.by = "Diabetes Status", 
        group.by = "celltype_qadir", 
        label = TRUE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.5,
        cols = c("dodgerblue3",      #beta
                 "turquoise2",       #beta+alpha
                 "lightseagreen",    #alpha
                 "darkseagreen2",    #cycling_endo
                 "khaki2",           #epsilon 
                 "springgreen4",     #gamma
                 "chartreuse3",      #delta
                 "burlywood3",       #beta+delta
                 "darkorange",       #ductal
                 "salmon3",          #acinar
                 "orangered",        #activated_setallate
                 "salmon",           #quiescent_stellate
                 "red",              #endothelial
                 "magenta3",         #macrophages
                 "orchid1",          #lymphocytes
                 "red4",             #mast
                 "grey30"            #schwann
        )
)


# Create a metadata slot for celltype_sex, celltype_sex_ancestry and celltype_sex_ancestry_disease
processed_rna < subset_clust
Idents(processed_rna) <- "celltype_qadir"
processed_rna$celltype_sex <- paste(Idents(processed_rna), processed_rna$Sex, sep = "_")
Idents(processed_rna) <- "celltype_sex"
processed_rna$celltype_sex_ancestry <- paste(Idents(processed_rna), processed_rna$ancestry, sep = "_")
Idents(processed_rna) <- "celltype_sex_ancestry"
processed_rna$celltype_sex_ancestry_disease <- paste(Idents(processed_rna), processed_rna$'Diabetes Status', sep = "_")
Idents(processed_rna) <- "celltype_sex"
processed_rna$celltype_sex_disease <- paste(Idents(processed_rna), processed_rna$'Diabetes Status', sep = "_")
table(processed_rna$celltype_qadir)
table(processed_rna$celltype_sex)
table(processed_rna$celltype_sex_ancestry)
table(processed_rna$celltype_sex_ancestry_disease)
table(processed_rna$celltype_sex_disease)

# Save data
#qsave(processed_rna, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\processed_rna.qs)")
#qsave(processed_rna, file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
})



############################ STAGE ############################
############################   5   ############################

###Step 1: Make Pseudobulk Matrices
#Read in final Seurat object
system.time({
adata <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
Idents(adata) <- adata@meta.data$celltype_qadir
samples <- unique(adata@meta.data$Library)

#Pull out list of all cell types
unique_cell_types <- unique(adata$celltype_qadir)

DefaultAssay(adata) <- 'RNA'
#Get counts data
gex.counts <- GetAssayData(adata,slot='counts')

dim(gex.counts)
head(gex.counts)
adata_matrices <- adata

##Pull out barcodes
sample_bcs <- list()
for (sample in samples){
  sample_bcs[[sample]] <- row.names(adata[[]][adata[[]]$Library == sample,])
}

lengths(sample_bcs)
head(sample_bcs[[1]])

#Looping through cell types by making ^ into a function
get_per_sample_gex_SUMS <- function(cell.type, mtx.fp){
  print(paste(cell.type))
  
  #pull out rows of gex.counts where BC Ident matches cell.type
  bcs <- names(Idents(adata_matrices)[Idents(adata_matrices) == cell.type])
  counts <- gex.counts[,colnames(gex.counts) %in% bcs]
  print(dim(counts))
  
  #initialize the matrix of sample gex
  counts.df <- as.data.frame(rep(0,length(row.names(gex.counts))))
  row.names(counts.df) <- row.names(gex.counts)
  colnames(counts.df) <- c('temp')
  
  #go through samples and calculate sum of gex values
  for (sample in samples){
    sample_cols <- colnames(counts) %in% sample_bcs[[sample]]
    counts.cut <- counts[,sample_cols]
    
    #if only one bc, this becomes a vector which is an issue
    if (typeof(counts.cut) == 'double'){
      mean.counts <- counts.cut
      #if there are NO bcs, this will return NA (just return 0 for everything)
    } else if(length(colnames(counts.cut)) == 0){
      mean.counts <- rep(0,length(row.names(counts)))
    } else {
      mean.counts <- rowSums(counts.cut)
    }
    counts.df <- cbind(counts.df,as.data.frame(mean.counts))
  }
  fin.counts.df <- counts.df[,-c(1)]
  colnames(fin.counts.df) <- samples
  head(fin.counts.df)
  
  #export df
  mtx.fp <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/analysis/pseudobulk_counts/%s_sample_gex_total_counts.txt',cell.type) # change to save dir
  write.table(fin.counts.df,mtx.fp,sep='\t',quote=FALSE)
}

#Run function to make matrices
unique_cell_types <- unique(adata$celltype_qadir)
for (cell.type in unique_cell_types){
  fp = sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/analysis/pseudobulk_counts/%s_pseudobulk.txt',cell.type) # change to save dir as above
  get_per_sample_gex_SUMS(cell.type, fp)
}

###Step 2: Make TPM Matrices
#Pull out gene exon info and calculate effective length
gene_annotations_gtf_fp <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.annotation.gtf'
suppressMessages(txdb <- makeTxDbFromGFF(gene_annotations_gtf_fp,format="gtf"))
exons.list.per.gene <- exonsBy(txdb,by="gene") #Collect the exons per gene_id
#Reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum them
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

#checking <- gene.info
gene.info <- rtracklayer::import('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.basic.annotation.gtf')
gene.info <- as.data.frame(gene.info)
gene.info <- gene.info %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
head(gene.info)

# Calculate entire gene boundaries
gene.info.start <- gene.info %>% group_by(gene_id) %>% slice_min(order_by = start)
gene.info.start <- gene.info.start[!duplicated(gene.info.start$gene_id),]
gene.info.start <- gene.info.start %>% select(gene_id, gene_name, gene_type, seqnames, end, strand, source, level)

gene.info.end <- gene.info %>% group_by(gene_id) %>% slice_max(order_by = end)
gene.info.end <- gene.info.end[!duplicated(gene.info.end$gene_id),]
gene.info.end <- gene.info.end %>% select(gene_id, start)

gene.info.comp <- merge(gene.info.end, gene.info.start, by = "gene_id")
gene.info.comp <- gene.info.comp %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
gene.info.comp$check <- ifelse(gene.info.comp$end > gene.info.comp$start, 'TRUE',
                               ifelse(gene.info.comp$end < gene.info.comp$start, 'FALSE'))

unique(gene.info.comp$check)
gene.info <- gene.info.comp

#Add the effective lengths to the original gene.info dataframe
temp_df <- gene.info
rownames(temp_df) <- gene.info$gene_id
temp_df2 <- as.data.frame(exonic.gene.sizes)
temp_df2$gene_id <- rownames(temp_df2)

new_df <- merge(temp_df,temp_df2, by='gene_id', all=TRUE)
#Remove duplicate rows from gene info df
fin.gene.info <- new_df[!duplicated(new_df$gene_name),]

#Read in psedobulk matrices from above 
dir = 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/analysis/pseudobulk_counts/'

files = list.files(dir, pattern =".txt")
cells = gsub("_sample_gex_total_counts.txt","", files)

make_tpm = function(raw_counts, gene_sizes){
  rpk <- raw_counts / gene_sizes
  tpm <- rpk
  for (i in 1:ncol(rpk)){
    tpm[,i] <- rpk[,i]/(sum(rpk[,i])/1e6)
    
  }
  return(tpm)
}

#Output dir for TPM matrices
outdir = "C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/analysis/TPM/"


for (FILE in files){
  cell <- cells[which(files == FILE)]
  raw_counts <- read.table(paste0(dir, FILE), row.names=1)
  raw_counts <- subset(raw_counts ,rownames(raw_counts) %in% fin.gene.info$gene_name)
  gene_sizes <- fin.gene.info$exonic.gene.sizes[match(rownames(raw_counts), fin.gene.info$gene_name )]
  
  tpm_mat <- make_tpm(raw_counts, gene_sizes)
  write.table(tpm_mat, paste0(outdir,  cell, "_TPM_per_sample.txt"), sep="\t", quote=F)
}

###Step 3: DESeq
#Create a metadata table
meta <- adata@meta.data[,c('Library', 'Sex', 'Tissue Source', 'Chemistry', 'ancestry', 'Diabetes Status')]
colnames(meta) <- c('Library', 'Sex', 'Tissue_Source', 'Chemistry', 'ancestry', 'Diabetes_Status')
rownames(meta) <- NULL
meta <- meta[!duplicated(meta),]
meta$sex_ancestry_diabetes <- paste0(meta$Sex, '_', meta$ancestry, '_', meta$Diabetes_Status)
meta$sex_diabetes <- paste0(meta$Sex, '_', meta$Diabetes_Status)

# #Create all combinations of tests comparing 2 conditions using the meta$sex_ancestry_diabetes column we created
# combinations <- combn(meta$sex_ancestry_diabetes, 2)
# combinations <- t(combinations)
# combinations <- as.data.frame(combinations)
# combinations <- combinations[which(combinations$V1 != combinations$V2),]
# split1 <- str_split_fixed(combinations$V1, pattern='_',n=3)
# split1 <- as.data.frame(split1)
# split2 <- str_split_fixed(combinations$V2, pattern='_',n=3)
# split2 <- as.data.frame(split2)
# combinations <- cbind(combinations, split1, split2)
# colnames(combinations) <- c('test1', 'test2', 'sex1', 'ancestry1', 'diabetes1','sex2', 'ancestry2', 'diabetes2')
# 
# #These next 4 lines will create the combinations where 2 variables remain the same and only one changes
# ## ie. would remove M_black_T2D vs. F_white_T2D since there is only a single common variable
# keep1 <- combinations[which(combinations$sex1 == combinations$sex2 & combinations$ancestry1 == combinations$ancestry2),]
# keep2 <- combinations[which(combinations$sex1 == combinations$sex2 & combinations$diabetes1 == combinations$diabetes2),]
# keep3 <- combinations[which(combinations$ancestry1 == combinations$ancestry2 & combinations$diabetes1 == combinations$diabetes2),]
# keep <- rbind(keep1,keep2,keep3)
# head(keep)

#Create all combinations of tests comparing 2 conditions using the meta$sex_ancestry_diabetes column we created
combinations <- as.data.frame(t(combn(meta$sex_ancestry_diabetes, 2)))
combinations <- combinations[which(combinations$V1 != combinations$V2),]
split1 <- as.data.frame(str_split_fixed(combinations$V1, pattern='_',n=3))
split2 <- as.data.frame(str_split_fixed(combinations$V2, pattern='_',n=3))
combinations <- cbind(combinations, split1, split2)
colnames(combinations) <- c('test1', 'test2', 'sex1', 'ancestry1', 'diabetes1','sex2', 'ancestry2', 'diabetes2')

#These next 4 lines will create the combinations where 2 variables remain the same and only one changes
## ie. would remove M_black_T2D vs. F_white_T2D since there is only a single common variable
keep1 <- combinations[which(combinations$sex1 == combinations$sex2 & combinations$ancestry1 == combinations$ancestry2),]
keep2 <- combinations[which(combinations$sex1 == combinations$sex2 & combinations$diabetes1 == combinations$diabetes2),]
keep3 <- combinations[which(combinations$ancestry1 == combinations$ancestry2 & combinations$diabetes1 == combinations$diabetes2),]
keep <- rbind(keep1,keep2,keep3)

keep$check <- paste0(keep$test1, '_', keep$test2) #should be testing combinations
keep <- keep[!duplicated(keep$check),] #remove duplicate testing combinations

#Pseudobulk matrices directory
dir <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/analysis/pseudobulk_counts/'
#Create outdir for results
outdir <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/analysis/DE_testing/'
#dir.create(outdir) #works like mkdir
files <- list.files('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/analysis/pseudobulk_counts/', pattern='gex')

##Create matrices for results
sumres <- matrix(nrow=length(cells), ncol = 3)
rownames(sumres) <- cells

# testing for sex_ancestry_diabetes
for (FILE in files) {
  cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
  raw_counts <- read.table(paste0(dir, FILE), row.names=1)
  sample_names <- unique(adata@meta.data$Library)
  sample_names <- gsub('-','.', sample_names)
  raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
  raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
  meta$Library2 <- gsub('-', '.', meta$Library)
  meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
  
  if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
    print(cell)
    print('Data for 2 sex present, however not all data may be present will check this at a later step')
    
    genes_to_keep <- c()
    for (i in 1:nrow(raw_counts)) {
      if (sum(raw_counts[i, ] >= 5) >= 2) {
        genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
      }
    }
    counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
    #counts <- raw_counts[rowSums(raw_counts) >= 10,] #Light pre-filtering
    
    if (length(unique(meta2$Chemistry)) > 1) {
      my_design <- as.formula ('~sex_ancestry_diabetes + Chemistry + Tissue_Source')
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    } else {
      my_design <- as.formula ('~sex_ancestry_diabetes + Tissue_Source')
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    }
    
    tests1 <- c('M_white_ND', 'M_white_ND', 'M_black_ND', 'F_white_ND', 'F_white_ND', 'F_black_ND', 'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'M_white_T2D', 
                'M_white_T2D', 'M_black_T2D', 'F_white_T2D', 'F_white_T2D', 'F_black_T2D', 'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 
                'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D')
    
    tests2 <- c('M_hispanic_ND', 'M_black_ND', 'M_hispanic_ND', 'F_hispanic_ND', 'F_black_ND', 'F_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND', 'M_hispanic_T2D', 
                'M_black_T2D', 'M_hispanic_T2D', 'F_hispanic_T2D', 'F_black_T2D', 'F_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D', 
                'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND')
    
    print('Preparing to run DESeq2')
    
    for (x in 1:length(tests1)){
      t1 <- tests1[[x]]
      t2 <- tests2[[x]]
      test <- c('sex_ancestry_diabetes', tests1[[x]],tests2[[x]])
      numoft1 <- length(which(meta2$sex_ancestry_diabetes==t1))
      numoft2 <- length(which(meta2$sex_ancestry_diabetes==t2))
      
      # if (numoft1 > 2 & numoft2 > 2) {
      #   print(paste(t1, 'and', t2, 'are present in the dataset', sep = ' '))
      #   print(paste("Data copied here:", outdir, sep = " "))
      # }
      
      if (numoft1 < 3) {
        message(paste("!!WARNING!!"))
        message(paste(t1, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste('####'))
        message(paste('####'))
      } else if (numoft2 < 3) {
        message(paste("!!WARNING!!"))
        message(paste(t2, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste("####"))
        message(paste("####"))
      } else if (numoft1 > 2 & numoft2 > 2) {
        #sprintf("%s and %s are present in the dataset", t1, t2)
        #sprintf("Find data here: %s", outdir)
        res <- results(dds, contrast=c(test))
        res <- as.data.frame(res)
        res <- res[order(res$pvalue),]
        outfile <- paste0(cell, '.deseq.WaldTest.', tests1[[x]],'.vs.',tests2[[x]],'.tsv')
        write.table(res,paste0(outdir, outfile) , sep='\t', quote=F)
        #print(paste(t1, 'and', t2, 'are present in dataset metadata', sep = " "))
        print(sprintf('%s and %s are present in the dataset metadata', t1, t2)) #just because I wanted to understand using sprintf
        print(paste("Data copied here:", outdir, sep = " "))
        print(paste('####'))
        print(paste('####'))
      }
      
    }
  }
}
        
# testing for sex_diabetes
for (FILE in files) {
  cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
  raw_counts <- read.table(paste0(dir, FILE), row.names=1)
  sample_names <- unique(adata@meta.data$Library)
  sample_names <- gsub('-','.', sample_names)
  raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
  raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
  meta$Library2 <- gsub('-', '.', meta$Library)
  meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
  
  if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
    print(cell)
    print('Data for 2 sex present, however not all data may be present will check this at a later step')
    
    genes_to_keep <- c()
    for (i in 1:nrow(raw_counts)) {
      if (sum(raw_counts[i, ] >= 5) >= 2) {
        genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
      }
    }
    counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
    #counts <- raw_counts[rowSums(raw_counts) >= 10,] #Light pre-filtering
    
    if (length(which(meta2$sex_diabetes == 'M_ND')) > 1 && 
        length(which(meta2$sex_diabetes == 'M_T2D')) > 1 && 
        length(which(meta2$sex_diabetes == 'F_ND')) > 1 && 
        length(which(meta2$sex_diabetes == 'F_T2D')) > 1) {
      my_design <- as.formula ('~sex_diabetes + Chemistry + Tissue_Source') # design for sex_diabetes
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    } else {
      print(sprintf('%s does not have sufficient diabetes samples to test, bypassing to test ND only', cell))
      meta2 <- subset(meta2, Diabetes_Status == 'ND') # it is possible some T2D are present so eliminate them from your dataset since you are restricted to sex
      counts <- counts[,meta2$Library2]
      my_design <- as.formula ('~Sex +  Tissue_Source') # design for sex_diabetes
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    }
    
    if (length(which(meta2$sex_diabetes == 'M_ND')) > 1 && 
        length(which(meta2$sex_diabetes == 'M_T2D')) > 1 && 
        length(which(meta2$sex_diabetes == 'F_ND')) > 1 && 
        length(which(meta2$sex_diabetes == 'F_T2D')) > 1) {
      # This is now the next test. Your samples need to be > 1  
      tests1 <- c('M_ND', 'M_ND', 'F_ND', 'M_T2D')
      tests2 <- c('F_ND', 'M_T2D', 'F_T2D', 'F_T2D')
    } else {
      tests1 <- c("M")
      tests2 <- c("F")
    }
    
    print('Preparing to run DESeq2')
    
    for (x in 1:length(tests1)){
      t1 <- tests1[[x]]
      t2 <- tests2[[x]]
      if (length(which(meta2$sex_diabetes == 'M_ND')) > 1 && 
          length(which(meta2$sex_diabetes == 'M_T2D')) > 1 && 
          length(which(meta2$sex_diabetes == 'F_ND')) > 1 && 
          length(which(meta2$sex_diabetes == 'F_T2D')) > 1) {
        test <- c('sex_diabetes', tests1[[x]],tests2[[x]]) # For Sex_diabetes
        numoft1 <- length(which(meta2$sex_diabetes==t1))
        numoft2 <- length(which(meta2$sex_diabetes==t2))
      } else {
        test <- c('Sex', tests1[[x]],tests2[[x]]) # For sex only (for example Schwann cells)
        numoft1 <- length(which(meta2$Sex==t1)) # For sex diabetes
        numoft2 <- length(which(meta2$Sex==t2))
      }
      
      # if (numoft1 > 2 & numoft2 > 2) {
      #   print(paste(t1, 'and', t2, 'are present in the dataset', sep = ' '))
      #   print(paste("Data copied here:", outdir, sep = " "))
      # }
      
      if (numoft1 < 3) {
        message(paste("!!WARNING CHECK METADATA!!"))
        message(paste(t1, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste('####'))
        message(paste('####'))
      } else if (numoft2 < 3) {
        message(paste("!!WARNING CHECK METADATA!!"))
        message(paste(t2, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste("####"))
        message(paste("####"))
      } else if (numoft1 > 2 & numoft2 > 2) {
        #sprintf("%s and %s are present in the dataset", t1, t2)
        #sprintf("Find data here: %s", outdir)
        res <- results(dds, contrast=c(test))
        res <- as.data.frame(res)
        res <- res[order(res$pvalue),]
        outfile <- paste0(cell, '.deseq.WaldTest.', tests1[[x]],'.vs.',tests2[[x]],'.tsv')
        write.table(res,paste0(outdir, outfile) , sep='\t', quote=F)
        #print(paste(t1, 'and', t2, 'are present in dataset metadata', sep = " "))
        print(sprintf('%s cells and %s cells are present in the dataset metadata', t1, t2)) #just because I wanted to understand using sprintf
        print(paste("Data copied here:", outdir, sep = " "))
        print(paste('####'))
        print(paste('####'))
      }
      
    }
  }
}
}) # System time

# #Create pseudobulk matrix from all cell types
# gene_x_cell <- adata@assays[['SCT']]@counts #gene by cell matrix
# 
# bulk <- data.frame(matrix(nrow = 26275, ncol = 0)) 
# for (x in samples){
#   sample_subset <- gene_x_cell[,grep(x, colnames(gene_x_cell))]
#   sample_means <- Matrix::rowSums(sample_subset)
#   rownames(bulk) <- rownames(sample_subset)
#   bulk[[x]] <- sample_means
#   
# }
# 
# dds_bulk <- DESeqDataSetFromMatrix(
#   countData = round(bulk),
#   meta,
#   design= ~sex_ancestry_diabetes + Chemistry + Tissue_Source)
# 
# dds_bulk <- estimateSizeFactors(dds_bulk)
# dds_bulk <- estimateDispersions(dds_bulk)    
# vsd_bulk <- varianceStabilizingTransformation(dds_bulk)
# 
# #Make a PCA
# options(repr.plot.height = 7, repr.plot.width = 14)
# pcaData <- plotPCA(vsd_bulk, intgroup=c('ancestry', 'Sex', 'Diabetes_Status', 'Tissue_Source'), returnData=TRUE, ntop=5000)
# percentVar <- round(100 * attr(pcaData, 'percentVar'))
# ggplot(pcaData, aes(PC1, PC2, color=ancestry, shape=Sex)) +
#   geom_point(size=3) +
#   xlab(paste0('PC1: ',percentVar[1],'% variance')) +
#   ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
#   coord_fixed() + theme_minimal()
# 
# #Make a heatmap
# options(repr.plot.height = 15, repr.plot.width = 20)
# sampleDists <- dist(t(assay(vsd_bulk)))
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- vsd_bulk$Library
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)

############################ STAGE ############################
############################   6   ############################
qsave(processed_rna, file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
DimPlot(subset_clust, 
        #split.by = "Diabetes Status", 
        group.by = "celltype_qadir", 
        label = TRUE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.5,
        cols = c("dodgerblue3",      #beta
                 "turquoise2",       #beta+alpha
                 "lightseagreen",    #alpha
                 "darkseagreen2",    #cycling_endo
                 "khaki2",           #epsilon 
                 "springgreen4",     #gamma
                 "chartreuse3",      #delta
                 "burlywood3",       #beta+delta
                 "darkorange",       #ductal
                 "salmon3",          #acinar
                 "orangered",        #activated_setallate
                 "salmon",           #quiescent_stellate
                 "red",              #endothelial
                 "magenta3",         #macrophages
                 "orchid1",          #lymphocytes
                 "red4",             #mast
                 "grey30"            #schwann
        )
)

# Gene ontology analysis Rapid Gene ontology Auto Loader (Rapid GOAL)
# Create a list of all files in directory
dgelist <- list.files(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\DE Testing\Windows\data_hpap)", 
                      all.files = FALSE, 
                      full.names = FALSE, 
                      pattern = "*.csv")

# Point towards WD using a function
for (sample in dgelist){
  wd <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/DE Testing/Windows/data_hpap/%s', dgelist)
}

# Run iterative function to perform GO on all data
for (x in wd) {
  sample_name <- str_split_fixed(x, "/", n=12)[12]
  datfile <- read.csv(file.path(x), row.names = 1)
  
  # Gene list of genes going UP
  sig_df_up <- dplyr::filter(datfile, p_val < 0.05 & avg_log2FC > 0.26303) # >1.2x
  sig_genes_up <- rownames(sig_df_up)
  
  # Gene list of genes going UP
  sig_df_down <- dplyr::filter(datfile, p_val < 0.05 & avg_log2FC < -0.32192) # <0.8x
  sig_genes_down <- rownames(sig_df_down)
  
  # All genes
  all_genes <- rownames(datfile)
  
  # Run GO enrichment analysis genes up
  GO.up <- enrichGO(gene = sig_genes_up, 
                    universe = all_genes, 
                    keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                    OrgDb = org.Hs.eg.db, 
                    ont = c("ALL"), 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 1, 
                    qvalueCutoff = 1, #if not set default is at 0.05
                    readable = TRUE)
  
  # Run GO enrichment analysis genes down
  GO.down <- enrichGO(gene = sig_genes_down, 
                      universe = all_genes, 
                      keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                      OrgDb = org.Hs.eg.db, 
                      ont = c("ALL"), 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1, 
                      qvalueCutoff = 1, #if not set default is at 0.05
                      readable = TRUE)
  
  go_data_up <- data.frame(GO.up)
  go_data_down <- data.frame(GO.down)
  
  go_data_up <- dplyr::filter(go_data_up, pvalue < 0.05)
  go_data_down <- dplyr::filter(go_data_down, pvalue < 0.05)
  
  # Save outputs
  write.csv(go_data_up, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/!FAHD/4. Sex and Race Based Study Project/Sequencing_Data/scRNAseq/hpap_combined/ORA/UP/%s", sample_name), row.names = FALSE)
  write.csv(go_data_down, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/!FAHD/4. Sex and Race Based Study Project/Sequencing_Data/scRNAseq/hpap_combined/ORA/DOWN/%s", sample_name), row.names = FALSE)
}



############################ END ############################
############################ END ############################

# 
#
#
#
#
#
#
#
#
#
#
#
#
# View clustering
DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)
pancreas_rna <- FindClusters(object = pancreas_rna, algorithm=4, resolution = 0, method = 'igraph')


Idents(pancreas_rna) <- "Chemistry"

# QC for Ruth split on chemistry
v2 <- subset(pancreas_rna, idents = c("10Xv2"))
v3 <- subset(pancreas_rna, idents = c("10Xv3"))

Idents(pancreas_rna) <- "Chemistry"
DimPlot(pancreas_rna, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)

p1 <- DimPlot(v2, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)
p2 <- DimPlot(v3, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)
p1|p2

RunUMAP(pancreas_rna, reduction = "harmony", dims = 1:20, return.model = TRUE, reduction.name = 'umap')
DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)
DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)














