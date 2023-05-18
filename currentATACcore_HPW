# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
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
BiocManager::install("scDblFinder")
BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")
BiocManager::install("motifmatchr")
BiocManager::install("GreenleafLab/chromVAR")

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

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')

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
    library(scDblFinder)
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
    library(JASPAR2020)
    library(TFBSTools)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(motifmatchr)
    library(chromVAR)
    library(SeuratWrappers)
    library(cicero)
  }
)


# Set global environment parameter par-proc
#options(future.globals.maxSize = 8000 * 1024^2)
set.seed(1234)

# Python env
if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/mqadir/AppData/Local/r-miniconda/envs/r-reticulate",Sys.getenv()["PATH"],sep=";"))
py_config()

# WD
setwd(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\WD)")
(WD <- getwd())


# Check package versions
packageVersion("clusterProfiler")
packageVersion("dittoSeq")
packageVersion("escape")
packageVersion("seurat")
packageVersion("signac")
packageVersion("EnrichmentBrowser")

# Calculation of Doublets https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scATAC.html
# Provide GRanges of repeat elements for exclusion:
suppressPackageStartupMessages(library(GenomicRanges))
repeats <- GRanges("chr6", IRanges(1000,2000))

# Combine with mitochondrial and sex chromosomes
otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))

# Combining them:
toExclude <- suppressWarnings(c(repeats, otherChroms))

# Running amulet method
{
  fragfile.HP2022801 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.SAMN15877725 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2024001 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2031401 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2105501 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2106201 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2107001 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2107901 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2108601 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2108901 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2110001 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2121601 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2123201 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2132801 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2202101 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\fragments.tsv.gz)", regionsToExclude=toExclude)
}

# Adding sampleID to barcode nomenclature
head(combined_atac) # shows that config of cell nomenclature is SAMPLEID_barcode [yeah I know, this shouldnt but hey live dangerous.]
{
  rownames(fragfile.HP2022801) <- paste0('HP2022801_', rownames(fragfile.HP2022801))
  rownames(fragfile.SAMN15877725) <- paste0('SAMN15877725_', rownames(fragfile.SAMN15877725))
  rownames(fragfile.HP2024001) <- paste0('HP2024001_', rownames(fragfile.HP2024001))
  rownames(fragfile.HP2031401) <- paste0('HP2031401_', rownames(fragfile.HP2031401))
  rownames(fragfile.HP2105501) <- paste0('HP2105501_', rownames(fragfile.HP2105501))
  rownames(fragfile.HP2106201) <- paste0('HP2106201_', rownames(fragfile.HP2106201))
  rownames(fragfile.HP2107001) <- paste0('HP2107001_', rownames(fragfile.HP2107001))
  rownames(fragfile.HP2107901) <- paste0('HP2107901_', rownames(fragfile.HP2107901))
  rownames(fragfile.HP2108601) <- paste0('HP2108601_', rownames(fragfile.HP2108601))
  rownames(fragfile.HP2108901) <- paste0('HP2108901_', rownames(fragfile.HP2108901))
  rownames(fragfile.HP2110001) <- paste0('HP2110001_', rownames(fragfile.HP2110001))
  rownames(fragfile.HP2121601) <- paste0('HP2121601_', rownames(fragfile.HP2121601))
  rownames(fragfile.HP2123201) <- paste0('HP2123201_', rownames(fragfile.HP2123201))
  rownames(fragfile.HP2132801) <- paste0('HP2132801_', rownames(fragfile.HP2132801))
  rownames(fragfile.HP2202101) <- paste0('HP2202101_', rownames(fragfile.HP2202101))
}

combined_atac_doublet <- do.call('rbind', list(fragfile.HP2022801, fragfile.SAMN15877725, fragfile.HP2024001, fragfile.HP2031401,
                                               fragfile.HP2105501, fragfile.HP2106201, fragfile.HP2107001, fragfile.HP2107901,
                                               fragfile.HP2108601, fragfile.HP2108901, fragfile.HP2110001, fragfile.HP2121601,
                                               fragfile.HP2123201, fragfile.HP2132801, fragfile.HP2202101))

#Save file
# saveRDS(fragfile.HP2022801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2022801.rds)")
# saveRDS(fragfile.SAMN15877725, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.SAMN15877725.rds)")
# saveRDS(fragfile.HP2024001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2024001.rds)")
# saveRDS(fragfile.HP2031401, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2031401.rds)")
# saveRDS(fragfile.HP2105501, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2105501.rds)")
# saveRDS(fragfile.HP2106201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2106201.rds)")
# saveRDS(fragfile.HP2107001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2107001.rds)")
# saveRDS(fragfile.HP2107901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2107901.rds)")
# saveRDS(fragfile.HP2108601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2108601.rds)")
# saveRDS(fragfile.HP2108901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2108901.rds)")
# saveRDS(fragfile.HP2110001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2110001.rds)")
# saveRDS(fragfile.HP2121601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2121601.rds)")
# saveRDS(fragfile.HP2123201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2123201.rds)")
# saveRDS(fragfile.HP2132801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2132801.rds)")
# saveRDS(fragfile.HP2202101, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2202101.rds)")
# saveRDS(combined_atac_doublet, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\combined_atac_doublet.rds)")

# Load
# fragfile.HP2022801 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2022801.rds)")
# fragfile.SAMN15877725 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.SAMN15877725.rds)")
# fragfile.HP2024001 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2024001.rds)")
# fragfile.HP2031401 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2031401.rds)")
# fragfile.HP2105501 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2105501.rds)")
# fragfile.HP2106201 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2106201.rds)")
# fragfile.HP2107001 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2107001.rds)")
# fragfile.HP2107901 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2107901.rds)")
# fragfile.HP2108601 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2108601.rds)")
# fragfile.HP2108901 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2108901.rds)")
# fragfile.HP2110001 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2110001.rds)")
# fragfile.HP2121601 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2121601.rds)")
# fragfile.HP2123201 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2123201.rds)")
# fragfile.HP2132801 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2132801.rds)")
# fragfile.HP2202101 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2202101.rds)")
combined_atac_doublet <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac_doublet.rds)")

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data #### https://stuartlab.org/signac/articles/merging.html
{
  # read in peak sets
  peaks.HP2022801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.SAMN15877725 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2024001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2031401 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2105501 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2106201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2107001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2107901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2108601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2108901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2110001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2121601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2123201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2132801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2202101 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
}

# Conversion of peaks to genomic Ranges
{
  gr.HP2022801 <- makeGRangesFromDataFrame(peaks.HP2022801)
  gr.SAMN15877725 <- makeGRangesFromDataFrame(peaks.SAMN15877725)
  gr.HP2024001 <- makeGRangesFromDataFrame(peaks.HP2024001)
  gr.HP2031401 <- makeGRangesFromDataFrame(peaks.HP2031401)
  gr.HP2105501 <- makeGRangesFromDataFrame(peaks.HP2105501)
  gr.HP2106201 <- makeGRangesFromDataFrame(peaks.HP2106201)
  gr.HP2107001 <- makeGRangesFromDataFrame(peaks.HP2107001)
  gr.HP2107901 <- makeGRangesFromDataFrame(peaks.HP2107901)
  gr.HP2108601 <- makeGRangesFromDataFrame(peaks.HP2108601)
  gr.HP2108901 <- makeGRangesFromDataFrame(peaks.HP2108901)
  gr.HP2110001 <- makeGRangesFromDataFrame(peaks.HP2110001)
  gr.HP2121601 <- makeGRangesFromDataFrame(peaks.HP2121601)
  gr.HP2123201 <- makeGRangesFromDataFrame(peaks.HP2123201)
  gr.HP2132801 <- makeGRangesFromDataFrame(peaks.HP2132801)
  gr.HP2202101 <- makeGRangesFromDataFrame(peaks.HP2202101)
  
  # Create a unified set of peaks to quantify in each dataset
  combined.peaks <- reduce(x = c(gr.HP2022801, gr.SAMN15877725, gr.HP2024001, gr.HP2031401,
                                 gr.HP2105501, gr.HP2106201, gr.HP2107001, gr.HP2107901,
                                 gr.HP2108601, gr.HP2108901, gr.HP2110001, gr.HP2121601,
                                 gr.HP2123201, gr.HP2132801, gr.HP2202101))
  
  # Filter out bad peaks based on length
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  combined.peaks
}

# Create Fragment objects
# Load metadata
{
  md.HP2022801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.SAMN15877725 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2024001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2031401 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2105501 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2106201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2107001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2107901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2108601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2108901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2110001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2121601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2123201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2132801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2202101 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
}

# perform an initial filtering of low count cells de-trash the data
md.HP2022801 <- md.HP2022801[md.HP2022801$passed_filters > 500, ]
md.SAMN15877725 <- md.SAMN15877725[md.SAMN15877725$passed_filters > 500, ]
md.HP2024001 <- md.HP2024001[md.HP2024001$passed_filters > 500, ]
md.HP2031401 <- md.HP2031401[md.HP2031401$passed_filters > 500, ]
md.HP2105501 <- md.HP2105501[md.HP2105501$passed_filters > 500, ]
md.HP2106201 <- md.HP2106201[md.HP2106201$passed_filters > 500, ]
md.HP2107001 <- md.HP2107001[md.HP2107001$passed_filters > 500, ]
md.HP2107901 <- md.HP2107901[md.HP2107901$passed_filters > 500, ]
md.HP2108601 <- md.HP2108601[md.HP2108601$passed_filters > 500, ]
md.HP2108901 <- md.HP2108901[md.HP2108901$passed_filters > 500, ]
md.HP2110001 <- md.HP2110001[md.HP2110001$passed_filters > 500, ]
md.HP2121601 <- md.HP2121601[md.HP2121601$passed_filters > 500, ]
md.HP2123201 <- md.HP2123201[md.HP2123201$passed_filters > 500, ]
md.HP2132801 <- md.HP2132801[md.HP2132801$passed_filters > 500, ]
md.HP2202101 <- md.HP2202101[md.HP2202101$passed_filters > 500, ]

# create fragment objects
{
  frags.HP2022801 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\fragments.tsv.gz)",
    cells = rownames(md.HP2022801)
  )
  
  frags.SAMN15877725 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\fragments.tsv.gz)",
    cells = rownames(md.SAMN15877725)
  )
  
  frags.HP2024001 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\fragments.tsv.gz)",
    cells = rownames(md.HP2024001)
  )
  
  frags.HP2031401 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\fragments.tsv.gz)",
    cells = rownames(md.HP2031401)
  )
  
  frags.HP2105501 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\fragments.tsv.gz)",
    cells = rownames(md.HP2105501)
  )
  
  frags.HP2106201 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\fragments.tsv.gz)",
    cells = rownames(md.HP2106201)
  )
  
  frags.HP2107001 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\fragments.tsv.gz)",
    cells = rownames(md.HP2107001)
  )
  
  frags.HP2107901 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\fragments.tsv.gz)",
    cells = rownames(md.HP2107901)
  )
  
  frags.HP2108601 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\fragments.tsv.gz)",
    cells = rownames(md.HP2108601)
  )
  
  frags.HP2108901 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\fragments.tsv.gz)",
    cells = rownames(md.HP2108901)
  )
  
  frags.HP2110001 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\fragments.tsv.gz)",
    cells = rownames(md.HP2110001)
  )
  
  frags.HP2121601 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\fragments.tsv.gz)",
    cells = rownames(md.HP2121601)
  )
  
  frags.HP2123201 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\fragments.tsv.gz)",
    cells = rownames(md.HP2123201)
  )
  
  frags.HP2132801 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\fragments.tsv.gz)",
    cells = rownames(md.HP2132801)
  )
  
  frags.HP2202101 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\fragments.tsv.gz)",
    cells = rownames(md.HP2202101)
  )
}

# Quantify Peaks
{
  HP2022801.counts <- FeatureMatrix(
    fragments = frags.HP2022801,
    features = combined.peaks,
    cells = rownames(md.HP2022801)
  )
  
  SAMN15877725.counts <- FeatureMatrix(
    fragments = frags.SAMN15877725,
    features = combined.peaks,
    cells = rownames(md.SAMN15877725)
  )
  
  HP2024001.counts <- FeatureMatrix(
    fragments = frags.HP2024001,
    features = combined.peaks,
    cells = rownames(md.HP2024001)
  )
  
  HP2031401.counts <- FeatureMatrix(
    fragments = frags.HP2031401,
    features = combined.peaks,
    cells = rownames(md.HP2031401)
  )
  
  HP2105501.counts <- FeatureMatrix(
    fragments = frags.HP2105501,
    features = combined.peaks,
    cells = rownames(md.HP2105501)
  )
  
  HP2106201.counts <- FeatureMatrix(
    fragments = frags.HP2106201,
    features = combined.peaks,
    cells = rownames(md.HP2106201)
  )
  
  HP2107001.counts <- FeatureMatrix(
    fragments = frags.HP2107001,
    features = combined.peaks,
    cells = rownames(md.HP2107001)
  )
  
  HP2107901.counts <- FeatureMatrix(
    fragments = frags.HP2107901,
    features = combined.peaks,
    cells = rownames(md.HP2107901)
  )
  
  HP2108601.counts <- FeatureMatrix(
    fragments = frags.HP2108601,
    features = combined.peaks,
    cells = rownames(md.HP2108601)
  )
  
  HP2108901.counts <- FeatureMatrix(
    fragments = frags.HP2108901,
    features = combined.peaks,
    cells = rownames(md.HP2108901)
  )
  
  HP2110001.counts <- FeatureMatrix(
    fragments = frags.HP2110001,
    features = combined.peaks,
    cells = rownames(md.HP2110001)
  )
  
  HP2121601.counts <- FeatureMatrix(
    fragments = frags.HP2121601,
    features = combined.peaks,
    cells = rownames(md.HP2121601)
  )
  
  HP2123201.counts <- FeatureMatrix(
    fragments = frags.HP2123201,
    features = combined.peaks,
    cells = rownames(md.HP2123201)
  )
  
  HP2132801.counts <- FeatureMatrix(
    fragments = frags.HP2132801,
    features = combined.peaks,
    cells = rownames(md.HP2132801)
  )
  
  HP2202101.counts <- FeatureMatrix(
    fragments = frags.HP2202101,
    features = combined.peaks,
    cells = rownames(md.HP2202101)
  )
}

# STEP 2: Create Seurat objects ####
{
  HP2022801_assay <- CreateChromatinAssay(HP2022801.counts, fragments = frags.HP2022801, min.features = 100)
  HP2022801_atac <- CreateSeuratObject(HP2022801_assay, assay = "ATAC", meta.data=md.HP2022801)
  
  SAMN15877725_assay <- CreateChromatinAssay(SAMN15877725.counts, fragments = frags.SAMN15877725, min.features = 100)
  SAMN15877725_atac <- CreateSeuratObject(SAMN15877725_assay, assay = "ATAC", meta.data=md.SAMN15877725)
  
  HP2024001_assay <- CreateChromatinAssay(HP2024001.counts, fragments = frags.HP2024001, min.features = 100)
  HP2024001_atac <- CreateSeuratObject(HP2024001_assay, assay = "ATAC", meta.data=md.HP2024001)
  
  HP2031401_assay <- CreateChromatinAssay(HP2031401.counts, fragments = frags.HP2031401, min.features = 100)
  HP2031401_atac <- CreateSeuratObject(HP2031401_assay, assay = "ATAC", meta.data=md.HP2031401)
  
  HP2105501_assay <- CreateChromatinAssay(HP2105501.counts, fragments = frags.HP2105501, min.features = 100)
  HP2105501_atac <- CreateSeuratObject(HP2105501_assay, assay = "ATAC", meta.data=md.HP2105501)
  
  HP2106201_assay <- CreateChromatinAssay(HP2106201.counts, fragments = frags.HP2106201, min.features = 100)
  HP2106201_atac <- CreateSeuratObject(HP2106201_assay, assay = "ATAC", meta.data=md.HP2106201)
  
  HP2107001_assay <- CreateChromatinAssay(HP2107001.counts, fragments = frags.HP2107001, min.features = 100)
  HP2107001_atac <- CreateSeuratObject(HP2107001_assay, assay = "ATAC", meta.data=md.HP2107001)
  
  HP2107901_assay <- CreateChromatinAssay(HP2107901.counts, fragments = frags.HP2107901, min.features = 100)
  HP2107901_atac <- CreateSeuratObject(HP2107901_assay, assay = "ATAC", meta.data=md.HP2107901)
  
  HP2108601_assay <- CreateChromatinAssay(HP2108601.counts, fragments = frags.HP2108601, min.features = 100)
  HP2108601_atac <- CreateSeuratObject(HP2108601_assay, assay = "ATAC", meta.data=md.HP2108601)
  
  HP2108901_assay <- CreateChromatinAssay(HP2108901.counts, fragments = frags.HP2108901, min.features = 100)
  HP2108901_atac <- CreateSeuratObject(HP2108901_assay, assay = "ATAC", meta.data=md.HP2108901)
  
  HP2110001_assay <- CreateChromatinAssay(HP2110001.counts, fragments = frags.HP2110001, min.features = 100)
  HP2110001_atac <- CreateSeuratObject(HP2110001_assay, assay = "ATAC", meta.data=md.HP2110001)
  
  HP2121601_assay <- CreateChromatinAssay(HP2121601.counts, fragments = frags.HP2121601, min.features = 100)
  HP2121601_atac <- CreateSeuratObject(HP2121601_assay, assay = "ATAC", meta.data=md.HP2121601)
  
  HP2123201_assay <- CreateChromatinAssay(HP2123201.counts, fragments = frags.HP2123201, min.features = 100)
  HP2123201_atac <- CreateSeuratObject(HP2123201_assay, assay = "ATAC", meta.data=md.HP2123201)
  
  HP2132801_assay <- CreateChromatinAssay(HP2132801.counts, fragments = frags.HP2132801, min.features = 100)
  HP2132801_atac <- CreateSeuratObject(HP2132801_assay, assay = "ATAC", meta.data=md.HP2132801)
  
  HP2202101_assay <- CreateChromatinAssay(HP2202101.counts, fragments = frags.HP2202101, min.features = 100)
  HP2202101_atac <- CreateSeuratObject(HP2202101_assay, assay = "ATAC", meta.data=md.HP2202101)
}

# Sample specific Metadata addition
{
  HP2022801_atac$sample <- "HP2022801"
  SAMN15877725_atac$sample <- "SAMN15877725"
  HP2024001_atac$sample <- "HP2024001"
  HP2031401_atac$sample <- "HP2031401"
  HP2105501_atac$sample <- "HP2105501"
  HP2106201_atac$sample <- "HP2106201"
  HP2107001_atac$sample <- "HP2107001"
  HP2107901_atac$sample <- "HP2107901"
  HP2108601_atac$sample <- "HP2108601"
  HP2108901_atac$sample <- "HP2108901"
  HP2110001_atac$sample <- "HP2110001"
  HP2121601_atac$sample <- "HP2121601"
  HP2123201_atac$sample <- "HP2123201"
  HP2132801_atac$sample <- "HP2132801"
  HP2202101_atac$sample <- "HP2202101"
  
  # Sex specific Metadata addition
  HP2022801_atac$sex <- "female"
  SAMN15877725_atac$sex <- "male"
  HP2024001_atac$sex <- "female"
  HP2031401_atac$sex <- "male"
  HP2105501_atac$sex <- "female"
  HP2106201_atac$sex <- "female"
  HP2107001_atac$sex <- "male"
  HP2107901_atac$sex <- "male"
  HP2108601_atac$sex <- "female"
  HP2108901_atac$sex <- "female"
  HP2110001_atac$sex <- "male"
  HP2121601_atac$sex <- "female"
  HP2123201_atac$sex <- "male"
  HP2132801_atac$sex <- "female"
  HP2202101_atac$sex <- "female"
  
  # Ancestry specific Metadata addition
  HP2022801_atac$ancestry <- "white"
  SAMN15877725_atac$ancestry <- "white"
  HP2024001_atac$ancestry <- "white"
  HP2031401_atac$ancestry <- "black"
  HP2105501_atac$ancestry <- "white"
  HP2106201_atac$ancestry <- "black"
  HP2107001_atac$ancestry <- "white"
  HP2107901_atac$ancestry <- "white"
  HP2108601_atac$ancestry <- "white"
  HP2108901_atac$ancestry <- "white"
  HP2110001_atac$ancestry <- "black"
  HP2121601_atac$ancestry <- "black"
  HP2123201_atac$ancestry <- "black"
  HP2132801_atac$ancestry <- "black"
  HP2202101_atac$ancestry <- "black"
            
  # Ancestry and sex specific Metadata addition
  HP2022801_atac$ancestry_sex <- "white_female"
  SAMN15877725_atac$ancestry_sex <- "white_male"
  HP2024001_atac$ancestry_sex <- "white_female"
  HP2031401_atac$ancestry_sex <- "black_male"
  HP2105501_atac$ancestry_sex <- "white_female"
  HP2106201_atac$ancestry_sex <- "black_female"
  HP2107001_atac$ancestry_sex <- "white_male"
  HP2107901_atac$ancestry_sex <- "white_male"
  HP2108601_atac$ancestry_sex <- "white_female"
  HP2108901_atac$ancestry_sex <- "white_female"
  HP2110001_atac$ancestry_sex <- "black_male"
  HP2121601_atac$ancestry_sex <- "black_female"
  HP2123201_atac$ancestry_sex <- "black_male"
  HP2132801_atac$ancestry_sex <- "black_female"
  HP2202101_atac$ancestry_sex <- "black_female"
                            
  # Ancestry, sex and assay specific Metadata addition
  HP2022801_atac$ancestry_sex_atac <- "white_female_atac"
  SAMN15877725_atac$ancestry_sex_atac <- "white_male_atac"
  HP2024001_atac$ancestry_sex_atac <- "white_female_atac"
  HP2031401_atac$ancestry_sex_atac <- "black_male_atac"
  HP2105501_atac$ancestry_sex_atac <- "white_female_atac"
  HP2106201_atac$ancestry_sex_atac <- "black_female_atac"
  HP2107001_atac$ancestry_sex_atac <- "white_male_atac"
  HP2107901_atac$ancestry_sex_atac <- "white_male_atac"
  HP2108601_atac$ancestry_sex_atac <- "white_female_atac"
  HP2108901_atac$ancestry_sex_atac <- "white_female_atac"
  HP2110001_atac$ancestry_sex_atac <- "black_male_atac"
  HP2121601_atac$ancestry_sex_atac <- "black_female_atac"
  HP2123201_atac$ancestry_sex_atac <- "black_male_atac"
  HP2132801_atac$ancestry_sex_atac <- "black_female_atac"
  HP2202101_atac$ancestry_sex_atac <- "black_female_atac"
}

# Add annotations
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'

# Gen annotation
genome(annotations) <- "hg38"

# add gene information to the object
{
  Annotation(HP2022801_atac) <- annotations
  Annotation(SAMN15877725_atac) <- annotations
  Annotation(HP2024001_atac) <- annotations
  Annotation(HP2031401_atac) <- annotations
  Annotation(HP2105501_atac) <- annotations
  Annotation(HP2106201_atac) <- annotations
  Annotation(HP2107001_atac) <- annotations
  Annotation(HP2107901_atac) <- annotations
  Annotation(HP2108601_atac) <- annotations
  Annotation(HP2108901_atac) <- annotations
  Annotation(HP2110001_atac) <- annotations
  Annotation(HP2121601_atac) <- annotations
  Annotation(HP2123201_atac) <- annotations
  Annotation(HP2132801_atac) <- annotations
  Annotation(HP2202101_atac) <- annotations
  
  # compute nucleosome signal score per cell
  HP2022801_atac <- NucleosomeSignal(object = HP2022801_atac)
  SAMN15877725_atac <- NucleosomeSignal(object = SAMN15877725_atac)
  HP2024001_atac <- NucleosomeSignal(object = HP2024001_atac)
  HP2031401_atac <- NucleosomeSignal(object = HP2031401_atac)
  HP2105501_atac <- NucleosomeSignal(object = HP2105501_atac)
  HP2106201_atac <- NucleosomeSignal(object = HP2106201_atac)
  HP2107001_atac <- NucleosomeSignal(object = HP2107001_atac)
  HP2107901_atac <- NucleosomeSignal(object = HP2107901_atac)
  HP2108601_atac <- NucleosomeSignal(object = HP2108601_atac)
  HP2108901_atac <- NucleosomeSignal(object = HP2108901_atac)
  HP2110001_atac <- NucleosomeSignal(object = HP2110001_atac)
  HP2121601_atac <- NucleosomeSignal(object = HP2121601_atac)
  HP2123201_atac <- NucleosomeSignal(object = HP2123201_atac)
  HP2132801_atac <- NucleosomeSignal(object = HP2132801_atac)
  HP2202101_atac <- NucleosomeSignal(object = HP2202101_atac)
  
  # compute TSS enrichment score per cell
  HP2022801_atac <- TSSEnrichment(object = HP2022801_atac, fast = FALSE)
  SAMN15877725_atac <- TSSEnrichment(object = SAMN15877725_atac, fast = FALSE)
  HP2024001_atac <- TSSEnrichment(object = HP2024001_atac, fast = FALSE)
  HP2031401_atac <- TSSEnrichment(object = HP2031401_atac, fast = FALSE)
  HP2105501_atac <- TSSEnrichment(object = HP2105501_atac, fast = FALSE)
  HP2106201_atac <- TSSEnrichment(object = HP2106201_atac, fast = FALSE)
  HP2107001_atac <- TSSEnrichment(object = HP2107001_atac, fast = FALSE)
  HP2107901_atac <- TSSEnrichment(object = HP2107901_atac, fast = FALSE)
  HP2108601_atac <- TSSEnrichment(object = HP2108601_atac, fast = FALSE)
  HP2108901_atac <- TSSEnrichment(object = HP2108901_atac, fast = FALSE)
  HP2110001_atac <- TSSEnrichment(object = HP2110001_atac, fast = FALSE)
  HP2121601_atac <- TSSEnrichment(object = HP2121601_atac, fast = FALSE)
  HP2123201_atac <- TSSEnrichment(object = HP2123201_atac, fast = FALSE)
  HP2132801_atac <- TSSEnrichment(object = HP2132801_atac, fast = FALSE)
  HP2202101_atac <- TSSEnrichment(object = HP2202101_atac, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  HP2022801_atac$pct_reads_in_peaks <- HP2022801_atac$peak_region_fragments / HP2022801_atac$passed_filters * 100
  HP2022801_atac$blacklist_ratio <- HP2022801_atac$blacklist_region_fragments / HP2022801_atac$peak_region_fragments
  
  SAMN15877725_atac$pct_reads_in_peaks <- SAMN15877725_atac$peak_region_fragments / SAMN15877725_atac$passed_filters * 100
  SAMN15877725_atac$blacklist_ratio <- SAMN15877725_atac$blacklist_region_fragments / SAMN15877725_atac$peak_region_fragments
  
  HP2024001_atac$pct_reads_in_peaks <- HP2024001_atac$peak_region_fragments / HP2024001_atac$passed_filters * 100
  HP2024001_atac$blacklist_ratio <- HP2024001_atac$blacklist_region_fragments / HP2024001_atac$peak_region_fragments
  
  HP2031401_atac$pct_reads_in_peaks <- HP2031401_atac$peak_region_fragments / HP2031401_atac$passed_filters * 100
  HP2031401_atac$blacklist_ratio <- HP2031401_atac$blacklist_region_fragments / HP2031401_atac$peak_region_fragments
  
  HP2105501_atac$pct_reads_in_peaks <- HP2105501_atac$peak_region_fragments / HP2105501_atac$passed_filters * 100
  HP2105501_atac$blacklist_ratio <- HP2105501_atac$blacklist_region_fragments / HP2105501_atac$peak_region_fragments
  
  HP2106201_atac$pct_reads_in_peaks <- HP2106201_atac$peak_region_fragments / HP2106201_atac$passed_filters * 100
  HP2106201_atac$blacklist_ratio <- HP2106201_atac$blacklist_region_fragments / HP2106201_atac$peak_region_fragments
  
  HP2107001_atac$pct_reads_in_peaks <- HP2107001_atac$peak_region_fragments / HP2107001_atac$passed_filters * 100
  HP2107001_atac$blacklist_ratio <- HP2107001_atac$blacklist_region_fragments / HP2107001_atac$peak_region_fragments
  
  HP2107901_atac$pct_reads_in_peaks <- HP2107901_atac$peak_region_fragments / HP2107901_atac$passed_filters * 100
  HP2107901_atac$blacklist_ratio <- HP2107901_atac$blacklist_region_fragments / HP2107901_atac$peak_region_fragments
  
  HP2108601_atac$pct_reads_in_peaks <- HP2108601_atac$peak_region_fragments / HP2108601_atac$passed_filters * 100
  HP2108601_atac$blacklist_ratio <- HP2108601_atac$blacklist_region_fragments / HP2108601_atac$peak_region_fragments
  
  HP2108901_atac$pct_reads_in_peaks <- HP2108901_atac$peak_region_fragments / HP2108901_atac$passed_filters * 100
  HP2108901_atac$blacklist_ratio <- HP2108901_atac$blacklist_region_fragments / HP2108901_atac$peak_region_fragments
  
  HP2110001_atac$pct_reads_in_peaks <- HP2110001_atac$peak_region_fragments / HP2110001_atac$passed_filters * 100
  HP2110001_atac$blacklist_ratio <- HP2110001_atac$blacklist_region_fragments / HP2110001_atac$peak_region_fragments
  
  HP2121601_atac$pct_reads_in_peaks <- HP2121601_atac$peak_region_fragments / HP2121601_atac$passed_filters * 100
  HP2121601_atac$blacklist_ratio <- HP2121601_atac$blacklist_region_fragments / HP2121601_atac$peak_region_fragments
  
  HP2123201_atac$pct_reads_in_peaks <- HP2123201_atac$peak_region_fragments / HP2123201_atac$passed_filters * 100
  HP2123201_atac$blacklist_ratio <- HP2123201_atac$blacklist_region_fragments / HP2123201_atac$peak_region_fragments
  
  HP2132801_atac$pct_reads_in_peaks <- HP2132801_atac$peak_region_fragments / HP2132801_atac$passed_filters * 100
  HP2132801_atac$blacklist_ratio <- HP2132801_atac$blacklist_region_fragments / HP2132801_atac$peak_region_fragments
  
  HP2202101_atac$pct_reads_in_peaks <- HP2202101_atac$peak_region_fragments / HP2202101_atac$passed_filters * 100
  HP2202101_atac$blacklist_ratio <- HP2202101_atac$blacklist_region_fragments / HP2202101_atac$peak_region_fragments
}

# view QC
HP2022801_atac$high.tss <- ifelse(HP2022801_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2022801_atac, group.by = 'high.tss') + NoLegend()

SAMN15877725_atac$high.tss <- ifelse(SAMN15877725_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(SAMN15877725_atac, group.by = 'high.tss') + NoLegend()

HP2024001_atac$high.tss <- ifelse(HP2024001_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2024001_atac, group.by = 'high.tss') + NoLegend()

HP2031401_atac$high.tss <- ifelse(HP2031401_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2031401_atac, group.by = 'high.tss') + NoLegend()

HP2105501_atac$high.tss <- ifelse(HP2105501_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2105501_atac, group.by = 'high.tss') + NoLegend()

HP2106201_atac$high.tss <- ifelse(HP2106201_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2106201_atac, group.by = 'high.tss') + NoLegend()

HP2107001_atac$high.tss <- ifelse(HP2107001_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2107001_atac, group.by = 'high.tss') + NoLegend()

HP2107901_atac$high.tss <- ifelse(HP2107901_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2107901_atac, group.by = 'high.tss') + NoLegend()

HP2108601_atac$high.tss <- ifelse(HP2108601_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2108601_atac, group.by = 'high.tss') + NoLegend()

HP2108901_atac$high.tss <- ifelse(HP2108901_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2108901_atac, group.by = 'high.tss') + NoLegend()

HP2110001_atac$high.tss <- ifelse(HP2110001_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2110001_atac, group.by = 'high.tss') + NoLegend()

HP2121601_atac$high.tss <- ifelse(HP2121601_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2121601_atac, group.by = 'high.tss') + NoLegend()

HP2123201_atac$high.tss <- ifelse(HP2123201_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2123201_atac, group.by = 'high.tss') + NoLegend()

HP2132801_atac$high.tss <- ifelse(HP2132801_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2132801_atac, group.by = 'high.tss') + NoLegend()

HP2202101_atac$high.tss <- ifelse(HP2202101_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2202101_atac, group.by = 'high.tss') + NoLegend()

# View Nucleosome signal  
HP2022801_atac$nucleosome_group <- ifelse(HP2022801_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2022801_atac, group.by = 'nucleosome_group')

SAMN15877725_atac$nucleosome_group <- ifelse(SAMN15877725_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = SAMN15877725_atac, group.by = 'nucleosome_group')

HP2024001_atac$nucleosome_group <- ifelse(HP2024001_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2024001_atac, group.by = 'nucleosome_group')

HP2031401_atac$nucleosome_group <- ifelse(HP2031401_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2031401_atac, group.by = 'nucleosome_group')

HP2105501_atac$nucleosome_group <- ifelse(HP2105501_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2105501_atac, group.by = 'nucleosome_group')

HP2106201_atac$nucleosome_group <- ifelse(HP2106201_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2106201_atac, group.by = 'nucleosome_group')

HP2107001_atac$nucleosome_group <- ifelse(HP2107001_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2107001_atac, group.by = 'nucleosome_group')

HP2107901_atac$nucleosome_group <- ifelse(HP2107901_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2107901_atac, group.by = 'nucleosome_group')

HP2108601_atac$nucleosome_group <- ifelse(HP2108601_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2108601_atac, group.by = 'nucleosome_group')

HP2108901_atac$nucleosome_group <- ifelse(HP2108901_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2108901_atac, group.by = 'nucleosome_group')

HP2110001_atac$nucleosome_group <- ifelse(HP2110001_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2110001_atac, group.by = 'nucleosome_group')

HP2121601_atac$nucleosome_group <- ifelse(HP2121601_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2121601_atac, group.by = 'nucleosome_group')

HP2123201_atac$nucleosome_group <- ifelse(HP2123201_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2123201_atac, group.by = 'nucleosome_group')

HP2132801_atac$nucleosome_group <- ifelse(HP2132801_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2132801_atac, group.by = 'nucleosome_group')

HP2202101_atac$nucleosome_group <- ifelse(HP2202101_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2202101_atac, group.by = 'nucleosome_group')

# Count fragments
HP2022801_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
SAMN15877725_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2024001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2031401_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2105501_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2106201_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2107001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2107901_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2108601_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2108901_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2110001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2121601_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2123201_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2132801_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2202101_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)

rownames(HP2022801_atac_fraginfo) <- HP2022801_atac_fraginfo$CB
rownames(SAMN15877725_atac_fraginfo) <- SAMN15877725_atac_fraginfo$CB
rownames(HP2024001_atac_fraginfo) <- HP2024001_atac_fraginfo$CB
rownames(HP2031401_atac_fraginfo) <- HP2031401_atac_fraginfo$CB
rownames(HP2105501_atac_fraginfo) <- HP2105501_atac_fraginfo$CB
rownames(HP2106201_atac_fraginfo) <- HP2106201_atac_fraginfo$CB
rownames(HP2107001_atac_fraginfo) <- HP2107001_atac_fraginfo$CB
rownames(HP2107901_atac_fraginfo) <- HP2107901_atac_fraginfo$CB
rownames(HP2108601_atac_fraginfo) <- HP2108601_atac_fraginfo$CB
rownames(HP2108901_atac_fraginfo) <- HP2108901_atac_fraginfo$CB
rownames(HP2110001_atac_fraginfo) <- HP2110001_atac_fraginfo$CB
rownames(HP2121601_atac_fraginfo) <- HP2121601_atac_fraginfo$CB
rownames(HP2123201_atac_fraginfo) <- HP2123201_atac_fraginfo$CB
rownames(HP2132801_atac_fraginfo) <- HP2132801_atac_fraginfo$CB
rownames(HP2202101_atac_fraginfo) <- HP2202101_atac_fraginfo$CB

HP2022801_atac$fragments <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "frequency_count"]
HP2022801_atac$mononucleosomal <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "mononucleosomal"]
HP2022801_atac$nucleosome_free <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "nucleosome_free"]
HP2022801_atac$reads_count <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "reads_count"]

SAMN15877725_atac$fragments <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "frequency_count"]
SAMN15877725_atac$mononucleosomal <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "mononucleosomal"]
SAMN15877725_atac$nucleosome_free <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "nucleosome_free"]
SAMN15877725_atac$reads_count <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "reads_count"]

HP2024001_atac$fragments <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "frequency_count"]
HP2024001_atac$mononucleosomal <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "mononucleosomal"]
HP2024001_atac$nucleosome_free <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "nucleosome_free"]
HP2024001_atac$reads_count <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "reads_count"]

HP2031401_atac$fragments <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "frequency_count"]
HP2031401_atac$mononucleosomal <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "mononucleosomal"]
HP2031401_atac$nucleosome_free <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "nucleosome_free"]
HP2031401_atac$reads_count <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "reads_count"]

HP2105501_atac$fragments <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "frequency_count"]
HP2105501_atac$mononucleosomal <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "mononucleosomal"]
HP2105501_atac$nucleosome_free <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "nucleosome_free"]
HP2105501_atac$reads_count <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "reads_count"]

HP2106201_atac$fragments <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "frequency_count"]
HP2106201_atac$mononucleosomal <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "mononucleosomal"]
HP2106201_atac$nucleosome_free <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "nucleosome_free"]
HP2106201_atac$reads_count <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "reads_count"]

HP2107001_atac$fragments <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "frequency_count"]
HP2107001_atac$mononucleosomal <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "mononucleosomal"]
HP2107001_atac$nucleosome_free <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "nucleosome_free"]
HP2107001_atac$reads_count <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "reads_count"]

HP2107901_atac$fragments <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "frequency_count"]
HP2107901_atac$mononucleosomal <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "mononucleosomal"]
HP2107901_atac$nucleosome_free <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "nucleosome_free"]
HP2107901_atac$reads_count <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "reads_count"]

HP2108601_atac$fragments <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "frequency_count"]
HP2108601_atac$mononucleosomal <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "mononucleosomal"]
HP2108601_atac$nucleosome_free <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "nucleosome_free"]
HP2108601_atac$reads_count <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "reads_count"]

HP2108901_atac$fragments <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "frequency_count"]
HP2108901_atac$mononucleosomal <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "mononucleosomal"]
HP2108901_atac$nucleosome_free <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "nucleosome_free"]
HP2108901_atac$reads_count <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "reads_count"]

HP2110001_atac$fragments <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "frequency_count"]
HP2110001_atac$mononucleosomal <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "mononucleosomal"]
HP2110001_atac$nucleosome_free <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "nucleosome_free"]
HP2110001_atac$reads_count <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "reads_count"]

HP2121601_atac$fragments <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "frequency_count"]
HP2121601_atac$mononucleosomal <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "mononucleosomal"]
HP2121601_atac$nucleosome_free <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "nucleosome_free"]
HP2121601_atac$reads_count <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "reads_count"]

HP2123201_atac$fragments <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "frequency_count"]
HP2123201_atac$mononucleosomal <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "mononucleosomal"]
HP2123201_atac$nucleosome_free <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "nucleosome_free"]
HP2123201_atac$reads_count <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "reads_count"]

HP2132801_atac$fragments <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "frequency_count"]
HP2132801_atac$mononucleosomal <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "mononucleosomal"]
HP2132801_atac$nucleosome_free <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "nucleosome_free"]
HP2132801_atac$reads_count <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "reads_count"]

HP2202101_atac$fragments <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "frequency_count"]
HP2202101_atac$mononucleosomal <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "mononucleosomal"]
HP2202101_atac$nucleosome_free <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "nucleosome_free"]
HP2202101_atac$reads_count <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "reads_count"]


HP2022801_atac <- FRiP(
  object = HP2022801_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

SAMN15877725_atac <- FRiP(
  object = SAMN15877725_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2024001_atac <- FRiP(
  object = HP2024001_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2031401_atac <- FRiP(
  object = HP2031401_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2105501_atac <- FRiP(
  object = HP2105501_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2106201_atac <- FRiP(
  object = HP2106201_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2107001_atac <- FRiP(
  object = HP2107001_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2107901_atac <- FRiP(
  object = HP2107901_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2108601_atac <- FRiP(
  object = HP2108601_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2108901_atac <- FRiP(
  object = HP2108901_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2110001_atac <- FRiP(
  object = HP2110001_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2121601_atac <- FRiP(
  object = HP2121601_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2123201_atac <- FRiP(
  object = HP2123201_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2132801_atac <- FRiP(
  object = HP2132801_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2202101_atac <- FRiP(
  object = HP2202101_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

# Add unique cell names otherwise integration will give errors
{
  HP2022801_atac <- RenameCells(object = HP2022801_atac, add.cell.id = "HP2022801")
  SAMN15877725_atac <- RenameCells(object = SAMN15877725_atac, add.cell.id = "SAMN15877725")
  HP2024001_atac <- RenameCells(object = HP2024001_atac, add.cell.id = "HP2024001")
  HP2031401_atac <- RenameCells(object = HP2031401_atac, add.cell.id = "HP2031401")
  HP2105501_atac <- RenameCells(object = HP2105501_atac, add.cell.id = "HP2105501")
  HP2106201_atac <- RenameCells(object = HP2106201_atac, add.cell.id = "HP2106201")
  HP2107001_atac <- RenameCells(object = HP2107001_atac, add.cell.id = "HP2107001")
  HP2107901_atac <- RenameCells(object = HP2107901_atac, add.cell.id = "HP2107901")
  HP2108601_atac <- RenameCells(object = HP2108601_atac, add.cell.id = "HP2108601")
  HP2108901_atac <- RenameCells(object = HP2108901_atac, add.cell.id = "HP2108901")
  HP2110001_atac <- RenameCells(object = HP2110001_atac, add.cell.id = "HP2110001")
  HP2121601_atac <- RenameCells(object = HP2121601_atac, add.cell.id = "HP2121601")
  HP2123201_atac <- RenameCells(object = HP2123201_atac, add.cell.id = "HP2123201")
  HP2132801_atac <- RenameCells(object = HP2132801_atac, add.cell.id = "HP2132801")
  HP2202101_atac <- RenameCells(object = HP2202101_atac, add.cell.id = "HP2202101")
}

# Visualize QC
VlnPlot(
  object = HP2022801_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = SAMN15877725_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2024001_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2031401_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2105501_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2106201_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2107001_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2107901_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2108601_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2108901_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2110001_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2121601_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2123201_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2132801_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2202101_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

# Add Gene activity matrix
{
  gene.activities.HP2022801 <- GeneActivity(HP2022801_atac)
  gene.activities.SAMN15877725 <- GeneActivity(SAMN15877725_atac)
  gene.activities.HP2024001 <- GeneActivity(HP2024001_atac)
  gene.activities.HP2031401 <- GeneActivity(HP2031401_atac)
  gene.activities.HP2105501 <- GeneActivity(HP2105501_atac)
  gene.activities.HP2106201 <- GeneActivity(HP2106201_atac)
  gene.activities.HP2107001 <- GeneActivity(HP2107001_atac)
  gene.activities.HP2107901 <- GeneActivity(HP2107901_atac)
  gene.activities.HP2108601 <- GeneActivity(HP2108601_atac)
  gene.activities.HP2108901 <- GeneActivity(HP2108901_atac)
  gene.activities.HP2110001 <- GeneActivity(HP2110001_atac)
  gene.activities.HP2121601 <- GeneActivity(HP2121601_atac)
  gene.activities.HP2123201 <- GeneActivity(HP2123201_atac)
  gene.activities.HP2132801 <- GeneActivity(HP2132801_atac)
  gene.activities.HP2202101 <- GeneActivity(HP2202101_atac)
}

# add the gene activity matrix to the Seurat object as a new assay and normalize it
{
  HP2022801_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2022801)
  HP2022801_atac <- NormalizeData(
    object = HP2022801_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2022801_atac$nCount_RNA)
  )
  
  SAMN15877725_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.SAMN15877725)
  SAMN15877725_atac <- NormalizeData(
    object = SAMN15877725_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(SAMN15877725_atac$nCount_RNA)
  )
  
  HP2024001_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2024001)
  HP2024001_atac <- NormalizeData(
    object = HP2024001_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2024001_atac$nCount_RNA)
  )
  
  HP2031401_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2031401)
  HP2031401_atac <- NormalizeData(
    object = HP2031401_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2031401_atac$nCount_RNA)
  )
  
  HP2105501_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2105501)
  HP2105501_atac <- NormalizeData(
    object = HP2105501_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2105501_atac$nCount_RNA)
  )
  
  HP2106201_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2106201)
  HP2106201_atac <- NormalizeData(
    object = HP2106201_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2106201_atac$nCount_RNA)
  )
  
  HP2107001_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2107001)
  HP2107001_atac <- NormalizeData(
    object = HP2107001_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2107001_atac$nCount_RNA)
  )
  
  HP2107901_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2107901)
  HP2107901_atac <- NormalizeData(
    object = HP2107901_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2107901_atac$nCount_RNA)
  )
  
  HP2108601_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2108601)
  HP2108601_atac <- NormalizeData(
    object = HP2108601_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2108601_atac$nCount_RNA)
  )
  
  HP2108901_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2108901)
  HP2108901_atac <- NormalizeData(
    object = HP2108901_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2108901_atac$nCount_RNA)
  )
  
  HP2110001_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2110001)
  HP2110001_atac <- NormalizeData(
    object = HP2110001_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2110001_atac$nCount_RNA)
  )
  
  HP2121601_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2121601)
  HP2121601_atac <- NormalizeData(
    object = HP2121601_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2121601_atac$nCount_RNA)
  )
  
  HP2123201_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2123201)
  HP2123201_atac <- NormalizeData(
    object = HP2123201_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2123201_atac$nCount_RNA)
  )
  
  HP2132801_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2132801)
  HP2132801_atac <- NormalizeData(
    object = HP2132801_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2132801_atac$nCount_RNA)
  )
  
  HP2202101_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2202101)
  HP2202101_atac <- NormalizeData(
    object = HP2202101_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2202101_atac$nCount_RNA)
  )
}


# merge all datasets, adding a cell ID to make sure cell names are unique
combined_atac <- merge(
  x = HP2022801_atac,
  y = list(SAMN15877725_atac, HP2024001_atac, HP2031401_atac, HP2105501_atac,
           HP2106201_atac, HP2107001_atac, HP2107901_atac, HP2108601_atac, 
           HP2108901_atac, HP2110001_atac, HP2121601_atac, HP2123201_atac,
           HP2132801_atac, HP2202101_atac))

VlnPlot(
  object = combined_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6, raster = FALSE
)

# QC Cleanup
combined_atac <- subset(
  x = combined_atac, 
  subset = peak_region_fragments > 2000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 30 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    FRiP > 0.2
)
combined_atac

#Save file
#qsave(combined_atac, file = r"(E:\2.SexbasedStudyCurrent\QS files\combined_atac.qs)")
#qsave(combined_atac_doublet, file = r"(E:\2.SexbasedStudyCurrent\QS files\combined_atac_doublet.qs)")
combined_atac <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\snATACfiles_earlierpartsofworkflow\combined_atac.qs)")
combined_atac_doublet <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\snATACfiles_earlierpartsofworkflow\combined_atac_doublet.qs)")

# Run TFDIF  
combined_atac <- FindTopFeatures(combined_atac, min.cutoff = 20)
combined_atac <- RunTFIDF(combined_atac)
combined_atac <- RunSVD(combined_atac)
combined_atac <- RunUMAP(combined_atac, dims = 2:30, reduction = 'lsi')
DimPlot(combined_atac, group.by = 'ancestry_sex', pt.size = 0.1)

# Batch correction using Harmony
hm.integrated <- RunHarmony(object = combined_atac, group.by.vars = 'sample', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated, group.by = 'ancestry_sex', pt.size = 0.1)

# We need to subset out all of those cell IDs in our doublet calculated data which are present in actual ATAC data
#create cell names as metadata colum
atac_cellnames <- colnames(combined_atac)

#subset the doublet file to have only those barcodes that have been identified as high quality cells
combined_atac_doublet <- subset(combined_atac_doublet, rownames(combined_atac_doublet) %in% atac_cellnames)

# Adding a column to define multiplet vs singlet
combined_atac_doublet$doublets <- ifelse(combined_atac_doublet$p.value >= 0.05, "singlet", "multiplet" )

# Check number of singlets and multiplets
sum(combined_atac_doublet$doublets == "singlet")
sum(combined_atac_doublet$doublets == "multiplet")

#Some more checking because AddMetadata was driving me INSANE
#cells_seuratobj <- Cells(combined_atac)
#cells_netadatobj <- rownames(combined_atac_doublet)
#setequal(cells_seuratobj, cells_netadatobj)

# subset out the columns you wish to add only, here we are only adding the doublets column
doubletdat <- subset(combined_atac_doublet, select = c("doublets"))

# Add doublet metadata to seurat object
hm.integrated <- AddMetaData(object = hm.integrated, metadata = doubletdat, col.name = 'doublets')
table(hm.integrated$doublets) #party time

# Check doublets
#Idents(hm.integrated) <- "doublets"
DimPlot(hm.integrated, group.by = 'doublets', pt.size = 0.1)

# Remove doublets
Idents(hm.integrated) <- "doublets" 
hm.integrated.dfree <- subset(hm.integrated, idents = c("singlet"))
hm.integrated.dfree

# Re-Run TFDIF  
hm.integrated.dfree <- FindTopFeatures(hm.integrated.dfree, min.cutoff = 20)
hm.integrated.dfree <- RunTFIDF(hm.integrated.dfree)
hm.integrated.dfree <- RunSVD(hm.integrated.dfree)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'lsi')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Batch correction using Harmony
hm.integrated.dfree <- RunHarmony(object = hm.integrated.dfree, group.by.vars = 'sample', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Normalize gene activities
DefaultAssay(hm.integrated.dfree) <- "RNA"
hm.integrated.dfree <- NormalizeData(hm.integrated.dfree)
hm.integrated.dfree <- ScaleData(hm.integrated.dfree, features = rownames(hm.integrated.dfree))

# Clustering
DefaultAssay(hm.integrated.dfree) <- "ATAC"
hm.integrated.dfree <- FindNeighbors(object = hm.integrated.dfree, reduction = 'harmony', dims = 2:30)
hm.integrated.dfree <- FindClusters(object = hm.integrated.dfree, verbose = FALSE, algorithm = 3)
DimPlot(object = hm.integrated.dfree, label = TRUE) + NoLegend()

#Subset out poor quality cells
Idents(hm.integrated.dfree) <- "seurat_clusters"
hm.integrated.dfree <- subset(hm.integrated.dfree, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                                              #"11", 
                                                              "12", "13", #"14", "15", 
                                                              "16", #"17", 
                                                              "18", "19", "20", "21"))

# Re-Run TFDIF  
hm.integrated.dfree <- FindTopFeatures(hm.integrated.dfree, min.cutoff = 20)
hm.integrated.dfree <- RunTFIDF(hm.integrated.dfree)
hm.integrated.dfree <- RunSVD(hm.integrated.dfree)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'lsi')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Normalize gene activities
DefaultAssay(hm.integrated.dfree) <- "RNA"
hm.integrated.dfree <- NormalizeData(hm.integrated.dfree)
hm.integrated.dfree <- ScaleData(hm.integrated.dfree, features = rownames(hm.integrated.dfree))

# Batch correction using Harmony
hm.integrated.dfree <- RunHarmony(object = hm.integrated.dfree, group.by.vars = 'sample', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Clusterning UMAP was created on basis of ATAC profile
DefaultAssay(hm.integrated.dfree) <- "ATAC"
hm.integrated.dfree <- FindNeighbors(object = hm.integrated.dfree, reduction = 'harmony', dims = 2:30)
hm.integrated.dfree <- FindClusters(object = hm.integrated.dfree, verbose = FALSE, algorithm = 3)
DimPlot(object = hm.integrated.dfree, label = TRUE) #+ NoLegend()

#Save file
#qsave(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
#combined_atac_doublet <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac_doublet.rds)")

# Open necessary scRNAseq and snATAC data
#hm.integrated.dfree <- readRDS(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
#pancreas.combined.h.s <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")

# Identify anchors
# Read in RNA data
DefaultAssay(pancreas.combined.h.s) <- 'RNA'
system.time({
#  user      system  elapsed 
#  25456.72  2710.84 28241.15 ~7.8hrs
transfer.anchors <- FindTransferAnchors(reference = pancreas.combined.h.s, 
                                        query = hm.integrated.dfree, 
                                        features = pancreas.combined.h.s@assays$RNA@var.features,
                                        reference.assay = "RNA", 
                                        query.assay = "RNA", 
                                        reduction = "cca")
})

# Save
#qsave(transfer.anchors, file = r"(E:\2.SexbasedStudyCurrent\QS files\transfer.anchors.qs)")

#### START VERY EARLY FROM HERE ####
####                            ####
transfer.anchors <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\transfer.anchors.qs)")
pancreas.combined.h.s <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
hm.integrated.dfree <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

# Annotation of scATAC cells via label transfer  
# map query onto the reference dataset
DefaultAssay(pancreas.combined.h.s) <- "RNA"
# you dont need to run this because return.model = T when generating UMAP for scRNAseq
#DimPlot(pancreas.combined.h.s, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
#pancreas.combined.h.s <- RunUMAP(object = pancreas.combined.h.s, assay = "SCT", reduction = "harmony", dims = 1:30, return.model = TRUE) # return model = TRUE

# View clustering
DefaultAssay(hm.integrated.dfree) <- "ATAC"
DimPlot(object = hm.integrated.dfree, label = TRUE) + NoLegend()

# Query mapping
hm.integrated.dfree <- MapQuery(
  anchorset = transfer.anchors,
  reference = pancreas.combined.h.s,
  query = hm.integrated.dfree,
  refdata = pancreas.combined.h.s$celltype_qadir,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

# Clusterning UMAP was created on basis of ATAC profile
# Add predicted ID data to the new Signac seurat object
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = pancreas.combined.h.s$celltype_qadir,
                                     weight.reduction = hm.integrated.dfree[["harmony"]], 
                                     dims = 2:30)

hm.integrated.dfree <- AddMetaData(hm.integrated.dfree, metadata = celltype.predictions)

# View UMAP
DimPlot(hm.integrated.dfree, group.by = "predicted.id", reduction = "umap", label = TRUE) + ggtitle("Predicted annotation")
DimPlot(hm.integrated.dfree, group.by = "ATAC_snn_res.0.8", reduction = "umap", label = TRUE) + ggtitle("Res")# + nolegend()
DimPlot(pancreas.combined.h.s, group.by = "celltype_qadir", reduction = "umap", label = TRUE) + ggtitle("Celltype Classification")

# rename Idents ans save as celltype
Idents(hm.integrated.dfree) <- "ATAC_snn_res.0.8"
hm.integrated.dfree <- RenameIdents(hm.integrated.dfree,
                                    "0" = "alpha", 
                                    "1" = "beta",
                                    "2" = "alpha", 
                                    "3" = "alpha",
                                    "4" = "beta", 
                                    "5" = "activated_stellate", #(PDGFRA+)
                                    "6" = "ductal", 
                                    "7" = "delta",
                                    "8" = "acinar", 
                                    "9" = "alpha",
                                    "10" = "beta", 
                                    "11" = "gamma",
                                    "12" = "endothelial",
                                    "13" = "endothelial", 
                                    "14" = "quiescent_stellate", #(RGS5+)
                                    "15" = "macrophage",
                                    "16" = "lymphocyte"
)
table(Idents(hm.integrated.dfree))
DimPlot(hm.integrated.dfree, reduction = "umap", label = TRUE)

# Saving this information in the metadata slot
table(Idents(hm.integrated.dfree))
hm.integrated.dfree$celltype <- Idents(hm.integrated.dfree)
head(hm.integrated.dfree@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "alpha", "delta", "gamma",
               "ductal", "acinar", 
               "activated_stellate", "quiescent_stellate", "endothelial",
               "macrophage", "lymphocyte"
)
# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree@meta.data$celltype <- factor(x = hm.integrated.dfree@meta.data$celltype, levels = my_levels)
Idents(hm.integrated.dfree) <- "celltype"
table(Idents(hm.integrated.dfree))

# Observing cells
Idents(hm.integrated.dfree) <- "celltype"
DimPlot(hm.integrated.dfree, 
        #split.by = "ancestry_sex", 
        #group.by = "celltype", 
        label = FALSE, 
        ncol = 2,  
        cols = c("dodgerblue3", #"beta"
                 "lightseagreen", #"alpha"
                 "chartreuse3", #"delta"
                 "springgreen4", #"gamma"
                 "darkorange2", #"ductal"
                 "salmon3", #"acinar"
                 "orange", #"activated-stellate"
                 "salmon", #"quiescent-stellate"
                 "red", #"endothelial"
                 "magenta3", #"macrophages"
                 "orchid1" #"lymphocyte"
                 )
) + ggtitle("Celltype Classification")


# Plotting
DefaultAssay(hm.integrated.dfree) <- "RNA"
FeaturePlot(
  object = hm.integrated.dfree,
  features = c('SST', 'INS'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  reduction = "umap",
  ncol = 3
)

# Add metadata for DT
Idents(object = hm.integrated.dfree) <- "celltype"
hm.integrated.dfree$celltype_sex <- paste(Idents(hm.integrated.dfree), hm.integrated.dfree$sex, sep = "_")
hm.integrated.dfree$celltype_sex_ancestry <- paste(Idents(hm.integrated.dfree), hm.integrated.dfree$sex, hm.integrated.dfree$ancestry, sep = "_")
table(hm.integrated.dfree@meta.data$celltype_sex)
table(hm.integrated.dfree@meta.data$celltype_sex_ancestry)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("beta_female_black", "beta_male_black", "beta_female_white", "beta_male_white", 
                "alpha_female_black", "alpha_male_black", "alpha_female_white", "alpha_male_white", 
                "delta_female_black", "delta_male_black", "delta_female_white", "delta_male_white",
                "gamma_female_black", "gamma_male_black", "gamma_female_white", "gamma_male_white",
                "ductal_female_black", "ductal_male_black", "ductal_female_white", "ductal_male_white",
                "acinar_female_black", "acinar_male_black", "acinar_female_white", "acinar_male_white",
                "activated_stellate_female_black", "activated_stellate_male_black", "activated_stellate_female_white", "activated_stellate_male_white", 
                "quiescent_stellate_female_black", "quiescent_stellate_male_black", "quiescent_stellate_female_white", "quiescent_stellate_male_white", 
                "endothelial_female_black", "endothelial_male_black", "endothelial_female_white", "endothelial_male_white", 
                "macrophage_female_black", "macrophage_male_black", "macrophage_female_white", "macrophage_male_white",
                "lymphocyte_female_black", "lymphocyte_male_black", "lymphocyte_female_white", "lymphocyte_male_white"
)

my_levels3 <- c("beta_female", "beta_male", 
                "alpha_female", "alpha_male", 
                "delta_female", "delta_male",
                "gamma_female", "gamma_male",
                "ductal_female", "ductal_male",
                "acinar_female", "acinar_male",
                "activated_stellate_female", "activated_stellate_male", 
                "quiescent_stellate_female", "quiescent_stellate_male", 
                "endothelial_female", "endothelial_male", 
                "macrophage_female", "macrophage_male",
                "lymphocyte_female", "lymphocyte_male"
)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree@meta.data$celltype_sex_ancestry <- factor(x = hm.integrated.dfree@meta.data$celltype_sex_ancestry, levels = my_levels2)
hm.integrated.dfree@meta.data$celltype_sex <- factor(x = hm.integrated.dfree@meta.data$celltype_sex, levels = my_levels3)
table(hm.integrated.dfree@meta.data$celltype_sex_ancestry)
table(hm.integrated.dfree@meta.data$celltype_sex)

#Peak calling
CoveragePlot(
  object = hm.integrated.dfree,
  group.by = "celltype",
  region = c("INS", 
             "GCG",
             "SST",
             "PPY"
             ),
  show.bulk = TRUE,
  #ymax = 80,
  ranges.title = "Ranges",
  window = 200,
  extend.upstream = 25000,
  extend.downstream = 10000,
) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) & 
  scale_fill_manual(values = c("dodgerblue3", #"beta"
                               "lightseagreen", #"alpha"
                               "chartreuse3", #"delta"
                               "springgreen4", #"gamma"
                               "darkorange2", #"ductal"
                               "salmon3", #"acinar"
                               "orange", #"activated-stellate"
                               "salmon", #"quiescent-stellate"
                               "red", #"endothelial"
                               "magenta3", #"macrophages"
                               "orchid1", #"lymphocyte"
                               "grey"
)) 

# Bask in the beauty of your optimised seurat object
hm.integrated.dfree[['ATAC']]
granges(hm.integrated.dfree)
table(hm.integrated.dfree@assays[["ATAC"]]@annotation@seqinfo@genome)

# Change back to peak data
DefaultAssay(hm.integrated.dfree) <- "ATAC"
Idents(hm.integrated.dfree) <- "celltype"

# Adding predicted RNA gene expression values
rna <- TransferData(
  anchorset = transfer.anchors,
  refdata = GetAssayData(pancreas.combined.h.s, assay = "RNA", slot = "data"),
  weight.reduction = hm.integrated.dfree[["lsi"]],
  dims = 2:30
)

# add predicted values as a new assay
hm.integrated.dfree[["predicted"]] <- rna

# Adding Motifs
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
hm.integrated.dfree <- AddMotifs(
  object = hm.integrated.dfree,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

#Compute motif activity
hm.integrated.dfree <- RunChromVAR(
  object = hm.integrated.dfree,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(hm.integrated.dfree) <- 'chromvar'

# look at the activity of Mef2c
# FeaturePlot(
#   object = hm.integrated.dfree,
#   features = "MA0497.1",
#   min.cutoff = 'q10',
#   max.cutoff = 'q90',
#   pt.size = 0.1
# )

#### START FROM HERE ####
####                 ####
#Save file
#qsave(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
hm.integrated.dfree <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

# Pseudobulk
Idents(hm.integrated.dfree) <- "celltype"
hm.integrated.dfree$celltype_sex_ancestry_lib <- paste(hm.integrated.dfree$celltype, hm.integrated.dfree$sex, hm.integrated.dfree$ancestry, hm.integrated.dfree$sample, sep = "_")
table(hm.integrated.dfree@meta.data$celltype_sex_ancestry_lib)
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry_lib"
DefaultAssay(hm.integrated.dfree) <- "ATAC"
combined_processed_atac <- Seurat:::PseudobulkExpression(object = hm.integrated.dfree, 
                                                         pb.method = 'aggregate', 
                                                         return.seurat = TRUE,
                                                         slot = 'counts')

DefaultAssay(combined_processed_atac) <- "ATAC"
combined_processed_atac <- RunTFIDF(combined_processed_atac, assay = "ATAC")

{
  combined_processed_atac$celltype_sex_ancestry_lib <- combined_processed_atac@active.ident
  Idents(combined_processed_atac) <- 'celltype_sex_ancestry_lib'
  combined_processed_atac$celltype <- combined_processed_atac$orig.ident
  metadat <- combined_processed_atac@meta.data
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "activated", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "quiescent", "quiescent-stellate"))
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex_ancestry_lib, "_", -3)
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -2)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -1)
  combined_processed_atac@meta.data = metadat
}

table(combined_processed_atac@meta.data[["celltype"]])
table(combined_processed_atac@meta.data[["sex"]])
table(combined_processed_atac@meta.data[["ancestry"]])
table(combined_processed_atac@meta.data[["lib"]])

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("delta", "beta", "alpha", "gamma",
               "ductal", "acinar",
               "activated", "quiescent", "endothelial",
               "lymphocyte", "macrophage") 

table(combined_processed_atac$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
table(combined_processed_atac$celltype)
Idents(combined_processed_atac) <- "celltype"

# DE testing to determine celltype specificity
DefaultAssay(combined_processed_atac) <- "ATAC"
Idents(combined_processed_atac) <- "celltype"
table(Idents(combined_processed_atac))

#combined_processed_atac[["ATAC"]]@counts<-as.matrix(combined_processed_rna[["RNA"]]@counts)+1
# Runtests
plan(strategy = "multicore", workers = 80)
system.time({
beta.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "beta", 
                                      min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                      only.pos = TRUE, verbose = TRUE)

alpha.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "alpha", 
                                      min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                      only.pos = TRUE, verbose = TRUE)

delta.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "delta", 
                                       min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                       only.pos = TRUE,verbose = TRUE)

gamma.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "gamma", 
                                       min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                       only.pos = TRUE, verbose = TRUE)

ductal.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "ductal", 
                                       min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                       only.pos = TRUE, verbose = TRUE)

acinar.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "acinar", 
                                        min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                        only.pos = TRUE, verbose = TRUE)

activated.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "activated", 
                                           min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                           only.pos = TRUE, verbose = TRUE)

quiescent.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "quiescent", 
                                           min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                           only.pos = TRUE,verbose = TRUE)

endothelial.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "endothelial", 
                                             min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                             only.pos = TRUE, verbose = TRUE)

lymphocyte.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "lymphocyte", 
                                             min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                             only.pos = TRUE, verbose = TRUE)

macrophage.conserved.markers <- FindMarkers(combined_processed_atac, min.pct = 1, logfc.threshold = 1, assay = "ATAC", slot = "counts", ident.1 = "macrophage", 
                                            min.cells.group = 1, test.use = "LR", densify = TRUE, latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                            only.pos = TRUE, verbose = TRUE)
# Save data
write.csv(beta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\beta.conserved.markers.csv)")
write.csv(alpha.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\alpha.conserved.markers.csv)")
write.csv(delta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\delta.conserved.markers.csv)")
write.csv(gamma.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\gamma.conserved.markers.csv)")
write.csv(ductal.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\ductal.conserved.markers.csv)")
write.csv(acinar.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\acinar.conserved.markers.csv)")
write.csv(activated.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\activated.conserved.markers.csv)")
write.csv(quiescent.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\quiescent.conserved.markers.csv)")
write.csv(endothelial.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\endothelial.conserved.markers.csv)")
write.csv(lymphocyte.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\lymphocyte.conserved.markers.csv)")
write.csv(macrophage.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\conserved\macrophages.conserved.markers.csv)")
})

# Testing across sex
# Advanced Differential Gene testing
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = hm.integrated.dfree) <- "celltype"
combined_processed_atac@meta.data$celltype_sex <- paste(combined_processed_atac$celltype, combined_processed_atac$sex, sep = "_")
table(combined_processed_atac$celltype_sex)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("beta_female", "beta_male",
                "alpha_female", "alpha_male",
                "delta_female", "delta_male",
                "gamma_female", "gamma_male",
                "ductal_female", "ductal_male",
                "acinar_female", "acinar_male",
                "activated_female", "activated_male",
                "quiescent_female", "quiescent_male",
                "endothelial_female", "endothelial_male",
                "macrophage_female", "macrophage_male",
                "lymphocyte_female", "lymphocyte_male"
)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_atac@meta.data$celltype_sex <- factor(x = combined_processed_atac@meta.data$celltype_sex, levels = my_levels2)
table(combined_processed_atac@meta.data$celltype_sex)

# Change back to peak data
DefaultAssay(combined_processed_atac) <- "ATAC"
Idents(combined_processed_atac) <- "celltype_sex"
table(Idents(combined_processed_atac))

# DE testing accessible sites sex_celltype
plan(strategy = "multicore", workers = 80)
system.time({
# 1.Beta-cells ####
# MALE VS. FEMALE
beta.mvsf <- FindMarkers(combined_processed_atac, 
                         ident.1 = "beta_male", ident.2 = "beta_female", 
                         slot = "counts",
                         assay = "ATAC",
                         test.use = "LR",
                         latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                         min.cells.group = 1,
                         min.pct = 1,
                         logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                         densify = TRUE,
                         verbose = TRUE)

# 2.Alpha-cells ####
# MALE VS. FEMALE
alpha.mvsf <- FindMarkers(combined_processed_atac, 
                          ident.1 = "alpha_male", ident.2 = "alpha_female", 
                          slot = "counts",
                          assay = "ATAC",
                          test.use = "LR",
                          latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                          min.cells.group = 1,
                          min.pct = 1,
                          logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                          densify = TRUE,
                          verbose = TRUE)

# 3.Delta-cells ####
# MALE VS. FEMALE
delta.mvsf <- FindMarkers(combined_processed_atac, 
                          ident.1 = "delta_male", ident.2 = "delta_female", 
                          slot = "counts",
                          assay = "ATAC",
                          test.use = "LR",
                          latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                          min.cells.group = 1,
                          min.pct = 1,
                          logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                          densify = TRUE,
                          verbose = TRUE)

# 4.Gamma-cells ####
# MALE VS. FEMALE
gamma.mvsf <- FindMarkers(combined_processed_atac, 
                          ident.1 = "gamma_male", ident.2 = "gamma_female", 
                          slot = "counts",
                          assay = "ATAC",
                          test.use = "LR",
                          latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                          min.cells.group = 1,
                          min.pct = 1,
                          logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                          densify = TRUE,
                          verbose = TRUE)

# 5.Ductal-cells ####
# MALE VS. FEMALE
ductal.mvsf <- FindMarkers(combined_processed_atac, 
                           ident.1 = "ductal_male", ident.2 = "ductal_female", 
                           slot = "counts",
                           assay = "ATAC",
                           test.use = "LR",
                           latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                           min.cells.group = 1,
                           min.pct = 1,
                           logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                           densify = TRUE,
                           verbose = TRUE)

# 6.Acinar-cells ####
# MALE VS. FEMALE
acinar.mvsf <- FindMarkers(combined_processed_atac, 
                           ident.1 = "acinar_male", ident.2 = "acinar_female", 
                           slot = "counts",
                           assay = "ATAC",
                           test.use = "LR",
                           latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                           min.cells.group = 1,
                           min.pct = 1,
                           logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                           densify = TRUE,
                           verbose = TRUE)

# 7.Activated-cells ####
# MALE VS. FEMALE
activated.mvsf <- FindMarkers(combined_processed_atac, 
                              ident.1 = "activated_male", ident.2 = "activated_female", 
                              slot = "counts",
                              assay = "ATAC",
                              test.use = "LR",
                              latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                              min.cells.group = 1,
                              min.pct = 1,
                              logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                              densify = TRUE,
                              verbose = TRUE)

# 8.Quiescent-cells ####
# MALE VS. FEMALE
quiescent.mvsf <- FindMarkers(combined_processed_atac, 
                              ident.1 = "quiescent_male", ident.2 = "quiescent_female", 
                              slot = "counts",
                              assay = "ATAC",
                              test.use = "LR",
                              latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                              min.cells.group = 1,
                              min.pct = 1,
                              logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                              densify = TRUE,
                              verbose = TRUE)

# 9.Endothelial-cells ####
# MALE VS. FEMALE
endothelial.mvsf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "endothelial_male", ident.2 = "endothelial_female", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# 10.Macrophage-cells ####
# MALE VS. FEMALE
macrophage.mvsf <- FindMarkers(combined_processed_atac, 
                               ident.1 = "macrophage_male", ident.2 = "macrophage_female", 
                               slot = "counts",
                               assay = "ATAC",
                               test.use = "LR",
                               latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                               min.cells.group = 1,
                               min.pct = 1,
                               logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                               densify = TRUE,
                               verbose = TRUE)

# 11.Lymphocyte-cells ####
# MALE VS. FEMALE
lymphocyte.mvsf <- FindMarkers(combined_processed_atac, 
                               ident.1 = "lymphocyte_male", ident.2 = "lymphocyte_female", 
                               slot = "counts",
                               assay = "ATAC",
                               test.use = "LR",
                               latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                               min.cells.group = 1,
                               min.pct = 1,
                               logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                               densify = TRUE,
                               verbose = TRUE)

# Save data
write.csv(beta.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\beta.mvsf.csv)")
write.csv(alpha.mvsfs, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\alpha.mvsf.csv)")
write.csv(delta.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\delta.mvsf.csv)")
write.csv(gamma.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\gamma.mvsf.csv)")
write.csv(ductal.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\ductal.mvsf.csv)")
write.csv(acinar.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\acinar.mvsf.csv)")
write.csv(activated.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\activated.mvsf.csv)")
write.csv(quiescent.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\quiescent.mvsf.csv)")
write.csv(endothelial.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\endothelial.mvsf.csv)")
write.csv(lymphocyte.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\lymphocyte.mvsf.csv)")
write.csv(macrophages.mvsf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\macrophages.mvsf.csv)")
})


# Advanced Differential Gene testing
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = combined_processed_atac) <- "celltype"
combined_processed_atac$celltype_sex_ancestry <- paste(combined_processed_atac$celltype, combined_processed_atac$sex, combined_processed_atac$ancestry, sep = "_")
table(combined_processed_atac@meta.data$celltype_sex_ancestry)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("beta_female_white", "beta_male_white", "beta_female_black", "beta_male_black",
                "alpha_female_white", "alpha_male_white", "alpha_female_black", "alpha_male_black",
                "delta_female_white", "delta_male_white", "delta_female_black", "delta_male_black",
                "gamma_female_white", "gamma_male_white", "gamma_female_black", "gamma_male_black",
                "ductal_female_white", "ductal_male_white", "ductal_female_black", "ductal_male_black",
                "acinar_female_white", "acinar_male_white", "acinar_female_black", "acinar_male_black",
                "activated_female_white", "activated_male_white", "activated_female_black", "activated_male_black",
                "quiescent_female_white", "quiescent_male_white", "quiescent_female_black", "quiescent_male_black",
                "endothelial_female_white", "endothelial_male_white", "endothelial_female_black", "endothelial_male_black",
                "macrophage_female_white", "macrophage_male_white", "macrophage_female_black", "macrophage_male_black",
                "lymphocyte_female_white", "lymphocyte_male_white", "lymphocyte_female_black", "lymphocyte_male_black"
)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_atac@meta.data$celltype_sex_ancestry <- factor(x = combined_processed_atac@meta.data$celltype_sex_ancestry, levels = my_levels2)
table(combined_processed_atac@meta.data$celltype_sex_ancestry)

# Change back to peak data
DefaultAssay(combined_processed_atac) <- "ATAC"
Idents(combined_processed_atac) <- "celltype_sex_ancestry"
table(Idents(combined_processed_atac))

plan(strategy = "multicore", workers = 80)
system.time({
# DE testing accessible sites sex_ancestry_celltype
# 1.Beta-cells ####
# WHITE MALE VS. WHITE FEMALE
beta.wmvswf <- FindMarkers(combined_processed_atac, 
                           ident.1 = "beta_male_white", ident.2 = "beta_female_white", 
                           slot = "counts",
                           assay = "ATAC",
                           test.use = "LR",
                           latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                           min.cells.group = 1,
                           min.pct = 1,
                           logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                           densify = TRUE,
                           verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
beta.bmvsbf <- FindMarkers(combined_processed_atac, 
                           ident.1 = "beta_male_black", ident.2 = "beta_female_black", 
                           slot = "counts",
                           assay = "ATAC",
                           test.use = "LR",
                           latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                           min.cells.group = 1,
                           min.pct = 1,
                           logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                           densify = TRUE,
                           verbose = TRUE)

# WHITE MALE VS. BLACK MALE
beta.wmvsbm <- FindMarkers(combined_processed_atac, 
                           ident.1 = "beta_male_white", ident.2 = "beta_male_black", 
                           slot = "counts",
                           assay = "ATAC",
                           test.use = "LR",
                           latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                           min.cells.group = 1,
                           min.pct = 1,
                           logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                           densify = TRUE,
                           verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
beta.wfvsbf <- FindMarkers(combined_processed_atac, 
                           ident.1 = "beta_female_white", ident.2 = "beta_female_black", 
                           slot = "counts",
                           assay = "ATAC",
                           test.use = "LR",
                           latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                           min.cells.group = 1,
                           min.pct = 1,
                           logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                           densify = TRUE,
                           verbose = TRUE)

# 2.Alpha-cells ####
# WHITE MALE VS. WHITE FEMALE
alpha.wmvswf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "alpha_male_white", ident.2 = "alpha_female_white", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
alpha.bmvsbf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "alpha_male_black", ident.2 = "alpha_female_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# WHITE MALE VS. BLACK MALE
alpha.wmvsbm <- FindMarkers(combined_processed_atac, 
                            ident.1 = "alpha_male_white", ident.2 = "alpha_male_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
alpha.wfvsbf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "alpha_female_white", ident.2 = "alpha_female_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# 3.Delta-cells ####
# WHITE MALE VS. WHITE FEMALE
delta.wmvswf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "delta_male_white", ident.2 = "delta_female_white", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
delta.bmvsbf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "delta_male_black", ident.2 = "delta_female_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# WHITE MALE VS. BLACK MALE
delta.wmvsbm <- FindMarkers(combined_processed_atac, 
                            ident.1 = "delta_male_white", ident.2 = "delta_male_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
delta.wfvsbf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "delta_female_white", ident.2 = "delta_female_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# 4.Gamma-cells ####
# WHITE MALE VS. WHITE FEMALE
gamma.wmvswf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "gamma_male_white", ident.2 = "gamma_female_white", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
gamma.bmvsbf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "gamma_male_black", ident.2 = "gamma_female_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# WHITE MALE VS. BLACK MALE
gamma.wmvsbm <- FindMarkers(combined_processed_atac, 
                            ident.1 = "gamma_male_white", ident.2 = "gamma_male_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
gamma.wfvsbf <- FindMarkers(combined_processed_atac, 
                            ident.1 = "gamma_female_white", ident.2 = "gamma_female_black", 
                            slot = "counts",
                            assay = "ATAC",
                            test.use = "LR",
                            latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                            min.cells.group = 1,
                            min.pct = 1,
                            logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                            densify = TRUE,
                            verbose = TRUE)

# 5.Ductal-cells ####
# WHITE MALE VS. WHITE FEMALE
ductal.wmvswf <- FindMarkers(combined_processed_atac, 
                             ident.1 = "ductal_male_white", ident.2 = "ductal_female_white", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
ductal.bmvsbf <- FindMarkers(combined_processed_atac, 
                             ident.1 = "ductal_male_black", ident.2 = "ductal_female_black", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)

# WHITE MALE VS. BLACK MALE
ductal.wmvsbm <- FindMarkers(combined_processed_atac, 
                             ident.1 = "ductal_male_white", ident.2 = "ductal_male_black", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
ductal.wfvsbf <- FindMarkers(combined_processed_atac, 
                             ident.1 = "ductal_female_white", ident.2 = "ductal_female_black", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)

# 6.Acinar-cells ####
# WHITE MALE VS. WHITE FEMALE
acinar.wmvswf <- FindMarkers(combined_processed_atac, 
                             ident.1 = "acinar_male_white", ident.2 = "acinar_female_white", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
acinar.bmvsbf <- FindMarkers(combined_processed_atac, 
                             ident.1 = "acinar_male_black", ident.2 = "acinar_female_black", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)

# WHITE MALE VS. BLACK MALE
acinar.wmvsbm <- FindMarkers(combined_processed_atac, 
                             ident.1 = "acinar_male_white", ident.2 = "acinar_male_black", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
acinar.wfvsbf <- FindMarkers(combined_processed_atac, 
                             ident.1 = "acinar_female_white", ident.2 = "acinar_female_black", 
                             slot = "counts",
                             assay = "ATAC",
                             test.use = "LR",
                             latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                             min.cells.group = 1,
                             min.pct = 1,
                             logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                             densify = TRUE,
                             verbose = TRUE)
    
# 7.Activated-cells ####
# WHITE MALE VS. WHITE FEMALE
activated.wmvswf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "activated_male_white", ident.2 = "activated_female_white", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
activated.bmvsbf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "activated_male_black", ident.2 = "activated_female_black", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# WHITE MALE VS. BLACK MALE
activated.wmvsbm <- FindMarkers(combined_processed_atac, 
                                ident.1 = "activated_male_white", ident.2 = "activated_male_black", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
activated.wfvsbf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "activated_female_white", ident.2 = "activated_female_black", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE) 

# 8.Quiescent-cells ####
# WHITE MALE VS. WHITE FEMALE
quiescent.wmvswf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "quiescent_male_white", ident.2 = "quiescent_female_white", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
quiescent.bmvsbf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "quiescent_male_black", ident.2 = "quiescent_female_black", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# WHITE MALE VS. BLACK MALE
quiescent.wmvsbm <- FindMarkers(combined_processed_atac, 
                                ident.1 = "quiescent_male_white", ident.2 = "quiescent_male_black", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
quiescent.wfvsbf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "quiescent_female_white", ident.2 = "quiescent_female_black", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE) 

# 9.Endothelial-cells ####
# WHITE MALE VS. WHITE FEMALE
endothelial.wmvswf <- FindMarkers(combined_processed_atac, 
                                  ident.1 = "endothelial_male_white", ident.2 = "endothelial_female_white", 
                                  slot = "counts",
                                  assay = "ATAC",
                                  test.use = "LR",
                                  latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                  min.cells.group = 1,
                                  min.pct = 1,
                                  logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                  densify = TRUE,
                                  verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
endothelial.bmvsbf <- FindMarkers(combined_processed_atac, 
                                ident.1 = "endothelial_male_black", ident.2 = "endothelial_female_black", 
                                slot = "counts",
                                assay = "ATAC",
                                test.use = "LR",
                                latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                min.cells.group = 1,
                                min.pct = 1,
                                logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                densify = TRUE,
                                verbose = TRUE)

# WHITE MALE VS. BLACK MALE
endothelial.wmvsbm <- FindMarkers(combined_processed_atac, 
                                  ident.1 = "endothelial_male_white", ident.2 = "endothelial_male_black", 
                                  slot = "counts",
                                  assay = "ATAC",
                                  test.use = "LR",
                                  latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                  min.cells.group = 1,
                                  min.pct = 1,
                                  logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                  densify = TRUE,
                                  verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
endothelial.wfvsbf <- FindMarkers(combined_processed_atac, 
                                  ident.1 = "endothelial_female_white", ident.2 = "endothelial_female_black", 
                                  slot = "counts",
                                  assay = "ATAC",
                                  test.use = "LR",
                                  latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                  min.cells.group = 1,
                                  min.pct = 1,
                                  logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                  densify = TRUE,
                                  verbose = TRUE) 

# 10.Macrophage-cells ####
# WHITE MALE VS. WHITE FEMALE
macrophage.wmvswf <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "macrophage_male_white", ident.2 = "macrophage_female_white", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
macrophage.bmvsbf <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "macrophage_male_black", ident.2 = "macrophage_female_black", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE)

# WHITE MALE VS. BLACK MALE
macrophage.wmvsbm <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "macrophage_male_white", ident.2 = "macrophage_male_black", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
macrophage.wfvsbf <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "macrophage_female_white", ident.2 = "macrophage_female_black", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE) 

# 11.Lymphocyte-cells ####
# WHITE MALE VS. WHITE FEMALE
lymphocyte.wmvswf <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "lymphocyte_male_white", ident.2 = "lymphocyte_female_white", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE)

# BLACK MALE VS. BLACK FEMALE
lymphocyte.bmvsbf <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "lymphocyte_male_black", ident.2 = "lymphocyte_female_black", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE)

# WHITE MALE VS. BLACK MALE
lymphocyte.wmvsbm <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "lymphocyte_male_white", ident.2 = "lymphocyte_male_black", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE)

# WHITE FEMALE VS. BLACK FEMALE
lymphocyte.wfvsbf <- FindMarkers(combined_processed_atac, 
                                 ident.1 = "lymphocyte_female_white", ident.2 = "lymphocyte_female_black", 
                                 slot = "counts",
                                 assay = "ATAC",
                                 test.use = "LR",
                                 latent.vars = c('lib', 'celltype_sex_ancestry_lib'),
                                 min.cells.group = 1,
                                 min.pct = 1,
                                 logfc.threshold = 1, # based on output log2 so 1 is 2 FC
                                 densify = TRUE,
                                 verbose = TRUE) 

# Save data
write.csv(beta.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\beta.wmvswf.csv)")
write.csv(beta.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\beta.bmvsbf.csv)")
write.csv(beta.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\beta.wmvsbm)")
write.csv(beta.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\beta.wfvsbf.csv)")

write.csv(alpha.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\alpha.wmvswf.csv)")
write.csv(alpha.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\alpha.bmvsbf.csv)")
write.csv(alpha.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\alpha.wmvsbm)")
write.csv(alpha.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\alpha.wfvsbf.csv)")

write.csv(delta.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\delta.wmvswf.csv)")
write.csv(delta.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\delta.bmvsbf.csv)")
write.csv(delta.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\delta.wmvsbm)")
write.csv(delta.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\delta.wfvsbf.csv)")

write.csv(gamma.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\gamma.wmvswf.csv)")
write.csv(gamma.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\gamma.bmvsbf.csv)")
write.csv(gamma.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\gamma.wmvsbm)")
write.csv(gamma.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\gamma.wfvsbf.csv)")

write.csv(acinar.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\acinar.wmvswf.csv)")
write.csv(acinar.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\acinar.bmvsbf.csv)")
write.csv(acinar.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\acinar.wmvsbm)")
write.csv(acinar.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\acinar.wfvsbf.csv)")

write.csv(ductal.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\ductal.wmvswf.csv)")
write.csv(ductal.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\ductal.bmvsbf.csv)")
write.csv(ductal.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\ductal.wmvsbm)")
write.csv(ductal.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\ductal.wfvsbf.csv)")

write.csv(activated.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\activated.wmvswf.csv)")
write.csv(activated.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\activated.bmvsbf.csv)")
write.csv(activated.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\activated.wmvsbm)")
write.csv(activated.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\activated.wfvsbf.csv)")

write.csv(quiescent.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\quiescent.wmvswf.csv)")
write.csv(quiescent.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\quiescent.bmvsbf.csv)")
write.csv(quiescent.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\quiescent.wmvsbm)")
write.csv(quiescent.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\quiescent.wfvsbf.csv)")

write.csv(endothelial.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\endothelial.wmvswf.csv)")
write.csv(endothelial.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\endothelial.bmvsbf.csv)")
write.csv(endothelial.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\endothelial.wmvsbm)")
write.csv(endothelial.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\endothelial.wfvsbf.csv)")

write.csv(lymphocyte.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\lymphocyte.wmvswf.csv)")
write.csv(lymphocyte.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\lymphocyte.bmvsbf.csv)")
write.csv(lymphocyte.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\lymphocyte.wmvsbm)")
write.csv(lymphocyte.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\lymphocyte.wfvsbf.csv)")

write.csv(macrophages.wmvswf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\macrophages.wmvswf.csv)")
write.csv(macrophages.bmvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\macrophages.bmvsbf.csv)")
write.csv(macrophages.wmvsbm, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\macrophages.wmvsbm)")
write.csv(macrophages.wfvsbf, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\macrophages.wfvsbf.csv)")
})

















 
