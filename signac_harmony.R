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
HP2022801_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\1_220628 Fahd_snATAC1_HP-20228-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
SAMN15877725_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\2_220701 Fahd_snATAC2_SAMN15877725\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2024001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\3_220701 Fahd_snATAC3_HP-20240-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2031401_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\4_220630 Fahd_snATAC4_HP-20314-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2105501_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\5_220303_snATAC_F52_HP-21055-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2106201_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\6_210401 snATAC_F62_HP-21062-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2107001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\7_210401 snATAC_F7a_HP-21070-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2107901_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\9_210628 snATAC_F9a_HP-21079-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2108601_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\10_210628 snATAC_F10a_HP-21086-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2108901_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\11_210714 snATAC_F11a_HP-21089-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2110001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\12_210714 snATAC_F12a_HP-21100-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2121601_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\13_211208_snATAC_F13_HP-21216-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2123201_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\14_211208_snATAC_F14_HP-21232-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2132801_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2202101_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw data\snATACseq\16_220630 Fahd_snATAC16_HP-22021-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)

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
#saveRDS(combined_atac, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac.rds)")
#combined_atac <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac.rds)")
#combined_atac_doublet <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac_doublet.rds)")

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

#Save file
#saveRDS(hm.integrated.dfree, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\hm.integrated.dfree.rds)")
#combined_atac_doublet <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac_doublet.rds)")

# Open necessary scRNAseq and snATAC data
hm.integrated.dfree <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\hm.integrated.dfree.rds)")
pancreas.combined.h.s <- readRDS(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\scRNAseq\CoreArch_CurrentR\processed_seurat\pancreas.combined.h.s.rds)")

# Identify anchors
# Read in RNA data
DefaultAssay(pancreas.combined.h.s) <- 'RNA'
transfer.anchors <- FindTransferAnchors(reference = pancreas.combined.h.s, 
                                        query = hm.integrated.dfree, 
                                        features = pancreas.combined.h.s@assays$RNA@var.features,
                                        reference.assay = "RNA", 
                                        query.assay = "RNA", 
                                        reduction = "cca")

# Save
#saveRDS(transfer.anchors, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\transfer.anchors.rds)")
transfer.anchors <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\transfer.anchors.rds)")

# Annotation of scATAC cells via label transfer  
# map query onto the reference dataset
DefaultAssay(pancreas.combined.h.s) <- "SCT"
DimPlot(pancreas.combined.h.s, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
pancreas.combined.h.s <- RunUMAP(object = pancreas.combined.h.s, assay = "SCT", reduction = "harmony", dims = 1:30, return.model = TRUE) # return model = TRUE

# Clusterning
hm.integrated.dfree <- FindNeighbors(object = hm.integrated.dfree, reduction = 'lsi', dims = 2:30)
hm.integrated.dfree <- FindClusters(object = hm.integrated.dfree, verbose = FALSE, algorithm = 3)
DimPlot(object = hm.integrated.dfree, label = TRUE) + NoLegend()


hm.integrated.dfree <- MapQuery(
  anchorset = transfer.anchors,
  reference = pancreas.combined.h.s,
  query = hm.integrated.dfree,
  refdata = pancreas.combined.h.s$celltype,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = pancreas.combined.h.s$celltype,
                                     weight.reduction = hm.integrated.dfree[["harmony"]], 
                                     dims = 2:30)

hm.integrated.dfree <- AddMetaData(hm.integrated.dfree, metadata = celltype.predictions)

DimPlot(hm.integrated.dfree, group.by = "predicted.id", reduction = "umap", label = TRUE) + ggtitle("Predicted annotation") # + nolegend()
DimPlot(pancreas.combined.h.s, group.by = "celltype", reduction = "umap", label = TRUE) + ggtitle("Celltype Classification")
p1+p2

# Clusterning
hm.integrated.dfree <- FindNeighbors(object = hm.integrated.dfree, reduction = 'lsi', dims = 2:30)
hm.integrated.dfree <- FindClusters(object = hm.integrated.dfree, verbose = FALSE, algorithm = 3)
DimPlot(object = hm.integrated.dfree, label = TRUE) + NoLegend()

#Save file
#saveRDS(combined_atac, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac.rds)")
#combined_atac <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac.rds)")




# # QC Cleanup
#   HP2022801_atac <- subset(
#     x = HP2022801_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2022801_atac
#   
#   SAMN15877725_atac <- subset(
#     x = SAMN15877725_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   SAMN15877725_atac
#   
#   HP2024001_atac <- subset(
#     x = HP2024001_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2024001_atac
#   
#   HP2031401_atac <- subset(
#     x = HP2031401_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2031401_atac
#   
#   HP2105501_atac <- subset(
#     x = HP2105501_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2105501_atac
#   
#   HP2106201_atac <- subset(
#     x = HP2106201_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2106201_atac
#   
#   HP2107001_atac <- subset(
#     x = HP2107001_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2107001_atac
#   
#   HP2107901_atac <- subset(
#     x = HP2107901_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2107901_atac
#   
#   HP2108601_atac <- subset(
#     x = HP2108601_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2108601_atac
#   
#   HP2108901_atac <- subset(
#     x = HP2108901_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2108901_atac
#   
#   HP2110001_atac <- subset(
#     x = HP2110001_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2110001_atac
#   
#   HP2121601_atac <- subset(
#     x = HP2121601_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2121601_atac
#   
#   HP2123201_atac <- subset(
#     x = HP2123201_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2123201_atac
#   
#   HP2132801_atac <- subset(
#     x = HP2132801_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2132801_atac
#   
#   HP2202101_atac <- subset(
#     x = HP2202101_atac,
#     subset = peak_region_fragments > 2000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 15 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
#   HP2202101_atac

# Normalization
# First compute LSI
# compute LSI
HP2022801_atac <- FindTopFeatures(HP2022801_atac, min.cutoff = 10)
HP2022801_atac <- RunTFIDF(HP2022801_atac)
HP2022801_atac <- RunSVD(HP2022801_atac)

# compute LSI
SAMN15877725_atac <- FindTopFeatures(SAMN15877725_atac, min.cutoff = 10)
SAMN15877725_atac <- RunTFIDF(SAMN15877725_atac)
SAMN15877725_atac <- RunSVD(SAMN15877725_atac)

# compute LSI
HP2024001_atac <- FindTopFeatures(HP2024001_atac, min.cutoff = 10)
HP2024001_atac <- RunTFIDF(HP2024001_atac)
HP2024001_atac <- RunSVD(HP2024001_atac)

# compute LSI
HP2031401_atac <- FindTopFeatures(HP2031401_atac, min.cutoff = 10)
HP2031401_atac <- RunTFIDF(HP2031401_atac)
HP2031401_atac <- RunSVD(HP2031401_atac)

# compute LSI
HP2105501_atac <- FindTopFeatures(HP2105501_atac, min.cutoff = 10)
HP2105501_atac <- RunTFIDF(HP2105501_atac)
HP2105501_atac <- RunSVD(HP2105501_atac)

# compute LSI
HP2106201_atac <- FindTopFeatures(HP2106201_atac, min.cutoff = 10)
HP2106201_atac <- RunTFIDF(HP2106201_atac)
HP2106201_atac <- RunSVD(HP2106201_atac)

# compute LSI
HP2107001_atac <- FindTopFeatures(HP2107001_atac, min.cutoff = 10)
HP2107001_atac <- RunTFIDF(HP2107001_atac)
HP2107001_atac <- RunSVD(HP2107001_atac)

# compute LSI
HP2107901_atac <- FindTopFeatures(HP2107901_atac, min.cutoff = 10)
HP2107901_atac <- RunTFIDF(HP2107901_atac)
HP2107901_atac <- RunSVD(HP2107901_atac)

# compute LSI
HP2108601_atac <- FindTopFeatures(HP2108601_atac, min.cutoff = 10)
HP2108601_atac <- RunTFIDF(HP2108601_atac)
HP2108601_atac <- RunSVD(HP2108601_atac)

# compute LSI
HP2108901_atac <- FindTopFeatures(HP2108901_atac, min.cutoff = 10)
HP2108901_atac <- RunTFIDF(HP2108901_atac)
HP2108901_atac <- RunSVD(HP2108901_atac)

# compute LSI
HP2110001_atac <- FindTopFeatures(HP2110001_atac, min.cutoff = 10)
HP2110001_atac <- RunTFIDF(HP2110001_atac)
HP2110001_atac <- RunSVD(HP2110001_atac)

# compute LSI
HP2121601_atac <- FindTopFeatures(HP2121601_atac, min.cutoff = 10)
HP2121601_atac <- RunTFIDF(HP2121601_atac)
HP2121601_atac <- RunSVD(HP2121601_atac)

# compute LSI
HP2123201_atac <- FindTopFeatures(HP2123201_atac, min.cutoff = 10)
HP2123201_atac <- RunTFIDF(HP2123201_atac)
HP2123201_atac <- RunSVD(HP2123201_atac)

# compute LSI
HP2132801_atac <- FindTopFeatures(HP2132801_atac, min.cutoff = 10)
HP2132801_atac <- RunTFIDF(HP2132801_atac)
HP2132801_atac <- RunSVD(HP2132801_atac)

# compute LSI
HP2202101_atac <- FindTopFeatures(HP2202101_atac, min.cutoff = 10)
HP2202101_atac <- RunTFIDF(HP2202101_atac)
HP2202101_atac <- RunSVD(HP2202101_atac)

# Add unique cell names otherwise integration will give errors
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

# head(x = colnames(x = HP2022801))
# head(x = colnames(x = SAMN15877725))
# head(x = colnames(x = HP2024001))
# head(x = colnames(x = HP2031401))
# head(x = colnames(x = HP2105501))
# head(x = colnames(x = HP2106201))
# head(x = colnames(x = HP2107001))
# head(x = colnames(x = HP2107901))
# head(x = colnames(x = HP2108601))
# head(x = colnames(x = HP2108901))
# head(x = colnames(x = HP2110001))
# head(x = colnames(x = HP2121601))
# head(x = colnames(x = HP2123201))
# head(x = colnames(x = HP2132801))
# head(x = colnames(x = HP2202101))

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
# combined_atac <- merge(HP2022801_atac, y = c(SAMN15877725_atac, HP2024001_atac, HP2031401_atac, HP2105501_atac,
#                        HP2106201_atac, HP2107001_atac, HP2107901_atac, HP2108601_atac, 
#                        HP2108901_atac, HP2110001_atac, HP2121601_atac, HP2123201_atac,
#                        HP2132801_atac, HP2202101_atac))

combined_atac <- merge(
  x = HP2022801_atac,
  y = list(SAMN15877725_atac, HP2024001_atac, HP2031401_atac, HP2105501_atac,
           HP2106201_atac, HP2107001_atac, HP2107901_atac, HP2108601_atac, 
           HP2108901_atac, HP2110001_atac, HP2121601_atac, HP2123201_atac,
           HP2132801_atac, HP2202101_atac)
  #   add.cell.ids = c("HP2022801", "SAMN15877725", "HP2024001", "HP2031401", "HP2105501",
  #                    "HP2106201", "HP2107001", "HP2107901", "HP2108601", "HP2108901",
  #                    "HP2110001", "HP2121601", "HP2123201", "HP2132801", "HP2202101")
)
combined_atac[["ATAC"]]


# Run TFDIF  
combined_atac <- FindTopFeatures(combined_atac, min.cutoff = 20)
combined_atac <- RunTFIDF(combined_atac)
combined_atac <- RunSVD(combined_atac)
combined_atac <- RunUMAP(combined_atac, dims = 2:30, reduction = 'lsi')
DimPlot(combined_atac, group.by = 'ancestry_sex', pt.size = 0.1)

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(HP2022801_atac, SAMN15877725_atac, HP2024001_atac, HP2031401_atac, 
                     HP2105501_atac, HP2106201_atac, HP2107001_atac, HP2107901_atac, 
                     HP2108601_atac, HP2108901_atac, HP2110001_atac, HP2121601_atac, 
                     HP2123201_atac, HP2132801_atac, HP2202101_atac),
  anchor.features = rownames(HP2022801_atac),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated_atac <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined_atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# We exclude the first dimension as this is typically correlated with sequencing depth
integrated_atac <- RunTFIDF(integrated_atac)
integrated_atac <- FindTopFeatures(integrated_atac, min.cutoff = "q0")
integrated_atac <- RunSVD(integrated_atac)

# create a new UMAP using the integrated embeddings
integrated_atac <- RunUMAP(integrated_atac, reduction = "integrated_lsi", 
                           dims = 2:30, 
                           return.model = TRUE,
                           seed.use = 42,
                           negative.sample.rate = 5,
                           a = NULL,
                           b = NULL)
DimPlot(integrated_atac, group.by = "ancestry_sex")

# Normalize gene activities
DefaultAssay(integrated_atac) <- "RNA"
integrated_atac <- NormalizeData(integrated_atac)
integrated_atac <- ScaleData(integrated_atac, features = rownames(integrated_atac))

# Using scRNAseq dataset to map anchors
pancreas.combined <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.combinedcorrectedSCT.rds)") # RNA dataset integrated

# Identify anchors
DefaultAssay(pancreas.combined) <- 'RNA'
transfer.anchors <- FindTransferAnchors(reference = pancreas.combined, 
                                        query = integrated_atac, 
                                        features = pancreas.combined@assays$integrated@var.features,
                                        reference.assay = "RNA", 
                                        query.assay = "RNA", 
                                        reduction = "cca")

# Save
saveRDS(transfer.anchors, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\transfer.anchors.rds)")
transfer.anchors <- readRDS(r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\transfer.anchors.rds)")

# Annotation of scATAC cells via label transfer  
# map query onto the reference dataset
DefaultAssay(pancreas.combined) <- "integrated"
pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", dims = 1:30, return.model = TRUE)
integrated_atac <- MapQuery(
  anchorset = transfer.anchors,
  reference = pancreas.combined,
  query = integrated_atac,
  refdata = pancreas.combined$celltype,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = pancreas.combined$celltype,
                                     weight.reduction = integrated_atac[["integrated_lsi"]], 
                                     dims = 2:30)

integrated_atac <- AddMetaData(integrated_atac, metadata = celltype.predictions)

DimPlot(integrated_atac, group.by = "predicted.id", label = TRUE) + ggtitle("Predicted annotation") # + nolegend()
DimPlot(pancreas.combined, group.by = "celltype", reduction = "umap", label = TRUE) + ggtitle("Celltype Classification")
p1+p2


# Visualization
VlnPlot(
  object = integrated_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

DefaultAssay(integrated_atac) <- 'ATAC'

FeaturePlot(
  object = integrated_atac,
  features = c('INS', 'GCG', 'SST', 'GHRL', 'PPY'),
  pt.size = 0.1,
  max.cutoff = '0.5',
  ncol = 2
)

FeaturePlot(
  object = integrated_atac,
  features = c('SDS', 'VWF', 'KRT19', 'CELA2A'),
  pt.size = 0.1,
  max.cutoff = '0.5',
  ncol = 2
)

FeaturePlot(
  object = integrated_atac,
  features = c("RGS5", "CSRP2", "FABP4", "COL3A1", "FMOD", "PDGFRB"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 2
)

FeaturePlot(
  object = integrated_atac,
  features = c("DDIT3"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 2
)

# Clustering
integrated_atac <- FindNeighbors(object = integrated_atac, reduction = 'integrated_lsi', dims = 2:50)
integrated_atac <- FindClusters(object = integrated_atac, verbose = FALSE, resolution = 0.6, algorithm = 3)
DimPlot(object = integrated_atac, label = TRUE) + NoLegend()

# Rename Idents
integrated_atac <- RenameIdents(integrated_atac, 
                                "0" = "Alpha", 
                                "1" = "Beta",
                                "2" = "Beta", 
                                "3" = "Alpha",
                                "4" = "Alpha", 
                                "5" = "Delta",
                                "6" = "Alpha", 
                                "7" = "Activated-Stellate",
                                "8" = "Ductal", 
                                "9" = "Acinar",
                                "10" = "Alpha", 
                                "11" = "Beta",
                                "12" = "Gamma",
                                "13" = "Alpha",
                                "14" = "Endothelial",
                                "15" = "Transdifferentiating-Endo",
                                "16" = "Transdifferentiating-Endo",
                                "17" = "Quiescent-Stellate",
                                "18" = "Transdifferentiating-Exo",
                                "19" = "Macrophage",
                                "20" = "Endothelial",
                                "21" = "Alpha",
                                "22" = "Activated-Stellate",
                                "23" = "Alpha",
                                "24" = "T-cell"
)

# Saving this information in the metadata slot
table(Idents(integrated_atac))
integrated_atac$celltype <- Idents(integrated_atac)
summary(integrated_atac@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("Beta", "Transdifferentiating-Endo", "Alpha", "Delta", "Gamma",
               "Ductal", "Transdifferentiating-Exo", "Acinar", 
               "Quiescent-Stellate", "Activated-Stellate",
               "Macrophages", "T-cells",
               "Endothelial")
head(integrated_atac@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
integrated_atac@meta.data$celltype <- factor(x = integrated_atac@meta.data$celltype, levels = my_levels)
Idents(pancreas.combined) <- "celltype"
DimPlot(object = integrated_atac, label = TRUE) + NoLegend()

# Save file
# Save env
# save.image("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Envs/ATAC_integrated_correctGenome.RData")
# saveRDS(integration.anchors, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\integration.anchors.rds)")
# saveRDS(integrated_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\integrated_atac.rds)")
# saveRDS(combined_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\combined_atac.rds)")
# saveRDS(HP2022801_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2022801_atac.rds)")
# saveRDS(SAMN15877725_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\SAMN15877725_atac.rds)")
# saveRDS(HP2024001_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2024001_atac.rds)")
# saveRDS(HP2031401_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2031401_atac.rds)")
# saveRDS(HP2105501_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2105501_atac.rds)")
# saveRDS(HP2106201_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2106201_atac.rds)")
# saveRDS(HP2107001_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2107001_atac.rds)")
# saveRDS(HP2107901_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2107901_atac.rds)")
# saveRDS(HP2108601_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2108601_atac.rds)")
# saveRDS(HP2108901_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2108901_atac.rds)")
# saveRDS(HP2110001_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2110001_atac.rds)")
# saveRDS(HP2121601_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2121601_atac.rds)")
# saveRDS(HP2123201_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2123201_atac.rds)")
# saveRDS(HP2132801_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2132801_atac.rds)")
# saveRDS(HP2202101_atac, file = r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\HP2202101_atac.rds)")

# Load RDS files
HP2022801_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2022801_atac.rds)")
SAMN15877725_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\SAMN15877725_atac.rds)")
HP2024001_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2024001_atac.rds)")
HP2031401_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2031401_atac.rds)")
HP2105501_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2105501_atac.rds)")
HP2106201_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2106201_atac.rds)")
HP2107001_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2107001_atac.rds)")
HP2107901_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2107901_atac.rds)")
HP2108601_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2108601_atac.rds)")
HP2108901_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2108901_atac.rds)")
HP2110001_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2110001_atac.rds)")
HP2121601_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2121601_atac.rds)")
HP2123201_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2123201_atac.rds)")
HP2132801_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2132801_atac.rds)")
HP2202101_atac <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\HP2202101_atac.rds)")

combined_atac <- readRDS(r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\combined_atac.rds)") 
integration.anchors <- readRDS(r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\integration.anchors.rds)")
integrated_atac <- readRDS(r"(D:\3. Coding Scripts and Data\Sex Based study\RDS files\integrated_atac.rds)")
pancreas.combined <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.combinedcorrectedSCT.rds)") # RNA dataset integrated


# ###############








































# Reference Mapping
atac_test <- integrated_atac
rna_test <- pancreas.combined

# map query onto the reference dataset
# find transfer anchors
gene.activities <- GeneActivity(atac_test)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
atac_test[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac_test <- NormalizeData(
  object = atac_test,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac_test$nCount_RNA)
)

DefaultAssay(rna_test) <- "RNA"
common_features <- lapply(list(rna_test, atac_test), row.names) %>% Reduce(intersect, .) 
transfer.anchors <- FindTransferAnchors(
  reference = rna_test,
  query = atac_test,
  reference.reduction = "umap",
  reduction = "lsiproject",
  dims = 2:30
)

multmodal.atac <- MapQuery(
  anchorset = transfer.anchors,
  reference = pbmc.multi,
  query = pbmc.atac,
  refdata = pbmc.multi$predicted.id,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)









# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(integrated_atac) <- annotations

# QC
# compute nucleosome signal score per cell
integrated_atac <- NucleosomeSignal(object = integrated_atac)

# compute TSS enrichment score per cell
integrated_atac <- TSSEnrichment(object = integrated_atac, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
integrated_atac$pct_reads_in_peaks <- integrated_atac$peak_region_fragments / integrated_atac$passed_filters * 100
integrated_atac$blacklist_ratio <- integrated_atac$blacklist_region_fragments / integrated_atac$peak_region_fragments

# Check QC
VlnPlot(
  object = integrated_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

