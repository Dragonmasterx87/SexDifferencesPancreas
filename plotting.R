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

############################ STAGE ############################
############################   7   ############################
# plotting first load seurat object
# First Plot cell Based clustering
processed_rna <- qread(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\processed_rna.qs)")
processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")

# Add metadata
Idents(processed_rna) <- "celltype_sex_ancestry_disease"
processed_rna$celltype_sex_ancestry_disease_lib <- paste(Idents(processed_rna), processed_rna$'Library', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib)

Idents(processed_rna) <- "celltype_sex_ancestry_disease_lib"
processed_rna$celltype_sex_ancestry_disease_lib_source <- paste(Idents(processed_rna), processed_rna$'Tissue Source', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib_source)

Idents(processed_rna) <- "Diabetes Status"
processed_rna$disease_ancestry_lib_sex <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', sep = "_")
table(processed_rna$disease_ancestry_lib_sex)

Idents(processed_rna) <- "Diabetes Status"
processed_rna$disease_ancestry_lib_sex_source <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'Tissue Source', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source)

Idents(processed_rna) <- "Diabetes Status"
processed_rna$disease_ancestry_lib_sex_source_celltype <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'Tissue Source', processed_rna$'celltype_qadir', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source_celltype)

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M", "ND_black_HP2110001_M", "ND_black_HP2123201_M", "ND_black_HPAP-052_M", #Black M ND
               "ND_black_HP2106201_F", "ND_black_HP2121601_F", "ND_black_HP2132801_F", "ND_black_HP2202101_F", #Black F ND
               
               "ND_hispanic_HPAP-080_M", #Hispanic M ND
               "ND_hispanic_HPAP-099_F", "ND_hispanic_HPAP-101_F", "ND_hispanic_HPAP-105_F", #Hispanic F ND
               
               "ND_white_HP2107001_M", "ND_white_HP2107901_M", "ND_white_HPAP-026_M", "ND_white_HPAP-035_M", "ND_white_HPAP-040_M", "ND_white_HPAP-056_M", "ND_white_HPAP-059_M", "ND_white_HPAP-075_M", "ND_white_HPAP-077_M", "ND_white_HPAP-082_M", "ND_white_SAMN15877725_M", #White M ND
               "ND_white_HP2022801_F", "ND_white_HP2024001_F", "ND_white_HP2105501_F", "ND_white_HP2108601_F", "ND_white_HP2108901_F", "ND_white_HPAP-022_F", "ND_white_HPAP-036_F", "ND_white_HPAP-037_F", "ND_white_HPAP-053_F", "ND_white_HPAP-054_F", "ND_white_HPAP-063_F",  "ND_white_HPAP-074_F", "ND_white_HPAP-103_F", #White F ND  
               
               "T2D_black_HPAP-070_M", "T2D_black_HPAP-083_M", "T2D_black_HPAP-108_M", #Black M T2D  
               "T2D_black_HPAP-051_F", "T2D_black_HPAP-058_F", "T2D_black_HPAP-061_F", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F", "T2D_hispanic_HPAP-091_F", "T2D_hispanic_HPAP-109_F", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M", "T2D_white_HPAP-100_M", "T2D_white_HPAP-106_M", # White M T2D
               "T2D_white_HPAP-057_F", "T2D_white_HPAP-081_F", "T2D_white_HPAP-085_F") # White F T2D

table(processed_rna$disease_ancestry_lib_sex)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex <- factor(x = processed_rna$disease_ancestry_lib_sex, levels = my_levels)
table(processed_rna$disease_ancestry_lib_sex)

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M_Tulane", "ND_black_HP2110001_M_Tulane", "ND_black_HP2123201_M_Tulane", "ND_black_HPAP-052_M_UPenn", #Black M ND
               "ND_black_HP2106201_F_Tulane", "ND_black_HP2121601_F_Tulane", "ND_black_HP2132801_F_Tulane", "ND_black_HP2202101_F_Tulane", #Black F ND
               
               "ND_hispanic_HPAP-080_M_nPod", #Hispanic M ND
               "ND_hispanic_HPAP-099_F_UPenn", "ND_hispanic_HPAP-101_F_nPod", "ND_hispanic_HPAP-105_F_nPod", #Hispanic F ND
               
               "ND_white_HP2107001_M_Tulane", "ND_white_HP2107901_M_Tulane", "ND_white_HPAP-026_M_nPod", "ND_white_HPAP-035_M_UPenn", "ND_white_HPAP-040_M_UPenn", "ND_white_HPAP-056_M_UPenn", "ND_white_HPAP-059_M_UPenn", "ND_white_HPAP-075_M_UPenn", "ND_white_HPAP-077_M_UPenn", "ND_white_HPAP-082_M_nPod", "ND_white_SAMN15877725_M_Tulane", #White M ND
               "ND_white_HP2022801_F_Tulane", "ND_white_HP2024001_F_Tulane", "ND_white_HP2105501_F_Tulane", "ND_white_HP2108601_F_Tulane", "ND_white_HP2108901_F_Tulane", "ND_white_HPAP-022_F_UPenn", "ND_white_HPAP-036_F_nPod", "ND_white_HPAP-037_F_UPenn", "ND_white_HPAP-053_F_UPenn", "ND_white_HPAP-054_F_UPenn", "ND_white_HPAP-063_F_UPenn",  "ND_white_HPAP-074_F_UPenn", "ND_white_HPAP-103_F_UPenn", #White F ND  
               
               "T2D_black_HPAP-070_M_nPod", "T2D_black_HPAP-083_M_UPenn", "T2D_black_HPAP-108_M_nPod", #Black M T2D  
               "T2D_black_HPAP-051_F_UPenn", "T2D_black_HPAP-058_F_nPod", "T2D_black_HPAP-061_F_UPenn", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F_nPod", "T2D_hispanic_HPAP-091_F_nPod", "T2D_hispanic_HPAP-109_F_nPod", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M_nPod", "T2D_white_HPAP-100_M_nPod", "T2D_white_HPAP-106_M_UPenn", # White M T2D
               "T2D_white_HPAP-057_F_UPenn", "T2D_white_HPAP-081_F_nPod", "T2D_white_HPAP-085_F_UPenn") # White F T2D

table(processed_rna$disease_ancestry_lib_sex_source)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex_source <- factor(x = processed_rna$disease_ancestry_lib_sex_source, levels = my_levels)
table(processed_rna$disease_ancestry_lib_sex_source)

# Tulane
Idents(processed_rna) <- 'Tissue Source'
tulane_rna <- subset(processed_rna, idents = "Tulane")
hpap_rna <- subset(processed_rna, idents = c("nPod", "UPenn"))


DimPlot(processed_rna, #switch here to plot
        #split.by = "Diabetes Status", 
        group.by = "celltype_qadir", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.05,
        cols = c("dodgerblue3",      #beta
                 "turquoise2",       #beta+alpha
                 "lightseagreen",    #alpha
                 "darkseagreen2",    #cycling_endo
                 "khaki2",           #epsilon 
                 "springgreen4",     #gamma
                 "chartreuse3",      #delta
                 "burlywood3",       #beta+delta
                 "darkorange2",      #ductal
                 "salmon3",          #acinar
                 "orange",           #activated_setallate
                 "salmon",           #quiescent_stellate
                 "red",              #endothelial
                 "magenta3",         #macrophages
                 "orchid1",          #lymphocytes
                 "red4",             #mast
                 "grey30"            #schwann
        )
)


dittoBarPlot(processed_rna, "celltype_qadir", 
             retain.factor.levels = TRUE,
             #scale = "count",
             color.panel = c("dodgerblue3",      #beta
                             "turquoise2",       #beta+alpha
                             "lightseagreen",    #alpha
                             "darkseagreen2",    #cycling_endo
                             "khaki2",           #epsilon 
                             "springgreen4",     #gamma
                             "chartreuse3",      #delta
                             "burlywood3",       #beta+delta
                             "darkorange2",      #ductal
                             "salmon3",          #acinar
                             "orange",           #activated_setallate
                             "salmon",           #quiescent_stellate
                             "red",              #endothelial
                             "magenta3",         #macrophages
                             "orchid1",          #lymphocytes
                             "red4",             #mast
                             "grey30"),          #schwann), group.by = "tx")
                             group.by = "disease_ancestry_lib_sex_source") + coord_flip()


# Umap of Diabetes status
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", 
        group.by = "Diabetes Status", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("dodgerblue", #ND
                 "red2"         #T2D  
        ))
# Umap of sex
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", http://127.0.0.1:42565/graphics/plot_zoom_png?width=1160&height=900
        group.by = "Sex", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("red4", #ND
                 "deepskyblue3"         #T2D  
        ))
# Umap of ancestry_sex
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", 
        group.by = "ancestry", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("black", #ND
                 "darkorange",
                 "deepskyblue3"
                 ))
# Umap of tissue source
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", 
        group.by = "Tissue Source", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("dodgerblue",
                 "springgreen4",         
                 "red4"
        )
)

# Heatmap
# Make average seurat object
Idents(processed_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AverageExpression(processed_rna, return.seurat = TRUE, slot = 'data')

# subset Tulane
Idents(processed_rna) <- "Tissue Source"
tulane_rna <- subset(processed_rna, idents = "Tulane")
Idents(tulane_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AverageExpression(tulane_rna, return.seurat = TRUE, slot = 'data')

# Split Metadata and add columns

{
  combined_processed_rna$disease_ancestry_lib_sex_source_celltype <- combined_processed_rna@active.ident
  Idents(combined_processed_rna) <- 'disease_ancestry_lib_sex_source_celltype'
  combined_processed_rna$disease <- combined_processed_rna$orig.ident
  metadat <- combined_processed_rna@meta.data
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "quiescent_stellate", "quiescent-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "cycling_endo", "cycling-endo"))
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, "_", -5)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -4)
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -3)
  metadat$source <- metadat[c('source')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -2)
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -1)
  combined_processed_rna@meta.data = metadat
}

genes.to.plot <- c('INS', 'GCG')

genes.to.plot <- c("DDX3Y", "EIF1AY",
                   "KDM5D", "NLGN4Y",
                   "RPS4Y1","USP9Y", 
                   "UTY", "ZFY",
                   "XIST", "TSIX",
                   "ZFX", "KDM5C",
                   "SEPTIN6", "EIF1AX",
                   "KDM6A", "PUDP", "DDX3X")

dittoHeatmap(
  combined_processed_rna,
  genes = genes.to.plot,
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("ancestry", "sex", "source", "celltype", "disease"),
  #annot.by = c("lib", "sex", "source"),
  order.by = c("sex", "ancestry", "disease"),
  # main = NA,
  # cell.names.meta = NULL,
  # assay = .default_assay(object),
  # slot = .default_slot(object),
  # swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
  scaled.to.max = FALSE,
  # heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  # annot.colors = c(dittoColors(), dittoColors(1)[seq_len(7)]),
  # annotation_col = NULL,
  annotation_colors = list(celltype = c("acinar" = "salmon3",
                                        "activated-stellate" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "beta+alpha" = "turquoise2",
                                        "beta+delta" = "burlywood3",
                                        "cycling-endo" = "darkseagreen2",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "epsilon" = "khaki2",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophages" = "magenta3",
                                        "mast" = "red4",
                                        "quiescent-stellate" = "salmon",
                                        "schwann" = "grey30"),
                           disease = c("ND" = "dodgerblue",
                                       "T2D" = "red2"),
                           sex = c("F" = "red4",
                                   "M" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black",
                                        "hispanic" = "darkorange"),
                           source = c("nPod" = "dodgerblue",
                                      "Tulane" = "springgreen4",         
                                      "UPenn" = "red4")),
  # # data.out = FALSE,
  # highlight.features = NULL,
  # highlight.genes = NULL,
  # show_colnames = isBulk(object),
  # show_rownames = TRUE,
  # scale = "row",
  #cluster_cols = TRUE,
  # border_color = NA,
  # legend_breaks = NA,
  # drop_levels = FALSE,
  # breaks = NA,
  # complex = FALSE
)

# Violin plot
Idents(processed_rna) <- "ancestry_sex"
processed_rna$ancestry_sex_library <- paste(Idents(processed_rna), processed_rna$'Library', sep = "_")
table(processed_rna@meta.data[["ancestry_sex_library"]])

Idents(processed_rna) <- "ancestry_sex_library"
Stacked_VlnPlot(processed_rna, features = c("DDX3Y", "EIF1AY",
                                            "KDM5D", "NLGN4Y",
                                            "RPS4Y1","USP9Y", 
                                            "UTY", "ZFY"), x_lab_rotate = TRUE)

Stacked_VlnPlot(processed_rna, features = c("XIST", "TSIX",
                                            "ZFX", "KDM5C",
                                            "SEPTIN6", "EIF1AX",
                                            "KDM6A", "PUDP", "DDX3X"), x_lab_rotate = TRUE)

# QC
processed_rna[["percent.mt"]] <- PercentageFeatureSet(processed_rna, pattern = "^MT-")
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
VlnPlot(processed_rna, group.by = "Tissue Source", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
Idents(hpap_rna) <- "Tissue Source"
plot1 <- FeatureScatter(hpap_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hpap_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", plot.cor = TRUE)
plot1 + plot2

# Volcano plot
# Plotting DE genes ###
Idents(processed_rna) <- "celltype_qadir"
beta.cells <- subset(processed_rna, idents = "beta")
plots <- VlnPlot(beta.cells, features = c("INS", "DDIT3", "MIF", "DEPP1", "PLCG2", "IAPP"), group.by = "Sex", 
                 pt.size = 0, combine = TRUE)
plots <- VlnPlot(beta.cells, features = c("MT-CO3", "MT-ND1", "MT-ND4", "MT-ATP6", "MT-CO1", "MT-CYB"), group.by = "Sex", 
                 pt.size = 0, combine = TRUE)
wrap_plots(plots = plots, nrow = 1, ncol = 1)

# Load data
volcanodat <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\analysis\DE_testing\old_data\tulane_samples\no_cookscutoff\beta\beta.deseq.WaldTest.M_ND.vs.F_ND.tsv)",
                         header = TRUE, sep = '\t', row.names = 1)
# volcanodat <- volcanodat %>% 
#   mutate_at(c('pvalue'), ~replace_na(.,0.0000000000000001))
#read.table(file.path(x), sep = '\t', row.names = 1) 
#volcanodat <- epsilon.DHT.response

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1.1
keyvals[which(volcanodat$log2FoldChange > 0.263034 & volcanodat$pvalue < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$log2FoldChange > 0.263034 & volcanodat$pvalue < 0.05)] <- 'high'
  
# modify keyvals for variables with fold change < -1.1
keyvals[which(volcanodat$log2FoldChange < -0.32193 & volcanodat$pvalue < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$log2FoldChange < -0.32193 & volcanodat$pvalue < 0.05)] <- 'low'
    
unique(names(keyvals))
    
unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = FALSE,
                selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                #selectLab = c('XIST', 'PNLIP', 'TSIX', 'REG1B', 'SLC6A2', 'FOXC1',
                #              'ADGRG7', 'RPS4Y1', 'EIF1AY', 'USP9Y', 'UTY', 'KDM5D'), # use this for labelling genes on plot
                # encircle = c('XIST', 'PNLIP', 'TSIX', 'REG1B', 'SLC6A2', 'FOXC1',
                #              'ADGRG7', 'RPS4Y1', 'EIF1AY', 'USP9Y', 'UTY', 'KDM5D'),
                #boxedLabels = TRUE,
                xlim = c(-15,15),
                ylim = c(0,30),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = c(--0.32192, 0.26303), 
                pointSize = c(ifelse((volcanodat$log2FoldChange > 0.26303 & volcanodat$pvalue < 0.05) | (volcanodat$log2FoldChange < -0.32192 & volcanodat$pvalue < 0.05), 5, 2)),
                #pointSize = 2,
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
                drawConnectors = FALSE,
                arrowheads = FALSE,
                widthConnectors = 1,
                typeConnectors = 'closed', 
                max.overlaps = 200,
                #min.segment.length = 20,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black') + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))

#Gene Ontology plotting
# Load data
# Tulane UP
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Tulane DOWN
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Alldat UP
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Alldat DOWN
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Hpap UP
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Hpap DOWN
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Convert frac to dec
beta.F_white_ND.vs.F_black_ND <- beta.F_white_ND.vs.F_black_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_black_ND.vs.F_black_ND <- beta.M_black_ND.vs.F_black_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_ND.vs.F_ND <- beta.M_ND.vs.F_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_white_ND.vs.F_white_ND <- beta.M_white_ND.vs.F_white_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_white_ND.vs.M_black_ND <- beta.M_white_ND.vs.M_black_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

# Add col identifier
beta.F_white_ND.vs.F_black_ND["condition"] = 'beta.F_white_ND.vs.F_black_ND'
beta.M_black_ND.vs.F_black_ND["condition"] = 'beta.M_black_ND.vs.F_black_ND'
beta.M_ND.vs.F_ND["condition"] = 'beta.M_ND.vs.F_ND'
beta.M_white_ND.vs.F_white_ND["condition"] = 'beta.M_white_ND.vs.F_white_ND'
beta.M_white_ND.vs.M_black_ND["condition"] = 'beta.M_white_ND.vs.M_black_ND'

# Merge
merged.df <- do.call("rbind", list(beta.F_white_ND.vs.F_black_ND, beta.M_black_ND.vs.F_black_ND, beta.M_ND.vs.F_ND, beta.M_white_ND.vs.F_white_ND, beta.M_white_ND.vs.M_black_ND))
table(merged.df$condition)

# data to plot
data_mod <- merged.df %>%                                     
  arrange(desc(qvalue)) %>%
  group_by(condition) %>%
  slice(1:5)
print(data_mod)


ggplot(data = merged.df,
       aes(x = condition, y = Description, 
       color = `qvalue`, size = GeneRatio)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=8, face = "bold"), 
        legend.text=element_text(size=8, face = "bold")) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

#Gene Ontology plotting
# Load data
# Tulane UP
beta.deseq.WaldTest.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\tulane\beta.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
beta.deseq.WaldTest.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\tulane\beta.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
beta.deseq.WaldTest.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\tulane\beta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
beta.deseq.WaldTest.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\tulane\beta.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
beta.deseq.WaldTest.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\tulane\beta.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)

# Extract gene lists UP
fwvsfb <- dplyr::filter(beta.deseq.WaldTest.F_white_ND.vs.F_black_ND, padj < 0.1 & log2FoldChange > 0.000000000014)
mbvsfb <- dplyr::filter(beta.deseq.WaldTest.M_black_ND.vs.F_black_ND, padj < 0.1 & log2FoldChange > 0.000000000014)
mvsf <- dplyr::filter(beta.deseq.WaldTest.M_ND.vs.F_ND, padj < 0.1 & log2FoldChange > 0.000000000014)
mwvsfw <- dplyr::filter(beta.deseq.WaldTest.M_white_ND.vs.F_white_ND, padj < 0.1 & log2FoldChange > 0.000000000014)
mwvsmb <- dplyr::filter(beta.deseq.WaldTest.M_white_ND.vs.M_black_ND, padj < 0.1 & log2FoldChange > 0.000000000014)

# Extract gene lists DOWN
fwvsfb <- dplyr::filter(beta.deseq.WaldTest.F_white_ND.vs.F_black_ND, log2FoldChange < -0.000000000014)
mbvsfb <- dplyr::filter(beta.deseq.WaldTest.M_black_ND.vs.F_black_ND, log2FoldChange < -0.000000000014)
mvsf <- dplyr::filter(beta.deseq.WaldTest.M_ND.vs.F_ND, log2FoldChange < -0.000000000014)
mwvsfw <- dplyr::filter(beta.deseq.WaldTest.M_white_ND.vs.F_white_ND, log2FoldChange < -0.000000000014)
mwvsmb <- dplyr::filter(beta.deseq.WaldTest.M_white_ND.vs.M_black_ND, log2FoldChange < -0.000000000014)

fwvsfb <- rownames(fwvsfb)
mbvsfb <- rownames(mbvsfb)
mvsf <- rownames(mvsf)
mwvsfw <- rownames(mwvsfw)
mwvsmb <- rownames(mwvsmb)

# Make list
gene.list <- list("fwvsfb" = fwvsfb, "mbvsfb" = mbvsfb, "mvsf" = mvsf, "mwvsfw" = mwvsfw, "mwvsmb" = mwvsmb)
all_genes <- rownames(beta.deseq.WaldTest.F_white_ND.vs.F_black_ND)
# Compare
ck <- compareCluster(geneCluster = gene.list, 
                     fun = enrichGO, 
                     universe = all_genes, 
                     keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                     OrgDb = org.Hs.eg.db, 
                     ont = c("ALL"), 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1, #if not set default is at 0.05
                     readable = TRUE)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
head(ck) 
cluster_summary <- data.frame(ck.sub)
ck.sub <- ck[ck@compareClusterResult[["qvalue"]] < 0.1, asis=T]
dotplot(ck, showCategory = 20)
dotplot(ck.sub, showCategory = 14)
    
    ############################ END ############################
    ############################ END ############################
    
    Idents(processed_rna) <- "ancestry_sex"
    processed_rna$ancestry_sex_diabetes <- paste(Idents(processed_rna), processed_rna$'Diabetes Status', sep = "_")
    table(processed_rna@meta.data[["ancestry_sex_diabetes"]])
    
    Idents(processed_rna) <- "ancestry_sex_diabetes"
    processed_rna$ancestry_sex_diabetes_sample <- paste(Idents(processed_rna), processed_rna$'Library', sep = "_")
    table(processed_rna@meta.data[["ancestry_sex_diabetes_sample"]])
    
    Idents(processed_rna) <- "ancestry_sex_diabetes_sample"
    processed_rna$ancestry_sex_diabetes_sample_origin <- paste(Idents(processed_rna), processed_rna$'Tissue Source', sep = "_")
    table(processed_rna@meta.data[["ancestry_sex_diabetes_sample_origin"]])
    
    
    # Dotplots
    Idents(processed_rna) <- "Diabetes Status"
    ND <- subset(processed_rna, idents = c('ND'))
    T2D <- subset(processed_rna, idents = c('T2D'))
    ND<- NULL
    bulk <- NULL
    genes.to.plot <- c('XIST', 'SRY')
    
    DefaultAssay(ND) <- 'SCT'
    Idents(ND) <- 'celltype_sex'
    table(ND$celltype_sex)
    # Define an order of cluster identities remember after this step-
    # cluster re-assignment occurs, which re-assigns clustering in my_levels
    my_levels <- c("alpha_F", 
                   "beta+alpha_F",
                   "beta_F", 
                   "beta+delta_F",
                   "delta_F",
                   "gamma_F",
                   "cycling_endo_F",
                   "epsilon_F",
                   "ductal_F",
                   "acinar_F",
                   "quiescent_stellate_F",
                   "activated_stellate_F",
                   "schwann_F",
                   "mast_F",
                   "endothelial_F",
                   "lymphocyte_F",
                   "macrophages_F",
                   "alpha_M", 
                   "beta+alpha_M",
                   "beta_M", 
                   "beta+delta_M", 
                   "delta_M", 
                   "gamma_M", 
                   "cycling_endo_M",
                   "epsilon_M",
                   "ductal_M",
                   "acinar_M",
                   "quiescent_stellate_M",
                   "activated_stellate_M",
                   "schwann_M",
                   "mast_M",
                   "endothelial_M",
                   "lymphocyte_M",
                   "macrophages_M")
    head(ND@meta.data$celltype_sex)
    
    # Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
    ND@meta.data$celltype_sex <- factor(x = ND@meta.data$celltype_sex, levels = my_levels)
    Idents(ND) <- "celltype_sex"
    
    # Selected genes
    markers.to.plot <- c("DDX3Y", "EIF1AY",
                         "KDM5D", "NLGN4Y",
                         "RPS4Y1","USP9Y", 
                         "UTY", "ZFY",
                         "XIST", "TSIX",
                         "ZFX", "KDM5C",
                         "SEPTIN6", "EIF1AX",
                         "KDM6A", "PUDP", "DDX3X")
    
    # Dotplot
    DotPlot(ND,  
            dot.scale = 8,
            col.min = -1, #minimum level
            col.max = 1,  #maximum level
            features = rev(markers.to.plot),
            scale = TRUE) + 
      geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
      theme_light() +
      #facet_wrap(~??? what metadata should be here??)
      #coord_flip() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.title=element_text(size=12, face = "bold"), 
            legend.text=element_text(size=12, face = "bold")) +
      scale_colour_gradient2(low =c("dodgerblue"), high =c("red3")) +
      guides(color = guide_colorbar(title = 'Average Expression')) + coord_flip() 
    
    # Dotplot
    DotPlot(ND,  
            dot.scale = 1,
            col.min = -1, #minimum level
            col.max = 1,  #maximum level
            features = rev(x.chrom),
            scale = TRUE) + 
      geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
      theme_light() +
      #facet_wrap(~??? what metadata should be here??)
      #coord_flip() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.title=element_text(size=12, face = "bold"), 
            legend.text=element_text(size=12, face = "bold")) +
      scale_colour_gradient2(low =c("dodgerblue"), high =c("red3")) +
      guides(color = guide_colorbar(title = 'Average Expression')) + coord_flip() 
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
    x.chrom <- c("ABCB7",
                 "ABCD1",
                 "ACE2",
                 "ACOT9",
                 "ACSL4",
                 "ACTRT1",
                 "ADGRG2",
                 "ADGRG4",
                 "AFF2",
                 "AGTR2",
                 "AIFM1",
                 "AKAP4",
                 "AKAP14",
                 "AKAP17A",
                 "ALAS2",
                 "ALG13",
                 "AMELX",
                 "AMER1",
                 "AMMECR1",
                 "AMOT",
                 "ANOS1",
                 "AP1S2",
                 "APEX2",
                 "APLN",
                 "APOO",
                 "APOOL",
                 "AR",
                 "ARAF",
                 "ARHGAP4",
                 "ARHGAP6",
                 "ARHGAP36",
                 "ARHGEF6",
                 "ARHGEF9",
                 "ARL13A",
                 "ARMCX1",
                 "ARMCX2",
                 "ARMCX3",
                 "ARMCX4",
                 "ARMCX5",
                 "ARMCX6",
                 "ARR3",
                 "ARSD",
                 "ARSF",
                 "ARSH",
                 "ARSL",
                 "ARX",
                 "ASB9",
                 "ASB11",
                 "ASB12",
                 "ASMT",
                 "ASMTL",
                 "ATG4A",
                 "ATP1B4",
                 "ATP2B3",
                 "ATP6AP1",
                 "ATP6AP2",
                 "ATP7A",
                 "ATP11C",
                 "ATRX",
                 "ATXN3L",
                 "AVPR2",
                 "AWAT1",
                 "AWAT2",
                 "BCAP31",
                 "BCLAF3",
                 "BCOR",
                 "BCORL1",
                 "BEND2",
                 "BEX1",
                 "BEX2",
                 "BEX3",
                 "BEX4",
                 "BEX5",
                 "BGN",
                 "BMP15",
                 "BMX",
                 "BRCC3",
                 "BRS3",
                 "BRWD3",
                 "BTK",
                 "C1GALT1C1",
                 "CA5B",
                 "CACNA1F",
                 "CAPN6",
                 "CASK",
                 "CBLL2",
                 "CCDC22",
                 "CCDC120",
                 "CCDC160",
                 "CCNB3",
                 "CCNQ",
                 "CD40LG",
                 "CD99",
                 "CD99L2",
                 "CDK16",
                 "CDKL5",
                 "CDX4",
                 "CENPI",
                 "CENPVL1",
                 "CENPVL2",
                 "CENPVL3",
                 "CETN2",
                 "CFAP47",
                 "CFP",
                 "CHIC1",
                 "CHM",
                 "CHRDL1",
                 "CHST7",
                 "CITED1",
                 "CLCN4",
                 "CLCN5",
                 "CLDN2",
                 "CLDN34",
                 "CLIC2",
                 "CLTRN",
                 "CMC4",
                 "CNGA2",
                 "CNKSR2",
                 "COL4A5",
                 "COL4A6",
                 "COX7B",
                 "CPXCR1",
                 "CRLF2",
                 "CSAG1",
                 "CSAG2",
                 "CSAG3",
                 "CSF2RA",
                 "CSTF2",
                 "CT45A1",
                 "CT45A2",
                 "CT45A3",
                 "CT45A5",
                 "CT45A6",
                 "CT45A7",
                 "CT45A8",
                 "CT45A9",
                 "CT45A10",
                 "CT47A1",
                 "CT47A2",
                 "CT47A3",
                 "CT47A4",
                 "CT47A5",
                 "CT47A6",
                 "CT47A7",
                 "CT47A8",
                 "CT47A9",
                 "CT47A10",
                 "CT47A11",
                 "CT47A12",
                 "CT47B1",
                 "CT47C1",
                 "CT55",
                 "CT83",
                 "CTAG1A",
                 "CTAG1B",
                 "CTAG2",
                 "CTPS2",
                 "CUL4B",
                 "CXCR3",
                 "CXorf38",
                 "CXorf49",
                 "CXorf49B",
                 "CXorf51A",
                 "CXorf51B",
                 "CXorf58",
                 "CXorf65",
                 "CXorf66",
                 "CYBB",
                 "CYLC1",
                 "CYSLTR1",
                 "DACH2",
                 "DCAF8L1",
                 "DCAF8L2",
                 "DCAF12L1",
                 "DCAF12L2",
                 "DCX",
                 "DDX3X",
                 "DDX53",
                 "DGAT2L6",
                 "DGKK",
                 "DHRSX",
                 "DIAPH2",
                 "DIPK2B",
                 "DKC1",
                 "DLG3",
                 "DMD",
                 "DMRTC1",
                 "DMRTC1B",
                 "DNAAF6",
                 "DNASE1L1",
                 "DOCK11",
                 "DRP2",
                 "DUSP9",
                 "DUSP21",
                 "DYNLT3",
                 "EBP",
                 "EDA",
                 "EDA2R",
                 "EFHC2",
                 "EFNB1",
                 "EGFL6",
                 "EIF1AX",
                 "EIF2S3",
                 "ELF4",
                 "ELK1",
                 "EMD",
                 "ENOX2",
                 "EOLA1",
                 "EOLA2",
                 "ERAS",
                 "ERCC6L",
                 "ESX1",
                 "ETDA",
                 "ETDB",
                 "ETDC",
                 "EZHIP",
                 "F8",
                 "F8A1",
                 "F8A2",
                 "F8A3",
                 "F9",
                 "FAAH2",
                 "FAM3A",
                 "FAM9A",
                 "FAM9B",
                 "FAM9C",
                 "FAM47A",
                 "FAM47B",
                 "FAM47C",
                 "FAM50A",
                 "FAM104B",
                 "FAM120C",
                 "FAM133A",
                 "FAM156A",
                 "FAM156B",
                 "FAM199X",
                 "FAM236A",
                 "FAM236B",
                 "FAM236C",
                 "FAM236D",
                 "FANCB",
                 "FATE1",
                 "FGD1",
                 "FGF13",
                 "FGF16",
                 "FHL1",
                 "FLNA",
                 "FMR1",
                 "FMR1NB",
                 "FOXO4",
                 "FOXP3",
                 "FOXR2",
                 "FRMD7",
                 "FRMPD3",
                 "FRMPD4",
                 "FTHL17",
                 "FTSJ1",
                 "FUNDC1",
                 "FUNDC2",
                 "G6PD",
                 "GAB3",
                 "GABRA3",
                 "GABRE",
                 "GABRQ",
                 "GAGE1",
                 "GAGE2A",
                 "GAGE2B",
                 "GAGE2C",
                 "GAGE2D",
                 "GAGE2E",
                 "GAGE4",
                 "GAGE5",
                 "GAGE6",
                 "GAGE7",
                 "GAGE8",
                 "GAGE10",
                 "GAGE12B",
                 "GAGE12C",
                 "GAGE12D",
                 "GAGE12E",
                 "GAGE12F",
                 "GAGE12G",
                 "GAGE12H",
                 "GAGE12I",
                 "GAGE12J",
                 "GAGE13",
                 "GATA1",
                 "GCNA",
                 "GDI1",
                 "GDPD2",
                 "GEMIN8",
                 "GJB1",
                 "GK",
                 "GLA",
                 "GLOD5",
                 "GLRA2",
                 "GLUD2",
                 "GNG5B",
                 "GNL3L",
                 "GPC3",
                 "GPC4",
                 "GPKOW",
                 "GPM6B",
                 "GPR34",
                 "GPR50",
                 "GPR82",
                 "GPR101",
                 "GPR119",
                 "GPR143",
                 "GPR173",
                 "GPR174",
                 "GPRASP1",
                 "GPRASP2",
                 "GPRASP3",
                 "GRIA3",
                 "GRIPAP1",
                 "GRPR",
                 "GSPT2",
                 "GTPBP6",
                 "GUCY2F",
                 "GYG2",
                 "H2AB1",
                 "H2AB2",
                 "H2AB3",
                 "H2AL3",
                 "H2AP",
                 "H2BW1",
                 "H2BW2",
                 "HAPSTR2",
                 "HAUS7",
                 "HCCS",
                 "HCFC1",
                 "HDAC6",
                 "HDAC8",
                 "HDX",
                 "HEPH",
                 "HMGB3",
                 "HMGN5",
                 "HNRNPH2",
                 "HPRT1",
                 "HS6ST2",
                 "HSD17B10",
                 "HSFX1",
                 "HSFX2",
                 "HSFX3",
                 "HSFX4",
                 "HTATSF1",
                 "HTR2C",
                 "HUWE1",
                 "IDH3G",
                 "IDS",
                 "IGBP1",
                 "IGSF1",
                 "IKBKG",
                 "IL1RAPL1",
                 "IL1RAPL2",
                 "IL2RG",
                 "IL3RA",
                 "IL9R",
                 "IL13RA1",
                 "IL13RA2",
                 "INTS6L",
                 "IQSEC2",
                 "IRAK1",
                 "IRS4",
                 "ITGB1BP2",
                 "ITIH6",
                 "ITM2A",
                 "JADE3",
                 "KANTR",
                 "KCND1",
                 "KCNE5",
                 "KDM5C",
                 "KDM6A",
                 "KIAA1210",
                 "KIF4A",
                 "KLF8",
                 "KLHL4",
                 "KLHL13",
                 "KLHL15",
                 "KLHL34",
                 "KRBOX4",
                 "L1CAM",
                 "LAGE3",
                 "LAMP2",
                 "LANCL3",
                 "LAS1L",
                 "LDOC1",
                 "LHFPL1",
                 "LONRF3",
                 "LPAR4",
                 "LRCH2",
                 "LUZP4",
                 "MAGEA1",
                 "MAGEA2",
                 "MAGEA2B",
                 "MAGEA3",
                 "MAGEA4",
                 "MAGEA6",
                 "MAGEA8",
                 "MAGEA9",
                 "MAGEA9B",
                 "MAGEA10",
                 "MAGEA11",
                 "MAGEA12",
                 "MAGEB1",
                 "MAGEB2",
                 "MAGEB3",
                 "MAGEB4",
                 "MAGEB5",
                 "MAGEB6",
                 "MAGEB6B",
                 "MAGEB10",
                 "MAGEB16",
                 "MAGEB17",
                 "MAGEB18",
                 "MAGEC1",
                 "MAGEC2",
                 "MAGEC3",
                 "MAGED1",
                 "MAGED2",
                 "MAGED4",
                 "MAGED4B",
                 "MAGEE1",
                 "MAGEE2",
                 "MAGEH1",
                 "MAGIX",
                 "MAGT1",
                 "MAMLD1",
                 "MAOA",
                 "MAOB",
                 "MAP3K15",
                 "MAP7D2",
                 "MAP7D3",
                 "MBNL3",
                 "MBTPS2",
                 "MCF2",
                 "MCTS1",
                 "MECP2",
                 "MED12",
                 "MED14",
                 "MED14OS",
                 "MID1",
                 "MID1IP1",
                 "MID2",
                 "MMGT1",
                 "MORC4",
                 "MORF4L2",
                 "MOSPD1",
                 "MOSPD2",
                 "MPC1L",
                 "MPP1",
                 "MSL3",
                 "MSN",
                 "MTCP1",
                 "MTM1",
                 "MTMR1",
                 "MTMR8",
                 "MXRA5",
                 "NAA10",
                 "NALF2",
                 "NAP1L2",
                 "NAP1L3",
                 "NBDY",
                 "NCBP2L",
                 "NDP",
                 "NDUFA1",
                 "NDUFB11",
                 "NEXMIF",
                 "NHS",
                 "NHSL2",
                 "NKAP",
                 "NKRF",
                 "NLGN3",
                 "NLGN4X",
                 "NLRP2B",
                 "NONO",
                 "NOX1",
                 "NR0B1",
                 "NRK",
                 "NSDHL",
                 "NUDT10",
                 "NUDT11",
                 "NUP62CL",
                 "NXF2",
                 "NXF2B",
                 "NXF3",
                 "NXF5",
                 "NXT2",
                 "NYX",
                 "OCRL",
                 "OFD1",
                 "OGT",
                 "OPHN1",
                 "OPN1LW",
                 "OPN1MW",
                 "OPN1MW2",
                 "OPN1MW3",
                 "OR13H1",
                 "OTC",
                 "OTUD5",
                 "OTUD6A",
                 "P2RY4",
                 "P2RY8",
                 "P2RY10",
                 "PABIR2",
                 "PABIR3",
                 "PABPC1L2A",
                 "PABPC1L2B",
                 "PABPC5",
                 "PAGE1",
                 "PAGE2",
                 "PAGE2B",
                 "PAGE3",
                 "PAGE4",
                 "PAGE5",
                 "PAK3",
                 "PASD1",
                 "PBDC1",
                 "PCDH11X",
                 "PCDH19",
                 "PCSK1N",
                 "PCYT1B",
                 "PDHA1",
                 "PDK3",
                 "PDZD4",
                 "PDZD11",
                 "PFKFB1",
                 "PGAM4",
                 "PGK1",
                 "PGRMC1",
                 "PHEX",
                 "PHF6",
                 "PHF8",
                 "PHKA1",
                 "PHKA2",
                 "PIGA",
                 "PIM2",
                 "PIN4",
                 "PIR",
                 "PJA1",
                 "PLAC1",
                 "PLCXD1",
                 "PLP1",
                 "PLP2",
                 "PLS3",
                 "PLXNA3",
                 "PLXNB3",
                 "PNCK",
                 "PNMA3",
                 "PNMA5",
                 "PNMA6A",
                 "PNMA6E",
                 "PNMA6F",
                 "PNPLA4",
                 "POF1B",
                 "POLA1",
                 "PORCN",
                 "POU3F4",
                 "PPEF1",
                 "PPP1R2C",
                 "PPP1R3F",
                 "PPP2R3B",
                 "PPP4R3C",
                 "PQBP1",
                 "PRAF2",
                 "PRDX4",
                 "PRICKLE3",
                 "PRKX",
                 "PRPS1",
                 "PRPS2",
                 "PRR32",
                 "PRRG1",
                 "PRRG3",
                 "PSMD10",
                 "PTCHD1",
                 "PUDP",
                 "PWWP3B",
                 "PWWP4",
                 "RAB9A",
                 "RAB9B",
                 "RAB33A",
                 "RAB39B",
                 "RAB40A",
                 "RAB40AL",
                 "RAB41",
                 "RADX",
                 "RAI2",
                 "RAP2C",
                 "RBBP7",
                 "RBM3",
                 "RBM10",
                 "RBM41",
                 "RBMX",
                 "RBMX2",
                 "RBMXL3",
                 "RENBP",
                 "REPS2",
                 "RGN",
                 "RHOXF1",
                 "RHOXF2",
                 "RHOXF2B",
                 "RIBC1",
                 "RIPPLY1",
                 "RLIM",
                 "RNF113A",
                 "RNF128",
                 "RP2",
                 "RPA4",
                 "RPGR",
                 "RPL10",
                 "RPL36A",
                 "RPL39",
                 "RPS4X",
                 "RPS6KA3",
                 "RPS6KA6",
                 "RRAGB",
                 "RS1",
                 "RTL3",
                 "RTL4",
                 "RTL5",
                 "RTL8A",
                 "RTL8B",
                 "RTL8C",
                 "RTL9",
                 "S100G",
                 "SAGE1",
                 "SASH3",
                 "SAT1",
                 "SATL1",
                 "SCML1",
                 "SCML2",
                 "SEPTIN6",
                 "SERPINA7",
                 "SERTM2",
                 "SH2D1A",
                 "SH3BGRL",
                 "SH3KBP1",
                 "SHOX",
                 "SHROOM2",
                 "SHROOM4",
                 "SLC6A8",
                 "SLC6A14",
                 "SLC7A3",
                 "SLC9A6",
                 "SLC9A7",
                 "SLC10A3",
                 "SLC16A2",
                 "SLC25A5",
                 "SLC25A6",
                 "SLC25A14",
                 "SLC25A43",
                 "SLC25A53",
                 "SLC35A2",
                 "SLC38A5",
                 "SLITRK2",
                 "SLITRK4",
                 "SMARCA1",
                 "SMC1A",
                 "SMIM9",
                 "SMIM10",
                 "SMIM10L2A",
                 "SMIM10L2B",
                 "SMPX",
                 "SMS",
                 "SNX12",
                 "SOWAHD",
                 "SOX3",
                 "SPACA5",
                 "SPACA5B",
                 "SPANXA1",
                 "SPANXA2",
                 "SPANXB1",
                 "SPANXC",
                 "SPANXD",
                 "SPANXN1",
                 "SPANXN2",
                 "SPANXN3",
                 "SPANXN4",
                 "SPANXN5",
                 "SPIN2A",
                 "SPIN2B",
                 "SPIN3",
                 "SPIN4",
                 "SPRY3",
                 "SRPK3",
                 "SRPX",
                 "SRPX2",
                 "SSR4",
                 "SSX1",
                 "SSX2",
                 "SSX2B",
                 "SSX3",
                 "SSX4",
                 "SSX4B",
                 "SSX5",
                 "SSX7",
                 "STAG2",
                 "STARD8",
                 "STEEP1",
                 "STK26",
                 "STS",
                 "SUPT20HL1",
                 "SUPT20HL2",
                 "SUV39H1",
                 "SYAP1",
                 "SYN1",
                 "SYP",
                 "SYTL4",
                 "SYTL5",
                 "TAB3",
                 "TAF1",
                 "TAF7L",
                 "TAF9B",
                 "TAFAZZIN",
                 "TASL",
                 "TBC1D8B",
                 "TBC1D25",
                 "TBL1X",
                 "TBX22",
                 "TCEAL1",
                 "TCEAL2",
                 "TCEAL3",
                 "TCEAL4",
                 "TCEAL5",
                 "TCEAL6",
                 "TCEAL7",
                 "TCEAL8",
                 "TCEAL9",
                 "TCEANC",
                 "TCP11X1",
                 "TCP11X2",
                 "TENM1",
                 "TENT5D",
                 "TEX11",
                 "TEX13A",
                 "TEX13B",
                 "TEX13C",
                 "TEX13D",
                 "TEX28",
                 "TFDP3",
                 "TFE3",
                 "TGIF2LX",
                 "THOC2",
                 "TIMM8A",
                 "TIMM17B",
                 "TIMP1",
                 "TKTL1",
                 "TLR7",
                 "TLR8",
                 "TMEM31",
                 "TMEM35A",
                 "TMEM47",
                 "TMEM164",
                 "TMEM185A",
                 "TMEM187",
                 "TMEM255A",
                 "TMLHE",
                 "TMSB4X",
                 "TMSB15A",
                 "TMSB15B",
                 "TMSB15C",
                 "TNMD",
                 "TRAPPC2",
                 "TREX2",
                 "TRMT2B",
                 "TRO",
                 "TRPC5",
                 "TRPC5OS",
                 "TSC22D3",
                 "TSPAN6",
                 "TSPAN7",
                 "TSPYL2",
                 "TSR2",
                 "TXLNG",
                 "UBA1",
                 "UBE2A",
                 "UBE2NL",
                 "UBL4A",
                 "UBQLN2",
                 "UPF3B",
                 "UPRT",
                 "USP9X",
                 "USP11",
                 "USP26",
                 "USP27X",
                 "USP51",
                 "UTP14A",
                 "UXT",
                 "VAMP7",
                 "VBP1",
                 "VCX",
                 "VCX2",
                 "VCX3A",
                 "VCX3B",
                 "VEGFD",
                 "VGLL1",
                 "VMA21",
                 "VSIG1",
                 "VSIG4",
                 "WAS",
                 "WDR13",
                 "WDR44",
                 "WDR45",
                 "WNK3",
                 "WWC3",
                 "XAGE1A",
                 "XAGE1B",
                 "XAGE2",
                 "XAGE3",
                 "XAGE5",
                 "XG",
                 "XIAP",
                 "XK",
                 "XKRX",
                 "XPNPEP2",
                 "YIPF6",
                 "YY2",
                 "ZBED1",
                 "ZBTB33",
                 "ZC3H12B",
                 "ZC4H2",
                 "ZCCHC12",
                 "ZCCHC13",
                 "ZCCHC18",
                 "ZDHHC9",
                 "ZDHHC15",
                 "ZFP92",
                 "ZFX",
                 "ZIC3",
                 "ZMAT1",
                 "ZMYM3",
                 "ZNF41",
                 "ZNF75D",
                 "ZNF81",
                 "ZNF157",
                 "ZNF182",
                 "ZNF185",
                 "ZNF275",
                 "ZNF280C",
                 "ZNF449",
                 "ZNF630",
                 "ZNF674",
                 "ZNF711",
                 "ZRSR2",
                 "ZXDA",
                 "ZXDB")
    
    BiocManager::install("scuttle", force = TRUE)
    library(scuttle)
    remotes::install_github("omnideconv/SimBu")
    library(SimBu)
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
