# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

# LOAD LIBRARIES ####
# Restart Rstudio or R
# make sure a Rtools version is installed that supports your current versin of R
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
install.packages("rlang")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'), force = TRUE)

BiocManager::install(c("EnhancedVolcano", "DoubletFinder", "glmGamPoi",
                       "GOSemSim", "org.Hs.eg.db", "AnnotationHub",
                       "GenomeInfoDb", "MeSHDbi", "clusterProfiler",
                       "dittoSeq", "escape", "ComplexHeatmap", "DropletUtils", 
                       "Nebulosa", "hdf5r", "scDblFinder", "JASPAR2020",
                       "TFBSTools", "motifmatchr", "GreenleafLab/chromVAR",
                       "EnrichmentBrowser"),
                     dependencies = T, force = TRUE)
BiocManager::install("BiocParallel")

# install Seurat from Github (automatically updates sctransform)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
install.packages("devtools")
devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("gaospecial/ggVennDiagram")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("harmony")
install.packages('SoupX')
install.packages('tidyverse')
install.packages("viridis")
install.packages("circlize")
install.packages("scCustomize")
install.packages("archive")
install.packages("R.utils")
install.packages("qs")
install.packages("ggseqlogo")

#devtools::install_github("stuart-lab/signac", ref = "develop", force = TRUE)

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')

# Run the following code once you have Seurat installed
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(leiden)
  library(stringr)
  library(hdf5r)
  library(SoupX)
  library(Rcpp)
  library(cowplot)
  library(Matrix)
  library(ggridges)
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
  library(BiocParallel)
  library(ggseqlogo)
  # Set global environment parameter par-proc
  # options(future.globals.maxSize = 8000 * 1024^2)
  set.seed(1234)
})
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
packageVersion("monocle3")
packageVersion("cicero")

# Load object
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
combined_atac <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\snATACfiles_earlierpartsofworkflow\combined_atac.qs)")

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma", 
               "acinar", "ductal",
               "quiescent_stellate", "activated_stellate", "endothelial",
               "lymphocyte", "macrophage")

table(hm.integrated.dfree$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree$celltype <- factor(x = hm.integrated.dfree$celltype, levels = my_levels)
table(unique(hm.integrated.dfree$celltype))

#QC
frac_reads_in_promoters <- combined_atac@meta.data$FRiP
quantile(frac_reads_in_promoters)
metadata <- combined_atac@meta.data
ggplot(data=metadata, mapping = aes(x=frac_reads_in_promoters)) +  geom_density(alpha = 0.2, fill= 'lightpink', color="pink") + 
  theme_linedraw() + geom_vline(xintercept=c(0.22,0.2,0.38), colour=c("blue", "red", "purple"),linetype = "longdash")

frac_reads_in_peaks <- combined_atac@meta.data$pct_reads_in_peaks
quantile(frac_reads_in_peaks)
ggplot(data=metadata, mapping = aes(x=frac_reads_in_peaks)) +  geom_density(alpha = 0.2, color="green", fill="lightgreen") + 
  theme_linedraw() + geom_vline(xintercept=c(30), colour=c("red"),linetype = "longdash")

TSS <- combined_atac@meta.data$TSS.enrichment
quantile(TSS)
ggplot(data=metadata, mapping = aes(x=TSS)) +  geom_density(alpha = 0.2, color="green", fill="lightgreen") + 
  theme_linedraw() + geom_vline(xintercept=c(2), colour=c("red"),linetype = "longdash")


DensityScatter(hm.integrated.dfree, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

Idents(hm.integrated.dfree) <- "sex"
Idents(combined_atac) <- "sex"
VlnPlot(
  object = hm.integrated.dfree,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks', "FRiP"),
  pt.size = 0,
  ncol = 3
)

VlnPlot(
  object = combined_atac,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks', "FRiP"),
  pt.size = 0,
  ncol = 3
)

# Dimplot
unique(hm.integrated.dfree@meta.data[["celltype"]])
DimPlot(hm.integrated.dfree, #switch here to plot
        #split.by = "Diabetes Status", 
        group.by = "celltype", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.05,
        cols = c("dodgerblue3",      #beta
                 "chartreuse3",      #delta
                 "lightseagreen",    #alpha
                 "springgreen4",     #gamma
                 "salmon3",          #acinar
                 "darkorange2",      #ductal
                 "salmon",           #quiescent_stellate
                 "orange",           #activated_setallate
                 "red",              #endothelial
                 "orchid1",          #lymphocytes
                 "magenta3"         #macrophages
        )
)

#
hm.integrated.dfree$sex_ancestry_sample <- paste(hm.integrated.dfree$sex, hm.integrated.dfree$ancestry, hm.integrated.dfree$sample, sep = "_")
dittoBarPlot(hm.integrated.dfree, "celltype", 
             retain.factor.levels = TRUE,
             scale = "percent",
             color.panel = c("dodgerblue3",      #beta
                             "chartreuse3",      #delta
                             "lightseagreen",    #alpha
                             "springgreen4",     #gamma
                             "salmon3",          #acinar
                             "darkorange2",      #ductal
                             "salmon",           #quiescent_stellate
                             "orange",           #activated_setallate
                             "red",              #endothelial
                             "orchid1",          #lymphocytes
                             "magenta3"),          
             group.by = "sex_ancestry_sample") + coord_flip()


# Umap of sex
DimPlot(hm.integrated.dfree, 
        #split.by = "Diabetes Status", http://127.0.0.1:42565/graphics/plot_zoom_png?width=1160&height=900
        group.by = "sex", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("red4", #ND
                 "deepskyblue3"         #T2D  
        ))
# Umap of ancestry_sex
DimPlot(hm.integrated.dfree, 
        #split.by = "Diabetes Status", 
        group.by = "ancestry", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("black", #ND
                 "deepskyblue3"
        ))

# Umap of Accessibility
Idents(hm.integrated.dfree) <- "celltype"
DefaultAssay(hm.integrated.dfree) <- "predicted"

FeaturePlot(
  object = hm.integrated.dfree,
  features = c("CD3D"),
  cols = c("lightgrey", "red4"),
  pt.size = 0.1,
  #max.cutoff = 1.6,
  ncol = 1,
  #order = TRUE
)

# #, "GCG", "SST", "PPY", "CFTR", "PRSS1", "ESM1", "SDS", "RGS5", "PDGFRA", "CD7"
# 
# # Corelation plot
# # # Pseudobulk
# # Idents(hm.integrated.dfree) <- "celltype"
# # hm.integrated.dfree$celltype_sex <- paste(hm.integrated.dfree$celltype, hm.integrated.dfree$sex, sep = "_")
# # table(hm.integrated.dfree@meta.data$celltype_sex)
# # Idents(hm.integrated.dfree) <- "celltype_sex"
# # DefaultAssay(hm.integrated.dfree) <- "ATAC"
# # combined_processed_atac <- Seurat:::PseudobulkExpression(object = hm.integrated.dfree, 
# #                                                          pb.method = 'aggregate', 
# #                                                          return.seurat = TRUE,
# #                                                          slot = 'counts')
# # 
# # DefaultAssay(combined_processed_atac) <- "ATAC"
# # combined_processed_atac <- RunTFIDF(combined_processed_atac, assay = "ATAC")
# # 
# # {
# #   combined_processed_atac$celltype_sex <- combined_processed_atac@active.ident
# #   Idents(combined_processed_atac) <- 'celltype_sex'
# #   combined_processed_atac$celltype <- combined_processed_atac$orig.ident
# #   metadat <- combined_processed_atac@meta.data
# #   metadat <- metadat %>% 
# #     mutate(celltype_sex = str_replace(celltype_sex, "activated", "activated-stellate"))
# #   metadat <- metadat %>% 
# #     mutate(celltype_sex = str_replace(celltype_sex, "quiescent", "quiescent-stellate"))
# #   metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex, "_", -1)
# #   combined_processed_atac@meta.data = metadat
# # }
# # 
# # table(combined_processed_atac@meta.data[["celltype"]])
# # table(combined_processed_atac@meta.data[["sex"]])
# # 
# # # cluster re-assignment occurs, which re-assigns clustering in my_levels
# # my_levels <- c("delta", "beta", "alpha", "gamma",
# #                "ductal", "acinar",
# #                "activated", "quiescent", "endothelial",
# #                "lymphocyte", "macrophage") 
# # 
# # table(combined_processed_atac$celltype)
# # 
# # # Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
# # combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
# # table(combined_processed_atac$celltype)
# # Idents(combined_processed_atac) <- "celltype"
# # DefaultAssay(combined_processed_atac) <- "ATAC"
# 
# # split object by cluster, take counts and aggregate by sum (pseudobulking)
# processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
# hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
# 
# # Idents(processed_rna) <- "Sex"
# # unique(processed_rna$Sex)
# # processed_rna$sex <- plyr::mapvalues(
# #   x= processed_rna$Sex,
# #   from = c("F",
# #            "M"),
# #   to = c("female",
# #          "male"))
# # table(processed_rna@meta.data[["Sex"]])
# # table(processed_rna@meta.data[["sex"]])
# # Idents(processed_rna) <- "celltype_qadir"
# # processed_rna$celltype_sex <- paste(processed_rna$celltype_qadir, processed_rna$sex, sep = "_")
# Idents(processed_rna) <- "tissue_source"
# processed_rna <- subset(processed_rna, idents = c("Tulane"))
# Idents(processed_rna) <- "celltype_qadir"
# gc()
# gc()
# 
# #processed_rna <- subset(processed_rna, features = rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features))
# pblist_seu <- lapply(unique(processed_rna$celltype_qadir), function(x) {
#   rowSums(SeuratObject::GetAssayData(processed_rna[,processed_rna$celltype_qadir == x], "data"))
# })
# names(pblist_seu) = unique(processed_rna$celltype_qadir)
# 
# # DefaultAssay(hm.integrated.dfree) <- "RNA"
# # Idents(hm.integrated.dfree) <- "celltype"
# # hm.integrated.dfree$celltype_sex <- paste(hm.integrated.dfree$celltype, hm.integrated.dfree$sex, sep = "_")
# DefaultAssay(hm.integrated.dfree) <- "predicted"
# Idents(hm.integrated.dfree) <- "celltype"
# #hm.integrated.dfree <- subset(hm.integrated.dfree, features = rownames(processed_rna@assays[["RNA"]]@meta.features))
# pblist_sig <- lapply(unique(hm.integrated.dfree$celltype), function(x) {
#   rowSums(SeuratObject::GetAssayData(hm.integrated.dfree[,hm.integrated.dfree$celltype == x], "data"))
# })
# names(pblist_sig) = unique(hm.integrated.dfree$celltype)
# 
# # Subset all lists
# # Create a intersect of gene/peak names
# common_genes <- intersect(intersect(rownames(processed_rna@assays[["RNA"]]@meta.features), rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features)), processed_rna@assays[["RNA"]]@var.features)
# common_genes <- intersect(rownames(processed_rna@assays[["RNA"]]@meta.features), rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features))
# common_genes <- c("INS", "GCG", "SST", "PPY")
# common_genes <- processed_rna@assays[["RNA"]]@var.features
# common_genes <- intersect(rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features), processed_rna@assays[["RNA"]]@var.features)
# 
# # Convert list of matrices to DFs
# ##### ATAC ###
# #####      ###
# pblist_sig <- lapply(
#   X = pblist_sig,
#   FUN = function(x){as.data.frame(x)})
# 
# # Initialize an empty list
# subset_dfs <- list()
# 
# # Iterate over each data frame
# for (i in 1:length(pblist_sig)) {
#   # Create subset based on common genes
#   subset_df <- pblist_sig[[i]][common_genes, ]
# 
#   # Append the subset data frame to the list
#   subset_dfs[[i]] <- subset_df
# }
# 
# # Add names
# names(subset_dfs) = unique(hm.integrated.dfree$celltype)
# pblist_sig <- subset_dfs
# 
# pblist_sig <- lapply(
#   X = pblist_sig,
#   FUN = function(x){as.matrix(x)})
# 
# 
# ##### scRNA ###
# #####       ###
# # Convert list of matrices to DFs
# pblist_seu <- lapply(
#   X = pblist_seu,
#   FUN = function(x){as.data.frame(x)})
# 
# # Initialize an empty list
# subset_dfs <- list()
# 
# # Iterate over each data frame
# for (i in 1:length(pblist_seu)) {
#   # Create subset based on common genes
#   subset_df <- pblist_seu[[i]][common_genes, ]
# 
#   # Append the subset data frame to the list
#   subset_dfs[[i]] <- subset_df
# }
# 
# # Add names
# names(subset_dfs) = unique(processed_rna$celltype_qadir)
# pblist_seu <- subset_dfs
# 
# pblist_seu <- lapply(
#   X = pblist_seu,
#   FUN = function(x){as.matrix(x)})
# 
# # create pairwise combinations
# ids = expand.grid(unique(processed_rna$celltype_qadir), unique(hm.integrated.dfree$celltype))
# 
# # iterate over combinations to extract correlations
# cors = apply(ids, 1, function(x) cor(pblist_seu[[x[1]]], pblist_sig[[x[2]]]))
# 
# pblist_sig[[2]]
# # fold the correlations vector into a matrix
# cmat = matrix(cors, nrow = length(unique(processed_rna$celltype_qadir)), byrow = TRUE)
# 
# # add cluster names. This implies that the names are matching, otherwise you need to add them separately for RNA and ATAC. This in turn depends on the order in which you expanded the grid above.
# rownames(cmat) = unique(processed_rna$celltype_qadir)
# colnames(cmat) = unique(hm.integrated.dfree$celltype)
# 
# # plot however you want
# pheatmap(cmat, clustering_distance_rows = "euclidean", 
#          cluster_rows = FALSE,cluster_cols = FALSE)
# 
# rna.seurat <- processed_rna[, sample(colnames(processed_rna), size =500, replace=F)]
# atac.signac <- hm.integrated.dfree[, sample(colnames(hm.integrated.dfree), size =500, replace=F)]
# qsave(rna.seurat, r"(E:\downsampling_collab\rna.seurat.qs)")
# qsave(atac.signac, r"(E:\downsampling_collab\atac.signac.qs)")



# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma",
               "acinar", "ductal",
               "quiescent_stellate", "activated_stellate", "endothelial",
               "lymphocyte", "macrophage")

table(hm.integrated.dfree$celltype)




# Confusion Matrix
predictions <- table(hm.integrated.dfree$celltype, hm.integrated.dfree$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
head(predictions)
predictions <- as.data.frame(predictions)
# predictions <- subset(predictions, Var2 == c("beta", "delta", "alpha", "gamma", 
#                                              "acinar", "ductal",
#                                              "quiescent_stellate", "activated_stellate", "endothelial",
#                                              "lymphocyte", "macrophage"))
head(predictions)
predicted_scale <- as.numeric(predictions$Freq)
predictions$Freq2 <- scale(predicted_scale, center = TRUE, scale = TRUE)

ggplot(predictions, aes(x=factor(Var1, level=c("delta", "beta", "alpha", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "macrophage", "endothelial", "lymphocyte")), 
                        y=factor(Var2, level=c("delta", "beta+delta", "beta", "beta+alpha", "alpha", "cycling_endo", "gamma", "epsilon", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "macrophages", "endothelial", "lymphocyte", "schwann", "mast")), 
                        fill = Freq2)) + 
  geom_tile() + scale_fill_gradient2(low = "dodgerblue3", high = "darkred", 
                                     midpoint = 0, limit = c(-5,5), space = "Lab", 
                                     name="Scaled Cell Contribution") +
  xlab("Cell type annotation (ATAC)") + ylab("Predicted cell type label (RNA)") +
  theme_minimal()+ 
  geom_text(aes(label = round(Freq, 1))) +
  ggtitle("Normalised Confusion Matrix") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))


# Plotting Chromatin accessibility
# when copying off the heatmap make sure to adjust dashes
DefaultAssay(hm.integrated.dfree) <- "ATAC"
Idents(hm.integrated.dfree) <- "celltype"
Idents(hm.integrated.dfree) <- "sex"
male_set <- subset(hm.integrated.dfree, idents = (c("male")))
female_set <- subset(hm.integrated.dfree, idents = (c("female")))
Idents(male_set) <- "celltype"
Idents(female_set) <- "celltype"
p1 <-  CoveragePlot(male_set, region = c("KDM5D"), 
                    assay = "ATAC",
                    window = 100,
                    ymax = 50,
             links = FALSE,
             #tile = TRUE,
             extend.upstream = 20000,
             extend.downstream = 20000
             ) & scale_fill_manual(values = c("dodgerblue3",      #beta
                                              "chartreuse3",      #delta
                                              "lightseagreen",    #alpha
                                              "springgreen4",     #gamma
                                              "salmon3",          #acinar
                                              "darkorange2",      #ductal
                                              "salmon",           #quiescent_stellate
                                              "orange",           #activated_setallate
                                              "red",              #endothelial
                                              "orchid1",          #lymphocytes
                                              "magenta3"))
p1
p1+p2

link_plot <- LinkPlot(
  object = hm.integrated.dfree,
  region = c("GCG"),
  extend.upstream = 50000,
  extend.downstream = 50000
)

expr_plot <- ExpressionPlot(
  object = hm.integrated.dfree,
  features = "GCG",
  assay = "predicted"
) & scale_fill_manual(values = c("dodgerblue3",      #beta
                                 "chartreuse3",      #delta
                                 "lightseagreen",    #alpha
                                 "springgreen4",     #gamma
                                 "salmon3",          #acinar
                                 "darkorange2",      #ductal
                                 "salmon",           #quiescent_stellate
                                 "orange",           #activated_setallate
                                 "red",              #endothelial
                                 "orchid1",          #lymphocytes
                                 "magenta3"))

CombineTracks(
  plotlist = list(cov_plot, link_plot),
  expression.plot = expr_plot,
  heights = c(10, 2),
  widths = c(10, 1)
)


# Plotting heatmap for all significant DA regions
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

# All genes
beta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\beta_peaks.csv)", sep = ",", row.names = 1)
alpha_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\alpha_peaks.csv)", sep = ",", row.names = 1)
delta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\delta_peaks.csv)", sep = ",", row.names = 1)
gamma_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\gamma_peaks.csv)", sep = ",", row.names = 1)
acinar_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\acinar_peaks.csv)", sep = ",", row.names = 1)
ductal_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\ductal_peaks.csv)", sep = ",", row.names = 1)
quiescentstellate_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\quiescentstellate_peaks.csv)", sep = ",", row.names = 1)
activatedstellate_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\activatedstellate_peaks.csv)", sep = ",", row.names = 1)
macrophage_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\macrophage_peaks.csv)", sep = ",", row.names = 1)
lymphocyte_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\lymphocyte_peaks.csv)", sep = ",", row.names = 1)
endothelial_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\endothelial_peaks.csv)", sep = ",", row.names = 1)
{
  beta_peaks$p_val_adj[beta_peaks$p_val_adj == 0] <- 2e-302
  beta_peaks <- dplyr::filter(beta_peaks, p_val_adj < 1e-5) 
  open_beta_peaks <- rownames(beta_peaks[beta_peaks$avg_log2FC > 1, ])
  
  alpha_peaks$p_val_adj[alpha_peaks$p_val_adj == 0] <- 2e-302
  alpha_peaks <- dplyr::filter(alpha_peaks, p_val_adj < 1e-5) 
  open_alpha_peaks <- rownames(alpha_peaks[alpha_peaks$avg_log2FC > 1, ])
  
  delta_peaks$p_val_adj[delta_peaks$p_val_adj == 0] <- 2e-302
  delta_peaks <- dplyr::filter(delta_peaks, p_val_adj < 1e-5) 
  open_delta_peaks <- rownames(delta_peaks[delta_peaks$avg_log2FC > 1, ])
  
  gamma_peaks$p_val_adj[gamma_peaks$p_val_adj == 0] <- 2e-302
  gamma_peaks <- dplyr::filter(gamma_peaks, p_val_adj < 1e-5) 
  open_gamma_peaks <- rownames(gamma_peaks[gamma_peaks$avg_log2FC > 1, ])
  
  acinar_peaks$p_val_adj[acinar_peaks$p_val_adj == 0] <- 2e-302
  acinar_peaks <- dplyr::filter(acinar_peaks, p_val_adj < 1e-5) 
  open_acinar_peaks <- rownames(acinar_peaks[acinar_peaks$avg_log2FC > 1, ])
  
  ductal_peaks$p_val_adj[ductal_peaks$p_val_adj == 0] <- 2e-302
  ductal_peaks <- dplyr::filter(ductal_peaks, p_val_adj < 1e-5) 
  open_ductal_peaks <- rownames(ductal_peaks[ductal_peaks$avg_log2FC > 1, ])
  
  quiescentstellate_peaks$p_val_adj[quiescentstellate_peaks$p_val_adj == 0] <- 2e-302
  quiescentstellate_peaks <- dplyr::filter(quiescentstellate_peaks, p_val_adj < 1e-5) 
  open_quiescentstellate_peaks <- rownames(quiescentstellate_peaks[quiescentstellate_peaks$avg_log2FC > 1, ])
  
  activatedstellate_peaks$p_val_adj[activatedstellate_peaks$p_val_adj == 0] <- 2e-302
  activatedstellate_peaks <- dplyr::filter(activatedstellate_peaks, p_val_adj < 1e-5) 
  open_activatedstellate_peaks <- rownames(activatedstellate_peaks[activatedstellate_peaks$avg_log2FC > 1, ])
  
  macrophage_peaks$p_val_adj[macrophage_peaks$p_val_adj == 0] <- 2e-302
  macrophage_peaks <- dplyr::filter(macrophage_peaks, p_val_adj < 1e-5) 
  open_macrophage_peaks <- rownames(macrophage_peaks[macrophage_peaks$avg_log2FC > 1, ])
  
  lymphocyte_peaks$p_val_adj[lymphocyte_peaks$p_val_adj == 0] <- 2e-302
  lymphocyte_peaks <- dplyr::filter(lymphocyte_peaks, p_val_adj < 1e-5) 
  open_lymphocyte_peaks <- rownames(lymphocyte_peaks[lymphocyte_peaks$avg_log2FC > 1, ])
  
  endothelial_peaks$p_val_adj[endothelial_peaks$p_val_adj == 0] <- 2e-302
  endothelial_peaks <- dplyr::filter(endothelial_peaks, p_val_adj < 1e-5) 
  open_endothelial_peaks <- rownames(endothelial_peaks[endothelial_peaks$avg_log2FC > 1, ])
  
  cgenes_beta <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks)
  cgenes_alpha <- ClosestFeature(hm.integrated.dfree, regions = open_alpha_peaks)
  cgenes_delta <- ClosestFeature(hm.integrated.dfree, regions = open_delta_peaks)
  cgenes_gamma <- ClosestFeature(hm.integrated.dfree, regions = open_gamma_peaks)
  cgenes_acinar <- ClosestFeature(hm.integrated.dfree, regions = open_acinar_peaks)
  cgenes_ductal <- ClosestFeature(hm.integrated.dfree, regions = open_ductal_peaks)
  cgenes_qstel <- ClosestFeature(hm.integrated.dfree, regions = open_quiescentstellate_peaks)
  cgenes_astel <- ClosestFeature(hm.integrated.dfree, regions = open_activatedstellate_peaks)
  cgenes_macro <- ClosestFeature(hm.integrated.dfree, regions = open_macrophage_peaks)
  cgenes_lympho <- ClosestFeature(hm.integrated.dfree, regions = open_lymphocyte_peaks)
  cgenes_endo <- ClosestFeature(hm.integrated.dfree, regions = open_endothelial_peaks)
  
  cgenes_beta <- dplyr::filter(cgenes_beta, distance < 100000) 
  cgenes_alpha <- dplyr::filter(cgenes_alpha, distance < 100000) 
  cgenes_delta <- dplyr::filter(cgenes_delta, distance < 100000) 
  cgenes_gamma <- dplyr::filter(cgenes_gamma, distance < 100000) 
  cgenes_acinar <- dplyr::filter(cgenes_acinar, distance < 100000) 
  cgenes_ductal <- dplyr::filter(cgenes_ductal, distance < 100000) 
  cgenes_qstel <- dplyr::filter(cgenes_qstel, distance < 100000) 
  cgenes_astel <- dplyr::filter(cgenes_astel, distance < 100000) 
  cgenes_macro <- dplyr::filter(cgenes_macro, distance < 100000) 
  cgenes_lympho <- dplyr::filter(cgenes_lympho, distance < 100000) 
  cgenes_endo <- dplyr::filter(cgenes_endo, distance < 100000)
}

allregions <- c(
  as.character(cgenes_beta$query_region),
  as.character(cgenes_delta$query_region),
  as.character(cgenes_alpha$query_region),
  as.character(cgenes_gamma$query_region),
  as.character(cgenes_ductal$query_region),
  as.character(cgenes_acinar$query_region),
  as.character(cgenes_astel$query_region),
  as.character(cgenes_qstel$query_region),
  as.character(cgenes_endo$query_region),
  as.character(cgenes_lympho$query_region),
  as.character(cgenes_macro$query_region)
)

cgenes_beta <- distinct(cgenes_beta, gene_name, .keep_all = TRUE)
cgenes_alpha <- distinct(cgenes_alpha, gene_name, .keep_all = TRUE)
cgenes_delta <- distinct(cgenes_delta, gene_name, .keep_all = TRUE)
cgenes_gamma <- distinct(cgenes_gamma, gene_name, .keep_all = TRUE)
cgenes_acinar <- distinct(cgenes_acinar, gene_name, .keep_all = TRUE)
cgenes_ductal <- distinct(cgenes_ductal, gene_name, .keep_all = TRUE)
cgenes_qstel <- distinct(cgenes_qstel, gene_name, .keep_all = TRUE)
cgenes_astel <- distinct(cgenes_astel, gene_name, .keep_all = TRUE)
cgenes_macro <- distinct(cgenes_macro, gene_name, .keep_all = TRUE)
cgenes_lympho <- distinct(cgenes_lympho, gene_name, .keep_all = TRUE)
cgenes_endo <- distinct(cgenes_endo, gene_name, .keep_all = TRUE)


allgenes <- c(
  as.character(cgenes_beta$gene_name),
  as.character(cgenes_delta$gene_name),
  as.character(cgenes_alpha$gene_name),
  as.character(cgenes_gamma$gene_name),
  as.character(cgenes_ductal$gene_name),
  as.character(cgenes_acinar$gene_name),
  as.character(cgenes_astel$gene_name),
  as.character(cgenes_qstel$gene_name),
  as.character(cgenes_endo$gene_name),
  as.character(cgenes_lympho$gene_name),
  as.character(cgenes_macro$gene_name)
)

# Sex differences
beta_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\beta.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
alpha_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\alpha.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
delta_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\delta.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
gamma_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\gamma.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
ductal_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\ductal.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
acinar_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\acinar.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
qstel_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\quiescent_stellate.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
astell_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\activated_stellate.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
macro_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\macrophage.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
lympho_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\lymphocyte.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)
endo_peaks <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\sexancestry\endothelial.deseq.WaldTest.male.vs.female.tsv)", row.names = 1)

{
beta_peaks$padj[beta_peaks$padj == 0] <- 2e-302
beta_peaks <- dplyr::filter(beta_peaks, padj < 1e-5) 
open_beta_peaks_male <- rownames(beta_peaks[beta_peaks$log2FoldChange > 1 & beta_peaks$padj < 0.1, ])
open_beta_peaks_female <- rownames(beta_peaks[beta_peaks$log2FoldChange < -1 & beta_peaks$padj < 0.1, ])

alpha_peaks$padj[alpha_peaks$padj == 0] <- 2e-302
alpha_peaks <- dplyr::filter(alpha_peaks, padj < 1e-5) 
open_alpha_peaks_male <- rownames(alpha_peaks[alpha_peaks$log2FoldChange > 1 & alpha_peaks$padj < 0.1, ])
open_alpha_peaks_female <- rownames(alpha_peaks[alpha_peaks$log2FoldChange < -1 & alpha_peaks$padj < 0.1, ])

delta_peaks$padj[delta_peaks$padj == 0] <- 2e-302
delta_peaks <- dplyr::filter(delta_peaks, padj < 1e-5) 
open_delta_peaks_male <- rownames(delta_peaks[delta_peaks$log2FoldChange > 1 & delta_peaks$padj < 0.1, ])
open_delta_peaks_female <- rownames(delta_peaks[delta_peaks$log2FoldChange < -1 & delta_peaks$padj < 0.1, ])

gamma_peaks$padj[gamma_peaks$padj == 0] <- 2e-302
gamma_peaks <- dplyr::filter(gamma_peaks, padj < 1e-5) 
open_gamma_peaks_male <- rownames(gamma_peaks[gamma_peaks$log2FoldChange > 1 & gamma_peaks$padj < 0.1, ])
open_gamma_peaks_female <- rownames(gamma_peaks[gamma_peaks$log2FoldChange < -1 & gamma_peaks$padj < 0.1, ])

ductal_peaks$padj[ductal_peaks$padj == 0] <- 2e-302
ductal_peaks <- dplyr::filter(ductal_peaks, padj < 1e-5) 
open_ductal_peaks_male <- rownames(ductal_peaks[ductal_peaks$log2FoldChange > 1 & ductal_peaks$padj < 0.1, ])
open_ductal_peaks_female <- rownames(ductal_peaks[ductal_peaks$log2FoldChange < -1 & ductal_peaks$padj < 0.1, ])

acinar_peaks$padj[acinar_peaks$padj == 0] <- 2e-302
acinar_peaks <- dplyr::filter(acinar_peaks, padj < 1e-5) 
open_acinar_peaks_male <- rownames(acinar_peaks[acinar_peaks$log2FoldChange > 1 & acinar_peaks$padj < 0.1, ])
open_acinar_peaks_female <- rownames(acinar_peaks[acinar_peaks$log2FoldChange < -1 & acinar_peaks$padj < 0.1, ])

qstel_peaks$padj[qstel_peaks$padj == 0] <- 2e-302
qstel_peaks <- dplyr::filter(qstel_peaks, padj < 1e-5) 
open_qstel_peaks_male <- rownames(qstel_peaks[qstel_peaks$log2FoldChange > 1 & qstel_peaks$padj < 0.1, ])
open_qstel_peaks_female <- rownames(qstel_peaks[qstel_peaks$log2FoldChange < -1 & qstel_peaks$padj < 0.1, ])

astell_peaks$padj[astell_peaks$padj == 0] <- 2e-302
astell_peaks <- dplyr::filter(astell_peaks, padj < 1e-5) 
open_astell_peaks_male <- rownames(astell_peaks[astell_peaks$log2FoldChange > 1 & astell_peaks$padj < 0.1, ])
open_astell_peaks_female <- rownames(astell_peaks[astell_peaks$log2FoldChange < -1 & astell_peaks$padj < 0.1, ])

macro_peaks$padj[macro_peaks$padj == 0] <- 2e-302
macro_peaks <- dplyr::filter(macro_peaks, padj < 1e-5) 
open_macro_peaks_male <- rownames(macro_peaks[macro_peaks$log2FoldChange > 1 & macro_peaks$padj < 0.1, ])
open_macro_peaks_female <- rownames(macro_peaks[macro_peaks$log2FoldChange < -1 & macro_peaks$padj < 0.1, ])

lympho_peaks$padj[lympho_peaks$padj == 0] <- 2e-302
lympho_peaks <- dplyr::filter(lympho_peaks, padj < 1e-5) 
open_lympho_peaks_male <- rownames(lympho_peaks[lympho_peaks$log2FoldChange > 1 & lympho_peaks$padj < 0.1, ])
open_lympho_peaks_female <- rownames(lympho_peaks[lympho_peaks$log2FoldChange < -1 & lympho_peaks$padj < 0.1, ])

endo_peaks$padj[endo_peaks$padj == 0] <- 2e-302
endo_peaks <- dplyr::filter(endo_peaks, padj < 1e-5) 
open_endo_peaks_male <- rownames(endo_peaks[endo_peaks$log2FoldChange > 1 & endo_peaks$padj < 0.1, ])
open_endo_peaks_female <- rownames(endo_peaks[endo_peaks$log2FoldChange < -1 & endo_peaks$padj < 0.1, ])

cgenes_beta_male <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks_male)
cgenes_beta_female <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks_female)
cgenes_alpha_male <- ClosestFeature(hm.integrated.dfree, regions = open_alpha_peaks_male)
cgenes_alpha_female <- ClosestFeature(hm.integrated.dfree, regions = open_alpha_peaks_female)
cgenes_delta_male <- ClosestFeature(hm.integrated.dfree, regions = open_delta_peaks_male)
cgenes_delta_female <- ClosestFeature(hm.integrated.dfree, regions = open_delta_peaks_female)
cgenes_gamma_male <- ClosestFeature(hm.integrated.dfree, regions = open_gamma_peaks_male)
cgenes_gamma_female <- ClosestFeature(hm.integrated.dfree, regions = open_gamma_peaks_female)
cgenes_ductal_male <- ClosestFeature(hm.integrated.dfree, regions = open_ductal_peaks_male)
cgenes_ductal_female <- ClosestFeature(hm.integrated.dfree, regions = open_ductal_peaks_female)
cgenes_acinar_male <- ClosestFeature(hm.integrated.dfree, regions = open_acinar_peaks_male)
cgenes_acinar_female <- ClosestFeature(hm.integrated.dfree, regions = open_acinar_peaks_female)
cgenes_qstell_male <- ClosestFeature(hm.integrated.dfree, regions = open_qstell_peaks_male)
cgenes_qstell_female <- ClosestFeature(hm.integrated.dfree, regions = open_qstell_peaks_female)
cgenes_astell_male <- ClosestFeature(hm.integrated.dfree, regions = open_astell_peaks_male)
cgenes_astell_female <- ClosestFeature(hm.integrated.dfree, regions = open_astell_peaks_female)
cgenes_macro_male <- ClosestFeature(hm.integrated.dfree, regions = open_macro_peaks_male)
cgenes_macro_female <- ClosestFeature(hm.integrated.dfree, regions = open_macro_peaks_female)
cgenes_lympho_male <- ClosestFeature(hm.integrated.dfree, regions = open_lympho_peaks_male)
cgenes_lympho_female <- ClosestFeature(hm.integrated.dfree, regions = open_lympho_peaks_female)
cgenes_endo_male <- ClosestFeature(hm.integrated.dfree, regions = open_endo_peaks_male)
cgenes_endo_female <- ClosestFeature(hm.integrated.dfree, regions = open_endo_peaks_female)


cgenes_beta_male <- dplyr::filter(cgenes_beta_male, distance < 100000) 
cgenes_beta_female <- dplyr::filter(cgenes_beta_female, distance < 100000) 
cgenes_alpha_male <- dplyr::filter(cgenes_alpha_male, distance < 100000) 
cgenes_alpha_female <- dplyr::filter(cgenes_alpha_female, distance < 100000) 
cgenes_delta_male <- dplyr::filter(cgenes_delta_male, distance < 100000) 
cgenes_delta_female <- dplyr::filter(cgenes_delta_female, distance < 100000)
cgenes_gamma_male <- dplyr::filter(cgenes_gamma_male, distance < 100000) 
cgenes_gamma_female <- dplyr::filter(cgenes_gamma_female, distance < 100000)
cgenes_ductal_male <- dplyr::filter(cgenes_ductal_male, distance < 100000) 
cgenes_ductal_female <- dplyr::filter(cgenes_ductal_female, distance < 100000)
cgenes_acinar_male <- dplyr::filter(cgenes_acinar_male, distance < 100000) 
cgenes_acinar_female <- dplyr::filter(cgenes_acinar_female, distance < 100000)
cgenes_qstell_male <- dplyr::filter(cgenes_qstell_male, distance < 100000) 
cgenes_qstell_female <- dplyr::filter(cgenes_qstell_female, distance < 100000)
cgenes_astell_male <- dplyr::filter(cgenes_astell_male, distance < 100000) 
cgenes_astell_female <- dplyr::filter(cgenes_astell_female, distance < 100000)
cgenes_macro_male <- dplyr::filter(cgenes_macro_male, distance < 100000) 
cgenes_macro_female <- dplyr::filter(cgenes_macro_female, distance < 100000)
cgenes_lympho_male <- dplyr::filter(cgenes_lympho_male, distance < 100000) 
cgenes_lympho_female <- dplyr::filter(cgenes_lympho_female, distance < 100000)
cgenes_endo_male <- dplyr::filter(cgenes_endo_male, distance < 100000) 
cgenes_endo_female <- dplyr::filter(cgenes_endo_female, distance < 100000)
}

allregions <- c(
  as.character(cgenes_beta_male$query_region),
  as.character(cgenes_beta_female$query_region),
  as.character(cgenes_alpha_male$query_region),
  as.character(cgenes_alpha_female$query_region),
  as.character(cgenes_delta_male$query_region),
  as.character(cgenes_delta_female$query_region),
  as.character(cgenes_gamma_male$query_region),
  as.character(cgenes_gamma_female$query_region),
  as.character(cgenes_ductal_male$query_region),
  as.character(cgenes_ductal_female$query_region),
  as.character(cgenes_acinar_male$query_region),
  as.character(cgenes_acinar_female$query_region),
  as.character(cgenes_endo_male$query_region),
  as.character(cgenes_endo_female$query_region)
)

cgenes_beta_male <- distinct(cgenes_beta_male, gene_name, .keep_all = TRUE)
cgenes_beta_female <- distinct(cgenes_beta_female, gene_name, .keep_all = TRUE)
cgenes_alpha_male <- distinct(cgenes_alpha_male, gene_name, .keep_all = TRUE)
cgenes_alpha_female <- distinct(cgenes_alpha_female, gene_name, .keep_all = TRUE)
cgenes_delta_male <- distinct(cgenes_delta_male, gene_name, .keep_all = TRUE)
cgenes_delta_female <- distinct(cgenes_delta_female, gene_name, .keep_all = TRUE)
cgenes_gamma_male <- distinct(cgenes_gamma_male, gene_name, .keep_all = TRUE)
cgenes_gamma_female <- distinct(cgenes_gamma_female, gene_name, .keep_all = TRUE)
cgenes_ductal_male <- distinct(cgenes_ductal_male, gene_name, .keep_all = TRUE)
cgenes_ductal_female <- distinct(cgenes_ductal_female, gene_name, .keep_all = TRUE)
cgenes_acinar_male <- distinct(cgenes_acinar_male, gene_name, .keep_all = TRUE)
cgenes_acinar_female <- distinct(cgenes_acinar_female, gene_name, .keep_all = TRUE)
cgenes_endo_male <- distinct(cgenes_endo_male, gene_name, .keep_all = TRUE)
cgenes_endo_female <- distinct(cgenes_endo_female, gene_name, .keep_all = TRUE)



allgenes <- c(
  as.character(cgenes_beta_male$gene_name),
  as.character(cgenes_beta_female$gene_name),
  as.character(cgenes_alpha_male$gene_name),
  as.character(cgenes_alpha_female$gene_name),
  as.character(cgenes_delta_male$gene_name),
  as.character(cgenes_delta_female$gene_name),
  as.character(cgenes_gamma_male$gene_name),
  as.character(cgenes_gamma_female$gene_name),
  as.character(cgenes_ductal_male$gene_name),
  as.character(cgenes_ductal_female$gene_name),
  as.character(cgenes_acinar_male$gene_name),
  as.character(cgenes_acinar_female$gene_name),
  as.character(cgenes_endo_male$gene_name),
  as.character(cgenes_endo_female$gene_name)
)

allregions
allgenes

# Concatenate and remove dupliates
uniquegenes <- unique(allgenes)
uniqueregions <- unique(allregions)

# Heatmap
# Pseudobulk
DefaultAssay(hm.integrated.dfree) <- "ATAC"
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry_lib"
combined_processed_atac <- AverageExpression(hm.integrated.dfree, 
                                               assays = c("ATAC", "RNA"), 
                                               features = NULL, return.seurat = TRUE,  
                                               group.by = "celltype_sex_ancestry_lib",
                                               slot = "counts", verbose = FALSE)

combined_processed_atac$celltype_sex_ancestry_lib <- Cells(combined_processed_atac) #6892 Seurat

{
  Idents(combined_processed_atac) <- 'celltype_sex_ancestry_lib'
  combined_processed_atac$celltype <- combined_processed_atac@meta.data[["orig.ident"]]
  metadat <- combined_processed_atac@meta.data
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "quiescent_stellate", "quiescent-stellate"))
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex_ancestry_lib, "_", -3)
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -2)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -1)
  combined_processed_atac@meta.data = metadat
}


table(combined_processed_atac@meta.data[["celltype"]])
table(combined_processed_atac@meta.data[["sex"]])

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma",
               "ductal", "acinar",
               "activated", "quiescent", "endothelial",
               "lymphocyte", "macrophage")

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
table(combined_processed_atac$celltype)
Idents(combined_processed_atac) <- "celltype"


# Plot heatmap
genes.to.plot <- uniquegenes
label_genes <- c("INS", "GCG", "SST") #uniquegenes
all_genes_inobj <- rownames(combined_processed_rna@assays[["RNA"]])
genes.to.plot.intr <- intersect(all_genes_inobj, genes.to.plot) # to get genes make a seurat object by combining data over counts

regions.to.plot <- uniqueregions
label_genes <- c("INS", "GCG", "SST") #uniquegenes


write.csv(regions.to.plot, R"(C:\Users\mqadir\Desktop\regions.csv)")
DefaultAssay(combined_processed_atac) <- "predicted"
pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 8,
    height = 6)
dittoHeatmap(
  object = combined_processed_atac,#(subset(combined_processed_rna, idents = c("alpha"))),
  genes = regions.to.plot, #this is a compete gene set
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("celltype", "sex", "ancestry"),
  order.by = c("sex"),
  # main = NA,
  # cell.names.meta = NULL,
  # assay = .default_assay(object),
  # slot = .default_slot(object),
  # swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("dodgerblue", "white", "red3"))(50),
  # scaled.to.max = FALSE,
  # heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  # annot.colors = c(dittoColors(), dittoColors(1)[seq_len(7)]),
  # annotation_col = NULL,
  annotation_colors = list(celltype = c("acinar" = "salmon3",
                                        "activated" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophage" = "magenta3",
                                        "quiescent" = "salmon"),
                           sex = c("female" = "red4",
                                   "male" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black")),
  # # data.out = FALSE,
  # highlight.features = c("INS", "MAFA", "IAPP", "ENTPD3", "NKX6-1", "PDX1", #beta
  #                        "GCG", "TTR", "IRX2", "ARX", "TM4SF4", "PCSK1N", #alpha
  #                        "SST", "RBP4", "HHEX", "LY6H", "F5", "BHLHE41", #delta
  #                        "PPY", #gamma
  #                        "GHRL", #epsilon
  #                        "TOP2A", "CCNB2", "HMGB2", "CDKN3", "MKI67", "CENPF", #cycendo
  #                        "SPP1", "TFPI2", "KRT19", "ONECUT1", "TM4SF1", #ductal
  #                        "CTRB1", "CTRB2", "PRSS2", "PRSS1", "PNLIP", "CELA2A", #acinar
  #                        "SFRP2", "VIM", "DCN", "COL1A1", "LUM", "PTGDS", #activated
  #                        "GADD45B", "HMGB1", "PDGFRB", "PRDX1", "PTMA", "RGS5", #quiescent
  #                        "PECAM1", "VWF", "SOX18", "FCN3", "CD59", "ESM1", #endo
  #                        "CCL5", "NKG7", "CD3E", "IL32", "TRAC", "HLA-B", #lymphocyte
  #                        "TPSAB1", "TPSB2", #mast
  #                        "CRYAB", "SOX10", "NGFR", "RUNX2", "BTC", "CDH19", #schwann
  #                        "SDS", "C1QB", "CD68", "APOE", "VMO1", "MS4A7"), #macrophages
  #right_annotations = rowAnnotation(foo = anno_mark(at = c(1), labels = c("HHEX"))),
  # show_colnames = isBulk(object),
  show_rownames = TRUE,
  # scale = "row",
  cluster_row = TRUE,
  cluster_cols = FALSE,
  # border_color = NA,
  # legend_breaks = NA,
  # drop_levels = FALSE,
  breaks=seq(-2, 2, length.out=50),
  complex = TRUE,
  #column_km = 1,
  use_raster = TRUE,
  raster_quality = 5,
  #column_split = combined_processed_rna$celltype,
  #border_color = "black",
  #gaps_col = c(15, 30, 45, 60, 75, 90, 105, 120, 135, 150),
  #gaps_row = c(1819, 3205, 5170, 5405, 10850, 14650, 21250, 23500, 25750, 26200)
  #gaps_row = c(2804, 4584, 7464, 7799, 17833, 23493, 36665, 39029, 41669, 42089)
) 
dev.off()
dev.off()

+ rowAnnotation(mark = anno_mark(at = match(label_genes, 
                                              rownames(combined_processed_atac[genes.to.plot,])), 
                                   labels = label_genes, 
                                   which = "row",
                                   labels_gp = list(cex=0.3),
                                   #link_width = unit(4, "mm"), link_height = unit(4, "mm"),
                                   padding = 0.1
)
)


# GO plot This uses allgenes that was derived from above
# Compare

# create gene list
# Sort out top 1000 accessible sites
cgenes_beta <- cgenes_beta %>% slice_max(cgenes_beta, n = 1000)
cgenes_delta <- cgenes_delta %>% slice_max(cgenes_delta, n = 1000)
cgenes_alpha <- cgenes_alpha %>% slice_max(cgenes_alpha, n = 1000)
cgenes_gamma <- cgenes_gamma %>% slice_max(cgenes_gamma, n = 1000)
cgenes_ductal <- cgenes_ductal %>% slice_max(cgenes_ductal, n = 1000)
cgenes_acinar <- cgenes_acinar %>% slice_max(cgenes_acinar, n = 1000)
cgenes_astel <- cgenes_astel %>% slice_max(cgenes_astel, n = 1000)
cgenes_qstel <- cgenes_qstel %>% slice_max(cgenes_qstel, n = 1000)
cgenes_endo <- cgenes_endo %>% slice_max(cgenes_endo, n = 1000)
cgenes_lympho <- cgenes_lympho %>% slice_max(cgenes_lympho, n = 1000)
cgenes_macro <- cgenes_macro %>% slice_max(cgenes_macro, n = 1000)

# For accessibility remmeber to run the code again generating the _bete files and subset 100kb window peaks
cgenes_beta <- as.character(cgenes_beta$gene_name)
cgenes_delta <- as.character(cgenes_delta$gene_name)
cgenes_alpha <- as.character(cgenes_alpha$gene_name)
cgenes_gamma <- as.character(cgenes_gamma$gene_name)
cgenes_ductal <- as.character(cgenes_ductal$gene_name)
cgenes_acinar <- as.character(cgenes_acinar$gene_name)
cgenes_astel <- as.character(cgenes_astel$gene_name)
cgenes_qstel <- as.character(cgenes_qstel$gene_name)
cgenes_endo <- as.character(cgenes_endo$gene_name)
cgenes_lympho <- as.character(cgenes_lympho$gene_name)
cgenes_macro <- as.character(cgenes_macro$gene_name)

#across sex
beta_male_ac <- as.character(cgenes_beta_male$gene_name)
beta_female_ac <- as.character(cgenes_beta_female$gene_name)
alpha_male_ac <- as.character(cgenes_alpha_male$gene_name)
alpha_female_ac <- as.character(cgenes_alpha_female$gene_name)
delta_male_ac <- as.character(cgenes_delta_male$gene_name)
delta_female_ac <- as.character(cgenes_delta_female$gene_name)
gamma_male_ac <- as.character(cgenes_gamma_male$gene_name)
gamma_female_ac <-as.character(cgenes_gamma_female$gene_name)
ductal_male_ac <- as.character(cgenes_ductal_male$gene_name)
ductal_female_ac <- as.character(cgenes_ductal_female$gene_name)
acinar_male_ac <- as.character(cgenes_acinar_male$gene_name)
acinar_female_ac <- as.character(cgenes_acinar_female$gene_name)
endo_male_ac <- as.character(cgenes_endo_male$gene_name)
endo_female_ac <- as.character(cgenes_endo_female$gene_name)


gene.list <- list("beta" = cgenes_beta, "alpha" = cgenes_alpha, "delta" = cgenes_delta, "gamma" = cgenes_gamma, 
                  "acinar" = cgenes_acinar, "ductal" = cgenes_ductal, "qstel" = cgenes_qstel, "astel" = cgenes_astel,
                  "macro" = cgenes_macro, "lympho" = cgenes_lympho, "endo" = cgenes_endo)

gene.list <- list("beta_male_ac" = beta_male_ac, "alpha_male_ac" = alpha_male_ac, "delta_male_ac" = delta_male_ac, "gamma_male_ac" = gamma_male_ac, "ductal_male_ac" = ductal_male_ac, "acinar_male_ac" = acinar_male_ac, "endo_male_ac" = endo_male_ac,
                  "beta_female_ac" = beta_female_ac,  "alpha_female_ac" = alpha_female_ac,  "delta_female_ac" = delta_female_ac,  "gamma_female_ac" = gamma_female_ac, "ductal_female_ac" = ductal_female_ac, "acinar_female_ac" = acinar_female_ac, "endo_female_ac" = endo_female_ac)

ck <- compareCluster(geneCluster = gene.list, # list of genes
                     fun = enrichGO, 
                     #universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                     keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                     OrgDb = org.Hs.eg.db, 
                     ont = c("ALL"), 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1, #if not set default is at 0.05
                     readable = TRUE) 
head(ck) 
cluster_summary <- data.frame(ck)
ck <- ck[ck@compareClusterResult[["qvalue"]] < 0.1, asis=T]
dotplot(ck, showCategory = 8) + #coord_flip() + 
  scale_x_discrete(limits=rev) + scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
dotplot(ck, showCategory = 1)
ck.save <- ck@compareClusterResult
write.csv(ck.save, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\ORA\sex_accessibility.save.csv)")

dotplot(ck, 
        # showCategory = c("insulin secretion", "calcium ion homeostasis", "protein localization to plasma membrane", "glucose homeostasis",
        #                      "regulation of G protein-coupled receptor signaling pathway", "pancreas development", "hormone metabolic process",
        #                      "glutamate receptor signaling pathway", "synaptic transmission, GABAergic", "peptide secretion",
        #                      "hormone transport",
        #                      "developmental growth involved in morphogenesis", "sodium ion transport", "regulation of actin filament-based process",
        #                      "regulation of Wnt signaling pathway", "integrin-mediated signaling pathway",
        #                      "muscle contraction", "establishment of endothelial barrier", "cellular response to fibroblast growth factor stimulus", "extracellular matrix organization",
        #                      "activation of immune response", "response to interferon-gamma", "cytokine-mediated signaling pathway", "leukocyte chemotaxis", "leukocyte degranulation", "interleukin-1 beta production", "phagocytosis",
        #                      "alpha-beta T cell activation",
        #                      "endothelial cell development", "cellular response to angiotensin"), 
        showCategory = c("histone lysine demethylation", "protein dealkylation", "histone modification", "sequestering of actin monomers", "gonad development", "sex determination", 
                         "androgen receptor signaling pathway", 
                         "dosage compensation", "regulation of gene expression, epigenetic", "chromatin remodeling"),
        font.size=14) + #coord_flip() + 
  scale_x_discrete(limits=rev) + scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

cnetplot(ck)

beta.alpha.delta <- list(
  delta_genes=as.character(delta_genes),
  beta_genes=as.character(beta_genes),
  alpha_genes=as.character(alpha_genes)
)


# Compare
ck.bad <- compareCluster(geneCluster = beta.alpha.delta, 
                         fun = enrichGO, 
                         universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                         keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                         OrgDb = org.Hs.eg.db, 
                         ont = c("ALL"), 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 0.1, #if not set default is at 0.05
                         readable = TRUE)
ck.bad <- setReadable(ck.bad, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
cnetplot(ck.bad,
         showCategory = c("gamma-aminobutyric acid signaling pathway", "hormone secretion",
                          "peptide transport", "peptide hormone secretion", "calcium-ion regulated exocytosis",
                          "neurotransmitter secretion", "Golgi to endosome transport", "potassium channel complex"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 2,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", "dodgerblue3", "lightseagreen"))

options(ggrepel.max.overlaps = Inf)
cnetplot(ck,
         showCategory = c("synapse organization", "gamma-aminobutyric acid signaling pathway",
                          "insulin secretion", "cilium assembly", "peptide hormone secretion", 
                          "cellular response to glucose starvation", "neurotransmitter secretion", "amide transport",
                          "neuropeptide signaling pathway", "protein secretion", "glucagon secretion",
                          "glucocorticoid secretion", "growth hormone secretion", "positive regulation of feeding behavior",
                          "nuclear division", "mitotic cell cycle phase transition", "organelle fission",
                          "epithelial cell proliferation", "digestive tract development", "water homeostasis", "organic anion transport", "SMAD protein signal transduction",
                          "digestion", "morphogenesis of a branching structure", "primary alcohol metabolic process",
                          "extracellular matrix organization", "collagen fibril organization",
                          "muscle contraction", "muscle cell differentiation", "regulation of systemic arterial blood pressure by hormone",
                          "regulation of angiogenesis", "blood vessel endothelial cell migration",
                          "T cell activation", "lymphocyte mediated immunity", "T cell selection",
                          "myeloid leukocyte activation", "antigen processing and presentation", "cell chemotaxis",
                          "immune response-regulating cell surface receptor signaling pathway", "mast cell activation", "activation of immune response",
                          "central nervous system myelination", "ensheathment of neurons", "axon development"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 10,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", #"delta" = 
                                                            "dodgerblue3", #"beta" = ,
                                                            "turquoise2", #"beta+alpha" =
                                                            "lightseagreen", #"alpha"= 
                                                            "springgreen4", #"gamma" =
                                                            "khaki2", #"epsilon" = 
                                                            "darkseagreen2", #"cycling-endo" = 
                                                            "darkorange2", #"ductal" =
                                                            "salmon3", #"acinar" = 
                                                            "orange", #"activated-stellate" = 
                                                            "salmon", #"quiescent-stellate" = 
                                                            "red", #"endothelial" = 
                                                            "orchid1", #"lymphocyte" = 
                                                            "magenta3", #"macrophages" = 
                                                            "red4", #"mast" = 
                                                            "grey30"#"schwann" = 
         ))

eg <- bitr(as.character(alpha_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
edox <- enrichDGN(as.character(eg$ENTREZID), readable = TRUE)
edox <- setReadable(edox, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
edox <- pairwise_termsim(edox)
emapplot(ck)
treeplot(edox)
mutate(edox, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")



# Motif analysis an plotting
# For accessibility remmeber to run the code again generating the _bete files and subset 100kb window peaks
allregions <- c(
  as.character(cgenes_beta$query_region),
  as.character(cgenes_delta$query_region),
  as.character(cgenes_alpha$query_region),
  as.character(cgenes_gamma$query_region),
  as.character(cgenes_ductal$query_region),
  as.character(cgenes_acinar$query_region),
  as.character(cgenes_astel$query_region),
  as.character(cgenes_qstel$query_region),
  as.character(cgenes_endo$query_region),
  as.character(cgenes_lympho$query_region),
  as.character(cgenes_macro$query_region)
)

uniqueregions <- unique(allregions)

peaks_beta <- as.character(cgenes_beta$query_region)
peaks_delta <- as.character(cgenes_delta$query_region)
peaks_alpha <- as.character(cgenes_alpha$query_region)
peaks_gamma <- as.character(cgenes_gamma$query_region)
peaks_ductal <- as.character(cgenes_ductal$query_region)
peaks_acinar <- as.character(cgenes_acinar$query_region)
peaks_astel <- as.character(cgenes_astel$query_region)
peaks_qstel <- as.character(cgenes_qstel$query_region)
peaks_endo <- as.character(cgenes_endo$query_region)
peaks_lympho <- as.character(cgenes_lympho$query_region)
peaks_macro <- as.character(cgenes_macro$query_region)

# Motif analysis
Idents(hm.integrated.dfree) <- "celltype"
open.peaks <- AccessiblePeaks(hm.integrated.dfree, idents = unique(as.character(hm.integrated.dfree$celltype)))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(hm.integrated.dfree, assay = "ATAC", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[uniqueregions, ],
  n = 50000
)

# Motif tests
# This is not super accurate its better to run motif analysis using chromvar
{
enriched.motifs.beta <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_beta)
enriched.motifs.delta <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_delta)
enriched.motifs.alpha <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_alpha)
enriched.motifs.gamma <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_gamma)
enriched.motifs.ductal <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_ductal)
enriched.motifs.acinar <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_acinar)
enriched.motifs.astel <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_astel)
enriched.motifs.qstel <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_qstel)
enriched.motifs.endo <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_endo)
enriched.motifs.lympho <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_lympho)
enriched.motifs.macro <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_macro)

write.csv(enriched.motifs.beta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.beta.csv)")
write.csv(enriched.motifs.delta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.delta.csv)")
write.csv(enriched.motifs.alpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.alpha.csv)")
write.csv(enriched.motifs.gamma, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.gamma.csv)")
write.csv(enriched.motifs.ductal, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.ductal.csv)")
write.csv(enriched.motifs.acinar, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.acinar.csv)")

write.csv(enriched.motifs.astel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.astel.csv)")
write.csv(enriched.motifs.qstel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.qstel.csv)")
write.csv(enriched.motifs.endo, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.endo.csv)")
write.csv(enriched.motifs.lympho, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.lympho.csv)")
write.csv(enriched.motifs.macro, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.macro.csv)")
}


# Add Motif scores to pseduo bulk and make a heat map
# Heatmap
# Pseudobulk
DefaultAssay(hm.integrated.dfree) <- "chromvar"
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry_lib"
combined_processed_atac <- AverageExpression(hm.integrated.dfree, 
                                             assays = c("chromvar", "ATAC", "RNA"), 
                                             features = NULL, return.seurat = TRUE,  
                                             group.by = "celltype_sex_ancestry_lib",
                                             slot = "data", verbose = FALSE)

combined_processed_atac$celltype_sex_ancestry_lib <- Cells(combined_processed_atac) #6892 Seurat

{
  Idents(combined_processed_atac) <- 'celltype_sex_ancestry_lib'
  combined_processed_atac$celltype <- combined_processed_atac@meta.data[["orig.ident"]]
  metadat <- combined_processed_atac@meta.data
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "quiescent_stellate", "quiescent-stellate"))
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex_ancestry_lib, "_", -3)
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -2)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -1)
  combined_processed_atac@meta.data = metadat
}


table(combined_processed_atac@meta.data[["celltype"]])
table(combined_processed_atac@meta.data[["sex"]])

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma",
               "ductal", "acinar",
               "activated", "quiescent", "endothelial",
               "lymphocyte", "macrophage")

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
table(combined_processed_atac$celltype)
Idents(combined_processed_atac) <- "celltype"


# # DE testing Motifs
# DefaultAssay(hm.integrated.dfree) <- "chromvar"
# Idents(hm.integrated.dfree) <- "celltype"
# 
# enriched.cvar.motifs.beta <- FindMarkers(
#   object = hm.integrated.dfree,
#   ident.1 = 'beta',
#   ident.2 = c('alpha', "delta", "gamma", "activated_stellate", "quiescent_stellate", "ductal", "acinar", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.2,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.alpha <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'alpha',
#   ident.2 = c('beta', "delta", "gamma", "activated", "quiescent", "ductal", "acinar", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.05,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.delta <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'delta',
#   ident.2 = c('alpha', "beta", "gamma", "activated", "quiescent", "ductal", "acinar", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.gamma <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'gamma',
#   ident.2 = c('alpha', "beta", "delta", "activated", "quiescent", "ductal", "acinar", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.astel <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'activated',
#   ident.2 = c('alpha', "beta", "delta", "gamma", "quiescent", "ductal", "acinar", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.qstel <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'quiescent',
#   ident.2 = c('alpha', "beta", "delta", "gamma", "activated", "ductal", "acinar", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.ductal <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'ductal',
#   ident.2 = c('alpha', "beta", "delta", "gamma", "activated", "quiescent", "acinar", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.acinar <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'acinar',
#   ident.2 = c('alpha', "beta", "delta", "gamma", "activated", "quiescent", "ductal", "macrophage", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.macro <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'macrophage',
#   ident.2 = c('alpha', "beta", "delta", "gamma", "activated", "quiescent", "ductal", "acinar", "lymphocyte", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.lympho <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'lymphocyte',
#   ident.2 = c('alpha', "beta", "delta", "gamma", "activated", "quiescent", "ductal", "acinar", "macrophage", "endothelial"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.endo <- FindMarkers(
#   object = combined_processed_atac,
#   ident.1 = 'endothelial',
#   ident.2 = c('alpha', "beta", "delta", "gamma", "activated", "quiescent", "ductal", "acinar", "macrophage", "lymphocyte"),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff",
#   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
#   min.pct = 0.1,
#   logfc.threshold = 0
# )
# 
# enriched.cvar.motifs.beta
# enriched.cvar.motifs.alpha
# enriched.cvar.motifs.delta
# enriched.cvar.motifs.gamma
# enriched.cvar.motifs.astel
# enriched.cvar.motifs.qstel
# enriched.cvar.motifs.ductal
# enriched.cvar.motifs.acinar
# enriched.cvar.motifs.macro
# enriched.cvar.motifs.lympho
# enriched.cvar.motifs.endo
# 
# # Translating Motifs
# motif_id <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\motif.ID.csv)", row.names = 1)
# enriched.motifs.beta <- merge(motif_id, enriched.cvar.motifs.beta, by = 'row.names')
# enriched.motifs.delta <- merge(motif_id, enriched.cvar.motifs.delta, by = 'row.names')
# enriched.motifs.alpha <- merge(motif_id, enriched.cvar.motifs.alpha, by = 'row.names')
# enriched.motifs.gamma <- merge(motif_id, enriched.cvar.motifs.gamma, by = 'row.names')
# enriched.motifs.astel <- merge(motif_id, enriched.cvar.motifs.astel, by = 'row.names')
# enriched.motifs.qstel <- merge(motif_id, enriched.cvar.motifs.qstel, by = 'row.names')
# enriched.motifs.ductal <- merge(motif_id, enriched.cvar.motifs.ductal, by = 'row.names')
# enriched.motifs.acinar <- merge(motif_id, enriched.cvar.motifs.acinar, by = 'row.names')
# enriched.motifs.macro <- merge(motif_id, enriched.cvar.motifs.macro, by = 'row.names')
# enriched.motifs.lympho <- merge(motif_id, enriched.cvar.motifs.lympho, by = 'row.names')
# enriched.motifs.endo <- merge(motif_id, enriched.cvar.motifs.endo, by = 'row.names')
# 
# enriched.motifs.beta <- enriched.motifs.beta %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.delta <- enriched.motifs.delta %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.alpha <- enriched.motifs.alpha %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.gamma <- enriched.motifs.gamma %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.astel <- enriched.motifs.astel %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.qstel <- enriched.motifs.qstel %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.ductal <- enriched.motifs.ductal %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.acinar <- enriched.motifs.acinar %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.macro <- enriched.motifs.macro %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.lympho <- enriched.motifs.lympho %>% remove_rownames %>% column_to_rownames(var="Row.names")
# enriched.motifs.endo <- enriched.motifs.endo %>% remove_rownames %>% column_to_rownames(var="Row.names")
# 
# write.csv(enriched.motifs.beta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.beta.csv)")
# write.csv(enriched.motifs.alpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.alpha.csv)")
# write.csv(enriched.motifs.delta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.delta.csv)")
# write.csv(enriched.motifs.gamma, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.gamma.csv)")
# write.csv(enriched.motifs.astel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.activated.stellate.csv)")
# write.csv(enriched.motifs.qstel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.quiescent.stellate.csv)")
# write.csv(enriched.motifs.ductal, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.ductal.csv)")
# write.csv(enriched.motifs.acinar, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.acinar.csv)")
# write.csv(enriched.motifs.macro, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.macro.csv)")
# write.csv(enriched.motifs.lympho, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.lympho.csv)")
# write.csv(enriched.motifs.endo, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.endo.csv)")

# Load in motifs From 
enriched.motifs.beta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.beta.csv)", row.names = 1)
enriched.motifs.alpha <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.alpha.csv)", row.names = 1)
enriched.motifs.delta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.delta.csv)", row.names = 1)
enriched.motifs.gamma <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.gamma.csv)", row.names = 1)
enriched.motifs.ductal <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.ductal.csv)", row.names = 1)
enriched.motifs.acinar <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.acinar.csv)", row.names = 1)
enriched.motifs.astel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.astel.csv)", row.names = 1)
enriched.motifs.qstel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.qstel.csv)", row.names = 1)
enriched.motifs.endo <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.endo.csv)", row.names = 1)
enriched.motifs.lympho <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.lympho.csv)", row.names = 1)
enriched.motifs.macro <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.macro.csv)", row.names = 1)

# Plotting heatmap from ->
# Adjust for incalculable pvals
enriched.motifs.beta$p.adjust[enriched.motifs.beta$p.adjust == 0] <- 2e-307
enriched.motifs.delta$p.adjust[enriched.motifs.delta$p.adjust == 0] <- 2e-307
enriched.motifs.alpha$p.adjust[enriched.motifs.alpha$p.adjust == 0] <- 2e-307
enriched.motifs.gamma$p.adjust[enriched.motifs.gamma$p.adjust == 0] <- 2e-307
enriched.motifs.astel$p.adjust[enriched.motifs.astel$p.adjust == 0] <- 2e-307
enriched.motifs.qstel$p.adjust[enriched.motifs.qstel$p.adjust == 0] <- 2e-307
enriched.motifs.ductal$p.adjust[enriched.motifs.ductal$p.adjust == 0] <- 2e-307
enriched.motifs.acinar$p.adjust[enriched.motifs.acinar$p.adjust == 0] <- 2e-307
enriched.motifs.macro$p.adjust[enriched.motifs.macro$p.adjust == 0] <- 2e-307
enriched.motifs.lympho$p.adjust[enriched.motifs.lympho$p.adjust == 0] <- 2e-307
enriched.motifs.endo$p.adjust[enriched.motifs.endo$p.adjust == 0] <- 2e-307

enriched.motifs.beta <- dplyr::filter(enriched.motifs.beta, p.adjust < 1e-10) 
enriched.motifs.delta <- dplyr::filter(enriched.motifs.delta, p.adjust < 1e-10) 
enriched.motifs.alpha <- dplyr::filter(enriched.motifs.alpha, p.adjust < 1e-10) 
enriched.motifs.gamma <- dplyr::filter(enriched.motifs.gamma, p.adjust < 1e-10) 
enriched.motifs.astel <- dplyr::filter(enriched.motifs.astel, p.adjust < 1e-10) 
enriched.motifs.qstel <- dplyr::filter(enriched.motifs.qstel, p.adjust < 1e-10) 
enriched.motifs.ductal <- dplyr::filter(enriched.motifs.ductal, p.adjust < 1e-10) 
enriched.motifs.acinar <- dplyr::filter(enriched.motifs.acinar, p.adjust < 1e-10) 
enriched.motifs.macro <- dplyr::filter(enriched.motifs.macro, p.adjust < 1e-10) 
enriched.motifs.lympho <- dplyr::filter(enriched.motifs.lympho, p.adjust < 1e-10)
enriched.motifs.endo <- dplyr::filter(enriched.motifs.endo, p.adjust < 1e-10)

# Order data #beta_motifs <- as.character(rownames(enriched.motifs.beta[enriched.motifs.beta$p_val_adj < 1e-50,]))
enriched.motifs.beta <- enriched.motifs.beta[order(enriched.motifs.beta$fold.enrichment),] 
enriched.motifs.delta <- enriched.motifs.delta[order(enriched.motifs.delta$fold.enrichment),] 
enriched.motifs.alpha <- enriched.motifs.alpha[order(enriched.motifs.alpha$fold.enrichment),] 
enriched.motifs.gamma <- enriched.motifs.gamma[order(enriched.motifs.gamma$fold.enrichment),]
enriched.motifs.astel <- enriched.motifs.astel[order(enriched.motifs.astel$fold.enrichment),] 
enriched.motifs.qstel <- enriched.motifs.qstel[order(enriched.motifs.qstel$fold.enrichment),]
enriched.motifs.ductal <- enriched.motifs.ductal[order(enriched.motifs.ductal$fold.enrichment),] 
enriched.motifs.acinar <- enriched.motifs.acinar[order(enriched.motifs.acinar$fold.enrichment),]
enriched.motifs.macro <- enriched.motifs.macro[order(enriched.motifs.macro$fold.enrichment),] 
enriched.motifs.lympho <- enriched.motifs.lympho[order(enriched.motifs.lympho$fold.enrichment),]
enriched.motifs.endo <- enriched.motifs.endo[order(enriched.motifs.endo$fold.enrichment),]

# -log10adjpval
enriched.motifs.beta$neglogpval <- -log10(enriched.motifs.beta$p.adjust)
enriched.motifs.delta$neglogpval <- -log10(enriched.motifs.delta$p.adjust)
enriched.motifs.alpha$neglogpval <- -log10(enriched.motifs.alpha$p.adjust)
enriched.motifs.gamma$neglogpval <- -log10(enriched.motifs.gamma$p.adjust)
enriched.motifs.astel$neglogpval <- -log10(enriched.motifs.astel$p.adjust)
enriched.motifs.qstel$neglogpval <- -log10(enriched.motifs.qstel$p.adjust)
enriched.motifs.ductal$neglogpval <- -log10(enriched.motifs.ductal$p.adjust)
enriched.motifs.acinar$neglogpval <- -log10(enriched.motifs.acinar$p.adjust)
enriched.motifs.macro$neglogpval <- -log10(enriched.motifs.macro$p.adjust)
enriched.motifs.lympho$neglogpval <- -log10(enriched.motifs.lympho$p.adjust)
enriched.motifs.endo$neglogpval <- -log10(enriched.motifs.endo$p.adjust)

# Plot
plot_celldata <- top_n(enriched.motifs.gamma, 75, fold.enrichment)
level_order <- plot_celldata$motif.name
mid <- mean(plot_celldata$fold.enrichment)
plot_celldata %>%
  arrange(fold.enrichment) %>%
  mutate(name=factor(motif.name, levels=motif.name)) %>%
  ggplot(aes(x = factor(motif.name, level = level_order),
             y = fold.enrichment, 
             color = neglogpval, 
             size = percent.observed)) + 
  geom_point(alpha=0.9) + 
  ylab("") + 
  xlab("") + 
  ggtitle("") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3, hjust=1, size =8, face = "plain", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =8, face = "plain", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "plain"),
        legend.title=element_text(size=8, face = "plain"), 
        legend.text=element_text(size=8, face = "plain")) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  coord_flip() +
  scale_size_continuous(
    breaks = c(0, 25, 50, 100),
    limits = c(0, 100)) +
  #scale_x_discrete(limits=rev) + 
  scale_color_gradient2(mid="white", high="springgreen4", space ="Lab" )

# TF Foot printing
# gather the footprinting information for sets of motifs
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
DefaultAssay(hm.integrated.dfree) <- "ATAC"
Idents(hm.integrated.dfree) <- "celltype"
hm.integrated.dfree <- Footprint(
  object = hm.integrated.dfree,
  motif.name = c("PDX1", "MAFA", "Arx", "NKX6-1", "Isl1", "GBX2", "GATA5", "GLIS3"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

PlotFootprint(hm.integrated.dfree, features = c("GATA5"))

# Corelation plots
# RNA
install.packages("ggcorrplot")
library(ggcorrplot)
seurobj <- processed_rna
cluster.id <- "beta"
gene.list <- c("INS", "GCG", "PDX1")
clusterCorPlot <- function(seurObj, cluster.id, gene.list, assay='RNA', slot='data') {
  gene.dat <- GetAssayData(
    subset(seurObj, idents = cluster.id), 
    assay=assay, 
    slot=slot)[gene.list,]
  ggcorrplot(as.data.frame(t(gene.dat)))
}

processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
Idents(processed_rna) <- "celltype_sex_ancestry_disease"
av.exp <- AverageExpression(processed_rna)$RNA
av.exp <- av.exp[rownames(av.exp) %in% processed_rna@assays[["RNA"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, unique(processed_rna@meta.data[["celltype_sex_ancestry_disease"]]))
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
    ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                    midpoint = 0.5, limit = c(0,1), space = "Lab", 
                                    name="Pearson\nCorrelation")
hm.integrated.dfree

# ATAC
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry"
av.exp <- AverageExpression(hm.integrated.dfree)$RNA
DefaultAssay(hm.integrated.dfree) <- "RNA"
hm.integrated.dfree <- FindVariableFeatures(hm.integrated.dfree, selection.method = "vst", nfeatures = 2000)
av.exp <- av.exp[rownames(av.exp) %in% hm.integrated.dfree@assays[["RNA"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
datanames_toplot <- unique(hm.integrated.dfree@meta.data[["celltype_sex_ancestry"]])
datanames_toplot <- as.character(datanames_toplot[!is.na(datanames_toplot)])
cor.df <- tidyr::gather(data = cor.exp, y, correlation, datanames_toplot)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
  ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                                                                                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                                                                                                         name="Pearson\nCorrelation")

av.exp <- AverageExpression(hm.integrated.dfree)$ATAC
DefaultAssay(hm.integrated.dfree) <- "ATAC"
av.exp <- av.exp[rownames(av.exp) %in% hm.integrated.dfree@assays[["ATAC"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
datanames_toplot <- unique(hm.integrated.dfree@meta.data[["celltype_sex_ancestry"]])
datanames_toplot <- as.character(datanames_toplot[!is.na(datanames_toplot)])
cor.df <- tidyr::gather(data = cor.exp, y, correlation, datanames_toplot)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
  ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                                                                                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                                                                                                         name="Pearson\nCorrelation")

av.exp <- AverageExpression(hm.integrated.dfree)$predicted
DefaultAssay(hm.integrated.dfree) <- "predicted"
hm.integrated.dfree <- FindVariableFeatures(hm.integrated.dfree, selection.method = "vst", nfeatures = 2000)
av.exp <- av.exp[rownames(av.exp) %in% hm.integrated.dfree@assays[["predicted"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
datanames_toplot <- unique(hm.integrated.dfree@meta.data[["celltype_sex_ancestry"]])
datanames_toplot <- as.character(datanames_toplot[!is.na(datanames_toplot)])
cor.df <- tidyr::gather(data = cor.exp, y, correlation, datanames_toplot)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
  ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                                                                                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                                                                                                         name="Pearson\nCorrelation")


# Load in motifs From 
enriched.motifs.beta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.beta.csv)")
enriched.motifs.delta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.alpha.csv)")
enriched.motifs.alpha <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.delta.csv)")
enriched.motifs.gamma <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.gamma.csv)")
enriched.motifs.ductal <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.ductal.csv)")
enriched.motifs.acinar <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.acinar.csv)")
enriched.motifs.astel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.astel.csv)")
enriched.motifs.qstel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.qstel.csv)")
enriched.motifs.endo <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.endo.csv)")
enriched.motifs.lympho <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.lympho.csv)")
enriched.motifs.macro <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\enriched.motifs.macro.csv)")

motifs.to.plot <- rownames(hm.integrated.dfree@assays[["chromvar"]])
combined_processed_atac <- FindVariableFeatures(combined_processed_atac, selection.method = "vst", nfeatures = 500)

DefaultAssay(combined_processed_atac) <- "ATAC"
pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 10,
    height = 9)
genes.to.plot <- combined_processed_atac@assays[["chromvar"]]@var.features #complete geneset
label_genes <- c("MA0132.2", #PDX1
                 "MA0148.4",  #FOXA1
                 "MA0874.1",  #Arx
                 "MA1608.1", #ISL1
                 "MA0674.1", #NKX6-1
                 "MA1645.1", #NKX2-2
                 "MA0117.2", #MAFB
                 "MA0077.1", #SOX9
                 "MA0068.2", #PAX4
                 "MA1618.1", #PTF1A
                 "MA0114.4", #HNF4A
                 "MA0046.2", #HNF1A
                 "MA0084.1", #SRY
                 "MA0007.3", #Ar
                 "MA0098.3", #ETS1
                 "MA1484.1" #ETS2
                 )
dittoHeatmap(
  object = combined_processed_atac,#(subset(combined_processed_rna, idents = c("alpha"))),
  genes = genes.to.plot, #this is a compete gene set
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("celltype", "sex", "ancestry"),
  order.by = c("celltype", "sex"),
  # main = NA,
  # cell.names.meta = NULL,
  # assay = .default_assay(object),
  # slot = .default_slot(object),
  # swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("dodgerblue", "white", "red3"))(50),
  # scaled.to.max = FALSE,
  # heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  # annot.colors = c(dittoColors(), dittoColors(1)[seq_len(7)]),
  # annotation_col = NULL,
  annotation_colors = list(celltype = c("acinar" = "salmon3",
                                        "activated" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophage" = "magenta3",
                                        "quiescent" = "salmon"),
                           sex = c("female" = "red4",
                                   "male" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black")),
  # # data.out = FALSE,
 #highlight.features = c("MA0132.2", #PDX1
  #                      "MA0148.4",  #FOXA1
   #                     "MA0874.1",  #Arx
    #                     "MA0674.1" #NKX6-1
                         #gamma
                         #ductal
                         #acinar
                         #activated
                         #quiescent
                         #endo
                         #lymphocyte
   #                      ), #macrophages
  #right_annotations = rowAnnotation(foo = anno_mark(at = c(1), labels = c("HHEX"))),
  # show_colnames = isBulk(object),
  show_rownames = FALSE,
  # scale = "row",
  cluster_row = TRUE,
  # cluster_cols = FALSE,
  # border_color = NA,
  # legend_breaks = NA,
  # drop_levels = FALSE,
  breaks=seq(-2, 2, length.out=50),
  complex = TRUE,
  #column_km = 1,
  use_raster = TRUE,
  raster_quality = 5,
  #column_split = combined_processed_rna$celltype,
  #border_color = "black",
  #gaps_col = c(15, 30, 45, 60, 75, 90, 105, 120, 135, 150),
  #gaps_row = c(1819, 3205, 5170, 5405, 10850, 14650, 21250, 23500, 25750, 26200)
  #gaps_row = c(2804, 4584, 7464, 7799, 17833, 23493, 36665, 39029, 41669, 42089)
) + rowAnnotation(mark = anno_mark(at = match(label_genes, 
                                              rownames(combined_processed_atac[genes.to.plot,])), 
                                   labels = label_genes, 
                                   which = "row",
                                   labels_gp = list(cex=0.3),
                                   #link_width = unit(4, "mm"), link_height = unit(4, "mm"),
                                   padding = 0.1))
dev.off()
dev.off()


# Plot motifs
# look at the activity of Mef2c
DefaultAssay(hm.integrated.dfree) <- "chromvar"
FeaturePlot(
  object = hm.integrated.dfree,
  features = "MA0077.1",
  min.cutoff = 0,
  max.cutoff = 1,
  cols = c("lightgrey", "red4"),
  #order = TRUE,
  raster = TRUE,
  pt.size = 1,
  raster.dpi = c(1024, 1024)
)

FeaturePlot(
  object = processed_rna,
  features = "SOX9",
  min.cutoff = 0,
  max.cutoff = 1,
  cols = c("lightgrey", "red4"),
  #order = TRUE,
  raster = TRUE,
  pt.size = 1,
  raster.dpi = c(1024, 1024)
)

MotifPlot(
  object = hm.integrated.dfree,
  motifs = "MA1608.1",
  assay = 'ATAC'
)

# Finding enriched motifs
enriched.motifs.xist <- FindMotifs(
  object = hm.integrated.dfree,
  features = c("chrX-73849537-73853078")
)

write.csv(enriched.motifs.xist, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\ORA\Motifs\xist.csv)")

# Plot
MotifPlot(
  object = hm.integrated.dfree,
  motifs = c("GLI3")
)

# Alternative to calculating 
bone <- Footprint(
  object = bone,
  motif.name = c("GATA2", "CEBPA", "EBF1"),
  genome = BSgenome.Hsapiens.UCSC.hg19
)


# Translating Motifs
motif_id <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\motif.ID.csv)", row.names = 1)
table1.df <- merge(motif_id, enriched.cvar.motifs.beta,
                   by = 'row.names', all = TRUE)
enriched.cvar.motifs.beta
# Stop



































      
system.time({
  ###Step 1: Make Pseudobulk Matrices
  #Read in final Seurat object
  adata <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
  Idents(adata) <- adata@meta.data$celltype
  samples <- unique(adata@meta.data$sample)
  
  #Pull out list of all cell types
  unique_cell_types <- unique(adata$celltype)
  DefaultAssay(adata) <- 'ATAC'
  
  #Get counts data
  da.counts <- GetAssayData(adata,slot='counts')
  
  dim(da.counts)
  head(da.counts)
  adata_matrices <- adata
  
  ##Pull out barcodes
  sample_bcs <- list()
  for (sample in samples){
    sample_bcs[[sample]] <- row.names(adata[[]][adata[[]]$sample == sample,])
  }
  
  lengths(sample_bcs)
  head(sample_bcs[[1]])
  
  #Looping through cell types by making ^ into a function
  get_per_sample_da_SUMS <- function(cell.type, mtx.fp){
    print(paste(cell.type))
    
    #pull out rows of da.counts where BC Ident matches cell.type
    bcs <- names(Idents(adata_matrices)[Idents(adata_matrices) == cell.type])
    counts <- da.counts[,colnames(da.counts) %in% bcs]
    print(dim(counts))
    
    #initialize the matrix of sample da
    counts.df <- as.data.frame(rep(0,length(row.names(da.counts))))
    row.names(counts.df) <- row.names(da.counts)
    colnames(counts.df) <- c('temp')
    
    #go through samples and calculate sum of da values
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
    mtx.fp <- sprintf('E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/%s_sample_da_total_counts.txt',cell.type) # change to save dir
    write.table(fin.counts.df,mtx.fp,sep='\t',quote=FALSE)
  }
  
  #Run function to make matrices
  unique_cell_types <- unique(adata$celltype)
  for (cell.type in unique_cell_types){
    fp = sprintf('E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/%s_pseudobulk.txt',cell.type) # change to save dir as above
    get_per_sample_da_SUMS(cell.type, fp)
  }
  
  
  ###Step 3: DESeq
  #Pseudobulk matrices directory
  dir = 'E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/'
  files = list.files(dir, pattern =".txt")
  cells = gsub("_sample_da_total_counts.txt","", files)
  
  #Create outdir for results
  outdir <- 'E:/2.SexbasedStudyCurrent/test_env/DEtesting/cells/' #changes based on analysis
  
  #Create a metadata table
  meta <- adata@meta.data[,c('sample', 'sex', 'ancestry')]
  colnames(meta) <- c('sample', 'sex', 'ancestry')
  rownames(meta) <- NULL
  meta <- meta[!duplicated(meta),]
  meta$sex_ancestry <- paste(meta$sex, meta$ancestry, sep = '_')
  
  # list of pseudobulk files
  files <- list.files('E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/', pattern='da')
  
  ##Create matrices for results
  sumres <- matrix(nrow=length(cells), ncol = 3)
  rownames(sumres) <- cells
  
  # testing for sex_ancestry_diabetes
  for (FILE in files) {
    cell <- gsub('_sample_da_total_counts.txt', '', FILE)
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    sample_names <- unique(adata@meta.data$sample)
    sample_names <- gsub('-','.', sample_names)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
    meta$Library2 <- gsub('-', '.', meta$sample)
    meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
    
    if ('male' %in% meta2$sex && 'female' %in% meta2$sex){
      print(cell)
      print('Data for 2 sex present, however not all data may be present will check this at a later step')
      
      accessibilities_to_keep <- c()
      for (i in 1:nrow(raw_counts)) {
        if (sum(raw_counts[i, ] >= 5) >= 2) {
          accessibilities_to_keep <- c(accessibilities_to_keep, rownames(raw_counts[i, ]))
        }
      }
      counts <- raw_counts[which(rownames(raw_counts) %in% accessibilities_to_keep),] 
      
      if ((length(which(meta2$sex == 'male')) > 1) | (length(which(meta2$sex == 'female')) > 1)) { #fix
        my_design <- as.formula ('~sex_ancestry') # alldata
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      }
      
      # No need for conditional formatting here
      # Specifying test combinations
      tests1 <- c('male_white', 'female_white', 'male_white', 'male_black')
      tests2 <- c('male_black', 'female_black', 'female_white', 'female_black')
      
      print('Preparing to run DESeq2')
      
      for (x in 1:length(tests1)){
        # No need for conditional formatting here
        t1 <- tests1[[x]]
        t2 <- tests2[[x]]
        test <- c('sex_ancestry', tests1[[x]],tests2[[x]]) # This should not change when you test subsetted data
        numoft1 <- length(which(meta2$sex_ancestry==t1))
        numoft2 <- length(which(meta2$sex_ancestry==t2))
        
        if (numoft1 < 3 | numoft2 < 3) {
          message(paste("!!WARNING!!"))
          message(paste(t1, "INSUFFICIENT N", sep= " "))
          message(paste('####'))
          message(paste('####'))
        } else if (numoft1 > 2 & numoft2 > 2) {
          #sprintf("%s and %s are present in the dataset", t1, t2)
          #sprintf("Find data here: %s", outdir)
          res <- results(dds, contrast=c(test), cooksCutoff=FALSE, independentFiltering=FALSE) #cooksCutoff = FALSE see here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#outlier
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














link_plot + scale_colour_viridis_d()

+ scale_fill_gradient2(low = "red",
                                 mid = "white",
                                 high = "blue",
                                 na.value = "grey50")





p1[["data"]][["Var2"]])
my_levels2 <- c("beta", "delta", "alpha", "gamma", 
                "acinar", "ductal",
                "quiescent_stellate", "activated_stellate", "endothelial",
                "lymphocyte", "macrophage")
processed_rna$disease_ancestry_lib_sex_source <- factor(x = processed_rna$disease_ancestry_lib_sex_source, levels = my_levels)

p1

predictions <- round(cor(predictions),2)
head(predictions)
melted_predictions <- melt(predictions)
head(melted_predictions)
ggplot(data = melted_predictions, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()



predictions <- t(predictions)
predictions <- as.data.frame(predictions)




predictions <- as.data.frame(cor(predictions, method = "spearman"))


predictions$x <- rownames(predictions)
predictions <- tidyr::gather(data = predictions, y, correlation, c('beta', 'alpha', 'delta', 'gamma', 'acinar', 'ductal', 'quiescent_stellate', 'activated_stellate', 'endothelial', 'macrophage'))
ggplot(predictions, aes(x, y, fill = correlation)) +
  geom_tile()

ggplot(predictions, aes(x= Var1, y = Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                                                                                  labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2






































