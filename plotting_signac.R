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

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')

# Run the following code once you have Seurat installed
suppressWarnings({
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

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma", 
               "acinar", "ductal",
               "quiescent_stellate", "activated_stellate", "endothelial",
               "lymphocyte", "macrophage")

table(hm.integrated.dfree$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree$celltype <- factor(x = hm.integrated.dfree$celltype, levels = my_levels)
table(unique(hm.integrated.dfree$celltype))

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


# Corelation plot
# # Pseudobulk
# Idents(hm.integrated.dfree) <- "celltype"
# hm.integrated.dfree$celltype_sex <- paste(hm.integrated.dfree$celltype, hm.integrated.dfree$sex, sep = "_")
# table(hm.integrated.dfree@meta.data$celltype_sex)
# Idents(hm.integrated.dfree) <- "celltype_sex"
# DefaultAssay(hm.integrated.dfree) <- "ATAC"
# combined_processed_atac <- Seurat:::PseudobulkExpression(object = hm.integrated.dfree, 
#                                                          pb.method = 'aggregate', 
#                                                          return.seurat = TRUE,
#                                                          slot = 'counts')
# 
# DefaultAssay(combined_processed_atac) <- "ATAC"
# combined_processed_atac <- RunTFIDF(combined_processed_atac, assay = "ATAC")
# 
# {
#   combined_processed_atac$celltype_sex <- combined_processed_atac@active.ident
#   Idents(combined_processed_atac) <- 'celltype_sex'
#   combined_processed_atac$celltype <- combined_processed_atac$orig.ident
#   metadat <- combined_processed_atac@meta.data
#   metadat <- metadat %>% 
#     mutate(celltype_sex = str_replace(celltype_sex, "activated", "activated-stellate"))
#   metadat <- metadat %>% 
#     mutate(celltype_sex = str_replace(celltype_sex, "quiescent", "quiescent-stellate"))
#   metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex, "_", -1)
#   combined_processed_atac@meta.data = metadat
# }
# 
# table(combined_processed_atac@meta.data[["celltype"]])
# table(combined_processed_atac@meta.data[["sex"]])
# 
# # cluster re-assignment occurs, which re-assigns clustering in my_levels
# my_levels <- c("delta", "beta", "alpha", "gamma",
#                "ductal", "acinar",
#                "activated", "quiescent", "endothelial",
#                "lymphocyte", "macrophage") 
# 
# table(combined_processed_atac$celltype)
# 
# # Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
# combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
# table(combined_processed_atac$celltype)
# Idents(combined_processed_atac) <- "celltype"
# DefaultAssay(combined_processed_atac) <- "ATAC"

# split object by cluster, take counts and aggregate by sum (pseudobulking)
processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

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
# DefaultAssay(hm.integrated.dfree) <- "RNA"
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
# 
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
# pheatmap(cmat, clustering_distance_rows = "euclidean", cluster_rows = FALSE,cluster_cols = FALSE)
# 
# rna.seurat <- processed_rna[, sample(colnames(processed_rna), size =500, replace=F)]
# atac.signac <- hm.integrated.dfree[, sample(colnames(hm.integrated.dfree), size =500, replace=F)]
# qsave(rna.seurat, r"(E:\downsampling_collab\rna.seurat.qs)")
# qsave(atac.signac, r"(E:\downsampling_collab\atac.signac.qs)")
# 
# 
# 
# # cluster re-assignment occurs, which re-assigns clustering in my_levels
# my_levels <- c("beta", "delta", "alpha", "gamma", 
#                "acinar", "ductal",
#                "quiescent_stellate", "activated_stellate", "endothelial",
#                "lymphocyte", "macrophage")
# 
# table(hm.integrated.dfree$celltype)




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


# Plotting Chromin accessibility
DefaultAssay(hm.integrated.dfree) <- "ATAC"
p1 <-  CoveragePlot(hm.integrated.dfree, region = c("INS"), 
                          window = 100,
                          ymax = 200,
             links = FALSE,
             #tile = TRUE,
             extend.upstream = 20000,
             extend.downstream = 20000) & scale_fill_manual(values = c("dodgerblue3",      #beta
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

regions.to.plot <- uniqueregions
label_genes <- c("INS", "GCG", "SST") #uniquegenes



DefaultAssay(combined_processed_atac) <- "ATAC"
pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 8,
    height = 6)
dittoHeatmap(
  object = combined_processed_atac,#(subset(combined_processed_rna, idents = c("alpha"))),
  genes = regions.to.plot, #this is a compete gene set
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
  show_rownames = FALSE,
  # scale = "row",
  cluster_row = FALSE,
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
  gaps_col = c(15, 30, 45, 60, 75, 90, 105, 120, 135, 150),
  gaps_row = c(1819, 3205, 5170, 5405, 10850, 14650, 21250, 23500, 25750, 26200)
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

c(1819, 3205, 5170, 5405, 10850, 14650, 21250, 23500, 25750, 26200)











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






































