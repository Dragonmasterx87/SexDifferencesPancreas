#Plotting data

# Cell population
cluster_order_y <- match(levels(pancreas.combined.h.s@meta.data[["celltype"]]), metaLevels("celltype", pancreas.combined.h.s))
dittoBarPlot(pancreas.combined.h.s, var = "celltype", group.by = "sample", var.labels.reorder = cluster_order_y, x.reorder = c(1,15,6,7,
                                                                                                                               2,4,8,9,
                                                                                                                               3,10,13,14,
                                                                                                                               5,11,12)) +
  scale_fill_manual(values = c("darkturquoise",
                               "lightgreen",
                               "springgreen4",
                               "lightgoldenrod3",
                               "green3",
                               "grey56",
                               "grey80",
                               "deeppink",
                               "violet",
                               "purple",
                               "coral2",
                               "magenta",
                               "red4",
                               "black",
                               "red")) + coord_flip()
                               
# Cell population
cluster_order_y <- match(levels(hm.integrated.dfree@meta.data[["celltype"]]), metaLevels("celltype", hm.integrated.dfree))
dittoBarPlot(hm.integrated.dfree, var = "celltype", group.by = "sample", var.labels.reorder = cluster_order_y, x.reorder = c(1,15,6,7,
                                                                                                                               2,4,8,9,
                                                                                                                               3,10,13,14,
                                                                                                                               5,11,12)) +
  scale_fill_manual(values = c("darkturquoise",
                               "lightgreen",
                               "springgreen4",
                               "lightgoldenrod3",
                               "grey56",
                               "grey80",
                               "deeppink",
                               "violet",
                               "coral2",
                               "magenta",
                               "red4")) + coord_flip()


# Plotting venn diagrams
beta.bmvsbf <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\bmvbf\alldata\beta.bmvsbf.csv)", row.names = 1)
beta.wfvsbf <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\wfvbf\alldata\beta.wfvsbf.csv)", row.names = 1)
beta.wmvsbm <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\wmvsbm\alldata\beta.wmvsbm.csv)", row.names = 1)
beta.wmvswf <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\wmvwf\All data\beta.wmvswf.csv)", row.names = 1)

sig_df_up <- dplyr::filter(beta.bmvsbf, p_val < 0.05 & avg_log2FC > 0.26303) # >1.2x
beta.bmvsbf_up <- rownames(sig_df_up)

sig_df_up <- dplyr::filter(beta.wfvsbf, p_val < 0.05 & avg_log2FC > 0.26303) # >1.2x
beta.wfvsbf_up <- rownames(sig_df_up)

sig_df_up <- dplyr::filter(beta.wmvsbm, p_val < 0.05 & avg_log2FC > 0.26303) # >1.2x
beta.wmvsbm_up <- rownames(sig_df_up)

sig_df_up <- dplyr::filter(beta.wmvswf, p_val < 0.05 & avg_log2FC > 0.26303) # >1.2x
beta.wmvswf_up <- rownames(sig_df_up)

x <- list(
  beta.bmvsbf = beta.bmvsbf_up,
  beta.wfvsb = beta.wfvsbf_up,
  beta.wmvsbm = beta.wmvsbm_up,
  beta.wmvswf = beta.wmvswf_up
  )

venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "darkturquoise")+
  scale_color_manual(values = c("beta.bmvsbf" = "black",
                                "beta.wfvsb" ="black",
                                "beta.wmvsbm" = 'black', 
                                "beta.wmvswf" = 'black'),
                     labels = c('D' = 'D = bdiv_human'))+
  theme_void()


# Plotting venn diagrams
beta.bmvsbf <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\bmvbf\alldata\beta.bmvsbf.csv)", row.names = 1)
beta.wfvsbf <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\wfvbf\alldata\beta.wfvsbf.csv)", row.names = 1)
beta.wmvsbm <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\wmvsbm\alldata\beta.wmvsbm.csv)", row.names = 1)
beta.wmvswf <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\Sequencing_Data\scRNAseq\updated analysis\DGE\wmvwf\All data\beta.wmvswf.csv)", row.names = 1)

sig_df_up <- dplyr::filter(beta.bmvsbf, p_val < 0.05 & avg_log2FC < -0.32192) # >1.2x
beta.bmvsbf_up <- rownames(sig_df_up)

sig_df_up <- dplyr::filter(beta.wfvsbf, p_val < 0.05 & avg_log2FC < -0.32192) # >1.2x
beta.wfvsbf_up <- rownames(sig_df_up)

sig_df_up <- dplyr::filter(beta.wmvsbm, p_val < 0.05 & avg_log2FC < -0.32192) # >1.2x
beta.wmvsbm_up <- rownames(sig_df_up)

sig_df_up <- dplyr::filter(beta.wmvswf, p_val < 0.05 & avg_log2FC < -0.32192) # >1.2x
beta.wmvswf_up <- rownames(sig_df_up)

x <- list(
  beta.bmvsbf = beta.bmvsbf_up,
  beta.wfvsb = beta.wfvsbf_up,
  beta.wmvsbm = beta.wmvsbm_up,
  beta.wmvswf = beta.wmvswf_up
)

venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "darkturquoise")+
  scale_color_manual(values = c("beta.bmvsbf" = "black",
                                "beta.wfvsb" ="black",
                                "beta.wmvsbm" = 'black', 
                                "beta.wmvswf" = 'black'),
                     labels = c('D' = 'D = bdiv_human'))+
  theme_void()












                        
