# This script is for Seurat integration of the M1 - M5 mouse data.

# Load libraries
library(Seurat) # v3.2.2
library(dplyr)
library(ggplot2)

# Load data
rds <- readRDS('data/RDSobjects/AllMice_integrated.rds')

# Load metadata that has the Dendritic cells removed
meta <- read.csv('data/CSVfiles/IntegratedMouse_CD8s_Only_withDendriticRemoved_meta.csv',
                 row.names = 1)

# Filter cells
cells_to_keep <- rds$seurat_clusters != 5
rds <- rds[, cells_to_keep]


#############################
#                           #
#       Run Integration     #
#                           #
#############################

# Split by sample and run SCTransform
# 'so' is 'Seurat object'
so <- SplitObject(rds, split.by = 'sample')
so <- lapply(so, FUN = SCTransform, return.only.var.genes = FALSE)

print(names(so))
## [1] "18-520-1" "18-520-4" "19-603-1" "19-603-2" "19-603-3" "19-603-5" "19-603-6" "19-603-7"

# Find integration features, run PrepSCTIntegration
features <- SelectIntegrationFeatures(so, nfeatures = Inf)
so <- PrepSCTIntegration(so, anchor.features = features, verbose = TRUE)

# Run PCA
so <- lapply(X = so, FUN = RunPCA, verbose = TRUE, features = features)

anchors <- FindIntegrationAnchors(object.list = so,
                                  normalization.method = "SCT",
                                  reduction = "rpca",
                                  anchor.features = features,
                                  verbose = TRUE)

# Integrate
integrated <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT",
                            verbose = TRUE)

# Print the data
print(integrated)
# An object of class Seurat 
# 25889 features across 10811 samples within 3 assays 
# Active assay: integrated (6083 features, 6083 variable features)
# 2 other assays present: RNA, SCT


##################################
#                                #
#      Run UMAP and Analysis     #
#                                #
##################################

DefaultAssay(integrated) <- "integrated"

# Run PCA and UMAP
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

# Find neighbors, clustering, and UMAP
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.4)

# Add the UMAP coordinates to the metadata
integrated$UMAP_1 <- integrated@reductions$umap@cell.embeddings[, 1]
integrated$UMAP_2 <- integrated@reductions$umap@cell.embeddings[, 2]

# Print the data
print(integrated)
# An object of class Seurat 
# 30772 features across 10889 samples within 3 assays 
# Active assay: integrated (6034 features, 6034 variable features)
# 2 other assays present: RNA, SCT
# 2 dimensional reductions calculated: pca, umap

DimPlot(integrated, split.by = 'sample', ncol = 3)

# Plot the UMAP
pdf('results/UMAPs/Integrated_ClustersUMAP.pdf')
DimPlot(integrated, label = TRUE, label.size = 6) +
  ggtitle("All Samples (Integrated)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Plot genes
genes_to_plot <- c("CD4", "CD3D", "HAVCR2", "IFNG",
                   "CD3E", "CD3G", "PDCD1", "CTLA4",
                   "GZMA", "GZMB", "MKI67", "SELL",
                   "TCF7")


pdf("results/UMAPs/Integrated_genes.pdf", width = 9, height = 8)
DefaultAssay(integrated) <- 'SCT'
FeaturePlot(integrated, features = genes_to_plot[9:16], ncol = 4, slot = 'data')
dev.off()

# Plot clone sizes
pdf('results/UMAPs/Integrated_CloneSizeUMAP.pdf')
FeaturePlot(integrated, features = 'CloneSize') +
  ggtitle("All Samples (Integrated) Clone Size") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed()
dev.off()


# Plot AB status
pdf('results/UMAPs/Integrated_ABstatus.pdf')
DimPlot(integrated, split.by = 'AB_status', ncol = 2) +
  ggtitle("All Samples (Integrated) AB Status") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed()
dev.off()

# Plot by sample
pdf('results/UMAPs/Integrated_bySample.pdf')
DimPlot(integrated, split.by = 'sample', ncol = 3) +
  ggtitle("All Samples (Integrated)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed()
dev.off()

# Bar plot of number of samples per cluster
pdf('results/UMAPs/Integrated_ClusterSizesBarplot.pdf')
p <- integrated@meta.data %>% group_by(seurat_clusters) %>% tally() %>%
  ggplot(aes(x = seurat_clusters, y = n, fill = seurat_clusters)) +
  geom_bar(stat = "identity") + geom_text(aes(label = n)) +
  theme_bw() + ggtitle("Cluster Sizes")
print(p)
dev.off()

# Bar plot of number of samples in data
pdf('results/UMAPs/Integrated_SampleSizesBarplot.pdf')
p <- integrated@meta.data %>% group_by(sample) %>% tally() %>%
  ggplot(aes(x = sample, y = n, fill = sample)) +
  geom_bar(stat = "identity") + geom_text(aes(label = n)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Sample Sizes")
print(p)
dev.off()

# Find DE genes per cluster
DefaultAssay(integrated) <- "RNA"
DEgenes <- FindAllMarkers(integrated, max.cells.per.ident = 1000)
write.csv(DEgenes, file = 'results/DEgenes_Integrated.csv')

# Save Seurat object
saveRDS(integrated, 'data/RDSobjects/Integrated.rds')

sessionInfo()

integrated <- readRDS('data/RDSobjects/Integrated.rds')

# R version 4.0.5 (2021-03-31)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.3      dplyr_1.0.5        Seurat_3.2.2       SeuratObject_4.0.0
# [5] renv_0.13.2       
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-152         matrixStats_0.58.0   RcppAnnoy_0.0.18     RColorBrewer_1.1-2  
# [5] httr_1.4.2           sctransform_0.3.2    tools_4.0.5          utf8_1.2.1          
# [9] R6_2.5.0             irlba_2.3.3          rpart_4.1-15         KernSmooth_2.23-18  
# [13] uwot_0.1.10          mgcv_1.8-34          DBI_1.1.1            lazyeval_0.2.2      
# [17] colorspace_2.0-0     withr_2.4.1          tidyselect_1.1.0     gridExtra_2.3       
# [21] compiler_4.0.5       cli_2.3.1            plotly_4.9.3         labeling_0.4.2      
# [25] scales_1.1.1         lmtest_0.9-38        spatstat.data_2.1-0  ggridges_0.5.3      
# [29] pbapply_1.4-3        spatstat_1.64-1      goftest_1.2-2        stringr_1.4.0       
# [33] digest_0.6.27        spatstat.utils_2.1-0 pkgconfig_2.0.3      htmltools_0.5.1.1   
# [37] parallelly_1.24.0    limma_3.44.3         fastmap_1.1.0        htmlwidgets_1.5.3   
# [41] rlang_0.4.10         rstudioapi_0.13      shiny_1.6.0          farver_2.1.0        
# [45] generics_0.1.0       zoo_1.8-9            jsonlite_1.7.2       ica_1.0-2           
# [49] magrittr_2.0.1       patchwork_1.1.1      Matrix_1.3-2         Rcpp_1.0.6          
# [53] munsell_0.5.0        fansi_0.4.2          abind_1.4-5          reticulate_1.18     
# [57] lifecycle_1.0.0      stringi_1.5.3        MASS_7.3-53.1        Rtsne_0.15          
# [61] plyr_1.8.6           grid_4.0.5           parallel_4.0.5       listenv_0.8.0       
# [65] promises_1.2.0.1     ggrepel_0.9.1        crayon_1.4.1         deldir_0.2-10       
# [69] miniUI_0.1.1.1       lattice_0.20-41      cowplot_1.1.1        splines_4.0.5       
# [73] tensor_1.5           pillar_1.5.1         igraph_1.2.6         future.apply_1.7.0  
# [77] reshape2_1.4.4       codetools_0.2-18     leiden_0.3.7         glue_1.4.2          
# [81] data.table_1.14.0    vctrs_0.3.7          png_0.1-7            httpuv_1.5.5        
# [85] gtable_0.3.0         RANN_2.6.1           purrr_0.3.4          polyclip_1.10-0     
# [89] tidyr_1.1.3          future_1.21.0        assertthat_0.2.1     rsvd_1.0.3          
# [93] mime_0.10            xtable_1.8-4         RSpectra_0.16-0      later_1.1.0.1       
# [97] survival_3.2-10      viridisLite_0.3.0    tibble_3.1.0         cluster_2.1.1       
# [101] globals_0.14.0       fitdistrplus_1.1-3   ellipsis_0.3.1       ROCR_1.0-11       
