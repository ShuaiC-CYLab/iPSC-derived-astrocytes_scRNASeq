library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(future)
library(Cairo)
set.seed(123)
setwd('/media/black/Data4/YiRan_2nd/Matrix/Figure2')

## Read 10X Data and remove previous cluster information----
# S1D21P
P <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S1D21P/Single Sample/S1D21P.rds')
P@meta.data <- P@meta.data[,-c(10:15)]
DefaultAssay(P) <- 'RNA'

# S1D21E
E <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S1D21E/Single Sample/S1D21E.rds')
E@meta.data <- E@meta.data[,-c(10:15)]
DefaultAssay(E) <- 'RNA'

# S2D21M
M <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S2D21M/Single Sample/S2D21M.rds')
M@meta.data <- M@meta.data[,-c(10:15)]
DefaultAssay(M) <- 'RNA'


## Merge Data
seurat_list <- list(M = M, P = P, E = E)
seurat_list_Standard <- variableFeatureSelection(seurat_list, method = 'Standard', nfeatures = 3000)
saveRDS(seurat_list_Standard, file = "Standard.rds")
seurat_list_SCT <- variableFeatureSelection(seurat_list, method = "SCT", nfeatures = 3000)
saveRDS(seurat_list_SCT, file = "SCT.rds")

# assay=SCT
data_merge <- merge(seurat_list_SCT[[1]], y = seurat_list_SCT[2:length(seurat_list_SCT)], project = "Time Course")
DefaultAssay(data_merge) <- "SCT"
seurat_features_SCT <- SelectIntegrationFeatures(seurat_list_SCT, nfeatures = 3000)
VariableFeatures(data_merge) <- seurat_features_SCT

# assay=RNA
seurat_features_RNA <- SelectIntegrationFeatures(seurat_list_Standard, nfeatures = 3000)
DefaultAssay(data_merge) <- "RNA"
VariableFeatures(data_merge) <- seurat_features_RNA
data_merge <- NormalizeData(data_merge, verbose = FALSE)
data_merge <- ScaleData(data_merge, verbose = FALSE, vars.to.regress = c("mito_pct"), features = rownames(data_merge@assays$RNA@data))
DefaultAssay(data_merge) <- "SCT"
saveRDS(data_merge, file = "data_merge.rds")


## Cell Cycle
S.Gene <- cc.genes.updated.2019$s.genes
G2M.Gene <- cc.genes.updated.2019$g2m.genes
data_merge <- CellCycleScoring(data_merge, s.features = S.Gene, g2m.features = G2M.Gene, set.ident = F)
data_merge <- RunPCA(data_merge, features = c(S.Gene, G2M.Gene))
data_merge$CC.Diff <- data_merge$S.Score - data_merge$G2M.Score
CairoPDF('CellCycle.pdf')
DimPlot(data_merge, dims = c(1, 2), reduction = 'pca', group.by = 'Phase')
DimPlot(data_merge, dims = c(1, 3), reduction = 'pca', group.by = 'Phase')
DimPlot(data_merge, dims = c(2, 3), reduction = 'pca', group.by = 'Phase')
dev.off()


set_resolutions <- seq(0.1, 1, by = 0.1)
## PCA Test
a <- data_merge
a <- RunPCA(a, npcs = 100, verbose = T)
CairoPDF('PCA Resolution.pdf')
ElbowPlot(object = a, ndims = 100)
pct <- a[['pca']]@stdev/sum(a[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
pcs <- ifelse(pcs %% 10 < 5, pcs - pcs %% 10 + 5, pcs + 10 - pcs %% 10)
a <- FindNeighbors(a, dims = 1:pcs, verbose = T)
a <- FindClusters(object = a, resolution = set_resolutions, verbose = T) 
clustree(a)
a <- RunUMAP(a, dims = 1:pcs)
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = 'orig.ident')
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = 'Phase')
merge_res <- sapply(set_resolutions, function(x){
  p <- DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = paste0('SCT_snn_res.', x))
  print(p)
})
dev.off()


## CCA
seurat_list.CCA <- PrepSCTIntegration(seurat_list_SCT, anchor.features = seurat_features_SCT)
CCA_Anchors <-  FindIntegrationAnchors(seurat_list.CCA, normalization.method = "SCT", anchor.features = seurat_features_SCT, verbose = FALSE)
data_merge.CCA <- IntegrateData(CCA_Anchors, normalization.method = "SCT", verbose = FALSE)
data_merge.CCA <- ScaleData(data_merge.CCA)
data_merge.CCA <- RunPCA(data_merge.CCA, npcs = 100)
pct <- data_merge.CCA[['pca']]@stdev/sum(data_merge.CCA[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
pcs <- ifelse(pcs %% 10 < 5, pcs - pcs %% 10 + 5, pcs + 10 - pcs %% 10)
data_merge.CCA <- RunUMAP(data_merge.CCA, reduction = "pca", dims = 1:pcs, verbose = FALSE)
data_merge.CCA <- RunTSNE(data_merge.CCA, reduction = "pca", dims = 1:pcs, verbose = FALSE)
data_merge.CCA <- FindNeighbors(data_merge.CCA, reduction = "pca", dims = 1:pcs, verbose = FALSE)
data_merge.CCA <- FindClusters(data_merge.CCA, resolution = 0.1, verbose = FALSE)
saveRDS(data_merge.CCA, file = "data_merge.CCA.rds")
