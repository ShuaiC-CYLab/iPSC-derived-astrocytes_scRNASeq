library(Seurat)
library(ggpubr)
library(clustree)
library(Cairo)
setwd('/media/black/Data4/YiRan_2nd/Matrix/Figure1/2023.09.15/')


## Read 10X Data and remove previous cluster information----
# S2_iPSCs
s2d0 <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S2_iPSCs/Single Sample/s2ipscM.rds')
s2d0@meta.data <- s2d0@meta.data[,-c(10:15)]
DefaultAssay(s2d0) <- 'RNA'
# S2D1
s2d1 <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S2D1/Single Sample/S2D1M.rds')
s2d1@meta.data <- s2d1@meta.data[,-c(10:15)]
DefaultAssay(s2d1) <- 'RNA'
# S2D3
s2d3 <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S2D3/Single Sample/S2D3M.rds')
s2d3@meta.data <- s2d3@meta.data[,-c(10:15)]
DefaultAssay(s2d3) <- 'RNA'
# S1D8
s1d8 <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S1D8/Single Sample/S1D8M.rds')
s1d8@meta.data <- s1d8@meta.data[,-c(10:15)]
DefaultAssay(s1d8) <- 'RNA'
# S2D14
s2d14 <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S2D14/Single Sample/S2D14M.rds')
s2d14@meta.data <- s2d14@meta.data[,-c(10:15)]
DefaultAssay(s2d14) <- 'RNA'
# S2D21
s2d21 <- readRDS('/media/black/Data4/YiRan_2nd/Matrix/S2D21M/Single Sample/S2D21M.rds')
s2d21@meta.data <- s2d21@meta.data[,-c(10:15)]
DefaultAssay(s2d21) <- 'RNA'


## Merge Data----
seurat_list <- list(S2D0 = s2d0, S2D1 = s2d1, S2D3 = s2d3, S1D8 = s1d8, S2D14 = s2d14, S2D21 = s2d21)
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

## Statistics----
DefaultAssay(data_merge) <- "RNA"
dir.create("./Statistics")
CairoPNG("./Statistics/Filtered Statistics.png", height = 2000, width = 3500, res = 300)
VlnPlot(data_merge, features = c("nFeature_RNA", "nCount_RNA", "mito_pct"), ncol = 3, group.by = "Day", cols = Palettes$group_pal[1:length(unique(data_merge@meta.data$Day))], pt.size = 0) +
  FeatureScatter(data_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Day") +
  FeatureScatter(data_merge, feature1 = "nCount_RNA", feature2 = "mito_pct", group.by = "Day")
dev.off()

# Number of Samples
cell_number <- as.data.frame(table(data_merge$Day))
CairoPNG("./Statistics/Cell Number.png", height = 4000, width = 3000, res = 300)
ggbarplot(cell_number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = Palettes$group_pal[1:length(unique(data_merge@meta.data$Day))],
          sort.by.groups=FALSE, 
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()


## Assess the cell cycle effect
S.Genes <- cc.genes.updated.2019$s.genes
G2M.Genes <- cc.genes.updated.2019$g2m.genes
data_merge <- CellCycleScoring(data_merge, s.features = S.Genes, g2m.features = G2M.Genes, set.ident = F)
data_merge <- RunPCA(data_merge, features = c(S.Genes, G2M.Genes))
CairoPNG("./Statistics/CellCycle1.png", height = 1500, width = 4000, res = 300)
DimPlot(data_merge, dims = c(1, 2), reduction = 'pca', group.by = 'Phase') +
  DimPlot(data_merge, dims = c(1, 3), reduction = 'pca', group.by = 'Phase') +
  DimPlot(data_merge, dims = c(2, 3), reduction = 'pca', group.by = 'Phase')
dev.off()
CairoPNG("./Statistics/CellCycle2.png", height = 3000, width = 4000, res = 300)
DimPlot(data_merge, dims = c(1, 2), reduction = 'pca', group.by = 'Phase', split.by = 'Day') +
  DimPlot(data_merge, dims = c(1, 3), reduction = 'pca', group.by = 'Phase', split.by = 'Day') +
  DimPlot(data_merge, dims = c(2, 3), reduction = 'pca', group.by = 'Phase', split.by = 'Day')
dev.off()
data_merge$CC.Diff <- data_merge$S.Score - data_merge$G2M.Score


set_resolutions <- seq(0.1, 1.0, by = 0.1)
## PCA Test
a <- data_merge
a <- RunPCA(a, npcs = 100, verbose = T)
pdf('Merge Cluster.pdf')
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
merge_res <- sapply(set_resolutions, function(x){
  p <- DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = paste0('RNA_snn_res.', x))
  print(p)
})
dev.off()

Idents(a) <- a$RNA_snn_res.0.1
pdf('./Statistics/Merge Observe Batch.pdf', width = 20)
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = 'orig.ident', split.by = 'orig.ident')
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = 'orig.ident', split.by = 'Dish')
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = 'orig.ident', split.by = 'Day')
dev.off()


## Resolution=0.1, without Harmony
data_merge <- RunPCA(data_merge, npcs = 100, verbose = T)
ElbowPlot(data_merge, ndims = 100)
pct <- data_merge[['pca']]@stdev/sum(data_merge[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
data_merge <- FindNeighbors(data_merge, dims = 1:pcs, verbose = T)
data_merge <- FindClusters(data_merge, resolution = 0.1, verbose = T) 
data_merge <- RunUMAP(data_merge, dims = 1:pcs)


## CSS----
library(simspec)
data_merge.CSS <- cluster_sim_spectrum(data_merge, label_tag = 'Day')
data_merge.CSS <- RunUMAP(data_merge.CSS, reduction = 'css', dims = 1:ncol(Embeddings(data_merge.CSS, 'css')))
data_merge.CSS <- FindNeighbors(data_merge.CSS, reduction = 'css', dims = 1:ncol(Embeddings(data_merge.CSS, 'css')))
data_merge.CSS <- FindClusters(data_merge.CSS, resolution = 1)
data_merge.CSS <- RunTSNE(data_merge.CSS, dims = 1:ncol(Embeddings(data_merge.CSS, 'css')))
saveRDS(data_merge.CSS, file = 'data_merge.CSS.rds')


## UMAP Plots----
Idents(data_merge.CSS) <- data_merge.CSS$Day
CairoPNG("UMAP1.png", height = 1600, width = 2000, res = 300)
DimPlot(data_merge.CSS, label = T) + NoLegend()
dev.off()


## Vln Plots----
CairoPNG("nFeature_RNA.png", height = 1600, width = 2000, res = 300)
VlnPlot(data_merge.CSS, features = "nFeature_RNA", pt.size = 0) + NoLegend()
dev.off()
CairoPNG("nCount_RNA.png", height = 1600, width = 2000, res = 300)
VlnPlot(data_merge.CSS, features = "nCount_RNA", pt.size = 0) + NoLegend()
dev.off()
