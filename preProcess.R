# ' @description: proProcess and cluster
# Load Packages
library(Seurat)
library(Cairo)
library(ggplot2)
library(ggpubr)
library(DoubletFinder)

## iPSCs_Day0
s2ipsc.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S2_iPSCs/Matrix/')
colnames(s2ipsc.data) <- paste('S2_iPSC', colnames(s2ipsc.data), sep = '_')
s2ipsc <- CreateSeuratObject(s2ipsc.data, project = 'S2_iPSCs', min.cells = 5, min.features = 300)
s2ipsc[['mito_pct']] <- PercentageFeatureSet(s2ipsc, pattern = '^MT-')
s2ipsc[['Day']] <- 'Day0'
s2ipsc[['CellLine']] <- 'M'
s2ipsc[['Dish']] <- 'S2'
s2ipsc <- subset(s2ipsc, subset = nCount_RNA > quantile(s2ipsc$nCount_RNA, c(0.025)) & 
                 nCount_RNA < quantile(s2ipsc$nCount_RNA, c(0.975)) & 
                 nFeature_RNA > quantile(s2ipsc$nFeature_RNA, c(0.025)) & 
                 nFeature_RNA < quantile(s2ipsc$nFeature_RNA, c(0.975)) & 
                 mito_pct < quantile(s2ipsc$mito_pct, c(0.975)))
s2ipsc.pro <- SCTransform(s2ipsc, vars.to.regress = c('mito_pct'), verbose = F)
s2ipsc.pro <- RunPCA(s2ipsc.pro, npcs = 50, verbose = F)
ElbowPlot(object = s2ipsc.pro, ndims = 50)
pct <- s2ipsc.pro[['pca']]@stdev/sum(s2ipsc.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s2ipsc.pro <- FindNeighbors(s2ipsc.pro, dims = 1:pcs, verbose = F)
s2ipsc.pro <- FindClusters(s2ipsc.pro, resolution = 0.1, verbose = F)
s2ipsc.pro  <- RunUMAP(s2ipsc.pro, dims = 1:pcs)
s2ipsc.pro <- doubletDetect(s2ipsc.pro, PCs = 1:pcs, doublet.rate = 0.116, annotation = 'SCT_snn_res.0.1', sct = T)
s2ipsc.pro.singlet <- subset(s2ipsc.pro, subset = Doublet == 'Singlet')
s2ipsc.pro.singlet <- RunTSNE(s2ipsc.pro.singlet, dims = 1:pcs)


## Day1_Mono
s2d1.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S2D1/')
colnames(s2d1.data) <- paste('S2_D1', colnames(s2d1.data), sep = '_')
s2d1 <- CreateSeuratObject(s2d1.data, project = 'S2D1M', min.cells = 5, min.features = 300)
s2d1[['mito_pct']] <- PercentageFeatureSet(s2d1, pattern = '^MT-')
s2d1[['Day']] <- 'Day1'
s2d1[['CellLine']] <- 'M'
s2d1[['Dish']] <- 'S2'
s2d1 <- subset(s2d1, subset = nCount_RNA > quantile(s2d1$nCount_RNA, c(0.025)) & 
                 nCount_RNA < quantile(s2d1$nCount_RNA, c(0.975)) & 
                 nFeature_RNA > quantile(s2d1$nFeature_RNA, c(0.025)) & 
                 nFeature_RNA < quantile(s2d1$nFeature_RNA, c(0.975)) & 
                 mito_pct < quantile(s2d1$mito_pct, c(0.975)))
s2d1.pro <- SCTransform(s2d1, vars.to.regress = c('mito_pct'), verbose = F)
s2d1.pro <- RunPCA(s2d1.pro, npcs = 50, verbose = F)
ElbowPlot(object = s2d1.pro, ndims = 50)
pct <- s2d1.pro[['pca']]@stdev/sum(s2d1.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s2d1.pro <- FindNeighbors(s2d1.pro, dims = 1:pcs, verbose = F)
s2d1.pro <- FindClusters(s2d1.pro, resolution = 0.1, verbose = F)
s2d1.pro  <- RunUMAP(s2d1.pro, dims = 1:pcs)
s2d1.pro <- doubletDetect(s2d1.pro, PCs = 1:pcs, doublet.rate = 0.087, annotation = 'SCT_snn_res.0.1', sct = T)
s2d1.pro.singlet <- subset(s2d1.pro, subset = Doublet == 'Singlet')
s2d1.pro.singlet <- RunTSNE(s2d1.pro.singlet, dims = 1:pcs)

## Day3_Mono
s2d3.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S2D3/')
colnames(s2d3.data) <- paste('S2_D3', colnames(s2d3.data), sep = '_')
s2d3 <- CreateSeuratObject(s2d3.data, project = 'S2D3M', min.cells = 5, min.features = 300)
s2d3[['mito_pct']] <- PercentageFeatureSet(s2d3, pattern = '^MT-')
s2d3[['Day']] <- 'Day3'
s2d3[['CellLine']] <- 'M'
s2d3[['Dish']] <- 'S2'
s2d3 <- subset(s2d3, subset = nCount_RNA > quantile(s2d3$nCount_RNA, c(0.025)) & 
                 nCount_RNA < quantile(s2d3$nCount_RNA, c(0.975)) & 
                 nFeature_RNA > quantile(s2d3$nFeature_RNA, c(0.025)) & 
                 nFeature_RNA < quantile(s2d3$nFeature_RNA, c(0.975)) & 
                 mito_pct < quantile(s2d3$mito_pct, c(0.975)))
s2d3.pro <- SCTransform(s2d3, vars.to.regress = c('mito_pct'), verbose = F)
s2d3.pro <- RunPCA(s2d3.pro, npcs = 50, verbose = F)
ElbowPlot(object = s2d3.pro, ndims = 50)
pct <- s2d3.pro[['pca']]@stdev/sum(s2d3.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s2d3.pro <- FindNeighbors(s2d3.pro, dims = 1:pcs, verbose = F)
s2d3.pro <- FindClusters(s2d3.pro, resolution = 0.1, verbose = F)
s2d3.pro  <- RunUMAP(s2d3.pro, dims = 1:pcs)
s2d3.pro <- doubletDetect(s2d3.pro, PCs = 1:pcs, doublet.rate = 0.099, annotation = 'SCT_snn_res.0.1', sct = T)
s2d3.pro.singlet <- subset(s2d3.pro, subset = Doublet == 'Singlet')
s2d3.pro.singlet <- RunTSNE(s2d3.pro.singlet, dims = 1:pcs)

## Day8_Mono
s1d8.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S1D8/')
colnames(s1d8.data) <- paste('S1_D8', colnames(s1d8.data), sep = '_')
s1d8 <- CreateSeuratObject(s1d8.data, project = 'S1D8M', min.cells = 5, min.features = 300)
s1d8[['mito_pct']] <- PercentageFeatureSet(s1d8, pattern = '^MT-')
s1d8[['Day']] <- 'Day8'
s1d8[['CellLine']] <- 'M'
s1d8[['Dish']] <- 'S1'
s1d8 <- subset(s1d8, subset = nCount_RNA > quantile(s1d8$nCount_RNA, c(0.025)) & 
                 nCount_RNA < quantile(s1d8$nCount_RNA, c(0.975)) & 
                 nFeature_RNA > quantile(s1d8$nFeature_RNA, c(0.025)) & 
                 nFeature_RNA < quantile(s1d8$nFeature_RNA, c(0.975)) & 
                 mito_pct < quantile(s1d8$mito_pct, c(0.975)))
s1d8.pro <- SCTransform(s1d8, vars.to.regress = c('mito_pct'), verbose = F)
s1d8.pro <- RunPCA(s1d8.pro, npcs = 50, verbose = F)
ElbowPlot(object = s1d8.pro, ndims = 50)
pct <- s1d8.pro[['pca']]@stdev/sum(s1d8.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s1d8.pro <- FindNeighbors(s1d8.pro, dims = 1:pcs, verbose = F)
s1d8.pro <- FindClusters(s1d8.pro, resolution = 0.3, verbose = F)
s1d8.pro  <- RunUMAP(s1d8.pro, dims = 1:pcs)
s1d8.pro <- doubletDetect(s1d8.pro, PCs = 1:pcs, doublet.rate = 0.072, annotation = 'SCT_snn_res.0.3', sct = T)
s1d8.pro.singlet <- subset(s1d8.pro, subset = Doublet == 'Singlet')
s1d8.pro.singlet <- RunTSNE(s1d8.pro.singlet, dims = 1:pcs)

## Day14_Mono
s2d14.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S2D14/')
colnames(s2d14.data) <- paste('S2_D14', colnames(s2d14.data), sep = '_')
s2d14 <- CreateSeuratObject(s2d14.data, project = 'S2D14M', min.cells = 5, min.features = 300)
s2d14[['mito_pct']] <- PercentageFeatureSet(s2d14, pattern = '^MT-')
s2d14[['Day']] <- 'Day14'
s2d14[['CellLine']] <- 'M'
s2d14[['Dish']] <- 'S2'
s2d14 <- subset(s2d14, subset = nCount_RNA > quantile(s2d14$nCount_RNA, c(0.025)) & 
                  nCount_RNA < quantile(s2d14$nCount_RNA, c(0.975)) & 
                  nFeature_RNA > quantile(s2d14$nFeature_RNA, c(0.025)) & 
                  nFeature_RNA < quantile(s2d14$nFeature_RNA, c(0.975)) & 
                  mito_pct < quantile(s2d14$mito_pct, c(0.975)))
s2d14.pro <- SCTransform(s2d14, vars.to.regress = c('mito_pct'), verbose = F)
s2d14.pro <- RunPCA(s2d14.pro, npcs = 50, verbose = F)
ElbowPlot(object = s2d14.pro, ndims = 50)
pct <- s2d14.pro[['pca']]@stdev/sum(s2d14.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s2d14.pro <- FindNeighbors(s2d14.pro, dims = 1:pcs, verbose = F)
s2d14.pro <- FindClusters(s2d14.pro, resolution = 0.1, verbose = F)
s2d14.pro  <- RunUMAP(s2d14.pro, dims = 1:pcs)
s2d14.pro <- doubletDetect(s2d14.pro, PCs = 1:pcs, doublet.rate = 0.074, annotation = 'SCT_snn_res.0.1', sct = T)
s2d14.pro.singlet <- subset(s2d14.pro, subset = Doublet == 'Singlet')
s2d14.pro.singlet <- RunTSNE(s2d14.pro.singlet, dims = 1:pcs)

## Day21_Mono
s2d21m.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S2D21M/')
colnames(s2d21m.data) <- paste('S2_D21', colnames(s2d21m.data), sep = '_')
s2d21m <- CreateSeuratObject(s2d21m.data, project = 'S2D21M', min.cells = 5, min.features = 300)
s2d21m[['mito_pct']] <- PercentageFeatureSet(s2d21m, pattern = '^MT-')
s2d21m[['Day']] <- 'Day21'
s2d21m[['CellLine']] <- 'M'
s2d21m[['Dish']] <- 'S2'
s2d21m <- subset(s2d21m, subset = nCount_RNA > quantile(s2d21m$nCount_RNA, c(0.025)) & 
                   nCount_RNA < quantile(s2d21m$nCount_RNA, c(0.975)) & 
                   nFeature_RNA > quantile(s2d21m$nFeature_RNA, c(0.025)) & 
                   nFeature_RNA < quantile(s2d21m$nFeature_RNA, c(0.975)) & 
                   mito_pct < quantile(s2d21m$mito_pct, c(0.975)))
s2d21m.pro <- SCTransform(s2d21m, vars.to.regress = c('mito_pct'), verbose = F)
s2d21m.pro <- RunPCA(s2d21m.pro, npcs = 50, verbose = F)
ElbowPlot(object = s2d21m.pro, ndims = 50)
pct <- s2d21m.pro[['pca']]@stdev/sum(s2d21m.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s2d21m.pro <- FindNeighbors(s2d21m.pro, dims = 1:pcs, verbose = F)
s2d21m.pro <- FindClusters(s2d21m.pro, resolution = 0.1, verbose = F)
s2d21m.pro  <- RunUMAP(s2d21m.pro, dims = 1:pcs)
s2d21m.pro <- doubletDetect(s2d21m.pro, PCs = 1:pcs, doublet.rate = 0.055, annotation = 'SCT_snn_res.0.1', sct = T)
s2d21m.pro.singlet <- subset(s2d21m.pro, subset = Doublet == 'Singlet')
s2d21m.pro.singlet <- RunTSNE(s2d21m.pro.singlet, dims = 1:pcs)

## Day21_E
s1d21e.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S1D21E/')
colnames(s1d21e.data) <- paste('E', colnames(s1d21e.data), sep = '_')
s1d21e <- CreateSeuratObject(s1d21e.data, project = 'S1D21E', min.cells = 5, min.features = 300)
s1d21e[['mito_pct']] <- PercentageFeatureSet(s1d21e, pattern = '^MT-')
s1d21e[['Day']] <- 'Day21'
s1d21e[['CellLine']] <- 'E'
s1d21e[['Dish']] <- 'S1'
s1d21e <- subset(s1d21e, subset = nCount_RNA > quantile(s1d21e$nCount_RNA, c(0.025)) & 
                   nCount_RNA < quantile(s1d21e$nCount_RNA, c(0.975)) & 
                   nFeature_RNA > quantile(s1d21e$nFeature_RNA, c(0.025)) & 
                   nFeature_RNA < quantile(s1d21e$nFeature_RNA, c(0.975)) & 
                   mito_pct < quantile(s1d21e$mito_pct, c(0.975)))
s1d21e.pro <- SCTransform(s1d21e, vars.to.regress = c('mito_pct'), verbose = F)
s1d21e.pro <- RunPCA(s1d21e.pro, npcs = 50, verbose = F)
ElbowPlot(object = s1d21e.pro, ndims = 50)
pct <- s1d21e.pro[['pca']]@stdev/sum(s1d21e.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s1d21e.pro <- FindNeighbors(s1d21e.pro, dims = 1:pcs, verbose = F)
s1d21e.pro <- FindClusters(s1d21e.pro, resolution = 0.1, verbose = F)
s1d21e.pro  <- RunUMAP(s1d21e.pro, dims = 1:pcs)
s1d21e.pro <- doubletDetect(s1d21e.pro, PCs = 1:pcs, doublet.rate = 0.107, annotation = 'SCT_snn_res.0.1', sct = T)
s1d21e.pro.singlet <- subset(s1d21e.pro, subset = Doublet == 'Singlet')
s1d21e.pro.singlet <- RunTSNE(s1d21e.pro.singlet, dims = 1:pcs)

## Day21_P
s1d21p.data <- Read10X('/media/black/Data4/YiRan_2nd/Matrix/S1D21P/')
colnames(s1d21p.data) <- paste('P', colnames(s1d21p.data), sep = '_')
s1d21p <- CreateSeuratObject(s1d21p.data, project = 'S1D21P', min.cells = 5, min.features = 300)
s1d21p[['mito_pct']] <- PercentageFeatureSet(s1d21p, pattern = '^MT-')
s1d21p[['Day']] <- 'Day21'
s1d21p[['CellLine']] <- 'P'
s1d21p[['Dish']] <- 'S1'
s1d21p <- subset(s1d21p, subset = nCount_RNA > quantile(s1d21p$nCount_RNA, c(0.025)) & 
                   nCount_RNA < quantile(s1d21p$nCount_RNA, c(0.975)) & 
                   nFeature_RNA > quantile(s1d21p$nFeature_RNA, c(0.025)) & 
                   nFeature_RNA < quantile(s1d21p$nFeature_RNA, c(0.975)) & 
                   mito_pct < quantile(s1d21p$mito_pct, c(0.975)))
s1d21p.pro <- SCTransform(s1d21p, vars.to.regress = c('mito_pct'), verbose = F)
s1d21p.pro <- RunPCA(s1d21p.pro, npcs = 50, verbose = F)
ElbowPlot(object = s1d21p.pro, ndims = 50)
pct <- s1d21p.pro[['pca']]@stdev/sum(s1d21p.pro[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs <- min(co1, co2)
pcs
s1d21p.pro <- FindNeighbors(s1d21p.pro, dims = 1:pcs, verbose = F)
s1d21p.pro <- FindClusters(s1d21p.pro, resolution = 0.2, verbose = F)
s1d21p.pro  <- RunUMAP(s1d21p.pro, dims = 1:pcs)
s1d21p.pro <- doubletDetect(s1d21p.pro, PCs = 1:pcs, doublet.rate = 0.089, annotation = 'SCT_snn_res.0.2', sct = T)
s1d21p.pro.singlet <- subset(s1d21p.pro, subset = Doublet == 'Singlet')
s1d21p.pro.singlet <- RunTSNE(s1d21p.pro.singlet, dims = 1:pcs)