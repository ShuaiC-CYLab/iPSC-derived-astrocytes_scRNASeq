data_merge <- merge(s2ipsc, y=c(s2d1,s2d3,s1d8,s2d14,s2d21m,s1d21p,s1d21e), project = "Unfiltered")
library(Seurat)
library(Cairo)
dir.create("./Unfiltered nFeature nCount Mito")

## nFeature_RNA----
Idents(data_merge) <- data_merge$Day
my_pal <- c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3","#F564E3","#F564E3")
CairoPNG("./Unfiltered nFeature nCount Mito/nFeautre_RNA.png", height = 1500, width = 2500, res = 300)
VlnPlot(data_merge, features = "nFeature_RNA", pt.size = 0, cols = my_pal)
dev.off()

## nCount_RNA----
CairoPNG("./Unfiltered nFeature nCount Mito/nCount_RNA.png", height = 1500, width = 2500, res = 300)
VlnPlot(data_merge, features = "nCount_RNA", pt.size = 0, cols = my_pal)
dev.off()

## mito_pct----
CairoPNG("./Unfiltered nFeature nCount Mito/mito_pct.png", height = 1500, width = 2500, res = 300)
VlnPlot(data_merge, features = "mito_pct", pt.size = 0, cols = my_pal)
dev.off()

