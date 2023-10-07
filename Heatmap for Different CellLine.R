library(Seurat)
library(Cairo)
library(do)
library(dplyr)
data_merge.CCA <- readRDS("/media/black/Data4/YiRan_2nd/Matrix/Figure2/data_merge.CCA.rds")
data_merge.CCA$CellLine <- Replace(data_merge.CCA$CellLine, from = "M", to = "DYR0100(M)")
data_merge.CCA$CellLine <- Replace(data_merge.CCA$CellLine, from = "P", to = "DYR0100")
data_merge.CCA$CellLine <- Replace(data_merge.CCA$CellLine, from = "E", to = "BIONi037-A")

## Markers----
Idents(data_merge.CCA) <- data_merge.CCA$seurat_clusters
data_merge.CCA_Markers <- FindAllMarkers(data_merge.CCA, min.pct = 0.25, logfc.threshold = 0.25,
                                         only.pos = T, assay = "RNA")
dir.create("./DEGs Heatmap")
top10 <- data_merge.CCA_Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
data4Heatmap <- subset(data_merge.CCA, downsample = 300)
CairoPNG('./DEGs Heatmap/Heatmap1 Bar.png', height = 5000, width = 3000, res = 300)
DoHeatmap(data4Heatmap, features = top10$gene) + theme(text = element_text(size = 16, face = 'bold'))
dev.off()
