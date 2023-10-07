library(Seurat)
library(WriteXLS)
library(scCustomize)
library(dplyr)
setwd("/media/black/Data4/YiRan_2nd/Matrix/Figure1")
data_merge.CSS <- readRDS("/media/black/Data4/YiRan_2nd/Matrix/Figure1/data_merge.CSS.rds")
Idents(data_merge.CSS) <- data_merge.CSS$Day
CSS_Markers <- FindAllMarkers(data_merge.CSS, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", 
                              only.pos = T)
WriteXLS(CSS_Markers, ExcelFileName = "CSS Markers Day.xlsx")

## DEGs Heatmap----
top10 <- CSS_Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
data4Heatmap <- subset(data_merge.CSS, downsample = 300)
dir.create("./DEGs Heatmap")
CairoPNG('./DEGs Heatmap/Heatmap1.png', height = 5000, width = 3000, res = 300)
DoHeatmap(data4Heatmap, features = top10$gene) + NoLegend() + theme(text = element_text(size = 16, face = 'bold'))
dev.off()
