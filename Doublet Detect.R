doubletDetect <- function(Seurat.object, PCs, doublet.rate = 0.076, annotation = "SCT_snn_res.1.1", pN_value = 0.25, GT = FALSE, sct = FALSE){
  
  library(DoubletFinder) #Require cleanup of low-quality cells in advance
  object.sample <- paramSweep_v3(Seurat.object, PCs = PCs, sct = sct)
  object.gt.calls <- summarizeSweep(object.sample, GT = GT)
  
  #Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value tested during doubletFinder_ParamSweep. 
  #Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions. 
  #If ground-truth doublet classifications are available, BCmvn is plotted along with mean ROC AUC for each pK.
  sweep.object <- find.pK(object.gt.calls)
  pK_value <- as.numeric(as.character(sweep.object$pK[sweep.object$BCmetric == max(sweep.object$BCmetric)])) #计算最优pK
  
  par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  plot(x = as.numeric(as.character(sweep.object$pK)), y = sweep.object$BCmetric, pch = 16,type="b", col = "blue",lty=1)
  abline(v=pK_value,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_value,max(sweep.object$BCmetric),as.character(pK_value),pos = 4,col = "red")
  
  #potential doublet rate，
  nExp_poi <- round(doublet.rate*nrow(Seurat.object@meta.data))
    annotations <- Seurat.object@meta.data[, annotation]
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    pANN_value <- paste0("pANN_", pN_value, "_", pK_value, '_', nExp_poi)
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = sct)
    label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi.adj)
  Seurat.object@meta.data$Doublet <- Seurat.object@meta.data[, label]
  
  return(Seurat.object)
}
