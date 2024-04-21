library(iSC.MEB)
library(Seurat)


run_ISC_MEB_multi <- function(samples, cnt_paths, xy_paths, genes_paths, cluster_num) {
    obj_list <- list()
    XY <- data.frame()
    N <- length(samples)
    for (i in 1:N){
      xy <- read.csv(xy_paths[[i]], row.names = 1)
      xy$sample_name <- samples[[i]]
      XY <- rbind(XY, xy)
      cnt_mtx <- t(readMM(cnt_paths[[i]]))
      rownames(cnt_mtx) <- read.csv(genes_paths[[i]])$genes
      colnames(cnt_mtx) <- rownames(xy)
      SeuratObj = CreateSeuratObject(counts = cnt_mtx)
      SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
      SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)
      SeuratObj@meta.data=merge(SeuratObj@meta.data,xy[, c('row', 'col')],by = 'row.names')
      rownames(SeuratObj@meta.data)=SeuratObj@meta.data$Row.names
      obj_list <- append(obj_list, SeuratObj)
    }

    iSCMEBObj <- CreateiSCMEBObject(seuList = obj_list, verbose = FALSE, gene.number=2000, premin.spots = 0, postmin.spots = 0,selectGenesMethod = c("SPARK-X"))
    iSCMEBObj <- CreateNeighbors(iSCMEBObj, platform = "Visium")
    iSCMEBObj <- runPCA(iSCMEBObj, npcs = 15, pca.method = "APCA")
    iSCMEBObj <- SetModelParameters(iSCMEBObj, verbose = TRUE, maxIter = 12)
    iSCMEBObj <- iSCMEB(iSCMEBObj, K = cluster_num)

    results <- list()
    for (i in 1:N){
      df <- subset(XY[, c('row', 'col', 'sample_name')], sample_name == samples[[i]])
      df <- df[rownames(iSCMEBObj@seulist[[i]]@meta.data), ]
      df <- cbind(C=idents(iSCMEBObj)[[i]], df)
      df$sample_name <- NULL
      results[[samples[[i]]]] <- df
    }
    return(results)
}