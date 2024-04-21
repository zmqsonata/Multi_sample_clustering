library(Seurat)
library(SeuratData)
library(ggplot2)
library(dplyr)
library(iSC.MEB)

run_ISC_MEB_Simple <- function(dir, file_name, cluster_num) {
options(future.globals.maxSize = 8000 * 1024^2)

data <- Seurat::Load10X_Spatial(
      data.dir = dir,
      filename = file_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE,
      to.upper = FALSE
    )

xy <- data@images[["slice1"]]@coordinates[,c("row","col")]
brain1 = CreateSeuratObject(counts = data@assays$Spatial$counts)
brain1 <- NormalizeData(brain1, normalization.method = "LogNormalize", scale.factor = 10000)
brain1 <- NormalizeData(brain1)
brain1 <- FindVariableFeatures(brain1, selection.method = "vst", nfeatures = 2000)

brain1@meta.data=merge(brain1@meta.data,xy,by = 'row.names')
rownames(brain1@meta.data)=brain1@meta.data$Row.names


brain.merge <- list(brain1)

iSCMEBObj <- CreateiSCMEBObject(seuList = brain.merge, verbose = FALSE,gene.number=2000, premin.spots = 0, postmin.spots = 0,selectGenesMethod = c("SPARK-X"))

iSCMEBObj <- CreateNeighbors(iSCMEBObj, platform = "Visium")

iSCMEBObj <- runPCA(iSCMEBObj, npcs = 15, pca.method = "APCA")

iSCMEBObj <- SetModelParameters(iSCMEBObj, verbose = TRUE,maxIter = 12)


iSCMEBObj <- iSCMEB(iSCMEBObj, K = cluster_num)

clusters <- data.frame(cbind(RowNames = rownames(iSCMEBObj@seulist[[1]]@meta.data),
                             C=idents(iSCMEBObj)[[1]]))
rownames(clusters) <- clusters$RowNames
pos <- data.frame(iSCMEBObj@resList@posList)
rownames(pos) <- clusters$RowNames
clusters$RowNames <- NULL
results <- merge(clusters, pos, by="row.names", all=FALSE)
rownames(results) <- results$Row.names
results$Row.names <- NULL
return(results)
}