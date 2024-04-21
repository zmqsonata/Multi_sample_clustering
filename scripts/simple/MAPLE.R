library(Seurat)
library(SeuratData)
library(spruce)
library(mclust)
library(ggplot2)
library(maple)

suppressPackageStartupMessages(require(optparse))

run_MAPLE_Simple <- function(dir, file_name, cluster_num){

  data <- Seurat::Load10X_Spatial(
      data.dir = dir,
      filename = file_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE,
      to.upper = FALSE
    )

  data.t <- SCTransform(data, assay = "Spatial", verbose = FALSE)
  
  d1 <- RunPCA(data.t)

  d1_fit_PCs <- fit_spruce(d1, K=cluster_num, emb = "PCs", MCAR = FALSE, CAR = FALSE)

  xy <- data@images[["slice1"]]@coordinates[,c("row","col")]
  clusters <- data.frame(cbind(RowNames = rownames(xy), C=d1_fit_PCs$z))
  rownames(clusters) <- clusters$RowNames
  clusters$RowNames <- NULL
  results <- merge(clusters, xy, by="row.names", all=FALSE)
  rownames(results) <- results$Row.names
  results$Row.names <- NULL
  return(results)
}