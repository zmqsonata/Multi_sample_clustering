library("Seurat")
library("stardust")
library(aricode)

run_Stardust_Simple <- function(dir, file_name, cluster_num){

data <- Seurat::Load10X_Spatial(
      data.dir = dir,
      filename = file_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE,
      to.upper = FALSE
    )

xy <- data@images[["slice1"]]@coordinates[,c("row","col")]
spotPositions <- data@images[["slice1"]]@coordinates
spotPositions=spotPositions[spotPositions$tissue=='1',]
countMatrix <- ReadMtx(
  mtx = file.path(dir, "filtered_feature_bc_matrix/matrix.mtx.gz"), 
  features = file.path(dir, "filtered_feature_bc_matrix/features.tsv.gz"),
  cells = file.path(dir, "filtered_feature_bc_matrix/barcodes.tsv.gz")
)
countMatrix=as.data.frame(countMatrix)


output <- autoStardust(countMatrix = countMatrix, spotPositions = spotPositions,
                      pcaDimensions=10, res=0.8)

clusters <- data.frame(cbind(RowNames = rownames(output@meta.data), C=output@active.ident))
rownames(clusters) <- clusters$RowNames
clusters$RowNames <- NULL
results <- merge(clusters, xy, by="row.names", all=FALSE)
rownames(results) <- results$Row.names
results$Row.names <- NULL

return(results)
}