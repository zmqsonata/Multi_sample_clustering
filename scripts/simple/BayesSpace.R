library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(tidyverse)
library(rhdf5)
library(Matrix)
library(Seurat)
library(assertthat)
library(mclust)

run_BayesSpace_Simple <- function(dir, file_name, cluster_num) {
    set.seed(1)

    data <- Seurat::Load10X_Spatial(
      data.dir = dir,
      filename = file_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE,
      to.upper = FALSE
    )

    data = Seurat::NormalizeData(data)
    data = Seurat::FindVariableFeatures(data)
    data = Seurat::ScaleData(data)
    data = Seurat::RunPCA(data)

    diet.data = Seurat::DietSeurat(data, graphs = "pca") #slim down Seurat obj prior to conversion
    sce = as.SingleCellExperiment(diet.data) #convert seurat to SCE
    colData(sce) = cbind(colData(sce), data@images$slice1@coordinates)

    q <- cluster_num
    d <- 15

    sce = spatialPreprocess(sce, platform = "Visium", skip.PCA = T, log.normalize = F)
    sce <- spatialCluster(sce, q=q, d=d, platform='Visium', nrep=50000, gamma=3, save.chain=TRUE, model = "t")
    labels <- sce$spatial.cluster
    
    xy <- data@images[["slice1"]]@coordinates[,c("row","col")]
    clusters <- data.frame(cbind(RowNames = rownames(xy), C=sce$spatial.cluster))
    rownames(clusters) <- clusters$RowNames
    clusters$RowNames <- NULL
    results <- merge(clusters, xy, by="row.names", all=FALSE)
    rownames(results) <- results$Row.names
    results$Row.names <- NULL
    return(results)
}