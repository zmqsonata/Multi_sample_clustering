library(dplyr)
library(Seurat)
library(patchwork)


run_Seurat_Simple <- function(dir, file_name, cluster_num){

data <- Seurat::Load10X_Spatial(
      data.dir = dir,
      filename = file_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE,
      to.upper = FALSE
    )

xy <- data@images[["slice1"]]@coordinates[,c("row","col")]
pbmc <- CreateSeuratObject(counts = data@assays$Spatial$counts, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes

# top10 <- head(VariableFeatures(pbmc), 10)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#
# DimPlot(pbmc, reduction = "pca")
#
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

#knn
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

clusters <- data.frame(cbind(RowNames = rownames(pbmc@meta.data), C=pbmc@meta.data$seurat_clusters))
rownames(clusters) <- clusters$RowNames
clusters$RowNames <- NULL
results <- merge(clusters, xy, by="row.names", all=FALSE)
rownames(results) <- results$Row.names
results$Row.names <- NULL
return(results)
}