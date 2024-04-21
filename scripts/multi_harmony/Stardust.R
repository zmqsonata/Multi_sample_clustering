library(Seurat)
# library(stardust)
library(aricode)
library(stringr)
library(harmony)


autoStardust <- function(countMatrix,
                         spotPositions,
                         pcaDimensions=10,
                         res=0.8){
  sample_names <- spotPositions$sample_name
  spotPositions <- spotPositions[, c('row', 'col')]
  countMatrix <- countMatrix[,sort(colnames(countMatrix))]
  d <- dim(spotPositions)[2]
  spotPositions <- spotPositions[,(d-1):d]
  spotPositions <- spotPositions[sort(rownames(spotPositions)),]
  
  seuratObject <- Seurat::CreateSeuratObject(countMatrix)
  SeuratObj@meta.data$orig.ident <- sample_names
  seuratObject <- suppressWarnings(Seurat::SCTransform(seuratObject, assay = "RNA", verbose = FALSE))
  seuratObject <- Seurat::RunPCA(seuratObject, assay = "SCT", verbose = FALSE)
  seuratObject <- RunHarmony(SeuratObj, group.by.vars = "orig.ident", verbose = F)
  if(pcaDimensions <= 2){
    pcaDimensions = 2
  }
  m <- seuratObject@reductions[["harmony"]]@cell.embeddings[,1:pcaDimensions]
  distPCA = dist(m,method="minkowski",p=2)
  distCoord <- dist(spotPositions,method="minkowski",p=2)
  distCoord <- distCoord*(max(distPCA)/max(distCoord))
  expr_norm <- (distPCA - min(distPCA)) / (max(distPCA) - min(distPCA))
  distCoord <- (distCoord)*(as.double(as.vector(expr_norm)))
  finalDistance <- as.matrix(distPCA + distCoord)
  neighbors <- suppressWarnings(Seurat::FindNeighbors(finalDistance))
  neighbors <- list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
  seuratObject@graphs <- neighbors
  seuratObject <- suppressWarnings(Seurat::FindClusters(seuratObject,  resolution = res,
                                                        verbose = FALSE, graph.name = "neighbors_snn"))
  seuratObject
}



run_Sardust_multi_harmony <- function(samples, cnt_paths, xy_paths, genes_paths, cluster_num) {
    genes_list <- list()
    cnt_cb <- c()
    XY <- data.frame()
    N <- length(samples)
    for (i in 1:N){
      xy <- read.csv(xy_paths[[i]], row.names = 1)
      xy$sample_name <- samples[[i]]
      XY <- rbind(XY, xy)
      genes_list[[i]] <- read.csv(genes_paths[[i]])$genes
    }
    common_gene <- Reduce(intersect, genes_list)
    for (i in 1:N){
      cnt_mtx <- t(readMM(cnt_paths[[i]]))
      rownames(cnt_mtx) <- genes_list[[i]]
      cnt_cb <- cbind(cnt_cb, cnt_mtx[rownames(cnt_mtx) %in% common_gene, ])
    }

    output <- autoStardust(countMatrix = data.frame(cnt_cb), spotPositions = XY[, c('row', 'col', 'sample_name')],
                          pcaDimensions=10, res=0.8)

    results <- list()
    clusters <- cbind(C=output@active.ident, XY[, c('row', 'col', 'sample_name')])
    for (i in samples){
      df <- subset(clusters, sample_name == i)
      df$sample_name <- NULL
      results[[i]] <- df
    }
    return(results)
}