library(Seurat)
library(SeuratData)
library(spruce)
library(mclust)
library(ggplot2)
library(data.table)
library(Matrix)


run_MAPLE_multi <- function(samples, cnt_paths, xy_paths, genes_paths, cluster_num) {
    genes_list <- list()
    cnt_cb <- c()
    XY <- data.frame()
    N <- length(samples)
    for (i in 1:N){
      xy <- read.csv(xy_paths[[i]], row.names = 1)
      xy$sample_name <- samples[[i]]
      XY <- rbind(XY, xy)
      genes_list[[i]] <- read.csv(genes_paths[i])$genes
    }
    common_gene <- Reduce(intersect, genes_list)
    for (i in 1:N){
      cnt_mtx <- t(readMM(cnt_paths[[i]]))
      rownames(cnt_mtx) <- genes_list[[i]]
      cnt_cb <- cbind(cnt_cb, cnt_mtx[rownames(cnt_mtx) %in% common_gene, ])
    }
    colnames(cnt_cb) <- rownames(XY)

    data = CreateSeuratObject(counts = cnt_cb, assay="Spatial")
    data@images$image =  new(
      Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = XY[, c('row', 'col')]
    )

    data <- SCTransform(data, assay = "Spatial", verbose = FALSE)
    d1 <- RunPCA(data)
    d1_fit_PCs <- fit_spruce(d1, K=cluster_num, emb = "PCs", MCAR = FALSE, CAR = FALSE)

    results <- list()
    clusters <- cbind(C=d1_fit_PCs$z, XY[, c('row', 'col', 'sample_name')])
    for (i in samples){
      df <- subset(clusters, sample_name == i)
      df$sample_name <- NULL
      results[[i]] <- df
    }
    return(results)
}