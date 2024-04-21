library(Seurat)
library(scater)
library(Matrix)
library(aricode)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)


run_Seurat_multi_harmony <- function(samples, cnt_paths, xy_paths, genes_paths, cluster_num) {
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
    colnames(cnt_cb) <- rownames(XY)

    SeuratObj = CreateSeuratObject(counts = cnt_cb, assay="Spatial")
    SeuratObj@images$image =  new(
      Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = XY[, c('row', 'col')]
    )
    SeuratObj@meta.data$orig.ident <- XY$sample_name
    SeuratObj <- SCTransform(SeuratObj, assay = "Spatial", verbose = FALSE)
    SeuratObj <- RunPCA(SeuratObj, assay = "SCT", verbose = FALSE,npcs = 30)
    SeuratObj <- RunHarmony(SeuratObj, group.by.vars = "orig.ident", verbose = T)
    SeuratObj <- FindNeighbors(SeuratObj, reduction = "harmony", dims = 1:20)
    SeuratObj <- FindClusters(SeuratObj, verbose = FALSE)

    results <- list()
    clusters <- cbind(C=SeuratObj@meta.data$seurat_clusters, XY[, c('row', 'col', 'sample_name')])
    for (i in samples){
      df <- subset(clusters, sample_name == i)
      df$sample_name <- NULL
      results[[i]] <- df
    }
    return(results)
}