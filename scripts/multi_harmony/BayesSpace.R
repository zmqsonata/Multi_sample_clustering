library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(tidyverse)
library(rhdf5)
library(Matrix)
library(Seurat)
library(data.table)
library(mclust)
library(harmony)
library(stringr)
library(scater)


run_BayesSpace_multi_harmony <- function(samples, cnt_paths, xy_paths, genes_paths, cluster_num) {
    sces <- c()
    XY <- data.frame()
    xy_list <- list()
    genes_list <- list()
    N <- length(samples)
    for (i in 1:N){
    genes_list[[i]] <- read.csv(genes_paths[[i]])$genes
    }
    common_gene <- Reduce(intersect, genes_list)
    for (i in 1:N){
    cnt_mtx <- t(readMM(cnt_paths[[i]]))
    rownames(cnt_mtx) <- genes_list[[i]]
    cnt <- list(counts=cnt_mtx[rownames(cnt_mtx) %in% common_gene, ])
    sce <- SingleCellExperiment(assays = cnt)
    xy <- read.csv(xy_paths[[i]], row.names = 1)
    xy$sample_name <- samples[[i]]
    XY <- rbind(XY, xy)
    colData(sce) = cbind(colData(sce), xy[, c('sample_name','row', 'col')])
    sces <- c(sces, sce)
    }

    sce.combined <- do.call(cbind, c(sces, deparse.level = 1))
    sce.combined <- scuttle::logNormCounts(sce.combined)
    sce.combined <- scater::runPCA(sce.combined)
    sce.combined <- spatialPreprocess(sce.combined, platform = "Visium",n.PCs = 50, skip.PCA = TRUE, log.normalize = FALSE)
    sce.combined <- RunHarmony(sce.combined, "sample_name", verbose = F)
    sce.combined <- spatialCluster(sce.combined, use.dimred = "HARMONY", q = cluster_num, nrep = 10000)

    results <- list()
    clusters <- cbind(C=sce.combined$spatial.cluster, XY[, c('row', 'col', 'sample_name')])
    for (i in samples){
    df <- subset(clusters, sample_name == i)
    df$sample_name <- NULL
    results[[i]] <- df
    }
    return(results)
}