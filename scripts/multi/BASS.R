library(BASS)
library(Matrix)
library(Seurat)
library(data.table)
library(tidyverse)


run_Bass_multi <- function(samples, cnt_paths, xy_paths, genes_paths, cluster_num) {

    cnt_list <- list()
    genes_list <- list()
    xy_list <- list()
    XY <- data.frame()
    N <- length(samples)

    for (i in 1:N){
      xy <- read.csv(xy_paths[[i]], row.names = 1)
      xy$sample_name <- samples[[i]]
      XY <- rbind(XY, xy)
      xy_list[[i]] <- xy[, c('row', 'col')]
      cnt_mtx <- t(readMM(cnt_paths[[i]]))
      rownames(cnt_mtx) <- read.csv(genes_paths[[i]])$genes
      colnames(cnt_mtx) <- rownames(xy)
      cnt_list[[i]] <- cnt_mtx
    }

    BASS <- createBASSObject(cnt_list, xy_list, C = 15, R = cluster_num, beta_method = "SW", burnin = 5000)
    BASS <- BASS.preprocess(BASS,doLogNormalize = TRUE, doPCA = TRUE, scaleFeature = TRUE, nPC = 20, doBatchCorrect = FALSE)
    BASS <- BASS.run(BASS)
    BASS <- BASS.postprocess(BASS)

    results <- list()
    for (i in 1:N){
      df <- subset(XY[, c('row', 'col', 'sample_name')], sample_name == samples[[i]])
      df <- cbind(C=BASS@results$z[[i]], df)
      df$sample_name <- NULL
      results[[samples[[i]]]] <- df
    }
    return(results)
}
