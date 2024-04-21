library(Seurat)
library(stardust)
library(aricode)
library(stringr)


run_Sardust_multi <- function(samples, cnt_paths, xy_paths, genes_paths, cluster_num) {
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

    output <- autoStardust(countMatrix = data.frame(cnt_cb), spotPositions = XY[, c('row', 'col')],
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