library(BASS)
library(Matrix)
library(Seurat)

run_Bass_Simple <- function(dir, file_name, cluster_num) {
    set.seed(1)

    data <- Seurat::Load10X_Spatial(
      data.dir = dir,
      filename = file_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE,
      to.upper = FALSE
    )

    xy <- list(as.matrix(data@images[["slice1"]]@coordinates[,c("row","col")]))
    cnts <- list(as.matrix(data@assays$Spatial$counts))

    BASS <- createBASSObject(cnts, xy, C = 15, R = cluster_num, beta_method = "SW")
    BASS <- BASS.preprocess(BASS,doLogNormalize = TRUE, doPCA = TRUE, scaleFeature = TRUE, nPC = 20)
    BASS <- BASS.run(BASS)
    BASS <- BASS.postprocess(BASS)
    
    xy <- data@images[["slice1"]]@coordinates[,c("row","col")]
    clusters <- data.frame(cbind(RowNames = rownames(xy), C=BASS@results$z[[1]]))
    rownames(clusters) <- clusters$RowNames
    clusters$RowNames <- NULL
    results <- merge(clusters, xy, by="row.names", all=FALSE)
    rownames(results) <- results$Row.names
    results$Row.names <- NULL

    return(results)
}
