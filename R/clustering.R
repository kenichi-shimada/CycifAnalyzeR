#_ -------------------------------------------------------

# fun: LdClustering matrix, Cycif, CycifStack ----

#' Find clusters using Seurat functions
#'
#' @param x An matrix, Cycif, or CycifStack object containing data to cluster.
#' @param ld_name A character scalar indicating the name of the layout to use.
#' @param k.param An integer specifying the number of nearest neighbors.
#' @param initial.membership A vector specifying the initial cluster membership of each cell.
#' @param node.sizes A vector specifying the size of each node in the clustering tree.
#' @param resolution A numeric scalar specifying the granularity of the clustering.
#' @param algorithm An integer specifying the clustering algorithm to use.
#' @param with.labels Logical scalar indicating whether to include labels.
#' @param method Community-detection method passed to \code{Seurat::FindClusters}. Default "igraph".
#' @param ... Additional arguments passed to the clustering function.
#'
#' @return
#' A modified object with cluster assignments added to the specified layout (for a `matrix` object).
#' An object of the same class as `x` with cluster assignments added (for a `Cycif` or `CycifStack` object).
#'
#' @importFrom Seurat FindNeighbors FindClusters
#'
#' @export
setGeneric("LdClustering", function(x,...) standardGeneric("LdClustering"))

#' @rdname LdClustering
#' @export
setMethod("LdClustering", "matrix",
          function(x, k.param = 20, initial.membership = NULL, node.sizes = NULL,
                   resolution = 0.8, algorithm = 1, with.labels = FALSE, method="igraph", ...) {
            ## fin neighbors
            g <- Seurat::FindNeighbors(object = x, k.param = k.param)

            cls <- Seurat::FindClusters(object = g$snn, initial.membership = initial.membership,
                                        node.sizes = node.sizes, resolution = resolution, algorithm = algorithm,
                                        method = method)[[1]] #,...

            return(cls)
          }
)

#' @rdname LdClustering
#' @export
setMethod("LdClustering", "Cycif",
          function(x,ld_name,k.param = 20,
                   initial.membership,node.sizes,resolution=0.8,algorithm=1,method="igraph", ...){
            .Defunct(msg=paste(
              "LdClustering(Cycif) has had zero call sites anywhere across",
              "CycifAnalyzeR or any downstream project repo (EP, TALAVE, HR-APM).",
              "Use subsetCells() + clusterCells() on a CellFeatures object instead",
              "(clusterCells() itself still uses LdClustering(\"matrix\", ...) as",
              "its engine, which remains live)."
            ))
          }
)

#' @rdname LdClustering
#' @export
setMethod("LdClustering", "CycifStack",
          function(x,ld_name,k.param = 20,
                   initial.membership,node.sizes,resolution=0.8,algorithm=1,method="igraph",...){
            .Defunct(msg=paste(
              "LdClustering(CycifStack) has had zero call sites anywhere across",
              "CycifAnalyzeR or any downstream project repo (EP, TALAVE, HR-APM).",
              "Use subsetCells() + clusterCells() on a CellFeatures object instead",
              "(clusterCells() itself still uses LdClustering(\"matrix\", ...) as",
              "its engine, which remains live)."
            ))
          }
)

#_ -------------------------------------------------------

