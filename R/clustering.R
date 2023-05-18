#_ -------------------------------------------------------

# fun: FindClusters matrix, Cycif, CycifStack ----
#' Find clusters using Seurat functions
#'
#' @param x an object containing data to cluster
#' @param ... additional arguments passed to FindClusters
#'
#' @return an object of the same class with x
#'
#' @export
#' @importFrom Seurat FindNeighbors FindClusters
#'
#' @export
setGeneric("FindClusters", function(x,...) standardGeneric("FindClusters"))

#' Find clusters in a matrix
#'
#' This function finds clusters in a given matrix using Seurat::FindNeighbors and Seurat::FindClusters functions.
#'
#' @rdname FindClusters
#'
#' @param x A matrix object.
#' @param k.param Integer scalar. The number of nearest neighbors.
#' @param initial.membership Integer vector. The initial cluster assignment for each cell.
#' @param node.sizes Integer vector. The number of cells in each cluster.
#' @param resolution Numeric scalar. The resolution parameter for clustering.
#' @param algorithm Integer scalar. The algorithm to use for clustering.
#' @param with.labels Logical scalar. Whether to include labels.
#' @param ... Additional arguments to pass to Seurat::FindNeighbors and Seurat::FindClusters functions.
#'
#' @return An integer vector with cluster assignments for each cell.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(100), nrow = 10)
#' FindClusters(mat, k.param = 5)
#' }
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @export
setMethod("FindClusters", "matrix",
          function(x, k.param = 20, initial.membership = NULL, node.sizes = NULL,
                   resolution = 0.8, algorithm = 1, with.labels = FALSE, ...) {
            ## fin neighbors
            g <- Seurat::FindNeighbors(object = x, k.param = k.param)

            cls <- Seurat::FindClusters(object = g$snn, initial.membership = initial.membership,
                                        node.sizes = node.sizes, resolution = resolution, algorithm = algorithm)[[1]]

            return(cls)
          }
)

#' Find clusters in a Cycif object.
#'
#' @rdname FindClusters
#'
#' @param x A \code{\link{Cycif}} object.
#' @param ld_name A character scalar indicating the name of the layout to use.
#' @param k.param An integer scalar indicating the number of nearest neighbors to use for clustering.
#' @param initial.membership A vector specifying the initial cluster membership of each cell.
#' @param node.sizes A vector specifying the size of each node in the clustering tree.
#' @param resolution A numeric scalar specifying the granularity of the clustering.
#' @param algorithm An integer specifying the clustering algorithm to use.
#' @param ... Additional arguments passed to \code{\link{Seurat::FindClusters}}.
#'
#' @return A modified \code{\link{Cycif}} object with cluster assignments added to the specified layout.
#'
#' @examples
#' \dontrun{
#' data("exampleCycif")
#' clusters <- FindClusters(exampleCycif, "Layout1", k.param = 10)
#' }
#'
#' @importFrom Seurat FindClusters
#' @export
#' @method FindClusters Cycif
setMethod("FindClusters", "Cycif",
          function(x,ld_name,k.param = 20,
                   initial.membership=NULL,node.sizes=NULL,resolution=0.8,algorithm=1,...){
            call1 <- sys.calls()[[1]]
            if(missing(ld_name)){
              stop("'ld_name' should be specified.")
            }else if(!ld_name %in% ld_names(x)){
              stop("Specified 'ld_name' does not exist.")
            }

            ## subsetting the expression matrix
            ld <- ld_coords(x,ld_name)
            used.cts <- ld@used.cts
            this.abs <- ld@used.abs
            is.used <- ld@is_used

            e <- exprs(x,type="logTh_normalized")
            e1 <- data.matrix(e[is.used,this.abs])
            cls <- FindClusters(e1,
                                k.param = k.param,
                                initial.membership = initial.membership,
                                node.sizes = node.sizes,
                                resolution = resolution,
                                algorithm = algorithm)
            x@ld_coords[[ld_name]]@clusters <- cls
            x@ld_coords[[ld_name]]@clust_call <- call1
            return(x)
          }
)

#' Find Clusters
#'
#' @param x Object of class "CycifStack"
#' @param ld_name Name of the label to use for clustering
#' @param k.param Integer, number of nearest neighbors
#' @param initial.membership Vector of initial cluster assignments
#' @param node.sizes List of node sizes for cluster centers
#' @param resolution Double, resolution of the clustering algorithm
#' @param algorithm Integer, clustering algorithm to use
#' @param ... Additional arguments to be passed
#'
#' @return The input CycifStack object with the cluster assignments added
#'
#' @examples
#' \dontrun{
#' FindClusters(x, ld_name = "label_name")
#' }
#' @rdname FindClusters
#' @export
setMethod("FindClusters", "CycifStack",
          function(x,ld_name,k.param = 20,
                   initial.membership,node.sizes,resolution=0.8,algorithm=1,...){
            call1 <- sys.calls()[[1]]
            if(missing(ld_name)){
              stop("'ld_name' should be specified.")
            }else if(!ld_name %in% ld_names(x)){
              stop("Specified 'ld_name' does not exist.")
            }

            ## subsetting the expression matrix
            ld <- ld_coords(x,ld_name)
            used.cts <- ld@used.cts
            this.abs <- ld@used.abs
            is.used <- ld@is_used

            e <- exprs(x,type="logTh_normalized")
            e1 <- data.matrix(e[is.used,this.abs])
            cls <- FindClusters(e1,
                                k.param = k.param,
                                initial.membership = initial.membership,
                                node.sizes = node.sizes,
                                resolution = resolution,
                                algorithm = algorithm)
            x@ld_coords[[ld_name]]@clusters <- cls
            x@ld_coords[[ld_name]]@clust_call <- call1
            return(x)
          }
)

#_ -------------------------------------------------------

