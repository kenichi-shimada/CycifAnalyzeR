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
#' @param ... Additional arguments passed to the clustering function.
#'
#' @return
#' A modified object with cluster assignments added to the specified layout (for a `matrix` object).
#' An object of the same class as `x` with cluster assignments added (for a `Cycif` or `CycifStack` object).
#'
#' @export
setGeneric("LdClustering", function(x,...) standardGeneric("LdClustering"))

#' @rdname LdClustering
#' @export
setMethod("LdClustering", "matrix",
          function(x, k.param = 20, initial.membership = NULL, node.sizes = NULL,
                   resolution = 0.8, algorithm = 1, with.labels = FALSE, ...) {
            ## fin neighbors
            g <- Seurat::FindNeighbors(object = x, k.param = k.param)

            cls <- Seurat::FindClusters(object = g$snn, initial.membership = initial.membership,
                                        node.sizes = node.sizes, resolution = resolution, algorithm = algorithm)[[1]]

            return(cls)
          }
)

#' @rdname LdClustering
#' @export
setMethod("LdClustering", "Cycif",
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
            norm_type <- ld@norm_type

            e <- exprs(x,type=norm_type)
            e1 <- data.matrix(e[is.used,this.abs])
            cls <- LdClustering(e1,
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

#' @rdname LdClustering
#' @export
setMethod("LdClustering", "CycifStack",
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
            norm_type <- ld@norm_type

            e <- exprs(x,type=norm_type)
            e1 <- data.matrix(e[is.used,this.abs])
            cls <- LdClustering(e1,
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

