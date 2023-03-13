#' @export
setGeneric("FindClusters", function(x,...) standardGeneric("FindClusters"))

#' @rdname FindClusters
#' @export
setMethod("FindClusters", "matrix",
          function(x,k.param=20,
                   initial.membership=NULL,node.sizes=NULL,resolution=0.8,algorithm=1,
                   with.labels=FALSE,...){
  ## fin neighbors
  g <- Seurat::FindNeighbors(
    object = x,
    k.param = k.param)

  cls <- Seurat::FindClusters(
    object = g$snn,
    initial.membership = initial.membership,
    node.sizes = node.sizes,
    resolution = resolution,
    algorithm = algorithm)[[1]]

  return(cls)
  }
)

#' @rdname FindClusters
#' @export
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

