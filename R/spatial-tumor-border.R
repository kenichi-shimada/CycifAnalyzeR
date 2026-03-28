#_ -------------------------------------------------------

# fun: defineTumorBorder Cycif ----

#' @title Define tumor bed regions within each Cycif dataset.
#'
#' @description This function identifies tumor bed regions within each Cycif dataset and assigns cells to these regions based on their spatial coordinates and cell types.
#'
#' @param x A Cycif or CycifStack object.
#' @param strict Logical, specifying whether to perform strict cell type calling (see \code{defineCellTypes}).
#' @param dth Numeric, the distance threshold for identifying adjacent cells.
#' @param parallel Logical, indicating whether to use parallel processing for edge point calculation in \code{concaveman}.
#' @param min.pts Integer, the minimum number of points required to form a cluster.
#' @param n.cores Integer, the number of CPU cores to use for parallel processing.
#' @param ct_name Character, the name of the cell type column (used in \code{cell_types}).
#' @param verbose Logical, indicating whether to display progress messages.
#'
#' @return A list containing information about tumor bed regions, including distance threshold, cell coordinates, regions within tumor bed, edge points, cell types, clusters, and neighbors.
#'
#' @export
#' @importFrom dbscan dbscan
#' @importFrom concaveman concaveman
#' @importFrom sp point.in.polygon
#' @importFrom parallel mclapply
#' @importFrom RColorBrewer brewer.pal
#'
#' @details
#' To define tumor bed regions within CyCIF data, the following steps are performed:
#' 1. If the `dth` parameter is not provided, it calculates an optimal distance threshold based on the MajorAxisLength of cells.
#' 2. Adjusts the Y-coordinate values to ensure proper orientation.
#' 3. Identifies cell types within the dataset, allowing for strict or non-strict matching.
#' 4. Separates tumor cells from non-tumor cells based on their cell types.
#' 5. Performs DBSCAN clustering on tumor cells to identify clusters of adjacent cells.
#' 6. Determines neighbors using NN (Fast Regular Nearest Neighbors) only for tumor cells.
#' 7. Calculates edge points for the tumor clusters using the concaveman algorithm.
#' 8. Determines whether cells are within the tumor bed regions based on edge point polygons.
#' 9. Returns a list containing various information about the tumor bed regions and assigned cells.
#'
#' @seealso
#' \code{\link{cyApply}} \code{\link{defineCellTypes}} \code{\link{concaveman}}
#'
#' @export
setGeneric("defineTumorBorder", function(x,...) standardGeneric("defineTumorBorder"))

#' @rdname defineTumorBorder
#' @export
setMethod("defineTumorBorder", "Cycif",
  function(x,strict=FALSE,dth,parallel=TRUE,min.pts=3,n.cores=1,ct_name="default",verbose=TRUE,cancer.cts = "Tumor"){
    require(dbscan)
    require(RColorBrewer)
    require(concaveman)
    require(sp)
    if(verbose){
      cat("processing defineTumorBorder() for ",names(x),"...\n")
    }

    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid

    this.cts <- cell_types(x,strict=strict,ct_name=ct_name)$cell_types
    if(any(levels(this.cts)!="NA")){
      this.cts <- factor(this.cts,levels=c(levels(this.cts),"NA"))
    }
    this.cts[is.na(this.cts)] <- "NA"

    is.tumor <- !is.na(this.cts) & this.cts %in% cancer.cts
    xy.tumor <- xy[is.tumor,]

    if(missing(dth)){
      majl <- x@segment_property$MajorAxisLength
      m1sd <- mean(majl) + sd(majl)*2
      dth <- m1sd * 2 # 35
    }

    dbs1 <- dbscan::dbscan(xy.tumor, eps=dth, minPts = min.pts, weights = NULL, borderPoints = TRUE)
    cls <- dbs1$cluster

    frt <- dbscan::frNN(xy.tumor,eps=dth)
    fr.ids <- sapply(seq(frt$id),function(i){
      tmp <- frt$id[[i]]
      tmp <- tmp[tmp > i]
    })

    is.in <- edge.pts <- list()

    nls <- unique(cls)
    nls <- sort(nls[nls!=0])

    if(parallel){
      edge.pts <- parallel::mclapply(seq(nls),function(i){
        xyt <- xy.tumor[cls==nls[i],]
        edge.pt <- as.data.frame(concaveman(as.matrix(xyt), concavity = .8, length_threshold = 0))
        return(edge.pt)
      },mc.cores=n.cores)
      is.in <- parallel::mclapply(seq(edge.pts),function(i){
        is.in <- point.in.polygon(xy$X,xy$Y,edge.pts[[i]]$V1,edge.pts[[i]]$V2)==1
        return(is.in)
      },mc.cores=n.cores)
    }else{
      for(i in seq(nls)){
        xyt <- xy.tumor[cls==nls[i],]
        edge.pts[[i]] <- as.data.frame(concaveman(as.matrix(xyt), concavity = .8, length_threshold = 0))
        is.in[[i]] <- point.in.polygon(xy2$X,xy2$Y,edge.pts[[i]]$V1,edge.pts[[i]]$V2)==1
      }
    }

    is.in.1 <- rowSums(do.call(cbind,is.in))>0

    lst <- list(
      dist.th = dth,
      min.pts = min.pts,
      xy = xy,
      within.tumor.bed = is.in.1,
      edge.pts = edge.pts,
      cell.types = this.cts,
      clusters = cls,
      neighbors = fr.ids
    )
    return(lst)
  }
)

#' @rdname defineTumorBorder
#' @export
setMethod("defineTumorBorder", "CycifStack",
  function(x,strict=FALSE,dth,parallel=TRUE,min.pts=3,n.cores=1,ct_name="default",verbose=TRUE,cancer.cts = "Tumor"){
    if(missing(dth)){
      dths <- cyApply(x,function(cy){
        majl <- cy@segment_property$MajorAxisLength
        m1sd <- mean(majl) + sd(majl)
        dth <- m1sd * 2 # 35
        return(dth)
      },simplify=T)
      dth <- median(dths)
    }
    lst <- cyApply(x,function(cy)defineTumorBorder(cy,strict=strict,dth=dth,parallel=parallel,min.pts=min.pts,n.cores=n.cores,ct_name=ct_name,verbose=verbose,cancer.cts=cancer.cts))
    return(lst)
  }
)
