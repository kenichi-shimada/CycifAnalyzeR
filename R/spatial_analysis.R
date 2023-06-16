#_ -------------------------------------------------------

# fun: defineTumorBed Cycif ----
#' Identify Tumor bed based on dbscan clusters
#' @export
setGeneric("defineTumorBed", function(x,...) standardGeneric("defineTumorBed"))

#' @rdname defineTumorBed
#'
#' @param x A Cycif obj
#' @param strict logical. if cell_types should be strict.
#' @param dth radius for dbscan (eps)
#' @param min.pts minimum points for dbscan
#' @param parallel logical. The default is TRUE.
#' @param n.cores number of cores.
#'
#' @importFrom dbscan dbscan
#' @importFrom dbscan frNN
#' @importFrom parallel mclapply
#' @importFrom concaveman concaveman
#' @importFrom sp point.in.polygon
#' @export
setMethod("defineTumorBed", "Cycif",
  function(x,strict=FALSE,dth,parallel=TRUE,min.pts=2,n.cores=1,ct_name){

    ## coordinates on slides
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    seg.prop <- x@segment_property



    ## cell types
    this.cts <- cell_types(x,strict=strict,ct_name=ct_name)$cell_types
    if(any(levels(this.cts)!="NA")){
      this.cts <- factor(this.cts,levels=c(levels(this.cts),"NA"))
    }
    this.cts[is.na(this.cts)] <- "NA"

    ## identify Tumor and non-Tumor cells within ROIs
    is.tumor <- !is.na(this.cts) & this.cts == "Tumor"
    xy.tumor <- xy[is.tumor,]

    is.imm <- !this.cts %in% c("Tumor","others","unknown","NA") # currently not used

    ### find optimal eps (dth) for adjacent cells - from all the cells (not only tumors)
    if(missing(dth)){
      # fr <- dbscan::frNN(xy,eps=100)
      # has.ns <- sapply(fr$dist,length)>0
      # min.ds <- sapply(fr$dist[has.ns],function(x)x[1])
      # # xh <- hist(min.ds,breaks=200,plot=TRUE)
      # dth <- quantile(min.ds,.99)
      dth <- median(seg.prop$MajorAxisLength)
    }

    ## dbscan clustering - only tumor cells
    dbs1 <- dbscan::dbscan(xy.tumor, eps=dth, minPts = min.pts, weights = NULL, borderPoints = TRUE)
    cls <- dbs1$cluster

    ## Find neighbors using frNN
    frt <- dbscan::frNN(xy.tumor,eps=dth)
    fr.ids <- sapply(seq(frt$id),function(i){
      tmp <- frt$id[[i]]
      tmp <- tmp[tmp > i] ## this gaurantees unique pairs
    })

    ## compute the size of the tumor chunks
    ## multinomial distribution? - check 'connectivity' of graph, or 'degree distribution'

    table(table(cls))

    ## Find within tumors
    is.in <- edge.pts <- list()

    nls <- unique(cls)
    nls <- sort(nls[nls!=0])

    if(parallel){
      ## edge.points
      edge.pts <- parallel::mclapply(seq(nls),function(i){
        xyt <- xy.tumor[cls==nls[i],]
        edge.pt <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = 0))
        return(edge.pt)
      },mc.cores=n.cores)
      is.in <- parallel::mclapply(seq(edge.pts),function(i){
        is.in <- sp::point.in.polygon(xy$X,xy$Y,edge.pts[[i]]$V1,edge.pts[[i]]$V2)==1
        return(is.in)
      },mc.cores=n.cores)
    }else{
      for(i in seq(nls)){
        xyt <- xy.tumor[cls==nls[i],]
        edge.pts[[i]] <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = 0))
        is.in[[i]] <- sp::point.in.polygon(xy2$X,xy2$Y,edge.pts[[i]]$V1,edge.pts[[i]]$V2)==1
      }
    }

    is.in.1 <- rowSums(do.call(cbind,is.in))>0

    lst <- list(
      dist.th = dth,
      min.pts = 2,
      xy = xy, ## all cells
      within.tumor.bed = is.in.1, ## all cells
      edge.pts = edge.pts, ## for non-empty clusters
      cell.types = this.cts, ## all cells - tumors can be retrieved here
      clusters = cls, # only tumors
      neighbors = fr.ids # only tumors
    )
    return(lst)
  })

#_ -------------------------------------------------------

# scimap ----

