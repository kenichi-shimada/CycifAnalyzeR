# plan
# - ratioTumorStroma(): change to 'generic function'
# - make defineTumorBorder() available for CycifStack obj
# - make slots for output of defineTumorBorder() in Cycif object
# - make slots for output of ratioTumorStroma() so this function doesn't have to run defineTumorBorder()
# spatialStat obj (temtative)
# for cycif
# - lst (should have a better name)
# for cycifstack
# - ratio.ts

#_ -------------------------------------------------------

# fun: distAdjacentCells Cycif ----

#' Calculate the Distance Between Adjacent Cells in CyCIF Data.
#' @export
setGeneric("distAdjacentCells", function(x,...) standardGeneric("distAdjacentCells"))

#' @rdname distAdjacentCells
#' @export
setMethod("distAdjacentCells", "Cycif",
  function(x){
    majl <- x@segment_property$MajorAxisLength
    m1sd <- mean(majl) + sd(majl)
    dth <- m1sd * 2 # 35
    return(dth)
  })

#' @rdname distAdjacentCells
#' @export
setMethod("distAdjacentCells", "CycifStack",
  function(x){
    dths <- cyApply(x,function(cy){
      majl <- cy@segment_property$MajorAxisLength
      m1sd <- mean(majl) + sd(majl)
      dth <- m1sd * 2 # 35
      return(dth)
    },simplify=T)
  dth <- mean(dths) # 34.86 px = 22.65 um
  return(dth)
  })

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
  function(x,strict=FALSE,dth,parallel=TRUE,min.pts=3,n.cores=1,ct_name="default",verbose=TRUE){
    require(dbscan)
    require(RColorBrewer)
    require(concaveman)
    require(sp)
    if(verbose){
      cat("processing defineTumorBorder() for ",names(x),"...\n")
    }

    ## coordinates on slides
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid

    ## cell types
    this.cts <- cell_types(x,strict=strict,ct_name=ct_name)$cell_types
    if(any(levels(this.cts)!="NA")){
      this.cts <- factor(this.cts,levels=c(levels(this.cts),"NA"))
    }
    this.cts[is.na(this.cts)] <- "NA"

    ## identify Tumor and non-Tumor cells within ROIs
    is.tumor <- !is.na(this.cts) & this.cts == "Tumor"
    xy.tumor <- xy[is.tumor,]

    ### find optimal eps (dth) for adjacent cells - from all the cells (not only tumors)
    if(missing(dth)){
      majl <- x@segment_property$MajorAxisLength
      m1sd <- mean(majl) + sd(majl)*2
      dth <- m1sd * 2 # 35
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

    ## Find within tumors
    is.in <- edge.pts <- list()

    nls <- unique(cls)
    nls <- sort(nls[nls!=0])

    if(parallel){
      ## edge.points
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
      xy = xy, ## all cells
      within.tumor.bed = is.in.1, ## all cells
      edge.pts = edge.pts, ## for non-empty clusters
      cell.types = this.cts, ## all cells - tumors can be retrieved here
      clusters = cls, # only tumors
      neighbors = fr.ids # only tumors
    )
    return(lst)
  }
)

#' @rdname defineTumorBorder
#' @export
setMethod("defineTumorBorder", "CycifStack",
  function(x,strict=FALSE,dth,parallel=TRUE,min.pts=3,n.cores=1,ct_name="default",verbose=TRUE){
    if(missing(dth)){
      dths <- cyApply(x,function(cy){
        majl <- cy@segment_property$MajorAxisLength
        m1sd <- mean(majl) + sd(majl)
        dth <- m1sd * 2 # 35
        return(dth)
      },simplify=T)
      dth <- median(dths) # roughly around 35 ptx = 23 um
    }
    lst <- cyApply(x,function(cy)defineTumorBorder(cy,strict=strict,dth=dth,parallel=parallel,min.pts=min.pts,n.cores=n.cores,ct_name=ct_name,verbose=verbose))
    return(lst)
  }
)

#_ -------------------------------------------------------
# fun: ratioTumorStroma Cycif,CycifStack ----

#' Calculate the Tumor-to-Stroma Ratio in CyCIF data.
#'
#' This function calculates the Tumor-to-Stroma Ratio (TSR) within CyCIF data for either a single CyCIF object or a CyCIFStack object. The TSR is calculated based on the presence of tumor cells within defined tumor bed regions and the specified stroma cell types.
#'
#' @param x A CyCIF or CyCIFStack object.
#' @param stroma.cts Character vector, specifying the stroma cell types to consider for TSR calculation. If not provided, stroma cell types are automatically determined.
#' @param ... Additional arguments to be passed to the \code{\link{defineTumorBorder}} function.
#'
#' @return A numeric vector or data frame containing the Tumor-to-Stroma Ratios (TSR) for each cell type in the dataset.
#'
#' @importFrom concaveman concaveman
#' @importFrom sp point.in.polygon
#' @importFrom dbscan dbscan
#' @importFrom parallel mclapply
#' @importFrom RColorBrewer brewer.pal
#'
#' @seealso
#' \code{\link{defineTumorBorder}}
#'
#' @details
#' The Tumor-to-Stroma Ratio (TSR) is calculated based on the presence of tumor cells within defined tumor bed regions and the specified stroma cell types. The function performs the following steps:
#' 1. If the input is a CyCIF object, it calculates the TSR for that object. If the input is a CyCIFStack object, it calculates the TSR for each CyCIF object in the stack.
#' 2. It determines the tumor bed regions within the CyCIF data using the \code{\link{defineTumorBorder}} function.
#' 3. It counts the number of cells of each cell type within the tumor bed regions.
#' 4. If the `stroma.cts` parameter is not provided, stroma cell types are automatically determined based on unique cell types in the dataset, excluding "Tumor" and "outOfROI" cell types.
#' 5. It calculates the overall TSR (r.all) as the ratio of stroma cells to all cells within the tumor bed.
#' 6. It calculates the TSR for each specified stroma cell type, considering the ratio of that cell type to all cells within the tumor bed.
#' 7. For CyCIFStack objects, the function returns a data frame containing TSR values for each CyCIF object in the stack.
#'
#' @rdname ratioTumorStroma
#' @export
setGeneric("ratioTumorStroma", function(x,...) standardGeneric("ratioTumorStroma"))

#' @rdname ratioTumorStroma
#' @export
setMethod("ratioTumorStroma", "Cycif",
          function(x,stroma.cts,...){
    n <- names(x)
    cat("processing ratioTumorStroma() for ",n,"...\n")
    lst0 <- defineTumorBorder(x,...)

    tab <- table(cell_types(x)$cell_types,lst0$within.tumor.bed)

    if(missing(stroma.cts)){
      uniq.cts <- rownames(tab)
      stroma.cts <- uniq.cts[!uniq.cts %in% c("Tumor","outOfROI")]
    }

    tab1 <- tab[stroma.cts,]
    sum.str <- colSums(tab1)
    r.all <- sum.str["TRUE"]/sum(sum.str)

    r.each <- tab1[,"TRUE"]/rowSums(tab1)
    r.ts <- r.each/r.all

    return(r.ts)
  }
)

#' @rdname ratioTumorStroma
#' @export
setMethod("ratioTumorStroma", "CycifStack",
          function(x,stroma.cts,...){
    ratio <- cyApply(x,ratioTumorStroma,stroma.cts=stroma.cts,...)
    r.ts <- data.frame(do.call(rbind,ratio))
    return(r.ts)
  }
)

#_ -------------------------------------------------------
# fun: dist2tumorBorder Cycif ----

#' Calculate the Distance to Tumor Border in CyCIF data.
#'
#' This function calculates the distance of each cell in CyCIF data to the nearest tumor border. The tumor border is defined based on the spatial arrangement of tumor cells within the dataset.
#'
#' @param x A CyCIF object.
#' @param n.cores Integer, number of CPU cores to use for parallel processing.
#' @param dth Numeric, the distance threshold used for tumor border detection.
#' @param minPts Integer, the minimum number of points required to form a cluster in DBSCAN.
#' @param concavity a relative measure of concavity. 1 results in a relatively detailed shape, Infinity results in a convex hull. You can use values lower than 1, but they can produce pretty crazy shapes (\code{concaveman}).
#' @param ... Additional arguments to be passed.
#'
#' @return A list containing the calculated distances and tumor border polygons.
#'
#' @importFrom concaveman concaveman
#' @importFrom sp point.in.polygon
#' @importFrom dbscan dbscan
#' @importFrom parallel mclapply
#' @importFrom rgeos gDistance
#'
#' @details
#' The `dist2tumorBorder` function calculates the distance of each cell in the CyCIF dataset to the nearest tumor border. It performs the following steps:
#'
#' 1. It extracts the coordinates of cells from the CyCIF data, considering the spatial arrangement of tumor cells.
#' 2. It identifies tumor cells and their clusters using DBSCAN clustering based on the specified distance threshold (`dth`).
#' 3. It computes the tumor border polygons by using the Concaveman algorithm on each tumor cell cluster.
#' 4. It determines if each non-tumor cell is inside or outside of the tumor border polygons using the `point.in.polygon` function.
#' 5. It calculates the distance of each cell to the nearest tumor border, considering whether the cell is inside or outside the tumor border.
#' 6. The function returns a list containing the calculated distances and the tumor border polygons.
#'
#' @seealso
#' \code{\link{concaveman}}, \code{\link{point.in.polygon}}, \code{\link{dbscan}} \code{\link{gDistance}}
#'
#' @rdname dist2tumorBorder
#' @export
setGeneric("dist2tumorBorder", function(x,...) standardGeneric("dist2tumorBorder"))

#' @rdname dist2tumorBorder
#' @export
setMethod("dist2tumorBorder","Cycif",
          function(x, n.cores = 7, dth, th, minPts = 3, concavity = 0.8, plot = FALSE,...) {
    n <- names(x)
    cat("Processing ",n," ... \n")

    this.cts <- cell_types(x)$cell_types # 89174

    ## coordinates: x => xy1, xy.sp1
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    # xy.sp <- sp::SpatialPoints(xy)
    xy.sf <- st_as_sf(xy, coords = c("X_centroid", "Y_centroid"), crs = NA)

    wr <- x@within_rois
    xy.sf1 <- xy.sf[wr,]
    this.cts1 <- this.cts[wr]
    is.tumor <- this.cts1 == "Tumor"

    xy.sf2 <- xy.sf1[is.tumor,]
    xy2 <- st_coordinates(xy.sf2)

    dbs1 <- dbscan::dbscan(xy2, eps=dth, minPts = minPts, weights = NULL, borderPoints = TRUE) # min 2 cells together
    cls3 <- dbs1$cluster
    nls3 <- unique(cls3) # unique cluster labels
    nls3 <- sort(nls3[nls3!=0])
    names(nls3) <- nls3

    idx <- cls3 >0
    cls3[idx] <- nls3[cls3[idx]] # rename clusters

    cat(" Computing borders ... \n")

    borders <- parallel::mclapply(seq(nls3),function(i){
      xyt <- xy.sf2[cls3==nls3[i],]
      conc <- concaveman::concaveman(xyt, concavity = concavity, length_threshold = dth)
      return(conc)
    },mc.cores=n.cores)
    names(borders) <- nls3

    ## rewrite the above using sf_within or sf_intersects. Note borders are a list of POLYGON type of sf objects

    cat(" Find out overlaps of borders ... \n")

    b <- do.call(rbind,borders)
    int <- st_intersects(b,sparse=FALSE)
    diag(int) <- FALSE
    int[lower.tri(int)] <- FALSE
    if(any(int)){
      ol.idx <- which(int,arr.ind=TRUE)[rev(seq(sum(int))),]

      used.borders <- rep(1,length(borders))
      for(i in seq(nrow(ol.idx))){
        borders[[ol.idx[i,1]]] <-
                 st_union(borders[[ol.idx[i,1]]],borders[[ol.idx[i,2]]])
        used.borders[ol.idx[i,2]] <- 0
      }
      borders.1 <- borders[as.logical(used.borders)]
    }else{
      borders.1 <- borders
    }

    b1 <- do.call(rbind,borders.1)
    b2 <- st_boundary(b1)

    cat(" Compute the distance bewteen each point and border of polygons ... \n")

    ## bounding box for each border and point
    bb <- st_make_grid(b1,cellsize=1000)
    if(0){
      plot(bb)
      plot(b1,add=T,col=2)
    }

    is.in.b <- st_intersects(b1,bb,sparse=FALSE)

    ## min bb index for each polygon (tumor region)
    min.bb.idx <- apply(is.in.b,1,function(i)min(which(i)))

    ## list of polygon indices for each bb
    bb.idx <- tapply(seq(min.bb.idx),min.bb.idx,identity)

    ## names of used bounding boxes
    used.bb.idx <- names(bb.idx)

    ## nearest polygon for each point
    near.feature <- st_nearest_feature(xy.sf1,b1)

    ## list of points per polygon
    near.idx <- tapply(seq(near.feature),near.feature,identity)

    ## distance to the nearest polygon
    sgn.d <- rep(NA,length(near.feature))
    for(bbi in used.bb.idx){
      # cat(bbi,"\n")
      border.i <- bb.idx[[bbi]]
      point.i <- unlist(near.idx[border.i])
      if(0){
        ## example - "cy29", bbi=15
        plot(xy.sf1[point.i,][[1]],col=uniq.cols[as.character(d1[point.i])],pch=20)
        plot(b1[border.i,],border=1,add=T,lwd=4)
      }
      d.i1 <- st_distance(b1[border.i,],xy.sf1[point.i,])
      tmp.sgn.d <- round(apply(d.i1,2,min))
      if(any(tmp.sgn.d==0)){
        is0 <- which(tmp.sgn.d==0)
        d.i2 <- st_distance(b2[border.i,],xy.sf1[point.i[is0],])
        tmp.sgn.d[is0] <- -round(apply(d.i2,2,min))
      }
      sgn.d[point.i] <- tmp.sgn.d
    }

    if(0){
      ## tolerance for distance (upper bound)
      tor <- round(dth*3)

      d1 <- sgn.d
      d1[d1 > tor] <- tor
      d1[d1 < -tor] <- -tor

      uniq.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(tor*2+1)
      names(uniq.cols) <- seq(-tor,tor)
      plot(xy.sf1[[1]],col=uniq.cols[as.character(d1)],pch=".")
      plot(b2,border=1,add=T)
      plot(bb,add=T)
    }

    if(plot){
      if(missing(th)){
        th <- round(dth * 3)
      }
      sgn.d1 <- sgn.d
      sgn.d1[sgn.d1 > th] <- th
      sgn.d1[sgn.d1 < -th] <- -th

      cols <- colorRampPalette((brewer.pal(11,"Spectral")))(th * 2 + 1)
      names(cols) <- seq(-th,th)

      plot(xy.sf1[[1]],pch=".",col=cols[as.character(round(sgn.d1))])
      plot(b1,add=T)
    }

    return(list(dist=sgn.d,points=near.idx,borders=b1))
  }
)

#_ -------------------------------------------------------
# fun: defineTumorStroma Cycif ----

#' Define tumor bed regions within each Cycif dataset.
#'
#' This function identifies tumor bed regions within each Cycif dataset and assigns cells to these regions based on their spatial coordinates and cell types.
#'
#' @param x A CyCIF object.
#' @param strict Logical, specifying whether to perform strict cell type calling (see \code{defineCellTypes}).
#' @param dth Numeric, the distance threshold for identifying adjacent cells.
#' @param n.cores Integer, the number of CPU cores to use for parallel processing.
#' @param n.cells.per.seg Integer, the minimum number of cells per segment.
#' @param n.cells.per.tumor.core Integer, the minimum number of cells per tumor region.
#' @param concavity a relative measure of concavity. 1 results in a relatively detailed shape, Infinity results in a convex hull. You can use values lower than 1, but they can produce pretty crazy shapes (\code{concaveman}).
#' @param plot Logical, indicating whether to display plots.
#'
#' @return A list containing information about tumor bed regions, including distance threshold, cell coordinates, regions within tumor bed, edge points, cell types, clusters, and neighbors.
#'
#' @importFrom concaveman concaveman
#' @importFrom sf st_as_sf st_make_grid st_nearest_feature st_distance st_intersects st_union st_boundary
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @importFrom RColorBrewer brewer.pal
#'
#' @details
#' The `dist2tumorBorder` function calculates the distance of each cell in the CyCIF dataset to the nearest tumor border. It performs the following steps:
#'
#' 1. It extracts the coordinates of cells from the CyCIF data, converting them into an sf object.
#' 2. It first identifies all the cells within ROIs and compute the boundary for each tissue segment. Note that a segment containing fewer than `n.cells.per.seg` cells is ignored.
#' 3. For each segment, it identifies tumor cells and their clusters using DBSCAN clustering based on the specified distance threshold (`dth`).
#' 4. It computes the tumor border polygons by using the DBSCAN and concaveman algorithm on each tumor cell cluster.
#' 5. It determines if each non-tumor cell is inside or outside of the tumor border polygons using the `st_distance` function. For the cells inside any polygons are assigned to the tumor bed, and their distance to the nearest border is calculated (expressed as negative value)
#' 6. The function returns a list containing the calculated distances and the tumor border polygons.
#'
#'
#' @seealso
#' \code{\link{concaveman}}, \code{\link{dbscan}}, \code{\link{frNN}}, \code{\link{st_distance}}, \code{\link{st_intersects}}, \code{\link{st_union}}, \code{\link{st_boundary}}
#'
#' @rdname defineTumorStroma
#' @export
setGeneric("defineTumorStroma", function(x,...) standardGeneric("defineTumorStroma"))

#' @rdname defineTumorStroma
#' @export
setMethod("defineTumorStroma","Cycif",
  function(x, n.cores = 7, n.cells.per.seg = 100,  n.cells.per.tumor.core = 10,
           dth, concavity = 0.8, plot = FALSE,...) {
    this.cts <- cell_types(x)$cell_types # 89174

    ## define connected blocks
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    # xy.sp <- sp::SpatialPoints(xy)
    xy.sf <- st_as_sf(xy, coords = c("X_centroid", "Y_centroid"), crs = NA)

    wr <- x@within_rois
    xy.sf1 <- xy.sf[wr,]
    this.cts1 <- this.cts[wr]
    xy1 <- st_coordinates(xy.sf1)

    cat("Identifying tissue segments ",n," ... ")
    dbs2 <- dbscan::dbscan(xy1, eps=dth*5, minPts = 3, weights = NULL, borderPoints = TRUE) # min 2 cells together
    cls2 <- dbs2$cluster
    nls2 <- unique(cls2) # unique cluster labels
    nls2 <- sort(nls2[nls2!=0])
    names(nls2) <- nls2
    n.cls2 <- table(cls2)
    used.cls2 <- names(which(n.cls2 > n.cells.per.seg))
    used.cls2 <- used.cls2[used.cls2 != "0"]
    idx2 <- cls2 %in% used.cls2
    # idx2 <- cls2 > 0

    cls2.1 <- cls2
    cls2.1[!idx2] <- 0
    cls2.1[idx2] <- as.numeric(factor(cls2.1[idx2]))
    nls2.1 <- unique(cls2.1)
    nls2.1 <- nls2.1[nls2.1!=0]

    cat(max(nls2.1)," segments identified\n")
    cat(" Computing borders for each segment ... \n")

    borders2 <- parallel::mclapply(seq(nls2.1),function(i){
      xyt <- xy.sf1[cls2.1==nls2.1[i],]
      conc <- concaveman::concaveman(xyt, concavity = concavity, length_threshold = dth*5)
      return(conc)
    },mc.cores=n.cores) # this gives an error when called within function -then run concaveman in console first
    if(class(borders2)=="try-error"){
      stop("Error in concaveman::concaveman")
    }

    names(borders2) <- nls2.1

    ## rewrite the above using sf_within or sf_intersects. Note borders2 are a list of POLYGON type of sf objects

    cat(" Find out overlaps of the segment borders ... \n")

    b2 <- do.call(rbind,borders2)
    int <- st_intersects(b2,sparse=FALSE)
    diag(int) <- FALSE
    int[lower.tri(int)] <- FALSE
    if(any(int)){
      ol.idx <- which(int,arr.ind=TRUE)[rev(seq(sum(int))),]

      used.borders2 <- rep(1,length(borders2))
      for(i in seq(nrow(ol.idx))){
        borders2[[ol.idx[i,1]]] <-
          st_union(borders2[[ol.idx[i,1]]],borders2[[ol.idx[i,2]]])
        used.borders2[ol.idx[i,2]] <- 0
      }
      borders2.1 <- borders2[as.logical(used.borders2)]
    }else{
      borders2.1 <- borders2
    }

    if(0){
      b2.1 <- do.call(rbind,borders2.1)
      plot(xy.sf1[[1]],col=cls2.1+1,pch=".")
      plot(b2.1,border=1,add=T)
    }

    ## define tumor regions
    cat("Identifying borders for tumor regions for each segment ... \n")

    ## following
    lst.segs <- list()
    for (seg.i in nls2.1){
      cat("Segment",seg.i,":\n")
      # subsetting cells for each segment
      this.seg <- cls2.1==nls2.1[seg.i]

      xy.sf2 <- xy.sf1[this.seg,]
      this.cts2 <- this.cts1[this.seg]
      this.cts2[is.na(this.cts2)] <- "all_other"

      is.tumor <- !is.na(this.cts2) & this.cts2 == "Cancer"

      if(sum(is.tumor)==0){
        cat("No tumor region found for seg",seg.i,"\n")
        lst.segs[[as.character(seg.i)]] <-
          list(this.seg = this.seg,
               dist = rep(NA,nrow(xy.sf2)),
               closest.border.idx = NA,
               borders=NA)
        next
      }

      xy.sf3 <- xy.sf2[is.tumor,]
      xy3 <- st_coordinates(xy.sf3)

      dbs1 <- dbscan::dbscan(xy3, eps=dth, minPts = 3, weights = NULL, borderPoints = TRUE) # min 2 cells together
      cls3 <- dbs1$cluster
      nls3 <- unique(cls3) # unique cluster labels
      nls3 <- sort(nls3[nls3!=0])
      names(nls3) <- nls3

      n.cls3 <- table(cls3)
      tum.reg.cls3 <- names(which(n.cls3 > n.cells.per.tumor.core)) ## tumor region
      tum.reg.cls3 <- tum.reg.cls3[tum.reg.cls3 != "0"] ## tumor bud
      tum.bud.cls3 <- as.character(nls3)[!nls3 %in% tum.reg.cls3]
      idx3 <- cls3 %in% tum.reg.cls3

      ## tumor regions (excluding tumor buds)
      cls3.1 <- cls3
      cls3.1[!idx3] <- 0
      cls3.1[idx3] <- as.numeric(factor(cls3.1[idx3]))
      nls3.1 <- unique(cls3.1)
      nls3.1 <- nls3.1[nls3.1!=0]

      cat(" Computing borders ... \n")

      borders3 <- parallel::mclapply(seq(nls3.1),function(i){
        xyt <- xy.sf3[cls3.1==nls3.1[i],]
        conc <- concaveman::concaveman(xyt, concavity = concavity, length_threshold = dth)
        return(conc)
      },mc.cores=n.cores)
      names(borders3) <- nls3.1

      if(length(borders3)==0){
        cat("No tumor region found for seg",seg.i,"\n")
        lst.segs[[as.character(seg.i)]] <-
          list(this.seg = this.seg,
               dist = rep(NA,nrow(xy.sf2)),
               closest.border.idx = NA,
               borders=NA)
        next
      }

      ## rewrite the above using sf_within or sf_intersects. Note borders are a list of POLYGON type of sf objects

      cat(" Find out overlaps of borders ... \n")

      b3 <- do.call(rbind,borders3)
      int3 <- st_intersects(b3,sparse=FALSE)
      diag(int3) <- FALSE
      int3[lower.tri(int3)] <- FALSE
      if(any(int3)){
        ol.idx <- which(int3,arr.ind=TRUE)[rev(seq(sum(int3))),]

        used.borders3 <- rep(1,length(borders3))
        for(i in seq(nrow(ol.idx))){
          borders3[[ol.idx[i,1]]] <-
            st_union(borders3[[ol.idx[i,1]]],borders3[[ol.idx[i,2]]])
          used.borders3[ol.idx[i,2]] <- 0
        }
        borders3.1 <- borders3[as.logical(used.borders3)]
      }else{
        borders3.1 <- borders3
      }

      b3.1 <- do.call(rbind,borders3.1)
      b3.2 <- st_boundary(b3.1)

      cat(" Compute the distance bewteen each point and border of polygons ... \n")

      ## bounding box for each border and point
      bb3 <- st_make_grid(b3.1,cellsize=1000)
      if(0){
        plot(bb3)
        plot(b3.1,add=T,col=2)
      }

      is.in.b3 <- st_intersects(b3.1,bb3,sparse=FALSE)

      min.bb.idx <- apply(is.in.b3,1,function(i)min(which(i))) ## min bb index for each polygon (tumor region)
      bb.idx <- tapply(seq(min.bb.idx),min.bb.idx,identity) ## list of polygon indices for each bb
      used.bb.idx <- names(bb.idx) ## names of used bounding boxes
      near.feature <- st_nearest_feature(xy.sf2,b3.1) ## nearest polygon for each point
      near.idx <- tapply(seq(near.feature),near.feature,identity) ## list of points per polygon

      sgn.d <- rep(NA,length(near.feature)) ## distance to the nearest polygon
      for(bbi in used.bb.idx){
        # cat(bbi,"\n")
        border.i <- bb.idx[[bbi]]
        point.i <- unlist(near.idx[border.i])
        d.i1 <- st_distance(b3.1[border.i,],xy.sf2[point.i,])
        tmp.sgn.d <- round(apply(d.i1,2,min))
        if(any(tmp.sgn.d==0)){
          is0 <- which(tmp.sgn.d==0)
          d.i2 <- st_distance(b3.2[border.i,],xy.sf2[point.i[is0],])
          tmp.sgn.d[is0] <- -round(apply(d.i2,2,min))
        }
        sgn.d[point.i] <- tmp.sgn.d
      }

      lst.segs[[as.character(seg.i)]] <-
        list(this.seg = this.seg,
             dist = sgn.d,
             closest.border.idx = near.feature,
             borders=b3.1)

    }

    closest.border.idxs <- dists <- rep(NA,nrow(xy.sf1))
    borders.all <- list()
    for(seg.i in names(lst.segs)){
      this.borders <- lst.segs[[seg.i]]$borders
      if(!is(this.borders,"sf") && is.na(this.borders)){
        next
      }
      this.seg <- lst.segs[[seg.i]]$this.seg
      dists[this.seg] <- lst.segs[[seg.i]]$dist
      closest.border.idxs[this.seg] <- paste0(seg.i,".",lst.segs[[seg.i]]$closest.border.idx)
      borders.all[[seg.i]] <- this.borders
    }

    lst <- list(dist=data.frame(seg=cls2.1,dist=dists,points=closest.border.idxs),borders=borders.all)

    if(0){
      ## tolerance for distance (upper bound)
      tor <- round(dth*3)

      d1 <- lst$dist$dist
      d1[d1 > tor] <- tor
      d1[d1 < -tor] <- -tor
      d1[is.na(d1)] <- tor+1

      uniq.cols <- c(colorRampPalette(brewer.pal(11,"Spectral"))(tor*2+1),"grey")
      names(uniq.cols) <- seq(-tor,tor+1)
      plot(xy.sf1[[1]],col=uniq.cols[as.character(d1)],pch=".")
      plot(lst$borders[[1]],border=1,add=T)
      # plot(bb3,add=T)
    }

    return(lst)
  }
)




# Print the result
#_ -------------------------------

# fun: computeArea (for density) Cycif ----

#' Compute the Area of Tumor Regions in CyCIF Data.
#'
#' This function calculates the total area of tumor regions within a CyCIF dataset based on the specified distance threshold for tumor border detection.
#'
#' @param x A CyCIF object.
#' @param dth Numeric, the distance threshold used for tumor border detection.
#' @param unit Character, the unit of the computed area in the output. Default is "mm2" (square millimeters).
#' @param plot Logical, whether to plot the tumor regions. Default is TRUE.
#' @param strict Logical, whether to use strict cell type filtering. Default is FALSE.
#' @param ct_name Character, the name of the cell type for tumor identification. Default is "default".
#' @param fn Character, the filename for saving the plot. Ignored if plot is FALSE.
#' @param minPts Integer, the minimum number of points required to form a cluster in DBSCAN.
#' @param eps Numeric, the maximum distance between two samples for one to be considered as in the neighborhood of the other in DBSCAN.
#'
#' @return A numeric value representing the computed area of tumor regions in the specified unit.
#'
#' @details
#' This function calculates the total area of tumor regions within a CyCIF dataset based on the specified distance threshold for tumor border detection. It uses a combination of DBSCAN clustering to identify tumor cell clusters and concave hull computation to estimate the tumor regions' boundaries. The area is then computed based on the identified tumor regions.
#' The process involves the following steps:
#' 1. Cell type filtering: If `strict` is set to TRUE, only cells with the specified `ct_name` (cell type name) will be considered as tumor cells; otherwise, all non-NA cell types will be considered as tumor cells.
#' 2. DBSCAN clustering: DBSCAN (Density-Based Spatial Clustering of Applications with Noise) is used to cluster the identified tumor cells into groups based on their spatial proximity. The parameters `minPts` (minimum number of points required to form a cluster) and `eps` (maximum distance between two points to be considered part of the same cluster) can be customized.
#' 3. Cluster merging: Overlapping clusters are merged to create distinct tumor regions.
#' 4. Concave hull computation: For each tumor region, a concave hull is computed using the `concaveman` package to approximate its boundary.
#' 5. Area calculation: The area of each tumor region is computed using the `sp::Polygon` function, and the areas of all tumor regions are summed to obtain the total area.
#' The computed area is returned as a numeric value, and the unit of measurement can be specified using the `unit` parameter (e.g., "mm2" for square millimeters).
#' If `plot` is set to TRUE, a plot displaying the tumor regions will be generated and saved to the specified `fn` (filename).
#'
#' @importFrom concaveman concaveman
#' @importFrom sp point.in.polygon
#' @importFrom dbscan dbscan
#' @importFrom parallel mclapply detectCores
#'
#' @seealso \code{\link{defineTumorBorder}} for defining tumor regions, \code{\link{concaveman::concaveman}} for concave hull computation, \code{\link{dbscan::dbscan}} for DBSCAN clustering.
#'
#' @rdname computeArea
#' @export
setGeneric("computeArea", function(x,...) standardGeneric("computeArea"))

#' @rdname computeArea
#' @export
setMethod("computeArea", "Cycif",
          function(x,dth,unit=c("mm2"),plot=TRUE,strict=FALSE,
                         ct_name="default",fn){

    ## coordinates
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    xy.sp <- sp::SpatialPoints(xy)

    ##
    ## cell types
    this.cts <- cell_types(x,strict=strict,ct_name=ct_name)$cell_types
    levs <- levels(this.cts)
    levs <- levs[!levs %in% c("NA","outOfROI")]
    this.cts <- factor(this.cts,levels=levs)

    ## find cells within rois
    wr <- x@within_rois
    xy1 <- xy[wr,]
    xy.sp1 <- xy.sp[wr,]
    this.cts1 <- this.cts[wr]

    ## define clusters of tumor chunk using dbscan
    dbs1 <- dbscan::dbscan(xy1, eps=dth*5, minPts = 3, weights = NULL, borderPoints = TRUE) # min 2 cells together
    cls3 <- dbs1$cluster

    ## Merge overlapped clustersFind tumor borders
    nls3 <- unique(cls3) # unique cluster label
    nls3 <- sort(nls3[nls3!=0])

    ncores <- parallel::detectCores(logical = TRUE)-1
    mcc <- min(ncores,length(nls3))

    ##
    i=1
    xyt <- xy1[cls3==nls3[i],]
    conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
    names(conc) <- c("X","Y")
    rm(i,xyt,conc)

    concs3 <- parallel::mclapply(seq(nls3),function(i){
    # for(i in seq(nls3)){
      xyt <- xy1[cls3==nls3[i],]
      conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
      names(conc) <- c("X","Y")
      return(conc)
    # }
    },mc.cores=mcc)

    if(length(concs3)>1){
      # Identify overlapping clusters and merge their points
      is_point_inside_polygon <- function(xyt1, xyt2){
        if(!is.data.frame(xyt1) | !is.data.frame(xyt2)){
          stop("xyt1 and xyt2 must be data frames")
        }
        is.inside <- sp::point.in.polygon(xyt2$X, xyt2$Y, xyt1$X, xyt1$Y)
        return(is.inside)
      }
      overlaps <- sapply(concs3,function(cluster_i){
        sapply(concs3,function(cluster_j){
          any(is_point_inside_polygon(cluster_j, cluster_i))
        })
      })
      diag(overlaps) <- 0

      # update clusters
      ol <- as.data.frame(which(overlaps==1,arr.ind=T))
      names(ol) <- c("child","parent")
      ol1 <- ol

      nls3.updated <- nls3
      names(nls3.updated) <- nls3

      while(nrow(ol1)){
        pa <- ol1$pa[1]
        ch <- ol1$ch[1]
        nls3.updated[ch] <- pa
        ol1$parent[ol1$parent==ch] <- pa
        ol1$child[ol1$child==ch] <- pa
        ol1 <- ol1[-1,]
      }

      idx <- cls3 >0
      cls3.updated <- cls3
      cls3.updated[idx] <- nls3.updated[cls3[idx]]
      nls3.updated <- unique(nls3.updated)

      concs3.updated <- parallel::mclapply(seq(nls3.updated),function(i){
        xyt <- xy1[cls3.updated==nls3.updated[i],]
        conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
        names(conc) <- c("X","Y")
        return(conc)
      },mc.cores=7)
      names(concs3.updated) <- nls3.updated
    }else{
      cls3.updated <- cls3
      nls3.updated <- nls3
      concs3.updated <- parallel::mclapply(seq(nls3.updated),function(i){
        xyt <- xy1[cls3.updated==nls3.updated[i],]
        conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
        names(conc) <- c("X","Y")
        return(conc)
      },mc.cores=7)
      names(concs3.updated) <- nls3.updated
    }

    if(plot){
      png(fn,width=700,height=700)

      par(mar=c(5,5,3,3))
      sp::plot(xy.sp1,col=cls3.updated+1,pch=c(2,20)[(cls3.updated!=0)+1],
           # asp=1,xlim=range(l1$x),ylim=range(l1$y),
           cex=.5,
           main="tumors, clustered (dbscan)")
      sapply(concs3.updated,function(conc)lines(conc))
      dev.off()
    }

    ## compute area
    areas <- sapply(concs3.updated,function(coords){
      polygon <- sp::Polygon(coords)
      return(polygon@area)
    })

    sum.area <- sum(areas) * (0.65 * (10^-3))^2 # unit: mm^2
    return(sum.area)
  }
)

#_ -------------------------------

# fun: computeCN ----

#' Compute Recurrent Cell Neighbors (RCN) for Cycif or CycifStack Objects
#'
#' This function computes the Recurrent Cell Neighbors (RCN) for Cycif or CycifStack objects.
#' RCN measures the relative frequencies of neighboring cell types around each cell within specified
#' radius 'r'. The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
#'
#' @param x A Cycif or CycifStack object.
#' @param r The radius within which neighboring cells are considered (in 'unit').
#' @param unit The unit of measurement for the radius ('pixel' or 'um'). If 'um' is specified,
#' the radius 'r' will be converted to pixels based on the assumed resolution (0.65 um per pixel).
#' @param cts.in.center A character vector specifying the cell types around which RCN is computed. If not is specified, all available cts are used.
#' @param cts.in.rcn A character vector specifying the cell types to consider when computing RCN values. If not is specified, all avaialble cts are used.
#' @param n.sampling The number of cells to randomly sample for RCN analysis.
#' @param seed The random seed for reproducibility.
#'
#' @return A CellNeighborhood object containing the following components:
#' - 'within.rois': A logical vector indicating whether each cell is within a region of interest (ROI). The length is the same as the number of cells in the dataset.
#' - 'cts.in.rcn': A character vector specifying the cell types considered when computing cell neighbors, cell type frequency, and expressions.
#' - 'n.cells.selected': The number of cells selected for RCN analysis, which is the smaller of `n.sampling` and the number of cells within ROIs.
#' - 'frnn': A list with recurrent Neighborhood information, including 'dist', 'id', 'eps', and 'sort'. The length of 'dist' and 'id' is the same as 'n.cells.selected'.
#' - 'cn_exp': A data frame containing expression data for selected cells.
#' - 'is.selected': A logical vector indicating whether each cell is selected for RCN analysis. The sum of the vector is the same as 'n.cells.selected'.
#' - 'rcn.count': A data frame containing the counts of neighboring cell types.
#' - 'rcn.freq': A data frame containing the relative frequencies of neighboring cell types.
#'
#' @details The RCN analysis is performed as follows:
#' 1. The function first identifies the cells that are within the ROIs.
#' 2. It then computes the Recurrent Neighborhood (frNN) for the selected cells using the 'dbscan::frNN' function.
#' 3. It then computes the RCN values for each cell type based on the relative frequencies of neighboring cell types.
#' 4. The function returns a list containing the RCN values for each cell type.
#' The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
#' If the input is a CycifStack object, the RCN analysis is performed on each Cycif object in the stack.
#'
#' @seealso \code{\link{cyApply}}, \code{\link{dbscan::frNN}}
#'
#' @importFrom dbscan frNN
#' @importFrom data.table := rbindlist setDT melt
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr na_if
#'
#' @rdname computeCN
#' @export
setGeneric("computeCN", function(x,...) standardGeneric("computeCN"))

#' @rdname computeCN
#' @export
setMethod("computeCN", "Cycif",
  function(x,r_um=20,k=20,
           type=c("knn","frnn"),
           used.cts,
           n.sampling=1000,
           seed=123){ # only for Cycif, and CycifStack
    ## find cells within rois - roi is a circle with a fixed radius (r_um)
    if(missing(type)){
      type <- "frnn"
    }else if(!type %in% c("knn","frnn")){
      stop("type must be either 'knn' or 'frnn'")
    }

    if(type=="frnn"){
      cat("compute frnn\n")
      r <- r_um/0.65 # unit converted to pixel
    }else{
      cat("compute knn\n")
    }

    smpl <- names(x)
    cts <- cell_types(x)

    ## coordinates: x => xy
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    # xy.sf <- st_as_sf(xy, coords = c("X_centroid", "Y_centroid"), crs = NA)

    wr <- within_rois(x)
    ## ---- cell types to consider (used.cts) - "OutOfROI" is removed for this subset of cells ----
    if(missing(used.cts)){
      lev.cts <- levels(cts$cell_types)
      lev.cts <- lev.cts[lev.cts != "outOfROI"]
      used.cts <- lev.cts
    }

    is.used <- cts$cell_types %in% used.cts

    # ---- expression matrix per cell type per CellNeighborhood ----
    xy1 <- xy[is.used,]
    cts1 <- cts[is.used,]
    df1 <- exprs(x,type="log")[is.used,]

    ## convert xy1, cts1, df1 to data.table
    data.table::setDT(xy1)
    data.table::setDT(cts1)
    data.table::setDT(df1)

    ## ---- compute frNN for all cells within ROIs and convert them to frNN object.  ----
    if(type=="frnn"){
      nn <- dbscan::frNN(xy1,eps=r,bucketSize=10)
      nn.ids <- nn$id <- lapply(seq_along(nn$id), function(i) c(i, nn$id[[i]]))
    }else if(type=="knn"){
      nn <- dbscan::kNN(xy1,k=k,bucketSize=10)
      nnids <- cbind(seq(nrow(nn$id)),nn$id)
      colnames(nnids) <- 0:k
      nn.ids <- nn$id <- apply(nnids,1,function(x)x,simplify=FALSE)
      nn$dist <- apply(nn$dist,1,function(x)c(0,x),simplify=FALSE)
    }

    # Add a row number column to xy1 and cts1 for joining
    xy1[, rn := .I]
    cts1[, rn := .I]
    df1[, rn := .I]

    ## find out which cell types express which cell state markers
    csts <- x@cell_types$default@cell_state_def
    if(colnames(csts)[1] != "cell_types"){
      # stop("The first column of the cell state marker definition table must be 'cell_types'")
      csts <- csts %>% rownames_to_column("cell_types")
    }
    csts <- csts %>% filter(cell_types %in% used.cts)
    protein_columns <- names(csts)
    protein_columns <- protein_columns[protein_columns != "cell_types"]

    # Convert csts to data.table and melt it to long format
    csts_dt <- data.table::melt(data.table::setDT(csts), id.vars = "cell_types", variable.name = "ab", value.name = "expression")

    # Filter for cell types that do not express each antibody (i.e., expression is NA)
    csts_dt <- csts_dt[is.na(expression)]

    # Join xy1, cts1, and df1 by row number
    combined_df <- xy1[cts1, on = "rn"][df1, on = "rn"]

    # Loop through each antibody and set expression to NA for cell types that do not express it
    if(0){ ## decided not to turn them as NAs as I want to compute expression per CN, irrespective of cell type composition.
      for (ab1 in as.character(unique(csts_dt$ab))) {
        non_expressing_cts <- csts_dt[ab == ab1, cell_types]
        combined_df[.(non_expressing_cts), (ab1) := NA, on = .(cell_types)]
      }
    }
    # tapply(combined_df$PD1,combined_df$cell_types,mean,na.rm=T) # sanity check

    # Create a neighborhood data table from nn
    nmat <- data.table::data.table(
      cell_id = rep(seq_along(nn$id), lengths(nn$id)),
      neighbor_id = unlist(nn$id, use.names = FALSE)
    )

    # Join with combined_df to get cell type and expression data for each neighbor
    neighborhood <- nmat[
      combined_df, on = .(neighbor_id = rn), nomatch = 0
    ]

    # Calculate the average expression per cell type in each neighborhood
    # Assuming protein columns in combined_df are named "protein1", "protein2", etc.

    exp_per_ct_cn <- neighborhood[, lapply(.SD, mean, na.rm = TRUE), by = .(cell_id, cell_types), .SDcols = protein_columns]
    exp_per_cn <- neighborhood[, lapply(.SD, mean, na.rm = TRUE), by = .(cell_id), .SDcols = protein_columns]

    if(type=="knn"){
      nn$eps <- -1
    }else if (type=="frnn"){
      nn$k <- -1
    }

    ##  ---- convert nn to nn object (not necessary?)  ----
    nn1 <- new("NN",
               type = type,
               dist = nn$dist,
               id = nn$id,
               k = nn$k,
               eps = nn$eps,
               sort = nn$sort)

    ##  ---- within positive ROIs + has neighbors ----
    n.nn <- lengths(nn1@id) # number of neighbors - the same as sum(wr)

    selected.ids1 <- which(n.nn>0 &
                           sapply(nn.ids,function(ids){
                             any(cts1$cell_types[ids] %in% used.cts)
                           })) # cell indices are after ROI filter

    n.this <- length(selected.ids1)
    n.cts <- min(n.sampling,n.this)

    ## sampling
    set.seed(seed)

    # Randomly sample n.cts cells from selected.ids1 (idx after ROI)
    selected.ids2 <- sample(selected.ids1,n.cts)

    ## selected ids among available focused celltypes (eg tumor cells)
    is.selected <- seq(nn.ids) %in% selected.ids2

    ## cell type count and frequency in each RCN
    rcn.freq <- t(sapply(nn1@id,function(ids){
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      ntab <- tab/sum(tab)
      return(ntab)
    }))

    rcn.count <- t(sapply(nn1@id,function(ids){
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      return(tab)
    }))

    if(type=="knn"){
      r <- sapply(nn$dist,function(x)x[k+1])
    }
    rcn.dens <- round(rcn.count/(pi*r^2),2)

    cn <- new("CellNeighborhood",
              within.rois=wr, # logical, within ROIs
              used.cts=used.cts, # character, cell types to consider
              n.cells.selected=n.cts, # integer, number of cells selected
              is.selected = is.selected,
              smpls = smpl,
              nn=nn1,
              n.neighbors = n.nn,
              exp.per.ct.cn=exp_per_ct_cn,
              exp.per.cn=exp_per_cn,
              rcn.count=rcn.count,
              rcn.dens=rcn.dens,
              rcn.freq=rcn.freq)
    return(cn)
})

#' @rdname computeCN
#' @export
setMethod("computeCN", "CycifStack",
  function(x,r_um = 20,
           used.cts,
           n.sampling,
           seed=123){ # only for Cycif, and CycifStack

    cat("Get neighbors ...\n")
    nn1 <- cyApply(x,function(cy){
      cat(names(cy),"\n")
      computeCN(x=cy,r=r,unit=c("pixel","um"),
        used.cts=used.cts,
        n.sampling=n.sampling,
        seed=seed)
    })

    cat("Restructure data ...\n")

    ## within.rois
    within.rois <- unlist(lapply(nn1,function(fr)fr@within.rois)) # same as nCells() for each sample

    ## n.cells.selected
    n.cells.selected <- sapply(nn1,function(fr)fr@n.cells.selected)

    # is.selected
    is.selected <- unlist(sapply(nn1, function(fr)fr@is.selected))

    if(sum(is.selected) != sum(n.cells.selected)){
      stop("is.selected and n.cells.selected are not consistent")
    }

    ## used.cts
    used.cts <- nn1[[1]]@used.cts

    ## smpls
    n.smpls <- sapply(nn1,function(fr)sum(fr@within.rois))
    smpls <- rep(names(n.smpls),n.smpls)

    ## nn
    lst.nn <- lapply(nn1,function(fr)fr@nn)

    ### nn@dists
    nn.dists <- do.call(c,lapply(lst.nn,function(fr)fr@dist))

    ### nn@id - combine indices so they can specify selected cells in the entire dataset
    n.nns <- sapply(lst.nn,function(nn)length(nn@id)) # 1325874, all cells, excluding outOfROIs
    n.nns.pre <- c(0,cumsum(n.nns)[-length(n.nns)])
    names(n.nns.pre) <- names(x)

    ## nn.ids, nn.ids1, nn.tum.ids - list of neighboring cells ids for tumors used for the nn analysis
    nn.ids <- lapply(names(x),function(nm){
      x <- lst.nn[[nm]]
      n.prior <- n.nns.pre[nm]
      this.ids <- x@id
      new.ids <- lapply(this.ids,function(id){
        new.id <- id + n.prior
        return(new.id)
      })
      return(new.ids)
    })
    nn.ids1 <- do.call(c,nn.ids) ## 1325874, now all data are combined - and the indices are after excluding outOfROIs

    ### nn@eps
    ### nn@sort
    eps <- unique(sapply(lst.nn,function(nn)nn@eps))
    sort <- unique(sapply(lst.nn,function(nn)nn@sort))

    ### assemble nn
    nn1 <- new("NN",
               type = type,
               dist = nn$dist,
               id = nn$id,
               k = nn$k,
               eps = nn$eps,
               sort = nn$sort)

    # n.neighbors
    n.neighbors <- lengths(nn@id)

    # exp.per.ct.cn
    exp.per.ct.cn <- data.table::rbindlist(lapply(nn1,function(nn)nn@exp.per.ct.cn))

    # exp.per.cn
    exp.per.cn <- data.table::rbindlist(lapply(nn1,function(nn)nn@exp.per.cn))

    ## rcn.count
    rcn.count <- as.matrix(data.table::rbindlist(lapply(nn1,function(fr)as.data.frame(fr@rcn.count))))

    ## rcn.freq
    rcn.freq <- as.matrix(data.table::rbindlist(lapply(nn1,function(fr)as.data.frame(fr@rcn.freq))))

    ## is.selected
    mclustda <- list()
    mclustda$sele <- mclustda$all <- list()

    cn <- new("CellNeighborhood",
              within.rois=within.rois,
              n.cells.selected=n.cells.selected,
              is.selected=is.selected,
              used.cts=used.cts,
              smpls=smpls,
              nn=nn,
              exp.per.ct.cn=exp.per.ct.cn,
              exp.per.cn=exp.per.cn,
              rcn.count=rcn.count,
              rcn.freq=rcn.freq,
              mclustda = mclustda)

    return(cn)
  }
)

#_ -------------------------------

# fun: setDist ----

#' set distance to tumorBorder for NN objects
#'
#' @param x A NN object.
#' @param value A numeric vector specifying the distance to tumor border for each cell.
#'
#' @export
setGeneric("setDist", function(x,...) standardGeneric("setDist"))

#' @rdname setDist
#' @export
setMethod("setDist", "CellNeighborhood",
  function(x,value){
    x@dist2tumorBorder <- value # only for Cycif, and CycifStack
    return(x)
  }
)

#_ -------------------------------
# fun: tcnClust ----

#' Cluster and Sort Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack Objects
#'
#' This function clusters and sorts the Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack objects.
#' It clusters cells based on their RCN profiles, sorts clusters based on the specified cell type,
#' and optionally extrapolates the clustering to the entire dataset.
#'
#' @param nn An object containing RCN information, typically obtained from 'computeCN'.
#' @param g The number of clusters to create.
#' @param seed The random seed for reproducibility.
#' @param sort.by The cell type to sort clusters by (e.g., "CD8T").
#' @param sort.type The type of data to use for sorting: "freq" (relative frequencies) or "count" (counts).
#' @param sort.smpls The subset of samples to use for sorting: "all" (entire dataset) or "selected" (selected cells).
#' @param data.type The type of data to use for clustering: "ct_exp" (cell type and expression data) or "ct" (cell type data only).
#' @param extrapolate Whether to extrapolate clusters to the entire dataset (TRUE) or not (FALSE).
#' @param mc.cores The number of CPU cores to use for parallel processing.
#'
#' @return An updated 'nn' object with clustering and sorting information.
#'
#' @details
#' The `tcnClust` function uses the provided `nn` object to perform clustering and classification of cells based on their neighborhood relationships. It allows you to specify the number of clusters (`g`), the cell type to sort clusters by (`sort.by`), and other clustering parameters.
#' The clustering process results in the classification of cells into distinct clusters, and the function provides information about these clusters, including the mean frequencies, counts, and more.
#' By specifying different options for `sort.by`, `sort.type`, and `sort.smpls`, you can customize the sorting behavior of clusters based on cell types and data types.
#' Additionally, you can choose to extrapolate clusters to the entire dataset using the `extrapolate` argument, which can be helpful for analyzing the overall dataset.
#'
#' @seealso \code{\link{computeCN}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom dbscan  frNN kNN
#' @importFrom data.table as.data.table rbindlist
#' @importFrom parallel mclapply
#' @export
setGeneric("tcnClust", function(nn,...) standardGeneric("tcnClust"))

#' @rdname tcnClust
#' @export
setMethod("tcnClust","data.frame",
  function(nn,
           g=50,
           seed=123,
           sort.by="CD8T",
           sort.type=c("freq","count"),
           sort.smpls=c("all","selected"),
           data.type=c("ct_exp","ct"),
           extrapolate=FALSE,
           mc.cores=1){
  mclustda <- nn@mclustda

  exps <- as.matrix(nn@exp)[,-1]
  exps.imp <- imputeData(exps)

  this.cts <- cts.in.rcn

  mat.count.all <- nn$rcn.count[,this.cts]
  mat.freq.all <- t(apply(nn$rcn.freq[,this.cts],1,function(x)x/sum(x)))

  is.selected <- mclustda$sele$is.used

  if(missing(data.type)){
    data.type="ct_exp"
  }

  if(data.type=="ct"){
    mat.count.sele <- mat.count.all[is.selected,]
    mat.freq.sele <- mat.freq.all[is.selected,]
  }else if(data.type=="ct_exp"){
    mat.count.sele <- cbind(mat.count.all,exps.imp)[is.selected,]
    mat.freq.sele <- cbind(mat.freq.all,exps.imp)[is.selected,]
  }

  ## clustering & classification
  set.seed(seed)

  cat("Clustering with Mclust ...\n")
  mem.ori <- mclust::Mclust(data=mat.freq.sele,G=g,modelNames="EII")$classification

  cat("Training MclustDA ...\n")
  mc1 <- mclust::MclustDA(data=mat.freq.sele,class=mem.ori,
                          G=g,modelNames="EII",modelType = "EDDA") # 100, EII
  mem.sele <- factor(predict(mc1)$classification)
  g1 <- nlevels(mem.sele)

  ## mean counts & frequencies

  mean.count.sele <- sapply(seq(g1),function(i)colMeans(mat.count.sele[mem.sele==i,]))
  mean.freq.sele <- sapply(seq(g1),function(i)colMeans(mat.freq.sele[mem.sele==i,]))
  colnames(mean.count.sele) <- colnames(mean.freq.sele) <- seq(g1)

  if(extrapolate){
    ## extrapolate clusters
    is.available <- !apply(mat.freq.all,1,function(x)any(is.na(x))) # 8587
    mat.count.all1 <- mat.count.all[is.available,]
    mat.freq.all1 <- mat.freq.all[is.available,]

    cat("Applying MclustDA to the entire data ...\n")
    mc.idx <- sort(rep(seq(mc.cores),length=nrow(mat.freq.all1)))

    mem.all <- parallel::mclapply(seq(mc.cores),function(i){
      predict(mc1,newdata=mat.freq.all1[mc.idx==i,nn$cts.in.rcn])$classification
    },mc.cores=mc.cores)
    mem.all <- do.call(c,mem.all)

    ## update labels
    mem.all <- factor(mem.all)
    g1 <- nlevels(mem.all)
    levels(mem.all) <- seq(g1)

    ## mean counts & frequencies
    mean.count.all <- sapply(seq(g1),function(i)colMeans(mat.count.all1[mem.all==i,]))
    mean.freq.all <- sapply(seq(g1),function(i)colMeans(mat.freq.all1[mem.all==i,]))

    colnames(mean.count.all) <- colnames(mean.freq.all) <- seq(g1)
  }

  ## Sort clusters based on frequency of a cell type (CD8T by default)
  if(missing(sort.type)){
    sort.type <- "freq"
  }

  if(missing(sort.smpls)){
    if(extrapolate){
      sort.smpls <- "all"
    }else{
      sort.smpls <- "selected"
    }
  }

  if(sort.type == "freq" & sort.smpls == "all"){
    mean.freq <- mean.freq.all
  }else if(sort.type == "freq" & sort.smpls == "selected"){
    mean.freq <- mean.freq.sele
  }else if(sort.type == "count" & sort.smpls == "all"){
    mean.freq <- mean.count.all
  }else if(sort.type == "count" & sort.smpls == "selected"){
    mean.freq <- mean.count.sele
  }

  if(sort.by=="CD8T"){
    o <- order(mean.freq["CD8T",],decreasing=T)
    mean.freq.sele <-  mean.freq.sele[,o]
    mean.count.sele <-  mean.count.sele[,o]
    mem.sele <- as.numeric(factor(as.character(mem.sele),levels=as.character(seq(ncol(mean.freq))[o])))
    colnames(mean.count.sele) <- colnames(mean.freq.sele) <- paste0("Rcn",seq(g1))

    mclustda$sele <- list(is.used = mclustda$sele$is.used,
                          mem = mem.sele,
                          mean.freq = mean.freq.sele,
                          mean.count = mean.count.sele)

    if(extrapolate){
      mean.freq.all <-  mean.freq.all[,o]
      mean.count.all <-  mean.count.all[,o]
      mem.all <- as.numeric(factor(as.character(mem.all),levels=as.character(seq(ncol(mean.freq))[o]))) ## mem redefined
      colnames(mean.count.all) <- colnames(mean.freq.all) <- paste0("RCN",seq(g1))

      mclustda$all <- list(is.used = is.available,
                           mem = mem.all,
                           mean.freq = mean.freq.all,
                           mean.count = mean.count.all)
    }
  }else{
    stop("Define 'sort.by' first")
  }

  ##
  mclustda$g <- g1

  nn$mclustda <- mclustda

  return(nn)
})


#_ -------------------------------

# fun: rcnClust ----

#' Cluster and Sort Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack Objects
#'
#' This function clusters and sorts the Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack objects.
#' It clusters cells based on their RCN profiles, sorts clusters based on the specified cell type,
#' and optionally extrapolates the clustering to the entire dataset.
#'
#' @param cn An object containing CN information, typically obtained from 'computeCN'.
#' @param g The number of clusters to create.
#' @param seed The random seed for reproducibility.
#' @param sort.by The cell type to sort clusters by (e.g., "CD8T").
#' @param sort.type The type of data to use for sorting: "freq" (relative frequencies) or "count" (counts).
#' @param sort.smpls The subset of samples to use for sorting: "all" (entire dataset) or "selected" (selected cells).
#' @param data.type The type of data to use for clustering: "ct_exp" (cell type and expression data) or "ct" (cell type data only).
#' @param extrapolate Whether to extrapolate clusters to the entire dataset (TRUE) or not (FALSE).
#' @param mc.cores The number of CPU cores to use for parallel processing.
#'
#' @return An updated 'nn' object with clustering and sorting information.
#'
#' @seealso \code{\link{computeCN}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom data.table as.data.table rbindlist
#' @importFrom parallel mclapply
#' @export
setGeneric("rcnClust", function(cn,...) standardGeneric("rcnClust"))

#' @rdname rcnClust
#' @export
setMethod("rcnClust","CellNeighborhood",
          function(cn,
                   g=50,
                   seed=123,
                   sort.by="dist",
                   sort.type=c("freq","count"),
                   sort.smpls=c("all","selected"),
                   data.type=c("ct_exp","ct"),
                   extrapolate=FALSE,
                   mc.cores=1){

            mclustda <- cn@mclustda
            is.selected <- cn@is.selected # n1 = 1325874
            smpls <- cn@smpls # n1

            ## exps
            exps <- as.matrix(cn@exp.per.cn)[,-1] # nrow = n1
            exps.imp <- imputeData(exps)

            ## rcn.freq
            rcn.freq <- cn@rcn.freq[,cts.in.rcn] # nrow = n1
            rf <- t(apply(rcn.freq,1,function(x)x/sum(x)))

            ## dist
            dist <- cn@dist2tumorBorder # length = n1

            ## data type
            if(missing(data.type)){
              data.type="ct_exp"
            }

            if(data.type=="ct"){
              df <- rcn.freq
            }else if(data.type=="ct_exp"){
              ## combine the data
              # df <- cbind(exps.imp,rcn.freq,dist)
              df <- cbind(exps.imp,rcn.freq)
            }

            df.sele <- df[is.selected,]
            norm.df.sele <- scale(df.sele,center=TRUE,scale=TRUE)

            ## clustering & classification
            set.seed(seed)

            cat("Clustering with Mclust ...\n")
            mem.ori <- mclust::Mclust(data=norm.df.sele,G=g,modelNames="EII")$classification

            if(extrapolate){
              cat("Training MclustDA ...\n")
              mc1 <- mclust::MclustDA(data=norm.df.sele,class=mem.ori,
                                      G=g,modelNames="EII",modelType = "EDDA") # 100, EII
              mem.sele <- factor(predict(mc1)$classification)
              g1 <- nlevels(mem.sele)
            }else{
              mem.sele <- mem.ori
              g1 <- nlevels(mem.sele)
            }

            ## mean counts & frequencies
            df.sele1 <- cbind(dist[is.selected],df.sele)
            colnames(df.sele1)[1] <- "dist"
            mean.df.sele <- sapply(seq(g1),function(i)colMeans(df.sele1[mem.sele==i,]))
            o <- order(mean.df.sele["dist",],decreasing=T)

            ## update labels
            mem.sele1 <- factor(as.numeric(factor(mem.sele,levels=levels(mem.sele)[o])))
            mean.df.sele1 <- sapply(seq(g1),function(i)colMeans(df.sele1[mem.sele1==i,]))

            boxplot(df.sele1[,"dist"] ~ mem.sele1,pch=NA)
            boxplot(df.sele1[,"pTBK1"] ~ mem.sele1,pch=NA)

            par(mfrow=c(2,1))
            boxplot(df.sele1[,"cCaspase3"] ~ mem.sele1,pch=NA)
            boxplot(df.sele1[,"CD8T"] ~ mem.sele1,pch=NA)
            # boxplot(df.sele1[,"pTBK1"] ~ mem.sele1,pch=NA)
            par(mfrow=c(2,1))
            boxplot(df.sele1[,"cCaspase3"] ~ mem.sele1,pch=NA)
            boxplot(df.sele1[,"BCLXL"] ~ mem.sele1,pch=NA)

            plot(t(mean.df.sele1)[,c("cCasepase3","BCLXL")])

            colnames(mean.df.sele1) <- seq(g1)

            heatmap3(mean.df.sele1,Rowv=NA,Colv=NA,scale="row",balanceColor = TRUE)
            heatmap3(mean.df.sele1[1:9,],Rowv=NA,Colv=NA,scale="row",balanceColor = TRUE)
            heatmap3(mean.df.sele1[10:18,],Rowv=NA,Colv=NA,scale="none",balanceColor = FALSE)

            if(extrapolate){
              ## extrapolate clusters
              is.available <- !apply(df,1,function(x)any(is.na(x))) # 8587
              df1 <- df[is.available,]

              cat("Applying MclustDA to the entire data ...\n")
              mc.idx <- sort(rep(seq(mc.cores),length=nrow(df1)))

              mem.all <- parallel::mclapply(seq(mc.cores),function(i){
                predict(mc1,newdata=df1[mc.idx==i,cn@cts.in.rcn])$classification
              },mc.cores=mc.cores)
              mem.all <- do.call(c,mem.all)

              ## update labels
              mem.all <- factor(mem.all)
              g1 <- nlevels(mem.all)
              levels(mem.all) <- seq(g1)

              ## mean counts & frequencies
              mean.df.all <- sapply(seq(g1),function(i)colMeans(df1[mem.all==i,]))
              colnames(mean.freq.all) <- seq(g1)
            }

            ## Sort clusters based on frequency of a cell type (CD8T by default)
            if(extrapolate){
              sort.smpls <- "all"
            }else{
              sort.smpls <- "selected"
            }

            if(sort.type == "freq" & sort.smpls == "all"){
              mean.freq <- mean.freq.all
            }else if(sort.type == "freq" & sort.smpls == "selected"){
              mean.freq <- mean.freq.sele
            }

            if(sort.by=="dist"){
              o <- order(mean.freq[sort.by,],decreasing=T)
              mean.freq.sele <-  mean.freq.sele[,o]
              mean.count.sele <-  mean.count.sele[,o]
              mem.sele <- as.numeric(factor(as.character(mem.sele),levels=as.character(seq(ncol(mean.freq))[o])))
              colnames(mean.count.sele) <- colnames(mean.freq.sele) <- paste0("Rcn",seq(g1))

              mclustda$sele <- list(is.used = mclustda$sele$is.used,
                                    mem = mem.sele,
                                    mean.freq = mean.freq.sele,
                                    mean.count = mean.count.sele)

              if(extrapolate){
                mean.freq.all <-  mean.freq.all[,o]
                mean.count.all <-  mean.count.all[,o]
                mem.all <- as.numeric(factor(as.character(mem.all),levels=as.character(seq(ncol(mean.freq))[o]))) ## mem redefined
                colnames(mean.count.all) <- colnames(mean.freq.all) <- paste0("RCN",seq(g1))

                mclustda$all <- list(is.used = is.available,
                                     mem = mem.all,
                                     mean.freq = mean.freq.all,
                                     mean.count = mean.count.all)
              }
            }

            ##
            mclustda$g <- g1

            cn$mclustda <- mclustda

            return(cn)
          })
#_ -------------------------------
# fun: meanExpRCN ----

#' @title Compute Mean Expression Profiles per RCN Cluster
#'
#' @description
#' This function computes the mean expression profiles for specified cell types
#' or features within Recurrent Cell Neighborhood (RCN) clusters. It takes a data frame,
#' typically containing expression data, and computes the mean expression values
#' for each RCN cluster based on the provided RCN information from 'computeCN'.
#' The function allows you to focus on specific cell types and features and can
#' extrapolate the clustering results to the entire dataset if needed.
#'
#' @param x A data frame containing expression data, typically from a CyCIF or similar dataset.
#' @param nn An object containing RCN information, typically obtained from 'computeCN'.
#' @param cts.in.center A character vector specifying the cell typesaround which RCN was computed (e.g., "Tumor").
#' @param cts.in.rcn A character vector specifying the cell types to include in the RCN analysis.
#' @param per.ct A logical value indicating whether to compute mean expression profiles per RCN cluster (TRUE) or for the entire dataset (FALSE).
#' @param extrapolate A logical value indicating whether to extrapolate clustering results to the entire dataset (TRUE) or not (FALSE).
#'
#' @return A list of data frames containing mean expression profiles for specified cell types or features within RCN clusters.
#'
#' @details
#' The 'meanExpRCN' function calculates the mean expression profiles for the specified cell types or features within RCN clusters. The function works as follows:
#' - It takes a data frame 'x' containing expression data and an object 'nn' containing RCN information obtained from the 'computeCN' function.
#' - You can specify the 'cts.in.center' argument to select specific cell types to focus on during the analysis.
#' - The 'cts.in.rcn' argument allows you to specify the cell types to include in the RCN analysis.
#' - If 'per.ct' is set to TRUE, the function computes mean expression profiles per RCN cluster; otherwise, it computes mean expression profiles for the entire dataset.
#' - The 'extrapolate' argument determines whether to extrapolate clustering results to the entire dataset based on RCN information.
#' The function returns a list of data frames containing mean expression profiles for the specified cell types or features within RCN clusters.
#'
#' @seealso \code{\link{computeCN}}, \code{\link{tcnClust}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr %>% arrange left_join summarize_at mutate_at mutate group_by summarize
#' @importFrom tidyr spread
#' @importFrom ggplot2 ggplot aes geom_point geom_line sym syms
#'
#' @export
setGeneric("meanExpRCN", function(x,...) standardGeneric("meanExpRCN"))

#' @rdname meanExpRCN
#' @export
setMethod("meanExpRCN","data.frame",
          function(x,
                   nn,
                   cts.in.center="Tumor",
                   cts.in.rcn=levels(cell_types(x)$cell_types)[1:10],
                   per.ct=TRUE,
                   extrapolate=TRUE){

    g <-   nn$mclustda$g

    lin.abs <- names(x@cell_types$default@cell_lineage_def[-(1:2)])
    cst.abs <- names(x@cell_types$default@cell_state_def)
    all.abs <- unique(c(lin.abs,cst.abs))

    if(extrapolate){
      mclustda <- nn$mclustda$all
    }else{
      mclustda <- nn$mclustda$sele
    }

    # clusts <- factor(mclustda$mem)
    clusts <- factor(paste0("RCN",mclustda$mem),labels=paste0("RCN",seq(g)))
    df1 <- cell_types(x) %>%
      filter(nn$within.rois) %>%
      tibble::rowid_to_column("idx") %>% ## idx is numbered within 'within.rois'
      dplyr::left_join(
        data.frame(tcn = clusts) %>%
          mutate(idx=which(mclustda$is.used)),by="idx") %>%
      dplyr::arrange(idx) %>%
      dplyr::left_join(exprs(x,type="log") %>%
                  filter(nn$within.rois) %>%
                  tibble::rowid_to_column("idx"),by="idx") %>%
      rename(idx.all="idx")## 1325874, within ROIs

    idx.nonna.tum <-!is.na(df1$tcn) & df1$cell_types %in% cts.in.center # 471910 / 1325874
    df.tum <- df1 %>%
      filter(idx.nonna.tum) %>%
      tibble::rowid_to_column("idx.tum") # 471910

    ## lst.nn: convert ids in each sample to ids in all samples
    nn.ids <- nn$nn$id[!is.na(df1$tcn)] # 1325311: 563 don't have proper RCNs.
    nn.tum.ids <- nn$nn$id[idx.nonna.tum] # 471910

    ## all unique ids per tcn cluster
    if(per.ct){
      lst.mean.exp <- tapply(nn.tum.ids,df.tum$tcn,function(this.ids){
        this.ids1 <- unique(sort(unlist(this.ids)))
        df.tmp <- df1[this.ids1,] %>%
          group_by(cell_types) %>%
          summarize_at(vars(!!!syms(all.abs)),mean,na.rm=TRUE)
        return(df.tmp)
      })
      df.exp <- as.data.frame(data.table::rbindlist(lapply(seq(lst.mean.exp),function(i){
        lst.mean.exp[[i]] %>%
          mutate(tcn=paste0("RCN",i))
      }))) %>%
        mutate(tcn = factor(tcn,levels=paste0("RCN",seq(g))))
      mes <- lapply(all.abs,function(ab){
        df.exp %>%
          select(cell_types,tcn,!!sym(ab)) %>%
          tidyr::spread(tcn,!!sym(ab))
      })
      names(mes) <- all.abs
      return(mes)
    }else{
      lst.mean.exp <- tapply(nn.tum.ids,df.tum$tcn,function(this.ids){
        this.ids1 <- unique(sort(unlist(this.ids)))
        df.tmp <- df1[this.ids1,] %>%
          group_by() %>%
          summarize_at(vars(!!!syms(all.abs)),mean,na.rm=TRUE)
        return(df.tmp)
      })

      mes <- as.data.frame(do.call(rbind,lst.mean.exp))
      rownames(mes) <- names(lst.mean.exp)
      colnames(mes) <- all.abs
      return(mes)
    }
  }
)

#_ -------------------------------
# fun: unused (non-exported) ----
dist2tgt <- function(this.cts, xy, seed.tum=1,n.tumors=100,summary=TRUE){
  if(sum(this.cts=="Tumor",na.rm=T)==0){
    stop("some of this.cts should be 'Tumor'")
  }

  has.cts <- !is.na(this.cts)
  uniq.cts <- levels(this.cts)

  lst.xys <- lapply(uniq.cts,function(ct)xy[!is.na(this.cts) & this.cts==ct,])
  names(lst.xys) <- uniq.cts

  ## sampling tumor cellsa
  set.seed(seed.tum)
  n.tumors.all <- nrow(lst.xys$Tumor)
  n.tumors <- min(n.tumors,n.tumors.all)
  idx1 <- sample(n.tumors.all,n.tumors,replace=FALSE)

  xys.tumors <- lst.xys$Tumor[idx1,]

  tgts <- names(lst.xys)
  tgts <- tgts[tgts != "Tumor"]
  dist2tumor <- sapply(tgts,function(tgt){
    if(nrow(lst.xys[[tgt]]) == 0){
      return(c(nn.dist=Inf))
    }
    lst.gds <- RANN::nn2(data=lst.xys[[tgt]],query=xys.tumors,k=1,searchtype="standard") #distance from every tumor to cd8t
    if(summary){
      med.d2tum <- median(lst.gds$nn.dists)
      return(med.d2tum)
    }else{
      return(lst.gds$nn.dists)
    }
  })
  return(dist2tumor)
}

# fun: stat.dist2tgt (immature) ----

stat.dist2tgt <- function(this.cts,xy,
                          seed.tum=1,n.tumors=100,summary=TRUE,
                          n.iter=1e4,seed.iter=123){

  i1 <- which(this.cts!="outOfROI")
  this.cts <- this.cts[i1]
  xy <- xy[i1,]

  ## distance - original
  cts.orig <- this.cts

  orig.d <- dist2tgt(cts.orig,xy,n.tumors=n.tumors,seed.tum=seed.tum,summary=summary)

  ## sampling
  ds <- do.call(rbind,parallel::mclapply(seq(n.iter),function(i){
    set.seed(i)
    cts <- sample(this.cts)
    ds <- dist2tgt(cts,xy,n.tumors=n.tumors,seed.tum=seed.tum,summary=summary)
    return(ds)
  },mc.cores=15))

  med.ds <- apply(ds,2,median)

  folds <- orig.d/med.ds
  return(folds)
}

