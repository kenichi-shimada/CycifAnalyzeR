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
    xy.sf <- st_as_sf(xy, coords = c("X_centroid", "Y_centroid"), crs = NA)

    wr <- x@within_rois
    xy.sf1 <- xy.sf[wr,]
    this.cts1 <- this.cts[wr]
    is.tumor <- this.cts1 == "Tumor"

    xy.sf2 <- xy.sf1[is.tumor,]
    xy2 <- st_coordinates(xy.sf2)

    dbs1 <- dbscan::dbscan(xy2, eps=dth, minPts = minPts, weights = NULL, borderPoints = TRUE)
    cls3 <- dbs1$cluster
    nls3 <- unique(cls3)
    nls3 <- sort(nls3[nls3!=0])
    names(nls3) <- nls3

    idx <- cls3 >0
    cls3[idx] <- nls3[cls3[idx]]

    cat(" Computing borders ... \n")

    borders <- parallel::mclapply(seq(nls3),function(i){
      xyt <- xy.sf2[cls3==nls3[i],]
      conc <- concaveman::concaveman(xyt, concavity = concavity, length_threshold = dth)
      return(conc)
    },mc.cores=n.cores)
    names(borders) <- nls3

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

    bb <- st_make_grid(b1,cellsize=1000)
    if(0){
      plot(bb)
      plot(b1,add=T,col=2)
    }

    is.in.b <- st_intersects(b1,bb,sparse=FALSE)

    min.bb.idx <- apply(is.in.b,1,function(i)min(which(i)))
    bb.idx <- tapply(seq(min.bb.idx),min.bb.idx,identity)
    used.bb.idx <- names(bb.idx)
    near.feature <- st_nearest_feature(xy.sf1,b1)
    near.idx <- tapply(seq(near.feature),near.feature,identity)

    sgn.d <- rep(NA,length(near.feature))
    for(bbi in used.bb.idx){
      border.i <- bb.idx[[bbi]]
      point.i <- unlist(near.idx[border.i])
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
