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
           cancer.cts = c("Cancer"),ct_name="default",
           dth, concavity = 0.8, plot = FALSE,...) {
    this.cts <- cell_types(x,ct_name=ct_name)$cell_types # 89174

    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    xy.sf <- st_as_sf(xy, coords = c("X_centroid", "Y_centroid"), crs = NA)

    wr <- x@within_rois
    xy.sf1 <- xy.sf[wr,]
    this.cts1 <- this.cts[wr]
    xy1 <- sf::st_coordinates(xy.sf1)

    cat("Identifying tissue segments  ... \n")
    dbs2 <- dbscan::dbscan(xy1, eps=dth*5, minPts = 3, weights = NULL, borderPoints = TRUE)
    cls2 <- dbs2$cluster
    nls2 <- unique(cls2)
    nls2 <- sort(nls2[nls2!=0])
    names(nls2) <- nls2
    n.cls2 <- table(cls2)
    used.cls2 <- names(which(n.cls2 > n.cells.per.seg))
    used.cls2 <- used.cls2[used.cls2 != "0"]
    idx2 <- cls2 %in% used.cls2

    cls2.1 <- cls2
    cls2.1[!idx2] <- 0
    cls2.1[idx2] <- as.numeric(factor(cls2.1[idx2]))
    nls2.1 <- unique(cls2.1)
    nls2.1 <- nls2.1[nls2.1!=0]

    cat(max(nls2.1)," segments identified\n")
    cat(" Computing borders for each segment ... \n")

    borders2 <- parallel::mclapply(seq(nls2.1),function(i){
      requireNamespace("sf", quietly = TRUE)
      requireNamespace("concaveman", quietly = TRUE)

      xyt <- xy.sf1[cls2.1==nls2.1[i],]
      conc <- concaveman::concaveman(xyt, concavity = concavity, length_threshold = dth*5)
      return(conc)
    },mc.cores=n.cores)
    if(class(borders2)=="try-error"){
      stop("Error in concaveman::concaveman")
    }
    names(borders2) <- nls2.1

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

    cat("Identifying borders for tumor regions for each segment ... \n")

    lst.segs <- list()
    for (seg.i in nls2.1){
      cat("Segment",seg.i,":\n")
      this.seg <- cls2.1==nls2.1[seg.i]

      xy.sf2 <- xy.sf1[this.seg,]
      this.cts2 <- this.cts1[this.seg]
      this.cts2[is.na(this.cts2)] <- "outOfROI"

      is.tumor <- !is.na(this.cts2) & this.cts2 %in% cancer.cts

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
      xy3 <- sf::st_coordinates(xy.sf3)

      dbs1 <- dbscan::dbscan(xy3, eps=dth, minPts = 3, weights = NULL, borderPoints = TRUE)
      cls3 <- dbs1$cluster
      nls3 <- unique(cls3)
      nls3 <- sort(nls3[nls3!=0])
      names(nls3) <- nls3

      n.cls3 <- table(cls3)
      tum.reg.cls3 <- names(which(n.cls3 > n.cells.per.tumor.core))
      tum.reg.cls3 <- tum.reg.cls3[tum.reg.cls3 != "0"]
      tum.bud.cls3 <- as.character(nls3)[!nls3 %in% tum.reg.cls3]
      idx3 <- cls3 %in% tum.reg.cls3

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

      bb3 <- st_make_grid(b3.1,cellsize=1000)
      if(0){
        plot(bb3)
        plot(b3.1,add=T,col=2)
      }

      is.in.b3 <- st_intersects(b3.1,bb3,sparse=FALSE)

      min.bb.idx <- apply(is.in.b3,1,function(i)min(which(i)))
      bb.idx <- tapply(seq(min.bb.idx),min.bb.idx,identity)
      used.bb.idx <- names(bb.idx)
      near.feature <- st_nearest_feature(xy.sf2,b3.1)
      near.idx <- tapply(seq(near.feature),near.feature,identity)

      sgn.d <- rep(NA,length(near.feature))
      for(bbi in used.bb.idx){
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
      tor <- round(dth*3)

      d1 <- lst$dist$dist
      d1[d1 > tor] <- tor
      d1[d1 < -tor] <- -tor
      d1[is.na(d1)] <- tor+1

      uniq.cols <- c(colorRampPalette(brewer.pal(11,"Spectral"))(tor*2+1),"grey")
      names(uniq.cols) <- seq(-tor,tor+1)
      plot(xy.sf1[[1]],col=uniq.cols[as.character(d1)],pch=".")
      plot(lst$borders[[1]],border=1,add=T)
    }

    return(lst)
  }
)
