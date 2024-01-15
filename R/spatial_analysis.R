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
#' 6. Determines neighbors using frNN (Fast Regular Nearest Neighbors) only for tumor cells.
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
    dbs1 <- dbscan(xy.tumor, eps=dth, minPts = min.pts, weights = NULL, borderPoints = TRUE)
    cls <- dbs1$cluster

    ## Find neighbors using frNN
    frt <- frNN(xy.tumor,eps=dth)
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
          function(x, n.cores = 7, dth, minPts = 3, concavity = 0.8, ...) {
    cat("Processing ",n," ... \n")

    this.cts <- cell_types(x)$cell_types # 89174

    ## coordinates: x => xy1, xy.sp1
    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    xy.sp <- sp::SpatialPoints(xy)

    wr <- x@within_rois
    xy1 <- xy[wr,]
    xy.sp1 <- xy.sp[wr,]

    this.cts1 <- this.cts[wr]
    is.tumor <- this.cts1 == "Tumor"

    xy2 <- xy1[is.tumor,]
    dbs1 <- dbscan::dbscan(xy2, eps=dth, minPts = minPts, weights = NULL, borderPoints = TRUE) # min 2 cells together
    cls3 <- dbs1$cluster
    nls3 <- unique(cls3) # unique cluster labels
    nls3 <- sort(nls3[nls3!=0])
    names(nls3) <- nls3

    idx <- cls3 >0
    cls3[idx] <- nls3[cls3[idx]] # rename clusters

    cat(" Computing borders ... \n")

    borders <- parallel::mclapply(seq(nls3),function(i){
      xyt <- xy2[cls3==nls3[i],]
      conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = concavity, length_threshold = dth))
      names(conc) <- c("X","Y")
      return(conc)
    },mc.cores=n.cores)
    names(borders) <- nls3

    # Identify overlapping clusters and merge their points
    overlaps <- sapply(borders,function(cluster_i){
      sapply(borders,function(cluster_j){
        any(is_point_inside_polygon(cluster_j, cluster_i))
      })
    })
    diag(overlaps) <- 0

    # update clusters
    ol <- as.data.frame(which(overlaps==1,arr.ind=T))
    names(ol) <- c("child","parent")
    ol1 <- ol

    is.ol <- !seq(max(nls3)) %in% ol1[,1]
    borders <- borders[is.ol]
    nls3 <- names(borders)

    # xy.sp1, (is.in.1), nls3, borders1
    cat(" Computing is.in.1 ... \n")

    is.in <- parallel::mclapply(nls3,function(i){
      is.in <- sp::point.in.polygon(xy1$X,xy1$Y,
                                    borders[[i]]$X,
                                    borders[[i]]$Y)==1
      return(is.in)
    },mc.cores=n.cores)
    is.in.1 <- rowSums(do.call(cbind,is.in))>0

    ## distance from each cell to tumor_border
    cat(" Computing sp.polys ... \n")

    lst.polys <- lapply(as.character(nls3),
                        function(x)sp::Polygons(list(sp::Polygon(borders[[x]])),ID=x))
    sp.polys <- sp::SpatialPolygons(lst.polys)

    # sp.polys => sp.points
    cat(" Computing sp.points ... \n")

    vertices <- lapply(sp.polys@polygons, function(p) p@Polygons[[1]]@coords)
    combined_vertices <- do.call(rbind, vertices)
    sp.points <- sp::SpatialPoints(combined_vertices)

    # distance from polygon to each point
    cat(" Computing gDistance ... \n")

    gds <- parallel::mclapply(seq(xy.sp1),function(i)
      if(is.in.1[i]){
        -rgeos::gDistance(xy.sp1[i],sp.points) # inside polygon gives 0
      }else{
        rgeos::gDistance(xy.sp1[i],sp.polys)
      },mc.cores=n.cores)

    gds <- unlist(gds)

    return(list(dist=gds,sp.polys=sp.polys))
  }
)

# Function to check if a point is inside a polygon using 'point.in.polygon()'
# not exported

is_point_inside_polygon <- function(point, polygon) {
  sp::point.in.polygon(point$X, point$Y, polygon$X, polygon$Y)
}

if(0){
  is_point_inside_polygon_sf <- function(point, polygon) {
    # Create an sf point object
    point_sf <- sf::st_sfc(sf::st_point(c(point$X, point$Y)), crs = NA_crs_)

    # Create an sf polygon object
    polygon_coords <- matrix(c(polygon$X, polygon$Y), ncol = 2, byrow = FALSE)
    polygon_sf <- sf::st_sfc(st_polygon(list(polygon_coords)), crs = NA_crs_)

    # Check if the point is inside the polygon
    sf::st_intersects(point_sf, polygon_sf, sparse = FALSE)[1, 1]
  }

  #example
  point <- list(X = 2, Y = 2)

  # Define a polygon
  # Note: The first and last coordinates should be the same to close the polygon
  polygon <- list(X = c(1, 1, 3, 3, 1), Y = c(1, 3, 3, 1, 1))

  # Call the function
  system.time({
    for(i in 1:1e4){
      is_inside <- is_point_inside_polygon_sf(point, polygon)
    }# 6.729 sec
  })
  system.time({
    for(i in 1:1e4){
      is_inside <- is_point_inside_polygon(point, polygon)
    } # 0.023 sec
  })
  print(is_inside)

}
# Print the result
#_ -------------------------------

# fun: computeArea (for density) ----

#' Compute the Area of Tumor Regions in CyCIF Data.
#'
#' This function calculates the total area of tumor regions within a CyCIF dataset based on the specified distance threshold for tumor border detection.
#'
#' @param x A CyCIF object.
#' @param dth Numeric, the distance threshold used for tumor border detection.
#' @param unit Character, the unit of the computed area. Default is "mm2" (square millimeters).
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
#' @importFrom parallel mclapply
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
    dbs1 <- dbscan::dbscan(xy1, eps=dth, minPts = 3, weights = NULL, borderPoints = TRUE) # min 2 cells together
    cls3 <- dbs1$cluster

    ## Merge overlapped clustersFind tumor borders
    nls3 <- unique(cls3) # unique cluster label
    nls3 <- sort(nls3[nls3!=0])

    mcc <- min(15,length(nls3))
    concs3 <- parallel::mclapply(seq(nls3),function(i){
      xyt <- xy1[cls3==nls3[i],]
      conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = .8, length_threshold = dth))
      names(conc) <- c("X","Y")
      return(conc)
    },mc.cores=mcc)

    if(length(concs3)>1){
      # Identify overlapping clusters and merge their points
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

# fun: computeRCN ----

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
#' @return A list containing the following components:
#' - 'within.rois': A logical vector indicating whether each cell is within a region of interest (ROI).
#' - 'cts.in.center': A character vector specifying the cell types around which RCN is computed.
#' - 'cts.in.rcn': A character vector specifying the cell types considered when computing RCN values.
#' - 'n.cells.selected': The number of cells selected for RCN analysis.
#' - 'frnn': A list with recurrent Neighborhood information, including 'dist', 'id', 'eps', and 'sort'.
#' - 'exp': A data frame containing expression data for selected cells.
#' - 'is.selected': A logical vector indicating whether each cell is selected for RCN analysis.
#' - 'rcn.count': A data frame containing the counts of neighboring cell types.
#' - 'rcn.freq': A data frame containing the relative frequencies of neighboring cell types.
#'
#' @seealso \code{\link{cyApply}}, \code{\link{dbscan::frNN}}
#'
#' @importFrom dbscan frNN
#' @importFrom data.table := rbindlist
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble rowid_to_column
#'
#' @rdname computeRCN
#' @export
setGeneric("computeRCN", function(x,...) standardGeneric("computeRCN"))

#' @rdname computeRCN
#' @export
setMethod("computeRCN", "Cycif",
          function(x,r,unit=c("pixel","um"),
                     cts.in.center,
                     cts.in.rcn,
                     n.sampling,
                     seed=123){ # only for Cycif, and CycifStack
    if(missing(r)){
      stop("radius (r) should be specified")
    }else if(unit[1]=="um"){
      r <- r/0.65
      unit <- "pixel"
    }

    xy <- xys(x)
    ymax <- max(xy$Y_centroid)
    xy$Y_centroid <- ymax - xy$Y_centroid
    cts <- cell_types(x)

    wr <- x@within_rois
    xy1 <- xy[wr,]
    cts1 <- cts[wr,]
    lev.cts <- levels(cts1$cell_types)
    lev.cts <- lev.cts[lev.cts != "outOfROI"]

    if(missing(cts.in.center)){
      cts.in.center <- lev.cts
    }

    if(missing(cts.in.rcn)){
      cts.in.rcn <- lev.cts
    }

    ## frnn object - all within_rois samples
    frnn <- dbscan::frNN(xy1,eps=r,bucketSize=10) # 265793 cells
    frnn <- new("frNN",
                dist = frnn$dist,
                id = frnn$id,
                eps = frnn$eps,
                sort = frnn$sort)

    frnn.ids <- frnn$id
    n.frnn <- sapply(frnn.ids,length)
    n0 <- which(n.frnn==0)
    idx.df <- data.frame(center_id=rep(seq(frnn.ids),n.frnn),
                         neighbor_id=.Internal(unlist(frnn.ids, FALSE, FALSE))) %>%
      data.table::as.data.table()

    ##
    df1 <- cell_types(x) %>%
      cbind(exprs(x,type="log")) %>%
      filter(cell_types != "outOfROI") %>%
      # tibble::rowid_to_column("center_id") %>%
      data.table::data.table()

    rn <- levels(df1$cell_types)
    csts <- x@cell_types$default@cell_state_def[rn,]

    lst.csts <- lapply(csts,function(x1){
      rn[!is.na(x1) & x1 == "CAN"]
    })
    cst.abs <- names(csts)

    for(ab in cst.abs){
      this.cts <- lst.csts[[ab]]
      df1 <- df1[!cell_types %in% this.cts ,(ab) := NA]
    }

    ##
    df1 <- cbind(idx.df,df1[idx.df$neighbor_id,])
    dt <-  df1[ ,lapply(.SD, mean, na.rm=T), by=center_id,.SDcols=cst.abs]
    dt0 <- data.table::as.data.table(cbind(n0,array(NA,c(length(n0),length(cst.abs)))))
    names(dt0) <- names(dt)

    dt <- data.table::rbindlist(list(dt,dt0))[order(center_id),]
    has.neighbors <- !seq(nrow(dt)) %in% n0

    ## within positive ROIs + has neighbors
    selected.ids1 <- which(cts1$cell_types %in% cts.in.center &
                          sapply(frnn$id,length)>0 &
                          sapply(frnn$id,function(ids){
                            any(cts1$cell_types[ids] %in% cts.in.rcn)
                          }) &
                          has.neighbors
                        )

    ## sampling
    set.seed(seed)
    n.this <- length(selected.ids1)
    n.cts <- min(n.sampling,n.this)

    selected.ids2 <- sample(selected.ids1,n.cts) # idx after ROI filter

    ## selected ids among avialble focused celltypes (eg tumor cells)
    sid <- seq(frnn$id) %in% selected.ids2

    ## cell type frequency (not absolute count)
    rcn.freq <- t(sapply(frnn$id,function(ids){ # [c(83538,101787)]
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      ntab <- tab/sum(tab)
      return(ntab)
    }))
    rcn.count <- t(sapply(frnn$id,function(ids){ # [c(83538,101787)]
      ct <- cts1$cell_types[ids]
      tab <- table(ct)
      return(tab)
    }))

    return(list(within.rois=wr,
                cts.in.center=cts.in.center,
                cts.in.rcn=cts.in.rcn,
                n.cells.selected=n.cts,
                frnn=frnn,
                exp=dt,
                is.selected=sid,
                rcn.count=rcn.count,
                rcn.freq=rcn.freq))

})

#' @rdname computeRCN
#' @export
setMethod("computeRCN", "CycifStack",
  function(x,r,unit=c("pixel","um"),
           cts.in.center,
           cts.in.rcn,
           n.sampling,
           seed=123){ # only for Cycif, and CycifStack

    cat("Get neighbors ...\n")
    frnn1 <- cyApply(x,function(cy){
      cat(names(cy),"\n")
      computeRCN(x=cy,r=r,unit=c("pixel","um"),
        cts.in.center=cts.in.center,
        cts.in.rcn=cts.in.rcn,
        n.sampling=n.sampling,
        seed=seed)
    })

    cat("Restructure data ...\n")
    ## assemble frnn
    lst.frnn <- lapply(frnn1,function(fr)fr$frnn)

    ### dists
    dists <- do.call(c,lapply(lst.frnn,function(frnn)frnn$dist))

    ### ids
    n.frnns <- sapply(lst.frnn,function(frnn)length(frnn$id)) # 1325874, all cells, excluding outOfROIs
    n.frnns.pre <- c(0,cumsum(n.frnns)[-length(n.frnns)])
    names(n.frnns.pre) <- names(x)

    ## frnn.ids, frnn.ids1, frnn.tum.ids - list of neighboring cells ids for tumors used for the frnn analysis
    frnn.ids <- lapply(names(x),function(nm){
      x <- lst.frnn[[nm]]
      n.prior <- n.frnns.pre[nm]
      this.ids <- x$id
      new.ids <- lapply(this.ids,function(id){
        new.id <- id + n.prior
        return(new.id)
      })
      return(new.ids)
    })
    frnn.ids1 <- do.call(c,frnn.ids) ## 1325874, now all data are combined - and the indices are after excluding outOfROIs

    ### eps, sort
    eps <- unique(sapply(lst.frnn,function(frnn)frnn$eps))
    sort <- unique(sapply(lst.frnn,function(frnn)frnn$sort))

    frnn <- list(dist=dists,id=frnn.ids1,eps=eps,sort=sort)

    ## within.rois
    within.rois <- do.call(c,lapply(frnn1,function(fr)fr$within.rois))

    ## cts.in.center
    cts.in.center <- unique(as.vector(sapply(frnn1,function(fr)fr$cts.in.center)))

    ## cts.in.rcn
    cts.in.rcn <- frnn1[[1]]$cts.in.rcn

    ## n.cells.selected
    n.cells.selected <- sapply(frnn1,function(fr)fr$n.cells.selected)
    v.smpls <- rep(names(n.cells.selected),n.cells.selected)

    ## is.selected
    mclustda <- list()
    mclustda$sele <- mclustda$all <- list()
    mclustda$sele$is.used <- do.call(c,lapply(frnn1,function(fr)fr$is.selected)) #67119 -> 66904, how did I do that?

    ##
    exps <- data.table::rbindlist(lapply(frnn1,function(fr)as.data.frame(fr$exp)))
    rcn.count <- data.table::rbindlist(lapply(frnn1,function(fr)as.data.frame(fr$rcn.count)))
    rcn.freq <- data.table::rbindlist(lapply(frnn1,function(fr)as.data.frame(fr$rcn.freq)))

    return(list(within.rois=within.rois, ## all cells (including outOfROIs)
                cts.in.center=cts.in.center, ## 1 or a few
                cts.in.rcn=cts.in.rcn, ## 1 or a few
                n.cells.selected=n.cells.selected,
                v.smpls=v.smpls,
                frnn=frnn,
                mclustda=mclustda, #is.selected=is.selected, <= included in mclustda
                exp=exps,
                rcn.count=rcn.count,
                rcn.freq=rcn.freq))
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
#' @param frnn An object containing RCN information, typically obtained from 'computeRCN'.
#' @param g The number of clusters to create.
#' @param seed The random seed for reproducibility.
#' @param sort.by The cell type to sort clusters by (e.g., "CD8T").
#' @param sort.type The type of data to use for sorting: "freq" (relative frequencies) or "count" (counts).
#' @param sort.smpls The subset of samples to use for sorting: "all" (entire dataset) or "selected" (selected cells).
#' @param data.type The type of data to use for clustering: "ct_exp" (cell type and expression data) or "ct" (cell type data only).
#' @param extrapolate Whether to extrapolate clusters to the entire dataset (TRUE) or not (FALSE).
#' @param mc.cores The number of CPU cores to use for parallel processing.
#'
#' @return An updated 'frnn' object with clustering and sorting information.
#'
#' @details
#' The `tcnClust` function uses the provided `frnn` object to perform clustering and classification of cells based on their neighborhood relationships. It allows you to specify the number of clusters (`g`), the cell type to sort clusters by (`sort.by`), and other clustering parameters.
#' The clustering process results in the classification of cells into distinct clusters, and the function provides information about these clusters, including the mean frequencies, counts, and more.
#' By specifying different options for `sort.by`, `sort.type`, and `sort.smpls`, you can customize the sorting behavior of clusters based on cell types and data types.
#' Additionally, you can choose to extrapolate clusters to the entire dataset using the `extrapolate` argument, which can be helpful for analyzing the overall dataset.
#'
#' @seealso \code{\link{computeRCN}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom dbscan frNN
#' @importFrom data.table as.data.table rbindlist
#' @importFrom parallel mclapply
#' @export
setGeneric("tcnClust", function(frnn,...) standardGeneric("tcnClust"))

#' @rdname tcnClust
#' @export
setMethod("tcnClust","data.frame",
  function(frnn,
           g=50,
           seed=123,
           sort.by="CD8T",
           sort.type=c("freq","count"),
           sort.smpls=c("all","selected"),
           data.type=c("ct_exp","ct"),
           extrapolate=FALSE,
           mc.cores=1){
  mclustda <- frnn$mclustda

  exps <- as.matrix(frnn$exp)[,-1]
  exps.imp <- imputeData(exps)

  this.cts <- cts.in.rcn

  mat.count.all <- as.matrix(frnn$rcn.count)[,this.cts]
  mat.freq.all <- t(apply(as.matrix(frnn$rcn.freq)[,this.cts],1,function(x)x/sum(x)))

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
      predict(mc1,newdata=mat.freq.all1[mc.idx==i,frnn$cts.in.rcn])$classification
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

  frnn$mclustda <- mclustda

  return(frnn)
})


# fun: rcnClust ----

#' Cluster and Sort Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack Objects
#'
#' This function clusters and sorts the Recurrent Cell Neighbors (RCN) for CyCIF or CyCIFStack objects.
#' It clusters cells based on their RCN profiles, sorts clusters based on the specified cell type,
#' and optionally extrapolates the clustering to the entire dataset.
#'
#' @param frnn An object containing RCN information, typically obtained from 'computeRCN'.
#' @param g The number of clusters to create.
#' @param seed The random seed for reproducibility.
#' @param sort.by The cell type to sort clusters by (e.g., "CD8T").
#' @param sort.type The type of data to use for sorting: "freq" (relative frequencies) or "count" (counts).
#' @param sort.smpls The subset of samples to use for sorting: "all" (entire dataset) or "selected" (selected cells).
#' @param data.type The type of data to use for clustering: "ct_exp" (cell type and expression data) or "ct" (cell type data only).
#' @param extrapolate Whether to extrapolate clusters to the entire dataset (TRUE) or not (FALSE).
#' @param mc.cores The number of CPU cores to use for parallel processing.
#'
#' @return An updated 'frnn' object with clustering and sorting information.
#'
#' @seealso \code{\link{computeRCN}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom dbscan frNN
#' @importFrom data.table as.data.table rbindlist
#' @importFrom parallel mclapply
#' @export
setGeneric("rcnClust", function(frnn,...) standardGeneric("rcnClust"))

#' @rdname rcnClust
#' @export
setMethod("rcnClust","data.frame",
          function(frnn,
                   g=50,
                   seed=123,
                   sort.by="CD8T",
                   sort.type=c("freq","count"),
                   sort.smpls=c("all","selected"),
                   data.type=c("ct_exp","ct"),
                   extrapolate=FALSE,
                   mc.cores=1){
            mclustda <- frnn$mclustda

            exps <- as.matrix(frnn$exp)[,-1]
            exps.imp <- imputeData(exps)

            this.cts <- cts.in.rcn

            mat.count.all <- as.matrix(frnn$rcn.count)[,this.cts]
            mat.freq.all <- t(apply(as.matrix(frnn$rcn.freq)[,this.cts],1,function(x)x/sum(x)))

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
                predict(mc1,newdata=mat.freq.all1[mc.idx==i,frnn$cts.in.rcn])$classification
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

            frnn$mclustda <- mclustda

            return(frnn)
          })
#_ -------------------------------
# fun: meanExpRCN ----

#' @title Compute Mean Expression Profiles per RCN Cluster
#'
#' @description
#' This function computes the mean expression profiles for specified cell types
#' or features within Recurrent Cell Neighborhood (RCN) clusters. It takes a data frame,
#' typically containing expression data, and computes the mean expression values
#' for each RCN cluster based on the provided RCN information from 'computeRCN'.
#' The function allows you to focus on specific cell types and features and can
#' extrapolate the clustering results to the entire dataset if needed.
#'
#' @param x A data frame containing expression data, typically from a CyCIF or similar dataset.
#' @param frnn An object containing RCN information, typically obtained from 'computeRCN'.
#' @param cts.in.center A character vector specifying the cell typesaround which RCN was computed (e.g., "Tumor").
#' @param cts.in.rcn A character vector specifying the cell types to include in the RCN analysis.
#' @param per.ct A logical value indicating whether to compute mean expression profiles per RCN cluster (TRUE) or for the entire dataset (FALSE).
#' @param extrapolate A logical value indicating whether to extrapolate clustering results to the entire dataset (TRUE) or not (FALSE).
#'
#' @return A list of data frames containing mean expression profiles for specified cell types or features within RCN clusters.
#'
#' @details
#' The 'meanExpRCN' function calculates the mean expression profiles for the specified cell types or features within RCN clusters. The function works as follows:
#' - It takes a data frame 'x' containing expression data and an object 'frnn' containing RCN information obtained from the 'computeRCN' function.
#' - You can specify the 'cts.in.center' argument to select specific cell types to focus on during the analysis.
#' - The 'cts.in.rcn' argument allows you to specify the cell types to include in the RCN analysis.
#' - If 'per.ct' is set to TRUE, the function computes mean expression profiles per RCN cluster; otherwise, it computes mean expression profiles for the entire dataset.
#' - The 'extrapolate' argument determines whether to extrapolate clustering results to the entire dataset based on RCN information.
#' The function returns a list of data frames containing mean expression profiles for the specified cell types or features within RCN clusters.
#'
#' @seealso \code{\link{computeRCN}}, \code{\link{tcnClust}}
#'
#' @importFrom mclust Mclust
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr %>% arrange left_join summarize_at mutate_at mutate group_by summarize
#' @export
setGeneric("meanExpRCN", function(x,...) standardGeneric("meanExpRCN"))

#' @rdname meanExpRCN
#' @export
setMethod("meanExpRCN","data.frame",
          function(x,
                   frnn,
                   cts.in.center="Tumor",
                   cts.in.rcn=levels(cell_types(x)$cell_types)[1:10],
                   per.ct=TRUE,
                   extrapolate=TRUE){

    g <-   frnn$mclustda$g

    lin.abs <- names(x@cell_types$default@cell_lineage_def[-(1:2)])
    cst.abs <- names(x@cell_types$default@cell_state_def)
    all.abs <- unique(c(lin.abs,cst.abs))

    if(extrapolate){
      mclustda <- frnn$mclustda$all
    }else{
      mclustda <- frnn$mclustda$sele
    }

    # clusts <- factor(mclustda$mem)
    clusts <- factor(paste0("RCN",mclustda$mem),labels=paste0("RCN",seq(g)))
    df1 <- cell_types(x) %>%
      filter(frnn$within.rois) %>%
      tibble::rowid_to_column("idx") %>% ## idx is numbered within 'within.rois'
      dplyr::left_join(
        data.frame(tcn = clusts) %>%
          mutate(idx=which(mclustda$is.used)),by="idx") %>%
      dplyr::arrange(idx) %>%
      dplyr::left_join(exprs(x,type="log") %>%
                  filter(frnn$within.rois) %>%
                  tibble::rowid_to_column("idx"),by="idx") %>%
      rename(idx.all="idx")## 1325874, within ROIs

    idx.nonna.tum <-!is.na(df1$tcn) & df1$cell_types %in% cts.in.center # 471910 / 1325874
    df.tum <- df1 %>%
      filter(idx.nonna.tum) %>%
      tibble::rowid_to_column("idx.tum") # 471910

    ## lst.frnn: convert ids in each sample to ids in all samples
    frnn.ids <- frnn$frnn$id[!is.na(df1$tcn)] # 1325311: 563 don't have proper RCNs.
    frnn.tum.ids <- frnn$frnn$id[idx.nonna.tum] # 471910

    ## all unique ids per tcn cluster
    if(per.ct){
      lst.mean.exp <- tapply(frnn.tum.ids,df.tum$tcn,function(this.ids){
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
      lst.mean.exp <- tapply(frnn.tum.ids,df.tum$tcn,function(this.ids){
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

