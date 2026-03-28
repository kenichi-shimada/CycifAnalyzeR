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
          function(x,dth,unit=c("mm2"),plot=TRUE,strict=FALSE,concavity=.8,
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

    ## need to run this becaus concaveman gives an error when run for the 1st time
    i=1
    xyt <- xy1[cls3==nls3[i],]
    conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = concavity, length_threshold = dth))
    names(conc) <- c("X","Y")
    rm(i,xyt,conc)

    concs3 <- parallel::mclapply(seq(nls3),function(i){
      xyt <- xy1[cls3==nls3[i],]
      conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = concavity, length_threshold = dth))
      names(conc) <- c("X","Y")
      return(conc)
    # }
    },mc.cores=mcc)
    names(concs3) <- nls3

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
    }else{
      cls3.updated <- cls3
      nls3.updated <- nls3
    }

    concs3.updated <- parallel::mclapply(seq(nls3.updated),function(i){
      xyt <- xy1[cls3.updated==nls3.updated[i],]
      conc <- as.data.frame(concaveman::concaveman(as.matrix(xyt), concavity = concavity, length_threshold = dth))
      names(conc) <- c("X","Y")
      return(conc)
    },mc.cores=7)
    names(concs3.updated) <- nls3.updated

    ## Interactive modifications
    satisfied <- FALSE
    while (!satisfied) {
      ## Plot the current state
      par(mar = c(5, 5, 3, 3))
      sp::plot(xy.sp1, col = cls3 + 1, pch = c(2, 20)[(cls3 != 0) + 1], cex = .5, main = "Tumors, clustered (dbscan)")
      sapply(concs3.updated, function(conc) lines(conc))

      ## Ask if user wants to zoom in
      zoom_in <- readline("Zoom in? [Y/N]: ")
      if (grepl("^[yY]", zoom_in)) {
        cat("Click bottom-left and top-right corners for zooming.\n")
        x1 <- locator(2)
        x2 <- list(x = c(x1$x[c(1, 2, 2, 1)]), y = c(x1$y[c(1, 1, 2, 2)]))

        ## Replot with zoomed-in area
        sp::plot(xy.sp1, col = cls3 + 1, pch = c(2, 20)[(cls3 != 0) + 1],
                 xlim = range(x1$x), ylim = range(x1$y), cex = .5,
                 main = "Tumors, clustered (dbscan)")

        is.in.x1 <- sapply(concs3.updated, function(cx) any(sp::point.in.polygon(cx$X, cx$Y, x2$x, x2$y) == 1))
        concs3.updated <- concs3.updated[is.in.x1]
        sapply(concs3.updated, function(conc) lines(conc))
      }

      ## Ask if user wants to remove specific points
      remove_points <- readline("Remove points based on drawn polygon? [Y/N]: ")
      if (grepl("^[yY]", remove_points)) {
        cat("Draw polygon to define points for removal. Click to complete.\n")
        rem_poly <- locator(type = "p")
        if (length(rem_poly$x) > 2) {
          rem_x <- rem_poly$x
          rem_y <- rem_poly$y
          concs3.updated <- lapply(concs3.updated, function(coords) {
            inside <- sp::point.in.polygon(coords$X, coords$Y, rem_x, rem_y) == 1
            coords[!inside, ]  # Keep points that are not inside removal region
          })
        }
      }

      concs3.updated <- Filter(function(coords) nrow(coords) > 0, concs3.updated)

      ## Initialize area status for each region
      area_status <- setNames(rep("positive", length(concs3.updated)), names(concs3.updated))

      ## Ask if user wants to remove blank spaces from drawn polygons
      remove_blank <- readline("Remove blank spaces from drawn polygon? [Y/N]: ")
      if (grepl("^[yY]", remove_blank)) {
        cat("Draw a polygon to define blank space to remove.\n")
        blank_poly <- locator(type = "p")

        if (length(blank_poly$x) > 2) {
          ## Create a valid polygon from user input
          blank_x <- c(blank_poly$x, blank_poly$x[1])
          blank_y <- c(blank_poly$y, blank_poly$y[1])
          drawn_polygon_sf <- sf::st_make_valid(sf::st_polygon(list(cbind(blank_x, blank_y))))

          ## Convert SpatialPoints to sf
          xy.sp1_sf <- st_as_sf(xy.sp1, coords = c("X", "Y"), crs = NA)

          ## Compute distances and select closest points
          distances <- sf::st_distance(xy.sp1_sf, drawn_polygon_sf)
          threshold <- quantile(distances, 0.2)
          closest_points <- xy.sp1_sf[distances <= threshold, ]

          ## Generate concave hull using alpha shape
          if (nrow(closest_points) > 2) {
            alpha_shape <- alphahull::ashape(sf::st_coordinates(closest_points)[,1],
                                             sf::st_coordinates(closest_points)[,2],
                                             alpha = 2)

            ## Extract edges safely
            if (nrow(alpha_shape$edges) > 0) {
              boundary_points <- unique(cbind(alpha_shape$edges$x1, alpha_shape$edges$y1))
              boundary_points <- rbind(boundary_points, boundary_points[1, ])  # Close polygon

              ## Convert boundary points to sf polygon
              enclosing_polygon_sf <- st_sfc(st_polygon(list(boundary_points)), crs = NA)

              ## Subtract from tumor regions
              concs3 <- lapply(names(concs3), function(region_name) {
                coords <- concs3[[region_name]]
                polygon_sf <- st_sfc(st_polygon(list(as.matrix(coords))), crs = NA)

                ## Ensure polygons are valid
                polygon_sf <- sf::st_make_valid(polygon_sf)
                enclosing_polygon_sf <- sf::st_make_valid(enclosing_polygon_sf)

                ## Perform subtraction if meaningful
                if (sf::st_area(enclosing_polygon_sf) < sf::st_area(polygon_sf)) {
                  new_coords <- sf::st_difference(polygon_sf, enclosing_polygon_sf)
                  if (!sf::st_is_empty(new_coords)) {
                    area_status[region_name] <- "negative"
                    return(as.data.frame(sf::st_coordinates(new_coords)))
                  }
                }
                return(coords)
              })

              ## Remove empty polygons
              concs3 <- Filter(function(coords) nrow(coords) > 0, concs3)

              ## Ask for user confirmation
              satisfied <- grepl("^[yY]", readline("Are the current modifications satisfactory? [Y/N]: "))
            } else {
              print("No valid edges found. Adjust alpha or include more points.")
            }
          } else {
            print("Not enough points to define a boundary. Try increasing the selection threshold.")
          }
        }
      }
      remove_blank <- readline("Remove blank spaces from drawn polygon? [Y/N]: ")
      if (grepl("^[yY]", remove_blank)) {
        cat("Draw a polygon to define blank space to remove.\n")
        blank_poly <- locator(type = "p")

        if (length(blank_poly$x) > 2) {
          ## Close the polygon by repeating the first coordinate at the end
          blank_x <- c(blank_poly$x, blank_poly$x[1])
          blank_y <- c(blank_poly$y, blank_poly$y[1])
          drawn_polygon_sf <- sf::st_polygon(list(cbind(blank_x, blank_y)))

          ## Ensure the polygon is valid
          if (!sf::st_is_valid(drawn_polygon_sf)) {
            drawn_polygon_sf <- sf::st_make_valid(drawn_polygon_sf)
          }

          # Convert SpatialPoints to sf
          xy.sp1_sf <- st_as_sf(xy.sp1, coords = c("X", "Y"), crs = NA)  # Cartesian CRS (no lat/lon)

          # Compute distances
          distances <- sf::st_distance(xy.sp1_sf, drawn_polygon_sf)

          threshold <- quantile(distances, 0.2)  # Set a threshold dynamically
          closest_points <- xy.sp1[distances <= threshold, ]  # Select only closest points

          # Generate a concave hull using alpha shape
          alpha_shape <- alphahull::ashape(closest_points$X_centroid,
                                           closest_points$Y_centroid,
                                           alpha = 2)  # Adjust alpha for tightness

          # Extract edges safely
          if (nrow(alpha_shape$edges) > 0) {
            edges <- alpha_shape$edges  # Retrieve edges
            edges_df <- data.frame(
              x1 = edges[, "x1"],
              y1 = edges[, "y1"],
              x2 = edges[, "x2"],
              y2 = edges[, "y2"]
            )

            # Plot edges
            plot(closest_points$X_centroid, closest_points$Y_centroid, col = "blue", pch = 20, main = "Alpha Shape Edges")
            # segments(edges_df$x1, edges_df$y1, edges_df$x2, edges_df$y2, col = "red", lwd = 2)
            # draw polygon drawn_polygon_sf on top of the closest_points

          } else {
            print("No edges created! Adjust alpha or use more points.")
          }
          # Extract unique boundary points
          boundary_points <- unique(cbind(edges$x1, edges$y1))  # Extract edge start points
          boundary_points <- rbind(boundary_points, boundary_points[1, ])  # Close the polygon

          # Convert to sf polygon
          alpha_polygon <- st_polygon(list(boundary_points))
          alpha_polygon_sf <- st_sfc(alpha_polygon, crs = NA)  # No CRS for Cartesian data

          # Ensure interior points are excluded
          boundary_points_sf <- sf::st_as_sf(as.data.frame(boundary_points), coords = c("X", "Y"))
          inside <- sf::st_within(boundary_points_sf, drawn_polygon_sf, sparse = FALSE)

          # Keep only the true edge points (not inside the drawn polygon)
          final_boundary_points <- boundary_points_sf[!inside, ]
          ## Convert all tumor cell coordinates to an sf object
          all_coords <- as.data.frame(xy.sp1@coords)
          colnames(all_coords) <- c("X", "Y")
          all_coords_sf <- sf::st_as_sf(all_coords, coords = c("X", "Y"), crs = 4326)
          all_coords_sf <- sf::st_transform(all_coords_sf, projected_crs)  ## Project to same CRS

          ## Compute distances from each point to the drawn polygon
          distances <- sf::st_distance(all_coords_sf, drawn_polygon_sf)

          ## Select points that are **just outside** the drawn polygon
          threshold_distance <- quantile(distances, 0.2)  # Get the 20th percentile of distance
          near_points <- all_coords[as.numeric(distances) > 0 & as.numeric(distances) <= threshold_distance, ]

          ## Compute the **alpha shape** instead of convex hull
          if (nrow(near_points) > 2) {
            alpha_shape <- alphahull::ashape(near_points$X, near_points$Y, alpha = 1)  # Alpha shape
            alpha_polygon <- alphahull::ashape2poly(alpha_shape)  # Convert to polygon
            enclosing_polygon <- as.data.frame(alpha_polygon@polygons[[1]]@Polygons[[1]]@coords)
            names(enclosing_polygon) <- c("X", "Y")
          } else {
            enclosing_polygon <- near_points  # Fallback if not enough points
          }

          ## Convert to sf polygon
          enclosing_polygon_sf <- sf::st_polygon(list(as.matrix(rbind(enclosing_polygon, enclosing_polygon[1, ]))))
          enclosing_polygon_sf <- sf::st_sfc(enclosing_polygon_sf, crs = sf::st_crs(4326))

          ## Subtract this enclosing polygon from tumor regions
          concs3 <- lapply(names(concs3), function(region_name) {
            coords <- concs3[[region_name]]
            polygon_sf <- sf::st_polygon(list(as.matrix(coords)))
            polygon_sf <- sf::st_sfc(polygon_sf, crs = sf::st_crs(4326))

            ## Ensure polygons are valid before computing difference
            if (!sf::st_is_valid(polygon_sf)) {
              polygon_sf <- sf::st_make_valid(polygon_sf)
            }
            if (!sf::st_is_valid(enclosing_polygon_sf)) {
              enclosing_polygon_sf <- sf::st_make_valid(enclosing_polygon_sf)
            }

            ## Subtract only if the enclosing polygon is meaningful
            if (sf::st_area(enclosing_polygon_sf) < sf::st_area(polygon_sf)) {
              new_coords <- sf::st_difference(polygon_sf, enclosing_polygon_sf)

              if (!sf::st_is_empty(new_coords)) {
                ## Mark region as "negative" (excluded from area)
                area_status[region_name] <- "negative"
                return(as.data.frame(sf::st_coordinates(new_coords)))
              }
            }
            return(coords)
          })

          ## Remove any empty polygons after blank space removal
          concs3 <- Filter(function(coords) nrow(coords) > 0, concs3)
      ## Ask if user is satisfied with the modifications
      satisfied <- grepl("^[yY]", readline("Are the current modifications satisfactory? [Y/N]: "))
    }

    ## Compute area
    areas <- sapply(concs3, function(coords) {
      polygon <- sp::Polygon(coords)
      return(polygon@area)
    })

    ## Convert area to mm^2
    sum.area <- sum(areas) * (0.65 * (10^-3))^2
    return(sum.area)
    }
  }
})
