#_ -------------------------------

# fun: computeCN ----

#' Compute Recurrent Cell Neighbors (RCN) for Cycif or CycifStack Objects
#'
#' This function computes the Recurrent Cell Neighbors (RCN) for Cycif or CycifStack objects.
#' RCN measures the relative frequencies of neighboring cell types around each cell within specified
#' radius 'r'. The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
#'
#' @param x A Cycif or CycifStack object.
#' @param r_um The radius within which neighboring cells are considered, in micrometers.
#' @param k The number of nearest neighbors to use when \code{type="knn"}.
#' @param type Neighbor-graph algorithm: "knn" (\code{\link[dbscan]{kNN}}) or "frnn"
#' (\code{\link[dbscan]{frNN}}, radius-based). Default "frnn".
#' @param used.cts A character vector of cell types to consider, both as centers and as
#' neighbors. If not specified, all available cell types except "outOfROI" are used.
#' @param ct_name Which cell-typing to use (default "default").
#' @param off_target A named list of cell type -> antibody vector, marking off-target/
#' non-expressing antibody-celltype pairs to mask out of the neighbor-averaged expression
#' (see \code{off_target_mode}). Optional; if NULL (default), no masking is applied.
#' @param off_target_mode "na" (default) sets off-target values to NA before averaging;
#' "zero" sets them to 0.
#' @param mc.cores (CycifStack method only) number of cores for \code{parallel::mclapply}
#' across samples. Default 1.
#' @param ... Additional arguments (unused).
#'
#' @return A CellNeighborhood object. Every per-cell slot is computed for EVERY within-ROI
#' cell belonging to `used.cts` -- computeCN() does not subsample. `is_used` records which
#' cells that is (eligibility, not a random sample); if you need a smaller working set for
#' a specific downstream analysis (clustering, plotting), use computeSelection()/
#' subsetCells() afterward -- that is a separate, explicit, reusable step, not something
#' computeCN() decides internally.
#'
#' @details The RCN analysis is performed as follows:
#' 1. The function first identifies the cells that are within the ROIs and belong to `used.cts`.
#' 2. It then computes the neighbor graph via 'dbscan::frNN' (or 'dbscan::kNN').
#' 3. It computes neighbor counts/frequencies/density and neighbor-averaged expression
#' for every such cell.
#' The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
#' If the input is a CycifStack object, the RCN analysis is performed on each Cycif object in the stack.
#'
#' @seealso \code{\link{cyApply}}, \code{\link[dbscan]{frNN}}
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
           ct_name="default",
           off_target = NULL, #list(Tumor_panCK = c("GrB","PD1"), Fibro = "gH2AX")
           off_target_mode = c("na","zero")
           ){ # only for Cycif, and CycifStack
    off_target_mode <- match.arg(off_target_mode)

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
    cts <- cell_types(x,ct_name=ct_name)

    ## coordinates: x => xy
    xy <- xys(x)
    ## NO Y-flip: used only for pairwise-distance-based frNN/kNN below, which is invariant
    ## under a uniform reflection anyway -- but max(Y) here was never a real image height,
    ## just this sample's own cell-Y max (see setWithinROIs() fix, preprocessing.R), so
    ## dropping it for consistency with the topleft-origin convention used elsewhere now.
    # ymax <- max(xy$Y_centroid)
    # xy$Y_centroid <- ymax - xy$Y_centroid
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
    csts <- x@cell_types[[ct_name]]@cell_state_def
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

    # [discontinued] Loop through each antibody and set expression to NA for cell types that do not express it
    if(0){ ## decided not to turn them as NAs as I want to compute expression per CN, irrespective of cell type composition.
      for (ab1 in as.character(unique(csts_dt$ab))) {
        non_expressing_cts <- csts_dt[ab == ab1, cell_types]
        combined_df[.(non_expressing_cts), (ab1) := NA, on = .(cell_types)]
      }
    }

    ## ### NEW: Drop off-target signals BEFORE neighborhood aggregation
    if (!is.null(off_target) && length(off_target) > 0) {
      # Normalize names to character vectors
      # Keep only valid (existing) cell types and protein columns
      valid_cts <- unique(combined_df$cell_types)
      valid_abs <- protein_columns

      # Build a two-column data.table of (cell_types, ab) pairs to zero/NA-out
      ot_pairs <- rbindlist(lapply(names(off_target), function(ct) {
        if (!(ct %in% valid_cts)) return(NULL)
        abs_here <- intersect(as.character(off_target[[ct]]), valid_abs)
        if (length(abs_here) == 0) return(NULL)
        data.table::data.table(cell_types = ct, ab = abs_here)
      }), use.names = TRUE, fill = TRUE)

      if (nrow(ot_pairs) > 0) {
        # Efficiently set each specified antibody to NA (or 0) for rows of that cell type
        # Use data.table::set for speed in a loop over antibodies (column-wise)
        if (off_target_mode == "na") {
          for (ab1 in unique(ot_pairs$ab)) {
            ct_to_null <- ot_pairs[cell_types %in% valid_cts & ab == ab1, unique(cell_types)]
            if (length(ct_to_null) == 0) next
            idx <- combined_df$cell_types %in% ct_to_null
            # set NA
            data.table::set(combined_df, which(idx), ab1, NA_real_)
          }
        } else { # "zero"
          for (ab1 in unique(ot_pairs$ab)) {
            ct_to_zero <- ot_pairs[cell_types %in% valid_cts & ab == ab1, unique(cell_types)]
            if (length(ct_to_zero) == 0) next
            idx <- combined_df$cell_types %in% ct_to_zero
            data.table::set(combined_df, which(idx), ab1, 0.0)
          }
        }
      }
    }

    # tapply(combined_df$PD1,combined_df$cell_types,mean,na.rm=T) # sanity check

    if(type=="knn"){
      nn$eps <- -1
    }else if (type=="frnn"){
      nn$k <- -1
    }

    ## Every row of xy1/cts1/df1 already belongs to used.cts + within ROI (the
    ## `is.used` filter above) -- that IS eligibility. computeCN() computes CN
    ## stats for every one of them; it does not additionally subsample. (Every
    ## cell always includes itself as its own neighbor, so lengths(nn$id) is
    ## always >=1 -- there's no further "has a neighbor" filter needed either.)
    nn1 <- new("NN",
               type = type,
               dist = nn$dist,
               id = nn$id,
               k = nn$k,
               eps = nn$eps,
               sort = nn$sort)

    n.neighbors <- lengths(nn1@id)

    # Create a neighborhood data table from every eligible cell's neighbor list
    nmat <- data.table::data.table(
      cell_id = rep(seq_along(nn1@id), lengths(nn1@id)),
      neighbor_id = unlist(nn1@id, use.names = FALSE)
    )

    # Join with combined_df to get cell type and expression data for each neighbor
    # (off-target ab x celltype combinations were already set to NA/0 above, so
    # both aggregates below inherit that masking)
    neighborhood <- nmat[
      combined_df, on = .(neighbor_id = rn), nomatch = 0
    ]

    # mean_expr_per_ct_nbhd: mean expression among neighbors, per neighbor cell type
    mean_expr_per_ct_nbhd <- neighborhood[, lapply(.SD, mean, na.rm = TRUE), by = .(cell_id, cell_types), .SDcols = protein_columns]
    # mean_expr_per_nbhd: mean expression pooled across all neighbors, any type
    mean_expr_per_nbhd <- neighborhood[, lapply(.SD, mean, na.rm = TRUE), by = .(cell_id), .SDcols = protein_columns]
    data.table::setorder(mean_expr_per_nbhd, cell_id)

    stopifnot(nrow(mean_expr_per_nbhd) == length(nn1@id)) # one row per selected cell, guaranteed by construction

    ## cell type count and frequency in each RCN -- reuses `neighborhood` (the
    ## same cell_id/cell_types pairs used for mean_expr_per_nbhd above), so
    ## rows line up by construction. This is a vectorized grouped count instead
    ## of calling table() once per cell in a loop: table() has real per-call
    ## overhead (factor coercion + hashing), so for a large eligible population
    ## one grouped count + reshape is roughly an order of magnitude faster than
    ## the equivalent sapply(nn1@id, function(ids) table(...)) loop.
    count_dt <- neighborhood[, .(cell_id, cell_types = factor(cell_types, levels = used.cts))][
      , .N, by = .(cell_id, cell_types)]
    ## drop = c(FALSE, FALSE): keep every used.cts column even if some cell type
    ## never actually occurs as a neighbor for any cell in this sample (small
    ## samples especially) -- otherwise dcast silently drops that column and the
    ## dimnames<- assignment below fails since used.cts no longer matches ncol.
    count_wide <- data.table::dcast(count_dt, cell_id ~ cell_types,
                                     value.var = "N", fill = 0, drop = c(FALSE, FALSE))
    data.table::setorder(count_wide, cell_id)
    stopifnot(nrow(count_wide) == length(nn1@id)) # every eligible cell has >=1 neighbor (itself)

    neighbor_counts <- as.matrix(count_wide[, -1, with = FALSE])
    dimnames(neighbor_counts) <- list(seq_len(nrow(neighbor_counts)), used.cts)
    neighbor_freq <- neighbor_counts / rowSums(neighbor_counts)

    if(type=="knn"){
      r <- sapply(nn$dist,function(x)x[k+1])
    }
    neighbor_density <- neighbor_counts/(pi*r^2)

    ## is_used: eligibility (within ROI + used.cts), already full-length
    ## (nCells(x)) and already exactly this population -- no further
    ## composition needed, since computeCN() doesn't subsample.
    is_used <- is.used

    stopifnot(sum(is_used) == length(nn1@id))

    smpls <- rep(smpl, length(nn1@id))

    cn <- new("CellNeighborhood",
              is_used = is_used,
              used.cts = used.cts,
              smpls = smpls,
              neighbor_graph = setNames(list(nn1), smpl),
              n.neighbors = n.neighbors,
              mean_expr_per_ct_nbhd = mean_expr_per_ct_nbhd,
              mean_expr_per_nbhd = mean_expr_per_nbhd,
              neighbor_counts = neighbor_counts,
              neighbor_density = neighbor_density,
              neighbor_freq = neighbor_freq)
    return(cn)
})

#' @rdname computeCN
#' @export
setMethod("computeCN", "CycifStack",
  function(x,r_um=20,k=20,
           type=c("knn","frnn"),
           used.cts,
           ct_name="default",
           off_target=NULL,
           off_target_mode=c("na","zero"),
           mc.cores=1){

    if(missing(type)) type <- "frnn"
    off_target_mode <- match.arg(off_target_mode)

    ## Each sample's computeCN() is fully independent (neighbor search never
    ## crosses sample boundaries), so this is embarrassingly parallel. mclapply
    ## forks per sample when mc.cores>1 (Unix/macOS only, per parallel::mclapply);
    ## on Windows this silently runs sequentially regardless of mc.cores.
    cat("Get neighbors ...\n")
    lst.cn <- parallel::mclapply(names(x),function(nm){
      cat(nm,"\n")
      computeCN(x=x[[nm]],r_um=r_um,k=k,type=type,
        used.cts=used.cts,
        ct_name=ct_name,
        off_target=off_target,
        off_target_mode=off_target_mode)
    }, mc.cores=mc.cores)

    is.err <- sapply(lst.cn,function(cn) inherits(cn,"try-error") || is(cn,"condition"))
    if(any(is.err)){
      stop("computeCN() failed for sample(s): ", paste(names(x)[is.err], collapse=", "),
           "\nFirst error: ", conditionMessage(lst.cn[[which(is.err)[1]]]))
    }
    names(lst.cn) <- names(x)

    cat("Restructure data ...\n")

    ## is_used: each per-sample cn@is_used has length nCells(that sample); the
    ## samples are concatenated in names(x) order, matching how exprs(x)/
    ## cell_types(x) stack cells for a CycifStack -- so c()'ing these gives a
    ## single mask of length nCells(x), directly indexable against exprs(x).
    is_used <- unlist(lapply(lst.cn,function(cn)cn@is_used), use.names=FALSE)
    stopifnot(length(is_used) == sum(nCells(x))) # nCells(CycifStack) is a per-sample vector

    used.cts <- lst.cn[[1]]@used.cts
    smpls <- unlist(lapply(lst.cn,function(cn)cn@smpls), use.names=FALSE)
    n.neighbors <- unlist(lapply(lst.cn,function(cn)cn@n.neighbors), use.names=FALSE)

    ## neighbor relationships are sample-local; keep one NN per sample rather
    ## than attempting to flatten indices into a single cross-sample NN object
    neighbor_graph <- do.call(c,lapply(lst.cn,function(cn)cn@neighbor_graph))

    mean_expr_per_ct_nbhd <- data.table::rbindlist(lapply(lst.cn,function(cn)cn@mean_expr_per_ct_nbhd), fill=TRUE)
    mean_expr_per_nbhd <- data.table::rbindlist(lapply(lst.cn,function(cn)cn@mean_expr_per_nbhd), fill=TRUE)

    neighbor_counts <- as.matrix(data.table::rbindlist(lapply(lst.cn,function(cn)as.data.frame(cn@neighbor_counts)), fill=TRUE))
    neighbor_density <- as.matrix(data.table::rbindlist(lapply(lst.cn,function(cn)as.data.frame(cn@neighbor_density)), fill=TRUE))
    neighbor_freq <- as.matrix(data.table::rbindlist(lapply(lst.cn,function(cn)as.data.frame(cn@neighbor_freq)), fill=TRUE))
    neighbor_counts[is.na(neighbor_counts)] <- 0
    neighbor_density[is.na(neighbor_density)] <- 0
    neighbor_freq[is.na(neighbor_freq)] <- 0

    stopifnot(sum(is_used) == length(smpls),
              sum(is_used) == nrow(mean_expr_per_nbhd),
              sum(is_used) == nrow(neighbor_counts))

    cn <- new("CellNeighborhood",
              is_used=is_used,
              used.cts=used.cts,
              smpls=smpls,
              neighbor_graph=neighbor_graph,
              n.neighbors=n.neighbors,
              mean_expr_per_ct_nbhd=mean_expr_per_ct_nbhd,
              mean_expr_per_nbhd=mean_expr_per_nbhd,
              neighbor_counts=neighbor_counts,
              neighbor_density=neighbor_density,
              neighbor_freq=neighbor_freq)

    x@cell_neighborhood <- cn
    return(x)
  }
)

#_ -------------------------------

# fun: setDist ----

#' set distance to tumorBorder for NN objects
#'
#' @param x A NN object.
#' @param value A numeric vector specifying the distance to tumor border for each cell.
#' @param ... Additional arguments (unused).
#'
#' @export
setGeneric("setDist", function(x,...) standardGeneric("setDist"))

#' @rdname setDist
#' @export
setMethod("setDist", "CellNeighborhood",
  function(x,value){
    x@dist_to_tumor_border <- value # only for Cycif, and CycifStack
    return(x)
  }
)

#_ -------------------------------

# fun: cell_neighborhood accessors ----

#' Get/set the CellNeighborhood object on a Cycif or CycifStack
#'
#' @param x A Cycif or CycifStack object.
#' @param value A CellNeighborhood object (for the setter).
#' @param ... Additional arguments (unused).
#'
#' @rdname cell_neighborhood
#' @export
setGeneric("cell_neighborhood", function(x,...) standardGeneric("cell_neighborhood"))

#' @rdname cell_neighborhood
#' @export
setMethod("cell_neighborhood", "Cycif", function(x) x@cell_neighborhood)

#' @rdname cell_neighborhood
#' @export
setMethod("cell_neighborhood", "CycifStack", function(x) x@cell_neighborhood)

#' @rdname cell_neighborhood
#' @export
setGeneric("cell_neighborhood<-", function(x,value) standardGeneric("cell_neighborhood<-"))

#' @rdname cell_neighborhood
#' @export
setReplaceMethod("cell_neighborhood", "Cycif", function(x,value){
  stopifnot(is(value,"CellNeighborhood"))
  x@cell_neighborhood <- value
  x
})

#' @rdname cell_neighborhood
#' @export
setReplaceMethod("cell_neighborhood", "CycifStack", function(x,value){
  stopifnot(is(value,"CellNeighborhood"))
  x@cell_neighborhood <- value
  x
})

# fun: CellNeighborhood field accessors ----

#' Convenience accessors for CellNeighborhood fields on a Cycif or CycifStack
#'
#' Each of these is a thin wrapper around \code{cell_neighborhood(x)@<slot>}.
#' All are sized \code{sum(is_used(cell_neighborhood(x)))} rows, in the same
#' order -- see \code{\linkS4class{CellNeighborhood}}.
#'
#' @param x A Cycif or CycifStack object with a CellNeighborhood already computed.
#'
#' @rdname cell_neighborhood_fields
#' @export
setGeneric("neighborCounts", function(x) standardGeneric("neighborCounts"))
#' @rdname cell_neighborhood_fields
#' @export
setMethod("neighborCounts", "ANY", function(x) cell_neighborhood(x)@neighbor_counts)

#' @rdname cell_neighborhood_fields
#' @export
setGeneric("neighborFreq", function(x) standardGeneric("neighborFreq"))
#' @rdname cell_neighborhood_fields
#' @export
setMethod("neighborFreq", "ANY", function(x) cell_neighborhood(x)@neighbor_freq)

#' @rdname cell_neighborhood_fields
#' @export
setGeneric("neighborDensity", function(x) standardGeneric("neighborDensity"))
#' @rdname cell_neighborhood_fields
#' @export
setMethod("neighborDensity", "ANY", function(x) cell_neighborhood(x)@neighbor_density)

#' @rdname cell_neighborhood_fields
#' @export
setGeneric("meanExprPerNbhd", function(x) standardGeneric("meanExprPerNbhd"))
#' @rdname cell_neighborhood_fields
#' @export
setMethod("meanExprPerNbhd", "ANY", function(x) cell_neighborhood(x)@mean_expr_per_nbhd)

#' @rdname cell_neighborhood_fields
#' @export
setGeneric("meanExprPerCtNbhd", function(x) standardGeneric("meanExprPerCtNbhd"))
#' @rdname cell_neighborhood_fields
#' @export
setMethod("meanExprPerCtNbhd", "ANY", function(x) cell_neighborhood(x)@mean_expr_per_ct_nbhd)

#' @rdname cell_neighborhood_fields
#' @export
setGeneric("distToTumorBorder", function(x) standardGeneric("distToTumorBorder"))
#' @rdname cell_neighborhood_fields
#' @export
setMethod("distToTumorBorder", "ANY", function(x) cell_neighborhood(x)@dist_to_tumor_border)

#' @rdname cell_neighborhood_fields
#' @export
setGeneric("neighborGraph", function(x) standardGeneric("neighborGraph"))
#' @rdname cell_neighborhood_fields
#' @export
setMethod("neighborGraph", "ANY", function(x) cell_neighborhood(x)@neighbor_graph)

#' Center-cell expression for the cells stored in a CellNeighborhood
#'
#' Not a stored slot -- \code{cell_neighborhood(x)} only ever holds
#' neighbor-averaged quantities. This reads the center cells' own expression
#' straight from \code{exprs(x)}, subset to the same rows and row order as
#' every other CellNeighborhood field (\code{which(is_used)}), so it can be
#' cbind'ed directly against neighborCounts(x)/meanExprPerNbhd(x)/etc.
#'
#' @param x A Cycif or CycifStack object with a CellNeighborhood already computed.
#' @param norm_type Which expression slot to read: "log" or "logTh".
#' @param ... Additional arguments (unused).
#'
#' @rdname centerCellExpr
#' @export
setGeneric("centerCellExpr", function(x,...) standardGeneric("centerCellExpr"))
#' @rdname centerCellExpr
#' @export
setMethod("centerCellExpr", "ANY", function(x,norm_type=c("log","logTh")){
  norm_type <- match.arg(norm_type)
  is_used <- cell_neighborhood(x)@is_used
  exprs(x,type=norm_type)[is_used,,drop=FALSE]
})

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
#' @param ... Additional arguments (unused).
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
  .Defunct(msg=paste(
    "tcnClust() referenced undefined variables (cts.in.rcn) even before the v2",
    "CellNeighborhood slot rename and was not reachable via any working path.",
    "Use computeCN() + your own clustering (e.g. clusterCells()) on",
    "neighborCounts(x)/neighborFreq(x) instead."
  ))
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
#' @param ... Additional arguments (unused).
#'
#' @return An updated 'nn' object with clustering and sorting information.
#'
#' @seealso \code{\link{computeCN}}
#'
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
            .Defunct(msg=paste(
              "rcnClust() referenced undefined variables (cts.in.rcn) and",
              "hardcoded study-specific marker names in inline boxplot()/heatmap3()",
              "calls even before the v2 CellNeighborhood slot rename -- it was not",
              "reachable via any working path. Use computeCN() + your own",
              "clustering (e.g. clusterCells()) on neighborCounts(x)/neighborFreq(x)",
              "instead."
            ))
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
#' @param ... Additional arguments (unused).
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
    .Defunct(msg=paste(
      "meanExpRCN() depended on the rcnClust()/tcnClust() mclustda pathway",
      "(already .Defunct), used $ instead of @ on an S4 nn object, and default-",
      "evaluated cts.in.rcn against a bare data.frame with no cell_types() method",
      "-- it was not reachable via any working path. Use",
      "meanExprPerNbhd(x)/meanExprPerCtNbhd(x) instead."
    ))
  }
)
