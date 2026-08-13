#_ -------------------------------------------------------

# fun: runUMAP ----

#' Run UMAP on a CellFeatures table
#'
#' Runs \code{uwot::umap()} directly on the requested columns of a
#' \code{CellFeatures}' table and writes the two resulting coordinates back as
#' new columns. Unlike the older \code{LdRunUMAP(Cycif/CycifStack, ...)}, this
#' does no subsampling or eligibility filtering itself -- \code{CellFeatures}
#' is already exactly the rows you want to embed (that's \code{computeSelection()}/
#' \code{subsetCells()}'s job), so this function is a thin, single-purpose wrapper.
#'
#' @param x A CellFeatures object.
#' @param used.abs Character vector of column names in \code{x@@data} to embed
#' (e.g. antibody expression columns).
#' @param ld_name Character string used as the column-name prefix for the
#' resulting coordinates (\code{<ld_name>_x}, \code{<ld_name>_y}).
#' @param n_neighbors Passed to \code{uwot::umap()}.
#' @param umap.seed Random seed, for reproducibility. Defaults to 12345, so
#' calling \code{runUMAP()} with no seed argument is still fully reproducible.
#' @param n_threads Passed to \code{uwot::umap()}. Defaults to 1: \code{uwot}'s
#' approximate-nearest-neighbor search runs its own internal RNG when
#' multi-threaded, which \code{set.seed()} does not control -- per uwot's own
#' documentation, \code{n_threads = 1} (together with \code{n_sgd_threads = 1},
#' also defaulted here) is required for run-to-run reproducibility, not just a
#' fixed seed. Pass a higher value only if you don't need reproducibility and
#' want the speedup.
#' @param n_sgd_threads Passed to \code{uwot::umap()}. See \code{n_threads}.
#' @param ... Additional arguments passed to \code{uwot::umap()}.
#'
#' @return \code{x}, with two new columns added to \code{x@@data}.
#'
#' @seealso \code{\link{subsetCells}}, \code{\link{clusterCells}}
#'
#' @rdname runUMAP
#' @export
setGeneric("runUMAP", function(x,...) standardGeneric("runUMAP"))

#' @rdname runUMAP
#' @export
setMethod("runUMAP", "CellFeatures",
  function(x, used.abs, ld_name="umap", n_neighbors=20, umap.seed=12345,
           n_threads=1, n_sgd_threads=1, ...){
    if (missing(used.abs)) {
      stop("'used.abs' should be specified: column names in x@data to embed")
    }
    if (any(!used.abs %in% names(x@data))) {
      stop("Some 'used.abs' are not columns of x@data: ",
           paste(setdiff(used.abs, names(x@data)), collapse=", "))
    }

    mat <- data.matrix(x@data[, used.abs, with=FALSE])
    if (anyNA(mat)) {
      stop("NA values found in x@data[, used.abs] -- clean these before running UMAP")
    }

    set.seed(umap.seed)
    ru <- uwot::umap(mat, n_neighbors=n_neighbors, scale=TRUE,
                      n_threads=n_threads, n_sgd_threads=n_sgd_threads, ...)

    x@data[[paste0(ld_name, "_x")]] <- ru[, 1]
    x@data[[paste0(ld_name, "_y")]] <- ru[, 2]
    x
  }
)

#_ -------------------------------------------------------

# fun: clusterCells ----

#' Cluster the cells in a CellFeatures table
#'
#' Runs the existing \code{LdClustering("matrix", ...)} Seurat-based clustering
#' (\code{FindNeighbors}/\code{FindClusters}) on the requested columns of a
#' \code{CellFeatures}' table, and writes the cluster assignment back as a new
#' column. Deliberately takes the same \code{used.abs} you passed to
#' \code{runUMAP()} (or a different set, e.g. clustering on markers not used
#' for the embedding) rather than remembering it internally -- \code{CellFeatures}
#' is a flat table, not a nested object, so nothing here is implicit.
#'
#' @param x A CellFeatures object.
#' @param used.abs Character vector of column names in \code{x@@data} to cluster on.
#' @param cluster_name Column name to store the resulting cluster assignment under.
#' @param k.param,resolution,algorithm,method,... Passed to \code{LdClustering("matrix", ...)}.
#'
#' @return \code{x}, with one new column added to \code{x@@data}.
#'
#' @seealso \code{\link{runUMAP}}, \code{\link{LdClustering}}
#'
#' @rdname clusterCells
#' @export
setGeneric("clusterCells", function(x,...) standardGeneric("clusterCells"))

#' @rdname clusterCells
#' @export
setMethod("clusterCells", "CellFeatures",
  function(x, used.abs, cluster_name="cluster",
           k.param=20, resolution=0.8, algorithm=1, method="igraph", ...){
    if (missing(used.abs)) {
      stop("'used.abs' should be specified: column names in x@data to cluster on")
    }
    if (any(!used.abs %in% names(x@data))) {
      stop("Some 'used.abs' are not columns of x@data: ",
           paste(setdiff(used.abs, names(x@data)), collapse=", "))
    }

    mat <- data.matrix(x@data[, used.abs, with=FALSE])
    if (anyNA(mat)) {
      stop("NA values found in x@data[, used.abs] -- clean these before clustering")
    }
    rownames(mat) <- as.character(seq_len(nrow(mat))) # Seurat::FindNeighbors requires rownames

    cls <- LdClustering(mat, k.param=k.param, resolution=resolution,
                         algorithm=algorithm, method=method, ...)

    x@data[[cluster_name]] <- cls
    x
  }
)
