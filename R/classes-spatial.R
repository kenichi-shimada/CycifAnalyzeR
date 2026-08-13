#_ -------------------------------------------------------

# class: NN ----

#' Class "NN" - output of dbscan::frNN() or dbscan::kNN()
#'
#' @slot type A character vector specifying the type of nearest neighbors. Possible values are 'frnn' or 'knn'.
#' @slot id A list of integer vector (length of the number of cell neighborhoods). Each vector contains the ids of the fixed radius nearest neighbors.
#' @slot dist A list with distances (same structure as id)
#' @slot k A numeric vector (scalar) containing the number of neighbors k that was used.
#' @slot eps A numeric vector (scalar) containing neighborhood radius eps that was used.
#' @slot sort A logical value indicating whether the distances are sorted.
#'
#' @seealso \code{\link[dbscan]{frNN}} \code{\link[dbscan]{kNN}}
#'
#' @rdname NN
#'
#' @export
setClass("NN",
         slots = c(
           type = "character",
           dist = "list",
           id = "list",
           k = "numeric",
           eps = "numeric",
           sort = "logical"
         )
)

# class: CellNeighborhood ----
#' Class "CellNeighborhood" - Cell Neighborhood Object
#'
#' An object of class "CellNeighborhood" represents the results of the computeCN function
#' for a Cycif or CycifStack object, providing information about the cell neighborhood analysis.
#'
#' Every per-cell slot below is stored *compact*: only the cells for which
#' \code{is_used} is \code{TRUE} have a row. \code{is_used} itself is the one
#' full-length mask (length = \code{nCells(x)}, matching the parent Cycif/
#' CycifStack) that records which original cells those rows correspond to --
#' \code{which(is_used)} gives their original row indices. No slot is ever
#' padded with NA rows for unused cells; this is the same pattern
#' \code{LdRunUMAP()}/\code{LDCoords} already use for \code{is_used}/coordinates.
#'
#' @slot is_used A logical vector, length = nCells(the source object). TRUE for
#'   cells that were selected (within ROI, of a used cell type, has >=1
#'   neighbor, and survived the n.sampling cap). All other slots below have
#'   \code{sum(is_used)} rows, in the same order as \code{which(is_used)}.
#' @slot used.cts A character vector of cell types considered in the neighborhood analysis.
#' @slot smpls A character vector, length sum(is_used): which sample each stored row belongs to.
#' @slot neighbor_graph A named list of \code{NN} objects, one per sample (neighbor
#'   relationships are always sample-local, so there is no single cross-sample NN).
#' @slot n.neighbors A numeric vector, length sum(is_used): number of neighbors (including self).
#' @slot mean_expr_per_ct_nbhd A data.table, sum(is_used) rows: mean marker expression
#'   among each cell's neighbors, grouped by neighbor cell type. Off-target
#'   ab x celltype combinations (per the \code{off_target} argument to
#'   \code{computeCN}) are NA, not silently averaged in.
#' @slot mean_expr_per_nbhd A data.table, sum(is_used) rows: mean marker expression
#'   pooled across all neighbors regardless of type. Uses the same off-target
#'   masking as mean_expr_per_ct_nbhd.
#' @slot neighbor_counts A matrix, sum(is_used) rows: counts of each neighboring cell type.
#' @slot neighbor_density A matrix, sum(is_used) rows: neighbor_counts normalized by search area.
#' @slot neighbor_freq A matrix, sum(is_used) rows: relative frequency of each neighboring cell type.
#' @slot dist_to_tumor_border A numeric vector, sum(is_used) rows: distance of each cell to the tumor border.
#' @slot mclustda A list; legacy slot used by the (currently unsupported) rcnClust()/tcnClust() path.
#'
#' @seealso \code{\link{computeCN}}
#'
#' @rdname CellNeighborhood
#' @export
setClass("CellNeighborhood",
         slots = c(
           is_used = "logical",
           used.cts = "character",
           smpls = "character",
           neighbor_graph = "list",
           n.neighbors = "numeric",
           mean_expr_per_ct_nbhd = "data.table",
           mean_expr_per_nbhd = "data.table",
           neighbor_counts = "matrix",
           neighbor_density = "matrix",
           neighbor_freq = "matrix",
           dist_to_tumor_border = "numeric",
           mclustda = "list"
         ))
