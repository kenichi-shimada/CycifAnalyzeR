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
#' @seealso \code{\link{dbscan::frNN}} \code{\link{dbscan::kNN}}
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
#' @slot within.rois A logical vector indicating whether each cell is within a region of interest (ROI).
#' @slot used.cts A character vector of cell types considered in the neighborhood analysis.
#' @slot n.cells.selected An integer indicating the number of cells selected for the analysis.
#' @slot smpls A character vector containing sample names.
#' @slot nn A frNN object containing information about cell neighbors.
#' @slot n.neighbors An integer indicating the number of neighbors (including self).
#' @slot exp.per.ct.cn A data.table containing expression data for each cell type in the neighborhood.
#' @slot exp.per.cn A data.table containing expression data for each cell neighborhood.
#' @slot is.selected A logical vector indicating whether each cell is selected.
#' @slot rcn.count A data.frame containing the counts of cell types in the neighborhood.
#' @slot rcn.dens A data.frame containing the density of cell types in the neighborhood.
#' @slot rcn.freq A data.frame containing the frequency of cell types in the neighborhood.
#' @slot dist2tumorBorder A numeric vector containing the distance of each cell to tumor border.
#' @slot mclustda A list containing the results of the mclustDA function.
#'
#' @seealso \code{\link{computeCN}}
#'
#' @rdname CellNeighborhood
#' @export
setClass("CellNeighborhood",
         slots = c(
           within.rois = "logical",
           n.cells.selected = "numeric",
           is.selected = "logical",
           used.cts = "character",
           smpls = "character",
           nn = "NN",
           n.neighbors = "numeric",
           exp.per.ct.cn = "data.table",
           exp.per.cn = "data.table",
           rcn.count = "matrix",
           rcn.dens = "matrix",
           rcn.freq = "matrix",
           dist2tumorBorder = "numeric",
           mclustda = "list"
         ))
