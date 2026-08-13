#_ -------------------------------------------------------

# class: CellTypes ----

#' @title CellTypes Class
#'
#' @description
#' This class represents the cell types and their characteristics derived from CyCIF data.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{cell_lineage_def}}{A data frame containing cell lineage definitions.}
#'   \item{\code{cell_state_def}}{A data frame containing cell state definitions.}
#'   \item{\code{markers}}{A data frame specifying markers associated with cell types.}
#'   \item{\code{n_cycles}}{A numeric value indicating the number of cycles in the CyCIF data.}
#'   \item{\code{sample_names}}{A character vector with sample names associated with the data.}
#'   \item{\code{cell_types}}{A factor representing the generated cell types per individual cells.}
#'   \item{\code{is_strict}}{A logical value indicating strictness in cell type calling.}
#' }
#'
#' @details
#' The \code{CellTypes} class stores information about cell lineage, cell state, markers, and sample-specific details from CyCIF data analysis.
#'
#' @param object A CellTypes object (for the \code{show} method).
#'
#' @seealso
#' Other classes: \code{\link{Cycif}}, \code{\link{CycifStack}}
#'
#' @rdname CellTypes
#' @export
setClass("CellTypes",
         slots = c(
           cell_lineage_def = "data.frame",
           cell_state_def = "data.frame",
           markers = "data.frame",
           n_cycles = "numeric",
           sample_names = "character",
           cell_types = "factor",
           is_strict = "logical"
         )
)

#_ -------------------------------------------------------

# class: LDCoords ----

#' @title LDCoords Class
#'
#' @description This class represents the coordinates and clustering results from dimensionality reduction techniques for CyCIF data.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{ld_type}}{A character vector specifying the dimensionality reduction technique used: 'PCA', 't-SNE', 'UMAP'.}
#'   \item{\code{norm_type}}{A character vector specifying the normalization type: 'log', or 'logTh'}
#'   \item{\code{smpls}}{A character vector containing sample names.}
#'   \item{\code{used.abs}}{A character vector containing antibody names used for the analysis.}
#'   \item{\code{used.cts}}{A character vector containing cell type names used for the analysis.}
#'   \item{\code{n_cells_per_smpl}}{A numeric vector representing the max number of cells per sample selected.}
#'   \item{\code{n_cells_total}}{A numeric value representing the total number of cells used for the analysis.}
#'   \item{\code{ld_coords}}{A data frame containing the coordinates from dimensionality reduction.}
#'   \item{\code{clusters}}{A factor vector containing cluster assignments.}
#'   \item{\code{is_used}}{A logical vector indicating whether each cell is used. The length of the vector is the same as the length of the total cells. The sum of TRUE's is the same as \code{n_cells_total}. }
#'   \item{\code{cts_params}}{A list containing cell type parameters.}
#'   \item{\code{ld_params}}{A list containing dimensionality reduction parameters.}
#'   \item{\code{ld_call}}{A call object representing the dimensionality reduction function call.}
#'   \item{\code{clust_call}}{A call object representing the clustering function call.}
#' }
#'
#' @details
#' The \code{LDCoords} class represents the coordinates and clustering results from dimensionality reduction techniques applied to CyCIF data.
#' It provides information about the samples, used abs, cell types, and various parameters of the cells used in the analysis so the plot can be color-coded accordingly.
#'
#' @rdname LDCoords
#' @export
setClass("LDCoords",
   slots = c(
     ld_type = "character",
     norm_type = "character",
     smpls = "character",
     used.abs = "character",
     used.cts = "character",
     n_cells_per_smpl = "numeric",
     n_cells_total = "numeric",
     ld_coords = "data.frame",
     clusters = "factor",
     is_used = "logical",
     cts_params = "list",
     ld_params = "list",
     ld_call = "call",
     clust_call = "call"
   )
)

#_ -------------------------------------------------------

# class: CellSelection ----

#' @title CellSelection Class
#'
#' @description
#' Records an explicit, reproducible choice of which cells to work with for a
#' specific downstream analysis (e.g. clustering, a UMAP). Not meant to be
#' constructed directly -- built by \code{computeSelection()}, and consumed by
#' \code{subsetCells()}. Centralizes the subsampling logic that used to be
#' duplicated ad hoc in individual analyses (e.g. computeCN()'s old
#' n.sampling, or hand-written per-sample sample() calls).
#'
#' @section Slots:
#' \describe{
#'   \item{\code{used.cts}}{A character vector of eligible cell types.}
#'   \item{\code{n_per_sample}}{A named numeric vector, one entry per sample, giving the
#'   cap applied to that sample (after resolving a scalar/percentile default to per-sample values).}
#'   \item{\code{seed}}{The random seed used.}
#'   \item{\code{is_used}}{A logical vector, length = nCells(the source CycifStack): TRUE
#'   for selected cells. This is the actual mask consumed by \code{subsetCells()}.}
#'   \item{\code{source_samples}}{Character vector of sample names this selection was computed against.}
#'   \item{\code{call}}{The \code{computeSelection()} call that produced this object, for provenance.}
#' }
#'
#' @seealso \code{\link{computeSelection}}, \code{\link{selectionSummary}}
#'
#' @rdname CellSelection
#' @export
setClass("CellSelection",
   slots = c(
     used.cts = "character",
     n_per_sample = "numeric",
     seed = "numeric",
     is_used = "logical",
     source_samples = "character",
     call = "call"
   )
)

#_ -------------------------------------------------------

# class: CellFeatures ----

#' @title CellFeatures Class
#'
#' @description
#' A thin, analysis-ready table: one row per selected cell, columns for
#' center-cell expression, xy coordinates, neighborhood composition/expression
#' (if available), and clinical metadata -- built once by \code{subsetCells()}
#' so downstream steps (UMAP, clustering) work off a flat \code{data.table}
#' rather than re-deriving row alignment against the CycifStack each time.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{data}}{A data.table, one row per selected cell.}
#'   \item{\code{selection}}{The CellSelection this table was built from, for provenance.}
#' }
#'
#' @seealso \code{\link{subsetCells}}, \code{\linkS4class{CellSelection}}
#'
#' @rdname CellFeatures
#' @export
setClass("CellFeatures",
   slots = c(
     data = "data.table",
     selection = "CellSelection"
   )
)
