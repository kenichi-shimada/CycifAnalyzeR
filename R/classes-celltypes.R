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
