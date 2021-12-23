#' @export
setGeneric("lineages", function(x) standardGeneric("lineages"))
setMethod("lineages", "CellTypeDef", function(x)x@lineages)

#' @export
setGeneric("markers", function(x) standardGeneric("markers"))
setMethod("markers", "CellTypeDef", function(x)x@markers)

#' @export
setGeneric("cell_lineages", function(x) standardGeneric("cell_lineages"))
setMethod("cell_lineages", "CellTypeDef", function(x)rownames(x@cell_lineage))

#' @export
setGeneric("cell_lineage_df", function(x) standardGeneric("cell_lineage_df"))
setMethod("cell_lineage_df", "CellTypeDef", function(x)x@cell_lineage_df)

#' @export
setGeneric("cell_state_df", function(x) standardGeneric("cell_state_df"))
setMethod("cell_state_df", "CellTypeDef", function(x)x@cell_state_df)
