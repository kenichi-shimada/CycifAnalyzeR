#' Accessing slots in CellType* class
#'
#' @param x A CellTypeDef object
#'
#' @rdname CellTypeDef-slots
#' @export
setGeneric("lineages", function(x) standardGeneric("lineages"))
setMethod("lineages", "CellTypeDef", function(x)x@lineages)

#' @rdname CellTypeDef-slots
#' @export
setGeneric("markers", function(x) standardGeneric("markers"))
setMethod("markers", "CellTypeDef", function(x)x@markers)

#' @rdname CellTypeDef-slots
#' @export
setGeneric("cell_lineages", function(x) standardGeneric("cell_lineages"))
setMethod("cell_lineages", "CellTypeDef", function(x)rownames(x@cell_lineage))

#' @rdname CellTypeDef-slots
#' @export
setGeneric("cell_lineage_df", function(x) standardGeneric("cell_lineage_df"))
setMethod("cell_lineage_df", "CellTypeDef", function(x)x@cell_lineage_df)

#' @rdname CellTypeDef-slots
#' @export
setGeneric("cell_state_df", function(x) standardGeneric("cell_state_df"))
setMethod("cell_state_df", "CellTypeDef", function(x)x@cell_state_df)
