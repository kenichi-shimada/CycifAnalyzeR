#' Accessing slots in CellType* class
#'
#' @param x A CellTypeDefault object
#'
#' @rdname CellTypeDefault-slots
#' @export
setGeneric("cellTypes", function(x) standardGeneric("cellTypes"))
setMethod("cellTypes", "CellTypeCycif", function(x)x@CellTypeCycif)

#' @rdname CellTypeDefault-slots
#' @export
setGeneric("markers", function(x) standardGeneric("markers"))
setMethod("markers", "CellTypeDefault", function(x)x@markers)

#' @rdname CellTypeDefault-slots
#' @export
setGeneric("cell_lineages", function(x) standardGeneric("cell_lineages"))
setMethod("cell_lineages", "CellTypeDefault", function(x)rownames(x@cell_lineage))

#' @rdname CellTypeDefault-slots
#' @export
setGeneric("cell_lineage_df", function(x) standardGeneric("cell_lineage_df"))
setMethod("cell_lineage_df", "CellTypeDefault", function(x)x@cell_lineage_df)

#' @rdname CellTypeDefault-slots
#' @export
setGeneric("cell_state_df", function(x) standardGeneric("cell_state_df"))
setMethod("cell_state_df", "CellTypeDefault", function(x)x@cell_state_df)
