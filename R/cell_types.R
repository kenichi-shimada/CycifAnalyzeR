#' @rdname CellTypeDefault-slots
#' @export
setGeneric("cell_types", function(x,...) standardGeneric("cell_types"))

#' @export
setMethod("cell_types", "CellTypeCycif", function(x,full=TRUE,leaves.only=TRUE){
  if(full){
    cts <- x@cell_types_full
    ctype <- x@expanded_lineage_df
  }else{
    cts <- x@cell_types
    ctype <- x@cell_lineage_df
  }
  if(leaves.only){
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    cts <- factor(cts,levels=leaves)
  }
  return(cts)
})

#' @export
setMethod("cell_types", "Cycif", function(x,full=TRUE,leaves.only=TRUE,within.rois=TRUE){
  cts <- cell_types(x@cell_type,full=full,leaves.only=leaves.only)
  if(within.rois){
    is.rois <- x@within_rois
  }
  cts[!is.rois] <- NA
  return(cts)
})

#' @export
setMethod("cell_types", "CycifStack", function(x,full=TRUE,leaves.only=TRUE,within.rois=TRUE){
  if(full){
    cts <- x@cell_type@cell_types_full
  }else{
    cts <- x@cell_type@cell_types
  }
  if(within.rois){
    is.rois <- within_rois(x)
  }
  cts[is.rois] <- NA
  return(cts)
})

#' #' @rdname CellTypeDefault-slots
#' #' @export
#' setGeneric("markers", function(x) standardGeneric("markers"))
#' setMethod("markers", "CellTypeDefault", function(x)x@markers)
#'
#' #' @rdname CellTypeDefault-slots
#' #' @export
#' setGeneric("cell_lineages", function(x) standardGeneric("cell_lineages"))
#' setMethod("cell_lineages", "CellTypeDefault", function(x)rownames(x@cell_lineage))
#'
#' #' @rdname CellTypeDefault-slots
#' #' @export
#' setGeneric("cell_lineage_df", function(x) standardGeneric("cell_lineage_df"))
#' setMethod("cell_lineage_df", "CellTypeDefault", function(x)x@cell_lineage_df)
#'
#' #' @rdname CellTypeDefault-slots
#' #' @export
#' setGeneric("cell_state_df", function(x) standardGeneric("cell_state_df"))
#' setMethod("cell_state_df", "CellTypeDefault", function(x)x@cell_state_df)
