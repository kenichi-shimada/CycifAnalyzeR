#' @rdname CellTypeDefault-slots
#' @export
setGeneric("cell_types", function(x,...) standardGeneric("cell_types"))

#' @export
setMethod("cell_types", "CellTypeCycif", function(x,ctype.full=TRUE,leaves.only=TRUE){
  if(ctype.full){
    cts <- x@cell_types_full
    ctype <- x@expanded_lineage_def
  }else{
    cts <- x@cell_types
    ctype <- x@cell_lineage_def
  }
  if(leaves.only){
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    cts <- factor(cts,levels=leaves)
  }
  return(cts)
})

#' @export
setMethod("cell_types", "Cycif", function(x,ctype.full=TRUE,leaves.only=TRUE,within.rois=TRUE){
  cts <- cell_types(x@cell_type,ctype.full=ctype.full,leaves.only=leaves.only)
  if(within.rois){
    is.rois <- x@within_rois
  }
  cts[!is.rois] <- NA
  return(cts)
})

#' @export
setMethod("cell_types", "CycifStack", function(x,ctype.full=TRUE,leaves.only=TRUE,within.rois=TRUE){
  if(ctype.full){
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
#' setGeneric("cell_lineage_def", function(x) standardGeneric("cell_lineage_def"))
#' setMethod("cell_lineage_def", "CellTypeDefault", function(x)x@cell_lineage_def)
#'
#' #' @rdname CellTypeDefault-slots
#' #' @export
#' setGeneric("cell_state_def", function(x) standardGeneric("cell_state_def"))
#' setMethod("cell_state_def", "CellTypeDefault", function(x)x@cell_state_def)
