#' @rdname CellTypeDefault-slots
#' @export
setGeneric("cell_types", function(x,...) standardGeneric("cell_types"))

#' @export
setMethod("cell_types", "CellTypeCycif", function(x,leaves.only=TRUE,strict=FALSE){
  cts <- x@cell_types
  ctype <- x@cell_lineage_def
  if(leaves.only){
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    leaves <- c(leaves[!grepl("unknown",leaves)],"unknown")
    cts <- factor(cts,levels=leaves)
  }
  if(strict){
    is.strict <- x@is_strict
    cts[!is.strict] <- NA
  }
  return(cts)
})

#' @export
setMethod("cell_types","Cycif",function(x,ctype.full=FALSE,leaves.only=TRUE,strict=FALSE,within.rois=TRUE){
  if(ctype.full){
    cts <- cell_types(x@cell_type_full,leaves.only=leaves.only,strict=strict)
  }else{
    cts <- cell_types(x@cell_type,leaves.only=leaves.only,strict=strict)
  }
  if(within.rois){
    is.rois <- x@within_rois
    cts[!is.rois] <- NA
  }
  return(cts)
})

#' @export
setMethod("cell_types", "CycifStack", function(x,ctype.full=FALSE,leaves.only=TRUE,strict=FALSE,within.rois=TRUE){
  if(ctype.full){
    ctd <- x@cell_type_full
  }else{
    ctd <- x@cell_type
  }
  cts <- ctd@cell_types

  if(leaves.only){
    ctype <- ctd@cell_lineage_def
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    cts <- factor(cts,levels=leaves)
  }

  if(strict){
    is.strict <- ctd@is_strict
    cts[!is.strict] <- NA
  }

  if(within.rois){
    is.rois <- within_rois(x)
    cts[!is.rois] <- NA
  }
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
