#' Perform a cell type calling function and set cell types in a Cycif or CycifStack object
#'
#' @rdname defineCellTypes
#' @export
setGeneric("defineCellTypes", function(x,...) standardGeneric("defineCellTypes"))
setMethod("defineCellTypes", "Cycif", function(x,ctype,cstate,gates,p_thres=0.5,...){
  # ctype and cstate is input
  ## run CellTypeCycif
  x@cell_type  <- CellTypeCycif(x,ctype,cstate,gates,ctype.full=FALSE)
  x@cell_type_full  <- CellTypeCycif(x,ctype,cstate,gates,ctype.full=TRUE)

  ## normalize - should be done within CellTypeCycif
  x <- normalize(x,method="logTh",p_thres=p_thres)

  ## set CellTypeCycif object in the cell_type slot
  cts.full <- CellTypeCalling(x,p_thres=p_thres,ctype.full=TRUE)
  cts.short <- CellTypeCalling(x,p_thres=p_thres,ctype.full=FALSE)

  x@cell_type_full@cell_types <- cts.full$cell_type
  x@cell_type_full@is_strict <- cts.full$is_strict
  x@cell_type@cell_types <- cts.short$cell_type
  x@cell_type@is_strict <- cts.short$is_strict

  return(x)
})

setMethod("defineCellTypes", "CycifStack", function(x,ctype,cstate,gates,p_thres=0.5,...){
  # x <- cyApply(x,defineCellTypes,ctype=ctype,cstate=cstate,gates=gates,p_thres=p_thres)
  x@cell_type <- CellTypeCycifStack(x,ctype,cstate,gates)
  x@cell_type@cell_types <- unlist(cyApply(x,cell_types,full=FALSE,leaves.only=TRUE,within.rois=TRUE))
  x@cell_type@cell_types_full <- unlist(cyApply(x,cell_types,full=TRUE,leaves.only=TRUE,within.rois=TRUE))

  return(x)
})
