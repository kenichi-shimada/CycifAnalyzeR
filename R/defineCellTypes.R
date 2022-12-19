#' Perform a cell type calling function and set cell types in a Cycif or CycifStack object
#'
#' @rdname defineCellTypes
#' @export
setGeneric("defineCellTypes", function(x,...) standardGeneric("defineCellTypes"))
setMethod("defineCellTypes", "Cycif", function(x,ctype,cstate,gates,p_thres=0.5,...){
  # ctype and cstate is input
  ## run CellTypeCycif
  x@cell_type  <- CellTypeCycif(x,ctype,cstate,gates)
  ## normalize - should be done within CellTypeCycif
  x <- normalize(x,method="logTh",p_thres=p_thres)

  ## set CellTypeCycif object in the cell_type slot
  cts.full <- CellTypeCalling(x,p_thres=p_thres,strict=FALSE,expanded_df=TRUE)
  cts.short <- CellTypeCalling(x,p_thres=p_thres,strict=FALSE,expanded_df=FALSE)
  x@cell_type@cell_types_full <- cts.full
  x@cell_type@cell_types <- cts.short
  return(x)
})

setMethod("defineCellTypes", "CycifStack", function(x,ctype,cstate,gates,p_thres=0.5,...){
  cys <- cyApply(x,defineCellTypes,ctype=ctype,cstate=cstate,gates=gates,p_thres=p_thres)
  ctcs <- CellTypeCycifStack(x,ctype,cstate,gates)
  cys@cell_type <-ctcs

  return(x)
})
