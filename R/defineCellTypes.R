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
  CellTypeCalling(x,)
  return(x)
})

setMethod("defineCellTypes", "CycifStack", function(x,ctype,cstate,gates,p_thres=0.5,...){
  cys <- cyApply(x,defineCellTypes,ctype=ctype,cstate=cstate)
  ctcs <- CellTypeCycifStack(x,ctype,cstate,gates)
  cys@cell_type <-ctcs

  return(x)
})
