#' Perform a cell type calling function and set cell types in a Cycif or CycifStack object
#'
#' @rdname findCellTypes
#' @export
setGeneric("findCellTypes", function(x,...) standardGeneric("findCellTypes"))
setMethod("findCellTypes", "Cycif", function(x,ctype,cstate,...){
  # ctype and cstate is input
  ## normalize
    normalize("logTh")
  ## run CellTypeCalling
  ## set CellTypeCycif object into  cell_type slot
  return(x)
})

setMethod("findCellTypes", "CycifStack", function(x,ctype,cstate,...){
  cys <- cyApply(x,findCellTypes,ctype=ctype,cstate=cstate)

  return(cys)
})
