#' Perform a cell type calling function and set cell types in a Cycif or CycifStack object
#'
#' @rdname findCellTypes
#' @export
setGeneric("findCellTypes", function(x,...) standardGeneric("findCellTypes"))
setMethod("findCellTypes", "Cycif", function(x,ctype,cstate,gates,...){
  # ctype and cstate is input
  ## run CellTypeCycif
  ctc <- CellTypeCycif(x[[1]],ctype,cstate,gates)
  ## normalize - should be done within CellTypeCycif


  ## set CellTypeCycif object in the cell_type slot
  x@cell_type <- ctc
  return(x)
})

setMethod("findCellTypes", "CycifStack", function(x,ctype,cstate,gates,...){
  cys <- cyApply(x,findCellTypes,ctype=ctype,cstate=cstate)
  ctcs <- CellTypeCycifStack(x,ctype,cstate,gates)
  cys@cell_type <-ctcs

  return(x)
})
