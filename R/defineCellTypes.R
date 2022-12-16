#' Perform a cell type calling function and set cell types in a Cycif or CycifStack object
#'
#' @rdname defineCellTypes
#' @export
setGeneric("defineCellTypes", function(x,...) standardGeneric("defineCellTypes"))
setMethod("defineCellTypes", "Cycif", function(x,ctype,cstate,gates,...){
  # ctype and cstate is input
  ## run CellTypeCycif
  ctc <- CellTypeCycif(x,ctype,cstate,gates)
  ## normalize - should be done within CellTypeCycif

  ## set CellTypeCycif object in the cell_type slot
  x@cell_type <- ctc
  return(x)
})

setMethod("defineCellTypes", "CycifStack", function(x,ctype,cstate,gates,...){
  cys <- cyApply(x,defineCellTypes,ctype=ctype,cstate=cstate)
  ctcs <- CellTypeCycifStack(x,ctype,cstate,gates)
  cys@cell_type <-ctcs

  return(x)
})
