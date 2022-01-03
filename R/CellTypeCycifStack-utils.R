#' Accessing slots of a CellTypeCycifStack object
#'
#' @include CellType-class.R
#' @include CycifStack-utils.R
#'
#' @param x A CellTypeCycifStack object
#' @param value A character vector containing antibodies as in column names of
#' the protein expression matrix.
#'
#' @rdname CellTypeCycifStack-slots
#' @export
setMethod("nCycles", "CellTypeCycifStack", function(x)x@n_cycles)

#' @rdname CellTypeCycifStack-slots
#' @export
setGeneric("uniq_abs<-", function(x,...) standardGeneric("uniq_abs<-"))
setMethod("uniq_abs<-", "CellTypeCycifStack", function(x,value){
  x@uniq_abs <- value
  return(x)
})
