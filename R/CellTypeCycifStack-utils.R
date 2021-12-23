#' @include CellType-class.R
#' @include CycifStack-utils.R
#' @export
setMethod("nCycles", "CellTypeCycifStack", function(x)x@n_cycles)

#' @export
setGeneric("uniq_abs<-", function(x,...) standardGeneric("uniq_abs<-"))
setMethod("uniq_abs<-", "CellTypeCycifStack", function(x,value){
  x@uniq_abs <- value
  return(x)
})
