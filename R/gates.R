#' Show gates in a Cycif/CycifStack/CellTypeCycif/CellTypeCicyfStack object
#'
#' @include Cycif-class.R
#' @include CycifStack-class.R
#' @include CellType-class.R
#'
#' @rdname gates
#' @export
setGeneric("gates", function(x) standardGeneric("gates"))
setMethod("gates", "Cycif", function(x)x@cell_type@gates)
setMethod("gates", "CycifStack", function(x)x@cell_type@gates)
setMethod("gates", "CellTypeCycif", function(x)x@gates)
setMethod("gates", "CellTypeCycifStack", function(x)x@gates)
