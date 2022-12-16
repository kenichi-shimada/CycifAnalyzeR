#' Accessing slots in a CellTypeCycifStack object
#'
#' @include Cycif-class.R
#' @include Cycif-utils.R
#' @include CycifStack-class.R
#'
#' @param x A CycifStack object
#'
#' @rdname CycifStack-slots
#' @export
setGeneric("nSamples", function(x)standardGeneric("nSamples"))
setMethod("nSamples", "CycifStack", function(x) x@n_samples)

#' @rdname CycifStack-slots
#' @export
setMethod("names", "CycifStack", function(x) x@names)
setMethod("names<-", "CycifStack", function(x,value){
  ori.names <- x@names
  for(i in seq(ori.names)){
    x@samples[[i]]@name <- as.character(value[i])
  }
  x@names <- value
  # validObject(x)
  return(x)
})

#' @rdname CycifStack-slots
#' @export
setMethod("abs_list", "CycifStack", function(x) x@abs_list)

#' @rdname CycifStack-slots
#' @export
setGeneric("maxCycles", function(x)standardGeneric("maxCycles"))
setMethod("maxCycles", "CycifStack", function(x) x@max_cycles)

#' @rdname CycifStack-slots
#' @export
setMethod("nCells", "CycifStack", function(x) x@n_cells)

#' @rdname CycifStack-slots
#' @export
setMethod("length","CycifStack",function(x)length(x@samples))




