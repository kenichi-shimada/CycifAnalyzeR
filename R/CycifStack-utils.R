#' @export
setGeneric("nSamples", function(x)standardGeneric("nSamples"))
setMethod("nSamples", "CycifStack", function(x) x@n_samples)

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

#' @export
setGeneric("uniq_abs", function(x)standardGeneric("uniq_abs"))
setMethod("uniq_abs", "CycifStack", function(x) x@uniq_abs)

#' @export
setMethod("nCycles", "CycifStack", function(x) x@n_cycles)

#' @export
setGeneric("maxCycles", function(x)standardGeneric("maxCycles"))
setMethod("maxCycles", "CycifStack", function(x) x@max_cycles)

#' @export
setMethod("nCells", "CycifStack", function(x) x@n_cells)

#' @export
setMethod("length","CycifStack",function(x)length(x@samples))

