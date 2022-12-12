#' @rdname Cycif-slots
#' @export
setGeneric("nCycles", function(x) standardGeneric("nCycles"))
setMethod("nCycles", "Cycif", function(x) x@n_cycles)

#' @rdname CycifStack-slots
#' @export
setMethod("nCycles", "CycifStack", function(x) x@n_cycles)

#' @rdname Cycif-slots
#' @export
setGeneric("nCycles<-", function(x,...) standardGeneric("nCycles<-"))
setMethod("nCycles<-", "Cycif", function(x,value){
  if(!is.numeric(value)){
    stop("value should be numeric")
  }
  max_n <- max(x@abs_list$cycle)
  if(value > max_n){
    stop(paste0("value should be not greater than  max_n = ",max_n))
  }
  x@n_cycles <- value

})
