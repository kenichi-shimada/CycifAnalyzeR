#' @export
setMethod("nCycles", "CellTypeCycifStack", function(x)x@n_cycles)

#' @export
setMethod("threshold", "CellTypeCycifStack", function(x)x@threshold)
setMethod("threshold<-", "CellTypeCycifStack", function(x,value,...){
  x@threshold <- value
  ## depending on the objec (Cycif vs CycifStack),
  ## value can be vector or data.frame
  return(x)
})

#' @export
setMethod("used_abs", "CellTypeCycifStack", function(x)x@used_abs)
setMethod("used_abs<-", "CellTypeCycifStack", function(x,value,...){
  all.abs <- levels(abs_list(x)$ab)
  if(!all(value %in% all.abs)){
    if(strict){
      unused.value <- value[!value %in% all.abs]
      stop("following channels not found in @abs_list:\n",
           paste(unused.value,collapse=","))
    }
    ## non-strict is useful, e.g. when the commonly used antibodies are
    ## specified but one particular sample was not stained with some antibodies
    ## but you want to keep the samples for the analysis anyway.
  }
  used.abs <- value[value %in% all.abs]
  cat(names(x),": ",length(used.abs)," out of ",length(all.abs)," channels are used.\n")

  x@used_abs <- value
  return(x)
})

#' @export
setGeneric("uniq_abs<-", function(x,...) standardGeneric("uniq_abs<-"))
setMethod("uniq_abs<-", "CellTypeCycifStack", function(x,value){
  x@uniq_abs <- value
  return(x)
})
