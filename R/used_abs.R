#' @export
setGeneric("used_abs", function(x) standardGeneric("used_abs"))
setMethod("used_abs", "Cycif", function(x) x@used_abs)
setMethod("used_abs", "CellTypeCycifStack", function(x)x@used_abs)

#' @export
setGeneric("used_abs<-", function(x,...) standardGeneric("used_abs<-"))
setMethod("used_abs<-", "Cycif", function(x,value,strict=TRUE){
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
  x@used_abs <- used.abs
  return(x)
})

#' @export
setMethod("used_abs<-", "CycifStack", function(x,value,strict=TRUE){
  if(length(x)==1){
    xs <- x[[1]]
    used_abs(xs,strict=strict) <- value
    xs <- as.CycifStack(xs)
  }else{
    xs <- cyApply(x,function(y){
      used_abs(y,strict=strict) <- value
      return(y)
    })
  }

  # stopifnot(validObject(xs))
  return(xs)
})

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
