#' @export
setMethod("names", "Cycif", function(x) x@name)

#' @export
setGeneric("abs_list", function(x) standardGeneric("abs_list"))
setMethod("abs_list", "Cycif", function(x) x@abs_list)

#' @export
setGeneric("used_abs", function(x) standardGeneric("used_abs"))
setMethod("used_abs", "Cycif", function(x) x@used_abs)

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
setGeneric("abs_params", function(x) standardGeneric("abs_params"))
setMethod("abs_params", "Cycif", function(x) x@abs_params)

#' @export
setGeneric("abs_params<-", function(x,...) standardGeneric("abs_params<-"))
setMethod("abs_params<-", "Cycif", function(x,value){
  stopifnot(!missing(value))
  x@abs_params <- value
  validObject(x)
  return(x)
})

#' @export
setGeneric("exprs", function(x,...) standardGeneric("exprs"))
setMethod("exprs", "Cycif", function(x,type=c("raw","normalized"),na.rm=TRUE,silent=TRUE){
  if(missing(type)){
    # cat("Error: specify type: \"raw\" or \"normalized\"\n")
    # stop()
    type <- "normalized"
  }

  if(type=="normalized" && !("normalized" %in% slotNames(x))){
    stop("Error: slot @normalized not found.\n")
  }

  if(type == "raw"){
    return(x@raw)
  }else if(type == "normalized"){
    if(!silent){
      cat(x@normalize.method,"\n")
    }
    return(x@normalized)
  }
})
