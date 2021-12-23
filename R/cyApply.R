#' @export
setGeneric("cyApply", function(x,...) standardGeneric("cyApply"))
setMethod("cyApply", "CycifStack", function(x,fun,simplify=FALSE,as.CycifStack=TRUE,...){
  if(simplify){
    out <- sapply(x@samples,fun,...)
  }else{
    out <- lapply(x@samples,fun,...)
    if(all(sapply(out,class)=="Cycif") && as.CycifStack){
      out <- as.CycifStack(out)
    }
  }
  return(out)
})

