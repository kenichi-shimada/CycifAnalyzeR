#' Apply-like loop function for a CycifStack object
#' @param x A CycifStack object
#' @param fun function to apply
#' @param simplfy logical (FALSE by default). If TRUE, the function returns a vector or a matrix rather
#'   than a list, if applicable. It uses 'sapply' instead of 'lapply' internally.
#' @param as.CycifStack logical (TRUE by default). If TRUE, the function attempts to convert the
#'   output into a CycifStack object.
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

