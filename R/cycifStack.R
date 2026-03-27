#_ -------------------------------------------------------

# fun: cyApply CycifStack ----
#' Apply-like loop function for a CycifStack object
#' @param x A CycifStack object
#' @param fun function to apply to each Cycif object
#' @param simplify logical (FALSE by default). If TRUE, the function returns a vector or a matrix rather
#'   than a list, if applicable. It uses 'sapply' instead of 'lapply' internally.
#' @param as.CycifStack logical (TRUE by default). If TRUE, the function attempts to convert the
#'   output into a CycifStack object.
#' @param ... Additional arguments passed to sapply/lapply functions.
#'
#' @usage
#' cyApply(x,fun,simplify=FALSE,as.CycifStack=TRUE,...)
#'
#' @export
setGeneric("cyApply", function(x,...) standardGeneric("cyApply"))
setMethod("cyApply", "CycifStack", function(x,fun,simplify=FALSE,as.CycifStack=TRUE,...){
  if(simplify){
    out <- sapply(x@samples,fun,...)
  }else{
    out <- lapply(x@samples,fun,...)
  }
  if(all(sapply(out,is,"Cycif")) && as.CycifStack){
    out <- list2CycifStack(out)

    if(length(x@cell_types)>0){
      out@cell_types <- x@cell_types
    }
    if(length(x@ld_coords)>0){
      out@ld_coords <- x@ld_coords
    }
    if(nrow(x@phenoData)>0){
      pData(out) <- x@phenoData
    }
    if(length(x@calls)>0){
      out@calls <- x@calls
    }
  }

  return(out)
})

#_ -------------------------------------------------------
