#' Show a raw or normalized protein expression matrix
#'
#' @param x A Cycif or CycifStack object.

#' @param type type. If "raw", a raw expression matrix is shown. If "normalized",
#'   a normalized expression matrix will be shown. The method of normalization is
#'   specified by "log" or "LogTh" in the normalize function provided previously.
#' @param silent logical. If FALSE, the method of normalization was shown in
#'  the error prompt.
#'
#' @export
setGeneric("exprs", function(x,...) standardGeneric("exprs"))

#' @export
setMethod("exprs", "Cycif", function(x,type=c("raw","normalized"),silent=TRUE){
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

#' @export
setMethod("exprs", "CycifStack",
          function(x,type="normalized"){
            return(x@normalized)
          }
)
