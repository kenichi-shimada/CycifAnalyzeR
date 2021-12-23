#' @export
setGeneric("exprs", function(x,...) standardGeneric("exprs"))

#' @export
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

#' @export
setMethod("exprs", "CycifStack",
          function(x,type="normalized"){
            return(x@normalized)
          }
)
