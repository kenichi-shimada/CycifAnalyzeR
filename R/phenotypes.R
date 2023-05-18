#_ -------------------------------------------------------

# fun: pData CycifStack ----

#' Show sample metadata
#'
#' @param x A Cycif object
#' @export
setGeneric("pData", function(x) standardGeneric("pData"))
setMethod("pData", "CycifStack", function(x) x@phenoData)

# fun: pData<- CycifStack ----

#' Insert sample metadata into a cycif object
#'
#' @param x A Cycif object
#' @param value input value of a data frame object
#' @param by the column name in the value. The variable should contain the sample names, which can be assessed by `names(x)`. The column name should be set as The variable is set as the 'id'.
#' @export
setGeneric("pData<-", function(x,...,value) standardGeneric("pData<-"))
setMethod("pData<-", "CycifStack", function(x,value){
  if(!is(value,"data.frame")){
    stop("Input phenoData should be a data frame\n")
  }
  smpls <- names(x)
  names(smpls) <- c()
  keys <- as.character(value$id)
  if(!identical(smpls,keys)){
    stop("'value$id' should be identical to 'names(x)'")
  }
  x@phenoData <- value
  return(x)
})
