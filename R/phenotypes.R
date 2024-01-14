#_ -------------------------------------------------------

# fun: pData CycifStack ----

#' Show or Set Sample metadata
#'
#' These functions allow you to show or set sample metadata for a Cycif or CycifStack object.
#'
#' @param x A Cycif or CycifStack object.
#' @param value For `pData<-`, a data frame containing sample metadata.
#' @param by (For `pData<-`) The column name in the value data frame that contains sample names.
#'        This column name should be set as the 'id' for sample identification.
#'
#' @details
#' - `pData` retrieves and displays the sample metadata associated with the input Cycif or CycifStack object.
#' - `pData<-` allows you to insert or replace sample metadata in the Cycif or CycifStack object with a provided data frame.
#'   Make sure the 'id' column in the data frame matches the sample names accessed by `names(x)`.
#'
#' @rdname pData
#' @export
setGeneric("pData", function(x) standardGeneric("pData"))

#' @rdname pData
#' @export
setMethod("pData", "CycifStack", function(x) x@phenoData)

# fun: pData<- CycifStack ----

#' @rdname pData
#' @export
setGeneric("pData<-", function(x,...,value) standardGeneric("pData<-"))

#' @rdname pData
#' @export
setMethod("pData<-", "CycifStack", function(x,value){
  if(!is(value,"data.frame")){
    stop("Input phenoData should be a data frame\n")
  }
  smpls <- names(x)
  keys <- as.character(value$id)
  if(!all(smpls==keys)){
    stop("'value$id' should be identical to 'names(x)'")
  }
  x@phenoData <- value
  return(x)
})
