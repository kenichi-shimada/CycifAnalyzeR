#' Show sample metadata
#'
#' @param x A Cycif object
#' @export
setGeneric("pData", function(x) standardGeneric("pData"))
setMethod("pData", "CycifStack", function(x) x@phenoData)

#' Insert sample metadata into a cycif object
#'
#' @param x A Cycif object
#' @param value input value of a data frame object
#' @param by the column name in the value. The variable should contain the sample names, which can be assessed by `names(x)`. The column name should be set as The variable is set as the 'id'.
#' @export
setGeneric("pData<-", function(x,...) standardGeneric("pData<-"))
setMethod("pData<-", "CycifStack", function(x,value,by){
  if(class(value)!="data.frame"){
    stop("Input phenoData should be a data frame\n")
  }
  if(missing(by)){
    stop("A variable name should be specified ('by')")
  }
  if(!any(names(value)==by)){
    stop("A variable should exist in the input data")
  }
  smpls <- names(x)
  df <- data.frame(id=smpls)
  value <- value %>%
    dplyr::mutate(id=!!sym(by)) %>%
    dplyr::relocate(id)
  df <- df %>% left_join(value,by="id")

  x@phenoData <- df
  return(x)
})
