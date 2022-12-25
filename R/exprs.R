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
setMethod("exprs", "Cycif", function(x,type=c("raw","log_normalized","logTh_normalized")){
  if(missing(type)){
    # cat("Error: specify type: \"raw\" or \"normalized\"\n")
    stop("'type'argument should be specified\n")
  }

  if(type=="log_normalized" && !("log_normalized" %in% slotNames(x))){
    stop("Error: slot @log_normalized not found.\n")
  }

  if(type=="logTh_normalized" && !("logTh_normalized" %in% slotNames(x))){
    stop("Error: slot @logTh_normalized not found.\n")
  }

  if(type == "raw"){
    return(x@raw)
  }else if(type == "log_normalized"){
    return(x@log_normalized)
  }else if(type == "logTh_normalized"){
    return(x@logTh_normalized)
  }
})

#' @export
setMethod("exprs", "CycifStack",function(x,type=c("raw","log_normalized","logTh_normalized")){
  if(missing(type)){
    stop("'type' should be specified\n")
  }

  all.abs <- abs_list(x)$ab
  exps <- do.call(rbind,cyApply(x,function(cy){
    tmp <- exprs(cy,type=type)
    if(any(!all.abs %in% names(tmp))){
      unused <- all.abs[!all.abs %in% names(tmp)]
      tmp1 <- data.frame(array(NA,dim=c(nrow(tmp),length(unused))))
      names(tmp1) <- unused
      tmp <- cbind(tmp,tmp1)[all.abs]
    }
    return(tmp)
  }))
  return(exps)
})
