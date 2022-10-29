#' Show or set a threshold in a Cycif or CycifStack object
#'
#' @include Cycif-class.R
#' @include CycifStack-class.R
#' @include CellType-class.R
#'
#' @rdname threshold
#' @export
setGeneric("threshold", function(x) standardGeneric("threshold"))
setMethod("threshold", "Cycif", function(x)x@threshold)
setMethod("threshold", "CycifStack", function(x)x@threshold)
setMethod("threshold", "CellTypeCycifStack", function(x)x@threshold)

#' @rdname threshold
#' @export
setGeneric("threshold<-", function(x,...) standardGeneric("threshold<-"))
setMethod("threshold<-", "Cycif", function(x,value,strict=FALSE,...){
  thres <- value

  smpl <- names(x)

  ## abs
  abs.list <- levels(abs_list(x)$ab)
  used.abs <- names(thres)
  used.abs <- used.abs[used.abs %in% abs.list]

  x@threshold <- thres

  used_abs(x,strict=strict) <- used.abs # x@used_abs
  return(x)
})

#' @rdname threshold
#' @export
setMethod("threshold<-", "CycifStack", function(x,value,strict=FALSE,...){
  thres <- value
  if(class(thres)!="data.frame"){
    stop("'threshold' should be a 'data.frame' object.")
  }

  ## samples
  smpls <- colnames(thres)
  all.smpls <- names(x)
  n.matched <- sum(smpls %in% all.smpls)
  # stopifnot(all(smpls %in% all.smpls))
  cat("all smpls: ",length(all.smpls),", thres: ",length(smpls),", matched:",n.matched,"\n")

  ## abs
  used.abs <- rownames(thres)
  all.abs <- levels(uniq_abs(x)$ab)
  stopifnot(all(used.abs %in% all.abs))

  thres <- as.data.frame(thres)
  x@threshold <- thres

  for(smpl in smpls){
    th <- thres[[smpl]]
    names(th) <- used.abs
    # used.abs1 <- used_abs(x[[smpl]])
    # used.abs1 <- intersect(used.abs,used.abs1)
    # th <- th[used.abs1]
    threshold(x[[smpl]],strict=strict) <- th
    # stop("here")
  }
  return(x)
})

#' @rdname threshold
#' @export
setMethod("threshold<-", "CellTypeCycifStack", function(x,value,...){
  x@threshold <- value
  ## depending on the objec (Cycif vs CycifStack),
  ## value can be vector or data.frame
  return(x)
})

