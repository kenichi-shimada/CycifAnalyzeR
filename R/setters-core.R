#_ -------------------------------------------------------
# nCycles<- Cycif, CycifStack ----

#' @rdname nCycles
#' @export
setGeneric("nCycles<-", function(x,...,value) standardGeneric("nCycles<-"))

#' @rdname nCycles
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
setMethod("nCycles<-", "Cycif", function(x,value){

  if(!is.numeric(value)){
    stop("value should be numeric")
  }
  value <- as.integer(value)

  current.nc <- nCycles(x)
  if(any(value > current.nc)){
    stop("New max n_cycles should be less than or equal to the current max n_cycles.")
  }

  x@n_cycles <- value

  x@abs_list <- abs_list(x) %>% filter(cycle <= value)
  this.abs <- as.character(abs_list(x)$ab)

  x@raw <- x@raw[this.abs]
  if(all(this.abs %in% names(x@log_normalized))){
    x@log_normalized <- x@log_normalized[this.abs]
  }
  if(all(this.abs %in% names(x@logTh_normalized))){
    x@logTh_normalized <- x@logTh_normalized[this.abs]
  }

  x@dna <- x@dna[seq(value)]

  if(nrow(x@used_cells)>0){
    x@used_cells <- x@used_cells[,seq(value)]
  }

  if(nrow(x@dna_thres)>0){
    x@dna_thres <- x@dna_thres[seq(value),]
  }

  if(length(x@rois)>0){
    ncys <- sapply(x@rois,function(y)y$cycle)
    x@rois <- x@rois[ncys <= value]
  }
  return(x)
})

#' @rdname nCycles
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
setMethod("nCycles<-", "CycifStack", function(x,value){
  if(!is.numeric(value) | length(value) != 1){
    stop("CycifStack@n_cycles will only accept a numeric scalar value")
  }

  ncys <- cyApply(x,nCycles,simplify=T)
  is.used.smpls <- ncys >= value

  if(any(!is.used.smpls)){
    stop("Some samples do not reach the specified cycle value. Please check the samples.\n",
         "You can use the 'nCycles' function to check the number of cycles in each sample.\n",
         "Use 'cyApply(x,nCycles,simplify=T)' to get the number of cycles in each sample.")
  }

  x <- cyApply(x,function(this){
    nCycles(this) <- value
    return(this)
  })

  x@n_cycles <- value
  x@max_cycles <- value

  x@abs_list <- abs_list(x[[1]])

  x@names <- smpls <- names(x)
  x@n_samples <- length(x)
  x@samples <- x@samples[smpls]

  if(nrow(x@phenoData)>0){
    x@phenoData <- x@phenoData %>%
      filter(id %in% smpls) %>%
      arrange(match(smpls, id))
  }

  return(x)
})

#_ -------------------------------------------------------
# names<- CycifStack ----

#' @export
setMethod("names<-", "CycifStack", function(x,value){
  ori.names <- x@names
  for(i in seq(ori.names)){
    x@samples[[i]]@name <- as.character(value[i])
  }
  names(x@samples) <- value
  x@names <- value
  return(x)
})

#_ -------------------------------------------------------
# set_abs CycifStack ----

#' @title Check and set antibodies for a Cycif or CycifStack object
#'
#' @description
#' The `check_abs` function checks whether all samples in a CycifStack object have the same number
#' of antibodies and whether they have identical antibody names. If not, it provides
#' suggestions on how to address the issue.
#'
#' The `set_abs` function allows you to set the
#' antibodies for a CycifStack object to a specified set of antibodies.
#'
#' @param x A CycifStack object.
#'
#' @param value A character vector containing the new set of antibodies (for `set_abs<-`).
#'
#' @return
#' For `check_abs`, a warning message is displayed if the samples do not have the same
#' number of antibodies or if the antibody names are not identical.
#'
#' For `set_abs<-`, a modified CycifStack object with the specified antibodies is returned.
#'
#' @export
setGeneric("check_abs", function(x)standardGeneric("check_abs"))

#' @rdname check_abs
#' @export
setMethod("check_abs", "CycifStack", function(x){
  abs <- cyApply(x,function(x){
    ab <- abs_list(x)$ab
  })
  n.abs <- sapply(abs,length)
  if(length(unique(n.abs))>1 | n.abs[1]==0){
    warning("not all the samples have the same # of abs;\nrun 'barplot(cyApply(x,function(cy)nrow(abs_list(cy)),simplify=T),las=2)'\n")
    warning("Run nCycles() to set a fixed # of cycles and/or drop samples that don't meet the criteria.")
  }else if(length(unique(n.abs))>1 | n.abs[1]==0){
    warning("not all the samples have the same Ab names;\nrun 'apply(do.call(cbind,cyApply(x,function(cy)abs_list(cy)$ab)),1,unique)'\n")
    warning("Run set_abs() to set the identical names in the Abs.")
  }
})

#' @rdname check_abs
#' @export
setGeneric("set_abs<-", function(x,...,value)standardGeneric("set_abs<-"))

#' @rdname check_abs
#' @export
setMethod("set_abs<-", "Cycif", function(x,value){
  new.abs <- value
  abs <- abs_list(x)$ab

  if(length(abs)!=length(new.abs)){
    stop("The number of 'new.abs' should be the same as the number of abs in abs_list(x)")
  }else if(!all(names(new.abs)==abs)){
    warning(names(new.abs))
    warning(paste(abs,collapse="\n"))
    stop("The names of the 'new.abs' should be the identical to abs_list(x)$ab")
  }

  x@abs_list$ab <- new.abs
  names(x@raw) <- new.abs

  if(nrow(x@log_normalized)>0){
    names(x@log_normalized) <- new.abs
  }

  if(nrow(x@logTh_normalized)>0){
    names(x@logTh_normalized) <- new.abs
  }

  return(x)
})

#' @rdname check_abs
#' @export
setMethod("set_abs<-", "CycifStack", function(x,value){
  new.abs <- value
  abs <- abs_list(x)$ab
  if(length(abs)!=length(new.abs)){
    stop("The number of 'new.abs' should be the same as the number of abs in abs_list(x)")
  }
  if(!all(names(new.abs)==abs)){
    ("The names of the 'new.abs' should be the identical to abs_list(x)$ab")
  }
  x@abs_list$ab <- new.abs
  x <- cyApply(x,function(cy){
    set_abs(cy) <- new.abs
    return(cy)
  },)
  return(x)
})
