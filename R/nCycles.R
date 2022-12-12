#' @rdname Cycif-slots
#' @export
setGeneric("nCycles", function(x) standardGeneric("nCycles"))
setMethod("nCycles", "Cycif", function(x) x@n_cycles)

#' @rdname CycifStack-slots
#' @export
setMethod("nCycles", "CycifStack", function(x){
  if(!is.null(x@n_cycles)){
    return(x@n_cycles)
  }else{
    return(cyApply(x,nCycles))
  }
})

#' @rdname Cycif-slots
#' @export
setGeneric("nCycles<-", function(x,...) standardGeneric("nCycles<-"))
setMethod("nCycles<-", "Cycif", function(x,value){
  require(dplyr)
  if(!is.numeric(value)){
    stop("value should be numeric")
  }
  value <- as.integer(value)

  current.nc <- nCycles(x)
  if(any(value > current.nc)){
    stop("New max n_cycles should be less than or equal to the current max n_cycles.")
  }

  x@abs_list <- abs_list(x) %>% filter(cycle <= value)
  this.abs <- abs_list(x)$ab

  x@raw <- x@raw[this.abs]
  if(all(this.abs %in% names(x@log_normalized))){
    x@log_normalized <- x@log_normalized[this.abs]
  }
  if(all(this.abs %in% names(x@logTh_normalized))){
    x@logTh_normalized <- x@logTh_normalized[this.abs]
  }

  x@dna <- x@dna[paste0("DNA",seq(value))]

  x@n_cycles <- value

  if(all(this.abs %in% names(x@threshold))){
    x@threshold <- x@threshold[this.abs]
  }

  if(nrow(x@used_cells)>0){
    x@used_cells <- x@used_cells[,seq(value)]
  }

  x@n_cycles <- value
  return(x)
})

setMethod("nCycles<-", "CycifStack", function(x,value){
  if(!is.numeric(value)){
    stop("value should be numeric")
  }
  smpls <- names(x)
  names(smpls) <- NULL
  if(!identical(smpls,names(value))){
    stop("Input should be a vector whose names are identical to sample names. If fewer samples are selected for analysis, first subset samples with [] then use this func.")
  }

  x <- cyApply(x,function(this){
    # return(class(this))
    smpl <- names(this)
    nCycles(this) <- value[smpl]
    return(this)
  })

  x@n_cycles <- value
  x@max_cycles <- max(value)

  if(0){
    ###### "Validate if samples with max_cycles use uniq_abs."
  }else{
    idx <- which(x@n_cycles == x@max_cycles)[1]
  }

  x@uniq_abs <- abs_list(x[[idx]])
  this.abs <- uniq_abs(cs1)$ab

  if(nrow(x@raw)>0){
    x@raw <- x@raw[this.abs]
  }

  if(all(this.abs %in% names(x@log_normalized))){
    x@log_normalized <- x@log_normalized[this.abs]
  }

  if(all(this.abs %in% names(x@logTh_normalized))){
    x@logTh_normalized <- x@logTh_normalized[this.abs]
  }

  return(x)
})

