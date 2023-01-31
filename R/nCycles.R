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

#' Set nCycles - note this should run before cell type calling
#'
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

  ## n_cycles, abs_list, raw, log_normalized, logTh_normalized, dna,
  ## dna_thres
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

  x@dna <- x@dna[paste0("DNA",seq(value))]

  if(nrow(x@used_cells)>0){
    x@used_cells <- x@used_cells[,seq(value)]
  }

  if(nrow(x@dna_thres)>0){
    x@dna_thres <- x@dna_thres[paste0("DNA",seq(value)),]
  }

  if(length(x@rois)>0){
    ncys <- sapply(x@rois,function(y)y$cycle)
    x@rois <- x@rois[ncys <= value]
  }
  return(x)
})

setMethod("nCycles<-", "CycifStack", function(x,value){
  if(!is.numeric(value) | length(value) != 1){
    stop("CycifStack@nCycles should be a single numeric value now (not accepting a vector with length of x)")
  }

  ## subsetting samples - trimming ones not reaching the cycle specified
  ncys <- cyApply(x,nCycles,simplify=T)
  is.used.smpls <- ncys >= value
  x <- x[is.used.smpls]

  x <- cyApply(x,function(this){
    # return(class(this))
    smpl <- names(this)
    nCycles(this) <- value
    return(this)
  })

  x@n_cycles <- value
  x@max_cycles <- value # will be discontinued

  x@abs_list <- abs_list(x[[1]]) ## all the samples should have same # of cycles
  this.abs <- as.character(abs_list(x)$ab)

  if(all(this.abs %in% names(x@log_normalized))){
    x@log_normalized <- x@log_normalized[this.abs]
  }

  if(all(this.abs %in% names(x@logTh_normalized))){
    x@logTh_normalized <- x@logTh_normalized[this.abs]
  }

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

