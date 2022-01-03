#' Converting a list to a CycifStack object
#'
#' @param x A list or a Cycif object (if Cycif, the length should be one)
#'
#' @return A CycifStack object
#' @export
#'
#' @export
setGeneric("as.CycifStack", function(x) standardGeneric("as.CycifStack"))

setMethod("as.CycifStack", "list",function(x){
  stopifnot(all(sapply(x,class) %in% c("Cycif","CycifStack")))
  n.samples <- sum(sapply(x,length))
  xs <- do.call(c,lapply(x,function(cs){
    if(class(cs)=="CycifStack"){
      out <- cs@samples
    }else if(class(cs)=="Cycif"){
      out <- cs
    }
  }))

  suffix <- xs[[1]]@suffix
  nms <- 	names(xs) <- sapply(xs,names)
  n_cells <- 	sapply(xs,nCells)
  n_cycles <- sapply(xs,nCycles)
  max_cycles <- max(n_cycles)

  idx <- which(n_cycles == max_cycles)[1]
  uniq_abs <- abs_list(xs[[idx]])
  thres <- as.data.frame(sapply(xs,function(y){
    if(.hasSlot(y,"threshold")){
      thres <- y@threshold
    }
  }))
  new("CycifStack",
      n_samples = n.samples,
      names = nms,
      uniq_abs = uniq_abs,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = xs,
      suffix = suffix,
      threshold = thres
  )
})

setMethod("as.CycifStack", "Cycif",function(x){
  as.CycifStack(list(x))
})
