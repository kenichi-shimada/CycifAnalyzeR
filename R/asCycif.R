#' Converting a list to a Cycif object
#'
#' @param x A list object
#'
#' @return A Cycif object
#' @export
#'
#' @export
setGeneric("as.Cycif", function(x) standardGeneric("as.Cycif"))
setMethod("as.Cycif", "list",function(x){
  slots.in.Cycif <- slotNames(getClassDef("Cycif"))
  stopifnot(all(names(x) %in% slots.in.Cycif))
  n.samples <- sum(sapply(x,length))
  xs <- do.call(c,lapply(x,function(cs){
    if(class(cs)=="Cycif"){
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
  new("Cycif",
      n_samples = n.samples,
      names = nms,
      uniq_abs = uniq_abs,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = xs,
      suffix = suffix
  )
})
