#_ -------------------------------------------------------
# fun: list2Cycif list ----

#' Converting a list to a Cycif object
#'
#' @param x A list object
#' @return A Cycif object
#'
#' @seealso
#' \code{\link{list2CycifStack}}, \code{\link{Cycif2CycifStack}},
#' \code{\link{Cycif-class}}, \code{\link{CycifStack-class}}
#'
#' @export
setGeneric("list2Cycif", function(x) standardGeneric("list2Cycif"))

#' @rdname list2Cycif
#' @export
setMethod("list2Cycif", "list",function(x){
  slots.in.Cycif <- slotNames(getClassDef("Cycif"))
  stopifnot(all(names(x) %in% slots.in.Cycif))
  n.samples <- sum(sapply(x,length))
  xs <- do.call(c,lapply(x,function(cs){
    if(is(cs,"Cycif")){
      out <- cs@samples
    }else if(is(cs,"Cycif")){
      out <- cs
    }
  }))

  mask_type <- xs[[1]]@mask_type
  nms <- names(xs) <- sapply(xs,names)
  n_cells <- sapply(xs,nCells)
  n_cycles <- sapply(xs,nCycles)
  max_cycles <- max(n_cycles)

  idx <- which(n_cycles == max_cycles)[1]
  abs_list <- abs_list(xs[[idx]])
  new("Cycif",
      n_samples = n.samples,
      names = nms,
      abs_list = abs_list,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = xs,
      mask_type = mask_type
  )
})
