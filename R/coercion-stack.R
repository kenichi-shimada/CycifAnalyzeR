#_ -------------------------------------------------------
# fun: list2CycifStack list ----
# fun: Cycif2CycifStack Cycif ----

#' @title Convert a list of Cycif or CycifStack objects to a CycifStack object
#'
#' @description
#' The `list2CycifStack` and `Cycif2CycifStack` functions convert a list of Cycif or CycifStack objects into a single CycifStack object.
#'
#' @param x A list containing Cycif or CycifStack objects.
#'
#' @return
#' A CycifStack object containing information from the input list of Cycif or CycifStack objects.
#'
#' @details
#' The `list2CycifStack` function takes a list of Cycif or CycifStack objects and combines them into a single CycifStack.
#' It extracts information such as the number of samples, names, antibodies, cycles, and cell counts from each input object,
#' and creates a new CycifStack object in which each element is a Cycif object.
#' The `Cycif2CycifStack` function handles a special case where a single Cycif object is converted to a CycifStack object.
#'
#' @seealso
#' \code{\link{list2Cycif}}, \code{\link{list2CycifStack}}, \code{\link{Cycif2CycifStack}}, \code{\link{Cycif}}, \code{\link{CycifStack}}, \code{\link{CellTypes}}
setGeneric("list2CycifStack", function(x) standardGeneric("list2CycifStack"))

#' @rdname list2CycifStack
#'@export
setMethod("list2CycifStack", "list",function(x){
  stopifnot(all(sapply(x,is,"Cycif")) | all(sapply(x,is,"CycifStack")))
  n.samples <- sum(sapply(x,length))
  xs <- do.call(c,lapply(x,function(cs){
    if(is(cs,"CycifStack")){
      out <- cs@samples
    }else if(is(cs,"Cycif")){
      out <- cs
    }
  }))

  mask_type <- xs[[1]]@mask_type
  nms <- names(xs) <- sapply(xs,names)
  n_cells <- sapply(xs,nCells)
  n_cycles <- sapply(xs,nCycles)
  if(length(unique(n_cycles))==1){
    n_cycles <- unique(n_cycles)
  }
  max_cycles <- max(n_cycles)

  idx <- which(n_cycles == max_cycles)[1]
  al <- abs_list(xs[[idx]])[1:2]
  for(i in seq(xs)){
    al <- al %>% dplyr::left_join(abs_list(xs[[i]]),by=c("ab","cycle"))
  }

  new("CycifStack",
      n_samples = n.samples,
      names = nms,
      abs_list = al,
      n_cycles = n_cycles,
      max_cycles = max_cycles,
      n_cells = n_cells,
      samples = xs,
      mask_type = mask_type
  )
})

#' @rdname list2CycifStack
#' @export
setGeneric("Cycif2CycifStack", function(x) standardGeneric("Cycif2CycifStack"))

#' @rdname list2CycifStack
#' @export
setMethod("Cycif2CycifStack", "Cycif",function(x){
  list2CycifStack(list(x))
})
