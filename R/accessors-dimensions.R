#_ -------------------------------------------------------
# nCycles Cycif, CycifStack ----

#' @title Get or Set the number of cycles in a CycifStack object
#'
#' @description
#' The `nCycles` function returns the number of cells in each sample of a CycifStack object.
#' The `nCycles<-` function sets the number of cycles in a CycifStack object. It trims samples that do not reach the specified cycle value.
#'
#' @param x A CycifStack object.
#' @param value A numeric scalar indicating the desired number of cycles.
#'
#' @return
#' The `nCycles` function returns a numeric vector containing the number of cells in each sample of the CycifStack object.
#' The `nCycles<-` function returns a modified CycifStack object with the updated number of cycles.
#'
#' @details
#' This `nCycles` function calculates the number of cells in each sample of a CycifStack object and returns the results as a numeric vector.
#' The `nCycles<-` function sets the number of cycles in a CycifStack object to the specified value. It trims samples that do not reach the specified cycle value, ensuring that all samples have the same number of cycles.
#'
#' @note
#' This function is designed to accept a numeric scalar (`value`) as the desired number of cycles. In the future, it will not accept a vector, but only a scalar value.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#' @export
setGeneric("nCycles", function(x) standardGeneric("nCycles"))

#' @rdname nCycles
#' @export
setMethod("nCycles", "Cycif", function(x) x@n_cycles)

#' @rdname nCycles
#' @export
setMethod("nCycles", "CycifStack", function(x){
  n.cycles <- unique(cyApply(x,nCycles,simplify=T))
  if(!is.null(x@n_cycles)){
    if(length(unique(n.cycles))!=1 | x@n_cycles != n.cycles){
      return(n.cycles)
      stop("The number of cycles in the samples are not consistent.")
    }
  }
  return(n.cycles)
})

#_ -------------------------------------------------------
# nSamples CycifStack ----

#' @title Get the number of samples in a CycifStack object
#'
#' @description
#' The `nSamples` function returns the number of samples in a CycifStack object.
#'
#' @param x A CycifStack object.
#'
#' @return
#' The number of samples in the CycifStack object.
#'
#' @seealso
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#'
#' @export
setGeneric("nSamples", function(x)standardGeneric("nSamples"))

#' @rdname nSamples
#' @export
setMethod("nSamples", "CycifStack", function(x) x@n_samples)

#' @rdname nSamples
#' @export
setMethod("length","CycifStack",function(x)x@n_samples)

#_ -------------------------------------------------------
# maxCycles CycifStack ----

#' @title Get the maximum number of cycles used in a CycifStack object
#'
#' @description
#' The `maxCycles` function returns the maximum number of cycles in a CycifStack object.
#'
#' @param x A CycifStack object.
#'
#' @return
#' The maximum number of cycles in the CycifStack object.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#'
#' @export
setGeneric("maxCycles", function(x)standardGeneric("maxCycles"))
setMethod("maxCycles", "CycifStack", function(x) x@max_cycles)

#_ -------------------------------------------------------
# nCells Cycif, CycifStack ----

#' @title Get the number of cells in each sample of a CycifStack object
#'
#' @description
#' The `nCells` function returns the number of cells in each sample of a CycifStack object.
#' Note that the cells are not assessed for their QC status. For the number of cells dropped
#' through DNAFiter or ROI selection, use other methods.
#'
#' @param x A CycifStack object.
#'
#' @return
#' A numeric vector containing the number of cells in each sample of the CycifStack object.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
#' @export
setGeneric("nCells", function(x) standardGeneric("nCells"))

#' @rdname nCells
#' @export
setMethod("nCells", "Cycif", function(x) x@n_cells)

#' @rdname nCells
#' @export
setMethod("nCells", "CycifStack", function(x){
  cyApply(x,nCells,simplify=T)
})
