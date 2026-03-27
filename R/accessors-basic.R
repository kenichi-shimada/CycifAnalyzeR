#_ -------------------------------------------------------

# names Cycif, CycifStack ----

#' @export
setMethod("names", "Cycif", function(x) x@name)

#' @export
setMethod("names", "CycifStack", function(x){
  ori.names <- x@names
  names(ori.names) <- c()
  return(ori.names)
})

#_ -------------------------------------------------------
# abs_list Cycif, CycifStack ----

#' @title Get the list of antibodies used in a CycifStack object
#'
#' @description
#' The `abs_list` function returns the absolute list of markers in a CycifStack object.
#'
#' @param x A CycifStack object.
#'
#' @return
#' A data frame containing the absolute list of markers with columns "ab" and "cycle."
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
#' @export
setGeneric("abs_list", function(x) standardGeneric("abs_list"))
setMethod("abs_list", "Cycif", function(x) x@abs_list)
setMethod("abs_list", "CycifStack", function(x) x@abs_list)

#_ -------------------------------------------------------

# dna Cycif ----

# create helpfile for the dna function below
#' @title Get the DNA channels in a Cycif object
#'
#' @description
#' The `dna` function returns the DNA channels in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A character vector containing the names of DNA channels in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
#' @examples
#' \dontrun{
#' dna(cycif_object)
#' }
#'
setGeneric("dna", function(x) standardGeneric("dna"))
setMethod("dna", "Cycif", function(x) x@dna)

#_ -------------------------------------------------------

# xys Cycif ----

#' @title Get the XY coordinates in a Cycif object
#'
#' @description
#' The `xys` function returns the XY coordinates in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A data frame containing the XY coordinates in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("xys", function(x) standardGeneric("xys"))
setMethod("xys", "Cycif", function(x) x@xy_coords)

#_ -------------------------------------------------------
# segProp Cycif ----

## roxygen2 help file for segProp()
#'
#' @title Get the segmentation property in a Cycif object
#'
#' @description
#' The `segProp` function returns the segmentation property in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A character vector containing the segmentation property in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("segProp", function(x) standardGeneric("segProp"))
setMethod("segProp", "Cycif", function(x) x@segment_property)

## roxygen2 help file for dna_thres()
#' @title Get the DNA threshold in a Cycif object
#'
#' @description
#' The `dna_thres` function returns the DNA threshold in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A data frame containing the DNA threshold in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("dna_thres", function(x) standardGeneric("dna_thres"))
setMethod("dna_thres", "Cycif", function(x) x@dna_thres)

## roxygen2 help file for used_cells()
#' @title Get the used cells in a Cycif object
#'
#' @description
#' The `used_cells` function returns the used cells in a Cycif object.
#'
#' @param x A Cycif object.
#'
#' @return
#' A data frame containing the used cells in the Cycif object.
#'
#' @export
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#' \code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
#' \code{\link{length}}: Equivalent to `nSamples` for CycifStack objects.
#'
setGeneric("used_cells", function(x) standardGeneric("used_cells"))
setMethod("used_cells", "Cycif", function(x) rowSums(x@used_cells==1)==ncol(x@used_cells))

#' @title Check if samples in a CycifStack are within ROIs
#'
#' @description
#' The `within_rois` function checks if samples in a CycifStack object are within regions of interest (ROIs).
#'
#' @param x A CycifStack object.
#'
#' @return
#' A logical vector indicating whether each sample is within an ROI.
#'
#' @seealso
#' \code{\link{nSamples}}: Get the number of samples in a CycifStack object.
#' \code{\link{names}}: Get the names of samples in a CycifStack object.
#' \code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
#'
#' @export
setGeneric("within_rois", function(x) standardGeneric("within_rois"))
setMethod("within_rois", "Cycif", function(x) x@within_rois)
