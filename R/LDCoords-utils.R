#'@export
#' @export
setGeneric("ld_names", function(x)standardGeneric("ld_names"))

#' @export
setMethod("ld_names", "Cycif",function(x)names(x@ld_coords))

#' @export
setMethod("ld_names", "CycifStack",function(x)names(x@ld_coords))