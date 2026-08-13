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
  .Defunct(msg=paste(
    "list2Cycif() has had zero call sites anywhere across CycifAnalyzeR or any",
    "downstream project repo (EP, TALAVE, HR-APM) -- it was also broken",
    "(constructed a 'Cycif' object using CycifStack-only slots like n_samples/",
    "samples, and had a dead is(cs,'Cycif')/is(cs,'Cycif') duplicate branch).",
    "Use list2CycifStack() instead."
  ))
})
