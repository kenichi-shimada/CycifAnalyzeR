#_ -------------------------------------------------------

# fun: distAdjacentCells Cycif ----

#' Calculate the Distance Between Adjacent Cells in CyCIF Data.
#'
#' @param x A Cycif or CycifStack object.
#' @param ... Additional arguments (unused).
#'
#' @return A numeric distance threshold (mean cell diameter + 2 SD, in pixels).
#'
#' @export
setGeneric("distAdjacentCells", function(x,...) standardGeneric("distAdjacentCells"))

#' @rdname distAdjacentCells
#' @export
setMethod("distAdjacentCells", "Cycif",
  function(x){
    majl <- x@segment_property$MajorAxisLength
    m1sd <- mean(majl) + sd(majl)
    dth <- m1sd * 2 # 35
    return(dth)
  })

#' @rdname distAdjacentCells
#' @export
setMethod("distAdjacentCells", "CycifStack",
  function(x){
    dths <- cyApply(x,function(cy){
      majl <- cy@segment_property$MajorAxisLength
      m1sd <- mean(majl) + sd(majl)
      dth <- m1sd * 2 # 35
      return(dth)
    },simplify=T)
  dth <- mean(dths) # 34.86 px = 22.65 um
  return(dth)
  })

#_ -------------------------------------------------------
# fun: dist2tumorBorder Cycif ----

#' Calculate the Distance to Tumor Border in CyCIF data.
#'
#' This function calculates the distance of each cell in CyCIF data to the nearest tumor border. The tumor border is defined based on the spatial arrangement of tumor cells within the dataset.
#'
#' @param x A CyCIF object.
#' @param n.cores Integer, number of CPU cores to use for parallel processing.
#' @param dth Numeric, the distance threshold used for tumor border detection.
#' @param th Numeric, distance threshold used when assigning cells to the border/interior split.
#' @param minPts Integer, the minimum number of points required to form a cluster in DBSCAN.
#' @param concavity a relative measure of concavity. 1 results in a relatively detailed shape, Infinity results in a convex hull. You can use values lower than 1, but they can produce pretty crazy shapes (\code{concaveman}).
#' @param plot Logical, whether to plot the tumor border polygons. Default FALSE.
#' @param ... Additional arguments to be passed.
#'
#' @return A list containing the calculated distances and tumor border polygons.
#'
#' @rdname dist2tumorBorder
#' @export
setGeneric("dist2tumorBorder", function(x,...) standardGeneric("dist2tumorBorder"))

#' @rdname dist2tumorBorder
#' @export
setMethod("dist2tumorBorder","Cycif",
          function(x, n.cores = 7, dth, th, minPts = 3, concavity = 0.8, plot = FALSE,...) {
    .Defunct(msg=paste(
      "dist2tumorBorder() has had zero call sites anywhere across CycifAnalyzeR",
      "or any downstream project repo (EP, TALAVE, HR-APM) -- defineTumorStroma()",
      "is what actually computes tumor-border distances in practice. Use",
      "defineTumorStroma() instead."
    ))
  }
)
