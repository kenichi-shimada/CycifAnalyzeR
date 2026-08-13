#_ -------------------------------------------------------

# fun: LdPlot LDCoords, Cycif, CycifStack ----

#' Create Low-Dimensional Plot
#'
#' This function generates low-dimensional plots for protein expression data from a `Cycif` or `CycifStack` object.
#'
#' @param x A `Cycif` or `CycifStack` object containing protein expression data.
#' @param ld_name The name of the low-dimensional representation (e.g., UMAP or clustering).
#' @param plot_type The type of low-dimensional plot to create, one of "cell_type", "clusters", "exp", or "smpl".
#' @param ab The name of the protein to be used for coloring points (required when `plot_type` is "exp").
#' @param meta.var The metadata column to color by (required when `plot_type` is "meta").
#' @param vals A vector of values to color by (required when `plot_type` is "custom").
#' @param uniq.cols Vector of unique colors to use for plotting points (optional).
#' @param with.labels Logical, indicating whether to label points (default is TRUE).
#' @param used.smpls Vector of samples to include in the plot.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param used.cts Vector of cell types to include in the plot.
#' @param cex.main Size of the main title text (default is 2).
#' @param balanceColor Logical, whether to center the color scale at 0 (default is FALSE).
#' @param main Title for the plot (default is auto-generated based on parameters).
#' @param p_thres Threshold for plotting (default is 0.5).
#' @param mar Margins of the plot.
#' @param cex.labs Size of label text (default is 1).
#' @param cex.leg Size of the legend text (default is 0.5).
#' @param cex Size of data point labels (default is 0.3).
#' @param ... Additional arguments (unused).
#'
#' @details
#' - The `LdPlot` function creates low-dimensional plots for protein expression data based on the specified `Cycif` or `CycifStack` object.
#' - Users can choose from different types of low-dimensional plots, including those based on cell types, clusters, protein expression, or samples.
#' - Various customization options are available for labeling, coloring, and scaling the plot.
#'
#' @rdname LdPlot
#' @export
setGeneric("LdPlot", function(x,...) standardGeneric("LdPlot"))

#' @rdname LdPlot
#' @export
setMethod("LdPlot", "LDCoords",
          function(x,
                   ld_name,
                   plot_type=c("cell_type","clusters","exp","smpl","meta","custom"),
                   ab, # exp
                   meta.var, # meta
                   vals, # custom

                   uniq.cols,with.labels=TRUE,
                   used.smpls,xlab="",ylab="",
                   used.cts,cex.main=2,
                   balanceColor=FALSE,
                   main,p_thres=0.5,mar,cex.labs=1,cex.leg=.5,
                   cex=.3,...){
            .Defunct(msg=paste(
              "LdPlot() has had zero call sites anywhere across CycifAnalyzeR or",
              "any downstream project repo (EP, TALAVE, HR-APM) -- superseded by",
              "plotting directly off a CellFeatures@data table (subsetCells() +",
              "runUMAP()/clusterCells())."
            ))
          }
)

#' @rdname LdPlot
#' @export
setMethod("LdPlot", "Cycif",
  function(x,
           ld_name,
           plot_type=c("cell_type","clusters","exp","smpl","meta","custom"),
           ab, # exp
           meta.var, # meta
           vals, # custom

           uniq.cols,with.labels=TRUE,
           used.smpls,xlab="",ylab="",
           used.cts,cex.main=2,
           main,mar,cex.labs=1,cex.leg=.5,
           cex=.3,...){
    .Defunct(msg=paste(
      "LdPlot() has had zero call sites anywhere across CycifAnalyzeR or any",
      "downstream project repo (EP, TALAVE, HR-APM) -- superseded by plotting",
      "directly off a CellFeatures@data table (subsetCells() +",
      "runUMAP()/clusterCells())."
    ))
  }
)

#' @rdname LdPlot
#' @export
setMethod("LdPlot", "CycifStack",
    function(x,
             ld_name,
             plot_type=c("cell_type","clusters","exp","smpl","meta","custom"),
             ab, # exp
             meta.var, # meta
             vals, # custom

             uniq.cols,with.labels=TRUE,
             used.smpls,xlab="",ylab="",
             used.cts,cex.main=2,
             main,p_thres=0.5,mar,cex.labs=1,cex.leg=.5,
             cex=.3,...){
    .Defunct(msg=paste(
      "LdPlot() has had zero call sites anywhere across CycifAnalyzeR or any",
      "downstream project repo (EP, TALAVE, HR-APM) -- superseded by plotting",
      "directly off a CellFeatures@data table (subsetCells() +",
      "runUMAP()/clusterCells())."
    ))
  }
)
