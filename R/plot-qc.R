#_ -------------------------------------------------------

# fun: plotUsedCellRatio Cycif,CycifStack ----

#' Plot Used Cell Ratio
#'
#' This function creates a plot of the ratio of used cells across channels in Cycif data.
#'
#' @param x A CycifStack object or a list of Cycif objects.
#' @param cumulative Should the cumulative ratio be plotted? (default is TRUE)
#' @param ncycle The number of cycles to consider for plotting (default is determined by the data).
#' @param mar Margins for the plot (default is c(5, 5, 4, 10)).
#' @param main Main title for the plot (default is "# cells attached on slide").
#' @param smpl.cols Colors for sample labels (default is generated using a color palette).
#' @param ncol Number of columns for the legend (default is 1).
#' @param leg.cex Size of the legend text (default is 0.8).
#' @param use_rois Should regions of interest (ROIs) be considered in the analysis? (default is TRUE)
#' @param ... Additional graphical parameters to customize the plot.
#'
#' @details
#' - The `plotUsedCellRatio` function creates a plot that visualizes the ratio of used cells in Cycif data.
#' - It calculates and plots the ratio of cells that are considered "used" based on the presence of regions of interest (ROIs).
#' - You can choose whether to display a cumulative ratio and specify the number of cycles to consider.
#' - The function also allows customization of various graphical parameters for the plot.
#'
#' @rdname plotUsedCellRatio
#' @export
setGeneric("plotUsedCellRatio",function(x,...) standardGeneric("plotUsedCellRatio"))

#' @rdname plotUsedCellRatio
#' @export
setMethod("plotUsedCellRatio", "Cycif", function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                                                 main="# cells attached on slide",...){
  .Defunct(msg=paste(
    "plotUsedCellRatio() plots statUsedCells()'s per-cycle retention ratio, which is retired --",
    "see ?statUsedCells. Under @used_cells set uniformly across cycles (e.g. ROI membership),",
    "this would just be a flat line. Summarize x@used_cells directly per sample instead."
  ))
})

#' @rdname plotUsedCellRatio
#' @export
setMethod("plotUsedCellRatio", "CycifStack",
          function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                   main="# cells attached on slide",smpl.cols,ncol=1,
                   leg.cex=0.8,use_rois=TRUE,...){
  .Defunct(msg=paste(
    "plotUsedCellRatio() plots statUsedCells()'s per-cycle retention ratio, which is retired --",
    "see ?statUsedCells. Under @used_cells set uniformly across cycles (e.g. ROI membership),",
    "this would just be a flat line. Summarize x@used_cells directly per sample instead."
  ))
})

#_ -------------------------------------------------------

# fun: plotAvailCellOnSlide Cycif ----
#' Plot Available Cells on Slide
#'
#' This function creates a plot that visualizes the available cells on a slide in Cycif data.
#'
#' @param x A Cycif object.
#' @param upside.down Should the plot be upside-down? (default is TRUE)
#' @param ncycle The number of cycles to consider for plotting (default is determined by the data).
#' @param mfrow Number of rows and columns in the grid of plots (default is c(3, 3)).
#' @param mar Margins for each individual plot (default is c(0, 0, 4, 0)).
#' @param legend Should a legend be included in the plot? (default is TRUE)
#' @param main Main title for the plot (default is the names of the Cycif object).
#' @param cex.title Size of the main title text (default is 1).
#' @param uniq.cols Colors for different cell states (default colors).
#' @param legend.cex Size of the legend text (default is 2).
#' @param xlab Label for the x-axis (default is empty).
#' @param ylab Label for the y-axis (default is empty).
#' @param ... Additional graphical parameters to customize the individual plots.
#'
#' @details
#' - The `plotAvailCellOnSlide` function creates a grid of plots, each representing the available cells on a slide for a specific cycle.
#' - The plots show the spatial distribution of available cells in different colors, with an optional legend.
#' - You can customize the appearance of the plots and the legend using various graphical parameters.
#'
#' @export
setGeneric("plotAvailCellOnSlide",function(x,...) standardGeneric("plotAvailCellOnSlide"))

#' @rdname plotAvailCellOnSlide
#' @export
setMethod("plotAvailCellOnSlide", "Cycif",
          function(x,upside.down=TRUE,ncycle,mfrow=c(3,3),mar=c(0,0,4,0),legend=TRUE,main=names(x),cex.title=1,
                   uniq.cols=c(lost="grey90",dropped="blue",available="black",bunched="red"),legend.cex=2,
                   xlab="",ylab="",...){
            .Defunct(msg=paste(
              "plotAvailCellOnSlide() draws one spatial panel per cycle (lost/dropped/available/",
              "bunched), and labels each panel using statUsedCells(), which is retired -- see",
              "?statUsedCells. Under @used_cells set uniformly across cycles (e.g. ROI membership),",
              "every panel would be identical. Plot xys(x) colored by x@used_cells[,1] directly",
              "for a single ROI-membership map instead."
            ))
          })

#_ -------------------------------------------------------


#_ -------------------------------------------------------

# fun: hist_1d numeric ----

#' 1D Histogram Plot
#'
#' This function creates a 1D histogram plot for a numeric vector.
#'
#' @param x A numeric vector for which the histogram is to be plotted.
#' @param n The number of bins or breaks for the histogram (default is 1000).
#' @param ths A vector of threshold values to mark on the plot (default is NULL).
#' @param mar A numerical vector of length 4 specifying the margin sizes (default is c(3, 4, 4, 2) + 0.1).
#' @param brks1 A vector of pre-defined breaks for the histogram (default is NULL).
#' @param ttl1 The title for the histogram plot (default is NULL).
#'
#' @details
#' - The `hist_1d` function creates a 1D histogram plot to visualize the distribution of a numeric vector.
#' - You can specify the number of bins with the `n` parameter or provide custom break points with the `brks1` parameter.
#' - Threshold values can be added to the plot using the `ths` parameter.
#' - Additional customization of the plot can be achieved using various graphical parameters.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats loess predict
#' @export
hist_1d <- function(x,n=1000,ths,mar=c(3,4,4,2)+.1,brks1,ttl1){
  omar <- par()$mar
  par(mar=mar)

  n.ab <- trim_fun(x,trim_th=1e-2)
  min.i <- which.min(abs(brks1-min(n.ab)))
  max.i <- which.min(abs(brks1-max(n.ab)))
  # stop(list(min.i,max.i))
  uniq.cols <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(max.i-min.i+1)
  # stop(length(uniq.cols))
  cols <- c(rep(uniq.cols[1],min.i-1),uniq.cols,rep(rev(uniq.cols)[1],n-max.i+1))
  a <- hist(x,breaks=brks1,main=ttl1,freq=FALSE,xlab="",col=cols,border=NA)

  ## smoothening the trail of histogram
  loessMod <- loess(a$density[seq(n)] ~ brks1[seq(n)], span=0.02)
  smoothened <- predict(loessMod)
  lines(smoothened, x=brks1[seq(n)], col=1,lwd=2)

  # cat("Showing current dna_thres\n")
  abline(v=ths[1],col=4,lty=2,lwd=2)
  abline(v=ths[2],col=2,lty=2,lwd=2)
  par(mar=omar)
  invisible(list(a=a,smoothened=smoothened))
}
