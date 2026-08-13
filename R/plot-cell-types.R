#_ -------------------------------------------------------

# fun: vlnPlot CycifStack ----
#' Violin Plots to Show Protein Expressions
#'
#' This function creates violin plots to visualize protein expressions in Cycif data.
#'
#' @param x A CycifStack object.
#' @param strat.by The strategy for stratifying the violin plots. Choose from "cell_type" or "smpl" (default is "cell_type").
#' @param ab The antibody or protein to plot.
#' @param use.pdata Should sample metadata be used for additional information? (default is FALSE).
#' @param fill.var The variable to use for filling the violin plots (default is "sample").
#' @param draw_thres Should the threshold be drawn on the plot? (default is FALSE).
#' @param type The type of data to use for plotting. Choose from "raw", "log", or "logTh" (default is "log").
#' @param strict Should strict cell type matching be enforced? (default is FALSE).
#' @param ct_name The name of the cell type column (default is "default").
#' @param ttl The title for the plot (default is determined based on inputs).
#' @param uniq.cts Unique cell types to include in the plot (default is all unique cell types).
#' @param uniq.smpls Unique samples to include in the plot (default is all samples).
#' @param ... Additional arguments (unused).
#'
#' @details
#' - The `vlnPlot` function creates violin plots to visualize the protein expressions in Cycif data.
#' - You can stratify the plots by either cell types or samples using the `strat.by` parameter.
#' - Additional customization of the plot can be achieved using various graphical parameters.
#'
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes geom_violin position_dodge ggtitle coord_fixed theme_void
#'
#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @rdname vlnPlot
#' @export
setMethod("vlnPlot", "CycifStack",
          function(x,strat.by=c("cell_type","smpl"), ab="PDL1",
                   use.pdata=FALSE,fill.var,draw_thres=FALSE,
                   type=c("raw","log","logTh"),
                   strict=FALSE,ct_name="default",ttl,
                   uniq.cts,uniq.smpls){
            .Defunct(msg=paste(
              "vlnPlot() has had zero call sites anywhere across CycifAnalyzeR",
              "or any downstream project repo (EP, TALAVE, HR-APM)."
            ))
          })
