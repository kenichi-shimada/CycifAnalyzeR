% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{vlnPlot}
\alias{vlnPlot}
\alias{vlnPlot,CycifStack-method}
\title{Violin Plots to Show Protein Expressions}
\usage{
vlnPlot(x, ...)

\S4method{vlnPlot}{CycifStack}(
  x,
  strat.by = c("cell_type", "smpl"),
  ab = "PDL1",
  use.pdata = FALSE,
  fill.var,
  draw_thres = FALSE,
  type = c("raw", "log", "logTh"),
  strict = FALSE,
  ct_name = "default",
  ttl,
  uniq.cts,
  uniq.smpls
)
}
\arguments{
\item{x}{A CycifStack object.}

\item{strat.by}{The strategy for stratifying the violin plots. Choose from "cell_type" or "smpl" (default is "cell_type").}

\item{ab}{The antibody or protein to plot.}

\item{use.pdata}{Should sample metadata be used for additional information? (default is FALSE).}

\item{fill.var}{The variable to use for filling the violin plots (default is "sample").}

\item{draw_thres}{Should the threshold be drawn on the plot? (default is FALSE).}

\item{type}{The type of data to use for plotting. Choose from "raw", "log", or "logTh" (default is "log").}

\item{strict}{Should strict cell type matching be enforced? (default is FALSE).}

\item{ct_name}{The name of the cell type column (default is "default").}

\item{ttl}{The title for the plot (default is determined based on inputs).}

\item{uniq.cts}{Unique cell types to include in the plot (default is all unique cell types).}

\item{uniq.smpls}{Unique samples to include in the plot (default is all samples).}
}
\description{
This function creates violin plots to visualize protein expressions in Cycif data.
}
\details{
\itemize{
\item The \code{vlnPlot} function creates violin plots to visualize the protein expressions in Cycif data.
\item You can stratify the plots by either cell types or samples using the \code{strat.by} parameter.
\item Additional customization of the plot can be achieved using various graphical parameters.
}
}
