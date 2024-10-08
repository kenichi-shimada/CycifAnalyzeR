% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_analysis.R
\name{computeCN}
\alias{computeCN}
\alias{computeCN,Cycif-method}
\alias{computeCN,CycifStack-method}
\title{Compute Recurrent Cell Neighbors (RCN) for Cycif or CycifStack Objects}
\usage{
computeCN(x, ...)

\S4method{computeCN}{Cycif}(
  x,
  r_um = 20,
  k = 20,
  type = c("knn", "frnn"),
  used.cts,
  n.sampling = 1000,
  seed = 123
)

\S4method{computeCN}{CycifStack}(x, r_um = 20, used.cts, n.sampling, seed = 123)
}
\arguments{
\item{x}{A Cycif or CycifStack object.}

\item{n.sampling}{The number of cells to randomly sample for RCN analysis.}

\item{seed}{The random seed for reproducibility.}

\item{r}{The radius within which neighboring cells are considered (in 'unit').}

\item{unit}{The unit of measurement for the radius ('pixel' or 'um'). If 'um' is specified,
the radius 'r' will be converted to pixels based on the assumed resolution (0.65 um per pixel).}

\item{cts.in.center}{A character vector specifying the cell types around which RCN is computed. If not is specified, all available cts are used.}

\item{cts.in.rcn}{A character vector specifying the cell types to consider when computing RCN values. If not is specified, all avaialble cts are used.}
}
\value{
A CellNeighborhood object containing the following components:
\itemize{
\item 'within.rois': A logical vector indicating whether each cell is within a region of interest (ROI). The length is the same as the number of cells in the dataset.
\item 'cts.in.rcn': A character vector specifying the cell types considered when computing cell neighbors, cell type frequency, and expressions.
\item 'n.cells.selected': The number of cells selected for RCN analysis, which is the smaller of \code{n.sampling} and the number of cells within ROIs.
\item 'frnn': A list with recurrent Neighborhood information, including 'dist', 'id', 'eps', and 'sort'. The length of 'dist' and 'id' is the same as 'n.cells.selected'.
\item 'cn_exp': A data frame containing expression data for selected cells.
\item 'is.selected': A logical vector indicating whether each cell is selected for RCN analysis. The sum of the vector is the same as 'n.cells.selected'.
\item 'rcn.count': A data frame containing the counts of neighboring cell types.
\item 'rcn.freq': A data frame containing the relative frequencies of neighboring cell types.
}
}
\description{
This function computes the Recurrent Cell Neighbors (RCN) for Cycif or CycifStack objects.
RCN measures the relative frequencies of neighboring cell types around each cell within specified
radius 'r'. The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
}
\details{
The RCN analysis is performed as follows:
\enumerate{
\item The function first identifies the cells that are within the ROIs.
\item It then computes the Recurrent Neighborhood (frNN) for the selected cells using the 'dbscan::frNN' function.
\item It then computes the RCN values for each cell type based on the relative frequencies of neighboring cell types.
\item The function returns a list containing the RCN values for each cell type.
The RCN analysis can be performed on a single Cycif object or across a CycifStack object.
If the input is a CycifStack object, the RCN analysis is performed on each Cycif object in the stack.
}
}
\seealso{
\code{\link{cyApply}}, \code{\link{dbscan::frNN}}
}
