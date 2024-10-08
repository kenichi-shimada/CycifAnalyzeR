% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_analysis.R
\name{computeArea}
\alias{computeArea}
\alias{computeArea,Cycif-method}
\title{Compute the Area of Tumor Regions in CyCIF Data.}
\usage{
computeArea(x, ...)

\S4method{computeArea}{Cycif}(
  x,
  dth,
  unit = c("mm2"),
  plot = TRUE,
  strict = FALSE,
  ct_name = "default",
  fn
)
}
\arguments{
\item{x}{A CyCIF object.}

\item{dth}{Numeric, the distance threshold used for tumor border detection.}

\item{unit}{Character, the unit of the computed area in the output. Default is "mm2" (square millimeters).}

\item{plot}{Logical, whether to plot the tumor regions. Default is TRUE.}

\item{strict}{Logical, whether to use strict cell type filtering. Default is FALSE.}

\item{ct_name}{Character, the name of the cell type for tumor identification. Default is "default".}

\item{fn}{Character, the filename for saving the plot. Ignored if plot is FALSE.}

\item{minPts}{Integer, the minimum number of points required to form a cluster in DBSCAN.}

\item{eps}{Numeric, the maximum distance between two samples for one to be considered as in the neighborhood of the other in DBSCAN.}
}
\value{
A numeric value representing the computed area of tumor regions in the specified unit.
}
\description{
This function calculates the total area of tumor regions within a CyCIF dataset based on the specified distance threshold for tumor border detection.
}
\details{
This function calculates the total area of tumor regions within a CyCIF dataset based on the specified distance threshold for tumor border detection. It uses a combination of DBSCAN clustering to identify tumor cell clusters and concave hull computation to estimate the tumor regions' boundaries. The area is then computed based on the identified tumor regions.
The process involves the following steps:
\enumerate{
\item Cell type filtering: If \code{strict} is set to TRUE, only cells with the specified \code{ct_name} (cell type name) will be considered as tumor cells; otherwise, all non-NA cell types will be considered as tumor cells.
\item DBSCAN clustering: DBSCAN (Density-Based Spatial Clustering of Applications with Noise) is used to cluster the identified tumor cells into groups based on their spatial proximity. The parameters \code{minPts} (minimum number of points required to form a cluster) and \code{eps} (maximum distance between two points to be considered part of the same cluster) can be customized.
\item Cluster merging: Overlapping clusters are merged to create distinct tumor regions.
\item Concave hull computation: For each tumor region, a concave hull is computed using the \code{concaveman} package to approximate its boundary.
\item Area calculation: The area of each tumor region is computed using the \code{sp::Polygon} function, and the areas of all tumor regions are summed to obtain the total area.
The computed area is returned as a numeric value, and the unit of measurement can be specified using the \code{unit} parameter (e.g., "mm2" for square millimeters).
If \code{plot} is set to TRUE, a plot displaying the tumor regions will be generated and saved to the specified \code{fn} (filename).
}
}
\seealso{
\code{\link{defineTumorBorder}} for defining tumor regions, \code{\link{concaveman::concaveman}} for concave hull computation, \code{\link{dbscan::dbscan}} for DBSCAN clustering.
}
