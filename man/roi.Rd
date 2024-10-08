% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\docType{class}
\name{roi-class}
\alias{roi-class}
\title{ROI (Region of Interest) class}
\description{
This class represents a Region of Interest (ROI) used in the context of CyCIF data.
}
\details{
The \code{roi} class is used to store information about ROIs, including their direction, cycle, type, and associated coordinates.
}
\section{Slots}{

\describe{
\item{\code{direction}}{A character vector specifying the direction of the ROI.
The direction is either "positive" or "negative", which correspond to the ROIs that are included or excluded from the analysis, respectively.}
\item{\code{cycle}}{A numeric vector representing the cycle in the CyCIF used in the analysis.}
\item{\code{roi_type}}{A character vector containing one element, that describes the type of the ROI. Currently accepted roi_types: Polygon, Ellipse, and Rectangle}
\item{\code{coords}}{A data frame containing coordinates associated with the ROI.}
}
}

\examples{
# Creating an roi object
roi_object <- new("roi",
                  direction = "positive",
                  cycle = 5,
                  roi_type = "Polygon",
                  coords = data.frame(x = c(1, 2, 3), y = c(4, 5, 6)))


}
\seealso{
Other classes: \code{\link{Cycif}}, \code{\link{CycifStack}}, \code{\link{CellTypes}}
}
