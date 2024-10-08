% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{xys}
\alias{xys}
\title{Get the XY coordinates in a Cycif object}
\usage{
xys(x)
}
\arguments{
\item{x}{A Cycif object.}
}
\value{
A data frame containing the XY coordinates in the Cycif object.
}
\description{
The \code{xys} function returns the XY coordinates in a Cycif object.
}
\seealso{
\code{\link{nSamples}}: Get the number of samples in a CycifStack object.
\code{\link{names}}: Get the names of samples in a CycifStack object.
\code{\link{maxCycles}}: Get the maximum number of cycles in a CycifStack object.
\code{\link{within_rois}}: Check if samples in a CycifStack are within ROIs.
\code{\link{length}}: Equivalent to \code{nSamples} for CycifStack objects.
}
