% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R, R/cycifStack.R
\docType{class}
\name{CycifStack-class}
\alias{CycifStack-class}
\alias{CycifStack}
\alias{show,CycifStack-method}
\title{CycifStack Class}
\usage{
CycifStack(
  ft_filenames,
  path = ".",
  mask_type = c("cellRing", "cell"),
  mcmicro = FALSE,
  use_scimap = FALSE
)

\S4method{show}{CycifStack}(object)
}
\arguments{
\item{ft_filenames}{Character vector of feature table filenames.}

\item{path}{The directory path where feature tables are located. Defaults to the current directory.}

\item{mask_type}{Character vector specifying the mask types to use for sample loading.
Default is c("cellRing", "cell").}

\item{mcmicro}{Logical. Whether the input files are from the MCMicro platform. Defaults to FALSE.}

\item{use_scimap}{Logical. Whether to use scimap for additional processing. Defaults to FALSE.}

\item{object}{a CycifStack object}
}
\value{
A CycifStack object containing information from the loaded Cycif or CycifStack samples.
}
\description{
This class represents a collection of CyCIF samples as a stack, allowing manipulation and analysis of multiple samples together.

The \code{CycifStack} function loads Cycif or CycifStack samples from provided feature table filenames
and creates a CycifStack object containing information from these samples.
}
\details{
The \code{CycifStack} class represents a collection of CyCIF samples organized as a stack. It allows users to analyze and manipulate multiple samples simultaneously, making it useful for high-throughput analysis of CyCIF data.
}
\section{Slots}{

\describe{
\item{\code{samples}}{A list of \code{Cycif} objects representing individual CyCIF samples.}
\item{\code{names}}{A character vector containing names for each sample.}
\item{\code{mask_type}}{A character vector specifying the type of mask used, e.g., 'cellRing'.}
\item{\code{n_samples}}{A numeric value indicating the total number of CyCIF samples in the stack.}
\item{\code{n_cycles}}{A numeric value indicating the maximum number of imaging cycles among all samples in the stack.}
\item{\code{max_cycles}}{A numeric value indicating the maximum number of imaging cycles among all samples in the stack.}
\item{\code{n_cells}}{A numeric value representing the total number of cells in the stack (including dropped cells).}
\item{\code{abs_list}}{A data frame containing information about antibodies and the cycles they were used imaging cycles.}
\item{\code{cell_types}}{A list of \code{CellType} objects representing cell types in the stack.}
\item{\code{ld_coords}}{A list of \code{LDCoords} objects.}
\item{\code{cell_neighborhood}}{An object of class \code{CellNeighborhood} representing cell neighborhood analysis results.}
\item{\code{phenoData}}{A data frame containing phenotypic information.}
\item{\code{calls}}{A list of functions called.}
}
}

\examples{
# Example usage:
ft_filenames <- c("unmicst-sample1_cellRing.csv", "unmicst-sample2_cell.csv")
stack <- CycifStack(ft_filenames)

}
\seealso{
Other classes: \code{\link{Cycif}}, \code{\link{CellTypes}}, \code{\link{LDCoords}}
}
