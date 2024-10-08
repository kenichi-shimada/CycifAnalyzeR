% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R, R/dimensionality_reduction.R
\docType{class}
\name{LDCoords-class}
\alias{LDCoords-class}
\alias{show,LDCoords-method}
\title{LDCoords Class}
\usage{
\S4method{show}{LDCoords}(object)
}
\description{
This class represents the coordinates and clustering results from dimensionality reduction techniques for CyCIF data.
}
\details{
The \code{LDCoords} class represents the coordinates and clustering results from dimensionality reduction techniques applied to CyCIF data.
It provides information about the samples, used abs, cell types, and various parameters of the cells used in the analysis so the plot can be color-coded accordingly.
}
\section{Slots}{

\describe{
\item{\code{ld_type}}{A character vector specifying the dimensionality reduction technique used: 'PCA', 't-SNE', 'UMAP'.}
\item{\code{norm_type}}{A character vector specifying the normalization type: 'log', or 'logTh'}
\item{\code{smpls}}{A character vector containing sample names.}
\item{\code{used.abs}}{A character vector containing antibody names used for the analysis.}
\item{\code{used.cts}}{A character vector containing cell type names used for the analysis.}
\item{\code{n_cells_per_smpl}}{A numeric vector representing the max number of cells per sample selected.}
\item{\code{n_cells_total}}{A numeric value representing the total number of cells used for the analysis.}
\item{\code{ld_coords}}{A data frame containing the coordinates from dimensionality reduction.}
\item{\code{clusters}}{A factor vector containing cluster assignments.}
\item{\code{is_used}}{A logical vector indicating whether each cell is used. The length of the vector is the same as the length of the total cells. The sum of TRUE's is the same as \code{n_cells_total}. }
\item{\code{cts_params}}{A list containing cell type parameters.}
\item{\code{ld_params}}{A list containing dimensionality reduction parameters.}
\item{\code{ld_call}}{A call object representing the dimensionality reduction function call.}
\item{\code{clust_call}}{A call object representing the clustering function call.}
}
}

