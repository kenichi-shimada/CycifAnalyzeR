% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{barplotCTS}
\alias{barplotCTS}
\title{Create a Barplot for Cell Type Composition (TALAVE)}
\usage{
barplotCTS(
  n.sh,
  anno.smpls,
  uniq.cols,
  ylim = c(0, 1),
  ylab = "CellType Composition"
)
}
\arguments{
\item{n.sh}{A matrix or data frame containing cell type composition data, with cell types as rows and samples as columns.}

\item{anno.smpls}{A data frame containing sample annotations, with sample names matching the column names of \code{n.sh}.}

\item{uniq.cols}{A vector of unique colors to use for plotting cell types.}

\item{ylim}{A numeric vector specifying the limits of the y-axis (default is c(0, 1)).}

\item{ylab}{Label for the y-axis (default is "CellType Composition").}
}
\description{
This function generates a composite barplot to visualize the composition of cell types across samples.
}
\details{
The \code{barplotCTS} function creates a composite barplot with the following features:
\itemize{
\item Visualizes the composition of various cell types across multiple samples.
\item Each bar in the plot represents a sample, and the height of the bar is proportional to the composition of cell types within that sample.
\item The cell types are color-coded using unique colors specified in the \code{uniq.cols} parameter for easy identification.
\item Sample annotations, such as patient ID, disease subtype, time point, and others, are displayed below the main barplot to provide additional context.
\item The y-axis represents the composition proportion of cell types, and you can customize the y-axis label using the \code{ylab} parameter.
\item The \code{ylim} parameter allows you to set specific limits for the y-axis, controlling the range of the composition proportions displayed.
}

This function is useful for visualizing and comparing cell type compositions across different samples, making it particularly valuable in biological and medical research for understanding the distribution of cell types in complex datasets.
}
