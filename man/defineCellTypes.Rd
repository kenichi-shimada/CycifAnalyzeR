% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/celltype.R
\name{defineCellTypes}
\alias{defineCellTypes}
\alias{defineCellTypes,data.frame-method}
\alias{defineCellTypes,Cycif-method}
\alias{defineCellTypes,CycifStack-method}
\title{Define cell types}
\usage{
defineCellTypes(x, ...)

defineCellTypes(x,ctype,cstate,gates,p_thres=0.5,...)

defineCellTypes(x,ctype,cstate,gates,p_thres=0.5,...)

\S4method{defineCellTypes}{CycifStack}(
  x,
  ctype,
  cstate,
  ct_name = "default",
  p_thres = 0.5,
  overwrite = FALSE,
  ...
)
}
\arguments{
\item{x}{a CycifStack object}

\item{...}{additional arguments (currently unused)}

\item{ctype}{a data.frame containing cell type definition}

\item{cstate}{a data.frame containing cell state definition}

\item{gates}{a data.frame containing gates (n.samples x n.proteins)}

\item{p_thres}{numerical between 0 and 1. A probability that corresponds to a threshold intensity}

\item{ct_anme}{name of the cell types}
}
\value{
a modified CycifStack object with updated cell type information
}
\description{
Perform a cell type calling function and set cell types in a Cycif or CycifStack object
}
