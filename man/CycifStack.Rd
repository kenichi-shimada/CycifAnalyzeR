% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cycifstack.R
\name{CycifStack}
\alias{CycifStack}
\alias{show,CycifStack-method}
\title{Instantiate and show a CycifStack object}
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
\item{ft_filenames}{quantification file (.csv)}

\item{path}{path to the dir containign filename}

\item{mask_type}{a common mask_type of each channel. If applicable, the mask_type is
removed from the channel names.}

\item{object}{a CycifStack object}
}
\description{
Instantiate and show a CycifStack object
}
