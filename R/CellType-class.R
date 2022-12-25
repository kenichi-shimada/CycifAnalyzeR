#' A class that defines cell types and antibodies in a CyCIF project (defunct)
#'
#'
setClass("CellTypeDefault",
   slots = c(
     markers = "data.frame",
     cell_lineage_def = "data.frame",
     cell_state_def = "data.frame"
   )
)

#' A class that defines cell types in Cycif class, extending CellTypeDefault class
#'
#' @export
setClass("CellTypeCycif",contains="CellTypeDefault",
   slots = c(
     name = "character",
     n_cycles = "numeric",
     gates = "numeric",
     cell_types = "factor", # generated from gates and lineage
     is_strict = "character"
   )
)

#' A class that defines cell types in CycifStack class, extending CellTypeDefault class
#'
#' @export
setClass("CellTypeCycifStack",contains="CellTypeDefault",
   slots = c(
     n_samples = "numeric",
     max_cycles = "numeric",
     gates = "data.frame",
     cell_types = "factor",
     is_strict = "character"
   )
)

