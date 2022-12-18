#' A class that defines cell types and antibodies in a CyCIF project (defunct)
#'
#'
setClass("CellTypeDefault",
   slots = c(
     markers = "data.frame",
     cell_lineage_df = "data.frame",
     cell_state_df = "data.frame",
     expanded_lineage_df = "data.frame",
     expanded_state_df = "data.frame"
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
     cell_types_full = "factor"
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
     cell_types_full = "factor"
   )
)

