#' A class that defines cell types and antibodies in a CyCIF project (defunct)
#'
#'
setClass("CellTypeDef",
   slots = c(
     cell_lineage_df = "data.frame",
     cell_state_mks = "character"
   )
)

#' A class that defines cell types in Cycif class, extending CellTypeDef class
#'
#' @export
setClass("CellTypeCycif",contains="CellTypeDef",
   slots = c(
     # cell_lineage_df = "data.frame",
     # cell_state_mks = "character",

     gates = "numeric",
     cell_types = "character" # generated from gates and lineage
   )
)

#' A class that defines cell types in CycifStack class, extending CellTypeDef class
#'
#' @export
setClass("CellTypeCycifStack",contains="CellTypeDef",
   slots = c(
     # cell_lineage_df = "data.frame",
     # cell_state_mks = "character",

     gates = "data.frame",
     cell_types = "character"
   )
)

