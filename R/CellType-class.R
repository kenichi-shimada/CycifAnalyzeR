#' A class that defines cell types and antibodies in a CyCIF project (defunct)
#'
#'
setClass("CellTypeDef",
         slots = c(
           cell_lineage_df = "data.frame",
           cell_state_df = "data.frame",
           markers = "list"
         )
)

#' A class that defines cell types in Cycif class, extending CellTypeDef class
#'
#' @export
setClass("CellTypeCycif",contains="CellTypeDef",
       slots = c(
         cell_lineage_df = "data.frame",
         cell_state_df = "data.frame",

         threshold = "data.frame",

         used_abs = "character", # generated from threshold and lineage
         cell_types = "character" # generaed from threshold and lineage
       )
)

#' A class that defines cell types in CycifStack class, extending CellTypeDef class
#'
#' @export
setClass("CellTypeCycifStack",contains="CellTypeDef",
         slots = c(
           cell_lineage_df = "data.frame",
           cell_state_df = "data.frame",

           threshold = "data.frame",

           used_abs = "character", # user-defined value
           cell_types = "character"
         )
)

