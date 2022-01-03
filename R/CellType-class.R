#' A class that defines cell types and antibodies in a CyCIF project
#'
#' @export
setClass("CellTypeDef",
         slots = c(
           cell_lineage_df = "data.frame",
           cell_state_df = "data.frame",
           lineages = "character",
           markers = "list"
         )
)

#' A class that defines cell types in Cycif class, extending CellTypeDef class
#'
#' @export
setClass("CellTypeCycif",contains="CellTypeDef",
         slots = c(
           threshold = "data.frame",

           used_abs = "character", # user-defined value
           n_cycles = "numeric", # user-defined value, can be used to subset used_abs

           uniq_abs = "data.frame", # import in a separate function
           meta = "character"
         )
)

#' A class that defines cell types in CycifStack class, extending CellTypeDef class
#'
#' @export
setClass("CellTypeCycifStack",contains="CellTypeDef",
         slots = c(
           threshold = "data.frame",

           used_abs = "character", # user-defined value
           n_cycles = "numeric", # user-defined value, can be used to subset used_abs

           uniq_abs = "data.frame", # import in a separate function
           meta = "character"
         )
)

