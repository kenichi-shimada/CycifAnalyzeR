#' A class that defines ld_coords
#'
#'@export
setClass("LDCoords",
   slots = c(
     type = "character", # PCA, TSNE, UMAP

     smpls = "character",
     used.abs = "character",
     used.cts = "character",

     n_cells_per_smpl = "numeric",
     n_cells_total = "numeric",
     ld_coords = "data.frame",
     is_used = "logical",
     ld_params = "list"
   )
)
