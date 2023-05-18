#_ -------------------------------------------------------

# class: roi ----
#' @export
setClass("roi",
   slots = c(
     direction = "character",
     cycle = "numeric",
     roi_type = "data.frame",
     coords = "data.frame"))

# class: file_paths ----
#' @export
setClass("file_paths",
   slots = c(
     sample_name = "character",
     root_dir = "character",
     mask_type = "character",
     registration = "character",
     segmentation = "character",
     quantification = "character")
)

#_ -------------------------------------------------------

# class: CellTypeDefault ----

#' A class that defines cell types and antibodies in a CyCIF project (defunct)
#'
#' @export
setClass("CellTypeDefault",
         slots = c(
           markers = "data.frame",
           cell_lineage_def = "data.frame",
           cell_state_def = "data.frame"
         )
)

# class: CellTypeCycif ----

#' A class that defines cell types in Cycif class, extending CellTypeDefault class
#'
#' @export
setClass("CellTypeCycif",contains="CellTypeDefault",
         slots = c(
           name = "character",
           n_cycles = "numeric",
           gates = "data.frame",
           cell_types = "factor", # generated from gates and lineage
           is_strict = "logical"
         )
)

# class: CellTypeCycifStack ----

#' A class that defines cell types in CycifStack class, extending CellTypeDefault class
#'
#' @export
setClass("CellTypeCycifStack",contains="CellTypeCycif",
         slots = c(
           n_samples = "numeric"
         )
)


#_ -------------------------------------------------------

# class: LDCoords ----

#' Coordinates for dimensionality reduction
#'
#' A class that contains outputs of dimensionality reduction algos
#'
#'@export
setClass("LDCoords",
   slots = c(
     # [ld_type, nom_type]
     ld_type = "character", # PCA, TSNE, UMAP
     norm_type = "character",

     # [used smpls, abs, celltypes]
     smpls = "character",
     used.abs = "character",
     used.cts = "character",

     # [n_cells]
     n_cells_per_smpl = "numeric",
     n_cells_total = "numeric",

     # [ld_coords, clusters]
     ld_coords = "data.frame",
     clusters = "factor",

     # [is_used, cts_params,ld_params]
     is_used = "logical",
     cts_params = "list",
     ld_params = "list",

     # [call]
     ld_call = "call",
     clust_call = "call"
   )
)

#_ -------------------------------------------------------

# class: neighborhood ----

# this class will be used by spatial_count and spatial_cluster

#' neighborhoods
#'
#' A class that contains outputs of neighborhoods
#'
#'@export
setClass("neighborhood",
         slots = c(
           # [neighbor_tyoe]
           neighbor_type = "character", # spatial_count, spatial_lda

           # [used smpls, abs, celltypes]
           smpls = "character",
           used.abs = "character",
           used.cts = "character",

           # [n_cells]
           n_cells_per_smpl = "numeric",
           n_cells_total = "numeric",

           # [ld_coords, clusters]
           neighborhood = "factor",
           ne_clusters = "factor",

           # [is_used, cts_params,ld_params]
           is_used = "logical",
           cts_params = "list",
           ld_params = "list",

           # [call]
           ld_call = "call",
           clust_call = "call"
         )
)

#_ -------------------------------------------------------

# class: Cycif ----

#' A class that contains data of one sample in a CyCIF project.
#'
#' @export
setClass("Cycif",
   slots = c(
     # [sample basic info]
     name = "character", # '13676'
     file_paths = "file_paths",
     mask_type = "character", # 'cellRing' => should be mask_type

     # [n_cycles, n_cells]
     n_cycles = "numeric", # n_cycles
     n_cells = "numeric", # total cells

     # [abs_list]
     abs_list = "data.frame", # cycle <= n_cycles (nrow = 3*n_cycles)

     # [quantification]
     dna = "data.frame", # ncol = n_cycles; nrow = n_cells
     raw = "data.frame", # ncol = n_cycles; nrow = n_cells

     # [coordinates]
     xy_coords = "data.frame", # X_centroid, Y_centroid, nrow = n_cells; not trimmed by qc metric

     # [properties]
     segment_property = "data.frame", # Area, MajorAxisLength, MinorAxisLength, Eccentricity, Solidity, Extent, Orientation

     # [qc]
     used_cells = "matrix", # numeric matrix, 0's and 1's, ncol = # of cycles
     rois= "list", # defined in Omero
     within_rois = "logical", # logical vector - chose column only for that @n_cycles
     dna_thres = "data.frame", # nrow = n_cycle, ncol =2 (low, high)

     # [normalized]
     log_normalized = "data.frame",
     logTh_normalized = "data.frame",

     # [cell types]
     cell_type = "list", # CellTypeCycif

     # [ld_coords, clusters]
     ld_coords = "list", # list of ld_coords object
     clusters = "numeric", # numeric

     # [scimap,spatial]
     adata = "anndata._core.anndata.AnnData",

     neighbors = "list", #

     # [call]
     calls = "list" # list of functions called
   )
)

# class: CycifStack ----
#' A class that contains data of multiple samples from a CyCIF project.
#'
#' @export
setClass("CycifStack",
   slots = c(
     # [sample basic info]
     names = "character",
     mask_type = "character",

     # [n_cycles, n_cells]
     n_samples = "numeric",
     n_cycles = "numeric",
     max_cycles = "numeric",
     n_cells = "numeric",

     # [abs_list]
     abs_list = "data.frame",

     # [normalized]
     log_normalized = "data.frame",
     logTh_normalized = "data.frame",

     # [cell types]
     cell_type = "list", # CellTypeCycifStack

     # [ld_coords, clusters]
     ld_coords = "list",
     samples = "list",

     # [ranges]
     # qt_gates = "vector", # quantification gates - for cell type calling
     # ch_ranges = "vector", # channel gates

     # [phenotypes]
     phenoData = "data.frame",

     # [call]
     calls = "list"
   )
)


#_ -------------------------------------------------------
