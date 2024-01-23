#_ -------------------------------------------------------

# class: roi ----

#' @title ROI (Region of Interest) class
#'
#' @description
#' This class represents a Region of Interest (ROI) used in the context of CyCIF data.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{direction}}{A character vector specifying the direction of the ROI.
#'   The direction is either "positive" or "negative", which correspond to the ROIs that are included or excluded from the analysis, respectively.}
#'   \item{\code{cycle}}{A numeric vector representing the cycle in the CyCIF used in the analysis.}
#'   \item{\code{roi_type}}{A character vector containing one element, that describes the type of the ROI. Currently accepted roi_types: Polygon, Ellipse, and Rectangle}
#'   \item{\code{coords}}{A data frame containing coordinates associated with the ROI.}
#' }
#'
#' @details
#' The \code{roi} class is used to store information about ROIs, including their direction, cycle, type, and associated coordinates.
#'
#' @seealso
#' Other classes: \code{\link{Cycif}}, \code{\link{CycifStack}}, \code{\link{CellTypes}}
#'
#' @examples
#' # Creating an roi object
#' roi_object <- new("roi",
#'                   direction = "positive",
#'                   cycle = 5,
#'                   roi_type = "Polygon",
#'                   coords = data.frame(x = c(1, 2, 3), y = c(4, 5, 6)))
#'
#'
#' @rdname roi
#' @export
setClass("roi",
   slots = c(
     direction = "character",
     cycle = "numeric",
     roi_type = "character",
     coords = "data.frame"))

# class: file_paths ----

#' @title
#' File Paths Class
#'
#' @description
#' This class represents file paths associated with CyCIF data processed with MCMICRO for a single sample.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{sample_name}}{A character vector specifying the sample name.}
#'   \item{\code{root_dir}}{A character vector indicating the root directory where data files are stored.}
#'   \item{\code{mask_type}}{A character vector representing the mask type used for processing only 'cell' and 'cellRing' are accepted.}
#'   \item{\code{registration}}{A character vector indicating the file path for registration data (ome.tif file).}
#'   \item{\code{segmentation}}{A character vector indicating the file path for segmentation data (ome.tif file).}
#'   \item{\code{quantification}}{A character vector indicating the file path for quantification data, also known as spatial feature table (csv file).}
#' }
#'
#' @details
#' The \code{file_paths} class is used to store information about the file paths for different steps of CyCIF data processing for a single sample.
#'
#' @seealso
#' Other classes: \code{\link{Cycif}}, \code{\link{CycifStack}}, \code{\link{CellTypes}}
#'
#' @examples
#' # Creating a file_paths object
#' file_paths_object <- new("file_paths",
#'                          sample_name = "Sample1",
#'                          root_dir = "/path/to/data/",
#'                          mask_type = "cell",
#'                          registration = "/path/to/registration/file",
#'                          segmentation = "/path/to/segmentation/file",
#'                          quantification = "/path/to/quantification/file")
#'
#' # Accessing slots
#' sample_name <- file_paths_object@sample_name
#' root_dir <- file_paths_object@root_dir
#' mask_type <- file_paths_object@mask_type
#' registration <- file_paths_object@registration
#' segmentation <- file_paths_object@segmentation
#' quantification <- file_paths_object@quantification
#'
#' @rdname file_paths
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

# class: CellTypes ----

#' @title CellTypes Class
#'
#' @description
#' This class represents the cell types and their characteristics derived from CyCIF data.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{cell_lineage_def}}{A data frame containing cell lineage definitions.}
#'   \item{\code{cell_state_def}}{A data frame containing cell state definitions.}
#'   \item{\code{markers}}{A data frame specifying markers associated with cell types.}
#'   \item{\code{n_cycles}}{A numeric value indicating the number of cycles in the CyCIF data.}
#'   \item{\code{sample_names}}{A character vector with sample names associated with the data.}
#'   \item{\code{cell_types}}{A factor representing the generated cell types per individual cells.}
#'   \item{\code{is_strict}}{A logical value indicating strictness in cell type calling.}
#' }
#'
#' @details
#' The \code{CellTypes} class stores information about cell lineage, cell state, markers, and sample-specific details from CyCIF data analysis.
#'
#' @seealso
#' Other classes: \code{\link{Cycif}}, \code{\link{CycifStack}}
#'
#' @rdname CellTypes
#' @export
setClass("CellTypes",# switch this to 'CellTypes' class, remove the others
         slots = c(
           cell_lineage_def = "data.frame",
           cell_state_def = "data.frame",
           markers = "data.frame",
           n_cycles = "numeric",
           sample_names = "character",
           cell_types = "factor", # generated from gates and lineage
           is_strict = "logical"
         )
)
# stop("CellTypes now has sample_names. defineCellTypes and cell_types() should change accordingly (should have the same # of cells with @cell_types)",
#      "then change LDCoords and spatial info so we can make the computation tomorrow")

#_ -------------------------------------------------------

# class: LDCoords ----

#' @title LDCoords Class
#'
#' @description This class represents the coordinates and clustering results from dimensionality reduction techniques for CyCIF data.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{ld_type}}{A character vector specifying the dimensionality reduction technique used: 'PCA', 't-SNE', 'UMAP'.}
#'   \item{\code{norm_type}}{A character vector specifying the normalization type: 'log', or 'logTh'}
#'   \item{\code{smpls}}{A character vector containing sample names.}
#'   \item{\code{used.abs}}{A character vector containing antibody names used for the analysis.}
#'   \item{\code{used.cts}}{A character vector containing cell type names used for the analysis.}
#'   \item{\code{n_cells_per_smpl}}{A numeric vector representing the max number of cells per sample selected.}
#'   \item{\code{n_cells_total}}{A numeric value representing the total number of cells used for the analysis.}
#'   \item{\code{ld_coords}}{A data frame containing the coordinates from dimensionality reduction.}
#'   \item{\code{clusters}}{A factor vector containing cluster assignments.}
#'   \item{\code{is_used}}{A logical vector indicating whether each cell is used. The length of the vector is the same as the length of the total cells. The sum of TRUE's is the same as \code{n_cells_total}. }
#'   \item{\code{cts_params}}{A list containing cell type parameters.}
#'   \item{\code{ld_params}}{A list containing dimensionality reduction parameters.}
#'   \item{\code{ld_call}}{A call object representing the dimensionality reduction function call.}
#'   \item{\code{clust_call}}{A call object representing the clustering function call.}
#' }
#'
#' @details
#' The \code{LDCoords} class represents the coordinates and clustering results from dimensionality reduction techniques applied to CyCIF data.
#' It provides information about the samples, used abs, cell types, and various parameters of the cells used in the analysis so the plot can be color-coded accordingly.
#'
#' @rdname LDCoords
#' @export
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
     clusters = "factor", # cluster should be contained within ld_coords

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

# class: frnn ----

#' Class "frNN" - output of dbscan::frNN()
#'
#' @slot id A list of integer vector (length of the number of cell neighborhoods). Each vector contains the ids of the fixed radius nearest neighbors.
#' @slot dist A list with distances (same structure as id)
#' @slot eps A numeric vector (scalar) containing neighborhood radius eps that was used.
#' @slot sort A logical value indicating whether the distances are sorted.
#'
#' @seealso \code{\link{dbscan::frNN}}
#'
#' @rdname frNN
#'
#' @export
setClass("frNN",
         slots = c(
           dist = "list",
           id = "list",
           eps = "numeric",
           sort = "logical"
         )
)

# class: CellNeighborhood ----

#' Class "CellNeighborhood" - Cell Neighborhood Object
#'
#' An object of class "CellNeighborhood" represents the results of the computeCN function
#' for a Cycif or CycifStack object, providing information about the cell neighborhood analysis.
#'
#' @slot within.rois A logical vector indicating whether each cell is within a region of interest (ROI).
#' @slot cts.in.rcn A character vector of cell types considered in the neighborhood analysis.
#' @slot n.cells.selected An integer indicating the number of cells selected for the analysis.
#' @slot smpls A character vector containing sample names.
#' @slot frnn A frNN object containing information about cell neighbors.
#' @slot n.neighbors An integer indicating the number of neighbors (including self).
#' @slot exp.per.ct.cn A data.table containing expression data for each cell type in the neighborhood.
#' @slot exp.per.cn A data.table containing expression data for each cell neighborhood.
#' @slot is.selected A logical vector indicating whether each cell is selected.
#' @slot rcn.count A data.frame containing the counts of cell types in the neighborhood.
#' @slot rcn.freq A data.frame containing the frequency of cell types in the neighborhood.
#' @slot dist2tumorBorder A numeric vector containing the distance of each cell to tumor border.
#' @slot mclustda A list containing the results of the mclustDA function.
#'
#' @seealso \code{\link{computeCN}}
#'
#' @rdname CellNeighborhood
#' @export
setClass("CellNeighborhood",
         slots = c(
           within.rois = "logical",
           n.cells.selected = "numeric",
           is.selected = "logical",
           cts.in.rcn = "character",
           smpls = "character",
           frnn = "frNN",
           n.neighbors = "numeric",
           exp.per.ct.cn = "data.table",
           exp.per.cn = "data.table",
           rcn.count = "matrix",
           rcn.freq = "matrix",
           dist2tumorBorder = "numeric",
           mclustda = "list"
         ))

#_ -------------------------------------------------------
# class: Cycif ----

#' Cycif Class
#'
#' This class represents CyCIF data and its associated metadata.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{name}}{A character vector representing the sample name correspoding to the CyCIF dataset.}
#'   \item{\code{file_paths}}{An object of class \code{file_paths} containing file paths for the data (not implemented).}
#'   \item{\code{mask_type}}{A character vector specifying the type of mask used, e.g., 'cell', 'cellRing'.}
#'   \item{\code{n_cycles}}{A numeric value indicating the total number of imaging cycles.}
#'   \item{\code{n_cells}}{A numeric value representing the total number of cells in the dataset (including the ones out of ROIs and dropped from dnaFilter).}
#'   \item{\code{abs_list}}{A data frame containing information about antibodies and the cycles they were used imaging cycles.}
#'   \item{\code{dna}}{A data frame representing raw DNA data per cell with columns corresponding to imaging cycles.}
#'   \item{\code{raw}}{A data frame representing raw protein data per cell with columns corresponding to imaging cycles.}
#'   \item{\code{xy_coords}}{A data frame containing X and Y coordinates of cell centroids.}
#'   \item{\code{segment_property}}{A data frame containing per-cell properties of image segmentations.}
#'   \item{\code{used_cells}}{A numeric matrix indicating the used cells with values, where 0 or 1 corresponds to not available or available for each cycle. The number of columns corresponds to the number of cycles.}
#'   \item{\code{rois}}{A list containing information about Regions of Interest (ROIs). Each element is an object of the roi class.}
#'   \item{\code{within_rois}}{A logical vector specifying whether cells are within ROIs for each imaging cycle.}
#'   \item{\code{dna_thres}}{A data frame representing DNA threshold values. The number of rows correspond to the number of cycles, and there are two columns that correspond to lower- and upper-boundary.}
#'   \item{\code{log_normalized}}{A data frame representing log-transformed normalized data, which are computed by \code{normalize} function.}
#'   \item{\code{logTh_normalized}}{A data frame representing log-transformed and thresholded normalized data, which are computed by \code{normalize} function.}
#'   \item{\code{cell_types}}{A named list of \code{CellTypes} objects representing cell types.}
#'   \item{\code{ld_coords}}{A list of \code{LDCoords} objects.}
#'   \item{\code{cell_neighborhood}}{An object of class \code{CellNeighborhood} representing cell neighborhood analysis results.}
#'   \item{\code{calls}}{A list of functions called.}
#' }
#'
#' @details
#' The \code{Cycif} class represents CyCIF data and its associated metadata, including sample name, file paths, imaging cycles, cell coordinates, cell properties, and more. This class is used to store and manipulate CyCIF data within R.
#'
#' @seealso
#' Other classes: \code{\link{CycifStack}}, \code{\link{CellTypes}}, \code{\link{LDCoords}}
#'
#' @rdname Cycif
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

           # [abs_list] - adding gates and dynamic range
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
           cell_types = "list", # CellTypeCycif

           # [ld_coords, clusters]
           ld_coords = "list", # list of ld_coords object
           # clusters = "numeric", # numeric

           # # [cellneighborhood]
           cell_neighborhood = "CellNeighborhood",

           # [call]
           calls = "list" # list of functions called
         )
)

# class: CycifStack ----
#' CycifStack Class
#'
#' This class represents a collection of CyCIF samples as a stack, allowing manipulation and analysis of multiple samples together.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{samples}}{A list of \code{Cycif} objects representing individual CyCIF samples.}
#'   \item{\code{names}}{A character vector containing names for each sample.}
#'   \item{\code{mask_type}}{A character vector specifying the type of mask used, e.g., 'cellRing'.}
#'   \item{\code{n_samples}}{A numeric value indicating the total number of CyCIF samples in the stack.}
#'   \item{\code{n_cycles}}{A numeric value indicating the maximum number of imaging cycles among all samples in the stack.}
#'   \item{\code{max_cycles}}{A numeric value indicating the maximum number of imaging cycles among all samples in the stack.}
#'   \item{\code{n_cells}}{A numeric value representing the total number of cells in the stack (including dropped cells).}
#'   \item{\code{abs_list}}{A data frame containing information about antibodies and the cycles they were used imaging cycles.}
#'   \item{\code{cell_types}}{A list of \code{CellType} objects representing cell types in the stack.}
#'   \item{\code{ld_coords}}{A list of \code{LDCoords} objects.}
#'   \item{\code{cell_neighborhood}}{An object of class \code{CellNeighborhood} representing cell neighborhood analysis results.}
#'   \item{\code{phenoData}}{A data frame containing phenotypic information.}
#'   \item{\code{calls}}{A list of functions called.}
#' }
#'
#' @details
#' The \code{CycifStack} class represents a collection of CyCIF samples organized as a stack. It allows users to analyze and manipulate multiple samples simultaneously, making it useful for high-throughput analysis of CyCIF data.
#'
#' @seealso
#' Other classes: \code{\link{Cycif}}, \code{\link{CellTypes}}, \code{\link{LDCoords}}
#'
#' @rdname CycifStack
#' @export
setClass("CycifStack",
         slots = c(
           # [list of Cycif objs]
           samples = "list",

           # [sample basic info]
           names = "character",
           mask_type = "character",

           # [n_cycles, n_cells]
           n_samples = "numeric",
           n_cycles = "numeric",
           max_cycles = "numeric",
           n_cells = "numeric",

           # [abs_list] - should contain gates and dynamic ranges
           abs_list = "data.frame",

           # [normalized] - should contain log_normalized and logTh_normalized
           # log_normalized = "data.frame",
           # logTh_normalized = "data.frame",

           # [cell types]
           cell_types = "list", # CellTypeCycifStack

           # [ld_coords]
           ld_coords = "list",

           # # [cellneighborhood]
           cell_neighborhood = "CellNeighborhood",

           # [phenotypes]
           phenoData = "data.frame",

           # [call]
           calls = "list"
         )
)
