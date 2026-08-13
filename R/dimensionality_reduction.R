#_ -------------------------------------------------------

# utils LDCoords ----

## ld_names ----

#' Get Names of UMAP or Clustering Results
#'
#' This function retrieves the names of UMAP embeddings or clustering results associated with a 'Cycif' or 'CycifStack' object.
#'
#' @param x A 'Cycif' or 'CycifStack' object containing UMAP embeddings or clustering results.
#' @param simplify Logical, if TRUE (default) print just each name; if FALSE, print the full \code{show()} of each LDCoords object.
#' @param show Logical, whether to print the names (default TRUE); if FALSE, return them invisibly with no printing.
#' @param ... Additional arguments (unused).
#'
#' @return A character vector containing the names of UMAP embeddings or clustering results.
#'
#' @details
#' The 'ld_names' function allows you to obtain the names of UMAP embeddings or clustering results stored within a 'Cycif' or 'CycifStack' object. These names can be used to reference specific UMAP embeddings or clustering outcomes for further analysis or visualization.
#'
#' @rdname ld_names
#' @export
setGeneric("ld_names", function(x,...)standardGeneric("ld_names"))

#' @rdname ld_names
#' @export
setMethod("ld_names", "Cycif",function(x,simplify=TRUE,show=T){
  .Defunct(msg=paste(
    "ld_names() has had zero call sites anywhere across CycifAnalyzeR or any",
    "downstream project repo (EP, TALAVE, HR-APM) -- superseded by CellFeatures",
    "(subsetCells() + runUMAP()/clusterCells()), which needs no named ld_coords",
    "registry."
  ))
  }
)

#' @rdname ld_names
#' @export
setMethod("ld_names", "CycifStack",function(x,simplify=TRUE,show=T){
  .Defunct(msg=paste(
    "ld_names() has had zero call sites anywhere across CycifAnalyzeR or any",
    "downstream project repo (EP, TALAVE, HR-APM) -- superseded by CellFeatures",
    "(subsetCells() + runUMAP()/clusterCells()), which needs no named ld_coords",
    "registry."
  ))
  }
)

## ld_coords ----

#' Get UMAP Coordinates or Clustering Results
#'
#' This function retrieves UMAP coordinates or clustering results associated with a specified name from a 'Cycif' or 'CycifStack' object.
#'
#' @param x A 'Cycif' or 'CycifStack' object containing UMAP embeddings or clustering results.
#' @param ld_name The name of the UMAP embedding or clustering result to retrieve.
#' @param ... Additional arguments (unused).
#'
#' @return UMAP coordinates or clustering results specified by 'ld_name'.
#'
#' @details
#' The 'ld_coords' function allows you to obtain UMAP coordinates or clustering outcomes stored within a 'Cycif' or 'CycifStack' object. You need to specify the name of the UMAP embedding or clustering result you want to retrieve using the 'ld_name' parameter. If the specified 'ld_name' does not exist in the object, an error will be raised.
#'
#' @rdname ld_coords
#' @export
setGeneric("ld_coords", function(x,...)standardGeneric("ld_coords"))

#' @rdname ld_coords
#' @export
setMethod("ld_coords", "Cycif",function(x,ld_name){
  .Defunct(msg=paste(
    "ld_coords() has had zero call sites anywhere across CycifAnalyzeR or any",
    "downstream project repo (EP, TALAVE, HR-APM) -- superseded by CellFeatures",
    "(subsetCells() + runUMAP()/clusterCells())."
  ))
})

#' @rdname ld_coords
#' @export
setMethod("ld_coords", "CycifStack",function(x,ld_name){
  .Defunct(msg=paste(
    "ld_coords() has had zero call sites anywhere across CycifAnalyzeR or any",
    "downstream project repo (EP, TALAVE, HR-APM) -- superseded by CellFeatures",
    "(subsetCells() + runUMAP()/clusterCells())."
  ))
})

#_ -------------------------------------------------------

# fun: constructor LDCoords ----

#' Construct an LDCoords object
#'
#' @param ld_type One of "PCA", "tSNE", "UMAP".
#' @param norm_type Normalization type used for the input expression matrix ("log" or "logTh").
#' @param smpls Character vector of sample names.
#' @param used.abs Character vector of antibodies used for the analysis.
#' @param used.cts Character vector of cell types used for the analysis.
#' @param n_cells_per_smpl Numeric, max number of cells per sample selected.
#' @param n_cells_total Numeric, total number of cells used for the analysis.
#' @param ld_coords A data.frame of the computed coordinates.
#' @param is_used A logical vector, length = total number of cells, TRUE for cells used.
#' @param cts_params A list of cell-type parameters.
#' @param ld_params A list of dimensionality-reduction parameters.
#' @param ld_call The call used to compute the coordinates.
#' @param clust_call The call used to compute the clustering.
#'
#' @return An LDCoords object.
#'
#' @rdname LDCoords
#' @export
LDCoords <- function(ld_type,norm_type,smpls,used.abs,used.cts,
                     n_cells_per_smpl,n_cells_total,
                     ld_coords,is_used,
                     cts_params,ld_params,
                     ld_call,clust_call){
  if(!ld_type %in% c("PCA","tSNE","UMAP")){
    stop("ld_type should be PCA, tSNE, or UMAP.")
  }
  new("LDCoords",
      ld_type = ld_type,
      norm_type = norm_type,

      smpls = smpls,
      used.abs = used.abs,
      used.cts = used.cts,

      n_cells_per_smpl = n_cells_per_smpl,
      n_cells_total = n_cells_total,

      ld_coords = ld_coords,
      is_used = is_used,
      cts_params = cts_params,
      ld_params = ld_params,
      ld_call=ld_call,
      clust_call=clust_call
  )
}

# fun: show LDCoords ----
#' @param object An LDCoords object (for the \code{show} method).
#' @rdname LDCoords
#' @export
setMethod("show", "LDCoords", function(object) {
  used.abs <- object@used.abs
  n.abs <- length(used.abs)
  abs <- paste(used.abs,collapse=", ")
  used.cts <- object@used.cts
  n.cts <- length(object@used.cts)
  cts <- paste(used.cts,collapse=", ")

  n.cells.per.smpl <- max(object@n_cells_per_smpl)
  n.smpls <- length(object@n_cells_per_smpl)
  mean.n.cells <- round(mean(object@n_cells_per_smpl),2)
  n.cells.total <- object@n_cells_total

  # aaa
  # is clustering performed yet?
  # ld_params?
  # clust_params?

  cat("[",is(object)[[1]], " object]\n\n",
      "Type: ", object@ld_type, "\n",
      "Normalization: ",object@norm_type, "\n\n",
      "Cell types (",n.cts,") : \n ",cts,"\n",
      "Abs (",n.abs,") : \n ",abs,"\n\n",
      "# samples : ", n.smpls,"\n",
      "# cells per smpl : ", mean.n.cells, "\t (max: ", n.cells.per.smpl, ")\n",
      "# cells in total : ", sum(n.cells.total), "\n\n",sep="")
  if(length(n.cells.total)>1){
    print(n.cells.total)
  }
})

#_ -------------------------------------------------------

# fun: LdRunUMAP Cycif, CycifStack ----

#' Run visualization with dimensionality reduction (t-SNE or U-MAP) on CyCIF data.
#'
#' This function performs dimensionality reduction (t-SNE or U-MAP) on CyCIF data, either for a single sample or across multiple samples (CycifStack). It generates UMAP or t-SNE coordinates for visualization purposes.
#'
#' @param x A Cycif or CycifStack object.
#' @param norm_type Character string specifying the type of expression values to use for dimensionality reduction. It should be one of "raw," "log," or "logTh" indicating whether to use raw or normalized data.
#' @param ld_name Character string specifying the name to assign to the dimensionality reduction results. This name will be used to retrieve the data later.
#' @param ct_name Character string specifying the cell type name to be used for analysis. Default is "default."
#' @param ncells.per.smpl Numeric value specifying the number of cells per sample. If NULL, the maximum number of cells per sample will be used.
#' @param used.abs Character vector containing a set of antibodies used for UMAP computation.
#' @param used.cts Character vector containing a set of cell types used for UMAP computation.
#' @param strict Logical indicating whether strict filtering of cell types should be applied. Default is TRUE.
#' @param n_neighbors Numeric specifying the number of neighbors to consider during dimensionality reduction.
#' @param smpl.seed Numeric specifying the random seed used for per-sample cell subsampling. Default 123.
#' @param umap.seed Numeric specifying the initial seed for reproducible UMAP results.
#' @param ... Additional arguments passed to uwot::umap().
#'
#' @return A Cycif or CycifStack object with added UMAP or t-SNE coordinates and related metadata.
#'
#' @details
#' The `LdRunUMAP` function performs dimensionality reduction on CyCIF data using UMAP or t-SNE algorithms, depending on the specified parameters. It generates UMAP or t-SNE coordinates for visualization purposes and adds them to the provided Cycif or CycifStack object. The resulting object contains the computed coordinates, along with metadata such as the type of dimensionality reduction, normalization type, used antibodies, and used cell types.
#'
#' @seealso \code{\linkS4class{LDCoords}}, \code{\link{Cycif}}, \code{\link{CycifStack}}
#'
#' @importFrom uwot umap
#'
#' @export
setGeneric("LdRunUMAP", function(x,...) standardGeneric("LdRunUMAP"))

#' @rdname LdRunUMAP
#' @export
setMethod("LdRunUMAP", "Cycif",
          function(x,
                   norm_type=c("raw","log","logTh"),
                   ld_name,
                   ct_name="default",
                   ncells.per.smpl=NULL,
                   used.abs,
                   used.cts,
                   strict=TRUE,
                   n_neighbors=20,
                   smpl.seed=123,
                   umap.seed=12345,...){
            .Defunct(msg=paste(
              "LdRunUMAP() has had zero call sites anywhere across CycifAnalyzeR",
              "or any downstream project repo (EP, TALAVE, HR-APM). Use",
              "subsetCells() + runUMAP() on a CellFeatures object instead."
            ))
          })

#' @rdname LdRunUMAP
#' @export
setMethod("LdRunUMAP", "CycifStack",
          function(x,
                   norm_type=c("raw","log","logTh"),
                   ld_name,
                   ct_name="default",
                   ncells.per.smpl=NULL,
                   used.abs,
                   used.cts,
                   strict=TRUE,
                   n_neighbors=20,
                   smpl.seed=123,
                   umap.seed=12345,...){
            .Defunct(msg=paste(
              "LdRunUMAP() has had zero call sites anywhere across CycifAnalyzeR",
              "or any downstream project repo (EP, TALAVE, HR-APM). Use",
              "subsetCells() + runUMAP() on a CellFeatures object instead."
            ))
          })

#_ -------------------------------------------------------
