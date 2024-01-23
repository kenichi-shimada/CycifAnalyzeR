#_ -------------------------------------------------------

# utils LDCoords ----

#' Get Names of UMAP or Clustering Results
#'
#' This function retrieves the names of UMAP embeddings or clustering results associated with a 'Cycif' or 'CycifStack' object.
#'
#' @param x A 'Cycif' or 'CycifStack' object containing UMAP embeddings or clustering results.
#'
#' @return A character vector containing the names of UMAP embeddings or clustering results.
#'
#' @details
#' The 'ld_names' function allows you to obtain the names of UMAP embeddings or clustering results stored within a 'Cycif' or 'CycifStack' object. These names can be used to reference specific UMAP embeddings or clustering outcomes for further analysis or visualization.
#'
#' @rdname ld_names
#' @export
setGeneric("ld_names", function(x)standardGeneric("ld_names"))

#' @rdname ld_names
#' @export
setMethod("ld_names", "Cycif",function(x)names(x@ld_coords))

#' @rdname ld_names
#' @export
setMethod("ld_names", "CycifStack",function(x)names(x@ld_coords))

#' Get UMAP Coordinates or Clustering Results
#'
#' This function retrieves UMAP coordinates or clustering results associated with a specified name from a 'Cycif' or 'CycifStack' object.
#'
#' @param x A 'Cycif' or 'CycifStack' object containing UMAP embeddings or clustering results.
#' @param ld_name The name of the UMAP embedding or clustering result to retrieve.
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
  if(missing(ld_name)){
    stop("'ld_name' should be specified to retrieve ld_coords")
  }else if(!ld_name %in% ld_names(x)){
    stop("Specified 'ld_name' doesn't exist")
  }
  return(x@ld_coords[[ld_name]])
})

#' @rdname ld_coords
#' @export
setMethod("ld_coords", "CycifStack",function(x,ld_name){
  if(missing(ld_name)){
    stop("'ld_name' should be specified to retrieve ld_coords")
  }else if(!ld_name %in% ld_names(x)){
    stop("Specified 'ld_name' doesn't exist")
  }
  return(x@ld_coords[[ld_name]])
})

#_ -------------------------------------------------------

# fun: constructor LDCoords ----
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
#' @export
setMethod("show", "LDCoords", function(object) {
  n.smpls <- length(object@smpls)
  used.abs <- object@used.abs
  n.abs <- length(used.abs)
  abs <- paste(used.abs,collapse=",")
  used.cts <- object@used.cts
  n.cts <- length(object@used.cts)
  cts <- paste(used.cts,collapse=",")

  n_cells_per_smpl <- object@n_cells_per_smpl
  n_cells_total <- object@n_cells_total

  cat("[",is(object)[[1]], " object]\n\n",
      "Type: ", object@ld_type, "\n",
      "Normalization: ",object@norm_type, "\n\n",
      "cts (",n.cts,") : ",cts,"\n",
      "abs (",n.abs,") : ",abs,"\n\n",
      "# samples : ", n.smpls,"\n",
      "# cells per smpl :\t", n_cells_per_smpl, "\n",
      "# cells in total :\t", sum(n_cells_total), "\n\n",sep="")
  if(length(n_cells_total)>1){
    print(n_cells_total)
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
#' @param smpls Character vector specifying the names of samples to be used for analysis.
#' @param ncells.per.smpl Numeric value specifying the number of cells per sample. If NULL, the maximum number of cells per sample will be used.
#' @param used.abs Character vector containing a set of antibodies used for UMAP computation.
#' @param used.cts Character vector containing a set of cell types used for UMAP computation.
#' @param strict Logical indicating whether strict filtering of cell types should be applied. Default is TRUE.
#' @param n_neighbors Numeric specifying the number of neighbors to consider during dimensionality reduction.
#' @param umap.seed Numeric specifying the initial seed for reproducible UMAP results.
#' @param save.coords Logical indicating whether to save the computed UMAP coordinates in the object. Default is FALSE.
#' @param ... Additional arguments passed to uwot::umap().
#'
#' @return A Cycif or CycifStack object with added UMAP or t-SNE coordinates and related metadata.
#'
#' @details
#' The `LdRunUMAP` function performs dimensionality reduction on CyCIF data using UMAP or t-SNE algorithms, depending on the specified parameters. It generates UMAP or t-SNE coordinates for visualization purposes and adds them to the provided Cycif or CycifStack object. The resulting object contains the computed coordinates, along with metadata such as the type of dimensionality reduction, normalization type, used antibodies, and used cell types.
#'
#' @seealso \code{\link{LDCoords}}, \code{\link{Cycif}}, \code{\link{CycifStack}}
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
            call1 <- sys.calls()[[1]]
            if(missing(ld_name)){
              stop("'ld_name' should be specified (it's used to retrieve the data later)")
            }

            ## norm_type - by default, should use normalized value ----
            if(missing(norm_type)){
              norm_type <- "log"
            }

            ## used.abs ----
            if(missing(used.abs)){
              stop("'used.abs' should be defined first")
            }else if(any(!used.abs %in% abs_list(x)$ab)){
              stop("Some abs in 'used.abs' are not available")
            }

            ## used.cts ----
            cts.df <- cell_types(x, ct_name = ct_name, strict = strict)
            cts <- cts.df$cell_type
            levels(cts) <- sub(",.+", "", levels(cts))
            if (missing(used.cts)) {
              used.cts <- levels(cts)
            }

            used.cts <- used.cts[used.cts != "outOfROI"]

            ## exprs ----
            mat <- exprs(x,type=norm_type)
            has.na <- apply(is.na(mat),1,any)
            if(length(has.na) != length(cts)){
              stop("exprs(x) and cell_types(x) should have the same number of rows")
            }

            ## is.used ----
            is.used <- cts %in% used.cts & !has.na

            is.used <- cts %in% used.cts
            n.used <- sum(is.used)

            ### Select 'ncells.per.smpl' cells ----
            smpl <- names(x)

            if(missing(ncells.per.smpl)){
              stop("`ncells.per.smpl` should be specified.")
            }else{
              ncells.per.smpl <- min(n.used,ncells.per.smpl)
            }

            set.seed(smpl.seed)
            used.idx <- sample(which(is.used),ncells.per.smpl)
            is.used.1 <- seq(is.used) %in% used.idx
            idx.per.smpl <- which(is.used.1)

            ## subsetting exprs ----
            mat <- mat[is.used.1,used.abs,drop=F]

            set.seed(umap.seed)
            ru <- uwot::umap(mat,n_neighbors=n_neighbors,...)

            ru <- data.frame(ru)
            rownames(ru) <- which(is.used.1)
            names(ru) <- c("x","y")

            ru <- ru %>% cbind(cts.df[is.used.1,])
            ru <- ru %>% mutate(idx = idx.per.smpl)

            ##
            ld <- LDCoords(
              ld_type = "UMAP",
              norm_type = norm_type,
              smpls = smpl,

              used.abs = used.abs,
              used.cts = used.cts,

              n_cells_per_smpl = ncells.per.smpl,
              n_cells_total = ncells.per.smpl,
              ld_coords = ru,
              is_used = is.used.1,
              cts_params = list(
                strict = strict,
                leaves.only= TRUE,
                within.rois = TRUE
              ),
              ld_params= list(
                umap.seed = umap.seed,
                n_neighbors = n_neighbors
              ),
              ld_call=call1)

            x@ld_coords[[ld_name]] <- ld
            return(x)
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
            call1 <- sys.calls()[[1]]
            if(missing(ld_name)){
              stop("'ld_name' should be specified (it's used to retrieve the data later)")
            }

            ## norm_type - by default, should use normalized value
            if(missing(norm_type)){
              norm_type <- "logTh"
            }

            ## used.abs
            if(missing(used.abs)){
              stop("'used.abs' should be defined first")
            }else if(any(!used.abs %in% abs_list(x)$ab)){
              stop("Some abs in 'used.abs' are not available")
            }

            ## used.cts ----
            cts.df <- cell_types(x, ct_name = ct_name, strict = strict)
            cts <- cts.all$cell_type
            smpls <- cts.all$sample
            levels(cts) <- sub(",.+", "", levels(cts))
            if (missing(used.cts)) {
              used.cts <- levels(cts)
            }

            used.cts <- used.cts[used.cts != "outOfROI"]

            ## exprs ----
            mat <- exprs(x,type=norm_type)
            has.na <- apply(is.na(mat),1,any)
            if(length(has.na) != length(cts)){
              stop("exprs(x) and cell_types(x) should have the same number of rows")
            }

            ## is.used ----
            is.used <- cts %in% used.cts & !has.na
            n.used <- table(is.used,smpls)["TRUE",]

            ### Select 'ncells.per.smpl' cells ----
            if(missing(ncells.per.smpl)){
              stop("`ncells.per.smpl` should be specified.")
            }else{
              ncells.per.smpl <- pmin(n.used,ncells.per.smpl)
            }

            ## list of is.used based on ncells.per.smpl and is.used for each sample
            set.seed(smpl.seed)
            uniq.smpls <- names(x)
            is.used.df <- lapply(uniq.smpls,function(smpl){
              is.used.1 <- is.used[smpls==smpl]
              if(ncells.per.smpl[smpl] < sum(is.used.1)){
                used.idx <- sample(which(is.used.1),ncells.per.smpl[smpl])
                is.used.2 <- seq(is.used.1) %in% used.idx
              }else{
                is.used.2 <- is.used.1
              }
              idx <- which(is.used.2)
              return(list(is.used=is.used.2,idx=idx))
            }) ## this value is used for is.used in LDCoords

            is.used.2 <- unlist(is.used.df,function(x)x$is.used)
            idx.per.smpl <- unlist(is.used.df,function(x)x$idx)

            ## subsetting exprs ----
            mat <- mat[is.used.2,used.abs,drop=F]

            ##
            set.seed(umap.seed)
            ru <- uwot::umap(mat,n_neighbors=n_neighbors,...)

            ru <- data.frame(ru)
            rownames(ru) <- rownames(mat)
            names(ru) <- c("x","y")

            ru <- ru %>% cbind(cts.df[is.used.2,])
            ru <- ru %>% mutate(idx = idx.per.smpl)

            ld <- LDCoords(
              ld_type = "UMAP",
              norm_type = norm_type,
              smpls = as.character(ru$sample),
              used.abs = used.abs,
              used.cts = used.cts,
              n_cells_per_smpl = ncells.per.smpl,
              n_cells_total = nrow(ru),
              ld_coords = ru,
              is_used = is.used.3,
              cts_params = list(
                strict = strict,
                leaves.only= TRUE,
                within.rois = TRUE
              ),
              ld_params= list(
                umap.seed = umap.seed,
                n_neighbors = n_neighbors
              ),
              ld_call=call1,
              clust_call=call("function"))

            x@ld_coords[[ld_name]] <- ld
            return(x)
            # validObject(x)
          })

#_ -------------------------------------------------------
