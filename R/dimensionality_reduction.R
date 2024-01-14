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

# fun: RunUMAP Cycif, CycifStack ----

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
#' @param init.seed Numeric specifying the initial seed for reproducible UMAP results.
#' @param save.coords Logical indicating whether to save the computed UMAP coordinates in the object. Default is FALSE.
#' @param ... Additional arguments passed to uwot::umap().
#'
#' @return A Cycif or CycifStack object with added UMAP or t-SNE coordinates and related metadata.
#'
#' @details
#' The `RunUMAP` function performs dimensionality reduction on CyCIF data using UMAP or t-SNE algorithms, depending on the specified parameters. It generates UMAP or t-SNE coordinates for visualization purposes and adds them to the provided Cycif or CycifStack object. The resulting object contains the computed coordinates, along with metadata such as the type of dimensionality reduction, normalization type, used antibodies, and used cell types.
#'
#' @seealso \code{\link{LDCoords}}, \code{\link{Cycif}}, \code{\link{CycifStack}}
#'
#' @export
setGeneric("RunUMAP", function(x,...) standardGeneric("RunUMAP"))

#' @rdname RunUMAP
#' @export
setMethod("RunUMAP", "Cycif",
          function(x,
                   norm_type=c("raw","log","logTh"),
                   ld_name,
                   ct_name="default",
                   ncells.per.smpl=NULL,
                   used.abs,
                   used.cts,
                   strict=TRUE,
                   n_neighbors=20,
                   init.seed=12345,...){
            call1 <- sys.calls()[[1]]
            if(missing(ld_name)){
              stop("'ld_name' should be specified (it's used to retrieve the data later)")
            }

            ## norm_type - by default, should use normalized value
            if(missing(norm_type)){
              norm_type <- "logTh_normalized"
            }

            ## used.abs
            if(missing(used.abs)){
              stop("'used.abs' should be defined first (have you set threshold?).")
            }

            ## is.used
            cts <- cell_types(x,ct_name=ct_name,strict=strict)$cell_type
            levels(cts) <- sub(",.+","",levels(cts))

            if(missing(used.cts)){
              used.cts <- levels(cts)
            }
            used.cts <- used.cts[!used.cts %in% c("unknown","outOfROI")]
            is.used <- cts %in% used.cts

            ## Select 'ncells.per.smpl' cells from available.
            smpl <- names(x)
            n.used <- sum(is.used)
            if(!is.null(ncells.per.smpl) & missing(ncells.per.smpl)){
              ncells.per.smpl <- n.used
            }else if(n.used < ncells.per.smpl){
              ncells.per.smpl <- n.used
              warning(smpl, ": try sampling ",ncells.per.smpl," cells but only ",n.used," cells available.\n")
            }else{
              cat(smpl, ": sampling ",ncells.per.smpl," out of ",n.used," cells.\n")
            }

            set.seed(123)
            used.idx <- sample(which(is.used),ncells.per.smpl)
            is.used.1 <- seq(is.used) %in% used.idx

            ## exprs matrix
            mat <- exprs(x,type=norm_type)
            mat <- mat[is.used.1,used.abs]

            set.seed(init.seed)
            ru <- uwot::umap(mat,n_neighbors=n_neighbors,...)

            ru <- data.frame(ru)
            rownames(ru) <- which(is.used.1)
            names(ru) <- c("x","y")

            ld <- LDCoords(
              ld_type = "UMAP",
              norm_type = norm_type,
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
                init.seed = init.seed,
                n_neighbors = n_neighbors
              ),
              ld_call=call1)

            x@ld_coords[[ld_name]] <- ld
            return(x)
          })

#' @rdname RunUMAP
#' @export
setMethod("RunUMAP", "CycifStack",
          function(x,
                   norm_type=c("raw","log","logTh"),
                   ld_name,
                   ct_name="default",
                   smpls,
                   used.abs,
                   used.cts,
                   strict=TRUE,
                   ncells.per.smpl=NULL,
                   n_neighbors=20,
                   init.seed=12345,
                   save.coords=FALSE,...){
            call1 <- sys.calls()[[1]]
            if(missing(ld_name)){
              stop("'ld_name' should be specified (it's used to retrieve the data later)")
            }

            ## norm_type - by default, should use normalized value
            if(missing(norm_type)){
              norm_type <- "logTh_normalized"
            }

            ## used.abs
            if(missing(used.abs)){
              stop("'used.abs' should be defined first (have you set threshold?).")
            }else if(any(!used.abs %in% abs_list(x)$ab)){
              missing.abs <- used.abs[(!used.abs %in% abs_list(x)$ab)]
              stop("missing abs in 'used.abs': ",paste(missing.abs,collapse=", "))
            }

            ## samples
            if(missing(smpls)){
              smpls <- names(x)
            }else if(any(!smpls %in% names(x))){
              missing.smpls <- smpls[!smpls %in% names(x)]
              stop("missing smpls in 'smpls': ",paste(missing.smpls,collapse=", "))
            }

            ## celltypes
            cts.df <- cell_types(x,ct_name=ct_name)
            cts <- cts.df$cell_types
            # ncts <- cyApply(x,function(y)length(cell_types(y,ct_name=ct_name)),simplify=T)
            levels(cts) <- sub(",.+","",levels(cts))
            if(missing(used.cts)){
              used.cts <- levels(cts)
              used.cts <- used.cts[!used.cts %in% c("unknown","outOfROI")]
            }else if(!all(used.cts %in% levels(cts))){
              stop("all 'used.cts' should be observed cell types")
            }
            is.used <- cts %in% used.cts

            ## ncells.per.smpl
            if(is.null(ncells.per.smpl) || missing(ncells.per.smpl)){
              smpls <- rep(names(x),nCells(x))
              max.per.smpls <- max(table(is.used,smpls)["TRUE",])
              ncells.per.smpl <- max.per.smpls
            }


            ## set 'n.cells'
            ncells <- nCells(x)
            v.ncells <- rep(names(x),ncells) # sample names
            n.used <- table(factor(is.used,levels=c("TRUE","FALSE")),v.ncells)["TRUE",]

            ## sample 'n.cells' cells per sample
            idx.list.mat <- cyApply(x,function(y){
              smpl <- names(y)

              set.seed(123)
              is.used.1 <- is.used[v.ncells==smpl]
              if(ncells.per.smpl < sum(is.used.1)){
                used.idx <- sample(which(is.used.1),ncells.per.smpl)
                is.used.2 <- seq(is.used.1) %in% used.idx
              }else{
                is.used.2 <- is.used.1
              }
            },as.CycifStack=FALSE)

            is.used.2 <- unlist(idx.list.mat)

            mat <- exprs(x,type=norm_type)[is.used.2,used.abs,drop=F]
            has.na <- apply(is.na(mat),1,any)

            mat <- mat[!has.na,]

            is.used.3 <- is.used.2
            is.used.3[which(is.used.2)[has.na]] <- FALSE

            ##
            set.seed(init.seed)
            ru <- uwot::umap(mat,n_neighbors=n_neighbors,...)
            ru <- data.frame(ru)
            rownames(ru) <- rownames(mat)
            names(ru) <- c("x","y")

            ru <- ru %>% cbind(cts.df[is.used.3,])
            n.smpls <- table(ru$sample)
            idx.per.smpl <- unlist(lapply(n.smpls,seq))
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
                init.seed = init.seed,
                n_neighbors = n_neighbors
              ),
              ld_call=call1,
              clust_call=call("function"))

            x@ld_coords[[ld_name]] <- ld
            return(x)
            # validObject(x)
          })

#_ -------------------------------------------------------
