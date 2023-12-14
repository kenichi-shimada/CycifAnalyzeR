#_ -------------------------------------------------------

# utils LDCoords ----
#' @export
setGeneric("ld_names", function(x)standardGeneric("ld_names"))
setMethod("ld_names", "Cycif",function(x)names(x@ld_coords))

#' @export
setMethod("ld_names", "CycifStack",function(x)names(x@ld_coords))

#' @export
setGeneric("ld_coords", function(x,...)standardGeneric("ld_coords"))
setMethod("ld_coords", "Cycif",function(x,ld_name){
  if(missing(ld_name)){
    stop("'ld_name' should be specified to retrieve ld_coords")
  }else if(!ld_name %in% ld_names(x)){
    stop("Specified 'ld_name' doesn't exist")
  }
  return(x@ld_coords[[ld_name]])
})

#'@export
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
#' @rdname LDCoords
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
#' Run visualization with dimensionality reduction (t-SNE and U-MAP) on CyCIF data.
#'
#' @param x A Cycif or CycifStack object.
#' @param type character. It should be "raw" or "normalized", indicating which expression
#'  values to use for the U-MAP.
#' @param n.cells numeric. The number of cells sampled from each CyCIF object.
#' @param n_neighbors numeric. The number of neighbors, passed on to dimensionality reduction function.
#' @param used.abs A character vector containing a set of antibodies used for UMAP computation.
#' @param init.seed initial seed to set for computing UMAP.
#' @param smpls Character vector containign a set of samples to be included inthe UMAP.
#' @param ncells.per.smpl A numeric scholar. The number of cells per sample to be set when a CycifStack is run.
#' @param ... arguments passed to uwot::umap().
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
