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
#' @param max.ncells.per.smpl A numeric scholar. The number of cells per sample to be set when a CycifStack is run.
#' @param ... arguments passed to uwot::umap().
#' @export
setGeneric("RunUMAP", function(x,...) standardGeneric("RunUMAP"))

#' @rdname RunUMAP
#' @export
setMethod("RunUMAP", "Cycif",
  function(x,type=c("raw","log_normalized","logTh_normalized"),
           ld_name,n.cells,used.abs,used.cts,ctype.full=FALSE,strict=TRUE,
           n_neighbors=20,init.seed=12345,...){
    call1 <- sys.calls()[[1]]
    if(missing(ld_name)){
      stop("'ld_name' should be specified (it's used to retrieve the data later)")
    }

    ## type - by default, should use normalized value
    if(missing(type)){
      type <- "logTh_normalized"
    }

    ## used.abs
    if(missing(used.abs)){
        stop("'used.abs' should be defined first (have you set threshold?).")
    }

    ## is.used
    cts <- cell_types(x,ctype.full=ctype.full,strict=strict)

    if(missing(used.cts)){
      used.cts <- levels(cts)
      used.cts <- used.cts[used.cts != "unknown"]
    }
    is.used <- cts %in% used.cts

    ## Select 'n.cells' cells from available.
    n.used <- sum(is.used)
    if(missing(n.cells)){
      n.cells <- n.used
    }

    smpl <- names(x)

    if(n.used < n.cells){
      warning(smpl, ": try sampling ",n.cells," cells but only ",n.used," cells available.\n")
    }else{
      cat(smpl, ": sampling ",n.cells," out of ",n.used," cells.\n")
    }

    set.seed(123)
    used.idx <- sample(which(is.used),n.cells)
    is.used.1 <- seq(is.used) %in% used.idx

    ## exprs matrix
    mat <- exprs(x,type=type)
    mat <- mat[is.used.1,used.abs]

    set.seed(init.seed)
    ru <- uwot::umap(mat,n_neighbors=n_neighbors,...)

    ru <- data.frame(ru)
    rownames(ru) <- which(is.used.1)
    names(ru) <- c("x","y")

    ld <- LDCoords(
      type = "UMAP",
      smpls = smpl,
      used.abs = used.abs,
      used.cts = used.cts,
      sn_cells_per_smpl = n.cells,
      n_cells_total = n.cells,
      ld_coords = ru,
      is_used = is.used.1,
      ld_params=call1)

    x@ld_coords[[ld_name]] <- ld
    return(x)
  })

#' @rdname RunUMAP
#' @export
setMethod("RunUMAP", "CycifStack",
  function(x,type=c("raw","log_normalized","logTh_normalized"),ld_name,
           smpls,used.abs,used.cts,ctype.full=FALSE,strict=TRUE,
           max.ncells.per.smpl,n_neighbors=20,init.seed=12345,
           save.coords=FALSE,...){
    call1 <- sys.calls()[[1]]
    if(missing(ld_name)){
      stop("'ld_name' should be specified (it's used to retrieve the data later)")
    }

    ## type - by default, should use normalized value
    if(missing(type)){
      type <- "logTh_normalized"
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
      x <- x[smpls] # redundant - no subsetting
    }else if(any(!smpls %in% names(x))){
      missing.smpls <- smpls[!smpls %in% names(x)]
      stop("missing smpls in 'smpls': ",paste(missing.smpls,collapse=", "))
    }

    ## celltypes
    cts <- cell_types(x,ctype.full=ctype.full,strict=strict)
    if(missing(used.cts)){
      used.cts <- levels(cts)
      used.cts <- used.cts[used.cts != "unknown"]
    }else if(!all(used.cts %in% levels(cts))){
      stop("all 'used.cts' should be observed cell types")
    }
    is.used <- cts %in% used.cts

    ## set 'n.cells'
    ncells <- nCells(x)
    v.ncells <- rep(names(ncells),ncells)
    n.used <- table(factor(is.used,levels=c("TRUE","FALSE")),v.ncells)["TRUE",]
    smpls <- names(which(n.used>0))

    ##
    idx.list.mat <- cyApply(x,function(y){
      smpl <- names(y)

      set.seed(123)
      is.used.1 <- is.used[v.ncells==smpl]
      if(max.ncells.per.smpl < sum(is.used.1)){
        used.idx <- sample(which(is.used.1),max.ncells.per.smpl)
        is.used.2 <- seq(is.used.1) %in% used.idx
      }else{
        is.used.2 <- is.used.1
      }
    },as.CycifStack=FALSE)

    is.used.2 <- unlist(idx.list.mat)
    mat <- exprs(x,type=type)[is.used.2,used.abs,drop=F]
    has.na <- apply(is.na(mat),1,any)
    mat <- mat[!has.na,]

    ##
    set.seed(init.seed)
    ru <- uwot::umap(mat,n_neighbors=n_neighbors,...)

    ru <- data.frame(ru)
    rownames(ru) <- rownames(mat)
    names(ru) <- c("x","y")

    tmp <- table(sub("\\..+","",rownames(mat)))
    n.used.cells <- as.vector(tmp)
    names(n.used.cells) <- names(tmp)


    ld <- LDCoords(
      type = "UMAP",
      smpls = smpls,
      used.abs = used.abs,
      used.cts = used.cts,
      n_cells_per_smpl = max.ncells.per.smpl,
      n_cells_total = n.used.cells,
      ld_coords = ru,
      is_used = is.used.2,
      ctype.full,ctype.full,
      ld_params = call1)

    x@ld_coords[[ld_name]] <- ld
    return(x)
    # validObject(x)
  })
