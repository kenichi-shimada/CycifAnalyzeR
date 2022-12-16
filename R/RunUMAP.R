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
#' @param n.cells.per.smpl A numeric scholar. The number of cells per sample to be set when a CycifStack is run.
#' @param ... arguments passed to uwot::umap().
#' @export
setGeneric("RunUMAP", function(x,...) standardGeneric("RunUMAP"))

#' @rdname RunUMAP
#' @export
setMethod("RunUMAP", "Cycif",
  function(x,type=c("raw","normalized"),n.cells,
           celltype,
           n_neighbors=20,used.abs,init.seed=12345,...){

    ## type - by default, should use normalized value
    if(missing(type)){
      type <- "normalized"
    }

    ## used.abs
    if(missing(used.abs)){
        stop("'used.abs' should be defined first (have you set threshold?).")
    }

    ## how many cycles to be used
    ab.cycle <- max((abs_list(x) %>% filter(ab %in% used.abs))$cycle) + 1 # note we use non-zero origin

    ## is.used
    if(!missing(celltype)){
      idx.used <- which(celltype != "Others")
      is.used <- rep(FALSE,length(celltype))
      is.used[idx.used] <- TRUE
    }else if(!.hasSlot(x,"used_cells")){
      stop("'used_cells' should be defined first.\n")
    }else{
      is.used <- cumUsedCells(x)[,ab.cycle]
    }

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

    x@ld_coords <- ru

    # validObject(x)
    return(x)
  })

#' @rdname RunUMAP
#' @export
setMethod("RunUMAP", "CycifStack",
          function(x,type=c("raw","normalized"),smpls,used.abs,
                   celltype,removed_ct="",
                   n.cells.per.smpl,n_neighbors=20,init.seed=12345,...){

            ## type - by default, should use normalized value
            if(missing(type)){
              type <- "normalized"
            }

            ## samples
            if(missing(smpls)){
              smpls <- names(x)
            }
            if(type=="normalized"){
              is.normalized <- all(cyApply(x,function(y).hasSlot(y,type),simplify=TRUE))
              if(!is.normalized){
                stop("Some samples aren't normalized yet.\n")
              }
            }

            x <- x[smpls] # redundant - no subsetting

            ## used.abs
            if(missing(used.abs)){
              list.used.abs <- cyApply(x,function(y)abs_list(y))
              tab.used.abs <- table(unlist(list.used.abs))
              used.abs <- names(which(tab.used.abs == length(smpls)))
            }

            ## how many cycles to be used
            ab.cycle <- max((abs_list(x) %>% filter(ab %in% used.abs))$cycle) + 1 # note we use non-zero origin

            ## is.used
            if(!missing(celltype)){
              list.is.used <- lapply(celltype,function(ct){
                  is.used <- !is.na(ct) & !(ct %in% removed_ct)
                  return(is.used)
              })
            }else{
              list.is.used <- cyApply(x,function(y){
                if(!.hasSlot(y,"used_cells")){
                  stop("'used_cells' should be defined first.\n")
                }
                is.used <- cumUsedCells(y)[,ab.cycle]
              })
            }

            n.used <- sapply(list.is.used,sum)

            ## Select 'n.cells' cells from available.
            if(missing(n.cells.per.smpl)){
              n.cells.per.smpl <- min(n.used)
            }

            ##
            if(min(n.used) < n.cells.per.smpl){
              warning("Tried to use ",n.cells.per.smpl," cells but some samples only have ",min(n.used)," cells.\n")
              n.cells.per.smpl <- min(n.used)
            }

            list.mat <- cyApply(x,function(y){
              smpl <- names(y)

              set.seed(123)
              is.used <- list.is.used[[smpl]]
              used.idx <- sample(which(is.used),n.cells.per.smpl) #
              is.used.1 <- seq(is.used) %in% used.idx

              mat <- exprs(y,type=type)
              mat <- mat[is.used.1,used.abs]
              rownames(mat) <- which(is.used.1)

              return(mat)
            },as.CycifStack=FALSE)

            mat <- do.call(rbind,list.mat)

            if(0){
              ## consistency with cts
              rn <- rownames(mat)
              idx <- which(dc[rn]=="Tumor")
              mat[head(idx),]
            }

            ##
            set.seed(init.seed)
            ru <- uwot::umap(mat,n_neighbors=n_neighbors,...)

            ru <- data.frame(ru)
            rownames(ru) <- rownames(mat)
            names(ru) <- c("x","y")

            x@ld_coords <- ru

            # validObject(x)
            return(x)
          })
