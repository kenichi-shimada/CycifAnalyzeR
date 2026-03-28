#_ -------------------------------------------------------

# fun: CellTypeSkeleton Character,Cycif,CycifStack ----

#' Generate CellTypes Object Based on Lineage and State Definitions
#'
#' @title Generate CellTypes Object
#'
#' @details
#' This function generates a \code{CellTypes} object based on provided lineage and state definitions. It accepts different input types:
#' - For a character vector `x` containing antibody names, `ctype` (cell lineage), and `cstate` (cell state) data frames, it creates a CellTypes object with lineage and state definitions.
#' - For a Cycif object `x`, it uses the antibodies from the object's abs_list to create a CellTypes object.
#' - For a CycifStack object `x`, it uses the antibodies from the object's abs_list to create a CellTypes object.
#'
#' @param x A character vector (for character input) or a Cycif/CycifStack object.
#' @param ctype A data frame defining cell lineage.
#' @param cstate A data frame defining cell state.
#' @param ctype.full Logical indicating whether to include the full lineage definition (default is FALSE).
#'
#' @return A CellTypes object containing cell lineage and state definitions.
#'
#' @seealso
#' \code{\link{CellTypes}}
#'
#' @rdname CellTypeSkeleton
#' @export
setGeneric("CellTypeSkeleton", function(x,...)standardGeneric("CellTypeSkeleton"))

#' @rdname CellTypeSkeleton
#' @export
setMethod("CellTypeSkeleton", "character",function(x,ctype,cstate,ctype.full=FALSE){
  if(missing(ctype) || missing(cstate)){
    stop("both lineage and state definitions should be provided")
  }
  if(!is(ctype,"data.frame")){
    stop("cell lineage definition should be a data.frame")
  }
  if(!is(cstate,"data.frame")){
    stop("cell state definition should be a data.frame")
  }
  if(!all(as.character(ctype$Child)==rownames(cstate))){
    stop("ctype$Child and cell types in cstate should be identical")
  }

  lmks <- names(ctype)[-c(1:2)]
  smks <- colnames(cstate)

  mks <- unique(c(lmks,smks))

  if(!any(mks %in% x)){
    stop("1st argument should be a character vector containing used abs")
  }

  if(!all(mks %in% x)){
    ## here ctype and cstate should be subsetted based on available mks
    ## Subsetting ctype and cstate so only used antibodies exist in the experiment
    used.abs1 <- lmks[lmks %in% x]
    unused.abs1 <- lmks[!lmks %in% x]
    used.ctype1 <- ctype[,used.abs1,drop=F]
    unused.ctype1 <- ctype[,unused.abs1,drop=F]
    is.used.ct <- !apply(unused.ctype1=="AND",1,any) #

    used.abs2 <- smks[smks %in% x]

    ctype <- ctype[is.used.ct,c("Parent","Child",used.abs1)]
    cstate <- cstate[is.used.ct,used.abs2]

  }

  elin <- expandLineageDef(ctype=ctype,cstate=cstate,ctype.full=ctype.full)

  ctype <- elin$ctype
  cstate <- elin$cstate

  lmks1 <- names(ctype)[-c(1:2)]
  smks1 <- colnames(cstate)

  mks1 <- unique(c(lmks1,smks1))

  mks.info <- data.frame(ab = mks1,
                         lineage = mks1 %in% lmks1,
                         state = mks1 %in% smks1)

  new("CellTypes",
      cell_lineage_def = elin$ctype,
      cell_state_def = elin$cstate,
      markers = mks.info
  )
})

#' @rdname CellTypeSkeleton
#' @export
setMethod("CellTypeSkeleton", "Cycif",
          function(x,ctype,cstate,ctype.full=FALSE){
  abs <- abs_list(x)$ab

  ## redefine ctype and cstate
  ctd <- CellTypeSkeleton(abs,ctype=ctype,cstate=cstate,ctype.full=ctype.full)

  new("CellTypes",
      sample_names = rep(names(x),nCells(x)),
      n_cycles = nCycles(x),
      cell_lineage_def = ctd@cell_lineage_def,
      cell_state_def = ctd@cell_state_def,
      markers = ctd@markers
  )
})

#' @rdname CellTypeSkeleton
#' @export
setMethod("CellTypeSkeleton", "CycifStack",function(x,ctype,cstate,ctype.full=FALSE){
  abs <- abs_list(x)$ab

  ## redefine ctype and cstate
  ctd <- CellTypeSkeleton(abs,ctype=ctype,cstate=cstate,ctype.full=ctype.full)

  new("CellTypes",
      sample_names = rep(names(x),nCells(x)),
      n_cycles = x@max_cycles,
      cell_lineage_def = ctd@cell_lineage_def,
      cell_state_def = ctd@cell_state_def,
      markers = ctd@markers
  )
})
#_ -------------------------------------------------------
# fun: defineCellTypes data.frame, Cycif, CycifStack ----

#' Define Cell Types and set them in a Cycif or CycifStack object
#'
#' This function performs a cell type calling operation and sets cell types in a Cycif or CycifStack object.
#'
#' @param x A Cycif or CycifStack object.
#' @param ctype A data.frame containing cell type definition.
#' @param cstate A data.frame containing cell state definition.
#' @param ct_name Name of the cell types.
#' @param p_thres Numerical value between 0 and 1. A probability that corresponds to a threshold intensity. Default is 0.5.
#' @param mc.cores Number of CPU cores to use for parallel processing. Default is 4.
#' @param overwrite Logical, if TRUE, overwrite existing cell type definitions with the same name. Default is FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The `defineCellTypes` function performs a cell type calling operation and sets cell types in a Cycif or CycifStack object. It takes as input a `ctype` data.frame, a `cstate` data.frame, and optional parameters for thresholding and parallel processing. Cell types can be defined at different hierarchy levels, and the results can be stored under a specified `ct_name`.
#'
#' If `ct_name` already exists and `overwrite` is set to FALSE, the function will not overwrite the existing cell type definitions and will issue a warning.
#'
#' @return A modified Cycif or CycifStack object with updated cell type information.
#'
#' @export
#'
#' @seealso \code{\link{CellTypeSkeleton}}, \code{\link{cyApply}}
setGeneric("defineCellTypes", function(x,...) standardGeneric("defineCellTypes"))

#' @rdname defineCellTypes
#' @export
setMethod("defineCellTypes", "data.frame",
          function(x,
                   ctype,
                   cstate,
                   p_thres=0.5,
                   mc.cores=4,
                   prioritized.celltypes=NULL,
                   ...){
# return a character vector containing 'cell_types'
  if (!is(x,"data.frame")){
    stop("input should be a logTh_normalized expression (a data frame)")
  }
  if(nrow(x)==0) {
    stop("run normalize(method=\"logTh\") before CellTypeCalling()")
  }

  ctlevs <- CellTypeGraph(ctype,plot=F)

  cell.types <- rep("all",nrow(x))

  is.strict <- rep(TRUE,nrow(x))  # return a character vector containing 'cell_types'

  ### Prioritize Cell Types ###
  if (!is.null(prioritized.celltypes)) {
#    stop(class(prioritized.celltypes))
    if(class(prioritized.celltypes) != "list"){
      stop("prioritized.celltypes should be a named list, whose names are cell types and values are lineage markers")
    }else if(!all(unlist(prioritized.celltypes) %in% names(x))){
      stop("prioritized.celltypes should contain only the markers that are in the expression data")
    }else if(!all(names(prioritized.celltypes) %in% unlist(ctlevs))){
      stop("prioritized.celltypes should contain only the cell types that are in the cell type definition")
    }

    ctlevs <- lapply(ctlevs,function(cts0){
      cts0[!cts0 %in% names(prioritized.celltypes)]
    })
    ctlevs <- ctlevs[sapply(ctlevs,length)>0]

    prioritized.prob <- parallel::mclapply(prioritized.celltypes, function(tmp) {
      if(length(unlist(tmp))==0){
        return(rep(NA,nrow(x)))
      }

      abs.and <- tmp$AND
      abs.or <- tmp$OR
      abs.not <- tmp$NOT

      suppressWarnings({
        if(length(abs.and)>0 & length(abs.or)>0){
          a <- apply(x[abs.and],1,min,na.rm=F)
          a[a==Inf] <- NA
          b <- apply(x[abs.or],1,max,na.rm=F)
          b[b==-Inf] <- NA
          pos.out <- pmin(a,b)
        }else if(length(abs.and)>0){
          pos.out <- apply(x[abs.and],1,min,na.rm=F)
          pos.out[pos.out==Inf] <- NA
        }else if(length(abs.or)>0){
          pos.out <- apply(x[abs.or],1,max,na.rm=F)
          pos.out[pos.out==-Inf] <- NA
        }else{
          # pos.out <- rep(p_thres,nrow(x)) # changed from 1 to p_thres
          pos.out <- rep(NA,nrow(x))
        }
        this.ct <- pos.out
        if(length(abs.not)>0){
          neg.out <- apply(x[abs.not],1,max,na.rm=T)
          neg.out[neg.out==-Inf] <- NA
          this.ct[which(neg.out > p_thres)] <- 0
          this.ct[is.na(neg.out)] <- NA
        }
      })
      return(this.ct)
    }, mc.cores=min(mc.cores, length(prioritized.celltypes)))

    prioritized.prob <- do.call(cbind, prioritized.prob)
    colnames(prioritized.prob) <- names(prioritized.celltypes)

    # is.used <- !apply(prioritized.prob,2,function(x)all(is.na(x)))
    # prioritized.prob <- prioritized.prob[,is.used,drop=F]

    cell.types0 <- apply(prioritized.prob,1,function(pr){
      ind <- which(pr > p_thres)
      if(length(ind)>0){
        return(names(pr[ind[1]]))
      }else{
        min(NA)
      }
    })
  } else{
    cell.types0 <- rep(NA,nrow(x))
  }

  ### End Prioritize Cell Types ###

  cat("ct_level=")
  for(l in seq(length(ctlevs)-1)){
    cat(l,"...",sep="")
    pas <- ctlevs[[l]]
    chs1 <- chs <- ctlevs[[l+1]]

    ct <- ctype %>% dplyr::filter(Parent %in% pas & Child %in% chs1)
    uniq.pas <- unique(ct$Parent)

    ## compute probability of being each child node in the following block
    mcc <- min(mc.cores,length(chs1))

    prs <- parallel::mclapply(chs1,function(ch){
    # prs <- sapply(chs1,function(ch){
      tmp <- unlist(ct %>% dplyr::filter(Child == ch))
      pa <- tmp[1]
      tmp <- tmp[-(1:2)]
      abs.and <- names(which(tmp=="AND"))
      abs.or <- names(which(tmp=="OR"))
      abs.not <- names(which(tmp=="NOT"))
      suppressWarnings({
        if(length(abs.and)>0 & length(abs.or)>0){
          a <- apply(x[abs.and],1,min,na.rm=F)
          a[a==Inf] <- NA
          b <- apply(x[abs.or],1,max,na.rm=T)
          b[b==-Inf] <- NA
          pos.out <- pmin(a,b)
        }else if(length(abs.and)>0){
          pos.out <- apply(x[abs.and],1,min,na.rm=F)
          pos.out[pos.out==Inf] <- NA
        }else if(length(abs.or)>0){
          pos.out <- apply(x[abs.or],1,max,na.rm=T)
          pos.out[pos.out==-Inf] <- NA
        }else{
          pos.out <- rep(p_thres,nrow(x)) # changed from 1 to p_thres
        }

        this.ct <- pos.out
        if(length(abs.not)>0){
          neg.out <- apply(x[abs.not],1,max,na.rm=T)
          neg.out[neg.out==-Inf] <- NA
          this.ct[which(!is.na(this.ct) & neg.out > p_thres)] <- 0
          this.ct[is.na(neg.out)] <- NA
        }
      })
      return(this.ct)
    # })
    },mc.cores=mcc)
    prs <- do.call(cbind,prs)
    colnames(prs) <- chs1

    ## warning("finish computing probs\n")
    ## for each parent, show which child is more likely to be the cell type each cell is
    for(pa in uniq.pas){
      # warning(pa,"\n")
      this.ind <- cell.types==pa
      this.chs <- (ct %>% dplyr::filter(Parent==pa))$Child
      this.chs <- colnames(prs)[colnames(prs) %in% this.chs]
      i.other <- grep("other",this.chs)
      if(!any(this.ind)){
        next
      }
      prs1 <- prs[this.ind,this.chs,drop=F]
      cell.types[this.ind] <- apply(prs1,1,function(pr){
        if(pa=="all"){
          pr1 <- pr[-i.other]
          if(all(is.na(pr1))){
            return("all_other")
          }
        }
        suppressWarnings({
          ind <- which(pr==max(pr,na.rm=T))
        })

        if(length(ind)>1){
          if(length(ind)==2 & length(i.other)==1 & any(ind==i.other)){
            ind <- ind[ind != i.other]
          }else{
            ind <- ind[1] # return("inc")
          }
        }
        if(length(ind)!=1){
          return(pa)
        }else{
          if(pr[ind] >= p_thres){
            return(this.chs[ind])
          }
        }
      })
      this.strict <- rowSums(prs1 > p_thres,na.rm=F) < 2
      this.strict[is.na(this.strict)] <- FALSE
      is.strict[this.ind] <- this.strict & is.strict[this.ind]
    }
    # warning("finish computing probs\n")
    # cat(paste0(sum(is.na(cell.types)),"..."))
  }

  ## convert to factor
  uniq.cts <- "all"
  for(l in seq(length(ctlevs)-1)){
    pas <- ctlevs[[l]]
    chs <- ctlevs[[l+1]]
    ct <- ctype %>% dplyr::filter(Parent %in% pas & Child %in% chs)
    tct <- tapply(ct$Child,ct$Parent,identity)
    for(i in seq(tct)){
      pa1 <- names(tct)[i]
      ch1 <- tct[[i]]
      uniq.cts <- unlist(replace(uniq.cts,which(uniq.cts==pa1),list(ch1)))
    }
  }
  if(!is.null(prioritized.celltypes)){
    uniq.cts <- c(names(prioritized.celltypes),uniq.cts)
  }

  leaves.cts <- ctype$Child[!ctype$Child %in% ctype$Parent]
  uniq.cts <- uniq.cts[uniq.cts %in% leaves.cts]
  # cts <- factor(cell.types,levels=c(uniq.cts,"unknown"))

  ## checking the prioritized.celltypes
  # table(cell.types,cell.types0))
  cell.types[!is.na(cell.types0)] <- cell.types0[!is.na(cell.types0)]
  cts <- factor(cell.types,levels=uniq.cts)
  # stop(cts)
  return(data.frame(cell_types=cts,is_strict=is.strict))
})

#' @rdname defineCellTypes
#' @export
setMethod("defineCellTypes", "Cycif",
          function(x,
                   ctype,
                   cstate,
                   ct_name="default",
                   p_thres=0.5,
                   mc.cores=4,
                   overwrite=FALSE,
                   prioritized.celltypes=NULL,
                   ...){

            if(missing(ct_name)){
              ct_name <- "default"
            }

            ## ct_name exists?
            if(ct_name %in% names(x@cell_types) && !overwrite){
              warning("cell type named '",ct_name,"' already exists and 'overwrite=FALSE'")
              return(x)
            }

            ## get within.rois
            within.rois <- within_rois(x)

            ## x@logTh_normalized exists?
            ex <- exprs(x,type="logTh")[within.rois,] ## expression only within rois
            if(!is(ex,"data.frame") && nrow(ex)>0){
              stop('normalize(method="logTh") should run first')
            }

            # load ctype, cstate, gates in Cycif obj (both stratification markers unexpanded and expanded)
            ctc  <- CellTypeSkeleton(x,ctype=ctype,cstate=cstate,ctype.full=FALSE)

            ctype <- ctc@cell_lineage_def
            cstate <- ctc@cell_state_def

            ## set plut info into defineCellTypes(df)
            nc <- nCells(x)
            cts1 <- defineCellTypes(ex,
                                   ctype=ctype,
                                   cstate=cstate,
                                   mc.cores=mc.cores,
                                   p_thres=p_thres,
                                   prioritized.celltypes=prioritized.celltypes)


            cts <- data.frame(
              cell_types = rep("outOfROI",nc),
              is_strict = rep(FALSE,nc)) %>%
              mutate(cell_types = factor(cell_types,levels=c(levels(cts1$cell_types),"outOfROI")))

            cts$cell_types[within.rois] <- cts1$cell_types
            cts$is_strict[within.rois] <- cts1$is_strict

            ctc@cell_types <- cts$cell_types
            ctc@is_strict <- cts$is_strict
            ctc@sample_names <- rep(names(x),nCells(x))

            x@cell_types[[ct_name]] <- ctc

            return(x)
          }
)

#' @rdname defineCellTypes
#' @export
setMethod("defineCellTypes", "CycifStack",
  function(x,
           ctype,
           cstate,
           ct_name="default",
           p_thres=0.5,
           mc.cores=4,
           overwrite=FALSE,
           prioritized.celltypes=NULL,
           ...){
    if(missing(ct_name)){
      ct_name <- "default"
    }

    ## ct_name exists?
    if(ct_name %in% names(x@cell_types) && !overwrite){
      stop("cell type named '",ct_name,"' already exists and 'overwrite=FALSE'")
    }else{
      cat(paste0("Compute cell_types, and save the result under ct_name='",ct_name,"'\n"))
    }

    ## create celltypeskeleton
    ctc  <- CellTypeSkeleton(x,ctype=ctype,cstate=cstate,ctype.full=FALSE)

    ctype <- ctc@cell_lineage_def
    cstate <- ctc@cell_state_def

    ## defineCellTypes for each sample (Cycif)
    for(nm in names(x)){
      cy <- x[[nm]]
      cat(paste0("Processing ",names(cy),"..."))

      ## define cell types
      cy <- defineCellTypes(cy,
                            ct_name=ct_name,
                            ctype=ctype,
                            cstate=cstate,
                            p_thres=p_thres,
                            mc.cores=mc.cores,
                            overwrite=overwrite,
                            prioritized.celltypes=prioritized.celltypes)

      x[[nm]] <- cy
      cat("done\n")
    }
    cat("Aggregating 'cell_types'...\n")
    ctc@cell_types <- unlist(cyApply(x,function(cy)cy@cell_types[[ct_name]]@cell_types))
    cat("Aggregating 'is_strict' flag...\n")
    ctc@is_strict <- unlist(cyApply(x,function(cy)cy@cell_types[[ct_name]]@is_strict))
    cat("Aggregating 'samples'...\n")
    ctc@sample_names <- unlist(cyApply(x,function(cy)cy@cell_types[[ct_name]]@sample_names))
    x@cell_types[[ct_name]] <- ctc

    return(x)
})
