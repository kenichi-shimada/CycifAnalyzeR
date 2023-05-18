#_ -------------------------------------------------------

# utils CellType* ----
#' @rdname gates
#' @export
setMethod("gates", "CellTypeCycif", function(x)x@gates)
setMethod("gates", "CellTypeCycifStack", function(x)x@gates)

#_ -------------------------------------------------------

# fun: expandLineageDef ctype,cstate ----
#'@export
expandLineageDef <- function(ctype,cstate,ctype.full=TRUE){
  if(!is.data.frame(ctype)){
    stop("ctype should be data.frame")
  }
  if(!all(ctype$Child==rownames(cstate))){
    stop("Input data: ctype$Child and rownames(cstate) should be identical")
  }
  uniq.cts <- ctype$Child
  uniq.abs <- names(ctype)[-c(1:2)]
  is.str <- sapply(ctype[uniq.abs],function(x)any(x=="CAN")) # markers used to stratify cell_types
  if(!any(is.str)){
    return(list(ctype=ctype,cstate=cstate))
  }else if(any(is.str) && !ctype.full){
    chs1 <- ctype$Child
    abs1 <- uniq.abs[!is.str]
    ctype.sub <- ctype[,c("Parent","Child",abs1)]

    str.abs <- uniq.abs[!uniq.abs %in% abs1]
    cstate.ext <- cbind(ctype[str.abs],cstate)
    rownames(cstate.ext) <- ctype.sub$Child
    return(list(ctype=ctype.sub,cstate=cstate.ext))
  }else{# if(any(is.str) & ctype.full){
    lin.abs <- uniq.abs[!is.str]
    str.abs <- uniq.abs[is.str]

    ## remove stratification at non-leaf nodes
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    for(ab in str.abs){
      this.cts <- ctype$Child[which(ctype[[ab]]=="CAN")]
      non.leaf <- this.cts[!this.cts %in% leaves]
      if(length(non.leaf)>0){
        ctype[[ab]][ctype$Child %in% non.leaf] <- ""
      }
    }

    ## expand
    ct <- ctype
    pas <- ctype$Parent
    chs <- ctype$Child
    for(ab in str.abs){
      is.can <- (ct[[ab]]=="CAN")
      i <- rep(which(is.can),each=2)
      ct1 <- ct[i,]

      ## replace "AND","NOT","OR" in lin.abs in ct1
      tmp.lin <- ct1[lin.abs]
      for(lab in c("AND","OR","NOT")){
        tmp.lin[tmp.lin==lab] <- ""
      }
      ct1[lin.abs] <- tmp.lin

      ## update Parent and Child cell types
      ct1$Parent <- ct1$Child
      ct1$Child <- paste0(ct1$Child,",",ab,c("+","-"))
      ct1[[ab]] <- rep(c("AND","NOT"),sum(is.can))

      ## replace "CAN" with ""
      tmp.can <- ct[is.can,]
      tmp.can[tmp.can=="CAN"] <- ""
      ct[is.can,] <- tmp.can

      ct <- rbind(ct,ct1)
    }

    ## remove non-leaf nodes that are generated in the earlier process
    leaf1 <- ct$Child[!ct$Child %in% ct$Parent] ## keep
    non.leaf1 <- ct$Child[ct$Child %in% ct$Parent]
    non.leaf1 <- non.leaf1[!non.leaf1 %in% chs] ## remove non-leaf that did't exist originally
    while(length(non.leaf1)>0){
      nls.i <- which(ct$Child %in% leaf1 & ct$Parent %in% non.leaf1)
      tmp <- ct[nls.i,]
      nls <- tmp$Parent
      pa.i <- match(nls,ct$Child)
      nls.pa <- ct$Parent[pa.i]
      ct$Parent[nls.i] <- nls.pa
      ct <- ct[-pa.i,]

      non.leaf1 <- ct$Child[ct$Child %in% ct$Parent]
      non.leaf1 <- non.leaf1[!non.leaf1 %in% chs] ## remove non-leaf that did't exist originally
    }

    # CellTypeGraph(ct,plot=T,transpose=T)
    # CellTypeGraph(ctype,plot=T,transpose = T)
    pas1 <- ct$Parent
    chs1 <- ct$Child
    pas1[chs1 %in% chs] <- chs1[chs1 %in% chs]
    if(!all(pas1 %in% chs)){
      stop("not all the parents are previously seen children")
    }
    idx <- match(pas1,chs)
    cstate <- cstate[idx,]
    cstate.ext <- cbind(ct[str.abs],cstate)
    rownames(cstate.ext) <- chs1

    return(list(ctype=ct,cstate=cstate.ext))
  }
}

#_ -------------------------------------------------------

# fun: constructor CellTypeDefault Cycif,CycifStack ----
#' @export
CellTypeDefault <- function(x,ctype,cstate,ctype.full=TRUE){
  if(is(x,"Cycif") | is(x,"CycifStack")){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be either a Cycif or CycifStack object")
  }
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
  if(!all(mks %in% used.abs)){
    ## here ctype and cstate should be subsetted based on available mks
    ## Subsetting ctype and cstate so only used antibodies exist in the experiment
    used.abs1 <- lmks[lmks %in% abs$ab]
    unused.abs1 <- lmks[!lmks %in% abs$ab]
    used.ctype1 <- ctype[,used.abs1,drop=F]
    unused.ctype1 <- ctype[,unused.abs1,drop=F]
    is.used.ct <- !apply(unused.ctype1=="AND",1,any) #

    used.abs2 <- smks[smks %in% abs$ab]

    ctype <- ctype[is.used.ct,c("Parent","Child",used.abs1)]
    cstate <- cstate[is.used.ct,used.abs2]

    ## up to here
    unknown <- mks[!mks %in% used.abs]
  }
  mks.info <- abs %>% dplyr::filter(ab %in% mks)

  elin <- expandLineageDef(ctype=ctype,cstate=cstate,ctype.full=ctype.full)

  new("CellTypeDefault",
      cell_lineage_def = elin$ctype,
      cell_state_def = elin$cstate,
      markers = mks.info
  )
}

# fun: show CellTypeDefault ----
#' @export
setMethod("show", "CellTypeDefault", function(object){
  nmk <- length(object@markers$ab)
  cts <- object@cell_lineage_def$Child
  cts <- cts[!cts == "unknown"]
  nct <- length(cts)
  mty <- apply(object@cell_lineage_def[-c(1:2)],2,function(x){
    if(all(x %in% c("AND","OR","NOT","NOR",""))){
      return("lin")
    }else if(all(x %in% c("CAN",""))){
      return("str")
    }else{
      return("others")
    }
  })
  nty <- tapply(names(mty),mty,identity)
  mst <- colnames(object@cell_state_def)

  nlin <- length(nty$lin)
  nstr <- length(nty$str)
  nst <- length(mst)

  # cat(m)
  cat("[",is(object)[[1]],"]\n",
      "# cell types:\t",nct,"\n",
      paste(cts,collapse=", "),"\n\n",

      "# markers in total:\t", nmk,"\n",
      "# cell lineage markers:",nlin,"\n",
      "# stratifying markers:\t",nstr,"\n",
      "# cell state markers:\t",nst,"\n\n",
      "lineage markers:",
      paste(nty$lin,collapse=", "),"\n",
      "stratifying markers:",
      paste(nty$str,collapse=", "),"\n",
      "state markers:",
      paste(mst,collapse=", ")
  )
})

#_ -------------------------------------------------------

# fun: constructor CellTypeCycif ----
#'@export
CellTypeCycif <- function(x,ctype,cstate,gates.df,ctype.full=FALSE){
  if(is(x,"Cycif")){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a Cycif object")
  }

  if(is.matrix(gates.df)){
    gates.df <- as.data.frame(gates.df)
  }
  ## redefine ctype and cstate
  ctd <- CellTypeDefault(x,ctype,cstate,ctype.full=ctype.full)
  ctype <- ctd@cell_lineage_def
  cstate <- ctd@cell_state_def
  mks <- ctd@markers$ab

  ##
  g <- rep(NA,length(mks))
  names(g) <- mks

  smpl <- names(x)
  g[] <- gates.df[mks,smpl]

  abs <- abs %>% dplyr::filter(ab %in% mks)

  new("CellTypeCycif",
      name = x@name,
      n_cycles = x@n_cycles,
      cell_lineage_def = ctype,
      cell_state_def = cstate,
      markers = abs,
      gates = g
  )
}

# fun: show CellTypeCycif ----
#' @export
setMethod("show", "CellTypeCycif", function(object){
  nmk <- length(object@markers$ab)
  cts <- object@cell_lineage_def$Child
  cts <- cts[!cts == "unknown"]
  nct <- length(cts)
  mty <- apply(object@cell_lineage_def[-c(1:2)],2,function(x){
    if(all(x %in% c("AND","OR","NOT","NOR",""))){
      return("lin")
    }else if(all(x %in% c("CAN",""))){
      return("str")
    }else{
      return("others")
    }
  })
  nty <- tapply(names(mty),mty,identity)
  mst <- colnames(object@cell_state_def)

  nlin <- length(nty$lin)
  nstr <- length(nty$str)
  nst <- length(mst)

  # cat(m)
  cat("[",is(object)[[1]],"]\n",
      "Sample name:\t",object@name,"\n",
      "# cycles:\t",object@n_cycles,"\n\n",
      "# cell types:\t",nct,"\n",
      paste(cts,collapse=", "),"\n\n",

      "# markers in total:\t", nmk,"\n",
      "# cell lineage markers:",nlin,"\n",
      "# stratifying markers:\t",nstr,"\n",
      "# cell state markers:\t",nst,"\n\n",

      "lineage markers:",
      paste(nty$lin,collapse=", "),"\n",
      "stratifying markers:",
      paste(nty$str,collapse=", "),"\n",
      "state markers:",
      paste(mst,collapse=", ")
  )
})

#_ -------------------------------------------------------

# fun: constructor CellTypeCycifStack ----
#' @export
CellTypeCycifStack <- function(x,ctype.full=FALSE){
  if(is(x,"CycifStack")){
    abs <- abs_list(x)
    used.abs <- as.character(abs$ab)
  }else{
    stop("1st argument should be a CycifStack object")
  }

  nc <- nCycles(x)
  min.i <- min(which(nc==max(nc)))
  smpl <- names(x)[min.i]
  x1 <- x[[smpl]]

  if(ctype.full){
    ctc <- x1@cell_type_full
  }else{
    ctc <- x1@cell_type
  }

  ## load ctype and cstate
  ctype <- ctc@cell_lineage_def
  cstate <- ctc@cell_state_def

  ## marker abs
  lmks <- colnames(ctype)[-c(1:2)]
  smks <- colnames(cstate)
  mks <- c(lmks,smks)

  ## compile gates
  gates.df <- as.data.frame(cyApply(x,function(cy){
    if(ctype.full){
      ctc <- cy@cell_type_full
    }else{
      ctc <- cy@cell_type
    }
    ctc@gates[mks]
  },simplify=TRUE))

  gates.smpls <- names(gates.df)
  smpls <- names(x)
  if(!all(smpls %in% gates.smpls)){
    stop("Check gates; All the samples should be gated before cell type calling")
  }

  is.ungated <- sapply(gates.df,function(g)all(is.na(g)))

  if(any(is.ungated)){
    warning("some samples are not gated:",paste(smpls[is.ungated],collapse=","))
  }

  abs <- abs %>% dplyr::filter(ab %in% mks)

  new("CellTypeCycifStack",
      n_samples = x@n_samples,
      max_cycles = x@max_cycles,
      cell_lineage_def = ctype,
      cell_state_def = cstate,
      markers = abs,
      gates = gates.df
  )
}

# fun: show CellTypeCycifStack ----
#' @export
setMethod("show", "CellTypeCycifStack", function(object){
  nmk <- length(object@markers$ab)
  cts <- object@cell_lineage_def$Child
  cts <- cts[!cts == "unknown"]
  nct <- length(cts)
  mty <- apply(object@cell_lineage_def[-c(1:2)],2,function(x){
    if(all(x %in% c("AND","OR","NOT","NOR",""))){
      return("lin")
    }else if(all(x %in% c("CAN",""))){
      return("str")
    }else{
      return("others")
    }
  })
  nty <- tapply(names(mty),mty,identity)
  mst <- colnames(object@cell_state_def)

  nlin <- length(nty$lin)
  nstr <- length(nty$str)
  nst <- length(mst)

  is.gated <- sapply(object@gates,function(g){
    any(!is.na(g))
  })

  # cat(m)
  cat("[",is(object)[[1]],"]\n",
      "# samples:\t",object@n_samples,
      paste0("(",sum(is.gated)," gated)"),"\n",
      "# max cycles:\t",object@max_cycles,"\n\n",
      "# cell types:\t",nct,"\n",
      paste(cts,collapse=", "),"\n\n",

      "# markers in total:\t", nmk,"\n",
      "# cell lineage markers:",nlin,"\n",
      "# stratifying markers:\t",nstr,"\n",
      "# cell state markers:\t",nst,"\n\n",

      "lineage markers:",
      paste(nty$lin,collapse=", "),"\n",
      "stratifying markers:",
      paste(nty$str,collapse=", "),"\n",
      "state markers:",
      paste(mst,collapse=", ")
  )
})

#_ -------------------------------------------------------

# fun: cell_types CellTypeCycif, Cycif, CellTypeCycifStack, CycifStack ----
#' @export
setGeneric("cell_types", function(x,...) standardGeneric("cell_types"))

#' @export
setMethod("cell_types", "CellTypeCycif", function(x,leaves.only=TRUE,strict=FALSE){
  cts <- x@cell_types
  ctype <- x@cell_lineage_def
  if(leaves.only){
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    leaves <- c(leaves[!grepl("unknown",leaves)],"unknown")
    cts <- factor(cts,levels=leaves)
  }
  if(strict){
    is.strict <- x@is_strict
    cts[!is.strict] <- NA
  }
  return(cts)
})

#' @export
setMethod("cell_types", "CellTypeCycifStack", function(x,leaves.only=TRUE,strict=FALSE){
  cts <- x@cell_types
  ctype <- x@cell_lineage_def
  if(leaves.only){
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    leaves <- c(leaves[!grepl("unknown",leaves)],"unknown")
    cts <- factor(cts,levels=leaves)
  }
  if(strict){
    is.strict <- x@is_strict
    cts[!is.strict] <- NA
  }
  return(cts)
})

#' @export
setMethod("cell_types","Cycif",function(x,ctype.full=FALSE,leaves.only=TRUE,strict=FALSE,within.rois=TRUE){
  if(ctype.full){
    cts <- cell_types(x@cell_type_full,leaves.only=leaves.only,strict=strict)
  }else{
    cts <- cell_types(x@cell_type,leaves.only=leaves.only,strict=strict)
  }
  if(within.rois){
    is.rois <- x@within_rois
    cts[!is.rois] <- NA
  }
  return(cts)
})

#' @export
setMethod("cell_types", "CycifStack", function(x,ctype.full=FALSE,leaves.only=TRUE,strict=FALSE,within.rois=TRUE){
  if(ctype.full){
    ctd <- x@cell_type_full
  }else{
    ctd <- x@cell_type
  }
  cts <- ctd@cell_types

  if(leaves.only){
    ctype <- ctd@cell_lineage_def
    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    cts <- factor(cts,levels=leaves)
  }

  if(strict){
    is.strict <- ctd@is_strict
    cts[!is.strict] <- NA
  }

  if(within.rois){
    is.rois <- within_rois(x)
    cts[!is.rois] <- NA
  }
  return(cts)
})

#_ -------------------------------------------------------

# fun: CellTypeCalling Cycif ----
#' Cell type calling function
#' @param cy A cycif obj.
#' @param p_thres A numerical between 0 and 1. A fixed probability that corresponds to the threshold intensity.
#' @param ctype.full logical. if cell types with fuller definition (e.g., stratification with PD-1/PD-L1 status) should be used.
#'
#' @usage
#' CellTypeCalling(cy,p_thres=0.5,ctype.full=TRUE)
#' @export
CellTypeCalling <- function(cy,p_thres=0.5,ctype.full=TRUE){
  # return a character vector containing 'cell_types'
  # cy <- x[[1]]
  lth <- exprs(cy,type="logTh_normalized")
  if (nrow(lth)==0) {
    stop("run normalize(method=\"logTh\") before CellTypeCalling()")
  }

  if(ctype.full){
    ctc <- cy@cell_type_full
  }else{
    ctc <- cy@cell_type
  }
  ctype <- ctc@cell_lineage_def

  ctlevs <- CellTypeGraph(ctype,plot=F)

  cell.types <- rep("all",nrow(lth))
  is.strict <- rep(TRUE,nrow(lth))

  for(l in seq(length(ctlevs)-1)){
    pas <- ctlevs[[l]]
    chs <- ctlevs[[l+1]]
    i.others <- chs %in% "unknown" | grepl("_other$",chs)
    chs1 <- chs[!i.others]

    ct <- ctype %>% dplyr::filter(Parent %in% pas & Child %in% chs1)
    uniq.pas <- unique(ct$Parent)

    ## compute probability of being each child node in the following block
    ## (not considering the parent node)
    prs <- sapply(chs1,function(ch){
      tmp <- unlist(ct %>% dplyr::filter(Child == ch))
      pa <- tmp[1]
      tmp <- tmp[-(1:2)]
      abs.and <- names(which(tmp=="AND"))
      abs.or <- names(which(tmp=="OR"))
      abs.not <- names(which(tmp=="NOT"))
      suppressWarnings({
        if(length(abs.and)>0 & length(abs.or)>0){
          a <- apply(lth[abs.and],1,min,na.rm=T)
          a[a==-Inf] <- NA
          b <- apply(lth[abs.or],1,max,na.rm=T)
          b[b==-Inf] <- NA
          pos.out <- pmin(a,b)
        }else if(length(abs.and)>0){
          pos.out <- apply(lth[abs.and],1,min,na.rm=T)
          pos.out[pos.out==Inf] <- NA
        }else if(length(abs.or)>0){
          pos.out <- apply(lth[abs.or],1,max,na.rm=T)
          pos.out[pos.out==-Inf] <- NA
        }else{
          pos.out <- rep(1,nrow(lth))
        }

        this.ct <- pos.out
        if(length(abs.not)>0){
          neg.out <- apply(lth[abs.not],1,max,na.rm=T)
          neg.out[neg.out==-Inf] <- NA
          this.ct[which(neg.out > p_thres)] <- 0
          this.ct[is.na(neg.out)] <- NA
        }
      })
      return(this.ct)
    })

    ## for each parent, show which child is more likely to be the cell type each cell is
    for(pa in uniq.pas){
      this.ind <- cell.types==pa
      this.chs <- (ct %>% dplyr::filter(Parent==pa))$Child
      this.chs <- colnames(prs)[colnames(prs) %in% this.chs]
      prs1 <- prs[this.ind,this.chs,drop=F]
      cell.types[this.ind] <- apply(prs1,1,function(pr){
        if(all(is.na(pr))){
          return(pa)
        }
        ind <- which(pr==max(pr,na.rm=T))
        if(length(ind)>1){
          return(pa)
        }
        if(pr[ind] > p_thres){
          return(this.chs[ind])
        }else{
          ch1 <- paste0(pa,"_other")
          if(ch1 == "all_other"){
            ch1 <- "unknown"
          }
          return(ch1)
        }
      })
      this.strict <- rowSums(prs1 > p_thres,na.rm=F) < 2
      this.strict[is.na(this.strict)] <- FALSE
      is.strict[this.ind] <- this.strict & is.strict[this.ind]
    }
  }

  ## convert to factor: Q: how to order cell types from ctype df?
  uniq.cts <- c("all",ctype$Child)
  leaves <- uniq.cts[!uniq.cts %in% ctype$Parent]

  uniq.cts <- c(uniq.cts[uniq.cts!="unknown"],"unknown")
  cts <- factor(cell.types,levels=uniq.cts)

  return(list(cell_type=cts,is_strict=is.strict))
}

#_ -------------------------------------------------------

# fun: modifyCellTypes ctype ----
#'@export
modifyCellTypes <- function(ctype,uniq.cts,check=TRUE){
  pas <- ctype$Parent
  chs <- ctype$Child
  cts <- unique(c(pas,chs))

  if(!all(uniq.cts %in% cts)){
    non.exist <- uniq.cts[!uniq.cts %in% cts]
    stop("Some cell types are not defined in `ctype': ",paste0(non.exist,collapse=","))
  }

  ch.cts <- lapply(cts,function(x)x)
  names(ch.cts) <- cts
  ctlevs <- rev(CellTypeGraph(ctype,plot=F))

  for(lev in seq(ctlevs)){
    this.lev <- ctlevs[[lev]]
    ctype1 <- ctype %>% dplyr::filter(Child %in% this.lev) %>% dplyr::select(Parent,Child)
    lst1 <- tapply(ctype1$Child,ctype1$Parent,identity)
    for(pa in names(lst1)){
      lst1.ch <- unlist(ch.cts[c(lst1[[pa]],pa)])
      names(lst1.ch) <- c()
      ch.cts[[pa]] <- lst1.ch
    }
  }

  ##
  idx <- sapply(uniq.cts,function(nu){
    sapply(uniq.cts,function(de){
      is.under <- all(ch.cts[[nu]] %in% ch.cts[[de]])
    })
  })
  diag(idx) <- FALSE
  tab.idx <- which(idx,arr.ind=T)
  colnames(tab.idx) <- c("Parent","Child")
  rownames(tab.idx) <- c()
  hierarchy.cts <- as.data.frame(apply(tab.idx,2,function(i)uniq.cts[i]))
  if(any(idx)){
    cat("Provided uniq.cts include redundant cell types:\n")
    print(hierarchy.cts)
    stop("Provide non-redundant cell types.")
  }
  ch.cts1 <- ch.cts[uniq.cts]
  return(ch.cts1)
}

#_ -------------------------------------------------------

# fun: defineCellTypes Cycif, CycifStack ----
#' Define cell types
#'
#'Perform a cell type calling function and set cell types in a Cycif or CycifStack object
#'
#' @export
setGeneric("defineCellTypes", function(x,...) standardGeneric("defineCellTypes"))

#' @rdname defineCellTypes
#'
#' @param x A Cycif object.
#' @param ctype a data.frame containing cell type definition
#' @param cstate a data.frame containing cell state definition
#' @param gates a data.frame containing gates (n.samples x n.proteins)
#' @param p_thres numerical between 0 and 1. A probability that corresponds to a threshold intensity
#'
#' @usage
#' defineCellTypes(x,ctype,cstate,gates,p_thres=0.5,...)
#'
#' @export
setMethod("defineCellTypes", "Cycif", function(x,ctype,cstate,gates,p_thres=0.5,...){
  # load ctype, cstate, gates in Cycif obj (both stratification markers unexpanded and expanded)
  x@cell_type  <- CellTypeCycif(x,ctype,cstate,gates,ctype.full=FALSE)
  x@cell_type_full  <- CellTypeCycif(x,ctype,cstate,gates,ctype.full=TRUE)

  ## normalize - should be done within CellTypeCycif
  x <- normalize(x,method="logTh",p_thres=p_thres)

  ## set CellTypeCycif object in the cell_type slot
  cts.full <- CellTypeCalling(x,p_thres=p_thres,ctype.full=TRUE)
  cts.short <- CellTypeCalling(x,p_thres=p_thres,ctype.full=FALSE)

  x@cell_type_full@cell_types <- cts.full$cell_type
  x@cell_type_full@is_strict <- cts.full$is_strict
  x@cell_type@cell_types <- cts.short$cell_type
  x@cell_type@is_strict <- cts.short$is_strict

  return(x)
})

#'#' Define cell types for a CycifStack object
#'
#' @param x a CycifStack object
#' @param ... additional arguments (currently unused)
#'
#' @return a modified CycifStack object with updated cell type information
#'
#' @export
#' @rdname defineCellTypes
#' @method defineCellTypes CycifStack
setMethod("defineCellTypes", "CycifStack", function(x,...){
  # x <- cyApply(x,defineCellTypes,ctype=ctype,cstate=cstate,gates=gates,p_thres=p_thres)
  nct <- cyApply(x,function(y)length(y@cell_type@cell_types),simplify=T)
  if(any(nct == 0)){
    stop("run defineCellTypes() for each Cycif object")
  }
  x@cell_type <- CellTypeCycifStack(x,ctype.full=FALSE)
  x@cell_type_full <- CellTypeCycifStack(x,ctype.full=TRUE)

  x@cell_type@cell_types <- unlist(cyApply(x,function(cy)cy@cell_type@cell_types))
  x@cell_type@is_strict <- unlist(cyApply(x,function(cy)cy@cell_type@is_strict))
  x@cell_type_full@cell_types <- unlist(cyApply(x,function(cy)cy@cell_type_full@cell_types))
  x@cell_type_full@is_strict <- unlist(cyApply(x,function(cy)cy@cell_type_full@is_strict))

  return(x)
})



#_ -------------------------------------------------------

# fun: barplotPosCells Cycif, CycifStack ----

#' Get stats for positive cells after gating each abs
#' @export
setGeneric("barplotPosCells", function(x,...) standardGeneric("barplotPosCells"))

#' @rdname barplotPosCells
#' @importFrom RColorBrewer brewer.pal
#' @export
setMethod("barplotPosCells", "Cycif",function(x,type=c("one","two"),mar,...){
  smpl <- names(x)
  lth <- exprs(x,type="logTh_normalized")
  ct <-  x@cell_type
  used.abs <- c(names(ct@cell_lineage_def)[-(1:2)],names(ct@cell_state_def))
  lin.abs <- colnames(ct@cell_lineage_def)[-(1:2)]

  pos <- as.data.frame(lapply(used.abs,function(ab){
    tmp <- lth[[ab]] > 0.5
    tmp <- factor(tmp,levels=c("TRUE","FALSE"))
    return(tmp)
  }))
  names(pos) <- used.abs

  is.str <- ct@is_strict

  opar <- par()
  if(type=="one"){
    diag <- sapply(used.abs,function(ab){
      pos1 <- pos[[ab]]
      if(ab %in% lin.abs){
        no <- sum(pos1=="FALSE",na.rm=T) # false
        yes.relax <- sum(pos1=="TRUE" & !is.str, na.rm=T) # false
        yes.str <- sum(pos1=="TRUE" & is.str, na.rm=T) # false
        ntab <- c(0,yes.str,yes.relax,no)
        n <- ntab/sum(ntab)
      }else{
        ntab <- table(pos1)
        ntab <- c(ntab[1],0,0,ntab[2])
        n <- ntab/sum(ntab)
      }
      return(n)
    })

    if(missing(mar)){
      par(mar=c(7,5,4,7))
    }else{
      par(mar=mar)
    }

    barplot(diag,beside=F,las=2,col=c(3,2,4,"grey80"),main=smpl,border=FALSE,ylab="Ratio of positive cells",...)
    par(xpd=T)
    legend(par()$usr[2],par()$usr[4],c("strict","non-strict","cell-state"),fill=c(2,4,3))
    suppressWarnings(par(opar))
    invisible(diag)
  }else if(type=="two"){
    ## compute overlap of positive cells
    rs <- sapply(used.abs,function(ab1){
      sapply(used.abs,function(ab2){
        a <- pos[[ab1]]
        b <- pos[[ab2]]
        idx <- !is.na(a) & !is.na(b)
        n <- sum(a[idx]=="TRUE" & b[idx]=="TRUE",na.rm=T)
        d <- sum(a[idx]=="TRUE",na.rm=T)
        r <- n/d
        return(r)
      })
    })

    dimnames(rs) <- list(
      numerator=rownames(rs),
      denominator=colnames(rs))

    rs <- round(rs,2)

    ##  heatmap
    if(missing(mar)){
      par(mar=c(1,10,10,1))
    }else{
      par(mar=mar)
    }

    uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(50)


    par(fg="white")
    image(seq(nrow(rs)),seq(ncol(rs)),t(rs[rev(seq(ncol(rs))),]),col=NA,xlab="",ylab="",axes=F,zlim=c(0,1))
    axis(2,at=seq(nrow(rs)),labels=rev(rownames(rs)),las=1)
    axis(3,at=seq(nrow(rs)),labels=rownames(rs),las=2)
    par(fg="black")

    image(seq(nrow(rs)),seq(ncol(rs)),rs[,rev(seq(ncol(rs)))],col=uniq.cols,xlab="",ylab="",axes=F,add=T,
          main="",zlim=c(0,1),...)
    title(main=smpl,line=7)
    mtext("Numerator",side=3,line=5)
    mtext("Denominator",side=2,line=5)
    suppressWarnings(par(opar))
    invisible(rs)
  }
})

#' @export
setMethod("barplotPosCells", "CycifStack",function(x,ab,mar,...){
  ct <-  x@cell_type
  used.abs <- c(names(ct@cell_lineage_def)[-(1:2)],names(ct@cell_state_def))
  lin.abs <- names(ct@cell_lineage_def)[-(1:2)]

  if(missing(ab)){
    stop("ab should be specified")
  }else if(!ab %in% used.abs){
    stop("ab should be one of the used abs")
  }

  lth <- exprs(x,type="logTh_normalized")
  sms <- sub("\\..+","",rownames(lth))

  tmp <- lth[[ab]] > 0.5
  lab <- rep(NA,length(tmp))
  if(ab %in% lin.abs){
    is.str <- ct@is_strict

    lab[!tmp] <- "no"
    lab[tmp & !is.str] <- "pos.relax"
    lab[tmp & is.str] <- "pos.str"

    lab <- factor(lab, levels=c("pos.str","pos.relax","no"))

    leg <- c("strict","non-strict","negative")
    cols <- c(2,4,"grey80")
    mar <- c(7,5,7,7)
  }else{
    lab[!tmp] <- "neg"
    lab[tmp] <- "pos"
    lab <- factor(lab,levels=c("pos","neg"))

    leg <- c("positive","negative")
    cols <- c(3,"grey80")
    mar <- c(7,5,7,7)
  }
  ntab <- table(lab,sms)
  diag <- apply(ntab,2,function(tmp)tmp/sum(tmp))

  opar <- par()

  par(mar=mar)

  xs <- barplot(diag,beside=F,las=2,main="",ylab="Ratio of positive cells",col=cols,...)
  title(main=ab,line=4)

  par(xpd=T)
  legend(par()$usr[2],par()$usr[4],leg,fill=cols)
  par(xpd=F)

  # suppressWarnings(par(opar))
  invisible(xs)

})

#_ -------------------------------------------------------

# fun: barplotCellTypes n_cts ----
#' @export
barplotCellTypes <- function(n.sh,anno.smpls,uniq.cols,ylim=c(0,1),ylab="CellType Composition"){
  nf <- layout(as.matrix(c(6:1,7)),widths=1,heights=c(1,1,1,1,1,1,16))
  if(!all(colnames(n.sh)==anno.smpls$id)){
    stop("colnames(n.sh) and anno.smpls$id should be the sample names with identical order.")
  }

  ##
  xs <- seq(ncol(n.sh))
  range.xs <- range(xs) + c(-1,1)*1.75

  ## 1: plot BOR
  par(mar=c(0,7,.5,10))
  par(tcl=0)
  image(xs,1,as.matrix(as.numeric(anno.smpls$BOR)),col=brewer.pal(3,"Set1")[c(1,3,2)],axes=F,
        xlab="",ylab="",xlim=range.xs)
  axis(2,at=1,labels=c("BOR"),las=1)

  ## 2: plot PFS
  par(mar=c(0,7,.5,10))
  pfs <- anno.smpls$PFS_days
  pfs.idx <- pfs - min(pfs,na.rm=T) + 1
  pfs.uniq.cols <- colorRampPalette(brewer.pal(9,"Blues"))(max(pfs.idx,na.rm=T))
  pfs.col <- c("grey80",pfs.uniq.cols[unique(sort(pfs.idx))])
  pfs.idx[is.na(pfs.idx)] <- -1
  pfs.idx <- as.numeric(factor(pfs.idx))

  par(tcl=0)
  image(xs,1,as.matrix(pfs.idx),col=pfs.col,axes=F,
        xlab="",ylab="",xlim=range.xs)
  axis(2,at=1,labels=c("PFS"),las=1)

  ## 3: plot Time point
  par(mar=c(0,7,.5,10))
  image(xs,1,as.matrix(as.numeric(anno.smpls$TimePoint)),col=brewer.pal(4,"Set2"),axes=F,
        xlab="",ylab="",xlim=range.xs)
  axis(2,at=1,labels=c("Time point"),las=1)

  ## 4: plot Subtype
  par(mar=c(0,7,.5,10))
  image(xs,1,as.matrix(as.numeric(anno.smpls$Subtype)),col=brewer.pal(3,"Dark2")[1:2],axes=F,
        xlab="",ylab="",xlim=range.xs)
  axis(2,at=1,labels=c("Subtype"),las=1)

  ## 5: plot gBRCA
  par(mar=c(0,7,.5,10))
  image(xs,1,as.matrix(as.numeric(factor(anno.smpls$gBRCA.status))),col=c("grey20","grey80"),axes=F,
        xlab="",ylab="",xlim=range.xs)
  axis(2,at=1,labels=c("gBRCA.status"),las=1)

  ## 6: plot Patient
  par(mar=c(0,7,.5,10))
  par(tcl=0)
  image(xs,1,as.matrix(as.numeric(anno.smpls$Patient)),col=c(brewer.pal(12,"Set3"),brewer.pal(8,"Pastel1")),
        axes=F,
        xlab="",ylab="",xlim=range.xs)
  axis(2,at=1,labels=c("Patient"),las=1)

  ## 7: main barplot
  par(mar=c(7,7,1,10))
  par(tcl=-0.5)

  ymax <- ylim[2]
  if(ymax!=1){
    if(!"others" %in% rownames(n.sh)){
      stop("Unless 'others' is in cell types, ymax should be 1")
    }else{
      n.sh["others",] <- n.sh["others",] - (1-ymax)
    }
  }

  xx <- barplot(n.sh,beside=F,border=F,col=uniq.cols,las=2,
                main="",names.arg=rep("",ncol(n.sh)),
                ylab=ylab)
  med.x <- tapply(seq(nrow(sub.pd)),
                  factor(sub.pd$Patient.ID,levels=unique(sub.pd$Patient.ID)),median)
  xidx <- sapply(med.x,function(x){
    if(x %% 1 == 0.5){
      median(xx[floor(x)+ (0:1)])
    }else{
      xx[x]
    }
  })

  par(xpd=T)
  legend(par()$usr[2],par()$usr[4],levels(anno.smpls$BOR),fill=brewer.pal(3,"Set1")[c(1,3,2)],title="BOR")
  legend(par()$usr[2],.8*ymax,levels(anno.smpls$TimePoint),fill=brewer.pal(4,"Set2"),title="TimePoint")
  legend(par()$usr[2],.55*ymax,c("MUT","WT"),fill=c("grey20","grey80"),title="gBRCA")
  legend(par()$usr[2],.375*ymax,levels(anno.smpls$Subtype),fill=brewer.pal(3,"Dark2"),title="Subtype")
  legend(par()$usr[2],0.2*ymax,sub("CD68_CD163_","DP_",rev(rownames(n.sh))),fill=rev(uniq.cols),title="Cell types")
  par(xpd=F)
  par(fg=NA)
  axis(1,at=xidx,labels=names(med.x),las=2,cex.axis=.8)
  par(fg="black")

  d <- diff(xx[1:2])
  par(xpd=T)
  pts <- factor(sub.pd$Patient.ID,levels=unique(sub.pd$Patient.ID))
  vs <- c(xx[1]-d/2,xx[cumsum(table(pts))] + d/2)
  sapply(vs,function(x1){
    lines(rep(x1,2),c(-0.05,1.4),lty=2)
  })
  par(xpd=F)

}


