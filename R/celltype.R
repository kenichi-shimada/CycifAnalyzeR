#_ -------------------------------------------------------

# utils CellType* ----

#' @export
setGeneric("ct_names", function(x)standardGeneric("ct_names"))

#' @rdname ct_names
#' @export
setMethod("ct_names", "Cycif", function(x)names(x@cell_types))

#' @rdname ct_names
#' @export
setMethod("ct_names", "CycifStack", function(x)names(x@cell_types))


#_ -------------------------------------------------------

# fun: getGates Cycif, CycifStack ----

#' @rdname gates
#' @export
setGeneric("getGates", function(x)standardGeneric("getGates"))

#' @rdname gates
#' @export
setMethod("getGates", "Cycif", function(x){
  gate.names <- paste0("gates_",names(x))
  if(!all(gate.names %in% names(abs_list(x)))){
    stop("Run 'setGates()' to insert gates in a CycifStack object")
  }else{
    out <- abs_list(x) %>% dplyr::select(any_of(c("ab",gate.names)))
  }
  return(out)
})

#' @rdname gates
#' @export
setMethod("getGates", "CycifStack", function(x){
  gate.names <- paste0("gates_",names(x))
  if(!all(gate.names %in% names(abs_list(x)))){
    stop("Run 'setGates()' to insert gates in a CycifStack object")
  }else{
    out <- abs_list(x) %>% dplyr::select(any_of(c("ab",gate.names)))
  }
  return(out)
})

# fun: setGates Cycif, CycifStack ----

#' @rdname gates
#' @export
setGeneric("setGates", function(x,...)standardGeneric("setGates"))

#' @rdname gates
#' @export
setMethod("setGates", "Cycif", function(x,gates.df,run_normalize=TRUE,p_thres=0.5,trim=1e-3){
  if(!is(gates.df,"data.frame")){
    "input should be in the data.frame"
  }

  if(all(rownames(gates.df) == as.character(seq(nrow(gates.df))))){
    if(!"ab" %in% names(gates.df)){
      stop("input obj should be a data.frame which contains two columns named 'ab', and 'gates_{smpl}'")
    }
  }else if(any(rownames(gates.df) %in% abs_list(x)$ab)){
    gates.df <- gates.df %>% tibble::rownames_to_column("ab")
  }

  this.col <- paste0("gates_",names(x))
  if(length(this.col)!=1 || !this.col %in% names(gates.df)){
    stop("the name of the column should be 'gates_*', where * are sample names")
  }else if(this.col %in% names(abs_list(x))){
    stop("Gates, ", this.col," already exists")
  }else{
    gates.df <- gates.df %>% dplyr::select(any_of(c("ab",this.col)))
  }

  x@abs_list <- abs_list(x) %>% dplyr::left_join(gates.df,by="ab")

  if(run_normalize){
    x <- normalize(x,method="logTh",p_thres=p_thres,trim=trim)
  }

  return(x)
})

#' @rdname gates
#' @export
setMethod("setGates", "CycifStack", function(x,gates.df,run_normalize=TRUE,p_thres=0.5,trim=1e-3){
  if(!is(gates.df,"data.frame")){
    "input should be in the data.frame"
  }

  if(all(rownames(gates.df) == as.character(seq(nrow(gates.df))))){
    if(!"ab" %in% names(gates.df)){
      stop("either rownames should be abs or one of the column's names should be 'ab'")
    }
  }else if(any(rownames(gates.df) %in% abs_list(x)$ab)){
    gates.df <- gates.df %>% tibble::rownames_to_column("ab")
  }

  this.col <- paste0("gates_",names(x))
  gates.df <- gates.df %>% dplyr::select(any_of(c("ab",this.col)))

  if(any(which(!this.col %in% names(gates.df)))){
    non.existing <- sub("gates_","",this.col[!this.col %in% names(gates.df)])
    stop("There are non-existing samples in the gates: ",non.existing)
  }


  x <- cyApply(x,function(cy){
    nm <- names(cy)
    cy <- setGates(cy,gates.df,run_normalize=run_normalize,p_thres=p_thres,trim=trim)
    return(cy)
  },as.CycifStack=TRUE)

  return(x)
})

#_ -------------------------------------------------------

# fun: expandLineageDef ctype,cstate ----

#' internally used only within CellTypeSkeleton, so don't export
#'
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

# fun: show CellTypes ----
#' @export
setMethod("show", "CellTypes", function(object){
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
  nst <- length(mst)

  # cat(m)
  cat("[",is(object)[[1]],"]\n",
      "Sample name:\t",object@name,"\n",
      "# cycles:\t",object@n_cycles,"\n\n",
      "# cell types:\t",nct,"\n",
      paste(cts,collapse=", "),"\n\n",

      "# markers in total:\t", nmk,"\n",
      "# cell lineage markers:",nlin,"\n",
      "# cell state markers:\t",nst,"\n\n",

      "lineage markers:",
      paste(nty$lin,collapse=", "),"\n",
      "state markers:",
      paste(mst,collapse=", ")
  )
})

#_ -------------------------------------------------------

# fun: CellTypeSkeleton Character,Cycif,CycifStack ----
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

# fun: CellTypeGraph ctype ----
#'@export
CellTypeGraph <- function(ctype,plot=F,transpose=T,...){
  uniq.cts <- c("all",ctype$Child)
  ctgraph <- ctype[c("Parent","Child")]
  ctgraph$Parent <- factor(ctgraph$Parent,levels=uniq.cts)
  ctgraph$Child <- factor(ctgraph$Child,levels=uniq.cts)
  g <- igraph::graph_from_data_frame(ctgraph)
  l <- igraph::layout_as_tree(g)
  levs <- l[,2]
  names(levs) <- names(igraph::V(g))
  ctlevs <- rev(tapply(names(levs),levs,identity))
  if(plot){
    if(transpose){
      l[,2] <- max(l[,2]) - l[,2]
      l <- l[,2:1]
    }
    igraph::V(g)$shape <- "rectangle"
    igraph::V(g)$label.family <- "Helvetica"

    plot(g, layout=l,
         edge.arrow.size=.5, vertex.color="gold", vertex.size=40,
         vertex.frame.color=NA, vertex.label.color="black",
         vertex.label.cex=0.8, vertex.label.dist=0, edge.curved=0
         ,...)
  }
  return(ctlevs)
}

#_ -------------------------------------------------------
# fun: defineCellTypes data.frame, Cycif, CycifStack ----
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
#' @param ct_anme name of the cell types
#' @param gates a data.frame containing gates (n.samples x n.proteins)
#' @param p_thres numerical between 0 and 1. A probability that corresponds to a threshold intensity
#'
#' @usage
#' defineCellTypes(x,ctype,cstate,gates,p_thres=0.5,...)
#'
#' @export
setMethod("defineCellTypes", "data.frame",
          function(x,
                   ctype,
                   cstate,
                   gates,
                   p_thres=0.5,...){
  # return a character vector containing 'cell_types'
  # cy <- x[[1]]
  if (!is(x,"data.frame")){
    stop("input should be a logTh_normalized expression (a data frame)")
  }
  if(nrow(x)==0) {
    stop("run normalize(method=\"logTh\") before CellTypeCalling()")
  }

  ctlevs <- CellTypeGraph(ctype,plot=F)

  cell.types <- rep("all",nrow(x))
  is.strict <- rep(TRUE,nrow(x))

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
          a <- apply(x[abs.and],1,min,na.rm=T)
          a[a==-Inf] <- NA
          b <- apply(x[abs.or],1,max,na.rm=T)
          b[b==-Inf] <- NA
          pos.out <- pmin(a,b)
        }else if(length(abs.and)>0){
          pos.out <- apply(x[abs.and],1,min,na.rm=T)
          pos.out[pos.out==Inf] <- NA
        }else if(length(abs.or)>0){
          pos.out <- apply(x[abs.or],1,max,na.rm=T)
          pos.out[pos.out==-Inf] <- NA
        }else{
          pos.out <- rep(1,nrow(x))
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

  return(data.frame(cell_types=cts,is_strict=is.strict))
})

#' @rdname defineCellTypes
#'
#' @param x A Cycif object.
#' @param ctype a data.frame containing cell type definition
#' @param cstate a data.frame containing cell state definition
#' @param ct_anme name of the cell types
#' @param gates a data.frame containing gates (n.samples x n.proteins)
#' @param p_thres numerical between 0 and 1. A probability that corresponds to a threshold intensity
#'
#' @usage
#' defineCellTypes(x,ctype,cstate,gates,p_thres=0.5,...)
#'
#' @export
setMethod("defineCellTypes", "Cycif",
          function(x,
                   ctype,
                   cstate,
                   ct_name=c("default","full"),
                   p_thres=0.5,
                   ctype.full=FALSE,
                   overwrite=FALSE,...){
            if(missing(ct_name)){
              if(ctype.full){
                ct_name <- "full"
              }else{
                ct_name <- "default"
              }
            }

            ## ct_name exists?
            if(ct_name %in% names(x@cell_types) && !overwrite){
              stop("cell type named '",ct_name,"' already exists and 'overwrite=FALSE'")
            }

            ## gates
            gates <- getGates(x)

            ## x@logTh_normalized exists?
            ex <- exprs(x,type="logTh_normalized")
            if(!is(ex,"data.frame") && nrow(ex)>0){
              stop('normalize(method="logTh") should run first')
            }

            # load ctype, cstate, gates in Cycif obj (both stratification markers unexpanded and expanded)
            ctc  <- CellTypeSkeleton(x,ctype=ctype,cstate=cstate,ctype.full=ctype.full)

            ctype <- ctc@cell_lineage_def
            cstate <- ctc@cell_state_def

            ## set plut info into defineCellTypes(df)
            cts <- defineCellTypes(ex,
                                   ctype=ctype,
                                   cstate=cstate,
                                   gates=gates,
                                   p_thres=p_thres)

            ctc@cell_types <- cts$cell_types
            ctc@is_strict <- cts$is_strict

            ctc@sample_names <- rep(names(x),nCells(x))

            x@cell_types[[ct_name]] <- ctc

            return(x)
          }
)

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
setMethod("defineCellTypes", "CycifStack",
  function(x,
           ctype,
           cstate,
           ct_name=c("default","full"),
           p_thres=0.5,
           ctype.full=FALSE,
           overwrite=FALSE,...){
    if(missing(ct_name)){
      if(ctype.full){
        ct_name <- "full"
      }else{
        ct_name <- "default"
      }
    }

    ## ct_name exists?
    if(ct_name %in% names(x@cell_types) && !overwrite){
      stop("cell type named '",ct_name,"' already exists and 'overwrite=FALSE'")
    }else{
      cat(paste0("Compute cell_types, and save the result under ct_name='",ct_name,"'\n"))
    }

    ## create celltypeskeleton
    ctc  <- CellTypeSkeleton(x,ctype=ctype,cstate=cstate,ctype.full=ctype.full)

    ctype <- ctc@cell_lineage_def
    cstate <- ctc@cell_state_def

    ## defineCellTypes for each sample (Cycif)
    for(nm in names(x)){
      cy <- x[[nm]]
      cat(paste0("Processing ",names(cy),"...\n"))
      cy <- defineCellTypes(cy,
                            ct_name=ct_name,
                            ctype=ctype,
                            cstate=cstate,
                            p_thres=p_thres,
                            ctype.full=ctype.full,
                            overwrite=overwrite)
      x[[nm]] <- cy
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

#_ -------------------------------------------------------

# fun: cell_types Cycif, CycifStack, CellTypeCycif, CellTypeCycifStack ----
#' @export
setGeneric("cell_types", function(x,...) standardGeneric("cell_types"))

#' @export
setMethod("cell_types", "CellTypes",
          function(x,
                   leaves.only=TRUE,
                   strict=FALSE,
                   uniq_cts){
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

  if(!missing(uniq_cts) && !is.null(uniq_cts)){
    # stop("uniq_cts not null")
    do.exist <- uniq_cts %in% levels(cts)
    if(!all(do.exist)){
      stop("non-existing cell_types specified in 'uniq_cts':",uniq_cts[!do.exist])
    }else{
      cts <- factor(cts,levels=uniq_cts)
    }
  }
  return(cts)
})

#' @export
setMethod("cell_types","Cycif",
          function(x,
                   ct_name="default",
                   leaves.only=TRUE,
                   strict=FALSE,
                   use_rois=TRUE,
                   uniq_cts=NULL){
  if(!ct_name %in% ct_names(x)){
    stop("ct_name doesn't exist: ",ct_name)
  }
  cts <- cell_types(x=x@cell_types[[ct_name]],
                    leaves.only=leaves.only,
                    strict=strict,
                    uniq_cts=uniq_cts)

  if(use_rois){
    is.rois <- x@within_rois
    cts[!is.rois] <- NA
  }

  smpls <- x@cell_types[[ct_name]]@sample_names
  df <- data.frame(sample=smpls,cell_types=cts)

  return(df)
})

#' @export
setMethod("cell_types", "CycifStack",
          function(x,
                   ct_name="default",
                   leaves.only=TRUE,
                   strict=FALSE,
                   use_rois=TRUE,
                   uniq_cts=NULL){
  if(!ct_name %in% ct_names(x)){
    stop("ct_name doesn't exist: ",ct_name)
  }
  cts <- cyApply(x,function(cy){
    cell_types(x=cy,
               ct_name=ct_name,
               leaves.only=leaves.only,
               strict=strict,
               use_rois=use_rois,
               uniq_cts=uniq_cts)
  })
  df <- do.call(rbind,cts)

  return(df)
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
  ct <-  x@cell_types
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
  ct <-  x@cell_types
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


