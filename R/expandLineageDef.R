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
