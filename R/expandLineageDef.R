#'@export
expandLineageDef <- function(ctype,ctype.full=TRUE){
  if(!is.data.frame(ctype)){
    stop("ctype should be data.frame")
  }
  uniq.cts <- ctype$Child
  uniq.abs <- colnames(ctype)[-c(1:2)]
  is.str <- sapply(ctype[uniq.abs],function(x)any(x=="CAN"))
  if(all(is.str!="CAN")){
    chs1 <- ctype$Child
    cts.conv <- data.frame(idx=seq(chs1),original=chs1,expanded=chs1)
    return(list(ctype=ctype,names=cts.conv))
  }else if(all(is.str!="CAN") & !ctype.full){
    chs1 <- ctype$Child
    uniq.abs <- uniq.abs[!is.str]
    ctype <- ctype[,c("Parent","Child",uniq.abs)]
    cts.conv <- data.frame(idx=seq(chs1),original=chs1,expanded=chs1)
    return(list(ctype=ctype,names=cts.conv))
  }

  lin.abs <- uniq.abs[!is.str]
  str.abs <- uniq.abs[is.str]
  is.pd <- any(str.abs %in% c("PDL1","PD1"))
  str.abs <- c(str.abs[!is.pd],str.abs[is.pd])
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
  cts.conv <- data.frame(idx=idx,original=pas1,expanded=chs1)

  return(list(ctype=ct,names=cts.conv))
}
