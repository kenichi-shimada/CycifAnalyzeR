# ctype v2
def2num <- function(ctype){
  is.and <- ctype == "AND"
  is.can <- ctype == "CAN"
  num.can <- num.and <- array(0,dim(is.and))

  num.and[is.and] <- 1
  num.can[is.can] <- 1

  uniq.bins <- 2^(seq(ncol(ctype))-1)
  idx <- lapply(seq(nrow(ctype)),function(i){
    and <- (num.and %*% uniq.bins)[i]
    idx <- which(num.can[i,]==1)
    can.mat <- cbind(0,uniq.bins[idx])
    cans <- rowSums(expand.grid(lapply(seq(nrow(can.mat)),function(j){
      can.mat[j,]
    })))
    all <- and + cans
    return(all)
  })
  names(idx) <- rownames(ctype)
  vec.idx <- unlist(idx)
  cell.names <- rep(names(idx),sapply(idx,length))
  df <- data.frame(idx=vec.idx,cell.name=cell.names,stringsAsFactors=F)
  df <- rbind(c(0,"others"),df)
  rownames(df) <- c()
  return(df)
}

# v3 and later
defCellType <- function(cy,ctype,p_thres=0.5,strict=FALSE){
  norm <- exprs(cy,type="normalized")[colnames(ctype)]

  used.cells <- apply(norm,1,function(x)all(!is.na(x)))
  n.cells <- sum(used.cells)

  ctdef <- apply(ctype,1,function(x){
    thisdef <- x[!is.na(x)] =="AND"
  })
  if(0){
    used.abs <- used_abs(cy)
    is.def <- sapply(ctdef,function(x){
      all(names(x) %in% used.abs)
    })
  }
  def1 <- sapply(ctdef,function(def){
    apply(norm >= p_thres,1,function(x){
      all(x[names(def)] == def)
    })
  })

  cts <- rep(NA,nrow(def1))
  i0 <- which(rowSums(def1)==0)
  i1 <- which(rowSums(def1)==1)
  ct1 <- apply(def1[i1,],1,function(x)names(which(x)))
  cts[i0] <- "unknown"
  cts[i1] <- ct1

  inc <- which(rowSums(def1)>1)
  if(strict){
    cts[inc] <- "inconc"
    uniq.cts <- c(rownames(ctype),"unknown","inconc")
  }else{
    ct.abs <- unique(unlist(sapply(ctdef,function(x)names(which(x)))))
    inc.ct <- sapply(inc,function(i){
      tc <- this.ct <- ctdef[names(which(def1[i,]))]
      names(tc) <- c()
      this.mkrs <- names(which(unlist(tc)))
      this.mkrs <- this.mkrs[this.mkrs %in% ct.abs]
      max.mkr <- names(which.max(norm[i,this.mkrs]))
      ct1 <- names(this.ct)[which(sapply(this.ct,function(x)any(names(x)==max.mkr)))]
      return(ct1)
    })
    cts[inc] <- inc.ct
    uniq.cts <- c(rownames(ctype),"unknown")
  }

  cts <- factor(cts,levels=uniq.cts)
  is.def <- !is.na(rowSums(def1))
  names(cts) <- rownames(norm)
  return(cts)
}


