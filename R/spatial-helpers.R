#_ -------------------------------
# fun: unused (non-exported) ----
dist2tgt <- function(this.cts, xy, seed.tum=1,n.tumors=100,summary=TRUE){
  if(sum(this.cts=="Tumor",na.rm=T)==0){
    stop("some of this.cts should be 'Tumor'")
  }

  has.cts <- !is.na(this.cts)
  uniq.cts <- levels(this.cts)

  lst.xys <- lapply(uniq.cts,function(ct)xy[!is.na(this.cts) & this.cts==ct,])
  names(lst.xys) <- uniq.cts

  ## sampling tumor cellsa
  set.seed(seed.tum)
  n.tumors.all <- nrow(lst.xys$Tumor)
  n.tumors <- min(n.tumors,n.tumors.all)
  idx1 <- sample(n.tumors.all,n.tumors,replace=FALSE)

  xys.tumors <- lst.xys$Tumor[idx1,]

  tgts <- names(lst.xys)
  tgts <- tgts[tgts != "Tumor"]
  dist2tumor <- sapply(tgts,function(tgt){
    if(nrow(lst.xys[[tgt]]) == 0){
      return(c(nn.dist=Inf))
    }
    lst.gds <- RANN::nn2(data=lst.xys[[tgt]],query=xys.tumors,k=1,searchtype="standard") #distance from every tumor to cd8t
    if(summary){
      med.d2tum <- median(lst.gds$nn.dists)
      return(med.d2tum)
    }else{
      return(lst.gds$nn.dists)
    }
  })
  return(dist2tumor)
}

# fun: stat.dist2tgt (immature) ----

stat.dist2tgt <- function(this.cts,xy,
                          seed.tum=1,n.tumors=100,summary=TRUE,
                          n.iter=1e4,seed.iter=123){

  i1 <- which(this.cts!="outOfROI")
  this.cts <- this.cts[i1]
  xy <- xy[i1,]

  ## distance - original
  cts.orig <- this.cts

  orig.d <- dist2tgt(cts.orig,xy,n.tumors=n.tumors,seed.tum=seed.tum,summary=summary)

  ## sampling
  ds <- do.call(rbind,parallel::mclapply(seq(n.iter),function(i){
    set.seed(i)
    cts <- sample(this.cts)
    ds <- dist2tgt(cts,xy,n.tumors=n.tumors,seed.tum=seed.tum,summary=summary)
    return(ds)
  },mc.cores=15))

  med.ds <- apply(ds,2,median)

  folds <- orig.d/med.ds
  return(folds)
}
