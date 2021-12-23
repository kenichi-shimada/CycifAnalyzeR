setGeneric("segmentHist", function(x,...) standardGeneric("segmentHist"))
setMethod("segmentHist", "Cycif", function(x,ch="DNA0",cycle=0,...){
  ## used cells
  ab.cycle <- (x@abs_list %>% filter(ab == ch))$cycle
  this.cycle <- max(cycle,ab.cycle)

  if(nrow(x@used_cells)>0){
    stopifnot(nrow(x@used_cells)==nrow(x@raw))
    this.dna <- paste0("DNA",this.cycle)
    idx <- match(this.dna,names(x@used_cells))
    is.used <- apply(x@used_cells[seq(idx)],1,all)
  }else{
    is.used <- rep(TRUE,nrow(xy))
  }

  ## data
  if(grepl("^DNA[0-9]+$",ch)){
    int <- x@dna[[ch]]
  }else if(ch %in% x@abs_list$ab){
    int <- x@raw[[ch]]
  }else{
    stop("ab should be either an exact DNA or Ab name registered (see abs_list(x))")
  }
  int <- int[is.used]

  ##
  col.palette <- scaled_color(x,subset=is.used,ch=ch,palette=TRUE)
  a <- hist(int,breaks=1000)

  counts <- log10(a$counts+1)
  breaks <- a$breaks[-1]-diff(a$breaks)/2
  lo <- loess(counts ~ breaks,span=0.1)
  pred.lo <- predict(lo)
  ylim <- range(pred.lo)

  # par(mar=c(2,2,0,1))
  plot(breaks[c(1,seq(breaks),length(breaks))],c(0,predict(lo),0),
       type="n",axes=F,xlab="",ylab="",xlim=c(0,65536),ylim=ylim)
  axis(1)
  mtext("log(frequency)",2,0)
  polygon(breaks[c(1,seq(breaks),length(breaks))],
          c(ylim[1],predict(lo),ylim[1]),col="grey80")
  return(ylim)
})
