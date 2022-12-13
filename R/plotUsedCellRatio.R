#' @export
setGeneric("plotUsedCellRatio",function(x,...) standardGeneric("plotUsedCellRatio"))
setMethod("plotUsedCellRatio", "Cycif", function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                                             main="# cells attached on slide",...){
  x <- as.CycifStack(list(x))
  ret <- plotUsedCellRatio(x,cumulative=cumulative,ncycle=ncycle,mar=mar,main=main,...)
  return(invisible(ret))
})

#' @export
setMethod("plotUsedCellRatio", "CycifStack", function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                                                  main="# cells attached on slide",smpl.cols,ncol=1,
                                                  leg.cex=0.8,roi.selected=NULL,...){
  stopifnot(all(cyApply(x,class)=="Cycif"))
  stopifnot(all(unlist(cyApply(x,function(cy)nrow(cy@used_cells))>0)))
  ncycles <- unlist(cyApply(x,nCycles))
  if(missing(ncycle)){
    ncycle <- max(ncycles)
  }

  used.ratio <- data.frame(sapply(names(x),function(n){
    y <- x[[n]]
    rs <- roi.selected[[n]]
    nc.ratio <- statUsedCells(y,roi.selected=rs,cumulative=TRUE,ratio=TRUE)
    if(length(nc.ratio) > ncycle){
      nc.ratio <- nc.ratio[seq(ncycle)]
    }
    return(nc.ratio)
  }))

  smpls <- names(x)
  if(missing(smpl.cols)){
    smpl.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(nSamples(x))
  }

  ##
  par(mar=mar)
  plot(c(1,ncycle),c(0,1),type="n",
       xlab="# cycles",
       ylab="Relative # cells on slide",
       axes=F,main=main,...)
  box()
  axis(1,at=seq(ncycle),labels=seq(ncycle))
  axis(2)

  abline(h=seq(0.2,0.8,length=4),lty=1,col="grey90")
  abline(h=0:1,lty=1,col="grey80")
  abline(v=seq(ncycle),lty=1,col="grey90")

  for(i in seq(used.ratio)){
    ur <- used.ratio[[i]]
    lines(seq(ur),ur,col=smpl.cols[i],lwd=2)
    points(seq(ur),ur,col=smpl.cols[i],pch=20)
    points(seq(ur),ur,col=1,pch=1)
  }

  par(xpd=T)
  legend(par()$usr[2],par()$usr[4],smpls,lty=1,lwd=2,col=smpl.cols,ncol=ncol,cex=leg.cex)
  par(xpd=F)

  return(invisible(used.ratio))
})

##
