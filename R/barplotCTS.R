#' @export
barplotCTS <- function(n.sh,uniq.cols,ylim=c(0,1)){
  nf <- layout(as.matrix(c(6:1)),widths=1,heights=c(1,1,1,1,1,16))
  ##
  par(mar=c(7,7,1,10))

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
                ylab="CellType Composition")
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
  legend(par()$usr[2],.8*ymax,levels(anno.smpls$TimePoint),fill=brewer.pal(3,"Set2"),title="TimePoint")
  legend(par()$usr[2],.55*ymax,c("MUT","WT"),fill=c("grey80","black"),title="gBRCA")
  legend(par()$usr[2],.375*ymax,levels(anno.smpls$Subtype),fill=brewer.pal(3,"Dark2"),title="Subtype")
  legend(par()$usr[2],0.2*ymax,sub("CD68_CD163_","DP_",rev(rownames(n.sh))[-1]),fill=rev(uniq.cols)[-1],title="cell types")
  par(xpd=F)
  par(fg=NA)
  axis(1,at=xidx,labels=names(med.x),las=2,cex.axis=.8)
  par(fg="black")

  d <- diff(xx[1:2])
  par(xpd=T)
  pts <- (pData(cs1) %>% filter(id %in% names(cs1)))$Patient.ID
  vs <- c(0,xx[cumsum(table(pts))] + d/2)
  sapply(vs,function(x1){
   lines(rep(x1,2),c(-0.05,1.05),lty=2)
  })
  text(xx[cumsum(table(pts))],1.15,unique(pts),srt=45,cex=.5)
  # legend(par()$usr[2],par()$usr[4],rev(names(uniq.cols)),fill=rev(uniq.cols),cex=.7,border=NA,box.col=NA)
  par(xpd=F)

  ##
  par(mar=c(0,7,.5,10))
  par(tcl=0)
  image(1:30,1,as.matrix(as.numeric(anno.smpls$BOR)),col=brewer.pal(3,"Set1")[c(1,3,2)],axes=F,
       xlab="",ylab="",xlim=c(-0.5,31.5))
  axis(2,at=1,labels=c("BOR"),las=1)
  ##
  par(mar=c(0,7,.5,10))
  image(1:30,1,as.matrix(as.numeric(anno.smpls$TimePoint)),col=brewer.pal(4,"Set2"),axes=F,
       xlab="",ylab="",xlim=c(-0.5,31.5))
  axis(2,at=1,labels=c("Time point"),las=1)
  ##
  par(mar=c(0,7,.5,10))
  image(1:30,1,as.matrix(as.numeric(anno.smpls$Subtype)),col=brewer.pal(3,"Dark2"),axes=F,
       xlab="",ylab="",xlim=c(-0.5,31.5))
  axis(2,at=1,labels=c("Subtype"),las=1)

  ##
  par(mar=c(0,7,.5,10))
  image(1:30,1,as.matrix(as.numeric(factor(anno.smpls$gBRCA.status))),col=c("black","grey80"),axes=F,
       xlab="",ylab="",xlim=c(-0.5,31.5))
  axis(2,at=1,labels=c("gBRCA.status"),las=1)

  ##
  par(mar=c(0,7,.5,10))
  par(tcl=0)
  image(1:30,1,as.matrix(as.numeric(anno.smpls$Patient)),col=brewer.pal(11,"Set3"),axes=F,
       xlab="",ylab="",xlim=c(-0.5,31.5))
  axis(2,at=1,labels=c("Patient"),las=1)
}
