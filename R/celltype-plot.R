#_ -------------------------------------------------------

# fun: barplotPosCells Cycif, CycifStack ----

#' @title Create barplots to visualize positive cell ratios for each antibody in a Cycif object.
#'
#' @description This function generates barplots to visualize the ratios of positive cells for each antibody
#' in a Cycif object. You can choose between two types of barplots ('one' or 'two') and specify
#' additional graphical parameters.
#'
#' @param x A Cycif object.
#' @param type Type of barplot to create: 'one' (default) or 'two'.
#' @param mar A numeric vector of length 4 specifying the margin size for the plot.
#' @param ab For the CycifStack method, the antibody to plot.
#' @param ... Additional graphical parameters (unused).
#'
#' @details
#' The 'one' type of barplot displays the ratios of positive cells for each antibody separately.
#' It distinguishes between 'strict', 'non-strict', and 'cell-state' positive cells if applicable.
#' The 'two' type of barplot displays the pairwise overlap ratios of positive cells for each antibody.
#'
#' @return A barplot showing the ratios of positive cells for each antibody.
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @rdname barplotPosCells
#' @export
setGeneric("barplotPosCells", function(x,...) standardGeneric("barplotPosCells"))

#' @rdname barplotPosCells
#' @export
setMethod("barplotPosCells", "Cycif",function(x,type=c("one","two"),mar,...){
  .Defunct(msg=paste(
    "barplotPosCells() has had zero call sites anywhere across CycifAnalyzeR",
    "or any downstream project repo (EP, TALAVE, HR-APM)."
  ))
})
#' @rdname barplotPosCells
#' @export
setMethod("barplotPosCells", "CycifStack",function(x,ab,mar,...){
  .Defunct(msg=paste(
    "barplotPosCells() has had zero call sites anywhere across CycifAnalyzeR",
    "or any downstream project repo (EP, TALAVE, HR-APM)."
  ))
})

#_ -------------------------------------------------------

# fun: barplotCellTypes n_cts ----

#' Stacked barplot of cell-type composition per sample
#'
#' @param n.sh A matrix of cell-type frequencies, one column per sample (colnames matching `anno.smpls$id`).
#' @param anno.smpls A data.frame of sample annotations, with an `id` column matching `colnames(n.sh)`.
#' @param uniq.cols A vector of colors, one per cell type (row of `n.sh`).
#' @param ylim y-axis limits.
#' @param ylab y-axis label.
#'
#' @return Invisibly, the barplot midpoints.
#'
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
