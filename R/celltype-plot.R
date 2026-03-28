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
