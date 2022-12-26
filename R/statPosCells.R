#' Get stats for positive cells after gating each abs
#' @export
setGeneric("statPosCells", function(x,...) standardGeneric("statPosCells"))

#' @export
setMethod("statPosCells", "Cycif",function(x,type=c("one","two"),mar,...){
  smpl <- names(x)
  lth <- exprs(x,type="logTh_normalized")            
  ct <-  x@cell_type
  used.abs <- c(names(ct@cell_lineage_def)[-(1:2)],names(ct@cell_state_def))

  pos <- as.data.frame(lapply(used.abs,function(ab){
    tmp <- lth[[ab]] > 0.5
    tmp <- factor(tmp,levels=c("TRUE","FALSE"))
    return(tmp)
  }))
  names(pos) <- used.abs
  
  opar <- par()
  if(type=="one"){
    diag <- sapply(pos,function(x)table(x)/sum(table(x)))
    if(missing(mar)){
      par(mar=c(7,5,4,1))      
    }else{
      par(mar=mar)
    }

    barplot(diag,beside=F,las=2,main=smpl,ylab="Ratio of positive cells",...)
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
    image(seq(nrow(rs)),seq(ncol(rs)),t(rs[rev(seq(ncol(rs))),]),col=NA,xlab="",ylab="",axes=F)
    axis(2,at=seq(nrow(rs)),labels=rev(rownames(rs)),las=1)
    axis(3,at=seq(nrow(rs)),labels=rownames(rs),las=2)
    par(fg="black")

    image(seq(nrow(rs)),seq(ncol(rs)),rs[,rev(seq(ncol(rs)))],col=uniq.cols,xlab="",ylab="",axes=F,add=T,
          main="",...)
    title(main=smpl,line=7)
    mtext("Numerator",side=3,line=5)
    mtext("Denominator",side=2,line=5)  
    suppressWarnings(par(opar))
    invisible(rs)
  }
})

        
  