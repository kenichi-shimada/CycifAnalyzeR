#' Plot 2D-coordinates computed by UMAP or T-SNE
#'
#' @param x A CycifStack object

#' @param type character, indicating how cells are colored. 'cell_type' is cell types
#'   computed prior to the plotting, 'smpl' show the samples the cells derived from,
#'   and 'exp' is expression of one antibody. When 'cell_type' and 'exp', additional
#'   parameters needs to be provided.
#' @param celltype a factor with the length of cell number, indicating cell types computed
#'   elsewhere. This needs to be provided when type is 'cell_type'.
#' @param ab character, indicating one channel/antibody. This needs to be provided when
#'   type is 'exp'.
#' @param uniq.cols color sets used for the plot()
#' @param pch pch passed to plot()
#' @param main main passed to plot()
#' @param ... other arguments  passed to plot()
#'
#' @export
setGeneric("lowDimPlot", function(x,...) standardGeneric("lowDimPlot"))


#' @rdname lowDimPlot
#' @export
setMethod("lowDimPlot", "Cycif",
  function(x,ld_name,type=c("celltype","cluster","exp"),
           pch=".",main,leg=TRUE,p_thres=p_thres,...){
    if(missing(ld_name)){
      stop("'ld_name' should be specified (it's used to retrieve the data later)")
    }
    ld <- ld_coords(x,ld_name=ld_name)
    xys <- ld@ld_coords
    is.used <- ld@is_used
    if(type=="celltype"){
      dc <- celltype
      levs <- levels(dc)
      levs <- levs[levs != "Others"]
      nlev <- length(levs)
      uniq.cols <- c(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nlev))
      cols <- uniq.cols[dc[rn]]

      if(missing(main)){
        main <- "Cell type"
      }

      par(mar=c(4,4,4,15))
      plot(ld,col=cols,pch=pch,main=main,xlab="",ylab="",...)

      if(leg){
        par(xpd=T)
        legend(par()$usr[2],par()$usr[4],levs,pch=pch,col=uniq.cols,cex=.9)
        par(xpd=F)
      }
    }else if(type=="exp"){
      stopifnot(!missing(ab))
      n <- x@normalized
      nn <- n[rn,ab]

      # rng <- range(n,na.rm=T)
      # rng <- range(zlim)
      # nn <- (n - rng[1])/diff(rng)

      nc <- 100
      uniq.cols <- rev(RColorBrewer::brewer.pal(11,"RdBu"))
      idx <-round(transform(seq(nc+1),method="Th",th=round((nc+1)/2),p_thres=p_thres)*nc)+1
      adj.uniq.cols <- colorRampPalette(uniq.cols)(nc+1)[idx]
      cols <- adj.uniq.cols[round(nn*nc)+1]

      if(missing(main)){
      #  main <- paste0(ab, "(expression)")
        main <- ab
      }

      def.par <- par(no.readonly = TRUE)
      # x <- layout(matrix(c(1,2),1,2,byrow=F),widths=c(3,1),heights = 1)
      # par(mar=c(4,4,4,0))
      par(mar=c(4,4,4,4))
      plot(ld,col=cols,pch=pch,main=main,...)
      # if(leg){
      #   par(mar=c(4,1,4,7))
      #   image(1,1:100,t(c(seq(1,50,length=40),seq(51,100,length=60))),
      #         col=colorRampPalette(uniq.cols)(100),
      #         axes=F,xlab="",ylab="")
      #   box()
      #   axis(4,at=c(1,40,100),labels=c("low","threshold","high"))
      #   par(def.par)
      # }
    }
  }
)

#' @rdname lowDimPlot
#' @export
setMethod("lowDimPlot", "CycifStack",
  function(x,pch=".",type=c("cell_type","smpl","exp"),
           type=c("logTh_normalized","log_normalized"),
           ab,main,leg=TRUE,p_thres=0.5,...){
    ld <- x@ld_coords
    rn <- rownames(ld)
    if(type=="cell_type"){
      dc <- factor(celltype[rn])
      stop(names(celltype))
      if(is.null(names(celltype))|all(rn %in% names(celltype))){
        stop("names(celltype) can't be blank")
      }

      levs <- levels(dc)
      nlev <- length(levs)
      uniq.cols <- c(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nlev))
      cols <- uniq.cols[dc[rn]]

      if(missing(main)){
        main <- "Cell type"
      }

      par(mar=c(4,4,4,13))
      plot(ld,col=cols,pch=pch,main=main,xlab="",ylab="",...)

      if(leg){
        par(xpd=T)
        legend(par()$usr[2],par()$usr[4],levs,pch=20,col=uniq.cols,,pt.cex=1)
        par(xpd=F)
      }

   }else if(type=="smpl"){
      smpls <- factor(sub("\\..+","",rn),levels=names(x))
      # set.seed(12345)
      # uniq.cols <- sample(colorRampPalette(brewer.pal(11,"Spectral"))(nlevels(smpls)))
      uniq.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(nlevels(smpls))
      cols <- uniq.cols[as.numeric(smpls)]
      pts <- levels(smpls)
      nlev <- length(pts)

      if(missing(main)){
        main <- "Samples"
      }

      par(mar=c(4,4,4,13))
      plot(ld,col=cols,pch=pch,main=main,xlab="",ylab="",...)

      if(leg){
        par(xpd=T)
        legend(par()$usr[2],par()$usr[4],pts,pch=20,col=uniq.cols,pt.cex=1)
        par(xpd=F)
      }
    }else if(type=="exp"){
      stopifnot(!missing(ab))
      n <- x@normalized
      mth <- x@normalize.method
      nn <- n[rn,ab]

      if(mth=="log"){
        rng <- range(nn,na.rm=T)
        nn <- (nn - rng[1])/diff(rng)
      }else if(mth=="logTh"){
        nn <- nn # it's already scaled between 0 and 1
      }

      nc <- 100
      uniq.cols <- rev(RColorBrewer::brewer.pal(11,"RdBu"))
      adj.uniq.cols <- colorRampPalette(uniq.cols)(nc+1)
      cols <- adj.uniq.cols[round(nn*nc)+1]

      if(missing(main)){
        #  main <- paste0(ab, "(expression)")
        main <- ab
      }

      def.par <- par(no.readonly = TRUE)
#      x <- layout(matrix(c(1,2),1,2,byrow=F),widths=c(3,1),heights = 1)
#      par(mar=c(4,4,4,0))
      par(mar=c(4,4,4,4))
      plot(ld,col=cols,pch=pch,main=main,...)
      # if(leg){
      #   par(mar=c(4,1,4,7))
      #   image(1,1:100,t(c(seq(1,50,length=40),seq(51,100,length=60))),
      #         col=colorRampPalette(uniq.cols)(100),
      #         axes=F,xlab="",ylab="")
      #   box()
      #   axis(4,at=c(1,40,100),labels=c("low","threshold","high"),las=1)
      # }
#      par(def.par)
    }
  }
)
