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
  function(x,ld_name,plot.type=c("celltype","cluster","exp"),
           na.col = "grey80",uniq.cols,
           pch=".",main,leg=TRUE,p_thres=0.5,...){
    if(missing(plot.type)){
      stop("Which plot type? (plot.type = celltype, cluster, exp)")
    }
    
    if(missing(ld_name)){
      stop("'ld_name' should be specified (it's used to retrieve the data later)")
    }else if(!ld_name %in% ld_names(x)){
      stop("Specified 'ld_name' does not exist.")
    }else{
      ld <- ld_coords(x,ld_name=ld_name)
    }

    xys <- ld@ld_coords
    ld.type <- ld@ld_type
    norm.type <- ld@norm_type
    
    is.used <- ld@is_used
    used.abs <- ld@used.abs
    used.cts <- ld@used.cts
    cts.params <- ld@cts_params
    
    if(plot.type=="celltype"){
      dc <- cell_types(x,
                       ctype.full = cts.params$ctype.full,
                       leaves.only=cts.params$leaves.only,
                       strict=cts.params$strict)
      dc <- factor(dc,levels=used.cts)
      levs <- levels(dc)
      nlev <- length(levs)
      if(missing(uniq.cols)){
        if(nlev < 9){
          uniq.cols <- RColorBrewer::brewer.pal(max(3,nlev),"Set1")
        }else{
          uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nlev)
        }
        names(uniq.cols) <- used.cts
      }
      cols <- rep(na.col,length(dc))
      cols[!is.na(dc)] <- uniq.cols[as.character(dc[!is.na(dc)])]
      cols <- cols[is.used]
      
      if(missing(main)){
        main <- paste0(ld.type,", Cell types")
      }

      par(mar=c(4,4,4,15))
      plot(xys,col=cols,pch=pch,main=main,xlab="",ylab="",...)

      if(leg){
        par(xpd=T)
        legend(par()$usr[2],par()$usr[4],levs,pch=pch,col=uniq.cols,cex=.9)
        par(xpd=F)
      }
    }else if(plot.type=="exp"){
      if(missing(ab)){
        stop("'ab' should be specified when plot.type == 'exp'")
      }
      n <- exprs(x,type=norm.type)[is.used,ab] 

      nc <- 100
      uniq.cols <- rev(RColorBrewer::brewer.pal(11,"RdBu"))
      
      idx <-round(transform(seq(nc+1),method="Th",th=round((nc+1)/2),p_thres=p_thres,trim=0)*nc)+1
      adj.uniq.cols <- colorRampPalette(uniq.cols)(nc+1)[idx]
      cols <- adj.uniq.cols[round(n*nc)+1]

      if(missing(main)){
        main <- paste0(ld.type,", ",ld_name,", ",ab)
      }

      def.par <- par(no.readonly = TRUE)
      par(mar=c(4,4,4,4))
      plot(xys,col=cols,pch=pch,main=main,...)
      # if(leg){
      #   par(mar=c(4,1,4,7))
      #   image(1,1:100,t(c(seq(1,50,length=40),seq(51,100,length=60))),
      #         col=colorRampPalette(uniq.cols)(100),
      #         axes=F,xlab="",ylab="")
      #   box()
      #   axis(4,at=c(1,40,100),labels=c("low","threshold","high"))
      #   par(def.par)
      # }
    }else if(plot.type=="cluster"){
      ## TO BE DONE
    }
  }
)

#' @rdname lowDimPlot
#' @export
setMethod("lowDimPlot", "CycifStack",
  function(x,ld_name,plot.type=c("celltype","smpl","cluster","exp"),
           norm_type=c("logTh_normalized","log_normalized"),pch=".",
           uniq.cols,ab,main,leg=TRUE,p_thres=0.5,...){
    if(missing(plot.type)){
      stop("Which plot type? (plot.type = celltype, smpl, cluster, exp)")
    }
    if(missing(ld_name)){
      stop("'ld_name' should be specified (it's used to retrieve the data later)")
    }else if(!ld_name %in% ld_names(x)){
      stop("Specified 'ld_name' does not exist.")
    }else{
      ld <- ld_coords(x,ld_name=ld_name)
    }
    xys <- ld@ld_coords
    ld.type <- ld@ld_type
    norm.type <- ld@norm_type
    
    is.used <- ld@is_used
    used.abs <- ld@used.abs
    used.cts <- ld@used.cts
    cts.params <- ld@cts_params
    
    if(plot.type=="celltype"){
      dc <- cell_types(x,
                       ctype.full = cts.params$ctype.full,
                       leaves.only=cts.params$leaves.only,
                       strict=cts.params$strict)
      dc <- factor(dc,levels=used.cts)
      levs <- levels(dc)
      nlev <- length(levs)
      if(missing(uniq.cols)){
        if(nlev < 9){
          uniq.cols <- RColorBrewer::brewer.pal(max(3,nlev),"Set1")
        }else{
          uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nlev)
        }
        names(uniq.cols) <- used.cts
      }
      cols <- rep(na.col,length(dc))
      cols[!is.na(dc)] <- uniq.cols[as.character(dc[!is.na(dc)])]
      cols <- cols[is.used]
      
      if(missing(main)){
        main <- paste0(ld.type,", Cell types")
      }
      
      par(mar=c(4,4,4,15))
      plot(xys,col=cols,pch=pch,main=main,xlab="",ylab="",...)
      if(leg){
        par(xpd=T)
        legend(par()$usr[2],par()$usr[4],levs,pch=pch,col=uniq.cols,cex=.9)
        par(xpd=F)
      }
   }else if(plot.type=="smpl"){
      n.cells.smpls <- ld@n_cells_total
      uniq.smpls <- names(n.cells.smpls)
        
      smpls <- factor(rep(uniq.smpls,n.cells.smpls),levels=uniq.smpls)
      # set.seed(12345)
      # uniq.cols <- sample(colorRampPalette(brewer.pal(11,"Spectral"))(nlevels(smpls)))
      uniq.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(nlevels(smpls))
      cols <- uniq.cols[as.numeric(smpls)]
      pts <- levels(smpls)
      nlev <- length(pts)

      if(missing(main)){
        main <- paste0(ld.type,", ",ld_name,", Samples")
      }

      par(mar=c(4,4,4,13))
      plot(xys,col=cols,pch=pch,main=main,xlab="",ylab="",...)

      if(leg){
        par(xpd=T)
        legend(par()$usr[2],par()$usr[4],pts,pch=20,col=uniq.cols,pt.cex=1)
        par(xpd=F)
      }
    }else if(plot.type=="exp"){
      if(missing(ab)){
        stop("'ab' should be specified when plot.type == 'exp'")
      }
      n <- exprs(x,type=norm.type)[is.used,ab] 
      
      nc <- 100
      uniq.cols <- rev(RColorBrewer::brewer.pal(11,"RdBu"))
      
      idx <-round(transform(seq(nc+1),method="Th",th=round((nc+1)/2),p_thres=p_thres,trim=0)*nc)+1
      adj.uniq.cols <- colorRampPalette(uniq.cols)(nc+1)[idx]
      cols <- adj.uniq.cols[round(n*nc)+1]
      
      if(missing(main)){
        main <- paste0(ld.type,", ",ld_name,", ",ab)
      }
      
      def.par <- par(no.readonly = TRUE)
      par(mar=c(4,4,4,4))
      plot(xys,col=cols,pch=pch,main=main)
    }else if(plot.type=="cluster"){
      ## TO BE DONE
    }
  }
)
