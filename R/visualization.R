#_ -------------------------------------------------------

# fun: AbsSummary CycifStack ----
#' Split a string
#'
#' @param x An object of Cycif class.
#' @param show.cycles.in.row logical. FALSE by default. If TRUE, the output plot shows the number
#'     of cycles each sample was stained for.
#' @param ... Arguments passed to image() function.
#'
#' @export
AbsSummary <- function(x,show.cycles.in.row=FALSE,...){
  if(!is(x,"CycifStack")){
    stop("input should be a CycifStack obj.")
  }
  uniq.abs <- as.character(x@abs_list$ab)
  n1 <- do.call(rbind,lapply(x@samples,function(y){
    this.abs <- as.character(abs_list(y)$ab)
    tested <- uniq.abs %in% this.abs
    return(tested)
  }))

  idx <- rep(seq(ncol(n1)/3),each=3)
  n2 <- n1 * idx[col(n1)]
  n2[n2==0] <- NA
  pcols <- grDevices::grey(c(.8,.6))[(seq(max(n2,na.rm=T)) %% 2) + 1]

  if(show.cycles.in.row){
    smpl.labs <- paste0(x@names,"\n(",x@n_cycles, " cycles)")
  }else{
    smpl.labs <- x@names
  }
  col.labs <- paste0("Cycle\n",seq(x@max_cycles))
  graphics::image(seq(ncol(n2)),seq(nrow(n2)),t(n2[rev(seq(nrow(n2))),]),col=pcols,
                  xlab="",ylab="",axes=F,...)
  graphics::box()
  graphics::abline(h=c(0,seq(nrow(n1)))+.5,lwd=.5)
  graphics::abline(v=c(0,seq(ncol(n1)))+.5,lwd=.5)
  graphics::axis(2,at=rev(seq(nrow(n1))),labels=smpl.labs,las=1)
  graphics::axis(1,at=seq(ncol(n2)),labels=uniq.abs,las=2)
  graphics::par(tcl=0)
  graphics::axis(3,at=seq(x@max_cycles)*3-1,labels=col.labs)
  graphics::par(tcl=-.5)

  return(invisible(NULL))
}

#_ -------------------------------------------------------

# fun: slidePlot Cycif ----
#' @export
setGeneric("slidePlot", function(x,...) standardGeneric("slidePlot"))
setMethod("slidePlot", "Cycif",
          function(x,pch=20,cex=2,plot_type=c("dna","exp","cell_type","custom"),
                   custom_labs,ctype.full=FALSE,strict=FALSE,ttl,ab,
                   uniq.cts,uniq.cols,draw.roi=TRUE,
                   remove.unknown=TRUE,cell.order,
                   na.col="grey80",use.roi=TRUE,use.thres=TRUE,ncells=1e4,
                   contour=FALSE,cont_nlevs=3,
                   trim_th=1e-2,legend=FALSE, legend.pos="bottomright",mar=c(3,3,3,3),
                   roi.sq,...){
            if(missing(plot_type)){
              stop("need to specify 'plot_type' argument: dna, exp, cell_type, custom")
            }

            smpl <- names(x)

            if(!missing(cell.order)){
              if(length(cell.order) != nCells(x)){
                stop("'cell.order' should be the same size as nrow(x@dna)")
              }
            }

            if(plot_type=="dna"){
              if(missing(ab) || !ab %in% paste0("DNA",seq(nCycles(x)))){
                stop(paste0("If DNA stain, `ab' should be one of DNA1, ..., DNA", nCycles(x),")"))
              }
              mat <- x@dna

              mat <- cbind(log1p(mat[[1]]),as.data.frame(lapply(mat,function(x)log1p(x/mat[[1]])))[-1])
              names(mat) <- names(x@dna)

              n.ab <- trim_fun(mat[[ab]],trim_th=trim_th)
              rn <- range(n.ab,na.rm=T)

              is.na <- is.na(n.ab)
              if(use.roi){
                within.rois <- x@within_rois
                if(length(within.rois)){
                  stop("ROIs not defined yet; 'use.roi' should be FALSE")
                }
                is.na <- is.na | !within.rois
              }

              if(missing(uniq.cols)){
                nlev <- 50
                idx <- round((n.ab[!is.na] - rn[1])/diff(rn)*(nlev-1))+1
                uniq.cols <- grDevices::colorRampPalette(rev(brewer.pal(11,"Spectral")))(nlev)
                cols <- rep(na.col,length(n.ab))
                cols[!is.na] <- uniq.cols[idx]
              }else{
                idx <- n.ab[!is.na]
                cols <- rep(na.col,length(n.ab))
                cols[!is.na] <- uniq.cols[idx]
              }

              if(missing(ttl)){
                if(ab=="DNA1"){
                  ab1 <- ab
                }else{
                  ab1 <- paste0(ab,"/DNA1")
                }
                ttl <- paste0(smpl," (",ab1,")")
              }
            }else if(plot_type=="exp"){
              cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=use.roi,strict=strict)
              if(missing(uniq.cts)){
                uniq.cts <- levels(cts)
                if(any(uniq.cts=="unknown")){
                  if(remove.unknown){
                    uniq.cts <- uniq.cts[uniq.cts!="unknown"]
                  }else{
                    ui <- which(uniq.cts=="unknown")
                    uniq.cts <- c(uniq.cts[-ui],"unknown")
                  }
                }
              }
              cts <- factor(cts,levels=uniq.cts)

              if(missing(ab) || !ab %in% abs_list(x)$ab){
                stop(ab, " is not specified or available in the sample ", names(x))
              }

              n <- exprs(x,plot_type="log_normalized")

              n.ab <- trim_fun(n[[ab]],trim_th=trim_th)
              rn <- range(n.ab,na.rm=T)
              is.na <- is.na(n.ab)
              if(use.roi){
                within.rois <- x@within_rois
                is.na <- is.na | !within.rois
              }

              if(missing(uniq.cols)){
                nlev <- 50

                idx <- round((n.ab[!is.na] - rn[1])/diff(rn)*(nlev-1))+1
                uniq.cols <- grDevices::colorRampPalette(rev(brewer.pal(11,"Spectral")))(nlev)
                cols <- rep(na.col,length(n.ab))
                cols[!is.na] <- uniq.cols[idx]
              }else{
                idx <- n.ab[!is.na]
                cols <- rep(na.col,length(n.ab))
                cols[!is.na] <- uniq.cols[idx]
              }

              if(missing(ttl)){
                ttl <- paste0(smpl,", ",ab," expression")
              }
            }else if(plot_type=="cell_type"){
              cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=use.roi,strict=strict)

              if(missing(ttl)){
                if(missing(uniq.cts)){
                  ttl <- paste0(smpl,", cell-types")
                }else{
                  ttl <- paste0(smpl,", ",paste0(uniq.cts,collapse=","))
                }
              }
              if(missing(uniq.cts)){
                uniq.cts <- levels(cts)

                if(any(uniq.cts=="Immune_other")){
                  ui <- which(uniq.cts=="Immune_other")
                  uniq.cts <- c(uniq.cts[-ui],"Immune_other")
                }

                if(any(uniq.cts=="unknown")){
                  if(remove.unknown){
                    uniq.cts <- uniq.cts[uniq.cts!="unknown"]
                  }else{
                    ui <- which(uniq.cts=="unknown")
                    uniq.cts <- c(uniq.cts[-ui],"unknown")
                  }
                }
              }
              cts <- factor(cts,levels=uniq.cts)

              if(missing(uniq.cols)){
                nct <- nlevels(cts)
                if(nct>10){
                  uniq.cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[-6])(nct)
                }else{
                  # set.seed(12)
                  uniq.cols <- (RColorBrewer::brewer.pal(11,"Spectral"))[c(7,1,4,11:9,8)]
                }
                names(uniq.cols) <- uniq.cts
              }
              cols <- uniq.cols[cts]
              is.na <- is.na(cts)
              if(use.roi){
                within.rois <- x@within_rois
                is.na <- is.na | !within.rois
              }


            }else if(plot_type=="custom"){
              if(missing(custom_labs)){
                stop("when plot_type='custom', the argument 'custom_labs' should be specified.")
              }
              if(!is(custom_labs,"factor")){
                custom_labs <- factor(custom_labs)
              }
              if(missing(uniq.cols)){
                nct <- nlevels(custom_labs)
                if(nct>11){
                  uniq.cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nct)
                }else{
                  uniq.cols <- RColorBrewer::brewer.pal(nct,"Spectral")
                }
                names(uniq.cols) <- levels(custom_labs)
              }
              cols <- uniq.cols[custom_labs]
              is.na <- is.na(custom_labs)
              if(use.roi){
                within.rois <- x@within_rois
                is.na <- is.na | !within.rois
              }

              if(missing(cell.order)){
                cell.order <- order(custom_labs,decreasing=T)
              }
              if(length(pch)==1){
                pch <- rep(pch,length(custom_labs))
              }
              if(length(cex)==1){
                cex <- rep(cex,length(custom_labs))
              }
              if(missing(ttl)){
                ttl <- ""
              }
            }

            ## plot
            xy <- xys(x)

            max.y <- max(xy$Y)
            xy$Y_centroid <- max.y - xy$Y

            if(!is.na(ncells) && nrow(xy) > ncells){
              set.seed(123)
              is.used <- seq(nrow(xy)) %in% sort(sample(nrow(xy),ncells))
            }else{
              is.used <- rep(TRUE,nrow(xy))
            }

            omar <- graphics::par()$mar
            ##
            if(!missing(cell.order)){
              xy <- xy[cell.order,]
              if(exists("is.na") && length(is.na)==length(cell.order)){
                # stop("is.na")
                is.na <- is.na[cell.order]
              }
              if(exists("is.used") && length(is.used)==length(cell.order)){
                is.used <- is.used[cell.order]
              }
              if(exists("pch") && length(pch)==length(cell.order)){
                pch <- pch[cell.order]
              }
              if(exists("cts") && length(cts)==length(cell.order)){
                cts <- cts[cell.order]
              }
              if(exists("cols") && length(cols)==length(cell.order)){
                cols <- cols[cell.order]
              }
            }

            if(!missing(roi.sq)){
              roi.sq$y <- range(max.y - roi.sq$y)
              is.na <- is.na |
                xy$X < roi.sq$x[1] | xy$X > roi.sq$x[2] |
                xy$Y < roi.sq$y[1] | xy$Y > roi.sq$y[2]
              xlim <- roi.sq$x
              ylim <- roi.sq$y
            }else{
              xlim <- range(xy$X)
              ylim <- range(xy$Y)
            }
            graphics::par(mar=mar)
            plot(xy$X,xy$Y,main=ttl,asp=1,xlab="",ylab="",type="n",xlim=xlim,ylim=ylim,...)

            cex1 = 4.8 * 8/3 * (graphics::par()$pin[1]/graphics::par()$cin[2])/diff(graphics::par()$usr[1:2])

            graphics::points(xy$X[is.na & is.used],xy$Y[is.na & is.used],col=na.col,pch=pch,cex=cex1)

            if(plot_type=="cell_type"){
              graphics::points(xy$X[!is.na & cts=="unknown" & is.used],
                     xy$Y[!is.na & cts=="unknown" & is.used],
                     col=cols[!is.na & cts=="unknown" & is.used],
                     pch=pch,cex=cex1*cex)
              graphics::points(xy$X[!is.na & cts!="unknown" & is.used],
                     xy$Y[!is.na & cts!="unknown" & is.used],
                     col=cols[!is.na & cts!="unknown" & is.used],
                     pch=pch,cex=cex1*cex)
            }else if(plot_type=="custom"){
              graphics::points(xy$X,
                     xy$Y,
                     col=cols,
                     pch=pch,cex=cex1*cex)
            }else{
              graphics::points(xy$X[!is.na & is.used],
                     xy$Y[!is.na & is.used],
                     col=cols[!is.na & is.used],
                     pch=pch,cex=cex1*cex)
            }
            if(draw.roi){
              prs <- x@rois
              # ncys <- sapply(prs,function(x)x$cycle) - guaranteed to be relevant to this ncycle
              for(i in seq(prs)){
                pr <- prs[[i]]
                # graphics::points(pr,pch=sub(".+(.)$","\\1",as.character(seq(length(pr$x)))))
                if(pr$dir=="positive"){
                  col1 <- 2
                }else if(pr$dir=="negative"){
                  col1 <- 4
                }
                graphics::polygon(pr$coords,lty=1,border=col1)
              }
            }
            if(plot_type=="cell_type" && contour){
              this.idx <- !is.na & cts %in% uniq.cts & is.used
              f1 <- MASS::kde2d(xy$X[this.idx],xy$Y[this.idx], n = 1000, lims = graphics::par()$usr,
                          h = c(width.SJ(xy$X), width.SJ(xy$Y)))
              contour(f1, nlevels  = cont_nlevs, add=T, labels=rep("",3),col="grey50")
            }
            if(legend){
              legend(legend.pos,names(uniq.cols),col=uniq.cols,pch=20)
            }
            graphics::par(mar=omar)
          }
)

trim_fun <- function(x,trim_th = 1e-3){
  qts <- stats::quantile(x,c(trim_th,1-trim_th),na.rm=T)
  x[x < qts[1]] <- qts[1]
  x[x > qts[2]] <- qts[2]
  return(x)
}

#_ -------------------------------------------------------

# fun: plotUsedCellRatio Cycif,CycifStack ----

#' Plot the ratio of used cells
#'
#' This function creates a plot of the ratio of used cells across channels.
#'
#' @export
#'
setGeneric("plotUsedCellRatio",function(x,...) standardGeneric("plotUsedCellRatio"))

#' @rdname plotUsedCellRatio
#' @export
setMethod("plotUsedCellRatio", "Cycif", function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                                                 main="# cells attached on slide",...){
  x <- as.CycifStack(list(x))
  ret <- plotUsedCellRatio(x,cumulative=cumulative,ncycle=ncycle,mar=mar,main=main,...)
  return(invisible(ret))
})


#' @rdname plotUsedCellRatio
#' @export
setMethod("plotUsedCellRatio", "CycifStack", function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                                                      main="# cells attached on slide",smpl.cols,ncol=1,
                                                      leg.cex=0.8,within.rois=NULL,...){
  stopifnot(all(cyApply(x,is,"Cycif")))
  stopifnot(all(unlist(cyApply(x,function(cy)nrow(cy@used_cells))>0)))
  ncycles <- unlist(cyApply(x,nCycles))
  if(missing(ncycle)){
    ncycle <- max(ncycles)
  }

  used.ratio <- data.frame(sapply(names(x),function(n){
    y <- x[[n]]
    rs <- within.rois[[n]]
    # if(!is.null(within.rois)){
    #   rs <- within.rois[[n]]
    # }else{
    #   rs <- FALSE
    # }
    nc.ratio <- statUsedCells(y,within.rois=rs,cumulative=TRUE,ratio=TRUE)
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

#_ -------------------------------------------------------

# fun: plotAvailCellOnSlide Cycif ----
#' @export
setGeneric("plotAvailCellOnSlide",function(x,...) standardGeneric("plotAvailCellOnSlide"))

#' @export
setMethod("plotAvailCellOnSlide", "Cycif",
          function(x,upside.down=TRUE,ncycle,mfrow=c(3,3),mar=c(0,0,4,0),legend=TRUE,main=names(x),cex.title=1,
                   uniq.cols=c(lost="grey90",dropped="blue",available="black",bunched="red"),legend.cex=2,...){
            stopifnot(nrow(x@used_cells)>0)

            u <- x@used_cells # not to be replaced with cumUsedCells
            nchannels <- ncol(u)

            nc.ratio <- round(statUsedCells(x)*100,1)

            xy <- xys(x)
            if(upside.down){
              xy$Y_centroid <- max(xy$Y) - xy$Y
            }

            graphics::par(oma=c(2,2,8,2))
            graphics::par(mfrow=mfrow)
            v <- sapply(seq(nchannels),function(i){
              if(i==1){
                id <- rep(2,nrow(u))
              }else{
                id <- u[,i] + 1
                is.avail <- rowSums(u[,seq(i-1),drop=F]==1)==i-1
                id[!is.avail] <- 0
              }

              rs <- id + 1

              main.txt <- paste0("cycle ",i-1," (",nc.ratio[i],"% available)")
              graphics::par(mar=mar)
              plot(xy,pch='.',col=uniq.cols[rs],main=main.txt,axes=F,asp=1,
                   xlab="",ylab="",...)
              graphics::box()
            })
            if(legend){
              graphics::plot.new()
              graphics::par(mar=c(1,1,1,1))
              # graphics::box()
              legend("topleft",c("available","dropped","bunched","lost previously"),
                     fill=uniq.cols[c("available","dropped","bunched","lost")],
                     cex=legend.cex)
            }
            if(!missing(main)){
              graphics::mtext(main, side=3, line=2, cex=cex.title, outer=TRUE)
            }
          })

#_ -------------------------------------------------------

# fun: vlnPlot CycifStack ----
#' Violin plots to show protein expressions
#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @export
setMethod("vlnPlot", "CycifStack",
          function(x,strat.by=c("celltype","smpl"), ab="PDL1", ctype.full=FALSE,
                   # type=c("raw","log_normalized","logTh_normalized"),
                   strict=FALSE,
                   within.rois = TRUE,uniq.cts,uniq.smpls,is.used){
            if(1){ # missing(type)){
              type <- "logTh_normalized"
            }

            if(missing(ab)){
              stop("ab should be always specified")
            }

            ## gates
            ab_thres <- unlist(x@cell_type@gates[ab,])
            # if(type=="raw"){
            #   ths <- expm1(ab_thres)
            # }else if(type=="log_normalized"){
            #   ths <- ab_thres
            # }else if(type=="logTh_normalized"){
            #   ths <- rep(0.5,length(ab_thres))
            # }
            ths <- rep(0.5,length(ab_thres))

            ## cell types
            cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=within.rois)
            is.strict <- unlist(cyApply(x,function(cy)cy@cell_type@is_strict))
            # table(cts) ## all zeros!!!

            if(strict){
              used <- within_rois(x) & is.strict
            }else{
              used <- within_rois(x)
            }
            if(!missing(is.used)){
              used <- used & is.used
            }

            df <- exprs(x,type=type) %>%
              tibble::rownames_to_column(var="smpl") %>%
              dplyr::mutate(smpl = sub("\\.[0-9]+$","",smpl)) %>%
              dplyr::mutate(celltype=factor(cts)) %>%
              dplyr::filter(used) %>%
              dplyr::left_join(pData(x) %>% dplyr::rename(smpl=id),by="smpl") %>%
              dplyr::mutate(smpl = factor(smpl))

            if(strat.by=="smpl"){
              if(missing(uniq.cts) || length(uniq.cts) != 1){
                stop("one celltype or 'all' should be set in 'uniq.cts'")
              }else if(uniq.cts == "all"){
                ttl <- paste0(ab,", all cells")
                uniq.cts <- levels(df$celltype)
              }else{
                if(length(uniq.cts)>1){
                  uc <- paste0(length(uniq.cts)," samples")
                }else if(length(uniq.cts)>1){
                  uc <- uniq.cts
                }

                ttl <- paste0(ab,", ",uc)
              }
              if(missing(uniq.smpls)){
                uniq.smpls <- levels(df$smpl)
              }
              df1 <- df %>%
                dplyr::filter(smpl %in% uniq.smpls) %>%
                dplyr::filter(celltype %in% uniq.cts) %>%
                dplyr::filter(TimePoint %in% c("BS","BX2","BX3")) %>%
                dplyr::mutate(TimePoint = factor(TimePoint)) %>%
                dplyr::mutate(thres=ab_thres[as.character(smpl)])
              p <- ggplot(df1,aes_string(x="Patient.ID",y=ab,fill="TimePoint")) +
                geom_violin(position = position_dodge(width = 0.9)) +
                ylab(paste0(ab," (",sub("_"," ",type),")")) +
                ggtitle(ttl) +
                geom_point(data=unique(df1 %>% select(Patient.ID,thres,TimePoint)),
                           aes(x=Patient.ID,y=thres,color=as.factor(TimePoint)),
                           position = position_dodge(width = 0.9),
                           shape = 95, size=10, show.legend=F) +
                scale_color_manual(values=rep("black",3),drop=F) +
                theme_bw() +
                theme(legend.position = "right") +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
              print(p)
            }else if(strat.by=="celltype"){ # not so useful?
              # if(missing(uniq.smpls) || length(uniq.smpls)!= 1){
              #   stop("one sample should be set in 'uniq.smpls'")
              # }
              if(missing(uniq.smpls)){
                uniq.smpls <- levels(df$smpl)
              }
              if(missing(uniq.cts)){
                uniq.cts <- levels(df$celltype)
                uniq.cts <- c(uniq.cts[uniq.cts != "Immune_other"],"Immune_other")
                uniq.cts <- c(uniq.cts[uniq.cts != "unknown"],"unknown")
              }
              df1 <- df %>%
                dplyr::mutate(celltype=factor(celltype,levels=uniq.cts)) %>%
                dplyr::filter(smpl %in% uniq.smpls) %>%
                dplyr::filter(!is.na(celltype))

              if(length(uniq.smpls)==1){
                us <- uniq.smpls
              }else if(length(uniq.smpls)>1){
                us <- paste0(length(uniq.smpls), " samples")
              }
              ttl <- paste0(ab,", ",us)
              p <- ggplot(df1,aes_string("celltype",ab)) +
                geom_violin(aes(fill = celltype)) +
                ggtitle(ttl) +
                geom_hline(yintercept = 0.5,col="black") + #ab_thres[[uniq.smpls]],col="black") +
                theme_bw() +
                theme(legend.position = "none") +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
              print(p)
            }
          })

#_ -------------------------------------------------------

# fun: hist_1d numeric ----
#'@importFrom RColorBrewer brewer.pal
#'@importFrom stats loess predict
#'@export
hist_1d <- function(x,n=1000,ths,mar=c(3,4,4,2)+.1,brks1,ttl1){
  omar <- par()$mar
  par(mar=mar)

  n.ab <- trim_fun(x,trim_th=1e-2)
  min.i <- which.min(abs(brks1-min(n.ab)))
  max.i <- which.min(abs(brks1-max(n.ab)))
  # stop(list(min.i,max.i))
  uniq.cols <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(max.i-min.i+1)
  # stop(length(uniq.cols))
  cols <- c(rep(uniq.cols[1],min.i-1),uniq.cols,rep(rev(uniq.cols)[1],n-max.i+1))
  a <- hist(x,breaks=brks1,main=ttl1,freq=FALSE,xlab="",col=cols,border=NA)

  ## smoothening the trail of histogram
  loessMod <- loess(a$density[seq(n)] ~ brks1[seq(n)], span=0.02)
  smoothened <- predict(loessMod)
  lines(smoothened, x=brks1[seq(n)], col=1,lwd=2)

  # cat("Showing current dna_thres\n")
  abline(v=ths[1],col=4,lty=2,lwd=2)
  abline(v=ths[2],col=2,lty=2,lwd=2)
  par(mar=omar)
  invisible(list(a=a,smoothened=smoothened))
}

# fun: h3/heatmap3 matrix, data.frame ----
#' @export
h3 <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
                # distfun = function(x) as.dist(1 - cor(t(x), use = "pa")),
                distfun = dist, distfunC, distfunR, balanceColor = F, ColSideLabs, RowSideLabs,
                showColDendro = T, showRowDendro = T, col = colorRampPalette(c("navy","white", "firebrick3"))(1024),
                legendfun, method = "complete", ColAxisColors = 0, RowAxisColors = 0, hclustfun = hclust,
                reorderfun = function(d, w) reorder(d, w), add.expr, symm = FALSE,
                revC = identical(Colv, "Rowv"), scale = c("row", "column","none"),
                na.rm = TRUE, ColSideFun, ColSideAnn, ColSideWidth = 0.4,
                ColSideCut, colorCell, highlightCell, file = "heatmap3.pdf",
                topN = NA, filterFun = sd, returnDistMatrix = FALSE, margins = c(5,5),
                ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nrow(x)),
                cexCol = 0.2 + 1/log10(ncol(x)), lasRow = 2, lasCol = 2,
                labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, na.color = "grey",
                keep.dendro = FALSE, verbose = getOption("verbose"),
                useRaster = if (ncol(x) * nrow(x) >= 50000) TRUE else FALSE, ...)
{
  hcc <- NULL
  hcr <- NULL
  if (!all(is.na(topN))) {
    temp <- apply(x, 1, filterFun)
    pdf(file)
    for (n in topN) {
      xSub <- x[rev(order(temp))[1:n], , drop = F]
      if (!missing(RowSideColors)) {
        RowSideColorsBak <- RowSideColors
        RowSideColors <- RowSideColors[rev(order(temp))[1:n],
                                       , drop = F]
      }
      result[[paste0(n)]] <- heatmap3(xSub, Rowv = Rowv,
                                      Colv = Colv, distfun = distfun, balanceColor = balanceColor,
                                      ColSideLabs = ColSideLabs, RowSideLabs = RowSideLabs,
                                      showColDendro = showColDendro, showRowDendro = showRowDendro,
                                      col = col, legendfun = legendfun, method = "complete",
                                      ColAxisColors = 0, RowAxisColors = 0, hclustfun = hclust,
                                      reorderfun = reorderfun, add.expr = add.expr,
                                      symm = symm, revC = revC, scale = scale, na.rm = na.rm,
                                      ColSideFun = ColSideFun, ColSideAnn = ColSideAnn,
                                      ColSideWidth = ColSideWidth, ColSideCut = ColSideCut,
                                      margins = margins, ColSideColors = ColSideColors,
                                      RowSideColors = RowSideColors, cexRow = cexRow,
                                      cexCol = cexCol, labRow = labRow, labCol = labCol,
                                      main = paste0("top ", n), xlab = xlab, ylab = ylab,
                                      keep.dendro = keep.dendro, verbose = verbose,
                                      ...)
      if (!missing(RowSideColors)) {
        RowSideColors <- RowSideColorsBak
      }
    }
    temp <- dev.off()
    cat(paste0("The heatmaps were generated at ", file, "\n"))
    return(invisible(result))
  }
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  ## a ---

  if (!missing(ColSideColors)) {
    if (is.vector(ColSideColors)) {
      ColSideColors <- cbind(ColSideColors)
    }
  }
  if (!missing(RowSideColors)) {
    if (is.vector(RowSideColors)) {
      RowSideColors <- cbind(RowSideColors)
    }
  }
  if (missing(distfunC)) {
    distfunC <- distfun
  }
  if (missing(distfunR)) {
    distfunR <- distfun
  }
  distMatrixC = NULL
  distMatrixR = NULL
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  nr <- di[1L]
  nc <- di[2L]
  if (nr <= 1 || nc <= 1)
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2L)
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (!doRdend && identical(Colv, "Rowv"))
    doCdend <- FALSE
  if (is.null(Rowv))
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv))
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram"))
      ddr <- Rowv
    else {
      distMatrixR = distfunR(x)
      hcr <- hclustfun(distMatrixR, method = method)
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv)
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr)))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram"))
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      distMatrixC = distfunC(if (symm)
        x
        else t(x))
      hcc <- hclustfun(distMatrixC, method = method)
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv)
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc)))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc
  x <- x[rowInd, colInd]
  labRow <- if (is.null(labRow))
    if (is.null(rownames(x)))
      (1L:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol))
    if (is.null(colnames(x)))
      (1L:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(1, 4)
  lhei <- c(1 + if (!is.null(main)) 0.2 else 0, 4)
  if (!missing(ColSideFun)) {
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], ColSideWidth, lhei[2L])
  }
  else if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) & nrow(ColSideColors) !=
        nc)
      stop("'ColSideColors' must be a character vector or matrix of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2 * round(ncol(ColSideColors)/2 +
                                      0.1), lhei[2L])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || nrow(RowSideColors) !=
        nr)
      stop("'RowSideColors' must be a character vector or matrix of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2 * round(ncol(RowSideColors)/2 +
                                      0.1), lwid[2L])
  }
  lmat <- lmat + 1
  lmat[is.na(lmat)] <- 0
  lmat[1, 1] <- 1
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  if (balanceColor) {
    if (abs(max(x, na.rm = T)) >= abs(min(x, na.rm = T))) {
      cut.off <- round(quantile(1:length(col), probs = 1 -
                                  (abs(max(x, na.rm = T)) + abs(min(x, na.rm = T)))/(2 *
                                                                                       abs(max(x, na.rm = T)))))
      col <- col[cut.off:length(col)]
    }
    else {
      cut.off <- round(quantile(1:length(col), probs = (abs(max(x,
                                                                na.rm = T)) + abs(min(x, na.rm = T)))/(2 * abs(min(x,
                                                                                                                   na.rm = T)))))
      col <- col[1:cut.off]
    }
  }
  graphics::layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  # layout.show(l)
  # stop()
  if (!missing(legendfun)) {
    par(mar = c(0, 0, 0, 0))
    par(xpd = NA)
    legendfun()
  }
  else {
    par(mar = c(5, 1, 1, 0))
    dummy.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                   length = length(col))
    dummy.z <- matrix(dummy.x, ncol = 1)
    image(x = dummy.x, y = 1, z = dummy.z, yaxt = "n", col = col,
          cex.axis = cexCol, xlab = "")
  }
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1L], 0, 0, 0.5))
    if (revC) {
      rsc = RowSideColors[rev(rowInd), , drop = F]
    }
    else {
      rsc = RowSideColors[rowInd, , drop = F]
    }
    rsc.colors = matrix()
    rsc.names = names(table(rsc))
    rsc.i = 1
    for (rsc.name in rsc.names) {
      rsc.colors[rsc.i] = rsc.name
      rsc[rsc == rsc.name] = rsc.i
      rsc.i = rsc.i + 1
    }
    rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
    image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
    if (missing(RowSideLabs)) {
      if (ncol(RowSideColors) == 1 & colnames(RowSideColors)[1] ==
          "") {
        RowSideLabs <- ""
      }
      else {
        RowSideLabs <- colnames(RowSideColors)
      }
    }
    if (dim(rsc)[2] == 1) {
      axis(1, 0, RowSideLabs, las = 2, tick = FALSE)
    }
    else {
      axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), RowSideLabs,
           las = 2, tick = FALSE)
    }
  }
  if (!missing(ColSideCut)) {
    ColSideCutResult <- cut(ddc, ColSideCut)$lower
    ColSideCutResultSubIndList <- list()
    for (i in 1:length(ColSideCutResult)) {
      ColSideCutResultSubInd <- order.dendrogram(ColSideCutResult[[i]])
      ColSideCutResultSubIndList[[i]] <- ColSideCutResultSubInd
    }
    cutTable <- NULL
    if (verbose) {
      cat(paste0("The samples could be cut into ", length(ColSideCutResult),
                 " parts with height ", ColSideCut))
      cat("\n")
      if (!missing(ColSideAnn)) {
        for (i in 1:ncol(ColSideAnn)) {
          if (is.factor(ColSideAnn[, i])) {
            cutTable[[i]] <- sapply(ColSideCutResultSubIndList,
                                    function(x) table(ColSideAnn[x, i]))
            colnames(cutTable[[i]]) <- paste0("Cluster ",
                                              1:length(ColSideCutResult))
            names(cutTable)[i] <- colnames(ColSideAnn)[i]
            pvalue <- chisq.test(cutTable[[i]])$p.value
            cat(paste0("Differential distribution for ",
                       colnames(ColSideAnn)[i], ", p value by chi-squared test: ",
                       round(pvalue, 3), "\n"))
            cutTable[[i]] <- rbind(cutTable[[i]], round(cutTable[[i]][1,
            ]/colSums(cutTable[[i]]), 2))
            row.names(cutTable[[i]])[nrow(cutTable[[i]])] <- paste0(row.names(cutTable[[i]])[1],
                                                                    "_Percent")
            cutTable[[i]] <- cbind(cutTable[[i]], pValue = c(pvalue,
                                                             rep(NA, nrow(cutTable[[i]]) - 1)))
          }
          else {
            cutTable[[i]] <- sapply(split(ColSideAnn[unlist(ColSideCutResultSubIndList),
                                                     i], rep(1:length(ColSideCutResultSubIndList),
                                                             sapply(ColSideCutResultSubIndList, length))),
                                    function(x) summary(na.omit(x)))
            colnames(cutTable[[i]]) <- paste0("Cluster ",
                                              1:length(ColSideCutResult))
            names(cutTable)[i] <- colnames(ColSideAnn)[i]
            temp <- aov(ColSideAnn[unlist(ColSideCutResultSubIndList),
                                   i] ~ as.factor(rep(1:length(ColSideCutResultSubIndList),
                                                      sapply(ColSideCutResultSubIndList, length))))
            pvalue <- summary(temp)[[1]]$"Pr(>F)"[1]
            cat(paste0("Differential distribution for ",
                       colnames(ColSideAnn)[i], ", p value by ANOVA: ",
                       round(pvalue, 3), "\n"))
            cutTable[[i]] <- cbind(cutTable[[i]], pValue = c(pvalue,
                                                             rep(NA, 5)))
          }
        }
      }
    }
    ColSideCutResultCol <- rainbow(length(ColSideCutResult),
                                   alpha = 0.2)
    ColNumber <- (ncol(x) - 1)
  }
  if (!missing(ColSideFun)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    ColSideAnn <- ColSideAnn[colInd, , drop = F]
    ColAnnHeight <- ColSideFun(ColSideAnn)
    if (!exists("ColAnnHeight")) {
      ColAnnHeight <- par("usr")[3:4]
    }
    if (!missing(ColSideCut)) {
      rect(c(0 - 1/ColNumber/2, (0 - 1/ColNumber/2) + 1/ColNumber *
               cumsum(sapply(ColSideCutResult, function(x) length(unlist(x))))[-length(ColSideCutResult)]),
           ColAnnHeight[1], c((0 - 1/ColNumber/2) + 1/ColNumber *
                                cumsum(sapply(ColSideCutResult, function(x) length(unlist(x))))),
           ColAnnHeight[2], col = ColSideCutResultCol)
    }
  }
  else if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    csc = ColSideColors[colInd, , drop = F]
    csc.colors = matrix()
    csc.names = names(table(csc))
    csc.i = 1
    for (csc.name in csc.names) {
      csc.colors[csc.i] = csc.name
      csc[csc == csc.name] = csc.i
      csc.i = csc.i + 1
    }
    csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
    image(csc, col = as.vector(csc.colors), axes = FALSE)
    if (missing(ColSideLabs)) {
      if (ncol(ColSideColors) == 1 & colnames(ColSideColors)[1] ==
          "") {
        ColSideLabs <- ""
      }
      else {
        ColSideLabs <- colnames(ColSideColors)
      }
    }
    if (dim(csc)[2] == 1) {
      axis(4, 0, ColSideLabs, las = 2, tick = FALSE)
    }
    else {
      axis(4, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), ColSideLabs,
           las = 2, tick = FALSE)
    }
  }
  par(mar = c(margins[1L], 0, 0, margins[2L]))
  if (!symm || scale != "none")
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend)
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        useRaster = useRaster, ...)
  if (any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1L:nc, 1L:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, useRaster=useRaster, add = TRUE)
  }
  if (!missing(colorCell)) {
    colorCell[, 1] <- match(colorCell[, 1], rowInd)
    colorCell[, 2] <- match(colorCell[, 2], colInd)
    rect(colorCell[, 2] - 0.5, colorCell[, 1] - 0.5, colorCell[,
                                                               2] + 0.5, colorCell[, 1] + 0.5, col = as.character(colorCell[,
                                                                                                                            3]), border = NA)
  }
  if (!missing(highlightCell)) {
    if (ncol(highlightCell) == 3) {
      highlightCell$lwd <- 1
    }
    highlightCell[, 1] <- match(highlightCell[, 1], rowInd)
    highlightCell[, 2] <- match(highlightCell[, 2], colInd)
    rect(highlightCell[, 2] - 0.5, highlightCell[, 1] - 0.5,
         highlightCell[, 2] + 0.5, highlightCell[, 1] + 0.5,
         border = as.character(highlightCell[, 3]), lwd = as.integer(highlightCell[,
                                                                                   4]))
  }
  if (!missing(ColSideColors) & ColAxisColors != 0) {
    mtext(1, at = 1L:nc, text = labCol, las = lasCol, line = 0.5,
          cex = cexCol, col = ColSideColors[colInd, ColAxisColors])
  }
  else {
    axis(1, 1L:nc, labels = labCol, las = lasCol, line = -0.5,
         tick = 0, cex.axis = cexCol)
  }
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1L] - 1.25)
  if (!missing(RowSideColors) & RowAxisColors != 0) {
    mtext(4, at = iy, text = labRow, las = lasRow, line = 0.5,
          cex = cexRow, col = RowSideColors[rowInd, RowAxisColors])
  }
  else {
    axis(4, iy, labels = labRow, las = lasRow, line = -0.5,
         tick = 0, cex.axis = cexRow)
  }
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2L] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend & showRowDendro)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend & showColDendro) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    if (!missing(ColSideCut)) {
      rect(c(0.5, 0.5 + cumsum(sapply(ColSideCutResult,
                                      function(x) length(unlist(x))))[-length(ColSideCutResult)]),
           0, cumsum(sapply(ColSideCutResult, function(x) length(unlist(x)))) +
             0.5, ColSideCut, col = ColSideCutResultCol)
    }
  }
  else if (!is.null(main))
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
                                                              doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc,
                 cutTable = if (!missing(ColSideAnn) && !missing(ColSideCut)) cutTable,
                 cutColoumIndList = if (!missing(ColSideCut)) ColSideCutResultSubIndList,
                 DistMatrixC = if (returnDistMatrix) distMatrixC, DistMatrixR = if (returnDistMatrix) distMatrixR,
                 hcr = hcr, hcc = hcc))
}

# fun: heatmapSummAb (Cycif),CycifStack ----
#' Violin plots to show protein expressions
#' @export
setGeneric("heatmapSummAb", function(x,...) standardGeneric("heatmapSummAb"))

#' @export
setMethod("heatmapSummAb", "CycifStack",
          function(x,norm_type=c("log_normalized","logTh_normalized"), ab,
                   sum_type=c("freq","mean","median","x percentile"),
                   ctype.full=FALSE,rois = TRUE,
                   uniq.cts,uniq.smpls,strict=FALSE,scale="row",...){
            options(dplyr.summarise.inform = FALSE)

            if(missing(sum_type)){
              stop("'sum_type' should be specified")
            }
            if(missing(norm_type)){
              norm_type <- "log_normalized"
            }
            if(missing(ctype.full)){
              ctype.full <- FALSE
            }
            if(missing(ab)){
              stop("ab should be always specified")
            }

            if(grepl("percentile$",sum_type)){
              pct <- strsplit(sum_type," ")[[1]]
              p.ile <- as.numeric(pct[[1]])
              sum_fun=function(y,...)quantile(y,p.ile/100,na.rm=T)
            }else if(sum_type=="median"){
              sum_fun=function(y,...)quantile(y,.5,na.rm=T)
            }else if(sum_type=="mean"){
              sum_fun=function(y,...)mean(y,na.rm=T)
            }else if(sum_type=="freq"){
              thres <- unlist(x@cell_type@gates[ab,]) # identical whichever ctype.full
              norm_type="logTh_normalized"
              sum_fun=function(y,th){
                if(length(y)==0){
                  return(NA)
                }else{
                  return(sum(y > 0.5,na.rm=T)/sum(!is.na(y)))
                }
              }
            }

            if(sum_type=="freq" || norm_type=="log_normalized"){
              balanceColor <- FALSE
            }else if(norm_type=="logTh_normalized"){
              balanceColor <- TRUE
              scale <- "none"
            }

            ## cell types
            cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=rois,strict=strict)

            df <- exprs(x,type=norm_type) %>%
              tibble::rownames_to_column(var="smpl") %>%
              mutate(smpl = sub("\\.[0-9]+$","",smpl)) %>%
              mutate(celltype=factor(cts))

            if(missing(uniq.smpls)){
              uniq.smpls <- unique(df$smpl)
            }
            if(missing(uniq.cts)){
              uniq.cts <- levels(df$celltype)
              uniq.cts <- uniq.cts[!uniq.cts %in% c("Immune_other","unknown")]
            }

            pd <- pData(x) %>% rename(smpl=id)
            df <- df %>%
              filter(within_rois(x)) %>%
              filter(!is.na(celltype)) %>%
              filter(smpl %in% uniq.smpls) %>%
              filter(celltype %in% uniq.cts) %>%
              mutate(celltype=factor(celltype,levels=uniq.cts)) %>%
              rename(this_ab=as.symbol(ab)) %>%
              mutate(thres1 = thres[smpl]) %>%
              group_by(smpl,celltype) %>%
              summarise(sum_ab=sum_fun(this_ab,thres1)) %>%
              left_join(pd,by="smpl") %>%
              mutate(smpl = factor(smpl))

            ttl <- paste0(ab,",",sum_type,",",sub("_.+","",norm_type),",(scale:",scale,")")
            m1 <- as.matrix(with(df, sparseMatrix(as.integer(smpl), as.integer(celltype), x=sum_ab)))
            m1[apply(m1,1,function(x)any(is.nan(x))),] <- NA
            if(sum_type!="freq"){
              m1[m1==0] <- NA
            }
            rownames(m1) <- levels(df$smpl)
            colnames(m1) <- levels(df$celltype)
            m1 <- m1[rev(seq(nrow(m1))),]

            if(sum_type!="freq" && norm_type=="logTh_normalized"){
              m1 <- m1 - 0.5 ## hardcode!!!!
            }

            h3(m1,
               margins=c(10,5),
               main=ttl,
               cexRow=.7,
               na.rm = T,
               balanceColor=balanceColor,
               scale=scale,...)
          })

# fun: heatmapSCExp CycifStack ----
#' Violin plots to show protein expressions
#' @export
setGeneric("heatmapSCExp", function(x,...) standardGeneric("heatmapSCExp"))

#' @export
setMethod("heatmapSCExp", "CycifStack",
          function(x,ld_name,norm_type=c("logTh_normalized","log_normalized"),
                   summarize=FALSE,summarize.by=c("BOR","TimePoint"),
                   order="cBTgP",margins=c(7,5),used.abs,ttl,...){
            options(dplyr.summarise.inform = FALSE)
            if(missing(norm_type)){
              norm_type <- "logTh_normalized"
            }

            ##
            if(missing(ld_name)){
              stop("'ld_name' should be specified (it's used to retrieve the data later)")
            }else if(!ld_name %in% ld_names(x)){
              stop("Specified 'ld_name' does not exist.")
            }else{
              ld <- ld_coords(x,ld_name=ld_name)
            }

            if(summarize){
              meta <- createRSC(x,with.clusters=T,ld_name=ld_name)
            }else{
              meta <- createRSC(x,with.clusters=T,ld_name=ld_name,order=order)
            }

            pd <- meta$pd
            rsc <- meta$rsc

            xys <- ld@ld_coords
            ld.type <- ld@ld_type
            norm_type <- ld@norm_type

            is.used <- ld@is_used
            if(missing(used.abs)){
              used.abs <- ld@used.abs
            }
            used.cts <- ld@used.cts
            cts.params <- ld@cts_params
            cls <- ld@clusters

            n <- exprs(x,type=norm_type)[is.used,used.abs]

            if(summarize){
              if(!all(summarize.by %in% colnames(rsc))){
                na.cols <-summarize.by[!summarize.by %in% colnames(rsc)]
                stop("'columns in 'summarize.by' (",paste(na.cols,collapse=","),") are not specified in 'order' (",paste(colnames(rsc),collapse=","),")")
              }
              # n1 <- n %>%
              #   tibble::rownames_to_column("smpl") %>%
              #   mutate(smpl=sub("\\..+","",smpl)) %>%
              #   left_join(pd[c("smpl",summarize.by)],by="smpl") %>%
              #   group_by(!!!syms(summarize.by)) %>%
              #   summarize_at(used.abs,mean,na.rm=T)
              n1 <- cbind(n,pd) %>%
                group_by(!!!syms(summarize.by)) %>%
                summarize_at(used.abs,mean,na.rm=T)
              rn <- apply(as.matrix(n1[summarize.by]),1,paste,collapse=",")

              n1 <- as.matrix(n1[used.abs])
              rownames(n1) <- rn

              rsc1 <- unique(cbind(pd[summarize.by],rsc[,summarize.by]))
              names(rsc1) <- paste0(rep(c("pd",""),each=length(summarize.by)),rep(summarize.by,times=2))
              rsc1 <- rsc1 %>% arrange(!!!syms(paste0("pd",summarize.by))) %>% select(!!! syms(summarize.by))
              rsc <- as.matrix(rsc1)
              rsl <- colnames(rsc)

              ## reverse row-order
              row.o <- rev(seq(nrow(n1)))
              n1 <- n1[row.o,]
              rsc <- rsc1[row.o,]
            }else{
              row.o <- meta$row.order
              n1 <- as.matrix(n[row.o,])
              rownames(n1) <- rep("",nrow(n1))
              rsl <- colnames(rsc)
            }

            if(missing(ttl)){
              ttl <- paste0(ld_name,", ",sub("_.+","",norm_type)," (no scaling)")
            }
            h3(n1,Colv=NA,Rowv=NA,scale="none",RowSideColors=rsc,useRaster=F,main=ttl,margins=margins,RowSideLabs=rsl,...)
          })

#_ -------------------------------------------------------

# fun: CellTypeGraph ctype ----
#'@export
CellTypeGraph <- function(ctype,plot=F,transpose=T,...){
  uniq.cts <- c("all",ctype$Child)
  ctgraph <- ctype[c("Parent","Child")]
  ctgraph$Parent <- factor(ctgraph$Parent,levels=uniq.cts)
  ctgraph$Child <- factor(ctgraph$Child,levels=uniq.cts)
  g <- igraph::graph_from_data_frame(ctgraph)
  l <- igraph::layout_as_tree(g)
  levs <- l[,2]
  names(levs) <- names(igraph::V(g))
  ctlevs <- rev(tapply(names(levs),levs,identity))
  if(plot){
    if(transpose){
      l[,2] <- max(l[,2]) - l[,2]
      l <- l[,2:1]
    }
    igraph::V(g)$shape <- "rectangle"
    igraph::V(g)$label.family <- "Helvetica"

    plot(g, layout=l,
         edge.arrow.size=.5, vertex.color="gold", vertex.size=40,
         vertex.frame.color=NA, vertex.label.color="black",
         vertex.label.cex=0.8, vertex.label.dist=0, edge.curved=0
         ,...)
  }
  return(ctlevs)
}

#_ -------------------------------------------------------

# fun: lowDimPlot Cycif, CycifStack ----

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
          function(x,ld_name,plot_type=c("celltype","cluster","exp"),ab,
                   na.col = "grey80",uniq.cols,with.labels = TRUE,leg=TRUE,
                   pch=".",main,p_thres=0.5,mar,cex.labs=1,cex.leg=1,cex.main=2,...){
            if(missing(plot_type)){
              stop("Which plot type? (plot_type = celltype, cluster, exp)")
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
            cls <- ld@clusters

            if(plot_type=="exp"){
              if(missing(ab)){
                stop("'ab' should be specified when plot_type == 'exp'")
              }
              leg <- with.labels <- FALSE

              n <- exprs(x,type=norm.type)[is.used,ab]

              ## fix color range
              nc <- 100
              uniq.cols <- rev(RColorBrewer::brewer.pal(11,"RdBu"))

              idx <-round(transform(seq(nc+1),method="Th",th=round((nc+1)/2),p_thres=p_thres,trim=0)*nc)+1
              adj.uniq.cols <- colorRampPalette(uniq.cols)(nc+1)[idx]
              cols <- adj.uniq.cols[round(n*nc)+1]

              ## rename plot_type
              plot_type <- paste0(ab,"(exp)")

            }else if(plot_type=="celltype"){
              facs <- cell_types(x,
                                 ctype.full = cts.params$ctype.full,
                                 leaves.only = cts.params$leaves.only,
                                 strict = cts.params$strict)[is.used]
              facs <- factor(facs,levels = used.cts)
              plot_type <- "Cell types"
            }else if(plot_type=="cluster"){
              facs <- ld@clusters
              plot_type <- "Clusters"
            }

            ## main
            if(missing(main)){
              main <- paste(ld.type,ld_name,plot_type,sep=", ")
            }

            ## colors - when plot_type
            if(!grepl("exp",plot_type)){
              if(missing(uniq.cols)){
                ## levels
                levs <- levels(facs)
                nlev <- nlevels(facs)

                if(nlev < 9){
                  uniq.cols <- RColorBrewer::brewer.pal(8,"Dark2")[seq(nlev)]
                }else{
                  uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nlev)
                }
                names(uniq.cols) <- levs
              }

              if(with.labels){
                uniq.cols <- sapply(uniq.cols,function(uc){
                  colorRampPalette(c(uc,"white"))(3)[2]
                })
              }

              ## dealing NA
              if(any(is.na(facs))){
                facs <- as.numeric(facs)
                facs[is.na(facs)] <- max(facs,na.rm=T)+1
                facs <- factor(facs,labels=c(levs,"NA"))
              }

              cols <- uniq.cols[as.numeric(facs)]
            }

            if(missing(mar)){
              if(leg){
                mar <- c(4,4,4,10)
              }else{
                mar <- c(4,4,4,4)
              }
            }

            ## plot
            def.par <- par(no.readonly = TRUE)
            par(mar=mar)
            plot(xys,col=cols,pch=pch,main=main,xlab="",ylab="",axes=F,cex=cex,
                 cex.main=cex.main,...)
            box()


            if(leg){
              par(xpd=T)
              legend(par()$usr[2],par()$usr[4],levs,pch=pch,col=uniq.cols,cex=cex.leg)
              par(xpd=F)
            }
            if(with.labels){
              mids <- sapply(xys,function(x){
                tapply(x,facs,mean)
              })
              text(mids,rownames(mids),cex=cex.labs)
            }
          }
)

#' @rdname lowDimPlot
#' @export
setMethod("lowDimPlot", "CycifStack",
          function(x,ld_name,plot_type=c("celltype","cluster","exp"),ab,
                   na.col = "grey80",uniq.cols,with.labels = TRUE,leg=TRUE,
                   cts.names,pal.name,cex.main=2,
                   pch=".",main,p_thres=0.5,mar,cex.labs=1,cex.leg=.5,cex=cex,...){
            if(missing(plot_type)){
              stop("Which plot type? (plot_type = celltype, cluster, exp)")
            }

            if(missing(ld_name)){
              stop("'ld_name' should be specified (it's used to retrieve the specific UMAP/Clustgering)")
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
            cls <- ld@clusters

            if(plot_type=="exp"){
              if(missing(ab)){
                stop("'ab' should be specified when plot_type == 'exp'")
              }
              leg <- with.labels <- FALSE

              n <- exprs(x,type=norm.type)[is.used,ab]

              ## fix color range
              nc <- 100
              uniq.cols <- rev(RColorBrewer::brewer.pal(11,"RdBu"))

              idx <-round(transform(seq(nc+1),method="Th",th=round((nc+1)/2),p_thres=p_thres,trim=0)*nc)+1
              adj.uniq.cols <- colorRampPalette(uniq.cols)(nc+1)[idx]
              cols <- adj.uniq.cols[round(n*nc)+1]

              ## rename plot_type
              plot_type <- paste0(ab,"(exp)")

            }else if(plot_type=="celltype"){
              facs <- cell_types(x,
                                 ctype.full = cts.params$ctype.full,
                                 leaves.only=cts.params$leaves.only)[is.used]
              levels(facs) <- sub(",.+","",levels(facs))
              if(!missing(cts.names)){
                if(!all(names(cts.names) %in%  used.cts)){
                  stop("all names(cts.names) should be original cell type names")
                }
                facs <- cts.names[as.character(facs)]
                facs <- factor(facs,levels=unique(cts.names))
              }else{
                facs <- factor(facs,levels=used.cts)
              }
              plot_type <- "Cell types"
            }else if(plot_type=="cluster"){
              facs <- ld@clusters
              plot_type <- "Clusters"
            }else if(plot_type=="smpl"){
              n.cells.smpls <- ld@n_cells_total
              uniq.smpls <- names(n.cells.smpls)

              facs <- factor(rep(uniq.smpls,n.cells.smpls),levels=uniq.smpls)
              plot_type <- "Samples"
            }else{
              pd <- pData(x)
              if(!plot_type %in% names(pd)){
                stop("'plot_type' should be 'exp','celltype','smpl','cluster' or one of colnames(pd)" )
              }
              n.cells.smpls <- ld@n_cells_total
              uniq.smpls <- names(n.cells.smpls)

              facs <- factor(rep(uniq.smpls,n.cells.smpls),levels=uniq.smpls)
              if(plot_type %in% c("Cohort","gBRCA.status","Subtype","Patient.ID")){
                pd[[plot_type]] <- factor(pd[[plot_type]])
              }

              this.levs <- (pd %>%
                              filter(id %in% uniq.smpls) %>%
                              arrange(match(uniq.smpls,id)))[[plot_type]]

              levels(facs) <- this.levs
              facs <- factor(facs,levels=levels(this.levs))
            }

            ## main
            if(missing(main)){
              main <- paste(ld.type,ld_name,plot_type,sep=", ")
            }



            ## colors - when plot_type
            if(!grepl("exp",plot_type)){
              levs <- levels(facs)
              if(missing(uniq.cols)){
                ## levels
                nlev <- nlevels(facs)

                if(nlev < 9){
                  if(missing(pal.name)){
                    pal.name <- "Dark2"
                  }
                  uniq.cols <- RColorBrewer::brewer.pal(8,pal.name)[seq(nlev)]
                }else{
                  if(missing(pal.name)){
                    pal.name <- "Spectral"
                  }
                  uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(11,pal.name))(nlev)
                }
                names(uniq.cols) <- levs
              }

              if(with.labels){
                uniq.cols <- sapply(uniq.cols,function(uc){
                  colorRampPalette(c(uc,"white"))(3)[2]
                })
              }

              ## dealing NA
              if(any(is.na(facs))){
                facs <- as.numeric(facs)
                facs[is.na(facs)] <- max(facs,na.rm=T)+1
                facs <- factor(facs,labels=c(levs,"NA"))
              }

              cols <- uniq.cols[as.numeric(facs)]
            }

            if(missing(mar)){
              if(leg){
                mar <- c(4,4,4,6)
              }else{
                mar <- c(4,4,4,6)
              }
            }

            ## plot
            def.par <- par(no.readonly = TRUE)
            par(mar=mar)
            plot(xys,col=cols,pch=pch,main=main,xlab="",ylab="",cex=cex,axes=F,
                 cex.main=cex.main,...)
            box()

            if(leg){
              par(xpd=T)
              legend(par()$usr[2],par()$usr[4],levs,pch=pch,col=uniq.cols,cex=cex.leg)
              par(xpd=F)
            }
            if(with.labels){
              mids <- sapply(xys,function(x){
                tapply(x,facs,mean)
              })
              text(mids,rownames(mids),cex=cex.labs)
            }
          }
)

#_ -------------------------------------------------------
