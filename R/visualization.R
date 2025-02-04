#_ -------------------------------------------------------

# fun: AbsSummary CycifStack ----

#' A graphical summary of cycles and antibodies of a CycifStack class.
#'
#' Generates an absolute summary plot for a CycifStack object.
#'
#' @param x A CycifStack object.
#' @param show.cycles.in.row Logical, indicating whether to display cycle information in row labels.
#' FALSE by default. If TRUE, the output plot shows the number of cycles each sample was stained for.
#' @param ... Additional parameters to be passed to the \code{graphics::image} function.
#'
#' @return
#' The function generates a graphical absolute summary plot and returns \code{NULL}.
#'
#' @details
#' This function generates a plot that summarizes the antibodies used in a CycifStack object.
#' It visualizes the antibodies across different samples and cycles.
#'
#' @seealso
#' \code{\link{abs_list}}
#'
#' @importFrom graphics image box abline axis par
#' @importFrom grDevices grey
#'
#' @export
AbsSummary <- function(x,show.cycles.in.row=FALSE,...){
  if(!is(x,"CycifStack")){
    stop("input should be a CycifStack obj.")
  }
  uniq.abs <- as.character(x@abs_list$ab)
  n1 <- data.frame(do.call(rbind,lapply(x@samples,function(y){
    this.abs <- as.character(abs_list(y)$ab)
    tested <- uniq.abs %in% this.abs
    return(tested)
  })))

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
#
#' @title Create a Slide Plot
#'
#' @description This function generates a slide plot for visualizing different aspects of Cycif data, including DNA intensity,
#' protein expression, cell types, or custom annotations.
#'
#' @param x A Cycif object.
#' @param pch Plotting character for points (default is 20).
#' @param cex Plotting character size (default is 2).
#' @param bg.col Background color for points.
#' @param plot_type Type of slide plot to create: "dna" for DNA intensity, "exp" for protein expression,
#'   "cell_type" for cell types, or "custom" for custom annotations.
#' @param custom_labs (For plot_type="custom") A factor vector of custom annotations to be used as labels.
#' @param ct_name The name of the cell type variable in the data.
#' @param strict (For plot_type="cell_type") If TRUE, cell type labels must exactly match those in the data.
#' @param ttl Title for the slide plot.
#' @param leg.ttl Legend title.
#' @param ab The name of the channel for plotting (e.g., "DNA1", "CD45", etc.).
#' @param uniq.cts (For plot_type="cell_type") A vector of unique cell type labels to be displayed.
#' @param uniq.cols (For plot_type="cell_type") A vector of unique colors corresponding to unique cell types.
#' @param draw.roi Should regions of interest (ROIs) be drawn on the slide plot? (default is FALSE)
#' @param show.na Should cells with missing values be displayed? (default is TRUE)
#' @param na.col Color for cells with missing values.
#' @param use_rois Should ROIs be used in the analysis? (default is TRUE)
#' @param use.thres Should thresholding be applied to the data? (default is TRUE)
#' @param contour Should contour lines be added to the plot? (default is FALSE)
#' @param cont_nlevs The number of contour levels (default is 3).
#' @param ncells Maximum number of cells to display (default is 1e4).
#' @param trim_th Trimming threshold for outlier removal (default is 1e-2).
#' @param legend Should a legend be displayed on the slide plot? (default is FALSE)
#' @param legend.pos Position of the legend on the slide plot (default is automatic).
#' @param mar Margins for the plot (default is c(3, 3, 3, 3)).
#' @param roi.rec A list specifying the region of interest (ROI) rectangle, including 'x' and 'y' ranges.
#' @param ... Additional graphical parameters to customize the plot.
#'
#' @details
#' - The `slidePlot` function generates a slide plot for visualizing Cycif data.
#' - You can choose the type of slide plot using the `plot_type` parameter, which can be "dna" for DNA intensity,
#'   "exp" for protein expression, "cell_type" for cell types, or "custom" for custom annotations.
#' - The function allows you to display ROIs, customize color schemes, and control various plotting parameters.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_polygon guides guide_colorbar guide_legend scale_color_distiller scale_color_manual
#' @importFrom dplyr mutate select filter rename
#'
#' @export
#' @rdname slidePlot
setGeneric("slidePlot", function(x,...) standardGeneric("slidePlot"))

#' @export
#' @rdname slidePlot
setMethod("slidePlot", "Cycif",
          function(x,pch=20,cex=2,bg.col,
                   plot_type=c("dna","exp","cell_type","custom"),
                   custom_labs,ct_name="default",strict=FALSE,ttl,leg.ttl,ab,
                   uniq.cts,uniq.cols,draw.roi=FALSE,show.na=TRUE,
                   na.col="grey80",use_rois=TRUE,use.thres=TRUE,
                   contour=FALSE,cont_nlevs=3,ncells=1e4,
                   trim_th=1e-2,legend=FALSE, legend.pos,mar=c(3,3,3,3),
                   roi.rec,...){
            if(missing(plot_type)){
              stop("need to specify 'plot_type' argument: dna, exp, cell_type, custom")
            }
            smpl <- names(x)

            if(use_rois && length(x@within_rois)==0){
                stop("ROIs not defined yet; 'use_rois' should be FALSE")
            }

            ## coordinates
            xy <- xys(x)
            max.y <- max(xy$Y_centroid)
            xy <- xy %>%
              mutate(Y_centroid = max(Y_centroid) - Y_centroid)

            ## cell types - specified by `uniq.cts`
            if(ct_name %in% names(x@cell_types)){
              cts <- cell_types(x,strict=strict)

              if(missing(uniq.cts)){
                uniq.cts <- levels(cts$cell_types)
                if(any(uniq.cts %in% c("unknown","outOfROI"))){
                  uniq.cts <- uniq.cts[!uniq.cts %in% c("unknown","outOfROI")]
                }
              }

              cts <- cts %>% mutate(cell_types=factor(cell_types,levels=uniq.cts)) # some are NAs
            }

            ## values to highlight color
            ## plot_type == "dna" ----
            if(plot_type=="dna"){
              if(missing(ab) || !ab %in% paste0("DNA",seq(nCycles(x)))){
                stop(paste0("If DNA stain, `ab' should be one of DNA1, ..., DNA", nCycles(x),")"))
              }else if(names(x@dna)[1]!="DNA1"){
                stop("DNA channels should be named as DNA1, DNA2, ... (1-origin)")
              }

              if(missing(ttl)){
                if(ab=="DNA1"){
                  ab1 <- ab
                }else{
                  ab1 <- paste0(ab,"/DNA1")
                }
                ttl <- paste0(smpl,", ",ab1)
              }
              if(missing(leg.ttl)){
                leg.ttl <- ab
              }

              ## Intensity of the channel
              n <- log1p(x@dna)

              mat <- n[1] %>%
                cbind(sapply(n[-1],function(x)x-n[[1]])) %>%
                dplyr::mutate(ab=trim_fun(!!sym(ab),trim_th=trim_th)) %>%
                mutate(is.na = (is.na(!!sym(ab))|!x@within_rois)) %>%
                select(ab, is.na)

              ## combine the mat
              df <- xy %>% cbind(mat)

            }else if(plot_type=="exp"){
              ## plot_type == "exp" -----
              if(missing(ab) || !ab %in% abs_list(x)$ab){
                stop(ab, " is not specified or available in the sample ", names(x))
              }

              if(missing(ttl)){
                ttl <- paste0(smpl,", ",ab," expression")
              }
              if(missing(leg.ttl)){
                leg.ttl <- ab
              }

              ## use ggplot now
              mat <- exprs(x,type="log") %>%
                dplyr::mutate(ab=trim_fun(!!sym(ab),trim_th=trim_th)) %>%
                mutate(is.na = (is.na(!!sym(ab)) | !x@within_rois)) %>%
                select(ab, is.na)

              ## combine the mat
              df <- xy %>% cbind(mat)

            }else if(plot_type=="cell_type"){
              ## plot_type == "cell_type" -----
              if(!missing(ab)){
                stop("ab shouldn't be specified when plot_type='cell_type'")
              }

              if(missing(ttl)){
                ttl <- paste0(smpl,", cell types")
              }
              if(missing(leg.ttl)){
                leg.ttl <- "cell types"
              }

              ## is.na
              mat <- cts %>%
                mutate(is.na = (is.na(cell_types)) | !x@within_rois) %>%
                dplyr::rename(ab="cell_types") %>%
                select(ab,is.na)

              ## combine the mat
              df <- xy %>% cbind(mat)

            }else if(plot_type=="custom"){
              ## plot_type == "custom" -----

              if(missing(custom_labs)){
                stop("when plot_type='custom', the argument 'custom_labs' should be specified.")
              }else if(length(custom_labs)!=nrow(xy)){
                stop("length(custom_labs) should be the same as nrow(xy)")
              }else if(!is(custom_labs,"factor")){
                custom_labs <- factor(custom_labs)
              }

              if(missing(ttl)){
                ttl <- "custom"
              }
              if(missing(leg.ttl)){
                leg.ttl <- "custom"
              }

              df <- xy %>%
                cbind(data.frame(ab=custom_labs,is.na=is.na(custom_labs)))
            }else{
              stop("'plot_type' should be one of dna,exp,cell_type,custom")
            }

            ## guide_legend ----
            if(is(df$ab,"numeric")){
              # stop("numeric")
              my_guides <- function(txt)guides(color=guide_colorbar(title=txt))
              cols <- function(uc)ggplot2::scale_color_distiller(palette="Spectral",direction=-1)
            }else{
              my_guides <- function(txt)guides(color=guide_legend(title=txt,ncol=1,override.aes = list(size = 2)))
              cols <- function(uc)ggplot2::scale_color_manual(labels=names(uc),values=uc)
            }

            ## roi.rec
            if(missing(roi.rec)){
              roi.rec <-
                list(x=range(xy$X_centroid) + c(-0.1,0.1),
                     y=range(xy$Y_centroid) + c(-0.1,0.1))
            }

            df <- df %>%
              mutate(within.roi =
                X_centroid > roi.rec$x[1] &
                X_centroid < roi.rec$x[2] &
                Y_centroid > roi.rec$y[1] &
                Y_centroid < roi.rec$y[2])

            ## choose cells to sample - depends on whether to show NAs
            idx1 <- which(!df$is.na & df$within.roi)
            idx2 <- which(df$is.na & df$within.roi)

            if(!show.na){
              used.ratio <- min(ncells/length(idx1),1)

              if(used.ratio < 1){
                set.seed(123)
                sorted.idx1 <- sort(sample(idx1,ncells))
              }else{
                sorted.idx1 <- idx1
              }
              used1 <- seq(nrow(df)) %in% sorted.idx1

              df <- df %>% mutate(is.used1=used1)
            }else{
              used.ratio <- min(ncells/(length(idx1)+length(idx2)),1)
              ncells.pos <- round(length(idx1) * used.ratio)
              ncells.neg <- round(length(idx2) * used.ratio)

              if(used.ratio < 1){
                set.seed(123)
                sorted.idx1 <- sort(sample(idx1,ncells.pos))
                sorted.idx2 <- sort(sample(idx2,ncells.neg))
              }else{
                sorted.idx1 <- idx1
                sorted.idx2 <- idx2
              }

              used1 <- seq(nrow(df)) %in% sorted.idx1
              used2 <- seq(nrow(df)) %in% sorted.idx2

              df <- df %>%
                mutate(is.used1=used1) %>%
                mutate(is.used2=used2)
            }

            ## core plot function ----
            p <- ggplot(df %>% filter(is.used1) ,aes(x=X_centroid,y=Y_centroid,color=ab))
            if(show.na){
              p <- p +
                geom_point(data=df %>% filter(is.used2),aes(x=X_centroid,y=Y_centroid),col="grey90",size=.1)
            }

            p <- p + geom_point(size=.2)

            uniq.abs.tab <- table(df$ab)
            uniq.abs <- names(uniq.abs.tab[uniq.abs.tab >0])

            if( !plot_type %in% c("dna","exp")){
              if(!all(uniq.abs %in% names(uniq.cols))){
                stop("uniq.cols should have names that correspond to the labels")
              }else{
                uniq.cols <- uniq.cols[uniq.abs]
              }
            }


            p <- p + cols(uniq.cols)

            p <- p +
              ggtitle(ttl) +
              coord_fixed() +
              theme_void() +
              my_guides(txt=leg.ttl)

            if(draw.roi){
              prs <- x@rois
              polys <- lapply(seq(prs),function(i){
                roi <- prs[[i]]
                dir <- roi$dir
                coords <- roi$coords
                cbind(idx=i,dir=dir,coords)
              })
              polys1 <- as.data.frame(do.call(rbind,polys)) %>%
                mutate(dir = factor(dir,levels=c("positive","negative"))) %>%
                mutate(idx = factor(idx))

              p <- p +
                geom_polygon(data=polys1 %>% filter(dir == "positive"),aes(x=x,y=y,group=idx),color="red",fill=NA) +
                geom_polygon(data=polys1 %>% filter(dir == "negative"),aes(x=x,y=y,group=idx),color="blue",fill=NA)

            }

            print(p)
          }
)


#_ -------------------------------------------------------

# fun: plotUsedCellRatio Cycif,CycifStack ----

#' Plot Used Cell Ratio
#'
#' This function creates a plot of the ratio of used cells across channels in Cycif data.
#'
#' @param x A CycifStack object or a list of Cycif objects.
#' @param cumulative Should the cumulative ratio be plotted? (default is TRUE)
#' @param ncycle The number of cycles to consider for plotting (default is determined by the data).
#' @param mar Margins for the plot (default is c(5, 5, 4, 10)).
#' @param main Main title for the plot (default is "# cells attached on slide").
#' @param smpl.cols Colors for sample labels (default is generated using a color palette).
#' @param ncol Number of columns for the legend (default is 1).
#' @param leg.cex Size of the legend text (default is 0.8).
#' @param use_rois Should regions of interest (ROIs) be considered in the analysis? (default is TRUE)
#' @param ... Additional graphical parameters to customize the plot.
#'
#' @details
#' - The `plotUsedCellRatio` function creates a plot that visualizes the ratio of used cells in Cycif data.
#' - It calculates and plots the ratio of cells that are considered "used" based on the presence of regions of interest (ROIs).
#' - You can choose whether to display a cumulative ratio and specify the number of cycles to consider.
#' - The function also allows customization of various graphical parameters for the plot.
#'
#' @rdname plotUsedCellRatio
#' @export
setGeneric("plotUsedCellRatio",function(x,...) standardGeneric("plotUsedCellRatio"))

#' @rdname plotUsedCellRatio
#' @export
setMethod("plotUsedCellRatio", "Cycif", function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                                                 main="# cells attached on slide",...){
  x <- list2CycifStack(list(x))
  ret <- plotUsedCellRatio(x,cumulative=cumulative,ncycle=ncycle,mar=mar,main=main,...)
  return(invisible(ret))
})

#' @rdname plotUsedCellRatio
#' @export
setMethod("plotUsedCellRatio", "CycifStack",
          function(x,cumulative=TRUE,ncycle,mar=c(5,5,4,10),
                   main="# cells attached on slide",smpl.cols,ncol=1,
                   leg.cex=0.8,use_rois=TRUE,...){
  stopifnot(all(cyApply(x,function(x)is(x,"Cycif"),simplify=TRUE)))
  stopifnot(all(unlist(cyApply(x,function(cy)nrow(cy@used_cells))>0)))

  if(missing(ncycle)){
    ncycle <- unique(nCycles(x))
    if(length(ncycle)>1){
      stop("there are more than one value in 'ncycle'; provide the argument explicitly")
    }
  }

  used.ratio <- data.frame(sapply(names(x),function(n){
    y <- x[[n]]
    nc.ratio <- statUsedCells(y,use_rois=use_rois,cumulative=TRUE,ratio=TRUE)
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
#' Plot Available Cells on Slide
#'
#' This function creates a plot that visualizes the available cells on a slide in Cycif data.
#'
#' @param x A Cycif object.
#' @param upside.down Should the plot be upside-down? (default is TRUE)
#' @param ncycle The number of cycles to consider for plotting (default is determined by the data).
#' @param mfrow Number of rows and columns in the grid of plots (default is c(3, 3)).
#' @param mar Margins for each individual plot (default is c(0, 0, 4, 0)).
#' @param legend Should a legend be included in the plot? (default is TRUE)
#' @param main Main title for the plot (default is the names of the Cycif object).
#' @param cex.title Size of the main title text (default is 1).
#' @param uniq.cols Colors for different cell states (default colors).
#' @param legend.cex Size of the legend text (default is 2).
#' @param xlab Label for the x-axis (default is empty).
#' @param ylab Label for the y-axis (default is empty).
#' @param ... Additional graphical parameters to customize the individual plots.
#'
#' @details
#' - The `plotAvailCellOnSlide` function creates a grid of plots, each representing the available cells on a slide for a specific cycle.
#' - The plots show the spatial distribution of available cells in different colors, with an optional legend.
#' - You can customize the appearance of the plots and the legend using various graphical parameters.
#'
#' @export
setGeneric("plotAvailCellOnSlide",function(x,...) standardGeneric("plotAvailCellOnSlide"))

#' @rdname plotAvailCellOnSlide
#' @export
setMethod("plotAvailCellOnSlide", "Cycif",
          function(x,upside.down=TRUE,ncycle,mfrow=c(3,3),mar=c(0,0,4,0),legend=TRUE,main=names(x),cex.title=1,
                   uniq.cols=c(lost="grey90",dropped="blue",available="black",bunched="red"),legend.cex=2,
                   xlab="",ylab="",...){
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
                   xlab=xlab,ylab=ylab,...)
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
#' Violin Plots to Show Protein Expressions
#'
#' This function creates violin plots to visualize protein expressions in Cycif data.
#'
#' @param x A CycifStack object.
#' @param strat.by The strategy for stratifying the violin plots. Choose from "cell_type" or "smpl" (default is "cell_type").
#' @param ab The antibody or protein to plot.
#' @param use.pdata Should sample metadata be used for additional information? (default is FALSE).
#' @param fill.var The variable to use for filling the violin plots (default is "sample").
#' @param draw_thres Should the threshold be drawn on the plot? (default is FALSE).
#' @param type The type of data to use for plotting. Choose from "raw", "log", or "logTh" (default is "log").
#' @param strict Should strict cell type matching be enforced? (default is FALSE).
#' @param ct_name The name of the cell type column (default is "default").
#' @param ttl The title for the plot (default is determined based on inputs).
#' @param uniq.cts Unique cell types to include in the plot (default is all unique cell types).
#' @param uniq.smpls Unique samples to include in the plot (default is all samples).
#'
#' @details
#' - The `vlnPlot` function creates violin plots to visualize the protein expressions in Cycif data.
#' - You can stratify the plots by either cell types or samples using the `strat.by` parameter.
#' - Additional customization of the plot can be achieved using various graphical parameters.
#'
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes geom_violin position_dodge ggtitle coord_fixed theme_void
#'
#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @rdname vlnPlot
#' @export
setMethod("vlnPlot", "CycifStack",
          function(x,strat.by=c("cell_type","smpl"), ab="PDL1",
                   use.pdata=FALSE,fill.var,draw_thres=FALSE,
                   type=c("raw","log","logTh"),
                   strict=FALSE,ct_name="default",ttl,
                   uniq.cts,uniq.smpls){
            ## cell types
            cts <- cell_types(x,
                              strict=strict,
                              ct_name=ct_name)

            ucts <- levels(cts$cell_types)
            ucts <- ucts[ucts != "outOfROI"]

            cts <- cts %>%
              mutate(cell_types=factor(cell_types,levels=ucts))

            if(missing(type)){
              type <- "log"
            }

            if(missing(ab)){
              stop("ab should be always specified")
            }

            if(0){ # note gates are not plotted anymore
              ## gates
              if(type=="log"){
                gates.df <- cyApply(x,function(cy){
                  tmp <- abs_list(cy)
                  idx1 <- grep("gates",names(tmp))
                  if(length(idx1)>0){
                    tmp[[idx1]]
                  }else{
                    stop("no gates found on smpl:",names(cy))
                  }
                },simplify=T)
                rownames(gates.df) <- abs_list(x[[1]])$ab
              }else if(type=="logTh"){
                gates.df <- array(0.5,c(nrow(abs_list(x)),nSamples(x)))
                rownames(gates.df) <- abs_list(x)$ab
                colnames(gates.df) <- names(x)
              }

              ab_thres <- gates.df[ab,]
            }

            ## expression
            within.rois <- unlist(cyApply(x,function(cy)cy@within_rois))
            df <- exprs(x,type=type) %>%
              cbind(cts) %>%
              filter(within.rois)

            if(strat.by=="smpl"){
              ## if multiple samples should be shown, uniq.cts should be specified
              if(missing(uniq.cts)){
                uniq.cts <- ucts
                if(missing(ttl)){
                  ttl <- paste0(ab,", all cells")
                }
              }else{
                if(length(uniq.cts)>1){
                  uc <- paste0(length(uniq.cts)," cell types")
                }else if(length(uniq.cts) == 1){
                  uc <- uniq.cts
                }
                if(missing(ttl)){
                  ttl <- paste0(ab,", ",uc)
                }
              }
              ## uniq.smlps can be subsetted, or unspecified
              if(missing(uniq.smpls)){
                uniq.smpls <- names(x)
              }

              ## can add additional information
              if(use.pdata){
                pd <- pData(x) %>% dplyr::rename(smpl="id")
                df1 <- df %>% dplyr::left_join(pd,by="smpl")
              }else{
                df1 <- df
              }

              df1 <- df1 %>%
                dplyr::filter(sample %in% uniq.smpls) %>%
                dplyr::mutate(sample = factor(sample,levels=uniq.smpls)) %>%
                dplyr::filter(cell_types %in% uniq.cts) %>%
                # dplyr::mutate(thres=ab_thres[as.character(sample)]) %>%
                dplyr::filter(!is.na(!!!syms(ab)))
              return(df1)

              if(missing(fill.var)){
                fill.var <- "sample"
              }
              p <- ggplot(df1,aes(x=sample,y=!!sym(ab))) +
                geom_violin(position = position_dodge(width = 0.9),aes(fill=!!sym(fill.var)))
            }else if(strat.by=="cell_type"){ # not so useful?
              if(missing(uniq.smpls)||length(uniq.smpls)!=1){
                # uniq.smpls <- names(x)
                stop("specify one sample in 'uniq.smpls'")
              }else if(!uniq.smpls %in% names(x)){
                stop("'uniq.smpls' not defined in 'names(x)'")
              }
              if(missing(uniq.cts)){
                uniq.cts <- ucts
              }
              df1 <- df %>%
                dplyr::mutate(cell_types=factor(cell_types,levels=uniq.cts)) %>%
                dplyr::filter(sample %in% uniq.smpls) %>%
                dplyr::filter(!is.na(cell_types)) %>%
                dplyr::filter(!is.na(!!!syms(ab)))

              if(missing(ttl)){
                ttl <- paste0(ab,", ",uniq.smpls)
              }

              p <- ggplot(df1,aes(cell_types,!!sym(ab))) +
                geom_violin(aes(fill = cell_types)) +
                ggplot2::scale_fill_manual(values=ct.cols[1:14])
            }

            p <- p +
              ggtitle(ttl) +
              ylab(paste0("Expression (",type,")")) +
              theme_bw() +
              theme(legend.position = "none") +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            print(p)
          })

#_ -------------------------------------------------------

# fun: hist_1d numeric ----

#' 1D Histogram Plot
#'
#' This function creates a 1D histogram plot for a numeric vector.
#'
#' @param x A numeric vector for which the histogram is to be plotted.
#' @param n The number of bins or breaks for the histogram (default is 1000).
#' @param ths A vector of threshold values to mark on the plot (default is NULL).
#' @param mar A numerical vector of length 4 specifying the margin sizes (default is c(3, 4, 4, 2) + 0.1).
#' @param brks1 A vector of pre-defined breaks for the histogram (default is NULL).
#' @param ttl1 The title for the histogram plot (default is NULL).
#'
#' @details
#' - The `hist_1d` function creates a 1D histogram plot to visualize the distribution of a numeric vector.
#' - You can specify the number of bins with the `n` parameter or provide custom break points with the `brks1` parameter.
#' - Threshold values can be added to the plot using the `ths` parameter.
#' - Additional customization of the plot can be achieved using various graphical parameters.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats loess predict
#' @export
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

#_ -------------------------------------------------------

# fun: h3/heatmap3 matrix, data.frame ----

#' 3D Heatmap Plot (extension of heatmap3)
#'
#' This function creates a 3D heatmap plot to visualize a matrix of data, extending \code{heatmap3::heatmap3}.
#'
#' @param x A numeric matrix containing the data to be plotted.
#' @param Rowv Specifies the row dendrogram or clustering (default is NULL).
#' @param Colv Specifies the column dendrogram or clustering (default is NULL).
#' @param distfun A function to calculate the distance matrix for rows (default is dist).
#' @param distfunC A function to calculate the distance matrix for columns (default is NULL).
#' @param distfunR A function to calculate the distance matrix for rows (default is NULL).
#' @param balanceColor Logical, indicating whether to balance colors (default is FALSE).
#' @param ColSideLabs Labels for the columns' side (default is NULL).
#' @param RowSideLabs Labels for the rows' side (default is NULL).
#' @param showColDendro Logical, indicating whether to show the column dendrogram (default is TRUE).
#' @param showRowDendro Logical, indicating whether to show the row dendrogram (default is TRUE).
#' @param col A color palette for the heatmap (default is a gradient from navy to firebrick3).
#' @param legendfun A custom function to create the legend (default is NULL).
#'
#' @details
#' - The `h3` function creates a 3D heatmap plot to visualize a matrix of data, extending and customizing \code{heatmap3::heatmap3} function.
#' - It allows customization of various aspects of the plot, including colors, labels, and dendrograms.
#' - Additional functionalities like distance matrix calculations and custom legends are available.
#'
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
  # browser()
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

#_ -------------------------------------------------------

# fun: heatmapSummAb (Cycif),CycifStack ----

#' Create Heatmap Summary of Protein Expression
#'
#' This function generates a heatmap summarizing protein expression data from a `CycifStack` object.
#'
#' @param x A `CycifStack` object containing protein expression data.
#' @param norm_type Normalization type for the data, one of "log" or "logTh" (default is "log").
#' @param ab The name of the protein to be used for summarization.
#' @param sum_type The type of summarization to be performed, one of "freq", "mean", "median", or "x percentile" (e.g., "95 percentile").
#' @param ct_name The name of the cell type column (default is "default").
#' @param uniq_cts Vector of unique cell types to include in the heatmap.
#' @param uniq.smpls Vector of unique samples to include in the heatmap.
#' @param strict Logical, indicating strict filtering of cell types (default is FALSE).
#' @param scale Scaling method for the heatmap, one of "none", "row", or "column" (default is "none").
#'
#' @details
#' - The `heatmapSummAb` function creates a heatmap summarizing protein expression data for the specified protein.
#' - Users can choose from different normalization types and summarization methods.
#' - Cell types and samples can be filtered and customized for the heatmap.
#'
#' @export
#' @rdname heatmapSummAb
setGeneric("heatmapSummAb", function(x,...) standardGeneric("heatmapSummAb"))

#' @rdname heatmapSummAb
#' @export
setMethod("heatmapSummAb", "CycifStack",
          function(x,norm_type=c("log","logTh"),
                   summBy=c("ab","heatmap"),ab,
                   sum_type=c("freq","mean","median","x percentile"),
                   ct_name="default",
                   uniq_cts,
                   uniq.smpls,
                   strict=FALSE,
                   p_thres=0.5,
                   scale="none",...){
    options(dplyr.summarise.inform = FALSE)

    if(missing(sum_type)){
      stop("'sum_type' should be specified")
    }
    if(missing(ab)){
      stop("ab should be always specified")
    }

    # norm_type
    if(sum_type=="freq"){
      norm_type="logTh"
    }else if(missing(norm_type)){
      norm_type <- "log"
    }

    # sum_fun
    if(grepl("percentile$",sum_type)){
      pct <- strsplit(sum_type," ")[[1]]
      p.ile <- as.numeric(pct[[1]])
      sum_fun=function(y,...)quantile(y,p.ile/100,na.rm=T)
    }else if(sum_type=="median"){
      sum_fun=function(y,...)quantile(y,.5,na.rm=T)
    }else if(sum_type=="mean"){
      sum_fun=function(y,...)mean(y,na.rm=T)
    }else if(sum_type=="freq"){
      sum_fun=function(y,th){
        if(length(y)==0){
          return(NA)
        }else{
          return(sum(y > th, na.rm=T)/sum(!is.na(y)))
        }
      }
    }

    # balanceColor
    if(sum_type=="freq" || norm_type=="log"){
      balanceColor <- FALSE
    }else if(norm_type=="logTh"){
      balanceColor <- TRUE
    }

    ## cell types
    cts <- cell_types(x,ct_name=ct_name,strict=strict)

    df <- exprs(x,type=norm_type) %>%
      cbind(cts) %>%
      filter(cell_types != "outOfROI") %>%
      rename(smpl = sample) %>%
      rename(cell_type = cell_types)

    if(missing(uniq.smpls)){
      uniq.smpls <- levels(df$smpl)
    }

    if(missing(uniq_cts)){
      uniq_cts <- levels(df$cell_type)
    }

    df <- df %>%
      mutate(smpl = factor(smpl,levels=uniq.smpls)) %>%
      filter(smpl %in% uniq.smpls) %>%
      mutate(cell_type=factor(cell_type,levels=uniq_cts)) %>%
      mutate(cell_type=factor(cell_type)) %>%
      filter(!is.na(cell_type)) %>%
      dplyr::rename(this_ab=!!sym(ab)) %>%
      group_by(smpl,cell_type) %>%
      dplyr::summarise(sum_ab=sum_fun(this_ab,th=p_thres))

    ttl <- paste0(ab,",",sum_type,",",norm_type,",(scale:",scale,")")

    m1 <- as.matrix(with(df, Matrix::sparseMatrix(as.integer(smpl), as.integer(cell_type), x=sum_ab)))

    m1[apply(m1,1,function(x)any(is.nan(x))),] <- NA
    if(sum_type!="freq"){
      m1[m1==0] <- NA
    }

    rownames(m1) <- levels(df$smpl)
    colnames(m1) <- levels(df$cell_type)
    m1 <- m1[rev(seq(nrow(m1))),]

    if(sum_type!="freq" && norm_type=="logTh"){
      m1 <- m1 - 0.5 ## hardcode!!!!
    }

    ## cols
    if(sum_type=="mean"){
      cols <- viridis::viridis(1024)
    }else if(sum_type == "freq"){
      cols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(1024)
    }
    h3(m1,
       margins=c(10,5),
       main=ttl,
       cexRow=.7,
       na.rm = T,
       col = cols,
       balanceColor=balanceColor,
       scale=scale,...)

    invisible(df)
  })

#_ -------------------------------------------------------

# fun: LdPlot LDCoords, Cycif, CycifStack ----

#' Create Low-Dimensional Plot
#'
#' This function generates low-dimensional plots for protein expression data from a `Cycif` or `CycifStack` object.
#'
#' @param x A `Cycif` or `CycifStack` object containing protein expression data.
#' @param ld_name The name of the low-dimensional representation (e.g., UMAP or clustering).
#' @param plot_type The type of low-dimensional plot to create, one of "cell_type", "clusters", "exp", or "smpl".
#' @param ab The name of the protein to be used for coloring points (required when `plot_type` is "exp").
#' @param uniq.cols Vector of unique colors to use for plotting points (optional).
#' @param with.labels Logical, indicating whether to label points (default is TRUE).
#' @param ct_name The name of the cell type column (default is "default").
#' @param used.smpls Vector of samples to include in the plot.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param used.cts Vector of cell types to include in the plot.
#' @param cex.main Size of the main title text (default is 2).
#' @param main Title for the plot (default is auto-generated based on parameters).
#' @param p_thres Threshold for plotting (default is 0.5).
#' @param mar Margins of the plot.
#' @param cex.labs Size of label text (default is 1).
#' @param cex.leg Size of the legend text (default is 0.5).
#' @param cex Size of data point labels (default is 0.3).
#'
#' @details
#' - The `LdPlot` function creates low-dimensional plots for protein expression data based on the specified `Cycif` or `CycifStack` object.
#' - Users can choose from different types of low-dimensional plots, including those based on cell types, clusters, protein expression, or samples.
#' - Various customization options are available for labeling, coloring, and scaling the plot.
#'
#' @rdname LdPlot
#' @export
setGeneric("LdPlot", function(x,...) standardGeneric("LdPlot"))

#' @rdname LdPlot
#' @export
setMethod("LdPlot", "LDCoords",
          function(x,
                   ld_name,
                   plot_type=c("cell_type","clusters","exp","smpl","meta","custom"),
                   ab, # exp
                   meta.var, # meta
                   vals, # custom

                   uniq.cols,with.labels=TRUE,
                   used.smpls,xlab="",ylab="",
                   used.cts,cex.main=2,
                   balanceColor=FALSE,
                   main,p_thres=0.5,mar,cex.labs=1,cex.leg=.5,
                   cex=.3,...){
            if(missing(plot_type)){
              stop("'plot_type' should be specified. (one of cell_type, clusters, exp, smpl, custom)")
            }

            ## xy coords ----
            xys <- x@ld_coords %>%
              dplyr::rename(smpl=sample) %>%
              dplyr::rename(cell_type=cell_types)

            if(plot_type=="exp"){
              xys <- xys %>%mutate(exp=vals)
              ptype <- paste0(ab," exp")

            }else if(plot_type=="cell_type"){
              if(!missing(used.cts)){
                if(any(!used.cts %in% levels(xys$cell_type))){
                  stop("used.cts contain cell types not used")
                }else{
                  xys$cell_type <- factor(xys$cell_type,levels=used.cts)
                }
              }

              ptype <- "Cell types"

            }else if(plot_type=="clusters"){
              # stop(paste(slotNames(x),collapse="\n"))
              if(length(x@clusters)==0){
                stop("first run LdClustering() to plot clusters")
              }else if(length(x@clusters)!=nrow(x@ld_coords)){
                stop("x@clusters doesn't have the same # of clusters as the cells")
              }else{
                # stop(paste(names(xys),collapse=", "))
                # stop(paste(dim(xys),collapse=","))
                clusts <- x@clusters
                # stop(paste(head(clusts),collapse=","))
                xys <- xys %>% mutate(clusters=clusts)
              }

              ptype <- "Clusters"

            }else if(plot_type=="smpl"){
              xys$smpl <- factor(xys$smpl)
              if(!missing(used.smpls)){
                if(!any(used.smpls %in% levels(xys$smpl))){
                  stop("used.smpls contain samples that don't exist")
                }
              }else{
                used.smpls <- levels(xys$smpl)
              }

              xys$smpl  <- factor(xys$smpl,levels=used.smpls)

              ptype <- "Samples"

            }else if(plot_type=="meta"){
              xys <- xys %>% mutate(meta=vals)
              ptype <- meta.var
              # return(xys)
            }else if(plot_type=="custom"){
              if(missing(vals)){
                stop("'vals' should be specified when 'plot_type' is 'custom'")
              }else{
                if(length(vals)!=nrow(x@ld_coords)){
                  stop("vals doesn't have the same # of values as the cells")
                }
              }
              xys <- xys %>% mutate(custom=vals)

              # stop(names(xys))

              ptype <- "custom"
            }else{
              stop("plot_type should be one of cell_type, cluster, exp, smpl, meta, custom")
            }

            ## main
            if(missing(main)){
              main <- paste(ld_name,ptype,x@ld_type,sep=", ")
            }

            ## colors, with.labels
            if(is.null(xys[[plot_type]])){
              stop("xys[[plot_type]] is null")
            }else{
              vals <- xys[[plot_type]]
              if(is(vals,"numeric")){
                if(plot_type=="custom"){
                  if(missing(main)){
                    stop("specify plot title in 'main'")
                  }else{
                    ab <- main
                  }
                }
                if(missing(ab)){
                  ab <- meta.var
                }
                my_guides <- guides(color=guide_colorbar(title=ab))
                if(balanceColor==TRUE){
                  cols <- scale_color_distiller(palette="RdBu",direction=-1,
                                                limits = c(-1,1)*max(abs(vals)))
                }else{
                  rng <- range(vals,is.na=T)
                  if(prod(rng)<0){
                    cols <- scale_color_distiller(palette="RdBu",direction=-1)
                  }else{
                    cols <- scale_color_distiller(palette="YlGnBu",direction=1)
                  }
                }
                with.labels <- FALSE
              }else if(is(vals,"factor")){
                my_guides <- guides(color=guide_legend(title=plot_type,override.aes = list(size = cex.main)))
                if(!missing(uniq.cols)){
                  tab <- table(vals)
                  names(tab[tab>0])
                  uc <- uniq.cols[names(tab[tab>0])]
                  cols <- scale_color_manual(values=uc)
                }else{
                  ## viridis color palette for factor
                  # cols <- scale_color_viridis_d(direction=1)
                  # uniq.cols <- rainbow(nlevels(xys$clusters))
                  # names(uniq.cols) <- levels(xys$clusters)

                  num_classes <- length(unique(vals))
                  spectral_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(num_classes)
                  cols <- scale_color_manual(values = spectral_colors)
                  # stop("uniq.cols should be provided when values are a factor")
                }
              }
            }

            p <- ggplot(xys,aes(x=x,y=y)) +
              geom_point(aes(color=!!sym(plot_type)),size=cex) +
              ggtitle(main) +
              xlab(xlab)+
              ylab(ylab)+
              theme_void() +
              theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
              my_guides +
              cols

            if(with.labels){
              xys1<- xys %>%
                select(x,y,!!sym(plot_type)) %>%
                group_by(!!sym(plot_type)) %>%
                summarize_at(c("x","y"),median)

              p <- p + geom_text_repel(data=xys1,aes(x=x,y=y,label=!!sym(plot_type)))
            }
            # print(p)
            # return(xys)
            return(p)
          }
)


#' @rdname LdPlot
#' @export
setMethod("LdPlot", "Cycif",
  function(x,
           ld_name,
           plot_type=c("cell_type","clusters","exp","smpl","meta","custom"),
           ab, # exp
           meta.var, # meta
           vals, # custom

           uniq.cols,with.labels=TRUE,
           used.smpls,xlab="",ylab="",
           used.cts,cex.main=2,
           main,mar,cex.labs=1,cex.leg=.5,
           cex=.3,...){
    if(missing(plot_type)){
      stop("'plot_type' should be specified. (one of cell_type, cluster, exp, smpl, meta, custom)")
    }

    ld <- x@ld_coords[[ld_name]]
    if(plot_type=="meta"){
      stop("plot_type='meta' is not supported for Cycif object")
    }else if(plot_type=="exp"){
      if(missing(ab)){
        stop("'ab' should be specified when 'plot_type' is 'exp'")
      }

      n <- exprs(x,type=ld@norm_type)[ld@is_used,]

      if(!ab %in% names(n)){
        stop(paste0(ab," is not included in the exprs(x)"))
      }

      vals <- n[[ab]]
    }
    LdPlot(ld,
           ld_name=ld_name,
           plot_type=plot_type,
           ab=ab,
           meta.var=meta.var,
           vals=vals,
           uniq.cols=uniq.cols,
           with.labels=with.labels,
           used.smpls=used.smpls,
           xlab=xlab,ylab=ylab,
           used.cts=used.cts,
           cex.main=2,
           main=main,
           mar=mar,
           cex.labs=cex.labs,
           cex.leg=cex.leg,
           cex=cex,...)
  }
)

#' @rdname LdPlot
#' @export
setMethod("LdPlot", "CycifStack",
    function(x,
             ld_name,
             plot_type=c("cell_type","clusters","exp","smpl","meta","custom"),
             ab, # exp
             meta.var, # meta
             vals, # custom

             uniq.cols,with.labels=TRUE,
             used.smpls,xlab="",ylab="",
             used.cts,cex.main=2,
             main,p_thres=0.5,mar,cex.labs=1,cex.leg=.5,
             cex=.3,...){
    if(missing(plot_type)){
      stop("'plot_type' should be specified. (one of cell_type, cluster, exp, smpl, meta, custom)")
    }

    ld <- x@ld_coords[[ld_name]]
    if(plot_type == "meta"){
      pd <- pData(x)
      if(missing(meta.var)){
        stop("'meta.var' should be specified when 'plot_type' is 'meta'")
      }else{
        if(!meta.var %in% names(pd)){
          stop(paste0(meta.var," is not included in the pData(x)"))
        }else{
          vals <- pd[[meta.var]][match(ld@ld_coords$sample,pd$id)]
        }
      }
    }else if(plot_type=="exp"){
      if(missing(ab)){
        stop("'ab' should be specified when 'plot_type' is 'exp'")
      }

      n <- exprs(x,type=ld@norm_type)[ld@is_used,]

      if(!ab %in% names(n)){
        stop(paste0(ab," is not included in the exprs(x)"))
      }

      vals <- n[[ab]]
    }
    LdPlot(ld,
           ld_name=ld_name,
           plot_type=plot_type,
           ab=ab,
           meta.var=meta.var,
           vals=vals,
           uniq.cols=uniq.cols,
           with.labels=with.labels,
           used.smpls=used.smpls,
           xlab=xlab,ylab=ylab,
           used.cts=used.cts,
           cex.main=2,
           main=main,
           mar=mar,
           cex.labs=cex.labs,
           cex.leg=cex.leg,
           cex=cex,...)
  }
)
#_ -------------------------------------------------------

# Following functions are unique to TALAVE project ----

## fun: barplotCTS ----

#' Create a Barplot for Cell Type Composition (TALAVE)
#'
#' This function generates a composite barplot to visualize the composition of cell types across samples.
#'
#' @param n.sh A matrix or data frame containing cell type composition data, with cell types as rows and samples as columns.
#' @param anno.smpls A data frame containing sample annotations, with sample names matching the column names of `n.sh`.
#' @param uniq.cols A vector of unique colors to use for plotting cell types.
#' @param ylim A numeric vector specifying the limits of the y-axis (default is c(0, 1)).
#' @param ylab Label for the y-axis (default is "CellType Composition").
#'
#' @details
#' The `barplotCTS` function creates a composite barplot with the following features:
#'
#' - Visualizes the composition of various cell types across multiple samples.
#' - Each bar in the plot represents a sample, and the height of the bar is proportional to the composition of cell types within that sample.
#' - The cell types are color-coded using unique colors specified in the `uniq.cols` parameter for easy identification.
#' - Sample annotations, such as patient ID, disease subtype, time point, and others, are displayed below the main barplot to provide additional context.
#' - The y-axis represents the composition proportion of cell types, and you can customize the y-axis label using the `ylab` parameter.
#' - The `ylim` parameter allows you to set specific limits for the y-axis, controlling the range of the composition proportions displayed.
#'
#' This function is useful for visualizing and comparing cell type compositions across different samples, making it particularly valuable in biological and medical research for understanding the distribution of cell types in complex datasets.
#'
#' @export
barplotCTS <- function(n.sh,anno.smpls,uniq.cols,ylim=c(0,1),ylab="CellType Composition"){
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
  image(xs,1,as.matrix(as.numeric(anno.smpls$TimePoint)),col=brewer.pal(3,"Set2"),axes=F,
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
  med.x <- tapply(seq(nrow(anno.smpls)),
                  factor(anno.smpls$Patient.ID,levels=unique(anno.smpls$Patient.ID)),median)
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
  legend(par()$usr[2],.55*ymax,c("MUT","WT"),fill=c("grey20","grey80"),title="gBRCA")
  legend(par()$usr[2],.375*ymax,levels(anno.smpls$Subtype),fill=brewer.pal(3,"Dark2"),title="Subtype")
  legend(par()$usr[2],0.2*ymax,sub("CD68_CD163_","DP_",rownames(n.sh)),fill=uniq.cols,title="Cell types")
  par(xpd=F)
  par(fg=NA)
  axis(1,at=xidx,labels=names(med.x),las=2,cex.axis=.8)
  par(fg="black")

  d <- diff(xx[1:2])
  par(xpd=T)
  pts <- factor(anno.smpls$Patient.ID,levels=unique(anno.smpls$Patient.ID))
  vs <- c(xx[1]-d/2,xx[cumsum(table(pts))] + d/2)
  sapply(vs,function(x1){
    lines(rep(x1,2),c(-0.05,1.4),lty=2)
  })
  par(xpd=F)

}

## fun: createRSC ----

# compute RowSideColors of heatmap3(), or its modified version, h3().

#' Create RowSideColors for Heatmap Visualization (TALAVE)
#'
#' This function generates RowSideColors (RSC) for use in heatmap visualization, typically with the `heatmap3` or `h3` functions. RSC assigns colors to rows (e.g., samples or features) based on sample annotations, such as Patient ID, gBRCA mutation status, time points, BOR (type of breast cancer), and clusters (for dimension reduction analysis).
#'
#' @param cs A `CycifStack` object or other data structure containing sample annotations used for color assignment.
#' @param with.clusters Logical value indicating whether to include clustering information in the RowSideColors (default is TRUE). If TRUE, the function expects a `ld_name` to specify the clustering to be used.
#' @param ld_name A character string specifying the name of the UMAP or clustering to be used when `with.clusters` is TRUE. Required if `with.clusters` is TRUE.
#' @param order A character string specifying the desired order of the RowSideColors. The order can be a combination of the following letters: "B" for BOR, "c" for clusters, "T" for TimePoint, "g" for gBRCA.status, and "P" for Patient.ID. This parameter allows you to customize the order of the RowSideColors according to your preference (default is NULL).
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{pd}: A data frame with row names and the specified annotations.
#'   \item \code{rsc}: A matrix representing the RowSideColors for the heatmap.
#'   \item \code{row.order}: The order of rows in the RSC matrix if custom ordering is applied (only returned if the \code{order} parameter is provided).
#' }
#'
#' @details
#' The `createRSC` function generates RowSideColors for heatmap visualization based on sample annotations. The function assigns colors to rows (samples or features) in the heatmap, making it easier to identify patterns and relationships between samples.
#'
#' The available annotations for color assignment are:
#' - Patient ID
#' - gBRCA mutation status (MUT or WT)
#' - Time point
#' - BOR (type of breast cancer)
#' - Clusters (dimension reduction clusters, requires specifying `ld_name` when `with.clusters` is TRUE)
#'
#' You can customize the order of the RowSideColors using the \code{order} parameter, which allows you to specify the order of annotations to be displayed in the heatmap from left to right. For example, setting \code{order = "BPcgT"} will order the annotations as BOR, Patient ID, clusters, gBRCA status, and TimePoint.
#'
#' This function is especially useful for visualizing multi-dimensional data, such as cytometry data, in a heatmap, allowing you to explore relationships between samples and their associated annotations.
#'
#' @export
createRSC <- function(cs,with.clusters=T,ld_name,order){
  pd <- pData(cs) %>%
    rename(smpl=id) %>%
    filter(smpl %in% names(cs))

  col1 <- RColorBrewer::brewer.pal(length(unique(pd$Patient.ID)),"Set3")
  names(col1) <- unique(pd$Patient.ID)
  col2 <- c("grey20","grey80")
  names(col2) <- c("MUT","WT")
  col3 <- RColorBrewer::brewer.pal(length(unique(pd$TimePoint)),"Set2")
  names(col3) <- levels(pd$TimePoint)
  col4 <- RColorBrewer::brewer.pal(length(unique(pd$BOR)),"Set1")[c(1,3,2)]
  names(col4) <- levels(pd$BOR)
  col5 <- RColorBrewer::brewer.pal(length(unique(pd$BOR)),"Dark2")[1:2]
  names(col5) <- levels(pd$Subtype)

  ## compile rsc (RowSideColors) matrix
  if(with.clusters){
    if(missing(ld_name)){
      stop("when 'with.clusters' is TRUE, 'ld_name' shoud be provided")
    }
    ld <- ld_coords(cs,ld_name)
    clusts <- ld@clusters
    ncls <- nlevels(clusts)
    nsmpls <- ld@n_cells_total
    smpls <- rep(names(nsmpls),nsmpls)
    pd <- data.frame(
      clusters = clusts,
      smpl=smpls
    ) %>% left_join(pd,by="smpl")

    set.seed(123)
    if(ncls < 9){
      col5 <- RColorBrewer::brewer.pal(8,"Spectral")[seq(ncls)]
    }else{
      col5 <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(ncls)
    }

    names(col5) <- unique(pd$clusters)
    rsc <- cbind(
      Patient.ID = col1[pd$Patient.ID],
      gBRCA.status = col2[pd$gBRCA.status],
      TimePoint = col3[pd$TimePoint],
      BOR = col4[pd$BOR],
      clusters = col5[pd$clusters])
  }else{
    rsc <- cbind(
      Patient.ID = col1[pd$Patient.ID],
      gBRCA.status = col2[pd$gBRCA.status],
      TimePoint = col3[pd$TimePoint],
      BOR = col4[pd$BOR])
  }

  pd <- pd %>% tibble::rownames_to_column("idx")

  ## reorder - specific to TALAVE project
  if(!missing(order)){
    o1 <- strsplit(order,"")[[1]]
    o2 <- sapply(o1,function(x){
      if(x=="B"){
        x <- "BOR"
      }else if(x=="c"){
        x <- "clusters"
      }else if(x=="T"){
        x <- "TimePoint"
      }else if(x=="g"){
        x <- "gBRCA.status"
      }else if(x=="P"){
        x <- "Patient.ID"
      }
      return(x)
    })

    row.o <- rev(as.numeric((pd %>% arrange(!!!rlang::syms(o2)))$idx))
    col.o <- rev(match(o1,sub("^(.).+","\\1",colnames(rsc))))
    rsc1 <- rsc[row.o,col.o,drop=F]
    pd1 <- pd[row.o,,drop=F]
    return(list(pd=pd1,rsc=rsc1,row.order=row.o))
  }else{
    return(list(pd=pd,rsc=rsc))
  }
}
