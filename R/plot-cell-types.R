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
