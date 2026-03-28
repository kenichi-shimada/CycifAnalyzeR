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
