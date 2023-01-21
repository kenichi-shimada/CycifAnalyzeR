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
           pch=".",main,p_thres=0.5,mar,cex.labs=1,...){
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
    plot(xys,col=cols,pch=pch,main=main,xlab="",ylab="",...)


    if(leg){
      par(xpd=T)
      legend(par()$usr[2],par()$usr[4],levs,pch=pch,col=uniq.cols,cex=.9)
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
                   cts.names,pal.name,
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
                                 leaves.only=cts.params$leaves.only,
                                 strict=cts.params$strict)[is.used]
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
            plot(xys,col=cols,pch=pch,main=main,xlab="",ylab="",cex=cex,...)


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
