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
setGeneric("dimPlot", function(x,...) standardGeneric("dimPlot"))

#' @rdname dimPlot
#' @export
setMethod("dimPlot", "CycifStack",
          function(x,pch=".",type=c("cell_type","smpl","exp"),celltype=cts,ab,uniq.cols,main,...){
            require(RColorBrewer)
            ld <- x@ld_coords
            rn <- rownames(ld)
            if(missing(uniq.cols)){
              uniq.cols <- brewer.pal(11,"Spectral")
            }
            if(type=="cell_type"){
              dc <- do.call(c,cts)
              uniq.cols <- c(RColorBrewer::brewer.pal(11,"Spectral"),"grey80","black")
              cols <- uniq.cols[dc[rn]]
              ttl <- "Cell type"
            }else if(type=="smpl"){
              smpls <- factor(sub("\\..+","",rn))
              # stop(table(smpls))
              cols <- as.numeric(smpls)+1
              ttl <- "Samples"
              pts <- levels(smpls)
            }else if(type=="exp"){
              stopifnot(!missing(ab))
              n <- x@normalized
              n <- n[rn,ab]

              rng <- range(n,na.rm=T)
              nn <- (n - rng[1])/diff(rng)

              nc <- 100
              uniq.cols <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
              cols <- colorRampPalette(uniq.cols)(nc+1)[round(nn*nc)+1]
              ttl <- paste0(ab, "(expression)")
            }
            if(missing(main)){
              main <- ttl
            }

            par(mar=c(4,4,4,10))
            plot(ld,col=cols,pch=".",main=main)
            par(xpd=F)
            legend(par()$usr[2],par()$usr[4],)
            par(xpd=T)
          }
)
