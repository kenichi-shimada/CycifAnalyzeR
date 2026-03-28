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

            ## cell types - specified by `uniq.cts`
            if(ct_name %in% names(x@cell_types)){
              cts <- cell_types(x,strict=strict,ct_name=ct_name)

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
                cbind(data.frame(ab=custom_labs,is.na=is.na(custom_labs))) %>%
                arrange(rev(ab))
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
              my_guides(txt=leg.ttl) +
              scale_y_reverse()

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
