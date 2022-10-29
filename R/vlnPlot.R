#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' #' @export
#' setMethod("vlnPlot", "CycifStack",
#'   function(x,ab,type=c("raw","normalized"),...){
#'     require(tidyr)
#'     require(ggplot2)
#'
#'     if(missing(type)){
#'       type <- "normalized"
#'     }
#'
#'     nms <- names(x)
#'     n <- exprs(x,type=type)
#'     thres <- x@threshold[ab,nms]
#'     if(type=="raw"){
#'       p <- ggplot(n,aes_string(x="smpl",y=ab),...) +
#'         ggtitle(ab) +
#'         geom_violin(scale="width") +
#'         labs(y = paste0("expression (raw)")) + xlab("") +
#'         theme(plot.title = element_text(hjust = 0.5))
#'       for(nm in nms){
#'         xidx <- match(nm,nms) + c(-.4,.4)
#'         p <- p + annotate("segment", x = xidx[1], y = thres[[nm]], xend = xidx[2], yend = thres[[nm]],col=2)
#'       }
#'     }else if(type=="normalized"){
#'       mth <- x@normalize.method
#'       if(mth == "log"){
#'         m <- "log"
#'       }else if(mth == "logTh"){
#'         m <- "log + threshold"
#'       }
#'       p <- ggplot(n,aes_string(x="smpl",y=ab),...) +
#'         ggtitle(ab) +
#'         geom_violin(scale="width") +
#'         labs(y = paste0("expression (",m,")")) + xlab("") +
#'         theme(plot.title = element_text(hjust = 0.5))
#'       if(x@normalize.method=="log"){
#'         log.thres <- log1p(thres)
#'         for(nm in nms){
#'           xidx <- match(nm,nms) + c(-.4,.4)
#'           p <- p + annotate("segment", x = xidx[1], y = log.thres[[nm]], xend = xidx[2], yend = log.thres[[nm]],col=2)
#'         }
#'       }else if(mth=="logTh"){
#'         p <- p + geom_hline(yintercept = 0.5,col=2)
#'       }
#'     }
#'     p
#'   }
#' )

#' @export
setMethod("vlnPlot", "CycifStack",
  function(x,type=c("raw","normalized"), xaxis=c("celltype","smpl"),
           ab="PDL1", cts, ctype="all", sample="all",ttl){
    require(tidyr)
    require(ggplot2)

    # x = ln
    # type="normalized"
    # xaxis="celltype"
    # celltype=levels(dc2)
    # # cts=ctypes
    # cts = dc2
    # ab=used.abs[1]

    if(missing(type)){
      type <- "normalized"
    }

    nms <- names(x)
    df <- exprs(x,type=type) %>%
      mutate(celltype=factor(cts))
    df <- df %>%
      mutate(smpl = factor(sub("\\.[0-9]+$","",rownames(df))))

    ab_thres <- log1p(x@threshold[ab,nms])

    df1 <- df

    if(ctype[1] != "all"){
      df1 <- df1 %>%
        mutate(celltype=factor(celltype,levels=ctype)) %>%
        filter(!is.na(celltype))
      # stop(length(levels(df1$celltype)))
    }
    if(sample[1] != "all"){
      df1 <- df1 %>%
        filter(smpl %in% sample) %>%
        mutate(smpl=factor(smpl))
    }
    if(missing(ttl)){
      ttl <- paste0(ab)
    }
    if(xaxis=="smpl"){
      df_thres <- data.frame(smpl=names(ab_thres),th=unlist(ab_thres))
      p <- ggplot(df1,aes_string("smpl",paste0("`",ab,"`"))) +
        geom_violin(aes(fill = sample)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(ttl)
      # for(smpl in smpls){
      #   i <- match(smpl,smpls)
      #   amt <- 0.2
      #   th <- ab_thres[[smpl]]
      #   p <- p +  geom_segment(aes(x=i-amt/2,xend=i+amt/2,y=th,yend=th))
      # }

      # stop(names(df_thres))
      print(p)
    }else if(xaxis=="celltype"){ # not so useful?
      p <- ggplot(df1,aes_string("celltype",paste0("`",ab,"`"))) +
        geom_hline(yintercept = ab_thres[[smpl]]) +
        geom_violin(aes(fill = celltype)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(ttl) +
        scale_fill_discrete(drop=FALSE) +
        scale_x_discrete(drop=FALSE)
      print(p)
    }
})
