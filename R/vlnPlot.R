#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @export
setMethod("vlnPlot", "CycifStack",
  function(x,ab,type=c("raw","normalized"),...){
    require(tidyr)
    require(ggplot2)

    if(missing(type)){
      type <- "normalized"
    }

    nms <- names(x)
    n <- exprs(x,type=type)
    thres <- x@threshold[ab,nms]
    if(type=="raw"){
      p <- ggplot(n,aes_string(x="smpl",y=ab),...) +
        ggtitle(ab) +
        geom_violin(scale="width") +
        labs(y = paste0("expression (raw)")) + xlab("") +
        theme(plot.title = element_text(hjust = 0.5))
      for(nm in nms){
        xidx <- match(nm,nms) + c(-.4,.4)
        p <- p + annotate("segment", x = xidx[1], y = thres[[nm]], xend = xidx[2], yend = thres[[nm]],col=2)
      }
    }else if(type=="normalized"){
      mth <- x@normalize.method
      if(mth == "log"){
        m <- "log"
      }else if(mth == "logTh"){
        m <- "log + threshold"
      }
      p <- ggplot(n,aes_string(x="smpl",y=ab),...) +
        ggtitle(ab) +
        geom_violin(scale="width") +
        labs(y = paste0("expression (",m,")")) + xlab("") +
        theme(plot.title = element_text(hjust = 0.5))
      if(x@normalize.method=="log"){
        log.thres <- log1p(thres)
        for(nm in nms){
          xidx <- match(nm,nms) + c(-.4,.4)
          p <- p + annotate("segment", x = xidx[1], y = log.thres[[nm]], xend = xidx[2], yend = log.thres[[nm]],col=2)
        }
      }else if(mth=="logTh"){
        p <- p + geom_hline(yintercept = 0.5,col=2)
      }
    }
    p
  }
)

#' @export
setMethod("vlnPlot", "CycifStack",
  function(x,type=c("raw","normalized"), xaxis=c("celltype","smpl"),
           ab="PDL1", cts, cell_type="all", sample="all", hline=0.4,ttl){
    require(tidyr)
    require(ggplot2)

    # x = ln
    # type="normalized"
    # xaxis="celltype"
    # cell_type=levels(dc)[1:15]
    # # cts=ctypes
    # cts = dc
    # ab="HER2"

    if(missing(type)){
      type <- "normalized"
    }

    nms <- names(x)
    df <- exprs(x,type=type) %>% mutate(celltype=factor(cts))
    thres <- x@threshold[ab,nms]

    df1 <- df
    if(cell_type[1] != "all"){
      df1 <- df1 %>%
        filter(celltype %in% cell_type) %>%
        mutate(celltype=factor(celltype))
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
      p <- ggplot(df1,aes_string("smpl",ab)) +
        geom_hline(yintercept = p_thres) +
        geom_violin(aes(fill = sample)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(ttl)
      print(p)
    }else if(xaxis=="cell_type"){ # not so useful?
      p <- ggplot(df1,aes_string("celltype",ab)) +
        geom_hline(yintercept = p_thres) +
        geom_violin(aes(fill = celltype)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(ttl) +
        scale_fill_discrete(drop=FALSE) +
        scale_x_discrete(drop=FALSE)
      print(p)
    }
})
