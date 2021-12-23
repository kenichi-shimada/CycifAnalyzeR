#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @export
setMethod("vlnPlot", "Cycif",
  function(x,type=c("raw","normalized"),...){
    require(tidyr)
    require(ggplot2)
    if(missing(type)){
      type <- "raw"
    }

    smpl <- names(x)
    thres <- x@threshold
    used.abs <- names(thres)

    n <- exprs(x,type=type)[,used.abs]
    df <- n %>% gather(key="channel",value="exprs")
    df$channel <- factor(df$channel,levels=abs)

    if(type=="raw"){
      p <- ggplot(df,aes(channel,exprs)) + geom_violin(scale="width") +
        ggtitle(smpl) +
        labs(y = paste0("expression (raw)")) + xlab("") +
        theme(plot.title = element_text(hjust = 0.5))
      for(ab in abs){
        if(ab %in% names(thres)){
          xidx <- match(ab,abs) + c(-.4,.4)
          p <- p + annotate("segment", x = xidx[1], y = thres[[ab]], xend = xidx[2], yend = thres[[ab]],col=2)
        }
      }
    }else{
      mth <- x@normalize.method
      if(mth == "log"){
        m <- "log"
      }else if(mth == "logTh"){
        m <- "log + threshold"
      }
      p <- ggplot(df,aes(channel,exprs)) + geom_violin(scale="width") +
        ggtitle(smpl) +
        labs(y = paste0("expression (",m,")")) + xlab("") +
        theme(plot.title = element_text(hjust = 0.5))

      if(mth=="log"){
        log.thres <- log1p(thres)
        for(ab in abs){
          if(ab %in% used.abs){
            xidx <- match(ab,abs) + c(-.4,.4)
            p <- p + annotate("segment", x = xidx[1], y = log.thres[[ab]], xend = xidx[2], yend = log.thres[[ab]],col=2)
          }
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
