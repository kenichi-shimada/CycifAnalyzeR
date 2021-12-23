#' @export
setGeneric("nrPlot", function(x,...) standardGeneric("nrPlot"))
setMethod("nrPlot", "Cycif",
          function(x,ab,xlog=FALSE,p_thres=0.5,main,...){
            stopifnot(!missing(ab))
            smpl <- names(x)
            mth <- x@normalize.method
            uabs <- names(x@threshold)
            th <- x@threshold

            stopifnot(ab %in% uabs)

            r <- exprs(x,type="raw")[[ab]]
            n <- exprs(x,type="normalized")[[ab]]

            if(missing(main)){
              main <- paste0(smpl," (",mth,")")
            }

            if(mth=="logTh"){
              m <- "log + threshold"
            }else if(mth =="log"){
              m <- "log"
            }
            if(xlog){
              plot(r,n,pch=".",main=main,log="x",
                   xlab=paste0(ab," (log)"),ylab=paste0(ab," (",m,")"),...)
            }else{
              plot(r,n,pch=".",main=main,
                   xlab=paste0(ab," (raw)"),ylab=paste0(ab," (",m,")"),...)
            }
            abline(v=th[ab],col=2)
            abline(h=p_thres,lty=2,col=1)
          }
)
