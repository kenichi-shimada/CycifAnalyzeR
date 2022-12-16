#' Split a string
#'
#' @param x An object of Cycif class.
#' @param show.cycles.in.row logical. FALSE by default. If TRUE, the output plot shows the number
#'     of cycles each sample was stained for.
#' @param ... Arguments passed to image() function.
#'
#' @export
#'
#' @include CycifStack-class.R
setGeneric("AbsSummary", function(x,...) standardGeneric("AbsSummary"))

#' @rdname AbsSummary
#' @export
setMethod("AbsSummary", "CycifStack", function(x,show.cycles.in.row=FALSE,...){
  uniq.abs <- x@abs_list$ab
  n1 <- do.call(rbind,lapply(x@samples,function(y){
    this.abs <- abs_list(y)$ab
    tested <- uniq.abs %in% this.abs
    return(tested)
  }))

  idx <- rep(seq(ncol(n1)/3),each=3)
  n2 <- n1 * idx[col(n1)]
  n2[n2==0] <- NA
  pcols <- grey(c(.8,.6))[(seq(max(n2,na.rm=T)) %% 2) + 1]

  if(show.cycles.in.row){
    smpl.labs <- paste0(x@names,"\n(",x@n_cycles, " cycles)")
  }else{
    smpl.labs <- x@names
  }
  col.labs <- paste0("Cycle\n",seq(x@max_cycles))
  image(seq(ncol(n2)),seq(nrow(n2)),t(n2[rev(seq(nrow(n2))),]),col=pcols,
        xlab="",ylab="",axes=F,...)
  box()
  abline(h=c(0,seq(nrow(n1)))+.5,lwd=.5)
  abline(v=c(0,seq(ncol(n1)))+.5,lwd=.5)
  axis(2,at=rev(seq(nrow(n1))),labels=smpl.labs,las=1)
  axis(1,at=seq(ncol(n2)),labels=uniq.abs,las=2)
  par(tcl=0)
  axis(3,at=seq(x@max_cycles)*3-1,labels=col.labs)
  par(tcl=-.5)

  return(invisible(NULL))
})
