setGeneric("slidePlot", function(x,...) standardGeneric("slidePlot"))
setMethod("slidePlot", "Cycif",
	function(x,pch=".",cex=1,type=c("exp","cell_type","custom"),ttl,ab,
	         uniq.cols,na.col="grey90",cell_type,roi.selected,ncells=1e4,
	         legend=FALSE, legend.pos="bottomright",mar=c(3,3,3,3),...){
	  n <- exprs(x,type="normalized")
	  smpl <- names(x)
	  if(missing(type)){
	    stop("need to specify the color_code by 'type' argument")
	  }
	  if(type=="exp"){
	    if(missing(ab) || !ab %in% abs_list(x)$ab){
  	    stop(ab, " is not specified or available in the sample ", names(x))
  	  }else{
    	  n.ab <- trim_fun(n[[ab]],trim_th=1e-3)
    	  rn <- range(n.ab,na.rm=T)
    	  # stop(sum(is.na(rn)))
        is.na <- is.na(n.ab)
        if(!missing(roi.selected)){
          is.na <- is.na | !roi.selected
        }
  	  }
	    if(missing(uniq.cols)){
	      nlev <- 50

	      idx <- round((n.ab[!is.na] - rn[1])/diff(rn)*(nlev-1))+1
        uniq.cols <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(nlev)
      	cols <- rep(na.col,length(n.ab))
      	cols[!is.na] <- uniq.cols[idx]
	    }else{
        idx <- n.ab[!is.na]
        cols <- rep(na.col,length(n.ab))
        cols[!is.na] <- uniq.cols[idx]
	    }

  	 if(missing(ttl)){
  	   ttl <- paste0(smpl,", ",ab)
  	 }
	  }else if(type=="cell_type"){
	    if(missing(cell_type)){
	      stop("when type='cell_type', argument 'cell_type' can't be missing.")
	    }
	    if(missing(uniq.cols)){
	        nct <- nlevels(cell_type)
	        if(nct>11){
  	        uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nct)
	        }else{
  	        uniq.cols <- RColorBrewer::brewer.pal(nct,"Spectral")
	        }
	        names(uniq.cols) <- levels(cell_type)
      }
      cols <- uniq.cols[cell_type]
      is.na <- is.na(cell_type)
      if(!missing(roi.selected)){
        is.na <- is.na | !roi.selected
      }

	    if(missing(ttl)){
	      ttl <- paste0(smpl,", cell-types")
	    }
	  }

    xy <- xys(x)
  	xy$Y_centroid <- max(xy$Y) - xy$Y

  	if(!is.na(ncells) && nrow(xy) > ncells){
  	  set.seed(123)
  	  is.used <- seq(nrow(xy)) %in% sort(sample(nrow(xy),ncells))
  	}else{
  	  is.used <- rep(TRUE,nrow(xy))
  	}
  	par(mar=mar)
  	plot(xy$X,xy$Y,main=ttl,asp=1,xlab="",ylab="",type="n",...)
  	points(xy$X[is.na & is.used],xy$Y[is.na & is.used],col=na.col,pch=pch,cex=cex*.8)
  	if(type=="cell_type"){
  	  points(xy$X[!is.na & cell_type=="_" & is.used],
  	         xy$Y[!is.na & cell_type=="_" & is.used],
  	         col=cols[!is.na & cell_type=="_" & is.used],
  	         pch=pch,cex=cex)
  	  points(xy$X[!is.na & cell_type!="_" & is.used],
  	         xy$Y[!is.na & cell_type!="_" & is.used],
  	         col=cols[!is.na & cell_type!="_" & is.used],
  	         pch=pch,cex=cex)
  	}else{
  	  points(xy$X[!is.na & is.used],
  	         xy$Y[!is.na & is.used],
  	         col=cols[!is.na & is.used],
  	         pch=pch,cex=cex)
  	}

  	if(legend){
  	  legend(legend.pos,names(uniq.cols),col=uniq.cols,pch=20)
  	}
  }
)

trim_fun <- function(x,trim_th = 1e-3){
  qts <- quantile(x,c(trim_th,1-trim_th),na.rm=T)
  x[x < qts[1]] <- qts[1]
  x[x > qts[2]] <- qts[2]
  return(x)
}


#' @export
setGeneric("plotAvailCellOnSlide",function(x,...) standardGeneric("plotAvailCellOnSlide"))
setMethod("plotAvailCellOnSlide", "Cycif",
	function(x,upside.down=TRUE,ncycle,mfrow=c(3,3),mar=c(0,0,4,0),legend=TRUE,main=names(x),cex.title=1,
		uniq.cols=c(lost="grey90",dropped="blue",available="black",bunched="red"),legend.cex=2,...){
	stopifnot(nrow(x@used_cells)>0)

	u <- x@used_cells # not to be replaced with cumUsedCells
	nchannels <- ncol(u)

	nc.ratio <- round(statUsedCells(x)*100,1)

	xy <- xys(x)
	if(upside.down){
		xy$Y_centroid <- max(xy$Y) - xy$Y
	}

	par(oma=c(2,2,8,2))
	par(mfrow=mfrow)
	v <- sapply(seq(nchannels),function(i){
		if(i==1){
			id <- rep(2,nrow(u))
		}else{
			id <- u[,i] + 1
			is.avail <- rowSums(u[,seq(i-1),drop=F]==1)==i-1
			id[!is.avail] <- 0
		}

		rs <- id + 1

		main.txt <- paste0("cycle ",i-1," (",nc.ratio[i],"% available)")
		par(mar=mar)
		plot(xy,pch='.',col=uniq.cols[rs],main=main.txt,axes=F,asp=1,
			xlab="",ylab="",...)
		box()
	})
	if(legend){
		plot.new()
		par(mar=c(1,1,1,1))
		# box()
		legend("topleft",c("available","dropped","bunched","lost previously"),
			fill=uniq.cols[c("available","dropped","bunched","lost")],
			cex=legend.cex)
	}
	if(!missing(main)){
  	mtext(main, side=3, line=2, cex=cex.title, outer=TRUE)
	}
})
