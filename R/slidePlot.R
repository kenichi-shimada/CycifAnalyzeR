#' @export
setGeneric("slidePlot", function(x,...) standardGeneric("slidePlot"))
setMethod("slidePlot", "Cycif",
	function(x,pch=20,cex=2,plot_type=c("dna","exp","cell_type","filter"),
	         within_filter_rng,ctype.full=FALSE,strict=FALSE,ttl,ab,
	         uniq.cts,uniq.cols,draw.roi=TRUE,roi.cycle,
	         remove.unknown=TRUE,cell.order,
	         na.col="grey80",use.roi=TRUE,use.thres=TRUE,ncells=1e4,
	         contour=FALSE,cont_nlevs=3,
	         trim_th=1e-2,legend=FALSE, legend.pos="bottomright",mar=c(3,3,3,3),...){
	  if(missing(plot_type)){
	    stop("need to specify 'plot_type' argument: dna, exp, cell_type, filter")
	  }

	  smpl <- names(x)

	  if(!missing(cell.order)){
	    if(length(cell.order) != nCells(x)){
	      stop("'cell.order' should be the same size as nrow(x@dna)")
	    }
	  }

	  if(plot_type=="dna"){
	    if(missing(ab) || !ab %in% paste0("DNA",seq(nCycles(x)))){
  	    stop(paste0("If DNA stain, `ab' should be one of DNA1, ..., DNA", nCycles(x),")"))
  	  }
	    mat <- x@dna

	    mat <- cbind(log1p(mat[[1]]),as.data.frame(lapply(mat,function(x)log1p(x/mat[[1]])))[-1])
	    names(mat) <- names(x@dna)

  	  n.ab <- trim_fun(mat[[ab]],trim_th=trim_th)
  	  rn <- range(n.ab,na.rm=T)

      is.na <- is.na(n.ab)
      if(use.roi){
        within.rois <- x@within_rois
        if(length(within.rois)){
          stop("ROIs not defined yet; 'use.roi' should be FALSE")
        }
        is.na <- is.na | !within.rois
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
  	   if(ab=="DNA1"){
  	     ab1 <- ab
  	   }else{
  	     ab1 <- paste0(ab,"/DNA1")
  	   }
  	   ttl <- paste0(smpl," (",ab1,")")
  	 }
	  }else if(plot_type=="exp"){
	    cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=use.roi,strict=strict)
	    if(missing(uniq.cts)){
	      uniq.cts <- levels(cts)
	      if(any(uniq.cts=="unknown")){
	        if(remove.unknown){
	          uniq.cts <- uniq.cts[uniq.cts!="unknown"]
	        }else{
	          ui <- which(uniq.cts=="unknown")
	          uniq.cts <- c(uniq.cts[-ui],"unknown")
	        }
	      }
	    }
	    cts <- factor(cts,levels=uniq.cts)

      if(missing(ab) || !ab %in% abs_list(x)$ab){
        stop(ab, " is not specified or available in the sample ", names(x))
      }

	    n <- exprs(x,plot_type="log_normalized")

	    n.ab <- trim_fun(n[[ab]],trim_th=trim_th)
      rn <- range(n.ab,na.rm=T)
      is.na <- is.na(n.ab)
      if(use.roi){
        within.rois <- x@within_rois
        is.na <- is.na | !within.rois
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
        ttl <- paste0(smpl,", ",ab," expression")
      }
    }else if(plot_type=="cell_type"){
      cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=use.roi,strict=strict)

      if(missing(ttl)){
        if(missing(uniq.cts)){
          ttl <- paste0(smpl,", cell-types")
        }else{
          ttl <- paste0(smpl,", ",paste0(uniq.cts,collapse=","))
        }
      }
      if(missing(uniq.cts)){
        uniq.cts <- levels(cts)

        if(any(uniq.cts=="Immune_other")){
          ui <- which(uniq.cts=="Immune_other")
          uniq.cts <- c(uniq.cts[-ui],"Immune_other")
        }

        if(any(uniq.cts=="unknown")){
          if(remove.unknown){
            uniq.cts <- uniq.cts[uniq.cts!="unknown"]
          }else{
            ui <- which(uniq.cts=="unknown")
            uniq.cts <- c(uniq.cts[-ui],"unknown")
          }
        }
      }
      cts <- factor(cts,levels=uniq.cts)

	    if(missing(uniq.cols)){
	        nct <- nlevels(cts)
	        if(nct>10){
  	        uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[-6])(nct)
	        }else{
	          # set.seed(12)
  	        uniq.cols <- (RColorBrewer::brewer.pal(11,"Spectral"))[c(7,1,2,4,5,11:9,8)]
	        }
	        names(uniq.cols) <- uniq.cts
      }
      cols <- uniq.cols[cts]
      is.na <- is.na(cts)
      if(use.roi){
        within.rois <- x@within_rois
        is.na <- is.na | !within.rois
      }


    }else if(plot_type=="filter"){
      if(missing(within_filter_rng)){
        stop("when plot_type='filter', the argument 'within_filter_rng' should be specified.")
      }
      if(class(within_filter_rng)!="factor"){
        within_filter_rng <- factor(within_filter_rng)
      }
      if(missing(uniq.cols)){
        nct <- nlevels(within_filter_rng)
        if(nct>11){
          uniq.cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nct)
        }else{
          uniq.cols <- RColorBrewer::brewer.pal(nct,"Spectral")
        }
        names(uniq.cols) <- levels(within_filter_rng)
      }
      cols <- uniq.cols[within_filter_rng]
      is.na <- is.na(within_filter_rng)
      if(use.roi){
        within.rois <- x@within_rois
        is.na <- is.na | !within.rois
      }

      if(missing(cell.order)){
        cell.order <- order(within_filter_rng,decreasing=T)
      }
      if(length(pch)==1){
        pch <- rep(pch,length(within_filter_rng))
      }
      if(length(cex)==1){
        cex <- rep(cex,length(within_filter_rng))
      }
      ttl <- ""
    }

 	  ## plot
    xy <- xys(x)

    xy$Y_centroid <- max(xy$Y) - xy$Y
  	prs <- x@rois

  	if(!is.na(ncells) && nrow(xy) > ncells){
  	  set.seed(123)
  	  is.used <- seq(nrow(xy)) %in% sort(sample(nrow(xy),ncells))
  	}else{
  	  is.used <- rep(TRUE,nrow(xy))
  	}

  	omar <- par()$mar
  	##
  	if(!missing(cell.order)){
  	  xy <- xy[cell.order,]
  	  if(exists("is.na") && length(is.na)==length(cell.order)){
  	    # stop("is.na")
  	    is.na <- is.na[cell.order]
  	  }
  	  if(exists("is.used") && length(is.used)==length(cell.order)){
  	    is.used <- is.used[cell.order]
  	  }
  	  if(exists("pch") && length(pch)==length(cell.order)){
  	    pch <- pch[cell.order]
  	  }
  	  if(exists("cts") && length(cts)==length(cell.order)){
  	    cts <- cts[cell.order]
  	  }
  	  if(exists("cols") && length(cols)==length(cell.order)){
  	    cols <- cols[cell.order]
  	  }
  	}

  	par(mar=mar)
  	plot(xy$X,xy$Y,main=ttl,asp=1,xlab="",ylab="",type="n")#,...)

  	cex1 = 4.8 * 8/3 * (par()$pin[1]/par()$cin[2])/diff(par()$usr[1:2])

  	points(xy$X[is.na & is.used],xy$Y[is.na & is.used],col=na.col,pch=pch,cex=cex1)

  	if(plot_type=="cell_type"){
  	  points(xy$X[!is.na & cts=="unknown" & is.used],
  	         xy$Y[!is.na & cts=="unknown" & is.used],
  	         col=cols[!is.na & cts=="unknown" & is.used],
  	         pch=pch,cex=cex1*cex)
  	  points(xy$X[!is.na & cts!="unknown" & is.used],
  	         xy$Y[!is.na & cts!="unknown" & is.used],
  	         col=cols[!is.na & cts!="unknown" & is.used],
  	         pch=pch,cex=cex1*cex)
  	}else if(plot_type=="filter"){
  	  points(xy$X,
  	         xy$Y,
  	         col=cols,
  	         pch=pch,cex=cex1*cex)
  	}else{
  	  points(xy$X[!is.na & is.used],
  	         xy$Y[!is.na & is.used],
  	         col=cols[!is.na & is.used],
  	         pch=pch,cex=cex1*cex)
  	}
    if(draw.roi){
      ncycles <- sapply(prs,function(x)x$cycle)
      for(i in seq(prs)){
        pr <- prs[[i]]
        if(pr$cycle != roi.cycle){
          next;
        }
        # points(pr,pch=sub(".+(.)$","\\1",as.character(seq(length(pr$x)))))
        if(pr$dir=="positive"){
          col1 <- 2
        }else if(pr$dir=="negative"){
          col1 <- 4
        }
        polygon(pr$coords,lty=1,border=col1)
      }
    }
  	if(plot_type=="cell_type" && contour){
  	  this.idx <- !is.na & cts %in% uniq.cts & is.used
  	  f1 <- kde2d(xy$X[this.idx],xy$Y[this.idx], n = 1000, lims = par()$usr,
  	              h = c(width.SJ(xy$X), width.SJ(xy$Y)))
  	  contour(f1, nlevels  = cont_nlevs, add=T, labels=rep("",3),col="grey50")
  	}
  	if(legend){
  	  legend(legend.pos,names(uniq.cols),col=uniq.cols,pch=20)
  	}
  	par(mar=omar)
  }
)

#' @export
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
