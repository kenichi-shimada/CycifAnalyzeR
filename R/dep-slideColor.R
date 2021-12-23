#' @include Cycif-class.R

setGeneric("scaled_color", function(x,...) standardGeneric("scaled_color"))
setMethod("scaled_color", "Cycif",
	function(x,ch="DNA0",subset=NULL,palette=FALSE,trim=0.01,int.range=NULL,...){
	n.levels <- 50 #
	col.palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"Spectral")))(n.levels+1)
	if(palette){
		return(col.palette)
	}

	if(grepl("^DNA[0-9]+$",ch)){
		int <- x@dna[[ch]]
	}else if(ch %in% x@abs_list$ab){
		int <- x@raw[[ch]]
	}else{
		stop("ab should be either an exact DNA or Ab name registered (see abs_list(x))")
	}

	if(!is.null(trim) && trim >= 0 && trim <= 1){
		qt1 <- quantile(int,c(trim,1-trim))
		int[int < qt1[1]] <- qt1[1]
		int[int > qt1[2]] <- qt1[2]
	}

	if(!is.null(int.range)){
		int[int < int.range[1]] <- int.range[1]
		int[int > int.range[2]] <- int.range[2]
	}

	rng <- range(int,int.range)

	if(!is.null(subset) && is.logical(subset) && length(subset) %in% c(1,length(int))){
		int <- int[subset]
	}

	nn <- round((int-rng[1])/diff(rng)*n.levels)
	cols <- col.palette[nn+1]

	return(cols)
})
