#' Filtering cells based on the area of segmented cells
#'
#' @param x An object of Cycif class.
#' @param factor numeric. The default value is 10. A threshold for acceptable 'area' is defined by
#'          threshold_area = median_area +  MAD_area * factor
#'     where MAD_area is a median-absolute deviation. Area larger than threshold_area was rejected.
#' @param xlim xlim passed to histogram.
#' @param main argument passed to plot().
#' @param ... arguments passed to hist() function.
#'
#'
#' @export
setGeneric("areaFilter", function(x,...) standardGeneric("areaFilter"))
setMethod("areaFilter", "Cycif",
	function(x,factor=10,xlim=c(0,3000),...){
		if(.hasSlot(x,"segment_property")){
			a <- x@segment_property$Area
		}else if(.hasSlot(x,"nuclear_property")){
			a <- x@nuclear_property$Area
		}

		ha <- hist(a,breaks=1e4,xlim=c(0,3000),main="histogram of 'segmented cell area'",
			xlab="segmented cell area",...)
		mod.i <- which.max(ha$count)
		mod <- ha$mids[mod.i]
		abline(v=mod,col=2,lwd=2)

		counts <- ha$counts[-seq(mod.i)]
		mids <- ha$mids[-seq(mod.i)]
		med <- median(rep(mids-mod,counts))
		med2 <- mod + med
		medf <- med * factor + mod
		abline(v=med2,col=3)
		abline(v=medf,col=4,lwd=2)

		is.used <- a < mad

		x@area_filter <- is.used

		return(x)
})
