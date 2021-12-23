#' @export
setGeneric("defCellType", function(x,...) standardGeneric("defCellType"))
setMethod("defCellType", "Cycif",


          # this is not done yet.

	function(x,ctype.def,...){
		## default method is logTh
		if(missing(method)){
			method <- "logTh"
		}

		## 'used abs' set already?
		if("used_abs" %in% slotNames(x)){
			used.abs <- used_abs(x)
		}else{
			used.abs <- as.character(abs_list(x)$ab)
		}

		smpl <- names(x)
		raw <- x@raw
		is.used <- cumUsedCells(x)
		x@normalize.method <- method

		## treatment is different between methods
		if(method=="log"){
			norm <- as.data.frame(
				sapply(used.abs,function(ab){
					cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
					cycle <- cycle + 1 # 0-origin
					is.used.1 <- is.used[,cycle]

					r <- raw[[ab]]
					n <- rep(NA,length(r))

					n[is.used.1] <- transform(r[is.used.1],method="log",trim=trim)
					return(n)
				})
			)
		}else if(method=="logTh"){
			if(!.hasSlot(x,"threshold")){
				stop(smpl,": set threshold first.\n")
			}

			thres <- threshold(x)
			used.abs <- used.abs[!is.na(thres[used.abs])]

			norm <- as.data.frame(
				sapply(used.abs,function(ab){
					cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
					cycle <- cycle + 1 # 0-origin
					is.used.1 <- is.used[,cycle]
					r <- raw[[ab]]
					n <- rep(NA,length(r))
					th <- thres[ab]
					n[is.used.1] <- transform(r[is.used.1],method="logTh",th=th,trim=trim)
					return(n)
				})
			)
		}
		x@normalized <- norm
		return(x)
	}
)
