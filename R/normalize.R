#' Show a raw or normalized protein expression matrix
#'
#' @param r A numeric scholar indicating a raw expression value.
#' @param x A Cycif or CycifStack object.
#' @param method Either "log" or "logTh". "log" indicates that
#'   the raw values are transformed using log1p() function. "logTh" method
#'   further trim the outliers outside outside [trm,1-trm] quantiles,
#'   where trm is specified by trim parameter. Outlier-trimmed, log1p-transformed
#'   raw values are further transformed such that the value indicated by the threshold,
#'   th (in the case of transform() function) or threshold(x) (in the case of
#'   normalize() function) will be set to p_thres. Values below or above the threshold
#'   are linearly transformed to values between 0 and p_thres, or values between p_thres and 1,
#'   respectively.
#' @param th A numeric scholar. A user-provided threshold for a raw expression value
#'   for each protein for each sample.
#' @export
#' @rdname normalize
setGeneric("normalize", function(x,...) standardGeneric("normalize"))

#' @rdname normalize
#' @export
setMethod("normalize", "Cycif",
  function(x,method=c("log","logTh","exp"),trim=1e-3,p_thres=0.5){
    ## default method is logTh
    if(missing(method)){
      method <- "logTh"
    }

    used.abs <- as.character(abs_list(x)$ab)

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
      x@log_normalized <- norm
    }else if(method=="logTh"){
      if(!.hasSlot(x,"threshold")){
        stop(smpl,": set threshold first.\n")
      }

      thres <- threshold(x)
      used.abs <- names(thres) # used.abs[!is.na(thres[used.abs])]

      # warning("cycle should be switched to 1-origin")

      ## change here somehow to convert values of unused abs to NA below:

      norm <- as.data.frame(
        sapply(used.abs,function(ab){
          n <- rep(NA,nrow(raw))
          if(ab %in% used_abs(x)){
            cycle <- abs_list(x)$cycle[abs_list(x)$ab==ab]
            cycle <- cycle + 1 # 0-origin - confusing
            is.used.1 <- is.used[,cycle]
            r <- raw[[ab]]
            th <- thres[ab]
            n[is.used.1] <- transform(r[is.used.1],method="logTh",th=th,trim=trim,p_thres=p_thres)
          }
          return(n)
        })
      )
      x@log_normalized <- norm
    }else if(method=="exp"){
      norm <- x@log_normalized
      raw <- expm1(norm)
      x@raw <- raw
    }

    return(x)
  }
)

#' @rdname normalize
#' @export
setMethod("normalize", "CycifStack",
  function(x,method=c("log","logTh","exp"),trim=1e-3,p_thres=0.5,nCycle){
    if(missing(method)){
      method <- "logTh"
    }
    xs <- cyApply(x,function(y){
      normalize(y,method=method,trim=trim,p_thres=p_thres)
    })

    if(method %in% c("log")){
      if(missing(nCycle)){
        nCycle <- max(cyApply(x,nCycles,simplify=TRUE)) - 1
      }
      norm <- do.call(rbind,cyApply(xs,function(cy)exprs(cy,type="log_normalized")[,seq((nCycle+1)*3)],as.CycifStack=FALSE))  %>%
        mutate(smpl=rep(names(xs),cyApply(xs,nCells,as.CycifStack=FALSE,simplify=TRUE)))
      xs@log_normalized <- norm
    }else if(method %in% c("logTh")){
      norm <- do.call(rbind,cyApply(xs,function(cy)exprs(cy,type="logTh_normalized"),as.CycifStack=FALSE))  %>%
        mutate(smpl=rep(names(xs),cyApply(xs,nCells,as.CycifStack=FALSE,simplify=TRUE)))
      xs@logTh_normalized <- norm
    }

    return(xs)
  }
)


## transform - at some point I should switch to nls() to apply sigmoidal curve
#' @rdname normalize
#' @export
transform <- function(r,method=c("log","logTh","Th","exp"),th,p_thres=0.5,trim=1e-3){
  if(missing(method)){
    method <- "logTh"
  }
  # r - raw intensity value (in quantification/*.csv)

  ## log-transformation (to be replaced with log with other bases)
  if(method=="Th"){
    r1 <- lr <- r
    lth <- th
  }else if(method=="logTh"){ ## when method contains log-transformation
    r1 <- lr <- log1p(r)
    lth <- log1p(th)
  }else if(method=="log"){
    r1 <- lr <- log1p(r)
  }else if(method=="exp"){
    r1 <- expm1(r)
  }

  ## method: either 'log' or 'logTh'
  if(method %in% c("logTh","Th")){
    stopifnot(!missing(th)) # th is an essential input from user

    ## trimming outliers
    qt <- quantile(lr,c(trim,1-trim))

    ## if lth is outside [trim, 1-trim], move lth to the bound
    if(lth < qt[1]){
      qt[1] <- lth
    }
    if(lth > qt[2]){
      qt[2] <- lth
    }

    if(any(lr > lth)){
      r1[lr > lth] <- (lr[lr > lth] - lth)/((qt[2] - lth)/(1-p_thres)) + p_thres
    }
    if(any(lr < lth)){
      r1[lr <= lth] <- (lr[lr <= lth] - qt[1])/((lth - qt[1])/p_thres)
    }
    r1[r1 < 0] <- 0
    r1[r1 > 1] <- 1
  }

  return(r1)
}

#' @rdname normalize
#' @export
trim_fun <- function(r,trim_th=1e-3){
  qt <- quantile(r,c(trim_th,1-trim_th),na.rm=T)

  r[r > qt[2]] <- qt[2]
  r[r < qt[1]] <- qt[1]

  return(r)
}

