#' Violin plots to show protein expressions
#' @export
setGeneric("heatmapSummAb", function(x,...) standardGeneric("heatmapSummAb"))

#' @export
setMethod("heatmapSummAb", "CycifStack",
  function(x,norm_type=c("log_normalized","logTh_normalized"), ab,
           sum_type=c("freq","mean","median","x percentile"),
           ctype.full=FALSE,rois = TRUE,
           uniq.cts,uniq.smpls,strict=FALSE,scale="row",...){
    options(dplyr.summarise.inform = FALSE)

    if(missing(sum_type)){
      stop("'sum_type' should be specified")
    }
    if(missing(norm_type)){
      norm_type <- "log_normalized"
    }
    if(missing(ctype.full)){
      ctype.full <- FALSE
    }
    if(missing(ab)){
      stop("ab should be always specified")
    }

    if(norm_type=="log_normalized"){
      balanceColor <- FALSE
    }else if(norm_type=="logTh_normalized"){
      balanceColor <- TRUE
      scale <- "none"
    }
    ## cell types
    cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=rois,strict=strict)

    if(grepl("percentile$",sum_type)){
      pct <- strsplit(sum_type," ")[[1]]
      p.ile <- as.numeric(pct[[1]])
      sum_fun=function(y,...)quantile(y,p.ile/100,na.rm=T)
    }else if(sum_type=="median"){
      sum_fun=function(y,...)quantile(y,.5,na.rm=T)
    }else if(sum_type=="mean"){
      sum_fun=function(y,...)mean(y,na.rm=T)
    }else if(sum_type=="freq"){
      thres <- unlist(x@cell_type@gates[ab,]) # identical whichever ctype.full
      sum_fun=function(y,th)sum(y > th,na.rm=T)/sum(!is.na(y))
    }

    df <- exprs(x,type=norm_type) %>%
      tibble::rownames_to_column(var="smpl") %>%
      mutate(smpl = sub("\\.[0-9]+$","",smpl)) %>%
      mutate(celltype=factor(cts))


    if(missing(uniq.smpls)){
      uniq.smpls <- unique(df$smpl)
    }
    if(missing(uniq.cts)){
      uniq.cts <- levels(df$celltype)
    }

    pd <- pData(x) %>% rename(smpl=id)
    df <- df %>%
      filter(within_rois(x)) %>%
      filter(!is.na(celltype)) %>%
      filter(smpl %in% uniq.smpls) %>%
      filter(celltype %in% uniq.cts) %>%
      rename(this_ab=as.symbol(ab)) %>%
      mutate(thres1 = thres[smpl]) %>%
      group_by(smpl,celltype) %>%
      summarise(sum_ab=sum_fun(this_ab,thres1)) %>%
      left_join(pd,by="smpl") %>%
      mutate(smpl = factor(smpl))

    ttl <- paste0(ab,",",sum_type,",",sub("_.+","",norm_type),",(scale:",scale,")")
    m1 <- as.matrix(with(df, sparseMatrix(as.integer(smpl), as.integer(celltype), x=sum_ab)))
    m1[m1==0] <- NA
    rownames(m1) <- levels(df$smpl)
    colnames(m1) <- levels(df$celltype)
    m1 <- m1[rev(seq(nrow(m1))),]

    if(norm_type=="logTh_normalized"){
      m1 <- m1 - 0.5 ## hardcode!!!!
    }

    h3(m1,
       margins=c(10,5),
       main=ttl,
       cexRow=.7,
       na.rm = T,
       balanceColor=balanceColor,
       scale=scale,
       ...)
})

