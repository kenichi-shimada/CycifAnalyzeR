#' Violin plots to show protein expressions
#' @export
setGeneric("heatmapExp", function(x,...) standardGeneric("heatmapExp"))

#' @export
setMethod("heatmapExp", "CycifStack",
  function(x,type=c("log_normalized","logTh_normalized"),
           ab="PDL1", sum_fun=function(y)quantile(y,.5,na.rm=T),
           ctype.full=FALSE,rois = TRUE,
           uniq.cts,uniq.smpls,strict=FALSE,...){
    options(dplyr.summarise.inform = FALSE)
    if(missing(type)){
      type <- "log_normalized"
    }
    if(missing(ctype.full)){
      ctype.full <- FALSE
    }
    if(missing(ab)){
      stop("ab should be always specified")
    }
    if(missing(uniq.smpls)){
      uniq.smpls <- levels(df$smpl)
    }
    if(missing(uniq.cts)){
      uniq.cts <- levels(df$celltype)
    }

    ## cell types
    cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=rois,strict=strict)

    df <- exprs(x,type=type) %>%
      tibble::rownames_to_column(var="smpl") %>%
      mutate(smpl = sub("\\.[0-9]+$","",smpl)) %>%
      mutate(celltype=factor(cts)) %>%
      filter(within_rois(x)) %>%
      filter(!is.na(celltype)) %>%
      filter(smpl %in% uniq.smpls) %>%
      filter(celltype %in% uniq.cts) %>%
      rename(this_ab=as.symbol(ab)) %>%
      group_by(smpl,celltype) %>%
      summarise(sum_ab=sum_fun(this_ab)) %>%
      left_join(pd,by="smpl") %>%
      mutate(smpl = factor(smpl))

    ttl <- paste0(ab," (",type,")")
    m1 <- as.matrix(with(df, sparseMatrix(as.integer(smpl), as.integer(celltype), x=sum_ab)))
    m1[m1==0] <- NA
    rownames(m1) <- levels(df$smpl)
    colnames(m1) <- levels(df$celltype)
    m1 <- m1[rev(seq(nrow(m1))),]

    h3(m1,margins=c(10,5),main=ttl,cexRow=.7,na.rm = T,...)
})

