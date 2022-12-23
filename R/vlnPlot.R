#' Violin plots to show protein expressions
#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @export
setMethod("vlnPlot", "CycifStack",
  function(x,xaxis=c("celltype","smpl"),
           type=c("raw","log_normalized","logTh_normalized"),
           ab="PDL1", ctype.full=FALSE,
           rois = TRUE,
           uniq.cts,uniq.smpls){
    if(missing(type)){
      type <- "log_normalized"
    }
    if(missing(ctype.full)){
      ctype.full <- FALSE
    }
    if(missing(ab)){
      stop("ab should be always specified")
    }

    if(xaxis=="celltype"){

    }else{

    }
    if(xaxis=="smpl"){

    }else{

    }

    ## gates
    ab_thres <- unlist(x@cell_type@gates[ab,])
    if(type=="raw"){
      ths <- expm1(ab_thres)
    }else if(type=="log_normalized"){
      ths <- ab_thres
    }else if(type=="logTh_normalized"){
      ths <- rep(0.5,length(ab_thres))
    }

    ## cell types
    cts <- cell_types(x,full=ctype.full,leaves.only=TRUE,within.rois=rois)

    df <- exprs(x,type=type) %>%
      tibble::rownames_to_column(var="smpl") %>%
      mutate(smpl = sub("\\.[0-9]+$","",smpl)) %>%
      filter(within_rois(x)) %>%
      mutate(celltype=factor(cts)) %>%
      left_join(pData(cs1) %>% rename(smpl=id),by="smpl") %>%
      mutate(smpl = factor(smpl))

    if(xaxis=="smpl"){
      if(missing(uniq.cts) || length(uniq.cts) != 1){
        stop("one celltype should be set in 'uniq.cts'")
      }
      if(missing(uniq.smpls)){
        uniq.smpls <- levels(df$smpl)
      }
      df1 <- df %>%
        filter(smpl %in% uniq.smpls) %>%
        filter(celltype %in% uniq.cts) %>%
        mutate(thres=ab_thres[as.character(smpl)])
      ttl <- paste0(ab,",",uniq.cts)
      p <- ggplot(df1,aes_string(x="Patient.ID",y=ab,fill="TimePoint")) +
        geom_violin(position = position_dodge(width = 0.9)) +
        ylab(paste0(ab," (",sub("_"," ",type),")")) +
        ggtitle(ttl) +
        geom_point(data=unique(df1 %>% select(Patient.ID,thres,TimePoint)),
                   aes(x=Patient.ID,y=thres,color=TimePoint),
                   position = position_dodge(width = 0.9),
                   shape = 95, size=10) +
        theme_bw() +
        theme(legend.position = "right") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      print(p)
    }else if(xaxis=="celltype"){ # not so useful?
      if(missing(uniq.smpls) || length(uniq.smpls)!= 1){
        stop("one sample should be set in 'uniq.smpls'")
      }
      if(missing(uniq.cts)){
        uniq.cts <- levels(df$celltype)
      }
      df1 <- df %>%
        mutate(celltype=factor(celltype,levels=uniq.cts)) %>%
        filter(smpl %in% uniq.smpls) %>%
        filter(!is.na(celltype))

      ttl <- paste0(ab,",",uniq.smpls)
      p <- ggplot(df1,aes_string("celltype",ab)) +
        geom_violin(aes(fill = celltype)) +
        ggtitle(ttl) +
        scale_fill_discrete(drop=FALSE) +
        scale_x_discrete(drop=FALSE) +
        geom_hline(yintercept = ab_thres[[smpl]],col="black") +
        theme_bw() +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      print(p)
    }
})

