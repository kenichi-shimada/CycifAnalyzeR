#' Violin plots to show protein expressions
#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @export
setMethod("vlnPlot", "CycifStack",
  function(x,strat.by=c("celltype","smpl"), ab="PDL1", ctype.full=FALSE,
           # type=c("raw","log_normalized","logTh_normalized"),
           strict=FALSE,
           within.rois = TRUE,uniq.cts,uniq.smpls,is.used){
    if(1){ # missing(type)){
      type <- "logTh_normalized"
    }

    if(missing(ab)){
      stop("ab should be always specified")
    }

    ## gates
    ab_thres <- unlist(x@cell_type@gates[ab,])
    # if(type=="raw"){
    #   ths <- expm1(ab_thres)
    # }else if(type=="log_normalized"){
    #   ths <- ab_thres
    # }else if(type=="logTh_normalized"){
    #   ths <- rep(0.5,length(ab_thres))
    # }
    ths <- rep(0.5,length(ab_thres))

    ## cell types
    cts <- cell_types(x,ctype.full=ctype.full,leaves.only=TRUE,within.rois=within.rois)
    is.strict <- unlist(cyApply(x,function(cy)cy@cell_type@is_strict))
    # table(cts) ## all zeros!!!

    if(strict){
      used <- within_rois(x) & is.strict
    }else{
      used <- within_rois(x)
    }
    if(!missing(is.used)){
      used <- used & is.used
    }

    df <- exprs(x,type=type) %>%
      tibble::rownames_to_column(var="smpl") %>%
      mutate(smpl = sub("\\.[0-9]+$","",smpl)) %>%
      mutate(celltype=factor(cts)) %>%
      filter(used) %>%
      left_join(pData(x) %>% rename(smpl=id),by="smpl") %>%
      mutate(smpl = factor(smpl))

    if(strat.by=="smpl"){
      if(missing(uniq.cts) || length(uniq.cts) != 1){
        stop("one celltype or 'all' should be set in 'uniq.cts'")
      }else if(uniq.cts == "all"){
        ttl <- paste0(ab,", all cells")
        uniq.cts <- levels(df$celltype)
      }else{
        if(length(uniq.cts)>1){
          uc <- paste0(length(uniq.cts)," samples")
        }else if(length(uniq.cts)>1){
          uc <- uniq.cts
        }

        ttl <- paste0(ab,", ",uc)
      }
      if(missing(uniq.smpls)){
        uniq.smpls <- levels(df$smpl)
      }
      df1 <- df %>%
        filter(smpl %in% uniq.smpls) %>%
        filter(celltype %in% uniq.cts) %>%
        filter(TimePoint %in% c("BS","BX2","BX3")) %>%
        mutate(TimePoint = factor(TimePoint)) %>%
        mutate(thres=ab_thres[as.character(smpl)])
      p <- ggplot(df1,aes_string(x="Patient.ID",y=ab,fill="TimePoint")) +
        geom_violin(position = position_dodge(width = 0.9)) +
        ylab(paste0(ab," (",sub("_"," ",type),")")) +
        ggtitle(ttl) +
        geom_point(data=unique(df1 %>% select(Patient.ID,thres,TimePoint)),
                   aes(x=Patient.ID,y=thres,color=as.factor(TimePoint)),
                   position = position_dodge(width = 0.9),
                   shape = 95, size=10, show.legend=F) +
        scale_color_manual(values=rep("black",3),drop=F) +
        theme_bw() +
        theme(legend.position = "right") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      print(p)
    }else if(strat.by=="celltype"){ # not so useful?
      # if(missing(uniq.smpls) || length(uniq.smpls)!= 1){
      #   stop("one sample should be set in 'uniq.smpls'")
      # }
      if(missing(uniq.smpls)){
        uniq.smpls <- levels(df$smpl)
      }
      if(missing(uniq.cts)){
        uniq.cts <- levels(df$celltype)
        uniq.cts <- c(uniq.cts[uniq.cts != "Immune_other"],"Immune_other")
        uniq.cts <- c(uniq.cts[uniq.cts != "unknown"],"unknown")
      }
      df1 <- df %>%
        mutate(celltype=factor(celltype,levels=uniq.cts)) %>%
        filter(smpl %in% uniq.smpls) %>%
        filter(!is.na(celltype))

      if(length(uniq.smpls)==1){
        us <- uniq.smpls
      }else if(length(uniq.smpls)>1){
        us <- paste0(length(uniq.smpls), " samples")
      }
      ttl <- paste0(ab,", ",us)
      p <- ggplot(df1,aes_string("celltype",ab)) +
        geom_violin(aes(fill = celltype)) +
        ggtitle(ttl) +
        geom_hline(yintercept = 0.5,col="black") + #ab_thres[[uniq.smpls]],col="black") +
        theme_bw() +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      print(p)
    }
})

