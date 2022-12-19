#' Violin plots to show protein expressions
#' @export
setGeneric("vlnPlot", function(x,...) standardGeneric("vlnPlot"))

#' @export
setMethod("vlnPlot", "CycifStack",
  function(x,type=c("raw","log_normalized"), xaxis=c("celltype","smpl"),
           ab="PDL1", ctype=c("short","long"),uniq.cts,uniq.smpls){
    require(tidyr)
    require(ggplot2)

    if(missing(type)){
      type <- "log_normalized"
    }

    if(missing(ctype)){
      ctype <- "short"
    }

    if(ctype=="short"){
      cts <- cell_types(x,full=FALSE,leaves.only=TRUE)
    }else if(ctype=="long"){
      cts <- cell_types(x,full=TRUE,leaves.only=TRUE)
    }

    df <- exprs(x,type=type) %>%  mutate(celltype=factor(cts))
    df <- df %>%
      mutate(smpl = factor(sub("\\.[0-9]+$","",rownames(df))))
      # select(c(one_of(ab),celltype,smpl))

    ab_thres <- x@cell_type@gates[ab,]

    if(missing(uniq.smpls)){
      uniq.smpls <- levels(df$smpl)
    }
    if(missing(uniq.cts)){
      uniq.cts <- levels(df$celltype)
    }

    if(xaxis=="smpl"){
      ab <- "gH2AX"
      uniq.smpls <- levels(df$smpl)
      df1 <- df %>%
        filter(smpl %in% uniq.smpls) %>%
        filter(celltype %in% "Tumor") %>%
        mutate(smpl=factor(smpl))
      ttl <- paste0(ab,",Tumor")
      p <- ggplot(df1,aes_string("smpl",ab)) +
        geom_violin(aes(fill = smpl)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(ttl)
      print(p)

      df_thres <- data.frame(smpl=factor(names(ab_thres)),th=unlist(ab_thres))
      for(smpl in df_thres$smpl){
        i <- match(smpl,df_thres$smpl)
        amt <- 0.2
        th <- ab_thres[[smpl]]
        p <- p +  geom_segment(aes(x=i-amt/2,xend=i+amt/2,y=th,yend=th))
      }

      # stop(names(df_thres))
      print(p)
    }else if(xaxis=="celltype"){ # not so useful?
      ab <- "pTBK1"
      ab <- "gH2AX"
      ab <- "pERK"
      ab <- "pAKT"
      smpl <- uniq.smpls <- "13693"
      df1 <- df %>%
        mutate(celltype=factor(celltype,levels=uniq.cts)) %>%
        filter(smpl %in% uniq.smpls) %>%
        filter(!is.na(celltype))

      ttl <- smpl
      p <- ggplot(df1,aes_string("celltype",ab)) +
        # geom_hline(yintercept = ab_thres[[smpl]]) +
        geom_violin(aes(fill = celltype)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(ttl) +
        scale_fill_discrete(drop=FALSE) +
        scale_x_discrete(drop=FALSE)
      print(p)
    }
})
