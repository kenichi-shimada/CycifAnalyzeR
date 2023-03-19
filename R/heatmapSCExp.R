#' Violin plots to show protein expressions
#' @export
setGeneric("heatmapSCExp", function(x,...) standardGeneric("heatmapSCExp"))

#' @export
setMethod("heatmapSCExp", "CycifStack",
  function(x,ld_name,norm_type=c("logTh_normalized","log_normalized"),
           summarize=FALSE,summarize.by=c("BOR","TimePoint"),
           order="cBTgP",margins=c(7,5),used.abs,ttl,...){
    options(dplyr.summarise.inform = FALSE)
    if(missing(norm_type)){
      norm_type <- "logTh_normalized"
    }

    ##
    if(missing(ld_name)){
      stop("'ld_name' should be specified (it's used to retrieve the data later)")
    }else if(!ld_name %in% ld_names(x)){
      stop("Specified 'ld_name' does not exist.")
    }else{
      ld <- ld_coords(x,ld_name=ld_name)
    }

    if(summarize){
      meta <- createRSC(x,with.clusters=T,ld_name=ld_name)
    }else{
      meta <- createRSC(x,with.clusters=T,ld_name=ld_name,order=order)
    }

    pd <- meta$pd
    rsc <- meta$rsc

    xys <- ld@ld_coords
    ld.type <- ld@ld_type
    norm_type <- ld@norm_type

    is.used <- ld@is_used
    if(missing(used.abs)){
      used.abs <- ld@used.abs
    }
    used.cts <- ld@used.cts
    cts.params <- ld@cts_params
    cls <- ld@clusters

    n <- exprs(x,type=norm_type)[is.used,used.abs]

    if(summarize){
      if(!all(summarize.by %in% colnames(rsc))){
        na.cols <-summarize.by[!summarize.by %in% colnames(rsc)]
        stop("'columns in 'summarize.by' (",paste(na.cols,collapse=","),") are not specified in 'order' (",paste(colnames(rsc),collapse=","),")")
      }
      # n1 <- n %>%
      #   tibble::rownames_to_column("smpl") %>%
      #   mutate(smpl=sub("\\..+","",smpl)) %>%
      #   left_join(pd[c("smpl",summarize.by)],by="smpl") %>%
      #   group_by(!!!syms(summarize.by)) %>%
      #   summarize_at(used.abs,mean,na.rm=T)
      n1 <- cbind(n,pd) %>%
          group_by(!!!syms(summarize.by)) %>%
          summarize_at(used.abs,mean,na.rm=T)
      rn <- apply(as.matrix(n1[summarize.by]),1,paste,collapse=",")

      n1 <- as.matrix(n1[used.abs])
      rownames(n1) <- rn

      rsc1 <- unique(cbind(pd[summarize.by],rsc[,summarize.by]))
      names(rsc1) <- paste0(rep(c("pd",""),each=length(summarize.by)),rep(summarize.by,times=2))
      rsc1 <- rsc1 %>% arrange(!!!syms(paste0("pd",summarize.by))) %>% select(!!! syms(summarize.by))
      rsc <- as.matrix(rsc1)
      rsl <- colnames(rsc)

      ## reverse row-order
      row.o <- rev(seq(nrow(n1)))
      n1 <- n1[row.o,]
      rsc <- rsc1[row.o,]
    }else{
      row.o <- meta$row.order
      n1 <- as.matrix(n[row.o,])
      rownames(n1) <- rep("",nrow(n1))
      rsl <- colnames(rsc)
    }

    if(missing(ttl)){
      ttl <- paste0(ld_name,", ",sub("_.+","",norm_type)," (no scaling)")
    }
    h3(n1,Colv=NA,Rowv=NA,scale="none",RowSideColors=rsc,useRaster=F,main=ttl,margins=margins,RowSideLabs=rsl,...)
  })

