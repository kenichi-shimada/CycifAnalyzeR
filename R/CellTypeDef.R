proj_dir <- "~/Dropbox (HMS)/_projects_/L01 TALAVE/_CyCIF_"
fig_dir <- file.path(proj_dir,"figs")
obj_dir <- file.path(proj_dir,"robj")

ctype.file <- file.path(proj_dir,"cell-type-def","cell-type-definition_TALAVE_120122.xlsx")
ctype.sheets <- getSheetNames(ctype.file)
ctype <- readWorkbook(ctype.file,sheet=ctype.sheets[2],colNames=TRUE,rowNames=FALSE)
cstate <- c("cCaspase3","pERK","pAKT","MCL1","BCLXL","pTBK1","gH2AX")

## leaves
leaves <- function(ctype){
  cts <- ctype$Child
  h <- ctype[c("Parent","Child")]
  leaves <- h$C[!h$C %in% h$P]
  return(leaves)
}

## expand cell types - CAN will be converted to AND and NOT (XOR is not implemented)
expand_celltype <- function(ctype){
  lvs <- leaves(ctype)
  cs <- names(which(sapply(ctype,function(ct)any(ct=="CAN"))))
  if(any (!cs %in%  c("PD1","PDL1"))){
    i.nonpd <- which(!cs %in%  c("PD1","PDL1"))
    cs <- cs[c(i.nonpd,seq(cs)[-i.nonpd])]
  }
  tmp <- ctype[cs]
  tmp[is.na(tmp)] <- "NA"

  ri <- which(apply(as.matrix(ctype[cs]),1,function(x)any(x=="CAN")))
  list.ct <- list()
  for(i in ri){
    cs1 <- names(tmp)[which(tmp[i,]=="CAN")]
    css <- apply(cbind(ctype$Child[i],expand.grid(lapply(cs1,function(x)paste0(x,c("+","-"))))),
                 1,paste0,collapse=",")

    n <- length(css)
    ct1 <- ctype[i,]
    ct1 <- ct1[rep(1,n),]
    ct1$Parent <- ct1$Child
    ct1$Child <- css

    ct2 <- expand.grid(lapply(cs1,function(i)c("AND","NOT")))
    names(ct2) <- cs1
    for(n in names(ct2)){
      ct1[[n]] <- ct2[[n]]
    }

    list.ct[[match(i,ri)]] <- ct1
  }
  lct <- do.call(rbind,c(list(ctype),list.ct))
  return(lct)
}

ct1 <- expand_celltype(ctype)
find_depth(ct1)
graph_hc(ct1)

used.abs <- as.character(cs1@uniq_abs$ab)


## make hierarchy
find_depth <- function(ctype,levels){
  h <- ctype[c("Parent","Child")]
  ucts <- unique(c(h$P,h$C))
  d <- rep(NA,length(ucts))
  names(d) <- ucts
  i <- 0
  chs1 <- "all"
  d[chs1] <- 0
  while(any(is.na(d))){
    i <- i+1
    chs <- h %>% filter(Parent %in% chs1) %>% select(Child) %>% .$Child
    chs1 <- names(which(is.na(d[chs])))
    d[chs1] <- i
  }
  leaves <- h$C[!h$C %in% h$P]
  if(!missing(levels)){
    if(!identical(sort(levels),sort(leaves))){
      stop("'levels' should include all the cell types.\n")
    }
  }
  if(any(table(leaves)>1)){
    stop("Cell types (leaves) should be unique names.\n")
  }
  return(list(depth=d,celltypes=leaves))
}

## graph_hc
graph_hc <- function(ctype){
  require(igraph)
  ct <- as.matrix(ctype[c("Parent","Child")])
  g <- igraph::graph_from_edgelist(ct)
  plot(g,layout=layout_as_tree)
}

##
#' Instanciate and show a CellTypeDef object
#'
#' @param filename Excel sheet that contains two spreadsheets, named 'cell lineage' and 'cell state', respectively.
#' @param object A CellTypeDef object
#'
#' @export

CellTypeDef <- function(x,ctype,cstate){
  ct.mks <- names(ctype)[-(1:2)]
  cs.mks <- cstate
  mks <- unique(c(ct.mks,cs.mks))

  is.ab.used <- cyApply(x,function(cy){
    all(mks %in% levels(cy@abs_list$ab))
  })

  if(!all(mks %in% )){
    stop("Cell types in `ctype' and `cstate' shoud be identical.")
  }

  if(!all(c("cell lineage","cell state") %in% sheet.names)){
    stop("Cell type and cell state should be defined in two spreadsheets named 'cell type' and 'cell state', respectively.")
  }

  clineage <-  openxlsx::readWorkbook(filename,sheet="cell lineage",colNames=TRUE,rowNames=TRUE)
  cstate <-  openxlsx::readWorkbook(filename,sheet="cell state",colNames=TRUE,rowNames=TRUE)

  r1 <- rownames(clineage)
  r2 <- rownames(cstate)
  c1 <- colnames(clineage)
  c2 <- colnames(cstate)

  if(!identical(r1,r2)){
    stop("row.names in 'cell lineage' and 'cell state' should be identical.")
  }else{
    lineages <- r1
    lineages <- gsub(", ","_",lineages)
    rownames(clineage) <- rownames(cstate) <- lineages
  }

  c12 <- intersect(c1,c2)
  if(length(c12)>0){
    stop("Abs defined as both lineage markers and state markers: ",paste(c12,collapse=", "))
  }else{
    markers <- list(lineage=c1,state=c2)
  }

  new("CellTypeDef",
      cell_lineage_df = clineage,
      cell_state_df = cstate,
      lineages = lineages,
      markers = markers
  )
}

#' @rdname CellTypeDef
#' @export
setMethod("show", "CellTypeDef", function(object) {
  lins <- lineages(object)
  n.lins <- length(lins)

  mks <- markers(object)
  n.mks <- length(unlist(mks))
  n.lm <- length(mks$lineage)
  n.sm <- length(mks$state)

  used.abs <- used_abs(object)
  n.uab <- length(used.abs)

  abs.list <- uniq_abs(object)

  ## first three examples
  lin.txt <- ifelse(n.lins >0,paste0(" (",paste(lins[seq(min(n.lins,3))],collapse=", "),",...)"),"")
  lm.txt <- ifelse(n.lm >0,paste0(" (",paste(mks$lineage[seq(min(n.lm,3))],collapse=", "),",...)"),"")
  sm.txt <- ifelse(n.sm >0,paste0(" (",paste(mks$state[seq(min(n.sm,3))],collapse=", "),",...)"),"")
  uab.txt <- ifelse(n.uab >0,paste0(" (",paste(used.abs[seq(min(n.uab,3))],collapse=", "),",...)"),"")

  cat(is(object)[[1]], "\n",
      " lineages:  ", n.lins, lin.txt,"\n",
      " markers:   ", n.mks, "\n",
      "   lineage: ", n.lm, lm.txt,"\n",
      "   state:   ", n.sm, sm.txt,"\n",
      " used_abs:  ", n.uab, uab.txt,"\n",
      " uniq_abs:  ",
      sep="")
  if(nrow(abs.list)>0){
    cat("\n")
    print(abs.list)
  }else{
    cat("not set\n")
  }
})
