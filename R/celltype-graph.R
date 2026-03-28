# fun: CellTypeGraph ctype ----

#' Generating the levels of cell types and their hierarchy and plotting the cell type hirarchy
#'
#' @param ctype A data frame containing information about parent-child relationships between cell types.
#' The object should have at least two columns named Parent and Child, and each row contains the relationship between
#' direct parent and child cell types.
#' @param plot Logical, indicating whether to plot the graph.
#' @param transpose Logical, indicating whether to transpose the graph layout.
#' @param with.hierarchy Logical, indicating whether to include the cell type hierarchy information in output.
#' @param ... Additional parameters to be passed to the \code{plot} function.
#'
#' @return
#' If \code{with.hierarchy} is \code{TRUE}, returns a list with components \code{ctlevs} (cell type levels)
#' and \code{hierarchy} (cell type hierarchy data frame). If \code{with.hierarchy} is \code{FALSE},
#' returns only \code{ctlevs}.
#'
#' @details
#' This function generates a graph that represents the hierarchy of cell types based on parent-child relationships provided as a user input.
#' If the user-input is a Cycif or CycifStack object, the function extracts cell_lineage_def data frame and use it as the input.
#' If the input is a data frame, the function assumes that's the definition data frame. It returns the numerical level of each cell type within the hierarchy.
#' If \code{with.hierarchy} is \code{TRUE}, the hierarchy as a data frame is included in the output.
#' sers can also generate the cell type hierarchy as a graphical output when \code{plot} is \code{TRUE}.
#'
#' @seealso
#' \code{\link{graph_from_data_frame}}, \code{\link{layout_as_tree}}, \code{\link{distances}}
#'
#' @importFrom igraph graph_from_data_frame layout_as_tree distances
#'
#' @export
CellTypeGraph <- function(ctype,
                          cname="default",
                          plot=FALSE,
                          transpose=TRUE,
                          with.hierarchy=FALSE,...){
  if(is(ctype,"Cycif") | is(ctype,"CycifStack")){
    ctype <- x@cell_types[[cname]]@cell_lineage_def
  }else if(is(ctype,"data.frame")){
    if(!all(c("Parent","Child") %in% names(ctype))){
      stop("If 'ctype' argument is a data.frame, it should contain two columns named 'Parent' and 'Child'.")
    }
  }

  uniq.cts <- c("all",ctype$Child)
  ctgraph <- ctype[c("Parent","Child")]
  ctgraph$Parent <- factor(ctgraph$Parent,levels=uniq.cts)
  ctgraph$Child <- factor(ctgraph$Child,levels=uniq.cts)
  g <- igraph::graph_from_data_frame(ctgraph)
  l <- igraph::layout_as_tree(g)

  levs <- l[,2]
  names(levs) <- names(igraph::V(g))
  ctlevs <- rev(tapply(names(levs),levs,function(x)uniq.cts[uniq.cts %in% x]))

  if(plot){
    if(transpose){
      l[,2] <- max(l[,2]) - l[,2]
      l <- l[,2:1]
    }
    igraph::V(g)$shape <- "rectangle"
    igraph::V(g)$label.family <- "Helvetica"

    plot(g, layout=l,
         edge.arrow.size=.5, vertex.color="gold", vertex.size=40,
         vertex.frame.color=NA, vertex.label.color="black",
         vertex.label.cex=0.8, vertex.label.dist=0, edge.curved=0
         ,...)
  }
  if(with.hierarchy){
    d <- igraph::distances(g,mode="out")[uniq.cts,uniq.cts]

    cts.hierarchy <- as.data.frame(apply(which(d > 0 & d < Inf,arr.ind=T),2,function(idx)uniq.cts[idx]))
    names(cts.hierarchy) <- c("ancestor","descendant")

    leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
    non.leaves <- unique(ctype$Parent[!ctype$Parent %in% leaves])

    cts.hierarchy <- cts.hierarchy %>% filter(ancestor %in% non.leaves & descendant %in% leaves) ## 41 -> 33
    cts.hierarchy <- cts.hierarchy
    return(list(ctlevs=ctlevs,hierarchy=cts.hierarchy))
  }else{
    return(ctlevs)
  }
}

#_ -------------------------------------------------------
# fun: cellTypeFrequency CycifStack ----

#' @title Calculate Cell Type Frequency
#'
#' @description This function calculates the frequency of different cell types in a CycifStack object. It can return the frequency of individual cell types or a hierarchical structure of cell type frequencies.
#'
#' @param x A CycifStack object.
#' @param ct_name The name of the cell type definition to use. Default is "default."
#' @param simple Logical, if TRUE, return a simple table of cell type frequencies. If FALSE, return a hierarchical structure of cell type frequencies. Default is TRUE.
#' @param count Logical, if TRUE, return the count of cells per cell type. If FALSE, return the frequency of cells per cell type.
#'
#' @return If `simple` is TRUE (default), a matrix with row names representing samples and column names representing cell types (that are leaf nodes in the cell type definition tree), with the number of cells per cell type. If `simple` is FALSE, a list of 'non-leaf' cell type frequency matrices is also returned in the matrix.
#'
#' @details
#' The `cellTypeFrequency` function calculates the frequency of different cell types within a CycifStack object. It allows you to specify a cell type definition using the `ct_name` parameter. By default, it uses the "default" cell type definition. You can choose to return a simple table of cell type frequencies or a hierarchical structure of cell type frequencies.
#'
#' If `simple` is TRUE, the function returns a matrix with rows representing samples and columns representing cell types, along with their corresponding frequencies. If `simple` is FALSE, it returns a hierarchical structure where each level of the hierarchy represents a different cell type grouping.
#'
#' @export
#' @seealso \code{\link{CellTypeGraph}}

setGeneric("cellTypeFrequency", function(x,...) standardGeneric("cellTypeFrequency"))

#' @export
setMethod("cellTypeFrequency", "CycifStack",
  function(x,ct_name="default",simple=TRUE,count=FALSE,cts.hierarchy){
    tab <- table(cell_types(x,ct_name=ct_name))
    mat <- matrix(tab,nrow=nrow(tab),dimnames=list(sample=rownames(tab),cell_types=colnames(tab)))
    cts1 <- colnames(mat)
    cts1 <- cts1[cts1 != "outOfROI"]
    if(!count){
      nsh <- apply(mat[,cts1],1,function(x){
        x <- x
        x/sum(x)
      })
    }else{
      nsh <- t(mat[,cts1])
    }

    if(any(rownames(nsh)=="unknown")){
      stop("'unknown' as a cell type is discontinued; run defineCellTypes() again to update the cell types")
    }

    if(simple){
      return(nsh)
    }else{
      ctype <- x@cell_types[[ct_name]]@cell_lineage_def
      leaves <- ctype$Child[!ctype$Child %in% ctype$Parent]
      names(leaves) <- leaves

      if(missing(cts.hierarchy)){
        ct.graph <- CellTypeGraph(ctype,plot=F,main='Cell type definition',with.hierarchy=TRUE)
        cts.hierarchy <- ct.graph$hierarchy
      }


      list.hie <- tapply(cts.hierarchy$descendant,cts.hierarchy$ancestor,identity)
      list.hie1 <- c(list.hie,as.list(leaves))

      nsh1 <- t(sapply(list.hie1,function(cts){
        colSums(nsh[cts,,drop=F])
      })) # `all` shouldn't be included
      if(count){
        return(nsh1)
      }

      ##
      uniq.cts <- c("all",ctype$Child)
      ctgraph <- ctype[c("Parent","Child")]
      ctgraph$Parent <- factor(ctgraph$Parent,levels=uniq.cts)
      ctgraph$Child <- factor(ctgraph$Child,levels=uniq.cts)
      g <- igraph::graph_from_data_frame(ctgraph)
      vts <- igraph::V(g)$name


      list.hie2 <- lapply(names(list.hie),function(ct){
        cts <- vts[find_descendants(g,ct)]
      })
      names(list.hie2) <- names(list.hie)

      ## frequency over `all`
      # all.nsh <- nsh1
      # rownames(all.nsh) <- paste0(rownames(all.nsh),":all")

      ##
      list.nshs <- lapply(names(list.hie2),function(ct){
        vs <- list.hie2[[ct]]
        tmp <- nsh1[ct,,drop=F]
        # if(!all(vs %in% rownames(nsh1))){
        #   stop(ct)
        # }
        nsh2 <- nsh1[vs,,drop=F]/tmp[col(nsh1[vs,,drop=F])]
        nsh2[is.nan(nsh2)] <- 0
        rownames(nsh2) <- paste0(rownames(nsh2),":",ct)
        return(nsh2)
      })
      names(list.nshs) <- names(list.hie2)
      return(list.nshs)
    }
})

# fun: find_descendants graph, node ----
# this should be called internally; not exported
find_descendants <- function(graph, node) {
  neighbors <- igraph::neighbors(graph, node, mode = "out")
  descendants <- numeric(0)

  while (length(neighbors) > 0) {
    descendants <- unique(c(descendants, neighbors))
    new_neighbors <- numeric(0)

    for (neighbor in neighbors) {
      new_neighbors <- c(new_neighbors, igraph::neighbors(graph, neighbor, mode = "out"))
    }

    neighbors <- setdiff(new_neighbors, descendants)
  }

  return(descendants)
}
