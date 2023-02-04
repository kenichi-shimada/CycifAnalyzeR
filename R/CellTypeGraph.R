#'@export
CellTypeGraph <- function(ctype,plot=F,transpose=T,...){
  uniq.cts <- c("all",ctype$Child)
  ctgraph <- ctype[c("Parent","Child")]
  ctgraph$Parent <- factor(ctgraph$Parent,levels=uniq.cts)
  ctgraph$Child <- factor(ctgraph$Child,levels=uniq.cts)
  g <- igraph::graph_from_data_frame(ctgraph)
  l <- igraph::layout_as_tree(g)
  levs <- l[,2]
  names(levs) <- names(igraph::V(g))
  ctlevs <- rev(tapply(names(levs),levs,identity))
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
  return(ctlevs)
}


