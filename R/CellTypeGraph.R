CellTypeGraph <- function(ctype,plot=F,transpose=T){
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
    plot(g,layout=l)
  }
  return(ctlevs)
}


