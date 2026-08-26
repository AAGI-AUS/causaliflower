edge_keys <- function(dag) {
  edges <- as.data.frame(dagitty::edges(dag))
  if (nrow(edges) == 0L) {
    return(character())
  }
  edges <- edges[, c("v", "e", "w"), drop = FALSE]
  sort(paste(edges$v, edges$e, edges$w))
}
