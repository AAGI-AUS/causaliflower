edge_keys <- function(dag) {
  edges <- as.data.frame(dagitty::edges(dag))
  if (nrow(edges) == 0L) {
    return(character())
  }
  edges <- edges[, c("v", "e", "w"), drop = FALSE]
  sort(paste(edges$v, edges$e, edges$w))
}


test_role_positions <- function(dag, role) {
  causaliflower:::existing_role_specs(dag)[[role]]$positions
}


role_mapping_sources <- function(dag, role) {
  causaliflower:::existing_role_specs(dag)[[role]]$mapping_sources
}
