# Preserve the scientific fixtures' query inputs while exercising the public recipe route
recipe_test_query <- function(crop, columns = NULL, treatments = NULL, outcomes = NULL,
                              randomised_treatment = FALSE, time_unroll = NULL,
                              label = NULL, quiet = FALSE, params = dagby_query_params()) {
  arguments <- list(crop = crop, treatments = treatments, outcomes = outcomes,
    randomised_treatment = randomised_treatment, columns = columns,
    label = if (is.null(label)) "id" else label, params = params)
  if (!missing(time_unroll)) arguments["time_unroll"] <- list(time_unroll)
  recipe <- do.call(dagby, arguments)
  return(dagby_query(recipe, quiet = quiet))
}

# Small enumerable source with C retained only at the largest size unless requested
recipe_fixture_source <- function(crop, columns = NULL, treatments = NULL,
                                  outcomes = NULL, randomised_treatment = FALSE,
                                  time_unroll = NULL, label = NULL, quiet = FALSE,
                                  params = dagby_query_params()) {
  route <- if (identical(time_unroll, "all")) "T -> A_t0 -> A_t1 -> Y" else "T -> A -> Y"
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", route, "; T -> C; T -> Y }"))
  nodes <- names(graph)
  quantities <- sub("_t[01]$", "", nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = crop, columns = columns, treatments = treatments,
      outcomes = outcomes, randomised_treatment = randomised_treatment,
      time_unroll = time_unroll, label = label, params = params),
    time_unroll_defaulted = missing(time_unroll),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = quantities, member_id = quantities,
      occurrence = nodes),
    policy = data.frame(node_id = unique(quantities), recovered_node = "yes"))
  return(graph)
}
