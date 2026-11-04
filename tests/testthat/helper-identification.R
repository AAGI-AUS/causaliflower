identification_result <- function(dag, ...) {
  suppressMessages(get_identification(dag, ..., quiet = TRUE))
}

expect_identification_error <- function(code, message) {
  error <- testthat::capture_error(code)
  testthat::expect_true(inherits(error, "error"))
  testthat::expect_identical(conditionMessage(error), message)
}

identification_strategy <- function(result, pattern) {
  index <- grep(pattern, result$strategy, ignore.case = TRUE)
  testthat::expect_length(index, 1L)
  result[index]
}

identification_edge_keys <- function(dag) {
  edges <- as.data.frame(dagitty::edges(dag))
  if (nrow(edges) == 0L) {
    return(character())
  }
  sort(paste(edges$v, edges$e, edges$w))
}

identification_roots <- function(dag) {
  edges <- as.data.frame(dagitty::edges(dag))
  directed <- edges[edges$e == "->", , drop = FALSE]
  sort(setdiff(names(dag), as.character(directed$v)))
}

identification_is_district <- function(dag) {
  vertices <- names(dag)
  if (length(vertices) <= 1L) {
    return(TRUE)
  }

  edges <- as.data.frame(dagitty::edges(dag))
  bidirected <- edges[edges$e == "<->", , drop = FALSE]
  if (nrow(bidirected) == 0L) {
    return(FALSE)
  }

  adjacency <- stats::setNames(lapply(vertices, function(x) character()),
                               vertices)
  for (i in seq_len(nrow(bidirected))) {
    left <- as.character(bidirected$v[i])
    right <- as.character(bidirected$w[i])
    adjacency[[left]] <- c(adjacency[[left]], right)
    adjacency[[right]] <- c(adjacency[[right]], left)
  }

  reached <- vertices[1L]
  frontier <- reached
  while (length(frontier) > 0L) {
    neighbours <- unique(unlist(adjacency[frontier], use.names = FALSE))
    frontier <- setdiff(neighbours, reached)
    reached <- union(reached, frontier)
  }
  setequal(reached, vertices)
}

expect_identification_shape <- function(result) {
  testthat::expect_s3_class(result, "data.table")
  testthat::expect_identical(
    names(result),
    c("strategy", "feasible", "tier", "variable_set", "note")
  )
  testthat::expect_gt(nrow(result), 0L)
  testthat::expect_type(result$strategy, "character")
  testthat::expect_type(result$feasible, "logical")
  testthat::expect_type(result$tier, "character")
  testthat::expect_type(result$variable_set, "character")
  testthat::expect_type(result$note, "character")

  scalar_logical_attributes <- c(
    "identifiable", "observationally_identifiable", "z_identifiable",
    "conditional"
  )
  for (attribute in scalar_logical_attributes) {
    value <- attr(result, attribute)
    testthat::expect_type(value, "logical")
    testthat::expect_length(value, 1L)
  }

  testthat::expect_type(attr(result, "data_availability"), "character")
  testthat::expect_length(attr(result, "data_availability"), 1L)
  testthat::expect_true(is.null(attr(result, "hedge")) ||
                        is.list(attr(result, "hedge")))
  testthat::expect_s3_class(attr(result, "roles"), "data.table")
  testthat::expect_type(attr(result, "treatment"), "character")
  testthat::expect_type(attr(result, "outcome"), "character")
  testthat::expect_length(attr(result, "treatment"), 1L)
  testthat::expect_length(attr(result, "outcome"), 1L)
  testthat::expect_true(attr(result, "effect") %in% c("total", "direct"))
  testthat::expect_type(attr(result, "conditioning"), "character")
  testthat::expect_type(attr(result, "conditioning_intervened"), "character")
  testthat::expect_type(attr(result, "conditioning_remaining"), "character")
  testthat::expect_type(attr(result, "caveat"), "character")
  testthat::expect_length(attr(result, "caveat"), 1L)
}

expect_valid_identification_hedge <- function(result, source_dag) {
  hedge <- attr(result, "hedge")
  testthat::expect_type(hedge, "list")
  testthat::expect_named(hedge, c("F", "F_prime", "roots", "query"))
  testthat::expect_s3_class(hedge$F, "dagitty")
  testthat::expect_s3_class(hedge$F_prime, "dagitty")
  testthat::expect_type(hedge$roots, "character")
  testthat::expect_type(hedge$query, "list")

  testthat::expect_setequal(identification_roots(hedge$F), hedge$roots)
  testthat::expect_setequal(identification_roots(hedge$F_prime), hedge$roots)
  testthat::expect_true(identification_is_district(hedge$F))
  testthat::expect_true(identification_is_district(hedge$F_prime))
  testthat::expect_true(all(names(hedge$F_prime) %in% names(hedge$F)))

  for (forest in list(hedge$F, hedge$F_prime)) {
    edges <- as.data.frame(dagitty::edges(forest))
    directed <- edges[edges$e == "->", , drop = FALSE]
    child_count <- table(as.character(directed$v))
    testthat::expect_true(length(child_count) == 0L ||
                          all(child_count <= 1L))
  }

  f_edges <- as.data.frame(dagitty::edges(hedge$F))
  f_keys <- identification_edge_keys(hedge$F)
  fp_keys <- identification_edge_keys(hedge$F_prime)
  in_f_prime <- as.character(f_edges$v) %in% names(hedge$F_prime) &
    as.character(f_edges$w) %in% names(hedge$F_prime)
  expected_fp <- sort(paste(f_edges$v[in_f_prime], f_edges$e[in_f_prime],
                            f_edges$w[in_f_prime]))
  testthat::expect_setequal(fp_keys, expected_fp)
  testthat::expect_true(all(f_keys %in% identification_edge_keys(source_dag)))

  intervention <- hedge$query$intervention
  testthat::expect_type(intervention, "character")
  testthat::expect_true(any(intervention %in% names(hedge$F)))
  testthat::expect_false(any(intervention %in% names(hedge$F_prime)))
}
