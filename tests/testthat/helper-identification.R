################################################################################
# helper-identification.R


# Support private identification oracles and public report checks

# Computational checks use one private query at a time. Public behaviour is
# tested through get_identification().
identification_result <- function(dag, ...) { # Return the private record for one query
  return(causaliflower:::id_identification_query(dag, ..., quiet = TRUE))
}

expect_identification_error <- function(code, message) { # Check the complete validation error
  error <- testthat::capture_error(code)
  testthat::expect_true(inherits(error, "error"))
  testthat::expect_identical(conditionMessage(error), message)
  return(invisible(NULL))
}

identification_strategy <- function(result, pattern) { # Select one private strategy row
  index <- grep(pattern, result$strategy, ignore.case = TRUE)
  testthat::expect_length(index, 1L)
  return(result[index])
}

identification_edge_keys <- function(dag) { # Compare typed graph edges
  edges <- as.data.frame(dagitty::edges(dag))
  if (nrow(edges) == 0L) {
    return(character())
  }
  return(sort(paste(edges$v, edges$e, edges$w)))
}

identification_roots <- function(dag) { # Find nodes with no outgoing directed edge
  edges <- as.data.frame(dagitty::edges(dag))
  directed <- edges[edges$e == "->", , drop = FALSE]
  return(sort(setdiff(names(dag), as.character(directed$v))))
}

identification_is_district <- function(dag) { # Check bidirected connectivity for a hedge witness
  vertices <- names(dag)
  if (length(vertices) <= 1L) {
    return(TRUE)
  }

  edges <- as.data.frame(dagitty::edges(dag))
  bidirected <- edges[edges$e == "<->", , drop = FALSE]
  if (nrow(bidirected) == 0L) {
    return(FALSE)
  }

  adjacency <- stats::setNames(lapply(vertices, function(node) { return(character()) }),
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
  return(setequal(reached, vertices))
}

expect_identification_shape <- function(result) { # Check every private record column and attribute
  # Check the private worker's five columns
  testthat::expect_s3_class(result, "data.table")
  testthat::expect_false(inherits(result, "causaliflower_identification"))
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
    attribute_value <- attr(result, attribute)
    testthat::expect_type(attribute_value, "logical")
    testthat::expect_length(attribute_value, 1L)
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
  return(invisible(NULL))
}

expect_valid_identification_hedge <- function(result, source_dag) { # Validate both hedge forests against the source graph
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
  return(invisible(NULL))
}
