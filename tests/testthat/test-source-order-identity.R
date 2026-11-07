# ==============================================================================
# test-source-order-identity.R


# Invalidate saved searches when the implementation of their recorded ordering changes
# ==============================================================================


ordered_identity_graph <- function(){ # Make a small query with a complete recorded ordering and two adjustment choices
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> B -> T; A -> Y; T -> Y }")
  nodes <- names(graph)
  count <- length(nodes)
  source_order <- list(version = 1L, source = "apsim-slice-2026-11-05",
    nodes = data.frame(node_id = nodes, concept_id = nodes),
    ranks = data.frame(node_id = nodes, base_order = seq_along(nodes),
      at_order = count + seq_along(nodes), post_order = 2L * count + seq_along(nodes),
      pre_order = 3L * count + seq_along(nodes), time_order = 4L * count + seq_along(nodes)))
  attr(graph, "cf_source_order") <- source_order
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", params = list(water_model = "SoilWater", query_depth = Inf)),
    requests = list(columns = list(resolved = character())),
    node_map = data.frame(node = nodes, quantity = nodes, member_id = nodes, occurrence = nodes),
    policy = data.frame(node_id = nodes, recovered_node = "yes"), source_order = source_order)
  class(graph) <- c("dag_by_query", class(graph))
  return(graph)
}


test_that("ordering implementation belongs only to a scoped search identity", {
  graph <- ordered_identity_graph()
  identity <- .cd_search_identity_value(.cd_prepare(graph, NULL, NULL))
  expect_identical(identity$ordering_adapter, list(
    order_key = deparse(body(.query_order_key)), sort = deparse(body(.query_sort)),
    applies = deparse(body(.query_order_applies)), set_key = deparse(body(.query_set_key))))
  attr(graph, "cf_source_order") <- NULL
  attr(graph, "query_record")$source_order <- NULL
  generic <- .cd_search_identity_value(.cd_prepare(graph, NULL, NULL))
  expect_false("source_order" %in% names(generic))
  expect_false("ordering_adapter" %in% names(generic))
})


test_that("unchanged ordering helpers reuse completed searches", {
  graph <- ordered_identity_graph()
  previous <- construct_dag(graph, quiet = TRUE, search = "complete")
  search_calls <- 0L
  source_search <- .cd_search
  result <- testthat::with_mocked_bindings(
    .cd_reconstruct(graph, previous = previous, quiet = TRUE, search = "complete"),
    .cd_search = function(...){ # Count source enumeration without changing the fingerprinted engine
      search_calls <<- search_calls + 1L
      return(source_search(...))
    })
  expect_identical(search_calls, 0L)
  expect_identical(result, previous)
})


test_that("changed ordering helper code invalidates reuse even when its outputs agree", {
  graph <- ordered_identity_graph()
  previous <- construct_dag(graph, quiet = TRUE, search = "complete")
  previous_identity <- .gd_state(previous)$construction$search_identity
  source_sort <- .query_sort
  source_search <- .cd_search
  search_calls <- 0L
  testthat::local_mocked_bindings(
    .query_sort = function(nodes, source_order = NULL, node_ids = NULL){ # Stand in for a revised ordering implementation
      result <- source_sort(nodes, source_order, node_ids)
      return(result)
    },
    .cd_search = function(...){ # Count both source searches after identity invalidation
      search_calls <<- search_calls + 1L
      return(source_search(...))
    })
  current_identity <- .cd_search_identity_value(.cd_prepare(graph, NULL, NULL))
  expect_false(identical(current_identity, previous_identity))
  expect_identical(current_identity$ordering_adapter$sort, deparse(body(.query_sort)))
  result <- .cd_reconstruct(graph, previous = previous, quiet = TRUE, search = "complete")
  expect_identical(search_calls, 2L)
  expect_identical(result, construct_dag(graph, quiet = TRUE, search = "complete"))
})
