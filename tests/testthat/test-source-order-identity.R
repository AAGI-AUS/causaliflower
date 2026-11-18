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


test_that("changed ordering helper code changes the search identity even when its outputs agree", {
  graph <- ordered_identity_graph()
  previous_identity <- .cd_search_identity_value(.cd_prepare(graph, NULL, NULL))
  source_sort <- .query_sort
  testthat::local_mocked_bindings(
    .query_sort = function(nodes, source_order = NULL, node_ids = NULL){
      result <- source_sort(nodes, source_order, node_ids)
      return(result)
    }, .package = "causaliflower")
  current_identity <- .cd_search_identity_value(.cd_prepare(graph, NULL, NULL))
  expect_false(identical(current_identity, previous_identity))
  expect_identical(current_identity$ordering_adapter$sort, deparse(body(.query_sort)))
})

test_that("the recorded source order sorts quantities and the members dagby_construct keeps", {
  source_order <- list(version = 1L, ranks = data.frame(node_id = c("B", "T", "A", "Y"), base_order = 1:4,
    at_order = 1:4, post_order = 1:4, pre_order = 1:4, time_order = 1:4))
  expect_identical(.query_sort(c("A", "B", "T", "Y"), source_order), c("B", "T", "A", "Y"))
  expect_identical(.query_sort(c("A", "B", "T", "Y")), c("A", "B", "T", "Y"))
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> Y; B -> T; B -> Y; T -> Y }")
  nodes <- names(graph)
  attr(graph, "cf_source_order") <- source_order
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = character())),
    node_map = data.frame(node = nodes, quantity = nodes, member_id = nodes, occurrence = nodes),
    policy = data.frame(node_id = nodes, recovered_node = "yes"), source_order = source_order)
  expect_identical(dagby_options(dagby_construct(graph, quiet = TRUE, size = "largest"))$quantities$total$members, c("B", "A"))
})
