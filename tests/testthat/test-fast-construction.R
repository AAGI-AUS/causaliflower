fast_query <- function(edges, columns = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  nodes <- names(graph)
  ids <- .strip_slice(nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = columns, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = ids, member_id = ids, occurrence = nodes),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  graph
}


# The ranked route (quantity-ranked-v1) has its own coverage below.

test_that("construction collapses repeated occasions and keeps the exact source membership", {
  graph <- fast_query("T -> A_t1 -> Y; T -> A_t2 -> Y; T -> Y")
  result <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_setequal(names(result), c("T", "Y", "A"))
  expect_true(dagby_options(result)$summary$applied)
  checks <- dagby_checks(result)
  expect_identical(checks$direct$status, "preserved")
  expect_identical(dagby_options(result)$quantities$direct$members, "A")
  expect_identical(.gd_adjustment_sets(result, "direct", Inf)$sets, list("A"))
  search <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), "direct")
  expect_identical(search$ranked$certificate, "complete")
  expect_identical(search$sets, list(c("A_t1", "A_t2")))
  expect_true(all(.gd_validate_sets(.cd_prepare(graph, NULL, NULL)$selection, search$sets, "direct")))
})


test_that("a ranked search error stops the size route and never falls back", {
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  testthat::local_mocked_bindings(.cd_ranked_search = function(...) stop("ranked failure"), .package = "causaliflower")
  expect_error(dagby_construct(graph, 0, quiet = TRUE), "ranked failure")
})
