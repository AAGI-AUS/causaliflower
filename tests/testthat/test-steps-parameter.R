# The ranked search's explicit step budget: construct_dag(steps), get_dag's
# construct_dag list, the recorded budget and select_dag's explicit replay.
# ranked_query is copied from test-ranked-search.R, where it is local.
ranked_query <- function(edges, columns = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  nodes <- names(graph)
  ids <- .strip_slice(nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = columns, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = if (is.null(columns)) list() else stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = ids, member_id = ids, occurrence = nodes),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  graph
}

test_that("the ranked search defaults to 20000 steps and an explicit value beats the option", {
  old <- options(causaliflower.ranked_steps = NULL)
  on.exit(options(old), add = TRUE)
  expect_null(getOption("causaliflower.ranked_steps"))
  expect_identical(.cd_ranked_steps(), 20000L)
  options(causaliflower.ranked_steps = 9L)
  expect_identical(.cd_ranked_steps(), 9L)
})

test_that("the ranked budget in force is used and recorded by the ranked search and the size route", {
  graph <- ranked_query("A -> T; A -> B; B -> Y; C -> T; C -> D; D -> Y; T -> Y")
  limited <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), "total", steps = 1L)
  expect_false(limited$complete)
  expect_identical(limited$ranked$certificate, "budget")
  expect_identical(limited$ranked$step_budget, 1L)
  withr::local_options(causaliflower.ranked_steps = 1L)
  full <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), "total", steps = 20000L)
  expect_true(full$complete)
  expect_identical(full$ranked$step_budget, 20000L)
  expect_identical(.cd_ranked_search(.cd_prepare(graph, NULL, NULL), "total")$ranked$step_budget, 1L)
  small <- suppressWarnings(construct_dag(graph, 0, quiet = TRUE))
  expect_identical(dag_options(small)$arguments$steps, 1L)
})


test_that("an invalid ranked budget option is refused by the size route", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y")
  withr::local_options(causaliflower.ranked_steps = "a")
  expect_error(construct_dag(graph, 0, quiet = TRUE), "Option causaliflower.ranked_steps must be", fixed = TRUE)
})
