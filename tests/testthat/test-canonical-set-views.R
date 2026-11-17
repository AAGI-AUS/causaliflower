set_view_graph <- function(text, columns = character()) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", text, "}"))
  ids <- .strip_slice(names(graph))
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "barley", params = list(water_model = "SoilWater", query_depth = Inf)),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = names(graph), quantity = ids, member_id = ids,
                          occurrence = names(graph), stringsAsFactors = FALSE),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  graph
}


test_that("summaries retain distinct endpoint occasions", {
  source <- dagitty::dagitty("dag { X_t2 [exposure]; Y [outcome]; X_t1 -> X_t2 -> A -> Y; X_t1 -> A }")
  nodes <- names(source)
  attr(source, "cf_node_map") <- data.frame(node = nodes, quantity = .strip_slice(nodes),
    member_id = .strip_slice(nodes), occurrence = nodes)
  expect_error(construct_dag(source, keep_nodes = "X_t1", quiet = TRUE),
    "distinct occasions of an endpoint", class = "cd_static_unavailable")
  dag <- suppressWarnings(construct_dag(source, keep_nodes = "X_t1", quiet = TRUE, unroll = TRUE))
  expect_true(all(c("X_t1", "X_t2") %in% names(dag)))
  expect_false(dag_options(dag)$summary$applied)
  expect_match(dag_options(dag)$level$static_reason, "distinct occasions of an endpoint")
  expect_identical(.gd_signature(dag), .gd_signature(temporal_dag(dag)))
})


test_that("time summaries do not invent an endpoint path by merging disconnected occasions", {
  source <- set_view_graph("T -> A_t1; A_t2 -> Y")
  proposed <- .gd_summary(source, keep = c("A_t1", "A_t2"))
  expect_false(proposed$applied)
  expect_match(proposed$reason, "proper directed path from treatment to outcome")
  expect_identical(.gd_signature(proposed$graph), .gd_signature(source))
  expect_error(construct_dag(source, keep_nodes = c("A_t1", "A_t2"), quiet = TRUE),
    "proper directed path", class = "cd_static_unavailable")
  dag <- construct_dag(source, keep_nodes = c("A_t1", "A_t2"), quiet = TRUE, unroll = TRUE)
  expect_false(dag_options(dag)$summary$applied)
  expect_setequal(names(dag), c("T", "A_t1", "A_t2", "Y"))
})

test_that("proper endpoint paths exclude other treatments from the interior", {
  source <- dagitty::dagitty("dag { X1 [exposure]; X2 [exposure]; Y [outcome]; X1 -> A_t1; A_t2 -> Y; X1 -> X2 -> Y }")
  paths <- .gd_endpoint_paths(source)
  expect_false(unname(paths[["4:2:X11:Y"]]))
  expect_true(unname(paths[["4:2:X21:Y"]]))
  proposed <- .gd_summary(source, keep = c("A_t1", "A_t2"))
  expect_false(proposed$applied)
  expect_match(proposed$reason, "proper directed path from treatment to outcome")
})

test_that("repeated compression cannot remove the only proper endpoint path", {
  source <- set_view_graph("T -> A_t0; B_t0 -> A_t0; A_t0 -> B_t1; B_t1 -> A_t1; B_t1 -> Y")
  proposed <- .gd_summary(source, keep = c("A_t0", "A_t1", "B_t0", "B_t1"))
  expect_true(unname(.gd_endpoint_paths(source)[["3:1:T1:Y"]]))
  expect_false(proposed$applied)
  expect_match(proposed$reason, "proper directed path from treatment to outcome")
  expect_identical(.gd_signature(proposed$graph), .gd_signature(source))
})
