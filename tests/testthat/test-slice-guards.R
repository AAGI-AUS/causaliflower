test_that("only the complete matching slice identity enables source ordering", {
  nodes <- apsim_dag_nodes
  source_order <- attr(nodes, "cf_source_order", exact = TRUE)
  expect_identical(.query_source_order(nodes), source_order)
  expect_identical(.query_source_order(nodes[nrow(nodes):1L, ]), source_order)
  expect_null(.query_source_order(nodes[-1L, ]))
  changed <- nodes
  changed$concept_id[1:2] <- rev(changed$concept_id[1:2])
  expect_null(.query_source_order(changed))
  changed <- nodes
  attr(changed, "cf_source_order")$ranks$base_order[1L] <- NA_integer_
  expect_null(.query_source_order(changed))
  plain <- nodes
  attr(plain, "cf_source_order") <- NULL
  expect_null(.query_source_order(plain))
})

test_that("readable keys and source metadata are retained by the default query", {
  expect_null(formals(query_graph)$label)
  graph <- query_graph("barley", treatments = "N_application", outcomes = "grain_yield", quiet = TRUE)
  expect_true("nitrogen_applied" %in% names(graph))
  expect_identical(attr(graph, "cf_source_order"), attr(apsim_dag_nodes, "cf_source_order"))
  expect_identical(attr(graph, "query_record")$source_order, attr(graph, "cf_source_order"))
})

test_that("the slice refuses unavailable grouping and unresolved requested quantities", {
  expect_error(query_graph("barley", treatments = "nitrogen_applied", outcomes = "grain_yield",
    params = query_params(apsim_nodes_detail = "subnode"), quiet = TRUE), "readable variable keys")
  expect_error(query_graph("barley", treatments = "nitrogen_applied", outcomes = "grain_yield",
    params = query_params(supernode_overrides = c(nitrogen_applied = "node")), quiet = TRUE), "without supernode_overrides")
  expect_error(query_graph("barley", treatments = c("nitrogen_applied", "unknown_treatment"),
    outcomes = "grain_yield", quiet = TRUE), "resolve to no single node.*unknown_treatment")
  expect_error(query_graph("barley", treatments = "nitrogen_applied",
    outcomes = c("grain_yield", "unknown_outcome"), quiet = TRUE), "resolve to no single node.*unknown_outcome")
  expect_error(query_graph("barley", treatments = "nitrogen_applied", outcomes = "grain_yield",
    columns = c(measured = "unknown_column"), quiet = TRUE), "match no node.*unknown_column")
})

test_that("valid imported recovery metadata does not depend on the slice's seven keys", {
  nodes <- data.frame(node_id = c("first", "second", "third"), concept_id = c("c:one", "c:two", "c:three"),
    recovered_node = c("yes", "yes:user", "no:Not a numerical quantity"), stringsAsFactors = FALSE)
  before <- serialize(nodes, NULL)
  expect_identical(.cf_effective_nodes(nodes), nodes)
  expect_identical(serialize(nodes, NULL), before)
})
