policy_fold_graph <- function(edges, policy, columns = character()) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  ids <- names(graph)
  rules <- stats::setNames(rep("yes", length(ids)), ids)
  rules[names(policy)] <- policy
  attr(graph, "query_record") <- list(
    arguments = list(crop = "barley", params = list(query_depth = Inf)),
    policy_applied = TRUE,
    policy = data.frame(node_id = names(rules), recovered_node = unname(rules)),
    node_map = data.frame(node = ids, quantity = ids, member_id = ids, occurrence = ids),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))))
  graph
}

test_that("constructor fallback bypasses metadata exclusions without latent projection", {
  graph <- policy_fold_graph("T -> X -> Z -> Y; U -> T; U -> Y",
                             c(X = "yes:user", Z = "no:Interpretive node", U = "yes:user"))
  before <- serialize(graph, NULL)
  prepared <- .cd_prepare(graph, latents = NULL, keep_nodes = NULL)
  expect_setequal(names(prepared$analysis), c("T", "Y"))
  expect_setequal(edge_key(dag_edges(prepared$analysis)), edge_key(dag_edges(dagitty::dagitty("dag { T -> Y }"))))
  expect_setequal(names(prepared$selection), names(prepared$analysis))
  expect_length(dagitty::latents(prepared$selection), 0L)
  expect_setequal(names(prepared$reference), names(graph))
  expect_identical(serialize(graph, NULL), before)
  total <- .cd_search(prepared, "total", 100L)
  expect_identical(total$sets, list(character()))
  expect_warning(constructed <- construct_dag(graph, quiet = TRUE), "Recovery policy exclusions")
  expect_identical(serialize(original_query(constructed), NULL), before)
  expect_identical(dag_checks(constructed)$total$status, "preserved")
})

test_that("constructor policy folding preserves existing retained edge types and statuses", {
  graph <- policy_fold_graph("T -> X -> M -> Y; U -> M; U -> Y; T <-> Y",
                             c(X = "yes:user", U = "yes:user"))
  prepared <- .cd_prepare(graph, latents = "M", keep_nodes = NULL)
  edges <- dag_edges(prepared$analysis)
  expected <- dag_edges(dagitty::dagitty("dag { T -> M -> Y; T <-> Y }"))
  expect_setequal(edge_key(edges), edge_key(expected))
  expect_identical(dagitty::latents(prepared$analysis), "M")
  expect_identical(dagitty::latents(prepared$selection), "M")
  expect_identical(dagitty::exposures(prepared$analysis), "T")
  expect_identical(dagitty::outcomes(prepared$analysis), "Y")
})

test_that("constructor policy folding honours explicit requests and rejects unsupported incidents", {
  graph <- policy_fold_graph("T -> X -> Y", c(X = "yes:user", T = "yes:user"))
  expect_true("X" %in% names(.cd_prepare(graph, NULL, "X")$analysis))
  columns <- policy_fold_graph("T -> X -> Y", c(X = "yes:user"), columns = "X")
  expect_true("X" %in% names(.cd_prepare(columns, NULL, NULL)$analysis))
  forbidden <- policy_fold_graph("T -> X -> Y", c(X = "no:Interpretive node"))
  expect_error(.cd_prepare(forbidden, NULL, "X"), "saved recovery policy")
  mixed <- policy_fold_graph("T -> Y; U <-> Y", c(U = "yes:user"))
  expect_error(.cd_prepare(mixed, NULL, NULL), "edges that are not directed and are incident to excluded nodes")
})
