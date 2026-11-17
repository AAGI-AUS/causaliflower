candidate_member_query <- function(edges, members = NULL, quantities = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  nodes <- names(graph)
  if (is.null(members)) members <- stats::setNames(nodes, nodes)
  if (is.null(quantities)) quantities <- members
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = nodes, quantity = unname(quantities[nodes]),
      member_id = unname(members[nodes]), occurrence = nodes),
    policy = data.frame(node_id = unique(unname(members)), recovered_node = "yes"))
  graph
}

test_that("a shared quantity label does not add source members to a candidate", {
  graph <- candidate_member_query("P -> T; P -> Y; Q -> Y; T -> Y",
    quantities = c(P = "Alias", Q = "Alias", T = "T", Y = "Y"))
  for (effect in c("direct", "total")) {
    search <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), effect, 500L)
    expect_identical(search$table$resolved, list("P"))
    expect_identical(search$table$raw_sets, list(list("P")))
    expect_identical(search$table$n_resolved, 1L)
    expect_identical(search$table$concepts, list("Alias"))
    expect_identical(search$table$n_concepts, 1L)
  }
  testthat::local_mocked_bindings(query_graph = function(...) graph, .package = "causaliflower")
  dag <- get_dag("barley", "T", "Y", FALSE, quiet = TRUE)
  expect_setequal(names(temporal_dag(dag)), c("P", "T", "Y"))
  for (effect in c("direct", "total")) {
    expect_identical(dag_options(dag)$quantities[[effect]]$members, "P")
    expect_identical(dag_checks(dag)[[effect]]$status, "preserved")
  }
})

test_that("distinct member alternatives with one label stay separate and correctly ranked", {
  graph <- candidate_member_query("P -> T; P -> Q; Q -> Y; T -> Y",
    quantities = c(P = "Alias", Q = "Alias", T = "T", Y = "Y"))
  for (effect in c("direct", "total")) {
    search <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), effect, 500L)
    expect_identical(search$table$resolved, list("P", "Q"))
    expect_identical(search$table$raw_sets, list(list("P"), list("Q")))
    expect_identical(search$table$n_resolved, c(1L, 1L))
    expect_identical(search$table$label, c("Alias", "Alias"))
    expect_identical(search$table$option, 1:2)
  }
  testthat::local_mocked_bindings(query_graph = function(...) graph, .package = "causaliflower")
  dag <- get_dag("barley", "T", "Y", FALSE, quiet = TRUE)
  for (effect in c("direct", "total")) {
    expect_identical(dag_options(dag)$quantities[[effect]]$set_size, 1)
    expect_setequal(dag_options(dag)$quantities[[effect]]$members, c("P", "Q"))
    expect_identical(dag_checks(dag)[[effect]]$status, "preserved")
  }
})

test_that("one canonical member retains all alias and occasion witnesses", {
  graph <- candidate_member_query("T -> Y; P_t0; P_t1",
    members = c(P_t0 = "P", P_t1 = "P", T = "T", Y = "Y"),
    quantities = c(P_t0 = "first alias", P_t1 = "second alias", T = "T", Y = "Y"))
  raw <- list("P_t1", "P_t0", c("P_t0", "P_t1"))
  table <- .cd_candidates(raw, .cd_prepare(graph, NULL, NULL))
  expect_identical(table$resolved, list("P"))
  expect_identical(table$n_resolved, 1L)
  expect_identical(table$raw_sets, list(raw))
  expect_identical(table$concepts, list("second alias"))
})

test_that("recombined nodes expand their own rows without following shared group labels", {
  graph <- candidate_member_query("T -> G_t0 -> Y; T -> Y; G")
  record <- attr(graph, "query_record")
  record$node_map$quantity[record$node_map$node == "G_t0"] <- "G"
  record$node_map$member_id[record$node_map$node == "G_t0"] <- "A"
  extra <- record$node_map[record$node_map$node == "G_t0", ]
  extra$member_id <- "B"
  record$node_map <- rbind(record$node_map, extra)
  attr(graph, "query_record") <- record
  table <- .cd_candidates(list("G_t0", "G"), .cd_prepare(graph, NULL, NULL))
  expect_identical(table$resolved, list("G", c("A", "B")))
  expect_identical(table$n_concepts, c(1L, 1L))
  expect_identical(table$n_resolved, c(1L, 2L))
  expect_identical(table$raw_sets, list(list("G"), list("G_t0")))
  empty <- .gd_candidates(list(character()), graph)
  absent <- .gd_candidates(list(), graph)
  expect_identical(empty$resolved, list(character()))
  expect_identical(empty$n_resolved, 0L)
  expect_identical(nrow(absent), 0L)
})

test_that("member-set keys cannot collide through separators in IDs", {
  graph <- candidate_member_query("T -> Y; A; B; C",
    members = c(A = "a\rb", B = "a", C = "b", T = "T", Y = "Y"))
  table <- .gd_candidates(list("A", c("B", "C")), graph)
  expect_identical(nrow(table), 2L)
  expect_identical(table$n_resolved, c(1L, 2L))
  expect_identical(table$raw_sets, list(list("A"), list(c("B", "C"))))
})
