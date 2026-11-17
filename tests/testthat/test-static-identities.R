static_identity_graph <- function(edges, members = NULL, printed = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  nodes <- names(graph)
  if (is.null(members)) members <- stats::setNames(nodes, nodes)
  if (is.null(printed)) printed <- members
  attr(graph, "cf_node_map") <- data.frame(node = nodes,
    quantity = unname(printed[nodes]), member_id = unname(members[nodes]),
    occurrence = nodes, stringsAsFactors = FALSE)
  graph
}

check_static_source_sets <- function(display, source) {
  map <- attr(display, "cf_node_map", exact = TRUE)
  checked <- 0L
  for (effect in c("direct", "total")) {
    sets <- .gd_adjustment_sets(.gd_plain(display), effect, Inf)$sets
    for (set in sets) {
      occurrences <- unique(map$occurrence[map$node %in% set])
      expect_true(all(.gd_validate_sets(source, list(occurrences), effect)), info = effect)
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 0L)
}

test_that("explicit distinct variable identities keep literal time-like names apart", {
  graph <- static_identity_graph("P_t0 -> T; P_t0 -> Y; P_t1 -> T; P_t1 -> Y; T -> Y")
  result <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  expect_setequal(names(result), names(graph))
  actual <- .gd_signature(result)
  expected <- .gd_signature(graph)
  actual$labels <- expected$labels <- NULL
  expect_identical(actual, expected)
  expect_identical(dag_options(result)$level$variables, 4L)
  expect_setequal(attr(result, "cf_node_map")$member_id, names(graph))
  check_static_source_sets(result, graph)
})

test_that("distinct variables stay apart alongside a genuinely repeated variable", {
  graph <- static_identity_graph(paste("P_t0 -> T; P_t0 -> Y; P_t1 -> T; P_t1 -> Y; T -> Y;",
    "A_t0 -> A_t1; A_t1 -> T; A_t1 -> Y"),
    c(P_t0 = "P_t0", P_t1 = "P_t1", A_t0 = "A", A_t1 = "A", T = "T", Y = "Y"))
  result <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  expect_true(dag_options(result)$summary$applied)
  expect_setequal(names(result), c("P_t0", "P_t1", "A", "T", "Y"))
  map <- attr(result, "cf_node_map")
  expect_setequal(map$occurrence[map$node == "A"], c("A_t0", "A_t1"))
  expect_setequal(map$occurrence[map$node == "P_t0"], "P_t0")
  expect_setequal(map$occurrence[map$node == "P_t1"], "P_t1")
  expect_identical(dag_options(result)$level$variables, 5L)
  check_static_source_sets(result, graph)
})

test_that("aliases fold by canonical member identity despite different printed quantities", {
  graph <- static_identity_graph("left_t0 -> right_t1; right_t1 -> T; right_t1 -> Y; T -> Y",
    c(left_t0 = "P", right_t1 = "P", T = "T", Y = "Y"),
    c(left_t0 = "first label", right_t1 = "second label", T = "T", Y = "Y"))
  summary <- .gd_summary(graph, keep = c("left_t0", "right_t1"))
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("left", "T", "Y"))
  map <- attr(summary$graph, "cf_node_map")
  expect_setequal(map$occurrence[map$node == "left"], c("left_t0", "right_t1"))
  expect_setequal(map$quantity[map$node == "left"], c("first label", "second label"))
  expect_identical(.cd_variable_count(summary$graph), 3L)
  reordered <- graph
  attr(reordered, "cf_node_map") <- attr(graph, "cf_node_map")[rev(seq_len(nrow(map))), ]
  expect_identical(.gd_signature(.gd_summary(reordered, keep = c("left_t0", "right_t1"))$graph),
    .gd_signature(summary$graph))
  check_static_source_sets(summary$graph, graph)
  built <- suppressWarnings(construct_dag(graph, keep_nodes = c("left_t0", "right_t1"), quiet = TRUE))
  actual <- .gd_signature(built)
  expected <- .gd_signature(summary$graph)
  actual$labels <- expected$labels <- NULL
  expect_identical(actual, expected)
  expect_identical(dag_options(built)$level$variables, 3L)
  expect_setequal(attr(built, "cf_node_map")$occurrence, names(graph))
})

test_that("recombined nodes fold only with the same complete member set", {
  graph <- static_identity_graph("AB_t0 -> AB_t1; AB_t1 -> T; AB_t1 -> Y; A_t0 -> T; A_t0 -> Y; T -> Y",
    c(AB_t0 = "A", AB_t1 = "A", A_t0 = "A", T = "T", Y = "Y"),
    c(AB_t0 = "combined", AB_t1 = "combined", A_t0 = "combined", T = "T", Y = "Y"))
  map <- attr(graph, "cf_node_map")
  second <- map[map$node %in% c("AB_t0", "AB_t1"), ]
  second$member_id <- "B"
  attr(graph, "cf_node_map") <- rbind(map, second)
  summary <- .gd_summary(graph, keep = c("AB_t0", "AB_t1", "A_t0"))
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("AB", "A", "T", "Y"))
  shown <- attr(summary$graph, "cf_node_map")
  expect_setequal(shown$member_id[shown$node == "AB"], c("A", "B"))
  expect_setequal(shown$occurrence[shown$node == "AB"], c("AB_t0", "AB_t1"))
  expect_identical(nrow(shown[shown$node == "AB", ]), 4L)
  expect_identical(.cd_variable_count(summary$graph), 4L)
  expect_true(all(summary$edge_ledger$from[summary$edge_ledger$v == "A_t0"] == "A"))
  expect_true(all(summary$edge_ledger$from[summary$edge_ledger$v == "AB_t1"] == "AB"))
  check_static_source_sets(summary$graph, graph)
})

test_that("endpoint comparisons use exact identities and never a colliding display stem", {
  graph <- dagitty::dagitty("dag { X_t0 [exposure]; X_t1 [exposure]; Y [outcome]; X_t0 -> A -> Y; X_t1 -> Y }")
  nodes <- names(graph)
  attr(graph, "cf_node_map") <- data.frame(node = nodes, quantity = nodes,
    member_id = nodes, occurrence = nodes)
  summary <- .gd_summary(graph)
  expect_true(summary$applied)
  expect_setequal(dagitty::exposures(summary$graph), c("X_t0", "X_t1"))
  expect_identical(.gd_endpoint_paths(graph), .gd_endpoint_paths(summary$graph))
  expect_length(.gd_endpoint_paths(summary$graph), 2L)
  expect_true(all(.gd_endpoint_paths(summary$graph)))
  check_static_source_sets(summary$graph, graph)
})

test_that("a folded alias cannot take another variable's exact node name", {
  graph <- static_identity_graph("A_t0 -> A_t1; A_t1 -> T; A_t1 -> Y; A -> T; A -> Y; T -> Y",
    c(A_t0 = "first", A_t1 = "first", A = "second", T = "T", Y = "Y"))
  summary <- .gd_summary(graph, keep = c("A_t0", "A_t1", "A"))
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("A_t0", "A", "T", "Y"))
  map <- attr(summary$graph, "cf_node_map")
  expect_setequal(map$member_id[map$node == "A_t0"], "first")
  expect_setequal(map$member_id[map$node == "A"], "second")
  check_static_source_sets(summary$graph, graph)
})

test_that("colliding display stems do not merge different recorded identities", {
  graph <- static_identity_graph("P_t0 -> T; P_t0 -> Y; P_t1 -> T; P_t1 -> Y; T -> Y",
    c(P_t0 = "first", P_t1 = "second", T = "T", Y = "Y"),
    c(P_t0 = "printed", P_t1 = "printed", T = "T", Y = "Y"))
  summary <- .gd_summary(graph, keep = c("P_t0", "P_t1"))
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("P_t0", "P_t1", "T", "Y"))
  map <- attr(summary$graph, "cf_node_map")
  expect_identical(map$member_id[map$node == "P_t0"], "first")
  expect_identical(map$member_id[map$node == "P_t1"], "second")
  check_static_source_sets(summary$graph, graph)
})

test_that("saved selections refuse older or altered construction versions", {
  graph <- static_identity_graph("T -> Y")
  dag <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  saved <- dag_call(dag)
  expect_identical(saved$candidate$binding$version, .cd_construction_version)
  expect_identical(.gd_signature(suppressWarnings(eval(saved))), .gd_signature(dag))
  for (version in c("variables-v1-root-variable-v1-explicit-helpers-v1", "changed")) {
    replay <- saved
    replay$candidate$binding$version <- version
    expect_error(eval(replay), "unsupported construction version")
  }
})
