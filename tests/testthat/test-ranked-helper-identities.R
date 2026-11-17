helper_query <- function(graph) {
  nodes <- names(graph)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = NULL, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = nodes, quantity = nodes, member_id = nodes, occurrence = nodes),
    policy = data.frame(node_id = nodes, recovered_node = "yes"))
  graph
}

helper_minimal <- function(graph, sets, effect) {
  validate <- .gd_validator(graph)
  vapply(sets, function(set) {
    isTRUE(validate(list(set), effect)) &&
      (!length(set) || !any(validate(lapply(seq_along(set), function(i) set[-i]), effect)))
  }, logical(1))
}

helper_subsets <- function(nodes) {
  c(list(character()), unlist(lapply(seq_along(nodes), function(k) {
    combn(nodes, k, simplify = FALSE)
  }), recursive = FALSE))
}

test_that("real names and inserted parents have distinct exact identities", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y;
    .bidirected_1 -> C; .bidirected_1 -> Y; C <-> T;
    L1 -> T; L1 -> Y; C <-> Y; U [latent]; U -> C }")
  attrs <- attributes(graph)
  indexed <- .cd_ranked_graph(graph)
  expect_identical(attributes(graph), attrs)
  expect_identical(indexed$nodes[seq_along(names(graph))], names(graph))
  expect_identical(indexed$id, stats::setNames(seq_along(indexed$nodes), indexed$nodes))
  expect_false(anyDuplicated(indexed$nodes) > 0L)
  expect_length(indexed$synthetic, 2L)
  expect_length(intersect(indexed$nodes[indexed$synthetic], names(graph)), 0L)
  expect_identical(indexed$latent, c(unname(indexed$id["U"]), indexed$synthetic))
  expect_false(unname(indexed$id["L1"]) %in% indexed$synthetic)
  expect_setequal(lapply(indexed$synthetic, function(i) sort(indexed$nodes[indexed$children[[i]]])),
    list(c("C", "T"), c("C", "Y")))
  expect_identical(.cd_ranked_graph(graph), indexed)
  expect_identical(.cd_ranked_without(indexed, rep(TRUE, length(indexed$from)))$synthetic, indexed$synthetic)
  expect_identical(unserialize(serialize(indexed, NULL)), indexed)
})

test_that("source ordering survives canonicalization and explicit helper indices", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A <-> T; B <-> Y; A -> Y; T -> Y }")
  nodes <- c("Y", "B", "T", "A")
  attr(graph, "cf_source_order") <- list(version = 1L,
    ranks = data.frame(node_id = nodes, base_order = seq_along(nodes), at_order = seq_along(nodes),
      post_order = seq_along(nodes), pre_order = seq_along(nodes), time_order = seq_along(nodes)))
  indexed <- .cd_ranked_graph(graph)
  expect_identical(indexed$nodes[1:4], nodes)
  expect_identical(indexed$synthetic, 5:6)
  expect_identical(indexed$from, c(3L, 4L, 5L, 6L, 5L, 6L))
  expect_identical(indexed$to, c(1L, 1L, 1L, 3L, 2L, 4L))
})

test_that("minimality is invariant under legal real-node renaming for both effects", {
  for (name in c("A", ".bidirected_1", ".bidirected_other", "L1")) {
    graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome]; T -> Y;",
      name, "-> C;", name, "-> Y; C <-> T }"))
    sets <- helper_subsets(c(name, "C"))
    for (effect in c("total", "direct")) {
      expect_identical(.cd_maximal_minimal_sets(graph, sets, effect, .gd_validator(graph)),
        helper_minimal(graph, sets, effect), info = paste(name, effect))
      expect_identical(.cd_maximal_minimal_sets(graph, list(character()), effect, .gd_validator(graph)), TRUE,
        info = paste(name, effect))
    }
  }
})

test_that("prefix-only real members remain adjustable in both searches", {
  for (name in c("A", ".bidirected_1", ".bidirected_other", "L1")) {
    graph <- helper_query(dagitty::dagitty(paste("dag { T [exposure]; Y [outcome]; T -> Y;",
      name, "-> T;", name, "-> Y }")))
    prepared <- .cd_prepare(graph, NULL, NULL)
    for (effect in c("total", "direct")) {
      result <- .cd_ranked_search(prepared, effect)
      expect_identical(result$sets, list(name), info = paste(name, effect))
      expect_identical(result$ranked$certificate, "complete")
      expect_identical(.cd_maximal_setup(prepared, effect)$vocabulary, name)
    }
  }
})

test_that("multiple helpers preserve exact minimal sets and ordinary latent nodes", {
  graph <- dagitty::dagitty("dag { T1 [exposure]; T2 [exposure]; Y [outcome];
    T1 -> M -> Y; T2 -> Y; .bidirected_1 -> T1; .bidirected_1 -> Y;
    L1 -> T2; L1 <-> Y; M <-> Y; U [latent]; U -> L1 }")
  sets <- helper_subsets(c(".bidirected_1", "L1", "M"))
  for (effect in c("total", "direct")) {
    expect_identical(.cd_maximal_minimal_sets(graph, sets, effect, .gd_validator(graph)),
      helper_minimal(graph, sets, effect))
    region <- if (effect == "total") .cd_ranked_total_region(graph) else .cd_ranked_direct_region(graph)
    expect_length(region$graph$synthetic, 2L)
    expect_length(intersect(region$removable, region$graph$synthetic), 0L)
    expect_false(unname(region$graph$id["U"]) %in% region$graph$synthetic)
  }
})

test_that("bidirected edges introduced by construction do not change valid empty sets", {
  graph <- helper_query(dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y;
    .bidirected_1 -> C; .bidirected_1 -> Y; H -> C; H -> T }"))
  expect_no_warning(dag <- construct_dag(graph, 0, keep_nodes = c(".bidirected_1", "C"), quiet = TRUE))
  expect_setequal(names(dag), c(".bidirected_1", "C", "T", "Y"))
  expect_true(any(dagitty::edges(dag)$e == "<->"))
  for (effect in c("total", "direct")) {
    expect_true(dag_checks(dag)[[effect]]$sets_minimal)
    expect_identical(dag_checks(dag)[[effect]]$status, "preserved")
    exact <- dagitty::adjustmentSets(dag, effect = effect)
    expect_length(exact, 1L)
    expect_length(exact[[1L]], 0L)
  }
})

test_that("ranked labels use exact mapped identities and never infer time from a name", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; P_t0 -> T; P_t1 -> Y }")
  map <- data.frame(node = c("P_t0", "P_t1"), quantity = c("printed", "printed"),
    member_id = c("source0", "source1"))
  attr(graph, "input_labels") <- c(P = "printed")
  prepared <- list(selection = graph, map = map)
  expect_identical(.cd_ranked_labels(prepared, c("P_t0", "P_t1", "Q_t0", "Q_t1")),
    list("source0", "source1", "Q_t0", "Q_t1"))
  prepared$map <- rbind(map, data.frame(node = "P_t0", quantity = "printed", member_id = "source2"))
  expect_identical(.cd_ranked_labels(prepared, "P_t0"), list(c("source0", "source2")))
  prepared$map$member_id <- "P"
  expect_identical(.cd_ranked_labels(prepared, c("P_t0", "P_t1")), list("P", "P"))
  prepared <- .cd_prepare(helper_query(dagitty::dagitty(
    "dag { T [exposure]; Y [outcome]; A -> T; A -> Y; T -> Y }")), NULL, NULL)
  another <- prepared$map[prepared$map$node == "A", , drop = FALSE]
  another$member_id <- "another_source"
  prepared$map <- rbind(prepared$map, another)
  for (effect in c("total", "direct")) {
    expect_error(.cd_maximal_setup(prepared, effect), class = "cd_several_members")
  }
})
