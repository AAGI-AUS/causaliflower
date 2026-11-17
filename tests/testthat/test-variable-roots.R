root_query <- function(edges, identities = c(P_t0 = "P", P_t1 = "P"), treatments = "T") {
  graph <- dagitty::dagitty(paste("dag {", paste(treatments, "[exposure];", collapse = " "),
    "Y [outcome];", edges, "}"))
  nodes <- names(graph)
  members <- nodes
  known <- nodes %in% names(identities)
  members[known] <- unname(identities[nodes[known]])
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = nodes, quantity = members, member_id = members, occurrence = nodes),
    policy = data.frame(node_id = unique(members), recovered_node = "yes"))
  graph
}

test_that("a qualifying occasion retains its variable within the prepared query", {
  graph <- root_query("P_t0 -> P_t1 -> K -> J -> T -> Y; Q -> P_t1")
  before <- as.character(graph)
  rules <- .cd_rule_nodes(.cd_prepare(graph, NULL, NULL)$analysis)
  expect_setequal(rules[[5L]], "J")
  expect_setequal(rules[[6L]], c("P_t0", "P_t1", "Q", "J"))
  expect_identical(as.character(graph), before)
})

test_that("same-variable temporal parents do not veto a qualifying occasion", {
  graph <- root_query("Q -> P_t0 -> P_t1 -> K -> J -> T -> Y")
  rules <- .cd_rule_nodes(.cd_prepare(graph, NULL, NULL)$analysis)
  expect_setequal(rules[[6L]], c("P_t0", "P_t1", "Q", "J"))
  expect_true(any(dag_edges(graph)$v == "P_t0" & dag_edges(graph)$w == "P_t1"))
})

test_that("external parents on every relevant occasion prevent root qualification", {
  graph <- root_query("Q0 -> P_t0 -> P_t1 -> K -> J -> T -> Y; Q1 -> P_t1")
  rules <- .cd_rule_nodes(.cd_prepare(graph, NULL, NULL)$analysis)
  expect_false(any(c("P_t0", "P_t1") %in% rules[[6L]]))
  expect_true(all(c("Q0", "Q1", "J") %in% rules[[6L]]))
})

test_that("an irrelevant parentless occasion does not qualify a variable", {
  graph <- root_query("Q0 -> P_t0 -> P_t1 -> K -> J -> T -> Y; Q1 -> P_t1; P_t2 -> R",
    c(P_t0 = "P", P_t1 = "P", P_t2 = "P"))
  rules <- .cd_rule_nodes(.cd_prepare(graph, NULL, NULL)$analysis)
  expect_false(any(c("P_t0", "P_t1", "P_t2", "R") %in% rules[[6L]]))
})

test_that("root qualification uses recorded identities rather than spelling", {
  graph <- root_query("Q -> early -> late -> K -> J -> T -> Y", c(early = "P", late = "P"))
  rules <- .cd_rule_nodes(.cd_prepare(graph, NULL, NULL)$analysis)
  expect_true(all(c("early", "late") %in% rules[[6L]]))
  plain <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; Q -> P_t0 -> P_t1 -> K -> J -> T -> Y }")
  expect_false(any(c("P_t0", "P_t1") %in% .cd_rule_nodes(plain)[[6L]]))
})

test_that("root eligibility follows source membership when display quantities differ", {
  graph <- root_query("Q -> early -> late -> K -> J -> T -> Y", c(early = "P", late = "P"))
  record <- attr(graph, "query_record")
  record$node_map$quantity[match(c("early", "late"), record$node_map$node)] <- c("First label", "Second label")
  attr(graph, "query_record") <- record
  rules <- .cd_rule_nodes(.cd_prepare(graph, NULL, NULL)$analysis)
  expect_true(all(c("early", "late") %in% rules[[6L]]))
})

test_that("root retention preserves the policy boundary and included parent edges", {
  graph <- root_query("P_t0 -> P_t1 -> K -> J -> T -> Y; Q -> P_t1; P_t1 -> P_t2",
    c(P_t0 = "P", P_t1 = "P", P_t2 = "outside"))
  record <- attr(graph, "query_record")
  record$policy_applied <- TRUE
  record$policy$recovered_node[record$policy$node_id == "outside"] <- "no:fixture"
  attr(graph, "query_record") <- record
  dag <- construct_dag(graph, unroll = TRUE, quiet = TRUE)
  temporal <- temporal_dag(dag)
  expect_true(all(c("P_t0", "P_t1", "Q") %in% names(temporal)))
  expect_false("P_t2" %in% names(temporal))
  edges <- dag_edges(temporal)
  expect_true(any(edges$v == "Q" & edges$w == "P_t1" & edges$e == "->"))
  expect_true(any(edges$v == "P_t0" & edges$w == "P_t1" & edges$e == "->"))
  expect_identical(as.character(original_query(dag)), as.character(graph))
})

test_that("identification attribution excludes parents already kept by a size rule", {
  graph <- root_query("X1 -> L -> X2 -> M -> Y; L -> M", character(), c("X1", "X2"))
  dag <- construct_dag(graph, unroll = TRUE, quiet = TRUE)
  identification <- dag_options(dag)$identification
  expect_true(identification$applied)
  expect_true("L" %in% identification$parents)
  expect_false("L" %in% identification$added_nodes)
  expect_true("L" %in% names(temporal_dag(dag)))
})

test_that("identification attribution excludes direct-path preservation", {
  graph <- root_query("X1 -> L -> X2 -> M -> Y; L -> M", character(), c("X1", "X2"))
  observed <- .cd_prepare(graph, NULL, NULL)
  empty <- list(state = "no_sets", members = character(), undecided = character())
  pins <- .cd_level_pins(observed, list(direct = empty, total = empty))
  expect_true("L" %in% pins$mediating)
  expect_true("L" %in% pins$identification$parents)
  expect_false("L" %in% pins$identification$added_nodes)
  expect_true("L" %in% pins$nodes)
})

test_that("only otherwise unretained treatment parents are identification additions", {
  graph <- root_query("U -> X1; V -> X1; X1 -> X2; X1 -> L -> X2 -> M -> Y; L -> M; L -> Y",
    character(), c("X1", "X2"))
  direct <- dagitty::adjustmentSets(graph, effect = "direct")
  total <- dagitty::adjustmentSets(graph, effect = "total")
  expect_length(direct, 1L)
  expect_setequal(direct[[1L]], c("L", "M"))
  expect_length(total, 0L)
  observed <- .cd_prepare(graph, NULL, "V")
  found <- list(direct = list(state = "sets", members = as.character(direct[[1L]]), undecided = character()),
    total = list(state = "no_sets", members = character(), undecided = character()))
  pins <- .cd_level_pins(observed, found)
  expect_setequal(pins$identification$parents, c("L", "U", "V", "X1"))
  expect_identical(pins$identification$added_nodes, "U")
  expect_true(all(pins$identification$parents %in% pins$nodes))
  for (reverse in c(FALSE, TRUE)) {
    .cf_cache$clear()
    build <- function(full) if (full) construct_dag(graph, keep_nodes = "V", unroll = TRUE, quiet = TRUE) else
      .cd_sized(graph, NULL, "V", FALSE, TRUE, 0, unroll = TRUE, quick = TRUE)
    results <- lapply(if (reverse) c(TRUE, FALSE) else c(FALSE, TRUE), build)
    if (reverse) results <- rev(results)
    expect_identical(.gd_signature(temporal_dag(results[[1L]])), .gd_signature(temporal_dag(results[[2L]])))
    expect_identical(dag_options(results[[1L]])$identification$added_nodes, "U")
    expect_false("U" %in% dag_options(results[[2L]])$identification$added_nodes)
  }
})
