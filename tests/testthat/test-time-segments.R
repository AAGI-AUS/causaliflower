# Author's rule for the static view: fold all occasions. The
# quantities on a common cycle of the fold are split symmetrically into the
# fewest contiguous segments that restore acyclicity (fewest display edges,
# then earliest cut). Analysis stays on the graph unrolled over time.

as_segment_query <- function(graph, columns = NULL) {
  nodes <- names(graph)
  ids <- .strip_slice(nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = columns, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = if (is.null(columns)) list() else stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = ids, member_id = ids, occurrence = nodes),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  class(graph) <- c("dag_by_query", class(graph))
  graph
}

segment_query <- function(edges, columns = NULL)
  as_segment_query(dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}")), columns)

segment_graph <- function(edges) dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))

# Every adjustment set a display reports, read as the occurrences its nodes
# stand for, validated on a reference graph.
display_sets_valid <- function(result, reference) {
  state <- .gd_state(result)
  map <- state$summary$node_map
  for (effect in c("direct", "total")) for (set in .gd_adjustment_sets(.gd_plain(result), effect, Inf)$sets) {
    reading <- if (isTRUE(state$summary$applied) && "display" %in% names(map))
      as.character(unlist(lapply(set, function(node) map$node[map$display == node]))) else as.character(set)
    if (!isTRUE(.gd_validate_sets(reference, list(reading), effect))) return(FALSE)
  }
  TRUE
}

# Lagged feedback: A -> B at lags 1 (t0 -> t1, t1 -> t2) and B -> A at lag 2.
lagged_feedback <- paste("T -> A_t0; A_t0 -> B_t1; A_t1 -> B_t2; B_t0 -> A_t2; A_t2 -> Y; B_t2 -> Y; T -> Y;",
                         "A_t0 -> A_t1; A_t1 -> A_t2; B_t0 -> B_t1; B_t1 -> B_t2")
# The same feedback as a confounder of T and Y rather than a mediator.
confounding_feedback <- paste("A_t2 -> T; A_t2 -> Y; B_t2 -> Y; T -> Y; A_t0 -> B_t1; A_t1 -> B_t2; B_t0 -> A_t2;",
                              "A_t0 -> A_t1; A_t1 -> A_t2; B_t0 -> B_t1; B_t1 -> B_t2")

test_that("feedback groups are the strongly connected quantities of the fold", {
  expect_identical(.gd_feedback_groups(c("A", "B", "A", "T"), c("B", "A", "Y", "A")), list(c("A", "B")))
  expect_identical(.gd_feedback_groups(c("Z", "Y", "C", "D"), c("Y", "Z", "D", "C")), list(c("C", "D"), c("Y", "Z")))
  expect_identical(.gd_feedback_groups(c("A", "B"), c("B", "C")), list())
  expect_identical(.gd_feedback_groups(character(), character()), list())
  expect_identical(.gd_cut_candidates(c(0, 1, 2, 3), 1L), list(numeric()))
  expect_identical(.gd_cut_candidates(c(0, 1, 2, 3), 3L), list(c(0, 1), c(0, 2), c(1, 2)))
  expect_identical(.gd_cut_candidates(c(0, 1), 3L), list())
  expect_identical(.gd_cut_candidates(5, 2L), list())
})

test_that("a cyclic fold is split into the fewest segments, with one cut shared by the group", {
  graph <- segment_graph(lagged_feedback)
  summary <- .gd_summary(graph, keep = "B_t0")
  expect_true(summary$applied)
  expect_null(summary$reason)
  display <- summary$graph
  expect_setequal(names(display), c("T", "Y", "A_t0", "A_t1to2", "B_t0", "B_t1to2"))
  expect_true(dag_is_acyclic(display))
  # Folding the retained graph to one node per quantity gives a cycle: A -> B -> A
  shown <- summary$edge_ledger[summary$edge_ledger$e == "->" & summary$edge_ledger$from != summary$edge_ledger$to, ]
  expect_identical(.gd_feedback_groups(shown$from, shown$to), list(c("A", "B")))
  expect_identical(summary$segments, data.frame(
    quantity = c("A", "A", "B", "B"), segment = c("A_t0", "A_t1to2", "B_t0", "B_t1to2"),
    first = c(0, 1, 0, 1), last = c(0, 2, 0, 2),
    members = c("A_t0", "A_t1, A_t2", "B_t0", "B_t1, B_t2"), stringsAsFactors = FALSE))
  map <- summary$node_map
  expect_identical(map$display[match(c("A_t0", "A_t1", "A_t2", "B_t0", "B_t1", "B_t2", "T", "Y"), map$node)],
                   c("A_t0", "A_t1to2", "A_t1to2", "B_t0", "B_t1to2", "B_t1to2", "T", "Y"))
  expect_identical(map$quantity, .strip_slice(map$node))
  # Edges within a quantity join different segments of the same quantity only.
  edges <- with(dag_edges(display), paste(v, e, w))
  expect_true(all(c("A_t0 -> A_t1to2", "B_t0 -> B_t1to2", "A_t1to2 -> B_t1to2", "B_t0 -> A_t1to2") %in% edges))
  expect_false(any(c("A_t1to2 -> A_t1to2", "A_t0 -> A_t0") %in% edges))
  # The display node map has one row per represented temporal occurrence and
  # reads every segment as its quantity.
  display_map <- attr(display, "cf_node_map")
  expect_identical(nrow(display_map), length(names(summary$temporal)))
  expect_setequal(display_map$occurrence[display_map$node == "A_t1to2"], c("A_t1", "A_t2"))
  expect_identical(unname(.gd_node_ids(display)[c("A_t0", "A_t1to2", "B_t1to2", "T")]), c("A", "A", "B", "T"))
  # Analysis stays on the graph unrolled over time.
  expect_identical(.gd_signature(summary$temporal), .gd_signature(.gd_canonical(graph)))
  # Earliest cut on the tie: {t0} | {t1, t2} and {t0, t1} | {t2} both need 9 edges.
  cuts <- .gd_segment_cuts(summary$edge_ledger, map[c("node", "quantity", "node_id")], list(c("A", "B")))
  expect_identical(cuts, list(A = 0, B = 0))
})

test_that("every displayed set remains valid on the source and segments map by occurrence", {
  graph <- segment_query(lagged_feedback, columns = c("A", "B"))
  reference <- .cd_prepare(graph, NULL, NULL)$selection
  summary <- .gd_summary(.cd_prepare(graph, NULL, NULL)$analysis, keep = "B_t0",
                         source_check = .cd_summary_source_check(reference))
  expect_true(summary$applied)
  display <- summary$graph
  for (effect in c("direct", "total")) {
    sets <- .gd_adjustment_sets(display, effect, Inf)$sets
    expanded <- lapply(sets, function(set) as.character(unlist(lapply(set, function(node)
      summary$node_map$node[summary$node_map$display == node]))))
    expect_true(all(.gd_validate_sets(reference, expanded, effect)), info = effect)
  }
  expect_identical(.gd_map_set("A_t1", reference, display), "A_t1to2")
  expect_identical(.gd_map_set("A_t2", reference, display), "A_t1to2")
  expect_identical(.gd_map_set("A_t0", reference, display), "A_t0")
  expect_identical(.gd_map_set(c("A_t1", "B_t2"), reference, display), c("A_t1to2", "B_t1to2"))
  # An occurrence that no segment represents is not mapped.
  expect_null(.gd_map_set("A_t3", dagitty::dagitty("dag { A_t3 }"), display))
  # A segment spanning several occasions has no single source reading (requires review).
  expect_null(.gd_map_set("A_t1to2", display, reference))
  expect_identical(.gd_map_set("A_t0", display, reference), "A_t0")
})

test_that("construction and views keep the segments and their checks", {
  graph <- segment_query(confounding_feedback, columns = c("A", "B"))
  before <- serialize(graph, NULL)
  result <- construct_dag(graph, quiet = TRUE)
  state <- .gd_state(result)
  expect_true(state$summary$applied)
  expect_setequal(names(result), c("T", "Y", "A_t0", "A_t1to2", "B_t0", "B_t1to2"))
  expect_identical(state$summary$segments$segment, c("A_t0", "A_t1to2", "B_t0", "B_t1to2"))
  expect_setequal(names(temporal_dag(result)), c("T", "Y", paste0(rep(c("A", "B"), each = 3), "_t", 0:2)))
  checks <- dag_checks(result)
  expect_identical(checks$total$status, "preserved")
  expect_identical(checks$direct$status, "preserved")
  expect_identical(checks$total$required, list("A_t2"))
  expect_identical(checks$total$final_sets, list("A_t1to2"))
  reference <- .cd_prepare(graph, NULL, NULL)$selection
  expect_true(.gd_validate_sets(reference, list(c("A_t1", "A_t2")), "total"))
  view <- dag_sets_view(result, list(total = 1L, direct = 1L))
  view_state <- .gd_state(view)
  expect_true(view_state$summary$applied)
  expect_identical(view_state$summary$segments, state$summary$segments)
  expect_setequal(names(view), names(result))
  expect_identical(vapply(dag_checks(view), `[[`, character(1), "status"),
                   c(total = "preserved", direct = "preserved"))
  expect_identical(.gd_signature(temporal_dag(view)), .gd_signature(temporal_dag(result)))
  expect_identical(serialize(original_query(result), NULL), before)
})

test_that("a segmented alternative without a single source reading is flagged, never certified", {
  # The mediating feedback adds the display set {A_t1to2, B_t1to2}: it is
  # valid on the source when read as its occurrences, but it is not a source
  # option and has no single reading at a specific time, so review is required.
  graph <- segment_query(lagged_feedback, columns = c("A", "B"))
  warned <- character()
  result <- withCallingHandlers(construct_dag(graph, quiet = TRUE),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  state <- .gd_state(result)
  expect_true(state$summary$applied)
  expect_identical(state$summary$segments$segment, c("A_t0", "A_t1to2", "B_t0", "B_t1to2"))
  direct <- dag_checks(result)$direct
  expect_identical(direct$status, "changed")
  introduced <- Filter(function(x) setequal(x$set, c("A_t1to2", "B_t1to2")), direct$introduced)
  expect_length(introduced, 1L)
  expect_identical(introduced[[1L]]$status, "temporal interpretation requires review")
  expect_identical(dag_checks(result)$total$status, "preserved")
  expect_true(any(grepl("changed or were not fully verified", warned)))
})

test_that("three segments are used when no single cut restores acyclicity", {
  chain <- paste("T -> A_t0; A_t0 -> B_t1; B_t1 -> A_t2; A_t2 -> B_t3; B_t3 -> A_t4; A_t4 -> B_t5; B_t5 -> Y; T -> Y;",
                 "A_t0 -> A_t1; A_t1 -> A_t2; A_t2 -> A_t3; A_t3 -> A_t4; A_t4 -> A_t5;",
                 "B_t0 -> B_t1; B_t1 -> B_t2; B_t2 -> B_t3; B_t3 -> B_t4; B_t4 -> B_t5")
  summary <- .gd_summary(segment_graph(chain))
  expect_true(summary$applied)
  mapping <- summary$node_map[c("node", "quantity", "node_id")]
  links <- summary$edge_ledger
  # Every single cut leaves an A -> B -> A chain inside one segment.
  for (cut in 0:4) {
    display <- stats::setNames(.gd_segment_display(mapping, list(A = cut, B = cut)), mapping$node)
    expect_false(.gd_edges_acyclic(.gd_display_edges(links, display)), info = cut)
  }
  expect_identical(.gd_segment_cuts(links, mapping, list(c("A", "B"))), list(A = c(1, 3), B = c(1, 3)))
  expect_identical(summary$segments$segment, c("A_t0to1", "A_t2to3", "A_t4", "B_t1", "B_t2to3", "B_t4to5"))
  expect_true(dag_is_acyclic(summary$graph))
})

test_that("fewer display edges outrank an earlier cut, and a member with one segment keeps its name", {
  graph <- segment_graph("T -> A_t0; A_t0 -> A_t1; A_t1 -> A_t2; A_t0 -> B_t1; B_t1 -> A_t2; A_t2 -> Y; T -> Y; C -> A_t0; C -> A_t1")
  summary <- .gd_summary(graph, keep = "C")
  expect_true(summary$applied)
  # {t0} | {t1, t2} needs 8 display edges. {t0, t1} | {t2} needs 7.
  expect_setequal(names(summary$graph), c("T", "Y", "C", "B", "A_t0to1", "A_t2"))
  expect_identical(summary$segments$segment, c("A_t0to1", "A_t2"))
  expect_identical(nrow(dag_edges(summary$graph)), 7L)
  expect_identical(summary$node_map$display[summary$node_map$node == "B_t1"], "B")
})

test_that("several feedback groups are cut together, not one after another by name", {
  pair <- function(a, b) paste0(a, "_t0 -> ", b, "_t1; ", a, "_t1 -> ", b, "_t2; ", b, "_t0 -> ", a, "_t2; ",
    a, "_t0 -> ", a, "_t1; ", a, "_t1 -> ", a, "_t2; ", b, "_t0 -> ", b, "_t1; ", b, "_t1 -> ", b, "_t2")
  two_groups <- function(a, b, c, d) segment_graph(paste0("T -> ", a, "_t0; ", pair(a, b), "; ", pair(c, d), "; ",
    a, "_t0 -> ", c, "_t0; ", a, "_t1 -> ", c, "_t0; ", a, "_t2 -> ", c, "_t2; ",
    b, "_t2 -> Y; ", c, "_t2 -> Y; ", d, "_t2 -> Y; T -> Y"))
  graph <- two_groups("A", "B", "C", "D")
  summary <- .gd_summary(graph, keep = names(graph))
  expect_true(summary$applied)
  mapping <- summary$node_map[c("node", "quantity", "node_id")]
  # Settling A,B first with C,D still folded would cut both groups at t0 (18
  # display edges). The joint choice cuts A,B at t1 and C,D at t0 (17 edges),
  # the fewest over all four combinations.
  expect_identical(.gd_segment_cuts(summary$edge_ledger, mapping, list(c("A", "B"), c("C", "D"))),
                   list(A = 1, B = 1, C = 0, D = 0))
  expect_identical(nrow(dag_edges(summary$graph)), 17L)
  expect_identical(summary$segment_alternatives, 1L)
  counts <- apply(expand.grid(ab = 0:1, cd = 0:1), 1, function(cut) {
    display <- stats::setNames(.gd_segment_display(mapping,
      list(A = cut[["ab"]], B = cut[["ab"]], C = cut[["cd"]], D = cut[["cd"]])), mapping$node)
    nrow(.gd_display_edges(summary$edge_ledger, display))
  })
  expect_identical(min(counts), 17L)
  expect_identical(unname(counts[2L]), 17L)  # row 2 of the grid is ab = 1, cd = 0
  # Renaming the first group so that it sorts after the second changes nothing.
  renamed <- two_groups("M", "N", "C", "D")
  other <- .gd_summary(renamed, keep = names(renamed))
  expect_identical(sort(chartr("MN", "AB", other$segments$segment)), sort(summary$segments$segment))
  expect_identical(nrow(dag_edges(other$graph)), 17L)
  expect_identical(other$segment_alternatives, 1L)
})

test_that("a residual tie between feedback groups is decided by name, recorded and warned about", {
  # Reviewer fixture (stage-d/review-d-r2-ties.R): both combinations have 18
  # display edges and the cut positions (0, 1). The group that comes first by name takes
  # its earliest cut.
  pair <- function(a, b) paste0(a, "_t0 -> ", b, "_t1; ", a, "_t1 -> ", b, "_t2; ", b, "_t0 -> ", a, "_t2; ",
    a, "_t0 -> ", a, "_t1; ", a, "_t1 -> ", a, "_t2; ", b, "_t0 -> ", b, "_t1; ", b, "_t1 -> ", b, "_t2")
  tied <- function(a, b, c, d) segment_graph(paste0(pair(a, b), "; ", pair(c, d), "; ",
    b, "_t1 -> ", c, "_t1; ", a, "_t0 -> ", c, "_t1; ", b, "_t0 -> ", c, "_t2; ", b, "_t1 -> ", d, "_t2; ",
    a, "_t0 -> ", c, "_t0; T -> ", a, "_t0; ", b, "_t2 -> Y; ", d, "_t2 -> Y; T -> Y"))
  graph <- tied("A", "B", "C", "D")
  expect_silent(summary <- .gd_summary(graph, keep = names(graph)))
  expect_true(summary$applied)
  expect_identical(summary$segment_alternatives, 2L)
  mapping <- summary$node_map[c("node", "quantity", "node_id")]
  expect_identical(.gd_segment_cuts(summary$edge_ledger, mapping, list(c("A", "B"), c("C", "D"))),
                   list(A = 0, B = 0, C = 1, D = 1))
  expect_identical(nrow(dag_edges(summary$graph)), 18L)
  renamed <- tied("M", "N", "C", "D")
  other <- .gd_summary(renamed, keep = names(renamed))
  expect_identical(other$segment_alternatives, 2L)
  expect_identical(.gd_segment_cuts(other$edge_ledger, other$node_map[c("node", "quantity", "node_id")],
                                    list(c("C", "D"), c("M", "N"))),
                   list(C = 0, D = 0, M = 1, N = 1))
  expect_identical(nrow(dag_edges(other$graph)), 18L)
  # The tie is reported by the caller that finally shows the segments, and
  # only then: a discarded summary (check of the set, family checks) is silent.
  expect_warning(.gd_segment_tie_notice(summary), "equally ranked cut choices")
  refused <- .gd_summary(graph, keep = names(graph), source_check = function(display, mapping) "refused by the test")
  expect_false(refused$applied)
  expect_null(refused$segment_alternatives)
  expect_silent(.gd_segment_tie_notice(refused))
  collect <- function(expr) {
    warned <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, tie = any(grepl("equally ranked cut choices", warned)))
  }
  query <- segment_query(paste0(pair("A", "B"), "; ", pair("C", "D"), "; B_t1 -> C_t1; A_t0 -> C_t1; B_t0 -> C_t2;",
    " B_t1 -> D_t2; A_t0 -> C_t0; T -> A_t0; B_t2 -> Y; D_t2 -> Y; T -> Y"), columns = c("A", "B", "C", "D"))
  built <- collect(construct_dag(query, quiet = TRUE))
  state <- .gd_state(built$value)
  expect_identical(built$tie, isTRUE(state$summary$applied) && isTRUE(state$summary$segment_alternatives > 1L))
  for (selection in list(list(total = 1L), list(direct = 1L), list(total = 1L, direct = 1L))) {
    viewed <- collect(dag_sets_view(built$value, selection))
    view_state <- .gd_state(viewed$value)
    expect_identical(viewed$tie, isTRUE(view_state$summary$applied) && isTRUE(view_state$summary$segment_alternatives > 1L))
  }
})

test_that("merged sets are checked on segmented displays and the shown edges across repeated stages keep their sets valid", {
  # Reviewer probe seed 51 (stage-d/review-d-setcheck.R): under Stage D the
  # segmented fold passed every structural guard and only the merged set
  # check refused it (an omitted edge across repeated stages made a displayed set
  # invalid). Under Stage D.1 that edge is shown between its segments, every
  # displayed set is valid when read as its occurrences, and the check
  # accepts the display. A refusing check still returns the temporal graph.
  edges <- paste(
    "Q1_t0 -> Q2_t0; Q1_t0 -> Q3_t1; Q1_t1 -> Q1_t2; Q1_t1 -> Q2_t1; Q1_t1 -> Q4_t2; Q1_t1 -> Q5_t3; Q1_t2 -> Q2_t2;",
    "Q1_t3 -> Q1_t4; Q1_t3 -> Q2_t3; Q1_t3 -> Y; Q1_t4 -> Q2_t4; Q1_t4 -> Y; Q2_t0 -> Q2_t1; Q2_t0 -> Q5_t2;",
    "Q2_t1 -> Q2_t2; Q2_t1 -> Q5_t3; Q2_t2 -> Q1_t3; Q2_t3 -> Q1_t4; Q2_t3 -> Q2_t4; Q2_t3 -> Q3_t4; Q2_t3 -> Q5_t4;",
    "Q2_t4 -> Y; Q3_t0 -> Q1_t0; Q3_t0 -> Q1_t1; Q3_t0 -> Q1_t2; Q3_t0 -> Q2_t0; Q3_t0 -> Q3_t1; Q3_t0 -> Q3_t2;",
    "Q3_t0 -> Q4_t1; Q3_t1 -> Q1_t1; Q3_t1 -> Q1_t3; Q3_t1 -> Q2_t1; Q3_t1 -> Q3_t2; Q3_t1 -> Q4_t2; Q3_t2 -> Q1_t3;",
    "Q3_t2 -> Q2_t2; Q3_t2 -> Q3_t3; Q3_t2 -> Q3_t4; Q3_t3 -> Q1_t3; Q3_t3 -> Q2_t3; Q3_t4 -> Q2_t4; Q4_t0 -> Q1_t2;",
    "Q4_t0 -> Q2_t0; Q4_t0 -> Q4_t1; Q4_t1 -> Q2_t1; Q4_t1 -> Q4_t2; Q4_t1 -> Q4_t3; Q4_t2 -> Q2_t2; Q4_t2 -> Q3_t4;",
    "Q4_t2 -> Q4_t3; Q4_t3 -> Q2_t3; Q4_t3 -> Q4_t4; Q4_t4 -> Q2_t4; Q5_t0 -> Q4_t0; Q5_t0 -> Q5_t1; Q5_t1 -> Q3_t3;",
    "Q5_t1 -> Q4_t1; Q5_t1 -> Q5_t2; Q5_t1 -> Q5_t3; Q5_t2 -> Q4_t2; Q5_t2 -> Q4_t4; Q5_t2 -> Q5_t3; Q5_t3 -> Q1_t4;",
    "Q5_t3 -> Q2_t4; Q5_t3 -> Q4_t3; Q5_t3 -> Q4_t4; Q5_t3 -> Q5_t4; Q5_t4 -> Q4_t4; T -> Q1_t1; T -> Q5_t1; T -> Y", collapse = " ")
  graph <- segment_query(edges, columns = c("Q1", "Q3", "Q5"))
  reference <- .cd_prepare(graph, NULL, NULL)$selection
  open <- .gd_summary(reference)
  expect_true(open$applied)
  expect_false(is.null(open$segments))
  expect_true(any(open$edge_ledger$display_role == "between_segments"))
  invalid <- FALSE
  for (effect in c("direct", "total")) for (set in .gd_adjustment_sets(open$graph, effect, Inf)$sets) {
    reading <- as.character(unlist(lapply(set, function(node) open$node_map$node[open$node_map$display == node])))
    if (!isTRUE(.gd_validate_sets(reference, list(reading), effect))) invalid <- TRUE
  }
  expect_false(invalid)
  checked <- .gd_summary(reference, source_check = .cd_summary_source_check(reference))
  expect_true(checked$applied)
  expect_identical(checked$segments, open$segments)
  refused <- .gd_summary(reference, source_check = function(display, mapping) "refused by the test")
  expect_false(refused$applied)
  expect_setequal(names(refused$graph), names(reference))
})

test_that("a segment is unobserved only when every occurrence it stands for is", {
  graph <- segment_graph(paste(lagged_feedback, "; B_t0 [latent]; B_t1 [latent]; B_t2 [latent]"))
  summary <- .gd_summary(graph, keep = "B_t0")
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("T", "Y", "A_t0", "A_t1to2", "B_t0", "B_t1to2"))
  expect_setequal(dagitty::latents(summary$graph), c("B_t0", "B_t1to2"))
})

test_that("a bidirected edge inside a feedback group is shown between the segments", {
  graph <- segment_graph(paste(lagged_feedback, "; A_t0 <-> B_t2"))
  summary <- .gd_summary(graph, keep = "B_t0")
  expect_true(summary$applied)
  expect_identical(summary$segments$segment, c("A_t0", "A_t1to2", "B_t0", "B_t1to2"))
  edges <- dag_edges(summary$graph)
  expect_true(any(edges$e == "<->" & ((edges$v == "A_t0" & edges$w == "B_t1to2") | (edges$v == "B_t1to2" & edges$w == "A_t0"))))
})

test_that("a segment name that coincides with another quantity keeps the temporal display", {
  graph <- segment_graph(paste(lagged_feedback, "; A_t1to2 -> Y"))
  summary <- .gd_summary(graph, keep = c("B_t0", "A_t1to2"))
  expect_false(summary$applied)
  expect_match(summary$reason, "time segment name coincides")
  expect_setequal(names(summary$graph), names(graph))
})

test_that("a tie decided by name between two eligible edges across repeated stages is counted and reported", {
  # Reviewer fixture (stage-d1/review-d1-conflict.R): A_t0 -> B_t1 and
  # C_t0 -> D_t1 each keep the display acyclic alone but together close a
  # cycle through the unsegmented connectors U and V. Both have lag 1 and
  # source occasion 0, so the edge key decides. The graph is symmetric under
  # (A,B,M,N,U) <-> (C,D,P,R,V).
  tied <- function(a, b, c, d, m, n, p, r, u, v) segment_graph(paste0(
    a, "_t0 -> ", m, "_t1; ", m, "_t0 -> ", a, "_t1; ", b, "_t0 -> ", n, "_t1; ", n, "_t0 -> ", b, "_t1; ",
    c, "_t0 -> ", p, "_t1; ", p, "_t0 -> ", c, "_t1; ", d, "_t0 -> ", r, "_t1; ", r, "_t0 -> ", d, "_t1; ",
    b, "_t0 -> ", a, "_t0; ", b, "_t1 -> ", a, "_t1; ", a, "_t0 -> ", b, "_t1; ",
    d, "_t0 -> ", c, "_t0; ", d, "_t1 -> ", c, "_t1; ", c, "_t0 -> ", d, "_t1; ",
    b, "_t1 -> ", u, "_t1; ", u, "_t0 -> ", c, "_t0; ", d, "_t1 -> ", v, "_t1; ", v, "_t0 -> ", a, "_t0; T -> Y"))
  graph <- tied("A", "B", "C", "D", "M", "N", "P", "R", "U", "V")
  summary <- .gd_summary(graph, keep = names(graph))
  expect_true(summary$applied)
  ledger <- summary$edge_ledger
  roles <- c(ledger$display_role[ledger$v == "A_t0" & ledger$w == "B_t1"], ledger$display_role[ledger$v == "C_t0" & ledger$w == "D_t1"])
  expect_identical(roles, c("between_segments", "between_repetitions"))
  expect_identical(summary$segment_edge_ties, 1L)
  expect_identical(summary$segment_alternatives, 1L)
  expect_warning(.gd_segment_tie_notice(summary), "same lag and occasion came first by name")
  # Renaming so that C sorts first flips the shown edge: the recorded tie is the
  # only trace of that name dependence.
  renamed <- tied("X", "Z", "C", "D", "M", "N", "P", "R", "U", "V")
  other <- .gd_summary(renamed, keep = names(renamed))
  ledger <- other$edge_ledger
  expect_identical(ledger$display_role[ledger$v == "C_t0" & ledger$w == "D_t1"], "between_segments")
  expect_identical(ledger$display_role[ledger$v == "X_t0" & ledger$w == "Z_t1"], "between_repetitions")
  expect_identical(other$segment_edge_ties, 1L)
  # Edges of different tiers never count as a tie: with the C/D pattern shifted
  # one occasion later, C_t1 -> D_t2 (lag 1, occasion 1) is blocked by the
  # earlier tier's A_t0 -> B_t1 (lag 1, occasion 0). The outcome does not depend
  # on names and is not a tie (reviewer fixture, stage-d1/review-d1-r2-ties.R).
  cross <- function(a, b, c, d, m, n, p, r, u, v) segment_graph(paste0(
    a, "_t0 -> ", m, "_t1; ", m, "_t0 -> ", a, "_t1; ", b, "_t0 -> ", n, "_t1; ", n, "_t0 -> ", b, "_t1; ",
    c, "_t1 -> ", p, "_t2; ", p, "_t1 -> ", c, "_t2; ", d, "_t1 -> ", r, "_t2; ", r, "_t1 -> ", d, "_t2; ",
    b, "_t0 -> ", a, "_t0; ", b, "_t1 -> ", a, "_t1; ", a, "_t0 -> ", b, "_t1; ",
    d, "_t1 -> ", c, "_t1; ", d, "_t2 -> ", c, "_t2; ", c, "_t1 -> ", d, "_t2; ",
    b, "_t1 -> ", u, "_t1; ", u, "_t0 -> ", c, "_t1; ", d, "_t2 -> ", v, "_t2; ", v, "_t0 -> ", a, "_t0; T -> Y"))
  for (names in list(c("A", "B", "C", "D", "M", "N", "P", "R", "U", "V"), c("C", "D", "A", "B", "P", "R", "M", "N", "U", "V"))) {
    graph <- do.call(cross, as.list(names))
    summary <- .gd_summary(graph, keep = names(graph))
    expect_true(summary$applied)
    ledger <- summary$edge_ledger
    expect_identical(ledger$display_role[ledger$v == paste0(names[1], "_t0") & ledger$w == paste0(names[2], "_t1")], "between_segments")
    expect_identical(ledger$display_role[ledger$v == paste0(names[3], "_t1") & ledger$w == paste0(names[4], "_t2")], "between_repetitions")
    expect_identical(summary$segment_edge_ties, 0L)
    expect_silent(.gd_segment_tie_notice(summary))
  }
})

test_that("manual edits and the view without selected sets validate their displayed sets on the source", {
  # Reviewer fixture (stage-d1/review-d1-restore.R): a timed treatment on a
  # feedback cycle with a latent quantity. The segmented display would report
  # the empty set although the source has an open latent backdoor path.
  base <- dagitty::dagitty(paste("dag { Q1_t2 -> Q3_t2; Q1_t3 -> Q3_t3; Q2_t0 -> Q1_t2; Q3_t2 -> Q1_t3;",
                                 "Q3_t0 -> Q3_t1; Q3_t1 -> Q3_t2; Q3_t0 -> T_t1; T_t1 -> Q2_t3; Q1_t3 -> Y }"))
  dagitty::exposures(base) <- "T_t1"
  dagitty::outcomes(base) <- "Y"
  latents <- c("Q3_t0", "Q3_t1", "Q3_t2", "Q3_t3")
  query <- as_segment_query(base)
  reference <- .cd_prepare(query, latents, NULL)$selection
  # Without the check of merged sets, the segmented display reports {} for an
  # effect that the source does not identify by {}.
  open <- .gd_summary(reference, keep = names(reference))
  expect_true(open$applied)
  expect_false(is.null(open$segments))
  checked <- .gd_summary(reference, keep = names(reference), source_check = .cd_summary_source_check(reference))
  expect_false(checked$applied)
  expect_match(checked$reason, "represented occurrences")
  result <- suppressWarnings(construct_dag(query, latents = latents, quiet = TRUE))
  expect_true(display_sets_valid(result, reference))
  absent <- setdiff(names(reference), names(temporal_dag(result)))
  restored <- suppressWarnings(restore_nodes(result, absent))
  expect_true(display_sets_valid(restored, reference))
  state <- .gd_state(restored)
  if (isTRUE(state$summary$applied)) expect_false(is.null(state$summary$segments)) else
    expect_match(state$summary$reason, "represented occurrences|from treatment to outcome")
})

test_that("segmentation does not depend on node order or quantity names", {
  statements <- trimws(strsplit(lagged_feedback, ";", fixed = TRUE)[[1L]])
  forward <- .gd_summary(segment_graph(paste(statements, collapse = "; ")), keep = "B_t0")
  backward <- .gd_summary(segment_graph(paste(rev(statements), collapse = "; ")), keep = "B_t0")
  expect_identical(.gd_signature(forward$graph), .gd_signature(backward$graph))
  expect_identical(forward$segments, backward$segments)
  # Swapping the roles of the two quantities (nitrate first or ammonium first)
  # swaps their segments and nothing else: both share one cut.
  swap <- function(x) chartr("AB", "BA", x)
  swapped <- .gd_summary(segment_graph(swap(lagged_feedback)), keep = "A_t0")
  expect_setequal(names(swapped$graph), swap(names(forward$graph)))
  expect_setequal(edge_key(dag_edges(swapped$graph)), swap(edge_key(dag_edges(forward$graph))))
  expect_identical(sort(swapped$segments$segment), sort(swap(forward$segments$segment)))
  nitrate <- gsub("B", "NO3", gsub("A", "NH4", lagged_feedback, fixed = TRUE), fixed = TRUE)
  ammonium <- gsub("A", "NO3", gsub("B", "NH4", lagged_feedback, fixed = TRUE), fixed = TRUE)
  first <- .gd_summary(segment_graph(nitrate), keep = "NO3_t0")
  second <- .gd_summary(segment_graph(ammonium), keep = "NH4_t0")
  expect_identical(first$segments$first, second$segments$first)
  expect_identical(first$segments$last, second$segments$last)
})

test_that("graphs without feedback collapse each quantity to one node as before", {
  summary <- .gd_summary(segment_graph("T -> A_t1 -> Y; T -> A_t2 -> Y"))
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("T", "A", "Y"))
  expect_null(summary$segments)
  expect_identical(summary$node_map$display, summary$node_map$quantity)
  expect_false(any(summary$edge_ledger$display_role == "between_segments"))
  result <- construct_dag(segment_query("T -> A_t1 -> Y; T -> A_t2 -> Y; T -> Y"), quiet = TRUE)
  expect_setequal(names(result), c("T", "A", "Y"))
  expect_null(.gd_state(result)$summary$segments)
})

test_that("a lagged chain with a gap is split at its exact occasions", {
  # Before the segment rule this fold (A -> B -> A) kept the temporal display.
  graph <- segment_graph("T -> A_t0 -> B_t1 -> A_t2 -> Y; T -> Y")
  summary <- .gd_summary(graph, keep = "B")
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("T", "Y", "A_t0", "A_t2", "B"))
  expect_identical(summary$segments$segment, c("A_t0", "A_t2"))
  expect_identical(summary$segments$members, c("A_t0", "A_t2"))
  expect_setequal(with(dag_edges(summary$graph), paste(v, e, w)),
                  c("T -> A_t0", "A_t0 -> B", "B -> A_t2", "A_t2 -> Y", "T -> Y"))
  result <- construct_dag(segment_query("T -> A_t0 -> B_t1 -> A_t2 -> Y; T -> Y"), keep_nodes = "B", quiet = TRUE)
  expect_setequal(names(result), c("T", "Y", "A_t0", "A_t2", "B"))
  expect_identical(vapply(dag_checks(result), `[[`, character(1), "status")[c("direct", "total")],
                   c(direct = "preserved", total = "preserved"))
})

test_that("a segment lies on a route from treatment to outcome only through its own occasions", {
  # T -> A_t1 and A_t2 -> Y are unconnected in the source, yet the segment
  # A_t1to2 would join them into a display route. A_t0 is the only mediating
  # occurrence of A. The feedback with B forces the cut at t0 (cutting at t1
  # leaves A_t0to1 <-> B_t0to1 cyclic). Reading the segment as its quantity
  # would accept this display, because A does lie on a source route.
  graph <- segment_graph(paste("T -> A_t0; A_t0 -> Y; T -> A_t1; A_t2 -> Y; T -> Y;",
                               "A_t0 -> B_t1; B_t0 -> A_t1; A_t0 -> A_t1; B_t0 -> B_t1"))
  summary <- .gd_summary(graph, keep = c("A_t1", "A_t2", "B_t0", "B_t1"))
  expect_false(summary$applied)
  expect_match(summary$reason, "changes which quantities lie on a directed route from treatment to outcome")
  expect_setequal(names(summary$graph), names(graph))
  expect_null(summary$segments)
  # Joining the occurrences in the source makes the segment a genuine mediator.
  joined <- segment_graph(paste("T -> A_t0; A_t0 -> Y; T -> A_t1; A_t1 -> A_t2; A_t2 -> Y; T -> Y;",
                                "A_t0 -> B_t1; B_t0 -> A_t1; A_t0 -> A_t1; B_t0 -> B_t1"))
  summary <- .gd_summary(joined, keep = c("B_t0", "B_t1"))
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("T", "Y", "A_t0", "A_t1to2", "B_t0", "B_t1"))
})

test_that("an endpoint with one occasion on a feedback cycle stays uncut while its partners are segmented", {
  graph <- segment_graph("A_t0 -> T; T -> A_t1; A_t1 -> Y")
  summary <- .gd_summary(graph, keep = "A_t0")
  expect_true(summary$applied)
  expect_setequal(names(summary$graph), c("T", "Y", "A_t0", "A_t1"))
  expect_identical(summary$segments$segment, c("A_t0", "A_t1"))
  expect_setequal(with(dag_edges(summary$graph), paste(v, e, w)), c("A_t0 -> T", "T -> A_t1", "A_t1 -> Y"))
  expect_identical(unname(.gd_node_ids(summary$graph)["T"]), "T")
  expect_identical(dagitty::exposures(summary$graph), "T")
  # In the fixture previously used to check node counts, the treatment sits on
  # the fold's feedback cycle with A. A is split at its exact occasions and the
  # display is kept.
  result <- construct_dag(segment_query("A_t0 -> T; T -> A_t1 -> Y; T -> Y"), keep_nodes = "A", quiet = TRUE)
  state <- .gd_state(result)
  expect_true(state$summary$applied)
  expect_setequal(names(result), c("T", "Y", "A_t0", "A_t1"))
  expect_identical(state$summary$segments$segment, c("A_t0", "A_t1"))
  expect_identical(vapply(dag_checks(result), `[[`, character(1), "status")[c("direct", "total")],
                   c(direct = "preserved", total = "preserved"))
})

test_that("edges across repeated stages are shown between segments when the display stays acyclic", {
  # Reviewer fixture (stage-d/review-d-betweenrep-min.log, seed 70): Q1_t1 -> Q2_t2
  # is an edge across repeated stages (Q2 -> Q1 within both occasions) whose display edge
  # Q1_t0to1 -> Q2_t2to3 keeps the display acyclic.
  graph <- segment_graph(paste("Q2_t1 -> Q1_t1; Q2_t2 -> Q1_t2; Q2_t2 -> Q4_t2; Q2_t3 -> Q4_t3; Q1_t1 -> Q2_t2;",
                               "Q1_t0 -> Q5_t2; Q3_t1 -> Q2_t3; Q4_t2 -> Q3_t3; Q5_t2 -> Q3_t3"))
  summary <- .gd_summary(graph, keep = names(graph))
  expect_true(summary$applied)
  ledger <- summary$edge_ledger
  expect_identical(ledger$display_role[ledger$v == "Q1_t1" & ledger$w == "Q2_t2"], "between_segments")
  expect_true("Q1_t0to1 -> Q2_t2to3" %in% with(dag_edges(summary$graph), paste(v, e, w)))
  expect_true(dag_is_acyclic(summary$graph))
  expect_false(any(ledger$display_role == "between_repetitions"))
  # An edge across repeated stages that would close a cycle stays omitted: soil nitrate
  # feeds uptake within every occasion and uptake depletes the next occasion's
  # nitrate. With nitrate segmented by the ammonium feedback, U -> NO3_t1to2
  # would close NO3_t1to2 -> U -> NO3_t1to2.
  depletion <- segment_graph(paste("T -> NH4_t0; NH4_t0 -> NO3_t1; NH4_t1 -> NO3_t2; NO3_t0 -> NH4_t2; NO3_t2 -> Y; T -> Y;",
                                   "NO3_t0 -> U_t0; NO3_t1 -> U_t1; NO3_t2 -> U_t2; U_t0 -> NO3_t1; U_t1 -> NO3_t2; U_t2 -> Y"))
  summary <- .gd_summary(depletion, keep = names(depletion))
  expect_true(summary$applied)
  expect_identical(summary$segments$segment, c("NH4_t0", "NH4_t1to2", "NO3_t0", "NO3_t1to2"))
  ledger <- summary$edge_ledger
  expect_identical(ledger$display_role[ledger$from == "U" & ledger$to == "NO3"], rep("between_repetitions", 2L))
  edges <- with(dag_edges(summary$graph), paste(v, e, w))
  expect_false(any(grepl("^U -> NO3", edges)))
  expect_true(all(c("NO3_t0 -> U", "NO3_t1to2 -> U") %in% edges))
  # Without segments an edge across repeated stages is never shown again (B stays off every
  # route from treatment to outcome so that the route guard accepts the fold).
  plain <- segment_graph("T -> A_t0; A_t0 -> Y; B_t0 -> A_t0; B_t1 -> A_t1; A_t0 -> B_t1; T -> Y")
  summary <- .gd_summary(plain, keep = c("A_t1", "B_t0", "B_t1"))
  expect_true(summary$applied)
  expect_null(summary$segments)
  expect_setequal(names(summary$graph), c("T", "Y", "A", "B"))
  ledger <- summary$edge_ledger
  expect_identical(ledger$display_role[ledger$v == "A_t0" & ledger$w == "B_t1"], "between_repetitions")
  expect_setequal(with(dag_edges(summary$graph), paste(v, e, w)), c("T -> A", "B -> A", "A -> Y", "T -> Y"))
})
