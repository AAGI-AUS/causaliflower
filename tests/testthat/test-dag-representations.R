# DAG sizes: dagby_construct(size), dagby_prune(), dagby_grow(), the short record they
# keep, and the contract of the exported functions.
# ranked_query is copied from test-ranked-search.R, where it is local, with exposures added.
ranked_query <- function(edges, columns = NULL, exposures = "T") {
  graph <- dagitty::dagitty(paste("dag {", paste(exposures, "[exposure];", collapse = " "), "Y [outcome];", edges, "}"))
  nodes <- names(graph)
  ids <- .strip_slice(nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = columns, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = if (is.null(columns)) list() else stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = ids, member_id = ids, occurrence = nodes),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  graph
}

# The exact reference: the quantity sets of dagitty's minimal sets, proper supersets dropped.
catalogue_sets <- function(graph) {
  lapply(list(direct = "direct", total = "total"), function(effect) {
    raw <- lapply(dagitty::adjustmentSets(graph, type = "minimal", effect = effect), as.character)
    keys <- unique(lapply(raw, function(set) sort(unique(.strip_slice(set)))))
    keep <- !vapply(seq_along(keys), function(i) any(vapply(seq_along(keys), function(j)
      j != i && length(keys[[j]]) < length(keys[[i]]) && all(keys[[j]] %in% keys[[i]]), logical(1))), logical(1))
    keys[keep]
  })
}

# The quantities of those sets.
catalogue_quantities <- function(graph) {
  lapply(catalogue_sets(graph), function(sets) sort(unique(as.character(unlist(sets, use.names = FALSE)))))
}

set_keys <- function(sets) sort(vapply(sets, function(set) paste(sort(set), collapse = "+"), character(1)))

# Run an expression and return its value with the messages and warnings it raised.
# Progress lines are dropped.
heard <- function(expr) {
  messages <- character()
  warnings <- character()
  value <- withCallingHandlers(expr,
    message = function(m) {
      text <- sub("\n$", "", conditionMessage(m))
      if (!grepl("^(Query|Sizes|Construction) \\[", text)) {
        messages <<- c(messages, text)
      }
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, messages = messages, warnings = warnings)
}

# Every level of a DAG holds all the source's minimal sets up to its cap (each
# effect at least its smallest sets), and lists exactly the source sets inside it.
expect_level_holds_sets <- function(graph, dag, label = "") {
  record <- dagby_options(dag)$level
  for (effect in c("direct", "total")) {
    sets <- catalogue_sets(graph)[[effect]]
    expect_gt(length(sets), 0L)
    cap <- if (is.finite(record$cap)) max(record$cap, min(lengths(sets))) else Inf
    wanted <- sets[lengths(sets) <= cap]
    listed <- set_keys(.gd_adjustment_sets(dag, effect)$sets)
    expect_gt(length(wanted), 0L)
    expect_identical(intersect(set_keys(wanted), listed), set_keys(wanted), info = paste(label, effect, "holds"))
    inside <- sets[vapply(sets, function(set) all(set %in% names(dag)), logical(1))]
    expect_identical(listed, set_keys(inside), info = paste(label, effect, "exact"))
  }
}

# E reaches the endpoints only through A, and Z is a descendant of the outcome.
# Neither is used by any minimal adjustment set. One size only.
wide_edges <- "E -> A; A -> T; A -> B; B -> Y; C -> T; C -> D; D -> Y; T -> M -> Y; T -> Y; Y -> Z"
# Total and direct sets {C, P}, {C, R}, {A, B, P}, {A, B, R}: sizes of 5 and 7 nodes. W is unused.
ladder_edges <- "A -> T; A -> C; B -> T; B -> C; C -> Y; T -> Y; P -> T; P -> Q; Q -> R; P -> R; R -> Y; W -> T"
# Sets {C, G}, {A, B, G}, {C, D, E, F}, {A, B, D, E, F}: sizes of 4, 6 and 9 nodes.
tied_edges <- "A -> T; A -> C; B -> T; B -> C; C -> Y; D -> T; E -> T; F -> T; D -> G; E -> G; F -> G; G -> Y; T -> Y"
# Total sets {C, D}, {A, B, C}; direct sets {C, D, M}, {A, B, C, M}: sizes of 5 and 7 nodes.
mixed_edges <- "T -> M; M -> Y; T -> Y; C -> T; C -> Y; A -> T; A -> D; B -> T; B -> D; D -> Y"
# The ladder with a repeated quantity: the time summary merges C_t1 and C_t2.
repeated_edges <- "A -> T; A -> C_t1; B -> T; B -> C_t1; C_t1 -> C_t2; C_t2 -> Y; C_t1 -> Y; T -> Y; P -> T; P -> Q; Q -> R; P -> R; R -> Y"
# Sets {C} and {A, B}: sizes of 3 and 5 nodes.
fork_edges <- "A -> T; A -> C; B -> T; B -> C; C -> Y; T -> Y"
# One minimal set {A}: one size of 3 nodes.
single_edges <- "A -> T; A -> Y; T -> Y"
# Total sets {C}, {A, B}; direct sets of 5 and 6 members (C or A and B, with D and three mediators): sizes of 6 and 9 nodes.
# At 10 steps no direct set is found.
mediator_edges <- "A -> T; A -> C; B -> T; B -> C; C -> Y; T -> M1; M1 -> Y; T -> M2; M2 -> Y; T -> M3; M3 -> Y; T -> Y; D -> M1; D -> Y"
# Two exposures. The exact treatment parents Q4 and S1 are kept at every size, and no total set exists:
# the levels of set sizes 3 and Inf draw the same 9-node graph.
pinned_edges <- paste("Q1_t1 -> Y; Q2_t1 -> Q2_t2; Q2_t1 -> S2; Q2_t1 -> S3; Q2_t1 -> T_t1; Q2_t1 -> Y; Q3_t2 -> Q2_t2;",
  "Q4_t0 -> Q1_t1; Q4_t0 -> Q2_t1; Q4_t0 -> S2; Q4_t0 -> S3; Q4_t0 -> T_t1; S1 -> Q1_t1; S1 -> Q2_t1; S1 -> Q2_t2;",
  "S1 -> S2; S1 -> T_t0; S2 -> Q1_t1; S2 -> Q3_t2; S2 -> T_t1; S2 -> Y; T_t0 -> Q2_t1; T_t0 -> Q2_t2; T_t1 -> Q2_t2; T_t1 -> Y")
pinned_query <- function() ranked_query(pinned_edges, exposures = c("T_t0", "T_t1"))

# The DAG of one rule row, built directly: row 1 is rule 0, the member search alone.
rule_dag <- function(graph, row = 1L, latents = NULL, unroll = FALSE) {
  query <- .cd_size_query(graph, NULL)
  rules <- .cd_rule_rows(query, TRUE)
  .cd_maximal(graph, latents, NULL, FALSE, TRUE, list(unroll = unroll), query, rules, rules$rows[[row]])
}

test_that("dagby_construct can request the largest DAG explicitly", {
  expect_identical(formals(dagby_construct)$size, "smallest")
  .cf_cache$clear()
  graph <- ranked_query(wide_edges)
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_s3_class(dag, "cf_dag")
  expect_false("dag_by_query" %in% class(dag))
  expect_false("dag_by_query" %in% class(graph))
  expect_identical(.gd_state(dag)$construction$representation, "maximal")
  expect_identical(dagby_options(dag)$search_method, "maximal-quantities-v1")
  kept <- dagby_options(dag)$quantities
  truth <- catalogue_quantities(graph)
  expect_identical(sort(kept$direct$members), truth$direct)
  expect_identical(sort(kept$total$members), truth$total)
  expect_identical(truth$total, c("A", "B", "C", "D"))
  expect_identical(truth$direct, c("A", "B", "C", "D", "M"))
  expect_identical(kept$direct$undecided, character())
  # Rule 4 adds the dead end Z and rule 6 the parentless variable E.
  expect_setequal(names(dag), c("T", "Y", "A", "B", "C", "D", "M", "E", "Z"))
  expect_setequal(edge_keys(dag), edge_keys(graph))
  expect_identical(dagby_options(dag)$level[c("rule", "top")], list(rule = 6L, top = TRUE))
  expect_identical(dagby_options(dag)$kept_by_rule, c(dead_ends = 1L, exogenous_causes = 1L))
  expect_identical(.gd_signature(dag), .gd_signature(dagby_construct(graph, Inf, quiet = TRUE)))
  zero <- rule_dag(graph)
  expect_setequal(names(zero), c("T", "Y", "A", "B", "C", "D", "M"))
  expect_setequal(edge_keys(zero), setdiff(edge_keys(graph), c("E -> A", "Y -> Z")))
})

test_that("the largest DAG lists exactly the source's minimal sets", {
  graph <- ranked_query(wide_edges)
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  for (effect in c("direct", "total")) {
    expect_gt(length(.gd_adjustment_sets(graph, effect)$sets), 0L)
    expect_identical(set_keys(.gd_adjustment_sets(dag, effect)$sets), set_keys(.gd_adjustment_sets(graph, effect)$sets))
  }
})

test_that("each size holds every minimal set up to its cap and no other set", {
  fixtures <- list(ladder = list(ladder_edges, c(5L, 7L), c(2, Inf)),
    tied = list(tied_edges, c(4L, 6L, 9L), c(2, 3, Inf)),
    mixed = list(mixed_edges, c(5L, 7L), c(2, Inf)))
  for (name in names(fixtures)) {
    graph <- ranked_query(fixtures[[name]][[1L]])
    .cf_cache$clear()
    for (index in seq_along(fixtures[[name]][[2L]])) {
      dag <- suppressWarnings(dagby_construct(graph, fixtures[[name]][[2L]][[index]], quiet = TRUE))
      expect_identical(length(names(dag)), fixtures[[name]][[2L]][[index]], info = name)
      expect_identical(dagby_options(dag)$level$cap, fixtures[[name]][[3L]][[index]], info = name)
      expect_level_holds_sets(graph, dag, paste(name, index))
    }
  }
})

test_that("the smallest size holds every smallest set of each effect", {
  # Total sets have 2 members, direct sets 3 (the mediator): the direct effect keeps its own smallest sets.
  graph <- ranked_query(mixed_edges)
  smallest <- dagby_construct(graph, 0, quiet = TRUE)
  expect_setequal(names(smallest), c("T", "Y", "C", "D", "M"))
  expect_identical(dagby_options(smallest)$level$cap, 2)
  expect_true(dagby_options(smallest)$level$smallest)
  sets <- catalogue_sets(graph)
  expect_identical(set_keys(.gd_adjustment_sets(smallest, "direct")$sets), "C+D+M")
  expect_identical(set_keys(.gd_adjustment_sets(smallest, "total")$sets), "C+D")
  expect_gt(min(lengths(sets$direct)), min(lengths(sets$total)))
  # Two smallest sets of one effect stay together.
  ladder <- ranked_query(ladder_edges)
  both <- dagby_construct(ladder, 0, quiet = TRUE)
  expect_setequal(names(both), c("T", "Y", "C", "P", "R"))
  expect_identical(set_keys(.gd_adjustment_sets(both, "total")$sets), c("C+P", "C+R"))
  expect_identical(set_keys(.gd_adjustment_sets(both, "direct")$sets), c("C+P", "C+R"))
})

test_that("rule 0 holds the union of every minimal set and the largest DAG adds the rule nodes", {
  for (edges in c(ladder_edges, tied_edges, mixed_edges)) {
    graph <- ranked_query(edges)
    .cf_cache$clear()
    top <- suppressWarnings(dagby_construct(graph, Inf, quiet = TRUE))
    truth <- catalogue_quantities(graph)
    union <- sort(unique(c("T", "Y", unlist(truth, use.names = FALSE))))
    zero <- suppressWarnings(dagby_construct(graph, length(union), quiet = TRUE))
    expect_identical(sort(names(zero)), union, info = edges)
    expect_identical(dagby_options(zero)$level$rule, 0L)
    expect_identical(dagby_options(zero)$search_method, "maximal-quantities-v1")
    rules <- .cd_rule_nodes(.cd_prepare(graph, NULL, NULL)$analysis)
    expect_identical(sort(names(top)), sort(union(union, rules[[6L]])), info = edges)
    expect_true(dagby_options(top)$level$top)
    expect_identical(dagby_options(top)$level$cap, Inf)
  }
})

test_that("a size snaps to the closest available size", {
  graph <- ranked_query(tied_edges)
  nodes <- function(size) length(names(suppressWarnings(dagby_construct(graph, size, quiet = TRUE))))
  expect_identical(nodes(0), 4L)
  expect_identical(nodes(1), 4L)
  expect_identical(nodes(4), 4L)
  expect_identical(nodes(6), 6L)
  expect_identical(nodes(7), 6L)
  expect_identical(nodes(8), 9L)
  expect_identical(nodes(9), 9L)
  expect_identical(nodes(500), 9L)
  expect_identical(nodes(Inf), 9L)
})

test_that("a tie between sizes goes up and warns with the call for the other size", {
  graph <- ranked_query(tied_edges)
  tie <- heard(dagby_construct(graph, 5))
  expect_identical(tie$warnings, "Sizes 4 and 6 are equally close to 5 variables. Returned 6. Use dagby_prune(dag, 2) for 4.")
  expect_identical(length(names(tie$value)), 6L)
  ladder <- ranked_query(ladder_edges)
  tie <- heard(dagby_construct(ladder, 6))
  expect_identical(tie$warnings, "Sizes 5 and 7 are equally close to 6 variables. Returned 7. Use dagby_prune(dag, 2) for 5.")
  expect_identical(length(names(tie$value)), 7L)
  # No tie, no warning.
  expect_identical(heard(dagby_construct(graph, 7))$warnings, character())
  # The warning's call gives the other size.
  other <- suppressMessages(dagby_prune(heard(dagby_construct(graph, 5))$value, 2))
  expect_identical(length(names(other)), 4L)
})

test_that("messages name the size returned and its neighbours", {
  graph <- ranked_query(tied_edges)
  # Guidance is based on the resolved choices, independent of earlier requests.
  .cf_cache$clear()
  expect_identical(heard(dagby_construct(graph, 0))$messages, "DAG with 4 variables. This is the smallest size. Next larger size: 6 variables. Use dagby_grow(dag, 2).")
  .cf_cache$clear()
  expect_identical(heard(dagby_construct(graph, size = "largest"))$messages, "DAG with 9 variables. This is the largest size. Next smaller size: 6 variables. Use dagby_prune(dag, 3).")
  .cf_cache$clear()
  middle <- heard(dagby_construct(graph, 6))
  expect_identical(middle$messages, "DAG with 6 variables. Next smaller size: 4 variables. Use dagby_prune(dag, 2). Next larger size: 9 variables. Use dagby_grow(dag, 3).")
  # Once listed, every request names its neighbours.
  expect_identical(heard(dagby_construct(graph, 0))$messages, "DAG with 4 variables. This is the smallest size. Next larger size: 6 variables. Use dagby_grow(dag, 2).")
  expect_identical(heard(dagby_construct(graph, 9))$messages, "DAG with 9 variables. This is the largest size. Next smaller size: 6 variables. Use dagby_prune(dag, 3).")
  expect_identical(heard(dagby_construct(graph, size = "largest"))$messages, "DAG with 9 variables. This is the largest size. Next smaller size: 6 variables. Use dagby_prune(dag, 3).")
  expect_identical(heard(dagby_construct(graph, 6, quiet = TRUE))$messages, character())
  expect_identical(heard(dagby_construct(ranked_query(wide_edges), 3))$messages, "DAG with 7 variables. This is the smallest size. Requested 3 variables. Returned 7. Next larger size: 8 variables. Use dagby_grow(dag, 1).")
  .cf_cache$clear()
  expect_identical(heard(dagby_construct(ranked_query(wide_edges), size = "largest"))$messages, "DAG with 9 variables. This is the largest size. Next smaller size: 8 variables. Use dagby_prune(dag, 1).")
  # A DAG that knows its sizes names both neighbours.
  small <- heard(dagby_construct(dagby_construct(graph, 6, quiet = TRUE), 0))
  expect_identical(small$messages, "DAG with 4 variables. This is the smallest size. Next larger size: 6 variables. Use dagby_grow(dag, 2).")
})

test_that("sizes are validated", {
  graph <- ranked_query(ladder_edges)
  for (bad in list(-1, NA_real_, c(1, 2), NULL, TRUE)) {
    expect_error(dagby_construct(graph, bad, quiet = TRUE), "size must be one number of at least 0, or Inf.", fixed = TRUE)
  }
  expect_error(dagby_construct(graph, "a", quiet = TRUE),
    'size must be "smallest", "largest" or a number of at least 0.', fixed = TRUE)
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_error(dagby_prune(dag, -1, quiet = TRUE), "size must be one number of at least 0, or Inf.", fixed = TRUE)
  expect_error(dagby_grow(dag, "a", quiet = TRUE), "size must be one number of at least 0, or Inf.", fixed = TRUE)
  expect_error(dagby_construct(graph, 0, quiet = "yes"), "quiet")
})

test_that("dagby_construct accepts a constructed DAG and rebuilds it from the saved query", {
  graph <- ranked_query(ladder_edges)
  top <- dagby_construct(graph, quiet = TRUE, size = "largest")
  small <- dagby_construct(graph, 0, quiet = TRUE)
  resized <- dagby_construct(top, 0, quiet = TRUE)
  expect_identical(.gd_signature(resized), .gd_signature(small))
  expect_identical(dagby_source(resized), graph)
  expect_identical(.gd_signature(dagby_construct(small, quiet = TRUE, size = "largest")), .gd_signature(top))
  # The saved latents and keep_nodes are used unless supplied.
  kept_result <- heard(dagby_construct(graph, 0, latents = "C", keep_nodes = "W", quiet = TRUE))
  expect_identical(kept_result$warnings, c(
    "Every direct adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: C).",
    "Every total adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: C)."))
  kept <- kept_result$value
  expect_identical(dagby_options(dagby_construct(kept, Inf, quiet = TRUE))$arguments[c("latents", "keep_nodes")],
    list(latents = "C", keep_nodes = "W"))
  replaced <- suppressWarnings(dagby_construct(kept, 0, latents = "P", quiet = TRUE))
  expect_identical(dagby_options(replaced)$arguments[c("latents", "keep_nodes")], list(latents = "P", keep_nodes = "W"))
  cleared <- dagby_construct(kept, 0, latents = NULL, keep_nodes = NULL, quiet = TRUE)
  expect_identical(dagby_options(cleared)$arguments[c("latents", "keep_nodes")], list(latents = NULL, keep_nodes = NULL))
  expect_identical(.gd_signature(cleared), .gd_signature(small))
  # A plan is rebuilt too.
  plan <- dagby_construct(graph, 0, dry_run = TRUE, quiet = TRUE)
  expect_identical(.gd_signature(dagby_construct(plan, Inf, quiet = TRUE)), .gd_signature(top))
})

test_that("dagby_prune and dagby_grow work by amounts and round trip", {
  .cf_cache$clear()
  graph <- ranked_query(tied_edges)
  top <- dagby_construct(graph, quiet = TRUE, size = "largest")
  six <- suppressMessages(dagby_prune(top, 3))
  expect_identical(length(names(six)), 6L)
  expect_identical(.gd_signature(six), .gd_signature(dagby_construct(graph, 6, quiet = TRUE)))
  expect_identical(length(names(suppressMessages(dagby_prune(top, 5)))), 4L)
  expect_identical(length(names(suppressMessages(dagby_prune(top, 1)))), 6L)
  # A tie takes the larger change, with the call for the other size.
  tie <- heard(dagby_prune(top, 4))
  expect_identical(length(names(tie$value)), 4L)
  expect_identical(tie$warnings, "Sizes 4 and 6 are equally close to 5 variables. Returned 4. Use dagby_grow(dag, 2) for 6.")
  small <- dagby_construct(graph, 0, quiet = TRUE)
  expect_identical(length(names(suppressMessages(dagby_grow(small, 2)))), 6L)
  expect_identical(length(names(suppressMessages(dagby_grow(small, 1)))), 6L)
  expect_identical(length(names(suppressMessages(dagby_grow(small, 5)))), 9L)
  growth <- heard(dagby_grow(small, 1))
  expect_identical(growth$warnings, character())
  expect_identical(growth$messages, "DAG with 6 variables. 2 variables were added, not 1. Next smaller size: 4 variables. Use dagby_prune(dag, 2). Next larger size: 9 variables. Use dagby_grow(dag, 3).")
  # Round trips return the same DAG.
  expect_identical(.gd_signature(suppressMessages(dagby_grow(suppressMessages(dagby_prune(top, 3)), 3))), .gd_signature(top))
  expect_identical(.gd_signature(suppressMessages(dagby_prune(suppressMessages(dagby_grow(small, 2)), 2))), .gd_signature(small))
  expect_identical(.gd_signature(suppressMessages(dagby_grow(suppressMessages(dagby_prune(six, 2)), 2))), .gd_signature(six))
  # Inf goes to the smallest or the largest.
  expect_identical(.gd_signature(suppressMessages(dagby_prune(top))), .gd_signature(small))
  expect_identical(.gd_signature(suppressMessages(dagby_grow(small))), .gd_signature(top))
  expect_identical(dagby_source(six), graph)
})

test_that("a DAG at its smallest or largest size is returned unchanged", {
  .cf_cache$clear()
  graph <- ranked_query(tied_edges)
  top <- dagby_construct(graph, quiet = TRUE, size = "largest")
  small <- dagby_construct(graph, 0, quiet = TRUE)
  expect_identical(heard(dagby_prune(small))$messages, "No smaller DAG is available. The DAG was returned unchanged.")
  expect_identical(heard(dagby_grow(top))$messages, "No larger DAG is available. The DAG was returned unchanged.")
  expect_identical(heard(dagby_prune(small, 2))$messages, "No smaller DAG is available. The DAG was returned unchanged.")
  expect_identical(heard(dagby_grow(top, 2))$messages, "No larger DAG is available. The DAG was returned unchanged.")
  expect_identical(suppressMessages(dagby_prune(small, 3)), small)
  expect_identical(suppressMessages(dagby_grow(top, 3)), top)
  expect_identical(heard(dagby_prune(small, quiet = TRUE))$messages, character())
  # An amount of zero changes nothing and says nothing.
  middle <- suppressMessages(dagby_prune(top, 3))
  silent <- heard(dagby_prune(middle, 0))
  expect_identical(silent$value, middle)
  expect_identical(silent$messages, character())
  expect_identical(heard(dagby_grow(middle, 0))$value, middle)
  # Only smaller sizes are candidates for dagby_prune, so a middle DAG cannot grow by pruning.
  expect_identical(length(names(suppressMessages(dagby_prune(middle, 100)))), 4L)
  # The whole sizes of a small DAG are known from the larger one it came from.
  expect_identical(length(names(suppressMessages(dagby_grow(middle, 100)))), 9L)
  # A graph with one size reports it.
  single <- dagby_construct(ranked_query(single_edges), quiet = TRUE, size = "largest")
  expect_identical(heard(dagby_prune(single))$messages, "No smaller DAG is available. The DAG was returned unchanged.")
  expect_identical(heard(dagby_grow(single))$messages, "No larger DAG is available. The DAG was returned unchanged.")
})

test_that("a plain DAG is accepted with a warning and its latent tags are unmeasured", {
  plain <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U [latent]; A -> T; A -> C; B -> T; B -> C; C -> Y; T -> Y; U -> T; U -> Y }")
  for (request in c("prune", "grow")) {
    changed <- get(paste0("dagby_", request))
    result <- heard(changed(plain, quiet = TRUE))
    expect_true(paste0("dagby_", request, "() is intended for working with DAGs created by querying the process-based structural causal model. It may have unexpected results when used with knowledge-based DAGs.") %in%
      result$warnings)
  }
  pruned <- heard(dagby_prune(plain, quiet = TRUE))
  expect_s3_class(pruned$value, "cf_dag")
  expect_setequal(names(pruned$value), c("C", "T", "U", "Y"))
  expect_identical(as.character(dagitty::latents(pruned$value)), "U")
  expect_identical(dagby_options(pruned$value)$arguments$latents, "U")
  expect_identical(dagby_options(pruned$value)$level$nodes, 4L)
  expect_identical(.gd_signature(suppressWarnings(dagby_prune(plain, 1, quiet = TRUE))), .gd_signature(pruned$value))
  # The largest DAG of a plain graph is not larger than the graph.
  grown <- suppressWarnings(dagby_grow(plain, 1, quiet = TRUE))
  expect_identical(grown, plain)
  expect_error(dagby_prune("T -> Y", quiet = TRUE), "Supply a DAG.", fixed = TRUE)
  adjusted <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A [adjusted]; A -> T; A -> Y; T -> Y }")
  expect_error(suppressWarnings(dagby_prune(adjusted, quiet = TRUE)), "The DAG has adjusted or selected nodes. Remove those tags first.", fixed = TRUE)
})

test_that("repeated occasions count once and a quantity outside every minimal member set is collapsed", {
  graph <- ranked_query("A_t1 -> T; A_t1 -> A_t2; A_t2 -> Y; B -> A_t2; B -> Y; C_t1 -> T; C_t1 -> C_t2; C_t2 -> Y; D -> Y; T -> Y")
  # B is in the exact-occasion set {A_t2, B, C_t2}, but {A, C} already suffices.
  expect_true(.gd_validate_sets(graph, list(c("A_t2", "B", "C_t2")), "total"))
  dag <- rule_dag(graph)
  kept <- dagby_options(dag)$quantities
  expect_identical(kept$total$members, c("A", "C"))
  expect_identical(kept$total$members, catalogue_quantities(graph)$total)
  expect_setequal(names(dagby_temporal(dag)), c("T", "Y", "A_t1", "A_t2", "C_t1", "C_t2"))
  expect_true(dagby_options(dag)$summary$applied)
  expect_setequal(names(dag), c("T", "Y", "A", "C"))
  # B and D are outside causes, kept by rule 5.
  expect_setequal(names(dagby_construct(graph, quiet = TRUE, size = "largest")), c("T", "Y", "A", "B", "C", "D"))
})

test_that("a listing cut short is reported as not fully checked", {
  # Six independent confounding chains give 64 minimal sets, more than the 50 listed.
  chains <- paste(sprintf("A%d -> T; A%d -> B%d; B%d -> Y", 1:6, 1:6, 1:6, 1:6), collapse = "; ")
  graph <- ranked_query(paste(chains, "; T -> Y"))
  expect_length(.gd_adjustment_sets(graph, "total")$sets, 64L)
  dag <- expect_no_warning(dagby_construct(graph, quiet = TRUE, size = "largest"))
  checks <- dagby_checks(dag)
  expect_identical(checks$total$listed, 50L)
  expect_false(checks$total$listing_complete)
  expect_true(checks$total$sets_minimal && checks$total$listed_valid_on_source)
  expect_identical(checks$total$status, "not fully checked")
  expect_setequal(dagby_options(dag)$quantities$total$members, c(paste0("A", 1:6), paste0("B", 1:6)))
  expect_setequal(names(dag), names(graph))
})

test_that("a time summary that only renames nodes is told from one that merges occasions", {
  renamed <- dagby_construct(ranked_query("A_t1 -> T; A_t1 -> Y; T -> Y"), quiet = TRUE, size = "largest")
  expect_setequal(names(dagby_temporal(renamed)), c("T", "Y", "A_t1"))
  expect_setequal(names(renamed), c("T", "Y", "A"))
  expect_identical(dagby_options(renamed)$summary[c("applied", "merges_occasions")], list(applied = TRUE, merges_occasions = FALSE))
  merged <- dagby_construct(ranked_query("A_t1 -> A_t2; A_t2 -> T; A_t1 -> Y; A_t2 -> Y; T -> Y"), quiet = TRUE, size = "largest")
  expect_setequal(names(dagby_temporal(merged)), c("T", "Y", "A_t1", "A_t2"))
  expect_setequal(names(merged), c("T", "Y", "A"))
  record <- dagby_options(merged)$summary
  expect_true(record$applied && record$merges_occasions)
  # Every set the merged display lists was validated on the source.
  expect_identical(record$displayed_sets_checked,
    length(.gd_adjustment_sets(merged, "direct")$sets) + length(.gd_adjustment_sets(merged, "total")$sets))
})

test_that("a merging summary whose sets cannot all be checked requires explicit time-unrolled output", {
  # Ten chains give 1024 minimal sets, and Q has two occasions the summary would merge.
  chains <- paste(sprintf("A%d -> T; A%d -> B%d; B%d -> Y", 1:10, 1:10, 1:10, 1:10), collapse = "; ")
  graph <- ranked_query(paste(chains, "; Q_t1 -> T; Q_t1 -> Q_t2; Q_t2 -> Y; T -> Y"))
  expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "unroll = TRUE", fixed = TRUE)
  dag <- expect_no_warning(dagby_construct(graph, quiet = TRUE, unroll = TRUE, size = "largest"))
  record <- dagby_options(dag)$summary
  expect_false(record$applied)
  expect_match(dagby_options(dag)$level$static_reason, "merges occasions and lists more than 1000 direct adjustment sets", fixed = TRUE)
  expect_true(all(c("Q_t1", "Q_t2") %in% names(dag)))
  expect_identical(.gd_signature(dag), .gd_signature(dagby_temporal(dag)))
  expect_setequal(dagby_options(dag)$quantities$total$members, c(paste0("A", 1:10), paste0("B", 1:10), "Q"))
})

test_that("declared unmeasured variables that every set needs are reported", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y; B -> A; B -> Y")
  seen <- character()
  dag <- withCallingHandlers(dagby_construct(graph, latents = "A", quiet = TRUE, size = "largest"), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, c(
    "Every direct adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A).",
    "Every total adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A)."))
  checks <- dagby_checks(dag)
  expect_identical(checks$total$unmeasured, "A")
  expect_identical(checks$total[c("availability", "temporal_availability")],
    list(availability = "needs an unmeasured variable", temporal_availability = "needs an unmeasured variable"))
  expect_identical(checks$total$status, "preserved")
  expect_length(dagitty::adjustmentSets(dag, effect = "total"), 0L)
  # B is collapsed, so declaring it changes nothing that is reported.
  collapsed <- expect_no_warning(dagby_construct(graph, 0, latents = "B", quiet = TRUE))
  expect_identical(dagby_checks(collapsed)$total$unmeasured, character())
  expect_identical(dagby_checks(collapsed)$total$availability, "available now")
  # With a second sufficient quantity an adjustment set remains without A.
  chain <- ranked_query("A -> T; A -> B; B -> Y; T -> Y")
  spared <- expect_no_warning(dagby_construct(chain, latents = "A", quiet = TRUE, size = "largest"))
  expect_identical(dagby_checks(spared)$total$unmeasured, "A")
  expect_identical(dagby_checks(spared)$total$availability, "available now")
  expect_identical(.gd_adjustment_sets(spared, "total")$sets, list("B"))
  seen <- character()
  both <- withCallingHandlers(dagby_construct(chain, latents = c("A", "B"), quiet = TRUE, size = "largest"), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, c(
    "Every direct adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A, B).",
    "Every total adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A, B)."))
  expect_identical(dagby_checks(both)$total$availability, "needs an unmeasured variable")
})

test_that("availability is stated for the returned graph and for the temporal graph", {
  # A set without the declared occasion exists on the source, but it needs C,
  # which no minimal member set uses, so the maximal representation collapses it.
  graph <- ranked_query("B_t1 -> T; B_t1 -> Y; B_t1 -> B_t2; B_t2 -> T; C -> B_t2; C -> Y; T -> Y")
  prepared <- .cd_prepare(graph, NULL, NULL)
  expect_true(.gd_validate_sets(.cd_set_latents(prepared$selection, "B_t2"), list(c("B_t1", "C")), "total"))
  dag <- suppressWarnings(rule_dag(graph, latents = "B_t2"))
  expect_false("C" %in% names(dagby_temporal(dag)))
  expect_identical(dagby_checks(dag)$total[c("unmeasured", "availability", "temporal_availability")],
    list(unmeasured = "B_t2", availability = "needs an unmeasured variable",
    temporal_availability = "needs an unmeasured variable"))
  expect_identical(.gd_adjustment_sets(dag, "total")$sets, list())
  # The largest DAG keeps C, an outside cause, so its temporal graph lists a set without the declared occasion.
  top <- suppressWarnings(dagby_construct(graph, latents = "B_t2", quiet = TRUE, size = "largest"))
  expect_true("C" %in% names(top))
  expect_identical(dagby_checks(top)$total$temporal_availability, "available now")
  # A merged node is marked when one of its occasions is declared, so the
  # returned graph can lose a set that the temporal graph still lists.
  merged <- ranked_query("A_t1 -> A_t2; A_t2 -> T; A_t1 -> Y; A_t2 -> Y; T -> Y")
  seen <- character()
  spare <- withCallingHandlers(dagby_construct(merged, latents = "A_t1", quiet = TRUE, size = "largest"), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, c(
    "A stands for both observed and unmeasured occasions (unmeasured: A_t1, observed: A_t2). It is shown unmeasured.",
    "Every direct adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A_t1). dagby_temporal() lists a set that uses none.",
    "Every total adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A_t1). dagby_temporal() lists a set that uses none."))
  expect_identical(as.character(dagitty::latents(spare)), "A")
  expect_identical(dagby_checks(spare)$total[c("unmeasured", "availability", "temporal_availability")],
    list(unmeasured = "A_t1", availability = "needs an unmeasured variable", temporal_availability = "available now"))
  expect_identical(.gd_adjustment_sets(spare, "total")$sets, list())
  expect_identical(.gd_adjustment_sets(dagby_temporal(spare), "total")$sets, list("A_t2"))
  needed <- suppressWarnings(dagby_construct(merged, latents = "A_t2", quiet = TRUE, size = "largest"))
  expect_identical(dagby_checks(needed)$total$temporal_availability, "needs an unmeasured variable")
  # A declared node kept only on request belongs to neither effect's record.
  kept <- expect_no_warning(dagby_construct(ranked_query(wide_edges), latents = "E", keep_nodes = "E", quiet = TRUE, size = "largest"))
  expect_identical(as.character(dagitty::latents(kept)), "E")
  expect_identical(dagby_checks(kept)$total[c("unmeasured", "availability")], list(unmeasured = character(), availability = "available now"))
  absent <- suppressWarnings(dagby_construct(ranked_query("T -> Y; T <-> Y; A -> Y"), latents = "A", quiet = TRUE, size = "largest"))
  expect_identical(dagby_checks(absent)$total$availability, "no adjustment set")
})

test_that("a merging summary that hides a smaller set requires explicit time-unrolled output", {
  graph <- ranked_query("A_t1 -> B_t1; B_t1 -> A_t2; A_t2 -> B_t2; A_t1 -> T; B_t2 -> Y; A_t2 -> Y; T -> Y; B_t1 -> T")
  expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "unroll = TRUE", fixed = TRUE)
  dag <- dagby_construct(graph, quiet = TRUE, unroll = TRUE, size = "largest")
  record <- dagby_options(dag)$summary
  expect_false(record$applied)
  expect_match(dagby_options(dag)$level$static_reason, "is not a minimal valid set of the time summary", fixed = TRUE)
  expect_null(record$merges_occasions)
  expect_identical(.gd_signature(dag), .gd_signature(dagby_temporal(dag)))
})

test_that("the quantities kept match the complete catalogue on several fixtures", {
  fixtures <- c(wide_edges,
    "A_t1 -> T; A_t1 -> A_t2; A_t2 -> Y; B -> A_t2; B -> Y; C_t1 -> T; C_t1 -> C_t2; C_t2 -> Y; D -> Y; T -> Y",
    "A_t1 -> B_t2; B_t1 -> A_t2; A_t1 -> T; B_t1 -> T; A_t2 -> Y; B_t2 -> Y; T -> Y",
    "A_t1 -> B_t1; B_t1 -> A_t2; A_t2 -> B_t2; B_t2 -> A_t3; A_t1 -> T; A_t3 -> Y; B_t2 -> Y; T -> Y",
    "T -> M; M -> Y; M <-> Y; C -> T; C -> Y")
  for (edges in fixtures) {
    graph <- ranked_query(edges)
    dag <- suppressWarnings(dagby_construct(graph, quiet = TRUE, unroll = TRUE, size = "largest"))
    truth <- catalogue_quantities(graph)
    kept <- dagby_options(dag)$quantities
    for (effect in c("direct", "total")) {
      expect_identical(sort(kept[[effect]]$members), truth[[effect]], info = paste(edges, effect))
      expect_identical(kept[[effect]]$undecided, character(), info = paste(edges, effect))
      expect_identical(names(kept[[effect]]), c("state", "members", "undecided"))
    }
  }
})

# Every field name the saved state of a sized DAG may hold. A name that is not on it, such as sets or options, fails.
state_fields <- c("added_nodes", "adjusted", "applied", "arguments", "availability", "cap", "certified", "checks",
  "construction", "declared", "direct", "direct_path_nodes", "display_nodes", "displayed_sets_checked",
  "edges", "effect", "exposures", "identification", "keep_nodes", "kept_by_rule", "labels",
  "latent", "latents", "level", "limitations", "listed", "listed_valid_on_source", "listing_complete", "members",
  "merges_occasions", "no_set_kept", "nodes", "outcomes", "parents", "partial", "quantities", "reason",
  "representation", "rule", "scope", "search", "search_method", "segment_alternatives",
  "segment_edge_ties", "selected", "semantics", "set_size", "sets_found", "sets_minimal", "signature", "size",
  "sizes", "sizes_dropped", "sizes_steps", "smallest", "state", "status", "steps", "summary",
  "temporal_availability", "temporal_nodes", "top", "total", "type", "undecided", "unmeasured",
  "metric", "variables", "static_available", "static_reason", "unroll", "select", "candidate", "id",
  "alternatives", "description", "alternatives_complete", "binding", "version", "source", "member_steps",
  "ranked_steps", "sizes_complete", "alternatives_limit")

# Every name in a list, at any depth, including the columns of data frames.
field_names <- function(x) {
  if (!is.list(x)) {
    return(character())
  }
  return(c(names(x), unlist(lapply(x, field_names), use.names = FALSE)))
}

test_that("the short record stores exact choices and size completeness without adjustment sets", {
  graph <- ranked_query(ladder_edges)
  .cf_cache$clear()
  small <- dagby_construct(graph, 0, quiet = TRUE)
  largest <- dagby_construct(graph, quiet = TRUE, size = "largest")
  zero <- dagby_construct(graph, 7, quiet = TRUE)
  listed <- dagby_prune(largest, 3, quiet = TRUE)
  construction_fields <- c("search_method", "representation", "arguments", "quantities", "summary", "identification",
    "unmeasured", "limitations", "scope", "direct_path_nodes", "level", "metric", "sizes", "candidate",
    "alternatives", "alternatives_complete", "sizes_complete")
  expect_setequal(names(.gd_state(small)$construction), c(construction_fields, "sizes_steps"))
  expect_setequal(names(.gd_state(largest)$construction), c(construction_fields, "kept_by_rule", "sizes_steps"))
  expect_setequal(names(.gd_state(zero)$construction), c(construction_fields, "sizes_steps"))
  expect_setequal(names(.gd_state(listed)$construction), c(construction_fields, "sizes_steps"))
  for (dag in list(listed, largest, zero, small)) {
    state <- attr(dag, "cf_provenance")
    expect_setequal(names(state), c("original_query", "temporal", "construction", "checks", "signature"))
    names_saved <- field_names(state[setdiff(names(state), c("original_query", "temporal"))])
    expect_identical(setdiff(names_saved, state_fields), character())
    expect_length(intersect(names_saved, c("sets", "options", "catalogue", "graphs")), 0L)
    expect_identical(dagby_source(dag), graph)
    expect_identical(.gd_signature(dagby_temporal(dag)), .gd_signature(dag))
    options <- dagby_options(dag)
    expect_identical(options$metric, "variables-v1")
    expect_identical(options$level$variables, length(names(dag)))
    expect_identical(options$level$nodes, length(names(dag)))
    expect_identical(options$alternatives$variables, options$level$variables)
    expect_identical(options$alternatives$selected, TRUE)
    expect_identical(options$alternatives$id, options$candidate$id)
    expect_type(options$candidate$signature, "raw")
    expect_true(options$alternatives_complete)
    checks <- dagby_checks(dag)
    expect_true(attr(checks, "current"))
    expect_identical(unname(vapply(checks, `[[`, character(1), "status")), c("preserved", "preserved"))
    expect_true(all(vapply(checks, function(check) check$sets_minimal && check$listed_valid_on_source, logical(1))))
    expect_true(all(vapply(checks, function(check) check$listing_complete, logical(1))))
    expect_identical(unname(vapply(checks, `[[`, character(1), "availability")), c("available now", "available now"))
    expect_identical(options$summary[c("applied", "merges_occasions")], list(applied = TRUE, merges_occasions = FALSE))
  }
  expect_identical(dagby_options(largest)$arguments, list(latents = NULL, keep_nodes = NULL,
    search = "exact", steps = 256L, size = 8L, unroll = FALSE))
  expect_identical(dagby_options(zero)$arguments, list(latents = NULL, keep_nodes = NULL,
    search = "exact", steps = 256L, size = 7L, unroll = FALSE))
  expect_identical(dagby_options(listed)$arguments, list(latents = NULL, keep_nodes = NULL,
    search = "auto", steps = 20000L, size = 5L, unroll = FALSE))
  expect_identical(dagby_options(largest)$kept_by_rule, c(outside_causes = 1L))
  expect_identical(dagby_options(largest)$search_method, "maximal-quantities-v1")
  expect_identical(dagby_options(zero)$search_method, "maximal-quantities-v1")
  expect_identical(dagby_options(small)$search_method, "quantity-ranked-v1")
  expect_identical(.gd_state(small)$construction$representation, "level")
  expect_identical(.gd_state(zero)$construction$representation, "maximal")
  expect_identical(dagby_options(largest)$level[c("cap", "rule", "certified", "smallest", "top", "variables")],
    list(cap = Inf, rule = 5L, certified = TRUE, smallest = FALSE, top = TRUE, variables = 8L))
  expect_identical(dagby_options(small)$quantities$total$set_size, 2)
  expect_identical(names(dagby_options(largest)$quantities$total), c("state", "members", "undecided"))
  expect_identical(dagby_options(small)$sizes, data.frame(variables = c(5L, 7L)))
  expect_identical(dagby_options(largest)$sizes, data.frame(variables = c(5L, 7L, 8L)))
  expect_false(dagby_options(small)$sizes_complete)
  expect_true(dagby_options(largest)$sizes_complete)
  expect_true(dagby_options(small)$level$certified)
  for (size in union(dagby_options(small)$sizes$variables, dagby_options(largest)$sizes$variables)) {
    reached <- dagby_select(small, size = size, quiet = TRUE)
    expect_identical(dagby_options(reached)$level$variables, size)
    expect_true(dagby_options(reached)$alternatives_complete)
    expect_true(dagby_options(reached)$level$certified)
    expect_level_holds_sets(graph, reached)
  }
  expect_identical(dagby_options(largest)$sizes_steps, 20000L)
  expect_identical(dagby_options(listed)$sizes_steps, 20000L)
  expect_null(dagby_options(listed)$sizes_dropped)
  .cf_cache$clear()
  expect_identical(dagby_options(dagby_construct(listed, 0, quiet = TRUE)), dagby_options(small))
})

test_that("the saved call replays the exact candidate, declarations and output mode", {
  graph <- ranked_query(tied_edges)
  built <- list(small = dagby_construct(graph, 0, quiet = TRUE), middle = dagby_construct(graph, 6, quiet = TRUE),
    top = dagby_construct(graph, quiet = TRUE, size = "largest"))
  expect_identical(vapply(built, function(dag) dagby_call(dag)$arguments$size, integer(1)),
    c(small = 4L, middle = 6L, top = 9L))
  for (dag in built) {
    call <- dagby_call(dag)
    expect_identical(call[[1L]], quote(causaliflower:::.cd_replay))
    expect_setequal(names(call)[-1L], c("graph", "arguments", "candidate", "dry_run"))
    expect_identical(call$graph, graph)
    expect_identical(call$candidate, dagby_options(dag)$candidate)
    .cf_cache$clear()
    expect_identical(eval(call), dag)
  }
  pruned <- dagby_prune(built$top, 3, quiet = TRUE)
  expect_identical(dagby_call(pruned)$arguments$size, 6L)
  .cf_cache$clear()
  expect_identical(eval(dagby_call(pruned)), pruned)
  declared <- dagby_construct(graph, 0, latents = "A", keep_nodes = "D", quiet = TRUE)
  expect_identical(eval(dagby_call(declared)), declared)
  expect_identical(dagby_call(declared)$arguments[c("latents", "keep_nodes")], list(latents = "A", keep_nodes = "D"))
  plan <- dagby_construct(graph, 0, dry_run = TRUE, quiet = TRUE)
  expect_true(dagby_call(plan)$dry_run)
  expect_identical(dagby_call(plan)$arguments$size, 4L)
  expect_identical(eval(dagby_call(plan)), plan)
})
test_that("different DAGs at the same variable count remain selectable", {
  graph <- ranked_query("T -> Y; T -> P_t0; P_t0 -> P_t1")
  full <- dagby_construct(graph, quiet = TRUE, size = "largest")
  choices <- dagby_options(full)$alternatives
  expect_identical(choices$variables, c(3L, 3L))
  expect_identical(choices$select, 1:2)
  expect_false(anyDuplicated(choices$id) > 0L)
  expect_setequal(names(dagby_temporal(full)), c("T", "Y", "P_t0", "P_t1"))
  other <- dagby_select(full, select = 2, quiet = TRUE)
  expect_setequal(names(dagby_temporal(other)), c("T", "Y", "P_t0"))
  expect_false(identical(.gd_signature(dagby_temporal(other)), .gd_signature(dagby_temporal(full))))
  expect_identical(dagby_options(other)$level$variables, 3L)
  expect_identical(dagby_options(other)$level$select, 2L)
  expect_identical(dagby_grow(other, 1, quiet = TRUE), other)
  expect_identical(.gd_signature(dagby_temporal(dagby_select(other, quiet = TRUE))),
    .gd_signature(dagby_temporal(other)))
  expect_identical(.gd_signature(dagby_temporal(dagby_select(other, select = 1, quiet = TRUE))),
    .gd_signature(dagby_temporal(full)))
})

test_that("levels with the same temporal graph are one alternative", {
  .cf_cache$clear()
  query <- pinned_query()
  first <- suppressWarnings(dagby_construct(query, 0, unroll = TRUE, quiet = TRUE))
  expect_identical(dagby_options(first)$level$variables, 7L)
  expect_identical(length(names(first)), 9L)
  expect_identical(nrow(dagby_options(first)$alternatives), 1L)
  grown <- suppressWarnings(dagby_grow(first, 1, quiet = TRUE))
  expect_identical(dagby_options(grown)$level$variables, 9L)
  expect_identical(length(names(grown)), 11L)
  back <- suppressWarnings(dagby_prune(grown, 2, quiet = TRUE))
  expect_identical(.gd_signature(back), .gd_signature(first))
  expect_identical(dagby_prune(back, 1, quiet = TRUE), back)
  expect_false(anyDuplicated(dagby_options(first)$alternatives$id) > 0L)
  expect_identical(dagby_options(first)$sizes$variables, c(7L, 9L))
})

test_that("identical candidates are listed once and each size keeps its own checks", {
  .cf_cache$clear()
  graph <- ranked_query(tied_edges)
  middle <- dagby_construct(graph, 6, quiet = TRUE)
  top <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_identical(dagby_options(top)$sizes$variables, c(4L, 6L, 9L))
  expect_identical(nrow(dagby_options(top)$alternatives), 1L)
  expect_null(dagby_options(top)$sizes_dropped)
  for (effect in c("direct", "total")) {
    expect_identical(dagby_checks(top)[[effect]]$sets_found, 3L, info = effect)
    expect_identical(dagby_checks(middle)[[effect]]$sets_found, 2L, info = effect)
    expect_length(.gd_adjustment_sets(top, effect)$sets, 4L)
    expect_length(.gd_adjustment_sets(middle, effect)$sets, 2L)
  }
})
test_that("without a certified smallest size only the smallest level and the top are offered", {
  withr::local_options(causaliflower.ranked_steps = 10)
  .cf_cache$clear()
  mediator <- ranked_query(mediator_edges)
  definitions <- .cd_level_definitions(.cd_size_query(mediator, NULL), 10L)
  expect_identical(vapply(definitions$effects, `[[`, character(1), "state"), c(direct = "not_found", total = "sets"))
  expect_length(definitions$levels, 1L)
  expect_false(definitions$levels[[1L]]$certified)
  small <- heard(dagby_construct(mediator, 0))
  expect_identical(small$warnings, "No direct adjustment set was found within the step budget. This is not evidence that no set exists.")
  expect_null(dagby_options(small$value)$quantities$direct$set_size)
  expect_identical(dagby_options(small$value)$quantities$total$set_size, 1)
  top <- heard(dagby_construct(mediator, 9))
  expect_identical(dagby_options(top$value)$sizes, data.frame(variables = c(6L, 9L)))
  again <- heard(dagby_construct(mediator, 6))
  # A size at the smallest level is answered from the list alone.
  expect_identical(again$messages, "DAG with 6 variables. This is the smallest size. Next larger size: 9 variables. Use dagby_grow(dag, 3).")
  expect_identical(again$warnings, small$warnings)
  # Sets found within the budget, but the smallest size is not certified.
  .cf_cache$clear()
  ladder <- ranked_query(ladder_edges)
  definitions <- .cd_level_definitions(.cd_size_query(ladder, NULL), 10L)
  expect_identical(vapply(definitions$effects, `[[`, character(1), "state"), c(direct = "sets", total = "sets"))
  expect_true(all(vapply(definitions$effects, function(effect) effect$certified_to < effect$smallest, logical(1))))
  expect_length(definitions$levels, 1L)
  expect_false(definitions$levels[[1L]]$certified)
  small <- heard(dagby_construct(ladder, 0))
  ranked <- paste0(c("direct", "total"), " options are ranked and certified complete up to 1 distinct source members. ",
    "Larger options may exist. Alternatives with the same member count as the smallest shown may be missing (step budget 10).")
  expect_identical(small$warnings, ranked)
  expect_identical(dagby_options(small$value)$level$certified, FALSE)
  expect_null(dagby_options(small$value)$quantities$direct$set_size)
  expect_null(dagby_options(small$value)$quantities$total$set_size)
  top <- heard(dagby_construct(ladder, 7))
  expect_identical(top$messages, "DAG with 7 variables. Next smaller size: 4 variables. Use dagby_prune(dag, 3). Next larger size: 8 variables. Use dagby_grow(dag, 1).")
  expect_identical(dagby_options(top$value)$sizes, data.frame(variables = c(4L, 7L, 8L)))
  expect_identical(heard(dagby_prune(top$value, 3))$warnings, ranked)
})

# Count the ranked searches and member searches of the code that follows.
spy_searches <- function(env = parent.frame()) {
  calls <- new.env()
  calls$reset <- function() {
    calls$steps <- integer()
    calls$members <- 0L
  }
  calls$reset()
  ranked <- get(".cd_ranked_search", envir = asNamespace("causaliflower"))
  members <- get(".cd_maximal_members", envir = asNamespace("causaliflower"))
  testthat::local_mocked_bindings(
    .cd_ranked_search = function(...) {
      calls$steps <- c(calls$steps, list(...)$steps)
      return(ranked(...))
    },
    .cd_maximal_members = function(...) {
      calls$members <- calls$members + 1L
      return(members(...))
    }, .package = "causaliflower", .env = env)
  return(calls)
}

test_that("requested discovery reuses its searches and defers other sizes until needed", {
  graph <- ranked_query(ladder_edges)
  .cf_cache$clear()
  calls <- spy_searches()
  small <- dagby_construct(graph, 0, quiet = TRUE)
  expect_identical(calls$steps, rep(200L, 2L))
  expect_true(all(calls$steps <= .cd_ranked_steps()))
  expect_identical(calls$members, 0L)
  expect_true(dagby_options(small)$level$certified)
  expect_true(dagby_options(small)$alternatives_complete)
  expect_false(dagby_options(small)$sizes_complete)
  calls$reset()
  again <- dagby_construct(graph, 0, quiet = TRUE)
  selected <- dagby_select(small, quiet = TRUE)
  expect_identical(calls$steps, integer())
  expect_identical(calls$members, 0L)
  expect_identical(.gd_signature(again), .gd_signature(small))
  expect_identical(dagby_options(selected)$candidate, dagby_options(small)$candidate)
  expect_identical(dagby_options(selected)$alternatives, dagby_options(small)$alternatives)
  .cf_cache$clear()
  calls$reset()
  cold <- dagby_construct(graph, 0, quiet = TRUE)
  expect_identical(calls$steps, rep(200L, 2L))
  expect_identical(calls$members, 0L)
  expect_identical(dagby_options(cold)$candidate, dagby_options(small)$candidate)
  expect_identical(dagby_options(cold)$alternatives, dagby_options(small)$alternatives)
  calls$reset()
  top <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_identical(calls$steps, rep(.cd_ranked_steps(), 2L))
  expect_identical(calls$members, 2L)
  expect_true(dagby_options(top)$level$certified)
  expect_true(dagby_options(top)$alternatives_complete)
  expect_true(dagby_options(top)$sizes_complete)
  calls$reset()
  expect_identical(dagby_options(dagby_construct(graph, quiet = TRUE, size = "largest")), dagby_options(top))
  expect_identical(dagby_options(dagby_construct(graph, 0, quiet = TRUE)), dagby_options(small))
  expect_identical(.gd_signature(dagby_grow(cold, Inf, quiet = TRUE)), .gd_signature(top))
  expect_identical(calls$steps, integer())
  expect_identical(calls$members, 0L)
})
test_that("sizes carried from a DAG built at another budget leave out its uncertified size", {
  graph <- ranked_query(ladder_edges)
  .cf_cache$clear()
  built <- withr::with_options(list(causaliflower.ranked_steps = 1), suppressWarnings(dagby_construct(graph, 6, quiet = TRUE)))
  expect_identical(length(names(built)), 6L)
  expect_identical(dagby_options(built)$sizes, data.frame(variables = c(6L, 7L, 8L)))
  expect_false(dagby_options(built)$level$certified)
  .cf_cache$clear()
  again <- heard(dagby_construct(built, 0))
  expect_identical(again$messages, "DAG with 5 variables. This is the smallest size. Next larger size: 7 variables. Use dagby_grow(dag, 2).")
  expect_identical(dagby_options(again$value)$sizes, data.frame(variables = c(5L, 7L)))
  expect_true(dagby_options(again$value)$level$certified)
  expect_true(dagby_options(again$value)$alternatives_complete)
  expect_false(dagby_options(again$value)$sizes_complete)
  for (size in dagby_options(again$value)$sizes$variables) {
    reached <- dagby_select(again$value, size = size, quiet = TRUE)
    expect_identical(dagby_options(reached)$level$variables, size)
    expect_true(dagby_options(reached)$level$certified)
    expect_true(dagby_options(reached)$alternatives_complete)
    expect_level_holds_sets(graph, reached)
  }
})

test_that("a query with one size says so", {
  .cf_cache$clear()
  single <- heard(dagby_construct(ranked_query(single_edges), size = "largest"))
  expect_identical(single$messages, "DAG with 3 variables. This is the only size.")
  .cf_cache$clear()
  wide <- heard(dagby_construct(ranked_query(wide_edges), size = "largest"))
  expect_identical(wide$messages, "DAG with 9 variables. This is the largest size. Next smaller size: 8 variables. Use dagby_prune(dag, 1).")
  expect_identical(heard(dagby_prune(wide$value))$messages, "DAG with 7 variables. This is the smallest size. Next larger size: 8 variables. Use dagby_grow(dag, 1).")
  expect_identical(heard(dagby_grow(wide$value))$messages, "No larger DAG is available. The DAG was returned unchanged.")
})

test_that("the printed lines give the set size kept for each effect", {
  lines <- function(graph, size) {
    capture.output(invisible(suppressMessages(suppressWarnings(dagby_construct(graph, size, quiet = FALSE)))), type = "output")[2:4]
  }
  .cf_cache$clear()
  expect_identical(lines(ranked_query(mixed_edges), 0), c("DAG with 5 variables",
    "direct: every minimal adjustment set of at most 3 members (3 quantities)",
    "total: every minimal adjustment set of at most 2 members (2 quantities)"))
  expect_identical(lines(ranked_query(fork_edges), 0), c("DAG with 3 variables",
    "direct: every minimal adjustment set of at most 1 member (1 quantity)",
    "total: every minimal adjustment set of at most 1 member (1 quantity)"))
  withr::local_options(causaliflower.ranked_steps = 10)
  .cf_cache$clear()
  expect_identical(lines(ranked_query(ladder_edges), 0), c("DAG with 4 variables",
    "direct: the smallest adjustment sets found (2 quantities)",
    "total: the smallest adjustment sets found (2 quantities)"))
  # Each effect follows its own certificate.
  .cf_cache$clear()
  expect_identical(lines(ranked_query(mediator_edges), 0), c("DAG with 6 variables",
    "direct: no adjustment set found within the step budget",
    "total: every minimal adjustment set of at most 1 member (1 quantity)"))
})

test_that("the same scientific candidate has the same record across construction and resizing", {
  for (edges in c(ladder_edges, tied_edges, mixed_edges)) {
    graph <- ranked_query(edges)
    .cf_cache$clear()
    small <- dagby_construct(graph, 0, quiet = TRUE)
    sizes <- dagby_options(small)$sizes$variables
    for (size in c(0, sizes, Inf)) {
      fresh <- suppressWarnings(dagby_construct(graph, size, quiet = TRUE))
      via_dag <- suppressWarnings(dagby_construct(small, size, quiet = TRUE))
      expect_identical(.gd_signature(via_dag), .gd_signature(fresh), info = paste(edges, size))
      expect_identical(dagby_options(via_dag), dagby_options(fresh), info = paste(edges, size))
    }
    gap <- sizes[[2L]] - sizes[[1L]]
    round_trip <- dagby_prune(dagby_grow(small, gap, quiet = TRUE), gap, quiet = TRUE)
    expect_identical(.gd_signature(round_trip), .gd_signature(small), info = edges)
    expect_identical(dagby_options(round_trip)$candidate, dagby_options(small)$candidate, info = edges)
    .cf_cache$clear()
    expect_identical(dagby_options(dagby_construct(graph, 0, quiet = TRUE)), dagby_options(small))
  }
  withr::local_options(causaliflower.ranked_steps = 5000)
  .cf_cache$clear()
  graph <- ranked_query(ladder_edges)
  expect_identical(dagby_options(dagby_construct(graph, 0, quiet = TRUE))$arguments$steps, 5000L)
  expect_identical(dagby_options(dagby_prune(dagby_construct(graph, quiet = TRUE, size = "largest"), 3, quiet = TRUE))$arguments$steps, 5000L)
})
test_that("the amount asked is named when it is not the amount changed", {
  .cf_cache$clear()
  graph <- ranked_query(ladder_edges)
  top <- dagby_construct(graph, quiet = TRUE, size = "largest")
  small <- dagby_construct(graph, 0, quiet = TRUE)
  expect_identical(heard(dagby_prune(top, 4))$messages,
    "DAG with 5 variables. This is the smallest size. 3 variables were removed, not 4. Next larger size: 7 variables. Use dagby_grow(dag, 2).")
  expect_identical(heard(dagby_grow(small, 1))$messages,
    "DAG with 7 variables. 2 variables were added, not 1. Next smaller size: 5 variables. Use dagby_prune(dag, 2). Next larger size: 8 variables. Use dagby_grow(dag, 1).")
  # No sentence when the amount is met or the amount is Inf.
  expect_identical(heard(dagby_prune(top, 1))$messages, "DAG with 7 variables. Next smaller size: 5 variables. Use dagby_prune(dag, 2). Next larger size: 8 variables. Use dagby_grow(dag, 1).")
  expect_identical(heard(dagby_grow(small, 3))$messages, "DAG with 8 variables. This is the largest size. Next smaller size: 7 variables. Use dagby_prune(dag, 1).")
  expect_identical(heard(dagby_grow(small))$messages, "DAG with 8 variables. This is the largest size. Next smaller size: 7 variables. Use dagby_prune(dag, 1).")
  # One node is singular.
  expect_identical(.cd_size_message(small, "prune", 2, 6L, larger = 7L), "DAG with 5 variables. This is the smallest size. 1 variable was removed, not 2. Next larger size: 7 variables. Use dagby_grow(dag, 2).")
  expect_identical(.cd_size_message(top, "grow", 2, 7L, smaller = 7L), "DAG with 8 variables. This is the largest size. 1 variable was added, not 2. Next smaller size: 7 variables. Use dagby_prune(dag, 1).")
})


test_that("old drawn-node size tables cannot override the exact source-variable count", {
  graph <- ranked_query(tied_edges)
  middle <- dagby_construct(graph, 6, quiet = TRUE)
  carried <- middle
  attr(carried, "cf_provenance")$construction$sizes <- data.frame(nodes = c(4L, 500L), cap = c(2, Inf), certified = TRUE)
  .cf_cache$clear()
  expect_identical(.gd_signature(dagby_select(carried, quiet = TRUE)), .gd_signature(middle))
  expect_identical(dagby_options(dagby_select(carried, quiet = TRUE))$level$variables, 6L)
  expect_identical(.gd_signature(dagby_construct(carried, 9, quiet = TRUE)),
    .gd_signature(dagby_construct(graph, 9, quiet = TRUE)))
  expect_identical(dagby_options(dagby_construct(carried, 9, quiet = TRUE))$sizes$variables, c(4L, 6L, 9L))
})
test_that("a smallest DAG that is not certified can be pruned at another budget", {
  graph <- ranked_query(tied_edges)
  .cf_cache$clear()
  built <- withr::with_options(list(causaliflower.ranked_steps = 3), suppressWarnings(dagby_construct(graph, 6, quiet = TRUE)))
  expect_identical(length(names(built)), 7L)
  expect_identical(dagby_options(built)$level[c("certified", "smallest")], list(certified = FALSE, smallest = TRUE))
  same <- withr::with_options(list(causaliflower.ranked_steps = 3), heard(dagby_prune(built, 1)))
  expect_identical(same$messages, "No smaller DAG is available. The DAG was returned unchanged.")
  expect_identical(same$value, built)
  .cf_cache$clear()
  other <- heard(dagby_prune(built, 1))
  expect_identical(length(names(other$value)), 6L)
  expect_identical(other$messages, "DAG with 6 variables. Next smaller size: 4 variables. Use dagby_prune(dag, 2). Next larger size: 9 variables. Use dagby_grow(dag, 3).")
})

test_that("a very large size is not a tie", {
  graph <- ranked_query(tied_edges)
  .cf_cache$clear()
  top <- heard(dagby_construct(graph, 1e300))
  expect_identical(length(names(top$value)), 9L)
  expect_identical(top$warnings, character())
  small <- dagby_construct(graph, 0, quiet = TRUE)
  pruned <- heard(dagby_prune(top$value, 1e300))
  expect_identical(pruned$warnings, character())
  expect_identical(pruned$messages, "DAG with 4 variables. This is the smallest size. 5 variables were removed, not 1e+300. Next larger size: 6 variables. Use dagby_grow(dag, 2).")
  grown <- heard(dagby_grow(small, 1e300))
  expect_identical(grown$warnings, character())
  expect_identical(grown$messages, "DAG with 9 variables. This is the largest size. 5 variables were added, not 1e+300. Next smaller size: 6 variables. Use dagby_prune(dag, 3).")
  expect_identical(heard(dagby_prune(top$value, 100000))$messages,
    "DAG with 4 variables. This is the smallest size. 5 variables were removed, not 100000. Next larger size: 6 variables. Use dagby_grow(dag, 2).")
  expect_identical(heard(dagby_construct(graph, 5))$warnings,
    "Sizes 4 and 6 are equally close to 5 variables. Returned 6. Use dagby_prune(dag, 2) for 4.")
})

test_that("requests give the same DAGs with the cache off", {
  graph <- ranked_query(tied_edges)
  requests <- function() {
    small <- dagby_construct(graph, 0, quiet = TRUE)
    top <- dagby_construct(graph, quiet = TRUE, size = "largest")
    five <- heard(dagby_construct(graph, 5))
    list(small = small, five = five$value, six = dagby_construct(graph, 6, quiet = TRUE), nine = dagby_construct(graph, 9, quiet = TRUE),
      top = top, pruned = suppressMessages(dagby_prune(top, 3)), grown = suppressMessages(dagby_grow(small, 2)),
      warnings = five$warnings)
  }
  .cf_cache$clear()
  on <- requests()
  withr::local_options(causaliflower.accelerate = FALSE)
  .cf_cache$clear()
  off <- requests()
  expect_identical(.cf_cache$stats()$entries, 0L)
  expect_identical(off$warnings, "Sizes 4 and 6 are equally close to 5 variables. Returned 6. Use dagby_prune(dag, 2) for 4.")
  expect_identical(off$warnings, on$warnings)
  for (name in setdiff(names(on), "warnings")) {
    expect_identical(.gd_signature(off[[name]]), .gd_signature(on[[name]]), info = name)
  }
  expect_identical(vapply(off[setdiff(names(off), "warnings")], function(dag) length(names(dag)), integer(1)),
    c(small = 4L, five = 6L, six = 6L, nine = 9L, top = 9L, pruned = 6L, grown = 6L))
})

test_that("the session cache says whether an entry exists without computing it", {
  .cf_cache$clear()
  expect_false(.cf_cache$has("levels", list("entry")))
  expect_identical(.cf_cache$call("levels", list("entry"), function() 1L), 1L)
  before <- .cf_cache$stats()$counts
  expect_true(.cf_cache$has("levels", list("entry")))
  expect_false(.cf_cache$has("levels", list("other")))
  expect_false(.cf_cache$has("metadata", list("entry")))
  expect_false(.cf_cache$has("levels", list(new.env())))
  expect_identical(.cf_cache$stats()$counts, before)
  withr::local_options(causaliflower.accelerate = FALSE)
  expect_false(.cf_cache$has("levels", list("entry")))
  .cf_cache$clear()
})

test_that("the saved call of every kind of DAG replays in a cleared cache", {
  graph <- ranked_query(tied_edges)
  .cf_cache$clear()
  built <- list(smallest = dagby_construct(graph, 0, quiet = TRUE), largest = dagby_construct(graph, quiet = TRUE, size = "largest"))
  built$middle <- dagby_construct(graph, 6, quiet = TRUE)
  expect_identical(dagby_options(built$largest)$search_method, "maximal-quantities-v1")
  for (name in names(built)) {
    .cf_cache$clear()
    expect_identical(suppressMessages(eval(dagby_call(built[[name]]))), built[[name]], info = name)
  }
})

test_that("latents and keep_nodes apply at every size", {
  graph <- ranked_query(ladder_edges)
  for (size in c(0, 7)) {
    plain <- dagby_construct(graph, size, quiet = TRUE)
    with_latent <- suppressWarnings(dagby_construct(graph, size, latents = "A", quiet = TRUE))
    expect_setequal(names(with_latent), names(plain))
    expect_setequal(edge_keys(with_latent), edge_keys(plain))
    expect_identical(dagby_options(with_latent)$quantities, dagby_options(plain)$quantities)
    expect_identical(dagby_options(with_latent)$arguments$latents, "A")
  }
  # A is not kept at the smallest size, so declaring it marks nothing there.
  small <- dagby_construct(graph, 0, latents = "A", quiet = TRUE)
  expect_false("A" %in% names(small))
  expect_identical(as.character(dagitty::latents(small)), character())
  expect_identical(dagby_options(small)$unmeasured$display_nodes, character())
  large <- suppressWarnings(dagby_construct(graph, 7, latents = "A", quiet = TRUE))
  expect_identical(as.character(dagitty::latents(large)), "A")
  expect_identical(dagby_options(large)$unmeasured$display_nodes, "A")
  expect_identical(as.character(dagitty::latents(dagby_temporal(large))), "A")
  # A requested node is kept at every size and shifts the node count by one.
  for (size in c(0, 7)) {
    kept_result <- heard(dagby_construct(graph, size, keep_nodes = "W", quiet = TRUE))
    expect_identical(kept_result$warnings, if (size == 7) {
      "Sizes 6 and 8 are equally close to 7 variables. Returned 8. Use dagby_prune(dag, 2) for 6."
    } else {
      character()
    })
    kept <- kept_result$value
    expect_true("W" %in% names(kept))
    expect_true("W -> T" %in% edge_keys(kept))
    expect_identical(dagby_options(kept)$arguments$keep_nodes, "W")
    expect_identical(length(names(kept)), length(names(dagby_construct(graph, size, quiet = TRUE))) + 1L)
  }
  expect_identical(length(names(dagby_construct(graph, 5, keep_nodes = "W", quiet = TRUE))), 6L)
  expect_identical(dagby_options(dagby_construct(graph, 5, keep_nodes = "W", quiet = TRUE))$level$nodes, 6L)
  # Resizing keeps both declarations.
  both <- suppressWarnings(dagby_construct(graph, Inf, latents = "A", keep_nodes = "W", quiet = TRUE))
  down <- suppressMessages(dagby_prune(both, 2))
  expect_identical(dagby_options(down)$arguments[c("latents", "keep_nodes")], list(latents = "A", keep_nodes = "W"))
  expect_true("W" %in% names(down))
  expect_identical(length(names(down)), 6L)
  expect_identical(.gd_signature(suppressWarnings(suppressMessages(dagby_grow(down, 2)))), .gd_signature(both))
  expect_error(dagby_construct(graph, 0, latents = "T", quiet = TRUE), "cannot be declared latent")
})

test_that("a quantity every set needs but a declared latent is reported at the size it is kept", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y; B -> A; B -> Y")
  result <- heard(dagby_construct(graph, 0, latents = "A", quiet = TRUE))
  expect_identical(result$warnings, c(
    "Every direct adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A).",
    "Every total adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: A)."))
  expect_identical(dagby_checks(result$value)$total$unmeasured, "A")
})

test_that("a time summary that merges occasions is applied at each size and checked on the source", {
  graph <- ranked_query(repeated_edges)
  small <- dagby_construct(graph, 0, quiet = TRUE)
  top <- suppressMessages(dagby_construct(graph, Inf, quiet = TRUE))
  expect_setequal(names(small), c("T", "Y", "C", "P", "R"))
  expect_setequal(names(dagby_temporal(small)), c("T", "Y", "C_t1", "C_t2", "P", "R"))
  expect_setequal(names(top), c("T", "Y", "A", "B", "C", "P", "R"))
  for (dag in list(small, top)) {
    record <- dagby_options(dag)$summary
    expect_true(record$applied && record$merges_occasions)
    expect_identical(record$displayed_sets_checked,
      length(.gd_adjustment_sets(dag, "direct")$sets) + length(.gd_adjustment_sets(dag, "total")$sets))
    expect_identical(unname(vapply(dagby_checks(dag), `[[`, character(1), "status")), c("preserved", "preserved"))
  }
  expect_identical(set_keys(.gd_adjustment_sets(small, "total")$sets), c("C+P", "C+R"))
  expect_identical(set_keys(.gd_adjustment_sets(dagby_temporal(small), "total")$sets), c("C_t1+P", "C_t1+R"))
  # The time-unrolled graph holds the source's sets for the same quantities.
  source_sets <- set_keys(.gd_adjustment_sets(graph, "total")$sets)
  expect_true(all(c("C_t1+P", "C_t1+R") %in% source_sets))
  expect_identical(length(names(dagby_construct(graph, 5, quiet = TRUE))), 5L)
})

test_that("the sizes of a query without a direct set keep the mediating node", {
  graph <- ranked_query("T -> M; M -> Y; M <-> Y; C -> T; C -> Y; A -> T; A -> D; B -> T; B -> D; D -> Y")
  expect_identical(.gd_adjustment_sets(graph, "direct")$sets, list())
  for (size in c(0, 7)) {
    dag <- suppressWarnings(dagby_construct(graph, size, quiet = TRUE))
    expect_identical(length(names(dag)), if (size) 7L else 5L)
    expect_identical(dagby_options(dag)$quantities$direct$state, "no_sets")
    expect_identical(dagby_options(dag)$quantities$total$state, "sets")
    expect_identical(dagby_options(dag)$direct_path_nodes, "M")
    expect_true("M" %in% names(dag))
    expect_identical(.gd_adjustment_sets(dag, "direct")$sets, list())
    expect_identical(dagby_checks(dag)$direct$status, "no eligible reference set")
    expect_identical(set_keys(.gd_adjustment_sets(dag, "total")$sets), if (size) c("A+B+C", "C+D") else "C+D")
  }
  warned <- heard(dagby_construct(graph, 0, quiet = TRUE))$warnings
  expect_identical(warned, "No eligible direct adjustment set was found. Returning a diagnostic graph. This is not a valid empty set.")
})

test_that("no set, the empty set and a set stay distinct at the smallest size", {
  empty <- dagby_construct(ranked_query("T -> Y; A -> Y"), 0, quiet = TRUE)
  expect_identical(dagby_options(empty)$quantities$total$state, "empty_set")
  expect_identical(dagby_options(empty)$quantities$direct$state, "empty_set")
  expect_setequal(names(empty), c("T", "Y"))
  expect_identical(.gd_adjustment_sets(empty, "total")$sets, list(character()))
  expect_identical(dagby_options(empty)$level[c("smallest", "top")], list(smallest = TRUE, top = FALSE))
  seen <- heard(dagby_construct(ranked_query("T -> Y; T <-> Y"), 0, quiet = TRUE))
  expect_match(seen$warnings, "No eligible direct/total adjustment set was found", fixed = TRUE)
  expect_identical(dagby_options(seen$value)$quantities$total$state, "no_sets")
  expect_identical(dagby_checks(seen$value)$total$status, "no eligible reference set")
  expect_identical(.gd_adjustment_sets(seen$value, "total")$sets, list())
  expect_true("T <-> Y" %in% edge_keys(seen$value) || "Y <-> T" %in% edge_keys(seen$value))
  mixed <- dagby_construct(ranked_query("T -> Y; A -> T; A -> Y"), 0, quiet = TRUE)
  expect_identical(dagby_options(mixed)$quantities$total$state, "sets")
  expect_setequal(names(mixed), c("T", "Y", "A"))
})

test_that("undecided quantities appear only in the largest DAG", {
  graph <- ranked_query(ladder_edges)
  for (size in c(0, 7)) {
    kept <- dagby_options(dagby_construct(graph, size, quiet = TRUE))$quantities
    expect_identical(kept$direct$undecided, character())
    expect_identical(kept$total$undecided, character())
  }
  # The cached list would answer the largest without the member search.
  .cf_cache$clear()
  testthat::local_mocked_bindings(.cd_maximal_steps = function(steps = NULL) 1L, .package = "causaliflower")
  seen <- heard(dagby_construct(graph, quiet = TRUE, size = "largest"))
  top <- dagby_options(seen$value)$quantities
  expect_setequal(top$total$undecided, c("A", "B", "C", "P", "R"))
  expect_identical(top$total$members, character())
  expect_identical(seen$warnings, c(
    "5 direct quantities were not decided within the step budget (1) and were retained. This is not evidence that a set uses them.",
    "5 total quantities were not decided within the step budget (1) and were retained. This is not evidence that a set uses them."))
  # A ranked budget too small to find a set leaves the smallest size without one.
  old <- options(causaliflower.ranked_steps = 1)
  on.exit(options(old), add = TRUE)
  small <- heard(dagby_construct(graph, 0, quiet = TRUE))
  quantities <- dagby_options(small$value)$quantities
  expect_identical(quantities$total$state, "not_found")
  expect_identical(quantities$total$undecided, character())
  expect_identical(small$warnings, "No direct/total adjustment set was found within the step budget. This is not evidence that no set exists.")
  expect_identical(dagby_checks(small$value)$total$status, "not fully checked")
  printed <- capture.output(invisible(suppressMessages(suppressWarnings(dagby_construct(graph, 0, quiet = FALSE)))), type = "output")
  expect_true("direct: no adjustment set found within the step budget" %in% printed)
  expect_true("total: no adjustment set found within the step budget" %in% printed)
})


test_that("no set, the empty set and a set stay distinct", {
  empty <- dagby_construct(ranked_query("T -> Y; A -> Y"), quiet = TRUE, size = "largest")
  expect_identical(dagby_options(empty)$quantities$total$state, "empty_set")
  expect_identical(dagby_options(empty)$quantities$direct$state, "empty_set")
  # A is an outside cause, kept by rule 5.
  expect_setequal(names(empty), c("T", "Y", "A"))
  expect_false(dagby_options(empty)$level$smallest)
  expect_identical(heard(dagby_construct(ranked_query("T -> Y; A -> Y"), size = "largest"))$messages,
    "DAG with 3 variables. This is the largest size. Next smaller size: 2 variables. Use dagby_prune(dag, 1).")
  expect_setequal(names(suppressMessages(dagby_prune(empty, 1))), c("T", "Y"))
  expect_identical(.gd_adjustment_sets(empty, "total")$sets, list(character()))
  seen <- character()
  absent <- withCallingHandlers(dagby_construct(ranked_query("T -> Y; T <-> Y"), quiet = TRUE, size = "largest"), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_match(seen, "No eligible direct/total adjustment set was found", fixed = TRUE)
  expect_identical(dagby_options(absent)$quantities$total$state, "no_sets")
  expect_identical(dagby_checks(absent)$total$status, "no eligible reference set")
  expect_identical(.gd_adjustment_sets(absent, "total")$sets, list())
  expect_true("T <-> Y" %in% edge_keys(absent) || "Y <-> T" %in% edge_keys(absent))
})

test_that("a query without a direct set keeps its mediating node and gains no direct set", {
  graph <- ranked_query("T -> M; M -> Y; M <-> Y")
  expect_identical(.gd_adjustment_sets(graph, "direct")$sets, list())
  dag <- suppressWarnings(dagby_construct(graph, quiet = TRUE, size = "largest"))
  expect_identical(dagby_options(dag)$quantities$direct$state, "no_sets")
  expect_identical(dagby_options(dag)$direct_path_nodes, "M")
  expect_true("M" %in% names(dag))
  expect_identical(.gd_adjustment_sets(dag, "direct")$sets, list())
  expect_identical(dagby_options(dag)$quantities$total$state, "empty_set")
})

test_that("exact treatment parents are retained when no total set exists", {
  graph <- dagitty::dagitty("dag { X1 [exposure] X2 [exposure] Y [outcome]
    X1 -> L -> X2 -> M -> Y; L -> M }")
  dag <- suppressWarnings(dagby_construct(graph, quiet = TRUE, size = "largest"))
  guard <- dagby_options(dag)$identification
  expect_true(guard$applied)
  expect_identical(guard$parents, "L")
  expect_identical(guard$status, "preserved by retaining exact treatment parents")
  expect_true("L" %in% names(dag))
  expect_false(dagby_options(dag)$summary$applied)
  expect_match(dagby_options(dag)$level$static_reason, "joint total identification safeguard")
  expect_length(dagitty::adjustmentSets(dag, effect = "total"), 0L)
})

test_that("declared unmeasured variables are marked afterwards and change nothing else", {
  graph <- ranked_query(wide_edges)
  observed <- dagby_construct(graph, quiet = TRUE, size = "largest")
  marked <- dagby_construct(graph, latents = "A", quiet = TRUE, size = "largest")
  expect_identical(as.character(dagitty::latents(marked)), "A")
  expect_identical(dagby_options(marked)$quantities, dagby_options(observed)$quantities)
  expect_setequal(names(marked), names(observed))
  expect_setequal(edge_keys(marked), edge_keys(observed))
  record <- dagby_options(marked)$unmeasured
  expect_identical(record$semantics, "observed-as-if-v1")
  expect_identical(record$display_nodes, "A")
  expect_identical(as.character(dagitty::latents(dagby_temporal(marked))), "A")
  expect_identical(dagby_call(marked)$arguments$latents, "A")
  expect_error(dagby_construct(graph, latents = "T", quiet = TRUE, size = "largest"), "cannot be declared latent")
})

test_that("keep_nodes are retained", {
  graph <- ranked_query(wide_edges)
  dag <- dagby_construct(graph, keep_nodes = "E", quiet = TRUE, size = "largest")
  expect_true("E" %in% names(dag))
  expect_true("E -> A" %in% edge_keys(dag))
  expect_identical(dagby_options(dag)$quantities, dagby_options(dagby_construct(graph, quiet = TRUE, size = "largest"))$quantities)
})

test_that("a quantity left undecided is retained and reported", {
  graph <- ranked_query(wide_edges)
  seen <- character()
  .cf_cache$clear()
  testthat::local_mocked_bindings(.cd_maximal_steps = function(steps = NULL) 1L, .package = "causaliflower")
  limited <- withCallingHandlers(dagby_construct(graph, quiet = TRUE, size = "largest"), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("direct quantities were not decided within the step budget (1) and were retained. This is not evidence that a set uses them.", seen, fixed = TRUE)))
  kept <- dagby_options(limited)$quantities
  expect_true(length(kept$direct$undecided) > 0L)
  testthat::local_mocked_bindings(.cd_maximal_steps = function(steps = NULL) 2000L, .package = "causaliflower")
  full <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_true(all(names(full) %in% names(limited)))
  expect_true(all(dagby_options(full)$quantities$direct$members %in% c(kept$direct$members, kept$direct$undecided)))
})

test_that("a dry-run plan is accepted for resizing", {
  graph <- ranked_query(ladder_edges)
  top <- dagby_construct(graph, quiet = TRUE, size = "largest")
  small <- dagby_construct(graph, 0, quiet = TRUE)
  plan <- dagby_construct(graph, quiet = TRUE, dry_run = TRUE, size = "largest")
  expect_identical(suppressMessages(dagby_grow(plan)), top)
  expect_identical(.gd_signature(suppressMessages(dagby_prune(plan))), .gd_signature(small))
})


test_that("the largest route has its own step budget and refuses adjusted nodes", {
  graph <- ranked_query(wide_edges)
  expect_identical(.cd_maximal_steps(), 256L)
  # The ranked search's option does not apply to this route.
  old <- options(causaliflower.ranked_steps = "a")
  on.exit(options(old), add = TRUE)
  expect_identical(dagby_options(dagby_construct(graph, quiet = TRUE, size = "largest"))$arguments$steps, 256L)
  adjusted <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A [adjusted]; A -> T; A -> Y; T -> Y }")
  for (size in c(0, 5, Inf)) {
    expect_error(dagby_construct(adjusted, size, quiet = TRUE),
      "The DAG has adjusted or selected nodes. Remove those tags first.", fixed = TRUE)
  }
})

test_that("counts of one are worded in the singular", {
  .cf_cache$clear()
  single <- ranked_query(single_edges)
  printed <- capture.output(invisible(suppressMessages(dagby_construct(single, quiet = FALSE, size = "largest"))), type = "output")
  expect_true(any(printed == "total: 1 quantity used by minimal adjustment sets"))
  expect_identical(sum(printed == "Largest DAG"), 1L)
  plan <- capture.output(invisible(suppressMessages(dagby_construct(single, quiet = FALSE, dry_run = TRUE, size = "largest"))), type = "output")
  expect_identical(sum(plan == "Largest DAG"), 1L)
  level <- capture.output(invisible(suppressMessages(dagby_construct(ranked_query(ladder_edges), 0, quiet = FALSE))), type = "output")
  expect_identical(sum(level == "DAG with 5 variables"), 1L)
  expect_identical(sum(level == "Largest DAG"), 0L)
  # Found by random search: with a budget of two steps only Q6 stays undecided.
  hard <- ranked_query("Q1_t1 -> Q2_t1; Q1_t1 -> Q3_t2; Q1_t1 -> Y; Q1_t1 <-> Q5_t0; Q2_t1 -> Q4_t1; Q3_t0 -> Q1_t1; Q3_t0 -> Q3_t2; Q3_t0 -> Q4_t1; Q3_t0 -> T; Q3_t0 -> Y; Q4_t1 -> Q3_t2; Q4_t1 <-> Q6_t1; Q5_t0 -> Q1_t1; Q5_t0 -> Q2_t1; Q5_t0 -> Q3_t2; Q5_t0 -> Q4_t1; Q5_t0 -> Q6_t1; Q5_t0 -> T; Q6_t1 -> Q3_t2; Q6_t1 -> Q4_t1; Q6_t1 -> Y; T -> Q1_t1; T -> Q4_t1; T -> Y; Y -> Q3_t2")
  .cf_cache$clear()
  testthat::local_mocked_bindings(.cd_maximal_steps = function(steps = NULL) 2L, .package = "causaliflower")
  limited <- heard(dagby_construct(hard, quiet = TRUE, size = "largest"))
  expect_identical(dagby_options(limited$value)$quantities$direct$undecided, "Q6")
  expect_identical(dagby_options(limited$value)$quantities$total$undecided, character())
  expect_identical(limited$warnings, "1 direct quantity was not decided within the step budget (2) and was retained. This is not evidence that a set uses it.")
  expect_true("Q6_t1" %in% names(dagby_temporal(limited$value)))
})

test_that("dry_run returns a plan for a size", {
  graph <- ranked_query(wide_edges)
  plan <- dagby_construct(graph, dry_run = TRUE, quiet = TRUE, size = "largest")
  expect_s3_class(plan, "cf_dag_plan")
  expect_identical(plan$graph, dagby_construct(graph, quiet = TRUE, size = "largest"))
  expect_identical(unname(plan$status), c("sets", "sets"))
  expect_length(plan$changes$removed_nodes, 0L)
  printed <- capture.output(print(plan))
  expect_true(any(grepl("Largest DAG", printed, fixed = TRUE)))
  expect_true(any(grepl("total: 4 quantities used by minimal adjustment sets", printed, fixed = TRUE)))
  expect_true(any(grepl("Final checks: direct preserved, total preserved", printed, fixed = TRUE)))
  expect_identical(dagby_call(plan)[[1L]], quote(causaliflower:::.cd_replay))
  expect_identical(dagby_source(plan), graph)
  ladder <- ranked_query(ladder_edges)
  small <- dagby_construct(ladder, 0, dry_run = TRUE, quiet = TRUE)
  expect_s3_class(small, "cf_dag_plan")
  expect_identical(small$graph, dagby_construct(ladder, 0, quiet = TRUE))
  expect_setequal(small$changes$removed_nodes, c("A", "B", "Q", "W"))
})

# The public query stage executes this fixture in place of the source model
mock_query <- function(fixture, env = parent.frame()) {
  attr(fixture, "query_record")$arguments$params <- dagby_query_params()
  testthat::local_mocked_bindings(.dagby_query_source = function(crop, columns = NULL, treatments = NULL,
      outcomes = NULL, randomised_treatment = FALSE, time_unroll = "all", label = NULL, quiet = NULL,
      params = dagby_query_params()) {
    fixture
  }, .package = "causaliflower", .env = env)
  invisible(fixture)
}

test_that("construction inherits recipe sizes and preserves snapped-size guidance", {
  .cf_cache$clear()
  fixture <- mock_query(ranked_query(ladder_edges))
  query <- dagby("barley", "T", "Y", TRUE, size = "largest") |> dagby_query(quiet = TRUE)
  seen <- heard(dagby_construct(query))
  expected <- heard(dagby_construct(fixture, size = "largest"))
  expect_identical(seen$messages, expected$messages)
  expect_identical(seen$warnings, expected$warnings)
  dag <- seen$value
  expect_identical(.gd_signature(dag), .gd_signature(expected$value))
  expect_identical(length(names(dag)), 8L)
  expect_identical(dagby_options(dag)$level$top, TRUE)
  call <- dagby_call(dag)
  expect_identical(call[[1L]], quote(causaliflower:::.cd_replay))
  expect_identical(call$arguments$size, 8L)
  expect_identical(call$graph, query)
  expect_identical(.gd_signature(suppressMessages(eval(call))), .gd_signature(dag))
  sized <- heard(dagby("barley", "T", "Y", TRUE, size = 6) |>
    dagby_query(quiet = TRUE) |> dagby_construct())
  expect_identical(sized$warnings, "Sizes 5 and 7 are equally close to 6 variables. Returned 7. Use dagby_prune(dag, 2) for 5.")
  expect_identical(length(names(sized$value)), 7L)
  expect_identical(dagby_call(sized$value)$arguments$size, 7L)
  expect_identical(.gd_signature(suppressWarnings(suppressMessages(eval(dagby_call(sized$value))))), .gd_signature(sized$value))
  expect_identical(length(names(dagby("barley", "T", "Y", TRUE, size = Inf) |>
    dagby_query(quiet = TRUE) |> dagby_construct(quiet = TRUE))), 8L)
  expect_identical(length(names(dagby("barley", "T", "Y", TRUE, size = 0) |>
    dagby_query(quiet = TRUE) |> dagby_construct(quiet = TRUE))), 5L)
  expect_error(dagby("barley", "T", "Y", TRUE, size = -1),
    "size must be one number of at least 0, or Inf.", fixed = TRUE)
  expect_error(dagby("barley", "T", "Y", TRUE, quiet = TRUE), "unused argument")
  # Resizing preserves the recipe and its public query call
  smaller <- suppressMessages(dagby_prune(dag))
  expect_identical(length(names(smaller)), 5L)
  expect_identical(.gd_signature(suppressMessages(dagby_grow(smaller))), .gd_signature(dag))
  expect_identical(.gd_state(smaller)$input_record, .gd_state(dag)$input_record)
})


test_that("the exported functions take the specified arguments", {
  expect_identical(names(formals(dagby_construct)), c("graph", "size", "latents", "keep_nodes", "dry_run", "quiet", "unroll"))
  expect_identical(as.list(formals(dagby_construct))[c("size", "latents", "keep_nodes", "dry_run", "quiet")],
    list(size = "smallest", latents = NULL, keep_nodes = NULL, dry_run = FALSE, quiet = NULL))
  expect_identical(names(formals(dagby_prune)), c("dag", "size", "quiet", "unroll"))
  expect_identical(names(formals(dagby_grow)), c("dag", "size", "quiet", "unroll"))
  expect_identical(names(formals(dagby_select)), c("dag", "size", "select", "unroll", "quiet"))
  expect_identical(formals(dagby_construct)$unroll, FALSE)
  expect_null(formals(dagby_prune)$unroll)
  expect_null(formals(dagby_grow)$unroll)
  expect_identical(formals(dagby_prune)$size, Inf)
  expect_identical(formals(dagby_grow)$size, Inf)
  expect_null(formals(dagby_prune)$quiet)
  expect_null(formals(dagby_grow)$quiet)
  expect_identical(names(formals(dagby)), c("crop", "treatments", "outcomes", "randomised_treatment",
    "size", "keep_nodes", "latents", "columns", "time_unroll", "label", "unroll", "params"))
  expect_identical(names(formals(dagby_query)), c("recipe", "quiet"))
  graph <- ranked_query(ladder_edges)
  for (argument in list(list(search = "auto"), list(steps = 5), list(direct_select = c(1, 1)),
    list(total_select = c(1, 1)), list(witness = TRUE))) {
    expect_error(do.call(dagby_construct, c(list(graph, quiet = TRUE), argument)), "unused argument")
  }
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_error(dagby_prune(dag, steps = 5), "unused argument")
  expect_error(dagby_grow(dag, steps = 5), "unused argument")
})

test_that("dagby refuses the search and selector inputs", {
  mock_query(ranked_query(ladder_edges))
  refused <- list(search = "auto", steps = 5, direct_select = c(1, 1), total_select = c(1, 1),
    max_sets = 5)
  for (name in names(refused)) {
    expect_error(do.call(dagby, c(list("barley", "T", "Y", TRUE),
      stats::setNames(list(refused[[name]]), name))), "unused argument")
  }
  expect_error(dagby("barley", "T", "Y", TRUE, search = "exact", steps = 5), "unused argument")
  expect_error(dagby("barley", "T", "Y", TRUE, search = "witness"), "unused argument")
})

test_that("checks of an edited DAG are not current", {
  dag <- dagby_construct(ranked_query(mixed_edges), quiet = TRUE, size = "largest")
  expect_no_warning(checks <- dagby_checks(dag))
  expect_true(attr(checks, "current"))
  edited <- set_latent(dag, "C")
  expect_warning(checks <- dagby_checks(edited),
    "These checks describe the DAG before it was edited. They do not verify the edited DAG.", fixed = TRUE)
  expect_false(attr(checks, "current"))
})

test_that("a set the temporal graph lists that the query rejects marks the checks changed", {
  listing <- .gd_adjustment_sets
  testthat::local_mocked_bindings(.gd_adjustment_sets = function(dag, effect, max_sets = Inf) {
    found <- listing(dag, effect, max_sets)
    if (identical(effect, "total") && identical(max_sets, 50L)) found$sets <- list(character())
    found
  }, .package = "causaliflower")
  seen <- heard(dagby_construct(ranked_query(single_edges), quiet = TRUE, size = "largest"))
  check <- dagby_checks(seen$value)$total
  expect_false(check$listed_valid_on_source)
  expect_identical(check$status, "changed")
  expect_true("The temporal graph changed the total adjustment sets. Inspect dagby_checks() and dagby_temporal()." %in% seen$warnings)
  expect_identical(dagby_checks(seen$value)$direct$status, "preserved")
})

test_that("a listed set where the query has none is refused and reported", {
  listing <- .gd_adjustment_sets
  testthat::local_mocked_bindings(.gd_adjustment_sets = function(dag, effect, max_sets = Inf) {
    found <- listing(dag, effect, max_sets)
    if (identical(effect, "total") && identical(max_sets, 50L)) found$sets <- list(character())
    found
  }, .package = "causaliflower")
  seen <- heard(dagby_construct(ranked_query("T -> Y; T <-> Y"), quiet = TRUE, size = "largest"))
  expect_identical(dagby_options(seen$value)$quantities$total$state, "no_sets")
  expect_false(dagby_options(seen$value)$summary$applied)
  expect_identical(dagby_options(seen$value)$summary$reason,
    "The source has no total adjustment set, but the time summary lists one. The temporal display was retained.")
  expect_false(dagby_checks(seen$value)$total$no_set_kept)
  expect_true("The source has no total adjustment set, but the temporal graph lists one. The graph is not reliable for this effect." %in% seen$warnings)
})


test_that("the rules keep the route, dead ends, outside causes and relevant parentless variables", {
  rules <- function(text) .cd_rule_nodes(dagitty::dagitty(paste("dag {", text, "}")))
  # Two treatments and two outcomes, a dead end after an outcome and a parentless ancestor.
  several <- rules(paste("T1 [exposure]; T2 [exposure]; Y1 [outcome]; Y2 [outcome]; T1 -> M -> T2; T2 -> N -> K -> L -> Y1;",
    "Y1 -> Y2; T1 -> D; D -> E; Y2 -> F; P -> T1; P -> Q; Q -> N; R -> Q"))
  expected <- list(c("L", "M", "N"), c("K", "L", "M", "N"), c("D", "K", "L", "M", "N"), c("D", "E", "F", "K", "L", "M", "N"),
    c("D", "E", "F", "K", "L", "M", "N", "P", "Q"), c("D", "E", "F", "K", "L", "M", "N", "P", "Q", "R"))
  expect_identical(lapply(several, sort), expected)
  # No route, only a direct edge; a <-> edge carries no ancestry.
  direct <- rules("T [exposure]; Y [outcome]; T -> Y; T -> D; A -> T; A <-> Y; B <-> T; C -> Y")
  expect_identical(lapply(direct, sort), list(character(), character(), "D", "D", c("A", "C", "D"), c("A", "C", "D")))
  # B's only link upstream is a <-> edge, so it has no parent.
  upstream <- rules("T [exposure]; Y [outcome]; A <-> B; B -> C; C -> T; T -> Y")
  expect_identical(lapply(upstream, sort), list(character(), character(), character(), character(), "C", c("B", "C")))
})

test_that("equal rules are one row and the sizes are nested on their pinned nodes", {
  .cf_cache$clear()
  graph <- ranked_query(wide_edges)
  rows <- .cd_rule_rows(.cd_size_query(graph, NULL), TRUE)$rows
  expect_identical(vapply(rows, `[[`, integer(1), "rule"), c(0L, 4L, 6L))
  expect_identical(vapply(rows, `[[`, integer(1), "variables"), c(7L, 8L, 9L))
  expect_identical(vapply(rows, `[[`, logical(1), "top"), c(FALSE, FALSE, TRUE))
  expect_identical(lapply(rows, `[[`, "extra"), list(character(), "Z", c("E", "Z")))
  # dagby_prune() and dagby_grow() step through the rule sizes.
  walk <- list(dagby_construct(graph, quiet = TRUE, size = "largest"))
  for (step in 1:2) {
    walk[[step + 1L]] <- suppressMessages(dagby_prune(walk[[step]], 1))
  }
  expect_identical(vapply(walk, function(dag) length(names(dag)), integer(1)), c(9L, 8L, 7L))
  expect_identical(vapply(walk, function(dag) dagby_options(dag)$level$rule, integer(1)), c(6L, 4L, 0L))
  for (step in 1:2) {
    expect_true(all(names(dagby_temporal(walk[[step + 1L]])) %in% names(dagby_temporal(walk[[step]]))))
  }
  source <- .cd_prepare(graph, NULL, NULL)$selection
  for (dag in walk) {
    for (effect in c("direct", "total")) {
      expect_true(all(.gd_validate_sets(source, .gd_adjustment_sets(dag, effect)$sets, effect)))
    }
  }
  expect_identical(.gd_signature(suppressMessages(dagby_grow(suppressMessages(dagby_grow(walk[[3L]], 1)), 1))), .gd_signature(walk[[1L]]))
  printed <- capture.output(invisible(suppressMessages(dagby_construct(graph, quiet = FALSE, size = "largest"))), type = "output")
  expect_true("also kept: 1 dead end, 1 upstream node" %in% printed)
  printed <- capture.output(.cd_print_level(.cd_construction(walk[[2L]])))
  expect_identical(printed[c(1L, length(printed))], c("DAG with 8 variables", "also kept: 1 dead end"))
  printed <- capture.output(.cd_print_level(.cd_construction(walk[[3L]])))
  expect_identical(printed[[1L]], "DAG with 7 variables")
  expect_false(any(grepl("also kept", printed)))
})

test_that("rule sizes count source variables independently of the checked display", {
  graph <- ranked_query("A_t1 -> T; A_t1 -> A_t2; A_t2 -> Y; B -> A_t2; B -> Y; C_t1 -> T; C_t1 -> C_t2; C_t2 -> Y; D -> Y; T -> Y")
  rows <- .cd_rule_rows(.cd_size_query(graph, NULL), TRUE)$rows
  expect_identical(vapply(rows, `[[`, integer(1), "variables"), c(4L, 6L))
  for (i in seq_along(rows)) {
    static <- rule_dag(graph, i)
    temporal <- rule_dag(graph, i, unroll = TRUE)
    map <- attr(graph, "query_record")$node_map
    counted <- length(unique(map$member_id[map$node %in% names(temporal)]))
    expect_identical(counted, rows[[i]]$variables)
    expect_identical(dagby_options(static)$level$variables, counted)
    expect_identical(dagby_options(temporal)$level$variables, counted)
    expect_identical(length(names(temporal)), c(6L, 8L)[[i]])
    expect_identical(.gd_signature(dagby_temporal(static)), .gd_signature(temporal))
  }
  six <- dagby_construct(graph, 6, quiet = TRUE)
  expect_identical(dagby_options(six)$level$rule, 5L)
  expect_identical(dagby_options(six)$sizes, data.frame(variables = c(4L, 6L)))
  .cf_cache$clear()
  four <- dagby_prune(six, 1, quiet = TRUE)
  expect_identical(dagby_options(four)$level[c("rule", "variables")], list(rule = 0L, variables = 4L))
  expect_identical(.gd_signature(dagby_grow(four, 1, quiet = TRUE)), .gd_signature(six))
})
test_that("each distinct rule candidate remains selectable at its exact variable count", {
  graph <- ranked_query(paste("Q1_t2 -> Q4_t2; Q2_t0 -> Q2_t2; Q2_t0 -> Q4_t2; Q2_t2 -> Q4_t2; Q3_t0 -> Q1_t2; Q3_t0 -> Q2_t2;",
    "Q3_t0 -> Q4_t2; Q3_t0 -> Q6_t2; Q3_t0 -> Q7_t0; Q5_t1 -> Y; Q6_t1 -> Q6_t2; Q6_t1 -> Q7_t1; Q7_t0 -> Q1_t2;",
    "Q7_t0 -> Q2_t2; Q7_t0 -> Q6_t2; Q7_t0 -> Q7_t1; Q7_t0 -> Y; Q8_t1 -> Q2_t2; Q8_t1 -> Q6_t1; Q8_t1 -> Q7_t1;",
    "T -> Q6_t2; T -> Q7_t1; T -> Q8_t1; Y -> Q1_t2; Y -> Q2_t2; Y -> Q4_t2; Y -> Q6_t2"))
  equal <- ranked_query(paste("Q1_t1 -> Q3_t1; Q1_t1 -> T_t1; Q1_t1 -> Y; Q3_t0 -> Q4_t1; Q3_t0 -> Q5_t0; Q3_t0 -> T_t1;",
    "Q4_t0 -> Q1_t1; Q4_t0 -> Q3_t0; Q4_t0 -> Q3_t1; Q4_t0 -> Q4_t1; Q4_t0 -> T_t0; Q4_t0 -> Y; Q5_t0 -> Q2_t2;",
    "Q5_t0 -> Q3_t1; Q5_t0 -> Q4_t1; Q5_t0 <-> T_t1; T_t0 -> Q1_t1; T_t0 -> Q2_t2; T_t0 -> Q3_t0; T_t0 -> Q3_t1;",
    "T_t0 -> T_t1; T_t0 -> Y; T_t1 -> Y"), exposures = c("T_t0", "T_t1"))

  for (fixture in list(graph, equal)) {
    query <- .cd_size_query(fixture, NULL)
    rules <- .cd_rule_rows(query, TRUE)
    map <- attr(fixture, "query_record")$node_map
    for (row in rules$rows) {
      built <- suppressWarnings(.cd_maximal(fixture, NULL, NULL, FALSE, TRUE,
        list(unroll = TRUE), query, rules, row))
      counted <- length(unique(map$member_id[map$node %in% names(built)]))
      expect_identical(row$variables, counted)
      public <- suppressWarnings(dagby_construct(fixture, counted, unroll = TRUE, quiet = TRUE))
      choices <- dagby_options(public)$alternatives
      signatures <- lapply(choices$select, function(index) {
        selected <- suppressWarnings(dagby_select(public, select = index, quiet = TRUE))
        expect_identical(dagby_options(selected)$level$variables, counted)
        expect_identical(length(unique(map$member_id[map$node %in% names(selected)])), counted)
        .gd_signature(dagby_temporal(selected))
      })
      expect_true(any(vapply(signatures, identical, logical(1), y = .gd_signature(dagby_temporal(built)))))
      expect_false(anyDuplicated(signatures) > 0L)
      .cf_cache$clear()
      cold <- suppressWarnings(dagby_construct(fixture, counted, unroll = TRUE, quiet = TRUE))
      expect_identical(dagby_options(cold)$alternatives, choices)
      expect_identical(.gd_signature(cold), .gd_signature(public))
    }
  }
})
test_that("a ranked candidate identical to the full rule has one stable public identity", {
  graph <- ranked_query(tied_edges)
  .cf_cache$clear()
  listed <- dagby_construct(graph, 9, quiet = TRUE)
  .cf_cache$clear()
  full <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_identical(.gd_signature(listed), .gd_signature(full))
  expect_identical(dagby_options(listed)$candidate, dagby_options(full)$candidate)
  expect_identical(dagby_options(listed)$alternatives, dagby_options(full)$alternatives)
  expect_identical(nrow(dagby_options(full)$alternatives), 1L)
  expect_identical(dagby_options(full)$level[c("cap", "rule", "top")], list(cap = Inf, rule = 0L, top = TRUE))
})
test_that("the size returned does not depend on the searches kept in the session", {
  # A size at or below the smallest level returns that level; above it the rule sizes take part whenever the last level could be returned.
  fixtures <- list(
    # One level of 3 nodes, as many as rule 0; rule 3 adds the dead end D.
    single = list(graph = ranked_query("A -> T; A -> Y; T -> Y; T -> D"), sizes = c(3, 3),
      rows = list(c(1, 0), c(1, 0), c(1, 0), c(1, 0), c(Inf, 3))),
    tied = list(graph = ranked_query(tied_edges), sizes = c(4, 9),
      rows = list(c(2, 0), c(2, 0), c(Inf, 0), c(2, 0), c(3, 0))))
  for (name in names(fixtures)) {
    graph <- fixtures[[name]]$graph
    sizes <- fixtures[[name]]$sizes
    query <- .cd_size_query(graph, NULL)
    .cf_cache$clear()
    top <- dagby_construct(graph, quiet = TRUE, size = "largest")
    smallest <- dagby_construct(graph, 0, quiet = TRUE)
    for (kept in list(character(), "members", "list", c("members", "list"))) {
      rows <- lapply(1:5, function(request) {
        .cf_cache$clear()
        if ("members" %in% kept) {
          invisible(.cd_rule_rows(query, TRUE))
        }
        if ("list" %in% kept) {
          definitions <- .cd_level_definitions(query, .cd_ranked_steps())
        }
        dag <- switch(request, dagby_construct(graph, 0, quiet = TRUE), dagby_construct(graph, sizes[[1L]], quiet = TRUE),
          dagby_construct(graph, sizes[[2L]], quiet = TRUE), dagby_prune(top, Inf, quiet = TRUE), dagby_grow(smallest, 1, quiet = TRUE))
        as.numeric(unlist(dagby_options(dag)$level[c("cap", "rule")]))
      })
      expect_identical(rows, fixtures[[name]]$rows, info = paste(name, paste(kept, collapse = " ")))
    }
  }
})

test_that("declared unmeasured variables and requested nodes apply at the rule sizes", {
  graph <- ranked_query(wide_edges)
  marked <- dagby_construct(graph, latents = "E", quiet = TRUE, size = "largest")
  expect_identical(as.character(dagitty::latents(marked)), "E")
  expect_identical(dagby_checks(marked)$total$unmeasured, character())
  # W is requested, so no rule adds a node and the largest DAG is rule 0.
  kept <- dagby_construct(ranked_query(ladder_edges), keep_nodes = "W", quiet = TRUE, size = "largest")
  expect_identical(dagby_options(kept)$level[c("rule", "top")], list(rule = 0L, top = TRUE))
  expect_null(dagby_options(kept)$kept_by_rule)
})

test_that("a multi-member query keeps supported ranked choices and consistently refuses the full rule", {
  graph <- ranked_query("T -> M; M -> Y; C -> T; C -> Y; D -> C; D -> Y; T -> Y")
  record <- attr(graph, "query_record")
  map <- record$node_map
  map <- rbind(map, map[map$node == "C", ])
  map$member_id[map$node == "C"] <- c("C1", "C2")
  record$node_map <- map
  attr(graph, "query_record") <- record
  for (accelerate in c(TRUE, FALSE)) {
    .cf_cache$clear()
    withr::with_options(list(causaliflower.accelerate = accelerate), {
      expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "stands for several source members")
      small <- dagby_construct(graph, 0, quiet = TRUE)
      expect_identical(dagby_options(small)$level$variables, 5L)
      expect_setequal(names(small), c("T", "Y", "M", "C"))
      expect_identical(dagby_options(small)$sizes$variables, 5L)
      expect_false(dagby_options(small)$level$top)
      expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "stands for several source members")
      expect_identical(.gd_signature(dagby_construct(graph, 5, quiet = TRUE)), .gd_signature(small))
      expect_identical(.gd_signature(dagby_construct(graph, 100, quiet = TRUE)), .gd_signature(small))
      expect_identical(dagby_grow(small, 1, quiet = TRUE), small)
    })
  }
})

test_that("a largest DAG saved before the rule sizes grows to them", {
  graph <- ranked_query(wide_edges)
  .cf_cache$clear()
  old <- rule_dag(graph)
  state <- attr(old, "cf_provenance")
  state$construction$level <- list(cap = Inf, certified = TRUE, smallest = TRUE, top = TRUE, nodes = 7L)
  state$construction$sizes <- data.frame(nodes = 7L, cap = Inf, certified = TRUE)
  attr(old, "cf_provenance") <- state
  .cf_cache$clear()
  grown <- suppressMessages(dagby_grow(old, 1))
  expect_identical(dagby_options(grown)$level[c("rule", "nodes")], list(rule = 4L, nodes = 8L))
  expect_identical(.gd_signature(grown), .gd_signature(suppressMessages(dagby_grow(rule_dag(graph), 1))))
  expect_setequal(names(suppressMessages(dagby_grow(old, Inf))), c("T", "Y", "A", "B", "C", "D", "M", "E", "Z"))
  # Its finite sizes are carried.
  graph <- ranked_query(ladder_edges)
  .cf_cache$clear()
  old <- rule_dag(graph)
  state <- attr(old, "cf_provenance")
  state$construction$level <- list(cap = Inf, certified = TRUE, smallest = FALSE, top = TRUE, nodes = 7L)
  state$construction$sizes <- data.frame(nodes = c(5L, 7L), cap = c(2, Inf), certified = TRUE)
  attr(old, "cf_provenance") <- state
  .cf_cache$clear()
  expect_identical(dagby_options(suppressMessages(dagby_grow(old, Inf)))$sizes,
    data.frame(variables = c(5L, 7L, 8L)))
})

test_that("unexpected member-search errors are not converted into a complete choice list", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y; T -> D")
  .cf_cache$clear()
  testthat::local_mocked_bindings(.cd_maximal_setup = function(...) stop("other"), .package = "causaliflower")
  for (size in c(0, 3, Inf)) {
    expect_error(dagby_construct(graph, size, quiet = TRUE), "other", fixed = TRUE)
  }
})
test_that("missing saved sizes are recomputed and never change resize choices", {
  graph <- ranked_query("A_t1 -> T; A_t1 -> A_t2; A_t2 -> Y; B -> A_t2; B -> Y; C_t1 -> T; C_t1 -> C_t2; C_t2 -> Y; D -> Y; T -> Y")
  six <- dagby_construct(graph, 6, quiet = TRUE)
  without <- six
  attr(without, "cf_provenance")$construction$sizes <- NULL
  .cf_cache$clear()
  four <- dagby_prune(without, 1, quiet = TRUE)
  expect_identical(dagby_options(four)$level$variables, 4L)
  expect_identical(.gd_signature(four), .gd_signature(dagby_prune(six, 1, quiet = TRUE)))
  expect_identical(.gd_signature(dagby_grow(four, 1, quiet = TRUE)), .gd_signature(six))
  full <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_identical(dagby_grow(full, 1, quiet = TRUE), full)
})
test_that("a larger size whose time summary gives no more nodes returns the DAG unchanged", {
  # Without a direct set M_t1 is kept; the top adds M_t2, which the time summary merges with it.
  graph <- ranked_query("T -> M_t1; M_t1 -> M_t2; M_t2 -> Y; M_t1 <-> Y")
  .cf_cache$clear()
  zero <- suppressWarnings(rule_dag(graph))
  expect_setequal(names(zero), c("T", "M", "Y"))
  grown <- heard(dagby_grow(zero, Inf))
  expect_identical(grown$value, zero)
  expect_identical(grown$messages, "No larger DAG is available. The DAG was returned unchanged.")
})
