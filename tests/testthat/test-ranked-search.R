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

ranked_random_graph <- function(seed) {
  set.seed(seed)
  K <- sample(3:6, 1L)
  TT <- sample(2:3, 1L)
  S <- sample(0:2, 1L)
  qn <- paste0("Q", seq_len(K))
  sn <- if (S) paste0("S", seq_len(S)) else character()
  order <- sample(qn)
  at <- function(q, t) paste0(q, "_t", t)
  nodes <- c(sn, unlist(lapply(0:(TT - 1L), function(t) at(qn, t))))
  E <- list()
  for (t in 0:(TT - 1L)) for (i in seq_along(order)) for (j in seq_along(order))
    if (i < j && stats::runif(1) < 0.3) E[[length(E) + 1L]] <- c(at(order[i], t), at(order[j], t))
  for (t in seq_len(TT - 1L) - 1L) for (a in qn) for (b in qn)
    if (stats::runif(1) < if (a == b) 0.7 else 0.2) E[[length(E) + 1L]] <- c(at(a, t), at(b, t + 1L))
  for (s in sn) for (v in setdiff(nodes, sn)) if (stats::runif(1) < 0.15) E[[length(E) + 1L]] <- c(s, v)
  X <- if (S && stats::runif(1) < 0.3) sample(sn, 1L) else at(sample(qn, 1L), 0L)
  Y <- at(sample(setdiff(qn, sub("_t[0-9]+$", "", X)), 1L), TT - 1L)
  if (stats::runif(1) < 0.5) E[[length(E) + 1L]] <- c(X, Y)
  E <- unique(E)
  g <- dagitty::dagitty(paste0("dag { ", paste(nodes, collapse = " "), " ",
    paste(vapply(E, paste, "", collapse = " -> "), collapse = " ; "), " }"))
  dagitty::exposures(g) <- X
  dagitty::outcomes(g) <- Y
  pool <- setdiff(nodes, c(X, Y))
  dagitty::latents(g) <- pool[stats::runif(length(pool)) < 0.1]
  g
}

ranked_quantity <- function(x) sub("_t[0-9]+$", "", x)
ranked_key <- function(s) paste(sort(unique(ranked_quantity(s))), collapse = "+")
ranked_canonical <- function(sets) {
  keys <- unique(vapply(sets, ranked_key, ""))
  members <- lapply(keys, function(k) if (nzchar(k)) strsplit(k, "+", fixed = TRUE)[[1L]] else character())
  keep <- !vapply(seq_along(members), function(i) any(vapply(seq_along(members), function(j)
    j != i && length(members[[j]]) < length(members[[i]]) && all(members[[j]] %in% members[[i]]),
    logical(1))), logical(1))
  sort(keys[keep])
}

test_that("ranked options equal exhaustive minimal sets collapsed to quantities", {
  quantity <- ranked_quantity
  canonical <- ranked_canonical
  for (seed in 1:25) for (effect in c("total", "direct")) {
    g <- ranked_random_graph(seed)
    oracle <- canonical(as.list(dagitty::adjustmentSets(g, type = "minimal", effect = effect)))
    region <- if (effect == "total") .cd_ranked_total_region(g) else .cd_ranked_direct_region(g)
    q <- quantity(region$graph$nodes)
    levels <- unique(q)
    result <- .cd_ranked_enumerate(region, lapply(match(q, levels), identity), 1e6)
    mine <- sort(vapply(result$options, function(set) paste(sort(levels[set]), collapse = "+"), ""))
    expect_identical(result$certificate, "complete", info = paste(seed, effect))
    expect_identical(mine, oracle, info = paste(seed, effect))
    if (length(result$realised))
      expect_true(all(.gd_validate_sets(g, result$realised, effect)), info = paste(seed, effect))
  }
})

test_that("the ranked search ranks options with a complete certificate", {
  prepared <- .cd_prepare(ranked_query("A -> T; A -> B; B -> Y; T -> M -> Y; T -> Y"), NULL, NULL)
  total <- .cd_ranked_search(prepared, "total")
  expect_identical(total$ranked$method, "quantity-ranked-v1")
  expect_identical(total$ranked$certificate, "complete")
  expect_setequal(total$table$label, c("A", "B"))
})

test_that("a step budget is reported as a partial ranking, never as no set", {
  graph <- ranked_query("A -> T; A -> B; B -> Y; C -> T; C -> D; D -> Y; T -> Y")
  search <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), "total", steps = 1L)
  expect_false(search$complete)
  expect_identical(search$ranked$certificate, "budget")
  seen <- character()
  withr::local_options(causaliflower.ranked_steps = 1L)
  dag <- withCallingHandlers(dagby_construct(graph, 0, quiet = TRUE), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("No eligible", seen)))
  expect_false(dagby_options(dag)$level$certified)
})

test_that("an option that fails the source check is dropped and uncertifies the ranking", {
  # C makes the full candidate larger, so size 0 exercises ranked construction
  graph <- ranked_query("A -> T; A -> B; B -> Y; T -> M -> Y; T -> Y; T -> C")
  enumerate <- .cd_ranked_enumerate
  testthat::local_mocked_bindings(.cd_ranked_enumerate = function(region, labels, step_budget) {
    result <- enumerate(region, labels, step_budget)
    result$options <- c(result$options, list(integer()))
    result$realised <- c(result$realised, list(region$graph$nodes[region$y]))
    result
  }, .dagby_query_source = function(...) graph, .package = "causaliflower")
  reference <- .cd_prepare(graph, NULL, NULL)$selection
  for (effect in c("direct", "total")) {
    expect_warning(search <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), effect),
      "failed validation on the selection reference. The ranking is not certified.", fixed = TRUE)
    expect_false(search$complete)
    expect_identical(search$ranked$certificate, "validation")
    expect_identical(search$ranked$dropped, 1L)
    expect_identical(search$ranked$complete_up_to, -1L)
    expect_true(search$omitted_equal_smallest)
    expect_true(length(search$sets) > 0L)
    expect_true(all(.gd_validate_sets(reference, search$sets, effect)))
  }
  seen <- character()
  dag <- withCallingHandlers(dagby("barley", "T", "Y", TRUE, size = 0) |>
    dagby_query(quiet = TRUE) |> dagby_construct(quiet = TRUE), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("No eligible", seen)))
  expect_false(dagby_options(dag)$level$certified)
})

test_that("a partial certificate is exact: every option within k is listed and nothing listed is spurious", {
  partial <- 0L
  for (seed in 1:12) for (effect in c("total", "direct")) for (budget in c(2L, 3L, 5L, 8L)) {
    g <- ranked_random_graph(seed)
    oracle <- ranked_canonical(as.list(dagitty::adjustmentSets(g, type = "minimal", effect = effect)))
    oracle_sizes <- lengths(lapply(oracle, function(k) if (nzchar(k)) strsplit(k, "+", fixed = TRUE)[[1L]] else character()))
    region <- if (effect == "total") .cd_ranked_total_region(g) else .cd_ranked_direct_region(g)
    q <- ranked_quantity(region$graph$nodes)
    levels <- unique(q)
    result <- .cd_ranked_enumerate(region, lapply(match(q, levels), identity), budget)
    mine <- sort(vapply(result$options, function(set) paste(sort(levels[set]), collapse = "+"), ""))
    k <- result$complete_up_to
    expect_true(all(mine %in% oracle), info = paste(seed, effect, budget))
    if (identical(result$certificate, "complete")) expect_identical(mine, oracle, info = paste(seed, effect, budget))
    else {
      partial <- partial + 1L
      expect_true(all(oracle[oracle_sizes <= k] %in% mine), info = paste(seed, effect, budget, "k", k))
    }
    if (length(result$realised))
      expect_true(all(.gd_validate_sets(g, result$realised, effect)), info = paste(seed, effect, budget))
  }
  expect_true(partial > 0L)
})

test_that("merging occasions is refused when the display hides or invents a route from treatment to outcome", {
  graph <- ranked_query("T -> Q_t1; Q_t1 -> Y; Q_t1 -> Z_t1; T -> Z_t1; Z_t1 -> Q_t2; Q_t2 -> Y; Q_t2 -> Z_t2; T -> Y",
                        columns = "Z")
  expect_error(suppressWarnings(dagby_construct(graph, quiet = TRUE, size = "largest")), "unroll = TRUE", fixed = TRUE)
  dag <- suppressWarnings(dagby_construct(graph, quiet = TRUE, unroll = TRUE, size = "largest"))
  record <- dagby_options(dag)$summary
  expect_false(record$applied)
  expect_match(dagby_options(dag)$level$static_reason, "route from treatment to outcome|represented occurrences")
  expect_true(all(c("Z_t1", "Q_t2") %in% names(dag)))
})

test_that("merging occasions is refused when a merged node's adjustment set is invalid on the source", {
  graph <- ranked_query("Q_t1 -> T; Q_t1 -> Z_t1; T -> Z_t1; Z_t1 -> Q_t2; Q_t2 -> Y; Q_t2 -> Z_t2; T -> Y",
                        columns = "Z")
  reference <- .cd_prepare(graph, NULL, NULL)$selection
  expect_error(suppressWarnings(dagby_construct(graph, quiet = TRUE, size = "largest")), "unroll = TRUE", fixed = TRUE)
  dag <- suppressWarnings(dagby_construct(graph, quiet = TRUE, unroll = TRUE, size = "largest"))
  record <- dagby_options(dag)$summary
  expect_false(record$applied)
  expect_match(dagby_options(dag)$level$static_reason, "represented occurrences|route from treatment to outcome")
  expect_true(all(c("Q_t1", "Q_t2") %in% names(dag)))
  for (effect in c("direct", "total")) {
    listed <- .gd_adjustment_sets(dag, effect, Inf)$sets
    expect_gt(length(listed), 0L)
    expect_true(all(.gd_validate_sets(reference, listed, effect)), info = effect)
  }
  expect_false(.gd_validate_sets(reference, list(c("Q_t1", "Q_t2")), "total"))
})

test_that("a budget stop with no option is recorded as unresolved, not as absence", {
  # C keeps the compact result distinct from the full membership result
  graph <- ranked_query("A -> T; A -> Y; T -> Y; T -> C")
  search <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), "total", steps = 1L)
  expect_length(search$sets, 0L)
  expect_false(search$complete)
  withr::local_options(causaliflower.ranked_steps = 1L)
  testthat::local_mocked_bindings(.dagby_query_source = function(...) graph, .package = "causaliflower")
  seen <- character()
  dag <- withCallingHandlers(dagby("barley", "T", "Y", TRUE, size = 0) |>
    dagby_query(quiet = TRUE) |> dagby_construct(quiet = TRUE), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(dagby_options(dag)$quantities$total$state, "not_found")
  expect_true(any(grepl("not evidence that no set exists", seen)))
  expect_false(any(grepl("No eligible", seen)))
})


test_that("the ranked search refuses adjusted or selected nodes and unbounded budgets", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y; S -> T; S -> Y")
  dagitty::adjustedNodes(graph) <- "S"
  expect_error(dagby_construct(graph, 0, quiet = TRUE), "adjusted or selected")
  withr::local_options(causaliflower.ranked_steps = 1e12)
  expect_error(dagby_construct(ranked_query("A -> T; A -> Y; T -> Y"), 0, quiet = TRUE), "no larger than")
})

test_that("a display that keeps time lists only sets valid on the observed source", {
  nodes <- c("S1", "S2", paste0("Q", 1:4, "_t0"), paste0("Q", 1:4, "_t1"))
  edges <- paste("Q4_t0 -> Q2_t0; Q1_t0 -> Q3_t0; Q4_t1 -> Q1_t1; Q2_t1 -> Q3_t1; Q1_t1 -> Q3_t1;",
    "Q1_t0 -> Q3_t1; Q2_t0 -> Q1_t1; Q2_t0 -> Q2_t1; Q2_t0 -> Q3_t1; Q3_t0 -> Q2_t1; Q3_t0 -> Q3_t1;",
    "Q4_t0 -> Q3_t1; Q4_t0 -> Q4_t1; S1 -> Q3_t1; S2 -> Q4_t1")
  graph <- dagitty::dagitty(paste("dag {", paste(nodes, collapse = " "), edges, "}"))
  dagitty::exposures(graph) <- "Q2_t0"
  dagitty::outcomes(graph) <- "Q3_t1"
  ids <- .strip_slice(names(graph))
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = NULL, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = names(graph), quantity = ids, member_id = ids, occurrence = names(graph)),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  latents <- c("S1", "Q3_t0", "Q4_t1")
  expect_error(suppressWarnings(dagby_construct(graph, latents = latents, quiet = TRUE, size = "largest")), "unroll = TRUE", fixed = TRUE)
  ranked <- suppressWarnings(dagby_construct(graph, latents = latents, quiet = TRUE, unroll = TRUE, size = "largest"))
  expect_identical(dagby_options(ranked)$quantities$direct$state, "sets")
  expect_false(dagby_options(ranked)$summary$applied)
  observed <- .cd_prepare(graph, NULL, NULL)$selection
  as_observed <- .gd_plain(ranked)
  dagitty::latents(as_observed) <- character()
  for (effect in c("direct", "total")) {
    listed <- .gd_adjustment_sets(as_observed, effect, Inf)$sets
    expect_gt(length(listed), 0L)
    expect_true(all(.gd_validate_sets(observed, listed, effect)), info = effect)
  }
  expect_identical(vapply(dagby_checks(ranked), `[[`, character(1), "status"), c(direct = "preserved", total = "preserved"))
  expect_setequal(dagitty::latents(ranked), c("Q3_t0", "Q4_t1", "S1"))
})

test_that("checking merged sets alone refuses a display whose merged node set is invalid on the source", {
  # Random graph seed 214 of ranked_random_graph(), without its declared latent: the structural guards
  # accept the merge, and only reading a merged direct set as all its occurrences refuses it.
  graph <- ranked_query(paste("Q1_t0 -> Q1_t1; Q1_t1 -> Q2_t1; Q1_t1 -> Q2_t2; Q1_t1 -> Q5_t2; Q1_t2 -> Q2_t2;",
    "Q2_t0 -> Q1_t1; Q2_t0 -> Q2_t1; Q3_t0 -> Q3_t1; Q3_t0 -> Q4_t0; Q3_t0 -> Q4_t1; Q3_t1 -> Q1_t1;",
    "Q3_t1 -> Q2_t1; Q3_t1 -> Q3_t2; Q3_t1 -> Q4_t1; Q3_t2 -> Q1_t2; Q4_t0 -> Q1_t0; Q4_t0 -> Q1_t2;",
    "Q4_t0 -> Q2_t0; Q4_t0 -> Q5_t1; Q4_t1 -> Q2_t2; Q4_t1 -> Q4_t2; Q5_t0 -> Q2_t0; Q5_t0 -> Q4_t0;",
    "Q5_t1 -> Q1_t1; Q5_t1 -> Q3_t1; Q5_t1 -> Q3_t2; Q5_t1 -> Q5_t2; Q5_t2 -> Q2_t2; Q5_t2 -> Q3_t2;",
    "Q2_t1; Q2_t2; Q4_t2"))
  record <- attr(graph, "query_record")
  graph <- dagitty::dagitty(sub("T \\[exposure\\]; Y \\[outcome\\];", "", as.character(graph)))
  dagitty::exposures(graph) <- "Q4_t0"
  dagitty::outcomes(graph) <- "Q1_t2"
  attr(graph, "query_record") <- record
  expect_error(dagby_construct(graph, 0, quiet = TRUE), "merged nodes are read as all represented occurrences", fixed = TRUE)
  dag <- suppressWarnings(dagby_construct(graph, 0, quiet = TRUE, unroll = TRUE))
  record <- dagby_options(dag)$summary
  expect_false(record$applied)
  expect_identical(dagby_options(dag)$level$static_reason, paste("A temporally merged direct adjustment set is not valid on the source when its",
    "merged nodes are read as all represented occurrences. The temporal display was retained."))
  expect_identical(.gd_signature(dag), .gd_signature(dagby_temporal(dag)))
  temporal <- dagby_temporal(dag)
  open <- .gd_summary(temporal, keep = names(temporal))
  expect_true(open$applied)
  sets <- .gd_adjustment_sets(open$graph, "direct", Inf)$sets
  expect_gt(length(sets), 0L)
  read <- lapply(sets, function(set) as.character(unique(unlist(lapply(set, function(node) {
    occurrences <- open$node_map$node[open$node_map$display == node]
    if (length(occurrences)) occurrences else node
  })))))
  expect_false(all(.gd_validate_sets(graph, read, "direct")))
  listed <- .gd_adjustment_sets(dag, "direct", Inf)$sets
  expect_gt(length(listed), 0L)
  expect_true(all(.gd_validate_sets(graph, listed, "direct")))
  expect_identical(dagby_checks(dag)$direct$status, "preserved")
})

test_that("ranked candidate tables preserve literal names and harmless renaming", {
  for (labels in list(c("P_t0", "P_t1"), c("A", "B"), c("first occasion", "second occasion"))) {
    graph <- dagitty::dagitty(paste0("dag { T [exposure]; Y [outcome]; \"", labels[[1L]],
      "\" -> \"", labels[[2L]], "\"; \"", labels[[2L]], "\" -> T; \"", labels[[2L]], "\" -> Y; T -> Y }"))
    prepared <- .cd_prepare(graph, NULL, NULL)
    for (effect in c("direct", "total")) {
      search <- .cd_ranked_search(prepared, effect)
      expect_identical(search$sets, list(labels[[2L]]))
      expect_identical(search$table$resolved, list(labels[[2L]]))
      expect_identical(search$table$n_resolved, 1L)
    }
    dag <- suppressWarnings(dagby_construct(graph, 0, quiet = TRUE))
    expect_setequal(names(dag), c(labels[[2L]], "T", "Y"))
    expect_identical(dagby_options(dag)$level$variables, 3L)
    for (effect in c("direct", "total")) {
      sets <- .gd_adjustment_sets(dag, effect, Inf)$sets
      expect_identical(sets, list(labels[[2L]]))
      expect_true(all(.gd_validate_sets(graph, sets, effect)))
      expect_identical(dagby_checks(dag)[[effect]]$status, "preserved")
      expect_true(dagby_checks(dag)[[effect]]$sets_minimal)
      expect_true(dagby_checks(dag)[[effect]]$listed_valid_on_source)
    }
  }
})

test_that("ranked options are ordered by member count", {
  graph <- ranked_query("T -> M; M -> Y; T -> Y; C -> T; C -> Y; A -> T; A -> D; B -> T; B -> D; D -> Y")
  table <- .cd_ranked_search(.cd_prepare(graph, NULL, NULL), "total")$table
  expect_identical(table$label, c("C + D", "A + B + C"))
  expect_identical(table$n_resolved, c(2L, 3L))
})
