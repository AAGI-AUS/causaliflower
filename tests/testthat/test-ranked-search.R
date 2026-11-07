ranked_query <- function(edges, columns = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
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

test_that("the default route ranks options and select_dag uses the saved ranking", {
  dag <- construct_dag(ranked_query("A -> T; A -> B; B -> Y; T -> M -> Y; T -> Y"), quiet = TRUE)
  options <- dag_options(dag)
  expect_identical(options$search_method, "quantity-ranked-v1")
  expect_identical(options$search$total$ranked$certificate, "complete")
  expect_setequal(dag_set_catalogue(dag)$total$label, c("A", "B"))
  switched <- select_dag(dag, total_select = c(1L, 2L))
  expect_identical(dag_options(switched)$search_method, "quantity-ranked-v1")
  expect_identical(dag_options(switched)$selected$total, c(1L, 2L))
})

test_that("a step budget is reported as a partial ranking, never as no set", {
  old <- options(causaliflower.ranked_steps = 1L)
  on.exit(options(old), add = TRUE)
  graph <- ranked_query("A -> T; A -> B; B -> Y; C -> T; C -> D; D -> Y; T -> Y")
  seen <- character()
  dag <- withCallingHandlers(construct_dag(graph, quiet = TRUE), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("certified complete up to", seen)))
  expect_false(any(grepl("No eligible", seen)))
  search <- dag_options(dag)$search$total
  expect_false(search$complete)
  expect_identical(search$ranked$certificate, "budget")
})

test_that("the witness route is explicit and needs automatic selectors", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y")
  dag <- construct_dag(graph, quiet = TRUE, search = "witness")
  expect_identical(dag_options(dag)$search_method, "source-witness-v1")
  expect_error(construct_dag(graph, direct_select = NULL, search = "witness", quiet = TRUE),
               "requires automatic")
})

test_that("an option that fails the source check is dropped and uncertifies the ranking", {
  graph <- ranked_query("A -> T; A -> B; B -> Y; T -> M -> Y; T -> Y")
  enumerate <- .cd_ranked_enumerate
  testthat::local_mocked_bindings(.cd_ranked_enumerate = function(region, labels, step_budget) {
    result <- enumerate(region, labels, step_budget)
    # Add an invalid option realised as the outcome itself. It must fail the source check.
    result$options <- c(result$options, list(integer()))
    result$realised <- c(result$realised, list(region$graph$nodes[region$y]))
    result
  })
  seen <- character()
  dag <- withCallingHandlers(construct_dag(graph, quiet = TRUE), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("failed validation", seen)))
  expect_true(any(grepl("not certified", seen)))
  expect_false(any(grepl("No eligible", seen)))
  expect_false(any(grepl("certified complete up to", seen)))
  reference <- .cd_prepare(graph, NULL, NULL)$selection
  for (effect in c("direct", "total")) {
    search <- dag_options(dag)$search[[effect]]
    expect_false(search$complete)
    expect_identical(search$ranked$certificate, "validation")
    expect_identical(search$ranked$dropped, 1L)
    expect_identical(search$ranked$complete_up_to, -1L)
    expect_true(search$omitted_equal_smallest)
    expect_true(length(search$sets) > 0L)
    expect_true(all(.gd_validate_sets(reference, search$sets, effect)))
  }
  expect_true(any(grepl("not certified", capture.output(print(dag_options(dag))), fixed = TRUE)))
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
  # Z_t1 mediates T -> Y through Q_t2. Merging drops the edge across repeated stages Z_t1 -> Q_t2
  # and would show Z on no route although it is a mediator on the source.
  graph <- ranked_query("T -> Q_t1; Q_t1 -> Y; Q_t1 -> Z_t1; T -> Z_t1; Z_t1 -> Q_t2; Q_t2 -> Y; Q_t2 -> Z_t2; T -> Y",
                        columns = "Z")
  dag <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  state <- .gd_state(dag)
  expect_false(state$summary$applied)
  expect_match(state$summary$reason, "route from treatment to outcome|represented occurrences")
  expect_true(all(c("Z_t1", "Q_t2") %in% names(dag)))
  # A view of the direct family of the earlier graph must not lose the route T -> Z_t1 -> Q_t2 -> Y either.
  graph2 <- ranked_query("Q_t1 -> T; Q_t1 -> Z_t1; T -> Z_t1; Z_t1 -> Q_t2; Q_t2 -> Y; Q_t2 -> Z_t2; T -> Y",
                         columns = c("Z", "Q"))
  dag2 <- suppressWarnings(construct_dag(graph2, quiet = TRUE))
  view <- dag_sets_view(dag2, list(direct = 1L))
  expect_false(.gd_state(view)$set_view$summary$applied)
  expect_match(.gd_state(view)$set_view$summary$reason, "route from treatment to outcome|represented occurrences")
})

test_that("merging occasions is refused when a merged node's adjustment set is invalid on the source", {
  # Q_t2 mediates T -> Y through Z_t1. Merging Q_t1 and Q_t2 into Q would make the
  # displayed total set {Q} stand for {Q_t1, Q_t2}, which is not valid on the source.
  graph <- ranked_query("Q_t1 -> T; Q_t1 -> Z_t1; T -> Z_t1; Z_t1 -> Q_t2; Q_t2 -> Y; Q_t2 -> Z_t2; T -> Y",
                        columns = "Z")
  reference <- .cd_prepare(graph, NULL, NULL)$selection
  dag <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  state <- .gd_state(dag)
  expect_false(state$summary$applied)
  # Either display guard may refuse first: the route guard or the check of merged sets.
  expect_match(state$summary$reason, "represented occurrences|route from treatment to outcome")
  expect_true(all(c("Q_t1", "Q_t2") %in% names(dag)))
  for (effect in c("direct", "total")) {
    final <- dag_checks(dag)[[effect]]$final_sets
    if (length(final)) expect_true(all(.gd_validate_sets(reference, final, effect)), info = effect)
  }
  view <- dag_sets_view(dag, list(total = 1L))
  expect_false(.gd_state(view)$set_view$summary$applied)
  expect_match(.gd_state(view)$set_view$summary$reason, "represented occurrences|route from treatment to outcome")
  # The merged reading itself is invalid, which is what the refusal protects against.
  expect_false(.gd_validate_sets(reference, list(c("Q_t1", "Q_t2")), "total"))
})

test_that("a budget stop with no option is recorded as unresolved, not as absence", {
  old <- options(causaliflower.ranked_steps = 1L)
  on.exit(options(old), add = TRUE)
  seen <- character()
  dag <- withCallingHandlers(construct_dag(ranked_query("A -> T; A -> Y; T -> Y"), quiet = TRUE),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  options <- dag_options(dag)
  expect_identical(options$status[["total"]], "no_sets")
  expect_false(options$search$total$complete)
  expect_match(options$choices$total$reason, "step budget")
  expect_match(options$identification$reason, "step budget")
  expect_true(any(grepl("not evidence that no set exists", seen)))
  expect_false(any(grepl("No eligible", seen)))
})

test_that("select_dag indexes the saved ranked catalogue even under a different step budget", {
  dag <- construct_dag(ranked_query("A -> T; A -> B; B -> Y; T -> M -> Y; T -> Y"), quiet = TRUE)
  saved <- dag_options(dag)$search$total$sets
  old <- options(causaliflower.ranked_steps = 1L)
  on.exit(options(old), add = TRUE)
  switched <- select_dag(dag, total_select = c(1L, 2L))
  expect_identical(dag_options(switched)$search$total$sets, saved)
  expect_identical(dag_options(switched)$search_identity$ranked_steps, dag_options(dag)$search_identity$ranked_steps)
  expect_identical(dag_options(switched)$selected$total, c(1L, 2L))
})

test_that("the ranked search refuses adjusted or selected nodes and unbounded budgets", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y; S -> T; S -> Y")
  dagitty::adjustedNodes(graph) <- "S"
  expect_error(construct_dag(graph, quiet = TRUE), "adjusted or selected")
  old <- options(causaliflower.ranked_steps = 1e12)
  on.exit(options(old), add = TRUE)
  expect_error(construct_dag(ranked_query("A -> T; A -> Y; T -> Y"), quiet = TRUE), "no larger than")
})

test_that("checking merged sets alone refuses a display whose merged node set is invalid on the source", {
  # Random probe seed 20 (stage-c/review3-random.R): the route guard accepts this
  # collapse. Only reading the merged direct set as all represented occurrences
  # refuses it. Nodes keep the probe's names. Exposure Q2_t0, outcome Q3_t1.
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
  class(graph) <- c("dag_by_query", class(graph))
  latents <- c("S1", "Q3_t0", "Q4_t1")
  # Declared unmeasured variables now change only the explicit complete and
  # witness searches, so the probe's declared refusal is reproduced on the
  # complete route, which applies the declarations before its search.
  warned <- character()
  dag <- withCallingHandlers(construct_dag(graph, latents = latents, quiet = TRUE, search = "complete"),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  state <- .gd_state(dag)
  expect_false(state$summary$applied)
  expect_match(state$summary$reason, "represented occurrences")
  # The direct family has no eligible set on this source (a diagnostic graph is
  # returned). Any set the retained temporal display reports that the source
  # rejects must be flagged "changed" and warned about, never shown silently.
  expect_identical(dag_options(dag)$status[["direct"]], "no_sets")
  reference <- .cd_prepare(graph, latents, NULL)$selection
  for (effect in c("direct", "total")) {
    check <- dag_checks(dag)[[effect]]
    if (!length(check$final_sets)) next
    valid <- .gd_validate_sets(reference, check$final_sets, effect)
    if (all(valid)) next
    expect_identical(check$status, "changed", info = effect)
    expect_true(any(vapply(check$introduced, function(x) identical(x$status, "not valid on reference"), logical(1))), info = effect)
  }
  expect_true(any(grepl("changed or were not fully verified", warned)))
  # The ranked default constructs as for a fully observed graph: the direct
  # family now has an option, the display keeps time for another reason, the
  # same invariant holds on the observed reference, and the declared nodes
  # retained by the selected sets are marked afterwards.
  warned <- character()
  ranked <- withCallingHandlers(construct_dag(graph, latents = latents, quiet = TRUE),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_identical(dag_options(ranked)$status[["direct"]], "selected")
  expect_false(.gd_state(ranked)$summary$applied)
  observed <- .cd_prepare(graph, NULL, NULL)$selection
  rejected <- FALSE
  for (effect in c("direct", "total")) {
    check <- dag_checks(ranked)[[effect]]
    if (!length(check$final_sets)) next
    valid <- .gd_validate_sets(observed, check$final_sets, effect)
    if (all(valid)) next
    rejected <- TRUE
    expect_identical(check$status, "changed", info = effect)
    expect_true(any(vapply(check$introduced, function(x) identical(x$status, "not valid on reference"), logical(1))), info = effect)
  }
  if (rejected) expect_true(any(grepl("changed or were not fully verified", warned)))
  expect_identical(vapply(dag_checks(ranked), `[[`, character(1), "status"), c(direct = "preserved", total = "preserved"))
  expect_setequal(dagitty::latents(ranked), c("Q3_t0", "Q4_t1"))
})
