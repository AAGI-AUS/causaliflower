fast_query <- function(edges, columns = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  nodes <- names(graph)
  ids <- .strip_slice(nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = columns, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = ids, member_id = ids, occurrence = nodes),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  class(graph) <- c("dag_by_query", class(graph))
  graph
}

test_that("explicit witness construction uses exact minimal witnesses without complete source work", {
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  before <- serialize(graph, NULL)
  oracle <- lapply(c("direct", "total"), function(effect)
    .gd_adjustment_sets(.cd_prepare(graph, NULL, NULL)$selection, effect, Inf)$sets)
  testthat::local_mocked_bindings(
    .cd_search = function(...) stop("complete source search called"),
    .cd_refine_search = function(...) stop("source refinement called"))
  result <- construct_dag(graph, quiet = TRUE, search = "witness")
  options <- dag_options(result)
  expect_identical(options$search_method, "source-witness-v1")
  expect_identical(options$enumeration, "witness")
  for (i in seq_along(oracle)) {
    search <- options$search[[i]]
    expect_identical(search$sets, oracle[[i]][1L])
    expect_false(search$complete)
    expect_false(search$exhaustive_requested)
    expect_true(all(.gd_validate_sets(.cd_prepare(graph, NULL, NULL)$selection,
      search$sets, names(options$search)[i])))
  }
  expect_identical(original_query(result), graph)
  expect_identical(serialize(graph, NULL), before)
  expect_true(all(vapply(dag_checks(result), `[[`, logical(1), "output_reference_complete")))
  expect_identical(dag_call(result)$direct_select, "auto")
  expect_identical(dag_call(result)$search, "witness")
  replay <- eval(dag_call(result))
  expect_identical(.gd_signature(replay), .gd_signature(result))
  expect_error(select_dag(result, direct_select = c(1, 1)), "source witnesses")
  catalogue <- dag_set_catalogue(result)
  expect_identical(catalogue$enumeration, "witness")
  expect_true(all(catalogue$incomplete))
})

test_that("exact temporal memberships replace ambiguous compact Output", {
  graph <- fast_query("T -> A_t1 -> Y; T -> A_t2 -> Y; T -> Y")
  result <- construct_dag(graph, quiet = TRUE, search = "witness")
  state <- .gd_state(result)
  expect_setequal(names(result), c("T", "Y", "A_t1", "A_t2"))
  expect_identical(dag_checks(result)$direct$required, list(c("A_t1", "A_t2")))
  expect_false(state$summary$applied)
  expect_match(state$summary$reason, "exact source temporal")
  expect_false(state$summary$attempted$checks$direct$output_reference_complete)
  expect_true(dag_checks(result)$direct$output_reference_complete)
  expect_identical(.gd_signature(temporal_dag(result)), .gd_signature(result))
  view <- dag_sets_view(result, list(direct = 1L, total = 1L))
  expect_setequal(names(view), names(result))
  expect_identical(.gd_state(view)$set_view$source_method, "source-witness-v1")
  expect_true(all(vapply(dag_checks(view), `[[`, logical(1), "output_reference_complete")))
})

test_that("all Output sets require exact occurrence evidence including known sets", {
  graph <- .cd_prepare(fast_query("A -> T; A -> Y; T -> Y"), NULL, NULL)$selection
  altered <- graph
  map <- attr(altered, "cf_node_map")
  row <- map[map$node == "A", ]
  row$occurrence <- "A_t2"
  attr(altered, "cf_node_map") <- rbind(map, row)
  expect_null(.cd_lift_exact("A", altered, graph))
  for (field in c("quantity", "member_id")) {
    changed <- graph
    mapping <- attr(changed, "cf_node_map")
    mapping[mapping$node == "A", field] <- "different-source-member"
    attr(changed, "cf_node_map") <- mapping
    expect_null(.cd_lift_exact("A", changed, graph))
  }
  expect_identical(.cd_lift_exact("A", graph, graph), "A")
  check <- list(status = "preserved", final_sets = list("A"), final_complete = TRUE)
  checked <- .cd_output_source_check(check, graph, altered, "total")
  expect_false(checked$output_reference_complete)
  expect_identical(checked$status, "changed")
  expect_match(checked$output_reference_status, "interpretation unavailable")
  # An invalid introduced empty set cannot pass just because the witness is known.
  checked <- .cd_output_source_check(list(status = "preserved",
    final_sets = list("A", character()), final_complete = TRUE), graph, graph, "total")
  expect_identical(checked$output_reference_valid, c(TRUE, FALSE))
  expect_false(checked$output_reference_complete)
})

test_that("empty and absent families remain distinct and errors never trigger full fallback", {
  empty <- construct_dag(fast_query("T -> Y"), quiet = TRUE, search = "witness")
  expect_identical(dag_options(empty)$search$direct$sets, list(character()))
  expect_true(dag_options(empty)$search$total$complete)
  expect_identical(dag_set_catalogue(empty)$enumeration, "witness")
  absent <- suppressWarnings(construct_dag(fast_query("U -> T; U -> Y; T -> Y"),
    latents = "U", quiet = TRUE, search = "witness"))
  expect_identical(dag_options(absent)$search$total$sets, list())
  expect_identical(dag_options(absent)$status[["total"]], "no_sets")
  expect_true(dag_options(absent)$search$total$complete)
  testthat::local_mocked_bindings(.gd_adjustment_engine = function(...) stop("witness failure"),
    .cd_search = function(...) stop("unexpected full fallback"))
  expect_error(construct_dag(fast_query("T -> Y"), quiet = TRUE, search = "witness"), "witness failure")
})

test_that("explicit and deferred construction preserve the complete source route and plans record their route", {
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  complete <- construct_dag(graph, quiet = TRUE, search = "complete")
  options <- construct_dag(graph, quiet = TRUE, direct_select = NULL, search = "complete")
  explicit <- construct_dag(graph, quiet = TRUE, direct_select = c(1, 1), search = "complete")
  witness_plan <- construct_dag(graph, quiet = TRUE, dry_run = TRUE, search = "witness")
  plan <- construct_dag(graph, quiet = TRUE, dry_run = TRUE)
  expect_identical(dag_options(complete)$search_method, "complete-source-sets-v2")
  expect_true(all(vapply(dag_options(complete)$search, `[[`, logical(1), "complete")))
  expect_s3_class(options, "cf_dag_options")
  expect_identical(options$search_method, "complete-source-sets-v2")
  expect_identical(options$direct, dag_options(complete)$direct)
  expect_identical(dag_options(explicit)$search_method, "complete-source-sets-v2")
  expect_identical(witness_plan$options$search_method, "source-witness-v1")
  expect_identical(dag_call(witness_plan)$search, "witness")
  expect_identical(plan$options$search_method, "quantity-ranked-v1")
  expect_identical(dag_call(plan)$search, "auto")
})

test_that("get_dag records the ranked default and explicit witness and complete modes", {
  testthat::local_mocked_bindings(query_graph = function(crop, treatments, outcomes,
      randomised_treatment, label = NULL, quiet = NULL, params = query_params())
    fast_query("A -> T; A -> B -> Y; T -> Y"))
  result <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE)
  expect_identical(dag_options(result)$search_method, "quantity-ranked-v1")
  expect_identical(dag_call(result)$construct_dag$search, "auto")
  # A ranked result records the completed numeric selection, as complete did.
  expect_identical(dag_call(result)$construct_dag$direct_select, c(1L, 1L))
  expect_identical(dag_call(result)$construct_dag$direct_select, .gd_state(result)$selected$direct)
  expect_identical(dag_call(result)$construct_dag$total_select, .gd_state(result)$selected$total)
  witness <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE,
    construct_dag = list(search = "witness"))
  expect_identical(dag_options(witness)$search_method, "source-witness-v1")
  expect_identical(dag_call(witness)$construct_dag$direct_select, "auto")
  expect_identical(dag_call(witness)$construct_dag$total_select, "auto")
  expect_identical(dag_call(witness)$construct_dag$search, "witness")
  full <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE,
    construct_dag = list(search = "complete"))
  expect_identical(dag_options(full)$search_method, "complete-source-sets-v2")
  expect_identical(dag_call(full)$construct_dag$search, "complete")
})

test_that("witness reconstruction reuses only matching intact witness receipts", {
  calls <- 0L
  witness_search <- .cd_witness_search
  testthat::local_mocked_bindings(.cd_witness_search = function(...) {
    calls <<- calls + 1L
    witness_search(...)
  })
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  base <- construct_dag(graph, quiet = TRUE, search = "witness")
  changed <- fast_query("A -> T; A -> B -> Y; T -> Y", "B")
  fresh <- construct_dag(changed, quiet = TRUE, search = "witness")
  calls <- 0L
  reused <- .cd_reconstruct(changed, previous = base, quiet = TRUE, search = "witness")
  expect_identical(calls, 0L)
  expect_identical(reused, fresh)
  expect_true(all(vapply(dag_checks(reused), `[[`, logical(1), "output_reference_complete")))
  state <- .gd_state(base)
  state$construction$search_identity <- NULL
  calls <- 0L
  .cd_reconstruct(changed, previous = .gd_attach(base, state), quiet = TRUE, search = "witness")
  expect_identical(calls, 2L)
  calls <- 0L
  .cd_reconstruct(changed, previous = base, latents = "A", quiet = TRUE, search = "witness")
  expect_identical(calls, 2L)
  withr::local_options(causaliflower.accelerate = FALSE)
  calls <- 0L
  .cd_reconstruct(changed, previous = base, quiet = TRUE, search = "witness")
  expect_identical(calls, 2L)
})

# The ranked default (quantity-ranked-v1) has its own coverage below.

test_that("the ranked default collapses repeated occasions and keeps the exact source membership", {
  graph <- fast_query("T -> A_t1 -> Y; T -> A_t2 -> Y; T -> Y")
  result <- construct_dag(graph, quiet = TRUE)
  options <- dag_options(result)
  expect_identical(options$search_method, "quantity-ranked-v1")
  expect_identical(options$search$direct$ranked$certificate, "complete")
  expect_setequal(names(result), c("T", "Y", "A"))
  expect_true(.gd_state(result)$summary$applied)
  checks <- dag_checks(result)
  expect_identical(checks$direct$status, "preserved")
  expect_identical(checks$direct$required, list(c("A_t1", "A_t2")))
  expect_identical(checks$direct$final_sets, list("A"))
  expect_identical(options$search$direct$sets, list(c("A_t1", "A_t2")))
  expect_true(all(.gd_validate_sets(.cd_prepare(graph, NULL, NULL)$selection,
    options$search$direct$sets, "direct")))
})

test_that("ranked reconstruction reuses only matching intact ranked receipts", {
  calls <- 0L
  ranked_search <- .cd_ranked_search
  testthat::local_mocked_bindings(.cd_ranked_search = function(...) {
    calls <<- calls + 1L
    ranked_search(...)
  })
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  base <- construct_dag(graph, quiet = TRUE)
  expect_identical(dag_options(base)$search_method, "quantity-ranked-v1")
  changed <- fast_query("A -> T; A -> B -> Y; T -> Y", "B")
  fresh <- construct_dag(changed, quiet = TRUE)
  calls <- 0L
  reused <- .cd_reconstruct(changed, previous = base, quiet = TRUE)
  expect_identical(calls, 0L)
  expect_identical(reused, fresh)
  state <- .gd_state(base)
  state$construction$search_identity <- NULL
  calls <- 0L
  .cd_reconstruct(changed, previous = .gd_attach(base, state), quiet = TRUE)
  expect_identical(calls, 2L)
  # Declared unmeasured variables do not change the ranked search, so a
  # declaration alone reuses it (the result is marked afterwards).
  calls <- 0L
  marked <- suppressWarnings(.cd_reconstruct(changed, previous = base, latents = "A", quiet = TRUE))
  expect_identical(calls, 0L)
  expect_identical(dagitty::latents(marked), "A")
  withr::local_options(causaliflower.accelerate = FALSE)
  calls <- 0L
  .cd_reconstruct(changed, previous = base, quiet = TRUE)
  expect_identical(calls, 2L)
})

test_that("a ranked search error never triggers the complete fallback", {
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  testthat::local_mocked_bindings(.cd_ranked_search = function(...) stop("ranked failure"),
    .cd_search = function(...) stop("unexpected full fallback"))
  expect_error(construct_dag(graph, quiet = TRUE), "ranked failure")
})

test_that("a ranked result records search 'auto' and its saved call replays", {
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  result <- construct_dag(graph, quiet = TRUE)
  expect_identical(dag_call(result)$search, "auto")
  replay <- eval(dag_call(result))
  expect_identical(dag_options(replay)$search_method, "quantity-ranked-v1")
  expect_identical(.gd_signature(replay), .gd_signature(result))
})

test_that("legacy inspection remains scalar safe and incomplete catalogues stay labelled", {
  graph <- fast_query("A -> T; A -> B -> Y; T -> Y")
  full <- construct_dag(graph, quiet = TRUE, search = "complete")
  state <- .gd_state(full)
  state$construction$search_method <- NULL
  legacy <- .gd_attach(full, state)
  expect_warning(call <- dag_call(legacy), "predates")
  expect_type(call, "language")
  state <- .gd_state(full)
  state$build$search$direct$complete <- FALSE
  state$construction$search <- state$build$search
  partial <- .gd_attach(full, state)
  catalogue <- dag_set_catalogue(partial)
  expect_identical(catalogue$enumeration, "incomplete")
  expect_identical(unname(catalogue$complete), c(FALSE, TRUE))
})
