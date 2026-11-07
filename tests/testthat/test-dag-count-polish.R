count_polish_graph <- function(edges) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  ids <- .strip_slice(names(graph))
  attr(graph, "query_record") <- list(
    arguments = list(crop = "barley", params = list(water_model = "SoilWater", query_depth = Inf)),
    node_map = data.frame(node = names(graph), quantity = ids, member_id = ids,
                          occurrence = names(graph)),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes")
  )
  graph
}

test_that("count receipts identify distinct nodes and explain the bounded prefix", {
  graph <- count_polish_graph("T -> A -> B -> Y; T -> Y")
  search <- .gd_adjustment_sets
  testthat::local_mocked_bindings(.gd_adjustment_sets = function(dag, effect, max_sets, ...)
    search(dag, effect, 1L), .package = "causaliflower")
  expect_warning(
    options <- construct_dag(graph, direct_select = NULL, quiet = TRUE, search = "complete"),
    "sets with exact occasions"
  )
  printed <- capture.output(print(options))
  expect_true(any(grepl("time repeats count once", printed, fixed = TRUE)))
  expect_true(any(grepl("distinct node", printed, fixed = TRUE)))
  expect_false(any(grepl("quantities", printed, fixed = TRUE)))
  expect_false(any(unlist(options$direct$raw_sets) %in% c("T", "Y")))
})

test_that("wrapper records complete construction without a hidden result limit", {
  graph <- count_polish_graph("T -> M -> Y; T -> Y")
  testthat::local_mocked_bindings(
    query_graph = function(crop, treatments, outcomes, randomised_treatment,
                           label = NULL, quiet = NULL, columns = NULL,
                           time_unroll = "all", params = NULL) graph,
    .package = "causaliflower"
  )
  result <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE,
    construct_dag = list(search = "complete"))
  expect_false("max_sets" %in% names(dag_options(result)$arguments))
  expect_false("max_sets" %in% names(dag_call(result)$construct_dag))
  expect_identical(dag_options(result)$search_method, "complete-source-sets-v2")
  expect_true(all(vapply(dag_options(result)$search, function(x) isTRUE(x$complete), logical(1))))
})

test_that("an identical temporal fallback reuses the same completed checks", {
  # Two retained occasions of the outcome quantity cannot be merged, so the
  # summary keeps the temporal display. (Feedback folds such as
  # T -> A_t0 -> B_t1 -> A_t2 -> Y or A_t0 -> T -> A_t1 are split into time
  # segments instead. See test-time-segments.R.)
  graph <- count_polish_graph("T -> A -> Y; Y_t0 -> A; T -> Y")
  original <- .cd_checks
  calls <- 0L
  testthat::local_mocked_bindings(
    .cd_checks = function(...) {
      calls <<- calls + 1L
      original(...)
    },
    .package = "causaliflower"
  )
  result <- construct_dag(graph, keep_nodes = "Y_t0", quiet = TRUE)
  expect_false(.gd_state(result)$summary$applied)
  expect_match(.gd_state(result)$summary$reason, "distinct occasions of an endpoint")
  expect_identical(.gd_state(result)$checks, .gd_state(result)$temporal_checks)
  expect_identical(calls, 1L)
  # One retained occasion of A declared unmeasured no longer refuses the fold as
  # mixed observation status: construction merges A_t0 and A_t1 as for a fully
  # observed graph, then marks A and records its observed occasion.
  calls <- 0L
  mixed <- suppressWarnings(construct_dag(count_polish_graph("T -> A_t0 -> Y; A_t1 -> A_t0; T -> Y"),
                                          latents = "A_t1", keep_nodes = "A_t1", quiet = TRUE))
  expect_true(.gd_state(mixed)$summary$applied)
  expect_identical(dagitty::latents(mixed), "A")
  expect_identical(dagitty::latents(temporal_dag(mixed)), "A_t1")
  expect_identical(dag_options(mixed)$unmeasured$partial, list(A = list(unmeasured = "A_t1", observed = "A_t0")))
  expect_identical(calls, 2L)
})

test_that("old catalogue indexes are not silently reinterpreted after a search change", {
  graph <- count_polish_graph("T -> A -> B -> Y; T -> Y")
  result <- construct_dag(graph, quiet = TRUE)
  for (method in list(NULL, "source-verified-refinement-v1")) {
    state <- .gd_state(result)
    state$construction$search_method <- method
    old <- .gd_attach(result, state)
    expect_identical(original_query(old), graph)
    expect_error(select_dag(old, direct_select = c(1L, 1L)), "predates the current search")
    expect_warning(saved <- dag_call(old), "before replaying")
    expect_true(is.call(saved))
  }
})
# Refinement must find a smaller set valid on the source, not just shrink a picture.
test_that("projection proposals reduce the number of distinct nodes on the original graph", {
  graph <- count_polish_graph("T -> A_t0 -> B_t0 -> Y; T -> A_t1 -> Y; T -> Y")
  prepared <- .cd_prepare(graph, NULL, NULL)
  initial <- setNames(lapply(c("direct", "total"), function(effect)
    .cd_search(prepared, effect, 100L)), c("direct", "total"))
  initial$direct$sets <- list(c("B_t0", "A_t1"))
  initial$direct$table <- .cd_candidates(initial$direct$sets, prepared)
  initial$direct$complete <- FALSE
  before <- serialize(prepared, NULL)
  result <- .cd_refine_search(prepared, initial, 100L)
  expect_identical(serialize(prepared, NULL), before)
  expect_identical(result$search$direct$source_sets, initial$direct$sets)
  expect_equal(result$search$direct$table$n_resolved[1L], 1L)
  expect_setequal(result$search$direct$table$raw_sets[[1L]][[1L]], c("A_t0", "A_t1"))
  expect_false(result$search$direct$complete)
  expect_false(result$refinement$global_minimum_certified)
  expect_lte(length(result$refinement$passes), 3L)
  for (set in result$search$direct$refinement_sets)
    expect_identical(.cd_certify_refinement(prepared$selection, set, "direct")$status, "accepted")
})

test_that("certification on the original reference rejects invalid and redundant proposals", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> M -> Y; Z; U [latent] }")
  before <- serialize(graph, NULL)
  expect_identical(.cd_certify_refinement(graph, "M", "direct")$status, "accepted")
  redundant <- .cd_certify_refinement(graph, c("M", "Z"), "direct")
  expect_true(redundant$sufficient)
  expect_false(redundant$minimal)
  expect_setequal(redundant$witness, "M")
  for (set in list(character(), "T", "U", "absent"))
    expect_identical(.cd_certify_refinement(graph, set, "direct")$status, "rejected")
  expect_identical(.cd_certify_refinement(graph, character(), "total")$status, "accepted")
  expect_identical(serialize(graph, NULL), before)
})

test_that("complete source enumeration is retained without refinement", {
  prepared <- .cd_prepare(count_polish_graph("T -> M -> Y"), NULL, NULL)
  initial <- setNames(lapply(c("direct", "total"), function(effect)
    .cd_search(prepared, effect, 100L)), c("direct", "total"))
  result <- .cd_refine_search(prepared, initial, 100L)
  expect_length(result$refinement$passes, 0L)
  expect_identical(result$search$direct$sets, initial$direct$sets)
  expect_identical(result$search$total$sets, initial$total$sets)
})
