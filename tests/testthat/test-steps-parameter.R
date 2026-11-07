# The ranked search's explicit step budget: construct_dag(steps), get_dag's
# construct_dag list, the recorded budget and select_dag's explicit replay.
# ranked_query is copied from test-ranked-search.R, where it is local.
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

test_that("the package default budget is 20000 steps and an explicit value beats the option", {
  old <- options(causaliflower.ranked_steps = NULL)
  on.exit(options(old), add = TRUE)
  expect_null(getOption("causaliflower.ranked_steps"))
  expect_identical(.cd_ranked_steps(), 20000L)
  expect_identical(.cd_ranked_steps(7), 7L)
  options(causaliflower.ranked_steps = 9L)
  expect_identical(.cd_ranked_steps(), 9L)
  expect_identical(.cd_ranked_steps(7L), 7L)
})

test_that("an explicit steps argument overrides the option in construct_dag", {
  old <- options(causaliflower.ranked_steps = NULL)
  on.exit(options(old), add = TRUE)
  graph <- ranked_query("A -> T; A -> B; B -> Y; C -> T; C -> D; D -> Y; T -> Y")
  seen <- character()
  limited <- withCallingHandlers(construct_dag(graph, quiet = TRUE, steps = 1L), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  search <- dag_options(limited)$search$total
  expect_false(search$complete)
  expect_identical(search$ranked$certificate, "budget")
  expect_identical(search$ranked$step_budget, 1L)
  expect_true(any(grepl("(step budget 1)", seen, fixed = TRUE)))
  expect_false(any(grepl("No eligible", seen)))
  options(causaliflower.ranked_steps = 1L)
  full <- construct_dag(graph, quiet = TRUE, steps = 20000L)
  for (effect in c("direct", "total")) {
    search <- dag_options(full)$search[[effect]]
    expect_true(search$complete, info = effect)
    expect_identical(search$ranked$certificate, "complete", info = effect)
    expect_identical(search$ranked$step_budget, 20000L, info = effect)
  }
  expect_identical(dag_options(full)$search_identity$ranked_steps, 20000L)
})

test_that("the budget used is recorded in the identity, the saved arguments and the call", {
  old <- options(causaliflower.ranked_steps = NULL)
  on.exit(options(old), add = TRUE)
  graph <- ranked_query("A -> T; A -> B; B -> Y; T -> M -> Y; T -> Y")
  cases <- list(list(dag = construct_dag(graph, quiet = TRUE), steps = 20000L),
    list(dag = construct_dag(graph, quiet = TRUE, steps = 500), steps = 500L),
    list(dag = construct_dag(graph, quiet = TRUE, direct_select = NULL, steps = 500L), steps = 500L),
    list(dag = construct_dag(graph, quiet = TRUE, dry_run = TRUE, steps = 500L), steps = 500L))
  expect_s3_class(cases[[3L]]$dag, "cf_dag_options")
  expect_s3_class(cases[[4L]]$dag, "cf_dag_plan")
  for (i in seq_along(cases)) {
    dag <- cases[[i]]$dag
    steps <- cases[[i]]$steps
    expect_identical(dag_options(dag)$search_identity$ranked_steps, steps, info = i)
    expect_identical(dag_options(dag)$arguments$steps, steps, info = i)
    expect_identical(dag_call(dag)$steps, steps, info = i)
    for (effect in c("direct", "total"))
      expect_identical(dag_options(dag)$search[[effect]]$ranked$step_budget, steps, info = paste(i, effect))
  }
  complete <- construct_dag(graph, quiet = TRUE, search = "complete")
  expect_false("steps" %in% names(dag_call(complete)))
  expect_false("steps" %in% names(dag_options(complete)$arguments))
  # A replayed call keeps its recorded budget whatever the option is now.
  options(causaliflower.ranked_steps = 1L)
  replay <- eval(dag_call(cases[[2L]]$dag))
  expect_identical(dag_options(replay)$search_identity$ranked_steps, 500L)
  expect_identical(dag_options(replay)$search$total$ranked$certificate, "complete")
  expect_identical(dag_options(replay)$search$total$sets, dag_options(cases[[2L]]$dag)$search$total$sets)
  expect_identical(.gd_signature(replay), .gd_signature(cases[[2L]]$dag))
})

test_that("get_dag forwards steps and its saved call records the budget used", {
  old <- options(causaliflower.ranked_steps = NULL)
  on.exit(options(old), add = TRUE)
  testthat::local_mocked_bindings(query_graph = function(crop, treatments, outcomes,
      randomised_treatment, label = NULL, quiet = NULL, params = query_params())
    ranked_query("A -> T; A -> B -> Y; T -> Y"))
  result <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE, construct_dag = list(steps = 7L))
  expect_identical(dag_call(result)$construct_dag$steps, 7L)
  expect_identical(dag_options(result)$search_identity$ranked_steps, 7L)
  expect_identical(dag_options(result)$arguments$steps, 7L)
  expect_identical(dag_options(result)$search$total$ranked$step_budget, 7L)
  defaulted <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE)
  expect_identical(dag_call(defaulted)$construct_dag$steps, 20000L)
  expect_identical(dag_options(defaulted)$search_identity$ranked_steps, 20000L)
  deferred <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE,
    construct_dag = list(direct_select = NULL, steps = 7L))
  expect_s3_class(deferred, "cf_dag_options")
  expect_identical(dag_call(deferred)$construct_dag$steps, 7L)
  plan <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE, dry_run = TRUE,
    construct_dag = list(steps = 7L))
  expect_s3_class(plan, "cf_dag_plan")
  expect_identical(dag_call(plan)$construct_dag$steps, 7L)
  witness <- get_dag("barley", "T", "Y", TRUE, quiet = TRUE, construct_dag = list(search = "witness"))
  expect_false("steps" %in% names(dag_call(witness)$construct_dag))
  expect_error(get_dag("barley", "T", "Y", TRUE, quiet = TRUE,
    construct_dag = list(search = "complete", steps = 7L)), "steps applies to the ranked search", fixed = TRUE)
})

test_that("select_dag replays a ranked catalogue under its recorded steps, not the option", {
  old <- options(causaliflower.ranked_steps = NULL)
  on.exit(options(old), add = TRUE)
  graph <- ranked_query("A -> T; A -> B; B -> Y; T -> M -> Y; T -> Y")
  dag <- construct_dag(graph, quiet = TRUE, steps = 500L)
  saved <- dag_options(dag)$search$total$sets
  options(causaliflower.ranked_steps = 1L)
  switched <- select_dag(dag, total_select = c(1L, 2L))
  expect_identical(getOption("causaliflower.ranked_steps"), 1L)
  expect_identical(dag_options(switched)$original_call$steps, 500L)
  expect_identical(dag_options(switched)$arguments$steps, dag_options(dag)$arguments$steps)
  expect_identical(dag_options(switched)$search_identity$ranked_steps, 500L)
  expect_identical(dag_options(switched)$search$total$ranked$certificate, "complete")
  expect_identical(dag_options(switched)$search$total$sets, saved)
  expect_identical(dag_options(switched)$selected$total, c(1L, 2L))
  # A receipt saved before steps was recorded replays under its search identity's budget.
  state <- .gd_state(dag)
  state$construction$arguments$steps <- NULL
  legacy <- .gd_attach(dag, state)
  expect_false("steps" %in% names(dag_options(legacy)$arguments))
  replayed <- select_dag(legacy, total_select = c(1L, 2L))
  expect_identical(dag_options(replayed)$arguments$steps, 500L)
  expect_identical(dag_options(replayed)$search_identity$ranked_steps, 500L)
  expect_identical(dag_options(replayed)$search$total$sets, saved)
})

test_that("invalid steps and steps outside the ranked route are refused", {
  graph <- ranked_query("A -> T; A -> Y; T -> Y")
  for (value in list(0, -1, 2.5, 1e12, "a", c(1, 2)))
    expect_error(construct_dag(graph, quiet = TRUE, steps = value),
      "Argument steps must be one positive whole number", fixed = TRUE, info = deparse(value))
  for (mode in c("complete", "witness"))
    expect_error(construct_dag(graph, quiet = TRUE, search = mode, steps = 5L),
      "steps applies to the ranked search (search = 'auto') only.", fixed = TRUE, info = mode)
  old <- options(causaliflower.ranked_steps = "a")
  on.exit(options(old), add = TRUE)
  expect_error(construct_dag(graph, quiet = TRUE), "Option causaliflower.ranked_steps must be", fixed = TRUE)
  expect_s3_class(construct_dag(graph, quiet = TRUE, steps = 20000L), "cf_dag")
})
