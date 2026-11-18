on_demand_query <- function() {
  branches <- unlist(lapply(1:15, function(parent) {
    paste0("U", parent, " -> U", c(2 * parent, 2 * parent + 1), ";")
  }))
  dagitty::dagitty(paste("dag { T [exposure]; Y [outcome]; U1 -> T;",
    "T -> M -> Y; T -> Y;", paste(branches, collapse = " "),
    paste0("U", 16:31, " -> Y;", collapse = " "), "}"))
}

expect_missing_query_record <- function(code) { # Assert the plain-graph warning and return the result
  warnings <- testthat::capture_warnings(value <- code)
  expect_identical(warnings, paste0(
    "No saved query record. The original requests and policy information are unavailable. ",
    "Size uses recorded variable identities or exact node names."))
  return(value)
}

test_that("a certified larger witness closes the requested group without global discovery", {
  withr::local_options(causaliflower.ranked_steps = 20000L)
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  graph <- on_demand_query()
  query <- .cd_size_query(graph, NULL)
  prefix <- .cd_choice_prefix(query, 20000L, 0, "size", NA_integer_)
  expect_false(is.null(prefix))
  expect_false(prefix$definitions$complete)
  expect_true(prefix$definitions$certified)
  expect_identical(prefix$definitions$steps, 20000L)
  expect_true(all(vapply(prefix$definitions$effects,
    function(effect) effect$ranked$step_budget == 200L, logical(1))))
  eager <- local({
    testthat::local_mocked_bindings(.cd_choice_prefix = function(...) NULL,
      .package = "causaliflower")
    expect_missing_query_record(dagby_construct(graph, 0, quiet = TRUE, unroll = TRUE))
  })
  expected <- .gd_signature(dagby_temporal(eager))
  expect_setequal(names(dagby_temporal(eager)), c("T", "Y", "M", "U1"))
  read_definitions <- .cd_level_definitions
  testthat::local_mocked_bindings(
    .cd_rule_rows = function(...) stop("The irrelevant rule family was searched."),
    .cd_level_definitions = function(query, steps) {
      if (steps != 200L) stop("The irrelevant larger ranked family was searched.")
      read_definitions(query, steps)
    }, .package = "causaliflower")
  for (mode in c("cold", "warm", "disabled")) {
    withr::local_options(causaliflower.accelerate = mode != "disabled")
    if (mode != "warm") .cf_cache$clear()
    dag <- expect_missing_query_record(dagby_construct(graph, 0, quiet = TRUE, unroll = TRUE))
    expect_identical(.gd_signature(dagby_temporal(dag)), expected)
    expect_identical(dagby_options(dag)$alternatives, dagby_options(eager)$alternatives)
    expect_identical(dagby_options(dag)$candidate, dagby_options(eager)$candidate)
    expect_identical(dagby_options(dag)$quantities, dagby_options(eager)$quantities)
    expect_identical(dagby_checks(dag), dagby_checks(eager))
    expect_true(dagby_options(dag)$alternatives_complete)
    expect_false(dagby_options(dag)$sizes_complete)
    expect_true(dagby_options(dag)$level$smallest)
    expect_false(dagby_options(dag)$level$top)
    expect_identical(dagby_options(dag)$arguments$steps, 20000L)
    expect_identical(.gd_signature(dagby_temporal(expect_missing_query_record(dagby_select(dag, quiet = TRUE)))), expected)
    expect_identical(.gd_signature(dagby_temporal(expect_missing_query_record(eval(dagby_call(dag))))), expected)
  }
})

test_that("prefix discovery preserves snapping and strict directional selection", {
  withr::local_options(causaliflower.ranked_steps = 20000L)
  graph <- on_demand_query()
  query <- .cd_size_query(graph, NULL)
  prefix <- .cd_choice_prefix(query, 20000L, 0, "size", NA_integer_)
  counts <- sort(unique(vapply(prefix$rows, `[[`, integer(1), "variables")))
  low <- counts[[1L]]
  high <- counts[[2L]]
  expect_identical(low, 4L)
  expect_true(high > low)
  expect_false(is.null(.cd_choice_prefix(query, 20000L, low + 0.1,
    "size", NA_integer_)))
  expect_null(.cd_choice_prefix(query, 20000L, (low + high) / 2,
    "prune", high + 1L))
  expect_false(is.null(.cd_choice_prefix(query, 20000L, low,
    "prune", high)))
  expect_false(is.null(.cd_choice_prefix(query, 20000L, low,
    "grow", low - 1L)))
  small <- expect_missing_query_record(dagby_construct(graph, 0, quiet = TRUE, unroll = TRUE))
  expect_message(suppressWarnings(dagby_construct(graph, low + 0.1, unroll = TRUE)),
    paste0("Requested ", low, "\\.1 variables\\. Returned ", low, ".*Next larger size: ", high))
  larger <- expect_missing_query_record(dagby_grow(small, high - low, quiet = TRUE))
  expect_identical(dagby_options(larger)$level$variables, high)
  expect_identical(.gd_signature(dagby_temporal(expect_missing_query_record(dagby_prune(larger, high - low, quiet = TRUE)))),
    .gd_signature(dagby_temporal(small)))
})

test_that("uncertain, diagnostic and same-count rule cases retain the full route", {
  graph <- on_demand_query()
  query <- .cd_size_query(graph, NULL)
  expect_null(.cd_choice_prefix(query, 200L, 0, "size", NA_integer_))
  expect_null(.cd_choice_prefix(query, 20000L, Inf, "size", NA_integer_))
  no_sets <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T <-> Y; T -> Y }")
  expect_null(.cd_choice_prefix(.cd_size_query(no_sets, NULL), 20000L,
    0, "size", NA_integer_))
  one_size <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; P -> T; P -> Y; T -> Y }")
  expect_null(.cd_choice_prefix(.cd_size_query(one_size, NULL), 20000L,
    0, "size", NA_integer_))
  testthat::local_mocked_bindings(.cd_maximal_minimal_sets = function(graph, sets, ...) {
    rep(FALSE, length(sets))
  }, .package = "causaliflower")
  expect_null(.cd_choice_prefix(query, 20000L, 0, "size", NA_integer_))
})

test_that("cardinality closure reuses the prefix while resolving same-size rule choices", {
  withr::local_options(causaliflower.ranked_steps = 20000L)
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome]; U -> T;",
    "T -> M -> Y; T -> Y;", paste0("U -> A", 1:20, " -> Y;", collapse = " "), "}"))
  query <- .cd_size_query(graph, NULL)
  prefix <- .cd_choice_prefix(query, 20000L, 0, "size", NA_integer_, retain = TRUE)
  expect_false(is.null(prefix))
  expect_false(prefix$closed)
  expect_false(prefix$definitions$complete)
  expect_length(prefix$rows, 1L)
  expect_true(all(.cd_prefix_bounds(query, prefix$definitions) > 4))
  rules <- .cd_rule_rows(query, TRUE)
  extended <- .cd_choice_extend(query, prefix, rules, 0, "size", NA_integer_)
  expect_true(extended$closed)
  expect_false(extended$complete)
  eager <- local({
    testthat::local_mocked_bindings(.cd_choice_prefix = function(...) NULL,
      .package = "causaliflower")
    expect_missing_query_record(dagby_construct(graph, 0, quiet = TRUE, unroll = TRUE))
  })
  read_definitions <- .cd_level_definitions
  testthat::local_mocked_bindings(.cd_level_definitions = function(query, steps) {
    if (steps != 200L) stop("The irrelevant larger ranked family was searched.")
    read_definitions(query, steps)
  }, .package = "causaliflower")
  for (mode in c("cold", "warm", "disabled")) {
    withr::local_options(causaliflower.accelerate = mode != "disabled")
    if (mode != "warm") .cf_cache$clear()
    dag <- expect_missing_query_record(dagby_construct(graph, 0, quiet = TRUE, unroll = TRUE))
    expect_setequal(names(dagby_temporal(dag)), c("T", "Y", "U", "M"))
    expect_identical(.gd_signature(dagby_temporal(dag)), .gd_signature(dagby_temporal(eager)))
    expect_identical(dagby_options(dag)$alternatives, dagby_options(eager)$alternatives)
    expect_identical(dagby_options(dag)$quantities, dagby_options(eager)$quantities)
    expect_identical(dagby_checks(dag), dagby_checks(eager))
    expect_true(dagby_options(dag)$alternatives_complete)
    expect_false(dagby_options(dag)$sizes_complete)
    expect_identical(dagby_options(dag)$candidate, dagby_options(eager)$candidate)
    expect_identical(.gd_signature(dagby_temporal(expect_missing_query_record(eval(dagby_call(dag))))), .gd_signature(dagby_temporal(eager)))
  }
})

test_that("known options omitted from levels still bound complete effects", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U -> T; U -> A; U -> B; A -> Y; B -> Y; T -> Y }")
  query <- .cd_size_query(graph, NULL)
  definitions <- .cd_level_definitions(query, 20000L)
  expect_true(definitions$complete)
  expect_length(definitions$levels, 2L)
  definitions$levels <- definitions$levels[1L]
  bounds <- .cd_prefix_bounds(query, definitions)
  expect_identical(unname(bounds), c(4, 4))
  expect_true(all(is.finite(bounds)))
})
