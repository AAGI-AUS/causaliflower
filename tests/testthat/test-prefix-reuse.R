prefix_reuse_input <- function() {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U -> T; U -> A; U -> B; A -> Y; B -> Y; T -> M -> Y; T -> Y; W -> T }")
  query <- .cd_size_query(graph, NULL)
  list(graph = graph, query = query, definitions = .cd_level_definitions(query, 200L),
    map = .cd_variable_map(query$observed$analysis))
}

test_that("pure prefix proofs reuse only bounded build storage", {
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  input <- prefix_reuse_input()
  proof <- function() .cd_prefix_proof(input$query, input$definitions, input$map)
  first <- proof()
  stored <- .cf_cache$stats()
  expect_identical(stored$entries_by_kind[["builds"]], 1L)
  expect_identical(proof(), first)
  expect_identical(.cf_cache$stats()$counts[["hits"]], stored$counts[["hits"]] + 1L)
  expect_identical(withr::with_options(list(causaliflower.accelerate = FALSE), proof()), first)
  for (i in seq_len(25L)) .cf_cache$call("builds", list("proof pressure", i), function() i)
  pressured <- .cf_cache$stats()
  expect_identical(pressured$entries_by_kind[["builds"]], pressured$maximum_entries[["builds"]])
  expect_identical(pressured$entries_by_kind[["levels"]], stored$entries_by_kind[["levels"]])
  expect_identical(proof(), first)
  expect_gt(.cf_cache$stats()$counts[["misses"]], pressured$counts[["misses"]])
  .cf_cache$clear()
  expect_identical(proof(), first)
})

test_that("proof keys bind exact definitions and current retention inputs", {
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  input <- prefix_reuse_input()
  first <- .cd_prefix_proof(input$query, input$definitions, input$map)
  changed <- input$query
  changed$observed$required <- union(changed$observed$required, "W")
  second <- .cd_prefix_proof(changed, input$definitions, input$map)
  expect_identical(second, withr::with_options(list(causaliflower.accelerate = FALSE),
    .cd_prefix_proof(changed, input$definitions, input$map)))
  expect_true(all(vapply(second$rows, function(row) "W" %in% row$pins, logical(1))))
  expect_false(identical(second, first))
  definitions <- input$definitions
  definitions$steps <- 201L
  before <- .cf_cache$stats()
  expect_identical(.cd_prefix_proof(input$query, definitions, input$map), first)
  expect_gt(.cf_cache$stats()$counts[["misses"]], before$counts[["misses"]])
  changed <- input$query
  attr(changed$observed$analysis, "input_labels") <- "private analysis label"
  before <- .cf_cache$stats()
  expect_identical(.cd_prefix_proof(changed, input$definitions, input$map), first)
  expect_gt(.cf_cache$stats()$counts[["misses"]], before$counts[["misses"]])
})

test_that("an external proof wrapper runs uncached without changing namespace lookup", {
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  input <- prefix_reuse_input()
  original <- .cd_prefix_proof
  first <- original(input$query, input$definitions, input$map)
  calls <- 0L
  testthat::local_mocked_bindings(.cd_prefix_proof = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "causaliflower")
  expect_identical(.cd_prefix_proof(input$query, input$definitions, input$map), first)
  expect_identical(.cd_prefix_proof(input$query, input$definitions, input$map), first)
  expect_identical(calls, 2L)
})

test_that("captured validator and minimality replacements bypass stored proofs", {
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  input <- prefix_reuse_input()
  proof <- function() .cd_prefix_proof(input$query, input$definitions, input$map)
  first <- proof()
  local({
    testthat::local_mocked_bindings(.gd_validator = function(dag) stop("validator probe"),
      .package = "causaliflower")
    expect_error(proof(), "validator probe")
    expect_error(proof(), "validator probe")
  })
  replacement <- function(answer) {
    force(answer)
    function(graph, sets, effect, validate) rep(answer, length(sets))
  }
  for (answer in c(FALSE, TRUE, FALSE)) local({
    testthat::local_mocked_bindings(.cd_maximal_minimal_sets = replacement(answer),
      .package = "causaliflower")
    if (answer) expect_identical(proof(), first) else expect_null(proof())
  })
  expect_identical(proof(), first)
})

test_that("warning-bearing native proof implementations are never stored", {
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  input <- prefix_reuse_input()
  original <- .cd_maximal_minimal_sets
  warned <- original
  body(warned) <- substitute({ warning("proof check", call. = FALSE); BODY }, list(BODY = body(original)))
  testthat::local_mocked_bindings(.cd_maximal_minimal_sets = warned, .package = "causaliflower")
  observe <- function() {
    warnings <- character()
    value <- withCallingHandlers(.cd_prefix_proof(input$query, input$definitions, input$map),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      })
    list(value = value, warnings = warnings)
  }
  first <- observe()
  expect_identical(first$warnings, rep("proof check", 2L))
  expect_identical(observe(), first)
  expect_identical(.cf_cache$stats()$entries_by_kind[["builds"]], 0L)
})
