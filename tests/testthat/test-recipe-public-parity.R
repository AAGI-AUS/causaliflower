expect_recipe_route_equal <- function(pipeline, direct) {
  pipeline_graph <- if (inherits(pipeline, "cf_dag_plan")) pipeline$graph else pipeline
  direct_graph <- if (inherits(direct, "cf_dag_plan")) direct$graph else direct
  expect_identical(class(pipeline), class(direct))
  expect_identical(as.character(pipeline_graph), as.character(direct_graph))
  expect_identical(dagby_temporal(pipeline), dagby_temporal(direct))
  source <- dagby_source(pipeline)
  # The public query adds its recipe record to the otherwise identical source
  expect_s3_class(attr(source, "dagby_record")$recipe, "dagby_recipe")
  attr(source, "dagby_record") <- NULL
  attr(source, "query_record")$call <- NULL
  expected_source <- dagby_source(direct)
  attr(expected_source, "query_record")$call <- NULL
  expect_identical(source, expected_source)
  expect_identical(dagby_options(pipeline), dagby_options(direct))
  expect_identical(dagby_checks(pipeline_graph), dagby_checks(direct_graph))
}

test_that("recipe sizes match explicit backend query and construction", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  graph <- recipe_fixture_source("barley", treatments = "T", outcomes = "Y", randomised_treatment = TRUE)
  choices <- list(list(), list(size = 0), list(size = 3), list(size = "largest"),
    list(size = Inf), list(size = "smallest", keep_nodes = "C"))
  expected_nodes <- list(c("A", "T", "Y"), c("A", "T", "Y"), c("A", "T", "Y"),
    c("A", "C", "T", "Y"), c("A", "C", "T", "Y"), c("A", "C", "T", "Y"))
  for (index in seq_along(choices)) {
    recipe <- do.call(dagby, c(list("barley", "T", "Y", TRUE), choices[[index]]))
    .cf_cache$clear()
    direct <- do.call(dagby_construct, c(list(graph = graph, quiet = TRUE), choices[[index]]))
    .cf_cache$clear()
    pipeline <- recipe |> dagby_query(quiet = TRUE) |> dagby_construct(quiet = TRUE)
    expect_setequal(names(pipeline), expected_nodes[[index]])
    expect_recipe_route_equal(pipeline, direct)
    expect_identical(dagby_call(pipeline)$input_record$recipe, recipe)
  }
})

test_that("recipes carry temporal output and declarations into DAGs and plans", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  graph <- recipe_fixture_source("barley", treatments = "T", outcomes = "Y",
    randomised_treatment = FALSE, time_unroll = "all")
  for (dry_run in c(FALSE, TRUE)) for (unroll in c(FALSE, TRUE)) {
    recipe <- dagby("barley", "T", "Y", FALSE, size = 0, keep_nodes = "C",
      latents = "C", time_unroll = "all", unroll = unroll)
    .cf_cache$clear()
    direct <- dagby_construct(graph, size = 0, keep_nodes = "C", latents = "C",
      unroll = unroll, dry_run = dry_run, quiet = TRUE)
    .cf_cache$clear()
    pipeline <- recipe |> dagby_query(quiet = TRUE) |>
      dagby_construct(dry_run = dry_run, quiet = TRUE)
    expect_recipe_route_equal(pipeline, direct)
    actual <- if (dry_run) pipeline$graph else pipeline
    expect_setequal(names(actual), if (unroll) c("A_t0", "A_t1", "C", "T", "Y") else c("A", "C", "T", "Y"))
    expect_identical(as.character(dagitty::latents(actual)), "C")
    expect_setequal(names(dagby_temporal(pipeline)), c("A_t0", "A_t1", "C", "T", "Y"))
    expect_identical(dagby_options(pipeline)$arguments$keep_nodes, "C")
  }
})

test_that("recipes preserve deliberate query restrictions and their scope warning", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  params <- dagby_query_params(query_depth = 1, show_hidden_nodes = FALSE)
  graph <- recipe_fixture_source("barley", treatments = "T", outcomes = "Y",
    randomised_treatment = FALSE, params = params)
  capture_warnings <- function(expr) {
    warnings <- character()
    result <- withCallingHandlers(expr, warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    return(list(value = result, warnings = warnings))
  }
  .cf_cache$clear()
  direct <- capture_warnings(dagby_construct(graph, quiet = TRUE))
  .cf_cache$clear()
  pipeline <- capture_warnings(dagby("barley", "T", "Y", FALSE, params = params) |>
    dagby_query(quiet = TRUE) |> dagby_construct(quiet = TRUE))
  expect_identical(pipeline$warnings, direct$warnings)
  expect_identical(pipeline$warnings,
    "Checks use the supplied query_depth=1 graph only. Excluded context was not substituted from a full snapshot.")
  expect_recipe_route_equal(pipeline$value, direct$value)
  expect_identical(attr(dagby_source(pipeline$value), "query_record")$arguments$params, params)
})
