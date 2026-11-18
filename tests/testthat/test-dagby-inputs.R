test_that("recipe creation validates inputs without querying or constructing", {
  testthat::local_mocked_bindings(
    .dagby_query_source = function(...) stop("Recipe creation queried the source."),
    .cd_sized = function(...) stop("Recipe creation constructed a graph."),
    .package = "causaliflower")
  recipe <- dagby("barley", "T", "Y", FALSE)
  expect_s3_class(recipe, "dagby_recipe")
  expect_named(recipe, c("query", "construction"))
  expect_identical(recipe$construction,
    list(size = "smallest", latents = NULL, keep_nodes = NULL, unroll = FALSE))
  expect_null(recipe$query$label)
  expect_false("time_unroll" %in% names(recipe$query))
  expect_identical(recipe$query$params, dagby_query_params())
  expect_output(print(recipe), "barley")
  expect_error(dagby(), "requires")
  expect_error(dagby("barley", "T", "Y"), "randomised_treatment")
  exploration <- dagby("barley", NULL, NULL, FALSE)
  expect_null(exploration$query$treatments)
  expect_null(exploration$query$outcomes)
  text <- 'stop("This is a node name, not executable input.")'
  expect_identical(dagby("barley", text, "Y", FALSE)$query$treatments, text)
})

test_that("invalid recipe settings fail before source work", {
  testthat::local_mocked_bindings(.dagby_query_source = function(...) stop("Source work occurred."),
    .package = "causaliflower")
  arguments <- list(crop = "barley", treatments = "T", outcomes = "Y", randomised_treatment = FALSE)
  invalid <- list(crop = "", treatments = 1, outcomes = NA_character_,
    randomised_treatment = NA, size = "medium", keep_nodes = 1, latents = NA_character_,
    columns = "", time_unroll = "sometimes", label = "pretty", unroll = NA)
  for (name in names(invalid)) {
    changed <- arguments
    changed[name] <- invalid[name]
    expect_error(do.call(dagby, changed), name, info = name)
  }
  expect_error(dagby("barley", "T", "Y", FALSE, size = -1), "size")
  expect_error(dagby("barley", "T", "Y", FALSE, params = list(forward_steps = 2)), "params")
  expect_error(dagby("barley", "T", "Y", FALSE, params = list(2)), "names")
  expect_error(dagby("barley", "T", "Y", FALSE, params = list(query_depth = 1, query_depth = 2)), "unique names")
  for (retired in c("query", "construct", "quiet", "dry_run", "focus")) {
    expect_error(do.call(dagby, c(arguments, stats::setNames(list(TRUE), retired))), "unused argument")
  }
  expect_error(dagby_query("barley"), "recipe")
  recipe <- dagby("barley", "T", "Y", FALSE)
  recipe$construction$size <- "medium"
  expect_error(dagby_query(recipe), "size")
})

test_that("recipes preserve labels, specialist settings and omitted temporal settings", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  omitted <- dagby("barley", "T", "Y", FALSE)
  explicit <- dagby("barley", "T", "Y", FALSE, time_unroll = NULL)
  expect_false("time_unroll" %in% names(omitted$query))
  expect_true("time_unroll" %in% names(explicit$query))
  expect_null(explicit$query$time_unroll)
  expect_false("time_unroll" %in% names(dagby_call(omitted)))
  expect_true("time_unroll" %in% names(dagby_call(explicit)))
  expect_true(attr(dagby_query(omitted, quiet = TRUE), "query_record")$time_unroll_defaulted)
  expect_false(attr(dagby_query(explicit, quiet = TRUE), "query_record")$time_unroll_defaulted)
  for (label in c("id", "input", "readable", "APSIM")) {
    recipe <- dagby("barley", "T", "Y", TRUE, label = label,
      params = list(query_depth = 1, show_hidden_nodes = FALSE))
    expected <- if (label == "id") NULL else label
    expect_identical(recipe$query$label, expected)
    expect_identical(recipe$query$params, dagby_query_params(query_depth = 1, show_hidden_nodes = FALSE))
    expect_identical(attr(dagby_query(recipe, quiet = TRUE), "query_record")$arguments$label, expected)
    expect_identical(dagby_call(recipe)$label, label)
    expect_identical(eval(dagby_call(recipe)), recipe)
  }
})

test_that("query recipes and saved calls use the public namespace", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  recipe <- dagby("barley", "T", "Y", FALSE, size = "largest", columns = "C")
  query <- dagby_query(recipe, quiet = TRUE)
  expect_identical(attr(query, "dagby_record")$recipe, recipe)
  saved <- dagby_call(query)
  expect_identical(saved[[1L]], quote(causaliflower::dagby_query))
  expect_identical(saved$quiet, TRUE)
  expect_identical(saved[[2L]][[1L]], quote(causaliflower::dagby))
  expect_identical(saved[[2L]]$size, "largest")
  expect_identical(saved[[2L]]$columns, "C")
  expect_identical(saved, attr(query, "dagby_record")$resolved_call)
  caller <- new.env(parent = globalenv())
  caller$dagby <- caller$dagby_query <- function(...) stop("A caller binding replaced a public function.")
  expect_identical(eval(dagby_call(recipe), caller), recipe)
  expect_identical(eval(saved, caller), query)
})

test_that("construction inherits recipe preferences and explicit overrides win", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  query <- dagby("barley", "T", "Y", FALSE, size = "largest", latents = "C",
    keep_nodes = "C", time_unroll = "all", unroll = TRUE) |> dagby_query(quiet = TRUE)
  inherited <- dagby_construct(query, quiet = TRUE)
  expect_setequal(names(inherited), c("A_t0", "A_t1", "C", "T", "Y"))
  expect_identical(as.character(dagitty::latents(inherited)), "C")
  expect_identical(dagby_options(inherited)$arguments$keep_nodes, "C")
  expect_true(dagby_options(inherited)$arguments$unroll)
  override <- dagby_construct(query, size = "smallest", latents = NULL,
    keep_nodes = NULL, unroll = FALSE, quiet = TRUE)
  expect_setequal(names(override), c("A", "T", "Y"))
  expect_length(dagitty::latents(override), 0L)
  expect_length(dagby_options(override)$arguments$keep_nodes, 0L)
  expect_false(dagby_options(override)$arguments$unroll)
  expect_identical(dagby_source(override), query)
  raw <- recipe_fixture_source("barley", treatments = "T", outcomes = "Y")
  compact <- dagby_construct(raw, quiet = TRUE)
  largest <- dagby_construct(raw, size = "largest", quiet = TRUE)
  expect_setequal(names(compact), c("A", "T", "Y"))
  expect_setequal(names(largest), c("A", "C", "T", "Y"))
  expect_true(dagby_options(compact)$level$smallest)
  expect_false(dagby_options(largest)$level$smallest)
})

test_that("default recipe guidance describes the smallest size without a zero request", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  messages <- character()
  dag <- withCallingHandlers(dagby("barley", "T", "Y", FALSE) |>
    dagby_query(quiet = TRUE) |> dagby_construct(), message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  expect_setequal(names(dag), c("A", "T", "Y"))
  guidance <- paste(messages, collapse = "\n")
  expect_match(guidance, "DAG with 3 variables. This is the smallest size.", fixed = TRUE)
  expect_match(guidance, "Next larger size: 4 variables.", fixed = TRUE)
  expect_false(grepl("Requested 0", guidance, fixed = TRUE))
})

test_that("selection, resizing and fresh construction retain recipes and exact replay", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  recipe <- dagby("barley", "T", "Y", FALSE, size = "largest", latents = "C", unroll = TRUE)
  query <- dagby_query(recipe, quiet = TRUE)
  dag <- dagby_construct(query, quiet = TRUE)
  inputs <- .gd_state(dag)$input_record
  expect_identical(inputs$recipe, recipe)
  expect_identical(inputs$resolved_call, dagby_call(query))
  testthat::local_mocked_bindings(.dagby_query_source = function(...) stop("Replay must not query again."),
    .package = "causaliflower")
  for (rebuilt in list(dagby_select(dag, quiet = TRUE), dagby_prune(dag, quiet = TRUE),
    dagby_grow(dag, quiet = TRUE), dagby_construct(dag, quiet = TRUE))) {
    expect_identical(.gd_state(rebuilt)$input_record, inputs)
    expect_identical(dagby_source(rebuilt), query)
    expect_true(dagby_options(rebuilt)$arguments$unroll)
    expect_identical(dagby_options(rebuilt)$arguments$latents, "C")
    replay <- eval(dagby_call(rebuilt))
    expect_identical(.gd_state(replay)$input_record, inputs)
    expect_identical(.gd_signature(replay), .gd_signature(rebuilt))
    expect_identical(dagby_checks(replay), dagby_checks(rebuilt))
  }
  retained <- dagby_construct(query, keep_nodes = "C", quiet = TRUE)
  fresh <- dagby_construct(retained, quiet = TRUE)
  expect_identical(dagby_options(fresh)$arguments$keep_nodes, "C")
  expect_identical(dagby_options(fresh)$arguments$latents, "C")
  expect_identical(as.character(dagitty::latents(fresh)), "C")
  expect_setequal(names(fresh), c("A", "C", "T", "Y"))
  expect_identical(.gd_state(fresh)$input_record, inputs)
})

test_that("plans retain recipe inputs through selection and exact replay", {
  testthat::local_mocked_bindings(.dagby_query_source = recipe_fixture_source, .package = "causaliflower")
  recipe <- dagby("barley", "T", "Y", TRUE, unroll = TRUE)
  plan <- recipe |> dagby_query(quiet = TRUE) |> dagby_construct(dry_run = TRUE, quiet = TRUE)
  inputs <- attr(plan, "dagby_record")
  expect_identical(inputs$recipe, recipe)
  expect_identical(dagby_call(plan)$input_record, inputs)
  testthat::local_mocked_bindings(.dagby_query_source = function(...) stop("Replay must not query again."),
    .package = "causaliflower")
  again <- eval(dagby_call(plan))
  expect_s3_class(again, "cf_dag_plan")
  expect_identical(attr(again, "dagby_record"), inputs)
  expect_null(.gd_state(again$graph)$input_record)
  expect_identical(dagby_source(again), dagby_source(plan))
  expect_identical(dagby_options(again)$candidate, dagby_options(plan)$candidate)
  expect_identical(.gd_state(dagby_select(plan, quiet = TRUE))$input_record, inputs)
  expect_identical(attr(dagby_construct(plan, dry_run = TRUE, quiet = TRUE), "dagby_record"), inputs)
  old_state <- .gd_state(plan$graph)
  old_state$construction$metric <- NULL
  old_state$construction$candidate <- NULL
  plan$graph <- .gd_attach(plan$graph, old_state)
  plan$options <- old_state$construction
  migrated <- eval(dagby_call(plan))
  expect_s3_class(migrated, "cf_dag_plan")
  expect_identical(attr(migrated, "dagby_record"), inputs)
  expect_identical(.gd_signature(dagby_temporal(migrated)), .gd_signature(dagby_temporal(plan)))
})

test_that("the README recipe produces a checked compact DAG from the unchanged source", {
  recipe <- dagby(crop = "barley", treatments = "nitrogen_applied", outcomes = "grain_protein",
    randomised_treatment = FALSE, size = "smallest")
  expect_no_warning(query <- dagby_query(recipe, quiet = TRUE))
  expect_identical(eval(dagby_call(query)), query)
  backend <- .dagby_query_source("barley", treatments = "nitrogen_applied", outcomes = "grain_protein",
    randomised_treatment = FALSE, quiet = TRUE)
  expect_identical(as.character(query), as.character(backend))
  metadata <- setdiff(names(attributes(backend)), "query_record")
  expect_identical(attributes(query)[metadata], attributes(backend)[metadata])
  public_record <- attr(query, "query_record")
  backend_record <- attr(backend, "query_record")
  fields <- setdiff(names(backend_record), "call")
  expect_identical(public_record[fields], backend_record[fields])
  expect_true(public_record$scope$full_query_defaults)
  expect_identical(public_record$arguments$params$query_depth, Inf)
  expect_no_warning(crop_dag <- query |> dagby_construct(quiet = TRUE))
  expect_identical(dagby_options(crop_dag)$level$variables, 6L)
  expect_length(names(crop_dag), 8L)
  expect_length(names(dagby_temporal(crop_dag)), 14L)
  expect_identical(vapply(dagby_checks(crop_dag), `[[`, "", "status"),
    c(direct = "preserved", total = "preserved"))
  expect_identical(dagby_source(crop_dag), query)
  expect_identical(dagby_call(crop_dag)$input_record$recipe, recipe)
})

test_that("dry_run returns a plan for an ordinary graph", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> B -> Y; T -> Y }")
  plan <- suppressWarnings(dagby_construct(q, dry_run = TRUE, quiet = TRUE, size = "largest"))
  expect_s3_class(plan, "cf_dag_plan")
  expect_s3_class(plan$graph, "cf_dag")
  expect_identical(dagby_source(plan), q)
  expect_identical(dagby_options(plan), plan$options)
  expect_s3_class(dagby_temporal(plan), "dagitty")
})
