# Small, explicitly synthetic query for wrapper/selection plumbing. These tests
# previously repeated the full Barley construction under its old default cap of
# 100 sets. The real full query boundary is checked separately below. Constructor
# completeness and scientific membership are tested on fully enumerable graphs.
input_wrapper_query <- function(crop, columns = NULL, treatments = NULL,
                                outcomes = NULL, randomised_treatment = FALSE,
                                time_unroll = "all", label = NULL, quiet = NULL,
                                params = query_params()) {
  text <- paste("dag {", paste0(treatments, " [exposure];", collapse = " "),
    paste0(outcomes, " [outcome];", collapse = " "),
    paste0(treatments, " -> A -> B;", collapse = " "),
    paste0("B -> ", outcomes, ";", collapse = " "),
    paste0(treatments[[1L]], " -> ", outcomes, ";", collapse = " "), "}")
  graph <- dagitty::dagitty(text)
  ids <- names(graph)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = crop, columns = columns, time_unroll = time_unroll,
      params = params, randomised_treatment = randomised_treatment),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = ids, quantity = ids, member_id = ids, occurrence = ids),
    policy = data.frame(node_id = ids, recovered_node = "yes"))
  attr(graph, "crop") <- crop
  graph
}

test_that("guided answers for query inputs are parsed without executing user text", {
  expect_identical(.gd_answer("randomised_treatment", "false"), FALSE)
  expect_error(.gd_answer("randomised_treatment", "maybe"), "TRUE or FALSE")
  expect_identical(.gd_answer("outcomes", "grain_N, grain_wt"), c("grain_N", "grain_wt"))
})

test_that("wrapper rejects retired or misplaced controls before querying", {
  args <- list(crop = "barley", treatments = "N_application", outcomes = "grain_N",
    randomised_treatment = TRUE)
  expect_error(get_dag(), "requires")
  expect_error(do.call(get_dag, c(args, list(focus = "soil_NO3"))), "unused argument")
  expect_error(do.call(get_dag, c(args, list(construct_dag = list(keep_nodes = "soil_NO3")))), "Unsupported")
  expect_error(do.call(get_dag, c(args, list(construct_dag = list(max_sets = 100L)))), "Unsupported")
  expect_error(construct_dag(max_sets = 100L), "unused argument")
  expect_error(do.call(get_dag, c(args, list(query_graph = list(data = data.frame(a = 1))))), "Unsupported")
  expect_error(do.call(get_dag, c(args, list(query_graph = list(params = list(query_depth = 1))))), "complete query")
  expect_error(do.call(get_dag, c(args, list(query_graph = list(params = list(forward_steps = 2))))), "Unsupported")
  expect_error(do.call(get_dag, c(args, list(construct_dag = list(size = 0, size = Inf)))), "unique names")
})

test_that("wrapper passes the full real query and metadata unchanged to construction", {
  q <- query_graph("barley", treatments = "N_application", outcomes = c("grain_N", "grain_wt"),
    randomised_treatment = TRUE, label = NULL, quiet = TRUE)
  received <- NULL
  testthat::local_mocked_bindings(.cd_sized = function(graph, latents, keep_nodes, dry_run, quiet, size,
    unroll, quick, requested_size) {
    received <<- list(graph = graph, latents = latents, keep_nodes = keep_nodes, dry_run = dry_run,
      quiet = quiet, size = size, unroll = unroll, quick = quick, requested_size = requested_size)
    graph
  }, .package = "causaliflower")
  wrapped <- get_dag("barley", "N_application", c("grain_N", "grain_wt"), TRUE, quiet = TRUE)
  expect_identical(as.character(received$graph), as.character(q))
  expect_identical(.gd_signature(received$graph), .gd_signature(q))
  metadata <- setdiff(names(attributes(q)), "query_record")
  expect_identical(attributes(received$graph)[metadata], attributes(q)[metadata])
  actual_record <- attr(received$graph, "query_record")
  expected_record <- attr(q, "query_record")
  # The recorded call correctly differs between direct and wrapper invocation.
  fields <- setdiff(names(expected_record), "call")
  expect_identical(actual_record[fields], expected_record[fields])
  expect_true(actual_record$scope$full_query_defaults)
  expect_identical(actual_record$arguments$params$query_depth, Inf)
  expect_identical(received[-1L], list(latents = NULL, keep_nodes = NULL, dry_run = FALSE, quiet = TRUE,
    size = 0, unroll = FALSE, quick = TRUE, requested_size = NA_real_))
})

test_that("default wrapper guidance does not invent a size request", {
  testthat::local_mocked_bindings(query_graph = input_wrapper_query, .package = "causaliflower")
  guidance <- function(choices) {
    messages <- character()
    dag <- withCallingHandlers(get_dag("barley", "N_application", "grain_N", TRUE,
      construct_dag = choices), message = function(condition) {
        text <- conditionMessage(condition)
        if (startsWith(text, "DAG with")) messages <<- c(messages, text)
        invokeRestart("muffleMessage")
      })
    list(dag = dag, messages = messages)
  }
  default <- guidance(list())
  explicit <- guidance(list(size = 0))
  expect_match(default$messages, "^DAG with 4 nodes")
  expect_false(grepl("Requested", default$messages, fixed = TRUE))
  expect_match(explicit$messages, "Requested 0 nodes. Returned 4.", fixed = TRUE)
  expect_identical(.gd_signature(default$dag), .gd_signature(explicit$dag))
})

test_that("wrapper is the shared query and construction in series", {
  testthat::local_mocked_bindings(query_graph = input_wrapper_query, .package = "causaliflower")
  q <- input_wrapper_query("barley", treatments = "N_application", outcomes = c("grain_N", "grain_wt"),
    randomised_treatment = TRUE, label = NULL, quiet = TRUE)
  for (size in list(Inf, 0)) {
    manual <- suppressWarnings(.cd_sized(q, NULL, NULL, FALSE, TRUE, size, quick = identical(size, 0)))
    wrapped <- suppressWarnings(get_dag("barley", "N_application", c("grain_N", "grain_wt"),
      TRUE, quiet = TRUE, construct_dag = list(size = size)))
    expect_identical(as.character(wrapped), as.character(manual))
    expect_identical(as.character(original_query(wrapped)), as.character(q))
    expect_identical(attr(original_query(wrapped), "query_record")$node_map, attr(q, "query_record")$node_map)
    expect_identical(dag_checks(wrapped), dag_checks(manual))
    saved <- dag_call(wrapped)
    resolved <- saved$input_record$resolved_call
    expect_equal(resolved$randomised_treatment, TRUE)
    expect_equal(resolved$query_graph$time_unroll, "all")
    expect_true("columns" %in% names(resolved$query_graph))
    expect_false("max_sets" %in% names(resolved$construct_dag))
    replay <- eval(saved)
    expect_identical(.gd_signature(replay), .gd_signature(wrapped))
    expect_identical(.gd_state(replay)$input_record, .gd_state(wrapped)$input_record)
  }
  expect_false("graph" %in% names(formals(get_dag)))
  expect_false("keep_nodes" %in% names(formals(get_dag)))
})


test_that("guided wrapper records answers to the query prompts", {
  answers <- c("barley", "N_application", "grain_N", "TRUE")
  testthat::local_mocked_bindings(query_graph = input_wrapper_query,
    .gd_interactive = function() TRUE,
    .gd_readline = function(prompt) {
      if (!length(answers)) stop("Unexpected extra prompt: ", prompt)
      value <- answers[1]
      answers <<- answers[-1]
      value
    }, .package = "causaliflower")
  expect_message(x <- suppressWarnings(get_dag(quiet = TRUE)), "responses")
  expect_s3_class(x, "cf_dag")
  expect_length(answers, 0L)
  expect_equal(nrow(.gd_state(x)$input_record$responses), 4L)
  saved <- dag_call(x)
  expect_identical(saved$input_record$resolved_call$crop, "barley")
  inputs <- .gd_state(x)$input_record
  expect_identical(saved$input_record, inputs)
  replay <- eval(saved)
  expect_identical(.gd_state(replay)$input_record, inputs)
  expect_identical(original_query(replay), original_query(x))
  expect_identical(.gd_signature(temporal_dag(replay)), .gd_signature(temporal_dag(x)))
  for (rebuilt in list(select_dag(x, quiet = TRUE), prune_dag(x, quiet = TRUE),
    grow_dag(x, quiet = TRUE), construct_dag(x, quiet = TRUE))) {
    expect_identical(.gd_state(rebuilt)$input_record, inputs)
    replay <- eval(dag_call(rebuilt))
    expect_identical(.gd_state(replay)$input_record, inputs)
  }
})

test_that("wrapper plans retain the original resolved inputs through exact replay", {
  testthat::local_mocked_bindings(query_graph = input_wrapper_query, .package = "causaliflower")
  plan <- suppressWarnings(get_dag("barley", "N_application", "grain_N", TRUE,
    dry_run = TRUE, quiet = TRUE, construct_dag = list(unroll = TRUE)))
  inputs <- attr(plan, "get_dag_record")
  expect_identical(inputs$resolved_call$construct_dag$unroll, TRUE)
  expect_identical(dag_call(plan)$input_record, inputs)
  testthat::local_mocked_bindings(query_graph = function(...) stop("Replay must not query again."),
    .package = "causaliflower")
  again <- eval(dag_call(plan))
  expect_s3_class(again, "cf_dag_plan")
  expect_identical(attr(again, "get_dag_record"), inputs)
  expect_null(.gd_state(again$graph)$input_record)
  expect_false("graph" %in% names(inputs$resolved_call))
  expect_identical(original_query(again), original_query(plan))
  expect_identical(dag_options(again)$candidate, dag_options(plan)$candidate)
  expect_identical(.gd_state(select_dag(plan, quiet = TRUE))$input_record, inputs)
  expect_identical(attr(construct_dag(plan, dry_run = TRUE, quiet = TRUE), "get_dag_record"), inputs)
  old_state <- .gd_state(plan$graph)
  old_state$construction$metric <- NULL
  old_state$construction$candidate <- NULL
  plan$graph <- .gd_attach(plan$graph, old_state)
  plan$options <- old_state$construction
  migrated <- eval(dag_call(plan))
  expect_s3_class(migrated, "cf_dag_plan")
  expect_identical(attr(migrated, "get_dag_record"), inputs)
  expect_identical(.gd_signature(temporal_dag(migrated)), .gd_signature(temporal_dag(plan)))
})

test_that("dry_run returns a plan", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> B -> Y; T -> Y }")
  plan <- suppressWarnings(construct_dag(q, dry_run = TRUE, quiet = TRUE))
  expect_s3_class(plan, "cf_dag_plan")
  expect_s3_class(plan$graph, "cf_dag")
  expect_identical(original_query(plan), q)
  expect_identical(dag_options(plan), plan$options)
  expect_s3_class(temporal_dag(plan), "dagitty")
})
