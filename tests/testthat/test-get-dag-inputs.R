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
  class(graph) <- c("dag_by_query", class(graph))
  graph
}

test_that("guided answers are parsed without executing user text", {
  expect_identical(.gd_answer("direct_select", "c(1, 2)"), c(1L, 2L))
  expect_identical(.gd_answer("total_select", "[2,1]"), c(2L, 1L))
  expect_identical(.gd_answer("direct_select", "auto"), "auto")
  expect_identical(.gd_answer("randomised_treatment", "false"), FALSE)
  expect_error(.gd_answer("direct_select", "system('anything')"), "Enter group")
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
  expect_error(do.call(get_dag, c(args, list(construct_dag = list(direct_select = "auto", direct_select = NULL)))), "unique names")
})

test_that("wrapper passes the full real query and metadata unchanged to construction", {
  q <- query_graph("barley", treatments = "N_application", outcomes = c("grain_N", "grain_wt"),
    randomised_treatment = TRUE, label = NULL, quiet = TRUE)
  received <- NULL
  testthat::local_mocked_bindings(construct_dag = function(graph, latents = NULL,
      direct_select = "auto", total_select = "auto", dry_run = FALSE, quiet = NULL) {
    received <<- list(graph = graph, latents = latents, direct_select = direct_select,
      total_select = total_select, dry_run = dry_run, quiet = quiet)
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
  expect_identical(received[-1L], list(latents = NULL, direct_select = "auto",
    total_select = "auto", dry_run = FALSE, quiet = TRUE))
  expect_identical(attr(wrapped, "get_dag_record")$query_record, actual_record)
})

test_that("wrapper is the shared query and complete construction in series", {
  testthat::local_mocked_bindings(query_graph = input_wrapper_query, .package = "causaliflower")
  q <- input_wrapper_query("barley", treatments = "N_application", outcomes = c("grain_N", "grain_wt"),
    randomised_treatment = TRUE, label = NULL, quiet = TRUE)
  manual <- suppressWarnings(construct_dag(q, quiet = TRUE, search = "complete"))
  wrapped <- suppressWarnings(get_dag("barley", "N_application", c("grain_N", "grain_wt"),
    TRUE, quiet = TRUE, construct_dag = list(search = "complete")))
  expect_identical(as.character(wrapped), as.character(manual))
  expect_identical(as.character(original_query(wrapped)), as.character(q))
  expect_identical(attr(original_query(wrapped), "query_record")$node_map,
                   attr(q, "query_record")$node_map)
  expect_equal(dag_options(wrapped)$selected, dag_options(manual)$selected)
  expect_identical(dag_checks(wrapped), dag_checks(manual))
  expect_equal(dag_call(wrapped)$randomised_treatment, TRUE)
  expect_equal(dag_call(wrapped)$query_graph$time_unroll, "all")
  expect_true("columns" %in% names(dag_call(wrapped)$query_graph))
  expect_false("max_sets" %in% names(dag_call(wrapped)$construct_dag))
  expect_true(all(vapply(dag_options(wrapped)$search, function(x) isTRUE(x$complete), logical(1))))
  expect_false("graph" %in% names(formals(get_dag)))
  expect_false("keep_nodes" %in% names(formals(get_dag)))
})

test_that("noninteractive NULL selectors return options without silently choosing", {
  testthat::local_mocked_bindings(query_graph = input_wrapper_query, .package = "causaliflower")
  x <- suppressWarnings(get_dag("barley", "N_application", "grain_N", TRUE,
    construct_dag = list(direct_select = NULL), quiet = TRUE))
  expect_s3_class(x, "cf_dag_options")
  expect_identical(x$status[["direct"]], "deferred")
  expect_null(dag_call(x)$construct_dag$direct_select)
  expect_s3_class(original_query(x), "dagitty")
  expect_identical(dag_options(x), x)
})

test_that("guided wrapper records answers and the actual choices", {
  answers <- c("barley", "N_application", "grain_N", "TRUE", "1,1")
  testthat::local_mocked_bindings(query_graph = input_wrapper_query,
    .gd_interactive = function() TRUE,
    .gd_readline = function(prompt) {
      if (!length(answers)) stop("Unexpected extra prompt: ", prompt)
      value <- answers[1]
      answers <<- answers[-1]
      value
    }, .package = "causaliflower")
  expect_message(x <- suppressWarnings(get_dag(
    construct_dag = list(direct_select = NULL), quiet = TRUE)), "responses")
  expect_s3_class(x, "cf_dag")
  expect_length(answers, 0L)
  expect_equal(dag_call(x)$construct_dag$direct_select, c(1L, 1L))
  expect_equal(nrow(.gd_state(x)$input_record$responses), 5L)
  expect_identical(dag_call(x)$crop, "barley")
})

test_that("constructor NULL and dry_run have distinct results", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> B -> Y; T -> Y }")
  choices <- suppressWarnings(construct_dag(q, direct_select = NULL, quiet = TRUE))
  plan <- suppressWarnings(construct_dag(q, dry_run = TRUE, quiet = TRUE))
  expect_s3_class(choices, "cf_dag_options")
  expect_s3_class(plan, "cf_dag_plan")
  expect_s3_class(plan$graph, "cf_dag")
  expect_identical(original_query(plan), q)
  expect_identical(dag_options(plan), plan$options)
  expect_s3_class(temporal_dag(plan), "dagitty")
})
