reuse_query <- function(columns = NULL, policy = NULL, depth = Inf) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];",
    "A -> B -> T; A -> M; T -> M -> Y; T -> Y; V -> Y }"))
  nodes <- names(graph)
  if (is.null(policy)) policy <- data.frame(node_id = nodes, recovered_node = "yes")
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "barley", columns = columns,
      params = list(water_model = "SoilWater", query_depth = depth)),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = nodes, member_id = nodes, occurrence = nodes),
    policy = policy)
  class(graph) <- c("dag_by_query", class(graph))
  graph
}

test_that("changing retained nodes reuses complete source searches and renews all results", {
  graph <- reuse_query()
  base <- construct_dag(graph, quiet = TRUE, search = "complete")
  before <- serialize(base, NULL)
  changed <- reuse_query("B")
  fresh <- construct_dag(changed, quiet = TRUE, search = "complete")
  source_calls <- check_calls <- 0L
  search <- .cd_search
  checks <- .cd_checks
  testthat::local_mocked_bindings(
    .cd_search = function(...) {
      source_calls <<- source_calls + 1L
      search(...)
    },
    .cd_checks = function(...) {
      check_calls <<- check_calls + 1L
      checks(...)
    })
  reused <- .cd_reconstruct(changed, previous = base, quiet = TRUE, search = "complete")
  expect_identical(source_calls, 0L)
  expect_gt(check_calls, 0L)
  expect_identical(reused, fresh)
  expect_identical(original_query(reused), changed)
  expect_identical(serialize(base, NULL), before)
  expect_identical(dag_checks(base)$total$required, list("A"))
  expect_identical(dag_checks(reused)$total$required, list("B"))
  expect_true(attr(dag_checks(reused), "current"))
})

test_that("reused raw families still honour explicit selection, deferred selection, dry run and node retention", {
  graph <- reuse_query()
  base <- construct_dag(graph, quiet = TRUE, search = "complete")
  source_calls <- 0L
  search <- .cd_search
  testthat::local_mocked_bindings(.cd_search = function(...) {
    source_calls <<- source_calls + 1L
    search(...)
  })
  cases <- list(list(keep_nodes = "B"), list(direct_select = NULL),
    list(total_select = c(1L, 2L)), list(keep_nodes = "V", dry_run = TRUE))
  for (control in cases) {
    args <- c(list(graph = graph, quiet = TRUE, search = "complete"), control)
    fresh <- do.call("construct_dag", args)
    source_calls <- 0L
    reused <- do.call(".cd_reconstruct", c(args, list(previous = base)))
    expect_identical(source_calls, 0L)
    expect_identical(reused, fresh)
  }
})

test_that("changed scientific identities and incomplete or damaged receipts rerun searches", {
  graph <- reuse_query()
  base <- construct_dag(graph, quiet = TRUE, search = "complete")
  source_calls <- 0L
  search <- .cd_search
  testthat::local_mocked_bindings(.cd_search = function(...) {
    source_calls <<- source_calls + 1L
    search(...)
  })
  changed_map <- graph
  attr(changed_map, "query_record")$node_map$occurrence[1L] <- "different exact occasion"
  changed_edge <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];",
    "A -> B -> T; A -> M; T -> M -> Y; T -> Y; V -> Y; V -> T }"))
  attr(changed_edge, "query_record") <- attr(graph, "query_record")
  class(changed_edge) <- class(graph)
  variants <- list(list(graph = graph, latents = "B"), list(graph = changed_map),
    list(graph = changed_edge), list(graph = reuse_query(depth = 2)))
  for (args in variants) {
    source_calls <- 0L
    result <- suppressWarnings(do.call(".cd_reconstruct",
      c(args, list(previous = base, quiet = TRUE, search = "complete"))))
    expect_identical(source_calls, 2L)
    expect_s3_class(result, "cf_dag")
  }
  for (field in c("incomplete", "missing", "version", "corrupt", "condition")) {
    state <- .gd_state(base)
    if (field == "incomplete") state$build$search$direct$complete <- FALSE
    if (field == "missing") state$construction$search_identity <- NULL
    if (field == "version") state$construction$search_identity$version <- 0L
    if (field == "corrupt") state$build$search$direct$table$label[1L] <- "wrong label"
    if (field == "condition") state$build$search$direct$condition_free <- FALSE
    state$construction$search <- state$build$search
    prior <- .gd_attach(base, state)
    source_calls <- 0L
    result <- .cd_reconstruct(graph, previous = prior, quiet = TRUE, search = "complete")
    expect_identical(source_calls, 2L, info = field)
    expect_identical(.gd_state(result)$checks, .gd_state(base)$checks)
  }
})

test_that("policy admission for a request cannot reuse a different prepared search", {
  graph <- reuse_query()
  policy <- attr(graph, "query_record")$policy
  policy$recovered_node[policy$node_id == "B"] <- "yes:user"
  graph <- reuse_query(policy = policy)
  base <- suppressWarnings(construct_dag(graph, quiet = TRUE, search = "complete"))
  changed <- reuse_query("B", policy = policy)
  fresh <- construct_dag(changed, quiet = TRUE, search = "complete")
  source_calls <- 0L
  search <- .cd_search
  testthat::local_mocked_bindings(.cd_search = function(...) {
    source_calls <<- source_calls + 1L
    search(...)
  })
  reused <- .cd_reconstruct(changed, previous = base, quiet = TRUE, search = "complete")
  expect_identical(source_calls, 2L)
  expect_identical(reused, fresh)
})

test_that("the shared private wrapper queries again and keeps current input provenance", {
  query_calls <- 0L
  synthetic_query <- function(crop, treatments, outcomes, randomised_treatment,
                              columns = NULL, label = NULL, quiet = NULL,
                              params = query_params()) {
    query_calls <<- query_calls + 1L
    graph <- reuse_query(columns)
    attr(graph, "query_record")$arguments$params <- params
    graph
  }
  testthat::local_mocked_bindings(query_graph = synthetic_query)
  args <- list(crop = "barley", treatments = "T", outcomes = "Y",
    randomised_treatment = TRUE, quiet = TRUE, construct_dag = list(search = "complete"))
  base <- do.call("get_dag", args)
  changed <- c(args, list(query_graph = list(columns = "B")))
  fresh <- do.call("get_dag", changed)
  source_calls <- 0L
  search <- .cd_search
  testthat::local_mocked_bindings(.cd_search = function(...) {
    source_calls <<- source_calls + 1L
    search(...)
  })
  reused <- do.call(".gd_get_dag", c(changed, list(previous = base)))
  expect_identical(query_calls, 3L)
  expect_identical(source_calls, 0L)
  expect_identical(reused, fresh)
  expect_identical(dag_call(reused)$query_graph$columns, "B")
  expect_false("previous" %in% names(dag_call(reused)))
  runtime <- do.call(getFromNamespace(".gd_get_dag", "causaliflower"),
    c(changed, list(previous = base)))
  expect_identical(runtime, fresh)
  expect_error(.gd_get_dag(previous = base), "requires")
  expect_error(do.call(".gd_get_dag", c(args, list(previous = base,
    query_graph = list(params = list(query_depth = 2))))), "complete query")
})

test_that("source search conditions remain observable and forbid reuse", {
  graph <- reuse_query()
  engine <- .gd_adjustment_sets
  testthat::local_mocked_bindings(.gd_adjustment_sets = function(dag, ...) {
    if (length(names(dag)) == length(names(graph))) warning("source diagnostic", call. = FALSE)
    engine(dag, ...)
  })
  conditions <- character()
  capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
    conditions <<- c(conditions, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  base <- capture(construct_dag(graph, quiet = TRUE, search = "complete"))
  expect_identical(conditions, rep("source diagnostic", 2L))
  expect_false(.gd_state(base)$build$search$direct$condition_free)
  conditions <- character()
  result <- capture(.cd_reconstruct(graph, previous = base, quiet = TRUE, search = "complete"))
  expect_identical(conditions, rep("source diagnostic", 2L))
  expect_identical(result, base)
})

test_that("disabling acceleration prevents search reuse across constructions", {
  graph <- reuse_query()
  base <- construct_dag(graph, quiet = TRUE, search = "complete")
  source_calls <- 0L
  search <- .cd_search
  testthat::local_mocked_bindings(.cd_search = function(...) {
    source_calls <<- source_calls + 1L
    search(...)
  })
  withr::local_options(causaliflower.accelerate = FALSE)
  result <- .cd_reconstruct(graph, previous = base, quiet = TRUE, search = "complete")
  expect_identical(source_calls, 2L)
  expect_identical(.gd_state(result)$checks, .gd_state(base)$checks)
  expect_identical(dag_set_catalogue(result), dag_set_catalogue(base))
})

test_that("CTYPE changes reuse only while prepared identities and existing guards agree", {
  ctype <- Sys.getlocale("LC_CTYPE")
  on.exit(invisible(Sys.setlocale("LC_CTYPE", ctype)), add = TRUE)
  if (!nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))))
    skip("A second CTYPE locale is unavailable")
  graphs <- list(reuse_query(), dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; T -> Y; Z -> T; Z -> z -> Y; ",
    "\"\u00c5\" -> T; \"\u00c5\" -> \"\u00e9\" -> Y }")))
  for (i in seq_along(graphs)) {
    graph <- graphs[[i]]
    invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
    base <- suppressWarnings(construct_dag(graph, quiet = TRUE, search = "complete"))
    before <- .cd_search_identity(.cd_prepare(graph, NULL, NULL))
    original_sets <- lapply(.gd_state(base)$build$search, `[[`, "source_sets")
    invisible(Sys.setlocale("LC_CTYPE", "C"))
    after <- .cd_search_identity(.cd_prepare(graph, NULL, NULL))
    if (i == 1L) expect_identical(after, before) else
      expect_false(identical(after, before))
    fresh <- suppressWarnings(construct_dag(graph, quiet = TRUE, search = "complete"))
    source_calls <- 0L
    search <- .cd_search
    reused <- testthat::with_mocked_bindings(
      suppressWarnings(.cd_reconstruct(graph, previous = base, quiet = TRUE, search = "complete")),
      .cd_search = function(...) {
        source_calls <<- source_calls + 1L
        search(...)
      })
    expect_identical(source_calls, if (i == 1L) 0L else 2L)
    expect_identical(reused, fresh)
    expect_identical(lapply(.gd_state(reused)$build$search, `[[`, "source_sets"),
      original_sets)
  }
})
