construction_parallel_paths <- function(n = 7L) {
  dagitty::dagitty(paste0("dag { T [exposure]; Y [outcome]; T -> Y; ",
    paste(sprintf("A%d -> T; A%d -> B%d -> Y;", seq_len(n), seq_len(n), seq_len(n)),
      collapse = " "), " }"))
}

test_that("construction searches the complete families beyond the former public cap", {
  graph <- construction_parallel_paths()
  expect_false("max_sets" %in% names(formals(construct_dag)))
  result <- suppressWarnings(construct_dag(graph, direct_select = NULL, quiet = TRUE,
    search = "complete"))
  for (effect in c("direct", "total")) {
    reference <- unname(lapply(dagitty::adjustmentSets(graph, effect = effect), as.character))
    # Dagitty's R adapter sorts names within each set. The private direct
    # adapter historically retains membership traversal order instead.
    expect_identical(lapply(result$search[[effect]]$sets, sort), lapply(reference, sort))
    expect_length(result$search[[effect]]$sets, 128L)
    expect_true(result$search[[effect]]$complete)
    expect_identical(result$search[[effect]]$limit, Inf)
  }
  expect_identical(result$refinement$proposal_limit, Inf)
  expect_length(result$refinement$passes, 0L)
  expect_false("max_sets" %in% names(result$arguments))
})

test_that("progress instrumentation preserves ordered results and finite witnesses", {
  for (n in c(2L, 7L, 10L)) for (effect in c("direct", "total")) {
    graph <- construction_parallel_paths(n)
    full <- .gd_adjustment_sets(graph, effect)
    shown <- .gd_adjustment_sets(graph, effect, progress = function(...) NULL)
    expect_identical(shown, full)
    expect_identical(length(full$sets), as.integer(2^n))
    bounded <- .gd_adjustment_sets(graph, effect, 1L)
    expect_identical(bounded$sets, full$sets[1L])
    expect_identical(bounded$probe, full$sets[[2L]])
    expect_false(bounded$complete)
  }
})

test_that("acceleration settings retain exact order and private witnesses", {
  graph <- construction_parallel_paths(3L)
  withr::local_options(causaliflower.accelerate = TRUE)
  for (effect in c("direct", "total")) for (limit in c(1L, Inf)) {
    accelerated <- .gd_adjustment_sets(graph, effect, limit)
    options(causaliflower.accelerate = FALSE)
    expect_identical(.gd_adjustment_sets(graph, effect, limit), accelerated)
    expect_identical(.gd_adjustment_sets(graph, effect, limit,
      progress = function(...) NULL), accelerated)
    options(causaliflower.accelerate = TRUE)
  }
})

test_that("both quiet and observed searches use the installed separator", {
  graph <- construction_parallel_paths(2L)
  withr::local_options(causaliflower.accelerate = TRUE)
  real_context <- V8::v8
  testthat::local_mocked_bindings(v8 = function(...) {
    context <- real_context(...)
    methods <- as.list(context)
    methods$assign <- function(name, value, ...) {
      assigned <- context$assign(name, value, ...)
      if (identical(name, "gd_graph_text")) {
        context$eval(paste(
          "DAGitty.GraphAnalyzer.listMinimalSeparators = function(graph, must, mustNot, max) {",
          "  throw new Error('installed route reached');",
          "};", sep = "\n"))
      }
      return(assigned)
    }
    return(methods)
  }, .package = "V8")
  for (effect in c("direct", "total")) {
    expect_error(.gd_adjustment_sets(graph, effect), "installed route reached", fixed = TRUE)
    expect_error(.gd_adjustment_sets(graph, effect, progress = function(...) {
      return(invisible(NULL))
    }), "installed route reached", fixed = TRUE)
  }
  expect_false(exists(".gd_separator_reference", envir = asNamespace("causaliflower"), inherits = FALSE))
  expect_false(exists(".gd_indexed_separator_source", envir = asNamespace("causaliflower"), inherits = FALSE))
})

test_that("construction enumeration retains the verified dagitty version", {
  withr::local_options(causaliflower.accelerate = TRUE)
  real_version <- utils::packageVersion
  testthat::local_mocked_bindings(packageVersion = function(package, ...) {
    if (identical(package, "dagitty")) {
      return(package_version("0.3.5"))
    }
    return(real_version(package, ...))
  }, .package = "utils")
  for (effect in c("direct", "total")) {
    expect_error(.gd_adjustment_sets(construction_parallel_paths(2L), effect),
                 "supports dagitty 0.3-4 only", fixed = TRUE)
  }
})

test_that("construction enumeration checks the installed separator API", {
  withr::local_options(causaliflower.accelerate = TRUE)
  real_context <- V8::v8
  testthat::local_mocked_bindings(v8 = function(...) {
    context <- real_context(...)
    methods <- as.list(context)
    methods$source <- function(...) {
      sourced <- context$source(...)
      context$eval("delete DAGitty.GraphAnalyzer.listMsasDirectEffect;")
      return(sourced)
    }
    return(methods)
  }, .package = "V8")
  for (effect in c("direct", "total")) {
    expect_error(.gd_adjustment_sets(construction_parallel_paths(2L), effect),
                 "lacks the verified separator API", fixed = TRUE)
  }
})

test_that("large enumeration transfers preserve complete ordered results in chunks", {
  graph <- construction_parallel_paths(11L)
  withr::local_options(causaliflower.accelerate = FALSE)
  reference <- lapply(c("direct", "total"), function(effect)
    .gd_adjustment_sets(graph, effect))
  options(causaliflower.accelerate = TRUE)
  chunks <- 0L
  real_context <- V8::v8
  testthat::local_mocked_bindings(v8 = function(...) {
    context <- real_context(...)
    methods <- as.list(context)
    methods$get <- function(name, ...) {
      if (identical(name, "gd_result")) stop("large result was not chunked")
      if (identical(name, "gd_chunk")) chunks <<- chunks + 1L
      context$get(name, ...)
    }
    methods
  }, .package = "V8")
  for (i in seq_along(reference)) {
    result <- .gd_adjustment_sets(graph, c("direct", "total")[[i]])
    expect_length(result$sets, 2048L)
    expect_identical(result, reference[[i]])
  }
  expect_identical(chunks, 6L)
})

test_that("complete search progress is concise and quiet construction keeps the same result", {
  graph <- construction_parallel_paths(2L)
  messages <- character()
  printed <- capture.output(shown <- withCallingHandlers(
    suppressWarnings(construct_dag(graph, search = "complete")), message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_length(printed, 0L)
  expect_true(any(grepl("Preparing the reference graph", messages, fixed = TRUE)))
  expect_true(any(grepl("direct search complete (4 alternatives)", messages, fixed = TRUE)))
  expect_true(any(grepl("total search complete (4 alternatives)", messages, fixed = TRUE)))
  expect_true(any(grepl("Complete (", messages, fixed = TRUE)))
  expect_false(any(grepl("remaining|ETA|%", messages)))
  messages <- character()
  suppressWarnings(quiet <- withCallingHandlers(construct_dag(graph, quiet = TRUE, search = "complete"),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_length(messages, 0L)
  expect_identical(as.character(shown), as.character(quiet))
  expect_identical(.gd_state(shown)$checks, .gd_state(quiet)$checks)
  expect_identical(dag_set_catalogue(shown), dag_set_catalogue(quiet))
})

test_that("enumeration heartbeats report existing counts without an R evaluation callback", {
  real_context <- V8::v8
  testthat::local_mocked_bindings(v8 = function(...) {
    context <- real_context(...)
    methods <- as.list(context)
    methods$eval <- function(src, ...) {
      if (grepl("var[[:space:]]+gd_steps[[:space:]]*=[[:space:]]*0;", src)) {
        context$eval("var clock=0; Date.now=function(){clock+=3000;return clock;};")
      }
      context$eval(src, ...)
    }
    methods
  }, .package = "V8")
  graph <- construction_parallel_paths(8L)
  printed <- capture.output(result <-
    .gd_adjustment_sets(graph, "total", progress = function(...) NULL))
  expect_length(result$sets, 256L)
  expect_true(any(grepl("total search \\[.* alternatives found", printed)))
  expect_false(any(grepl("remaining|ETA|%", printed)))
  counts <- as.integer(sub(".*: ([0-9]+) alternatives found.*", "\\1", printed))
  expect_true(all(counts >= 0L & counts <= length(result$sets)))
  expect_length(.gd_adjustment_sets(graph, "direct")$sets, 256L)
})

test_that("candidate grouping preserves temporal membership and the order first seen", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y; A_t0; A_t1; B; C }")
  raw <- list(c("B", "C"), "A_t1", c("C", "B"), "A_t0", character(), "B")
  before <- serialize(raw, NULL)
  result <- .gd_candidates(raw, graph, list())
  expect_identical(result$label, c("No adjustment nodes", "A", "B", "B + C"))
  expect_identical(result$raw_sets, list(list(character()), list("A_t1", "A_t0"),
    list("B"), list(c("B", "C"), c("C", "B"))))
  expect_identical(result$group, c(1L, 2L, 2L, 3L))
  expect_identical(result$option, c(1L, 1L, 2L, 1L))
  expect_identical(serialize(raw, NULL), before)
  expect_equal(nrow(.gd_candidates(list(), graph, list())), 0L)
})

test_that("verification still distinguishes exact, redundant, incomplete and invalid sets", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; Z -> T; Z -> Y; W; T -> Y }")
  cases <- list(
    list(base = "Z", final = list("Z"), complete = TRUE, status = "preserved"),
    list(base = c("Z", "W"), final = list("Z"), complete = FALSE, status = "valid but not minimal"),
    list(base = "Z", final = list(), complete = FALSE, status = "minimality not verified"),
    list(base = "Z", final = list(), complete = TRUE, status = "valid but not minimal"),
    list(base = "W", final = list("Z"), complete = TRUE, status = "not valid")
  )
  for (case in cases) {
    checked <- .gd_check_family(graph, graph, list(case$base), case$final,
      "total", case$complete, list(case$base))
    expect_identical(checked$records[[1L]]$status, case$status)
    expect_true(checked$records[[1L]]$chosen)
  }
})
