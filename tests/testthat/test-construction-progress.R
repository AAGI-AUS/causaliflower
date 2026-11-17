construction_parallel_paths <- function(n = 7L) {
  dagitty::dagitty(paste0("dag { T [exposure]; Y [outcome]; T -> Y; ",
    paste(sprintf("A%d -> T; A%d -> B%d -> Y;", seq_len(n), seq_len(n), seq_len(n)),
      collapse = " "), " }"))
}


test_that("searches preserve ordered results and finite witnesses", {
  for (n in c(2L, 7L, 10L)) for (effect in c("direct", "total")) {
    graph <- construction_parallel_paths(n)
    full <- .gd_adjustment_sets(graph, effect)
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
    options(causaliflower.accelerate = TRUE)
  }
})

test_that("searches use the installed separator", {
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

test_that("construction progress is concise and quiet construction keeps the same result", {
  graph <- construction_parallel_paths(2L)
  for (size in list(Inf, 0)) {
    messages <- character()
    printed <- capture.output(shown <- withCallingHandlers(
      suppressWarnings(construct_dag(graph, size)), message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }))
    expect_true(any(grepl("Preparing the reference graph", messages, fixed = TRUE)), info = size)
    expect_false(any(grepl("remaining|ETA|%", messages)), info = size)
    messages <- character()
    suppressWarnings(quiet <- withCallingHandlers(construct_dag(graph, size, quiet = TRUE),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }))
    expect_length(messages, 0L)
    expect_identical(as.character(shown), as.character(quiet))
    expect_identical(.gd_state(shown)$checks, .gd_state(quiet)$checks)
  }
})

test_that("candidate grouping preserves temporal membership and the order first seen", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y; A_t0; A_t1; B; C }")
  nodes <- names(graph)
  ids <- c(A_t0 = "A", A_t1 = "A", B = "B", C = "C", T = "T", Y = "Y")
  attr(graph, "cf_node_map") <- data.frame(node = nodes, quantity = unname(ids[nodes]),
    member_id = unname(ids[nodes]), occurrence = nodes)
  raw <- list(c("B", "C"), "A_t1", c("C", "B"), "A_t0", character(), "B")
  before <- serialize(raw, NULL)
  result <- .gd_candidates(raw, graph)
  expect_identical(result$label, c("No adjustment nodes", "A", "B", "B + C"))
  expect_identical(result$raw_sets, list(list(character()), list("A_t1", "A_t0"),
    list("B"), list(c("B", "C"), c("C", "B"))))
  expect_identical(result$group, c(1L, 2L, 2L, 3L))
  expect_identical(result$option, c(1L, 1L, 2L, 1L))
  expect_identical(serialize(raw, NULL), before)
  expect_equal(nrow(.gd_candidates(list(), graph)), 0L)
})
