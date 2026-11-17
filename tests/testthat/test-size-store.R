# Built sizes kept for the session, and sets validated with the graph parsed once.
store_query <- function() {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];",
    "A_t1 -> T ; A_t1 -> A_t2 ; A_t2 -> Y ; B -> T ; B -> C1 ; B -> C2 ; C1 -> Y ; C2 -> Y ;",
    "T -> M_t1 ; M_t1 -> M_t2 ; M_t2 -> Y ; E -> M_t1 ; E -> Y }"))
  nodes <- names(graph)
  ids <- .strip_slice(nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = NULL, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = nodes, quantity = ids, member_id = ids, occurrence = nodes),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  graph
}

heard <- function(expr) {
  messages <- character()
  warnings <- character()
  output <- utils::capture.output(value <- withCallingHandlers(expr,
    message = function(m) {
      text <- sub("\n$", "", conditionMessage(m))
      if (!grepl("^(Sizes|Construction) \\[", text)) {
        messages <<- c(messages, text)
      }
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }))
  return(list(value = value, output = output, messages = messages, warnings = warnings))
}

test_that("a second identical construction is answered from the session store", {
  query <- store_query()
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  first <- heard(construct_dag(query, quiet = FALSE))
  stored <- .cf_cache$stats()
  expect_gt(stored$entries_by_kind[["builds"]], 0L)
  second <- heard(construct_dag(query, quiet = FALSE))
  expect_identical(second, first)
  after <- .cf_cache$stats()
  expect_identical(after$entries_by_kind[["builds"]], stored$entries_by_kind[["builds"]])
  expect_gt(after$counts[["hits"]], stored$counts[["hits"]])
  grown <- heard(grow_dag(first$value, Inf))
  expect_identical(grown$value, first$value)
  expect_identical(.cf_cache$stats()$entries_by_kind[["builds"]], stored$entries_by_kind[["builds"]])
})

test_that("without acceleration nothing is kept and the results are the same", {
  query <- store_query()
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  kept <- list(heard(construct_dag(query)), heard(construct_dag(query, 0)), heard(construct_dag(query)))
  old <- options(causaliflower.accelerate = FALSE)
  on.exit(options(old), add = TRUE)
  .cf_cache$clear()
  plain <- list(heard(construct_dag(query)), heard(construct_dag(query, 0)), heard(construct_dag(query)))
  expect_identical(.cf_cache$stats()$entries, 0L)
  expect_identical(plain, kept)
})

test_that("a size counted and then returned with declared latents equals a fresh build", {
  query <- store_query()
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  largest <- suppressWarnings(suppressMessages(construct_dag(query, quiet = TRUE)))
  sizes <- dag_options(suppressWarnings(suppressMessages(prune_dag(largest, 1, quiet = TRUE))))$sizes$variables
  expect_gt(length(sizes), 1L)
  size <- sizes[[length(sizes) - 1L]]
  latents <- "B"
  .cf_cache$clear()
  counted <- heard(construct_dag(query, size, latents = latents))
  builds <- .cf_cache$stats()$entries_by_kind[["builds"]]
  expect_identical(heard(construct_dag(query, size))$value, suppressWarnings(suppressMessages(
    construct_dag(query, size, quiet = TRUE))))
  expect_identical(.cf_cache$stats()$entries_by_kind[["builds"]], builds)
  old <- options(causaliflower.accelerate = FALSE)
  on.exit(options(old), add = TRUE)
  expect_identical(counted, heard(construct_dag(query, size, latents = latents)))
})

test_that("total-effect validation gives dagitty's answer for every set", {
  graphs <- list(
    dagitty::dagitty("dag { x -> m -> y ; z -> x ; z -> y ; w -> m ; u -> z ; u -> y ; u [latent] ; x [exposure] ; y [outcome] }"),
    dagitty::dagitty("dag { x -> m -> y ; x <-> z ; z -> y ; w -> x ; w -> y ; v -> m ; v -> y ; x [exposure] ; y [outcome] }"),
    dagitty::dagitty(paste("dag { x1 -> m -> y1 ; x2 -> m ; z -> x1 ; z -> y2 ; m -> y2 ; c -> x2 ; c -> y1 ;",
      "x1 [exposure] ; x2 [exposure] ; y1 [outcome] ; y2 [outcome] }")),
    dagitty::dagitty("pdag { x -- z ; z -> y ; x -> y ; w -> y ; x [exposure] ; y [outcome] }"))
  for (graph in graphs) {
    nodes <- names(graph)
    excluded <- c(dagitty::exposures(graph), dagitty::outcomes(graph), dagitty::latents(graph))
    sets <- c(list(character()), as.list(nodes), utils::combn(nodes, 2L, simplify = FALSE),
      list(c(nodes[[1L]], nodes[[1L]]), "absent", c(nodes[[2L]], "absent")))
    expected <- vapply(sets, function(set) {
      return(all(set %in% nodes) && !any(set %in% excluded) && isTRUE(dagitty::isAdjustmentSet(graph, set)))
    }, logical(1))
    validate <- .gd_validator(graph)
    expect_identical(validate(sets, "total"), expected)
    expect_identical(validate(rev(sets), "total"), rev(expected))
    expect_identical(.gd_validate_sets(graph, stats::setNames(sets, seq_along(sets)), "total"), expected)
    if (dagitty::graphType(graph) == "dag") {
      expect_identical(validate(sets, "direct"), .gd_validator(graph)(sets, "direct"))
      expect_identical(validate(sets, "total"), expected)
    }
  }
  expect_identical(.gd_validate_sets(graphs[[1L]], list(), "total"), logical())
  expect_error(.gd_validator(graphs[[1L]])(list(NA_character_), "total"),
    "sets must be a list of character vectors without missing values.", fixed = TRUE)
  expect_error(.gd_validator(dagitty::dagitty("dag { x -> y }"))(list("x"), "total"),
    "Both treatment and outcome tags are required to validate adjustment sets.", fixed = TRUE)
  expect_error(.gd_validator(graphs[[4L]])(list("w"), "direct"), "direct-effect adjustment sets")
  expect_identical(.gd_validator(graphs[[4L]])(list("w", "z"), "total"),
    c(isTRUE(dagitty::isAdjustmentSet(graphs[[4L]], "w")), isTRUE(dagitty::isAdjustmentSet(graphs[[4L]], "z"))))
})

test_that("stored builds never push out entries of another kind", {
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  .cf_cache$call("levels", list("small"), function() {
    return(1L)
  })
  for (index in 1:20) {
    .cf_cache$call("builds", list("large", index), function() {
      return(raw(4e6))
    })
  }
  stats <- .cf_cache$stats()
  expect_true(.cf_cache$has("levels", list("small")))
  expect_lte(stats$accounted_bytes, stats$maximum_bytes)
  expect_lt(stats$entries_by_kind[["builds"]], 20L)
})

test_that("builds that differ only in kept nodes or column requests are kept apart", {
  columns <- function(name) {
    query <- store_query()
    record <- attr(query, "query_record")
    record$requests$columns$resolved <- stats::setNames(name, name)
    attr(query, "query_record") <- record
    return(query)
  }
  builds <- function() {
    return(list(heard(construct_dag(store_query(), 0, keep_nodes = "C1")),
      heard(construct_dag(store_query(), 0, keep_nodes = "C2")),
      heard(construct_dag(columns("C1"), 0)), heard(construct_dag(columns("C2"), 0))))
  }
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  kept <- builds()
  expect_true("C1" %in% names(kept[[1L]]$value) && !"C2" %in% names(kept[[1L]]$value))
  expect_true("C2" %in% names(kept[[2L]]$value) && !"C1" %in% names(kept[[2L]]$value))
  expect_true("C1" %in% names(kept[[3L]]$value) && !"C2" %in% names(kept[[3L]]$value))
  expect_true("C2" %in% names(kept[[4L]]$value) && !"C1" %in% names(kept[[4L]]$value))
  old <- options(causaliflower.accelerate = FALSE)
  on.exit(options(old), add = TRUE)
  expect_identical(builds(), kept)
})

test_that("a build step that signals a condition other than progress is not stored", {
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  expect_warning(value <- .cd_kept(list("probe"), list("core", 1), function() {
    warning("probe warning")
    return(1L)
  }), "probe warning", fixed = TRUE)
  expect_identical(value, 1L)
  expect_message(.cd_kept(list("probe"), list("core", 2), function() {
    message("probe message")
    return(2L)
  }), "probe message", fixed = TRUE)
  expect_identical(.cf_cache$stats()$entries_by_kind[["builds"]], 0L)
  expect_message(.cd_kept(list("probe"), list("core", 3), function() {
    message("Construction [probe]")
    return(3L)
  }), "Construction [probe]", fixed = TRUE)
  expect_identical(.cf_cache$stats()$entries_by_kind[["builds"]], 1L)
  query <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y ; T <-> Y ; A -> T ; A -> Y }")
  record <- attr(store_query(), "query_record")
  record$node_map <- data.frame(node = names(query), quantity = names(query), member_id = names(query),
    occurrence = names(query))
  record$policy <- data.frame(node_id = names(query), recovered_node = "yes")
  attr(query, "query_record") <- record
  first <- heard(construct_dag(query))
  expect_gt(.cf_cache$stats()$entries_by_kind[["builds"]], 1L)
  expect_match(first$warnings, "No eligible direct/total adjustment set was found.", fixed = TRUE, all = FALSE)
  expect_identical(heard(construct_dag(query)), first)
})
