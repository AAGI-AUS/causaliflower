acceleration_mode <- function(enabled, code) {
  old <- options(causaliflower.accelerate = enabled)
  on.exit(options(old))
  force(code)
}

test_that("compiled kernels are registered in the package namespace", {
  acceleration_mode(TRUE, {
    expect_s3_class(.cf_native_symbol("C_cf_project_edges"), "NativeSymbolInfo")
    expect_s3_class(.cf_native_symbol("C_cf_curve_contact"), "NativeSymbolInfo")
    expect_identical(environment(.cf_native_symbol), asNamespace("causaliflower"))
  })
  acceleration_mode(FALSE, expect_null(.cf_native_symbol("C_cf_project_edges")))
  local_mocked_bindings(.cf_native_symbol = function(name) NULL)
  edges <- data.frame(v = c("A", "B"), e = "->", w = c("B", "C"))
  expect_identical(collapse_project_edges(c("A", "B", "C"), edges, "B"),
    collapse_project_edges_R(c("A", "B", "C"), edges, "B"))
})

test_that("registered projection preserves complete typed edge tables", {
  set.seed(7499)
  acceleration_mode(TRUE, {
    for (n in c(0:10, 20L)) for (iteration in seq_len(8L)) {
      vertices <- if (n) paste0("n", seq_len(n)) else character()
      if (n > 2L) vertices[2:3] <- c("\u00e9", 'quoted " name')
      grid <- which(upper.tri(matrix(FALSE, n, n)), arr.ind = TRUE)
      directed <- grid[runif(nrow(grid)) < .2, , drop = FALSE]
      bidirected <- grid[runif(nrow(grid)) < .15, , drop = FALSE]
      edges <- rbind(data.frame(v = vertices[directed[, 1L]], e = rep("->", nrow(directed)),
        w = vertices[directed[, 2L]]),
        data.frame(v = vertices[bidirected[, 1L]], e = rep("<->", nrow(bidirected)),
          w = vertices[bidirected[, 2L]]))
      omitted <- vertices[runif(n) < .6]
      expect_identical(collapse_project_edges(vertices, edges, omitted),
        collapse_project_edges_R(vertices, edges, omitted))
    }
    duplicate <- c("A", "A", "B")
    edges <- data.frame(v = "A", e = "->", w = "B")
    expect_identical(collapse_project_edges(duplicate, edges, "A"),
      collapse_project_edges_R(duplicate, edges, "A"))
    invalid <- data.frame(v = "missing", e = "->", w = "B")
    message <- function(code) tryCatch({
      force(code)
      NULL
    }, error = conditionMessage)
    expect_identical(message(collapse_project_edges(c("A", "B"), invalid, "A")),
      message(collapse_project_edges_R(c("A", "B"), invalid, "A")))
  })
})

test_that("native contacts preserve exact curve plans and numerical fallback", {
  graphs <- c(
    'dag { A[pos="0,0"] B[pos="1,0"] C[pos="2,0"] D[pos="1,1"] A->C A->D B->D D->C A<->B }',
    'dag { A[pos="0,0"] B[pos="0,1"] C[pos="0,2"] D[pos="1,1"] A->C A->D B->D D->C A<->B }',
    'dag { A[pos="0,0"] B[pos="1,1"] C[pos="2,0"] D[pos="1,-1"] A->B B->C A->D D->C A->C B<->D }')
  for (text in graphs) for (mode in c("relative", "absolute"))
    for (position in c("default", "flipped", "over", "under")) {
      data <- tidy_ggdagitty(dagitty::dagitty(text))
      expected <- acceleration_mode(FALSE, curve_edge_plan(data, mode = mode, position = position))
      actual <- acceleration_mode(TRUE, curve_edge_plan(data, mode = mode, position = position))
      expect_identical(actual, expected)
    }
  symbol <- acceleration_mode(TRUE, .cf_native_symbol("C_cf_curve_contact"))
  expect_true(is.na(.Call(symbol, list(), list(), matrix(0, 2L, 2L))))
  poly <- function(points, from, to) list(from = from, to = to, lateral = TRUE,
    points = points, low = apply(points, 2, min), high = apply(points, 2, max),
    segment_low = pmin(points[-nrow(points), , drop = FALSE], points[-1L, , drop = FALSE]),
    segment_high = pmax(points[-nrow(points), , drop = FALSE], points[-1L, , drop = FALSE]))
  xy <- rbind(A = c(0, 0), B = c(1e308, 0), C = c(0, 1), D = c(1e308, 1))
  expect_true(is.na(.Call(symbol, poly(xy[c("A", "B"), ], "A", "B"),
    poly(xy[c("C", "D"), ], "C", "D"), xy)))
})

test_that("exact caches preserve all attributes, visibility and independent results", {
  acceleration_mode(TRUE, {
    .cf_cache$clear()
    on.exit(.cf_cache$clear(), add = TRUE)
    graph <- dagitty::dagitty("dag { X [exposure]; Y [outcome]; U -> X; U -> Y; X -> Y }")
    expected <- acceleration_mode(FALSE, .cd_prepare(graph, NULL, NULL))
    first <- .cd_prepare(graph, NULL, NULL)
    second <- .cd_prepare(graph, NULL, NULL)
    expect_identical(first, expected)
    expect_identical(second, expected)
    expect_gt(.cf_cache$stats()$counts[["hits"]], 0L)
    second$map$quantity[] <- "changed by caller"
    expect_identical(.cd_prepare(graph, NULL, NULL), expected)
    attr(graph, "source_revision") <- "new"
    expect_identical(.cd_prepare(graph, NULL, NULL),
      acceleration_mode(FALSE, .cd_prepare(graph, NULL, NULL)))
    attr(graph, "input_labels") <- c(U = "source_U")
    expect_identical(.gd_node_ids(graph), acceleration_mode(FALSE, .gd_node_ids(graph)))
    attr(graph, "cf_node_map") <- data.frame(node = names(graph),
      quantity = paste0("q_", names(graph)), occurrence = names(graph))
    expect_identical(.gd_node_ids(graph), acceleration_mode(FALSE, .gd_node_ids(graph)))
    for (i in 1:2) expect_false(withVisible(.cf_cache$call("metadata", list("invisible"),
      function() invisible(3L)))$visible)
    for (i in seq_len(12L)) .cf_cache$call("metadata", list(i), function() i)
    stats <- .cf_cache$stats()
    expect_lte(stats$entries_by_kind[["metadata"]], 8L)
    expect_lte(stats$entries_by_kind[["preparation"]], 2L)
    expect_lte(stats$accounted_bytes, stats$maximum_bytes)
  })
})

test_that("conditions and live references never disappear into a cache hit", {
  acceleration_mode(TRUE, {
    .cf_cache$clear()
    on.exit(.cf_cache$clear(), add = TRUE)
    calls <- 0L
    warn <- function() {
      calls <<- calls + 1L
      warning("cache warning")
      1L
    }
    for (i in 1:2) expect_warning(.cf_cache$call("metadata", list("warning"), warn), "cache warning")
    emit_message <- function() {
      calls <<- calls + 1L
      message("cache message")
      1L
    }
    for (i in 1:2) expect_message(.cf_cache$call("metadata", list("message"), emit_message), "cache message")
    for (i in 1:2) expect_error(.cf_cache$call("metadata", list("error"),
      function() {
        calls <<- calls + 1L
        stop("cache error")
      }), "cache error")
    for (i in 1:2) .cf_cache$call("metadata", list("condition"), function() {
      calls <<- calls + 1L
      signalCondition(simpleCondition("custom condition"))
      1L
    })
    expect_identical(calls, 8L)
    expect_identical(.cf_cache$stats()$entries, 0L)
    reference <- new.env(parent = emptyenv())
    reference$value <- 1L
    expect_identical(.cf_cache$call("metadata", list(reference), function() reference$value), 1L)
    reference$value <- 2L
    expect_identical(.cf_cache$call("metadata", list(reference), function() reference$value), 2L)
    expect_identical(.cf_cache$stats()$entries, 0L)
    graph <- dagitty::dagitty("dag { X [exposure]; Y [outcome]; X -> Y }")
    attr(graph, "live") <- reference
    expect_identical(.cd_prepare(graph, NULL, NULL),
      acceleration_mode(FALSE, .cd_prepare(graph, NULL, NULL)))
    class(graph) <- c("custom_graph", class(graph))
    expect_identical(.gd_node_ids(graph), acceleration_mode(FALSE, .gd_node_ids(graph)))
  })
})

test_that("mapping preserves quantity identity and exact temporal membership", {
  graphs <- lapply(c("dag { A -> B }", "dag { A_t1 -> B_t2 }", "dag { A_t2 -> B_t3 }"), dagitty::dagitty)
  mapped <- graphs[[1L]]
  attr(mapped, "cf_node_map") <- data.frame(node = c("A", "B"),
    quantity = c("A", "B"), occurrence = c("A_t2", "B_t3"))
  graphs[[4L]] <- mapped
  for (from in graphs) for (to in graphs)
    for (set in list(character(), "A", "A_t1", "A_t2", c("A", "B"), "unknown"))
      expect_identical(acceleration_mode(TRUE, .gd_map_set(set, from, to)),
        acceleration_mode(FALSE, .gd_map_set(set, from, to)))
})

test_that("an interrupted cache miss preserves retained entries and accounting", {
  acceleration_mode(TRUE, {
    .cf_cache$clear()
    on.exit(.cf_cache$clear(), add = TRUE)
    expect_identical(.cf_cache$call("metadata", list("retained"), function() 7L), 7L)
    before <- .cf_cache$stats()
    interrupted <- tryCatch(.cf_cache$call("metadata", list("interrupted"), function() {
      stop(structure(list(message = "cancelled", call = NULL), class = c("interrupt", "condition")))
    }), interrupt = function(e) TRUE)
    expect_true(interrupted)
    after <- .cf_cache$stats()
    expect_identical(after$entries, before$entries)
    expect_identical(after$accounted_bytes, before$accounted_bytes)
    expect_identical(.cf_cache$call("metadata", list("retained"),
      function() stop("retained entry lost")), 7L)
    .cf_cache$clear()
    expect_identical(.cf_cache$stats()$entries, 0L)
    expect_identical(.cf_cache$stats()$accounted_bytes, 0)
  })
})

test_that("locale changes invalidate preparation reuse", {
  acceleration_mode(TRUE, {
    .cf_cache$clear()
    on.exit(.cf_cache$clear(), add = TRUE)
    locale <- c(LC_CTYPE = "first encoding", LC_COLLATE = "first collation")
    expect_identical(.cf_cache_locale(), c(Sys.getlocale("LC_CTYPE"), Sys.getlocale("LC_COLLATE")))
    local_mocked_bindings(.cf_cache_locale = function() unname(locale))
    graph <- dagitty::dagitty("dag { X [exposure]; Y [outcome]; X -> Y }")
    first <- .cd_prepare(graph, NULL, NULL)
    expect_identical(.cd_prepare(graph, NULL, NULL), first)
    hits <- .cf_cache$stats()$counts[["hits"]]
    expect_gt(hits, 0L)
    locale[["LC_CTYPE"]] <- "changed encoding"
    misses <- .cf_cache$stats()$counts[["misses"]]
    expect_identical(.cd_prepare(graph, NULL, NULL), first)
    expect_gt(.cf_cache$stats()$counts[["misses"]], misses)
    locale[["LC_COLLATE"]] <- "changed collation"
    misses <- .cf_cache$stats()$counts[["misses"]]
    expect_identical(.cd_prepare(graph, NULL, NULL), first)
    expect_gt(.cf_cache$stats()$counts[["misses"]], misses)
  })
})

test_that("known query graphs reuse while custom nested classes bypass", {
  acceleration_mode(TRUE, {
    .cf_cache$clear()
    on.exit(.cf_cache$clear(), add = TRUE)
    graph <- dagitty::dagitty("dag { X [exposure]; Y [outcome]; X -> Y }")
    class(graph) <- c("dag_by_query", "dagitty")
    first <- .cd_prepare(graph, NULL, NULL)
    expect_identical(first, acceleration_mode(FALSE, .cd_prepare(graph, NULL, NULL)))
    expect_identical(.cd_prepare(graph, NULL, NULL), first)
    expect_gt(.cf_cache$stats()$counts[["hits"]], 0L)
    expect_identical(.cf_cache$stats()$entries_by_kind[["preparation"]], 1L)
    expect_false(.cf_standard_graph(structure(graph, class = c("custom", class(graph)))))
    .cf_cache$clear()
    labels <- structure(c(X = "source_X"), class = "external_labels")
    calls <- 0L
    read <- function() {
      calls <<- calls + 1L
      calls
    }
    for (i in 1:2) expect_identical(.cf_cache$call("metadata", list(list(labels = labels)), read), i)
    expect_identical(.cf_cache$stats()$entries, 0L)
    map <- structure(data.frame(node = "X", quantity = "source_X"),
      class = c("external_map", "data.frame"))
    for (i in 3:4) expect_identical(.cf_cache$call("metadata", list(list(map = map)), read), i)
    expect_identical(.cf_cache$stats()$entries, 0L)
    time <- structure(0, class = c("POSIXct", "POSIXt"))
    for (i in 5:6) expect_identical(.cf_cache$call("metadata", list(list(time = time)), read), i)
    expect_identical(.cf_cache$stats()$entries, 0L)
    attr(graph, "input_labels") <- labels
    expect_identical(.gd_node_ids(graph), acceleration_mode(FALSE, .gd_node_ids(graph)))
    expect_identical(.cf_cache$stats()$entries, 0L)
  })
})
