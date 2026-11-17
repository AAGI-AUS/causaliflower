query_signature_graph <- function() {
  dagitty::dagitty("dag { T [exposure]; Y [outcome]; P -> T; P -> Y; Q -> Y; T -> Y }")
}

test_that("request-local signatures preserve the exact saved binding and cache key", {
  graph <- query_signature_graph()
  query <- .cd_size_query(graph, NULL)
  expect_identical(query$key, list(.cd_construction_version,
    .cd_search_identity(query$observed), .gd_signature(query$observed$analysis), NULL))
  legacy <- query[c("observed", "key")]
  expected <- .cd_choice_binding(graph, legacy, "Q", "P", 20000L)
  read_signature <- .gd_signature
  calls <- 0L
  testthat::local_mocked_bindings(.gd_signature = function(graph) {
    calls <<- calls + 1L
    read_signature(graph)
  }, .package = "causaliflower")
  expect_identical(.cd_choice_binding(graph, query, "Q", "P", 20000L), expected)
  expect_identical(calls, 1L)
  calls <- 0L
  expect_identical(.cd_choice_binding(graph, legacy, "Q", "P", 20000L), expected)
  expect_identical(calls, 3L)
})

test_that("changed prepared graph inputs invalidate only their saved signature", {
  graph <- query_signature_graph()
  query <- .cd_size_query(graph, NULL)
  read_signature <- .gd_signature
  calls <- 0L
  testthat::local_mocked_bindings(.gd_signature = function(graph) {
    calls <<- calls + 1L
    read_signature(graph)
  }, .package = "causaliflower")
  for (name in c("analysis", "selection")) {
    for (field in c("text", "roles", "labels", "class")) {
      changed <- query
      dag <- changed$observed[[name]]
      if (field == "text") {
        dag[[1L]] <- as.character(dagitty::dagitty(
          "dag { T [exposure]; Y [outcome]; P -> T; Q -> Y; T -> Y }"))
      } else if (field == "roles") {
        dagitty::exposures(dag) <- "P"
      } else if (field == "labels") {
        attr(dag, "input_labels") <- c(P = "a different label")
      } else {
        class(dag) <- c("cf_dag", "dagitty")
      }
      changed$observed[[name]] <- dag
      expected <- .cd_choice_binding(graph, changed[c("observed", "key")], NULL, NULL, 20000L)
      calls <- 0L
      expect_identical(.cd_choice_binding(graph, changed, NULL, NULL, 20000L), expected,
        info = paste(name, field))
      expect_identical(calls, 2L, info = paste(name, field))
    }
  }
})

test_that("missing identities and changed locales retain fresh signature calculation", {
  graph <- query_signature_graph()
  query <- .cd_size_query(graph, NULL)
  expected <- .cd_choice_binding(graph, query[c("observed", "key")], NULL, NULL, 20000L)
  read_signature <- .gd_signature
  calls <- 0L
  testthat::local_mocked_bindings(.gd_signature = function(graph) {
    calls <<- calls + 1L
    read_signature(graph)
  }, .package = "causaliflower")
  for (field in seq_along(query$signatures$locale)) {
    calls <- 0L
    local({
      locale <- query$signatures$locale
      locale[[field]] <- paste0(locale[[field]], "-changed")
      testthat::local_mocked_bindings(.cf_cache_locale = function() locale,
        .package = "causaliflower")
      expect_identical(.cd_choice_binding(graph, query, NULL, NULL, 20000L), expected)
      expect_identical(calls, 3L)
    })
  }
  testthat::local_mocked_bindings(.cd_search_identity = function(...) NULL,
    .package = "causaliflower")
  missing <- .cd_size_query(graph, NULL)
  expect_null(missing$key)
  expect_null(missing$signatures)
  calls <- 0L
  expect_identical(.cd_choice_binding(graph, missing, NULL, NULL, 20000L), expected)
  expect_identical(calls, 3L)
})

test_that("signature reuse keeps original-source and variable-map checks fresh", {
  graph <- query_signature_graph()
  query <- .cd_size_query(graph, NULL)
  before <- .cd_choice_binding(graph, query, NULL, NULL, 20000L)
  changed <- graph
  dagitty::latents(changed) <- "Q"
  after <- .cd_choice_binding(changed, query, NULL, NULL, 20000L)
  expect_false(identical(after, before))
  expect_identical(after,
    .cd_choice_binding(changed, query[c("observed", "key")], NULL, NULL, 20000L))

  changed_query <- query
  map <- attr(changed_query$observed$analysis, "cf_node_map", exact = TRUE)
  map$member_id[map$node == "P"] <- "different-P"
  attr(changed_query$observed$analysis, "cf_node_map") <- map
  after <- .cd_choice_binding(graph, changed_query, NULL, NULL, 20000L)
  expect_false(identical(after, before))
  expect_identical(after, .cd_choice_binding(graph,
    changed_query[c("observed", "key")], NULL, NULL, 20000L))
  attr(changed_query$observed$analysis, "cf_node_map") <- map[FALSE, ]
  expect_error(.cd_choice_binding(graph, changed_query, NULL, NULL, 20000L),
    "source identity map has missing nodes or invalid values")
  attr(changed, "cf_node_map") <- data.frame(node = "P")
  expect_error(.cd_choice_binding(changed, query, NULL, NULL, 20000L),
    "source identity map is incomplete or invalid")
})
