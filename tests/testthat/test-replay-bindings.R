test_that("saved selections bind the original source, occurrences and requests", {
  query <- function(edges) {
    graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
    nodes <- names(graph)
    attr(graph, "query_record") <- list(version = 1L, policy_applied = TRUE,
      arguments = list(crop = "fixture", params = list(query_depth = Inf)),
      requests = list(columns = list(resolved = list())),
      node_map = data.frame(node = nodes, quantity = nodes, member_id = nodes, occurrence = nodes),
      policy = data.frame(node_id = nodes, recovered_node = ifelse(nodes == "B", "no:fixture", "yes")))
    graph
  }
  graph <- query("A -> B -> T; A -> T; A -> Y; T -> Y")
  record <- attr(graph, "query_record")
  dag <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  saved <- dag_call(dag)
  expect_identical(.gd_signature(suppressWarnings(eval(saved))), .gd_signature(dag))

  changed <- query("A -> B -> T; A -> Y; T -> Y")
  expect_false(identical(.gd_signature(graph), .gd_signature(changed)))
  expect_identical(.gd_signature(.cd_prepare(graph, NULL, NULL)$analysis),
    .gd_signature(.cd_prepare(changed, NULL, NULL)$analysis))
  replay <- saved
  replay$graph <- changed
  expect_error(eval(replay), "source or construction settings no longer match")

  for (field in c("occurrence", "requests")) {
    changed <- graph
    changed_record <- record
    if (field == "occurrence") {
      changed_record$node_map$occurrence[changed_record$node_map$node == "A"] <- "A_other"
    } else {
      changed_record$requests$columns$raw <- "a different original request"
    }
    attr(changed, "query_record") <- changed_record
    replay <- saved
    replay$graph <- changed
    expect_error(eval(replay), "source or construction settings no longer match")
    if (field == "occurrence") {
      state <- .gd_state(dag)
      state$original_query <- changed
      modified <- .gd_attach(dag, state)
      expect_error(suppressWarnings(select_dag(modified, quiet = TRUE)),
        "does not exactly match an available alternative")
    }
  }
})

test_that("saved selections bind effective attributes without a query receipt", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> Y; T -> Y }")
  nodes <- names(graph)
  ranks <- seq_along(nodes)
  source_order <- list(version = 1L, ranks = data.frame(node_id = nodes,
    base_order = ranks, at_order = ranks, post_order = ranks, pre_order = ranks,
    time_order = ranks))
  attr(graph, "cf_source_order") <- source_order
  attr(graph, "query_depth") <- 2
  dag <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  saved <- dag_call(dag)
  expect_identical(.gd_signature(suppressWarnings(eval(saved))), .gd_signature(dag))

  changed <- graph
  changed_order <- source_order
  changed_order$ranks[-1L] <- lapply(changed_order$ranks[-1L], rev)
  attr(changed, "cf_source_order") <- changed_order
  replay <- saved
  replay$graph <- changed
  expect_error(suppressWarnings(eval(replay)), "source or construction settings no longer match")

  changed <- graph
  attr(changed, "query_depth") <- 3
  replay <- saved
  replay$graph <- changed
  expect_error(suppressWarnings(eval(replay)), "source or construction settings no longer match")
})
