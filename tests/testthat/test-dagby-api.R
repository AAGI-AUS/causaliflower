dagby_api_query <- function(edges = "T -> A -> Y; T -> C; T -> Y", columns = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  nodes <- names(graph)
  quantities <- sub("_t[0-9]+$", "", nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = columns, params = dagby_query_params()),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = quantities, member_id = quantities,
      occurrence = nodes),
    policy = data.frame(node_id = unique(quantities), recovered_node = "yes"))
  return(graph)
}

test_that("the dagby family replaces the development names without aliases", {
  public <- c("dagby", "dagby_query", "dagby_construct", "dagby_query_params",
    "dagby_prune", "dagby_grow", "dagby_select", "dagby_options", "dagby_checks",
    "dagby_call", "dagby_source", "dagby_temporal")
  # These old spellings are intentional assertions of their removal
  retired <- c("get_dag", "query_graph", "construct_dag", "query_params", "prune_dag",
    "grow_dag", "select_dag", "dag_options", "dag_checks", "dag_call", "original_query", "temporal_dag")
  exported <- getNamespaceExports("causaliflower")
  expect_true(all(public %in% exported))
  expect_false(any(retired %in% exported))
  for (name in public) expect_true(is.function(getExportedValue("causaliflower", name)), info = name)
  for (name in retired) expect_false(exists(name, envir = asNamespace("causaliflower"),
    inherits = FALSE), info = name)
})

test_that("recipe execution resolves the query backend when the query stage runs", {
  queried <- 0L
  query <- function(crop, treatments = NULL, outcomes = NULL, randomised_treatment = FALSE,
                    columns = NULL, time_unroll = NULL, label = NULL, quiet = FALSE,
                    params = dagby_query_params()) {
    queried <<- queried + 1L
    return(dagby_api_query(columns = columns))
  }
  recipe <- dagby("fixture", "T", "Y", TRUE, columns = "C", size = 0, keep_nodes = "C")
  expect_identical(queried, 0L)
  testthat::local_mocked_bindings(.dagby_query_source = query, .package = "causaliflower")
  dag <- recipe |> dagby_query(quiet = TRUE) |> dagby_construct(quiet = TRUE)
  expect_identical(queried, 1L)
  expect_setequal(names(dag), c("T", "A", "C", "Y"))
  resolved <- dagby_call(dag)$input_record$resolved_call
  expect_identical(resolved[[1L]], quote(causaliflower::dagby_query))
  expect_identical(resolved[[2L]]$columns, "C")
  expect_identical(resolved[[2L]]$size, 0)
  expect_identical(resolved[[2L]]$keep_nodes, "C")
  expect_false(any(c("query_graph", "construct_dag") %in% names(resolved)))
  caller <- new.env(parent = globalenv())
  caller$dagby <- function(...) stop("A caller binding must not replace the recorded function.")
  reproduced <- dagby_construct(eval(resolved, caller), quiet = TRUE)
  expect_identical(.gd_signature(reproduced), .gd_signature(dag))
})

test_that("saved construction recipes use the public namespace and valid arguments", {
  graph <- dagby_api_query()
  saved <- .gd_saved_call(list(graph = graph, size = 0, quiet = TRUE, retired = "discard"))
  expect_identical(saved[[1L]], quote(causaliflower::dagby_construct))
  expect_false("retired" %in% names(saved))
  caller <- new.env(parent = globalenv())
  caller$dagby_construct <- function(...) stop("A caller binding must not replace the recorded function.")
  reproduced <- eval(saved, caller)
  expect_setequal(names(reproduced), c("T", "A", "Y"))
  expect_identical(dagby_source(reproduced), graph)
})

test_that("renamed exact replay retains source checks and rejects changed source", {
  graph <- dagby_api_query()
  dag <- dagby_construct(graph, 0, quiet = TRUE)
  saved <- dagby_call(dag)
  expect_identical(saved[[1L]], quote(causaliflower:::.cd_replay))
  reproduced <- eval(saved)
  expect_identical(dagby_source(reproduced), graph)
  expect_identical(dagby_temporal(reproduced), dagby_temporal(dag))
  expect_identical(dagby_checks(reproduced), dagby_checks(dag))
  saved$graph <- dagby_api_query("T -> A -> Y; T -> C; T -> Y; C -> Y")
  expect_error(eval(saved), "source or construction settings no longer match", fixed = TRUE)
})

test_that("size guidance and static retries use the renamed operations", {
  messages <- function(expr) {
    heard <- character()
    withCallingHandlers(expr, message = function(condition) {
      heard <<- c(heard, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
    return(paste(heard, collapse = "\n"))
  }
  graph <- dagby_api_query()
  expect_match(messages(dagby_construct(graph, 0)), "dagby_grow(dag, 1)", fixed = TRUE)
  expect_match(messages(dagby_construct(graph, Inf)), "dagby_prune(dag, 1)", fixed = TRUE)
  temporal <- dagby_api_query(paste("Q_t1 -> T; Q_t1 -> Z_t1; T -> Z_t1;",
    "Z_t1 -> Q_t2; Q_t2 -> Y; Q_t2 -> Z_t2; T -> Y"), columns = "Z")
  expect_error(suppressWarnings(dagby_construct(temporal, quiet = TRUE, size = "largest")),
    "Use dagby_construct(graph, size = ", fixed = TRUE, class = "cd_static_unavailable")
  dag <- suppressWarnings(dagby_construct(temporal, unroll = TRUE, quiet = TRUE, size = "largest"))
  expect_error(suppressWarnings(dagby_select(dag, unroll = FALSE, quiet = TRUE)),
    "Use dagby_select(dag, size = ", fixed = TRUE, class = "cd_static_unavailable")
})
