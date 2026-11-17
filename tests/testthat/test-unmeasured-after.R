# Declared unmeasured variables never change the ranked construction: it runs as
# for a fully observed graph and the declarations are marked on the finished DAG
# (author decision, 23 October 2026). ranked_query is copied from
# test-ranked-search.R, where it is local.
unmeasured_query <- function(edges, columns = NULL, node_map = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  nodes <- names(graph)
  ids <- .strip_slice(nodes)
  if (is.null(node_map)) node_map <- data.frame(node = nodes, quantity = ids, member_id = ids, occurrence = nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = columns, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = if (is.null(columns)) list() else stats::setNames(columns, columns))),
    node_map = node_map,
    policy = data.frame(node_id = unique(node_map$member_id), recovered_node = "yes"))
  graph
}

unmeasured_warnings <- function(expr) {
  warned <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warned)
}

test_that("a declared unmeasured confounder keeps the observed quantities and sizes and is marked", {
  graph <- unmeasured_query("U -> T; U -> Y; T -> Y")
  base <- construct_dag(graph, quiet = TRUE)
  run <- unmeasured_warnings(construct_dag(graph, latents = "U", quiet = TRUE))
  dag <- run$value
  expect_identical(dag_options(dag)$quantities, dag_options(base)$quantities)
  expect_identical(dag_options(dag)$sizes, dag_options(base)$sizes)
  expect_setequal(names(dag), names(base))
  expect_identical(dagitty::latents(dag), "U")
  expect_identical(dagitty::latents(temporal_dag(dag)), "U")
  expect_length(dagitty::latents(base), 0L)
  checks <- dag_checks(dag)
  expect_identical(checks$total$status, "preserved")
  expect_identical(checks$total$unmeasured, "U")
  expect_match(checks$total$availability, "^needs an unmeasured")
  expect_identical(dag_checks(base)$total$availability, "available now")
  expect_true(any(grepl("Every total adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: U).",
                        run$warnings, fixed = TRUE)))
  options <- dag_options(dag)
  expect_identical(options$unmeasured$semantics, "observed-as-if-v1")
  expect_identical(options$unmeasured$declared, "U")
  expect_identical(options$unmeasured$nodes, "U")
  expect_identical(options$unmeasured$display_nodes, "U")
  expect_true(attr(dag_checks(dag), "current"))
  expect_no_error(get_identification(dag, quiet = TRUE))
})

test_that("a declaration changes no size and the smallest size keeps the smaller set", {
  graph <- unmeasured_query("U -> T; U -> P; U -> Q; P -> Y; Q -> Y; T -> Y")
  for (size in c(0, Inf)) {
    base <- construct_dag(graph, size, quiet = TRUE)
    dag <- suppressWarnings(construct_dag(graph, size, latents = "U", quiet = TRUE))
    expect_setequal(names(dag), names(base))
    expect_identical(dag_options(dag)$quantities, dag_options(base)$quantities)
    expect_identical(dagitty::latents(dag), "U")
  }
  expect_setequal(names(construct_dag(graph, 0, quiet = TRUE)), c("T", "Y", "U"))
  expect_setequal(names(construct_dag(graph, quiet = TRUE)), c("T", "Y", "U", "P", "Q"))
})

test_that("a declared unmeasured mediator keeps its direct quantity and is marked", {
  graph <- unmeasured_query("T -> M -> Y; T -> Y")
  run <- unmeasured_warnings(construct_dag(graph, latents = "M", quiet = TRUE))
  dag <- run$value
  expect_identical(dag_options(dag)$quantities$direct$state, "sets")
  expect_identical(dag_options(dag)$quantities$direct$members, "M")
  expect_identical(dag_checks(dag)$direct$status, "preserved")
  expect_identical(dagitty::latents(dag), "M")
  expect_identical(dag_checks(dag)$direct$availability, "needs an unmeasured variable")
  expect_identical(dag_checks(dag)$total$availability, "available now")
  expect_identical(run$warnings,
    "Every direct adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: M).")
})

test_that("a folded display node standing for a declared occasion is marked and recorded as partial", {
  graph <- unmeasured_query("A_t0 -> T; A_t0 -> Y; A_t0 -> A_t1; A_t1 -> T; A_t1 -> Y; T -> Y")
  run <- unmeasured_warnings(construct_dag(graph, latents = "A_t0", quiet = TRUE))
  dag <- run$value
  expect_true(dag_options(dag)$summary$applied)
  expect_setequal(names(dag), c("T", "A", "Y"))
  expect_identical(dag_options(dag)$quantities$total$members, "A")
  expect_identical(dagitty::latents(dag), "A")
  expect_identical(dagitty::latents(temporal_dag(dag)), "A_t0")
  partial <- dag_options(dag)$unmeasured$partial
  expect_identical(partial, list(A = list(unmeasured = "A_t0", observed = "A_t1")))
  expect_true(any(grepl("A stands for both observed and unmeasured occasions (unmeasured: A_t0, observed: A_t1). It is shown unmeasured.",
                        run$warnings, fixed = TRUE)))
  expect_identical(dag_checks(dag)$total$availability, "needs an unmeasured variable")
  whole <- suppressWarnings(construct_dag(graph, latents = "A", quiet = TRUE))
  expect_identical(dagitty::latents(whole), "A")
  expect_setequal(dagitty::latents(temporal_dag(whole)), c("A_t0", "A_t1"))
  expect_length(dag_options(whole)$unmeasured$partial, 0L)
})

test_that("a fold with one declared occasion is marked and a declared member marks its group node", {
  feedback <- unmeasured_query(paste("A_t2 -> T; A_t2 -> Y; B_t2 -> Y; T -> Y; A_t0 -> B_t1; A_t1 -> B_t2;",
                                     "B_t0 -> A_t2; A_t0 -> A_t1; A_t1 -> A_t2; B_t0 -> B_t1; B_t1 -> B_t2"))
  run <- unmeasured_warnings(construct_dag(feedback, 0, latents = "A_t0", quiet = TRUE))
  dag <- run$value
  expect_identical(dagitty::latents(dag), "A")
  expect_identical(dag_options(dag)$unmeasured$partial, list(A = list(unmeasured = "A_t0", observed = c("A_t1", "A_t2"))))
  expect_identical(dag_checks(dag)$total[c("availability", "temporal_availability")],
    list(availability = "needs an unmeasured variable", temporal_availability = "available now"))
  expect_true(any(grepl("temporal_dag() lists a set that uses none", run$warnings, fixed = TRUE)))
  expect_length(.gd_adjustment_sets(dag, "total", Inf)$sets, 0L)
  observed <- dag
  dagitty::latents(observed) <- character()
  observed_sets <- .gd_adjustment_sets(observed, "total", Inf)$sets
  expect_gt(length(observed_sets), 0L)
  expect_true(all(vapply(observed_sets, function(set) "A" %in% set, logical(1))))
  temporal_sets <- .gd_adjustment_sets(temporal_dag(dag), "total", Inf)$sets
  expect_gt(length(temporal_sets), 0L)
  expect_true(any(vapply(temporal_sets, function(set) !"A_t0" %in% set, logical(1))))
  expect_true(all(.gd_validate_sets(.cd_prepare(feedback, NULL, NULL)$selection, temporal_sets, "total")))
  map <- data.frame(node = c("T", "Y", "G", "G", "C"), quantity = c("T", "Y", "G", "G", "C"),
                    member_id = c("T", "Y", "A", "B", "C"), occurrence = c("T", "Y", "G", "G", "C"))
  group <- unmeasured_query("G -> T; G -> Y; T -> Y; C -> T; C -> Y", node_map = map)
  run <- unmeasured_warnings(construct_dag(group, 0, latents = "A", quiet = TRUE))
  dag <- run$value
  expect_identical(dag_options(dag)$unmeasured$nodes, "G")
  expect_identical(dagitty::latents(dag), "G")
  expect_true(any(grepl("declared unmeasured variable", run$warnings, fixed = TRUE)))
  expect_error(construct_dag(group, latents = "A", quiet = TRUE), "stands for several source members")
})

test_that("a requested declared unmeasured variable is retained and marked", {
  graph <- unmeasured_query("T -> M -> Y; T -> Y; V -> Y")
  dag <- construct_dag(graph, latents = "V", keep_nodes = "V", quiet = TRUE)
  expect_true("V" %in% names(dag))
  expect_identical(dagitty::latents(dag), "V")
  expect_identical(vapply(dag_checks(dag), `[[`, character(1), "availability"),
                   c(direct = "available now", total = "available now"))
  unrequested <- construct_dag(graph, 0, latents = "V", quiet = TRUE)
  expect_false("V" %in% names(unrequested))
  expect_length(dagitty::latents(unrequested), 0L)
  expect_length(dag_options(unrequested)$unmeasured$display_nodes, 0L)
})


test_that("a family blocked by a bidirected edge still has no set when there are no declarations", {
  graph <- unmeasured_query("T <-> Y; T -> Y")
  dag <- NULL
  expect_warning(dag <- construct_dag(graph, quiet = TRUE), "No eligible")
  expect_identical(vapply(dag_options(dag)$quantities, `[[`, "", "state"), c(direct = "no_sets", total = "no_sets"))
  expect_identical(dag_checks(dag)$total$status, "no eligible reference set")
  expect_identical(vapply(dag_checks(dag), `[[`, character(1), "availability"),
                   c(direct = "no adjustment set", total = "no adjustment set"))
  expect_length(dag_options(dag)$unmeasured$nodes, 0L)
  empty <- construct_dag(unmeasured_query("T -> Y"), quiet = TRUE)
  expect_identical(dag_options(empty)$quantities$total$state, "empty_set")
  expect_identical(dag_checks(empty)$total$availability, "available now")
})

test_that("a budget stop without a set is never reported as no set", {
  withr::local_options(causaliflower.ranked_steps = 1L)
  graph <- unmeasured_query("A -> T; A -> Y; T -> Y")
  testthat::local_mocked_bindings(query_graph = function(...) graph, .package = "causaliflower")
  run <- unmeasured_warnings(get_dag("barley", "T", "Y", TRUE, quiet = TRUE))
  dag <- run$value
  expect_identical(vapply(dag_options(dag)$quantities, `[[`, "", "state"), c(direct = "not_found", total = "not_found"))
  expect_identical(vapply(dag_checks(dag), `[[`, "", "status"), c(direct = "not fully checked", total = "not fully checked"))
  expect_true(dag_options(dag)$identification$applied)
  listed <- .gd_adjustment_sets(dag, "total", Inf)$sets
  expect_gt(length(listed), 0L)
  expect_true(all(.gd_validate_sets(unmeasured_query("A -> T; A -> Y; T -> Y"), listed, "total")))
  expect_identical(dag_checks(dag)$total$availability, "available now")
  expect_true(any(grepl("not evidence that no set exists", run$warnings, fixed = TRUE)))
  expect_false(any(grepl("No eligible", run$warnings, fixed = TRUE)))
})
