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
  class(graph) <- c("dag_by_query", class(graph))
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

test_that("a declared unmeasured confounder keeps the observed catalogue and selection and is marked", {
  graph <- unmeasured_query("U -> T; U -> Y; T -> Y")
  base <- construct_dag(graph, quiet = TRUE)
  run <- unmeasured_warnings(construct_dag(graph, latents = "U", quiet = TRUE))
  dag <- run$value
  # Search, certificates, steps and catalogue equal the latents = NULL run.
  expect_identical(dag_options(dag)$search_identity, dag_options(base)$search_identity)
  for (effect in c("direct", "total")) {
    declared <- dag_options(dag)[[effect]]
    observed <- dag_options(base)[[effect]]
    expect_identical(declared$label, observed$label, info = effect)
    expect_identical(declared$raw_sets, observed$raw_sets, info = effect)
    expect_identical(declared[setdiff(names(declared), "unmeasured")],
                     observed[setdiff(names(observed), "unmeasured")], info = effect)
    expect_identical(dag_options(dag)$search[[effect]]$sets, dag_options(base)$search[[effect]]$sets, info = effect)
    expect_identical(dag_options(dag)$search[[effect]]$ranked, dag_options(base)$search[[effect]]$ranked, info = effect)
  }
  expect_identical(dag_options(dag)$selected, dag_options(base)$selected)
  expect_identical(dag_checks(dag)$total$required, list("U"))
  # The declaration is marked on the finished DAG and its temporal graph only.
  expect_identical(dagitty::latents(dag), "U")
  expect_identical(dagitty::latents(temporal_dag(dag)), "U")
  expect_length(dagitty::latents(base), 0L)
  checks <- dag_checks(dag)
  expect_identical(checks$total$status, "preserved")
  expect_identical(checks$total$unmeasured, "U")
  expect_match(checks$total$availability, "^needs unmeasured")
  expect_identical(checks$total$records[[1L]]$unmeasured, "U")
  expect_identical(dag_checks(base)$total$availability, "available now")
  expect_true(any(grepl("The selected total adjustment set relies on declared unmeasured variable(s) U",
                        run$warnings, fixed = TRUE)))
  options <- dag_options(dag)
  expect_identical(options$unmeasured$semantics, "observed-as-if-v1")
  expect_identical(options$unmeasured$declared, "U")
  expect_identical(options$unmeasured$nodes, "U")
  expect_identical(options$unmeasured$members, "U")
  expect_identical(options$unmeasured$display_nodes, "U")
  expect_identical(options$unmeasured$selected, list(direct = "U", total = "U"))
  expect_identical(options$total$unmeasured[[1L]], "U")
  expect_identical(.gd_state(dag)$build$availability$display_latent, "U")
  # Nothing is set before the search. Marking is the last step.
  history <- dag_history(dag)$operation
  expect_identical(tail(history, 1L), "mark_unmeasured")
  expect_false("set_latents" %in% history)
  mark <- tail(dag_history(dag, details = TRUE), 1L)[[1L]]
  expect_identical(mark$delta$statuses$latent$added, "U")
  expect_identical(mark$details$availability$total, "needs unmeasured: U")
  expect_identical(.gd_signature(mark$after), .gd_signature(dag))
  expect_true(attr(dag_checks(dag), "current"))
  expect_no_error(get_identification(dag, quiet = TRUE))
})

test_that("options that need an unmeasured variable stay listed and the automatic choice is unchanged", {
  graph <- unmeasured_query("U -> T; U -> P; U -> Q; P -> Y; Q -> Y; T -> Y")
  base <- construct_dag(graph, quiet = TRUE)
  dag <- suppressWarnings(construct_dag(graph, latents = "U", quiet = TRUE))
  for (effect in c("direct", "total")) {
    expect_identical(dag_options(dag)[[effect]]$label, c("U", "P + Q"), info = effect)
    expect_identical(dag_options(dag)[[effect]]$label, dag_options(base)[[effect]]$label, info = effect)
    expect_identical(dag_options(dag)[[effect]]$unmeasured, list("U", character()), info = effect)
  }
  expect_identical(dag_options(dag)$selected, dag_options(base)$selected)
  expect_identical(dag_options(dag)$selected$total, c(1L, 1L))
  # An explicit choice of {P, Q} in both families builds that DAG without U.
  pq <- construct_dag(graph, latents = "U", direct_select = c(2L, 1L), total_select = c(2L, 1L), quiet = TRUE)
  expect_false("U" %in% names(pq))
  expect_length(dagitty::latents(pq), 0L)
  expect_identical(dag_checks(pq)$total$status, "preserved")
  expect_identical(dag_checks(pq)$total$availability, "available now")
  expect_identical(dag_options(pq)$unmeasured$selected, list(direct = character(), total = character()))
  printed <- capture.output(print(dag_options(dag)))
  needs <- grep("needs unmeasured", printed, value = TRUE)
  expect_length(needs, 2L)
  expect_true(all(grepl("{U} needs unmeasured {U}", needs, fixed = TRUE)))
  expect_true(any(grepl("Declared unmeasured variables: U (marked on the finished DAG", printed, fixed = TRUE)))
  expect_false(any(grepl("unmeasured", capture.output(print(dag_options(base))), fixed = TRUE)))
  # Deferred options carry the same annotations before any DAG exists.
  deferred <- construct_dag(graph, latents = "U", direct_select = NULL, quiet = TRUE)
  expect_s3_class(deferred, "cf_dag_options")
  expect_identical(deferred$total$unmeasured, list("U", character()))
  expect_identical(deferred$unmeasured$nodes, "U")
  expect_length(deferred$unmeasured$display_nodes, 0L)
})

test_that("a declared unmeasured mediator keeps its direct option and is marked", {
  graph <- unmeasured_query("T -> M -> Y; T -> Y")
  run <- unmeasured_warnings(construct_dag(graph, latents = "M", quiet = TRUE))
  dag <- run$value
  expect_identical(dag_options(dag)$status[["direct"]], "selected")
  expect_identical(dag_options(dag)$direct$label, "M")
  expect_identical(dag_checks(dag)$direct$required, list("M"))
  expect_identical(dag_checks(dag)$direct$status, "preserved")
  expect_identical(dagitty::latents(dag), "M")
  expect_identical(dag_checks(dag)$direct$availability, "needs unmeasured: M")
  expect_identical(dag_checks(dag)$total$availability, "available now")
  expect_identical(run$warnings,
    "The selected direct adjustment set relies on declared unmeasured variable(s) M. They are marked on the DAG.")
})

test_that("a folded display node standing for a declared occasion is marked and recorded as partial", {
  graph <- unmeasured_query("A_t0 -> T; A_t0 -> Y; A_t0 -> A_t1; A_t1 -> T; A_t1 -> Y; T -> Y")
  run <- unmeasured_warnings(construct_dag(graph, latents = "A_t0", quiet = TRUE))
  dag <- run$value
  expect_true(.gd_state(dag)$summary$applied)
  expect_setequal(names(dag), c("T", "A", "Y"))
  expect_identical(dag_checks(dag)$total$required, list(c("A_t0", "A_t1")))
  expect_identical(dagitty::latents(dag), "A")
  expect_identical(dagitty::latents(temporal_dag(dag)), "A_t0")
  partial <- dag_options(dag)$unmeasured$partial
  expect_identical(partial, list(A = list(unmeasured = "A_t0", observed = "A_t1")))
  expect_identical(.gd_state(dag)$build$availability$partial, partial)
  expect_true(any(grepl("A stands for both observed and unmeasured occasions (unmeasured: A_t0, observed: A_t1). It is shown unmeasured.",
                        run$warnings, fixed = TRUE)))
  expect_identical(dag_checks(dag)$total$availability, "needs unmeasured: A")
  # The catalogue and the warning name the declared occasion the realised set
  # uses. The checks read the folded display node.
  expect_identical(dag_options(dag)$total$unmeasured, list("A_t0"))
  expect_identical(dag_options(dag)$unmeasured$selected, list(direct = "A_t0", total = "A_t0"))
  expect_true(any(grepl("total adjustment set relies on declared unmeasured variable(s) A_t0.", run$warnings, fixed = TRUE)))
  mark <- tail(dag_history(dag, details = TRUE), 1L)[[1L]]
  expect_identical(mark$details$nodes, "A")
  expect_identical(mark$details$temporal_nodes, "A_t0")
  expect_identical(mark$details$partial, partial)
  # A view of the same receipt applies the same rule.
  expect_warning(view <- dag_set_view(dag, "total", 1L), "stands for both observed and unmeasured occasions")
  expect_identical(dagitty::latents(view), "A")
  expect_identical(dag_checks(view)$total$unmeasured, "A")
  # Declaring the whole quantity marks it without a partial record.
  whole <- suppressWarnings(construct_dag(graph, latents = "A", quiet = TRUE))
  expect_identical(dagitty::latents(whole), "A")
  expect_setequal(dagitty::latents(temporal_dag(whole)), c("A_t0", "A_t1"))
  expect_length(dag_options(whole)$unmeasured$partial, 0L)
})

test_that("the catalogue and the warning name the declared nodes the realised sets use", {
  # The declared A_t0 is not used: the realised sets use the observed A_t2.
  feedback <- unmeasured_query(paste("A_t2 -> T; A_t2 -> Y; B_t2 -> Y; T -> Y; A_t0 -> B_t1; A_t1 -> B_t2;",
                                     "B_t0 -> A_t2; A_t0 -> A_t1; A_t1 -> A_t2; B_t0 -> B_t1; B_t1 -> B_t2"))
  run <- unmeasured_warnings(construct_dag(feedback, latents = "A_t0", quiet = TRUE))
  dag <- run$value
  expect_identical(dag_checks(dag)$total$required, list("A_t2"))
  expect_false(any(grepl("relies on declared unmeasured", run$warnings, fixed = TRUE)))
  expect_identical(dag_options(dag)$unmeasured$selected, list(direct = character(), total = character()))
  expect_identical(dag_options(dag)$total$unmeasured, list(character()))
  expect_false(any(grepl("needs unmeasured", capture.output(print(dag_options(dag))), fixed = TRUE)))
  expect_identical(vapply(dag_checks(dag), `[[`, character(1), "availability"),
                   c(direct = "available now", total = "available now"))
  # A declared member marks its recorded group node. The column names that node.
  map <- data.frame(node = c("T", "Y", "G", "G", "C"), quantity = c("T", "Y", "G", "G", "C"),
                    member_id = c("T", "Y", "A", "B", "C"), occurrence = c("T", "Y", "G", "G", "C"))
  group <- unmeasured_query("G -> T; G -> Y; T -> Y; C -> T; C -> Y", node_map = map)
  run <- unmeasured_warnings(construct_dag(group, latents = "A", quiet = TRUE))
  dag <- run$value
  expect_identical(dag_options(dag)$unmeasured$nodes, "G")
  expect_identical(dag_options(dag)$unmeasured$members, c("A", "B"))
  expect_identical(dag_options(dag)$total$unmeasured, list("G"))
  expect_identical(dag_options(dag)$unmeasured$selected, list(direct = "G", total = "G"))
  expect_identical(dagitty::latents(dag), "G")
  expect_true(any(grepl("total adjustment set relies on declared unmeasured variable(s) G.", run$warnings, fixed = TRUE)))
})

test_that("a requested declared unmeasured variable is retained and marked", {
  graph <- unmeasured_query("T -> M -> Y; T -> Y; V -> Y")
  dag <- construct_dag(graph, latents = "V", keep_nodes = "V", quiet = TRUE)
  expect_true("V" %in% names(dag))
  expect_identical(dagitty::latents(dag), "V")
  expect_identical(dag_options(dag)$unmeasured$selected, list(direct = character(), total = character()))
  expect_identical(vapply(dag_checks(dag), `[[`, character(1), "availability"),
                   c(direct = "available now", total = "available now"))
  unrequested <- construct_dag(graph, latents = "V", quiet = TRUE)
  expect_false("V" %in% names(unrequested))
  expect_length(dagitty::latents(unrequested), 0L)
  expect_identical(dag_options(unrequested)$unmeasured$nodes, "V")
  expect_length(dag_options(unrequested)$unmeasured$display_nodes, 0L)
})

test_that("declarations reuse the observed ranked search and older declared receipts are refused", {
  calls <- 0L
  ranked_search <- .cd_ranked_search
  testthat::local_mocked_bindings(.cd_ranked_search = function(...) {
    calls <<- calls + 1L
    ranked_search(...)
  })
  graph <- unmeasured_query("U -> T; U -> P; U -> Q; P -> Y; Q -> Y; T -> Y")
  base <- construct_dag(graph, quiet = TRUE)
  expect_identical(calls, 2L)
  calls <- 0L
  marked <- suppressWarnings(.cd_reconstruct(graph, previous = base, latents = "U", quiet = TRUE))
  expect_identical(calls, 0L)
  expect_identical(dagitty::latents(marked), "U")
  expect_identical(marked, suppressWarnings(construct_dag(graph, latents = "U", quiet = TRUE)))
  # The reverse direction reuses the declared receipt's search and marks nothing.
  calls <- 0L
  plain <- .cd_reconstruct(graph, previous = marked, quiet = TRUE)
  expect_identical(calls, 0L)
  expect_length(dagitty::latents(plain), 0L)
  expect_identical(.gd_signature(plain), .gd_signature(base))
  # select_dag replays the receipt's declarations on the same catalogue indexes.
  calls <- 0L
  switched <- suppressWarnings(select_dag(marked, total_select = c(2L, 1L)))
  expect_identical(calls, 0L)
  expect_identical(dag_options(switched)$selected$total, c(2L, 1L))
  expect_identical(dag_options(switched)$unmeasured$nodes, "U")
  expect_identical(dagitty::latents(switched), "U")
  expect_identical(dag_checks(switched)$total$availability, "available now")
  expect_identical(tail(dag_history(switched)$operation, 1L), "select_dag")
  # A ranked receipt with declarations saved before this change is refused.
  state <- .gd_state(marked)
  state$construction$unmeasured <- NULL
  old <- .gd_attach(marked, state)
  expect_error(select_dag(old, total_select = c(2L, 1L)),
               "saved before declared unmeasured variables stopped changing the search", fixed = TRUE)
  expect_warning(saved <- dag_call(old),
                 "predates the rule that declared unmeasured variables do not change the search", fixed = TRUE)
  expect_true(is.call(saved))
  expect_silent(dag_call(marked))
  # Receipts without declarations replay as before.
  expect_identical(dag_options(select_dag(base, total_select = c(2L, 1L)))$selected$total, c(2L, 1L))
})

test_that("set views of a declared receipt are built as observed and marked afterwards", {
  graph <- unmeasured_query("U -> T; U -> Y; T -> Y")
  dag <- suppressWarnings(construct_dag(graph, latents = "U", quiet = TRUE))
  catalogue <- dag_set_catalogue(dag)
  expect_identical(catalogue$total$unmeasured, list("U"))
  expect_identical(catalogue$direct$unmeasured, list("U"))
  expect_identical(catalogue$source_signature, dag_options(dag)$refinement$source_signature)
  combined <- dag_sets_view(dag, list(total = 1L))
  expect_identical(dagitty::latents(combined), "U")
  expect_identical(dagitty::latents(temporal_dag(combined)), "U")
  expect_identical(dag_checks(combined)$total$status, "preserved")
  expect_identical(dag_checks(combined)$total$unmeasured, "U")
  expect_identical(dag_checks(combined)$total$availability, "needs unmeasured: U")
  expect_true(attr(dag_checks(combined), "current"))
  single <- dag_set_view(dag, "total", 1L)
  expect_identical(dagitty::latents(single), "U")
  expect_identical(.gd_state(single)$set_view$raw_sets, list("U"))
  expect_identical(dag_checks(single)$total$unmeasured, "U")
  # Receipts without declarations keep unmarked views.
  plain <- dag_set_view(construct_dag(graph, quiet = TRUE), "total", 1L)
  expect_length(dagitty::latents(plain), 0L)
  expect_identical(dag_checks(plain)$total$availability, "available now")
})

test_that("restore_nodes keeps the observed fold and the marks of a declared receipt", {
  graph <- unmeasured_query("A_t0 -> T; A_t0 -> Y; A_t0 -> A_t1; A_t1 -> T; A_t1 -> Y; T -> Y; B -> Y")
  dag <- suppressWarnings(construct_dag(graph, latents = "A_t0", quiet = TRUE))
  expect_false("B" %in% names(dag))
  expect_true(.gd_state(dag)$summary$applied)
  restored <- NULL
  expect_warning(restored <- restore_nodes(dag, "B"), "stands for both observed and unmeasured occasions")
  expect_setequal(names(restored), c("T", "A", "Y", "B"))
  expect_true(.gd_state(restored)$summary$applied)
  expect_identical(dagitty::latents(restored), "A")
  expect_identical(dagitty::latents(temporal_dag(restored)), "A_t0")
  expect_identical(dagitty::latents(.gd_state(restored)$analysis), "A_t0")
  expect_identical(as.character(undo_dag(restored)), as.character(dag))
})

test_that("the explicit witness and complete routes keep applying declarations before their search", {
  graph <- unmeasured_query("U -> T; U -> Y; T -> Y")
  for (mode in c("witness", "complete")) {
    dag <- NULL
    expect_warning(dag <- construct_dag(graph, latents = "U", quiet = TRUE, search = mode), "diagnostic graph")
    expect_identical(dag_options(dag)$status[["total"]], "no_sets", info = mode)
    expect_null(dag_options(dag)$unmeasured, info = mode)
    expect_false("unmeasured" %in% names(dag_options(dag)$total), info = mode)
    expect_identical(dag_history(dag)$operation[[1L]], "set_latents", info = mode)
    expect_false("mark_unmeasured" %in% dag_history(dag)$operation, info = mode)
    expect_null(dag_checks(dag)$total$availability, info = mode)
  }
})

test_that("a family blocked by a bidirected edge still has no set when there are no declarations", {
  graph <- unmeasured_query("T <-> Y; T -> Y")
  dag <- NULL
  expect_warning(dag <- construct_dag(graph, quiet = TRUE), "diagnostic graph")
  expect_identical(dag_options(dag)$status, c(direct = "no_sets", total = "no_sets"))
  expect_true(dag_options(dag)$search$total$complete)
  expect_identical(dag_set_catalogue(dag)$source_sets$total, list())
  expect_identical(dag_checks(dag)$total$status, "no eligible reference set")
  expect_identical(vapply(dag_checks(dag), `[[`, character(1), "availability"),
                   c(direct = "no selected set", total = "no selected set"))
  expect_length(dag_options(dag)$unmeasured$nodes, 0L)
  # A valid empty set is different: it is selected and available now.
  empty <- construct_dag(unmeasured_query("T -> Y"), quiet = TRUE)
  expect_identical(dag_checks(empty)$total$required, list(character()))
  expect_identical(dag_checks(empty)$total$availability, "available now")
})

test_that("a budget stop without an option has no selected set, never an available one", {
  run <- unmeasured_warnings(construct_dag(unmeasured_query("A -> T; A -> Y; T -> Y"), quiet = TRUE, steps = 1L))
  dag <- run$value
  expect_identical(dag_options(dag)$status, c(direct = "no_sets", total = "no_sets"))
  expect_false(dag_options(dag)$search$total$complete)
  expect_identical(vapply(dag_checks(dag), `[[`, character(1), "availability"),
                   c(direct = "no selected set", total = "no selected set"))
  expect_true(any(grepl("not evidence that no set exists", run$warnings, fixed = TRUE)))
  expect_false(any(grepl("No eligible", run$warnings, fixed = TRUE)))
})
