set_view_graph <- function(text, columns = character()) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", text, "}"))
  ids <- .strip_slice(names(graph))
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "barley", params = list(water_model = "SoilWater", query_depth = Inf)),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = names(graph), quantity = ids, member_id = ids,
                          occurrence = names(graph), stringsAsFactors = FALSE),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  graph
}

test_that("canonical set families merge equal members and prune only proper supersets", {
  graph <- set_view_graph("A_t0; A_t1; A_t2; B_t0; C_t0; D_t0; G_t0; H_t0; T -> Y")
  prepared <- .cd_prepare(graph, NULL, NULL)
  prepared$map$member_id[prepared$map$quantity %in% c("G", "H")] <- "C"
  raw <- list("A_t0", c("A_t1", "B_t0"), "A_t2", c("C_t0", "D_t0"),
              c("G_t0", "D_t0"), c("H_t0", "D_t0"))
  table <- .cd_candidates(raw, prepared)
  before <- serialize(table, NULL)
  family <- .gd_set_family(table, "direct")
  expect_identical(family$table$variables, list("A", c("C", "D")))
  expect_identical(family$table$n_variables, c(1L, 2L))
  expect_identical(family$table$id, c("direct-set-1", "direct-set-2"))
  expect_identical(family$table$index, 1:2)
  expect_length(family$table$raw_sets[[1L]], 2L)
  expect_length(family$table$raw_sets[[2L]], 3L)
  expect_equal(nrow(family$table$source_options[[2L]]), 3L)
  expect_identical(family$supersets$variables, list(c("A", "B")))
  expect_identical(serialize(table, NULL), before)
  expect_identical(.gd_set_family(table, "direct"), family)
})

test_that("catalogue keeps valid empty, absent and incomplete families distinct", {
  graph <- set_view_graph("T -> A -> B -> Y; T -> Y")
  dag <- construct_dag(graph, quiet = TRUE, search = "complete")
  original_options <- serialize(dag_options(dag), NULL)
  cat <- dag_set_catalogue(dag)
  expect_identical(cat$direct$variables, list("A", "B"))
  expect_identical(cat$total$variables, list(character()))
  expect_identical(cat$total$raw_sets, list(list(character())))
  expect_true(all(cat$complete))
  expect_false(any(cat$incomplete))
  expect_identical(serialize(dag_options(dag), NULL), original_options)
  # Reproduce a saved bounded receipt through the private backend, rather than
  # restoring a public cap. Existing catalogues remain inspectable as incomplete.
  bounded <- local({
    search <- .gd_adjustment_sets
    local_mocked_bindings(.gd_adjustment_sets = function(dag, effect, max_sets, ...)
      search(dag, effect, 1L), .package = "causaliflower")
    saved <- suppressWarnings(construct_dag(graph, quiet = TRUE, search = "complete"))
    state <- .gd_state(saved)
    state$construction$search_method <- "source-verified-refinement-v1"
    state$construction$arguments$max_sets <- 1L
    .gd_attach(saved, state)
  })
  expect_false(dag_set_catalogue(bounded)$complete[["direct"]])
  expect_true(dag_set_catalogue(bounded)$incomplete[["direct"]])
  bounded_view <- suppressWarnings(dag_set_view(bounded, "direct", 1L))
  expect_false(.gd_state(bounded_view)$set_view$source_complete)
  expect_false(dag_checks(bounded_view)$direct$reference_complete)
  expect_true(dag_checks(bounded_view)$direct$final_complete)
  expect_error(select_dag(bounded, direct_select = c(1L, 1L)), "predates the current search")
  confounded <- suppressWarnings(construct_dag(set_view_graph("U -> T; U -> Y; T -> Y"),
                                               latents = "U", quiet = TRUE, search = "complete"))
  absent <- dag_set_catalogue(confounded)
  expect_equal(nrow(absent$total), 0L)
  expect_true(absent$complete[["total"]])
  expect_length(absent$source_sets$total, 0L)
  expect_error(dag_set_view(confounded, "total", 1L), "available row")
})

test_that("set views reuse cached source, preserve time and scope checks to one family", {
  source <- set_view_graph("T -> A_t1 -> A_t2 -> Y; T -> Y; V", columns = "V")
  source_before <- serialize(source, NULL)
  combined <- construct_dag(source, quiet = TRUE, search = "complete")
  combined_before <- serialize(combined, NULL)
  search <- .gd_adjustment_sets
  calls <- list()
  forbidden <- function(...) stop("Source query or construction was called")
  local_mocked_bindings(query_graph = forbidden, get_dag = forbidden,
    construct_dag = forbidden,
    .gd_adjustment_sets = function(dag, effect, max_sets) {
      calls[[length(calls) + 1L]] <<- list(graph = .gd_signature(dag), effect = effect)
      search(dag, effect, max_sets)
    })
  catalogue <- dag_set_catalogue(combined)
  expect_length(calls, 0L)
  expect_identical(catalogue$direct$variables, list("A"))
  direct <- dag_set_view(combined, "direct", 1L, simplify = FALSE)
  expect_length(calls, 1L)
  expect_identical(calls[[1L]]$effect, "direct")
  expect_identical(calls[[1L]]$graph, .gd_signature(direct))
  expect_setequal(names(direct), c("T", "Y", "A_t1", "A_t2", "V"))
  expect_s3_class(direct, "cf_dag")
  expect_identical(names(dag_checks(direct)), "direct")
  expect_identical(dag_checks(direct)$direct$status, "preserved")
  expect_length(.gd_state(direct)$set_view$mapped_sets, 2L)
  expect_identical(.gd_signature(temporal_dag(direct)), .gd_signature(direct))
  expect_false(.gd_state(direct)$summary$applied)
  expect_false(.gd_state(direct)$controls$simplify)
  expect_identical(tail(dag_history(direct)$operation, 1L), "dag_set_view")
  total <- dag_set_view(direct, "total", 1L, simplify = FALSE)
  expect_length(calls, 2L)
  expect_identical(calls[[2L]]$effect, "total")
  expect_setequal(names(total), c("T", "Y", "V"))
  expect_identical(dag_checks(total)$total$required, list(character()))
  expect_identical(dag_checks(total)$total$status, "preserved")
  expect_identical(serialize(original_query(total), NULL), source_before)
  expect_identical(serialize(combined, NULL), combined_before)
  expect_identical(serialize(source, NULL), source_before)
  expect_identical(.gd_signature(undo_dag(total)), .gd_signature(direct))
})

test_that("one effect view does not retain the other selected family by accident", {
  source <- set_view_graph("U -> T; U -> M; T -> M -> Y; T -> Y")
  combined <- construct_dag(source, quiet = TRUE)
  direct <- dag_set_view(combined, "direct", 1L)
  total <- dag_set_view(combined, "total", 1L)
  expect_setequal(.gd_state(direct)$set_view$variables, "M")
  expect_setequal(.gd_state(total)$set_view$variables, "U")
  expect_false("U" %in% names(direct))
  expect_false("M" %in% names(total))
  expect_true(any(dag_edges(direct)$e == "<->"))
  expect_identical(dag_checks(direct)$direct$status, "preserved")
  expect_identical(dag_checks(total)$total$status, "preserved")
})

test_that("views preserve explicit latent requests and retain parents when no total set exists", {
  source <- set_view_graph("T -> A_t1 -> A_t2 -> Y; T -> Y; U -> A_t1")
  combined <- suppressWarnings(construct_dag(source, latents = "U", keep_nodes = "U", quiet = TRUE))
  view <- dag_set_view(combined, "direct", 1L)
  expect_setequal(dagitty::latents(view), "U")
  expect_true("U" %in% names(view))
  joint <- dagitty::dagitty("dag { X1 [exposure]; X2 [exposure]; Y [outcome]; P -> X1 -> M -> X2 -> Y; M -> Y }")
  joint_dag <- suppressWarnings(construct_dag(joint, quiet = TRUE))
  joint_view <- dag_set_view(joint_dag, "direct", 1L)
  expect_true(.gd_state(joint_view)$set_view$identification$applied)
  expect_true(all(c("P", "M", "X1", "X2", "Y") %in% names(joint_view)))
  expect_match(.gd_state(joint_view)$set_view$identification$status, "preserved")
  expect_true(.gd_state(joint_view)$set_view$summary$requested)
  expect_false(.gd_state(joint_view)$set_view$summary$applied)
  expect_match(.gd_state(joint_view)$set_view$summary$reason, "joint total identification")
})

test_that("set view switches guard graph, source, preparation, methods and indexes", {
  source <- set_view_graph("T -> A -> B -> Y; T -> Y; V")
  combined <- construct_dag(source, quiet = TRUE, search = "complete")
  edited <- combined
  dagitty::latents(edited) <- "A"
  attr(edited, "cf_provenance") <- .gd_state(combined)
  expect_error(dag_set_catalogue(edited), "changed outside")
  restored <- restore_nodes(combined, "V")
  expect_error(dag_set_view(restored, "direct", 1), "Commands changed")
  old <- .gd_state(combined)
  old$construction$search_method <- "old"
  expect_error(dag_set_catalogue(.gd_attach(combined, old)), "unsupported search method")
  changed_source <- .gd_state(combined)
  attr(changed_source$original_query, "external_edit") <- TRUE
  expect_error(dag_set_catalogue(.gd_attach(combined, changed_source)), "source or search receipt")
  changed_reference <- .gd_state(combined)
  dagitty::latents(changed_reference$analysis) <- "B"
  expect_error(dag_set_catalogue(.gd_attach(combined, changed_reference)), "source or search receipt")
  changed_map <- .gd_state(combined)
  attr(changed_map$analysis, "cf_node_map")$member_id[1L] <- "Wrong"
  expect_error(dag_set_view(.gd_attach(combined, changed_map), "direct", 1), "source identities")
  moved <- combined
  dagitty::coordinates(moved) <- list(x = stats::setNames(seq_along(names(moved)), names(moved)),
                                     y = stats::setNames(seq_along(names(moved)), names(moved)))
  attr(moved, "input_labels") <- attr(combined, "input_labels")
  attr(moved, "cf_node_map") <- attr(combined, "cf_node_map")
  attr(moved, "cf_provenance") <- .gd_state(combined)
  expect_identical(dag_set_catalogue(moved), dag_set_catalogue(combined))
  for (index in list(0, -1, 1.5, NA_real_, Inf, "1", c(1, 2), 999))
    expect_error(dag_set_view(combined, "direct", index), "available row")
  expect_error(dag_set_view(combined, "other", 1), "arg")
  for (simplify in list(NA, NULL, 1, "TRUE", c(TRUE, FALSE)))
    expect_error(dag_set_view(combined, "direct", 1, simplify = simplify), "simplify")
})

test_that("default summaries preserve raw set witnesses and the full temporal projection", {
  source <- set_view_graph("T -> A_t1 -> A_t2 -> Y")
  combined <- construct_dag(source, quiet = TRUE, search = "complete")
  before <- serialize(combined, NULL)
  catalogue <- dag_set_catalogue(combined)
  exact <- dag_set_view(combined, "direct", 1L, simplify = FALSE)
  summary <- dag_set_view(combined, "direct", 1L)
  state <- .gd_state(summary)
  expect_setequal(names(summary), c("T", "A", "Y"))
  expect_setequal(names(temporal_dag(summary)), c("T", "A_t1", "A_t2", "Y"))
  expect_identical(.gd_signature(temporal_dag(summary)), .gd_signature(exact))
  expect_true(state$set_view$summary$requested)
  expect_true(state$set_view$summary$applied)
  expect_null(state$set_view$summary$reason)
  expect_identical(unname(state$set_view$summary$counts$displayed), c(3L, 2L))
  expect_identical(state$set_view$raw_sets, catalogue$direct$raw_sets[[1L]])
  expect_identical(state$set_view$mapped_sets, list("A", "A"))
  expect_identical(dag_checks(summary)$direct$status, "preserved")
  expect_identical(state$temporal_checks$direct$status, "preserved")
  expect_identical(dag_set_catalogue(summary), catalogue)
  expect_identical(serialize(original_query(summary), NULL), serialize(source, NULL))
  expect_identical(serialize(combined, NULL), before)
  expect_identical(.gd_signature(undo_dag(summary)), .gd_signature(combined))
  expect_true(tail(dag_history(summary, details = TRUE), 1L)[[1L]]$arguments$simplify)
})

test_that("summaries retain distinct endpoint occasions", {
  source <- dagitty::dagitty("dag { X_t2 [exposure]; Y [outcome]; X_t1 -> X_t2 -> A -> Y; X_t1 -> A }")
  combined <- suppressWarnings(construct_dag(source, keep_nodes = "X_t1", quiet = TRUE))
  view <- dag_set_view(combined, "direct", 1L)
  expect_true(all(c("X_t1", "X_t2") %in% names(view)))
  expect_false(.gd_state(view)$summary$applied)
  expect_match(.gd_state(view)$summary$reason, "distinct occasions of an endpoint")
  expect_identical(.gd_signature(view), .gd_signature(temporal_dag(view)))
})

test_that("changed or incompletely checked summaries fall back with the attempted receipt", {
  combined <- construct_dag(set_view_graph("T -> A_t1 -> A_t2 -> Y"), quiet = TRUE, search = "complete")
  search <- .gd_adjustment_sets
  local_mocked_bindings(.gd_adjustment_sets = function(dag, effect, max_sets) {
    result <- search(dag, effect, max_sets)
    if ("A" %in% names(dag)) result$complete <- FALSE
    result
  })
  view <- dag_set_view(combined, "direct", 1L)
  state <- .gd_state(view)
  expect_false(state$summary$applied)
  expect_match(state$summary$reason, "not fully checked")
  expect_true(state$summary$attempted$summary$applied)
  expect_identical(state$summary$attempted$checks$status, "not fully checked")
  expect_identical(dag_checks(view)$direct$status, "preserved")
  expect_setequal(names(view), c("T", "A_t1", "A_t2", "Y"))
})

test_that("time summaries do not invent an endpoint path by merging disconnected occasions", {
  source <- set_view_graph("T -> A_t1; A_t2 -> Y")
  proposed <- .gd_summary(source, keep = c("A_t1", "A_t2"))
  expect_false(proposed$applied)
  expect_match(proposed$reason, "proper directed path from treatment to outcome")
  expect_identical(.gd_signature(proposed$graph), .gd_signature(source))
  combined <- construct_dag(source, keep_nodes = c("A_t1", "A_t2"), quiet = TRUE)
  expect_false(.gd_state(combined)$summary$applied)
  view <- dag_set_view(combined, "total", 1L)
  expect_false(.gd_state(view)$summary$applied)
  expect_identical(.gd_state(view)$set_view$variables, character())
  expect_setequal(names(view), c("T", "A_t1", "A_t2", "Y"))
})

test_that("proper endpoint paths exclude other treatments from the interior", {
  source <- dagitty::dagitty("dag { X1 [exposure]; X2 [exposure]; Y [outcome]; X1 -> A_t1; A_t2 -> Y; X1 -> X2 -> Y }")
  paths <- .gd_endpoint_paths(source)
  expect_false(unname(paths[["X1\rY"]]))
  expect_true(unname(paths[["X2\rY"]]))
  proposed <- .gd_summary(source, keep = c("A_t1", "A_t2"))
  expect_false(proposed$applied)
  expect_match(proposed$reason, "proper directed path from treatment to outcome")
})

test_that("repeated compression cannot remove the only proper endpoint path", {
  source <- set_view_graph("T -> A_t0; B_t0 -> A_t0; A_t0 -> B_t1; B_t1 -> A_t1; B_t1 -> Y")
  proposed <- .gd_summary(source, keep = c("A_t0", "A_t1", "B_t0", "B_t1"))
  expect_true(unname(.gd_endpoint_paths(source)[["T\rY"]]))
  expect_false(proposed$applied)
  expect_match(proposed$reason, "proper directed path from treatment to outcome")
  expect_identical(.gd_signature(proposed$graph), .gd_signature(source))
})
