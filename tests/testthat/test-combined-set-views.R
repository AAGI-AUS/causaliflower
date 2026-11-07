combined_set_graph <- function(text) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", text, "}"))
  ids <- .strip_slice(names(graph))
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "barley", params = list(water_model = "SoilWater", query_depth = Inf)),
    requests = list(columns = list(resolved = character())),
    node_map = data.frame(node = names(graph), quantity = ids, member_id = ids,
                          occurrence = names(graph), stringsAsFactors = FALSE),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes"))
  graph
}

test_that("combined views preserve every alternative with one search per family", {
  source <- combined_set_graph("A -> T; A -> B -> Y; A -> C -> Y; T -> Y")
  combined <- construct_dag(source, quiet = TRUE, search = "complete")
  before <- serialize(combined, NULL)
  catalogue <- dag_set_catalogue(combined)
  expect_identical(catalogue$total$variables, list("A", c("B", "C")))
  expect_identical(catalogue$direct$variables, list("A", c("B", "C")))
  search <- .gd_adjustment_sets
  calls <- list()
  forbidden <- function(...) stop("A source query or construction was called")
  local_mocked_bindings(query_graph = forbidden, get_dag = forbidden,
    construct_dag = forbidden,
    .gd_adjustment_sets = function(graph, effect, max_sets) {
      calls[[length(calls) + 1L]] <<- list(effect = effect, nodes = names(graph))
      search(graph, effect, max_sets)
    })
  view <- dag_sets_view(combined, list(total = 1:2, direct = 1:2), simplify = FALSE)
  state <- .gd_state(view)
  expect_length(calls, 2L)
  expect_identical(vapply(calls, `[[`, character(1), "effect"), c("total", "direct"))
  expect_setequal(names(view), c("A", "B", "C", "T", "Y"))
  expect_identical(state$set_view$mode, "combined")
  expect_identical(state$set_view$selections, list(total = 1:2, direct = 1:2))
  for (effect in c("total", "direct")) {
    expect_identical(state$set_view$alternatives[[effect]], catalogue[[effect]])
    expect_identical(state$checks[[effect]]$required, list("A", c("B", "C")))
    expect_identical(state$checks[[effect]]$status, "preserved")
    expect_true(all(vapply(state$checks[[effect]]$records, function(x) x$status == "preserved", logical(1))))
  }
  expect_identical(serialize(combined, NULL), before)
  expect_identical(serialize(original_query(view), NULL), serialize(source, NULL))
  expect_identical(dag_options(view), dag_options(combined))
  expect_identical(dag_set_catalogue(view), catalogue)
  expect_identical(.gd_signature(undo_dag(view)), .gd_signature(combined))
  expect_identical(tail(dag_history(view)$operation, 1L), "dag_sets_view")
})

test_that("the union of context from both families is never validated as a total adjustment set", {
  combined <- construct_dag(combined_set_graph("T -> M -> Y; T -> Y"), quiet = TRUE)
  view <- dag_sets_view(combined, list(total = 1L, direct = 1L))
  state <- .gd_state(view)
  expect_identical(state$set_view$variables, "M")
  expect_identical(state$checks$total$required, list(character()))
  expect_identical(state$checks$direct$required, list("M"))
  expect_false(.gd_validate_sets(view, list("M"), "total")[[1L]])
  expect_true(.gd_validate_sets(view, list("M"), "direct")[[1L]])
  expect_true(all(vapply(state$checks, function(x) x$status == "preserved", logical(1))))
  direct <- dag_sets_view(combined, list(direct = 1L))
  expect_identical(names(dag_checks(direct)), "direct")
  expect_identical(names(.gd_state(direct)$set_view$source_status), "direct")
})

test_that("no returned sets differ from valid empty sets and disabled families", {
  # Declared unmeasured variables no longer remove options on the ranked route,
  # so the family without a set comes from a structural blocker without declarations.
  source <- combined_set_graph("T <-> Y; T -> Y")
  combined <- suppressWarnings(construct_dag(source, quiet = TRUE))
  catalogue <- dag_set_catalogue(combined)
  expect_equal(nrow(catalogue$total), 0L)
  expect_equal(nrow(catalogue$direct), 0L)
  view <- dag_sets_view(combined, list(total = integer(), direct = integer()))
  state <- .gd_state(view)
  expect_identical(state$set_view$variables, character())
  expect_identical(state$set_view$source_status, c(total = "no_sets", direct = "no_sets"))
  expect_identical(state$set_view$raw_sets, list(total = list(), direct = list()))
  expect_true(all(vapply(state$checks, function(x) x$status == "no eligible reference set", logical(1))))
  expect_true(any(dag_edges(view)$e == "<->"))
  expect_identical(.gd_signature(temporal_dag(view)), .gd_signature(view))
  expect_error(dag_sets_view(combined, list(total = 1L)), "available rows")
  ordinary <- construct_dag(combined_set_graph("T -> Y"), quiet = TRUE)
  expect_error(dag_sets_view(ordinary, list(total = integer())), "available rows")
  empty <- dag_sets_view(ordinary, list(total = 1L))
  expect_identical(.gd_state(empty)$set_view$raw_sets$total, list(character()))
  expect_identical(dag_checks(empty)$total$status, "preserved")
})

test_that("enabled families without options retain the exact parents needed for identification", {
  source <- dagitty::dagitty("dag { X1 [exposure]; X2 [exposure]; Y [outcome]; P -> X1 -> M -> X2 -> Y; M -> Y }")
  combined <- suppressWarnings(construct_dag(source, quiet = TRUE))
  expect_equal(nrow(dag_set_catalogue(combined)$total), 0L)
  view <- dag_sets_view(combined, list(total = integer()))
  state <- .gd_state(view)
  expect_true(state$set_view$identification$applied)
  expect_true(all(c("P", "M", "X1", "X2", "Y") %in% names(view)))
  expect_false(state$summary$applied)
  expect_match(state$summary$reason, "joint total identification")
  expect_identical(dag_checks(view)$total$status, "no eligible reference set")
})

test_that("combined summary acceptance requires every enabled family to pass", {
  combined <- construct_dag(combined_set_graph("T -> A_t1 -> A_t2 -> Y"), quiet = TRUE, search = "complete")
  search <- .gd_adjustment_sets
  local_mocked_bindings(.gd_adjustment_sets = function(graph, effect, max_sets) {
    result <- search(graph, effect, max_sets)
    if (effect == "direct" && "A" %in% names(graph)) result$complete <- FALSE
    result
  })
  view <- dag_sets_view(combined, list(total = 1L, direct = 1L))
  state <- .gd_state(view)
  expect_false(state$summary$applied)
  expect_match(state$summary$reason, "direct: not fully checked")
  expect_identical(state$summary$attempted$checks$total$status, "preserved")
  expect_identical(state$summary$attempted$checks$direct$status, "not fully checked")
  expect_true(all(vapply(state$checks, function(x) x$status == "preserved", logical(1))))
  expect_setequal(names(view), c("T", "A_t1", "A_t2", "Y"))
})

test_that("combined view selection inputs are explicit and validated", {
  combined <- construct_dag(combined_set_graph("T -> A -> B -> Y; T -> Y"), quiet = TRUE, search = "complete")
  bad <- list(list(), list(1L), list(wrong = 1L), list(total = 1L, total = 1L),
    list(total = NULL), list(total = NA_real_), list(total = c(1, 1)),
    list(total = 1.5), list(total = -1), list(total = Inf), list(total = 99),
    list(total = "1"), list(total = 1 + 1i), data.frame(total = 1L))
  for (selections in bad) expect_error(dag_sets_view(combined, selections), "selections")
  for (simplify in list(NA, NULL, 1, "TRUE", c(TRUE, FALSE)))
    expect_error(dag_sets_view(combined, list(total = 1L), simplify), "simplify")
  expect_identical(.gd_state(dag_sets_view(combined, list(direct = c(2L, 1L))))$set_view$selections,
                   list(direct = 1:2))
})

test_that("a family without options cannot hide an invalid introduced set", {
  combined <- suppressWarnings(construct_dag(
    combined_set_graph("U -> T; U -> Y; T -> Y"), latents = "U", quiet = TRUE, search = "complete"))
  search <- .gd_adjustment_sets
  local_mocked_bindings(.gd_adjustment_sets = function(graph, effect, max_sets) {
    result <- search(graph, effect, max_sets)
    result$sets <- list(character())
    result$complete <- TRUE
    result
  })
  expect_warning(view <- dag_sets_view(combined, list(total = integer())), "changed or were not fully verified")
  expect_identical(dag_checks(view)$total$status, "changed")
  expect_identical(dag_checks(view)$total$introduced[[1L]]$status, "not valid on reference")
  expect_false(.gd_state(view)$summary$applied)
  # The invalid introduced set is refused either by the display's check against the source
  # (merged nodes read as all represented occurrences) or by the later family check.
  expect_match(.gd_state(view)$summary$reason, "total: changed|not valid on the source")
})
