cd_test_graph <- function(edges, latent = character(), columns = character(), depth = Inf,
                          policy = NULL) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  dagitty::latents(graph) <- latent
  ids <- .strip_slice(names(graph))
  map <- data.frame(node = names(graph), quantity = ids, member_id = ids,
    occurrence = names(graph), stringsAsFactors = FALSE)
  if (is.null(policy)) policy <- data.frame(node_id = unique(ids), recovered_node = "yes")
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "barley", params = list(water_model = "SoilWater", query_depth = depth)),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = map, policy = policy)
  attr(graph, "crop") <- "barley"
  attr(graph, "water_model") <- "SoilWater"
  class(graph) <- c("dag_by_query", class(graph))
  graph
}

test_that("constructor retains direct and total targets separately and preserves source", {
  graph <- cd_test_graph("U -> T; U -> M; T -> M -> Y; T -> Y", latent = "U")
  before <- serialize(graph, NULL)
  dag <- construct_dag(graph, quiet = TRUE)
  expect_s3_class(dag, "cf_dag")
  expect_identical(serialize(original_query(dag), NULL), before)
  expect_setequal(names(dag), c("T", "Y", "M", "U"))
  expect_length(dagitty::latents(dag), 0L)
  checks <- dag_checks(dag)
  expect_identical(checks$direct$status, "preserved")
  expect_identical(checks$total$status, "preserved")
  expect_identical(checks$direct$required, list("M"))
  expect_identical(checks$total$required, list("U"))
  expect_false(dagitty::isAdjustmentSet(graph, c("M", "U")))
  expect_identical(serialize(graph, NULL), before)
})

test_that("valid empty total family differs from no eligible total sets", {
  randomised <- cd_test_graph("T -> M -> Y; T -> Y")
  dag <- construct_dag(randomised, quiet = TRUE)
  expect_identical(dag_checks(dag)$total$required, list(character()))
  expect_identical(.gd_state(dag)$construction$status[["total"]], "selected")
  confounded <- cd_test_graph("U -> T; U -> Y; T -> Y")
  # A declared unmeasured confounder no longer removes its option: {U} is
  # selected and marked on the DAG. No eligible set needs a structural blocker.
  marked <- suppressWarnings(construct_dag(confounded, latents = "U", quiet = TRUE))
  expect_identical(.gd_state(marked)$construction$status[["total"]], "selected")
  expect_identical(dag_checks(marked)$total$required, list("U"))
  expect_identical(dagitty::latents(marked), "U")
  expect_identical(dag_checks(marked)$total$availability, "needs unmeasured: U")
  blocked <- cd_test_graph("T <-> Y; T -> Y")
  expect_warning(diagnostic <- construct_dag(blocked, quiet = TRUE), "diagnostic graph")
  expect_identical(.gd_state(diagnostic)$construction$status[["total"]], "no_sets")
  expect_length(dag_checks(diagnostic)$total$required, 0L)
  expect_identical(dag_checks(diagnostic)$total$status, "no eligible reference set")
})

test_that("NULL defers without constructing and total waits for direct context", {
  graph <- cd_test_graph("T -> A -> B -> Y; T -> Y")
  options <- construct_dag(graph, direct_select = NULL, quiet = TRUE)
  expect_s3_class(options, "cf_dag_options")
  expect_identical(options$status, c(direct = "deferred", total = "deferred"))
  expect_match(options$choices$total$reason, "Pending direct")
  expect_identical(options$graph, graph)
  expect_equal(options$direct$group, c(1L, 1L))
  expect_equal(options$direct$option, c(1L, 2L))
  dag <- construct_dag(graph, direct_select = c(1, 2), quiet = TRUE)
  expect_identical(.gd_state(dag)$selected$direct, c(1L, 2L))
  expect_setequal(names(dag), c("T", "B", "Y"))
  expect_error(construct_dag(graph, direct_select = c(2, 1), quiet = TRUE), "absent")
})

test_that("total selection uses smallest size group before overlap and stable ties", {
  table <- data.frame(group = c(1L, 1L, 2L), option = c(1L, 2L, 1L), n_resolved = c(1L, 1L, 2L))
  table$resolved <- list("A", "B", c("B", "C"))
  expect_identical(.cd_select(table, "auto", "total", c("B", "C"))$selection, c(1L, 2L))
  choice <- .cd_select(table, "auto", "total", character())
  expect_identical(choice$selection, c(1L, 1L))
  expect_identical(choice$ties, c(1L, 2L))
  expect_identical(.cd_select(table, c(2, 1), "total")$selection, c(2L, 1L))
})

test_that("recorded requests supplement construction and influence total overlap", {
  graph <- cd_test_graph("A -> B -> T; A -> M; T -> M -> Y; T -> Y; V -> Y", columns = "B")
  dag <- construct_dag(graph, quiet = TRUE, search = "complete")
  expect_identical(dag_checks(dag)$total$required, list("B"))
  expect_true("B" %in% names(dag))
  expect_false("V" %in% names(dag))
  manual <- construct_dag(graph, keep_nodes = "V", quiet = TRUE, search = "complete")
  expect_true("V" %in% names(manual))
  expect_error(construct_dag(graph, keep_nodes = "absent", quiet = TRUE, search = "complete"), "absent from")
})

test_that("policy eligibility is separate from declared observation status", {
  graph <- cd_test_graph("A -> B -> T; A -> M; T -> M -> Y; T -> Y")
  attr(graph, "query_record")$policy$recovered_node[attr(graph, "query_record")$policy$node_id == "A"] <- "no: interpretive"
  expect_warning(dag <- construct_dag(graph, quiet = TRUE), "Recovery policy exclusions")
  expect_identical(dag_checks(dag)$total$required, list(character()))
  expect_length(dagitty::latents(.gd_state(dag)$analysis), 0L)
  expect_false("A" %in% names(.gd_state(dag)$analysis))
  expect_error(construct_dag(graph, keep_nodes = "A", quiet = TRUE), "saved recovery policy")
  attr(graph, "query_record")$policy$recovered_node[attr(graph, "query_record")$policy$node_id == "A"] <- "yes:user"
  allowed <- construct_dag(graph, keep_nodes = "A", quiet = TRUE)
  expect_true("A" %in% names(allowed))
})

test_that("temporal names and latent quantities retain exact status", {
  graph <- cd_test_graph("T -> A_t1 -> A_t2 -> Y; T -> Y")
  exact <- suppressWarnings(construct_dag(graph, latents = "A_t1", keep_nodes = "A", quiet = TRUE))
  expect_setequal(dagitty::latents(temporal_dag(exact)), "A_t1")
  # Construction folds as for a fully observed graph. The display node A stands
  # for the declared A_t1 and the observed A_t2, so it is marked and recorded.
  expect_setequal(names(exact), c("T", "A", "Y"))
  expect_identical(dagitty::latents(exact), "A")
  expect_identical(dag_options(exact)$unmeasured$partial,
                   list(A = list(unmeasured = "A_t1", observed = "A_t2")))
  # The explicit complete route still applies the declaration before its fold.
  complete <- suppressWarnings(construct_dag(graph, latents = "A_t1", keep_nodes = "A", quiet = TRUE,
                                             search = "complete"))
  expect_true(all(c("A_t1", "A_t2") %in% names(complete)))
  expect_match(.gd_state(complete)$summary$reason, "Mixed observation status")
  both <- .cd_prepare(graph, "A", NULL)
  expect_setequal(dagitty::latents(both$reference), c("A_t1", "A_t2"))
  expect_error(construct_dag(graph, latents = "T", quiet = TRUE), "cannot be declared latent")
})

test_that("counts expand saved groups without counting temporal occurrences", {
  graph <- cd_test_graph("T -> G_t1 -> G_t2 -> Y; T -> Y")
  record <- attr(graph, "query_record")
  map <- record$node_map
  added <- map[map$quantity == "G", ]
  added$member_id <- "H"
  record$node_map <- rbind(map, added)
  record$policy <- rbind(record$policy, data.frame(node_id = "H", recovered_node = "yes"))
  attr(graph, "query_record") <- record
  options <- construct_dag(graph, direct_select = NULL, total_select = NULL, quiet = TRUE,
    search = "complete")
  expect_equal(options$direct$n_resolved, 2L)
  expect_equal(options$direct$n_concepts, 1L)
  expect_length(options$direct$raw_sets[[1L]], 2L)
  expect_setequal(options$direct$resolved[[1L]], c("G", "H"))
  # The ranked default counts the same members and realises each option once.
  ranked <- construct_dag(graph, direct_select = NULL, total_select = NULL, quiet = TRUE)
  expect_identical(ranked$search_method, "quantity-ranked-v1")
  expect_equal(ranked$direct$n_resolved, 2L)
  expect_equal(ranked$direct$n_concepts, 1L)
  expect_setequal(ranked$direct$resolved[[1L]], c("G", "H"))
  expect_length(ranked$direct$raw_sets[[1L]], 1L)
  expect_true(all(lengths(ranked$direct$raw_sets) == 1L))
  expect_true(all(.gd_validate_sets(.cd_prepare(graph, NULL, NULL)$selection,
    unlist(ranked$direct$raw_sets, recursive = FALSE), "direct")))
})

test_that("cap probe records known omitted temporal alternatives of the same size", {
  graph <- cd_test_graph("T -> A -> B -> Y; T -> Y")
  bounded <- .cd_search(.cd_prepare(graph, NULL, NULL), "direct", 1L)
  expect_false(bounded$complete)
  expect_true(bounded$omitted_equal_smallest)
  expect_length(bounded$sets, 1L)
  expect_length(bounded$probe, 1L)
  expect_identical(bounded$probe_quantity_count, bounded$table$n_resolved[[1L]])
})

test_that("dry run previews complete construction and keeps narrower supplied reference", {
  graph <- cd_test_graph("T -> M -> Y; T -> Y", depth = 1L)
  attr(graph, "full_vertices") <- c(names(graph), "UnknownExtra")
  expect_warning(plan <- construct_dag(graph, dry_run = TRUE, quiet = TRUE), "supplied query_depth=1")
  expect_s3_class(plan, "cf_dag_plan")
  expect_s3_class(plan$graph, "cf_dag")
  expect_false("UnknownExtra" %in% names(.gd_state(plan$graph)$analysis))
  expect_identical(original_query(plan$graph), graph)
  actual <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  expect_identical(.gd_signature(plan$graph), .gd_signature(actual))
  expect_identical(plan$checks, .gd_state(actual)$checks)
})

test_that("constructor call replays without querying or changing the source", {
  graph <- cd_test_graph("T -> A -> B -> Y; T -> Y")
  dag <- construct_dag(graph, direct_select = c(1, 2), quiet = TRUE)
  call <- dag_call(dag)
  expect_identical(call[[1L]], as.name("construct_dag"))
  again <- eval(call)
  expect_identical(.gd_signature(again), .gd_signature(dag))
  expect_identical(original_query(again), graph)
  expect_gte(nrow(dag_history(dag)), 3L)
  alternative <- select_dag(dag, c(1, 1))
  expect_setequal(names(alternative), c("T", "A", "Y"))
  expect_identical(tail(dag_history(alternative)$operation, 1L), "select_dag")
  expect_identical(.gd_state(alternative)$input_history[[1L]], .gd_state(dag)$input_record)
  expect_identical(.gd_signature(undo_dag(alternative)), .gd_signature(dag))
  options <- construct_dag(graph, direct_select = NULL, quiet = TRUE)
  expect_s3_class(eval(dag_call(options)), "cf_dag_options")
  preview <- construct_dag(graph, dry_run = TRUE, quiet = TRUE)
  expect_s3_class(eval(dag_call(preview)), "cf_dag_plan")
})

test_that("new invalid alternatives do not pass simply because the chosen set survived", {
  graph <- cd_test_graph("T -> A -> B -> Y; T -> Y; C")
  prepared <- .cd_prepare(graph, NULL, NULL)
  search <- list(direct = .cd_search(prepared, "direct", 20L), total = .cd_search(prepared, "total", 20L))
  choices <- list(direct = .cd_select(search$direct$table, c(1, 1), "direct"),
    total = .cd_select(search$total$table, "auto", "total"))
  final <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> A -> C -> Y; T -> Y }")
  checks <- .cd_checks(prepared, final, search, choices, 20L)
  expect_identical(checks$direct$status, "changed")
  expect_true(any(vapply(checks$direct$introduced, function(x) x$status == "not valid on reference", logical(1))))
})

test_that("detailed receipts name query and genuinely redundant requests", {
  graph <- cd_test_graph("U -> T; U -> M; T -> M -> Y; T -> Y", columns = "U")
  expect_output(expect_message(dag <- construct_dag(graph, quiet = FALSE), "Requested U already retained: selected total membership"), "Causal query:.*treatment")
  expect_identical(.gd_state(dag)$request_notice$reason, "selected total membership")
  expect_silent(construct_dag(graph, quiet = TRUE))
})

test_that("recorded human labels and bare lone occasions retain source identities", {
  graph <- cd_test_graph("T -> A -> Y; T -> Y")
  record <- attr(graph, "query_record")
  record$node_map$quantity[record$node_map$node == "A"] <- "Canonical.A"
  record$node_map$member_id[record$node_map$node == "A"] <- "Canonical.A"
  record$node_map$occurrence[record$node_map$node == "A"] <- "Canonical.A_t2"
  record$policy$node_id[record$policy$node_id == "A"] <- "Canonical.A"
  attr(graph, "query_record") <- record
  dag <- construct_dag(graph, quiet = TRUE)
  expect_identical(.gd_state(dag)$construction$direct$concepts, list("Canonical.A"))
  expect_identical(.cd_match_nodes(graph, "Canonical.A_t2"), "A")
  from <- .cd_prepare(graph, NULL, NULL)$reference
  to <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> A_t3 -> Y; T -> Y }")
  attr(to, "cf_node_map") <- data.frame(node = "A_t3", quantity = "Canonical.A", member_id = "Canonical.A", occurrence = "Canonical.A_t3")
  expect_null(.gd_map_set("A", from, to))
})
