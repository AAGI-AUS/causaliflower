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
  graph
}

test_that("constructor retains direct and total targets separately and preserves source", {
  graph <- cd_test_graph("U -> T; U -> M; T -> M -> Y; T -> Y", latent = "U")
  before <- serialize(graph, NULL)
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_s3_class(dag, "cf_dag")
  expect_identical(serialize(dagby_source(dag), NULL), before)
  expect_setequal(names(dag), c("T", "Y", "M", "U"))
  expect_length(dagitty::latents(dag), 0L)
  checks <- dagby_checks(dag)
  expect_identical(checks$direct$status, "preserved")
  expect_identical(checks$total$status, "preserved")
  kept <- dagby_options(dag)$quantities
  expect_identical(kept$direct$members, "M")
  expect_identical(kept$total$members, "U")
  expect_false(dagitty::isAdjustmentSet(graph, c("M", "U")))
  expect_identical(serialize(graph, NULL), before)
})

test_that("valid empty total family differs from no eligible total sets", {
  randomised <- cd_test_graph("T -> M -> Y; T -> Y")
  dag <- dagby_construct(randomised, quiet = TRUE, size = "largest")
  expect_identical(dagby_options(dag)$quantities$total$state, "empty_set")
  expect_identical(dagby_checks(dag)$total$state, "empty_set")
  confounded <- cd_test_graph("U -> T; U -> Y; T -> Y")
  marked <- suppressWarnings(dagby_construct(confounded, latents = "U", quiet = TRUE, size = "largest"))
  expect_identical(dagby_options(marked)$quantities$total$state, "sets")
  expect_identical(dagby_options(marked)$quantities$total$members, "U")
  expect_identical(dagitty::latents(marked), "U")
  expect_identical(dagby_checks(marked)$total$availability, "needs an unmeasured variable")
  blocked <- cd_test_graph("T <-> Y; T -> Y")
  expect_warning(diagnostic <- dagby_construct(blocked, quiet = TRUE, size = "largest"), "No eligible")
  expect_identical(dagby_options(diagnostic)$quantities$total$state, "no_sets")
  expect_identical(dagby_checks(diagnostic)$total$status, "no eligible reference set")
})


test_that("recorded requests supplement construction", {
  graph <- cd_test_graph("A -> B -> T; A -> M; T -> M -> Y; T -> Y; V -> Y", columns = "B")
  dag <- dagby_construct(graph, 0, quiet = TRUE)
  expect_true("B" %in% names(dag))
  expect_false("V" %in% names(dag))
  manual <- dagby_construct(graph, 0, keep_nodes = "V", quiet = TRUE)
  expect_true("V" %in% names(manual))
  expect_error(dagby_construct(graph, keep_nodes = "absent", quiet = TRUE, size = "largest"), "absent from")
  expect_true("B" %in% names(dagby_construct(graph, 0, quiet = TRUE)))
})

test_that("policy eligibility is separate from declared observation status", {
  graph <- cd_test_graph("A -> B -> T; A -> M; T -> M -> Y; T -> Y")
  attr(graph, "query_record")$policy$recovered_node[attr(graph, "query_record")$policy$node_id == "A"] <- "no: interpretive"
  expect_warning(dag <- dagby_construct(graph, quiet = TRUE, size = "largest"), "Recovery policy exclusions")
  expect_identical(dagby_options(dag)$quantities$total$state, "empty_set")
  prepared <- suppressWarnings(.cd_prepare(graph, NULL, NULL))
  expect_length(dagitty::latents(prepared$analysis), 0L)
  expect_false("A" %in% names(prepared$analysis))
  expect_length(dagitty::latents(dagby_temporal(dag)), 0L)
  expect_false("A" %in% names(dagby_temporal(dag)))
  expect_error(dagby_construct(graph, keep_nodes = "A", quiet = TRUE, size = "largest"), "saved recovery policy")
  attr(graph, "query_record")$policy$recovered_node[attr(graph, "query_record")$policy$node_id == "A"] <- "yes:user"
  allowed <- dagby_construct(graph, keep_nodes = "A", quiet = TRUE, size = "largest")
  expect_true("A" %in% names(allowed))
})

test_that("temporal names and latent quantities retain exact status", {
  graph <- cd_test_graph("T -> A_t1 -> A_t2 -> Y; T -> Y")
  exact <- suppressWarnings(dagby_construct(graph, latents = "A_t1", keep_nodes = "A", quiet = TRUE, size = "largest"))
  expect_setequal(dagitty::latents(dagby_temporal(exact)), "A_t1")
  expect_setequal(names(exact), c("T", "A", "Y"))
  expect_identical(dagitty::latents(exact), "A")
  expect_identical(dagby_options(exact)$unmeasured$partial,
                   list(A = list(unmeasured = "A_t1", observed = "A_t2")))
  both <- .cd_prepare(graph, "A", NULL)
  expect_setequal(dagitty::latents(both$reference), c("A_t1", "A_t2"))
  expect_error(dagby_construct(graph, latents = "T", quiet = TRUE, size = "largest"), "cannot be declared latent")
})

test_that("the ranked search counts saved groups without counting temporal occurrences", {
  graph <- cd_test_graph("T -> G_t1 -> G_t2 -> Y; T -> Y")
  record <- attr(graph, "query_record")
  map <- record$node_map
  added <- map[map$quantity == "G", ]
  added$member_id <- "H"
  record$node_map <- rbind(map, added)
  record$policy <- rbind(record$policy, data.frame(node_id = "H", recovered_node = "yes"))
  attr(graph, "query_record") <- record
  prepared <- .cd_prepare(graph, NULL, NULL)
  ranked <- .cd_ranked_search(prepared, "direct")
  expect_identical(ranked$ranked$method, "quantity-ranked-v1")
  expect_equal(ranked$table$n_resolved, 2L)
  expect_equal(ranked$table$n_concepts, 1L)
  expect_setequal(ranked$table$resolved[[1L]], c("G", "H"))
  expect_length(ranked$table$raw_sets[[1L]], 1L)
  expect_true(all(lengths(ranked$table$raw_sets) == 1L))
  expect_true(all(.gd_validate_sets(prepared$selection, unlist(ranked$table$raw_sets, recursive = FALSE), "direct")))
  expect_true("G" %in% names(dagby_construct(graph, 0, quiet = TRUE)))
  expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "stands for several source members")
})


test_that("dry run previews complete construction and keeps narrower supplied reference", {
  graph <- cd_test_graph("T -> M -> Y; T -> Y", depth = 1L)
  attr(graph, "full_vertices") <- c(names(graph), "UnknownExtra")
  expect_warning(plan <- dagby_construct(graph, dry_run = TRUE, quiet = TRUE, size = "largest"), "supplied query_depth=1")
  expect_s3_class(plan, "cf_dag_plan")
  expect_s3_class(plan$graph, "cf_dag")
  expect_false("UnknownExtra" %in% names(suppressWarnings(.cd_prepare(graph, NULL, NULL))$analysis))
  expect_false("UnknownExtra" %in% names(dagby_temporal(plan$graph)))
  expect_identical(dagby_source(plan$graph), graph)
  actual <- suppressWarnings(dagby_construct(graph, quiet = TRUE, size = "largest"))
  expect_identical(.gd_signature(plan$graph), .gd_signature(actual))
  expect_identical(plan$checks, .gd_state(actual)$checks)
})

test_that("constructor call replays without querying or changing the source", {
  graph <- cd_test_graph("T -> A -> B -> Y; T -> Y")
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  call <- dagby_call(dag)
  expect_identical(call[[1L]], quote(causaliflower:::.cd_replay))
  expect_identical(call$graph, graph)
  expect_identical(call$candidate, dagby_options(dag)$candidate)
  expect_identical(call$arguments, dagby_options(dag)$arguments)
  again <- eval(call)
  expect_identical(.gd_signature(again), .gd_signature(dag))
  expect_identical(dagby_source(again), graph)
  preview <- dagby_construct(graph, dry_run = TRUE, quiet = TRUE, size = "largest")
  expect_s3_class(eval(dagby_call(preview)), "cf_dag_plan")
})


test_that("detailed receipts name the query", {
  graph <- cd_test_graph("U -> T; U -> M; T -> M -> Y; T -> Y", columns = "U")
  expect_output(dag <- dagby_construct(graph, quiet = FALSE, size = "largest"), "Causal query:.*treatment")
  expect_silent(dagby_construct(graph, quiet = TRUE, size = "largest"))
})

test_that("recorded human labels and bare lone occasions retain source identities", {
  graph <- cd_test_graph("T -> A -> Y; T -> Y")
  record <- attr(graph, "query_record")
  record$node_map$quantity[record$node_map$node == "A"] <- "Canonical.A"
  record$node_map$member_id[record$node_map$node == "A"] <- "Canonical.A"
  record$node_map$occurrence[record$node_map$node == "A"] <- "Canonical.A_t2"
  record$policy$node_id[record$policy$node_id == "A"] <- "Canonical.A"
  attr(graph, "query_record") <- record
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_identical(dagby_options(dag)$quantities$direct$members, "Canonical.A")
  expect_identical(.cd_match_nodes(graph, "Canonical.A_t2"), "A")
})
