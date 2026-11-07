test_that("constructor retains source and typed history and prints only the graph", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U -> T; U -> M; U -> Y; T -> M -> Y }")
  attr(q, "author_note") <- "preserve"
  before <- serialize(q, NULL)
  g <- suppressWarnings(construct_dag(q, quiet = TRUE))
  expect_identical(serialize(original_query(g), NULL), before)
  expect_identical(paste(capture.output(print(g)), collapse = "\n"), as.character(g))
  h <- dag_history(g, details = TRUE)
  # The ranked route sets nothing before its search and ends by marking the
  # declared unmeasured variables (none here, so that step changes nothing).
  expect_identical(vapply(h, `[[`, character(1), "operation"),
    c("recover_nodes", "refine_adjustment_search", "select_direct_and_total", "trim_time_routes", "summarise_time", "mark_unmeasured"))
  expect_identical(h[[2]]$details, dag_options(g)$refinement)
  expect_true(all(unlist(dag_history(g)[6L, -(1:2)]) == 0L))
  expect_identical(.gd_signature(h[[1]]$before), .gd_signature(q))
  expect_identical(.gd_signature(tail(h, 1)[[1]]$after), .gd_signature(g))
  expect_null(attr(.gd_state(g)$checkpoint, "cf_provenance"))
  expect_identical(serialize(q, NULL), before)
})

test_that("restores preserve declared states and the common cause measurement rule", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> M -> Y; U [latent]; U -> M; U -> Y; V -> M; V -> Y }")
  # Declarations no longer remove options: the direct option {M, U, V} is
  # selected, so both declared common causes are retained and marked.
  g <- suppressWarnings(construct_dag(q, latents = c("U", "V"), keep_nodes = "M", quiet = TRUE))
  expect_setequal(dagitty::latents(g), c("U", "V"))
  expect_identical(restore_nodes(g, "U"), g)
  # The view of the total family omits both, a focused graph for manual editing.
  focused <- dag_set_view(g, "total", 1L)
  expect_false(any(c("U", "V") %in% names(focused)))
  restored <- restore_nodes(focused, "U")
  expect_true("U" %in% dagitty::latents(restored))
  expect_false("V" %in% names(restored))
  expect_identical(original_query(restored), q)
  confounded <- add_confounders(focused, node_name = "M", include = "V")
  expect_true("V" %in% dagitty::latents(confounded))
  expect_identical(as.character(undo_dag(restored)), as.character(focused))
  expect_warning(dag_checks(restored), "before later edits")
  changed <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y }")
  attr(changed, "cf_provenance") <- attr(g, "cf_provenance")
  expect_error(restore_nodes(changed, "U"), "outside the traced")
})

test_that("manual construction uses a restricted emitted query as its reference", {
  q <- query_graph("barley", treatments = "N_application", outcomes = "grain_N",
    randomised_treatment = TRUE, quiet = TRUE, params = query_params(query_depth = 1))
  g <- suppressWarnings(construct_dag(q, quiet = TRUE))
  expect_true(any(grepl("supplied query_depth", dag_options(g)$limitations)))
  expect_identical(original_query(g), q)
  expect_true(all(names(temporal_dag(g)) %in% names(q)))
})
test_that("restoring a node keeps saved grouped members and aliases", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> A -> Y; T -> Y; U -> A; U -> B -> Y; C -> Y }")
  ids <- names(q)
  map <- data.frame(node = ids, quantity = ids, member_id = ids, occurrence = ids)
  extra <- map[map$node == "A", ]
  extra$member_id <- "A2"
  map$member_id[map$node == "A"] <- "A1"
  map <- rbind(map, extra)
  attr(q, "query_record") <- list(arguments = list(crop = "fixture", params = list(water_model = "SoilWater")),
    node_map = map, aliases = list(U = "special_u"),
    requests = list(columns = list(resolved = character())))
  g <- suppressWarnings(construct_dag(q, quiet = TRUE))
  restored <- restore_nodes(g, "C")
  kept <- .cd_node_map(.gd_state(restored)$analysis)
  expect_setequal(kept$member_id[kept$node == "A"], c("A1", "A2"))
  # When U is omitted from the selected graph, its saved alias remains usable.
  expect_false("U" %in% names(g))
  alias <- add_confounders(g, node_name = "A", include = "special_u")
  expect_true("U" %in% names(alias))
})
