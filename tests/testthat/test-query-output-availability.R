test_that("the bundled slice makes no claims about individual output locators", {
  barley <- .query_apsim_outputs("barley")
  expect_equal(nrow(barley), 0L)
  expect_named(barley, c("crop", "node_id", "source_expression", "status", "report_tested", "evidence", "note"))
  expect_false(anyDuplicated(barley$node_id) > 0L)
  expect_true(all(barley$node_id %in% apsim_dag_nodes$node_id))
  expect_false(any(c("cgr_actual", "n_translocation", "canopy_conductance") %in% barley$node_id))
  expect_equal(nrow(.query_apsim_outputs("maize")), 0L)
  expect_equal(nrow(.query_apsim_outputs("wheat")), 0L)
  expect_equal(nrow(.query_apsim_outputs("rice")), 0L)
  expect_error(query_params(observe_apsim_outputs = NA), "TRUE or FALSE")
  expect_error(query_params(observe_apsim_outputs = "TRUE"), "TRUE or FALSE")
})

test_that("the output register filters crop, expression and verification status", {
  fixture_path <- testthat::test_path("fixtures", "output-availability.csv")
  read_outputs <- .query_apsim_outputs
  environment(read_outputs) <- list2env(list(system.file = function(...) { # Supply the explicit synthetic register
    return(fixture_path)
  }), parent = environment(.query_apsim_outputs))
  barley <- read_outputs("barley")
  expect_setequal(barley$node_id, c("left_quantity", "right_quantity"))
  expect_false(anyDuplicated(barley$node_id) > 0L)
  expect_true(all(barley$status == "source_verified_public"))
  expect_setequal(read_outputs("maize")$node_id, "maize_quantity")
  expect_equal(nrow(read_outputs("rice")), 0L)
})

test_that("hidden output assumption never overrides supplied columns", {
  v <- c("Mapped_t0", "Mapped_t1", "Unknown", "Visible", "Y")
  hidden <- setNames(c(TRUE, TRUE, TRUE, FALSE, TRUE), v)
  assign <- function(columns, matched = character(), show = TRUE, outputs = "Mapped") {
    .assign_latent(v, hidden, columns, matched, show,
                   always_observed = "Y", apsim_observed_ids = outputs)
  }
  expect_identical(assign(FALSE), "Unknown")
  expect_setequal(assign(FALSE, outputs = character()), c("Mapped_t0", "Mapped_t1", "Unknown"))
  expect_setequal(assign(TRUE, "Visible"), c("Mapped_t0", "Mapped_t1", "Unknown"))
  expect_setequal(assign(TRUE, "Mapped"), c("Unknown", "Visible"))
  expect_length(assign(FALSE, show = FALSE), 0L)
})

test_that("barley query assumes ordinary nodes obtainable and records exceptions", {
  query <- function(observe = TRUE, columns = NULL, label = NULL) {
    query_graph("barley", columns = columns, treatments = "N_application",
      outcomes = c("grain_N", "grain_wt"), time_unroll = "selective",
      randomised_treatment = FALSE, label = label, quiet = TRUE,
      params = query_params(observe_apsim_outputs = observe))
  }
  g <- query()
  old <- query(FALSE)
  expect_identical(attr(g, "edges"), attr(old, "edges"))
  expect_identical(names(g), names(old))
  policy <- attr(g, "apsim_output_assumption")$policy
  expected <- attr(old, "latent")[.strip_slice(attr(old, "latent")) %in% policy$node_id[policy$recovered_node == "yes"]]
  expect_setequal(setdiff(attr(old, "latent"), attr(g, "latent")), expected)
  expect_gt(length(expected), 0L)
  expect_length(attr(g, "latent"), 0L)
  expect_false(any(c("cgr_actual", "lodging") %in% .strip_slice(names(g))))
  expect_true("soil_thermal_regime" %in% .strip_slice(names(g)))
  expect_true(attr(g, "apsim_output_assumption")$applied)
  expect_length(attr(g, "apsim_output_assumption")$conflicts$node_id, 0L)
  expect_true(attr(g, "query_record")$arguments$params$observe_apsim_outputs)
  expect_setequal(attr(g, "apsim_output_mappings")$node_id, .query_apsim_outputs("barley")$node_id)
  expect_true(length(unique(.strip_slice(expected))) > nrow(attr(g, "apsim_output_mappings")))
  expect_false(any(.strip_slice(attr(g, "full_latent")) %in% .query_apsim_outputs("barley")$node_id))
  expect_false(any(attr(g, "rolled_latent") %in% .query_apsim_outputs("barley")$node_id))
  cols <- c("N_application", "grain_N", "grain_wt")
  expect_identical(attr(query(TRUE, cols), "dag_text"), attr(query(FALSE, cols), "dag_text"))
  expect_false(attr(query(TRUE, cols), "apsim_output_assumption")$applied)
  # Display aliases must not change the canonical availability decision.
  labelled <- query(label = "readable")
  expect_equal(length(attr(labelled, "latent")), length(attr(g, "latent")))
})

test_that("policy conflicts cannot be overridden by broad APSIM locators", {
  tab <- data.frame(node_id = c("ordinary", "interpreted", "nonapsim", "prohibited"),
    recovered_node = c("yes", "yes:user", "yes", "no:Not a quantity"),
    provenance = c("internal", "output", "not_apsim", "not_apsim"),
    apsim_obtainable = c("", "fixture:interpreted", "fixture:nonapsim", ""))
  p <- .query_output_policy(tab, "barley")
  expect_identical(p$recovered_node, c("yes", "yes:user", "yes:user", "no:Not a quantity"))
  expect_setequal(p$node_id[p$metadata_conflict], c("interpreted", "nonapsim"))
  expect_identical(.query_output_policy(tab, "maize")$recovered_node,
                   c("yes", "yes", "yes:user", "yes:user"))
})

test_that("construction folds saved interpretations by directed bypass", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> M -> Y; T -> S -> Y; U -> M; U -> Y }")
  ids <- names(q)
  rule <- setNames(rep("yes", length(ids)), ids)
  rule[c("S", "U")] <- "yes:user"
  attr(q, "query_record") <- list(arguments = list(crop = "maize", params = list(water_model = "SoilWater")),
    policy_applied = TRUE,
    policy = data.frame(node_id = ids, recovered_node = unname(rule)),
    node_map = data.frame(node = ids, quantity = ids, member_id = ids, occurrence = ids),
    requests = list(columns = list(resolved = character())))
  before <- serialize(q, NULL)
  p <- .cd_prepare(q, latents = NULL, keep_nodes = NULL)
  expect_setequal(p$ledger$node[!p$ledger$eligible], c("S", "U"))
  expect_setequal(names(p$analysis), c("T", "M", "Y"))
  e <- dagitty::edges(p$analysis)
  expect_true(any(e$v == "T" & e$w == "Y" & e$e == "->"))
  expect_false(any(e$e == "<->"))
  expect_identical(serialize(q, NULL), before)
  expect_true("U" %in% names(.cd_prepare(q, latents = NULL, keep_nodes = "U")$analysis))
})
