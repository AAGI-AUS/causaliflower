# A synthetic subsystem tests temporal query mechanics using the source structure.
# Package tables are supplied through a private lexical environment.
# No namespace, global data or scientific helper is replaced.
phase_history_query <- function(treatment = "sowing_depth", outcome = "root_depth",
                                treatment_stage = "0", outcome_stage = "3",
                                randomised = FALSE, columns = NULL, depth = 2L,
                                window = "origin", recombine = FALSE, label = NULL,
                                confounded = FALSE, show = TRUE,
                                hidden_history = FALSE, observe_outputs = TRUE) {
  ids <- c("sowing_depth", "root_depth", "root_area", "sw_mm")
  nodes <- apsim_dag_nodes[match(ids, apsim_dag_nodes$node_id), , drop = FALSE]
  nodes$concept_id <- paste0("c:", ids)
  attr(nodes, "cf_source_order") <- NULL
  nodes$crop_type <- "barley"
  nodes$function_of <- c("", "c:root_depth@t-1", "c:root_depth", "")
  nodes$composition_of <- ""
  nodes$init_of <- c("", "c:sowing_depth", "", "")
  nodes$path_dependent <- c(FALSE, TRUE, FALSE, FALSE)
  if (hidden_history) nodes$hidden_node[2L] <- TRUE
  if (confounded) {
    nodes$function_of[c(1L, 2L)] <- c("c:sw_mm", "c:root_depth@t-1;c:sw_mm")
  }
  tables <- list(apsim_dag_nodes = nodes, node_aliases = node_aliases[ids],
                 management_defaults = list(core = "sowing_depth"),
                 crop_phases = crop_phases[crop_phases$crop == "barley" & crop_phases$phase_index <= 3L, ],
                 phase_of = phase_of[FALSE, ])
  q <- query_graph
  environment(q) <- list2env(list(getExportedValue = function(ns, name) {
    if (identical(ns, "causaliflower") && name %in% names(tables)) tables[[name]] else base::getExportedValue(ns, name)
  }), parent = environment(query_graph))
  q("barley", treatments = treatment, outcomes = outcome, columns = columns,
    randomised_treatment = randomised, label = label, quiet = TRUE,
    params = query_params(query_depth = depth, treatment_growth_stage = treatment_stage,
                         outcome_growth_stage = outcome_stage, phase_window = window,
                         node_recombination = recombine, show_hidden_nodes = show,
                         observe_apsim_outputs = observe_outputs))
}

test_that("a measured state retains its initialisation and earlier state chain", {
  g <- phase_history_query()
  expect_setequal(names(g), c("sowing_depth", paste0("root_depth_t", 0:3)))
  expect_setequal(paste(attr(g, "edges")$from, attr(g, "edges")$to),
                  c("sowing_depth root_depth_t0", paste(paste0("root_depth_t", 0:2), paste0("root_depth_t", 1:3))))
  expect_identical(dagitty::outcomes(g), "root_depth_t3")
  expect_identical(dagitty::exposures(g), "sowing_depth")
  expect_true(dagitty::isAcyclic(g))
  expect_setequal(attr(g, "query_record")$node_map$node, names(g))
})

test_that("earlier outcome states are not automatically observed measurements", {
  g <- phase_history_query(columns = "sowing_depth")
  history <- paste0("root_depth_t", 0:2)
  expect_setequal(attr(g, "latent"), history)
  expect_true(all(history %in% attr(g, "full_latent")))
  expect_false("root_depth_t3" %in% attr(g, "full_latent"))
  expect_identical(dagitty::outcomes(g), "root_depth_t3")
})

test_that("treatment and outcome with the same base name retain distinct roles and precise surgery", {
  observational <- phase_history_query(treatment = "root_depth", treatment_stage = "1", confounded = TRUE)
  randomised <- phase_history_query(treatment = "root_depth", treatment_stage = "1", confounded = TRUE, randomised = TRUE)
  for (g in list(observational, randomised)) {
    expect_identical(dagitty::exposures(g), "root_depth_t1")
    expect_identical(dagitty::outcomes(g), "root_depth_t3")
    expect_true(dagitty::isAcyclic(g))
  }
  before <- attr(observational, "full_edges")
  after <- attr(randomised, "full_edges")
  expect_true(any(before$to == "root_depth_t1"))
  expect_false(any(after$to == "root_depth_t1"))
  expect_true(any(after$from == "root_depth_t1" & after$to == "root_depth_t2"))
  expect_true(any(after$to == "root_depth_t0"))
  expect_true(any(after$to == "root_depth_t3" & grepl("sw_mm", after$from, fixed = TRUE)))
})

test_that("joint outcome and treatment anchors do not relabel intermediate states", {
  joint <- phase_history_query(outcome_stage = c("1", "3"))
  expect_setequal(dagitty::outcomes(joint), c("root_depth_t1", "root_depth_t3"))
  expect_false("root_depth_t2" %in% dagitty::outcomes(joint))
  treatment <- phase_history_query(treatment = "root_depth", treatment_stage = c("0", "1"), randomised = TRUE)
  expect_setequal(dagitty::exposures(treatment), c("root_depth_t0", "root_depth_t1"))
  expect_identical(dagitty::outcomes(treatment), "root_depth_t3")
  expect_true(any(attr(treatment, "edges")$from == "root_depth_t1" & attr(treatment, "edges")$to == "root_depth_t2"))
})

test_that("compact history respects its boundary and aliases preserve endpoint roles", {
  g <- phase_history_query(treatment = "root_depth", treatment_stage = "2", window = "compact")
  expect_false(any(grepl("_t0$", names(g))))
  expect_false(any(attr(g, "edges")$from == "sowing_depth"))
  plain <- phase_history_query(recombine = TRUE)
  labelled <- phase_history_query(recombine = TRUE, label = "readable")
  m <- attr(labelled, "query_record")$node_map
  expect_setequal(m$occurrence[m$node %in% dagitty::outcomes(labelled)], dagitty::outcomes(plain))
  expect_setequal(m$occurrence[m$node %in% dagitty::exposures(labelled)], dagitty::exposures(plain))
  expect_true(dagitty::isAcyclic(labelled))
})

test_that("an earlier outcome can remain a common cause without an invented bypass", {
  g <- phase_history_query(treatment = "root_depth", treatment_stage = "1",
                          columns = "root_area", confounded = TRUE)
  expect_true(all(c("root_depth_t0", "root_depth_t2") %in% attr(g, "latent")))
  expect_true(any(grepl("sw_mm", names(g), fixed = TRUE)))
  expect_true(all(dagitty::edges(g)$e == "->"))
  expect_identical(dagitty::outcomes(g), "root_depth_t3")
})

test_that("hidden history uses existing fold policy without losing the selected endpoint", {
  g <- phase_history_query(hidden_history = TRUE, show = FALSE, observe_outputs = FALSE)
  expect_identical(dagitty::outcomes(g), "root_depth")
  expect_setequal(names(g), c("sowing_depth", "root_depth"))
  expect_identical(paste(attr(g, "edges")$from, attr(g, "edges")$to), "sowing_depth root_depth")
  expect_setequal(attr(g, "recovery_folding")$hidden_removed, paste0("root_depth_t", 0:2))
})

test_that("the endpoint exemption does not observe hidden history without columns", {
  g <- phase_history_query(hidden_history = TRUE, observe_outputs = FALSE)
  expect_setequal(attr(g, "latent"), paste0("root_depth_t", 0:2))
  expect_true(all(paste0("root_depth_t", 0:2) %in% attr(g, "full_latent")))
  expect_false("root_depth_t3" %in% attr(g, "full_latent"))
})
