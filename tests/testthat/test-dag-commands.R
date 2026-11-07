command_fixture <- function(text, keep, observed = character(), simplify = FALSE) {
  original <- dagitty::dagitty(text)
  ids <- .gd_node_ids(original)
  omitted <- setdiff(names(original), keep)
  projected <- if (length(omitted)) collapse_nodes(original, omitted) else original
  surface <- .gd_recover(original, character(), apply_policy = FALSE)
  state <- list(original_query = original, analysis = original,
    temporal = projected, retained = names(projected), summary = NULL,
    focus = unique(unname(.gd_node_ids(projected))),
    controls = list(crop = "fixture", water_model = "SoilWater", disambiguate = NULL,
                    simplify = simplify, keep_common_causes = FALSE),
    build = list(recovery = surface, observed_ids = observed,
                 availability = list(mode = "scm", statement = "Fixture")),
    history = list(), undo = list(), view_id = "fixture")
  .gd_attach(projected, state)
}

test_that("restoring a common cause explains only its own projected edges", {
  g <- command_fixture("dag { T [exposure]; Y [outcome]; T -> M -> Y; U -> M; U -> Y; V -> M; V -> Y }",
                       c("T", "M", "Y"))
  before <- serialize(g, NULL)
  preview <- NULL
  expect_output(preview <- restore_nodes(g, "U", dry_run = TRUE), "U")
  expect_identical(serialize(g, NULL), before)
  restored <- restore_nodes(g, "U")
  expect_true("U" %in% names(restored))
  expect_false("U" %in% dagitty::latents(restored)) # ordinary restore preserves status
  edges <- dagitty::edges(restored)
  expect_true(any(edges$e == "<->")) # V is still omitted
  both <- restore_nodes(restored, "V")
  expect_false(any(dagitty::edges(both)$e == "<->"))
  expect_identical(original_query(both), original_query(g))
  back <- undo_dag(both)
  expect_setequal(names(back), names(restored))
  expect_setequal(edge_key(dagitty::edges(back)), edge_key(dagitty::edges(restored)))
  expect_setequal(names(undo_dag(back)), names(g))
  expect_error(undo_dag(g), "no editing command")
  selected <- add_confounders(g, node_name = "M", include = "U")
  expect_true("U" %in% names(selected))
  expect_false("V" %in% names(selected))
  expect_error(add_confounders(g, node_name = "M", include = "T"), "do not match")
})

test_that("common cause selection rejects a serial path and reports witnesses", {
  g <- command_fixture("dag { T [exposure]; Y [outcome]; T -> A -> Y; U -> A; W -> A; W -> Y }",
                       c("T", "A", "Y"))
  out <- add_confounders(g, node_name = "A")
  expect_true("W" %in% names(out))
  expect_false("U" %in% names(out)) # U -> A -> Y is not a fork
  expect_equal(tail(.gd_state(out)$history, 1)[[1]]$details$witnesses$cause, "W")
  by_role <- add_confounders(g, role = "outcome")
  expect_setequal(names(out), names(by_role))
  expect_error(add_confounders(g, role = "outcome", node_name = "A"), "Choose role")
  expect_error(add_confounders(g, role = "made_up"), "Unknown")
  expect_error(add_confounders(g, node_name = "missing"), "absent")
})

test_that("supplied columns determine status of restored common causes", {
  g <- command_fixture("dag { T [exposure]; Y [outcome]; T -> M -> Y; U -> M; U -> Y; V -> M; V -> Y }",
                       c("T", "M", "Y"), observed = "U")
  hidden <- add_confounders(g, node_name = "M")
  expect_true("V" %in% names(hidden))
  expect_false("U" %in% names(hidden))
  expect_true("V" %in% dagitty::latents(hidden))
  all <- add_confounders(g, node_name = "M", unmeasured = FALSE)
  expect_setequal(intersect(c("U", "V"), dagitty::latents(all)), "V")
  direct <- restore_nodes(g, "U")
  expect_false("U" %in% dagitty::latents(direct))
  direct_v <- restore_nodes(g, "V")
  expect_false("V" %in% dagitty::latents(direct_v)) # no measurement assumption
})

test_that("restoration respects explicit recovery policy", {
  g <- command_fixture("dag { nitrogen_applied [exposure]; grain_n [outcome]; nitrogen_applied -> canopy_senescence -> grain_n; unspecified_organ_live_n -> grain_n }",
                       c("nitrogen_applied", "grain_n"))
  state <- .gd_state(g)
  state$controls$crop <- "barley"
  state$build$recovery <- .gd_recover(original_query(g), character(), apply_policy = TRUE)
  state$analysis <- state$build$recovery$graph
  g <- .gd_attach(g, state)
  expect_true("canopy_senescence" %in% names(restore_nodes(g, "canopy_senescence")))
  expect_error(restore_nodes(g, "unspecified_organ_live_n"), "recovery policy")
})

test_that("time occurrences alone do not count as two distinct confounded quantities", {
  g <- command_fixture("dag { T [exposure]; Y [outcome]; T -> A_t2 -> A_t3 -> Y; U -> A_t2; U -> A_t3 }",
                       c("T", "A_t2", "A_t3", "Y"))
  expect_message(out <- add_confounders(g, node_name = "A_t2"), "No omitted")
  expect_false("U" %in% names(out))
})

test_that("dry run, repeat calls, history and protection from stale edits work", {
  g <- command_fixture("dag { T [exposure]; Y [outcome]; T -> M -> Y; U -> M; U -> Y }",
                       c("T", "M", "Y"))
  initial <- serialize(g, NULL)
  expect_output(p <- add_confounders(g, role = "outcome", dry_run = TRUE), "U")
  expect_identical(serialize(g, NULL), initial)
  edited <- add_confounders(g, role = "outcome")
  expect_s3_class(edited, "dagitty")
  expect_identical(restore_nodes(edited, "U"), edited)
  expect_output(undo_dag(edited, dry_run = TRUE), "U")
  expect_true(nrow(dag_history(edited)) >= 1L)
  # Native external edits may keep attributes. Such provenance must not be
  # silently used to overwrite that edit during reference reconstruction.
  changed <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y }")
  attr(changed, "cf_provenance") <- attr(edited, "cf_provenance")
  class(changed) <- c("cf_dag", "dagitty")
  expect_error(restore_nodes(changed, "U"), "changed|edited|match")
  expect_error(temporal_dag(changed), "changed|edited|match")
})

test_that("exact temporal targets and include requests remain exact", {
  g <- command_fixture("dag { T [exposure]; Y [outcome]; T -> A_t2 -> A_t3 -> Y; U -> A_t3; U -> Y }",
                       c("T", "A_t2", "A_t3", "Y"))
  expect_message(out <- add_confounders(g, node_name = "A_t2"), "No omitted")
  expect_false("U" %in% names(out))
  expect_true("U" %in% names(add_confounders(g, node_name = "A_t3")))
  g <- command_fixture("dag { T [exposure]; Y [outcome]; T -> A -> Y; U_t1 -> T; U_t2 -> A; U_t2 -> Y }",
                       c("T", "A", "Y"))
  expect_error(add_confounders(g, node_name = "A", include = c("U_t1", "U_t2")), "U_t1")
  expect_true("U_t2" %in% names(add_confounders(g, node_name = "A", include = "U")))
  phase <- command_fixture("dag { T [exposure]; Y [outcome]; T -> A_pre_treatment -> A_post_treatment -> Y; U -> A_post_treatment; U -> Y }",
                           c("T", "A_pre_treatment", "A_post_treatment", "Y"))
  expect_message(out <- add_confounders(phase, node_name = "A_pre_treatment"), "No omitted")
  expect_false("U" %in% names(out))
})
