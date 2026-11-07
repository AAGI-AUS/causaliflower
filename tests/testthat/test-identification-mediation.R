################################################################################
# test-identification-mediation.R


# Check mediation formulas and mediator validation

# Mediation formula checks use the private query records.
# Explicit get_identification() calls test public argument handling.
mediation_entries <- function(result) { # Read the private mediation strategy entries
  return(attr(result, "mediation"))
}

mediation_entry <- function(result, mediators) { # Select one mediator set without losing alternatives
  entries <- mediation_entries(result)
  matches <- vapply(entries, function(entry) {
    return(identical(entry$mediators, sort(unique(mediators))))
  }, logical(1))
  expect_equal(sum(matches), 1L)
  return(entries[[which(matches)]])
}

expect_mediation_formula <- function(dag, result, mediators, adjustment) { # Check both adjustment conditions for the mediation formula
  entry <- mediation_entry(result, mediators)
  expect_true(entry$feasible)
  expect_length(entry$adjustment_sets, 1L)
  expect_setequal(entry$adjustment_sets[[1L]], adjustment)
  # These are the two distinct adjustment queries of Shpitser and
  # VanderWeele (2011), Theorems 1-2, not a total-effect adjustment check.
  expect_true(dagitty::isAdjustmentSet(
    dag, entry$adjustment_sets[[1L]], exposure = "T", outcome = mediators
  ))
  expect_true(dagitty::isAdjustmentSet(
    dag, entry$adjustment_sets[[1L]], exposure = c("T", mediators),
    outcome = "Y"
  ))
  return(invisible(NULL))
}

test_that("direct queries report a separate natural-effects strategy", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> M M -> Y T -> Y }"
  )
  result <- identification_result(dag, effect = "direct")
  expect_identification_shape(result)
  expect_equal(nrow(result), 2L)
  expect_identical(result$strategy[[1L]], "general identification")
  expect_identical(result$variable_set[[1L]], "fixing set: {M}")
  expect_true(attr(result, "identifiable"))
  expect_true(attr(result, "observationally_identifiable"))

  row <- identification_strategy(result, "^mediation \\(natural effects\\)$")
  expect_true(row$feasible)
  expect_identical(row$tier, "assumption-dependent")
  expect_match(row$variable_set, "M", fixed = TRUE)
  expect_match(row$variable_set, "{}", fixed = TRUE)

  entries <- mediation_entries(result)
  expect_length(entries, 1L)
  entry <- entries[[1L]]
  expect_named(entry, c("mediators", "feasible", "adjustment_sets", "note",
                        "assumptions"))
  expect_identical(entry$mediators, "M")
  expect_type(entry$note, "character")
  expect_length(entry$note, 1L)
  expect_type(entry$assumptions, "character")
  expect_length(entry$assumptions, 1L)
  expect_match(entry$assumptions, "cross-world")
  expect_match(entry$assumptions, "consistency")
  expect_match(entry$assumptions, "positivity")
  expect_mediation_formula(dag, result, "M", character())
  expect_identical(
    identification_result(dag, effect = "direct", mediator = "M"), result
  )
})

test_that("mediation uses one covariate set for both adjustment queries", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome]",
    "T -> M M -> Y T -> Y",
    "X -> T X -> M Z -> M Z -> Y W -> T W -> Y }"
  ))
  result <- identification_result(dag, effect = "direct", mediator = "M")
  expect_mediation_formula(dag, result, "M", c("W", "X", "Z"))
  # Neither relation's own adjustment set suffices for the other relation.
  expect_false(dagitty::isAdjustmentSet(
    dag, c("W", "Z"), exposure = "T", outcome = "M"
  ))
  expect_false(dagitty::isAdjustmentSet(
    dag, "X", exposure = c("T", "M"), outcome = "Y"
  ))
})

test_that("latent confounding of each required relation prevents the formula", {
  relations <- c("T <-> M", "M <-> Y", "T <-> Y")
  for (edge in relations) {
    dag <- dagitty::dagitty(paste(
      "dag { T [exposure] Y [outcome] T -> M M -> Y T -> Y", edge, "}"
    ))
    result <- identification_result(dag, effect = "direct", mediator = "M")
    entry <- mediation_entry(result, "M")
    expect_false(entry$feasible, info = edge)
    expect_length(entry$adjustment_sets, 0L)
    row <- identification_strategy(result, "mediation")
    expect_false(row$feasible, info = edge)
    expect_identical(row$tier, "assumption-dependent")
  }

  explicit_latent <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent]",
    "T -> M M -> Y T -> Y U -> M U -> Y }"
  ))
  expect_false(mediation_entry(
    identification_result(explicit_latent, effect = "direct"), "M"
  )$feasible)
})

test_that("treatment-induced mediator-outcome confounding is target-specific", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> L L -> M M -> Y",
    "L -> Y T -> Y }"
  ))
  result <- identification_result(dag, effect = "direct")
  expect_false(mediation_entry(result, "M")$feasible)
  expect_length(mediation_entry(result, "M")$adjustment_sets, 0L)
  expect_mediation_formula(dag, result, "L", character())
  expect_mediation_formula(dag, result, c("L", "M"), character())

  # Without L -> Y, L no longer confounds M and Y: a treatment-induced
  # ancestor of the mediator is not by itself a reason to reject mediation.
  serial <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> L L -> M M -> Y T -> Y }"
  ))
  serial_result <- identification_result(serial, effect = "direct", mediator = "M")
  expect_mediation_formula(serial, serial_result, "M", character())
})

test_that("automatic queries include singleton and joint observed mediators", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> M1 M1 -> Y T -> M2 M2 -> Y",
    "T -> Y T -> D X -> Y I }"
  ))
  result <- identification_result(dag, effect = "direct")
  entries <- mediation_entries(result)
  expect_length(entries, 3L)
  expect_equal(sum(result$strategy == "mediation (natural effects)"), 3L)
  expect_mediation_formula(dag, result, "M1", character())
  expect_mediation_formula(dag, result, "M2", character())
  expect_mediation_formula(dag, result, c("M1", "M2"), character())
  expect_false(any(c("D", "X", "I") %in% unlist(lapply(
    entries, function(entry) { return(entry$mediators) }
  ), use.names = FALSE)))

  selected <- identification_result(
    dag, effect = "direct", mediator = c(second = "M2", first = "M1")
  )
  expect_length(mediation_entries(selected), 1L)
  expect_identical(mediation_entries(selected)[[1L]]$mediators, c("M1", "M2"))
  expect_mediation_formula(dag, selected, c("M1", "M2"), character())
})

test_that("confounding within a joint mediator vector is not endpoint confounding", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> M1 M1 -> Y T -> M2 M2 -> Y",
    "T -> Y M1 <-> M2 }"
  ))
  result <- identification_result(dag, effect = "direct")
  expect_false(mediation_entry(result, "M1")$feasible)
  expect_false(mediation_entry(result, "M2")$feasible)
  expect_mediation_formula(dag, result, c("M1", "M2"), character())
})

test_that("mediation reporting does not require a direct treatment-outcome arrow", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> M M -> Y }"
  )
  result <- identification_result(dag, effect = "direct")
  expect_true(result$feasible[[1L]])
  expect_true(attr(result, "identifiable"))
  expect_mediation_formula(dag, result, "M", character())

  # A front-door total effect can be identified when the standard mediation
  # adjustment formula is unavailable. The latter is not a verdict on ID.
  frontdoor <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> M M -> Y T <-> Y }"
  )
  expect_true(attr(identification_result(frontdoor), "identifiable"))
  expect_false(mediation_entry(
    identification_result(frontdoor, effect = "direct"), "M"
  )$feasible)
})

test_that("latent projection retains observed mediator choices only", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] L [latent]",
    "T -> L L -> M M -> Y T -> Y }"
  ))
  result <- identification_result(dag, effect = "direct")
  expect_length(mediation_entries(result), 1L)
  expect_mediation_formula(dag, result, "M", character())

  latent_only <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] L [latent] T -> L L -> Y }"
  )
  empty <- identification_result(latent_only, effect = "direct")
  expect_identical(mediation_entries(empty), list())
  expect_equal(nrow(empty), 1L)
})

test_that("natural-effect feasibility does not replace the controlled-effect verdict", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> M M -> Y T -> Y T <-> M }"
  ))
  result <- identification_result(dag, effect = "direct")
  expect_false(mediation_entry(result, "M")$feasible)
  expect_true(result$feasible[[1L]])
  expect_true(attr(result, "identifiable"))
  expect_true(attr(result, "observationally_identifiable"))
  expect_identical(attr(result, "z_identifiable"), NA)
  expect_null(attr(result, "hedge"))
  expect_identical(attr(result, "effect"), "direct")
  roles <- attr(result, "roles")
  expect_equal(nrow(roles), 1L)
  expect_identical(roles$variable, "M")
  expect_identical(roles$role, "outcome parent")

  # Conversely, intervening on every other outcome parent can fail even
  # though mediation through M has a valid observational formula. The
  # unrelated parent A is not an intervention in that mediation target.
  other_parent <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> M M -> Y T -> Y A -> Y A <-> Y }"
  ))
  other_result <- identification_result(other_parent, effect = "direct")
  expect_false(other_result$feasible[[1L]])
  expect_false(attr(other_result, "identifiable"))
  expect_false(attr(other_result, "observationally_identifiable"))
  expect_type(attr(other_result, "hedge"), "list")
  expect_mediation_formula(other_parent, other_result, "M", character())
})

test_that("mediation leaves total, conditional, and experimental output unchanged", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] X -> T X -> Y T -> M M -> Y T -> Y }"
  ))
  for (result in list(
    identification_result(dag),
    identification_result(dag, conditioning = "X"),
    identification_result(dag, experiments = "X")
  )) {
    expect_false(any(grepl("mediation", result$strategy)))
    expect_null(mediation_entries(result))
    expect_identical(attr(result, "effect"), "total")
  }
  no_mediator <- identification_result(
    dagitty::dagitty("dag { T [exposure] Y [outcome] T -> Y }"),
    effect = "direct"
  )
  expect_identical(mediation_entries(no_mediator), list())
  expect_equal(nrow(no_mediator), 1L)

  # Appending mediator must not reinterpret any existing positional input.
  positional <- suppressMessages(get_identification(
    dag, NULL, NULL, "direct", NULL, NULL, TRUE
  ))
  expect_identical(
    positional, get_identification(dag, effect = "direct", quiet = TRUE)
  )
})

test_that("explicit mediators must be observed internal causal-path nodes", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] L [latent] T -> L L -> M M -> Y",
    "T -> Y X -> T X -> Y T -> D I }"
  ))
  invalid <- list(character(), NA_character_, "", 1, list("M"),
                  c("M", NA_character_), "T", "Y", "L", "X", "D", "I",
                  "missing", c("M", "X"))
  for (value in invalid) {
    expect_error(
      get_identification(dag, effect = "direct", mediator = value, quiet = TRUE),
      "mediator"
    )
  }
  overview <- get_identification(dag, mediator = "M", quiet = TRUE)
  expect_identical(
    attr(overview, "analyses")$mediation,
    mediation_entries(identification_result(
      dag, effect = "direct", mediator = "M"
    ))
  )
  expect_error(get_identification(dag, effect = "total", mediator = "M",
                                  quiet = TRUE), "direct")
  expect_error(get_identification(dag, mediator = character(), quiet = TRUE),
               "mediator")
  expect_error(get_identification(dag, effect = "direct", mediator = "M",
                                  conditioning = "X", quiet = TRUE))
  expect_error(get_identification(dag, effect = "direct", mediator = "M",
                                  experiments = "X", quiet = TRUE))
})

test_that("mediator choices are stable under node order and duplicate inputs", {
  first <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> M1 M1 -> M2 M2 -> Y T -> Y }"
  ))
  second <- dagitty::dagitty(paste(
    "dag { Y [outcome] M2 M1 T [exposure] T -> Y M2 -> Y M1 -> M2 T -> M1 }"
  ))
  expect_identical(
    identification_result(first, effect = "direct"),
    identification_result(second, effect = "direct")
  )
  expect_identical(
    identification_result(first, effect = "direct", mediator = c("M2", "M1", "M2")),
    identification_result(first, effect = "direct", mediator = c("M1", "M2"))
  )
})

test_that("mediation row notes separate the result from its assumptions", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> M1 M1 -> Y T -> M2 M2 -> Y",
    "T -> Y M2 <-> Y }"
  ))
  result <- identification_result(dag, effect = "direct")
  rows <- result[result$strategy == "mediation (natural effects)"]
  entries <- mediation_entries(result)
  expect_equal(nrow(rows), length(entries))
  expect_identical(rows$feasible, c(TRUE, FALSE, FALSE))
  expect_match(rows$note[[1L]], "adjustment checks. Requires a", fixed = TRUE)
  expect_match(
    rows$note[[2L]],
    "no covariate set of observed treatment nondescendants satisfies both",
    fixed = TRUE
  )
  for (index in seq_along(entries)) {
    expect_identical(
      rows$note[[index]],
      paste0(entries[[index]]$note, ". ", entries[[index]]$assumptions)
    )
  }
})

test_that("a failed common-adjustment search stays polynomial", {
  # Each Ai -> Bi -> T, Ai -> M pair gives T -> M two blockers, so that query
  # has 2^12 minimal sets. The joint query needs K, which opens
  # T <- U1 -> K <- U2 -> M, so no common set exists.
  index <- seq_len(12L)
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U1 [latent] U2 [latent]",
    "T -> M M -> Y T -> Y U1 -> T U1 -> K U2 -> K U2 -> M K -> Y",
    paste(sprintf("A%d -> B%d B%d -> T A%d -> M A%d -> Y",
                  index, index, index, index, index), collapse = " "),
    "}"
  ))
  checks <- 0L
  is_adjustment_set <- dagitty::isAdjustmentSet
  testthat::local_mocked_bindings(
    isAdjustmentSet = function(...) {
      checks <<- checks + 1L
      is_adjustment_set(...)
    },
    .package = "dagitty"
  )
  result <- identification_result(dag, effect = "direct")
  entry <- mediation_entry(result, "M")
  expect_false(entry$feasible)
  expect_identical(entry$adjustment_sets, list())
  # Two candidates with two checks each; enumeration needed thousands.
  expect_lte(checks, 4L)

  # With no eligible covariates the empty set is the only candidate.
  checks <- 0L
  unadjusted <- identification_result(dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> M M -> Y T -> Y M <-> Y }"
  ), effect = "direct")
  expect_false(mediation_entry(unadjusted, "M")$feasible)
  expect_lte(checks, 2L)
})
