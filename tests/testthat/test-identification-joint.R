################################################################################
# test-identification-joint.R


# Check identification for joint endpoints and handling of unsupported methods

test_that("joint totals support one or several treatments and outcomes", {
  expect_joint_total <- function(dag, treatments, outcomes) { # Check one joint endpoint query
    before <- as.character(dag)
    result <- get_identification(dag, treatments, outcomes, quiet = TRUE)
    total <- attr(result, "analyses")$total
    expect_identical(names(attr(result, "analyses")), "total")
    expect_identical(attr(result, "scope"), "total")
    expect_identical(attr(result, "treatment"), treatments)
    expect_identical(attr(result, "outcome"), outcomes)
    expect_true(total$query$joint)
    expect_identical(total$query$intervention, treatments)
    expect_true(total$identifiable)
    expect_true(total$observationally_identifiable)
    expect_null(total$hedge)
    expect_identical(total$adjustment_sets, list(character()))
    expect_identical(total$canonical_adjustment, list(character()))
    expect_identical(names(result), "total")
    expect_null(names(result$total))
    expect_identical(vapply(result$total, `[[`, character(1L), "method"), "adjustment")
    expect_identical(attr(result, "verdicts")$total$id$id_status, "identifiable")
    expect_identical(as.character(dag), before)
    expect_identical(total$query$conditioning, character())
    expect_identical(total$query$fixing_set, character())
    return(invisible(NULL))
  }
  expect_joint_total(dagitty::dagitty("dag { X -> Y1 X -> Y2 }"), "X", c("Y1", "Y2"))
  expect_joint_total(dagitty::dagitty("dag { X1 -> Y X2 -> Y }"), c("X1", "X2"), "Y")
  expect_joint_total(dagitty::dagitty("dag { X1 -> Y1 X2 -> Y2 }"), c("X1", "X2"), c("Y1", "Y2"))
})

test_that("joint adjustment records complete sets", {
  dag <- dagitty::dagitty("dag {
    C -> X1 C -> X2 C -> Y1 C -> Y2 X1 -> Y1 X2 -> Y2
  }")
  result <- get_identification(dag, c("X1", "X2"), c("Y1", "Y2"), quiet = TRUE)
  total <- attr(result, "analyses")$total
  expect_true(total$identifiable)
  expect_identical(total$adjustment_sets, list("C"))
  expect_identical(total$canonical_adjustment, list("C"))
  expect_true(dagitty::isAdjustmentSet(dag, total$adjustment_sets[[1L]],
    exposure = c("X1", "X2"), outcome = c("Y1", "Y2")))
  expect_identical(total$roles$variable, "C")
  expect_identical(total$roles$role, "adjustment variable")
  rows <- as.data.frame(result)
  adjustment <- rows[rows$method == "adjustment", , drop = FALSE]
  expect_equal(nrow(adjustment), 1L)
  expect_identical(adjustment$variables[[1L]]$covariates, "C")
})

test_that("joint ID handles observed bidirected edges and explicit directed latents", {
  texts <- c(
    "dag { X -> Y1 X -> Y2 Y1 <-> Y2 }",
    "dag { X -> Y1 X -> Y2 U [latent] U -> Y1 U -> Y2 }"
  )
  for(text in texts){
    result <- get_identification(dagitty::dagitty(text), "X", c("Y1", "Y2"), quiet = TRUE)
    total <- attr(result, "analyses")$total
    expect_true(total$identifiable)
    expect_null(total$hedge)
    expect_identical(total$adjustment_sets, list(character()))
  }
  # ID can establish a joint distribution even though no adjustment set exists
  mediated <- dagitty::dagitty("dag {
    X -> M M -> Y1 M -> Y2 X <-> Y1 X <-> Y2
  }")
  result <- get_identification(mediated, "X", c("Y1", "Y2"), quiet = TRUE)
  total <- attr(result, "analyses")$total
  expect_true(total$identifiable)
  expect_identical(total$adjustment_sets, list())
  expect_identical(total$frontdoor_sets, list())
  expect_equal(nrow(as.data.frame(result)), 0L)
  expect_identical(result$total, list())
  expect_identical(attr(result, "verdicts")$total$id$id_status, "identifiable")
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "ID")
  expect_match(printed, "[Ii]dentifiable")
  expect_false(grepl("[Nn]ot identifiable", printed))
  full <- as.data.frame(summary(result))
  expect_identical(full$id_status[full$method == "front_door"], "not_assessed")
  expect_identical(full$id_status[full$method == "adjustment"], "not_established")
})

test_that("joint adjustment retains multiple complete alternatives", {
  dag <- dagitty::dagitty("dag {
    C1 -> C2 C2 -> X C1 -> Y1 C1 -> Y2 X -> Y1 X -> Y2
  }")
  result <- get_identification(dag, "X", c("Y1", "Y2"), quiet = TRUE)
  total <- attr(result, "analyses")$total
  expect_setequal(vapply(total$adjustment_sets, paste, character(1), collapse = ","), c("C1", "C2"))
  expect_equal(sum(as.data.frame(result)$method == "adjustment"), 2L)
  expect_true(all(vapply(total$adjustment_sets, function(adjustment_set) {
    return(dagitty::isAdjustmentSet(dag, adjustment_set,
      exposure = "X", outcome = c("Y1", "Y2")))
  }, logical(1))))
})

test_that("joint nonidentification retains a valid hedge", {
  expect_joint_hedge <- function(dag, treatments, outcomes) { # Check one nonidentified joint query
    result <- get_identification(dag, treatments, outcomes, quiet = TRUE)
    total <- attr(result, "analyses")$total
    expect_false(total$identifiable)
    expect_false(total$observationally_identifiable)
    expect_identical(total$adjustment_sets, list())
    expect_type(total$hedge, "list")
    expect_identical(result$total, list())
    expect_identical(attr(result, "verdicts")$total$id$id_status, "not_identifiable")
    expect_identical(attr(summary(result), "verdicts"), attr(result, "verdicts"))
    private <- identification_result(dag, treatments, outcomes)
    expect_valid_identification_hedge(private, dag)
    return(invisible(NULL))
  }
  expect_joint_hedge(dagitty::dagitty("dag { X -> Y1 X <-> Y1 X -> Y2 }"), "X", c("Y1", "Y2"))
  expect_joint_hedge(dagitty::dagitty("dag { X1 -> Y X1 <-> Y X2 -> Y }"), c("X1", "X2"), "Y")
  expect_joint_hedge(dagitty::dagitty("dag { X1 -> Y1 X1 <-> Y1 X2 -> Y2 }"), c("X1", "X2"), c("Y1", "Y2"))
})

test_that("joint summaries distinguish not assessed methods from unsuccessful checks", {
  result <- get_identification(dagitty::dagitty("dag { X -> Y1 X -> Y2 }"),
    "X", c("Y1", "Y2"), quiet = TRUE)
  total <- attr(result, "analyses")$total
  expect_identical(total$method_support$assessed, c("ID", "covariate adjustment"))
  expect_setequal(total$method_support$not_assessed,
    c("front-door", "instrumental variable", "direct effects", "mediation", "IDC", "z-ID"))
  full <- as.data.frame(summary(result))
  expect_identical(names(full), c(
    "effect", "method", "assumptions", "variables", "id_status", "notes"
  ))
  expect_true(all(c("front_door", "conditional_front_door", "instrumental_variable") %in% full$method))
  expect_true(all(full$id_status[
    full$method %in% c("front_door", "conditional_front_door", "instrumental_variable")
  ] == "not_assessed"))
  expect_false(any(full$effect %in% c("direct", "indirect")))
  expect_identical(total$frontdoor_sets, list())
  expect_identical(total$instruments, list())
  expect_true(is.na(total$z_identifiable))
  expect_match(paste(capture.output(print(result)), collapse = "\n"), "Y1, Y2", fixed = TRUE)
  expect_identical(causaliflower:::id_query_label(total$query), "P(Y1, Y2 | do(X))")
})

test_that("joint endpoint order and duplicates preserve the same query", {
  dag <- dagitty::dagitty("dag {
    X1 [exposure] X2 [exposure] Y1 [outcome] Y2 [outcome]
    X1 -> Y1 X2 -> Y2
  }")
  tagged <- get_identification(dag, quiet = TRUE)
  reordered <- get_identification(dag, c("X2", "X1", "X2"),
    c("Y2", "Y1", "Y1"), effect = "total", quiet = TRUE)
  expect_equal(attr(tagged, "analyses"), attr(reordered, "analyses"))
  expect_equal(as.data.frame(tagged), as.data.frame(reordered))
  expect_identical(attr(reordered, "treatment"), c("X1", "X2"))
  expect_identical(attr(reordered, "outcome"), c("Y1", "Y2"))
  expect_message(get_identification(dag), "X1, X2 -> Y1, Y2", fixed = TRUE)
})

test_that("joint endpoints reject latent missing invalid and overlapping selections", {
  dag <- dagitty::dagitty("dag { X1 -> Y1 X2 -> Y2 U [latent] U -> Y1 }")
  for(x in list(c("X1", NA_character_), c("X1", ""), c("X1", "missing"), c("X1", "U"))){
    expect_error(get_identification(dag, x, c("Y1", "Y2"), quiet = TRUE), "Invalid treatment set")
  }
  for(y in list(c("Y1", NA_character_), c("Y1", ""), c("Y1", "missing"), c("Y1", "U"))){
    expect_error(get_identification(dag, c("X1", "X2"), y, quiet = TRUE), "Invalid outcome set")
  }
  expect_error(get_identification(dag, character(), c("Y1", "Y2")), "No treatment detected")
  expect_error(get_identification(dag, c("X1", "X2"), character()), "No outcome detected")
  expect_error(get_identification(dag, factor(c("X1", "X2")), "Y1"), "Invalid treatment")
  expect_error(get_identification(dag, c("X1", "X2"), c("X2", "Y1")), "Overlapping")
  expect_error(get_identification(dag, c("X1", "X2"), c("X1", "X2")), "Overlapping")
})

test_that("explicit unsupported joint requests fail without changing the estimand", {
  dag <- dagitty::dagitty("dag { X -> M M -> Y1 M -> Y2 Z }")
  base <- list(dag = dag, treatment = "X", outcome = c("Y1", "Y2"), quiet = TRUE)
  for(extra in list(list(effect = "direct"), list(mediator = "M"),
      list(conditioning = "Z"), list(experiments = "Z"),
      list(conditioning = character()), list(experiments = character()),
      list(mediator = character()))){
    expect_error(do.call(get_identification, c(base, extra)), "Joint endpoint queries support only")
  }
  expect_error(do.call(get_identification, c(base, list(effect = "indirect"))),
    "Joint endpoints, conditioning and experiments support only total effect queries.",
    fixed = TRUE)
  expect_error(identification_result(dag, "X", c("Y1", "Y2"), effect = "direct"),
    "Joint endpoint queries support only")
})

test_that("joint roles do not mislabel paths interrupted by another intervention", {
  dag <- dagitty::dagitty("dag { X1 -> M M -> X2 X2 -> Y A }")
  result <- get_identification(dag, c("X1", "X2"), "Y", quiet = TRUE)
  total <- attr(result, "analyses")$total
  expect_true(total$identifiable)
  expect_identical(total$adjustment_sets, list(character()))
  expect_identical(total$roles$role[total$roles$variable == "M"], "other")
  expect_false(any(total$roles$role %in% c("mediator", "instrument", "outcome parent")))
  expect_setequal(total$roles$variable, c("M", "A"))
})

test_that("joint structural zero and the scalar contract remain distinct", {
  dag <- dagitty::dagitty("dag { X Y1 Y2 Y1 <-> Y2 }")
  result <- get_identification(dag, "X", c("Y1", "Y2"), quiet = TRUE)
  expect_identical(attr(result, "verdicts")$total$graph_structure$id_status, "zero_by_graph")
  expect_false("zero_by_graph" %in% as.data.frame(result)$id_status)
  expect_true(attr(result, "analyses")$total$identifiable)
  expect_identical(attr(result, "analyses")$total$adjustment_sets, list(character()))
  joint_rows <- as.data.frame(summary(result))
  expect_identical(joint_rows$id_status[joint_rows$method == "front_door"], "not_assessed")
  expect_identical(joint_rows$id_status[joint_rows$method == "instrumental_variable"], "not_assessed")
  scalar <- get_identification(dagitty::dagitty("dag { X -> Y }"), "X", "Y", quiet = TRUE)
  expect_null(attr(scalar, "analyses")$total$method_support)
  expect_null(attr(scalar, "analyses")$total$query$joint)
  expect_setequal(names(attr(scalar, "analyses")),
                  c("total", "controlled", "mediation", "direct", "conditional_frontdoor"))
  expect_identical(attr(scalar, "scope"), c("total", "direct", "indirect"))
  scalar_rows <- as.data.frame(summary(scalar))
  not_assessed <- scalar_rows[scalar_rows$id_status == "not_assessed", , drop = FALSE]
  expect_identical(unique(not_assessed$method), "mediation_formula")
  expect_setequal(not_assessed$effect, c("direct", "indirect"))
  expect_identical(names(attr(scalar, "verdicts")$total), "id")
})
