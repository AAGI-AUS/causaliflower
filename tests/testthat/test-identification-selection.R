################################################################################
# test-identification-selection.R


# Check public selection without changing the underlying causal query

test_that("strategy status selection preserves the separate identification verdicts", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] X -> M -> Y M <-> Y }"
  )
  ordinary <- get_identification(dag, quiet = TRUE)
  complete <- get_identification(dag, id_status = "all", quiet = TRUE)
  rows <- as.data.frame(ordinary)
  all_rows <- as.data.frame(complete)

  expect_true(all(rows$id_status == "identifiable"))
  expect_identical(attr(ordinary, "verdicts")$direct$graph_structure$id_status,
                   "zero_by_graph")
  expect_true("not_established" %in% all_rows$id_status)
  expect_identical(complete, get_identification(dag, id_status = NULL, quiet = TRUE))
  expect_identical(as.data.frame(summary(ordinary)), all_rows)
  expect_identical(attr(ordinary, "analyses"), attr(complete, "analyses"))
  expect_identical(attr(ordinary, "verdicts"), attr(complete, "verdicts"))

  zeros <- get_identification(dag, id_status = "zero_by_graph", quiet = TRUE)
  expect_equal(nrow(as.data.frame(zeros)), 0L)
  expect_identical(attr(zeros, "verdicts"), attr(complete, "verdicts"))
  failed <- get_identification(dag, id_status = "not_established", quiet = TRUE)
  expect_true(all(as.data.frame(failed)$id_status == "not_established"))
  expect_identical(attr(failed, "analyses"), attr(complete, "analyses"))
  expect_identical(attr(failed, "verdicts"), attr(complete, "verdicts"))

  mixed <- get_identification(
    dag, id_status = c("zero_by_graph", "not_established"), quiet = TRUE
  )
  expect_identical(unique(as.data.frame(mixed)$id_status), "not_established")
  expect_identical(attr(mixed, "verdicts"), attr(complete, "verdicts"))
})

test_that("empty strategies retain the query verdict and convert without invented rows", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] X -> Y X <-> Y }"
  )
  result <- get_identification(dag, effect = "total", quiet = TRUE)
  expect_s3_class(result, "causaliflower_identification")
  expect_identical(names(result), "total")
  expect_identical(result$total, list())
  expect_false(attr(result, "analyses")$total$identifiable)
  expect_identical(attr(result, "verdicts")$total$id$id_status, "not_identifiable")
  expect_equal(nrow(as.data.frame(result)), 0L)
  expect_equal(nrow(data.table::as.data.table(result)), 0L)
  expect_identical(names(as.data.frame(result)), c(
    "effect", "method", "assumptions", "variables", "id_status", "notes"
  ))
  expect_identical(attr(as.data.frame(result), "verdicts"), attr(result, "verdicts"))
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_gt(nchar(printed), 0L)
  expect_match(printed, "[Nn]ot identifiable")
  expect_false(grepl("\\[\\[|\\[1\\]|\\[,", printed))
  expect_gt(nrow(as.data.frame(summary(result))), 0L)
  expect_false("not_identifiable" %in% as.data.frame(summary(result))$id_status)
})

test_that("method selectors preserve complete alternatives within compatible effects", {
  dag <- dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] C1 -> C2 C2 -> X C1 -> Y",
    "X -> M -> Y X -> Y }"
  ))
  complete <- get_identification(dag, id_status = "all", quiet = TRUE)
  selected <- get_identification(
    dag, method = "single_door", id_status = "all", quiet = TRUE
  )
  expect_identical(names(selected), "direct")
  expect_null(names(selected$direct))
  expect_length(selected$direct, 2L)
  complete_direct <- Filter(function(entry) {
    return(identical(entry$method, "single_door"))
  }, complete$direct)
  expect_identical(selected$direct, complete_direct)
  expect_identical(as.data.frame(summary(selected)), as.data.frame(selected))

  natural <- get_identification(dag, method = "mediation_formula", quiet = TRUE)
  expect_identical(names(natural), c("direct", "indirect"))
  expect_setequal(as.data.frame(natural)$effect, c("direct", "indirect"))
  expect_true(all(as.data.frame(summary(natural))$method == "mediation_formula"))
  expect_false(any(as.data.frame(summary(natural))$method == "single_door"))

  methods <- get_identification(
    dag, effect = "total", method = c("adjustment", "front_door"),
    id_status = "all", quiet = TRUE
  )
  expect_setequal(as.data.frame(methods)$method, c("adjustment", "front_door"))
  expect_true(all(as.data.frame(summary(methods))$method %in% c("adjustment", "front_door")))
  expect_identical(names(summary(methods)), "total")
  expect_identical(attr(methods, "verdicts")$total, attr(complete, "verdicts")$total)
})

test_that("summary restores statuses within the selected causal scope", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] X -> M -> Y X -> Y M <-> Y }"
  )
  direct <- get_identification(
    dag, effect = "direct", method = "mediation_formula", quiet = TRUE
  )
  expect_identical(direct$direct, list())
  expanded <- summary(direct)
  expect_identical(names(expanded), "direct")
  expect_null(names(expanded$direct))
  rows <- as.data.frame(expanded)
  expect_identical(rows$method, "mediation_formula")
  expect_identical(rows$id_status, "not_established")
  expect_false(any(grepl("controlled", rows$method)))
})

test_that("incompatible and unknown selectors fail before returning a different question", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] X -> M -> Y X -> Y Z -> X }"
  )
  expect_error(get_identification(dag, effect = "total", method = "single_door"), "method")
  expect_error(get_identification(dag, effect = "direct", method = "adjustment"), "method")
  expect_error(get_identification(dag, effect = "indirect", method = "single_door"), "method")
  expect_error(get_identification(dag, conditioning = "M", method = "mediation_formula"))
  expect_error(get_identification(dag, experiments = "Z", method = "mediation_formula"))
  expect_error(get_identification(dag, mediator = "M", method = "single_door"), "mediator")
  expect_error(get_identification(dag, target = "natural_direct_effect"), "unused argument")

  for (value in list("", NA_character_, 1, list("adjustment"), "unknown")) {
    expect_error(get_identification(dag, method = value), "method")
    expect_error(get_identification(dag, id_status = value), "id_status")
  }
  for (value in c("id", "idc", "z_id", "graph_structure")) {
    expect_error(get_identification(dag, method = value), "method")
  }
  expect_error(get_identification(dag, method = c("adjustment", NA_character_)), "method")
  expect_error(get_identification(dag, id_status = c("identifiable", NA_character_)),
               "id_status")
})

test_that("mediator inputs define the same natural effects in both groups", {
  dag <- dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] X -> M1 -> Y",
    "X -> M2 -> Y X -> Y }"
  ))
  direct <- get_identification(
    dag, effect = "direct", method = "mediation_formula",
    mediator = c("M2", "M1", "M2"), quiet = TRUE
  )
  indirect <- get_identification(
    dag, effect = "indirect", mediator = c("M1", "M2"), quiet = TRUE
  )
  expect_identical(direct$direct[[1L]]$variables$mediators, c("M1", "M2"))
  expect_identical(indirect$indirect[[1L]]$variables$mediators, c("M1", "M2"))
  expect_identical(attr(direct, "analyses")$mediation,
                   attr(indirect, "analyses")$mediation)
})
