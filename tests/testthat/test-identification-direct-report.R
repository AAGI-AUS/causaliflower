################################################################################
# test-identification-direct-report.R


# Keep the original coefficient separate from the retained parent intervention query

expect_direct_report <- function(dag, expected_verdict) { # Check direct output and its legacy diagnostic
  result <- get_identification(dag, effect = "direct", quiet = TRUE)
  expanded <- summary(result)
  direct <- attr(result, "verdicts")$direct$graph_structure
  controlled <- attr(result, "analyses")$controlled

  expect_identical(names(result), "direct")
  expect_identical(result$direct, list())
  expect_identical(direct$id_status, "zero_by_graph")
  expect_identical(direct$method, "graph_structure")
  expect_identical(controlled$identifiable, expected_verdict)
  expect_identical(attr(expanded, "analyses")$controlled, controlled)
  expect_identical(attr(expanded, "verdicts")$direct$graph_structure, direct)
  expect_false(any(grepl("controlled", as.data.frame(result)$method)))
  expect_false(any(grepl("controlled", as.data.frame(expanded)$method)))

  ordinary <- paste(capture.output(print(result)), collapse = "\n")
  detailed <- paste(capture.output(print(result, details = TRUE)), collapse = "\n")
  expect_match(ordinary, "[Zz]ero")
  expect_false(grepl("[Cc]ontrolled intervention", ordinary))
  expect_false(grepl("[Cc]ontrolled intervention", detailed))

  selected <- result["direct"]
  expect_type(selected, "list")
  expect_false(inherits(selected, "causaliflower_identification"))
  expect_identical(selected$direct, result$direct)
  expect_null(attr(selected, "analyses"))
  expect_null(attr(selected, "verdicts"))
  selected$direct <- list(list(notes = "author edit"))
  expect_false(identical(selected$direct, result$direct))

  overview <- get_identification(dag, quiet = TRUE)
  explicit_overview <- get_identification(dag, effect = NULL, quiet = TRUE)
  expect_identical(overview, explicit_overview)
  expect_identical(attr(overview, "verdicts")$direct$graph_structure, direct)
  expect_identical(attr(overview, "analyses")$controlled, controlled)
  return(invisible(NULL))
}

test_that("direct output keeps the zero coefficient and an identifiable parent query distinct", {
  dag <- dagitty::dagitty(
    "dag { L [latent] T [exposure] Y [outcome] T -> L L -> Y }"
  )
  expect_direct_report(dag, TRUE)
})

test_that("direct output keeps the zero coefficient and a failed parent query distinct", {
  dag <- dagitty::dagitty(paste(
    "dag { L [latent] U [latent] T [exposure] Y [outcome]",
    "T -> L L -> Y U -> T U -> Y }"
  ))
  expect_direct_report(dag, FALSE)
})

test_that("direct reports retain every alternative coefficient adjustment set", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] C1 -> C2 C2 -> X C1 -> Y X -> Y }"
  )
  result <- get_identification(
    dag, effect = "direct", method = "single_door", quiet = TRUE
  )
  alternatives <- result$direct
  expect_length(alternatives, 2L)
  expect_identical(vapply(alternatives, `[[`, character(1L), "method"),
                   rep("single_door", 2L))
  covariates <- lapply(alternatives, function(alternative) {
    return(alternative$variables$covariates)
  })
  expect_setequal(unlist(covariates, use.names = FALSE), c("C1", "C2"))
  expect_true(all(lengths(covariates) == 1L))
  expect_identical(summary(result)$direct, alternatives)
  expect_false(any(grepl("controlled", as.data.frame(summary(result))$method)))
})
