################################################################################
# test-identification-overview.R


# Check effect groups and the scientific records behind the public report

overview_result <- function(dag, ...) { # Request a report without console messages
  return(get_identification(dag, ..., quiet = TRUE))
}

overview_set_signature <- function(sets) { # Compare sets independently of their order
  if (is.null(sets)) {
    return(character())
  }
  return(unname(sort(vapply(unclass(sets), function(set) {
    return(paste(sort(as.character(set)), collapse = "\r"))
  }, character(1L)))))
}

test_that("the default separates total direct and indirect effect groups", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] C -> X C -> Y X -> M -> Y X -> Y }"
  )
  result <- overview_result(dag)
  analyses <- attr(result, "analyses")
  rows <- as.data.frame(result)

  expect_s3_class(result, "causaliflower_identification")
  expect_type(result, "list")
  expect_false(inherits(result, "data.frame"))
  expect_identical(names(result), c("total", "direct", "indirect"))
  expect_null(names(result$total))
  expect_null(names(result$direct))
  expect_null(names(result$indirect))
  expect_identical(vapply(result$total, `[[`, character(1L), "method"), "adjustment")
  expect_identical(vapply(result$direct, `[[`, character(1L), "method"),
                   c("single_door", "mediation_formula"))
  expect_identical(vapply(result$indirect, `[[`, character(1L), "method"),
                   "mediation_formula")
  expect_identical(attr(result, "effect"), "overview")
  expect_identical(attr(result, "scope"), c("total", "direct", "indirect"))
  expect_identical(result, overview_result(dag, effect = NULL))
  expect_true(all(rows$id_status %in% c("identifiable", "zero_by_graph")))
  expect_false(any(grepl("controlled", rows$method)))

  expect_identical(analyses$total$query$treatment, "X")
  expect_identical(analyses$total$query$outcome, "Y")
  expect_identical(analyses$total$query$intervention, "X")
  expect_true(analyses$total$identifiable)
  expect_true(analyses$direct$direct_arrow)
  expect_identical(analyses$direct$treatment, "X")
  expect_identical(analyses$direct$outcome, "Y")
  expect_setequal(analyses$controlled$query$fixing_set, c("C", "M"))
  expect_identical(analyses$mediation[[1L]]$mediators, "M")

  expect_null(attr(result, "identifiable", exact = TRUE))
  expect_null(attr(result, "observationally_identifiable", exact = TRUE))
})

test_that("effect selection keeps only the requested effect group", {
  dag <- dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] C -> X C -> Y",
    "X -> M -> Y X -> Y Z -> X }"
  ))
  total <- overview_result(dag, effect = "total")
  direct <- overview_result(dag, effect = "direct")
  indirect <- overview_result(dag, effect = "indirect")

  expect_identical(names(total), "total")
  expect_identical(names(direct), "direct")
  expect_identical(names(indirect), "indirect")
  expect_null(names(indirect$indirect))
  expect_true(all(vapply(indirect$indirect, function(entry) {
    return(identical(entry$method, "mediation_formula"))
  }, logical(1L))))
  expect_identical(attr(total, "scope"), "total")
  expect_identical(attr(direct, "scope"), "direct")
  expect_identical(attr(indirect, "scope"), "indirect")
  expect_null(attr(total, "analyses")$controlled)
  expect_null(attr(direct, "analyses")$total)
  expect_identical(overview_result(dag, effect = c(input = "total")), total)
  expect_identical(overview_result(dag, effect = c(input = "direct")), direct)
  expect_identical(overview_result(dag, effect = c(input = "indirect")), indirect)

  for (result in list(total, direct, indirect)) {
    expect_null(attr(result, "identifiable", exact = TRUE))
    expect_null(attr(result, "observationally_identifiable", exact = TRUE))
  }

  conditional <- overview_result(dag, conditioning = "M")
  expect_identical(names(conditional), "total")
  expect_identical(attr(conditional, "scope"), "total")
  expect_identical(attr(conditional, "analyses")$total$query$conditioning, "M")
  expect_identical(conditional,
                   overview_result(dag, effect = "total", conditioning = "M"))
  expect_equal(nrow(as.data.frame(conditional)), 0L)
  expect_identical(attr(conditional, "verdicts")$total$idc$id_status, "identifiable")

  experimental <- overview_result(dag, experiments = "Z")
  expect_identical(names(experimental), "total")
  expect_identical(attr(experimental, "analyses")$total$experiments, "Z")
  expect_identical(experimental,
                   overview_result(dag, effect = "total", experiments = "Z"))

  expect_error(overview_result(dag, effect = "direct", conditioning = "M"))
  expect_error(overview_result(dag, effect = "direct", experiments = "Z"))
  expect_error(overview_result(dag, effect = "indirect", conditioning = "M"))
  expect_error(overview_result(dag, effect = "indirect", experiments = "Z"))
  expect_error(overview_result(dag, conditioning = "M", experiments = "Z"))
  expect_error(overview_result(dag, effect = "total", mediator = "M"))
  expect_error(overview_result(dag, conditioning = "C", mediator = "M"))
  expect_error(overview_result(dag, effect = c("total", "direct")))
  expect_error(overview_result(dag, effect = "unknown"))

  selected <- overview_result(dag, mediator = "M")
  expect_identical(attr(selected, "effect"), "overview")
  expect_identical(attr(selected, "analyses")$mediation[[1L]]$mediators, "M")
  expect_setequal(attr(selected, "analyses")$controlled$query$fixing_set,
                  c("C", "M"))
})

test_that("single door sets match dagitty and minimal_sets on the supplied graph", {
  graphs <- c(
    fully_observed = "dag { X [exposure] Y [outcome] C -> X C -> Y X -> M -> Y X -> Y }",
    latent_confounding = "dag { X [exposure] Y [outcome] U [latent] U -> X U -> Y X -> Y }",
    bidirected = "dag { X [exposure] Y [outcome] X <-> C C -> Y X -> Y }"
  )
  for (name in names(graphs)) {
    dag <- dagitty::dagitty(graphs[[name]])
    expected <- dagitty::adjustmentSets(
      dag, exposure = "X", outcome = "Y", effect = "direct", type = "minimal"
    )
    wrapped <- suppressMessages(minimal_sets(
      dag, treatment = "X", outcome = "Y", effect = "direct", num_sets = 100L
    ))
    result <- overview_result(dag, effect = "direct", id_status = "all")
    direct <- attr(result, "analyses")$direct
    rows <- as.data.frame(result)
    rows <- rows[rows$method == "single_door", , drop = FALSE]

    expect_true(direct$direct_arrow, info = name)
    expect_identical(overview_set_signature(direct$adjustment_sets),
                     overview_set_signature(expected), info = name)
    expect_identical(overview_set_signature(direct$adjustment_sets),
                     overview_set_signature(wrapped), info = name)
    expect_identical(direct$feasible, length(expected) > 0L, info = name)
    expected_status <- "not_found"
    expected_result <- "not_established"
    if (length(expected) > 0L) {
      expected_status <- "available"
      expected_result <- "identifiable"
    }
    expect_identical(direct$status, expected_status, info = name)
    expect_equal(nrow(rows), max(1L, length(expected)), info = name)
    expect_true(all(rows$id_status == expected_result), info = name)
    expect_false(any(rows$id_status == "zero_by_graph"), info = name)
  }
})

test_that("a missing original direct arrow remains a visible zero coefficient", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] X -> M -> Y M <-> Y }"
  )
  result <- overview_result(dag)
  analyses <- attr(result, "analyses")
  direct <- attr(result, "verdicts")$direct$graph_structure

  expect_true(analyses$total$identifiable)
  expect_false(analyses$direct$direct_arrow)
  expect_identical(analyses$direct$status, "structural_zero")
  expect_identical(direct$id_status, "zero_by_graph")
  expect_identical(direct$method, "graph_structure")
  expect_identical(analyses$direct$feasible, NA)
  expect_false(analyses$controlled$identifiable)
  expect_false(analyses$controlled$query$direct_arrow)
  expect_setequal(analyses$controlled$query$fixing_set, "M")
  expect_false(analyses$mediation[[1L]]$feasible)
  expect_false(any(as.data.frame(result)$effect == "direct" &
                   as.data.frame(result)$method == "mediation_formula"))
  expect_identical(result$indirect, list())

  expanded <- as.data.frame(summary(result))
  natural <- expanded[expanded$method == "mediation_formula", , drop = FALSE]
  expect_identical(natural$id_status, rep("not_established", 2L))
  expect_setequal(natural$effect, c("direct", "indirect"))
  expect_false(any(grepl("controlled", expanded$method)))
})

test_that("a latent mediator keeps the original coefficient distinct from projection", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] L [latent] X -> L -> Y }"
  )
  result <- overview_result(dag)
  analyses <- attr(result, "analyses")

  expect_false(analyses$direct$direct_arrow)
  expect_identical(analyses$direct$status, "structural_zero")
  expect_true(analyses$controlled$query$direct_arrow)
  expect_identical(analyses$controlled$query$fixing_set, character())
  expect_true(analyses$controlled$identifiable)
  expect_true(analyses$total$identifiable)
  expect_length(analyses$mediation, 0L)
  expect_false(any(as.data.frame(result)$effect == "direct" &
                   as.data.frame(result)$method == "mediation_formula"))
  expect_identical(result$indirect, list())
  expect_null(attr(result, "identifiable", exact = TRUE))
})

test_that("failed parent interventions remain diagnostic without changing natural results", {
  dag <- dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] X -> M -> Y X -> Y",
    "C -> Y C <-> Y }"
  ))
  result <- overview_result(dag, mediator = "M")
  analyses <- attr(result, "analyses")

  expect_true(analyses$total$identifiable)
  expect_false(analyses$controlled$identifiable)
  expect_setequal(analyses$controlled$query$intervention, c("X", "C", "M"))
  expect_true(analyses$mediation[[1L]]$feasible)
  expect_identical(analyses$mediation[[1L]]$mediators, "M")
  natural <- as.data.frame(result)
  natural <- natural[natural$method == "mediation_formula", , drop = FALSE]
  expect_identical(natural$effect, c("direct", "indirect"))
  expect_identical(natural$id_status, rep("identifiable", 2L))
  expect_false(any(grepl("controlled", as.data.frame(summary(result))$method)))
  expect_null(attr(result, "identifiable", exact = TRUE))
})

test_that("public queries preserve the supplied graph and its attributes", {
  dag <- dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] U [latent]",
    "U -> X U -> Y X -> M -> Y X -> Y }"
  ))
  attr(dag, "author_record") <- list(note = "preserve supplied graph", revision = 3L)
  before <- serialize(dag, NULL)
  overview_result(dag)
  overview_result(dag, effect = "total")
  overview_result(dag, effect = "direct", mediator = "M")
  overview_result(dag, effect = "indirect", mediator = "M")
  expect_identical(serialize(dag, NULL), before)
})
