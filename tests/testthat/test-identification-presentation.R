################################################################################
# test-identification-presentation.R


# Check nested reports, selections, conversions and escaped output

identification_print_text <- function(result, ...) { # Capture a printed report
  return(paste(capture.output(print(result, ...)), collapse = "\n"))
}

identification_html <- function(result, ...) { # Render through knitr's HTML method
  previous <- knitr::opts_knit$get("rmarkdown.pandoc.to")
  on.exit(knitr::opts_knit$set(rmarkdown.pandoc.to = previous))
  knitr::opts_knit$set(rmarkdown.pandoc.to = "html")
  printed <- capture.output(rendered <- knitr::knit_print(result, ...))
  if (inherits(rendered, "knit_asis")) {
    return(paste(as.character(rendered), collapse = "\n"))
  }
  return(paste(printed, collapse = "\n"))
}

identification_method <- function(result, method) { # Select flattened rows for one method
  rows <- as.data.frame(result)
  return(rows[rows$method == method, , drop = FALSE])
}

test_that("a zero coefficient does not display a verdict for the parent intervention", {
  for (confounding in c("", "M <-> Y")) {
    dag <- dagitty::dagitty(paste(
      "dag { X [exposure] Y [outcome] X -> M -> Y", confounding, "}"
    ))
    result <- get_identification(dag, effect = "direct", quiet = TRUE)
    analyses <- attr(result, "analyses")
    controlled <- analyses$controlled
    query <- controlled$query

    expect_identical(controlled$observationally_identifiable,
                     identical(confounding, ""))
    expect_false(query$direct_arrow)
    expect_true(query$direct_path)
    expect_setequal(query$intervention, c("X", "M"))
    expect_identical(query$fixing_set, "M")
    expect_identical(analyses$direct$status, "structural_zero")
    expect_identical(analyses$direct$feasible, NA)
    expect_null(attr(result, "identifiable", exact = TRUE))
    expect_false(any(grepl("controlled", as.data.frame(summary(result))$method)))

    printed <- identification_print_text(result)
    expect_match(printed, "[Zz]ero")
    expect_false(grepl("[Cc]ontrolled", printed))
    expect_false(grepl("P\\(", printed))
    expect_false(grepl("[Cc]ontrolled",
                       identification_print_text(result, details = TRUE)))
  }
})

test_that("IV candidates remain separate from nonparametric identification", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] Z -> X -> Y X <-> Y }"
  )
  result <- get_identification(dag, effect = "total", quiet = TRUE)
  total <- attr(result, "analyses")$total
  expanded <- summary(result)
  general <- attr(expanded, "verdicts")$total$id
  iv <- identification_method(expanded, "instrumental_variable")

  expect_identical(result$total, list())
  expect_identical(general$method, "id")
  expect_identical(general$id_status, "not_identifiable")
  expect_false(total$identifiable)
  expect_false("id" %in% as.data.frame(expanded)$method)
  expect_identical(iv$id_status, "candidate")
  expect_match(paste(iv$assumptions, collapse = " "), "[Ll]inear")
  expect_identical(total$instruments[[1L]]$instrument, "Z")
  expect_identical(total$instruments[[1L]]$conditioning, character())
  expect_identical(iv$variables[[1L]]$instrument, "Z")
  expect_identical(iv$variables[[1L]]$conditioning, character())

  printed <- identification_print_text(expanded)
  expect_match(printed, "[Nn]onparametric")
  expect_match(printed, "[Ll]inear")
  expect_null(attr(result, "identifiable", exact = TRUE))
})

test_that("conditional roles retain the requested query and IDC movement", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] X -> M -> Y X -> Y }"
  )
  result <- get_identification(dag, conditioning = "M", quiet = TRUE)
  total <- attr(result, "analyses")$total
  query <- total$query
  roles <- total$roles
  guidance <- paste(roles[variable == "M", conditioning_guidance], collapse = " ")
  rows <- as.data.frame(result)

  expect_true(total$identifiable)
  expect_equal(nrow(rows), 0L)
  expect_identical(attr(result, "verdicts")$total$idc$id_status, "identifiable")
  expect_identical(query$conditioning, "M")
  expect_identical(query$conditioning_intervened, "M")
  expect_identical(query$conditioning_remaining, character())
  expect_match(guidance, "[Cc]ondition")
  expect_match(guidance, "IDC|[Ii]nterven")
  expect_false(grepl("do not include.*adjustment set for the total effect", guidance))
  expect_identical(attr(result, "verdicts")$total$idc$variables$conditioning, "M")
  printed <- identification_print_text(result)
  expect_false(grepl("do not include.*adjustment set for the total effect", printed))
})

test_that("natural mediator selection does not relabel the parent intervention diagnostic", {
  dag <- dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] X -> M -> Y X -> Y",
    "C -> Y C <-> Y }"
  ))
  result <- get_identification(dag, effect = "direct", mediator = "M", quiet = TRUE)
  analyses <- attr(result, "analyses")
  query <- analyses$controlled$query
  natural <- analyses$mediation[[1L]]

  expect_false(analyses$controlled$identifiable)
  expect_setequal(query$fixing_set, c("C", "M"))
  expect_setequal(query$intervention, c("X", "C", "M"))
  expect_identical(natural$mediators, "M")
  expect_true(natural$feasible)
  expect_identical(natural$adjustment_sets, list(character()))
  expect_identical(names(result), "direct")
  expect_null(result$indirect)

  printed <- identification_print_text(result)
  expect_false(grepl("[Cc]ontrolled", printed))
  expect_match(printed, "[Nn]atural")
  expect_match(printed, "[Mm]ediator")
  expect_false(grepl("do(X, C, M)", printed, fixed = TRUE))
})

test_that("structured strategies retain complete sets independently of display", {
  frontdoor <- get_identification(dagitty::dagitty(
    "dag { X [exposure] Y [outcome] X -> M -> Y X <-> Y }"
  ), effect = "total", quiet = TRUE)
  total <- attr(frontdoor, "analyses")$total
  expect_identical(total$frontdoor_sets, list("M"))
  expect_identical(total$query$intervention, "X")
  expect_identical(total$query$fixing_set, character())

  experiment <- get_identification(dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] Z -> X -> Y",
    "X <-> Z Z <-> Y }"
  )), experiments = "Z", quiet = TRUE)
  total <- attr(experiment, "analyses")$total
  expect_false(total$observationally_identifiable)
  expect_true(total$identifiable)
  expect_identical(total$experiments, "Z")
  expect_identical(total$experiment_set, "Z")
  experimental <- attr(experiment, "verdicts")$total$z_id
  expect_identical(experimental$variables$experiments, "Z")
  expect_identical(experimental$variables$experiment_set, "Z")
})

test_that("nested alternatives and table conversions preserve typed variable roles", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] C -> X C -> Y X -> M -> Y X -> Y }"
  )
  result <- get_identification(dag, quiet = TRUE)
  before <- serialize(result, NULL)
  rows <- as.data.frame(result)
  table <- data.table::as.data.table(result)

  expect_identical(names(rows), c(
    "effect", "method", "assumptions", "variables", "id_status", "notes"
  ))
  expect_s3_class(rows, "data.frame")
  expect_s3_class(table, "data.table")
  expect_equal(as.data.frame(table), rows)
  expect_type(rows$variables, "list")
  expect_gt(nrow(rows), 0L)

  for (effect in names(result)) {
    alternatives <- result[[effect]]
    expect_type(alternatives, "list")
    expect_null(names(alternatives))
    expect_gt(length(alternatives), 0L)
    for (alternative in alternatives) {
      expect_identical(names(alternative),
                       c("method", "assumptions", "variables", "id_status", "notes"))
      expect_type(alternative$method, "character")
      expect_length(alternative$method, 1L)
      expect_type(alternative$assumptions, "character")
      expect_type(alternative$variables, "list")
      expect_true(all(vapply(alternative$variables, is.character, logical(1L))))
      expect_type(alternative$id_status, "character")
      expect_length(alternative$id_status, 1L)
      expect_type(alternative$notes, "character")
    }
  }
  expect_identical(attr(rows, "verdicts"), attr(result, "verdicts"))
  expect_false(any(rows$method %in% c("id", "idc", "z_id", "graph_structure")))

  adjustment <- identification_method(result, "adjustment")
  expect_identical(adjustment$variables[[1L]]$covariates, "C")
  natural <- identification_method(result, "mediation_formula")
  expect_identical(natural$variables[[1L]]$mediators, "M")
  expect_identical(natural$variables[[1L]]$covariates, "C")
  expect_identical(serialize(result, NULL), before)
})

test_that("printing and summary respect edits and independent effect subsets", {
  dag <- dagitty::dagitty(
    "dag { X [exposure] Y [outcome] C -> X C -> Y X -> Y }"
  )
  result <- get_identification(dag, quiet = TRUE)
  before <- serialize(result, NULL)
  selected <- result["direct"]
  expect_type(selected, "list")
  expect_false(inherits(selected, "causaliflower_identification"))
  expect_null(attr(selected, "analyses"))
  expect_identical(selected$direct, result$direct)
  selected$direct[[1L]]$variables$covariates <- "author choice"
  expect_identical(result$direct[[1L]]$variables$covariates, "C")

  edited <- result
  edited$total[[1L]]$notes <- "author edit"
  expect_match(identification_print_text(edited), "author edit", fixed = TRUE)
  expect_identical(summary(edited)$total, edited$total)
  expect_null(attr(as.data.frame(edited), "verdicts"))
  expect_identical(serialize(result, NULL), before)

  printed <- identification_print_text(result)
  expect_match(printed, "[Tt]otal")
  expect_match(printed, "[Dd]irect")
  expect_match(printed, "[Mm]ethod")
  expect_match(printed, "[Aa]ssumptions")
  expect_match(printed, "[Ss]tatus")
  expect_false(grepl("\\[\\[|\\[,", printed))
  expect_false(grepl("References:|https?://", printed))
  expect_identical(identification_print_text(result, details = TRUE),
                   identification_print_text(summary(result)))
  expect_length(testthat::capture_messages(get_identification(dag, quiet = TRUE)), 0L)
})

test_that("complete query failures remain distinct from unsuccessful restricted strategies", {
  expect_report <- function(dag, name, ...) { # Check reported statuses against retained analyses
    result <- get_identification(dag, ..., id_status = "all", quiet = TRUE)
    before <- serialize(result, NULL)
    rows <- as.data.frame(result)
    analyses <- attr(result, "analyses")
    expect_false(any(grepl("controlled", rows$method)))
    expect_false(any(grepl("P\\(", rows$method)))

    if (identical(analyses$direct$status, "structural_zero")) {
      direct <- attr(result, "verdicts")$direct$graph_structure
      expect_identical(direct$id_status, "zero_by_graph", info = name)
    }
    if (name == "frontdoor") {
      expect_true(analyses$total$identifiable)
      expect_identical(identification_method(result, "front_door")$id_status,
                       "identifiable")
      expect_identical(analyses$total$adjustment_sets, list())
      expect_identical(analyses$total$instruments, list())
    }
    if (name == "instrument") {
      expect_false(analyses$total$identifiable)
      expect_identical(attr(result, "verdicts")$total$id$id_status,
                       "not_identifiable")
      iv <- identification_method(result, "instrumental_variable")
      expect_identical(iv$id_status, "candidate")
      expect_match(paste(iv$assumptions, collapse = " "), "[Ll]inear")
      expect_identical(iv$variables[[1L]]$instrument, "Z")
    }
    if (name == "general_id") {
      expect_true(analyses$total$identifiable)
      expect_identical(analyses$total$frontdoor_sets, list())
      expect_identical(analyses$total$adjustment_sets, list())
      expect_identical(attr(result, "verdicts")$total$id$id_status, "identifiable")
      expect_false(any(rows$id_status == "identifiable"))
    }
    if (name == "failed_experiments") {
      expect_false(analyses$total$identifiable)
      expect_identical(attr(result, "verdicts")$total$z_id$id_status,
                       "not_identifiable")
    }
    if (name %in% c("experiments", "failed_experiments")) {
      experimental <- attr(result, "verdicts")$total$z_id
      expect_match(paste(experimental$notes, collapse = " "), "[Oo]bservational")
      expect_match(paste(experimental$notes, collapse = " "), "[Ee]xperimental")
    }
    if (name %in% c("failed_mediation", "zero_direct")) {
      natural <- identification_method(result, "mediation_formula")
      expect_gt(nrow(natural), 0L)
      expect_true(all(natural$id_status == "not_established"))
      expect_false(any(natural$id_status == "not_identifiable"))
    }
    if (name == "zero_direct") {
      adjustment <- identification_method(result, "adjustment")
      expect_identical(adjustment$id_status, "identifiable")
      expect_identical(adjustment$variables[[1L]]$covariates, character())
    }
    if (!is.null(analyses$total) && !length(analyses$total$query$conditioning)) {
      expect_true(all(c("adjustment", "front_door", "conditional_front_door",
                        "instrumental_variable") %in% rows$method))
      expect_false(any(rows$method %in% c("id", "idc", "z_id", "graph_structure")))
      restricted <- rows[rows$method %in% c("adjustment", "front_door", "conditional_front_door", "instrumental_variable"),
                         , drop = FALSE]
      expect_false(any(restricted$id_status == "not_identifiable"))
    }
    expect_identical(serialize(result, NULL), before)
    return(invisible(NULL))
  }

  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] C -> X C -> Y X -> M -> Y X -> Y }"),
                "overview")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] X -> M -> Y X <-> Y }"),
                "frontdoor", effect = "total")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] Z -> X -> Y X <-> Y }"),
                "instrument", effect = "total")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] X -> M -> Y X <-> M }"),
                "conditional", conditioning = "M")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] X -> Z1 Z1 -> Z2 Z1 -> Y Z2 -> Y X -> Y X <-> Z2 Z1 <-> Y }"),
                "general_id", effect = "total")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] Z -> X -> Y X <-> Z Z <-> Y }"),
                "experiments", experiments = "Z")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] X -> Y X <-> Y Z }"),
                "failed_experiments", experiments = "Z")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] X -> M -> Y X -> Y M <-> Y }"),
                "failed_mediation", effect = "direct")
  expect_report(dagitty::dagitty("dag { X [exposure] Y [outcome] X -> M -> Y M <-> Y }"),
                "zero_direct")
})

test_that("all adjustment alternatives remain separate records", {
  result <- get_identification(dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] X -> Y",
    "A -> B -> C -> D -> E -> F -> X A -> Y }"
  )), effect = "total", quiet = TRUE)
  before <- serialize(result, NULL)
  sets <- attr(result, "analyses")$total$adjustment_sets
  expect_identical(sets, as.list(LETTERS[1:6]))
  adjustment <- identification_method(result, "adjustment")
  expect_equal(nrow(adjustment), 6L)
  covariates <- lapply(adjustment$variables, `[[`, "covariates")
  expect_identical(covariates, as.list(LETTERS[1:6]))
  expect_identical(identification_method(summary(result), "adjustment"), adjustment)
  expect_identical(attr(result, "analyses")$total$adjustment_sets, sets)
  expect_identical(serialize(result, NULL), before)
})

test_that("each mediator choice has its own direct and indirect records", {
  dag <- dagitty::dagitty(paste0(
    "dag { X [exposure] Y [outcome] X -> Y ",
    paste(sprintf("X -> M%d -> Y", seq_len(6L)), collapse = " "), " }"
  ))
  result <- get_identification(dag, effect = "direct", quiet = TRUE)
  before <- serialize(result, NULL)
  entries <- attr(result, "analyses")$mediation
  expect_length(entries, 7L)
  natural <- identification_method(result, "mediation_formula")
  expect_equal(nrow(natural), 7L)
  expect_true(all(natural$effect == "direct"))
  expect_true(all(natural$id_status == "identifiable"))
  expect_length(unique(natural$variables), 7L)
  expect_true(all(grepl("[Cc]ross world|cross-world", natural$assumptions)))
  expect_identical(serialize(result, NULL), before)

  indirect <- get_identification(dag, effect = "indirect", quiet = TRUE)
  indirect_rows <- identification_method(indirect, "mediation_formula")
  expect_identical(indirect_rows$variables, natural$variables)
  expect_identical(indirect_rows$id_status, natural$id_status)
  expect_true(all(indirect_rows$effect == "indirect"))

  mixed <- get_identification(dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] X -> M1 -> Y",
    "X -> M2 -> Y X -> Y M2 <-> Y }"
  )), effect = "direct", quiet = TRUE)
  entries <- attr(mixed, "analyses")$mediation
  expect_identical(vapply(entries, `[[`, logical(1L), "feasible"),
                   c(TRUE, FALSE, FALSE))
  expect_identical(identification_method(mixed, "mediation_formula")$id_status,
                   "identifiable")
  all_natural <- identification_method(summary(mixed), "mediation_formula")
  expect_identical(all_natural$id_status,
                   c("identifiable", "not_established", "not_established"))
})

test_that("HTML escapes values and displays the same visible alternatives", {
  skip_if_not(requireNamespace("knitr", quietly = TRUE))
  dag <- dagitty::dagitty(paste(
    'dag { "T<unsafe&" [exposure] "Y&value" [outcome]',
    '"C<b>&" -> "T<unsafe&" "C<b>&" -> "Y&value"',
    '"T<unsafe&" -> "Y&value" }'
  ))
  source <- serialize(dag, NULL)
  result <- get_identification(dag, quiet = TRUE)
  before <- serialize(result, NULL)
  html <- identification_html(result)

  expect_match(html, "C&lt;b&gt;&amp;", fixed = TRUE)
  expect_match(html, "T&lt;unsafe&amp;", fixed = TRUE)
  expect_match(html, "Y&amp;value", fixed = TRUE)
  expect_false(grepl("<b>|T<unsafe&|Y&value", html))
  expect_false(grepl("&amp;lt;|&amp;amp;", html))
  expect_false(grepl("[Cc]ontrolled intervention", html))
  expect_false(grepl("References:|\\[\\[", html))

  edited <- result
  index <- 0L
  for (effect in names(edited)) {
    for (alternative in seq_along(edited[[effect]])) {
      index <- index + 1L
      edited[[effect]][[alternative]]$notes <-
        paste0("author <edit>& ", index)
    }
  }
  edited_html <- identification_html(edited)
  printed <- identification_print_text(edited)
  for (number in seq_len(index)) {
    expect_match(edited_html, paste0("author &lt;edit&gt;&amp; ", number), fixed = TRUE)
    expect_match(printed, paste0("author <edit>& ", number), fixed = TRUE)
  }
  expect_false(grepl("author <edit>", edited_html, fixed = TRUE))
  expect_identical(identification_html(result, details = TRUE),
                   identification_html(summary(result)))
  expect_identical(serialize(result, NULL), before)
  expect_identical(serialize(dag, NULL), source)
})

test_that("conditional front door reports paired sets separately from the ID verdict", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent]",
    "C -> T C -> M U -> T U -> Y T -> M M -> Y }"
  ))
  result <- get_identification(dag, effect = "total", quiet = TRUE)
  rows <- as.data.frame(result)
  verdict <- attr(result, "verdicts")$total$id
  before <- serialize(result, NULL)

  expect_identical(rows$method, "conditional_front_door")
  expect_identical(rows$id_status, "identifiable")
  expect_identical(rows$variables[[1L]]$mediators, "M")
  expect_identical(rows$variables[[1L]]$covariates, "C")
  expect_identical(verdict$id_status, "identifiable")
  expect_identical(verdict$query, attr(result, "analyses")$total$query)
  expect_identical(attr(rows, "verdicts"), attr(result, "verdicts"))

  selected <- get_identification(dag, method = "conditional_front_door", quiet = TRUE)
  expect_identical(selected$total, result$total)
  expect_identical(attr(selected, "scope"), "total")
  expect_identical(attr(selected, "verdicts"), attr(result, "verdicts"))
  adjustment <- get_identification(
    dag, effect = "total", method = "adjustment", quiet = TRUE
  )
  expect_identical(adjustment$total, list())
  expect_identical(attr(adjustment, "verdicts"), attr(result, "verdicts"))
  expect_identical(attr(adjustment, "analyses"), attr(result, "analyses"))

  printed <- identification_print_text(result)
  expect_match(printed, "Conditional front door adjustment", fixed = TRUE)
  expect_false(grepl("Method: ID|General causal identification", printed))
  expect_identical(serialize(result, NULL), before)
})

test_that("an identified query without a named strategy retains a positive verdict", {
  dag <- dagitty::dagitty(paste(
    "dag { X [exposure] Y [outcome] X -> Z1 Z1 -> Z2 Z1 -> Y",
    "Z2 -> Y X -> Y X <-> Z2 Z1 <-> Y }"
  ))
  result <- get_identification(dag, effect = "total", quiet = TRUE)
  expect_identical(result$total, list())
  expect_equal(nrow(as.data.frame(result)), 0L)
  expect_identical(attr(result, "verdicts")$total$id$id_status, "identifiable")
  expect_identical(attr(as.data.frame(result), "verdicts"), attr(result, "verdicts"))
  printed <- identification_print_text(result)
  expect_match(printed, "ID")
  expect_match(printed, "[Ii]dentifiable")
  expect_false(grepl("[Nn]ot identifiable", printed))
})

test_that("edited strategies do not display or export a stale verdict", {
  dag <- dagitty::dagitty("dag { X [exposure] Y [outcome] X -> Y }")
  result <- get_identification(dag, effect = "total", quiet = TRUE)
  edited <- result
  edited$total[[1L]]$notes <- "author revision"
  expect_match(identification_print_text(edited), "author revision", fixed = TRUE)
  expect_null(attr(as.data.frame(edited), "verdicts"))
  expect_null(attr(data.table::as.data.table(edited), "verdicts"))
  expect_identical(summary(edited)$total, edited$total)
  expect_identical(attr(result, "verdicts")$total$id$id_status, "identifiable")
})
