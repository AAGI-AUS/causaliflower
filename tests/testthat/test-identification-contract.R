test_that("successful results use one stable contract", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] X X -> T X -> Y T -> Y }"
  )

  result <- identification_result(dag)

  expect_identification_shape(result)
  expect_identical(attr(result, "treatment"), "T")
  expect_identical(attr(result, "outcome"), "Y")
  expect_identical(attr(result, "effect"), "total")
  expect_false(attr(result, "conditional"))
  expect_identical(attr(result, "conditioning"), character())
  expect_identical(attr(result, "conditioning_intervened"), character())
  expect_identical(attr(result, "conditioning_remaining"), character())
  expect_true(result$feasible[1L])

  directed_latent_chain <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] L1 [latent] L2 [latent] ",
      "T -> L1 L1 -> L2 L2 -> Y }"
    )
  )
  chain_result <- identification_result(directed_latent_chain)
  expect_identification_shape(chain_result)
  expect_true(attr(chain_result, "observationally_identifiable"))
  expect_true(attr(chain_result, "identifiable"))
  expect_null(attr(chain_result, "hedge"))

  projected_bow <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] U [latent] L [latent] ",
      "U -> T U -> Y T -> L L -> Y }"
    )
  )
  bow_result <- identification_result(projected_bow)
  expect_identification_shape(bow_result)
  expect_false(attr(bow_result, "observationally_identifiable"))
  expect_false(attr(bow_result, "identifiable"))
  expect_type(attr(bow_result, "hedge"), "list")
  expect_setequal(
    identification_edge_keys(attr(bow_result, "hedge")$F),
    c("T -> Y", "T <-> Y")
  )
})

test_that("strategy and role tables use the documented schema", {
  dag <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] X M ",
      "X -> T X -> Y T -> M M -> Y T -> Y }"
    )
  )

  result <- identification_result(dag)
  roles <- attr(result, "roles")

  expect_identical(
    unname(unlist(result[1L], use.names = FALSE)),
    c(
      "general identification", "TRUE", "nonparametric",
      "not applicable", "the ID algorithm identifies the total effect"
    )
  )

  adjustment <- identification_strategy(result, "^covariate adjustment$")
  expect_identical(adjustment$strategy, "covariate adjustment")
  expect_true(adjustment$feasible)
  expect_identical(adjustment$tier, "nonparametric")
  expect_identical(adjustment$variable_set, "{X}")
  expect_identical(adjustment$note,
                   "minimal adjustment sets for the total effect")

  expect_identical(
    names(roles),
    c("variable", "role", "conditioning_guidance")
  )
  confounder <- roles[variable == "X" & role == "adjustment variable"]
  expect_equal(nrow(confounder), 1L)
  expect_identical(
    confounder$conditioning_guidance,
    "member of at least one minimal adjustment set for the total effect; use a complete set"
  )
  mediator <- roles[variable == "M" & role == "mediator"]
  expect_equal(nrow(mediator), 1L)
  expect_identical(
    mediator$conditioning_guidance,
    "do not include in an adjustment set for the total effect"
  )

  direct <- identification_result(dag, effect = "direct")
  direct_roles <- attr(direct, "roles")
  outcome_parent <- direct_roles[
    variable == "M" & role == "outcome parent"
  ]
  expect_equal(nrow(outcome_parent), 1L)
  expect_identical(
    outcome_parent$conditioning_guidance,
    "held fixed in the direct effect query"
  )

  expect_identical(attr(result, "data_availability"),
                   "observational data")
  expect_identical(
    attr(result, "caveat"),
    paste(
      "Identification is a graphical result. Estimation also requires",
      "positivity and data that support the comparisons used by the",
      "reported strategy."
    )
  )
})

test_that("rendered sets distinguish node names that contain delimiters", {
  expect_identical(id_render_set("A, B"), "{\"A, B\"}")
  expect_identical(id_render_set(c("A", "B")), "{A, B}")
})

test_that("explicit untagged endpoints are accepted", {
  untagged <- dagitty::dagitty("dag { A -> B }")
  result <- identification_result(untagged, treatment = "A", outcome = "B")

  expect_true(attr(result, "identifiable"))
  expect_identical(attr(result, "treatment"), "A")
  expect_identical(attr(result, "outcome"), "B")

  multi <- dagitty::dagitty(
    "dag { T1 [exposure] T2 [exposure] Y1 [outcome] Y2 [outcome]
           T1 -> Y1 T2 -> Y2 }"
  )
  expect_error(get_identification(multi, quiet = TRUE))
  selected <- identification_result(
    multi, treatment = "T1", outcome = "Y1"
  )
  expect_true(attr(selected, "identifiable"))

  no_treatment <-
    "No treatment detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again."
  multiple_treatments <-
    "Multiple treatments detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again."
  invalid_treatment <-
    "Invalid treatment detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again."
  no_outcome <-
    "No outcome detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again."
  multiple_outcomes <-
    "Multiple outcomes detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again."
  invalid_outcome <-
    "Invalid outcome detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again."
  distinct_endpoints <-
    "Identical treatment and outcome detected: 'get_identification()' requires the named treatment and outcome to refer to different nodes in the DAG. Please ensure different treatment and outcome nodes are supplied and try again."

  expect_identification_error(
    get_identification(untagged, quiet = TRUE), no_treatment
  )
  expect_identification_error(
    get_identification(untagged, treatment = "A", quiet = TRUE), no_outcome
  )
  expect_identification_error(
    get_identification(untagged, treatment = c("A", "B"), outcome = "B",
                       quiet = TRUE),
    multiple_treatments
  )
  expect_identification_error(
    get_identification(untagged, treatment = "A", outcome = c("A", "B"),
                       quiet = TRUE),
    multiple_outcomes
  )
  expect_identification_error(
    get_identification(untagged, treatment = "A", outcome = "A",
                       quiet = TRUE),
    distinct_endpoints
  )
  expect_identification_error(
    get_identification(untagged, treatment = "missing", outcome = "B",
                       quiet = TRUE),
    invalid_treatment
  )
  expect_identification_error(
    get_identification(untagged, treatment = 1, outcome = "B",
                       quiet = TRUE),
    invalid_treatment
  )
  expect_identification_error(
    get_identification(untagged, treatment = "A", outcome = 1,
                       quiet = TRUE),
    invalid_outcome
  )

  latent <- dagitty::dagitty(
    "dag { L [latent] T Y L -> T L -> Y T -> Y }"
  )
  expect_identification_error(
    get_identification(latent, treatment = "L", outcome = "Y", quiet = TRUE),
    invalid_treatment
  )
})

test_that("get_identification validates the plain-DAG boundary", {
  dag <- dagitty::dagitty("dag { T [exposure] Y [outcome] T -> Y }")
  attr(dag, "full_dag_text") <-
    "dag { T [exposure] Y [outcome] T -> Y T <-> Y }"
  attr(dag, "full_treatments") <- "not-T"
  attr(dag, "full_outcomes") <- "not-Y"

  result <- identification_result(dag)

  expect_true(attr(result, "identifiable"))
  expect_identical(attr(result, "treatment"), "T")
  expect_identical(attr(result, "outcome"), "Y")

  query_carrier <- structure(list(), class = "dag_by_query")
  attr(query_carrier, "dag_text") <- as.character(dag)
  expect_error(get_identification(query_carrier, treatment = "T",
                                  outcome = "Y", quiet = TRUE))

  invalid_input <-
    "Invalid input detected: only dagitty objects may be supplied as 'dag' input. Please check the input and try again."
  cycle_error <-
    "Cycles detected: only acyclic dagitty objects are supported as 'dag' input. Please check the input and try again."
  undirected_error <-
    "Unsupported edge type detected: only directed or bidirected edges may be supplied as 'dag' input. Please supply a dagitty object containing only directed and/or bidirected edges."
  graph_type_error <-
    "Unsupported graph type detected: 'get_identification()' requires a 'dagitty' DAG. Please supply one and try again."
  status_error <-
    "Adjusted/selected node(s) detected: these are not supported by 'get_identification()' as 'dag' input. Please remove these from the supplied dagitty object and try again."
  latent_bidirected_error <-
    "Unsupported bidirected edge detected: a bidirected edge in the DAG connects to a latent node. 'get_identification()' currently supports bidirected edges only between observed nodes. Please revise the supplied DAG and try again."
  cyclic_boundary <- structure(
    "dag { T [exposure]; Y [outcome]; T -> Y; Y -> T }",
    class = "dagitty"
  )
  undirected_boundary <- structure(
    "dag { T [exposure]; Y [outcome]; T -- Y }",
    class = "dagitty"
  )

  expect_identification_error(
    get_identification("not a DAG", quiet = TRUE), invalid_input
  )
  expect_identification_error(
    get_identification(cyclic_boundary, quiet = TRUE),
    cycle_error
  )
  expect_identification_error(
    get_identification(undirected_boundary, quiet = TRUE),
    undirected_error
  )

  adjusted <- dagitty::dagitty(
    "dag { A [adjusted] T [exposure] Y [outcome] A -> T T -> Y }"
  )
  selected <- dagitty::dagitty(
    "dag { S [selected] T [exposure] Y [outcome] T -> S <- Y }"
  )
  for (dag in list(adjusted, selected)) {
    expect_identification_error(
      get_identification(dag, quiet = TRUE), status_error
    )
  }

  latent_bidirected <- dagitty::dagitty(
    "dag { L [latent] T [exposure] Y [outcome] L <-> T T -> Y }"
  )
  expect_identification_error(
    get_identification(latent_bidirected, quiet = TRUE),
    latent_bidirected_error
  )

  pdag <- dagitty::dagitty(
    "pdag { T [exposure]; Y [outcome]; T -- Y }"
  )
  expect_identification_error(
    get_identification(pdag, quiet = TRUE), undirected_error
  )

  mag <- dagitty::dagitty(
    "mag { T [exposure]; Y [outcome]; T -> Y }"
  )
  expect_identification_error(
    get_identification(mag, quiet = TRUE), graph_type_error
  )
})

test_that("unsupported combinations fail instead of being ignored", {
  dag <- dagitty::dagitty(
    "dag { L [latent] T [exposure] Y [outcome] Z W
           L -> T L -> Y T -> Z Z -> Y W -> Y }"
  )

  bad_conditioning <-
    "Invalid conditioning input detected: 'conditioning' requires observed DAG node names as input. Please ensure each supplied name corresponds to an observed node in the supplied DAG other than the treatment and outcome, and try again."
  bad_experiments <-
    "Invalid experiments input detected: 'experiments' requires observed DAG node names as input. Please ensure each supplied name corresponds to an observed node in the supplied DAG other than the treatment and outcome, and try again."
  conditioning_total_only <-
    "Unsupported input combination detected: 'conditioning' must be 'NULL' when 'effect = \"direct\"'. Please remove the supplied 'conditioning' input or use 'effect = \"total\"' and try again."
  experiments_total_only <-
    "Unsupported input combination detected: 'experiments' must be 'NULL' when 'effect = \"direct\"'. Please remove the supplied 'experiments' input or use 'effect = \"total\"' and try again."
  incompatible_data <-
    "Unsupported input combination detected: 'conditioning' and 'experiments' cannot both be supplied when 'effect = \"total\"'. Please remove either the supplied 'conditioning' input or the supplied 'experiments' input and try again."
  invalid_effect <-
    "Invalid input detected: only 'total' or 'direct' may be supplied as 'effect' input. Please check the input and try again."
  invalid_quiet <-
    "Invalid input detected: only TRUE or FALSE may be supplied as 'quiet' input."

  expect_identification_error(
    get_identification(dag, conditioning = "missing", quiet = TRUE),
    bad_conditioning
  )
  expect_identification_error(
    get_identification(dag, experiments = "missing", quiet = TRUE),
    bad_experiments
  )
  for (bad in list(NA_character_, "", "L", "T", "Y")) {
    expect_identification_error(
      get_identification(dag, conditioning = bad, quiet = TRUE),
      bad_conditioning
    )
    expect_identification_error(
      get_identification(dag, experiments = bad, quiet = TRUE),
      bad_experiments
    )
  }

  expect_identification_error(
    get_identification(dag, conditioning = "Z", experiments = "W",
                       quiet = TRUE),
    incompatible_data
  )
  expect_identification_error(
    get_identification(dag, effect = "direct", conditioning = "Z",
                       quiet = TRUE),
    conditioning_total_only
  )
  expect_identification_error(
    get_identification(dag, effect = "direct", experiments = "W",
                       quiet = TRUE),
    experiments_total_only
  )
  expect_identification_error(
    get_identification(dag, effect = "unknown", quiet = TRUE),
    invalid_effect
  )

  for (bad_quiet in list(NA, c(TRUE, FALSE), 1, "yes")) {
    expect_identification_error(
      get_identification(dag, quiet = bad_quiet), invalid_quiet
    )
  }
})

test_that("quiet suppresses the summary message", {
  expect_summary <- function(code, expected) {
    messages <- testthat::capture_messages(code)
    expect_identical(messages, paste0(expected, "\n"))
  }

  total_identified <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> Y }"
  )
  expect_summary(
    get_identification(total_identified),
    paste0(
      "The total effect of 'T' on 'Y' is identifiable from ",
      "observational data under this DAG."
    )
  )

  z_rescue <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z Z -> T T -> Y T <-> Z Z <-> Y }"
  )
  expect_summary(
    get_identification(z_rescue, experiments = "Z"),
    paste0(
      "The total effect of 'T' on 'Y' is not identifiable from ",
      "observational data alone, but is identifiable using surrogate ",
      "experimental data for {Z} declared through 'experiments'."
    )
  )

  bow <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> Y T <-> Y }"
  )
  expect_summary(
    get_identification(bow),
    paste0(
      "The total effect of 'T' on 'Y' is not identifiable from the ",
      "supplied data under this DAG."
    )
  )

  condition_helps <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z T -> Z Z -> Y T <-> Z }"
  )
  expect_summary(
    get_identification(condition_helps, conditioning = "Z"),
    paste0(
      "The effect of 'T' on 'Y' given {Z} is identifiable from ",
      "observational data under this DAG."
    )
  )

  condition_hurts <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z T -> Z Y -> Z T <-> Z }"
  )
  expect_summary(
    get_identification(condition_hurts, conditioning = "Z"),
    paste0(
      "The effect of 'T' on 'Y' given {Z} is not identifiable from ",
      "observational data under this DAG."
    )
  )

  direct_identified <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] M T -> M M -> Y T -> Y }"
  )
  expect_summary(
    get_identification(direct_identified, effect = "direct"),
    paste0(
      "The outcome distribution for the controlled direct effect of 'T' ",
      "on 'Y' is identifiable from observational data in the ",
      "latent-projected graph (with the outcome's other parents held fixed)."
    )
  )
  expect_summary(
    get_identification(bow, effect = "direct"),
    paste0(
      "The outcome distribution for the controlled direct effect of 'T' ",
      "on 'Y' is not identifiable from observational data in the ",
      "latent-projected graph (with the outcome's other parents held fixed)."
    )
  )

  no_total_path <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] }"
  )
  expect_summary(
    get_identification(no_total_path),
    paste0(
      "The DAG contains no directed path from 'T' to 'Y', so it encodes ",
      "no total causal effect of 'T' on 'Y'."
    )
  )

  no_direct_arrow <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] M T -> M M -> Y }"
  )
  expect_summary(
    get_identification(no_direct_arrow, effect = "direct"),
    paste0(
      "The latent-projected graph contains no direct arrow from 'T' to ",
      "'Y'. The verdict attributes refer to the intervention on the ",
      "treatment and fixing set."
    )
  )

  expect_identical(
    testthat::capture_messages(
      get_identification(total_identified, quiet = TRUE)
    ),
    character()
  )
})

test_that("names, isolates, and edge order do not change the verdict", {
  isolated <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] W }"
  )
  isolated_result <- identification_result(isolated)

  expect_identification_shape(isolated_result)
  expect_true(attr(isolated_result, "identifiable"))
  expect_true(attr(isolated_result, "observationally_identifiable"))
  expect_true(isolated_result$feasible[1L])

  named <- dagitty::dagitty(
    paste0(
      "dag { \"soil water\" [exposure] \"grain-yield\" [outcome] ",
      "\"field block\" ",
      "\"soil water\" -> \"field block\" ",
      "\"field block\" -> \"grain-yield\" ",
      "\"soil water\" <-> \"field block\" }"
    )
  )
  named_result <- identification_result(named, conditioning = "field block")

  expect_identification_shape(named_result)
  expect_true(attr(named_result, "identifiable"))
  expect_identical(attr(named_result, "treatment"), "soil water")
  expect_identical(attr(named_result, "outcome"), "grain-yield")
  expect_identical(attr(named_result, "conditioning_intervened"),
                   "field block")

  first <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z T -> Z Z -> Y T <-> Z }"
  )
  second <- dagitty::dagitty(
    "dag { Y [outcome] Z T [exposure] T <-> Z Z -> Y T -> Z }"
  )

  first_result <- identification_result(first, conditioning = "Z")
  second_result <- identification_result(second, conditioning = "Z")

  expect_identical(first_result, second_result)

  relabelled <- dagitty::dagitty(
    "dag { A [exposure] B [outcome] C A -> C C -> B A <-> C }"
  )
  relabelled_result <- identification_result(
    relabelled, conditioning = "C"
  )

  expect_identical(attr(first_result, "identifiable"),
                   attr(relabelled_result, "identifiable"))
  expect_identical(attr(first_result, "observationally_identifiable"),
                   attr(relabelled_result, "observationally_identifiable"))
  expect_identical(first_result$feasible, relabelled_result$feasible)
  expect_identical(first_result$tier, relabelled_result$tier)
})


test_that("named scalar inputs preserve the requested query", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] M T -> M M -> Y T -> Y }"
  )
  for (effect in c("total", "direct")) {
    plain <- identification_result(dag, effect = effect)
    named <- identification_result(
      dag, treatment = c(input = "T"), outcome = c(input = "Y"),
      effect = stats::setNames(effect, "query")
    )
    expect_identical(named, plain)
    expect_identical(
      testthat::capture_messages(get_identification(
        dag, effect = stats::setNames(effect, "query")
      )),
      testthat::capture_messages(get_identification(dag, effect = effect))
    )
  }
  expect_error(
    get_identification(dag, treatment = c(input = "T"), outcome = "T"),
    "Identical treatment and outcome detected", fixed = TRUE
  )
  expect_error(
    get_identification(dag, treatment = "T", outcome = c(input = "T")),
    "Identical treatment and outcome detected", fixed = TRUE
  )
  expect_error(
    get_identification(
      dag, treatment = c(left = "T"), outcome = c(right = "T")
    ),
    "Identical treatment and outcome detected", fixed = TRUE
  )
})


test_that("only NULL omits an incompatible optional argument", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] M T -> M M -> Y T -> Y }"
  )
  expect_error(
    get_identification(dag, effect = "direct", conditioning = character()),
    "Unsupported input combination detected", fixed = TRUE
  )
  expect_error(
    get_identification(dag, effect = "direct", experiments = character()),
    "Unsupported input combination detected", fixed = TRUE
  )
  expect_error(
    get_identification(dag, conditioning = character(), experiments = "M"),
    "Unsupported input combination detected", fixed = TRUE
  )
  expect_error(
    get_identification(dag, conditioning = "M", experiments = character()),
    "Unsupported input combination detected", fixed = TRUE
  )
  expect_error(
    get_identification(
      dag, conditioning = character(), experiments = character()
    ),
    "Unsupported input combination detected", fixed = TRUE
  )
  for (value in list(numeric(), logical(), list())) {
    expect_error(
      get_identification(dag, conditioning = value),
      "Invalid conditioning input detected", fixed = TRUE
    )
    expect_error(
      get_identification(dag, experiments = value),
      "Invalid experiments input detected", fixed = TRUE
    )
  }
  expect_identical(
    identification_result(dag, conditioning = character()),
    identification_result(dag)
  )
  expect_identical(
    identification_result(dag, experiments = character()),
    identification_result(dag)
  )
})


test_that("identification keeps node names with backslashes and double quotes", {
  backslash <- dagitty::dagitty('dag { "a\\\\b" [exposure] Y [outcome] "a\\\\b" -> Y }')
  expect_s3_class(identification_result(backslash), "data.frame")

  quoted <- dagitty::dagitty('dag { "a\\"b" [exposure] Y [outcome] "a\\"b" -> Y }')
  expect_s3_class(identification_result(quoted), "data.frame")
})
