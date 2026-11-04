test_that("legacy build_graph contract is preserved", {
  implicit <- build_graph("X", "T", "Y")
  explicit <- build_graph(variables = "X", treatments = "T", outcomes = "Y",
                          type = "full")
  modern <- build_graph(confounders = "X", treatments = "T", outcomes = "Y",
                        type = "full")

  expect_equal(edge_keys(implicit), edge_keys(explicit))
  expect_equal(edge_keys(modern), edge_keys(explicit))
  expect_identical(names(formals(build_graph)),
                   c("variables", "treatments", "outcomes", "mediators",
                     "latents", "instrumental_variables",
                     "mediator_outcome_confounders", "competing_causes",
                     "colliders", "type", "confounders", "latent_variables"))
  expect_identical(formals(build_graph)$type, "full")
})

test_that("build_graph() accepts the deprecated latent_variables alias", {
  canonical <- build_graph("X", "T", "Y", mediators = "M", latents = "U",
                           type = "ordered")
  expect_warning(
    aliased <- build_graph("X", "T", "Y", mediators = "M",
                           latent_variables = "U", type = "ordered"),
    paste0("build_graph(): 'latent_variables' is deprecated as of 0.0.2. ",
           "Use 'latents' instead."),
    fixed = TRUE
  )

  expect_identical(edge_keys(aliased), edge_keys(canonical))
  expect_identical(dagitty::latents(aliased), dagitty::latents(canonical))
  expect_identical(dagitty::latents(canonical), "U")

  expect_error(
    build_graph("X", "T", "Y", latents = "U", latent_variables = "U"),
    paste0("build_graph(): supply either 'latents' or its deprecated alias ",
           "'latent_variables', not both."),
    fixed = TRUE
  )

  expect_warning(
    sentinel <- build_graph("X", "T", "Y", mediators = "M", latents = NA,
                            latent_variables = "U", type = "ordered"),
    paste0("build_graph(): 'latent_variables' is deprecated as of 0.0.2. ",
           "Use 'latents' instead."),
    fixed = TRUE
  )
  expect_identical(edge_keys(sentinel), edge_keys(canonical))
  expect_identical(dagitty::latents(sentinel), "U")

  positional <- build_graph("X", "T", "Y", "M", "U", type = "ordered")
  expect_identical(edge_keys(positional), edge_keys(canonical))
})

test_that("standard graph equations retain their edge contract", {
  dag <- build_graph(
    variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
    mediators = "M", instrumental_variables = "IV", type = "saturated"
  )

  expect_identical(
    edge_keys(dag),
    sort(c("IV -> T", "M -> Y", "T -> M", "T -> Y", "X1 -> M",
           "X1 -> T", "X1 -> Y", "X1 -> X2", "X2 -> M", "X2 -> T",
           "X2 -> Y"))
  )
  expect_identical(get_roles(dag)$instrument, "IV")
  expect_identical(get_roles(dag)$mediator, "M")
  expect_identical(get_roles(dag)$confounder, c("X1", "X2"))
})

test_that("competing_causes is the single competing vocabulary", {
  dag <- build_graph(
    variables = "X", treatments = "T", outcomes = "Y",
    competing_causes = "C", type = "saturated"
  )

  expect_identical(competing_causes(dag), "C")
  expect_warning(
    aliased <- competing_exposures(dag),
    paste0("competing_exposures(): 'competing_exposures()' is deprecated as of ",
           "0.0.1. Use 'competing_causes()' instead."),
    fixed = TRUE, class = "deprecatedWarning"
  )
  expect_identical(aliased, "C")
  expect_identical(
    names(get_roles(dag))[1:10],
    c("outcome", "treatment", "confounder", "mediator",
      "mediator_outcome_confounder", "instrument", "competing_cause",
      "collider", "latent", "observed")
  )
  expect_false("competing_exposure" %in% names(get_roles(dag)))
  expect_identical(tail(names(get_roles(dag)), 1L), "undetermined")
  expect_error(get_edges(dag, "competing_exposure"), "not valid roles")
  expect_true("competing_cause" %in%
                c(get_edges(dag)$ancestor_role,
                  get_edges(dag)$descendant_role))
  expect_false("competing_exposure" %in%
                 c(get_edges(dag)$ancestor_role,
                   get_edges(dag)$descendant_role))
  expect_identical(get_nodes(dag)$C, "competing_cause")
  structure_roles <- unlist(lapply(get_structure(dag), function(node){
    node$role
  }), use.names = FALSE)
  expect_true("competing_cause" %in% structure_roles)
  expect_false("competing_exposure" %in% structure_roles)

  without_competing <- build_graph("X", "T", "Y", type = "saturated")
  role_diff <- get_diff_roles(dag, without_competing)
  expect_identical(role_diff$competing_cause, "C")
  expect_false("competing_exposure" %in% names(role_diff))
  expect_s3_class(
    plot_dagitty(dag, include_legend = "competing_cause"),
    "ggplot"
  )
  expect_error(
    build_graph("X", "T", "Y", competing_exposures = "C1"),
    "unused argument"
  )
})

test_that("construction preserves special node names exactly", {
  cases <- list(
    c("soil moisture", "air temp"),
    c("soil-moisture", "air-temp"),
    c("soil%", "air%")
  )

  for (variables in cases) {
    dag <- build_graph(variables, "T", "Y", type = "saturated")
    expect_setequal(names(dag), c(variables, "T", "Y"))
  }
})

test_that("graph-building entry points keep their established defaults", {
  dag <- build_graph(c("X1", "X2"), "T", "Y", type = "ordered")

  set.seed(1)
  added <- suppressMessages(
    add_nodes(dag, "W", role = "confounder", print_edges = FALSE)
  )
  connected <- suppressMessages(
    connect_nodes(dag, nodes = "X1", print_edges = FALSE)
  )

  expect_true("W" %in% names(added))
  expect_setequal(names(connected), names(dag))
  expect_identical(formals(add_nodes)$type, "saturated")
  expect_identical(formals(connect_nodes)$type, "full")
  expect_error(build_graph("X", "T", "Y", type = "invalid"),
               "Invalid type")
})

test_that("position additions have basic structural coverage", {
  dag <- build_graph(c("X1", "X2"), "T", "Y", type = "ordered")

  set.seed(1)
  positioned <- suppressMessages(
    add_nodes(dag, "W", role = "confounder", type = "ordered",
              position = first(), print_edges = FALSE)
  )
  expect_true("W -> X1" %in% edge_keys(positioned))

})

test_that("rebuilding an existing graph keeps latent nodes without children", {
  childless <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] U [latent] T -> Y }"
  )
  rebuilt <- build_graph(childless, type = "ordered")

  expect_setequal(dagitty::latents(rebuilt), "U")
  expect_setequal(edge_keys(rebuilt), c("U -> T", "U -> Y", "T -> Y"))

  mixed <- dagitty::dagitty(
    paste("dag { T [exposure] Y [outcome] X U1 [latent] U2 [latent]",
          "U1 -> T U1 -> Y X -> T X -> Y T -> Y }")
  )
  remixed <- build_graph(mixed, type = "ordered")

  expect_setequal(dagitty::latents(remixed), c("U1", "U2"))
  expect_setequal(edge_keys(remixed),
                  c("U1 -> T", "U1 -> Y", "U2 -> T", "U2 -> Y", "T -> Y",
                    "X -> T", "X -> Y"))
})


test_that("build_graph() rejects absent treatments or outcomes with its own message", {
  message <- "The 'treatments' and 'outcomes' inputs should be used"
  expect_error(build_graph(variables = "X", treatments = NULL, outcomes = "Y"),
               message, fixed = TRUE)
  expect_error(build_graph(variables = "X", treatments = "T", outcomes = character(0)),
               message, fixed = TRUE)
  expect_error(build_graph(variables = "X", treatments = list(), outcomes = "Y"),
               message, fixed = TRUE)
})


test_that("graph type inputs must be a single string", {
  expect_error(
    build_graph(treatments = "T", outcomes = "Y", type = NULL),
    "Invalid type: \"\". Must be one of: full, chain, first, last, ordered, saturated.",
    fixed = TRUE
  )
  expect_error(
    build_graph(treatments = "T", outcomes = "Y", type = c("full", "chain")),
    "Invalid type: \"full, chain\". Must be one of: full, chain, first, last, ordered, saturated.",
    fixed = TRUE
  )
  dag <- build_graph(variables = "X", treatments = "T", outcomes = "Y",
                     type = "saturated")
  expect_error(
    add_nodes(dag, nodes = "W", role = "confounder", type = character(0)),
    "Invalid type: \"\". Must be one of: full, chain, ordered, saturated.",
    fixed = TRUE
  )
  expect_error(
    connect_nodes(dag, type = c("full", "chain")),
    "Invalid type: \"full, chain\". Must be one of: full, chain, ordered, saturated.",
    fixed = TRUE
  )
})
