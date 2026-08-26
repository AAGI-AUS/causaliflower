test_that("legacy build_graph contract is preserved", {
  implicit <- build_graph("Z", "X", "Y")
  explicit <- build_graph(variables = "Z", treatments = "X", outcomes = "Y",
                          type = "full")
  modern <- build_graph(confounders = "Z", treatments = "X", outcomes = "Y",
                        type = "full")

  expect_equal(edge_keys(implicit), edge_keys(explicit))
  expect_equal(edge_keys(modern), edge_keys(explicit))
  expect_identical(names(formals(build_graph))[1:10],
                   c("variables", "treatments", "outcomes", "mediators",
                     "latent_variables", "instrumental_variables",
                     "mediator_outcome_confounders", "competing_causes",
                     "colliders", "type"))
  expect_identical(formals(build_graph)$type, "full")
})

test_that("standard graph equations retain their edge contract", {
  dag <- build_graph(
    variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
    mediators = "M", instrumental_variables = "IV", type = "saturated"
  )

  expect_identical(
    edge_keys(dag),
    sort(c("IV -> X", "M -> Y", "X -> M", "X -> Y", "Z1 -> X",
           "Z1 -> Y", "Z1 -> Z2", "Z2 -> X", "Z2 -> Y"))
  )
  expect_identical(get_roles(dag)$instrument, "IV")
  expect_identical(get_roles(dag)$mediator, "M")
  expect_identical(get_roles(dag)$confounder, c("Z1", "Z2"))
})

test_that("competing_causes is the single competing vocabulary", {
  dag <- build_graph(
    variables = "Z", treatments = "X", outcomes = "Y",
    competing_causes = "C", type = "saturated"
  )

  expect_identical(competing_causes(dag), "C")
  expect_error(competing_exposures(dag), "renamed competing_causes")
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
    c(node$ancestor$role, node$descendant$role)
  }), use.names = FALSE)
  expect_true("competing_cause" %in% structure_roles)
  expect_false("competing_exposure" %in% structure_roles)

  without_competing <- build_graph("Z", "X", "Y", type = "saturated")
  role_diff <- get_diff_roles(dag, without_competing)
  expect_identical(role_diff$competing_cause, "C")
  expect_false("competing_exposure" %in% names(role_diff))
  expect_s3_class(
    plot_dagitty(dag, include_legend = "competing_cause"),
    "ggplot"
  )
  expect_error(
    build_graph("Z", "X", "Y", competing_exposures = "C1"),
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
    dag <- build_graph(variables, "X", "Y", type = "saturated")
    expect_setequal(names(dag), c(variables, "X", "Y"))
  }
})

test_that("graph-building entry points keep their established defaults", {
  dag <- build_graph(c("Z1", "Z2"), "X", "Y", type = "ordered")

  set.seed(1)
  added <- suppressMessages(
    add_nodes(dag, "W", node_role = "confounder", print_edges = FALSE)
  )
  connected <- suppressMessages(
    connect_nodes(dag, nodes = "Z1", print_edges = FALSE)
  )

  expect_true("W" %in% names(added))
  expect_setequal(names(connected), names(dag))
  expect_identical(formals(add_nodes)$type, "saturated")
  expect_identical(formals(connect_nodes)$type, "full")
  expect_error(build_graph("Z", "X", "Y", type = "invalid"),
               "Invalid type")
})

test_that("position additions have basic structural coverage", {
  dag <- build_graph(c("Z1", "Z2"), "X", "Y", type = "ordered")

  set.seed(1)
  positioned <- suppressMessages(
    add_nodes(dag, "W", node_role = "confounder", type = "ordered",
              position = first(), print_edges = FALSE)
  )
  expect_true("W -> Z1" %in% edge_keys(positioned))

})
