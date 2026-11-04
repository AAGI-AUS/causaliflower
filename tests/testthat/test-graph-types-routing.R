role_edge <- function(dag, from, to, edge = "->") {
  keys <- edge_keys(dag)
  paste(from, edge, to) %in% keys
}


symmetric_role_edge <- function(dag, left, right) {
  keys <- edge_keys(dag)
  paste(left, "<->", right) %in% keys ||
    paste(right, "<->", left) %in% keys
}


test_that("chain and ordered respect temporal positions", {
  input <- list("M1", c("M2", "M3"), "M4")
  chain <- build_graph("X", "T", "Y", mediators = input, type = "chain")
  ordered <- build_graph("X", "T", "Y", mediators = input, type = "ordered")

  adjacent <- c("M1 -> M2", "M1 -> M3", "M2 -> M4", "M3 -> M4")
  expect_true(all(adjacent %in% edge_keys(chain)))
  expect_false("M1 -> M4" %in% edge_keys(chain))
  expect_false(any(c("M2 -> M3", "M3 -> M2") %in% edge_keys(chain)))

  expect_true(all(c(adjacent, "M1 -> M4") %in% edge_keys(ordered)))
  expect_false(any(c("M2 -> M3", "M3 -> M2") %in% edge_keys(ordered)))
})


test_that("chain adds no arrows between roles", {
  dag <- build_graph(
    c("X1", "X2"), "T", "Y", mediators = c("M1", "M2"),
    mediator_outcome_confounders = c("K1", "K2"),
    competing_causes = c("P1", "P2"), colliders = c("L1", "L2"),
    latents = c("U1", "U2"), type = "chain"
  )

  expect_setequal(edge_keys(dag), c(
    "X1 -> X2", "M1 -> M2", "K1 -> K2", "P1 -> P2", "L1 -> L2", "U1 -> U2"
  ))
})


test_that("chain keeps mapped latent and instrumental edges", {
  dag <- build_graph(
    "X", "T", "Y", mediators = c("M1", "M2"),
    latents = list(list("U", c("M1", "M2"))),
    instrumental_variables = "IV", type = "chain"
  )

  expect_setequal(edge_keys(dag), c(
    "M1 -> M2", "U -> M1", "U -> M2", "IV -> T"
  ))
})


test_that("first and last differ only in the sending position", {
  roles <- list(
    variables = c("X1", "X2", "X3"), treatments = "T", outcomes = "Y",
    mediators = c("M1", "M2", "M3")
  )
  within <- c("X1 -> X2", "X2 -> X3", "M1 -> M2", "M2 -> M3")

  from_last <- do.call(build_graph, c(roles, list(type = "last")))
  from_first <- do.call(build_graph, c(roles, list(type = "first")))

  expect_setequal(edge_keys(from_last), c(
    within, "X3 -> T", "X3 -> Y", "T -> M1", "T -> Y", "M3 -> Y"
  ))
  expect_setequal(edge_keys(from_first), c(
    within, "X1 -> T", "X1 -> Y", "T -> M1", "T -> Y", "M1 -> Y"
  ))
})


test_that("first and last keep direct arrows for source and collider roles", {
  for (type in c("first", "last")) {
    dag <- build_graph(
      "X", "T", "Y", mediators = "M",
      mediator_outcome_confounders = c("K1", "K2"),
      competing_causes = c("P1", "P2"), colliders = c("L1", "L2"),
      type = type
    )
    keys <- edge_keys(dag)

    expect_true(all(c(
      "K1 -> M", "K2 -> M", "K1 -> Y", "K2 -> Y",
      "P1 -> Y", "P2 -> Y",
      "T -> L1", "T -> L2", "Y -> L1", "Y -> L2"
    ) %in% keys), info = type)
  }
})


test_that("routing markers override the first sending position", {
  dag <- build_graph(
    "X", "T", "Y", mediators = c("M1", "M2", outgoing("M3")),
    type = "first"
  )
  keys <- edge_keys(dag)

  expect_true(all(c("M3 -> Y", "T -> M1") %in% keys))
  expect_false(any(c("M1 -> Y", "M2 -> Y") %in% keys))
})


test_that("build_graph() accepts first and last as graph types", {
  for (type in c("full", "chain", "first", "last", "ordered", "saturated")) {
    expect_s3_class(
      build_graph("X", "T", "Y", mediators = "M", type = type),
      "dagitty", exact = FALSE
    )
  }

  expect_error(build_graph("X", "T", "Y", type = "middle"), "Invalid type")
})


test_that("last connects the ends of successive role sequences", {
  dag <- build_graph(
    c("X1", "X2", "X3", "X4"), "T", "Y",
    mediators = c("M1", "M2", "M3", "M4"), type = "last"
  )

  expect_setequal(edge_keys(dag), c(
    "X1 -> X2", "X2 -> X3", "X3 -> X4",
    "X4 -> T", "X4 -> Y",
    "T -> M1", "T -> Y",
    "M1 -> M2", "M2 -> M3", "M3 -> M4", "M4 -> Y"
  ))
})


test_that("last includes every node in a parallel boundary position", {
  arrows <- function(from, to) {
    as.vector(outer(from, to, function(v, w) paste(v, "->", w)))
  }
  dag <- build_graph(
    variables = list(c("X1", "X2"), c("X3", "X4")),
    treatments = list(c("T1", "T2"), c("T3", "T4")),
    outcomes = list(c("Y1", "Y2"), c("Y3", "Y4")),
    mediators = list(c("M1", "M2"), c("M3", "M4")),
    type = "last"
  )

  expected <- c(
    arrows(c("X1", "X2"), c("X3", "X4")),
    arrows(c("T1", "T2"), c("T3", "T4")),
    arrows(c("M1", "M2"), c("M3", "M4")),
    arrows(c("Y1", "Y2"), c("Y3", "Y4")),
    arrows(c("X3", "X4"), c("T1", "T2")),
    arrows(c("X3", "X4"), c("Y1", "Y2")),
    arrows(c("T3", "T4"), c("M1", "M2")),
    arrows(c("T3", "T4"), c("Y1", "Y2")),
    arrows(c("M3", "M4"), c("Y1", "Y2"))
  )
  expect_setequal(edge_keys(dag), expected)
})


test_that("routing annotations override last boundary defaults", {
  dag <- build_graph(
    variables = c("X1", outgoing("X2"), "X3"),
    treatments = "T", outcomes = "Y",
    mediators = c(outgoing("M1"), incoming("M2"), "M3"),
    mediator_outcome_confounders = c("K1", outgoing("K2"), "K3"),
    competing_causes = c("C1", outgoing("C2"), "C3"),
    colliders = c("L1", incoming("L2"), "L3"),
    type = "last"
  )
  keys <- edge_keys(dag)

  expect_true(all(c(
    "X2 -> T", "X2 -> Y", "T -> M2", "M1 -> Y",
    "K2 -> M2", "K2 -> Y", "C2 -> Y", "T -> L2", "Y -> L2"
  ) %in% keys))
  expect_false(any(c(
    "X3 -> T", "X3 -> Y", "T -> M1", "T -> M3", "M3 -> Y",
    "K1 -> M2", "K3 -> Y", "C1 -> Y", "C3 -> Y",
    "T -> L1", "Y -> L3"
  ) %in% keys))
  expect_true(all(c(
    "X1 -> X2", "X2 -> X3", "M1 -> M2", "M2 -> M3",
    "K1 -> K2", "K2 -> K3", "C1 -> C2", "C2 -> C3",
    "L1 -> L2", "L2 -> L3"
  ) %in% keys))
})


test_that("last retains direct edges needed by structurally inferred roles", {
  dag <- build_graph(
    c("X1", "X2"), "T", "Y", mediators = c("M1", "M2"),
    mediator_outcome_confounders = c("K1", "K2"),
    competing_causes = c("C1", "C2"), colliders = c("L1", "L2"),
    type = "last"
  )
  keys <- edge_keys(dag)

  expect_true(all(c(
    "K1 -> M1", "K1 -> Y", "K2 -> M1", "K2 -> Y",
    "C1 -> Y", "C2 -> Y",
    "T -> L1", "Y -> L1", "T -> L2", "Y -> L2"
  ) %in% keys))
  roles <- get_roles(dag)
  expect_setequal(roles$confounder, c("X1", "X2"))
  expect_setequal(roles$mediator, c("M1", "M2"))
  expect_setequal(roles$mediator_outcome_confounder, c("K1", "K2"))
  expect_setequal(roles$competing_cause, c("C1", "C2"))
  expect_setequal(roles$collider, c("L1", "L2"))
})


test_that("every supported non-instrument role uses the same temporal grouping", {
  cases <- c(
    variables = "confounder",
    treatments = "treatment",
    outcomes = "outcome",
    mediators = "mediator",
    mediator_outcome_confounders = "mediator_outcome_confounder",
    competing_causes = "competing_cause",
    colliders = "collider",
    latents = "latent"
  )

  for (argument in names(cases)) {
    atomic_args <- list(
      variables = "X", treatments = "T", outcomes = "Y", type = "chain"
    )
    atomic_args[[argument]] <- c("A", "B")
    atomic <- do.call(build_graph, atomic_args)
    expect_true(
      role_edge(atomic, "A", "B"),
      info = paste("atomic order for", cases[[argument]])
    )

    parallel_args <- list(
      variables = "X", treatments = "T", outcomes = "Y", type = "chain"
    )
    parallel_args[[argument]] <- list(c("A", "B"))
    parallel <- do.call(build_graph, parallel_args)
    expect_false(
      any(c("A -> B", "B -> A") %in% edge_keys(parallel)),
      info = paste("parallel tier for", cases[[argument]])
    )
    expect_true(dagitty::isAcyclic(parallel))
  }
})


test_that("saturated adds every non-instrument edge allowed by the role order", {
  fully_ordered <- build_graph(
    "X", "T", "Y", mediators = c("M1", "M2", "M3"),
    type = "saturated"
  )
  expect_length(edge_keys(fully_ordered), choose(6L, 2L))
  expect_true(dagitty::isAcyclic(fully_ordered))

  dag <- build_graph(
    variables = "X", treatments = "T", outcomes = "Y",
    mediators = "M", instrumental_variables = "IV",
    mediator_outcome_confounders = "K", competing_causes = "C",
    colliders = "L", type = "saturated"
  )
  keys <- edge_keys(dag)

  expect_true(all(c(
    "X -> T", "X -> M", "X -> Y", "X -> L",
    "T -> M", "T -> Y", "T -> L",
    "K -> M", "K -> Y", "K -> L",
    "M -> Y", "M -> L",
    "C -> Y", "C -> L", "Y -> L"
  ) %in% keys))
  expect_false(any(c("X -> K", "T -> K", "C -> M") %in% keys))

  # IV inputs and generated edges remain unchanged.
  expect_identical(keys[grepl("^IV ", keys)], "IV -> T")
})


test_that("ordered and saturated route within and between roles", {
  arrows <- function(from, to) {
    as.vector(outer(from, to, function(v, w) paste(v, "->", w)))
  }
  ordered <- build_graph(
    c("X1", "X2"), "T", "Y", mediators = c("M1", "M2"),
    type = "ordered"
  )
  expected_ordered <- c(
    "X1 -> X2", "M1 -> M2",
    arrows(c("X1", "X2"), c("T", "Y")),
    arrows("T", c("M1", "M2", "Y")),
    arrows(c("M1", "M2"), "Y")
  )
  expect_setequal(edge_keys(ordered), expected_ordered)

  saturated <- build_graph(
    c("X1", "X2"), "T", "Y", mediators = c("M1", "M2"),
    type = "saturated"
  )
  expected_saturated <- c(
    expected_ordered,
    arrows(c("X1", "X2"), c("M1", "M2"))
  )
  expect_setequal(edge_keys(saturated), expected_saturated)
})


test_that("chain boundary routing leaves instrumental edges unchanged", {
  flat <- build_graph(
    c("X1", "X2"), c("T1", "T2"), "Y",
    instrumental_variables = c("IV1", "IV2"), type = "chain"
  )
  expect_setequal(
    edge_keys(flat)[grepl("^IV", edge_keys(flat))],
    c("IV1 -> T1", "IV2 -> T2")
  )

  mapped <- build_graph(
    c("X1", "X2"), c("T1", "T2"), "Y",
    instrumental_variables = list(list("IV1", c("T1", "T2"))),
    type = "chain"
  )
  expect_setequal(
    edge_keys(mapped)[grepl("^IV", edge_keys(mapped))],
    c("IV1 -> T1", "IV1 -> T2")
  )
})


test_that("last keeps mapped latents separate from temporal boundaries", {
  dag <- build_graph(
    "X", c("T1", "T2"), c("Y1", "Y2"),
    latents = list(
      "U1", "U2", list("U_map", c("T2", "Y2"))
    ),
    type = "last"
  )
  keys <- edge_keys(dag)

  expect_setequal(keys[grepl("^U", keys)], c(
    "U1 -> U2", "U2 -> T1", "U2 -> Y1",
    "U_map -> T2", "U_map -> Y2"
  ))
  expect_setequal(dagitty::latents(dag), c("U1", "U2", "U_map"))
})


test_that("routing annotations constrain arrows between roles but not temporal edges", {
  dag <- build_graph(
    "X", "T", "Y",
    mediators = c(incoming(c("M1", "M2")), outgoing("M3")),
    colliders = "C", type = "saturated"
  )
  keys <- edge_keys(dag)

  expect_true(all(c("T -> M1", "T -> M2", "X -> M1", "X -> M2") %in% keys))
  expect_false(any(c("T -> M3", "X -> M3") %in% keys))
  expect_true(all(c("M3 -> Y", "M3 -> C") %in% keys))
  expect_false(any(c("M1 -> Y", "M2 -> Y", "M1 -> C", "M2 -> C") %in% keys))
  expect_true(all(c("M1 -> M2", "M1 -> M3", "M2 -> M3") %in% keys))

  both <- build_graph(
    "X", "T", "Y", mediators = c(both = "B1", "B2"),
    type = "saturated"
  )
  both_keys <- edge_keys(both)
  expect_true(all(c("T -> B1", "X -> B1", "B1 -> Y", "B1 -> B2") %in%
                    both_keys))
  expect_false(any(c("T -> B2", "X -> B2", "B2 -> Y") %in% both_keys))
  expect_error(
    build_graph("X", "T", "Y",
                mediators = c(both = c("B1", "B2")),
                type = "saturated"),
    "scalar annotation"
  )
})


test_that("full applies routing markers to bidirected edges", {
  dag <- build_graph(
    "X", "T", "Y",
    mediators = c(incoming("M1"), outgoing("M2"), "M3"),
    type = "full"
  )

  expect_true(symmetric_role_edge(dag, "T", "M1"))
  expect_false(symmetric_role_edge(dag, "T", "M2"))
  expect_false(symmetric_role_edge(dag, "T", "M3"))
  expect_true(symmetric_role_edge(dag, "M2", "Y"))
  expect_false(symmetric_role_edge(dag, "M1", "Y"))
  expect_false(symmetric_role_edge(dag, "M3", "Y"))
})


test_that("build_graph() and add_nodes() with a role generate the same edges", {
  input <- list("M1", c("M2", "M3"), "M4")
  for (type in c("full", "chain", "ordered", "saturated")) {
    built <- build_graph("X", "T", "Y", mediators = input, type = type)
    base <- build_graph("X", "T", "Y", type = type)
    added <- add_nodes(
      base, nodes = input, role = "mediator", type = type,
      print_edges = FALSE
    )
    expect_identical(edge_keys(added), edge_keys(built), info = type)
  }
})


test_that("ADD and CONNECT keep their own graph types", {
  base <- build_graph("X", "T", "Y", mediators = "M1", type = "last")

  for (type in c("full", "chain", "ordered", "saturated")) {
    expect_s3_class(
      add_nodes(base, nodes = "N", role = "mediator", type = type,
                print_edges = FALSE),
      "dagitty", exact = FALSE
    )
  }

  expect_error(
    add_nodes(base, nodes = "N", role = "mediator", type = "middle",
              print_edges = FALSE),
    "Invalid type"
  )
  expect_error(
    connect_nodes(base, nodes = "M1", type = "middle", print_edges = FALSE),
    "Invalid type"
  )
  expect_warning(
    add_nodes(base, nodes = "N", role = "mediator", type = "first",
              print_edges = FALSE),
    paste0("add_nodes(): 'type = \"first\"' is deprecated as of 0.0.1. ",
           "Use 'position = first()' instead."),
    fixed = TRUE
  )
  expect_warning(
    connect_nodes(base, nodes = "M1", type = "last", print_edges = FALSE),
    paste0("connect_nodes(): 'type = \"last\"' is deprecated as of 0.0.1. ",
           "Use 'type = \"ordered\"' instead."),
    fixed = TRUE
  )
})


test_that("full connects the new or selected nodes to every other node", {
  base <- build_graph("X", "T", "Y", type = "full")
  added <- add_nodes(
    base, nodes = c("N1", "N2"), role = "mediator", type = "full",
    print_edges = FALSE
  )

  expect_setequal(edge_keys(added), c(
    "T <-> Y", "T <-> X", "X <-> Y",
    "N1 <-> N2", "N1 <-> T", "N1 <-> Y", "N1 <-> X",
    "N2 <-> T", "N2 <-> Y", "N2 <-> X"
  ))

  directed <- build_graph("X", "T", "Y", mediators = "M", type = "last")
  connected <- connect_nodes(
    directed, nodes = "X", type = "full", print_edges = FALSE
  )

  expect_true(all(edge_keys(directed) %in% edge_keys(connected)))
  expect_true(all(c("M <-> X", "T <-> X", "X <-> Y") %in%
                    edge_keys(connected)))
})


test_that("full leaves instrumental variable edges to the IV rules", {
  directed <- build_graph(
    "X", "T", "Y", mediators = "M", instrumental_variables = "IV",
    type = "last"
  )
  connected <- connect_nodes(
    directed, nodes = "M", type = "full", print_edges = FALSE
  )
  keys <- edge_keys(connected)

  expect_true("IV -> T" %in% keys)
  expect_false(any(grepl("IV <->", keys)))
})


test_that("ADD and CONNECT do not remove existing arrows when applying chain", {
  base <- build_graph(
    "X", "T", "Y", mediators = c("M1", "M2"), type = "last"
  )
  added <- add_nodes(
    base, nodes = "M3", role = "mediator", type = "chain",
    print_edges = FALSE
  )
  expect_true(all(edge_keys(base) %in% edge_keys(added)))
  expect_true(all(c("M2 -> Y", "M2 -> M3") %in% edge_keys(added)))

  disconnected <- remove_edges(base, "M1 -> M2")
  connected <- connect_nodes(
    disconnected, nodes = c("M1", "M2"), type = "chain", print_edges = FALSE
  )
  expect_true(all(edge_keys(disconnected) %in% edge_keys(connected)))
})


test_that("ADD preserves parallel tiers across operations", {
  base <- build_graph(
    "X", "T", "Y",
    mediators = list("M0", c("M1", "M2"), "M3"), type = "last"
  )
  added <- add_nodes(
    base, nodes = list(c("M4", "M5")), role = "mediator",
    type = "ordered", print_edges = FALSE
  )
  positions <- test_role_positions(added, "mediator")
  expect_identical(positions, list("M0", c("M1", "M2"), "M3",
                                   c("M4", "M5")))
  expect_true(all(c("M3 -> M4", "M3 -> M5") %in% edge_keys(added)))
  expect_false(any(c("M1 -> M2", "M2 -> M1", "M4 -> M5", "M5 -> M4") %in%
                     edge_keys(added)))

  expect_identical(test_role_positions(base, "mediator"),
                   list("M0", c("M1", "M2"), "M3"))
})


test_that("ADD uses routing annotations for supplied roles", {
  base <- build_graph("X", "T", "Y", mediators = "M0", type = "saturated")
  expect_warning(
    added <- add_nodes(
      base,
      nodes = c(incoming("M0"), "M1", outgoing("M2")),
      role = "mediator", type = "saturated", print_edges = FALSE
    ),
    "skipping"
  )
  added_keys <- edge_keys(added)
  expect_true(all(c("T -> M0", "X -> M0", "M2 -> Y") %in% added_keys))
  expect_false(any(c("T -> M1", "X -> M1") %in% added_keys))
})


test_that("routing marks apply to the operation that supplies them", {
  base <- build_graph(
    "X", "T", "Y", mediators = incoming("M0"), type = "saturated"
  )
  expect_true(all(c("T -> M0", "X -> M0") %in% edge_keys(base)))

  unmarked_add <- add_nodes(
    base, nodes = "M1", role = "mediator", type = "saturated",
    print_edges = FALSE
  )
  expect_true(all(c("T -> M1", "X -> M1", "M1 -> Y") %in%
                    edge_keys(unmarked_add)))

  marked_add <- add_nodes(
    base, nodes = c(incoming("M1"), outgoing("M2")), role = "mediator",
    type = "saturated", print_edges = FALSE
  )
  expect_true(all(c("T -> M1", "X -> M1", "M2 -> Y") %in%
                    edge_keys(marked_add)))
  expect_false(any(c("T -> M2", "X -> M2", "M1 -> Y") %in%
                     edge_keys(marked_add)))
})


test_that("connect_nodes() applies chain and saturated rules", {
  sequence <- build_graph(c("X1", "X2", "X3"), "T", "Y", type = "last")
  expect_identical(test_role_positions(sequence, "confounder"),
                   list("X1", "X2", "X3"))
  without_chain <- keep_edges(
    setdiff(edge_keys(sequence), c("X1 -> X2", "X2 -> X3")), sequence
  )
  expect_identical(names(without_chain), names(sequence))
  restored <- connect_nodes(without_chain, type = "chain", print_edges = FALSE)
  expect_identical(edge_keys(restored), edge_keys(without_chain))

  ordered <- build_graph("X", "T", "Y", mediators = "M", type = "ordered")
  saturated <- connect_nodes(ordered, type = "saturated", print_edges = FALSE)
  expect_true("X -> M" %in% edge_keys(saturated))
})


test_that("connect_nodes() warns for the arguments it no longer takes", {
  base <- build_graph("X", "T", "Y", mediators = "M", type = "ordered")
  plain <- connect_nodes(base, nodes = "X", type = "saturated",
                         print_edges = FALSE)

  expect_warning(
    with_role <- connect_nodes(base, nodes = "X", type = "saturated",
                               role = "confounder", print_edges = FALSE),
    "connect_nodes(): 'role' is deprecated as of 0.0.2 and was ignored.",
    fixed = TRUE
  )
  expect_identical(edge_keys(with_role), edge_keys(plain))
  expect_identical(dagitty::exposures(with_role), dagitty::exposures(base))
  expect_identical(dagitty::outcomes(with_role), dagitty::outcomes(base))
  expect_identical(dagitty::latents(with_role), dagitty::latents(base))

  expect_warning(
    with_alias <- connect_nodes(base, nodes = "X", type = "saturated",
                                node_role = "confounder", print_edges = FALSE),
    "connect_nodes(): 'node_role' is deprecated as of 0.0.2 and was ignored.",
    fixed = TRUE
  )
  expect_identical(edge_keys(with_alias), edge_keys(plain))
  expect_identical(dagitty::exposures(with_alias), dagitty::exposures(base))
  expect_identical(dagitty::outcomes(with_alias), dagitty::outcomes(base))
  expect_identical(dagitty::latents(with_alias), dagitty::latents(base))

  expect_warning(
    with_position <- connect_nodes(base, nodes = "X", type = "saturated",
                                   position = first(), print_edges = FALSE),
    "connect_nodes(): 'position' is deprecated as of 0.0.2 and was ignored.",
    fixed = TRUE
  )
  expect_identical(edge_keys(with_position), edge_keys(plain))
  expect_identical(dagitty::exposures(with_position),
                   dagitty::exposures(base))
  expect_identical(dagitty::outcomes(with_position), dagitty::outcomes(base))
  expect_identical(dagitty::latents(with_position), dagitty::latents(base))
})


test_that("a supplied role no longer sets a node status", {
  base <- build_graph("X", "T", "Y", mediators = "M", type = "ordered")

  for (deprecated_role in c("treatment", "outcome", "latent")) {
    expect_warning(
      connected <- connect_nodes(base, nodes = "X", type = "saturated",
                                 role = deprecated_role, print_edges = FALSE),
      "connect_nodes(): 'role' is deprecated as of 0.0.2 and was ignored.",
      fixed = TRUE
    )
    expect_identical(dagitty::exposures(connected), dagitty::exposures(base))
    expect_identical(dagitty::outcomes(connected), dagitty::outcomes(base))
    expect_identical(dagitty::latents(connected), dagitty::latents(base))
  }
})


test_that("connect_nodes() rejects arguments it does not know", {
  base <- build_graph("X", "T", "Y", type = "ordered")

  expect_error(
    connect_nodes(base, nodes = "X", nodes_role = "confounder",
                  print_edges = FALSE),
    "connect_nodes(): unused argument(s): nodes_role.",
    fixed = TRUE
  )
  expect_error(
    connect_nodes(base, "X", "full", FALSE, "confounder"),
    "connect_nodes(): unused argument(s): argument 1.",
    fixed = TRUE
  )
})


test_that("connect_nodes() reads the deprecated first and last types as ordered", {
  base <- build_graph("X", "T", "Y", mediators = "M1", type = "last")
  ordered <- connect_nodes(base, nodes = "M1", type = "ordered",
                           print_edges = FALSE)

  expect_warning(
    legacy_first <- connect_nodes(base, nodes = "M1", type = "first",
                                  print_edges = FALSE),
    paste0("connect_nodes(): 'type = \"first\"' is deprecated as of 0.0.1. ",
           "Use 'type = \"ordered\"' instead."),
    fixed = TRUE
  )
  expect_identical(edge_keys(legacy_first), edge_keys(ordered))

  expect_warning(
    legacy_last <- connect_nodes(base, nodes = "M1", type = "last",
                                 print_edges = FALSE),
    paste0("connect_nodes(): 'type = \"last\"' is deprecated as of 0.0.1. ",
           "Use 'type = \"ordered\"' instead."),
    fixed = TRUE
  )
  expect_identical(edge_keys(legacy_last), edge_keys(ordered))
})


test_that("deep latent target mappings remain distinct from temporal tiers", {
  mapped <- build_graph(
    "X", "T", "Y",
    latents = list(list("U", c("T", "Y"))),
    type = "saturated"
  )
  expect_true(all(c("U -> T", "U -> Y") %in% edge_keys(mapped)))
  expect_identical(dagitty::latents(mapped), "U")

  parallel <- build_graph(
    "X", "T", "Y",
    latents = list(c("U1", "U2")),
    type = "saturated"
  )
  expect_setequal(dagitty::latents(parallel), c("U1", "U2"))
  expect_false(any(c("U1 -> U2", "U2 -> U1") %in% edge_keys(parallel)))

  base <- build_graph("X", "T", "Y", type = "ordered")
  added <- add_nodes(
    base, nodes = list(list("U", c("T", "Y"))), role = "latent",
    type = "saturated", print_edges = FALSE
  )
  expect_true(all(c("U -> T", "U -> Y") %in% edge_keys(added)))
  expect_identical(dagitty::latents(added), "U")
  expect_error(
    add_nodes(
      base, nodes = list(list("U", c("T", "unknown"))), role = "latent",
      print_edges = FALSE
    ),
    "targets are not in the graph"
  )
})


test_that("mapped and temporal latent nodes retain separate edge rules", {
  mixed <- build_graph(
    "X", "T", "Y", mediators = "M",
    latents = list("U_flat", list("U_map", "M")),
    type = "saturated"
  )
  keys <- edge_keys(mixed)
  expect_true(all(c("U_flat -> X", "U_flat -> T", "U_flat -> Y",
                    "U_map -> M") %in% keys))
  expect_false(any(c("U_map -> X", "U_map -> T", "U_map -> Y",
                     "U_flat -> U_map", "U_map -> U_flat") %in% keys))
  expect_setequal(
    role_mapping_sources(mixed, "latent"), c("U_flat", "U_map")
  )

  mapped <- build_graph(
    "X", "T", "Y", mediators = "M",
    latents = list(list("U_map", "M")), type = "saturated"
  )
  extended <- add_nodes(
    mapped, nodes = "U_flat", role = "latent", type = "saturated",
    print_edges = FALSE
  )
  expect_true(all(c("U_map -> M", "U_flat -> X", "U_flat -> T",
                    "U_flat -> Y") %in% edge_keys(extended)))
  expect_false(any(c("U_map -> X", "U_map -> T", "U_map -> Y") %in%
                     edge_keys(extended)))
})


test_that("non-instrument changes preserve existing IV edges without adding new ones", {
  base <- build_graph(
    "X", "T", "Y", mediators = "M", instrumental_variables = "IV",
    type = "saturated"
  )
  added <- add_nodes(
    base, nodes = "T2", role = "treatment", type = "saturated",
    print_edges = FALSE
  )
  expect_identical(sort(edge_keys(added)[grepl("^IV ", edge_keys(added))]),
                   "IV -> T")

  role_free <- connect_nodes(
    base, nodes = names(base), type = "saturated", print_edges = FALSE
  )
  expect_identical(
    sort(edge_keys(role_free)[grepl("^IV ", edge_keys(role_free))]),
    "IV -> T"
  )
  expect_identical(get_roles(role_free)$instrument, "IV")
})


test_that("BUILD keeps legacy IV edges when treatments use temporal tiers", {
  treatment_inputs <- list(
    atomic = c("T1", "T2"),
    sequential = list("T1", "T2"),
    parallel = list(c("T1", "T2"))
  )
  iv_inputs <- list(
    flat = c("IV1", "IV2", "IV3"),
    mapped = list(list("IV1", c("T1", "T2")))
  )

  for (type in c("full", "ordered", "saturated")) {
    for (iv_name in names(iv_inputs)) {
      reference <- build_graph(
        "X", treatment_inputs$atomic, "Y",
        instrumental_variables = iv_inputs[[iv_name]], type = type
      )
      reference_iv <- sort(edge_keys(reference)[
        grepl("^IV", edge_keys(reference))
      ])
      for (treatment_name in c("sequential", "parallel")) {
        tiered <- build_graph(
          "X", treatment_inputs[[treatment_name]], "Y",
          instrumental_variables = iv_inputs[[iv_name]], type = type
        )
        tiered_iv <- sort(edge_keys(tiered)[grepl("^IV", edge_keys(tiered))])
        expect_identical(
          tiered_iv, reference_iv,
          info = paste(type, iv_name, treatment_name)
        )
      }
    }
  }
})


test_that("new argument names and deprecated aliases produce the same graph", {
  expect_identical(
    names(formals(add_nodes))[1:8],
    c("dag", "nodes", "parents", "children", "role", "type",
      "print_edges", "position")
  )
  expect_identical(
    names(formals(connect_nodes)),
    c("dag", "nodes", "type", "print_edges", "...")
  )

  base <- build_graph("X", "T", "Y", type = "ordered")
  canonical <- add_nodes(
    base, nodes = "A", parents = "X", children = "T",
    print_edges = FALSE
  )
  expect_warning(
    aliased <- add_nodes(base, new_nodes = "A", parents = "X", children = "T",
                         print_edges = FALSE),
    "add_nodes(): 'new_nodes' is deprecated as of 0.0.2. Use 'nodes' instead.",
    fixed = TRUE
  )
  expect_identical(edge_keys(aliased), edge_keys(canonical))

  expect_warning(
    role_alias <- add_nodes(base, nodes = "X2", node_role = "confounder",
                            print_edges = FALSE),
    "add_nodes(): 'node_role' is deprecated as of 0.0.2. Use 'role' instead.",
    fixed = TRUE
  )
  expect_true("X2 -> T" %in% edge_keys(role_alias))

  expect_warning(
    parent_alias <- add_nodes(base, nodes = "A", ancestors = "X", children = "T",
                              print_edges = FALSE),
    "add_nodes(): 'ancestors' is deprecated as of 0.0.2. Use 'parents' instead.",
    fixed = TRUE
  )
  expect_identical(edge_keys(parent_alias), edge_keys(canonical))
  expect_warning(
    child_alias <- add_nodes(base, nodes = "A", parents = "X", descendants = "T",
                             print_edges = FALSE),
    "add_nodes(): 'descendants' is deprecated as of 0.0.2. Use 'children' instead.",
    fixed = TRUE
  )
  expect_identical(edge_keys(child_alias), edge_keys(canonical))

  expect_error(
    add_nodes(base, nodes = "A", new_nodes = "B", print_edges = FALSE),
    "not both"
  )
  expect_error(
    add_nodes(base, nodes = "A", parents = "X", ancestors = "X",
              print_edges = FALSE),
    "not both"
  )
  expect_error(
    add_nodes(base, nodes = "A", children = "T", descendants = "T",
              print_edges = FALSE),
    "add_nodes(): supply either 'children' or its deprecated alias 'descendants', not both.",
    fixed = TRUE
  )
  expect_error(
    add_nodes(base, nodes = "A", role = "confounder",
              node_role = "confounder", print_edges = FALSE),
    "not both"
  )
})


test_that("positions work for non-instrument roles and observed nodes stay explicit", {
  collider_base <- build_graph(
    "X", "T", "Y", colliders = "C1", type = "last"
  )
  collider <- add_nodes(
    collider_base, nodes = "C0", role = "collider", type = "chain",
    position = first(), print_edges = FALSE
  )
  expect_true("C0 -> C1" %in% edge_keys(collider))

  observed_base <- dagitty::dagitty("dag { T [exposure] Y [outcome] O1 T -> Y }")
  observed <- add_nodes(
    observed_base, nodes = "O2", role = "observed", type = "ordered",
    position = first(), print_edges = FALSE
  )
  expect_false(any(c("O1 -> O2", "O2 -> O1") %in% edge_keys(observed)))
})


test_that("construction and connection warn on directed cycles", {
  expect_warning(
    cyclic_build <- build_graph(
      variables = NA, treatments = c("A", "B"), outcomes = "Y",
      mediators = c("B", "A"), type = "ordered"
    ),
    "build_graph(): the supplied roles and temporal positions create a directed cycle.",
    fixed = TRUE
  )
  expect_false(dagitty::isAcyclic(cyclic_build))

  chain <- string_to_dag("T -> Y")
  expect_warning(
    expect_warning(
      cyclic_add <- add_nodes(chain, nodes = "W", parents = "Y", children = "T",
                              print_edges = FALSE),
      "add_nodes(): the requested edges create a directed cycle.",
      fixed = TRUE
    ),
    "acyclic subset", fixed = TRUE
  )
  expect_false(dagitty::isAcyclic(cyclic_add))
  expect_true(all(c("Y -> W", "W -> T") %in% edge_keys(cyclic_add)))

  dag <- dagitty::dagitty(
    "dag { U [latent] T [exposure] Y [outcome] T -> Y T -> U }"
  )
  expect_warning(
    cyclic <- connect_nodes(dag, type = "saturated", print_edges = FALSE),
    "The connected graph contains a directed cycle.",
    fixed = TRUE
  )
  expect_true("U -> T" %in% edge_keys(cyclic))
  expect_false(dagitty::isAcyclic(cyclic))
  expect_true(dagitty::isAcyclic(dag))
})


test_that("role inputs accept incoming and outgoing as names", {
  base <- build_graph("X", "T", "Y", type = "saturated")

  helpers <- add_nodes(
    base, nodes = c("M1", incoming("M2"), outgoing("M3")),
    role = "mediator", type = "saturated", print_edges = FALSE
  )
  equations <- add_nodes(
    base, nodes = c("M1", incoming = "M2", outgoing = "M3"),
    role = "mediator", type = "saturated", print_edges = FALSE
  )
  expect_identical(edge_keys(equations), edge_keys(helpers))
  expect_true(all(c("T -> M2", "X -> M2", "M3 -> Y") %in%
                    edge_keys(equations)))
  expect_false(any(c("T -> M1", "T -> M3", "M1 -> Y") %in%
                     edge_keys(equations)))

  numbered <- add_nodes(
    base, nodes = c(incoming = c("M2", "M3")), role = "mediator",
    type = "saturated", print_edges = FALSE
  )
  expect_true(all(c("T -> M2", "X -> M2", "T -> M3", "X -> M3") %in%
                    edge_keys(numbered)))

  marked_list <- add_nodes(
    base, nodes = list(incoming = "M1", outgoing = "M2"), role = "mediator",
    type = "saturated", print_edges = FALSE
  )
  expect_identical(
    edge_keys(marked_list),
    edge_keys(add_nodes(
      base, nodes = list(incoming("M1"), outgoing("M2")), role = "mediator",
      type = "saturated", print_edges = FALSE
    ))
  )

  grouped <- add_nodes(
    base, nodes = list(incoming = c("M1", "M2")), role = "mediator",
    type = "saturated", print_edges = FALSE
  )
  expect_true(all(c("T -> M1", "X -> M1", "T -> M2", "X -> M2") %in%
                    edge_keys(grouped)))
  expect_false(any(c("M1 -> M2", "M2 -> M1") %in% edge_keys(grouped)))
})


test_that("position accepts first, last, before and after as values", {
  base <- build_graph(
    "X", "T", "Y", mediators = list("M1", "M2", "M3"), type = "ordered"
  )
  place <- function(position){
    edge_keys(add_nodes(
      base, nodes = "M9", role = "mediator", type = "ordered",
      position = position, print_edges = FALSE
    ))
  }
  expect_identical(place("first"), place(first()))
  expect_identical(place("last"), place(last()))
  expect_identical(place(c(before = "M2")), place(before("M2")))
  expect_identical(place(c(after = "M2")), place(after("M2")))
  expect_true(all(c("M1 -> M9", "M9 -> M2") %in% place(c(before = "M2"))))

  together <- add_nodes(
    base, nodes = list(c("M8", "M9")), role = "mediator", type = "ordered",
    position = "first", print_edges = FALSE
  )
  expect_true(all(c("T -> M8", "T -> M9", "M8 -> M1", "M9 -> M1") %in%
                    edge_keys(together)))

  expect_error(
    add_nodes(base, nodes = "M9", role = "mediator", position = "middle",
              print_edges = FALSE),
    "use first"
  )
})


test_that("add_nodes() rejects position for the instrument role", {
  dag <- build_graph(variables = "X", treatments = "T", outcomes = "Y",
                     type = "saturated")
  expect_error(
    add_nodes(dag, nodes = "IV", role = "instrument", position = first(),
              print_edges = FALSE),
    paste0("The position argument is not supported for role \"instrument\"; ",
           "it would be ignored. Remove position, or use a supported role."),
    fixed = TRUE
  )
})


test_that("add_nodes() connects new instruments to the treatments", {
  dag <- build_graph(confounders = "X", treatments = c("T1", "T2"),
                     outcomes = "Y", type = "saturated")
  instrument_edges <- function(type){
    added <- add_nodes(dag, c("IV1", "IV2"), role = "instrument", type = type)
    edges <- dagitty::edges(added)
    edges <- edges[edges$v %in% c("IV1", "IV2") | edges$w %in% c("IV1", "IV2"), ]
    sort(paste(edges$v, edges$e, edges$w))
  }
  expect_identical(instrument_edges("saturated"),
                   c("IV1 -> T1", "IV1 -> T2", "IV2 -> T1", "IV2 -> T2"))
  expect_identical(instrument_edges("full"),
                   c("IV1 <-> T1", "IV1 <-> T2", "IV2 <-> T1", "IV2 <-> T2"))
  expect_identical(instrument_edges("chain"), c("IV1 -> T1", "IV2 -> T2"))
  expect_identical(instrument_edges("ordered"), c("IV1 -> T1", "IV2 -> T2"))

  # Other roles are planned by role_operation_edges(), never by this helper.
  expect_error(add_nodes_helper(dag, "N", "confounder"),
               "only connects instrument nodes", fixed = TRUE)
})
