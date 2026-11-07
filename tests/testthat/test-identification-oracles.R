################################################################################
# test-identification-oracles.R


# Compare identification results with independent causal graph oracles

# Computational checks use the private worker through identification_result().
# Public overview behaviour has its own tests.
test_that("IDC matches published conditional-effect oracles", {
  # Shpitser and Pearl (2008), Figure 6(a), p. 1954.
  condition_helps <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z T -> Z Z -> Y T <-> Z }"
  )

  marginal_helps <- identification_result(condition_helps)
  conditional_helps <- identification_result(
    condition_helps, conditioning = "Z"
  )

  expect_false(attr(marginal_helps, "observationally_identifiable"))
  expect_false(attr(marginal_helps, "identifiable"))
  expect_valid_identification_hedge(marginal_helps, condition_helps)
  expect_true(attr(conditional_helps, "observationally_identifiable"))
  expect_true(attr(conditional_helps, "identifiable"))
  expect_true(attr(conditional_helps, "conditional"))
  expect_identical(attr(conditional_helps, "conditioning_intervened"), "Z")
  expect_identical(attr(conditional_helps, "conditioning_remaining"),
                   character())
  expect_null(attr(conditional_helps, "hedge"))

  # Shpitser and Pearl (2008), Figure 6(b), p. 1954.
  condition_hurts <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z T -> Z Y -> Z T <-> Z }"
  )

  marginal_hurts <- identification_result(condition_hurts)
  conditional_hurts <- identification_result(
    condition_hurts, conditioning = "Z"
  )

  expect_true(attr(marginal_hurts, "observationally_identifiable"))
  expect_true(attr(marginal_hurts, "identifiable"))
  expect_false(attr(conditional_hurts, "observationally_identifiable"))
  expect_false(attr(conditional_hurts, "identifiable"))
  expect_identical(attr(conditional_hurts, "conditioning_intervened"),
                   character())
  expect_identical(attr(conditional_hurts, "conditioning_remaining"), "Z")
  expect_valid_identification_hedge(conditional_hurts, condition_hurts)

  # Shpitser and Pearl (2006), Figure 1(c), PDF p. 2.
  dag <- dagitty::dagitty(
    "dag { T [exposure] W [outcome] Y Z
           T -> Z Y -> Z W -> Y T <-> Z }"
  )

  marginal <- identification_result(dag)
  conditional <- identification_result(dag, conditioning = "Z")

  expect_true(attr(marginal, "observationally_identifiable"))
  expect_false(attr(conditional, "observationally_identifiable"))
  expect_false(attr(conditional, "identifiable"))
  expect_identical(attr(conditional, "conditioning_intervened"), character())
  expect_identical(attr(conditional, "conditioning_remaining"), "Z")
  expect_valid_identification_hedge(conditional, dag)
})

test_that("ID matches published marginal-effect oracles", {
  # Shpitser and Pearl (2008), Figure 1(a), p. 1946.
  bow <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> Y T <-> Y }"
  )

  result <- identification_result(bow)

  expect_false(attr(result, "observationally_identifiable"))
  expect_false(attr(result, "identifiable"))
  expect_valid_identification_hedge(result, bow)
  hedge <- attr(result, "hedge")
  expect_setequal(names(hedge$F), c("T", "Y"))
  expect_setequal(identification_edge_keys(hedge$F),
                  c("T -> Y", "T <-> Y"))
  expect_identical(names(hedge$F_prime), "Y")
  expect_length(identification_edge_keys(hedge$F_prime), 0L)
  expect_identical(hedge$roots, "Y")
})

test_that("failed ID returns a validated edge-bearing hedge", {
  dag <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Q [outcome] W1 W2 Y1 Y2 ",
      "W1 -> T T -> Y1 W1 -> W2 W2 -> Y2 Y1 -> Q Y2 -> Q ",
      "W1 <-> Y1 W1 <-> W2 W1 <-> Y2 T <-> W2 }"
    )
  )

  result <- identification_result(dag)

  expect_false(attr(result, "observationally_identifiable"))
  expect_valid_identification_hedge(result, dag)
  hedge <- attr(result, "hedge")
  expect_setequal(names(hedge$F), c("W1", "T", "W2", "Y1", "Y2"))
  expect_setequal(
    identification_edge_keys(hedge$F),
    c(
      "W1 -> W2", "W2 -> Y2", "T -> Y1",
      "W1 <-> W2", "W1 <-> Y1", "W1 <-> Y2", "T <-> W2"
    )
  )
  expect_setequal(names(hedge$F_prime), c("W1", "W2", "Y1", "Y2"))
  expect_setequal(
    identification_edge_keys(hedge$F_prime),
    c("W1 -> W2", "W2 -> Y2", "W1 <-> W2", "W1 <-> Y1",
      "W1 <-> Y2")
  )
  expect_setequal(hedge$roots, c("Y1", "Y2"))
})

test_that("general ID succeeds beyond adjustment and front-door", {
  # Shpitser and Pearl (2008), Figure 2(f), p. 1948.
  dag <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] Z1 Z2 ",
      "T -> Z1 Z1 -> Z2 Z1 -> Y Z2 -> Y T -> Y ",
      "T <-> Z2 Z1 <-> Y }"
    )
  )

  result <- identification_result(dag)

  expect_true(attr(result, "observationally_identifiable"))
  expect_true(attr(result, "identifiable"))
  expect_true(result$feasible[1L])
  expect_false(identification_strategy(result, "adjust")$feasible)
  expect_false(identification_strategy(result, "front-door")$feasible)
  expect_null(attr(result, "hedge"))

  subset_dag <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] M N ",
      "T -> M M -> N N -> Y T <-> Y T <-> N }"
    )
  )
  subset_result <- identification_result(subset_dag)
  subset_frontdoor <- identification_strategy(subset_result, "front-door")

  expect_true(attr(subset_result, "observationally_identifiable"))
  expect_true(subset_frontdoor$feasible)
  expect_identical(subset_frontdoor$variable_set, "{M}")

  parallel_dag <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] M1 M2 ",
      "T -> M1 M1 -> Y T -> M2 M2 -> Y T <-> Y }"
    )
  )
  parallel_result <- identification_result(parallel_dag)
  parallel_frontdoor <- identification_strategy(
    parallel_result, "front-door"
  )

  expect_true(attr(parallel_result, "observationally_identifiable"))
  expect_true(parallel_frontdoor$feasible)
  expect_identical(parallel_frontdoor$variable_set, "{M1, M2}")
})

test_that("z-identification matches published surrogate-experiment oracles", {
  # Bareinboim and Pearl (2016), Figure 3A, pp. 7348--7349.
  graph_a <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z Z -> T T -> Y T <-> Z Z <-> Y }"
  )
  result_a <- identification_result(graph_a, experiments = "Z")
  expect_false(attr(result_a, "observationally_identifiable"))
  expect_true(attr(result_a, "z_identifiable"))
  expect_true(attr(result_a, "identifiable"))
  expect_true(identification_strategy(result_a, "z-ident|surrogate")$feasible)

  # Bareinboim and Pearl (2016), Figure 3B, pp. 7348--7349.
  graph_b <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] Z1 Z2 ",
      "Z1 -> T Z2 -> T T -> Y ",
      "T <-> Z1 Z1 <-> Z2 Z1 <-> Y Z2 <-> Y }"
    )
  )
  result_b <- identification_result(
    graph_b, experiments = "Z1"
  )
  expect_false(attr(result_b, "observationally_identifiable"))
  expect_true(attr(result_b, "z_identifiable"))
  expect_true(attr(result_b, "identifiable"))
  z_b <- identification_strategy(result_b, "z-ident|surrogate")
  expect_true(z_b$feasible)
  expect_identical(z_b$variable_set, "{Z1}")

  result_b_multiple <- identification_result(
    graph_b,
    experiments = c("Z1", "Z2")
  )
  z_b_multiple <- identification_strategy(
    result_b_multiple,
    "z-ident|surrogate"
  )
  expect_true(attr(result_b_multiple, "z_identifiable"))
  expect_true(z_b_multiple$feasible)
  expect_identical(z_b_multiple$variable_set, "{Z1}")

  # Bareinboim and Pearl (2016), Figure 3C, pp. 7348--7349.
  graph_c <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z Z -> T T -> Y T <-> Y Z <-> Y }"
  )
  result_c <- identification_result(graph_c, experiments = "Z")
  expect_false(attr(result_c, "observationally_identifiable"))
  expect_false(attr(result_c, "z_identifiable"))
  expect_false(attr(result_c, "identifiable"))
  expect_false(identification_strategy(result_c,
                                       "z-ident|surrogate")$feasible)

  # Bareinboim and Pearl (2016), Figure 3D, pp. 7348--7349.
  graph_d <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] Z T -> Y Z -> Y T <-> Z Z <-> Y }"
  )
  result_d <- identification_result(graph_d, experiments = "Z")
  expect_false(attr(result_d, "observationally_identifiable"))
  expect_false(attr(result_d, "z_identifiable"))
  expect_false(attr(result_d, "identifiable"))
  expect_false(identification_strategy(result_d,
                                       "z-ident|surrogate")$feasible)
})

test_that("direct effect identification uses the outcome parent intervention", {
  direct <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] M T -> M M -> Y T -> Y }"
  )
  direct_result <- identification_result(direct, effect = "direct")

  expect_identification_shape(direct_result)
  expect_identical(attr(direct_result, "effect"), "direct")
  expect_true(attr(direct_result, "observationally_identifiable"))
  expect_true(attr(direct_result, "identifiable"))
  expect_true(direct_result$feasible[1L])
  expect_identical(direct_result$variable_set[1L], "fixing set: {M}")
  direct_roles <- attr(direct_result, "roles")
  expect_equal(
    nrow(direct_roles[variable == "M" & role == "outcome parent"]),
    1L
  )

  projected_direct <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] Z L [latent] ",
      "T -> L Z -> L L -> Y }"
    )
  )
  projected_direct_result <- identification_result(
    projected_direct, effect = "direct"
  )
  projected_roles <- attr(projected_direct_result, "roles")
  expect_true(attr(projected_direct_result,
                   "observationally_identifiable"))
  expect_true(attr(projected_direct_result, "identifiable"))
  expect_equal(
    nrow(projected_roles[
      variable == "Z" & role == "outcome parent"
    ]),
    1L
  )

  projected_bow <- dagitty::dagitty(
    paste0(
      "dag { T [exposure] Y [outcome] U [latent] L [latent] ",
      "U -> T U -> Y T -> L L -> Y }"
    )
  )
  projected_bow_result <- identification_result(
    projected_bow, effect = "direct"
  )
  expect_false(attr(projected_bow_result,
                    "observationally_identifiable"))
  expect_false(attr(projected_bow_result, "identifiable"))
  expect_type(attr(projected_bow_result, "hedge"), "list")

  no_direct_edge <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] M T -> M M -> Y }"
  )
  absent_result <- identification_result(no_direct_edge, effect = "direct")
  expect_true(attr(absent_result, "observationally_identifiable"))
  expect_true(attr(absent_result, "identifiable"))
  expect_true(absent_result$feasible[1L])
  expect_identical(
    absent_result$note[1L],
    paste(
      "no direct arrow in the latent-projected graph;",
      "ID identifies the outcome distribution under intervention on",
      "the treatment and fixing set in the latent-projected graph"
    )
  )

  linear_single_door <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] A T -> Y A -> Y A <-> Y }"
  )
  nonparametric_result <- identification_result(
    linear_single_door,
    effect = "direct"
  )
  expect_false(attr(nonparametric_result, "identifiable"))
  expect_false(any(
    nonparametric_result$strategy == "single-door adjustment"
  ))

  confounded <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> Y T <-> Y }"
  )
  confounded_result <- identification_result(confounded, effect = "direct")
  expect_false(attr(confounded_result, "observationally_identifiable"))
  expect_false(attr(confounded_result, "identifiable"))
  expect_valid_identification_hedge(confounded_result, confounded)
})

test_that("front-door returns before direct-arrow subset enumeration", {
  graph <- new.env(parent = emptyenv())
  graph$D <- data.frame(from = "T", to = "Y")
  makeActiveBinding("V", function(value) {
    stop("front-door candidate enumeration started")
  }, graph)

  expect_identical(
    id_frontdoor(NULL, graph, "T", "Y"),
    list(feasible = FALSE, sets = list())
  )
})

test_that("named adjustment and instrument rows cross-check dagitty", {
  adjusted <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] X X -> T X -> Y T -> Y }"
  )
  adjustment_result <- identification_result(adjusted)
  adjustment <- identification_strategy(adjustment_result, "adjust")
  dagitty_adjustment <- unclass(dagitty::adjustmentSets(
    adjusted, exposure = "T", outcome = "Y", type = "canonical"
  ))

  expect_length(dagitty_adjustment, 1L)
  expect_setequal(as.character(dagitty_adjustment[[1L]]), "X")
  expect_true(adjustment$feasible)
  expect_match(adjustment$variable_set, "X")

  instrumented <- dagitty::dagitty(
    "dag { Z T [exposure] Y [outcome] Z -> T T -> Y T <-> Y }"
  )
  iv_result <- identification_result(instrumented)
  iv <- identification_strategy(iv_result, "instrument")
  dagitty_iv <- dagitty::instrumentalVariables(
    instrumented, exposure = "T", outcome = "Y"
  )

  expect_gt(length(dagitty_iv), 0L)
  expect_true(iv$feasible)
  expect_match(iv$variable_set, "Z")
  expect_false(attr(iv_result, "observationally_identifiable"))
  expect_false(attr(iv_result, "identifiable"))
})


test_that("direct identification rows describe the intervention without an arrow", {
  graphs <- c(
    identified = "dag { T [exposure] Y [outcome] M T -> M M -> Y }",
    unidentified = paste(
      "dag { T [exposure] Y [outcome] M U [latent]",
      "T -> M M -> Y U -> M U -> Y }"
    )
  )
  for( name in names(graphs) ){
    dag <- dagitty::dagitty(graphs[[name]])
    result <- identification_result(dag, effect = "direct")
    expected <- identical(name, "identified")
    expect_identical(result$feasible[1L], expected)
    expect_identical(attr(result, "identifiable"), expected)
    expect_identical(attr(result, "observationally_identifiable"), expected)
    expect_identical(result$variable_set[1L], "fixing set: {M}")
    if( expected ){
      expect_null(attr(result, "hedge"))
    }else{
      projected <- causaliflower:::id_as_dagitty(
        causaliflower:::id_admg(dag)
      )
      expect_valid_identification_hedge(result, projected)
    }
  }
})


test_that("direct fixing sets contain projected outcome parents", {
  observed <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] M X D",
    "T -> M M -> Y T -> Y X -> T X -> Y D -> Y }"
  ))
  observed_result <- identification_result(observed, effect = "direct")
  expect_identical(observed_result$variable_set[1L], "fixing set: {D, M, X}")
  expect_true(observed_result$feasible[1L])

  for( arrow in c("", "T -> Y") ){
    latent <- dagitty::dagitty(paste(
      "dag { T [exposure] Y [outcome] L [latent] T -> L L -> Y",
      arrow, "}"
    ))
    latent_result <- identification_result(latent, effect = "direct")
    expect_identical(latent_result$variable_set[1L], "fixing set: {}")
    expect_true(latent_result$feasible[1L])
    expect_true(attr(latent_result, "identifiable"))
    expect_equal(nrow(attr(latent_result, "roles")), 0L)
    expect_match(latent_result$note[1L], "latent-projected graph", fixed = TRUE)
  }

  projected <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] L [latent] Z",
    "T -> L Z -> L L -> Y }"
  ))
  projected_result <- identification_result(projected, effect = "direct")
  expect_identical(projected_result$variable_set[1L], "fixing set: {Z}")
  expect_true(projected_result$feasible[1L])

  confounded <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] M U [latent]",
    "T -> M M -> Y T -> Y U -> M U -> Y }"
  ))
  confounded_result <- identification_result(confounded, effect = "direct")
  expect_identical(confounded_result$variable_set[1L], "fixing set: {M}")
  expect_false(confounded_result$feasible[1L])
  expect_false(attr(confounded_result, "identifiable"))
})

test_that("front-door pruning stops when no candidate set can qualify", {
  # Each Mi is confounded with Y, so no front-door set exists; the old search
  # tried all 2^14 - 1 mediator subsets.
  index <- seq_len(14L)
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U1 [latent] U2 [latent]",
    "U1 -> T U1 -> Y U2 -> Y",
    paste(sprintf("T -> M%d M%d -> Y U2 -> M%d", index, index, index),
          collapse = " "),
    "}"
  ))
  G <- id_admg(dag)
  strategy_dag <- id_as_dagitty(G)
  calls <- 0L
  induced <- id_induced
  testthat::local_mocked_bindings(id_induced = function(...) {
    calls <<- calls + 1L
    induced(...)
  })
  expect_identical(
    id_frontdoor(strategy_dag, G, "T", "Y"),
    list(feasible = FALSE, sets = list())
  )
  expect_lte(calls, 4L * length(index))
  expect_false(
    identification_strategy(identification_result(dag), "front-door")$feasible
  )
})

test_that("front-door pruning keeps every minimum front-door set", {
  brute_force <- function(dag, G, treatment, outcome) {
    candidates <- Filter(function(node) {
      id_directed_reachable(G, treatment, node) &&
        id_directed_reachable(G, node, outcome)
    }, setdiff(G$V, c(treatment, outcome)))
    for (size in seq_along(candidates)) {
      sets <- Filter(function(mediators) {
        !id_directed_reachable(id_induced(G, setdiff(G$V, mediators)),
                               treatment, outcome) &&
          isTRUE(dagitty::isAdjustmentSet(dag, character(),
                                          exposure = treatment,
                                          outcome = mediators)) &&
          isTRUE(dagitty::isAdjustmentSet(dag, treatment,
                                          exposure = mediators,
                                          outcome = outcome))
      }, utils::combn(candidates, size, simplify = FALSE))
      if (length(sets) > 0L) {
        return(lapply(sets, sort))
      }
    }
    list()
  }
  set.seed(20261006)
  checked <- 0L
  while (checked < 60L) {
    nodes <- paste0("V", seq_len(7L))
    pairs <- utils::combn(7L, 2L)
    directed <- pairs[, stats::runif(ncol(pairs)) < 0.4, drop = FALSE]
    lines <- sprintf("%s -> %s", nodes[directed[1L, ]], nodes[directed[2L, ]])
    bidirected <- pairs[, stats::runif(ncol(pairs)) < 0.12, drop = FALSE]
    lines <- c(lines, sprintf("%s <-> %s", nodes[bidirected[1L, ]],
                              nodes[bidirected[2L, ]]))
    dag <- dagitty::dagitty(paste("dag {", paste(c(nodes, lines),
                                                 collapse = " "), "}"))
    G <- id_admg(dag)
    if (any(G$D$from == "V1" & G$D$to == "V7") ||
        !id_directed_reachable(G, "V1", "V7")) {
      next
    }
    checked <- checked + 1L
    strategy_dag <- id_as_dagitty(G)
    expect_identical(
      id_frontdoor(strategy_dag, G, "V1", "V7")$sets,
      brute_force(strategy_dag, G, "V1", "V7"),
      info = as.character(dag)
    )
  }
})
