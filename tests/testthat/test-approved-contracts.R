test_that("connect_nodes() restricts new edges to the selected nodes", {
  set.seed(1)
  dag <- suppressMessages(build_graph(c("Z1", "Z2", "Z3"), "X", "Y", type = "ordered"))

  base_keys <- edge_keys(dag)

  connected <- suppressWarnings(suppressMessages(
    connect_nodes(dag, nodes = "Z1", type = "saturated", print_edges = FALSE)
  ))

  added <- setdiff(edge_keys(connected), base_keys)
  expect_true(length(added) > 0L)
  added_ends <- strsplit(added, " ")
  expect_true(all(vapply(added_ends, function(x) "Z1" %in% x[c(1, 3)], logical(1))))
})

test_that("connect_nodes() selection default is the whole graph", {
  set.seed(1)
  dag <- suppressMessages(build_graph(c("Z1", "Z2"), "X", "Y", type = "ordered"))

  all_default <- suppressWarnings(suppressMessages(
    connect_nodes(dag, type = "full", print_edges = FALSE)
  ))
  all_explicit <- suppressWarnings(suppressMessages(
    connect_nodes(dag, nodes = names(dag), type = "full", print_edges = FALSE)
  ))

  expect_identical(edge_keys(all_default), edge_keys(all_explicit))
})

test_that("add_nodes() rejects a node_role combined with explicit targets", {
  set.seed(1)
  dag <- suppressMessages(build_graph(c("Z1", "Z2"), "X", "Y", type = "ordered"))

  expect_error(
    add_nodes(dag, "W", node_role = "confounder", ancestors = "X", print_edges = FALSE),
    "not both"
  )
  expect_error(
    add_nodes(dag, "W", node_role = "confounder", descendants = "Y", print_edges = FALSE),
    "not both"
  )
})

test_that("add_nodes() rejects position for unsupported node roles", {
  set.seed(1)
  dag <- suppressMessages(build_graph(c("Z1", "Z2"), "X", "Y", type = "ordered"))

  expect_error(
    add_nodes(dag, "K", node_role = "collider", position = first(), print_edges = FALSE),
    "supported for node_role"
  )
})

test_that("assess_edges() honours edges_to_assess = 'bidirectional' without the causal sequence", {
  dag <- string_to_dag(c("A -> B", "C <-> D", "X -> Y"))

  capture.output(
    res <- suppressMessages(assess_edges(dag, edges_to_assess = "bidirectional"))
  )

  expect_type(res, "list")
  expect_identical(nrow(res$edges_to_assess), 1L)
  expect_identical(unlist(res$edges_to_assess[, "e"]), c(e = "<->"))
  expect_setequal(unname(res$edges_to_keep), c("A -> B", "X -> Y"))
})

test_that("assess_edges() bidirectional mode assesses only symmetric edges in the sequence", {
  dag <- string_to_dag(c("A -> B", "C <-> D", "X -> Y"))
  criteria <- data.frame(
    name = "temporality", question = "q", description = "d",
    source = "s", required = "yes", stringsAsFactors = FALSE
  )

  testthat::with_mocked_bindings(
    .assessment_is_interactive = function() TRUE,
    .assessment_readline = function(prompt) "y",
    {
      capture.output(
        res <- suppressMessages(assess_edges(
          dag,
          edges_to_assess = "bidirectional",
          assess_causal_criteria = TRUE,
          save_answers = TRUE,
          causal_criteria = criteria
        ))
      )
    }
  )

  expect_identical(nrow(res$answers), 1L)
  expect_identical(unname(unlist(res$answers[, c("v", "e", "w")])), c("C", "<->", "D"))
})

test_that("minimal_sets() separates the empty set, no valid set, and invalid calls", {
  clean <- string_to_dag("X -> Y")
  dagitty::exposures(clean) <- "X"
  dagitty::outcomes(clean) <- "Y"
  expect_warning(res <- minimal_sets(clean), NA)
  expect_identical(length(res), 1L)
  expect_identical(names(res), "0")
  expect_identical(length(res[[1]]), 0L)

  latent <- dagitty::dagitty("dag { X [exposure] Y [outcome] U [latent] U -> X U -> Y X -> Y }")
  expect_message(none <- minimal_sets(latent), "assess edges using")
  expect_null(none)

  expect_error(minimal_sets("not a dag"), "Failed to compute")
})

test_that("minimal_sets() validates num_sets as a positive whole number", {
  dag <- dagitty::dagitty("dag { X [exposure] Y [outcome] Z Z -> X Z -> Y X -> Y }")

  for( bad in list(0, -1, NA, c(2, 3), 2.5, "3", Inf) ){
    expect_error(minimal_sets(dag, num_sets = bad), "positive integer")
  }
})

test_that("minimal_sets() warns once when a same-length set is excluded", {
  dag <- dagitty::dagitty("dag { X [exposure] Y [outcome] Z W Z -> X Z -> W W -> Y X -> Y }")

  expect_warning(res <- minimal_sets(dag, num_sets = 1),
                 "A same length adjustment set has been excluded")
  expect_identical(length(res), 1L)
})

test_that("ordered graphs connect every role to every outcome", {
  set.seed(1)
  dag <- suppressMessages(build_graph(
    variables = c("Z1", "Z2"), treatments = "X", outcomes = c("Y1", "Y2"),
    mediators = "M", latent_variables = "U", type = "ordered"
  ))

  keys <- edge_keys(dag)
  expect_true("X -> Y2" %in% keys)
  expect_true("M -> Y2" %in% keys)
  expect_true("U -> Y2" %in% keys)
  expect_true("Z1 -> Y2" %in% keys)
  expect_true("Y1 -> Y2" %in% keys)
})

test_that("join_graphs() warns on conflicting shared declarations and cyclic output", {
  dag1 <- dagitty::dagitty("dag { X [exposure] Y [outcome] X -> Y }")
  dag2 <- dagitty::dagitty("dag { X [outcome] Y [exposure] Y -> X }")

  expect_warning(
    expect_warning(joined <- join_graphs(dag1, dag2), "disagree on the exposure/outcome/latent status"),
    "contains cycles"
  )
  expect_identical(treatments(joined), "X")
  expect_identical(causaliflower:::.outcomes(joined), "Y")
  expect_false(dagitty::isAcyclic(joined))
})

test_that("join_graphs() warns when direction disagreement creates a cycle", {
  dag1 <- string_to_dag("A -> B")
  dag2 <- string_to_dag("B -> A")

  expect_warning(joined <- join_graphs(dag1, dag2), "contains cycles")
  expect_false(dagitty::isAcyclic(joined))
})

test_that("status accessors match dagitty's own readers", {
  dag <- dagitty::dagitty("dag { X [exposure] Y [outcome] U [latent] Z X -> Y U -> Y }")
  expect_setequal(treatments(dag), dagitty::exposures(dag))
  expect_setequal(causaliflower:::.outcomes(dag), dagitty::outcomes(dag))
  expect_setequal(unobserved(dag), dagitty::latents(dag))
  expect_setequal(observed(dag), setdiff(names(dag), dagitty::latents(dag)))

  short <- dagitty::dagitty("dag { X [e] Y [o] U [u] X -> Y }")
  expect_setequal(treatments(short), dagitty::exposures(short))
  expect_setequal(causaliflower:::.outcomes(short), dagitty::outcomes(short))
  expect_setequal(unobserved(short), dagitty::latents(short))

  dup <- dagitty::dagitty("dag { X [exposure] Y [outcome] X [pos=\"1,1\"] X -> Y }")
  expect_identical(treatments(dup), "X")
  expect_setequal(treatments(dup), dagitty::exposures(dup))

  upper <- dagitty::dagitty("dag { X [Exposure] Y [outcome] X -> Y }")
  expect_setequal(treatments(upper), dagitty::exposures(upper))

  empty <- string_to_dag("")
  expect_identical(treatments(empty), character(0))
  expect_identical(observed(empty), character(0))
})

test_that("outcomes() pairs with treatments() and keeps the dagitty replacement form", {
  dag <- dagitty::dagitty("dag { X [exposure] Y [outcome] X -> Y }")
  expect_identical(outcomes(dag), "Y")
  expect_setequal(outcomes(dag), dagitty::outcomes(dag))

  g <- string_to_dag("A -> B")
  outcomes(g) <- "B"
  expect_identical(outcomes(g), "B")
  expect_setequal(dagitty::outcomes(g), "B")
})

test_that("status accessors reject unsupported graph types", {
  pag <- dagitty::dagitty("pag { A @-> B }")

  expect_error(treatments(pag), "does not currently support")
  expect_error(observed(pag), "does not currently support")
  expect_error(unobserved(pag), "does not currently support")
})

test_that("add_coords() is deterministic and temporally ordered", {
  dag <- suppressMessages(build_graph(c("Z1", "Z2"), "X", "Y",
                                      mediators = "M", instrumental_variables = "IV",
                                      type = "saturated"))
  a <- add_coords(dag)
  set.seed(42)
  b <- add_coords(dag)
  expect_identical(dagitty::coordinates(a), dagitty::coordinates(b))

  co <- dagitty::coordinates(a)
  expect_false(anyNA(unlist(co)))

  ed <- dagitty::edges(a)
  dir <- ed[ed$e == "->", ]
  expect_true(all(co$x[as.character(dir$v)] < co$x[as.character(dir$w)]))

  for( col in split(names(co$x), co$x) ){
    if( length(col) > 1 ){
      ys <- sort(co$y[col])
      expect_true(all(diff(ys) >= 1 - 1e-9))
    }
  }
})

test_that("join_graphs() keeps the first graph's coordinates and places only new nodes", {
  dag <- add_coords(suppressMessages(build_graph(c("Z1", "Z2"), "X", "Y", type = "ordered")))
  before <- dagitty::coordinates(dag)

  joined <- suppressWarnings(join_graphs(dag, string_to_dag("Y -> K")))
  after <- dagitty::coordinates(joined)

  common <- names(before$x)
  expect_identical(after$x[common], before$x[common])
  expect_identical(after$y[common], before$y[common])
  expect_false(anyNA(c(after$x["K"], after$y["K"])))
})

test_that("cyclic graphs are laid out on an acyclic subset with a warning", {
  cyc <- string_to_dag(c("A -> B", "B -> C", "C -> A"))

  expect_warning(laid <- add_coords(cyc), "acyclic subset of the supplied")
  co <- dagitty::coordinates(laid)
  expect_false(anyNA(unlist(co)))
  expect_setequal(edge_keys(laid), c("A -> B", "B -> C", "C -> A"))
})

test_that("add_coords() validates its layout controls", {
  dag <- string_to_dag("A -> B")

  expect_error(add_coords(dag, x_step = 0), "positive number")
  expect_error(add_coords(dag, y_step = -1), "positive number")
  expect_error(add_coords(dag, layout_shape = 2), "between 0 and 1")
})

test_that("symmetric edges across ranks share a height after layout", {
  g <- dagitty::dagitty("dag { A -> B A <-> C B -> C }")
  laid <- add_coords(g)
  co <- dagitty::coordinates(laid)
  expect_identical(co$y[["A"]], co$y[["C"]])
})

test_that("connect_nodes() position requires a node_role", {
  set.seed(1)
  dag <- suppressMessages(build_graph(c("Z1", "Z2"), "X", "Y", type = "ordered"))

  expect_error(
    connect_nodes(dag, nodes = "Z1", position = first(), print_edges = FALSE),
    "requires a 'node_role'"
  )
})
