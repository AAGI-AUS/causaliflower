status_setter_graph <- function() {
  dag <- dagitty::dagitty(paste0(
    "dag { bb=\"0,0,1,1\" ",
    "U [latent,pos=\"0.10,0.10\"] T [exposure,pos=\"0.30,0.50\"] ",
    "Y [outcome,adjusted,pos=\"0.70,0.50\"] Z [pos=\"0.50,0.10\"] ",
    "T -> Y [pos=\"0.50,0.30\"] Z -> T U -> Z }"
  ))
  attr(dag, "note") <- "keep"
  dag
}


status_setter_native <- function() {
  structure(
    paste0(
      "dag rich {\n",
      "bb=\"0,0,1,1\"\n",
      "A [label=\"Alpha\",pos=\"0.123456789,0.5\"]\n",
      "B [outcome,adjusted,note=\"a note\",pos=\"0.9,0.25\"]\n",
      "C [selected,pos=\"0.4,0.66666666\"]\n",
      "D\n",
      "A -> B [pos=\"0.5,0.4\"]\n",
      "B -> C\n",
      "C -> D\n}"
    ),
    class = "dagitty"
  )
}


status_setters <- list(
  exposure = set_treatment,
  outcome = set_outcome,
  latent = set_latent
)


test_that("each setter adds only its own status", {
  dag <- status_setter_graph()

  treated <- set_treatment(dag, "Z")
  expect_setequal(treatments(treated), c("T", "Z"))
  expect_identical(outcomes(treated), "Y")
  expect_identical(unobserved(treated), "U")

  outcome <- set_outcome(dag, "Z")
  expect_setequal(outcomes(outcome), c("Y", "Z"))
  expect_identical(treatments(outcome), "T")
  expect_identical(unobserved(outcome), "U")

  latent <- set_latent(dag, "Z")
  expect_setequal(unobserved(latent), c("U", "Z"))
  expect_identical(treatments(latent), "T")
  expect_identical(outcomes(latent), "Y")
})


test_that("setters leave edges, coordinates, curves and attributes alone", {
  dag <- status_setter_graph()

  for (setter in list(set_treatment, set_outcome, set_latent)) {
    updated <- setter(dag, "Z")
    expect_identical(edge_keys(updated), edge_keys(dag))
    expect_identical(dagitty::coordinates(updated), dagitty::coordinates(dag))
    expect_identical(dagitty::adjustedNodes(updated), "Y")
    expect_identical(attr(updated, "note"), "keep")
    expect_true(grepl("bb=\"0,0,1,1\"", as.character(updated), fixed = TRUE))
    expect_true(grepl("T -> Y [pos=\"0.500,0.300\"]", as.character(updated),
                      fixed = TRUE))
  }
})


test_that("a setter inserts the keyword and no other byte on a built graph", {
  set.seed(1)
  dag <- build_graph(
    variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
    mediators = "M", type = "saturated"
  )
  dag <- add_coords(dag)
  dag <- add_curves(dag, height = 0.25)
  before <- as.character(dag)
  expect_match(before, "pos=\"[^\"]*\\.[0-9]{4}")

  for (keyword in names(status_setters)) {
    updated <- as.character(status_setters[[keyword]](dag, "X1"))
    expect_identical(
      sub(paste0("X1 [", keyword, ",pos"), "X1 [pos", updated, fixed = TRUE),
      before
    )
  }
})


test_that("a setter keeps a native graph's name, labels, notes and precision", {
  dag <- status_setter_native()
  before <- as.character(dag)

  for (keyword in names(status_setters)) {
    updated <- status_setters[[keyword]](dag, c("A", "D"))
    stripped <- sub(paste0("A [", keyword, ",label"), "A [label",
                    as.character(updated), fixed = TRUE)
    stripped <- sub(paste0("D [", keyword, "]"), "D", stripped, fixed = TRUE)
    expect_identical(stripped, before)
    expect_identical(dagitty::adjustedNodes(updated), "B")
    expect_identical(causaliflower:::selected_nodes(updated), "C")
    expect_identical(edge_keys(updated), edge_keys(dag))
  }

  expect_setequal(
    causaliflower:::get_outcome_names(set_outcome(dag, c("A", "D"))),
    c("A", "B", "D")
  )
  expect_identical(
    dagitty::coordinates(set_treatment(dag, "A"))$x[["A"]],
    0.123456789
  )
})


test_that("a node used only in an edge statement gains a declaration", {
  dag <- structure("dag { X -> Y\n}", class = "dagitty")

  updated <- set_latent(dag, "Y")
  expect_identical(unobserved(updated), "Y")
  expect_identical(as.character(updated), "dag { X -> Y\n\n\"Y\" [latent]\n}")
  expect_identical(edge_keys(updated), "X -> Y")
  expect_identical(names(updated), c("X", "Y"))
})


test_that("setting a status a node already carries is a no-op", {
  dag <- status_setter_graph()

  expect_identical(as.character(set_treatment(dag, "T")), as.character(dag))
  expect_identical(as.character(set_outcome(dag, "Y")), as.character(dag))
  expect_identical(as.character(set_latent(dag, "U")), as.character(dag))

  twice <- set_treatment(set_treatment(dag, "Z"), "Z")
  expect_identical(as.character(twice), as.character(set_treatment(dag, "Z")))
})


test_that("setters reject a non-graph and unknown node names", {
  dag <- status_setter_graph()

  expect_error(
    set_treatment("dag { T }", "T"),
    "set_treatment(): 'dag' must be a dagitty object.",
    fixed = TRUE
  )
  expect_error(
    set_outcome("dag { T }", "T"),
    "set_outcome(): 'dag' must be a dagitty object.",
    fixed = TRUE
  )
  expect_error(
    set_latent("dag { T }", "T"),
    "set_latent(): 'dag' must be a dagitty object.",
    fixed = TRUE
  )

  expect_error(
    set_treatment(dag, c("Q", "R")),
    "set_treatment(): 'nodes' must name nodes in the graph: Q, R.",
    fixed = TRUE
  )
  expect_error(
    set_outcome(dag, c("Q", "R")),
    "set_outcome(): 'nodes' must name nodes in the graph: Q, R.",
    fixed = TRUE
  )
  expect_error(
    set_latent(dag, c("Q", "R")),
    "set_latent(): 'nodes' must name nodes in the graph: Q, R.",
    fixed = TRUE
  )

  expect_error(
    set_treatment(dag, 42),
    "set_treatment(): 'nodes' must name nodes in the graph: 42.",
    fixed = TRUE
  )
})
