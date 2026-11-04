test_that("labels work for edgeless graphs and remain unique", {
  empty <- string_to_dag("")
  isolated <- string_to_dag(c("soil moisture", "solar mass", "soil matter"))

  expect_identical(get_labels(empty), stats::setNames(character(), character()))
  expect_identical(get_labels(string_to_dag("A")), c(A = "A"))
  expect_false(anyDuplicated(unname(get_labels(isolated, "initials"))) > 0L)
  expect_false(anyDuplicated(unname(get_labels(isolated, "short"))) > 0L)
})

test_that("instrument helpers have stable zero-length outputs", {
  dag <- string_to_dag("A -> B")

  expect_identical(
    causaliflower:::extract_instrumental_variables(
      dag, character(), character()
    ),
    character()
  )
  instruments <- causaliflower:::instrumental_variables_helper(
    dag, character(), character()
  )
  expect_s3_class(instruments, "data.table")
  expect_identical(names(instruments),
                   c("treatment", "outcome", "instrument",
                     "conditioning_set"))
  expect_equal(nrow(instruments), 0L)

  edges <- causaliflower:::draw_iv_edges(
    "saturated", "IV", character(), "->"
  )
  expect_identical(names(edges), c("ancestor", "edge", "descendant"))
  expect_equal(nrow(edges), 0L)
})

test_that("minimal_sets calls the current dagitty interface", {
  dag <- dagitty::dagitty(
    "dag { X -> T; X -> Y; T -> Y; T [exposure]; Y [outcome] }"
  )

  expect_identical(unname(minimal_sets(dag)[[1L]]), "X")
  expect_identical(
    unname(minimal_sets(dag, treatment = "T", outcome = "Y")[[1L]]),
    "X"
  )
})

test_that("isolated nodes are retained by structural helpers", {
  dag <- string_to_dag(c("A -> B", "isolated"))
  expect_true("isolated" %in% names(dag))
  expect_true("isolated" %in% names(get_nodes(dag)))

  kept <- keep_edges("A -> B", dag)
  expect_true("isolated" %in% names(kept))
  expect_identical(get_roles(string_to_dag("isolated"))$observed,
                   "isolated")
})


test_that("node statuses are not read from inside quoted node names", {
  bracketed <- build_graph(variables = "note A [latent] here", treatments = "A",
                           outcomes = "B", type = "saturated")
  expect_identical(unobserved(bracketed), character(0))
  expect_setequal(observed(bracketed), c("note A [latent] here", "A", "B"))

  quoted <- string_to_dag(c('"U [latent]" -> Z', "Z -> W"))
  expect_identical(unobserved(quoted), character(0))
  expect_setequal(names(get_structure(quoted)), c("U [latent]", "Z", "W"))
  expect_no_error(get_roles(quoted))

  declared <- dagitty::dagitty('dag { "soil N" [latent] X -> "soil N" }')
  expect_identical(unobserved(declared), "soil N")
})
