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
    "dag { Z -> X; Z -> Y; X -> Y; X [exposure]; Y [outcome] }"
  )

  expect_identical(unname(minimal_sets(dag)[[1L]]), "Z")
  expect_identical(
    unname(minimal_sets(dag, treatment = "X", outcome = "Y")[[1L]]),
    "Z"
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
