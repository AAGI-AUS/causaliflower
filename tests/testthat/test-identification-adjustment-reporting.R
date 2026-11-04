test_that("adjustment reporting separates minimal and canonical sets", {
  dag <- dagitty::dagitty(paste(
    "dag { M U T [exposure] Y [outcome] X",
    "M -> Y U -> T U -> X T -> M T -> Y X -> T X -> Y }"
  ))
  result <- identification_result(dag)
  adjustment <- identification_strategy(result, "^covariate adjustment$")

  expect_true(adjustment$feasible)
  expect_identical(adjustment$variable_set, "{X}")
  expect_identical(attr(result, "adjustment_sets"), list("X"))
  expect_identical(attr(result, "canonical_adjustment"), list(c("U", "X")))
  expect_true(attr(result, "identifiable"))
  expect_identical(
    attr(result, "roles")[role == "adjustment variable", variable],
    "X"
  )
})

test_that("the display limit preserves every minimal adjustment set", {
  edges <- vapply(seq_len(3L), function(i){
    edge <- paste0("A", i, " -> B", i, " -> T A", i, " -> Y")

    return(edge)
  }, character(1L))
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] T -> Y",
    paste(edges, collapse = " "), "}"
  ))
  result <- identification_result(dag)
  adjustment <- identification_strategy(result, "^covariate adjustment$")
  sets <- attr(result, "adjustment_sets")
  expected <- dagitty::adjustmentSets(dag, type = "minimal")
  expected <- lapply(unclass(expected), function(set){
    nodes <- sort(as.character(set))

    return(nodes)
  })

  expect_length(sets, 8L)
  expect_identical(sets, unname(expected))
  expect_identical(
    sum(strsplit(adjustment$variable_set, "", fixed = TRUE)[[1L]] == "{"),
    5L
  )
  expect_match(adjustment$variable_set, "; ... (3 more)", fixed = TRUE)
  expect_true(all(vapply(sets, function(set){
    valid <- dagitty::isAdjustmentSet(dag, Z = set)

    return(valid)
  }, logical(1L))))
  expect_setequal(
    attr(result, "roles")[role == "adjustment variable", variable],
    c("A1", "A2", "A3", "B1", "B2", "B3")
  )
  expect_false(grepl("more", causaliflower:::id_render_sets(sets), fixed = TRUE))
  expect_identical(
    sum(strsplit(causaliflower:::id_render_sets(sets), "", fixed = TRUE)[[1L]] == "{"),
    8L
  )
})

test_that("empty adjustment and unavailable adjustment remain distinct", {
  empty <- identification_result(dagitty::dagitty(
    "dag { T [exposure] Y [outcome] T -> Y }"
  ))
  unavailable <- identification_result(dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent]",
    "U -> T U -> Y T -> Y }"
  )))

  expect_identical(attr(empty, "adjustment_sets"), list(character()))
  expect_identical(attr(empty, "canonical_adjustment"), list(character()))
  expect_identical(
    identification_strategy(empty, "^covariate adjustment$")$variable_set,
    "{}"
  )
  expect_identical(attr(unavailable, "adjustment_sets"), list())
  expect_identical(attr(unavailable, "canonical_adjustment"), list())
  expect_identical(
    identification_strategy(unavailable, "^covariate adjustment$")$variable_set,
    "none"
  )
  expect_false(
    identification_strategy(unavailable, "^covariate adjustment$")$feasible
  )
})

test_that("adjustment attributes describe only assessed marginal total effects", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] X",
    "X -> T X -> Y T -> Y }"
  ))
  direct <- identification_result(dag, effect = "direct")
  conditional <- identification_result(dag, conditioning = "X")
  marginal <- identification_result(dag, conditioning = character())
  experimental <- identification_result(dag, experiments = "X")

  for(result in list(direct, conditional)){
    expect_null(attr(result, "adjustment_sets"))
    expect_null(attr(result, "canonical_adjustment"))
  }
  for(result in list(marginal, experimental)){
    expect_identical(attr(result, "adjustment_sets"), list("X"))
    expect_identical(attr(result, "canonical_adjustment"), list("X"))
  }
})

test_that("latent ancestors are not reported as adjustment variables", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] L [latent] X",
    "L -> X X -> T X -> Y T -> Y }"
  ))
  result <- identification_result(dag)

  expect_identical(attr(result, "adjustment_sets"), list("X"))
  expect_identical(attr(result, "canonical_adjustment"), list("X"))
  expect_false("L" %in% attr(result, "roles")$variable)
})
