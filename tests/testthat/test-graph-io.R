test_that("string_to_dag parses chains and every supported operator", {
  chain <- string_to_dag("A -> B -> C")
  mixed <- string_to_dag("A <- B; C <-> D; E -- F")

  expect_identical(edge_keys(chain), c("A -> B", "B -> C"))
  expect_true("B -> A" %in% edge_keys(mixed))
  expect_true("C <-> D" %in% edge_keys(mixed))
  expect_true("E -- F" %in% edge_keys(mixed))
})

test_that("string_to_dag handles special and isolated names", {
  cases <- c(
    '"soil moisture" -> "air temp"',
    "soil-moisture -> air-temp",
    "rate% -> yield%"
  )
  expected <- list(
    c("soil moisture", "air temp"),
    c("soil-moisture", "air-temp"),
    c("rate%", "yield%")
  )

  for (i in seq_along(cases)) {
    expect_setequal(names(string_to_dag(cases[i])), expected[[i]])
  }
  expect_identical(names(string_to_dag("isolated node")), "isolated node")
  expect_length(names(string_to_dag("")), 0L)

  operator_name <- string_to_dag('"A->B" -> C')
  separator_name <- string_to_dag('"A;B" -> C')
  expect_setequal(names(operator_name), c("A->B", "C"))
  expect_identical(edge_keys(operator_name), "A->B -> C")
  expect_setequal(names(separator_name), c("A;B", "C"))
  expect_identical(edge_keys(separator_name), "A;B -> C")
  expect_identical(names(string_to_dag('"A->B"')), "A->B")
})

test_that("string parsing rejects malformed edge text", {
  expect_error(string_to_dag("A ->"), "node names cannot be empty")
  expect_error(string_to_dag("A => B"), "could not parse")
  expect_error(string_to_dag('"A -> B'), "unmatched quote")
})

test_that("keep_edges accepts strings, tables, dagitty objects, and reverse arrows", {
  dag <- build_graph(c("Z1", "Z2"), "X", "Y", type = "saturated")

  from_strings <- keep_edges(c("Z1 -> X", "X -> Y"), dag)
  from_reverse <- keep_edges(c("X <- Z1", "X -> Y"), dag)
  from_table <- keep_edges(
    data.frame(v = c("Z1", "X"), e = "->", w = c("X", "Y")),
    dag
  )
  from_dag <- keep_edges(from_strings, dag)

  expected <- sort(c("Z1 -> X", "X -> Y"))
  expect_identical(edge_keys(from_strings), expected)
  expect_identical(edge_keys(from_reverse), expected)
  expect_identical(edge_keys(from_table), expected)
  expect_identical(edge_keys(from_dag), expected)
})

test_that("keep_edges preserves nodes, roles, coordinates, and symmetric edges", {
  dag <- dagitty::dagitty(
    'dag { "soil moisture" [exposure] "yield%" [outcome]
           "air-temp" -> "soil moisture"
           "soil moisture" -> "yield%" A <-> B }'
  )
  dagitty::coordinates(dag) <- list(
    x = stats::setNames(seq_along(names(dag)), names(dag)),
    y = stats::setNames(rep(0, length(names(dag))), names(dag))
  )
  original_coordinates <- dagitty::coordinates(dag)

  kept <- keep_edges(
    c('"air-temp" -> "soil moisture"', '"soil moisture" -> "yield%"',
      "B <-> A"),
    dag
  )

  expect_setequal(names(kept), names(dag))
  expect_identical(treatments(kept), "soil moisture")
  expect_identical(causaliflower:::.outcomes(kept), "yield%")
  expect_equal(dagitty::coordinates(kept), original_coordinates)
  expect_true("A <-> B" %in% edge_keys(kept))
})

test_that("keep_edges has deliberate empty and missing-edge contracts", {
  dag <- build_graph(c("Z1", "Z2"), "X", "Y", type = "saturated")

  expect_error(keep_edges(character(), dag), "no edges were supplied")
  expect_error(keep_edges("A => B", dag), "could not parse")
  expect_error(keep_edges(data.frame(a = 1, b = 2), dag),
               "v, e, and w")

  expect_warning(empty <- keep_edges("absent -> edge", dag),
                 "not present")
  expect_length(edge_keys(empty), 0L)
  expect_setequal(names(empty), names(dag))
  expect_identical(treatments(empty), "X")
  expect_identical(causaliflower:::.outcomes(empty), "Y")

  expect_warning(
    expect_error(
      keep_edges(data.frame(v = NA, e = NA, w = NA), dag),
      "no complete edges"
    ),
    "containing NA"
  )
})
