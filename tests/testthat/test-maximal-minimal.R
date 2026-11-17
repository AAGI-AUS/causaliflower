# The minimality test of the largest DAG against the full validity test.
old_minimal <- function(graph, sets, effect) {
  vapply(sets, function(set) {
    isTRUE(.gd_validate_sets(graph, list(set), effect)) &&
      (!length(set) || !any(.gd_validate_sets(graph, lapply(seq_along(set), function(i) set[-i]), effect)))
  }, logical(1))
}

counted <- function(graph) {
  validate <- .gd_validator(graph)
  calls <- 0L
  list(validate = function(sets, effect) {
    calls <<- calls + 1L
    validate(sets, effect)
  }, calls = function() calls)
}

expect_minimal <- function(graph, sets, effect, expected, validated = FALSE) {
  validator <- counted(graph)
  decided <- .cd_maximal_minimal_sets(graph, sets, effect, validator$validate)
  expect_identical(decided, old_minimal(graph, sets, effect))
  expect_identical(decided, expected)
  expect_identical(validator$calls() > 0L, validated)
}

test_that("minimal, extra, dropped and empty sets of a total effect match the full test", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> B; B -> Y; C -> T; C -> Y;
    T -> M; M -> Y; M -> D; F -> E; E -> Y; A -> W }")
  expect_minimal(graph, list(c("A", "C"), c("B", "C"), c("A", "C", "E"), c("A", "B", "C"), "A", "C", character()),
    "total", c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_true(.cd_maximal_minimal_on(graph, list(c("A", "C"), c("B", "C")), "total"))
  expect_false(.cd_maximal_minimal_on(graph, list(c("A", "C"), c("A", "C", "E")), "total"))
  expect_true(.cd_maximal_minimal_on(graph, list(), "total"))
})

test_that("the empty set is minimal exactly when it is valid", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y; A -> Y; B -> T }")
  expect_minimal(graph, list(character(), "A", "B"), "total", c(TRUE, FALSE, FALSE))
  expect_minimal(graph, list(character(), "A"), "direct", c(TRUE, FALSE))
})

test_that("members outside the ancestors or forbidden go to the full test", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> B; B -> Y; C -> T; C -> Y;
    T -> M; M -> Y; M -> D; F -> E; E -> Y; A -> W }")
  expect_minimal(graph, list(c("A", "C", "W")), "total", FALSE, validated = TRUE)
  expect_minimal(graph, list(c("A", "C", "D")), "total", FALSE, validated = TRUE)
  expect_minimal(graph, list(c("A", "C"), c("A", "C", "D")), "total", c(TRUE, FALSE), validated = TRUE)
})

test_that("direct effect sets match the full test", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> M; M -> Y; T -> Y; A -> T; A -> Y;
    U -> M; U -> Y; B -> A }")
  expect_minimal(graph, list(c("A", "M", "U"), c("B", "M", "U"), c("A", "B", "M", "U"), c("A", "M"), c("M", "U"),
    character()), "direct", c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  mediated <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> B; B -> Y; C -> T; C -> Y;
    T -> M; M -> Y; M -> D; F -> E; E -> Y; A -> W }")
  expect_minimal(mediated, list(c("A", "C", "M"), c("B", "C", "M"), c("A", "C", "M", "E"), c("A", "C")), "direct",
    c(TRUE, TRUE, FALSE, FALSE))
})

test_that("bidirected edges and several treatments match the full test", {
  bidirected <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y; A -> T; A <-> Y; B -> A; B -> Y }")
  expect_minimal(bidirected, list("A", c("A", "B"), "B", character()), "total", c(TRUE, FALSE, FALSE, FALSE))
  expect_minimal(bidirected, list("A", c("A", "B"), character()), "direct", c(TRUE, FALSE, FALSE))
  several <- dagitty::dagitty("dag { T1 [exposure]; T2 [exposure]; Y [outcome]; T1 -> Y; T2 -> Y;
    A -> T1; A -> Y; B -> T2; B -> Y; C -> Y }")
  expect_minimal(several, list(c("A", "B"), c("A", "B", "C"), "A", "B", character()), "total",
    c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_minimal(several, list(c("A", "B"), c("A", "B", "C"), "A"), "direct", c(TRUE, FALSE, FALSE))
})
