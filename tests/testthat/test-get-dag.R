test_that("concept resolution counts definitions, not equation inputs", {
  definitions <- list(
    Bundle = list(members = c("A", "B"), reason = "Two explicitly grouped quantities", source = "test definition"),
    Unknown = list(members = character(), reason = "Meaning unresolved", source = "test definition"),
    Loop = list(members = c("Loop", "A"), reason = "Repeated unresolved grouping", source = "test definition")
  )
  registry <- .gd_components(definitions)
  plain <- .gd_resolve_components(c("A", "B"), registry)
  expanded <- .gd_resolve_components(c("Bundle", "A"), registry)
  expect_identical(plain$members, expanded$members)
  expect_length(plain$unresolved, 0L)
  loop <- .gd_resolve_components("Loop", registry)
  expect_identical(loop$members, "A")
  expect_identical(loop$unresolved, "unresolved:Loop")
  expect_true("cycle" %in% loop$trace$status)
  expect_error(.gd_components(list(X = list(members = "A"))), "reason and a source")
})

test_that("condensed alternatives retain every temporal membership", {
  g <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A_t0; A_t1; B_t0; U_t0; T -> Y }")
  original <- serialize(g, NULL)
  raw <- list(c("A_t0", "A_t1"), "A_t0", c("B_t0", "U_t0"))
  definition <- list(U = list(members = character(), reason = "Unknown process", source = "test"))
  result <- .gd_candidates(raw, g, definition)
  expect_equal(nrow(result), 2L)
  expect_equal(result$n_concepts, c(1L, 2L))
  expect_length(result$raw_sets[[1]], 2L)
  expect_setequal(result$raw_sets[[1]][[1]], c("A_t0", "A_t1"))
  expect_equal(result$status[2], "needs definition")
  expect_true(is.na(result$n_resolved[2]))
  expect_identical(serialize(g, NULL), original)
  empty <- .gd_candidates(list(character()), g, list())
  absent <- .gd_candidates(list(), g, list())
  expect_equal(nrow(empty), 1L)
  expect_equal(empty$n_concepts, 0L)
  expect_equal(nrow(absent), 0L)
})

test_that("bounded searches match dagitty and distinguish no set from the empty set", {
  fixtures <- c(
    "dag { T [exposure]; Y [outcome]; T -> Y }",
    "dag { T [exposure]; Y [outcome]; U [latent]; U -> T; U -> Y; T -> Y }",
    "dag { T [exposure]; Y [outcome]; Z -> T; Z -> Y; T -> M -> Y; T -> Y }",
    "dag { T [exposure]; Y [outcome]; T -> M -> Y; M <-> Y; T -> Y }"
  )
  keys <- function(sets) sort(vapply(sets, function(x) paste(sort(as.character(unlist(x, use.names = FALSE))), collapse = "|"), character(1)))
  for (text in fixtures) for (effect in c("total", "direct")) {
    g <- dagitty::dagitty(text)
    ours <- .gd_adjustment_sets(g, effect, 10L)
    reference <- dagitty::adjustmentSets(g, effect = effect)
    expect_identical(keys(ours$sets), unname(keys(reference)))
    expect_true(ours$complete)
  }
  g <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A -> T; A -> B -> Y; C -> T; C -> D -> Y; T -> Y }")
  for (effect in c("total", "direct")) {
    bounded <- .gd_adjustment_sets(g, effect, 1L)
    expect_length(bounded$sets, 1L)
    expect_false(bounded$complete)
  }
})
