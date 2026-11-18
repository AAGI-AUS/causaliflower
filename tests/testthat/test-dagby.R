test_that("condensed alternatives retain every temporal membership", {
  g <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A_t0; A_t1; B_t0; U_t0; T -> Y }")
  nodes <- names(g)
  ids <- c(A_t0 = "A", A_t1 = "A", B_t0 = "B", U_t0 = "U", T = "T", Y = "Y")
  attr(g, "cf_node_map") <- data.frame(node = nodes, quantity = unname(ids[nodes]),
    member_id = unname(ids[nodes]), occurrence = nodes)
  original <- serialize(g, NULL)
  raw <- list(c("A_t0", "A_t1"), "A_t0", c("B_t0", "U_t0"))
  result <- .gd_candidates(raw, g)
  expect_equal(nrow(result), 2L)
  expect_equal(result$n_concepts, c(1L, 2L))
  expect_length(result$raw_sets[[1]], 2L)
  expect_setequal(result$raw_sets[[1]][[1]], c("A_t0", "A_t1"))
  expect_identical(serialize(g, NULL), original)
  empty <- .gd_candidates(list(character()), g)
  absent <- .gd_candidates(list(), g)
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
