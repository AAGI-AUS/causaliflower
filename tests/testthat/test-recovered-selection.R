test_that("an unobserved set is not available adjustment", {
  g <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U [latent]; U -> T; U -> Y; T -> Y }")
  expect_false(.gd_validate_sets(g, list("U"), "total"))
  expect_false(.gd_validate_sets(g, list("U"), "direct"))
})

test_that("a set with an endpoint, a latent node or an unknown name is not valid", {
  g <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; L [latent]; A -> T; A -> Y; T -> Y; L -> Y }")
  for (effect in c("total", "direct")) {
    expect_identical(.gd_validate_sets(g, list("A", c("A", "T"), c("A", "Y"), c("A", "L"), c("A", "Z")), effect),
      c(TRUE, FALSE, FALSE, FALSE, FALSE), info = effect)
  }
})
