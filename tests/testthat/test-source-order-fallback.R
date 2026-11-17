# ==============================================================================
# test-source-order-fallback.R


# Leave whole graph and catalogue order unchanged when recorded ranks are incomplete
# ==============================================================================


test_that("unknown endpoints prevent partial reordering in ranked search", {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome]; T -> Y;",
    "A -> T; A -> B -> Y; C -> T; C -> D -> Y; }"))
  prepared <- .cd_prepare(graph, NULL, NULL)
  expected <- .cd_ranked_search(prepared, "total", steps = 5000L)
  source_order <- list(version = 1L, ranks = data.frame(node_id = c("A", "B", "C", "D"),
    base_order = 4:1, at_order = 4:1, post_order = 4:1, pre_order = 4:1, time_order = 4:1))
  prepared$record$source_order <- source_order
  attr(prepared$selection, "cf_source_order") <- source_order
  actual <- .cd_ranked_search(prepared, "total", steps = 5000L)
  expect_identical(actual, expected)
})
