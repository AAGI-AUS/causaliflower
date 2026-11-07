# ==============================================================================
# test-source-order-fallback.R


# Leave whole graph and catalogue order unchanged when recorded ranks are incomplete
# ==============================================================================


test_that("an unmapped family retains the original tie order by display label", {
  table <- data.frame(id = c("first", "second"), group = c(1L, 2L), option = c(1L, 1L))
  table$resolved <- list(c("A", "Z"), c("A !", "B"))
  table$raw_sets <- list(list(c("A", "Z")), list(c("A !", "B")))
  table$concepts <- table$resolved
  source_order <- list(version = 1L, ranks = data.frame(node_id = "unrelated",
    base_order = 1L, at_order = 1L, post_order = 1L, pre_order = 1L, time_order = 1L))
  expected <- .gd_set_family(table, "direct")
  expect_identical(.gd_set_family(table, "direct", source_order), expected)
  expect_identical(expected$table$label, c("A ! + B", "A + Z"))
  source_order$ranks <- data.frame(node_id = c("A", "Z"), base_order = 2:1,
    at_order = 2:1, post_order = 2:1, pre_order = 2:1, time_order = 2:1)
  expect_identical(.gd_set_family(table, "direct", source_order), expected)
})


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
