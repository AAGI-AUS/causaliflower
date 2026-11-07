# ==============================================================================
# test-adjustment-source-order.R


# Keep complete enumeration and finite witnesses in recorded source order
# ==============================================================================


test_that("public node names preserve the order of complete separators and first witnesses", {
  source <- dagitty::dagitty(paste(
    'dag { T [exposure,pos="0,0"]; Y [outcome,pos="1,0"];',
    'T -> Y [pos="0.5,0.2"]; A1 -> T; A1 -> B1 -> Y;',
    'A2 -> T; A2 -> B2 -> Y; A3 -> T; A3 -> B3 -> Y; }'))
  source_nodes <- names(source)
  public_nodes <- paste0("public_", sprintf("%03d", rev(seq_along(source_nodes))))
  rename_map <- stats::setNames(public_nodes, source_nodes)
  attr(source, "labels") <- stats::setNames(paste("Label", source_nodes), source_nodes)
  graph <- copy_rename_native_nodes(source, rename_map, "source-order test")
  ranks <- data.frame(node_id = public_nodes, base_order = seq_along(public_nodes),
    at_order = seq_along(public_nodes), post_order = seq_along(public_nodes),
    pre_order = seq_along(public_nodes), time_order = seq_along(public_nodes))
  attr(graph, "cf_source_order") <- list(version = 1L, ranks = ranks)
  original_graph <- graph
  public_members <- function(value) { # Translate expected members without changing their shape or attributes
    if (is.list(value)) {
      for (index in seq_along(value)) {
        value[[index]] <- public_members(value[[index]])
      }
      return(value)
    }
    if (is.character(value)) {
      renamed <- unname(rename_map[value])
      attributes(renamed) <- attributes(value)
      return(renamed)
    }
    return(value)
  }
  for (enabled in c(FALSE, TRUE)) {
    withr::local_options(causaliflower.accelerate = enabled)
    for (effect in c("direct", "total")) {
      for (limit in c(1L, Inf)) {
        expected <- public_members(.gd_adjustment_engine(source, effect, limit))
        expect_identical(.gd_adjustment_engine(graph, effect, limit), expected)
        progress <- function(...) {
          return(invisible(NULL))
        }
        expected_observed <- public_members(.gd_adjustment_engine(source, effect, limit, progress))
        expect_identical(.gd_adjustment_engine(graph, effect, limit, progress), expected_observed)
      }
    }
  }
  expect_identical(graph, original_graph)
})


test_that("mixed graphs retain their ordinary name order", {
  graph <- dagitty::dagitty("dag { X [exposure]; Y [outcome]; U -> X; U -> Y; X -> Y; }")
  original <- graph
  ranks <- data.frame(node_id = "X", base_order = 1L, at_order = 1L,
                      post_order = 1L, pre_order = 1L, time_order = 1L)
  attr(graph, "cf_source_order") <- list(version = 1L, ranks = ranks)
  for (effect in c("direct", "total")) {
    expect_identical(.gd_adjustment_engine(graph, effect, Inf),
                     .gd_adjustment_engine(original, effect, Inf))
  }
})
