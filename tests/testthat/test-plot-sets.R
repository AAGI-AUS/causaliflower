################################################################################
# test-plot-sets.R


# Preserve complete adjustment alternatives, plot controls and the supplied graph

plot_set_nodes <- function(plot){ # Read the named groups from the plotted point layer
  layer <- Filter(function(layer){ return(inherits(layer$geom, "GeomPoint")) }, plot$layers)
  nodes <- layer[[1L]]$data
  return(stats::setNames(as.character(nodes$role), as.character(nodes$name)))
}


test_that("plot_sets highlights adjustment membership without changing role plots", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure]; Y [outcome]; U [latent];",
    "STATE -> LOCATION; STATE -> U; LOCATION -> T; LOCATION -> Y;",
    "U -> T; U -> VARIETY; VARIETY -> Y; T -> Y }"
  ))
  incoming <- serialize(dag, NULL)
  expect_no_warning(plot <- plot_sets(dag))
  nodes <- plot_set_nodes(plot)
  expect_setequal(names(nodes)[nodes == "adjustment set"],
                  c("LOCATION", "VARIETY"))
  expect_identical(unname(nodes[c("T", "Y", "U", "STATE")]),
                    c("treatment", "outcome", "latent", "observed"))
  built <- ggplot2::ggplot_build(plot)$plot
  expect_identical(built$scales$get_scales("colour")$map("adjustment set"),
                    "coral2")
  expect_identical(built$scales$get_scales("fill")$map("latent"), "white")

  ordinary <- plot_set_nodes(plot_dagitty(dag))
  expect_identical(unname(ordinary["STATE"]), "confounder")
  expect_false("adjustment set" %in% ordinary)
  expect_identical(serialize(dag, NULL), incoming)
})


test_that("valid empty adjustment is distinguished from unavailable adjustment", {
  empty <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y }")
  unavailable <- dagitty::dagitty(paste(
    "dag { T [exposure]; Y [outcome]; U [latent];",
    "U -> T; U -> Y; T -> Y }"
  ))
  expect_no_warning(empty_plot <- plot_sets(empty))
  expect_false("adjustment set" %in% plot_set_nodes(empty_plot))
  expect_warning(unavailable_plot <- plot_sets(unavailable),
                  "no adjustment set is available for the total effect")
  expect_false("adjustment set" %in% plot_set_nodes(unavailable_plot))
  expect_identical(unname(plot_set_nodes(unavailable_plot)["U"]), "latent")
})


test_that("multiple sets retain minimal_sets ordering and are not truncated", {
  edges <- vapply(seq_len(3L), function(i){
    return(paste0("A", i, " -> T; A", i, " -> B", i, " -> Y;"))
  }, character(1L))
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure]; Y [outcome]; T -> Y;",
    paste(edges, collapse = " "), "}"
  ))
  sets <- minimal_sets(dag, num_sets = .Machine$integer.max)
  expect_length(sets, 8L)
  expect_warning(plot <- plot_sets(dag),
                  "8 adjustment sets are available.*first in minimal_sets")
  nodes <- plot_set_nodes(plot)
  expect_setequal(names(nodes)[nodes == "adjustment set"],
                  as.character(unlist(sets[[1L]], use.names = FALSE)))
})


test_that("effect selection changes the highlighted set to the single door set", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure]; Y [outcome]; U -> T; U -> Y;",
    "T -> M -> Y; T -> Y }"
  ))
  total <- plot_set_nodes(plot_sets(dag))
  direct <- plot_set_nodes(plot_sets(dag, effect = "direct"))
  expect_identical(names(total)[total == "adjustment set"], "U")
  expect_setequal(names(direct)[direct == "adjustment set"], c("M", "U"))
  expect_identical(unname(total["M"]), "observed")
  expect_error(plot_sets(dag, effect = "natural"), "effect must be")
})


test_that("explicit endpoints control both calculation and displayed groups", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure]; Y [outcome]; T -> Y;",
    "X -> A; X -> B; A -> B }"
  ))
  incoming <- serialize(dag, NULL)
  nodes <- plot_set_nodes(plot_sets(dag, treatment = "A", outcome = "B"))
  expect_identical(unname(nodes[c("A", "B", "X", "T", "Y")]),
                    c("treatment", "outcome", "adjustment set", "observed", "observed"))
  expect_identical(serialize(dag, NULL), incoming)
  expect_error(plot_sets(dag, treatment = "missing"), "treatment must name")
  expect_error(plot_sets(dag, outcome = character()), "outcome must name")
  unmarked <- dagitty::dagitty("dag { A -> B }")
  expect_error(plot_sets(unmarked), "treatment must name")
})


test_that("joint queries preserve exact temporal names of adjustment nodes", {
  dag <- dagitty::dagitty(paste(
    "dag { T_t0 [exposure]; Y_t1 [outcome]; Y_t2 [outcome];",
    "C_t0 -> T_t0; C_t0 -> Y_t1; C_t1 -> T_t0; C_t1 -> Y_t2;",
    "T_t0 -> Y_t1; T_t0 -> Y_t2 }"
  ))
  nodes <- plot_set_nodes(plot_sets(dag))
  expect_setequal(names(nodes)[nodes == "adjustment set"], c("C_t0", "C_t1"))
  expect_setequal(names(nodes)[nodes == "outcome"], c("Y_t1", "Y_t2"))
  expect_setequal(names(nodes), names(dag))
  expect_true(dagitty::isAdjustmentSet(dag,
    names(nodes)[nodes == "adjustment set"],
    exposure = "T_t0", outcome = c("Y_t1", "Y_t2")))
})


test_that("set plotting bypasses causal role and node classification", {
  forbidden <- function(...){ # Fail if set plotting requests causal role classification
    stop("causal role classification was called")
  }
  local_mocked_bindings(classify_node_roles = forbidden,
    get_roles = forbidden, get_nodes = forbidden, .package = "causaliflower")
  dag <- dagitty::dagitty(
    "dag { T [exposure]; Y [outcome]; X -> T; X -> Y; T -> Y }"
  )
  expect_no_error(plot <- plot_sets(dag))
  expect_identical(unname(plot_set_nodes(plot)["X"]), "adjustment set")
  expect_error(plot_dagitty(dag), "causal role classification was called")
})


test_that("set plots retain labels, coordinates and legend controls", {
  dag <- dagitty::dagitty(
    "dag { T [exposure]; Y [outcome]; X -> T; X -> Y; T -> Y }"
  )
  dagitty::coordinates(dag) <- list(x = c(T = 1, X = 0, Y = 2),
                                    y = c(T = 0, X = 1, Y = 0))
  attr(dag, "source_note") <- "Keep supplied metadata"
  incoming <- serialize(dag, NULL)
  labels <- c(T = "Treatment", X = "Measured context", Y = "Outcome")
  plotted <- plot_sets(dag, labels = labels, curve = FALSE,
                       label_placement = "node", legend_position = "bottom-left")
  text <- Filter(function(layer){ return(inherits(layer$geom, "GeomText")) }, plotted$layers)[[1L]]
  expect_setequal(text$data$label, unname(labels))
  expect_identical(plotted$theme$legend.position, "bottom")
  expect_identical(plotted$theme$legend.justification, "left")
  hidden <- plot_sets(dag, include_legend = FALSE)
  expect_identical(hidden$theme$legend.position, "none")
  built_hidden <- ggplot2::ggplot_build(hidden)$plot
  expect_identical(built_hidden$scales$get_scales("colour")$map("adjustment set"),
                    "#111111")
  expect_identical(serialize(dag, NULL), incoming)
})


test_that("set plots share SVG and inferred PNG export with invisible returns", {
  dag <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y }")
  svg_file <- tempfile(fileext = ".svg")
  png_file <- tempfile(fileext = ".png")
  on.exit(unlink(c(svg_file, png_file)), add = TRUE)
  visible <- withVisible(plot_sets(dag))
  svg_result <- withVisible(plot_sets(dag, file = svg_file, width = 5, height = 3))
  png_result <- withVisible(plot_sets(dag, file = png_file, width = 5,
                                      height = 3, dpi = 72))
  expect_true(visible$visible)
  expect_false(svg_result$visible)
  expect_false(png_result$visible)
  expect_s3_class(svg_result$value, "ggplot")
  expect_s3_class(png_result$value, "ggplot")
  expect_match(paste(readLines(svg_file, warn = FALSE), collapse = ""), "<svg")
  expect_identical(readBin(png_file, "raw", n = 8L),
                    as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
  expect_error(plot_sets(dag, file = png_file, output_format = "svg"),
                 "does not match")
})
