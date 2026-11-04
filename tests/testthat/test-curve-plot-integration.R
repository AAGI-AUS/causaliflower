curve_test_paths <- function(plot){
  layers <- Filter(function(layer) inherits(layer$geom, "GeomEdgePath"), plot$layers)
  paths <- as.data.frame(data.table::rbindlist(lapply(layers, `[[`, "data"), fill = TRUE))
  paths <- paths[order(paths$edge_id, paths$index),
                 c("edge_id", "from", "to", "kind", "index", "x", "y"), drop = FALSE]
  rownames(paths) <- NULL
  paths
}

test_that("plot_dagitty() honours stored control points without changing the graph", {
  dag <- dagitty::dagitty('dag { T [exposure,pos="-1,0"] M [pos="0,0"] Y [outcome,pos="1,0"] T -> M -> Y; T -> Y }')
  manual <- causaliflower:::set_edge_controls(dag,
    data.frame(v = "T", e = "->", w = "Y", x = 0, y = 2.123456789012345))
  original <- manual
  saved_plot <- plot_dagitty(manual, curve = FALSE)
  saved_paths <- curve_test_paths(saved_plot)
  xy <- saved_paths[saved_paths$from == "T" & saved_paths$to == "Y", ]
  middle <- xy[which.min(abs(xy$index - 0.5)), ]
  expect_equal(middle$x, 0, tolerance = 1e-14)
  expect_equal(middle$y, -1.0617283945061725, tolerance = 1e-14)
  expect_identical(manual, original)

  # One row with unit column gaps: automatic curves use height 0.5.
  rerouted <- add_curves(manual, height = 0.5, mode = "absolute")
  expect_equal(curve_test_paths(plot_dagitty(manual)),
               curve_test_paths(plot_dagitty(rerouted, curve = FALSE)),
               tolerance = 1e-12)
  expect_identical(manual, original)
  expect_false(identical(causaliflower:::edge_controls(manual),
                         causaliflower:::edge_controls(rerouted)))
})

test_that("plotting honours stored non-directed control points", {
  for( type in c("<->", "--") ){
    graph_type <- if( type == "--" ) "pdag" else "dag"
    dag <- dagitty::dagitty(paste0(graph_type,
      ' { A [pos="-1,0"] B [pos="1,0"] A ', type, ' B }'))
    dag <- causaliflower:::set_edge_controls(dag,
      data.frame(v = "A", e = type, w = "B", x = 0, y = 2))
    original <- dag
    manual <- curve_test_paths(plot_dagitty(dag, curve = FALSE))
    automatic <- curve_test_paths(plot_dagitty(dag, curve = TRUE))
    expect_equal(manual, automatic, tolerance = 1e-12)
    expect_equal(manual$y[which.min(abs(manual$index - 0.5))], -1)
    expect_identical(dag, original)
  }
})

test_that("plot coordinates regenerate as a complete layout and remain local", {
  dag <- dagitty::dagitty('dag { A [pos="99,99"] A -> B -> C; A -> C }')
  original <- dag
  expected <- dagitty::coordinates(add_coords(dag, x_step = 4, y_step = 3))
  plot <- plot_dagitty(dag, x_step = 4, y_step = 3)
  nodes <- plot$data[!duplicated(plot$data$name), c("name", "x", "y")]
  expect_equal(nodes$x, unname(expected$x[as.character(nodes$name)]))
  expect_equal(nodes$y, -unname(expected$y[as.character(nodes$name)]))
  expect_identical(dag, original)
})

test_that("labels cannot influence graph routing", {
  dag <- dagitty::dagitty('dag { A [pos="0,0"] B [pos="1,0"] C [pos="2,0"] A -> B -> C; A -> C }')
  plain <- curve_test_paths(plot_dagitty(dag))
  labelled <- curve_test_paths(plot_dagitty(dag,
    labels = c(A = "A much longer treatment label", B = "Middle", C = "Result")))
  inside <- curve_test_paths(plot_dagitty(dag, label_placement = "node"))
  expect_equal(plain, labelled)
  expect_equal(plain, inside)
})
