plot_edge_paths <- function(plot){
  built <- ggplot2::ggplot_build(plot)
  edge_layers <- which(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomEdgePath"),
    logical(1L)
  ))

  paths <- list()
  for( layer_index in edge_layers ){
    layer_data <- built$data[[layer_index]]
    if( nrow(layer_data) == 0L ) next
    for( group in unique(layer_data$group) ){
      path <- layer_data[layer_data$group == group, , drop = FALSE]
      path <- path[order(path$index), , drop = FALSE]
      start <- c(path$x[1L], path$y[1L])
      end <- c(path$x[nrow(path)], path$y[nrow(path)])
      vector <- end - start
      length <- sqrt(sum(vector^2))
      offsets <- cbind(path$x - start[1L], path$y - start[2L])
      signed <- if( length > 0 ){
        (vector[1L] * offsets[, 2L] - vector[2L] * offsets[, 1L]) / length
      }else{
        rep(0, nrow(path))
      }
      paths[[length(paths) + 1L]] <- data.frame(
        x = start[1L],
        y = start[2L],
        xend = end[1L],
        yend = end[2L],
        deviation = max(abs(signed)),
        signed_deviation = signed[which.max(abs(signed))]
      )
    }
  }

  out <- do.call(rbind, paths)
  out[order(out$x, out$y, out$xend, out$yend), , drop = FALSE]
}

plot_layers_of <- function(plot, geom_class){
  Filter(function(layer) inherits(layer$geom, geom_class), plot$layers)
}


test_that("plot_dagitty() curves obstructed arrows and preserves supplied nodes", {
  dag <- dagitty::dagitty("dag { A -> B -> C; A -> C }")
  dagitty::coordinates(dag) <- list(
    x = c(A = 0, B = 1, C = 2),
    y = c(A = 0, B = 0, C = 0)
  )
  original_coordinates <- dagitty::coordinates(dag)

  for( placement in c("above", "node") ){
    default_plot <- plot_dagitty(dag, label_placement = placement)
    auto_plot <- plot_dagitty(dag, label_placement = placement, curve = TRUE)
    repeated_plot <- plot_dagitty(dag, label_placement = placement, curve = TRUE)
    none_plot <- plot_dagitty(dag, label_placement = placement, curve = FALSE)

    default_paths <- plot_edge_paths(default_plot)
    auto_paths <- plot_edge_paths(auto_plot)
    none_paths <- plot_edge_paths(none_plot)

    expect_equal(default_paths, auto_paths)
    expect_equal(auto_paths, plot_edge_paths(repeated_plot))
    expect_equal(auto_paths$deviation, c(0, auto_paths$deviation[2L], 0),
                 tolerance = 1e-10)
    expect_gt(auto_paths$deviation[2L], 0.075)
    expect_true(all(none_paths$deviation < 1e-10))

    built <- ggplot2::ggplot_build(auto_plot)
    point_layer <- which(vapply(
      auto_plot$layers,
      function(layer) inherits(layer$geom, "GeomPoint"),
      logical(1L)
    ))
    expect_length(point_layer, 1L)
    expect_setequal(built$data[[point_layer]]$x, original_coordinates$x)
    expect_setequal(built$data[[point_layer]]$y, -original_coordinates$y)
  }

  expect_equal(dagitty::coordinates(dag), original_coordinates)
  y_range <- ggplot2::ggplot_build(plot_dagitty(dag))$
    layout$panel_params[[1L]]$y.range
  expect_gt(diff(y_range), 0)
  expect_lte(y_range[1L], 0)
  expect_gte(y_range[2L], 0)
})


test_that("automatic curves match add_curves() with the documented height", {
  dag <- dagitty::dagitty("dag { A -> B -> C; A -> C; D -> C; A -> D }")
  dagitty::coordinates(dag) <- list(
    x = c(A = 0, B = 1, C = 2, D = 1),
    y = c(A = 0, B = 0, C = 0, D = 2)
  )
  # One column holds two nodes two units apart, so the height is 0.75 * 2.
  curved <- add_curves(dag, height = 1.5, mode = "absolute")
  expect_equal(plot_edge_paths(plot_dagitty(dag)),
               plot_edge_paths(plot_dagitty(curved, curve = FALSE)),
               tolerance = 1e-12)
})


test_that("labels sit a fixed physical distance above or below their nodes", {
  dag <- dagitty::dagitty("dag { A -> B -> C }")
  dagitty::coordinates(dag) <- list(
    x = c(A = 0, B = 1, C = 2),
    y = c(A = 0, B = 0, C = 0)
  )
  labels <- c(A = "alpha node", B = "beta node", C = "gamma node")

  above <- plot_layers_of(plot_dagitty(dag, labels = labels), "GeomDagNodeLabel")
  expect_length(above, 1L)
  expect_gt(above[[1L]]$geom_params$offset, 0)
  expect_identical(above[[1L]]$aes_params$vjust, 0)
  expect_setequal(above[[1L]]$data$display_label, unname(labels))

  below <- plot_layers_of(
    plot_dagitty(dag, labels = labels, label_placement = "below"),
    "GeomDagNodeLabel"
  )
  expect_length(below, 1L)
  expect_lt(below[[1L]]$geom_params$offset, 0)
  expect_identical(below[[1L]]$aes_params$vjust, 1)

  inside <- plot_dagitty(dag, labels = labels, label_placement = "node")
  expect_length(plot_layers_of(inside, "GeomDagNodeLabel"), 0L)
  expect_length(plot_layers_of(inside, "GeomText"), 1L)

  # The offset is in millimetres, so it does not depend on the data scale.
  wide <- dag
  dagitty::coordinates(wide) <- list(
    x = c(A = 0, B = 1000, C = 2000),
    y = c(A = 0, B = 0, C = 0)
  )
  wide_above <- plot_layers_of(plot_dagitty(wide, labels = labels), "GeomDagNodeLabel")
  expect_identical(wide_above[[1L]]$geom_params$offset,
                   above[[1L]]$geom_params$offset)
})


test_that("automatic label placement uses the side with more room", {
  nodes <- data.frame(
    x = c(0, 0, 0, 1),
    y = c(0, 1, 3, 0)
  )
  # The bottom node has unlimited room below it; the middle node has more
  # room above (2) than below (1).
  expect_identical(plot_label_sides(nodes, "auto"),
                   c("below", "above", "above", "above"))
  expect_identical(plot_label_sides(nodes, "below"), rep("below", 4L))
})


test_that("labels wrap at spaces only and margins make room for them", {
  dag <- dagitty::dagitty("dag { A -> B }")
  dagitty::coordinates(dag) <- list(x = c(A = 0, B = 1), y = c(A = 0, B = 0))
  long <- c(A = "first node with a long descriptive name",
            B = "Second.node_without.spaces_but.very.long")
  plot <- plot_dagitty(dag, labels = long)
  layer <- plot_layers_of(plot, "GeomDagNodeLabel")[[1L]]
  wrapped <- layer$data$display_label[layer$data$name == "A"]
  expect_match(wrapped, "\n", fixed = TRUE)
  expect_identical(gsub("\n", " ", wrapped, fixed = TRUE), long[["A"]])
  expect_identical(layer$data$display_label[layer$data$name == "B"], long[["B"]])

  margin <- as.numeric(plot$theme$plot.margin)
  short <- as.numeric(plot_dagitty(dag, labels = c(A = "A", B = "B"))$theme$plot.margin)
  expect_gt(margin[4L], short[4L])
  expect_gt(margin[2L], short[2L])
})


test_that("legend position is one lever with a centred fallback", {
  cases <- list(
    list("top", "top", "center"),
    list("top-left", "top", "left"),
    list("Top-Right", "top", "right"),
    list("bottom-centre", "bottom", "center"),
    list("bottom center", "bottom", "center"),
    list("left-top", "left", "top"),
    list("right-bottom", "right", "bottom"),
    list("left-left", "left", "center"),
    list("top-top", "top", "center"),
    list("right-middle", "right", "center"),
    list("left", "left", "center"),
    list("none", "none", "center")
  )
  for( case in cases ){
    expect_identical(plot_legend_position(case[[1L]]),
                     list(side = case[[2L]], justification = case[[3L]]),
                     info = case[[1L]])
  }
  expect_identical(plot_legend_position(NULL),
                   list(side = "top", justification = "center"))
  for( value in list("middle", 1, c("top", "left"), NA_character_) ){
    expect_error(plot_dagitty(string_to_dag("A -> B"), legend_position = value),
                 "legend_position")
  }

  dag <- string_to_dag("A -> B")
  plot <- plot_dagitty(dag, legend_position = "right-top")
  expect_identical(plot$theme$legend.position, "right")
  expect_identical(plot$theme$legend.justification, "top")
  expect_identical(plot_dagitty(dag)$theme$legend.position, "top")
  expect_identical(plot_dagitty(dag, legend_position = "none")$theme$legend.position,
                   "none")
})


test_that("every role present is coloured and latent nodes are hollow", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent] C -> T C -> Y T -> M M -> Y",
    "L -> M L -> Y U -> T U -> Y }"
  ))
  plot <- plot_dagitty(dag)
  point_data <- plot_layers_of(plot, "GeomPoint")[[1L]]$data
  present <- unique(as.character(point_data$role))
  breaks <- plot$scales$get_scales("colour")$breaks
  expect_setequal(breaks, present)
  expect_true("latent" %in% breaks)

  built <- ggplot2::ggplot_build(plot)
  point_index <- which(vapply(plot$layers, function(layer){
    inherits(layer$geom, "GeomPoint")
  }, logical(1L)))
  points <- built$data[[point_index]]
  roles <- as.character(point_data$role)
  expect_identical(unique(points$fill[roles == "latent"]), "white")
  expect_identical(unique(points$colour[roles == "latent"]), "black")
  expect_true(all(points$fill[roles != "latent"] == points$colour[roles != "latent"]))

  hidden <- plot_dagitty(dag, include_legend = FALSE)
  expect_length(hidden$scales$get_scales("colour")$breaks, 0L)
})


test_that("automatic routing leaves unpositioned bidirected arcs unchanged", {
  dag <- dagitty::dagitty(
    "dag { A [exposure] C [outcome] A -> B -> C; A -> C; D <-> E }"
  )
  dagitty::coordinates(dag) <- list(
    x = c(A = 0, B = 1, C = 2, D = 0, E = 2),
    y = c(A = 0, B = 0, C = 0, D = 1, E = 1)
  )
  original_coordinates <- dagitty::coordinates(dag)
  auto <- plot_dagitty(dag, curve = TRUE)
  none <- plot_dagitty(dag, curve = FALSE)

  bidirected_layer <- function(plot){
    layer <- which(vapply(plot$layers, function(candidate){
      arrow <- candidate$geom_params$arrow
      !is.null(arrow) && identical(arrow$ends, 3L)
    }, logical(1L)))
    expect_length(layer, 1L)
    layer
  }
  bidirected_path <- function(plot){
    built <- ggplot2::ggplot_build(plot)$data[[bidirected_layer(plot)]]
    built[order(built$group, built$index), c("x", "y"), drop = FALSE]
  }

  expect_equal(bidirected_path(auto), bidirected_path(none))
  expect_identical(
    auto$layers[[bidirected_layer(auto)]]$aes_params$edge_linetype,
    "dashed"
  )
  expect_equal(dagitty::coordinates(dag), original_coordinates)
})


test_that("automatic routing retains undirected edges", {
  dag <- dagitty::dagitty(
    "pdag { A -> B -> C; A -> C; D -- E }"
  )
  dagitty::coordinates(dag) <- list(
    x = c(A = 0, B = 1, C = 2, D = 0, E = 2),
    y = c(A = 0, B = 0, C = 0, D = 1, E = 1)
  )

  auto_paths <- plot_edge_paths(plot_dagitty(dag, curve = TRUE))
  none_paths <- plot_edge_paths(plot_dagitty(dag, curve = FALSE))

  undirected_auto <- auto_paths[
    abs(auto_paths$y) == 1 & abs(auto_paths$yend) == 1,
    , drop = FALSE
  ]
  undirected_none <- none_paths[
    abs(none_paths$y) == 1 & abs(none_paths$yend) == 1,
    , drop = FALSE
  ]
  expect_equal(nrow(undirected_auto), 1L)
  expect_equal(nrow(undirected_none), 1L)
  expect_equal(undirected_auto$deviation, 0, tolerance = 1e-10)
})


test_that("edges fade only on dense graphs", {
  sparse <- plot_layers_of(plot_dagitty(string_to_dag("A -> B")), "GeomEdgePath")
  expect_equal(sparse[[1L]]$aes_params$edge_alpha, 0.85)
  # 105 edges: opacity 3 / sqrt(edges) reaches its 0.3 floor.
  nodes <- paste0("N", 1:15)
  edges <- unlist(lapply(seq_len(14L), function(i){
    paste(nodes[i], "->", nodes[(i + 1L):15L])
  }))
  dense <- plot_layers_of(plot_dagitty(string_to_dag(edges)), "GeomEdgePath")
  expect_equal(dense[[1L]]$aes_params$edge_alpha, 0.3)
})


test_that("plot_dagitty() validates arguments without shifting old ones", {
  dag <- string_to_dag("A -> B")
  invalid <- list("auto", 1, c(TRUE, FALSE), NA, NULL)
  for( value in invalid ){
    expect_error(
      plot_dagitty(dag, curve = value),
      "plot_dagitty(): curve must be TRUE or FALSE.",
      fixed = TRUE
    )
  }

  # Arguments up to curve keep their 0.0.1 positions; new ones follow.
  expect_identical(
    names(formals(plot_dagitty))[1:10],
    c("dag", "seed", "x_step", "y_step", "layout_shape", "labels",
      "label_type", "label_placement", "include_legend", "curve")
  )
  expect_identical(formals(plot_dagitty)$curve, TRUE)

  expect_error(plot_dagitty(dag, label_placement = "side"), "label_placement")
  for( argument in c("node_size", "label_size", "arrow_size", "dpi") ){
    for( value in list(0, -1, Inf, "1", c(1, 2)) ){
      arguments <- stats::setNames(list(dag, value), c("dag", argument))
      expect_error(do.call(plot_dagitty, arguments), argument, fixed = TRUE)
    }
  }
  expect_error(plot_dagitty(dag, edge_gap = -1), "edge_gap")
  expect_error(plot_dagitty(dag, file = tempfile(), width = 0), "width")
  expect_error(plot_dagitty(dag, file = tempfile(), height = "4"), "height")
  expect_error(plot_dagitty(dag, file = c("a.svg", "b.svg")), "file")
})


test_that("seed and text_box are deprecated aliases", {
  dag <- string_to_dag("A -> B")
  expect_warning(seeded <- plot_dagitty(dag, seed = 1),
                 "'seed' is deprecated as of 0.0.2", fixed = TRUE)
  expect_warning(boxed <- plot_dagitty(dag, label_placement = "text_box"),
                 "deprecated as of 0.0.2. Use \"above\"", fixed = TRUE)
  plain <- plot_dagitty(dag)
  expect_equal(plot_edge_paths(seeded), plot_edge_paths(plain))
  expect_equal(plot_layers_of(boxed, "GeomDagNodeLabel")[[1L]]$geom_params,
               plot_layers_of(plain, "GeomDagNodeLabel")[[1L]]$geom_params)
})


test_that("plots are saved as SVG or PNG with a fitted canvas", {
  dag <- dagitty::dagitty("dag { X -> T X -> Y T -> M M -> Y T -> Y }")
  svg_file <- tempfile(fileext = ".svg")
  returned <- withVisible(plot_dagitty(dag, file = svg_file))
  expect_false(returned$visible)
  expect_s3_class(returned$value, "ggplot")
  expect_true(file.exists(svg_file))
  expect_match(paste(readLines(svg_file, n = 5L, warn = FALSE), collapse = ""),
               "<svg", fixed = TRUE)

  png_file <- tempfile(fileext = ".png")
  plot_dagitty(dag, file = png_file)
  expect_identical(readBin(png_file, "raw", 8L)[2:4], charToRaw("PNG"))

  untyped <- tempfile()
  plot_dagitty(dag, file = untyped, output_format = "PNG")
  expect_identical(readBin(untyped, "raw", 8L)[2:4], charToRaw("PNG"))
  expect_error(plot_dagitty(dag, file = tempfile(fileext = ".svg"),
                            output_format = "png"),
               "does not match")
  expect_error(plot_dagitty(dag, file = tempfile(), output_format = "pdf"),
               "output_format")

  # The fitted canvas grows with the graph and stays within bounds.
  small <- plot_fit_size(
    data.frame(x = c(0, 1), y = c(0, 0), display_label = c("A", "B")),
    c("above", "above"), 5, 7, 4.5,
    c(top = 15, right = 8, bottom = 8, left = 8), 2, "top", "outcome"
  )
  large <- plot_fit_size(
    data.frame(x = rep(0:6, each = 6), y = rep(seq(0, 7.5, by = 1.5), 7),
               display_label = "a label of moderate length"),
    rep("above", 42L), 5, 7, 4.5,
    c(top = 15, right = 40, bottom = 8, left = 40), 2, "top",
    c("outcome", "treatment", "confounder", "mediator", "instrument", "latent")
  )
  expect_true(all(small >= c(4, 3)))
  expect_true(all(large > small))
  expect_true(all(is.finite(large)))
})


test_that("saving falls back to base devices without svglite or ragg", {
  testthat::skip_if_not(isTRUE(capabilities("cairo")))
  dag <- string_to_dag("A -> B")
  testthat::local_mocked_bindings(plot_package_available = function(package) FALSE)
  svg_file <- tempfile(fileext = ".svg")
  plot_dagitty(dag, file = svg_file)
  expect_match(paste(readLines(svg_file, n = 5L, warn = FALSE), collapse = ""),
               "<svg", fixed = TRUE)
  png_file <- tempfile(fileext = ".png")
  plot_dagitty(dag, file = png_file)
  expect_identical(readBin(png_file, "raw", 8L)[2:4], charToRaw("PNG"))
})


test_that("SVG startup failures leave existing devices and files intact", {
  testthat::skip_if_not(isTRUE(capabilities("cairo")))
  testthat::local_mocked_bindings(plot_package_available = function(package) FALSE)
  testthat::local_mocked_bindings(svg = function(...) {
    warning("failed to load cairo DLL", call. = FALSE)
    invisible(NULL)
  }, .package = "grDevices")

  grDevices::pdf(NULL, width = 5, height = 3)
  device <- grDevices::dev.cur()
  devices <- grDevices::dev.list()
  size <- grDevices::dev.size("in")
  on.exit({
    if( device %in% grDevices::dev.list() ) grDevices::dev.off(device)
  }, add = TRUE)
  svg_file <- tempfile(fileext = ".svg")
  writeLines("existing SVG", svg_file)
  on.exit(unlink(svg_file), add = TRUE)

  expect_warning(
    expect_error(plot_dagitty(string_to_dag("A -> B"), file = svg_file),
                 "SVG device could not start.*svglite.*XQuartz"),
    "failed to load cairo DLL", fixed = TRUE
  )
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(grDevices::dev.list(), devices)
  expect_identical(grDevices::dev.size("in"), size)
  expect_identical(readLines(svg_file), "existing SVG")
})


test_that("every node keeps its label when labels repeat", {
  drawn_labels <- function(plot){
    built <- ggplot2::ggplot_build(plot)
    drawn <- NULL
    for( index in seq_along(plot$layers) ){
      if( inherits(plot$layers[[index]]$geom, "GeomEdgePath") ) next
      layer_data <- built$data[[index]]
      if( !"label" %in% names(layer_data) ) next
      keep <- !is.na(layer_data$label) & nzchar(as.character(layer_data$label))
      drawn <- rbind(drawn, data.frame(
        x = layer_data$x[keep],
        label = as.character(layer_data$label)[keep]
      ))
    }
    drawn <- drawn[order(drawn$x), , drop = FALSE]
    drawn$label
  }

  dag <- dagitty::dagitty("dag { A -> B -> C }")
  dagitty::coordinates(dag) <- list(
    x = c(A = 0, B = 1, C = 2),
    y = c(A = 0, B = 0, C = 0)
  )

  for( placement in c("node", "above", "below", "auto") ){
    repeated <- plot_dagitty(dag, labels = c(A = "Q", B = "Q", C = "R"),
                             label_placement = placement)
    expect_identical(drawn_labels(repeated), c("Q", "Q", "R"))
  }
})


test_that("unnamed labels are refused", {
  dag <- dagitty::dagitty("dag { A -> B -> C }")
  dagitty::coordinates(dag) <- list(
    x = c(A = 0, B = 1, C = 2),
    y = c(A = 0, B = 0, C = 0)
  )

  for( placement in c("node", "above") ){
    expect_error(
      plot_dagitty(dag, labels = c("Q", "R", "S"), label_placement = placement),
      "The supplied labels are not named after the nodes in the graph. Please name each label with its node and try again.",
      fixed = TRUE
    )
    expect_error(
      plot_dagitty(dag, labels = c(A = "Q", "R", C = "S"),
                   label_placement = placement),
      "The supplied labels are not named after the nodes in the graph. Please name each label with its node and try again.",
      fixed = TRUE
    )
    named <- plot_dagitty(dag, labels = c(A = "Q", B = "R", C = "S"),
                          label_placement = placement)
    expect_s3_class(named, "ggplot")
  }
})


test_that("free-form coordinates save on a bounded canvas", {
  expect_identical(plot_columns(c(0, 0.01, 1, 1.02, 2)), c(1L, 1L, 2L, 2L, 3L))
  # Many evenly spaced columns stay separate rather than chaining together.
  expect_identical(plot_columns(seq(0, 87, by = 3)), 1:30)
  dag <- dagitty::getExample("Shrier")
  svg_file <- tempfile(fileext = ".svg")
  expect_no_error(plot_dagitty(dag, file = svg_file))
  expect_true(file.exists(svg_file))
  header <- paste(readLines(svg_file, n = 5L, warn = FALSE), collapse = " ")
  dimension <- function(name){
    pattern <- paste0("(?<= ", name, "=')[0-9.]+")
    as.numeric(regmatches(header, regexpr(pattern, header, perl = TRUE)))
  }
  expect_lte(dimension("width"), 40 * 72)
  expect_lte(dimension("height"), 40 * 72)

  # User-supplied sizes beyond ggsave's 50-inch limit are allowed.
  wide <- tempfile(fileext = ".svg")
  expect_no_error(plot_dagitty(string_to_dag("A -> B"), file = wide,
                               width = 60, height = 3))
})


test_that("a graph without a legend keeps room for its labels", {
  dag <- dagitty::dagitty("dag { A [exposure] C [outcome] A -> B -> C }")
  dagitty::coordinates(dag) <- list(x = c(A = 0, B = 1, C = 2),
                                    y = c(A = 0, B = 1, C = 0))
  hidden <- plot_dagitty(dag, include_legend = FALSE)
  shown <- plot_dagitty(dag)
  expect_identical(hidden$theme$legend.position, "none")
  hidden_top <- as.numeric(hidden$theme$plot.margin)[1L]
  shown_room <- as.numeric(shown$theme$legend.box.spacing)
  expect_gt(hidden_top, 4)
  expect_equal(hidden_top, shown_room + 2)
  unmatched <- plot_dagitty(dag, include_legend = "collider")
  expect_identical(unmatched$theme$legend.position, "none")
})


test_that("supplied line breaks are kept and missing labels are not drawn", {
  dag <- dagitty::dagitty("dag { A -> B -> C }")
  dagitty::coordinates(dag) <- list(x = c(A = 0, B = 1, C = 2),
                                    y = c(A = 0, B = 0, C = 0))
  drawn <- function(plot){
    layers <- plot_layers_of(plot, "GeomDagNodeLabel")
    unlist(lapply(layers, function(layer) layer$data$display_label), use.names = FALSE)
  }
  broken <- plot_dagitty(dag, labels = c(A = "Line one\nLine two", B = "B", C = "C"))
  expect_true("Line one\nLine two" %in% drawn(broken))

  missing <- plot_dagitty(dag, labels = c(A = NA, B = "B", C = ""))
  expect_identical(drawn(missing), "B")
  misnamed <- plot_dagitty(dag, labels = c(A = "A", B = "B", Z = "Z"))
  expect_setequal(drawn(misnamed), c("A", "B"))
  inside <- plot_dagitty(dag, labels = c(A = NA, B = "B", C = "C"),
                         label_placement = "node")
  expect_setequal(plot_layers_of(inside, "GeomText")[[1L]]$data$label, c("B", "C"))
})


test_that("unsupported file extensions are refused", {
  dag <- string_to_dag("A -> B")
  for( file in c("figure.pdf", "figure.jpeg") ){
    expect_error(plot_dagitty(dag, file = file.path(tempdir(), file)),
                 "must end in .svg or .png", fixed = TRUE)
  }
})


test_that("the node label geom follows the installed GeomLabel parameters", {
  parameters <- GeomDagNodeLabel$parameters()
  expect_true(all(ggplot2::GeomLabel$parameters() %in% parameters))
  expect_true("offset" %in% parameters)
})
