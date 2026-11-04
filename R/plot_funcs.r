# ==============================================================================
# plot_funcs.r


# plot_dagitty() and its drawing helpers
# ==============================================================================


#' Plot a dagitty object
#'
#' Draws a dagitty object as a causal diagram: nodes coloured by causal role,
#' labels attached to their nodes, curved edges and a role legend. Node
#' coordinates and edge paths are added when needed, for plotting only. The
#' supplied dagitty object is unchanged. The plot can also be saved as an SVG
#' or PNG file sized to fit the graph.
#'
#' @importFrom dagitty coordinates
#' @importFrom ggplot2 ggplot scale_colour_manual labs guides guide_legend theme
#' @importFrom ggraph circle
#' @importFrom grid unit arrow
#' @param dag A dagitty object.
#' @param seed Deprecated as of 0.0.2 and ignored: label placement no longer
#'   breaks ties at random.
#' @param x_step Horizontal spacing between node columns when coordinates are
#'   generated. Single positive number.
#' @param y_step Minimum vertical gap between nodes in the same column when
#'   coordinates are generated. Single positive number.
#' @param layout_shape Vertical placement of nodes in the same column, between
#'   0 and 1: 0 places nodes nearest their parents' paths, 1 spreads the column
#'   evenly.
#' @param labels A named vector of labels, one per node, named after the
#'   nodes. Labels wrap at spaces only.
#' @param label_type Label style passed to [get_labels()] when `labels` is not
#'   supplied. Defaults to `"name"`. `"initials"` and `"short"` are also
#'   available.
#' @param label_placement Where node labels are drawn: `"above"` (the default),
#'   `"below"`, `"auto"` or `"node"`. `"auto"` chooses the side with more space
#'   to the nearest node in the same column, using above when tied. `"node"`
#'   places labels inside nodes and suits `label_type = "initials"`.
#'   `"text_box"` is a deprecated alias of `"above"` as of 0.0.2.
#' @param include_legend Node roles drawn in colour and listed in the legend.
#'   Defaults to all roles. Only roles present in the graph are listed. Roles
#'   left out are drawn in almost black. Use `FALSE` to draw every node in
#'   almost black without a legend.
#' @param curve By default, directed arrows automatically curve to avoid
#'   passing through the centre of another node and reduce edge crossings
#'   (`curve = TRUE`). Use `FALSE` to skip automatic curving. Any curves already
#'   stored in the graph are retained.
#' @param legend_position Legend side and alignment in one string. The side is
#'   `"top"` (the default), `"bottom"`, `"left"` or `"right"`, optionally
#'   followed by an alignment along that side: `"left"`, `"right"` or
#'   `"centre"` for the top and bottom, and `"top"`, `"bottom"` or `"centre"`
#'   for the left and right, as in `"top-left"` or `"right-bottom"`. A missing
#'   or incompatible alignment is centred. Use `"none"` to omit the legend.
#' @param node_size Node size, as a ggplot2 point size. Defaults to `13`.
#' @param label_size Label text size in millimetres. Defaults to `4.5`.
#' @param arrow_size Arrowhead length in points. Defaults to `4.5`.
#' @param edge_gap Extra space in millimetres between each edge end and its
#'   node. Defaults to `1`.
#' @param file Optional file path. When supplied, the plot is also saved to
#'   this file and returned invisibly.
#' @param output_format File format used with `file`: `"svg"` (the default) or
#'   `"png"`. A `file` ending in `.png` is saved as PNG unless `output_format`
#'   is supplied. Other file extensions are refused.
#' @param width,height Size of the saved file in inches. By default, the size
#'   is fitted to the graph's columns, rows, labels and legend.
#' @param dpi Resolution of a saved PNG file. Defaults to `300`, reduced if
#'   needed to keep the longer side within 10,000 pixels.
#'
#' @details
#' Node, label, arrow and legend sizes use physical units in both SVG and PNG
#' output. The plot margins are calculated to leave room for labels and the
#' legend. A small plotting window can still be crowded. By default, saving to
#' a file chooses dimensions to fit the graph.
#'
#' Automatic curves (`curve = TRUE`) match
#' `add_curves(dag, height = h, mode = "absolute")`, where `h` is three
#' quarters of the median vertical gap between nodes in the same column, or
#' half the median horizontal gap between columns when no column has two
#' nodes. Neighbouring nodes whose horizontal positions differ by less than 4%
#' of the graph's width, and by less than half the median gap between
#' positions, count as one column. To choose the curve height or position, call [add_curves()] first,
#' then plot with `curve = FALSE`. Bidirected edges are dashed with arrowheads
#' at both ends. Undirected edges have no arrowheads. Latent nodes are drawn
#' as hollow circles.
#'
#' SVG files are written with the svglite package and PNG files with the ragg
#' package when they are installed. Otherwise [grDevices::svg()] or
#' [grDevices::png()] is used.
#'
#' @return A ggplot object, returned invisibly when `file` is supplied.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' plot_dagitty(dag)
#'
#' plot_dagitty(dag, legend_position = "right-top", label_placement = "below")
#'
#' \donttest{
#' svg_file <- tempfile(fileext = ".svg")
#' plot_dagitty(dag, file = svg_file)
#' }
#'
#' @export
plot_dagitty <- function(dag,
                         seed = NULL,
                         x_step = 3,
                         y_step = 1.5,
                         layout_shape = 1,
                         labels = NULL,
                         label_type = "name",
                         label_placement = "above",
                         include_legend = c("outcome",
                                            "treatment",
                                            "confounder",
                                            "mediator",
                                            "mediator_outcome_confounder",
                                            "instrument",
                                            "competing_cause",
                                            "collider",
                                            "latent",
                                            "observed",
                                            "undetermined"),
                         curve = TRUE,
                         legend_position = "top",
                         node_size = 13,
                         label_size = 4.5,
                         arrow_size = 4.5,
                         edge_gap = 1,
                         file = NULL,
                         output_format = "svg",
                         width = NULL,
                         height = NULL,
                         dpi = 300
                         ){

  if( !is.logical(curve) || length(curve) != 1L ||
      is.na(curve) ){
    stop("plot_dagitty(): curve must be TRUE or FALSE.", call. = FALSE)
  }
  if( !is.null(seed) ){
    warning("plot_dagitty(): 'seed' is deprecated as of 0.0.2 and is ignored; ",
            "label placement no longer breaks ties at random.", call. = FALSE)
  }
  if( !is.character(label_placement) || length(label_placement) != 1L ||
      is.na(label_placement) ||
      !label_placement %in% c("above", "below", "auto", "node", "text_box") ){
    stop("plot_dagitty(): label_placement must be \"above\", \"below\", ",
         "\"auto\" or \"node\".", call. = FALSE)
  }
  if( label_placement == "text_box" ){
    warning("plot_dagitty(): label_placement = \"text_box\" is deprecated as of ",
            "0.0.2. Use \"above\" instead.", call. = FALSE)
    label_placement <- "above"
  }
  legend <- plot_legend_position(legend_position)
  for( argument in c("node_size", "label_size", "arrow_size", "dpi") ){
    value <- get(argument)
    if( !is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0 ){
      stop("plot_dagitty(): ", argument, " must be a single positive number.",
           call. = FALSE)
    }
  }
  if( !is.numeric(edge_gap) || length(edge_gap) != 1L ||
      !is.finite(edge_gap) || edge_gap < 0 ){
    stop("plot_dagitty(): edge_gap must be a single non-negative number.",
         call. = FALSE)
  }
  format <- plot_output_format(file, output_format,
                               missing(output_format))
  for( argument in c("width", "height") ){
    value <- get(argument)
    if( !is.null(value) && ( !is.numeric(value) || length(value) != 1L ||
                             !is.finite(value) || value <= 0 ) ){
      stop("plot_dagitty(): ", argument, " must be NULL or a single positive ",
           "number of inches.", call. = FALSE)
    }
  }

  # Share structural reads between the edge check and role extraction.
  .cache <- graph_cache()

  if( nrow( data.table::as.data.table( cached_edges(dag, .cache) ) ) == 0 ){
    stop("plot_dagitty() requires a graph with at least one edge. Roles and coordinates for graphs without edges remain available via get_roles() and add_coords().")
  }

  roles <- classify_node_roles(dag, .cache = .cache)

  controls <- edge_controls(dag, "plot_dagitty")
  dag <- dag_for_layout(dag, "plot_dagitty")
  if( nrow(controls) > 0L ){
    if( curve ) controls <- controls[controls$e != "->", , drop = FALSE]
    dag <- set_edge_controls(dag, controls, replace = dag_edges(dag, "plot_dagitty()"),
                              caller = "plot_dagitty")
  }
  existing_coordinates <- dagitty::coordinates(dag)

  if( anyNA(existing_coordinates$x) || anyNA(existing_coordinates$y) ){
    dag <- layout_coords(dag, x_step = x_step, y_step = y_step, layout_shape = layout_shape)
  }
  if( any(!is.finite(dagitty::coordinates(dag)$x)) ||
      any(!is.finite(dagitty::coordinates(dag)$y)) ){
    stop("plot_dagitty(): node coordinates must be finite.", call. = FALSE)
  }

  # aesthetics for the role legend
  col <- c("outcome"="deepskyblue",
           "treatment"="darkolivegreen2",
           "confounder"="coral2",
           "mediator"="darkorchid1",
           "mediator_outcome_confounder"="magenta4",
           "instrument"="deeppink1",
           "competing_cause"="darkseagreen4",
           "collider"="darkred",
           "latent"="black",
           "undetermined"="grey60",
           "observed"="#111111")
  col[!names(col) %in% include_legend] <- "#111111"
  # Latent nodes are hollow; other nodes are filled with their role colour.
  fill_col <- col
  fill_col[["latent"]] <- "white"
  # variable for legend order
  order_col <- c("outcome",
                 "treatment",
                 "confounder",
                 "mediator",
                 "mediator_outcome_confounder",
                 "instrument",
                 "competing_cause",
                 "collider",
                 "latent",
                 "undetermined",
                 "observed")
  order_col <- order_col[ order_col %in% include_legend ]

  if( is.null(labels) ){
    labels <- get_labels(dag, label_type)
  }

  if(!is.null(labels) & length(labels) != length(names(dag))){
    stop("The length of supplied labels does not equal the number of nodes in the graph. Please check labels input and try again.")
  }

  if(!is.null(labels) && (is.null(names(labels)) || any(!nzchar(names(labels))))){
    stop("The supplied labels are not named after the nodes in the graph. Please name each label with its node and try again.")
  }

  dag_df <- tidy_ggdagitty(dag, labels, roles = roles)
  dag_df_complete_cases <- dag_df[complete.cases(dag_df[, "role"]), ]

  if( nrow(dag_df) > nrow(dag_df_complete_cases) ){
    message("Nodes without arrows removed and not displayed in the plotted graph.")
    dag_df <- dag_df_complete_cases
  }

  dag_df[sapply(dag_df, is.character)] <- lapply( dag_df[sapply(dag_df, is.character)], as.factor )
  rendered_paths <- stored_edge_paths(dag_df, edge_controls(dag, "plot_dagitty"))
  if( curve ){
    other_paths <- rendered_paths[rendered_paths$kind != "directed", , drop = FALSE]
    edge_plan <- checked_curve_edge_plan(dag_df,
                                          height = plot_curve_height(dag_df),
                                          mode = "absolute",
                                          other_paths = other_paths,
                                          caller = "plot_dagitty")
    rendered_paths <- as.data.frame(data.table::rbindlist(
      list(edge_plan$paths, other_paths), use.names = TRUE, fill = TRUE
    ))
    rendered_paths$.edge_group <- match(rendered_paths$edge_id,
                                        unique(rendered_paths$edge_id))
  }else{
    edge_plan <- list(anchors = curve_plot_anchors(dag_df))
  }

  nodes <- dag_df[!duplicated(as.character(dag_df$name)), , drop = FALSE]
  # Missing or empty labels are not drawn. Lines wrap at spaces only, and line
  # breaks in supplied labels are kept.
  nodes$label <- vapply(seq_len(nrow(nodes)), function(i){
    value <- nodes$label[[i]]
    text <- if( is.null(value) || length(value) == 0L || is.na(value[1L]) ){
      NA_character_
    }else{
      as.character(value[1L])
    }
    if( !is.na(text) && !nzchar(trimws(text)) ) text <- NA_character_

    return(text)
  }, character(1L))
  nodes$display_label <- vapply(nodes$label, function(label){
    if( is.na(label) ) return(NA_character_)
    lines <- strsplit(label, "\n", fixed = TRUE)[[1L]]
    wrapped <- unlist(lapply(lines, function(line){
      parts <- strwrap(line, width = 26L)
      if( length(parts) == 0L ) parts <- ""

      return(parts)
    }), use.names = FALSE)
    text <- paste(wrapped, collapse = "\n")

    return(text)
  }, character(1L), USE.NAMES = FALSE)
  present <- order_col[order_col %in% as.character(nodes$role)]
  if( length(present) == 0L ){
    legend$side <- "none"
  }

  # Physical sizes in millimetres. A filled circle's radius is 0.375 of the
  # point's font size; edges stop just outside it.
  node_radius <- 0.375 * (node_size * ggplot2::.pt + ggplot2::.stroke / 2) /
    ggplot2::.pt
  label_offset <- node_radius + 2
  edge_cap <- ggraph::circle(node_radius + 0.8 + edge_gap, "mm")
  edge_count <- length(unique(rendered_paths$edge_id))
  edge_alpha <- min(0.85, max(0.3, 3 / sqrt(max(1L, edge_count))))

  edge_layers <- plot_dag_edge_layers(
    rendered_paths,
    edge_plan = edge_plan,
    arrow_directed = grid::arrow(angle = 25, length = grid::unit(arrow_size, "pt"),
                                 type = "closed"),
    arrow_bidirected = grid::arrow(angle = 25, length = grid::unit(arrow_size, "pt"),
                                   ends = "both", type = "closed"),
    start_cap = edge_cap,
    end_cap = edge_cap,
    edge_alpha = edge_alpha
  )

  ggdag <- ggplot2::ggplot(data = dag_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    edge_layers +
    ggplot2::geom_point(
      data = nodes,
      mapping = ggplot2::aes(x = x, y = y, colour = role, fill = role),
      shape = 21,
      size = node_size,
      stroke = 1,
      inherit.aes = FALSE
    )

  sides <- NULL
  if( label_placement == "node" ){
    labelled <- nodes[!is.na(nodes$label), , drop = FALSE]
    labelled$text_colour <- ifelse(as.character(labelled$role) == "latent",
                                   "black", "white")
    ggdag <- ggdag +
      ggplot2::geom_text(
        data = labelled,
        mapping = ggplot2::aes(x = x, y = y, label = label),
        colour = labelled$text_colour,
        size = label_size * 2 / 3,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }else{
    sides <- plot_label_sides(nodes, label_placement)
    for( side in c("above", "below") ){
      side_nodes <- nodes[sides == side & !is.na(nodes$display_label), , drop = FALSE]
      if( nrow(side_nodes) == 0L ) next
      ggdag <- ggdag + plot_node_label_layer(
        side_nodes,
        offset = if( side == "above" ) label_offset else -label_offset,
        vjust = if( side == "above" ) 0 else 1,
        label_size = label_size
      )
    }
  }

  room <- plot_label_room(nodes, sides, node_radius, label_offset, label_size)
  margin <- c(top = 4, right = 4, bottom = 4, left = 4) + room
  box_spacing <- 2
  if( legend$side != "none" ){
    box_spacing <- room[[legend$side]] + 2
    margin[[legend$side]] <- 4
  }

  ggdag <- ggdag +
    ggplot2::scale_colour_manual(values = col,
                                 name = "Group",
                                 breaks = present) +
    ggplot2::scale_fill_manual(values = fill_col,
                               name = "Group",
                               breaks = present) +
    ggplot2::coord_cartesian(clip = "off", expand = FALSE) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::guides(colour = ggplot2::guide_legend(
      nrow = if( legend$side %in% c("top", "bottom") && length(present) > 5L ) 2L else NULL,
      ncol = if( legend$side %in% c("left", "right") ) 1L else NULL,
      byrow = TRUE,
      override.aes = list(size = 7)
    )) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = legend$side,
      legend.justification = legend$justification,
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12, colour = "#27343E"),
      legend.key.width = grid::unit(19, "pt"),
      legend.spacing.x = grid::unit(8, "pt"),
      legend.box.spacing = grid::unit(box_spacing, "mm"),
      plot.margin = ggplot2::margin(margin[["top"]], margin[["right"]],
                                    margin[["bottom"]], margin[["left"]],
                                    unit = "mm"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )

  if( !is.null(file) ){
    size <- plot_fit_size(nodes, sides, node_radius, label_offset, label_size,
                          margin, box_spacing, legend$side, present)
    plot_save_file(ggdag, file, format,
                   width = if( is.null(width) ) size[["width"]] else width,
                   height = if( is.null(height) ) size[["height"]] else height,
                   dpi = dpi)

    return(invisible(ggdag))
  }

  return(ggdag)
}


#' Parse the single legend position argument
#'
#' @param legend_position A side, optionally followed by an alignment, such
#'   as `"top"`, `"top-left"` or `"right-bottom"`; or `"none"`.
#' @return A list with the ggplot2 legend `side` and `justification`.
#' @noRd
plot_legend_position <- function(legend_position){
  if( is.null(legend_position) ){
    legend_position <- "top"
  }
  if( !is.character(legend_position) || length(legend_position) != 1L ||
      is.na(legend_position) ){
    stop("plot_dagitty(): legend_position must be a single string such as ",
         "\"top\", \"top-left\" or \"none\".", call. = FALSE)
  }
  parts <- strsplit(tolower(trimws(legend_position)), "[-_ ]+")[[1L]]
  side <- if( length(parts) == 0L || !nzchar(parts[1L]) ) "top" else parts[1L]
  if( !side %in% c("top", "bottom", "left", "right", "none") ){
    stop("plot_dagitty(): legend_position must start with \"top\", ",
         "\"bottom\", \"left\", \"right\" or \"none\".", call. = FALSE)
  }
  alignment <- if( length(parts) >= 2L ) parts[2L] else "centre"
  allowed <- if( side %in% c("top", "bottom") ) c("left", "right") else c("top", "bottom")
  justification <- if( alignment %in% allowed ) alignment else "center"
  if( side == "none" ){
    justification <- "center"
  }
  position <- list(side = side, justification = justification)

  return(position)
}


#' Resolve the file format used to save a plot
#'
#' @param file Optional file path.
#' @param output_format Requested format.
#' @param format_missing Whether `output_format` was omitted.
#' @return `"svg"` or `"png"`, or `NULL` when no file is saved.
#' @noRd
plot_output_format <- function(file, output_format, format_missing){
  if( is.null(file) ){
    format <- NULL
    return(format)
  }
  if( !is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file) ){
    stop("plot_dagitty(): file must be NULL or a single file path.", call. = FALSE)
  }
  if( !is.character(output_format) || length(output_format) != 1L ||
      is.na(output_format) || !tolower(output_format) %in% c("svg", "png") ){
    stop("plot_dagitty(): output_format must be \"svg\" or \"png\".", call. = FALSE)
  }
  format <- tolower(output_format)
  base_name <- basename(file)
  extension <- if( grepl("[.][[:alnum:]]+$", base_name) ){
    tolower(sub("^.*[.]", "", base_name))
  }else{
    ""
  }
  if( format_missing && extension == "png" ){
    format <- "png"
  }
  if( nzchar(extension) && !extension %in% c("svg", "png") ){
    stop("plot_dagitty(): file must end in .svg or .png, or have no extension.",
         call. = FALSE)
  }
  if( extension %in% c("svg", "png") && extension != format ){
    stop("plot_dagitty(): the file extension '.", extension, "' does not match ",
         "output_format = \"", format, "\".", call. = FALSE)
  }

  return(format)
}


#' Group nodes into columns
#'
#' Neighbouring horizontal positions closer than 4% of the graph's width, and
#' closer than half the median gap between positions, form one column. Free-form
#' coordinates are grouped while evenly spaced columns stay separate.
#'
#' @param x Node horizontal coordinates.
#' @return An integer column number for each node, increasing from left to
#'   right.
#' @noRd
plot_columns <- function(x){
  span <- diff(range(x))
  gaps <- diff(sort(unique(x)))
  tolerance <- if( span > 0 ) min(0.04 * span, 0.5 * stats::median(gaps)) else 0
  ordering <- order(x)
  groups <- cumsum(c(TRUE, diff(x[ordering]) > tolerance))
  column <- integer(length(x))
  column[ordering] <- groups

  return(column)
}


#' Choose the height of automatic curves
#'
#' @param dag_df Node and edge table with node coordinates.
#' @return Three quarters of the median usable vertical gap within columns,
#'   otherwise half the median usable gap between columns. Returns `1` when
#'   neither gap is available.
#' @noRd
plot_curve_height <- function(dag_df){
  nodes <- dag_df[!duplicated(as.character(dag_df$name)), c("x", "y"), drop = FALSE]
  tolerance <- max(sqrt(.Machine$double.eps), 0.04 * diff(range(nodes$y)))
  column <- plot_columns(nodes$x)
  row_gaps <- unlist(lapply(split(nodes$y, column), function(values){
    gaps <- diff(sort(unique(values)))

    return(gaps)
  }), use.names = FALSE)
  row_gaps <- row_gaps[row_gaps > tolerance]
  if( length(row_gaps) > 0L ){
    height <- 0.75 * stats::median(row_gaps)
    return(height)
  }
  column_gaps <- diff(as.numeric(tapply(nodes$x, column, mean)))
  column_gaps <- column_gaps[column_gaps > sqrt(.Machine$double.eps)]
  if( length(column_gaps) > 0L ){
    height <- 0.5 * stats::median(column_gaps)
    return(height)
  }
  height <- 1

  return(height)
}


#' Choose the side of each node's label
#'
#' @param nodes One row per node with `x` and `y` coordinates.
#' @param placement `"above"`, `"below"` or `"auto"`.
#' @return A character vector of `"above"` or `"below"`, one per node.
#' @noRd
plot_label_sides <- function(nodes, placement){
  if( placement != "auto" ){
    sides <- rep(placement, nrow(nodes))
    return(sides)
  }
  column <- plot_columns(nodes$x)
  sides <- vapply(seq_len(nrow(nodes)), function(i){
    same <- column == column[i] & seq_len(nrow(nodes)) != i
    offsets <- nodes$y[same] - nodes$y[i]
    room_above <- if( any(offsets > 0) ) min(offsets[offsets > 0]) else Inf
    room_below <- if( any(offsets < 0) ) min(-offsets[offsets < 0]) else Inf
    side <- if( room_below > room_above ) "below" else "above"

    return(side)
  }, character(1L))

  return(sides)
}


#' Estimate label boxes in millimetres
#'
#' @param labels Displayed label text, possibly with line breaks.
#' @param label_size Label text size in millimetres.
#' @return A data frame with `width` and `height` in millimetres.
#' @noRd
plot_label_extent <- function(labels, label_size){
  missing <- is.na(labels)
  labels <- ifelse(missing, "", as.character(labels))
  lines <- strsplit(labels, "\n", fixed = TRUE)
  characters <- vapply(lines, function(text){
    value <- if( length(text) == 0L ) 0L else max(nchar(text))

    return(value)
  }, numeric(1L))
  line_count <- pmax(1L, lengths(lines))
  padding <- 2 * 0.18 * 1.2 * label_size
  extent <- data.frame(
    width = ifelse(missing, 0, 0.56 * label_size * characters + padding),
    height = ifelse(missing, 0, 1.2 * label_size * line_count + padding)
  )

  return(extent)
}


#' Find the space labels need beyond the outermost nodes
#'
#' @param nodes One row per node with `x`, `y` and `display_label`.
#' @param sides Label sides, or `NULL` when labels are inside the nodes.
#' @param node_radius,label_offset Node radius and label offset in
#'   millimetres.
#' @param label_size Label text size in millimetres.
#' @return Named margins in millimetres: `top`, `right`, `bottom`, `left`.
#' @noRd
plot_label_room <- function(nodes, sides, node_radius, label_offset, label_size){
  room <- c(top = node_radius, right = node_radius,
            bottom = node_radius, left = node_radius)
  if( is.null(sides) ){
    return(room)
  }
  extent <- plot_label_extent(nodes$display_label, label_size)
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(c(nodes$x, nodes$y)))
  top <- abs(nodes$y - max(nodes$y)) <= tolerance
  bottom <- abs(nodes$y - min(nodes$y)) <= tolerance
  left <- abs(nodes$x - min(nodes$x)) <= tolerance
  right <- abs(nodes$x - max(nodes$x)) <= tolerance
  above <- sides == "above"
  if( any(top & above) ){
    room[["top"]] <- max(room[["top"]], label_offset + extent$height[top & above])
  }
  if( any(bottom & !above) ){
    room[["bottom"]] <- max(room[["bottom"]],
                            label_offset + extent$height[bottom & !above])
  }
  room[["left"]] <- max(room[["left"]], extent$width[left] / 2)
  room[["right"]] <- max(room[["right"]], extent$width[right] / 2)

  return(room)
}


#' Fit a canvas to the plotted graph
#'
#' One scale, in millimetres per coordinate unit, keeps the layout's
#' proportions. Each pair of nodes needs enough horizontal room for their
#' labels or enough vertical room to stack a node and label, whichever needs
#' the smaller scale; the canvas uses the largest scale any pair needs.
#'
#' @param nodes One row per node with `x`, `y` and `display_label`.
#' @param sides Label sides, or `NULL` when labels are inside the nodes.
#' @param node_radius,label_offset Node radius and label offset in
#'   millimetres.
#' @param label_size Label text size in millimetres.
#' @param margin Plot margins in millimetres.
#' @param box_spacing Space between the legend and the graph in millimetres.
#' @param legend_side Legend side, or `"none"`.
#' @param legend_roles Roles listed in the legend.
#' @return Named `width` and `height` in inches, between 4 x 3 and 40 x 40.
#' @noRd
plot_fit_size <- function(nodes, sides, node_radius, label_offset, label_size,
                          margin, box_spacing, legend_side, legend_roles){
  extent <- plot_label_extent(nodes$display_label, label_size)
  if( is.null(sides) ){
    extent$width <- 0
    extent$height <- 0
  }
  node_width <- pmax(extent$width, 2 * node_radius)
  scale <- NA_real_
  if( nrow(nodes) > 1L ){
    pairs <- utils::combn(nrow(nodes), 2L)
    first <- pairs[1L, ]
    second <- pairs[2L, ]
    dx <- abs(nodes$x[first] - nodes$x[second])
    dy <- abs(nodes$y[first] - nodes$y[second])
    need_x <- (node_width[first] + node_width[second]) / 2 + 6
    need_y <- 2 * node_radius + (label_offset - node_radius) + 2 +
      pmax(extent$height[first], extent$height[second])
    scale_x <- ifelse(dx > 0, need_x / dx, Inf)
    scale_y <- ifelse(dy > 0, need_y / dy, Inf)
    pair_scale <- pmin(scale_x, scale_y)
    pair_scale <- pair_scale[is.finite(pair_scale)]
    if( length(pair_scale) > 0L ) scale <- max(pair_scale)
  }
  if( is.na(scale) ) scale <- 30
  # Legend text is 12 pt (4.2 mm); keys and spacing add about 16 mm.
  legend_width <- if( legend_side %in% c("left", "right") && length(legend_roles) > 0L ){
    0.56 * 4.2 * max(nchar(legend_roles)) + 16
  }else{
    0
  }
  legend_height <- if( legend_side %in% c("top", "bottom") && length(legend_roles) > 0L ){
    if( length(legend_roles) > 5L ) 20 else 11
  }else{
    0
  }
  spacing <- if( legend_side == "none" ) 0 else box_spacing
  extra_width <- margin[["left"]] + margin[["right"]] + legend_width +
    if( legend_side %in% c("left", "right") ) spacing else 0
  extra_height <- margin[["top"]] + margin[["bottom"]] + legend_height +
    if( legend_side %in% c("top", "bottom") ) spacing else 0
  x_span <- diff(range(nodes$x))
  y_span <- diff(range(nodes$y))
  # Keep the larger side within 40 inches by reducing the scale if needed.
  limit <- 40 * 25.4
  if( x_span > 0 ) scale <- min(scale, max(0, limit - extra_width) / x_span)
  if( y_span > 0 ) scale <- min(scale, max(0, limit - extra_height) / y_span)
  width <- x_span * scale + extra_width
  height <- y_span * scale + extra_height
  # A top or bottom legend needs its widest row; a side legend its full height.
  if( length(legend_roles) > 0L ){
    item_width <- 0.56 * 4.2 * max(nchar(legend_roles)) + 13
    per_row <- if( length(legend_roles) > 5L ) ceiling(length(legend_roles) / 2) else length(legend_roles)
    if( legend_side %in% c("top", "bottom") ){
      width <- max(width, per_row * item_width + 8)
    }
    if( legend_side %in% c("left", "right") ){
      height <- max(height, 9 * length(legend_roles) + 8)
    }
  }
  size <- c(width = min(40, max(4, width / 25.4)),
            height = min(40, max(3, height / 25.4)))

  return(size)
}


#' Save a plot as SVG or PNG
#'
#' @param plot A ggplot object.
#' @param file Output file path.
#' @param format `"svg"` or `"png"`.
#' @param width,height Size in inches.
#' @param dpi PNG resolution.
#' @return The file path, invisibly.
#' @noRd
plot_save_file <- function(plot, file, format, width, height, dpi){
  if( format == "svg" ){
    if( plot_package_available("svglite") ){
      device <- svglite::svglite
    }else if( isTRUE(capabilities("cairo")) ){
      device <- grDevices::svg
    }else{
      stop("plot_dagitty(): saving SVG files requires the svglite package or ",
           "an R installation with cairo support.", call. = FALSE)
    }
    ggplot2::ggsave(file, plot, device = device, width = width,
                    height = height, units = "in", bg = "white",
                    limitsize = FALSE)
  }else{
    dpi <- min(dpi, floor(10000 / max(width, height)))
    if( plot_package_available("ragg") ){
      ggplot2::ggsave(file, plot, device = ragg::agg_png, width = width,
                      height = height, units = "in", dpi = dpi, bg = "white",
                      limitsize = FALSE)
    }else{
      ggplot2::ggsave(file, plot, device = "png", width = width,
                      height = height, units = "in", dpi = dpi, bg = "white",
                      limitsize = FALSE)
    }
  }

  return(invisible(file))
}


#' Check whether an optional graphics package is installed
#'
#' @param package Package name.
#' @return `TRUE` if the package can be loaded.
#' @noRd
plot_package_available <- function(package){
  available <- requireNamespace(package, quietly = TRUE)

  return(available)
}


#' Label nodes a fixed physical distance above or below them
#'
#' A `GeomLabel` whose labels are shifted by `offset` millimetres, so their
#' distance from the nodes does not depend on the device size.
#'
#' @noRd
GeomDagNodeLabel <- ggplot2::ggproto(
  "GeomDagNodeLabel",
  ggplot2::GeomLabel,
  # Accept exactly the parameters of the installed ggplot2's GeomLabel, plus
  # the offset, so the label options track the ggplot2 version.
  parameters = function(self, extra = FALSE){
    params <- union(ggplot2::GeomLabel$parameters(extra), "offset")

    return(params)
  },
  draw_panel = function(self, data, panel_params, coord, ..., offset = 0){
    grob <- ggplot2::ggproto_parent(ggplot2::GeomLabel, self)$draw_panel(
      data, panel_params, coord, ...
    )
    grob$vp <- grid::viewport(y = grid::unit(0.5, "npc") + grid::unit(offset, "mm"))

    return(grob)
  }
)


#' Create a node label layer
#'
#' @param nodes Nodes to label, with `x`, `y` and `display_label`.
#' @param offset Vertical label shift in millimetres.
#' @param vjust `0` for labels above nodes or `1` for labels below them.
#' @param label_size Label text size in millimetres.
#' @return A ggplot2 layer.
#' @noRd
plot_node_label_layer <- function(nodes, offset, vjust, label_size){
  params <- list(
    colour = "black",
    fill = "white",
    size = label_size,
    vjust = vjust,
    hjust = 0.5,
    lineheight = 1.05,
    label.padding = grid::unit(0.18, "lines"),
    label.r = grid::unit(0, "lines"),
    offset = offset,
    na.rm = FALSE
  )
  # ggplot2 4.0 replaced label.size with the linewidth aesthetic.
  if( "label.size" %in% ggplot2::GeomLabel$parameters() ){
    params$label.size <- 0.25
  }else{
    params$linewidth <- 0.25
  }
  layer_out <- ggplot2::layer(
    data = nodes,
    mapping = ggplot2::aes(x = x, y = y, label = display_label),
    stat = "identity",
    geom = GeomDagNodeLabel,
    position = "identity",
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = params
  )

  return(layer_out)
}


#' Create identifiers for rendered graph edges
#'
#' @param edges Node-edge rows with `name` and `to` endpoint columns.
#' @param kind Edge kind: `"directed"`, `"bidirected"`, or `"undirected"`.
#' @return A character identifier for each edge.
#' @noRd
edge_path_ids <- function(edges, kind){
  if( nrow(edges) == 0L ) {
    path_ids <- character()
    return(path_ids)
  }
  path_ids <- paste(
    kind,
    as.character(edges$name),
    as.character(edges$to),
    sep = "\036"
  )

  return(path_ids)
}


#' Sample straight graph edge paths
#'
#' @param edges Edge rows with endpoint coordinates `x`, `y`, `xend`, and
#'   `yend`.
#' @param kind Edge kind: `"directed"`, `"bidirected"`, or `"undirected"`.
#' @return An edge-point data frame with path indices, identifiers, endpoints,
#'   and kinds.
#' @noRd
straight_edge_paths <- function(edges, kind){
  progress <- if( kind == "directed" ) seq(0, 1, length.out = 17L) else c(0, 1)
  paths <- lapply(seq_len(nrow(edges)), function(i){
    path <- edges[rep(i, length(progress)), , drop = FALSE]
    path$x <- edges$x[i] + progress * (edges$xend[i] - edges$x[i])
    path$y <- edges$y[i] + progress * (edges$yend[i] - edges$y[i])
    path$index <- progress
    path$edge_id <- edge_path_ids(edges[i, , drop = FALSE], kind)
    path$from <- as.character(edges$name[i])
    path$to <- as.character(edges$to[i])
    path$kind <- kind

    return(path)
  })
  if( length(paths) > 0L ){
    paths_out <- do.call(rbind, paths)
    return(paths_out)
  }

  paths <- edges
  paths$index <- numeric()
  paths$edge_id <- character()
  paths$from <- character()
  paths$to <- character()
  paths$kind <- character()

  return(paths)
}


#' Sample curved paths for bidirected edges
#'
#' @param edges Edge rows with endpoint coordinates `x`, `y`, `xend`, and
#'   `yend`.
#' @param strength Amount of rotation used to place the cubic control points.
#' @param n Number of sampled points per edge.
#' @return An edge-point data frame with cubic Bezier path coordinates.
#' @noRd
bidirected_edge_paths <- function(edges, strength = 0.3, n = 101L){
  progress <- seq(0, 1, length.out = n)
  paths <- lapply(seq_len(nrow(edges)), function(i){
    start <- c(as.numeric(edges$x[i]), as.numeric(edges$y[i]))
    end <- c(as.numeric(edges$xend[i]), as.numeric(edges$yend[i]))
    vector <- end - start
    half_length <- sqrt(sum(vector^2)) / 2
    edge_angle <- atan2(vector[2L], vector[1L])
    turn <- -pi / 2 * strength
    control_1 <- start + half_length * c(
      cos(edge_angle - turn), sin(edge_angle - turn)
    )
    control_2 <- end + half_length * c(
      cos(edge_angle - pi + turn), sin(edge_angle - pi + turn)
    )
    complement <- 1 - progress
    coordinates <-
      complement^3 %o% start +
      (3 * complement^2 * progress) %o% control_1 +
      (3 * complement * progress^2) %o% control_2 +
      progress^3 %o% end

    path <- edges[rep(i, n), , drop = FALSE]
    path$x <- coordinates[, 1L]
    path$y <- coordinates[, 2L]
    path$index <- progress
    path$edge_id <- edge_path_ids(edges[i, , drop = FALSE], "bidirected")
    path$from <- as.character(edges$name[i])
    path$to <- as.character(edges$to[i])
    path$kind <- "bidirected"

    return(path)
  })
  if( length(paths) > 0L ){
    paths_out <- do.call(rbind, paths)
    return(paths_out)
  }

  paths <- edges
  paths$index <- numeric()
  paths$edge_id <- character()
  paths$from <- character()
  paths$to <- character()
  paths$kind <- character()

  return(paths)
}


#' Render edge paths using stored control points
#'
#' @param dag_df Node and edge table with endpoint coordinates.
#' @param controls Native control-point table with `v`, `e`, `w`, `x`, and `y`
#'   columns.
#' @return An edge-point data frame ordered by edge and path position.
#' @noRd
stored_edge_paths <- function(dag_df, controls){
  direction <- as.character(dag_df$direction)
  directed <- dag_df[!is.na(direction) & direction == "->", , drop = FALSE]
  bidirected <- dag_df[!is.na(direction) & direction == "<->", , drop = FALSE]
  undirected <- dag_df[!is.na(direction) & direction == "--", , drop = FALSE]
  paths <- as.data.frame(data.table::rbindlist(list(
    straight_edge_paths(directed, "directed"),
    bidirected_edge_paths(bidirected),
    straight_edge_paths(undirected, "undirected")
  ), use.names = TRUE, fill = TRUE))
  if( nrow(paths) == 0L ) return(paths)

  edge_rows <- dag_df[!is.na(direction), , drop = FALSE]
  keys <- edge_control_keys(as.character(edge_rows$name),
                             as.character(edge_rows$direction),
                             as.character(edge_rows$to))
  control_keys <- edge_control_keys(controls$v, controls$e, controls$w)
  matched <- match(keys, control_keys)
  curved <- which(!is.na(matched))
  if( length(curved) > 0L ){
    native_paths <- lapply(curved, function(i){
      edge <- edge_rows[i, , drop = FALSE]
      control <- controls[matched[i], , drop = FALSE]
      kind <- switch(as.character(edge$direction),
                     "->" = "directed", "<->" = "bidirected", "--" = "undirected")
      t <- seq(0, 1, length.out = 201L)
      path <- edge[rep(1L, length(t)), , drop = FALSE]
      path$x <- (1 - t)^2 * edge$x + 2 * t * (1 - t) * control$x + t^2 * edge$xend
      path$y <- (1 - t)^2 * edge$y - 2 * t * (1 - t) * control$y + t^2 * edge$yend
      path$index <- t
      path$edge_id <- edge_path_ids(edge, kind)
      path$from <- as.character(edge$name)
      path$to <- as.character(edge$to)
      path$kind <- kind

      return(path)
    })
    replacements <- as.data.frame(data.table::rbindlist(native_paths, use.names = TRUE, fill = TRUE))
    paths <- as.data.frame(data.table::rbindlist(list(
      paths[!paths$edge_id %in% replacements$edge_id, , drop = FALSE], replacements
    ), use.names = TRUE, fill = TRUE))
  }
  edge_order <- edge_path_ids(directed, "directed")
  edge_order <- c(edge_order, edge_path_ids(bidirected, "bidirected"),
                 edge_path_ids(undirected, "undirected"))
  paths$.edge_group <- match(paths$edge_id, edge_order)
  paths_out <- paths[order(paths$.edge_group, paths$index), , drop = FALSE]

  return(paths_out)
}


#' Add plot anchors when all nodes share an x or y coordinate
#'
#' @param dag_df Node and edge table with endpoint coordinates.
#' @return A two-row `x`, `y` data frame, or `NULL` when both coordinate spans
#'   are positive.
#' @noRd
curve_plot_anchors <- function(dag_df){
  nodes <- dag_df[!duplicated(dag_df$name), c("x", "y"), drop = FALSE]
  x_bounds <- range(nodes$x)
  y_bounds <- range(nodes$y)
  x_span <- diff(x_bounds)
  y_span <- diff(y_bounds)
  if( x_span > 0 && y_span > 0 ) {
    out <- NULL
    return(out)
  }
  if( x_span == 0 ) x_bounds <- mean(x_bounds) + c(-1, 1) * curve_axis_gap(nodes$y)
  if( y_span == 0 ) y_bounds <- mean(y_bounds) + c(-0.25, 0.25) * curve_axis_gap(nodes$x)
  anchors <- data.frame(x = x_bounds, y = y_bounds)

  return(anchors)
}


#' Create a capped graph edge path layer
#'
#' @param data Sampled edge paths with grouping and path-index columns.
#' @param arrow Arrow specification for the path layer.
#' @param start_cap Shape used to cap the start of each path.
#' @param end_cap Shape used to cap the end of each path.
#' @param edge_alpha Edge opacity.
#' @param edge_linetype Edge line type.
#' @return A ggraph edge-path layer.
#' @noRd
dag_edge_path_layer <- function(data,
                                 arrow,
                                 start_cap,
                                 end_cap,
                                 edge_alpha = 1,
                                 edge_linetype = "solid"
                                 ){
  layer_out <- ggplot2::layer(
    data = data,
    mapping = ggplot2::aes(group = .edge_group, index = index),
    stat = "identity",
    geom = ggraph::GeomEdgePath,
    position = "identity",
    show.legend = NA,
    inherit.aes = TRUE,
    check.aes = FALSE,
    params = list(
      arrow = arrow,
      interpolate = FALSE,
      na.rm = TRUE,
      start_cap = start_cap,
      end_cap = end_cap,
      edge_width = 0.35,
      edge_colour = "#647484",
      edge_alpha = edge_alpha,
      edge_linetype = edge_linetype,
      lineend = "butt",
      linejoin = "round",
      linemitre = 1
    )
  )

  return(layer_out)
}


#' Create plot layers for all graph edge kinds
#'
#' @param rendered_paths Sampled edge paths with a `kind` column.
#' @param edge_plan List containing optional plot range anchors.
#' @param arrow_directed Arrow specification for directed edges.
#' @param arrow_bidirected Arrow specification for bidirected edges.
#' @param start_cap Shape used to cap the start of each path.
#' @param end_cap Shape used to cap the end of each path.
#' @param edge_alpha Edge opacity.
#' @return A list of edge-path layers and any required range-anchor layer.
#' @noRd
plot_dag_edge_layers <- function(rendered_paths,
                                  edge_plan,
                                  arrow_directed,
                                  arrow_bidirected,
                                  start_cap,
                                  end_cap,
                                  edge_alpha = 1
                                  ){
  layers <- list()
  directed <- rendered_paths[rendered_paths$kind == "directed", , drop = FALSE]
  if( nrow(directed) > 0L ){
    layers[[length(layers) + 1L]] <- dag_edge_path_layer(
      data = directed,
      arrow = arrow_directed,
      start_cap = start_cap,
      end_cap = end_cap,
      edge_alpha = edge_alpha
    )
  }
  bidirected <- rendered_paths[
    rendered_paths$kind == "bidirected", , drop = FALSE
  ]
  if( nrow(bidirected) > 0L ){
    layers[[length(layers) + 1L]] <- dag_edge_path_layer(
      data = bidirected,
      arrow = arrow_bidirected,
      start_cap = start_cap,
      end_cap = end_cap,
      edge_alpha = edge_alpha,
      edge_linetype = "dashed"
    )
  }
  undirected <- rendered_paths[
    rendered_paths$kind == "undirected", , drop = FALSE
  ]
  if( nrow(undirected) > 0L ){
    layers[[length(layers) + 1L]] <- dag_edge_path_layer(
      data = undirected,
      arrow = NULL,
      start_cap = start_cap,
      end_cap = end_cap,
      edge_alpha = edge_alpha
    )
  }
  if( !is.null(edge_plan$anchors) ){
    layers[[length(layers) + 1L]] <- ggplot2::geom_blank(
      data = edge_plan$anchors,
      mapping = ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE
    )
  }


  return(layers)
}


#' Prepare node and edge rows for plotting
#'
#' Converts the graph with `ggdag::tidy_dagitty()` and reverses its y
#' coordinates for plotting. When labels are supplied, `add_labels()`
#' adds the node labels and roles.
#'
#' @importFrom ggdag tidy_dagitty
#' @param dag dagitty object
#' @param labels vector of labels for nodes in dagitty object
#' @param roles Optional precomputed node roles. When labels are supplied,
#'   roles are computed if this is NULL.
#' @return A data frame of node and edge rows, with label and role columns
#'   when labels are supplied.
#' @noRd
tidy_ggdagitty <- function(dag, labels = NULL, roles = NULL){
  # Cleaning the dags and turning it into a data frame.
  dag_df <- data.frame(ggdag::tidy_dagitty(dag))

  # flip y axis for ggplot
  dag_df$y <- dag_df$y*-1
  dag_df$yend <- dag_df$yend*-1

  if( is.null(labels) ){

    return(dag_df)
  }

  dag_df <- add_labels(dag, dag_df, labels, roles = roles)

  return(dag_df)
}


#' Add node labels and causal roles to a plotting table
#'
#' Matches labels by node name and assigns the first matching role
#' from the role list.
#'
#' @param dag dagitty object
#' @param dag_df a dag object converted to data frame using the ggdag::tidy_dagitty() function, or similar.
#' @param labels Named node labels matched to the table's name column.
#' @param roles Optional precomputed get_roles(dag) output. Computed internally when NULL.
#' @return The input data frame with label and role columns.
#' @noRd
add_labels <- function(dag, dag_df, labels, roles = NULL){

  # Labeling variables
  dag_df$label <- sapply( seq_along( dag_df$name ),
                          function(x) if( dag_df$name[x] %in% attr(labels, "names")) {
                            value <- labels[[ dag_df$name[x] ]]

                            return(value)
                          }else{
                            value <- NULL
                            return(invisible(value))
                          } )

  if( is.null(roles) ){
    roles <- classify_node_roles(dag)
  }

  node_roles <- roles
  dag_df$role <- vapply(as.character(dag_df$name), function(node){
    matches <- names(node_roles)[vapply(
      node_roles, function(role_nodes) {
        matches <- node %in% role_nodes

        return(matches)
      }, logical(1L)
    )]
    if( length(matches) == 0L ) {
      value <- NA_character_

      return(value)
    } else {
      value <- matches[1L]

      return(value)
    }
  }, character(1L), USE.NAMES = FALSE)

  return(dag_df)
}
