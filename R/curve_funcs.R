# ==============================================================================
# curve_funcs.R


# add_curves() and the routing calculations used to choose control points
# ==============================================================================


#' Curve the edges in a dagitty object
#'
#' Curves directed arrows to avoid nodes and reduce edge crossings, storing
#' their control points in the dagitty object.
#'
#' @param dag A dagitty object.
#' @param height Single nonnegative finite number controlling arc height.
#'   Interpreted using `mode`. Defaults to `0.2`. Use `0` to straighten
#'   directed edges.
#' @param mode In `"relative"` mode (the default), the starting curve height is
#'   scaled to the edge and graph coordinates as described below. A taller curve
#'   may be used to clear nodes. In `"absolute"` mode, the curve departs from the
#'   straight edge by exactly `height` coordinate units and is not enlarged for
#'   clearance.
#' @param position Curve position for a layout running from left to right. `"default"`
#'   chooses positions automatically. `"flipped"` mirrors each default arc,
#'   keeping its height. `"over"` and `"under"` curve above or below the straight
#'   edge. Defaults to `"default"`.
#'
#' @details
#' Complete existing node coordinates are kept. If any coordinates are missing,
#' [add_coords()] regenerates the entire layout using its defaults. Call
#' [add_coords()] first to choose the layout settings.
#' This function does not change nodes, edges or node roles, and is not called
#' by [build_graph()].
#'
#' Existing bends on directed edges are replaced. Edges that need no bend stay
#' straight. Stored control points on other edge types are kept. Curve placement
#' takes nodes and edges into account, but does not consider labels or the lines
#' connecting labels to nodes. Specifying a curve position or a small absolute
#' height can leave overlaps.
#'
#' For an edge that is not vertical, height is the largest vertical departure from the
#' straight edge, reached halfway along a centred quadratic Bezier curve. For a
#' horizontal edge at y = 0, absolute height 1 reaches y = 1 or y = -1. The
#' control point is twice this distance from the straight edge. For a vertical
#' edge, height is the largest horizontal departure.
#' Height is not a physical distance or a fixed fraction of the final plot.
#' Plot limits, aspect ratio and device size affect its appearance.
#'
#' In relative mode, let `s = xspan / xref` for an edge that is not vertical. Its
#' starting height is
#' `height * yref * (s + 0.25 * s * (1 - s))`. A vertical edge uses
#' `s = yspan / yref` and the corresponding horizontal reference `xref`.
#' `xspan` and `yspan` are the edge's absolute horizontal and vertical spans.
#' `xref` and `yref` normally equal the graph's horizontal and vertical spans
#' across all node coordinates.
#' If the vertical span is zero, `yref` becomes half the median positive
#' horizontal gap between nodes. If the horizontal span is zero, `xref` becomes
#' twice the median positive vertical gap. Reversing an arrow does not change
#' its span.
#' In relative mode, the curve height can increase up to the larger of its
#' starting height and the reference span in the direction of the bend.
#'
#' The graph stores the curve control points, but not the `height`, `mode` or
#' `position` settings.
#' A later call to [add_coords()] moves curved edges with their endpoints and can
#' change their absolute height. Call this function again with `mode = "absolute"`
#' to set an exact height after moving nodes. Curving coincident endpoints, or
#' moving an edge's stored control point without usable endpoint coordinates,
#' raises an error.
#'
#' To display these curved edges without recalculating their paths, use
#' `plot_dagitty(dag, curve = FALSE)`. In dagitty 0.3-4, reconstructing
#' graph text or assigning coordinates outside causaliflower can round control
#' points or discard a control point whose x coordinate is zero.
#'
#' @return A dagitty object with edge control points stored in its native edge
#'   `pos` fields.
#' @seealso [add_coords()], [plot_dagitty()]
#' @examples
#' dag <- dagitty::dagitty("dag { T -> M -> Y; T -> Y }")
#' dag <- add_curves(dag, position = "flipped")
#' plot_dagitty(dag, curve = FALSE)
#'
#' dag <- add_curves(dag, height = 1, mode = "absolute")
#' plot_dagitty(dag, curve = FALSE)
#' @export

add_curves <- function(dag,
                       height = 0.2,
                       mode = "relative",
                       position = "default"
                       ){
  if( !inherits(dag, "dagitty") ){
    stop("add_curves(): dag must be a dagitty object.", call. = FALSE)
  }
  if( !is.numeric(height) || is.complex(height) || length(height) != 1L ||
      !is.finite(height) || height < 0 ){
    stop("add_curves(): height must be a single non-negative finite number.",
         call. = FALSE)
  }
  if( !is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% c("relative", "absolute") ){
    stop("add_curves(): mode must be 'relative' or 'absolute'.", call. = FALSE)
  }
  if( !is.character(position) || length(position) != 1L || is.na(position) ||
      !position %in% c("default", "flipped", "over", "under") ){
    stop("add_curves(): position must be 'default', 'flipped', 'over' or 'under'.",
         call. = FALSE)
  }
  height <- as.numeric(height)
  mode <- unname(mode)
  position <- unname(position)

  edges <- dag_edges(dag, "add_curves()")
  directed <- edges[edges$e == "->", , drop = FALSE]
  node_names <- names(dag)
  if( length(node_names) == 0L ) return(dag)

  coordinates <- dagitty::coordinates(dag)
  missing_coordinates <- !all(node_names %in% names(coordinates$x)) ||
    !all(node_names %in% names(coordinates$y)) ||
    anyNA(coordinates$x[node_names]) || anyNA(coordinates$y[node_names])
  empty_controls <- data.frame(v = character(), e = character(), w = character(),
                               x = numeric(), y = numeric())
  if( missing_coordinates ){
    # Remove directed control points before regenerating the layout.
    dag <- set_edge_controls(dag, empty_controls, replace = directed,
                              caller = "add_curves()")
    dag <- add_coords(dag)
  }
  coordinates <- dagitty::coordinates(dag)
  if( !all(node_names %in% names(coordinates$x)) ||
      !all(node_names %in% names(coordinates$y)) ||
      any(!is.finite(coordinates$x[node_names])) ||
      any(!is.finite(coordinates$y[node_names])) ){
    stop("add_curves(): node coordinates must be finite.", call. = FALSE)
  }
  if( nrow(directed) == 0L ) return(dag)
  if( height == 0 ){
    dag_out <- set_edge_controls(dag, empty_controls, replace = directed,
                                  caller = "add_curves()")
    return(dag_out)
  }

  dag_df <- tidy_ggdagitty(dag_for_layout(dag, "add_curves"))
  retained_controls <- edge_controls(dag, "add_curves()")
  retained_controls <- retained_controls[retained_controls$e != "->", , drop = FALSE]
  stored_paths <- stored_edge_paths(dag_df, retained_controls)
  other_paths <- stored_paths[stored_paths$kind != "directed", , drop = FALSE]
  plan <- checked_curve_edge_plan(dag_df, height, mode, position, other_paths)

  curved <- Filter(function(control) {
    matches <- isTRUE(control$curved)

    return(matches)
  }, plan$controls)
  controls <- if( length(curved) == 0L ) empty_controls else {
    do.call(rbind, lapply(curved, function(control){
      values <- data.frame(v = control$from, e = "->", w = control$to,
                           x = control$control_x, y = -control$control_y)

      return(values)
    }))
  }
  dag_out <- set_edge_controls(dag, controls, replace = directed,
                                caller = "add_curves()")

  return(dag_out)
}


#' Check coordinate geometry and construct a curve plan
#'
#' @param dag_df Node and edge table with endpoint coordinates.
#' @param height Positive starting curve height.
#' @param mode Curve height mode, `"relative"` or `"absolute"`.
#' @param position Curve position, `"default"`, `"flipped"`, `"over"`, or
#'   `"under"`.
#' @param other_paths Optional sampled paths for fixed non-directed edges.
#' @param caller Name of the calling function for error messages.
#' @return A curve plan with finite sampled paths and control points.
#' @noRd
checked_curve_edge_plan <- function(dag_df,
                                     height = 0.2,
                                     mode = "relative",
                                     position = "default",
                                     other_paths = NULL,
                                     caller = "add_curves"
                                     ){
  caller <- sub("\\(\\)$", "", caller)
  fail <- function(message){
    stop(caller, "(): ", message, call. = FALSE)
  }
  if( any(!is.finite(dag_df$x)) || any(!is.finite(dag_df$y)) ){
    fail("node coordinates must be finite.")
  }
  arrow <- dag_df[!is.na(dag_df$direction) & dag_df$direction == "->", , drop = FALSE]
  dx <- arrow$xend-arrow$x
  dy <- arrow$yend-arrow$y
  if( any(dx == 0 & dy == 0) ){
    fail("cannot curve an edge whose endpoints have identical coordinates.")
  }
  x_span <- diff(range(dag_df$x))
  y_span <- diff(range(dag_df$y))
  x_reference <- if(x_span > 0) x_span else 2*curve_axis_gap(dag_df$y)
  y_reference <- if(y_span > 0) y_span else curve_axis_gap(dag_df$x)/2
  if( any(!is.finite(c(x_reference, y_reference))) ||
      x_reference <= 0 || y_reference <= 0 ){
    fail("the requested height cannot be represented in this coordinate layout.")
  }
  lateral <- curve_is_lateral(arrow$x, arrow$xend)
  base_height <- if(mode == "absolute") rep(height, nrow(arrow)) else {
    curve_relative_base_height(dx, dy, x_reference, y_reference,
                                height, lateral)
  }
  normalised_height <- base_height/ifelse(lateral, x_reference, y_reference)
  if( any(!is.finite(4*base_height)) || any(base_height <= 0) ||
      any(!is.finite(4*normalised_height)) || any(normalised_height <= 0) ){
    fail("the requested height cannot be represented in this coordinate layout.")
  }
  plan <- curve_edge_plan(dag_df, height, mode, position, other_paths)
  finite_paths <- all(is.finite(plan$paths$x)) && all(is.finite(plan$paths$y))
  finite_controls <- all(vapply(plan$controls, function(control){
    if( !all(is.finite(c(control$control_x, control$control_y, control$height))) ){
      matches <- FALSE
      return(matches)
    }
    if( !control$curved ) {
      matches <- TRUE
      return(matches)
    }
    axis <- if(control$lateral) "x" else "y"
    endpoints <- plan$nodes[match(c(control$from, control$to), plan$nodes$name), axis]
    midpoint <- endpoints[1L]/2+endpoints[2L]/2
    matches <- control[[paste0("control_", axis)]] != midpoint

    return(matches)
  }, logical(1L)))
  if( length(plan$errors) > 0L || !finite_paths || !finite_controls ){
    fail("the requested height cannot be represented in this coordinate layout.")
  }
  return(plan)
}


#' Scale relative curve height by edge span
#'
#' @param span_proportion Edge span divided by the graph reference span.
#' @return The numeric span factors.
#' @noRd
curve_span_factor <- function(span_proportion){
  factor <- span_proportion + 0.25 * span_proportion * (1 - span_proportion)

  return(factor)
}


#' Identify edges with matching horizontal coordinates
#'
#' @param x Starting horizontal coordinates.
#' @param xend Ending horizontal coordinates.
#' @return Logical values identifying edges that bend horizontally.
#' @noRd
curve_is_lateral <- function(x, xend){
  lateral <- abs(xend - x) <= 1e-12 * pmax(1, abs(x), abs(xend))

  return(lateral)
}


#' Compute starting relative curve heights
#'
#' @param dx Horizontal endpoint differences.
#' @param dy Vertical endpoint differences.
#' @param xref Horizontal reference span.
#' @param yref Vertical reference span.
#' @param height Non-negative starting curve height.
#' @param lateral Logical values identifying edges that bend horizontally.
#' @return Numeric curve heights in graph coordinate units.
#' @noRd
curve_relative_base_height <- function(dx, dy, xref, yref, height, lateral){
  span_proportion <- ifelse(lateral, abs(dy) / yref, abs(dx) / xref)
  bend_reference <- ifelse(lateral, xref, yref)
  height_out <- height * bend_reference * curve_span_factor(span_proportion)

  return(height_out)
}


# Route selection uses graph coordinates and paths, independently of labels.
#' Choose directed edge paths that avoid nodes and reduce crossings
#'
#' @param df Node and edge table with endpoint coordinates.
#' @param height Non-negative starting curve height.
#' @param mode Curve height mode: `"relative"` or `"absolute"`.
#' @param position Curve position: `"default"`, `"flipped"`, `"over"`, or
#'   `"under"`.
#' @param other_paths Optional sampled paths for the fixed non-directed edges.
#' @return A curve plan containing sampled paths, control points, nodes, and
#'   reference spans.
#' @noRd
curve_edge_plan <- function(df, height = 0.2, mode = "relative",
                             position = "default", other_paths = NULL){
  initial <- initial_curve_plan(df, height, mode, "default")
  nodes <- initial$nodes
  origin <- c(min(nodes$x), min(nodes$y))
  scale <- c(initial$xref, initial$yref)
  xy <- sweep(sweep(as.matrix(nodes[, c("x", "y")]), 2, origin), 2, scale, "/")
  rownames(xy) <- nodes$name
  edges <- df[!is.na(df$direction) & as.character(df$direction) == "->", , drop = FALSE]
  n <- nrow(edges)
  if( !n ) return(initial)
  # Clearance is measured in normalised graph coordinates.
  nearest <- as.matrix(stats::dist(xy))
  diag(nearest) <- Inf
  margin <- min(.1, .35 * stats::median(apply(nearest, 1, min)))
  if( !is.finite(margin) || margin <= 0 ) margin <- .05
  states <- lapply(seq_len(n), function(i){
    control <- initial$controls[[i]]
    a <- xy[as.character(edges$name[i]),]
    b <- xy[as.character(edges$to[i]),]
    values <- list(a = unname(a), b = unname(b), from = as.character(edges$name[i]), to = as.character(edges$to[i]),
                   lateral = control$lateral, h = control$height / scale[if( control$lateral ) 1L else 2L],
                   sign = control$sign)

    return(values)
  })
  coordinates <- function(s, t = seq(0, 1, length.out = 33L)){ # Sample one candidate path in normalised coordinates
    out <- (1 - t) %o% s$a + t %o% s$b
    axis <- if( s$lateral ) 1L else 2L
    out[, axis] <- out[, axis] + (if( s$lateral ) -1 else 1) * s$sign * 4 * s$h * t * (1 - t)

    return(out)
  }
  next_geometry_id <- 0L # Unique IDs for candidate curves within this routing call
  cache_geometry <- function(s){ # Cache sampled points and bounding boxes for one candidate curve
    next_geometry_id <<- next_geometry_id + 1L
    s$geometry_id <- next_geometry_id
    s$points <- if( !is.null(s$points) ) s$points else coordinates(s, seq(0, 1, length.out = 33L))
    p <- s$points
    s$low <- c(min(p[, 1]), min(p[, 2]))
    s$high <- c(max(p[, 1]), max(p[, 2]))
    if( !is.null(s[["h"]]) && s[["h"]] > 0 ){
      axis <- if( s$lateral ) 1L else 2L
      bend <- (if( s$lateral ) -1 else 1) * s$sign * s$h
      t <- .5 + (s$b[axis] - s$a[axis]) / (8 * bend)
      if( t > 0 && t < 1 ){
        extreme <- coordinates(s, t)[1,]
        s$low <- pmin(s$low, extreme)
        s$high <- pmax(s$high, extreme)
      }
    }
    s$segment_low <- pmin(p[-nrow(p), , drop = FALSE], p[-1L, , drop = FALSE])
    s$segment_high <- pmax(p[-nrow(p), , drop = FALSE], p[-1L, , drop = FALSE])

    return(s)
  }
  node_score <- function(s){ # Score a path by nearby nodes, then by the squared clearance shortfall
    others <- xy[!rownames(xy) %in% c(s$from, s$to), , drop = FALSE]
    if( !nrow(others) ) {
      score <- c(0, 0)
      return(score)
    }
    path <- coordinates(s)
    a <- path[-nrow(path), , drop = FALSE]
    v <- diff(path)
    vv <- rowSums(v * v)
    distances <- apply(others, 1, function(point){
      w <- matrix(point, nrow(a), 2, byrow = TRUE) - a
      t <- pmax(0, pmin(1, rowSums(w * v) / pmax(vv, 1e-30)))
      value <- min(sqrt(rowSums((w - v * t)^2)))

      return(value)
    })
    score <- c(sum(distances < margin - 1e-8), sum((pmax(0, margin - distances) / margin)^2))

    return(score)
  }
  y_at <- function(s, x){
    t <- (x - s$a[1]) / (s$b[1] - s$a[1])
    y <- s$a[2] + t * (s$b[2] - s$a[2]) + s$sign * 4 * s$h * t * (1 - t)

    return(y)
  }
  shared_point <- function(s, r, p){
    shared <- intersect(c(s$from, s$to), c(r$from, r$to))
    shared_out <- length(shared) > 0 && any(vapply(shared, function(k) {
      matches <- sum((p - xy[k,])^2) < 1e-14

      return(matches)
    }, logical(1)))

    return(shared_out)
  }
  cross <- function(a, b) {
    cross_product <- a[1] * b[2] - a[2] * b[1]

    return(cross_product)
  }
  poly_contact <- function(s, r){
    p <- s$points
    q <- r$points
    for(i in seq_len(nrow(p) - 1L)){
      a <- p[i,]
      v <- p[i + 1L,] - a
      low <- pmin(p[i,], p[i + 1L,])
      high <- pmax(p[i,], p[i + 1L,])
      possible <- which(r$segment_low[, 1] <= high[1] + 1e-10 &
        r$segment_high[, 1] >= low[1] - 1e-10 &
        r$segment_low[, 2] <= high[2] + 1e-10 &
        r$segment_high[, 2] >= low[2] - 1e-10)
      for(j in possible){
        b <- q[j,]
        w <- q[j + 1L,] - b
        den <- cross(v, w)
        if( abs(den) > 1e-12 ){
          t <- cross(b - a, w) / den
          u <- cross(b - a, v) / den
          if( t >= -1e-9 && t <= 1 + 1e-9 && u >= -1e-9 && u <= 1 + 1e-9 && !shared_point(s, r, a + t * v) ) {
            contact <- TRUE
            return(contact)
          }
        }else if( sum(v * v) > 1e-24 && abs(cross(b - a, v)) < 1e-10 ){
          range_q <- sort(c(sum((b - a) * v), sum((b + w - a) * v)) / sum(v * v))
          lo <- max(0, range_q[1])
          hi <- min(1, range_q[2])
          if( hi - lo > 1e-8 || (hi >= lo && !shared_point(s, r, a + lo * v)) ) {
            contact <- TRUE
            return(contact)
          }
        }
      }
    }
    contact <- FALSE

    return(contact)
  }
  contact <- function(s, r){
    if( any(s$low > r$high + 1e-9) || any(r$low > s$high + 1e-9) ) {
      contact <- FALSE
      return(contact)
    }
    if( s$lateral || r$lateral ) {
      contact <- poly_contact(s, r)
      return(contact)
    }
    lo <- max(min(s$a[1], s$b[1]), min(r$a[1], r$b[1]))
    hi <- min(max(s$a[1], s$b[1]), max(r$a[1], r$b[1]))
    if( hi < lo - 1e-10 ) {
      contact <- FALSE
      return(contact)
    }
    if( hi - lo < 1e-10 ){
      contact <- abs(y_at(s, lo) - y_at(r, lo)) < 1e-9 && !shared_point(s, r, c(lo, y_at(s, lo)))
      return(contact)
    }
    x <- c(lo, (lo + hi) / 2, hi)
    f <- y_at(s, x) - y_at(r, x)
    a <- 2 * (f[3] + f[1] - 2 * f[2])
    b <- 4 * f[2] - f[3] - 3 * f[1]
    c <- f[1]
    if( max(abs(c(a, b, c))) < 1e-10 ) {
      contact <- TRUE
      return(contact)
    }
    if( abs(a) < 1e-10 ) roots <- if( abs(b) < 1e-10 ) numeric() else -c / b else{
      discriminant <- b * b - 4 * a * c
      tolerance <- 1e-12 * (b * b + abs(4 * a * c))
      roots <- if( discriminant < -tolerance ) numeric() else (-b + c(-1, 1) * sqrt(max(0, discriminant))) / (2 * a)
    }
    roots <- roots[roots >= -1e-8 & roots <= 1 + 1e-8]
    roots <- roots[abs(a * roots^2 + b * roots + c) <= 1e-9 * max(abs(c(a, b, c)), 1e-12)]
    contact <- any(vapply(roots, function(t){
      at <- lo + t * (hi - lo)
      matches <- !shared_point(s, r, c(at, y_at(s, at)))

      return(matches)
    }, logical(1)))

    return(contact)
  }
  contact_cache <- new.env(hash = TRUE, parent = emptyenv())
  cached_contact <- function(s, r){
    key <- if( s$lateral || r$lateral ) paste(s$geometry_id, r$geometry_id, sep = ":") else
      paste(min(s$geometry_id, r$geometry_id), max(s$geometry_id, r$geometry_id), sep = ":")
    result <- contact_cache[[key]]
    if( !is.null(result) ) return(result)
    result <- contact(s, r)
    contact_cache[[key]] <- result

    return(result)
  }
  if( is.null(other_paths) ){
    direction <- as.character(df$direction)
    bidirected <- df[!is.na(direction) & direction == "<->", , drop = FALSE]
    undirected <- df[!is.na(direction) & direction == "--", , drop = FALSE]
    other_paths <- as.data.frame(data.table::rbindlist(
      list(bidirected_edge_paths(bidirected),
           straight_edge_paths(undirected, "undirected")),
      use.names = TRUE, fill = TRUE
    ))
  }
  other <- other_paths
  fixed <- lapply(split(other, other$edge_id), function(path){
    points <- sweep(sweep(as.matrix(path[, c("x", "y")]), 2, origin), 2, scale, "/")
    value <- cache_geometry(list(from = path$from[1], to = path$to[1], lateral = TRUE, points = points))

    return(value)
  })
  contacts <- function(s, index) {
    contacts <- sum(vapply(c(states[-index], fixed), function(r) {
    value <- cached_contact(s, r)

    return(value)
  }, logical(1)))

    return(contacts)
  }
  candidates <- lapply(seq_len(n), function(i){
    state <- states[[i]]
    v <- state$b - state$a
    bend_axis <- if( state$lateral ) 1L else 2L
    delta <- v * scale
    base <- if( mode == "absolute" ) height / scale[bend_axis] else
      curve_relative_base_height(delta[1L], delta[2L], scale[1L], scale[2L],
                                  height, state$lateral) / scale[bend_axis]
    if( sum(v * v) < 1e-24 || base == 0 ) {
      values <- list(state)
      return(values)
    }
    options <- list(state)
    straight <- state
    straight$h <- 0
    options[[2]] <- straight
    for(sign in c(state$sign, -state$sign)){
      s <- state
      s$sign <- sign
      s$h <- base
      options[[length(options) + 1L]] <- s
      if( mode == "relative" ){
        axis <- if( state$lateral ) 2L else 1L
        dep <- 3L - axis
        t <- (xy[, axis] - state$a[axis]) / v[axis]
        inside <- t > 1e-8 & t < 1 - 1e-8
        if( any(inside) ){
          direction <- sign * (if( state$lateral ) -1 else 1)
          # Fit above or below the node margin, allowing for the arc's slope.
          angle <- seq(0, pi, length.out = 13L)
          obstacle_axis <- as.vector(outer(xy[inside, axis], margin * cos(angle), "+"))
          obstacle_dep <- as.vector(outer(direction * xy[inside, dep], margin * sin(angle), "+"))
          sample_t <- (obstacle_axis - state$a[axis]) / v[axis]
          valid <- sample_t > 1e-8 & sample_t < 1 - 1e-8
          sample_t <- sample_t[valid]
          line <- state$a[dep] + sample_t * v[dep]
          required <- max(base, (obstacle_dep[valid] - direction * line) / (4 * sample_t * (1 - sample_t)))
          if( required > base + 1e-8 && required * 1.015 <= max(base, 1) ){
            s$h <- required * 1.015
            options[[length(options) + 1L]] <- s
          }
        }
      }
    }
    key <- vapply(options, function(s) {
      text <- paste(signif(s$h, 12), s$sign)

      return(text)
    }, character(1))
    value <- options[!duplicated(key)]

    return(value)
  })
  candidates <- lapply(candidates, function(opts) {
    values <- lapply(opts, cache_geometry)

    return(values)
  })
  states <- lapply(candidates, function(opts) {
    value <- opts[[1L]]

    return(value)
  })
  candidate_node_scores <- lapply(candidates, function(opts) {
    values <- lapply(opts, node_score)

    return(values)
  })
  order_edges <- order(-abs(edges$xend - edges$x) / initial$xref, as.character(edges$name), as.character(edges$to)) # Process wider horizontal spans first
  changes <- 0L
  passes <- 0L
  for(pass in 1:2){ # Try up to two passes to improve the routes
    changed <- FALSE
    passes <- pass
    for(i in order_edges){
      current <- states[[i]]
      own_score <- node_score(current)
      current_score <- c(own_score, contacts(current, i))
      if( current_score[1] == 0 && current_score[3] == 0 ) next
      node_scores <- t(vapply(candidate_node_scores[[i]], identity, numeric(2)))
      node_scores[, 2] <- round(node_scores[, 2], 8)
      best <- do.call(order, as.data.frame(node_scores))[1]
      eligible <- which(node_scores[, 1] == node_scores[best, 1] & node_scores[, 2] == node_scores[best, 2])
      scores <- t(vapply(eligible, function(j){
        values <- c(candidate_node_scores[[i]][[j]], contacts(candidates[[i]][[j]], i), candidates[[i]][[j]]$h)

        return(values)
      }, numeric(4)))
      scores[, 2] <- round(scores[, 2], 8)
      ranked <- do.call(order, as.data.frame(scores))[1]
      choice <- eligible[ranked]
      proposal_score <- scores[ranked, 1:3]
      current_score[2] <- round(current_score[2], 8)
      delta <- which(abs(proposal_score - current_score) > 1e-8)
      if( length(delta) && proposal_score[delta[1]] < current_score[delta[1]] ){
        states[[i]] <- candidates[[i]][[choice]]
        changes <- changes + 1L
        changed <- TRUE
      }
    }
    if( !changed ) break
  }
  default_signs <- vapply(states, function(s) {
    value <- s$sign

    return(value)
  }, numeric(1))
  if( position == "flipped" ) for(i in seq_len(n)) states[[i]]$sign <- -states[[i]]$sign
  if( position %in% c("over", "under") ) for(i in seq_len(n)) states[[i]]$sign <- if( position == "over" ) 1 else -1
  paths <- lapply(seq_len(n), function(i){
    s <- states[[i]]
    # Intermediate straight-path samples retain short capped arrows.
    t <- seq(0, 1, length.out = if( s$h > 0 ) 201L else 17L)
    points <- sweep(sweep(coordinates(s, t), 2, scale, "*"), 2, origin, "+")
    p <- edges[rep(i, length(t)), , drop = FALSE]
    p$x <- points[, 1]
    p$y <- points[, 2]
    p$index <- t
    p$.edge_group <- i
    p$edge_id <- paste("directed", s$from, s$to, sep = "\036")
    p$from <- s$from
    p$to <- s$to
    p$kind <- "directed"
    p$curved <- s$h > 0

    return(p)
  })
  controls <- lapply(seq_len(n), function(i){
    s <- states[[i]]
    axis <- if( s$lateral ) 1L else 2L
    delta <- (s$b - s$a) * scale
    base <- if( mode == "absolute" ) height else
      curve_relative_base_height(delta[1L], delta[2L], scale[1L], scale[2L],
                                  height, s$lateral)
    cp <- (s$a + s$b) / 2
    cp[axis] <- cp[axis] + (if( s$lateral ) -1 else 1) * 2 * s$sign * s$h
    cp <- cp * scale + origin
    values <- list(edge_id = paste(s$from, s$to, sep = "->"), from = s$from, to = s$to, curved = s$h > 0,
                   height = s$h * scale[axis], base_height = base, height_multiplier = if( base > 0 ) s$h * scale[axis] / base else 0,
                   control_x = cp[1], control_y = cp[2], sign = s$sign, default_sign = default_signs[i], lateral = s$lateral)

    return(values)
  })
  initial$paths <- do.call(rbind, paths)
  initial$controls <- controls
  initial$specifications <- NULL
  initial$routing_metrics <- list(changes = changes, passes = passes, normalised_margin = margin,
    candidate_count = sum(lengths(candidates)), label_inputs = FALSE)

  return(initial)
}


#' Find the typical gap along a coordinate axis
#'
#' @param x Numeric coordinates along one axis.
#' @return The median positive gap, or `1` when no gap is available.
#' @noRd
curve_axis_gap <- function(x){
  gaps <- diff(sort(unique(x[is.finite(x)])))
  if( length(gaps) ) {
    gap <- stats::median(gaps)

    return(gap)
  } else {
    gap <- 1

    return(gap)
  }
}


#' Measure the distance from a point to a sampled path
#'
#' @param point Numeric point coordinates in x, y order.
#' @param path Data frame of ordered path points with `x` and `y` columns.
#' @return The shortest Euclidean distance to any path segment.
#' @noRd
curve_point_path_distance <- function(point, path){
  a <- as.matrix(path[-nrow(path), c("x", "y"), drop = FALSE])
  b <- as.matrix(path[-1, c("x", "y"), drop = FALSE])
  v <- b - a
  w <- matrix(point, nrow(a), 2, byrow = TRUE) - a
  vv <- rowSums(v * v)
  t <- pmax(0, pmin(1, rowSums(w * v) / ifelse(vv > 0, vv, 1)))
  distance <- min(sqrt(rowSums((w - v * t)^2)))

  return(distance)
}


#' Construct initial directed curves around intervening nodes
#'
#' @param df Node and edge table with endpoint coordinates.
#' @param height Non-negative starting curve height.
#' @param mode Curve height mode: `"relative"` or `"absolute"`.
#' @param position Curve position: `"default"`, `"flipped"`, `"over"`, or
#'   `"under"`.
#' @return A curve plan containing paths, control points, node coordinates,
#'   reference spans, and endpoint errors.
#' @noRd
initial_curve_plan <- function(df, height, mode, position){
  nodes <- df[!duplicated(df$name), c("name", "x", "y"), drop = FALSE]
  nodes$name <- as.character(nodes$name)
  X <- diff(range(nodes$x))
  Y <- diff(range(nodes$y))
  xref <- if( X > 0 ) X else 2 * curve_axis_gap(nodes$y)
  yref <- if( Y > 0 ) Y else curve_axis_gap(nodes$x) / 2
  edges <- df[!is.na(df$direction) & as.character(df$direction) == "->", , drop = FALSE]
  specifications <- vector("list", nrow(edges))
  errors <- character()
  for(i in seq_len(nrow(edges))){
    e <- edges[i, , drop = FALSE]
    a <- c(e$x, e$y)
    b <- c(e$xend, e$yend)
    v <- b - a
    key <- paste(as.character(e$name), as.character(e$to), sep = "->")
    lateral <- curve_is_lateral(e$x, e$xend)
    normalised <- v / c(xref, yref)
    others <- nodes[!nodes$name %in% c(as.character(e$name), as.character(e$to)),]
    offsets <- sweep(as.matrix(others[, c("x", "y"), drop = FALSE]), 2, a)
    offsets <- sweep(offsets, 2, c(xref, yref), "/")
    length_squared <- sum(normalised^2)
    if( length_squared <= 1e-24 ){
      obstruction <- FALSE
      errors <- c(errors, paste("Coincident directed endpoints:", key))
    }else{
      projections <- as.numeric(offsets %*% normalised) / length_squared
      perpendicular <- (normalised[1] * offsets[, 2] - normalised[2] * offsets[, 1]) / sqrt(length_squared)
      obstruction <- any(projections > 1e-10 & projections < 1 - 1e-10 & abs(perpendicular) < 1e-9)
    }
    h <- if( !obstruction ) 0 else if( mode == "absolute" ) height else
      curve_relative_base_height(v[1L], v[2L], xref, yref, height, lateral)
    specifications[[i]] <- list(index = i, key = key, from = as.character(e$name), to = as.character(e$to),
      a = a, b = b, v = v, lateral = lateral, obstruction = obstruction, h = h, length = sqrt(length_squared))
  }
  path_for <- function(spec, sign, n = 201L){
    t <- if( spec$h > 0 ) seq(0, 1, length.out = n) else c(0, 1)
    path <- edges[rep(spec$index, length(t)), , drop = FALSE]
    path$x <- spec$a[1] + t * spec$v[1]
    path$y <- spec$a[2] + t * spec$v[2]
    if( spec$lateral ) path$x <- path$x - sign * 4 * spec$h * t * (1 - t)
    else path$y <- path$y + sign * 4 * spec$h * t * (1 - t)
    path$index <- t
    path$.edge_group <- spec$index
    path$edge_id <- paste("directed", spec$from, spec$to, sep = "\036")
    path$from <- spec$from
    path$to <- spec$to
    path$kind <- "directed"
    path$curved <- spec$h > 0

    return(path)
  }
  signs <- rep(1L, length(specifications))
  tied <- rep(FALSE, length(signs))
  for(i in seq_along(specifications)){
    spec <- specifications[[i]]
    if( !spec$obstruction ) next
    others <- nodes[!nodes$name %in% c(spec$from, spec$to),]
    others$x <- others$x / xref
    others$y <- others$y / yref
    scores <- vapply(c(1L, -1L), function(sign){
      path <- path_for(spec, sign)
      path$x <- path$x / xref
      path$y <- path$y / yref
      value <- min(vapply(seq_len(nrow(others)), function(j) {
        value <- curve_point_path_distance(c(others$x[j], others$y[j]), path)

        return(value)
      }, numeric(1)))

      return(value)
    }, numeric(1))
    tied[i] <- abs(diff(scores)) <= 1e-9
    signs[i] <- if( scores[1] >= scores[2] ) 1L else -1L
  }
  # Retain deterministic alternation for tied, collinear overlapping routes.
  free <- which(tied)
  adjacency <- matrix(FALSE, length(signs), length(signs))
  if( length(free) > 1 ) for(aa in seq_len(length(free) - 1L)) for(bb in seq.int(aa + 1L, length(free))){
    i <- free[aa]
    j <- free[bb]
    first <- specifications[[i]]
    second <- specifications[[j]]
    a <- first$a / c(xref, yref)
    u <- first$v / c(xref, yref)
    u <- u / sqrt(sum(u * u))
    c0 <- second$a / c(xref, yref)
    d0 <- second$b / c(xref, yref)
    cross <- function(z) {
      cross_product <- u[1] * z[2] - u[2] * z[1]

      return(cross_product)
    }
    if( abs(cross(c0 - a)) < 1e-9 && abs(cross(d0 - a)) < 1e-9 ){
      interval <- sort(c(sum((c0 - a) * u), sum((d0 - a) * u)))
      if( min(first$length, interval[2]) - max(0, interval[1]) > 1e-9 ) adjacency[i, j] <- adjacency[j, i] <- TRUE
    }
  }
  unassigned <- free
  while(length(unassigned)){
    component <- unassigned[1]
    frontier <- component
    while(length(frontier)){
      neighbours <- setdiff(which(colSums(adjacency[frontier, , drop = FALSE]) > 0), component)
      component <- c(component, neighbours)
      frontier <- neighbours
    }
    unassigned <- setdiff(unassigned, component)
    if( length(component) > 1 ){
      ranking <- order(-vapply(specifications[component], function(s) {
        value <- s$length

        return(value)
      }, numeric(1)),
                       vapply(specifications[component], function(s) {
                         value <- s$key

                         return(value)
                       }, character(1)))
      signs[component[ranking]] <- rep(c(-1L, 1L), length.out = length(component))
    }
  }
  default_signs <- signs
  if( position == "flipped" ) signs <- -signs
  if( position == "over" ) signs[] <- 1L
  if( position == "under" ) signs[] <- -1L
  paths <- if( length(specifications) ) do.call(rbind, lapply(seq_along(specifications), function(i) {
    value <- path_for(specifications[[i]], signs[i])

    return(value)
  })) else NULL
  controls <- lapply(seq_along(specifications), function(i){
    spec <- specifications[[i]]
    control <- (spec$a + spec$b) / 2
    if( spec$lateral ) control[1] <- control[1] - 2 * signs[i] * spec$h
    else control[2] <- control[2] + 2 * signs[i] * spec$h
    values <- list(edge_id = spec$key, from = spec$from, to = spec$to, curved = spec$h > 0, height = spec$h,
                   control_x = control[1], control_y = control[2], sign = signs[i], default_sign = default_signs[i], lateral = spec$lateral)

    return(values)
  })
  anchors <- if( X == 0 || Y == 0 ) data.frame(
    x = if( X == 0 ) mean(range(nodes$x)) + c(-.5, .5) * xref else range(nodes$x),
    y = if( Y == 0 ) mean(range(nodes$y)) + c(-.5, .5) * yref else range(nodes$y)) else NULL
  plan <- list(paths = paths, controls = controls, nodes = nodes, anchors = anchors, errors = errors,
               xref = xref, yref = yref, X = X, Y = Y, specifications = specifications)

  return(plan)
}
