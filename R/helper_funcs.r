# ==============================================================================
# helper_funcs.r


# Shared graph parsing, editing and inspection helpers
# Shared helpers appear once under their main caller's file
# ==============================================================================


################################################################################
# coords_funcs.r
################################################################################
#' Helper for ranking nodes for layout and role ordering
#'
#' Roots receive rank 0. Every other node receives one rank more than its
#' highest ranked parent.
#'
#' @param parents Named list with one entry per node. Each entry is a character
#'   vector of that node's parent names, or character(0) for a root.
#' @return Named integer vector of ranks, one per node, in names(parents) order.
#' @noRd
topological_rank <- function(parents){

  if( !is.list( parents ) ){
    stop( "The 'parents' argument must be a named list mapping each node to a character vector of its parent names." )
  }

  if( length( parents ) == 0L ){
    rank_out <- integer( 0L )
    return(rank_out)
  }

  if( is.null( names( parents ) ) ){
    stop( "The 'parents' argument must be a named list mapping each node to a character vector of its parent names." )
  }

  nodes <- names( parents )
  n <- length( nodes )

  if( anyDuplicated( nodes ) > 0L ){
    stop( "The names of 'parents' must be unique; each node may appear only once." )
  }

  all_parents <- unique( unlist( parents, use.names = FALSE ) )
  unknown <- all_parents[ !all_parents %in% nodes ]
  if( length( unknown ) > 0L ){
    stop( "Every parent must also be a node in 'parents'; the following parent names are not present: ",
          paste( unknown, collapse = ", " ), "." )
  }

  ## build child adjacency and in-degree from the parent list ##
  children <- vector( "list", n )
  names( children ) <- nodes
  in_degree <- integer( n )
  names( in_degree ) <- nodes

  for( v in nodes ){
    node_parents <- parents[[v]]
    in_degree[[v]] <- length( node_parents )
    for( p in node_parents ){
      children[[p]] <- c( children[[p]], v )
    }
  }

  ## Kahn's algorithm: process parents before assigning each node's final rank ##
  rank <- integer( n )
  names( rank ) <- nodes
  remaining <- in_degree
  queue <- nodes[ in_degree == 0L ]
  processed <- 0L

  while( length( queue ) > 0L ){

    v <- queue[[1L]]
    queue <- queue[ -1L ]
    processed <- processed + 1L

    for( w in children[[v]] ){
      if( rank[[v]] + 1L > rank[[w]] ){
        rank[[w]] <- rank[[v]] + 1L
      }
      remaining[[w]] <- remaining[[w]] - 1L
      if( remaining[[w]] == 0L ){
        queue <- c( queue, w )
      }
    }
  }

  if( processed < n ){
    stop( "The graph contains a directed cycle. Coordinates can only be assigned to a directed acyclic graph." )
  }

  return( rank )
}


#' Helper for finding edges that close directed cycles
#'
#' Uses depth-first search to identify back edges, including self-loops.
#' Removing these edges makes the directed graph acyclic. Visits nodes in
#' `node_names` order and outgoing edges in their supplied order.
#'
#' @param node_names Character vector of all node names.
#' @param from,to Equal-length vectors of directed-edge source and target names.
#' @return Integer positions of the back edges in `from` and `to`.
#' @noRd
feedback_arcs <- function( node_names, from, to ){
  n <- length( node_names )
  if( n == 0L || length( from ) == 0L ){
    arcs <- integer(0)
    return(arcs)
  }
  idx_by_node <- split( seq_along(from), factor(from, levels = node_names) )
  color <- stats::setNames( rep(0L, n), node_names )   # 0 unseen, 1 on stack, 2 done
  back  <- integer(0)
  for( s in node_names ){
    if( color[[s]] != 0L ){ next }
    color[[s]] <- 1L
    node_stack <- s
    ptr_stack  <- 1L
    while( length( node_stack ) > 0L ){
      u   <- node_stack[[ length(node_stack) ]]
      ei  <- ptr_stack[[ length(ptr_stack) ]]
      out <- idx_by_node[[ u ]]
      if( ei <= length( out ) ){
        ptr_stack[[ length(ptr_stack) ]] <- ei + 1L
        e  <- out[[ ei ]]
        v  <- to[[ e ]]
        cv <- color[[ v ]]
        if( cv == 0L ){
          color[[ v ]] <- 1L
          node_stack <- c( node_stack, v )
          ptr_stack  <- c( ptr_stack, 1L )
        } else if( cv == 1L ){
          back <- c( back, e )   # v is on the current stack, so (u -> v) closes a cycle
        }
      } else {
        color[[ u ]] <- 2L
        node_stack <- node_stack[ -length(node_stack) ]
        ptr_stack  <- ptr_stack[ -length(ptr_stack) ]
      }
    }
  }
  return( back )
}


#' Shared helper for updating node and curve coordinates
#'
#' @param dag A dagitty object.
#' @param coordinates List of named x and y vectors in dagitty coordinates.
#' @param caller Function name to include in errors.
#' @param controls Stored edge controls to carry through the coordinate change.
#' @return A dagitty object with updated coordinates and transformed edge
#'   controls.
#' @noRd
set_coords_preserving_curves <- function(dag, coordinates, caller = "add_coords",
                                          controls = edge_controls(dag, caller)){
  before <- dagitty::coordinates(dag)
  result <- write_node_coordinates(dag, coordinates, caller)
  if( nrow(controls) > 0L ){
    controls <- transfer_edge_controls(controls, before, dagitty::coordinates(result), caller)
    result <- set_edge_controls(result, controls, caller = caller)
  }

  return(result)
}


#' Write node coordinates in native graph text
#'
#' @param dag A dagitty object.
#' @param coordinates List of named x and y vectors in dagitty coordinates.
#' @param caller Function name to include in errors.
#' @return A dagitty object with coordinates updated for the matching nodes.
#' @noRd
write_node_coordinates <- function(dag, coordinates, caller){
  nodes <- intersect(names(coordinates$x), names(coordinates$y))
  nodes <- intersect(nodes, names(dag))
  if( !length(nodes) ) {
    return(dag)
  }
  if( any(!is.finite(coordinates$x[nodes]) | !is.finite(coordinates$y[nodes])) ){
    curve_layout_error(caller)
  }
  position <- stats::setNames(paste0('pos="', native_layout_number(coordinates$x[nodes]),
                                     ",", native_layout_number(coordinates$y[nodes]), '"'), nodes)
  text <- as.character(dag)
  index <- native_layout_index(text, caller)
  edits <- list()
  positioned <- character()
  for( block in index$blocks ){
    if( block$kind != "node" || !block$node %in% nodes ) {
      next
    }
    edits <- c(edits, native_pos_edits(index, block))
    retained <- any(vapply(block$fields, function(field) {
      matches <- field$key != "pos"

      return(matches)
    }, logical(1L)))
    edits[[length(edits) + 1L]] <- list(start = index$tokens$start[block$last],
                                       end = index$tokens$start[block$last] - 1L,
                                       text = paste0(if( retained ) { "," } else { "" }, position[block$node]))
    positioned <- union(positioned, block$node)
  }
  for( declaration in index$declarations ){
    node <- declaration$node
    if( !node %in% setdiff(nodes, positioned) ) {
      next
    }
    last <- index$tokens$end[declaration$token]
    edits[[length(edits) + 1L]] <- list(start = last + 1L, end = last,
                                       text = paste0(" [", position[node], "]"))
    positioned <- union(positioned, node)
  }
  missing <- setdiff(nodes, positioned)
  if( length(missing) ){
    lines <- paste0(native_layout_quote(missing), " [", position[missing], "]")
    edits[[length(edits) + 1L]] <- list(start = index$closing, end = index$closing - 1L,
                                       text = paste0("\n", paste(lines, collapse = "\n"), "\n"))
  }
  result <- native_layout_value(dag, native_layout_edits(text, edits))
  updated <- tryCatch(dagitty::coordinates(result),
                     error = function(e) curve_layout_error(caller))
  if( !identical(as.numeric(updated$x[nodes]), as.numeric(coordinates$x[nodes])) ||
      !identical(as.numeric(updated$y[nodes]), as.numeric(coordinates$y[nodes])) ){
    curve_layout_error(caller)
  }

  return(result)
}


#' Format native coordinate values
#'
#' @param value Numeric coordinate values.
#' @return Character vector of decimal values without exponent notation.
#' @noRd
native_layout_number <- function(value){
  value <- as.numeric(value)
  text <- sprintf("%.17g", value)
  for( i in which(grepl("e", text, fixed = TRUE)) ){
    text[i] <- native_layout_decimal(value[i], as.integer(sub(".*e", "", text[i])))
  }

  return(text)
}


#' Write a number in decimal notation
#'
#' @param value A numeric value.
#' @param exponent The base-10 exponent from the formatted value.
#' @return A character string containing the decimal value.
#' @noRd
native_layout_decimal <- function(value, exponent){
  text <- sprintf("%.*f", max(0L, 16L - exponent), value)

  return(text)
}


#' Shared helper for moving control points with their endpoints
#'
#' @param controls Data frame of edge controls with columns v, e, w, x and y.
#' @param before List of named x and y endpoint coordinates before the change.
#' @param after List of named x and y endpoint coordinates after the change.
#' @param caller Function name to include in errors.
#' @return The edge control table with positions adjusted to the new
#'   coordinates.
#' @noRd
transfer_edge_controls <- function(controls, before, after, caller){
  for( i in seq_len(nrow(controls)) ){
    endpoints <- c(controls$v[i], controls$w[i])
    old <- cbind(before$x[endpoints], before$y[endpoints])
    new <- cbind(after$x[endpoints], after$y[endpoints])
    if( identical(as.numeric(old), as.numeric(new)) ) {
      next
    }
    if( any(!is.finite(old)) || any(!is.finite(new)) ) {
      curve_layout_error(caller, TRUE)
    }
    old_chord <- old[2L, ] - old[1L, ]
    new_chord <- new[2L, ] - new[1L, ]
    old_length <- sqrt(sum(old_chord^2))
    new_length <- sqrt(sum(new_chord^2))
    if( !is.finite(old_length) || !is.finite(new_length) ||
        old_length == 0 || new_length == 0 ) {
      curve_layout_error(caller, TRUE)
    }
    unit <- old_chord / old_length
    displacement <- c(controls$x[i], controls$y[i]) - old[1L, ]
    along <- sum(displacement * unit) / old_length
    normal <- sum(displacement * c(-unit[2L], unit[1L])) / old_length
    control <- new[1L, ] + along * new_chord + normal * c(-new_chord[2L], new_chord[1L])
    if( any(!is.finite(control)) ) {
      curve_layout_error(caller, TRUE)
    }
    controls$x[i] <- control[1L]
    controls$y[i] <- control[2L]
  }

  return(controls)
}


################################################################################
# copy_nodes.R
################################################################################
#' Helper for copy_nodes() and merge_nodes()
#'
#' @param dag A dagitty object.
#' @param rename_map Named character vector mapping old node names to new
#'   names.
#' @param caller Name of the calling function for messages.
#' @return A dagitty object with renamed nodes and label names.
#' @noRd
copy_rename_native_nodes <- function(dag, rename_map, caller){
  unchanged <- names(rename_map) == unname(rename_map)
  rename_map <- rename_map[!unchanged]
  if( length(rename_map) == 0L ) {
    return(dag)
  }

  index <- copy_native_node_index(dag, caller)
  raw_values <- index$tokens$value[index$node_tokens]
  decoded <- vapply(raw_values, function(value){
    if( startsWith(value, '"') ){
      value <- substr(value, 2L, nchar(value) - 1L)
      value <- gsub('\\\\(.)', '\\1', value, perl = TRUE)
    }

    return(value)
  }, character(1L), USE.NAMES = FALSE)
  selected <- which(decoded %in% names(rename_map))
  edits <- lapply(selected, function(i){
    token <- index$node_tokens[i]
    edit <- list(
      start = index$tokens$start[token],
      end = index$tokens$end[token],
      text = native_layout_quote(unname(rename_map[decoded[i]]))
    )

    return(edit)
  })
  result <- native_layout_value(
    dag, native_layout_edits(index$text, edits)
  )

  labels <- attr(result, "labels", exact = TRUE)
  if( !is.null(labels) && !is.null(names(labels)) ){
    names(labels) <- copy_remap_names(names(labels), rename_map)
    attr(result, "labels") <- labels
  }

  expected <- c(setdiff(names(dag), names(rename_map)), unname(rename_map))
  if( !setequal(names(result), expected) ){
    stop(caller, ": existing node declarations could not be renamed safely.",
         call. = FALSE)
  }

  return(result)
}


#' Helper for copy_rename_native_nodes()
#'
#' @param dag A dagitty object.
#' @param caller Name of the calling function for messages.
#' @return A list containing the graph `text`, its `tokens`, and the
#'   `node_tokens` indices.
#' @noRd
copy_native_node_index <- function(dag, caller){
  text <- as.character(dag)
  tokens <- native_layout_tokens(text, caller)
  values <- tokens$value
  position <- which(values == "{")[1L]
  if( is.na(position) ) {
    curve_layout_error(caller)
  }
  node_tokens <- integer()

  decode <- function(value){
    if( startsWith(value, '"') ){
      value <- substr(value, 2L, nchar(value) - 1L)
      value <- gsub('\\\\(.)', '\\1', value, perl = TRUE)
    }

    return(value)
  }
  record <- function(endpoint){
    if( !is.null(endpoint) ) {
      recorded_tokens <- (node_tokens <<- c(node_tokens, endpoint$token))

      return(invisible(recorded_tokens))
    }else{
      recorded_tokens <- NULL
      return(invisible(recorded_tokens))
    }
  }
  skip_attributes <- function(){
    position <<- position + 1L
    while( position <= length(values) && values[position] != "]" ){
      position <<- position + 1L
    }
    if( position > length(values) ) {
      curve_layout_error(caller)
    }
    next_position <- (position <<- position + 1L)

    return(invisible(next_position))
  }
  endpoint <- function(record_now = TRUE){
    if( values[position] == "{" ){
      body()
      out <- NULL
      return(out)
    }
    if( values[position] %in% c("}", "[", "]", "=", ",", ";") ){
      curve_layout_error(caller)
    }
    result <- list(node = decode(values[position]), token = position)
    position <<- position + 1L
    if( record_now ) {
      record(result)
    }

    return(result)
  }
  body <- function(){
    position <<- position + 1L
    while( position <= length(values) && values[position] != "}" ){
      if( values[position] %in% c(",", ";") ){
        position <<- position + 1L
        next
      }
      first <- endpoint(record_now = FALSE)
      if( position <= length(values) && values[position] == "=" ){
        position <<- position + 2L
        next
      }
      record(first)
      while( position <= length(values) &&
             values[position] %in% c("->", "<-", "<->", "--") ){
        position <<- position + 1L
        endpoint(record_now = TRUE)
      }
      while( position <= length(values) && values[position] == "[" ){
        skip_attributes()
      }
    }
    if( position > length(values) ) {
      curve_layout_error(caller)
    }
    next_position <- (position <<- position + 1L)

    return(invisible(next_position))
  }
  body()
  node_index <- list(text = text, tokens = tokens, node_tokens = unique(node_tokens))

  return(node_index)
}


#' Helper for replacing names during node copying and renaming
#'
#' @param values Node names to rename.
#' @param rename_map Named character vector mapping old node names to new
#'   names.
#' @return A character vector with matched names replaced.
#' @noRd
copy_remap_names <- function(values, rename_map){
  values <- as.character(values)
  matched <- match(values, names(rename_map))
  replace <- !is.na(matched)
  values[replace] <- unname(rename_map[matched[replace]])

  return(values)
}


#' Shared helper for checking whether a new node shares coordinates
#'
#' @param coordinates A list containing named `x` and `y` coordinate vectors.
#' @param new_nodes Names of nodes to add.
#' @return `TRUE` if a new node coincides with another node, otherwise `FALSE`.
#' @noRd
copy_coincident_copies <- function(coordinates, new_nodes){
  x <- coordinates$x
  y <- coordinates$y
  shared <- intersect(names(x), names(y))
  shared <- shared[is.finite(x[shared]) & is.finite(y[shared])]
  if( length(shared) < 2L ) {
    coincident <- FALSE
    return(coincident)
  }
  places <- paste(x[shared], y[shared])
  repeated <- unique(places[duplicated(places)])
  if( length(repeated) == 0L ) {
    coincident <- FALSE
    return(coincident)
  }
  coincident <- any(shared[places %in% repeated] %in% new_nodes)

  return(coincident)
}


################################################################################
# curve_funcs.R
################################################################################
#' Read stored edge control points
#'
#' @param dag A dagitty object.
#' @param caller Function name to include in errors.
#' @return Data frame of stored controls with columns v, e, w, x and y.
#' @noRd
edge_controls <- function(dag, caller = "add_curves"){
  edges <- tryCatch(dagitty::edges(dag),
                    error = function(e) curve_layout_error(caller))
  if( nrow(edges) == 0L ){
    controls <- data.frame(v = character(), e = character(), w = character(),
                           x = numeric(), y = numeric())
    return(controls)
  }
  edges <- edges[, c("v", "e", "w", "x", "y"), drop = FALSE]
  stored <- !is.na(edges$x) | !is.na(edges$y)
  if( any(stored & (!is.finite(edges$x) | !is.finite(edges$y))) ){
    curve_layout_error(caller)
  }
  controls <- edges[stored, , drop = FALSE]

  return(controls)
}


#' Helper for reporting graph editing and curve transfer errors
#'
#' @param caller Function name to include in the error.
#' @param transfer Whether the error concerns moving an existing curve.
#' @return Stops with an error.
#' @noRd
curve_layout_error <- function(caller, transfer = FALSE){
  caller <- sub("\\(\\)$", "", caller)
  suffix <- if( transfer ){
    "cannot preserve a curved edge's shape without finite, distinct endpoint coordinates before and after the layout change."
  }else{
    "cannot read or update native edge control points without changing the graph."
  }
  stop(caller, "(): ", suffix, call. = FALSE)
}


#' Replace stored edge control points
#'
#' @param dag A dagitty object.
#' @param controls Data frame with columns v, e, w, x and y. Two missing
#'   coordinates clear a control point.
#' @param replace Edge table identifying the control points to replace or clear.
#' @param caller Function name to include in errors.
#' @return A dagitty object with the specified edge controls replaced.
#' @noRd
set_edge_controls <- function(dag, controls,
                               replace = controls[, c("v", "e", "w"), drop = FALSE],
                               caller = "add_curves"){
  previous <- edge_controls(dag, caller)
  keys <- function(x) {
    keys <- edge_control_keys(x$v, x$e, x$w)

    return(keys)
  }
  before <- dagitty::edges(dag)
  if( any(!keys(replace) %in% keys(before)) || anyDuplicated(keys(controls)) ){
    curve_layout_error(caller)
  }
  clear <- is.na(controls$x) & is.na(controls$y)
  if( any(!clear & (!is.finite(controls$x) | !is.finite(controls$y))) ||
      any(!keys(controls) %in% keys(replace)) ) {
    curve_layout_error(caller)
  }
  controls <- controls[!clear, , drop = FALSE]
  controls <- rbind(previous[!keys(previous) %in% keys(replace), , drop = FALSE],
                    controls[, c("v", "e", "w", "x", "y"), drop = FALSE])
  text <- as.character(dag)
  index <- native_layout_index(text, caller)
  edits <- unlist(lapply(index$blocks, function(block){
    if( block$kind == "edge" ) {
      position_edits <- native_pos_edits(index, block)

      return(position_edits)
    } else {
      position_edits <- list()

      return(position_edits)
    }
  }), recursive = FALSE)
  lines <- character()
  for( i in seq_len(nrow(controls)) ){
    control <- controls[i, , drop = FALSE]
    pos <- paste0('pos="', native_layout_number(control$x), ",",
                  native_layout_number(control$y), '"')
    matches <- which(vapply(index$blocks, function(block){
      matches <- block$kind == "edge" && length(block$keys) == 1L && block$keys == keys(control)

      return(matches)
    }, logical(1L)))
    if( length(matches) ){
      block <- index$blocks[[tail(matches, 1L)]]
      retained <- any(vapply(block$fields, function(field) {
        matches <- field$key != "pos"

        return(matches)
      }, logical(1L)))
      edits[[length(edits) + 1L]] <- list(start = index$tokens$start[block$last],
                                         end = index$tokens$start[block$last] - 1L,
                                         text = paste0(if( retained ) { "," } else { "" }, pos))
    }else{
      lines <- c(lines, paste0(native_layout_quote(control$v), " ", control$e, " ",
                                native_layout_quote(control$w), " [", pos, "]"))
    }
  }
  if( length(lines) ){
    edits[[length(edits) + 1L]] <- list(start = index$closing,
                                       end = index$closing - 1L,
                                       text = paste0("\n", paste(lines, collapse = "\n"), "\n"))
  }
  result <- native_layout_value(dag, native_layout_edits(text, edits))
  after <- tryCatch(dagitty::edges(result),
                   error = function(e) curve_layout_error(caller))
  actual <- edge_controls(result, caller)
  matched <- match(keys(controls), keys(actual))
  if( !setequal(keys(before), keys(after)) || nrow(actual) != nrow(controls) ||
      anyNA(matched) ||
      !identical(as.numeric(actual$x[matched]), as.numeric(controls$x)) ||
      !identical(as.numeric(actual$y[matched]), as.numeric(controls$y)) ||
      !identical(dagitty::coordinates(dag), dagitty::coordinates(result)) ){
    curve_layout_error(caller)
  }

  return(result)
}


#' Create edge keys from endpoints and operators
#'
#' @param v,w Node names at each end of the edges.
#' @param e Edge operators.
#' @return Character vector of edge keys; reversing a symmetric edge gives the
#'   same key.
#' @noRd
edge_control_keys <- function(v, e, w){
  if( !length(v) ) {
    keys <- character()
    return(keys)
  }
  v <- as.character(v)
  e <- as.character(e)
  w <- as.character(w)
  reverse <- e == "<-"
  first <- v[reverse]
  v[reverse] <- w[reverse]
  w[reverse] <- first
  e[reverse] <- "->"
  swap <- e %in% c("<->", "--") & v > w
  first <- v[swap]
  v[swap] <- w[swap]
  w[swap] <- first
  keys <- paste0(nchar(v, type = "bytes"), ":", v,
                 nchar(e, type = "bytes"), ":", e,
                 nchar(w, type = "bytes"), ":", w)

  return(keys)
}


#' Remove position fields from an attribute block
#'
#' @param index The token and declaration index for the graph.
#' @param block A node or edge attribute block from the index.
#' @return List of text edits removing pos fields and adjacent commas.
#' @noRd
native_pos_edits <- function(index, block){
  values <- index$tokens$value
  removed <- integer()
  for( field in block$fields ){
    if( field$key != "pos" ) {
      next
    }
    removed <- c(removed, seq.int(field$first, field$last))
    if( values[field$last + 1L] == "," ){
      removed <- c(removed, field$last + 1L)
    }else if( values[field$first - 1L] == "," ){
      removed <- c(removed, field$first - 1L)
    }
  }
  edits <- lapply(unique(removed), function(i){
    edit <- list(start = index$tokens$start[i], end = index$tokens$end[i], text = "")

    return(edit)
  })

  return(edits)
}


################################################################################
# draw_funcs.r
################################################################################
#' Collect children in node and edge order
#'
#' @param edges An edge table.
#' @param nodes Parent node names in the required order.
#' @return Character vector of unique children, ordered by the supplied nodes
#'   and edges.
#' @noRd
children_edge_order <- function(edges, nodes){
  edges <- normalise_edge_frame(edges)
  directed <- edges[edges$e == "->", , drop = FALSE]
  children <- Reduce(function(found, node){
    children <- union(found, directed$w[directed$v == node])

    return(children)
  }, nodes, character())

  return(children)
}


#' Collect parents in node and edge order
#'
#' @param edges An edge table.
#' @param nodes Child node names in the required order.
#' @return Character vector of unique parents, ordered by the supplied nodes and
#'   edges.
#' @noRd
parents_edge_order <- function(edges, nodes){
  edges <- normalise_edge_frame(edges)
  directed <- edges[edges$e == "->", , drop = FALSE]
  parents <- Reduce(function(found, node){
    parents <- union(found, directed$v[directed$w == node])

    return(parents)
  }, nodes, character())

  return(parents)
}


################################################################################
# get_funcs.r
################################################################################
#' Create a graph read cache
#'
#' Use one cache for one graph during a public operation.
#'
#' @return An empty environment.
#' @noRd
graph_cache <- function() { # Create a cache for one operation on one unchanged graph
  cache <- new.env(parent = emptyenv())

  return(cache)
}


#' Read the cached edge table
#'
#' Computes the table on the first request. Copy it before modifying it
#' in place.
#'
#' @importFrom dagitty edges
#' @param dag A dagitty object.
#' @param .cache An environment from graph_cache(), or NULL to read directly.
#' @return A data frame of graph edges with character columns v, e and w.
#' @noRd
cached_edges <- function(dag, .cache = NULL) {
  edges <- cached_graph_value(.cache, "edges", function() {
    edges <- dag_edges(dag)

    return(edges)
  })

  return(edges)
}


#' Read or compute a cached graph value
#'
#' A NULL cache calls `compute()` directly.
#'
#' @param .cache An environment from graph_cache(), or NULL to bypass the cache.
#' @param key Character key for the cached value.
#' @param compute Function with no arguments that computes a missing value.
#' @return The cached or freshly computed value.
#' @noRd
cached_graph_value <- function(.cache, key, compute) {
  if( is.null(.cache) ){
    value_visibility <- withVisible(compute())
    value <- value_visibility$value
    if( !value_visibility$visible ){
      return(invisible(value))
    }
    return(value)
  }
  if( exists(key, envir = .cache, inherits = FALSE) ){
    value <- get(key, envir = .cache, inherits = FALSE)
    return(value)
  }
  val <- compute()
  assign(key, val, envir = .cache)

  return(val)
}


#' Read all node names in a dagitty object
#'
#' Includes isolated nodes.
#'
#' @param dag A dagitty object.
#' @returns Character vector of node names, or character(0) if none.
#' @noRd
get_node_names <- function(dag){

  s  <- dag_to_string(dag)
  # Remove the graph type and optional graph name.
  s  <- sub("^\\s*(dag|pdag|mag|pag)\\b[^{]*", " ", s, perl = TRUE)
  s  <- gsub('("(?:\\\\"|[^"])*")|\\[[^][]*\\]', "\\1 ", s, perl = TRUE)   # drop [attribute] blocks (and pos="x,y")
  # Remove graph-level assignments, including the bounding box.
  s  <- gsub('[[:alnum:]_.]+\\s*=\\s*"(?:\\\\"|[^"])*"', " ", s, perl = TRUE)
  qn <- regmatches(s, gregexpr('"(?:\\\\"|[^"])*"', s, perl = TRUE))[[1]]  # quoted names (now pos-free, may hold spaces)
  qn <- vapply(qn, unquote_node, character(1L))
  qn <- gsub('\\"', '"', qn, fixed = TRUE)
  s  <- gsub('"(?:\\\\"|[^"])*"', " ", s, perl = TRUE)                     # remove quoted names from the string
  s  <- gsub("[{}]", " ", s)                                    # drop braces
  s  <- gsub("<->|<-|->|--", " ", s)                            # drop edge operators -> leave bare endpoints
  bn <- strsplit(s, "[[:space:];]+", perl = TRUE)[[1]]          # bare node names
  bn <- bn[nzchar(bn)]
  vals <- unique(c(as.character(qn), as.character(bn)))

  if( length(vals) > 0L ){

    return(vals)
  }

  nodes <- character(0L)

  return(nodes)
}


#' Read the dagitty model string
#'
#' @param dag A dagitty object.
#' @returns A single character string containing the graph model.
#' @noRd
dag_to_string <- function(dag){

  cands <- list(tryCatch(as.character(dag), error = function(e) {
    out <- NULL

    return(out)
  }),
                tryCatch(unclass(dag),      error = function(e) {
                  out <- NULL

                  return(out)
                }),
                tryCatch(format(dag),       error = function(e) {
                  out <- NULL

                  return(out)
                }))
  for( s in cands ){
    if( is.null(s) ) {
      next
    }
    s <- paste(as.character(s), collapse = "\n")
    if( grepl("\\{|->|--|<->|\\[", s) ){
      graph_type <- sub("^\\s*([[:alpha:]]+).*", "\\1", s)
      if( !graph_type %in% c("dag", "pdag", "mag") ){
        stop("causaliflower does not currently support dagitty graph type '",
             graph_type, "'.", call. = FALSE)
      }
      return(s)
    }
  }

  stop("Could not read a dagitty model string from the supplied object. ",
       "Please check that a dagitty object is supplied and try again.", call. = FALSE)
}


#' Remove surrounding double quotes
#'
#' @param x Character vector.
#' @returns `x` with one pair of surrounding double quotes removed where
#'   present.
#' @noRd
unquote_node <- function(x){

  text <- sub('^"(.*)"$', "\\1", x)

  return(text)
}


#' Shared helper for directed reachability over an edge table
#'
#' directed_reachable() walks the directed ("->") edges of a v/e/w edge table,
#' returning every node with a directed path to (direction = "to") or from
#' (direction = "from") any node in 'nodes'. Nodes in 'avoiding' block traversal.
#' Bidirected ("<->") edges are excluded. Undirected ("--") edges are followed
#' in both directions only when `possible = TRUE`. The seed nodes themselves
#' are excluded from the result.
#'
#' @param edges A data frame or data table of edges with columns v, e, w.
#' @param nodes Character vector of seed node names.
#' @param direction "to" for ancestors of nodes, "from" for descendants.
#' @param avoiding Optional character vector of node names paths may not pass through.
#' @param possible Include paths through undirected edges. Defaults to `FALSE`.
#' @returns Character vector of reachable node names.
#' @noRd
directed_reachable <- function(edges,
                               nodes,
                               direction = "to",
                               avoiding = NULL,
                               possible = FALSE
                               ){
  arrow_edges <- edges[ unlist(edges[, "e"]) == "->", ]
  tails <- as.character(unlist(arrow_edges[, "v"]))
  heads <- as.character(unlist(arrow_edges[, "w"]))

  if( isTRUE(possible) ){
    # Undirected edges permit reachability in either direction (Meek 1995).
    undirected_edges <- edges[ unlist(edges[, "e"]) == "--", ]
    uv <- as.character(unlist(undirected_edges[, "v"]))
    uw <- as.character(unlist(undirected_edges[, "w"]))
    tails <- c(tails, uv, uw)
    heads <- c(heads, uw, uv)
  }

  if( direction == "to" ){
    step_from <- heads # ancestors: walk each edge head to tail
    step_to   <- tails
  }else{
    step_from <- tails # descendants: walk each edge tail to head
    step_to   <- heads
  }

  reached  <- nodes[ !nodes %in% avoiding ]
  frontier <- reached

  while( length(frontier) > 0 ){

    next_nodes <- unique( step_to[ step_from %in% frontier ] )
    frontier   <- next_nodes[ !next_nodes %in% reached & !next_nodes %in% avoiding ]
    reached    <- c(reached, frontier)

  }

  nodes_out <- reached[ !reached %in% nodes ]

  return(nodes_out)
}


#' Helper for extract_unique_node_roles() and extract_node_roles()
#'
#' Defaults to instrument detection in dagitty::instrumentalVariables() (see van der Zander, Textor & Liskiewicz (2015, IJCAI)).
#' With usable_instruments = TRUE, candidates that are descendants of their own exposure are dropped.
#' Classical instrumental conditions consider instruments as an exogenous, upstream
#' source of variation in the exposure (see Greenland 2000, Int J Epidemiol 29:722-729; Hernan & Robins 2006, Epidemiology 17:360-372; Pearl 2009, Causality 2nd ed., section 7.4).
#' A post-treatment candidate is essentially redundant as an instrument, since it is certified by the graphical criterion only via a conditioning set that
#' already identifies the effect by back-door adjustment.
#'
#' @importFrom dagitty instrumentalVariables descendants
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param treatments Vector of treatment node names.
#' @param outcomes Vector of outcome node names.
#' @param usable_instruments When FALSE (default), returns every instrument dagitty identifies (the faithful graphical set). When TRUE, restricts to usable instruments by dropping candidates that are descendants of their own exposure.
#' @returns Vector of instrumental variable names.
#' @noRd
extract_instrumental_variables <- function(dag,
                                           treatments,
                                           outcomes,
                                           usable_instruments = FALSE
                                           ){
  instrumental_vars <- list()

  if( length(treatments) > 0 & length(outcomes) > 0 ){
    instrumental_vars <- suppressWarnings( lapply(seq_along(treatments), function(x){

      per_treatment <- lapply(seq_along(outcomes), function(y){
        ivs <- tryCatch(dagitty::instrumentalVariables(dag, exposure = treatments[x], outcome = outcomes[y]),
                        error = function(e) {
                          instruments <- list()

                          return(instruments)
                        })
        instruments <- unlist( lapply(ivs, function(iv) {
          instrument <- iv$I

          return(instrument)
        }) )

        return(instruments)
      })

      per_treatment <- unique( unlist(per_treatment) )
      if( usable_instruments ){
        per_treatment <- per_treatment[ !per_treatment %in% dagitty::descendants(dag, treatments[x]) ] # drop post-treatment candidates
      }

      return(per_treatment)

    }) )
  }
  instrumental_vars <- as.character(unique(unlist(instrumental_vars)))
  if( is.null(instrumental_vars) ) {
    instrumental_vars <- character()
  }

  return(instrumental_vars)
}


################################################################################
# graph_funcs.r
################################################################################
#' Read the supported graph type
#'
#' @param dag A dagitty object.
#' @return The graph type: dag, pdag or mag.
#' @noRd
dag_graph_type <- function(dag){
  if( !dagitty::is.dagitty(dag) ){
    stop("Invalid input detected: only dagitty objects may be supplied as 'dag' input. Please check the input and try again.",
         call. = FALSE)
  }
  model <- as.character(dag)
  graph_type <- sub("^\\s*([[:alpha:]]+).*", "\\1", model)
  if( !graph_type %in% c("dag", "pdag", "mag") ){
    stop("causaliflower does not currently support dagitty graph type '",
         graph_type, "'.", call. = FALSE)
  }

  return(graph_type)
}


#' Construct a graph
#'
#' construct_graph() constructs a dagitty object from an edge table.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty
#' @param edges A data frame of edges.
#' @param node_names Vector of node names to declare, including isolated nodes, e.g., c("Z1", "Z2", "Z3").
#' @param treatments Exposure node names to declare in the graph.
#' @param outcomes Outcome node names to declare in the graph.
#' @param latent_vec Character or vector of additional or already supplied latent (unobserved) variable names, e.g. "U" or c("U1", "U2", "M1").
#' @param adjusted_vec Adjusted node names to declare in the graph.
#' @param selected_vec Selected node names to declare in the graph.
#' @param graph_type Native graph type: `"dag"`, `"pdag"`, or `"mag"`.
#' @returns A dagitty object
#' @noRd
construct_graph <- function(edges,
                            node_names,
                            treatments,
                            outcomes,
                            latent_vec,
                            adjusted_vec = NULL,
                            selected_vec = NULL,
                            graph_type = "dag"
                            ){
  if( !graph_type %in% c("dag", "pdag", "mag") ){
    stop("Unsupported graph type for construction: ", graph_type, ".",
         call. = FALSE)
  }
  edges <- normalise_edge_frame(edges, "construct_graph()")
  edges <- edges[complete.cases(edges), , drop = FALSE]

  treatments <- clean_dag_names(treatments)
  outcomes <- clean_dag_names(outcomes)
  latent_vec <- clean_dag_names(latent_vec)
  adjusted_vec <- clean_dag_names(adjusted_vec)
  selected_vec <- clean_dag_names(selected_vec)
  all_nodes <- unique(c(clean_dag_names(node_names), edges$v, edges$w,
                        treatments, outcomes, latent_vec, adjusted_vec,
                        selected_vec))
  declarations <- vapply(all_nodes, function(node){
    attributes <- c(
      if( node %in% treatments ) { "exposure" },
      if( node %in% outcomes ) { "outcome" },
      if( node %in% latent_vec ) { "latent" },
      if( node %in% adjusted_vec ) { "adjusted" },
      if( node %in% selected_vec ) { "selected" }
    )
    quoted <- quote_dag_names(node)
    if( length(attributes) ){
      text <- paste0(quoted, " [", paste(attributes, collapse = ","), "]")

      return(text)
    }else{

      return(quoted)
    }
  }, character(1L), USE.NAMES = FALSE)
  edge_lines <- if( nrow(edges) ){
    paste(quote_dag_names(edges$v), edges$e, quote_dag_names(edges$w))
  }else{
    character()
  }

  dag <- dagitty::dagitty(paste(c(paste(graph_type, "{"), declarations,
                                  edge_lines, "}"),
                                collapse = " "))

  return(dag)
}


#' Standardise an edge table
#'
#' @param edges An edge table, or NULL for no edges.
#' @param caller Function name to include in errors.
#' @return Data frame with character columns v, e and w, with reverse arrows
#'   written as forward arrows.
#' @noRd
normalise_edge_frame <- function(edges, caller = "edge input"){
  if( is.null(edges) ) {
    edges_out <- empty_edge_frame()
    return(edges_out)
  }
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  if( nrow(edges) == 0L && ncol(edges) == 0L ){
    edges_out <- empty_edge_frame()
    return(edges_out)
  }
  if( ncol(edges) < 3L ){
    stop(caller, ": edge tables must contain v, e, and w columns.", call. = FALSE)
  }
  if( all(c("v", "e", "w") %in% names(edges)) ){
    edges <- edges[, c("v", "e", "w"), drop = FALSE]
  }else{
    edges <- edges[, seq_len(3L), drop = FALSE]
    names(edges) <- c("v", "e", "w")
  }
  edges[] <- lapply(edges, as.character)
  reverse <- !is.na(edges$e) & edges$e == "<-"
  if( any(reverse) ){
    old_v <- edges$v[reverse]
    edges$v[reverse] <- edges$w[reverse]
    edges$w[reverse] <- old_v
    edges$e[reverse] <- "->"
  }
  invalid <- !is.na(edges$e) & !edges$e %in% c("->", "<->", "--")
  if( any(invalid) ){
    stop(caller, ": unsupported edge operator(s): ",
         paste(unique(edges$e[invalid]), collapse = ", "), ".", call. = FALSE)
  }

  return(edges)
}


#' Create an empty edge data frame
#'
#' @return A zero-row data frame with character columns v, e and w.
#' @noRd
empty_edge_frame <- function(){
  edges <- data.frame(v = character(), e = character(), w = character(),
                      stringsAsFactors = FALSE)

  return(edges)
}


#' Collect unique non-missing node names
#'
#' @param x Node names supplied as a vector or list.
#' @return Character vector of unique, non-missing, non-empty names.
#' @noRd
clean_dag_names <- function(x){
  x <- as.character(unlist(x, use.names = FALSE))
  nodes <- unique(x[!is.na(x) & nzchar(x)])

  return(nodes)
}


#' Validate and quote node names
#'
#' @param x Node names supplied as a vector or list.
#' @param caller Function name to include in errors.
#' @return Character vector of names enclosed in double quotes.
#' @noRd
quote_dag_names <- function(x, caller = "graph construction"){
  x <- as.character(unlist(x, use.names = FALSE))
  if( length(x) == 0L ) {
    text <- character()
    return(text)
  }
  if( anyNA(x) || any(!nzchar(x)) ){
    stop(caller, ": node names cannot be missing or empty.", call. = FALSE)
  }
  invalid <- grepl('"', x, fixed = TRUE) |
    grepl(intToUtf8(92L), x, fixed = TRUE) |
    grepl("\r", x, fixed = TRUE) |
    grepl("\n", x, fixed = TRUE)
  if( any(invalid) ){
    stop(caller, ": node names cannot contain quotes, backslashes, or line breaks: ",
         paste(x[invalid], collapse = ", "), ".", call. = FALSE)
  }
  text <- paste0('"', x, '"')

  return(text)
}


#' Check a dagitty graph for directed cycles
#'
#' Removes nodes without remaining directed parents until none are left
#' (Kahn's algorithm). `dagitty::isAcyclic()` instead follows every directed
#' path and can take minutes on dense graphs. As there, only directed edges can
#' form a cycle.
#'
#' @param dag A dagitty object.
#' @return `TRUE` if the directed edges contain no cycle, otherwise `FALSE`.
#' @noRd
dag_is_acyclic <- function(dag){
  graph_edges <- dagitty::edges(dag)
  directed <- graph_edges[graph_edges$e == "->", , drop = FALSE]
  if( nrow(directed) == 0L ){
    acyclic <- TRUE
    return(acyclic)
  }
  from <- as.character(directed$v)
  to <- as.character(directed$w)
  nodes <- unique(c(from, to))
  indegree <- stats::setNames(tabulate(match(to, nodes), length(nodes)), nodes)
  children <- split(to, factor(from, levels = nodes))
  frontier <- nodes[indegree == 0L]
  removed <- 0L
  while( length(frontier) > 0L ){
    node <- frontier[[1L]]
    frontier <- frontier[-1L]
    removed <- removed + 1L
    for(child in children[[node]]){
      indegree[[child]] <- indegree[[child]] - 1L
      if( indegree[[child]] == 0L ){
        frontier <- c(frontier, child)
      }
    }
  }
  acyclic <- removed == length(nodes)

  return(acyclic)
}


#' Read and standardise graph edges
#'
#' @param dag A dagitty object.
#' @param caller Function name to include in errors.
#' @return Data frame of edges with character columns v, e and w.
#' @noRd
dag_edges <- function(dag, caller = "graph"){
  dag_graph_type(dag)
  edges_out <- normalise_edge_frame(dagitty::edges(dag), caller)

  return(edges_out)
}


#' Append nodes and edges to the native dagitty text
#'
#' Declares the new nodes with their statuses and adds the new edges before
#' the closing brace, leaving the graph name, bounding box, coordinates,
#' labels, notes and stored control points unchanged. Returns NULL when the
#' result does not hold exactly the expected nodes and edges.
#'
#' @param dag A dagitty object.
#' @param new_nodes Names of the nodes to declare.
#' @param new_edges Edge table (v, e, w) of the edges to add.
#' @param statuses Named list of status keywords for the new nodes.
#' @param caller Function label for messages.
#' @returns A dagitty object, or NULL.
#' @noRd
native_append_graph <- function(dag, new_nodes, new_edges, statuses = list(),
                                 caller){
  valid_names <- tryCatch({
    quote_dag_names(new_nodes)
    TRUE
  }, error = function(error) {
    matches <- FALSE

    return(matches)
  })
  if( !valid_names ){
    out <- NULL
    return(out)
  }
  existing_edges <- dag_edges(dag, caller)
  existing_keys <- edge_key(existing_edges)
  new_edges <- normalise_edge_frame(new_edges, caller)
  new_edges <- new_edges[complete.cases(new_edges), , drop = FALSE]
  new_edges <- new_edges[!duplicated(edge_key(new_edges)), , drop = FALSE]
  new_edges <- new_edges[!edge_key(new_edges) %in% existing_keys, , drop = FALSE]
  declarations <- vapply(new_nodes, function(node){
    status <- statuses[[node]]
    text <- paste0(native_layout_quote(node),
                   if( length(status) ) { paste0(" [", paste(status, collapse = ","), "]") } else { "" })

    return(text)
  }, character(1L), USE.NAMES = FALSE)
  edge_lines <- if( nrow(new_edges) ){
    paste(native_layout_quote(new_edges$v), new_edges$e,
          native_layout_quote(new_edges$w))
  }else{
    character()
  }
  if( length(declarations) == 0L && length(edge_lines) == 0L ){
    return(dag)
  }
  candidate <- tryCatch({
    text <- as.character(dag)
    closing <- native_layout_index(text, caller)$closing
    addition <- paste0("\n", paste(c(declarations, edge_lines), collapse = "\n"),
                       "\n")
    native_layout_value(dag, native_layout_edits(text, list(list(
      start = closing, end = closing - 1L, text = addition
    ))))
  }, error = function(error) {
    out <- NULL

    return(out)
  })
  if( is.null(candidate) ){
    out <- NULL
    return(out)
  }
  actual_edges <- tryCatch(dag_edges(candidate, caller), error = function(error) {
    out <- NULL

    return(out)
  })
  candidate_nodes <- tryCatch(names(candidate), error = function(error) {
    out <- NULL

    return(out)
  })
  expected_nodes <- c(names(dag), new_nodes)
  expected_keys <- unique(c(existing_keys, edge_key(new_edges)))
  actual_keys <- if( is.null(actual_edges) ) { character() } else { edge_key(actual_edges) }
  if( is.null(candidate_nodes) ||
      !setequal(candidate_nodes, expected_nodes) ||
      length(candidate_nodes) != length(expected_nodes) ||
      !setequal(actual_keys, expected_keys) ||
      length(unique(actual_keys)) != length(expected_keys) ){
    out <- NULL
    return(out)
  }

  return(candidate)
}


#' Create comparison keys for an edge table
#'
#' @param edges An edge table.
#' @return Character vector of edge keys, with symmetric endpoints in a
#'   consistent order.
#' @noRd
edge_key <- function(edges){
  edges <- normalise_edge_frame(edges)
  symmetric <- edges$e %in% c("<->", "--")
  left <- edges$v
  right <- edges$w
  if( any(symmetric) ){
    left[symmetric] <- pmin(edges$v[symmetric], edges$w[symmetric])
    right[symmetric] <- pmax(edges$v[symmetric], edges$w[symmetric])
  }
  keys <- paste(left, edges$e, right, sep = "\r")

  return(keys)
}


#' Quote node names for native graph text
#'
#' @param value Character vector of node names.
#' @return Character vector of quoted names with backslashes and double quotes
#'   escaped.
#' @noRd
native_layout_quote <- function(value){
  value <- gsub("\\", "\\\\", value, fixed = TRUE)
  text <- paste0('"', gsub('"', '\\"', value, fixed = TRUE), '"')

  return(text)
}


#' Index node and edge declarations
#'
#' @param text A dagitty model string.
#' @param caller Function name to include in errors.
#' @return A list of tokens, attribute blocks, node declarations and arrows,
#'   with the closing-brace position.
#' @noRd
native_layout_index <- function(text, caller){ # Nested readers update this call's parser state with <<-
  tokens <- native_layout_tokens(text, caller)
  values <- tokens$value
  position <- which(values == "{")[1L]
  if( is.na(position) ) {
    curve_layout_error(caller)
  }
  blocks <- list()
  declarations <- list()
  arrows <- list()
  decode <- function(value){
    if( startsWith(value, '"') ){
      value <- substr(value, 2L, nchar(value) - 1L)
      value <- gsub('\\\\(.)', '\\1', value, perl = TRUE)
    }

    return(value)
  }
  attributes <- function(kind, node = NULL, keys = character()){
    start <- position
    position <<- position + 1L
    fields <- list()
    while( position <= length(values) && values[position] != "]" ){
      if( values[position] %in% c(",", ";") ){
        position <<- position + 1L
        next
      }
      first <- position
      key <- decode(values[position])
      position <<- position + 1L
      if( position <= length(values) && values[position] == "=" ){
        position <<- position + 2L
      }
      fields[[length(fields) + 1L]] <- list(key = key, first = first,
                                           last = position - 1L)
    }
    if( position > length(values) ) {
      curve_layout_error(caller)
    }
    blocks[[length(blocks) + 1L]] <<- list(kind = kind, node = node,
                                          first = start, last = position,
                                          fields = fields, keys = keys)
    next_position <- (position <<- position + 1L)

    return(invisible(next_position))
  }
  endpoint <- function(){
    if( values[position] == "{" ){
      body()
      out <- NULL
      return(out)
    }
    if( values[position] %in% c("}", "[", "]", "=", ",", ";") ){
      curve_layout_error(caller)
    }
    result <- list(node = decode(values[position]), token = position)
    position <<- position + 1L

    return(result)
  }
  body <- function(){
    position <<- position + 1L
    while( position <= length(values) && values[position] != "}" ){
      if( values[position] %in% c(",", ";") ){
        position <<- position + 1L
        next
      }
      first <- endpoint()
      if( position <= length(values) && values[position] == "=" ){
        position <<- position + 2L
        next
      }
      edge <- FALSE
      keys <- character()
      source <- first
      grouped <- FALSE
      while( position <= length(values) &&
             values[position] %in% c("->", "<-", "<->", "--") ){
        edge <- TRUE
        operator <- values[position]
        marker <- position
        position <<- position + 1L
        target <- endpoint()
        if( !is.null(source) && !is.null(target) ){
          key <- edge_control_keys(source$node, operator, target$node)
          keys <- c(keys, key)
          arrows[[length(arrows) + 1L]] <<- list(key = key, operator = marker,
                                                source = source,
                                                target = target)
        }else{
          grouped <- TRUE
        }
        source <- target
      }
      if( !edge && !is.null(first) ){
        declarations[[length(declarations) + 1L]] <<- first
      }
      while( position <= length(values) && values[position] == "[" ){
        attributes(if( edge ) { "edge" } else { "node" },
                   if( !edge && !is.null(first) ) { first$node } else { NULL },
                   if( grouped ) { character() } else { keys })
      }
    }
    if( position > length(values) ) {
      curve_layout_error(caller)
    }
    closing <- position
    position <<- position + 1L

    return(closing)
  }
  closing <- body()
  layout_index <- list(tokens = tokens, blocks = blocks, declarations = declarations,
                       arrows = arrows, closing = tokens$start[closing])

  return(layout_index)
}


#' Locate tokens in native graph text
#'
#' @param text A dagitty model string.
#' @param caller Function name to include in errors.
#' @return Data frame of token values and their start and end positions,
#'   excluding comments.
#' @noRd
native_layout_tokens <- function(text, caller){ # Keep character offsets so later edits can preserve untouched graph text
  pattern <- paste0('"(?:\\\\.|[^"\\\\])*"|',
                    '/\\*[^*]*\\*+(?:[^/*][^*]*\\*+)*/|//[^\\r\\n]*|',
                    '#[^\\r\\n]*|<->|->|<-|--|[{}\\[\\]=,;]|',
                    '(?:(?!<->|->|<-|--)[^\\s{}\\[\\]=,;"])+')
  starts <- gregexpr(pattern, text, perl = TRUE)[[1L]]
  if( starts[1L] < 0L ) {
    curve_layout_error(caller)
  }
  lengths <- attr(starts, "match.length")
  values <- substring(text, starts, starts + lengths - 1L)
  comment <- startsWith(values, "//") | startsWith(values, "/*") |
    startsWith(values, "#")
  tokens <- data.frame(value = values[!comment], start = starts[!comment],
                       end = (starts + lengths - 1L)[!comment])

  return(tokens)
}


#' Restore graph attributes on edited text
#'
#' @param dag The original dagitty object.
#' @param text The edited model string.
#' @return The supplied text with the attributes of dag.
#' @noRd
native_layout_value <- function(dag, text){
  attributes(text) <- attributes(dag)

  return(text)
}


#' Apply edits to native graph text
#'
#' Edits are applied from right to left so their positions remain valid.
#'
#' @param text A dagitty model string.
#' @param edits List of edits, each with start, end and replacement text.
#' @return The edited model string.
#' @noRd
native_layout_edits <- function(text, edits){
  if( length(edits) == 0L ) {
    return(text)
  }
  order <- order(vapply(edits, `[[`, numeric(1L), "start"), decreasing = TRUE)
  for( i in order ){
    edit <- edits[[i]]
    text <- paste0(substr(text, 1L, edit$start - 1L), edit$text,
                   substring(text, edit$end + 1L))
  }

  return(text)
}


#' Rebuild dag
#'
#' rebuild_dag() rebuilds a dag using a dagitty object and data frame of edges input.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges coordinates dagitty
#' @param dag A dagitty object.
#' @param edges An edge table with columns v, e and w.
#' @param treatments Optional character vector of exposure names to declare; defaults to the declarations in dag.
#' @param outcomes Optional character vector of outcome names to declare; defaults to the declarations in dag.
#' @param latent_vec Optional character vector of latent names to declare; defaults to the declarations in dag.
#' @param extra_nodes Optional character vector of additional node names to declare beyond names(dag) and the edge endpoints.
#' @param adjusted_vec Optional adjusted node names; defaults to those in dag.
#' @param selected_vec Optional selected node names; defaults to those in dag.
#' @param .curve_caller Function name to include in curve-control messages.
#' @returns A dagitty object
#' @noRd
rebuild_dag <- function(dag,
                        edges,
                        treatments = NULL,
                        outcomes = NULL,
                        latent_vec = NULL,
                        extra_nodes = NULL,
                        adjusted_vec = NULL,
                        selected_vec = NULL,
                        .curve_caller = "connect_nodes"
                        ){
  if( is.null(latent_vec) ) {
    latent_vec <- unobserved(dag)
  }
  if( is.null(treatments) ) {
    treatments <- treatments(dag)
  }
  if( is.null(outcomes) ) {
    outcomes <- get_outcome_names(dag)
  }
  if( is.null(adjusted_vec) ) {
    adjusted_vec <- dagitty::adjustedNodes(dag)
  }
  if( is.null(selected_vec) ) {
    selected_vec <- selected_nodes(dag)
  }
  graph_type <- dag_graph_type(dag)

  controls <- edge_controls(dag, .curve_caller)
  edges <- normalise_edge_frame(edges, "rebuild_dag()")
  coordinates <- dagitty::coordinates(dag)
  node_names <- unique(c(names(dag), extra_nodes, edges$v, edges$w))
  node_names <- setdiff(node_names, c(treatments, outcomes, latent_vec))

  rebuilt <- construct_graph(
    edges, node_names, treatments, outcomes, latent_vec,
    adjusted_vec = adjusted_vec,
    selected_vec = selected_vec,
    graph_type = graph_type
  )
  coordinate_names <- intersect(names(coordinates$x), names(coordinates$y))
  coordinate_names <- intersect(coordinate_names, names(rebuilt))
  positioned <- coordinate_names[
    !is.na(coordinates$x[coordinate_names]) &
      !is.na(coordinates$y[coordinate_names])
  ]
  if( length(positioned) > 0L ){
    rebuilt <- write_node_coordinates(rebuilt, list(
      x = coordinates$x[positioned],
      y = coordinates$y[positioned]
    ), .curve_caller)
  }
  retained <- edge_control_keys(controls$v, controls$e, controls$w) %in%
    edge_control_keys(edges$v, edges$e, edges$w)
  rebuilt <- set_edge_controls(rebuilt, controls[retained, , drop = FALSE],
                                caller = .curve_caller)
  attributes(rebuilt) <- attributes(dag)

  return(rebuilt)
}


#' Read selected-node declarations
#'
#' @param dag A dagitty object.
#' @return Character vector of node names carrying the selected status.
#' @noRd
selected_nodes <- function(dag){
  lines <- strsplit(as.character(dag), "\n", fixed = TRUE)[[1L]]
  declaration <- '^\\s*("(?:\\\\"|[^"])*"|[^"[:space:]]+)\\s+\\[(.*)\\]\\s*$'
  quoted_value <- '"(?:\\\\.|[^"\\\\])*"'
  attribute_lines <- lines[grepl(declaration, lines, perl = TRUE)]
  attribute_text <- sub(declaration, "\\2", attribute_lines, perl = TRUE)
  attribute_keys <- gsub(quoted_value, "", attribute_text, perl = TRUE)
  selected_lines <- attribute_lines[
    grepl("(^|[,[:space:]])selected([,[:space:]]|$)", attribute_keys)
  ]
  if( length(selected_lines) == 0L ){
    nodes <- character()
    return(nodes)
  }
  node_text <- sub(declaration, "\\1", selected_lines, perl = TRUE)
  node_text <- sub('^"(.*)"$', "\\1", node_text)
  nodes <- unique(gsub('\\"', '"', node_text, fixed = TRUE))

  return(nodes)
}


#' Shared helper for formatting graph edges
#'
#' @importFrom data.table as.data.table
#' @param new_edges Data frame / data table of edges.
#' @returns A list of formatted edge strings, grouped by the first endpoint.
#' @noRd
print_edges_helper <- function(new_edges){
  if( nrow(new_edges) == 0L ){
    edges <- list()
    return(edges)
  }
  # collapse new_edges to a vector and output grouped by nodes
  unique_ancestors <- unique( new_edges[[1]] )
  num_unique_ancestors <- length(unique_ancestors)
  new_edges_list <- list()

  # Group edges by their first endpoint
  new_edges_list <- suppressWarnings( lapply(seq_len(num_unique_ancestors), function(x){
    ancestor_edges <- new_edges[ unlist(new_edges[,1]) %in% unlist(unique_ancestors)[x], ]

    return(ancestor_edges)
  }) )

  # Format each group as quoted edges for the caller to display
  new_edges_list <- lapply(seq_len(num_unique_ancestors), function(x){
    new_edges_list <- noquote(
      paste0( paste0( "'", sapply(seq_len(nrow(new_edges_list[[x]])), function(y){
        edge_text <- (new_edges_list[x] <- noquote( paste( new_edges_list[[x]][y,], collapse=" "  ) ))

        return(invisible(edge_text))
        }), "'", collapse=", " ), sep = "")
      )

    return(invisible(new_edges_list))
  })

  return(new_edges_list)
}


#' Read edge statements into a data frame
#'
#' @param edge_strings Character vector of edge statements.
#' @param caller Function name to include in errors.
#' @return Data frame of edges with columns v, e and w.
#' @noRd
parse_edge_strings <- function(edge_strings, caller = "edge input"){
  edges_out <- parse_graph_text(edge_strings, caller, allow_isolated = FALSE,
                    allow_empty = FALSE)$edges

  return(edges_out)
}


#' Remove matching edges from a graph
#'
#' Uses native text removal where possible, otherwise rebuilds the retained
#' edges.
#'
#' @param dag A dagitty object.
#' @param edges Edge table identifying the edges to remove.
#' @param caller Function name to include in errors.
#' @return A dagitty object without the matching edges.
#' @noRd
native_remove_edges <- function(dag, edges, caller){
  existing <- dag_edges(dag, paste0(caller, "()"))
  existing_keys <- edge_control_keys(existing$v, existing$e, existing$w)
  edges <- normalise_edge_frame(edges, paste0(caller, "()"))
  removed <- unique(edge_control_keys(edges$v, edges$e, edges$w))
  removed <- removed[removed %in% existing_keys]
  if( length(removed) == 0L ) {
    return(dag)
  }
  result <- tryCatch(native_edge_removal(dag, existing_keys, removed, caller),
                     error = function(error) {
                       out <- NULL

                       return(out)
                     })
  if( !is.null(result) ) {
    return(result)
  }
  dag_out <- rebuild_dag(dag, existing[!existing_keys %in% removed, , drop = FALSE],
                         .curve_caller = caller)

  return(dag_out)
}


#' Remove edges from native graph text
#'
#' @param dag A dagitty object.
#' @param existing_keys Control-point keys for all existing edges.
#' @param removed Control-point keys for the edges to remove.
#' @param caller Function name to include in errors.
#' @return The edited dagitty object, or NULL if the removal cannot preserve the
#'   remaining graph.
#' @noRd
native_edge_removal <- function(dag, existing_keys, removed, caller){
  quote_dag_names(names(dag), caller)
  text <- as.character(dag)
  index <- native_layout_index(text, caller)
  values <- index$tokens$value
  operators <- c("->", "<-", "<->", "--")
  decode <- function(value){
    if( startsWith(value, '"') ){
      value <- substr(value, 2L, nchar(value) - 1L)
      value <- gsub('\\\\(.)', '\\1', value, perl = TRUE)
    }

    return(value)
  }
  edits <- list()
  matched <- character()
  cut <- integer()
  ends <- list()
  for( arrow in index$arrows ){
    if( !arrow$key %in% removed ) {
      next
    }
    matched <- c(matched, arrow$key)
    cut <- c(cut, arrow$operator)
    ends <- c(ends, list(arrow$source, arrow$target))
  }
  if( !setequal(matched, removed) ) {
    out <- NULL
    return(out)
  }
  for( block in index$blocks ){
    if( block$kind != "edge" || !any(block$keys %in% removed) ) {
      next
    }
    if( length(block$keys) != 1L ) {
      out <- NULL
      return(out)
    }
    edits[[length(edits) + 1L]] <- list(
      start = index$tokens$start[block$first],
      end = index$tokens$end[block$last],
      text = ""
    )
  }
  linked <- function(i){
    matches <- i >= 1L && i <= length(values) && values[i] %in% operators &&
      !(i %in% cut)

    return(matches)
  }
  present <- vapply(index$declarations, `[[`, character(1L), "node")
  for( i in which(values %in% operators) ){
    if( i %in% cut ) {
      next
    }
    present <- c(present, decode(values[i - 1L]), decode(values[i + 1L]))
  }
  orphans <- list()
  for( end in ends ){
    if( linked(end$token - 1L) || linked(end$token + 1L) ) {
      next
    }
    orphans <- c(orphans, list(end))
  }
  restore <- integer()
  for( orphan in orphans ){
    if( orphan$node %in% present ) {
      next
    }
    present <- c(present, orphan$node)
    restore <- c(restore, orphan$token)
  }
  gone <- sort(unique(c(cut, unlist(lapply(orphans, `[[`, "token")))))
  opening <- which(!(gone - 1L) %in% gone)
  closing <- which(!(gone + 1L) %in% gone)
  for( i in seq_along(opening) ){
    span <- gone[opening[i]:closing[i]]
    declared <- if( length(span) > 1L ){
      paste(values[span[span %in% restore]], collapse = "\n")
    }else{
      "\n"
    }
    edits[[length(edits) + 1L]] <- list(
      start = index$tokens$start[span[1L]],
      end = index$tokens$end[span[length(span)]],
      text = declared
    )
  }
  result <- native_layout_value(dag, native_layout_edits(text, edits))
  after <- dag_edges(result, caller)
  nodes <- names(dag)
  before_coordinates <- dagitty::coordinates(dag)
  after_coordinates <- dagitty::coordinates(result)
  controls <- edge_controls(dag, caller)
  controls <- controls[
    !edge_control_keys(controls$v, controls$e, controls$w) %in% removed, ,
    drop = FALSE
  ]
  actual <- edge_controls(result, caller)
  found <- match(edge_control_keys(controls$v, controls$e, controls$w),
                 edge_control_keys(actual$v, actual$e, actual$w))
  if( !setequal(edge_control_keys(after$v, after$e, after$w),
                setdiff(existing_keys, removed)) ||
      !setequal(names(result), nodes) ||
      length(names(result)) != length(nodes) ||
      !identical(after_coordinates$x[nodes], before_coordinates$x[nodes]) ||
      !identical(after_coordinates$y[nodes], before_coordinates$y[nodes]) ||
      nrow(actual) != nrow(controls) || anyNA(found) ||
      !identical(as.numeric(actual$x[found]), as.numeric(controls$x)) ||
      !identical(as.numeric(actual$y[found]), as.numeric(controls$y)) ||
      !identical(treatments(result), treatments(dag)) ||
      !identical(get_outcome_names(result), get_outcome_names(dag)) ||
      !identical(unobserved(result), unobserved(dag)) ||
      !identical(dagitty::adjustedNodes(result), dagitty::adjustedNodes(dag)) ||
      !identical(selected_nodes(result), selected_nodes(dag)) ||
      !identical(dag_graph_type(result), dag_graph_type(dag)) ){
    out <- NULL
    return(out)
  }

  return(result)
}


#' Remove nodes and their incident edges from native graph text
#'
#' @param dag A dagitty object.
#' @param nodes Existing node names to remove.
#' @param caller Function name to include in errors.
#' @return A dagitty object containing the remaining nodes and edges.
#' @noRd
native_remove_nodes <- function(dag, nodes, caller){
  existing <- dag_edges(dag, paste0(caller, "()"))
  incident <- existing$v %in% nodes | existing$w %in% nodes
  dag <- native_remove_edges(dag, existing[incident, , drop = FALSE], caller)
  remaining <- setdiff(names(dag), nodes)
  text <- as.character(dag)
  index <- native_layout_index(text, caller)
  edits <- list()
  for( declaration in index$declarations ){
    if( !declaration$node %in% nodes ){
      next
    }
    edits[[length(edits) + 1L]] <- list(
      start = index$tokens$start[declaration$token],
      end = index$tokens$end[declaration$token], text = ""
    )
  }
  for( block in index$blocks ){
    if( block$kind != "node" || !block$node %in% nodes ){
      next
    }
    edits[[length(edits) + 1L]] <- list(
      start = index$tokens$start[block$first],
      end = index$tokens$end[block$last], text = ""
    )
  }
  result <- native_layout_value(dag, native_layout_edits(text, edits))
  labels <- attr(result, "labels", exact = TRUE)
  if( !is.null(labels) && !is.null(names(labels)) ){
    attr(result, "labels") <- labels[!names(labels) %in% nodes]
  }

  actual_edges <- dag_edges(result, paste0(caller, "()"))
  before_coordinates <- dagitty::coordinates(dag)
  after_coordinates <- dagitty::coordinates(result)
  before_status <- dagitty::dagitty(text)
  after_status <- dagitty::dagitty(as.character(result))
  if( !setequal(names(result), remaining) ||
      !setequal(edge_key(actual_edges),
                edge_key(existing[!incident, , drop = FALSE])) ||
      !identical(as.numeric(after_coordinates$x[remaining]),
                 as.numeric(before_coordinates$x[remaining])) ||
      !identical(as.numeric(after_coordinates$y[remaining]),
                 as.numeric(before_coordinates$y[remaining])) ||
      !setequal(dagitty::exposures(after_status),
                setdiff(dagitty::exposures(before_status), nodes)) ||
      !setequal(dagitty::outcomes(after_status),
                setdiff(dagitty::outcomes(before_status), nodes)) ||
      !setequal(dagitty::latents(after_status),
                setdiff(dagitty::latents(before_status), nodes)) ||
      !setequal(dagitty::adjustedNodes(after_status),
                setdiff(dagitty::adjustedNodes(before_status), nodes)) ||
      !setequal(selected_nodes(after_status),
                setdiff(selected_nodes(before_status), nodes)) ||
      !identical(dag_graph_type(result), dag_graph_type(dag)) ){
    stop(caller, "(): the remaining nodes and edges could not be preserved.",
         call. = FALSE)
  }

  return(result)
}


################################################################################
# helper_funcs.r
################################################################################
#' Build a dagitty object from edge text
#'
#' `string_to_dag()` converts edge text to a dagitty object.
#'
#' Input may contain edges separated by newlines or semicolons, edge chains, and isolated node declarations. Supported operators are \code{->},
#' \code{<-}, \code{<->}, and \code{--}.
#' @param edges_string Character vector containing edge or node declarations.
#' @returns A dagitty object.
#' @examples
#' dag <- string_to_dag(c(
#'   "rainfall -> irrigation",
#'   "rainfall -> yield",
#'   "irrigation -> yield"
#' ))
#' get_edges(dag)
#'
#' @export
string_to_dag <- function(edges_string){
  parsed <- parse_graph_text(edges_string, "string_to_dag()",
                              allow_isolated = TRUE, allow_empty = TRUE)
  node_names <- unique(c(parsed$edges$v, parsed$edges$w, parsed$isolated))
  dag <- construct_graph(parsed$edges, node_names, NULL, NULL, NULL)

  return(dag)
}


#' Read edges and isolated nodes from text
#'
#' @param text Character vector of edge statements and optional isolated nodes.
#' @param caller Function name to include in errors.
#' @param allow_isolated Whether statements may contain a node without an edge.
#' @param allow_empty Whether an empty input is accepted.
#' @return A list containing an edges data frame and a character vector of
#'   isolated nodes.
#' @noRd
parse_graph_text <- function(text, caller, allow_isolated = FALSE,
                              allow_empty = FALSE){
  if( !is.character(text) ){
    stop(caller, ": edge input must be a character vector.", call. = FALSE)
  }
  if( anyNA(text) ){
    stop(caller, ": edge strings cannot contain NA.", call. = FALSE)
  }

  lines <- split_graph_statements(text, caller)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if( length(lines) == 0L ){
    if( allow_empty ){
      parsed <- list(edges = empty_edge_frame(), isolated = character())
      return(parsed)
    }
    stop(caller, ": no edges were supplied.", call. = FALSE)
  }

  edge_rows <- list()
  isolated <- character()
  for( line in lines ){
    matches <- find_graph_operators(line, caller)
    if( isTRUE(matches$invalid_marker) ){ # a stray "<" or ">" outside any recognised operator
      stop(caller, ": could not parse '", line,
           "'. Expected 'node op node'.", call. = FALSE)
    }
    if( length(matches$position) == 0L ){
      if( !allow_isolated ){
        stop(caller, ": could not parse '", line,
             "'. Expected 'node op node'.", call. = FALSE)
      }
      isolated <- c(isolated, strip_dag_quotes(line, caller))
      next
    }

    operators <- matches$operator
    starts <- c(1L, matches$position + matches$length)
    ends <- c(matches$position - 1L, nchar(line))
    nodes <- vapply(seq_along(starts), function(i){
      node <- strip_dag_quotes(substr(line, starts[i], ends[i]), caller)

      return(node)
    }, character(1L))

    for( i in seq_along(operators) ){
      v <- nodes[i]
      w <- nodes[i + 1L]
      e <- operators[i]
      if( e == "<-" ){
        edge_rows[[length(edge_rows) + 1L]] <- c(v = w, e = "->", w = v)
      }else{
        edge_rows[[length(edge_rows) + 1L]] <- c(v = v, e = e, w = w)
      }
    }
  }

  edges <- if( length(edge_rows) ){
    out <- as.data.frame(do.call(rbind, edge_rows), stringsAsFactors = FALSE)
    names(out) <- c("v", "e", "w")
    out
  }else{
    empty_edge_frame()
  }
  parsed <- list(edges = edges, isolated = unique(isolated))

  return(parsed)
}


#' Split graph text into statements
#'
#' @param text Character vector of graph statements.
#' @param caller Function name to include in errors.
#' @return Character vector split at semicolons and line breaks outside quoted
#'   names.
#' @noRd
split_graph_statements <- function(text, caller){
  chars <- strsplit(paste(text, collapse = "\n"), "", fixed = TRUE)[[1L]]
  if( length(chars) == 0L ) {
    statements_out <- character()
    return(statements_out)
  }

  statements <- character()
  current <- character()
  in_quote <- FALSE
  for( char in chars ){
    if( identical(char, '"') ){
      in_quote <- !in_quote
      current <- c(current, char)
    }else if( !in_quote && char %in% c(";", "\r", "\n") ){
      statements <- c(statements, paste0(current, collapse = ""))
      current <- character()
    }else{
      current <- c(current, char)
    }
  }
  if( in_quote ){
    stop(caller, ": unmatched quote in graph text.", call. = FALSE)
  }
  statements_out <- c(statements, paste0(current, collapse = ""))

  return(statements_out)
}


#' Locate edge operators outside quoted names
#'
#' @param line A single graph statement.
#' @param caller Function name to include in errors.
#' @return A list of operators, positions and lengths, with a flag for stray
#'   angle brackets.
#' @noRd
find_graph_operators <- function(line, caller){
  chars <- strsplit(line, "", fixed = TRUE)[[1L]]
  positions <- integer()
  lengths <- integer()
  operators <- character()
  invalid_marker <- FALSE
  in_quote <- FALSE
  i <- 1L
  while( i <= length(chars) ){
    if( identical(chars[i], '"') ){
      in_quote <- !in_quote
      i <- i + 1L
      next
    }
    if( !in_quote ){
      remaining <- paste0(chars[i:length(chars)], collapse = "")
      candidates <- c("<->", "->", "<-", "--")
      operator <- candidates[vapply(candidates, function(candidate){
        matches <- startsWith(remaining, candidate)

        return(matches)
      }, logical(1L))]
      if( length(operator) > 0L ){
        operator <- operator[1L]
        positions <- c(positions, i)
        lengths <- c(lengths, nchar(operator))
        operators <- c(operators, operator)
        i <- i + nchar(operator)
        next
      }
      if( chars[i] %in% c("<", ">") ) {
        invalid_marker <- TRUE
      }
    }
    i <- i + 1L
  }
  if( in_quote ){
    stop(caller, ": unmatched quote in graph text.", call. = FALSE)
  }
  operators_out <- list(position = positions, length = lengths, operator = operators,
                        invalid_marker = invalid_marker)

  return(operators_out)
}


#' Read a node name from edge text
#'
#' @param x A single node name, optionally enclosed in double quotes.
#' @param caller Function name to include in errors.
#' @return The node name with surrounding whitespace and quotes removed.
#' @noRd
strip_dag_quotes <- function(x, caller){
  x <- trimws(x)
  starts <- startsWith(x, '"')
  ends <- endsWith(x, '"')
  if( xor(starts, ends) ){
    stop(caller, ": unmatched quote in node name '", x, "'.", call. = FALSE)
  }
  if( starts && ends ) {
    x <- substr(x, 2L, nchar(x) - 1L)
  }
  if( !nzchar(x) ){
    stop(caller, ": node names cannot be empty.", call. = FALSE)
  }
  if( grepl('"', x, fixed = TRUE) ){
    stop(caller, ": embedded quotes are not supported in node names.", call. = FALSE)
  }

  return(x)
}


################################################################################
# utility_funcs.r
################################################################################
#' Read node names with a native dagitty status
#'
#' Reads status keywords from the model string.
#'
#' @param dag A dagitty object.
#' @param role One of "treatment", "outcome", or "unobserved".
#' @returns Character vector of node names, or character(0) if none.
#' @noRd
get_dagitty_role <- function(dag, role){

  tags <- switch(role,
                 treatment   = c("exposure", "e"),
                 outcome    = c("outcome", "o"),
                 unobserved = c("latent", "unobserved", "u"),
                 stop("Unknown role: ", role, ". Expected treatment, outcome, or unobserved."))
  lowercase_node_names <- node_status_attributes(dag_to_string(dag))
  matched_node_names   <- names(lowercase_node_names)[vapply(lowercase_node_names, function(a) {
    matches <- any(tags %in% a)

    return(matches)
  }, logical(1L))]

  if( length(matched_node_names) > 0L ){
    nodes <- as.character(matched_node_names)
    return(nodes)
  }

  nodes <- character(0L)

  return(nodes)
}


#' Read node status keywords from graph text
#'
#' Combines bare attribute keywords from repeated declarations of each node.
#' Node names and keywords retain their case.
#'
#' @param s A single character string containing the graph model.
#' @returns A list of keyword vectors indexed by node name.
#' @noRd
node_status_attributes <- function(s){

  quoted <- '"(?:\\\\"|[^"])*"'
  pat  <- paste0('^(', quoted, '|[A-Za-z0-9_.]+)\\s*\\[([^][]*)\\]$')
  scan <- paste0('(?:', quoted, '|[A-Za-z0-9_.]+)\\s*\\[[^][]*\\]|', quoted) # a quoted name without attributes is consumed whole
  m   <- gregexpr(scan, s, perl = TRUE)[[1]]
  if( m[1] == -1L ) {
    nodes <- list()
    return(nodes)
  }
  starts <- m; lens <- attr(m, "match.length")
  out <- list()
  for( i in seq_along(starts) ){
    chunk <- substr(s, starts[i], starts[i] + lens[i] - 1L)
    if( !grepl(pat, chunk, perl = TRUE) ) {
      next
    }
    nm   <- unquote_node(sub(pat, "\\1", chunk, perl = TRUE))
    att  <- sub(pat, "\\2", chunk, perl = TRUE)
    att  <- gsub('[A-Za-z0-9_.]+\\s*=\\s*"[^"]*"', "", att)     # drop key="value" pairs, e.g. pos="x,y"
    toks <- trimws(strsplit(att, "[,[:space:]]+", perl = TRUE)[[1]]) # tags are case-sensitive, matching dagitty's grammar
    out[[nm]] <- unique( c( out[[nm]], toks[nzchar(toks)] ) ) # merge tokens across repeated declarations of one node
  }

  return(out)
}


#' Shared helper for reading outcome names from native node statuses
#'
#' @param dag A dagitty object.
#' @returns Character vector of outcome names, or character(0) if none.
#' @noRd
get_outcome_names <- function(dag){

  nodes <- get_dagitty_role(dag, "outcome")

  return(nodes)
}


#' Helper for instruments()
#'
#' instrumental_variables_helper() delegates to dagitty::instrumentalVariables()
#' (van der Zander, Textor & Liskiewicz 2015, IJCAI), returning a data.table of
#' instruments and their conditioning sets (one row per treatment, outcome, instrument).
#' Conditional instruments are included. With
#' usable_instruments = TRUE, candidates that are descendants of their own exposure are
#' dropped, per exposure (Greenland 2000, Int J Epidemiol 29:722-729; Hernan & Robins
#' 2006, Epidemiology 17:360-372; Pearl 2009, section 7.4) - see
#' extract_instrumental_variables() for the full rationale.
#'
#' @importFrom dagitty instrumentalVariables descendants
#' @importFrom data.table as.data.table data.table setorder
#' @param dag A dagitty object
#' @param treatments A vector of treatment node names
#' @param outcomes A vector of outcome node names
#' @param usable_instruments When FALSE (default), returns every instrument dagitty identifies. When TRUE, drops candidates that are descendants of their own exposure.
#' @returns A data.table with one row per identified instrument (columns: treatment, outcome, instrument, conditioning_set; conditioning_set lists the variables to condition on, comma-separated, empty if unconditional).
#' @noRd
instrumental_variables_helper <- function(dag,
                                          treatments,
                                          outcomes,
                                          usable_instruments = FALSE
                                          ){
  instrumental_variables <- suppressWarnings( lapply(seq_along(treatments), function(x){

    per_outcome <- lapply(seq_along(outcomes), function(y){
      ivs <- tryCatch(dagitty::instrumentalVariables(dag, exposure = treatments[x], outcome = outcomes[y]),
                      error = function(e) {
                        instruments <- list()

                        return(instruments)
                      })
      instrument_rows <- lapply(ivs, function(iv){
        instrument_row <- c( treatment        = treatments[x],
                     outcome          = outcomes[y],
                     instrument       = iv$I,
                     conditioning_set = paste0(iv$Z, collapse = ", ") )

        return(instrument_row)
      })

      return(instrument_rows)
    })

    per_treatment <- unlist(per_outcome, recursive = FALSE)
    if( usable_instruments ){
      post_treatment <- dagitty::descendants(dag, treatments[x]) # drop post-treatment candidates
      per_treatment  <- Filter(function(iv) {
        matches <- !iv[["instrument"]] %in% post_treatment

        return(matches)
      }, per_treatment)
    }

    return(per_treatment)

  }) )

  instrumental_variables <- Filter(Negate(anyNA), unlist(instrumental_variables, recursive = FALSE))

  if( length(instrumental_variables) == 0 ){
    instruments <- data.table::data.table( treatment        = character(0),
                                           outcome          = character(0),
                                           instrument       = character(0),
                                           conditioning_set = character(0) )
    return(instruments)
  }

  instrumental_variables <- unique( data.table::as.data.table( do.call( rbind, instrumental_variables ) ) )
  instrumental_variables <- data.table::setorder( instrumental_variables, treatment, outcome, instrument )

  return(instrumental_variables)
}
