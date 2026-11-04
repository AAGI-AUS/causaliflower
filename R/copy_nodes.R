# ==============================================================================
# copy_nodes.R


# copy_nodes() and its internal helpers
# ==============================================================================


#' Copy nodes in a dagitty object
#'
#' `copy_nodes()` keeps each selected node, adds copies that inherit its direct
#' incoming and outgoing arrows and its bidirected and undirected edges, and
#' connects the original and copies in temporal positions. An anchor may
#' reroute inherited arrows as described below. `type = "chain"` places the
#' copies in a path instead of inheriting the original's edges.
#'
#' @param dag A dagitty object containing a directed acyclic graph.
#' @param nodes Existing node names to copy.
#' @param suffix Optional suffix for the existing nodes. When `copy_suffix` is
#'   `NULL`, numbering covers the renamed original and all copies.
#' @param copy_suffix Optional character vector or list of suffixes for copies.
#'   A character vector places copies in successive temporal positions.
#'   Suffixes grouped in one list element place copies at the same position.
#' @param anchor Optional existing node connected to one end of the copies. An
#'   anchor that is a descendant of the copied node, or is unrelated to it, is
#'   a parent of the first copy. An ancestor of the copied node is a child of
#'   the last copy.
#' @param repeats Optional number of copies, excluding each original.
#' @param type Optional `"chain"` to place the copies in the anchor's path
#'   instead of inheriting the original's edges. Defaults to `NULL`.
#' @param renew_coords Regenerate the coordinates of every node after copying.
#'   Defaults to `TRUE`. Stored curves are carried to the new coordinates. With
#'   `FALSE`, only the copies are placed and existing coordinates are kept. If
#'   a copy would sit on an existing node, every coordinate is regenerated.
#' @param print_edges Print added and rerouted edges in the console. Defaults to
#'   `FALSE`.
#'
#' @details
#' `suffix` renames the existing node. If `copy_suffix` is `NULL`, numbering
#' covers the renamed original and its copies. With `copy_suffix`, a character
#' vector gives successive copy positions and values grouped in one list element
#' are parallel. If neither suffix is supplied, one copy is added by default
#' with `_copy` appended to the node name. With `repeats` greater than one,
#' the copies use `_copy1`, `_copy2`, and so on. `repeats` counts copies,
#' excluding the original, and must agree with
#' the number of explicit copy suffixes.
#'
#' The copies are connected to each other in order, and the original joins one
#' end of that order. The copies follow an anchor that is a descendant of the
#' copied node or unrelated to it, and precede an anchor that is an ancestor.
#' An inherited arrow that conflicts with this placement is rerouted for the
#' copy that joins the anchor and is not added to the other copies.
#'
#' `type = "chain"` drops the inherited edges and places the copies in the
#' path leading out of the anchor, or into the anchor when it is an ancestor of
#' the copied node. Without an anchor the copied node itself is the anchor.
#' Each arrow the copies take over is removed from the anchor and reported as a
#' reroute. The graph name, bounding box, node statuses, labels, notes, and R
#' attributes are unchanged apart from those attached to the removed arrows.
#' Coordinates and stored curves follow `renew_coords`.
#'
#' @returns A dagitty object.
#' @examples
#' dag <- dagitty::dagitty("dag { T -> M -> Y }")
#'
#' sequential <- copy_nodes(
#'   dag, "M", copy_suffix = c("survey", "harvest"), print_edges = FALSE
#' )
#'
#' parallel <- copy_nodes(
#'   dag, "M", copy_suffix = list(c("survey", "sensor")),
#'   print_edges = FALSE
#' )
#'
#' chained <- copy_nodes(dag, "M", type = "chain", print_edges = FALSE)
#' @export
copy_nodes <- function(dag,
                       nodes,
                       suffix = NULL,
                       copy_suffix = NULL,
                       anchor = NULL,
                       repeats = NULL,
                       type = NULL,
                       renew_coords = TRUE,
                       print_edges = FALSE
                       ){
  caller <- "copy_nodes()"

  if( !dagitty::is.dagitty(dag) ){
    stop(caller, ": dag must be a dagitty object.", call. = FALSE)
  }
  if( dag_graph_type(dag) != "dag" ){
    stop(caller, ": dag must use the dag graph type.", call. = FALSE)
  }
  if( !isTRUE(dag_is_acyclic(dag)) ){
    stop(caller, ": dag must be acyclic before nodes are copied.",
         call. = FALSE)
  }
  if( !is.logical(print_edges) || length(print_edges) != 1L ||
      is.na(print_edges) ){
    stop(caller, ": print_edges must be TRUE or FALSE.", call. = FALSE)
  }
  if( !is.logical(renew_coords) || length(renew_coords) != 1L ||
      is.na(renew_coords) ){
    stop(caller, ": renew_coords must be TRUE or FALSE.", call. = FALSE)
  }
  if( !is.null(type) && !identical(type, "chain") ){
    stop(caller, ": type must be NULL or \"chain\".", call. = FALSE)
  }

  nodes <- copy_validate_nodes(nodes, names(dag), "nodes")
  if( anyDuplicated(nodes) ){
    stop(caller, ": nodes must not contain duplicates.", call. = FALSE)
  }

  if( !is.null(anchor) ){
    anchor <- copy_validate_nodes(anchor, names(dag), "anchor")
    if( length(anchor) != 1L ){
      stop(caller, ": anchor must name one existing node.", call. = FALSE)
    }
    if( anchor %in% nodes ){
      stop(caller, ": anchor cannot also be one of the copied nodes.",
           call. = FALSE)
    }
  }

  specification <- copy_name_specification(
    nodes = nodes,
    suffix = suffix,
    copy_suffix = copy_suffix,
    repeats = repeats
  )
  original_names <- specification$original_names
  copy_names <- specification$copy_names
  copy_groups <- specification$copy_groups
  all_copies <- unname(unlist(copy_names, recursive = TRUE,
                              use.names = FALSE))

  prospective_names <- c(setdiff(names(dag), nodes), original_names,
                         all_copies)
  duplicated_names <- unique(prospective_names[duplicated(prospective_names)])
  if( length(duplicated_names) > 0L ){
    stop(caller, ": generated node names collide with existing or generated names: ",
         paste(duplicated_names, collapse = ", "), ".", call. = FALSE)
  }

  rename_map <- stats::setNames(original_names, nodes)
  original_edges <- dag_edges(dag, caller)
  placements <- stats::setNames(rep("after", length(nodes)), nodes)
  if( !is.null(anchor) ){
    placements <- vapply(nodes, function(node){
      ancestors <- directed_reachable(original_edges, node, direction = "to")
      if( anchor %in% ancestors ) {
        value <- "before"

        return(value)
      } else {
        value <- "after"

        return(value)
      }
    }, character(1L), USE.NAMES = TRUE)
  }

  renamed <- copy_rename_native_nodes(dag, rename_map, caller)
  renamed_edges <- original_edges
  renamed_edges$v <- copy_remap_names(renamed_edges$v, rename_map)
  renamed_edges$w <- copy_remap_names(renamed_edges$w, rename_map)

  generated <- copy_generate_edges(
    original_edges = original_edges,
    nodes = nodes,
    renamed = rename_map,
    copy_names = copy_names,
    copy_groups = copy_groups,
    anchor = anchor,
    placements = placements,
    type = type
  )
  new_edges <- generated$edges
  removed_edges <- generated$removed
  rerouted <- generated$rerouted

  if( nrow(removed_edges) > 0L ){
    kept <- !edge_key(renamed_edges) %in% edge_key(removed_edges)
    renamed_edges <- renamed_edges[kept, , drop = FALSE]
    renamed <- native_remove_edges(renamed, removed_edges, "copy_nodes")
  }

  candidate <- copy_append_graph(
    renamed,
    new_nodes = all_copies,
    new_edges = new_edges,
    source_for_copy = specification$source_for_copy,
    source_dag = dag,
    caller = caller
  )

  expected_edges <- unique(rbind(renamed_edges, new_edges))
  actual_edges <- dag_edges(candidate, caller)
  if( !setequal(edge_key(expected_edges), edge_key(actual_edges)) ){
    stop(caller, ": dagitty could not represent the requested edges exactly.",
         call. = FALSE)
  }
  if( !isTRUE(dag_is_acyclic(candidate)) ){
    stop(caller, ": the proposed copies would create a directed cycle; ",
         "the graph was not changed.", call. = FALSE)
  }

  new_coordinates <- NULL
  if( !isTRUE(renew_coords) ){
    existing_coordinates <- dagitty::coordinates(renamed)
    candidate_coordinates <- dagitty::coordinates(candidate)
    new_coordinates <- renew_coords(
      candidate,
      new_node_names = all_copies,
      coordinates = candidate_coordinates
    )
    positioned <- intersect(names(existing_coordinates$x),
                            names(existing_coordinates$y))
    positioned <- positioned[
      is.finite(existing_coordinates$x[positioned]) &
        is.finite(existing_coordinates$y[positioned])
    ]
    new_coordinates$x[positioned] <- existing_coordinates$x[positioned]
    new_coordinates$y[positioned] <- existing_coordinates$y[positioned]
    if( copy_coincident_copies(new_coordinates, all_copies) ){
      warning(caller,
              ": the copies could not be placed apart from existing nodes; ",
              "coordinates were regenerated for every node.", call. = FALSE)
      new_coordinates <- NULL
    }
  }
  if( is.null(new_coordinates) ){
    new_coordinates <- renew_coords(candidate)
  }
  candidate <- set_coords_preserving_curves(
    candidate, new_coordinates, caller = "copy_nodes"
  )

  if( isTRUE(print_edges) ){
    copy_print_reroutes(rerouted)
    printed <- print_edges_helper(new_edges)
    cat(paste("c(", paste(unlist(printed), collapse = ",\n\n"), ")",
              sep = "\n", collapse = ""))
    message("\nPrinted new edges - copy and paste to use as a vector object.\n")
  }

  return(candidate)
}


## Copy inputs and names ##
#' Check node names supplied for copying
#'
#' @param value Node names to check.
#' @param available Node names already present in the graph.
#' @param argument Argument name for messages.
#' @return The supplied node-name vector.
#' @noRd
copy_validate_nodes <- function(value, available, argument){
  caller <- "copy_nodes()"
  if( !is.character(value) || length(value) == 0L || anyNA(value) ||
      any(!nzchar(value)) ){
    stop(caller, ": ", argument,
         " must contain non-missing, non-empty node names.", call. = FALSE)
  }
  missing <- unique(value[!value %in% available])
  if( length(missing) > 0L ){
    stop(caller, ": ", argument, " not in the graph: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }

  return(value)
}


#' Check the number of node copies
#'
#' @param repeats Optional positive whole number of copies per original node.
#' @return The copy count as an integer, or `NULL` when omitted.
#' @noRd
copy_validate_repeats <- function(repeats){
  if( is.null(repeats) ) {
    out <- NULL
    return(out)
  }
  if( !is.numeric(repeats) || length(repeats) != 1L ||
      !is.finite(repeats) || repeats < 1 || repeats != as.integer(repeats) ){
    stop("copy_nodes(): repeats must be one positive whole number.",
         call. = FALSE)
  }
  repeats_out <- as.integer(repeats)

  return(repeats_out)
}


#' Group copy suffixes by temporal position
#'
#' @param value Character suffixes or a list of suffix groups.
#' @param argument Argument name for messages.
#' @return A list of character suffix vectors, one per temporal position.
#' @noRd
copy_suffix_positions <- function(value, argument){
  if( is.character(value) ){
    if( length(value) == 0L || anyNA(value) || any(!nzchar(value)) ){
      stop("copy_nodes(): ", argument,
           " must contain non-missing, non-empty suffixes.", call. = FALSE)
    }
    positions_out <- lapply(as.character(value), function(x) {

      return(x)
    })
    return(positions_out)
  }
  if( !is.list(value) || length(value) == 0L ){
    stop("copy_nodes(): ", argument,
         " must be a character vector or a list of character positions.",
         call. = FALSE)
  }
  if( any(vapply(value, is.list, logical(1L))) ){
    stop("copy_nodes(): ", argument,
         " does not support lists nested more than one level.", call. = FALSE)
  }
  positions <- lapply(value, function(position){
    if( !is.character(position) || length(position) == 0L ||
        anyNA(position) || any(!nzchar(position)) ){
      stop("copy_nodes(): every ", argument,
           " list element must contain non-missing, non-empty character suffixes.",
           call. = FALSE)
    }
    value <- as.character(position)

    return(value)
  })
  flattened <- unlist(positions, use.names = FALSE)
  if( anyDuplicated(flattened) ){
    stop("copy_nodes(): ", argument, " suffixes must be unique.",
         call. = FALSE)
  }

  return(positions)
}


#' Check the suffix for an original node
#'
#' @param suffix Optional suffix for original nodes.
#' @return A list containing the suffix `value` and a `parallel` flag.
#' @noRd
copy_single_suffix <- function(suffix){
  if( is.null(suffix) ){
    suffix_spec <- list(value = NULL, parallel = FALSE)
    return(suffix_spec)
  }
  parallel <- is.list(suffix)
  if( parallel ){
    if( length(suffix) != 1L || is.list(suffix[[1L]]) ){
      stop("copy_nodes(): suffix must contain one character suffix.",
           call. = FALSE)
    }
    suffix <- suffix[[1L]]
  }
  if( !is.character(suffix) || length(suffix) != 1L || is.na(suffix) ||
      !nzchar(suffix) ){
    stop("copy_nodes(): suffix must contain one non-missing, non-empty character suffix.",
         call. = FALSE)
  }
  suffix_spec <- list(value = suffix, parallel = parallel)

  return(suffix_spec)
}


#' Construct original and copy node names
#'
#' @param nodes Character vector of node names.
#' @param suffix Optional suffix for original nodes.
#' @param copy_suffix Character suffixes or a list grouping parallel copies.
#' @param repeats Optional positive whole number of copies per original node.
#' @return A list containing `original_names`, `copy_names`, `copy_groups`, and
#'   `source_for_copy`.
#' @noRd
copy_name_specification <- function(nodes, suffix, copy_suffix, repeats){
  repeats <- copy_validate_repeats(repeats)
  original_suffix <- copy_single_suffix(suffix)

  if( !is.null(copy_suffix) ){
    positions <- copy_suffix_positions(copy_suffix, "copy_suffix")
    suffixes <- unlist(positions, use.names = FALSE)
    count <- length(suffixes)
    if( !is.null(repeats) ){
      if( repeats != count ){
        stop("copy_nodes(): repeats is ", repeats,
             ", but copy_suffix supplies ", count,
             " suffixes. Remove repeats or make the values agree.",
             call. = FALSE)
      }
      warning("copy_nodes(): repeats is redundant when copy_suffix supplies one suffix per copy.",
              call. = FALSE)
    }
    original_names <- if( is.null(original_suffix$value) ){
      nodes
    }else{
      paste0(nodes, "_", original_suffix$value)
    }
    groups <- positions
  }else if( !is.null(original_suffix$value) ){
    count <- if( is.null(repeats) ) 1L else repeats
    original_names <- paste0(nodes, "_", original_suffix$value, "1")
    suffixes <- paste0(original_suffix$value, seq.int(2L, count + 1L))
    groups <- if( original_suffix$parallel ) list(suffixes) else as.list(suffixes)
  }else{
    count <- if( is.null(repeats) ) 1L else repeats
    original_names <- nodes
    suffixes <- if( count == 1L ) "copy" else paste0("copy", seq_len(count))
    groups <- as.list(suffixes)
  }

  copy_names <- lapply(nodes, function(node){
    values <- lapply(groups, function(group) {
      text <- paste0(node, "_", group)

      return(text)
    })

    return(values)
  })
  names(copy_names) <- nodes
  all_copies <- unname(unlist(copy_names, recursive = TRUE,
                              use.names = FALSE))
  source_for_copy <- rep(nodes, vapply(copy_names, function(groups){
    value <- length(unlist(groups, recursive = TRUE, use.names = FALSE))

    return(value)
  }, integer(1L)))
  names(source_for_copy) <- all_copies

  name_spec <- list(
    original_names = original_names,
    copy_names = copy_names,
    copy_groups = groups,
    source_for_copy = source_for_copy
  )

  return(name_spec)
}


## Copied edges and graph updates ##
#' Connect each source node to each target node
#'
#' @param from Source node names.
#' @param to Target node names.
#' @return A `v`, `e`, `w` data frame of directed edges.
#' @noRd
copy_edges_between <- function(from, to){
  if( length(from) == 0L || length(to) == 0L ) {
    edges <- empty_edge_frame()
    return(edges)
  }
  pairs <- expand.grid(v = as.character(from), w = as.character(to),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  edges <- data.frame(v = pairs$v, e = "->", w = pairs$w,
                      stringsAsFactors = FALSE)

  return(edges)
}


#' Transfer an anchor path to its copies
#'
#' @param directed Directed edges with columns `v`, `e`, and `w`.
#' @param pivot Original name of the anchor node.
#' @param copies Copy names at the outer end of the copied path.
#' @param placement `"before"` or `"after"` the copied node.
#' @param renamed Named character vector mapping original node names to renamed
#'   nodes.
#' @return A list containing `taken` edges for the copies and `removed` edges
#'   for the anchor.
#' @noRd
copy_path_edges <- function(directed, pivot, copies, placement, renamed){
  point <- copy_remap_names(pivot, renamed)
  if( placement == "before" ){
    parents <- copy_remap_names(directed$v[directed$w == pivot], renamed)
    edge_changes <- list(taken = copy_edges_between(parents, copies),
                         removed = copy_edges_between(parents, point))
    return(edge_changes)
  }
  children <- copy_remap_names(directed$w[directed$v == pivot], renamed)
  edge_changes <- list(taken = copy_edges_between(copies, children),
                       removed = copy_edges_between(point, children))

  return(edge_changes)
}


#' Construct inherited and connecting edges for copies
#'
#' @param original_edges Original graph edges with columns `v`, `e`, and `w`.
#' @param nodes Character vector of node names.
#' @param renamed Named character vector mapping original node names to renamed
#'   nodes.
#' @param copy_names Nested copy names grouped by source node and temporal
#'   position.
#' @param copy_groups Suffix groups defining successive or parallel copy
#'   positions.
#' @param anchor Optional existing node that anchors the copies.
#' @param placements Named placement values, `"before"` or `"after"`, for each
#'   source node.
#' @param type `NULL` to inherit edges or `"chain"` to splice copies into a
#'   path.
#' @return A list containing added `edges`, `removed` edges, and a `rerouted`
#'   edge table.
#' @noRd
copy_generate_edges <- function(original_edges, nodes, renamed, copy_names,
                                 copy_groups, anchor, placements, type){
  directed <- original_edges[original_edges$e == "->", , drop = FALSE]
  inherited <- if( is.null(type) ){
    directed
  }else{
    empty_edge_frame()
  }
  selected_v <- inherited$v %in% nodes
  selected_w <- inherited$w %in% nodes
  pieces <- list()
  removals <- list()
  rerouted <- list()
  add_piece <- function(piece){
    if( nrow(piece) > 0L ) {
      value <- (pieces[[length(pieces) + 1L]] <<- piece)

      return(invisible(value))
    }else{
      value <- NULL
      return(invisible(value))
    }
  }
  add_removal <- function(piece){
    if( nrow(piece) > 0L ) {
      value <- (removals[[length(removals) + 1L]] <<- piece)

      return(invisible(value))
    }else{
      value <- NULL
      return(invisible(value))
    }
  }
  add_reroute <- function(from, to, replacement_from, replacement_to){
    if( !length(from) || !length(to) ) return()
    size <- max(length(from), length(to), length(replacement_from),
                length(replacement_to))
    value <- (rerouted[[length(rerouted) + 1L]] <<- data.frame(
      old_v = rep(from, length.out = size),
      old_w = rep(to, length.out = size),
      new_v = rep(replacement_from, length.out = size),
      new_w = rep(replacement_to, length.out = size),
      stringsAsFactors = FALSE
    ))

    return(invisible(value))
  }

  anchor_ancestors <- anchor_descendants <- character()
  if( !is.null(anchor) ){
    anchor_ancestors <- directed_reachable(original_edges, anchor,
                                           direction = "to")
    anchor_descendants <- directed_reachable(original_edges, anchor,
                                             direction = "from")
  }

  external_in <- inherited[!selected_v & selected_w, , drop = FALSE]
  for( i in seq_len(nrow(external_in)) ){
    source <- external_in$w[i]
    groups <- copy_names[[source]]
    copies <- unlist(groups, recursive = TRUE, use.names = FALSE)
    parent <- external_in$v[i]
    reroute <- !is.null(anchor) && placements[[source]] == "before" &&
      parent %in% c(anchor, anchor_descendants)
    if( reroute ){
      joining <- unlist(groups[[length(groups)]], recursive = TRUE,
                        use.names = FALSE)
      add_reroute(rep(parent, length(joining)), joining,
                  joining, rep(anchor, length(joining)))
    }else{
      add_piece(copy_edges_between(parent, copies))
    }
  }

  external_out <- inherited[selected_v & !selected_w, , drop = FALSE]
  for( i in seq_len(nrow(external_out)) ){
    source <- external_out$v[i]
    groups <- copy_names[[source]]
    copies <- unlist(groups, recursive = TRUE, use.names = FALSE)
    child <- external_out$w[i]
    anchor_after_source <- !is.null(anchor) &&
      placements[[source]] == "after" &&
      anchor %in% directed_reachable(original_edges, source,
                                     direction = "from")
    reroute <- anchor_after_source && child %in% c(anchor, anchor_ancestors)
    if( reroute ){
      joining <- unlist(groups[[1L]], recursive = TRUE, use.names = FALSE)
      add_reroute(joining, rep(child, length(joining)),
                  rep(anchor, length(joining)), joining)
    }else{
      add_piece(copy_edges_between(copies, child))
    }
  }

  internal <- inherited[selected_v & selected_w, , drop = FALSE]
  for( i in seq_len(nrow(internal)) ){
    from_groups <- copy_names[[internal$v[i]]]
    to_groups <- copy_names[[internal$w[i]]]
    for( group in seq_along(copy_groups) ){
      add_piece(data.frame(
        v = from_groups[[group]], e = "->", w = to_groups[[group]],
        stringsAsFactors = FALSE
      ))
    }
  }

  ## symmetric edges are inherited as they stand; an anchor never reroutes them ##
  symmetric <- if( is.null(type) ){
    original_edges[original_edges$e %in% c("<->", "--"), , drop = FALSE]
  }else{
    empty_edge_frame()
  }
  for( i in seq_len(nrow(symmetric)) ){
    lhs <- symmetric$v[i]
    rhs <- symmetric$w[i]
    operator <- symmetric$e[i]
    if( lhs %in% nodes && rhs %in% nodes ){
      from_groups <- copy_names[[lhs]]
      to_groups <- copy_names[[rhs]]
      for( group in seq_along(copy_groups) ){
        add_piece(data.frame(
          v = from_groups[[group]], e = operator, w = to_groups[[group]],
          stringsAsFactors = FALSE
        ))
      }
    }else if( lhs %in% nodes ){
      copies <- unlist(copy_names[[lhs]], recursive = TRUE, use.names = FALSE)
      add_piece(data.frame(v = copies, e = operator, w = rhs,
                           stringsAsFactors = FALSE))
    }else if( rhs %in% nodes ){
      copies <- unlist(copy_names[[rhs]], recursive = TRUE, use.names = FALSE)
      add_piece(data.frame(v = lhs, e = operator, w = copies,
                           stringsAsFactors = FALSE))
    }
  }

  for( source in nodes ){
    groups <- copy_names[[source]]
    if( placements[[source]] == "before" ){
      if( length(groups) > 1L ){
        for( group in seq_len(length(groups) - 1L) ){
          add_piece(copy_edges_between(groups[[group]], groups[[group + 1L]]))
        }
      }
      add_piece(copy_edges_between(groups[[length(groups)]], renamed[[source]]))
    }else{
      add_piece(copy_edges_between(renamed[[source]], groups[[1L]]))
      if( length(groups) > 1L ){
        for( group in seq_len(length(groups) - 1L) ){
          add_piece(copy_edges_between(groups[[group]], groups[[group + 1L]]))
        }
      }
    }
    if( !is.null(anchor) ){
      if( placements[[source]] == "before" ){
        add_piece(copy_edges_between(groups[[length(groups)]], anchor))
      }else{
        add_piece(copy_edges_between(anchor, groups[[1L]]))
      }
    }
    if( identical(type, "chain") ){
      pivot <- if( is.null(anchor) ) source else anchor
      end <- if( placements[[source]] == "before" ){
        groups[[1L]]
      }else{
        groups[[length(groups)]]
      }
      path <- copy_path_edges(directed, pivot, end, placements[[source]],
                               renamed)
      add_piece(path$taken)
      add_removal(path$removed)
      if( placements[[source]] == "before" ){
        add_reroute(path$taken$v, path$removed$w[1L], path$taken$v,
                    path$taken$w)
      }else{
        add_reroute(path$removed$v[1L], path$taken$w, path$taken$v,
                    path$taken$w)
      }
    }
  }

  edges <- if( length(pieces) ){
    unique(do.call(rbind, pieces))
  }else{
    empty_edge_frame()
  }
  rownames(edges) <- NULL
  removed <- if( length(removals) ){
    unique(do.call(rbind, removals))
  }else{
    empty_edge_frame()
  }
  rownames(removed) <- NULL
  rerouted <- if( length(rerouted) ){
    unique(do.call(rbind, rerouted))
  }else{
    data.frame(old_v = character(), old_w = character(),
               new_v = character(), new_w = character(),
               stringsAsFactors = FALSE)
  }
  rownames(rerouted) <- NULL
  edge_changes <- list(edges = edges, removed = removed, rerouted = rerouted)

  return(edge_changes)
}


#' Append copied nodes and their edges to native graph text
#'
#' @param dag A dagitty object.
#' @param new_nodes Names of nodes to add.
#' @param new_edges Edge table to append, with columns `v`, `e`, and `w`.
#' @param source_for_copy Named character vector mapping each copy to its
#'   source node.
#' @param source_dag Original dagitty object supplying node statuses.
#' @param caller Name of the calling function for messages.
#' @return A dagitty object containing the new nodes, their inherited treatment,
#'   outcome and latent statuses, and new edges.
#' @noRd
copy_append_graph <- function(dag, new_nodes, new_edges, source_for_copy,
                               source_dag, caller){
  exposures <- treatments(source_dag)
  outcomes <- get_outcome_names(source_dag)
  latents <- unobserved(source_dag)
  declarations <- vapply(new_nodes, function(node){
    source <- unname(source_for_copy[[node]])
    status <- c(
      if( source %in% exposures ) "exposure",
      if( source %in% outcomes ) "outcome",
      if( source %in% latents ) "latent"
    )
    text <- paste0(native_layout_quote(node),
                   if( length(status) ) paste0(" [", paste(status, collapse = ","), "]") else "")

    return(text)
  }, character(1L), USE.NAMES = FALSE)
  edge_lines <- if( nrow(new_edges) ){
    paste(native_layout_quote(new_edges$v), new_edges$e,
          native_layout_quote(new_edges$w))
  }else{
    character()
  }

  text <- as.character(dag)
  closing <- native_layout_index(text, caller)$closing
  addition <- paste0("\n", paste(c(declarations, edge_lines), collapse = "\n"),
                     "\n")
  result <- native_layout_value(
    dag,
    native_layout_edits(text, list(list(
      start = closing, end = closing - 1L, text = addition
    )))
  )
  expected_nodes <- c(names(dag), new_nodes)
  if( !setequal(names(result), expected_nodes) ){
    stop(caller, ": new node declarations could not be added safely.",
         call. = FALSE)
  }

  return(result)
}


#' Print the edges rerouted through node copies
#'
#' @param rerouted Table with `old_v`, `old_w`, `new_v`, and `new_w` endpoint
#'   columns.
#' @return Invisible `NULL`.
#' @noRd
copy_print_reroutes <- function(rerouted){
  if( nrow(rerouted) == 0L ) {
    out <- NULL
    return(invisible(out))
  }
  lines <- paste0("'", rerouted$old_v, " -> ", rerouted$old_w,
                  "' => '", rerouted$new_v, " -> ", rerouted$new_w, "'")
  message(paste(c("copy_nodes(): rerouted edges:", lines),
                collapse = "\n"))
  out <- NULL

  return(invisible(out))
}
