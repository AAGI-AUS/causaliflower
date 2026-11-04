# ==============================================================================
# graph_planner.R


# Edge planning and position handling for graph_funcs.r
# ==============================================================================


## Package constants ##
CF_INCOMING_MARKER <- "..causaliflower.incoming.." # Vector name used by incoming()


CF_OUTGOING_MARKER <- "..causaliflower.outgoing.." # Vector name used by outgoing()


CF_NON_IV_ROLES <- c( # Roles reconstructed from a graph, with instruments handled separately
  "confounder", "treatment", "outcome", "mediator",
  "mediator_outcome_confounder", "competing_cause", "collider",
  "latent", "observed"
)


CF_CHAIN_TYPES <- c("chain", "first", "last") # Types joining adjacent positions within a role


## Routing helpers ##
#' Choose where arrows enter or leave a role
#'
#' `incoming()` marks nodes that receive generated arrows from other roles.
#' `outgoing()` marks nodes that send generated arrows to other roles. The
#' names `incoming` and `outgoing` do the same inside `c()` or `list()`. Node
#' names are still supplied only once in the role input.
#'
#' @param nodes Character vector of node names.
#' @return A character vector for use in any role except instrumental variables with
#'   [build_graph()] or [add_nodes()].
#' @examples
#' routed <- build_graph(
#'   confounders = "X", treatments = "T", outcomes = "Y",
#'   mediators = c(incoming("M1"), "M2", outgoing("M3")),
#'   type = "saturated"
#' )
#' get_edges(routed)
#'
#' named <- build_graph(
#'   confounders = "X", treatments = "T", outcomes = "Y",
#'   mediators = c(incoming = "M1", "M2", outgoing = "M3"),
#'   type = "saturated"
#' )
#' get_edges(named)
#' @name routing_helpers
#' @export
incoming <- function(nodes){
  nodes_out <- routing_annotation(nodes, CF_INCOMING_MARKER, "incoming")

  return(nodes_out)
}


#' @rdname routing_helpers
#' @export
outgoing <- function(nodes){
  nodes_out <- routing_annotation(nodes, CF_OUTGOING_MARKER, "outgoing")

  return(nodes_out)
}


#' Mark nodes for incoming or outgoing edges
#'
#' @param nodes Atomic vector of node names.
#' @param marker The name used to mark the selected routing direction.
#' @param caller Function name to include in errors.
#' @return Character vector of node names carrying the routing marker as their
#'   names.
#' @noRd
routing_annotation <- function(nodes, marker, caller){
  if( is.list(nodes) ){
    stop(caller, "(): 'nodes' must be an atomic vector.", call. = FALSE)
  }
  nodes <- as.character(nodes)
  if( length(nodes) == 0L || anyNA(nodes) || any(!nzchar(nodes)) ){
    stop(caller, "(): supply non-missing, non-empty node names.",
         call. = FALSE)
  }
  nodes_out <- stats::setNames(nodes, rep(marker, length(nodes)))

  return(nodes_out)
}


#' Check a routing annotation name
#'
#' @param name A single annotation name.
#' @param side The routing direction: incoming or outgoing.
#' @return TRUE if the name matches the direction, with or without a numeric
#'   suffix.
#' @noRd
is_routing_name <- function(name, side){
  if( length(name) != 1L || is.na(name) ){
    matches <- FALSE
    return(matches)
  }
  matches <- grepl(paste0("^", side, "[0-9]*$"), name)

  return(matches)
}


#' Find routing annotations in vector names
#'
#' @param value_names Names attached to a role input vector.
#' @param marker The routing marker to recognise.
#' @param side The routing direction: incoming or outgoing.
#' @return Logical vector identifying names that match the marker or direction.
#' @noRd
routing_names <- function(value_names, marker, side){
  matches <- value_names == marker | grepl(paste0("^", side, "[0-9]*$"), value_names)

  return(matches)
}


#' Create an empty role edge table
#'
#' @return A zero-row data table with character columns v, e and w.
#' @noRd
empty_role_edges <- function(){
  edges <- data.table::data.table(v = character(), e = character(), w = character())

  return(edges)
}


#' Convert role input to node names
#'
#' @param x Node names to convert to character.
#' @return Character vector of node names, including character(0) for empty
#'   input.
#' @noRd
role_input_names <- function(x){
  if( length(x) == 0L ) {
    nodes <- character()

    return(nodes)
  } else {
    nodes <- as.character(x)

    return(nodes)
  }
}


#' Create edges between two sets of nodes
#'
#' @param from Source node names.
#' @param to Target node names.
#' @param edge The edge operator to assign.
#' @return Data table with one edge for each source-target pair and columns v, e
#'   and w.
#' @noRd
role_edge_frame <- function(from, to, edge = "->"){
  from <- as.character(from)
  to <- as.character(to)
  from <- from[!is.na(from) & nzchar(from)]
  to <- to[!is.na(to) & nzchar(to)]
  if( length(from) == 0L || length(to) == 0L ){
    edges <- empty_role_edges()
    return(edges)
  }
  pairs <- expand.grid(
    v = from, w = to, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  edges <- data.table::data.table(v = pairs$v, e = edge, w = pairs$w)

  return(edges)
}


#' Combine role edge tables
#'
#' @param ... Edge tables with columns v, e and w.
#' @return Data table of unique edges without missing endpoints or self-edges.
#' @noRd
bind_role_edges <- function(...){
  pieces <- Filter(function(x) {
    matches <- !is.null(x) && nrow(x) > 0L

    return(matches)
  }, list(...))
  if( length(pieces) == 0L ){
    edges_out <- empty_role_edges()
    return(edges_out)
  }
  edges <- data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
  edges <- edges[complete.cases(edges), c("v", "e", "w"), with = FALSE]
  edges <- edges[edges$v != edges$w,]

  if( any(edges$e %in% c("<->", "--")) ){
    symmetric <- edges$e %in% c("<->", "--")
    left <- pmin(edges$v[symmetric], edges$w[symmetric])
    right <- pmax(edges$v[symmetric], edges$w[symmetric])
    edges$v[symmetric] <- left
    edges$w[symmetric] <- right
  }
  edges_out <- unique(edges)

  return(edges_out)
}


#' Check for an absent role input
#'
#' @param x Role input supplied as a vector or list.
#' @return TRUE if the input is empty or contains only missing values.
#' @noRd
absent_role_input <- function(x){
  if( is.null(x) || length(x) == 0L ){
    absent <- TRUE
    return(absent)
  }
  flat <- unlist(x, recursive = TRUE, use.names = FALSE)
  absent <- length(flat) == 0L || all(is.na(flat))

  return(absent)
}


#' Check the shape of a latent mapping
#'
#' @param x A role input entry to inspect.
#' @return TRUE for a two-element list with one source and non-empty targets.
#' @noRd
is_latent_mapping_entry <- function(x){
  mapping <- is.list(x) && length(x) == 2L &&
    !is.list(x[[1L]]) &&
    length(unlist(x[[1L]], recursive = TRUE, use.names = FALSE)) == 1L &&
    length(unlist(x[[2L]], recursive = TRUE, use.names = FALSE)) > 0L

  return(mapping)
}


#' Read node names and routing annotations
#'
#' @param x An atomic vector of node names for one temporal position.
#' @param top_name The name of the containing list entry, or NULL.
#' @param caller Function name to include in errors.
#' @return A list of node names, incoming and outgoing selections, and flags
#'   indicating which sides were marked.
#' @noRd
role_leaf <- function(x, top_name = NULL, caller = "role input"){
  if( is.list(x) ){
    stop(caller, ": nested lists are not temporal positions.", call. = FALSE)
  }
  values <- as.character(x)
  value_names <- names(x)
  if( is.null(value_names) ){
    value_names <- rep("", length(values))
  }

  numbered_both <- grepl("^both[0-9]+$", value_names)
  if( any(numbered_both) ){
    stop(caller, ": 'both' is a scalar annotation; mark one node at a time.",
         call. = FALSE)
  }
  both <- value_names == "both"
  if( identical(top_name, "both") ){
    if( length(values) != 1L ){
      stop(caller, ": 'both' is a scalar annotation; mark one node at a time.",
           call. = FALSE)
    }
    both[] <- TRUE
  }
  if( sum(both) > 1L ){
    stop(caller, ": 'both' is a scalar annotation; mark one node at a time.",
         call. = FALSE)
  }

  incoming <- routing_names(value_names, CF_INCOMING_MARKER, "incoming")
  outgoing <- routing_names(value_names, CF_OUTGOING_MARKER, "outgoing")
  if( is_routing_name(top_name, "incoming") ){
    incoming[] <- TRUE
  }
  if( is_routing_name(top_name, "outgoing") ){
    outgoing[] <- TRUE
  }
  incoming <- incoming | both
  outgoing <- outgoing | both

  spec <- list(
    nodes = unname(values),
    incoming = unname(values[incoming]),
    outgoing = unname(values[outgoing]),
    incoming_marked = any(incoming),
    outgoing_marked = any(outgoing)
  )

  return(spec)
}


#' Create an empty role specification
#'
#' @param role The role name.
#' @return A list with empty node positions, routing selections and latent
#'   mappings.
#' @noRd
empty_role_spec <- function(role){
  spec <- list(
    role = role,
    positions = list(),
    nodes = character(),
    incoming = character(),
    outgoing = character(),
    incoming_marked = FALSE,
    outgoing_marked = FALSE,
    mapping = FALSE,
    mapping_sources = character(),
    mapping_edges = empty_role_edges(),
    same_role = !identical(role, "observed")
  )

  return(spec)
}


#' Read positions and routing from a role input
#'
#' @param x Role input supplied as a vector or list.
#' @param role The role name.
#' @param caller Function name to include in errors.
#' @return A role specification containing node positions, routing selections
#'   and any latent mappings.
#' @noRd
normalise_role_input <- function(x, role, caller = "role input"){
  spec <- empty_role_spec(role)
  if( absent_role_input(x) ){
    return(spec)
  }

  if( is.list(x) ){
    top_names <- names(x)
    if( is.null(top_names) ){
      top_names <- rep("", length(x))
    }
    leaves <- lapply(seq_along(x), function(i){
      entry <- x[[i]]
      if( identical(role, "latent") && is_latent_mapping_entry(entry) ){
        source <- role_leaf(entry[[1L]], caller = caller)
        targets <- as.character(unlist(entry[[2L]], recursive = TRUE,
                                       use.names = FALSE))
        if( length(targets) == 0L || anyNA(targets) || any(!nzchar(targets)) ){
          stop(caller, ": latent mapping targets must be non-missing node names.",
               call. = FALSE)
        }
        source$mapping_targets <- targets
        return(source)
      }
      if( is.list(entry) ){
        if( identical(role, "latent") ){
          stop(caller, ": nested latent entries must use one source and its targets.",
               call. = FALSE)
        }
        stop(caller, ": nested lists are not supported for this role.",
             call. = FALSE)
      }
      leaf <- role_leaf(entry, top_name = top_names[[i]], caller = caller)

      return(leaf)
    })
  }else{
    leaf <- role_leaf(x, caller = caller)
    leaves <- lapply(seq_along(leaf$nodes), function(i){
      position <- list(
        nodes = leaf$nodes[[i]],
        incoming = intersect(leaf$nodes[[i]], leaf$incoming),
        outgoing = intersect(leaf$nodes[[i]], leaf$outgoing),
        incoming_marked = leaf$nodes[[i]] %in% leaf$incoming,
        outgoing_marked = leaf$nodes[[i]] %in% leaf$outgoing
      )

      return(position)
    })
  }

  positions <- lapply(leaves, function(leaf) {
    position_nodes <- leaf$nodes

    return(position_nodes)
  })
  positions <- Filter(length, positions)
  nodes <- role_input_names(unlist(positions, use.names = FALSE))
  if( length(nodes) == 0L || anyNA(nodes) || any(!nzchar(nodes)) ){
    stop(caller, ": node names must be non-missing and non-empty.",
         call. = FALSE)
  }
  duplicates <- unique(nodes[duplicated(nodes)])
  if( length(duplicates) > 0L ){
    stop(caller, ": each node may occupy only one temporal position: ",
         paste(duplicates, collapse = ", "), ".", call. = FALSE)
  }

  spec$positions <- positions
  spec$nodes <- nodes
  spec$incoming <- unique(unlist(lapply(leaves, function(leaf){
    incoming_nodes <- leaf$incoming

    return(incoming_nodes)
  }), use.names = FALSE))
  spec$outgoing <- unique(unlist(lapply(leaves, function(leaf){
    outgoing_nodes <- leaf$outgoing

    return(outgoing_nodes)
  }), use.names = FALSE))
  spec$incoming_marked <- any(vapply(leaves, function(leaf){
    matches <- isTRUE(leaf$incoming_marked)

    return(matches)
  }, logical(1L)))
  spec$outgoing_marked <- any(vapply(leaves, function(leaf){
    matches <- isTRUE(leaf$outgoing_marked)

    return(matches)
  }, logical(1L)))
  mapped <- vapply(leaves, function(leaf){
    matches <- !is.null(leaf$mapping_targets)

    return(matches)
  }, logical(1L))
  spec$mapping_sources <- unique(role_input_names(unlist(lapply(leaves[mapped], function(leaf){
    source_nodes <- leaf$nodes

    return(source_nodes)
  }), use.names = FALSE)))
  spec$mapping_edges <- do.call(bind_role_edges, lapply(leaves[mapped], function(leaf){
    mapping_edges <- role_edge_frame(leaf$nodes, leaf$mapping_targets, "->")

    return(mapping_edges)
  }))
  spec$mapping <- length(spec$mapping_sources) > 0L

  return(spec)
}


#' Keep selected nodes in a role specification
#'
#' @param spec A role specification.
#' @param nodes Node names to retain.
#' @param preserve_routing Whether to retain the original routing selections and
#'   flags.
#' @return The role specification restricted to the selected nodes.
#' @noRd
subset_role_spec <- function(spec, nodes, preserve_routing = FALSE){
  if( is.null(spec$mapping_sources) ){
    spec$mapping_sources <- unique(as.character(spec$mapping_edges$v))
  }
  routing <- list(
    incoming = spec$incoming,
    outgoing = spec$outgoing,
    incoming_marked = spec$incoming_marked,
    outgoing_marked = spec$outgoing_marked
  )
  nodes <- intersect(spec$nodes, nodes)
  spec$positions <- Filter(length, lapply(spec$positions, intersect, y = nodes))
  spec$nodes <- role_input_names(unlist(spec$positions, use.names = FALSE))
  spec$incoming <- role_input_names(intersect(
    spec[["incoming", exact = TRUE]], spec$nodes
  ))
  spec$outgoing <- role_input_names(intersect(
    spec[["outgoing", exact = TRUE]], spec$nodes
  ))
  spec$incoming_marked <- length(spec$incoming) > 0L
  spec$outgoing_marked <- length(spec$outgoing) > 0L
  spec$mapping_sources <- role_input_names(intersect(
    spec[["mapping_sources", exact = TRUE]], spec$nodes
  ))
  spec$mapping_edges <- spec$mapping_edges[
    spec$mapping_edges$v %in% spec$mapping_sources,
  ]
  spec$mapping <- length(spec$mapping_sources) > 0L
  spec$same_role <- !identical(spec$role, "observed")
  if( isTRUE(preserve_routing) ){
    spec$incoming <- routing$incoming
    spec$outgoing <- routing$outgoing
    spec$incoming_marked <- routing$incoming_marked
    spec$outgoing_marked <- routing$outgoing_marked
  }

  return(spec)
}


## Edge planning ##
#' Choose role nodes at a chain boundary
#'
#' @param spec A role specification.
#' @param side The endpoint direction: incoming or outgoing.
#' @param source_position Which outgoing source position to use: first or last.
#' @return Character vector of eligible endpoint names.
#' @noRd
chain_role_side <- function(spec, side, source_position){
  positions <- Filter(length, lapply(
    spec$positions, setdiff, y = spec$mapping_sources
  ))
  if( length(positions) == 0L ){
    nodes <- character()
    return(nodes)
  }

  direct_source_roles <- c(
    "mediator_outcome_confounder", "competing_cause"
  )
  if( identical(side, "outgoing") &&
      spec$role %in% direct_source_roles ){
    nodes <- unique(role_input_names(unlist(positions, use.names = FALSE)))
    return(nodes)
  }
  if( identical(side, "incoming") &&
      identical(spec$role, "collider") ){
    nodes <- unique(role_input_names(unlist(positions, use.names = FALSE)))
    return(nodes)
  }

  position <- if( identical(side, "incoming") ){
    positions[[1L]]
  }else if( identical(source_position, "first") ){
    positions[[1L]]
  }else{
    positions[[length(positions)]]
  }
  nodes <- unique(role_input_names(position))

  return(nodes)
}


#' Choose endpoints for edges between roles
#'
#' @param spec A role specification.
#' @param side The endpoint direction: incoming or outgoing.
#' @param type The graph construction type.
#' @return Character vector of endpoint names, excluding explicit latent mapping
#'   sources.
#' @noRd
role_side <- function(spec, side, type){
  marked <- spec[[paste0(side, "_marked")]]
  nodes <- if( isTRUE(marked) ){
    spec[[side]]
  }else if( type %in% c("first", "last") ){
    chain_role_side(spec, side, source_position = type)
  }else{
    spec$nodes
  }
  nodes_out <- setdiff(nodes, spec$mapping_sources)

  return(nodes_out)
}


#' Create edges within a role
#'
#' @param spec A role specification.
#' @param type The graph construction type.
#' @param edge The edge operator to assign.
#' @return Data table of within-role edges with columns v, e and w.
#' @noRd
position_edges <- function(spec, type, edge){
  if( !isTRUE(spec$same_role) ){
    edges <- empty_role_edges()
    return(edges)
  }

  positions <- Filter(length, lapply(
    spec$positions, setdiff, y = spec$mapping_sources
  ))
  nodes <- unlist(positions, use.names = FALSE)
  if( length(nodes) < 2L ){
    edges <- empty_role_edges()
    return(edges)
  }

  if( identical(type, "full") ){
    full_roles <- c("confounder", "treatment", "outcome", "mediator", "latent")
    if( !spec$role %in% full_roles ){
      edges <- empty_role_edges()
      return(edges)
    }
    pairs <- utils::combn(nodes, 2L)
    edges <- data.table::data.table(
      v = pairs[1L,], e = edge, w = pairs[2L,]
    )
    return(edges)
  }

  if( length(positions) < 2L ){
    edges <- empty_role_edges()
    return(edges)
  }
  pairs <- if( type %in% CF_CHAIN_TYPES ){
    cbind(seq_len(length(positions) - 1L), seq_len(length(positions) - 1L) + 1L)
  }else{
    utils::combn(seq_along(positions), 2L)
  }
  if( type %in% CF_CHAIN_TYPES ){
    pieces <- lapply(seq_len(nrow(pairs)), function(i){
      edges <- role_edge_frame(positions[[pairs[i, 1L]]], positions[[pairs[i, 2L]]], edge)

      return(edges)
    })
  }else{
    pieces <- lapply(seq_len(ncol(pairs)), function(i){
      edges <- role_edge_frame(positions[[pairs[1L, i]]], positions[[pairs[2L, i]]], edge)

      return(edges)
    })
  }
  edges <- do.call(bind_role_edges, pieces)

  return(edges)
}


#' List the role pairs connected by a graph type
#'
#' @param type The graph construction type.
#' @return Two-column character matrix of source and target roles, one pair per
#'   row.
#' @noRd
cross_role_rules <- function(type){
  if( identical(type, "chain") ){
    role_pairs <- matrix(character(), nrow = 0L, ncol = 2L)
    return(role_pairs)
  }

  if( identical(type, "full") ){
    role_pairs <- rbind(
      c("confounder", "treatment"),
      c("confounder", "outcome"),
      c("confounder", "mediator"),
      c("confounder", "mediator_outcome_confounder"),
      c("confounder", "latent"),
      c("treatment", "outcome"),
      c("treatment", "mediator"),
      c("treatment", "collider"),
      c("treatment", "mediator_outcome_confounder"),
      c("latent", "treatment"),
      c("outcome", "collider"),
      c("mediator", "outcome"),
      c("mediator", "latent"),
      c("mediator_outcome_confounder", "outcome"),
      c("mediator_outcome_confounder", "mediator"),
      c("mediator_outcome_confounder", "latent"),
      c("competing_cause", "outcome"),
      c("latent", "outcome")
    )
    return(role_pairs)
  }

  if( identical(type, "saturated") ){
    role_pairs <- rbind(
      c("confounder", "treatment"),
      c("confounder", "mediator"),
      c("confounder", "outcome"),
      c("confounder", "collider"),
      c("treatment", "mediator"),
      c("treatment", "outcome"),
      c("treatment", "collider"),
      c("mediator_outcome_confounder", "mediator"),
      c("mediator_outcome_confounder", "outcome"),
      c("mediator_outcome_confounder", "collider"),
      c("mediator", "outcome"),
      c("mediator", "collider"),
      c("competing_cause", "outcome"),
      c("competing_cause", "collider"),
      c("outcome", "collider"),
      c("latent", "outcome"),
      c("latent", "treatment"),
      c("latent", "confounder")
    )
    return(role_pairs)
  }

  role_pairs <- rbind(
    c("confounder", "treatment"),
    c("confounder", "outcome"),
    c("treatment", "outcome"),
    c("treatment", "mediator"),
    c("treatment", "collider"),
    c("outcome", "collider"),
    c("mediator", "outcome"),
    c("mediator_outcome_confounder", "outcome"),
    c("mediator_outcome_confounder", "mediator"),
    c("competing_cause", "outcome"),
    c("latent", "outcome"),
    c("latent", "treatment")
  )

  return(role_pairs)
}


#' Plan graph edges from role specifications
#'
#' @param specs Named list of role specifications.
#' @param type The graph construction type.
#' @param instrumental_variables Instrument names or explicit instrument
#'   mappings.
#' @param existing_dag An existing dagitty object for connecting observed nodes,
#'   or NA.
#' @param include_iv Whether to include instrument edges.
#' @return Data table of planned edges with columns v, e and w.
#' @noRd
plan_graph_edges <- function(specs, type, instrumental_variables = NA,
                                 existing_dag = NA, include_iv = TRUE){
  edge <- if( identical(type, "full") ) { "<->" } else { "->" }
  same_role <- lapply(specs, position_edges, type = type, edge = edge)

  rules <- cross_role_rules(type)
  cross_role <- lapply(seq_len(nrow(rules)), function(i){
    source_role <- rules[i, 1L]
    target_role <- rules[i, 2L]
    source <- specs[[source_role]]
    target <- specs[[target_role]]
    edges <- role_edge_frame(
      role_side(source, "outgoing", type),
      role_side(target, "incoming", type),
      edge
    )

    return(edges)
  })

  latent_mapping <- specs$latent$mapping_edges
  if( nrow(latent_mapping) > 0L ){
    latent_mapping$e <- edge
  }

  iv_edges <- empty_role_edges()
  if( isTRUE(include_iv) ){
    iv_edges <- draw_iv_edges(
      type = type,
      instrumental_variables = instrumental_variables,
      treatments = specs$treatment$nodes,
      e = edge
    )
    if( nrow(iv_edges) > 0L ){
      names(iv_edges) <- c("v", "e", "w")
    }else{
      iv_edges <- empty_role_edges()
    }
  }

  observed_edges <- empty_role_edges()
  if( dagitty::is.dagitty(existing_dag) && length(specs$observed$nodes) > 0L ){
    observed_edges <- draw_observed_edges(
      observed = specs$observed$nodes,
      existing_dag = existing_dag,
      e = edge
    )
    if( nrow(observed_edges) > 0L ){
      names(observed_edges) <- c("v", "e", "w")
    }
  }

  edges <- do.call(bind_role_edges, c(same_role, cross_role,
                                     list(latent_mapping, iv_edges, observed_edges)))

  return(edges)
}


#' Collect the graph role specifications
#'
#' @param confounders Confounder node names or positions.
#' @param treatments Treatment node names or positions.
#' @param outcomes Outcome node names or positions.
#' @param mediators Mediator node names or positions.
#' @param mediator_outcome_confounders Mediator-outcome confounder names or
#'   positions.
#' @param competing_causes Competing cause names or positions.
#' @param colliders Collider node names or positions.
#' @param latent_variables Latent node names, positions or source-target
#'   mappings.
#' @param observed Names or positions of nodes without an assigned role.
#' @param caller Function name to include in errors.
#' @return Named list of role specifications for graph construction.
#' @noRd
role_specs <- function(confounders, treatments, outcomes, mediators,
                           mediator_outcome_confounders, competing_causes,
                           colliders, latent_variables, observed = NA,
                           caller = "build_graph()"){
  inputs <- list(
    confounder = confounders,
    treatment = treatments,
    outcome = outcomes,
    mediator = mediators,
    mediator_outcome_confounder = mediator_outcome_confounders,
    competing_cause = competing_causes,
    collider = colliders,
    latent = latent_variables,
    observed = observed
  )
  specs <- Map(function(input, role){
    spec <- normalise_role_input(
      input, role,
      caller = paste0(caller, " ", role, " input")
    )

    return(spec)
  }, inputs, names(inputs))
  names(specs) <- names(inputs)

  return(specs)
}


## Roles and positions in existing graphs ##
#' Infer temporal positions from directed edges
#'
#' @param dag A dagitty object.
#' @param nodes Node names belonging to one role.
#' @return List of node vectors grouped by directed rank within the role.
#' @noRd
role_positions <- function(dag, nodes){
  nodes <- unique(as.character(nodes))
  nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
  if( length(nodes) == 0L ){
    positions <- list()
    return(positions)
  }

  edges <- dag_edges(dag)
  directed <- edges[edges$e == "->",]
  from <- as.character(directed$v)
  to <- as.character(directed$w)
  within <- from %in% nodes & to %in% nodes
  from <- from[within]
  to <- to[within]

  back <- feedback_arcs(nodes, from, to)
  if( length(back) > 0L ){
    from <- from[-back]
    to <- to[-back]
  }

  parents <- lapply(nodes, function(node) {
    parent_nodes <- unique(from[to == node])

    return(parent_nodes)
  })
  names(parents) <- nodes
  rank <- topological_rank(parents)
  positions <- lapply(sort(unique(rank)), function(step) {
    position_nodes <- nodes[rank == step]

    return(position_nodes)
  })

  return(positions)
}


#' Read outgoing latent mappings
#'
#' @param dag A dagitty object.
#' @param nodes Latent node names whose outgoing edges are read.
#' @return Data table of directed edges from the selected nodes to nodes outside
#'   that set.
#' @noRd
latent_mapping_edges <- function(dag, nodes){
  if( length(nodes) == 0L ){
    edges_out <- empty_role_edges()
    return(edges_out)
  }
  edges <- dag_edges(dag)
  mapped <- edges[edges$e == "->" & edges$v %in% nodes & !edges$w %in% nodes,]
  if( nrow(mapped) == 0L ){
    edges_out <- empty_role_edges()
    return(edges_out)
  }
  edges_out <- data.table::data.table(
    v = as.character(mapped$v), e = "->", w = as.character(mapped$w)
  )

  return(edges_out)
}


#' Infer role specifications from a graph
#'
#' @param dag A dagitty object.
#' @return Named list of role specifications containing node positions and
#'   latent mappings.
#' @noRd
existing_role_specs <- function(dag){
  roles <- classify_node_roles(dag)
  node_names <- names(dag)

  specs <- lapply(CF_NON_IV_ROLES, function(role){
    nodes <- roles[[role]]
    if( identical(role, "observed") ){
      nodes <- unique(c(nodes, roles$undetermined))
    }
    nodes <- as.character(nodes)
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    nodes <- node_names[node_names %in% nodes]

    spec <- empty_role_spec(role)
    spec$positions <- role_positions(dag, nodes)
    spec$nodes <- role_input_names(unlist(spec$positions, use.names = FALSE))
    if( identical(role, "observed") ){
      spec$same_role <- FALSE
    }
    if( identical(role, "latent") ){
      spec$mapping_edges <- latent_mapping_edges(dag, spec$nodes)
      spec$mapping_sources <- unique(as.character(spec$mapping_edges$v))
      spec$mapping <- length(spec$mapping_sources) > 0L
    }

    return(spec)
  })
  names(specs) <- CF_NON_IV_ROLES

  return(specs)
}


#' Find nodes sharing selected temporal positions
#'
#' @param spec A role specification.
#' @param nodes Node names whose temporal positions are selected.
#' @return Character vector of all nodes at the selected positions.
#' @noRd
tier_members <- function(spec, nodes){
  nodes_out <- unique(role_input_names(unlist(Filter(function(position){
    matches <- any(position %in% nodes)

    return(matches)
  }, spec$positions), use.names = FALSE)))

  return(nodes_out)
}


#' Separate positions to move from positions to retain
#'
#' @param existing The existing role specification.
#' @param supplied The supplied role specification.
#' @param strict Whether selecting part of an existing position is an error.
#' @return A list containing blocks to place and the remaining positions as
#'   spine.
#' @noRd
position_blocks <- function(existing, supplied, strict = TRUE){
  selected_existing <- intersect(supplied$nodes, existing$nodes)
  moving_indices <- which(vapply(existing$positions, function(position){
    matches <- any(position %in% selected_existing)

    return(matches)
  }, logical(1L)))

  if( isTRUE(strict) && length(moving_indices) > 0L ){
    split <- vapply(existing$positions[moving_indices], function(position){
      matches <- !all(position %in% selected_existing)

      return(matches)
    }, logical(1L))
    if( any(split) ){
      stop("position: nodes supplied together at one temporal position must be moved together.",
           call. = FALSE)
    }
  }

  existing_blocks <- existing$positions[moving_indices]
  new_blocks <- Filter(length, lapply(supplied$positions, function(position){
    new_nodes <- setdiff(position, existing$nodes)

    return(new_nodes)
  }))
  position_blocks <- list(
    blocks = c(existing_blocks, new_blocks),
    spine = if( length(moving_indices) ){
      existing$positions[-moving_indices]
    }else{
      existing$positions
    }
  )

  return(position_blocks)
}


#' Read position helpers and their name/value form
#'
#' @param position A cf_position object or a character position specification.
#' @return List of placement clauses, each containing kind, anchor and nodes.
#' @noRd
position_clauses <- function(position){
  if( inherits(position, "cf_position") ){
    clauses <- unclass(position)
    return(clauses)
  }
  if( is.list(position) || !is.atomic(position) ){
    stop("position: use first(), last(), before(), or after().", call. = FALSE)
  }
  values <- as.character(position)
  labels <- names(position)
  if( is.null(labels) ){
    labels <- rep("", length(values))
  }
  clauses <- lapply(seq_along(values), function(i){
    kind <- sub("[0-9]+$", "", labels[[i]])
    anchor <- NULL
    if( kind %in% c("before", "after") ){
      anchor <- values[[i]]
    }else if( !nzchar(kind) ){
      kind <- values[[i]]
    }
    if( !kind %in% c("first", "last") ){
      if( is.null(anchor) || is.na(anchor) || !nzchar(anchor) ){
        stop("position: use first(), last(), before(), or after().",
             call. = FALSE)
      }
    }
    values_out <- list(kind = kind, anchor = anchor, nodes = NULL)

    return(values_out)
  })

  return(clauses)
}


#' Place complete temporal positions
#'
#' @param position Position helpers or their name/value form.
#' @param blocks List of node groups to place.
#' @param spine List of existing node groups that remain in order.
#' @param supplied_nodes Node names available for selection by the position
#'   helpers.
#' @return List of node groups in the resolved temporal order.
#' @noRd
resolve_position_blocks <- function(position, blocks, spine,
                                        supplied_nodes){
  clauses <- position_clauses(position)
  required <- c("kind", "anchor", "nodes")
  if( !is.list(clauses) || any(!vapply(clauses, function(clause){
    matches <- is.list(clause) && all(required %in% names(clause))

    return(matches)
  }, logical(1L))) ){
    stop("position: use first(), last(), before(), or after().",
         call. = FALSE)
  }

  empty <- vapply(clauses, function(clause) {
    matches <- is.null(clause$nodes)

    return(matches)
  }, logical(1L))
  if( sum(empty) > 1L ){
    stop("position: only one position helper may omit its nodes argument.",
         call. = FALSE)
  }
  claimed <- unlist(lapply(clauses, `[[`, "nodes"), use.names = FALSE)
  duplicates <- unique(claimed[duplicated(claimed)])
  if( length(duplicates) > 0L ){
    stop("position: each selected node may appear in only one position helper: ",
         paste(duplicates, collapse = ", "), ".", call. = FALSE)
  }
  unknown <- setdiff(claimed, supplied_nodes)
  if( length(unknown) > 0L ){
    stop("position: these nodes are not among the selected nodes: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }

  block_clause <- rep(NA_integer_, length(blocks))
  for(i in seq_along(blocks)){
    hits <- which(vapply(clauses, function(clause){
      matches <- length(intersect(blocks[[i]], clause$nodes)) > 0L

      return(matches)
    }, logical(1L)))
    if( length(hits) > 1L ||
        (length(hits) == 1L && !all(blocks[[i]] %in% clauses[[hits]]$nodes)) ){
      stop("position: nodes supplied together at one temporal position must be moved together.",
           call. = FALSE)
    }
    if( length(hits) == 1L ){
      block_clause[[i]] <- hits
    }
  }
  catch_all <- which(empty)
  if( length(catch_all) == 1L ){
    block_clause[is.na(block_clause)] <- catch_all
  }else if( anyNA(block_clause) ){
    clauses[[length(clauses) + 1L]] <- list(
      kind = "last", anchor = NULL, nodes = NULL
    )
    block_clause[is.na(block_clause)] <- length(clauses)
  }

  clause_blocks <- lapply(seq_along(clauses), function(i){
    selected_blocks <- blocks[which(block_clause == i)]

    return(selected_blocks)
  })
  order <- spine
  first_blocks <- list()
  last_blocks <- list()
  after_counts <- list()
  for(i in seq_along(clauses)){
    insert <- clause_blocks[[i]]
    if( length(insert) == 0L ){
      next
    }
    clause <- clauses[[i]]
    if( identical(clause$kind, "first") ){
      first_blocks <- c(first_blocks, insert)
    }else if( identical(clause$kind, "last") ){
      last_blocks <- c(last_blocks, insert)
    }else{
      anchor_index <- which(vapply(order, function(item){
        matches <- clause$anchor %in% item

        return(matches)
      }, logical(1L)))
      if( length(anchor_index) != 1L ){
        stop("position: ", clause$kind, "() anchor '", clause$anchor,
             "' must name an unselected node in this role.",
             call. = FALSE)
      }
      if( identical(clause$kind, "before") ){
        order <- append(order, insert, after = anchor_index - 1L)
      }else if( identical(clause$kind, "after") ){
        prior <- after_counts[[clause$anchor]]
        if( is.null(prior) ) {
          prior <- 0L
        }
        order <- append(order, insert, after = anchor_index + prior)
        after_counts[[clause$anchor]] <- prior + length(insert)
      }else{
        stop("position: use first(), last(), before(), or after().",
             call. = FALSE)
      }
    }
  }
  positions <- c(first_blocks, order, last_blocks)

  return(positions)
}


#' Merge supplied positions and routing into a role
#'
#' @param existing The existing role specification.
#' @param supplied The supplied role specification.
#' @param position Position helpers or their name/value form, or NULL to append.
#' @param legacy_type The legacy first or last placement, or NULL.
#' @return The combined role specification.
#' @noRd
merge_role_spec <- function(existing, supplied, position = NULL,
                                legacy_type = NULL){
  if( is.null(existing$mapping_sources) ){
    existing$mapping_sources <- unique(as.character(existing$mapping_edges$v))
  }
  if( is.null(supplied$mapping_sources) ){
    supplied$mapping_sources <- unique(as.character(supplied$mapping_edges$v))
  }

  if( !is.null(legacy_type) ){
    placement <- position_blocks(existing, supplied, strict = FALSE)
    positions <- if( identical(legacy_type, "first") ){
      c(placement$blocks, placement$spine)
    }else{
      c(placement$spine, placement$blocks)
    }
  }else if( !is.null(position) ){
    placement <- position_blocks(existing, supplied, strict = TRUE)
    positions <- resolve_position_blocks(
      position, placement$blocks, placement$spine, supplied$nodes
    )
  }else{
    new_positions <- Filter(length, lapply(supplied$positions, function(item){
      new_nodes <- setdiff(item, existing$nodes)

      return(new_nodes)
    }))
    positions <- c(existing$positions, new_positions)
  }
  existing$positions <- positions
  existing$nodes <- unique(role_input_names(unlist(positions, use.names = FALSE)))

  for(side in c("incoming", "outgoing")){
    flag <- paste0(side, "_marked")
    if( isTRUE(supplied[[flag]]) ){
      routed_nodes <- if( isTRUE(existing[[flag]]) ){
        union(existing[[side]], supplied[[side]])
      }else{
        supplied[[side]]
      }
      existing[[side]] <- intersect(routed_nodes, existing$nodes)
      existing[[flag]] <- length(existing[[side]]) > 0L
    }
  }

  retained_mapping <- existing$mapping_sources
  retained_edges <- existing$mapping_edges[
    existing$mapping_edges$v %in% retained_mapping,
  ]
  existing$mapping_sources <- union(retained_mapping,
                                    supplied$mapping_sources)
  existing$mapping_edges <- bind_role_edges(retained_edges,
                                            supplied$mapping_edges)
  existing$mapping <- length(existing$mapping_sources) > 0L
  existing$same_role <- !identical(existing$role, "observed")

  return(existing)
}


#' Update one role in the graph specifications
#'
#' @param dag A dagitty object.
#' @param supplied The supplied role specification.
#' @param role The role to update.
#' @param type The graph construction type.
#' @param position Position helpers or their name/value form, or NULL to append.
#' @return Named list of graph role specifications with the selected role
#'   updated.
#' @noRd
specs_after_role_operation <- function(dag, supplied, role, type,
                                           position = NULL){
  specs <- existing_role_specs(dag)
  specs[[role]] <- merge_role_spec(
    specs[[role]], supplied, position = position,
    legacy_type = if( type %in% c("first", "last") ) { type } else { NULL }
  )

  return(specs)
}


#' Connect selected nodes with bidirected edges
#'
#' @param dag A dagitty object.
#' @param nodes Node names to connect.
#' @return Data table of bidirected edges joining the selected nodes to one
#'   another and to all other graph nodes.
#' @noRd
full_pair_edges <- function(dag, nodes){
  nodes <- unique(role_input_names(nodes))
  others <- setdiff(names(dag), nodes)
  within <- empty_role_edges()
  if( length(nodes) > 1L ){
    pairs <- utils::combn(nodes, 2L)
    within <- data.table::data.table(
      v = pairs[1L,], e = "<->", w = pairs[2L,]
    )
  }
  edges <- bind_role_edges(role_edge_frame(nodes, others, "<->"), within)

  return(edges)
}


#' Plan edges affected by a role edit
#'
#' @param dag A dagitty object.
#' @param supplied The supplied role specification.
#' @param role The role being edited.
#' @param type The graph construction type.
#' @param position Position helpers or their name/value form, or NULL to append.
#' @param specs Updated role specifications, or NULL to compute them.
#' @return Data table of edges for the role operation, excluding new instrument
#'   edges.
#' @noRd
role_operation_edges <- function(dag, supplied, role, type,
                                     position = NULL,
                                     specs = NULL){
  before <- existing_role_specs(dag)
  if( is.null(specs) ){
    specs <- specs_after_role_operation(
      dag, supplied, role, type, position
    )
  }
  legacy <- type %in% c("first", "last")
  planner_type <- if( legacy ) { "ordered" } else { type }
  planning_specs <- specs

  boundary <- empty_role_edges()
  touch <- unique(c(
    supplied$nodes, supplied$incoming, supplied$outgoing,
    specs[[role]]$incoming, specs[[role]]$outgoing,
    tier_members(before[[role]], supplied$nodes)
  ))
  if( legacy ){
    planning_specs[[role]]$same_role <- FALSE
    moved <- tier_members(before[[role]], supplied$nodes)
    moved <- unique(c(moved, setdiff(supplied$nodes, before[[role]]$nodes)))
    spine <- Filter(function(item) {
      matches <- !any(item %in% moved)

      return(matches)
    },
                    before[[role]]$positions)
    if( length(spine) > 0L && length(moved) > 0L ){
      boundary <- if( identical(type, "first") ){
        role_edge_frame(moved, spine[[1L]], "->")
      }else{
        role_edge_frame(spine[[length(spine)]], moved, "->")
      }
    }
    touch <- union(touch, moved)
  }

  candidate <- if( identical(type, "full") ){
    full_pair_edges(dag, supplied$nodes)
  }else{
    plan_graph_edges(
      planning_specs,
      type = planner_type,
      existing_dag = dag,
      include_iv = FALSE
    )
  }
  instruments <- classify_node_roles(dag)$instrument
  instruments <- instruments[!is.na(instruments)]
  iv_source <- candidate$v %in% instruments |
    (candidate$e %in% c("<->", "--") & candidate$w %in% instruments)
  candidate <- candidate[
    (candidate$v %in% touch | candidate$w %in% touch) &
      !iv_source,
  ]
  boundary <- boundary[
    !boundary$v %in% instruments,
  ]
  edges <- bind_role_edges(candidate, boundary)

  return(edges)
}


## Position helpers ##
#' Position helpers for placing nodes
#'
#' These helpers specify where selected nodes sit in a temporal order when
#' passed to the \code{position} argument of \code{\link{add_nodes}}. Combine helpers with \code{c()}.
#'
#' Supplying \code{nodes} limits a helper to those new nodes. One helper may
#' omit \code{nodes} to place all remaining new nodes. If no helper does this,
#' the remaining nodes are added at the end of the temporal order.
#'
#' Attaching causaliflower masks \code{first()} and \code{last()} from
#' packages loaded earlier, such as data.table and dplyr. Use \code{::} to
#' call the versions from those packages.
#'
#' @param nodes Optional character vector of new node names to place. If
#'   omitted, the helper applies to all new nodes not named in another helper.
#' @param anchor The existing node to place new nodes before or after. If
#'   omitted, \code{before()} places them first and \code{after()} places them
#'   last.
#'
#' @return An object of class \code{cf_position}: a list of one placement
#'   clause. Combining helpers with \code{c()} returns one \code{cf_position}
#'   containing their clauses (see \code{\link{c.cf_position}}).
#'
#' @examples
#' # all new nodes before everything in their role
#' first()
#' # X3 first, X4 and X5 last, and all remaining nodes immediately after X2
#' c(first("X3"), last(c("X4", "X5")),
#'   after("X2"))
#'
#' @name position_helpers
NULL


#' Create one position clause
#'
#' @param kind Placement type: first, last, before or after.
#' @param anchor The reference node, or NULL for an end position.
#' @param nodes Node names to place, or NULL for all otherwise unplaced nodes.
#' @return A cf_position object containing one placement clause.
#' @noRd
position_clause <- function(kind, anchor, nodes) {
  if (!is.null(nodes))  {
    nodes  <- as.character(nodes)
  }
  if (!is.null(anchor)) {
    anchor <- as.character(anchor)[1]
  }
  clause <- structure(list(list(kind = kind, anchor = anchor, nodes = nodes)),
                      class = "cf_position")

  return(clause)
}


#' @rdname position_helpers
#' @export
first <- function(nodes = NULL) {
  clause <- position_clause("first", anchor = NULL, nodes = nodes)

  return(clause)
}


#' @rdname position_helpers
#' @export
last <- function(nodes = NULL) {
  clause <- position_clause("last", anchor = NULL, nodes = nodes)

  return(clause)
}


#' @rdname position_helpers
#' @export
before <- function(anchor = NULL, nodes = NULL) {
  if (is.null(anchor)) {
    clause <- position_clause("first", anchor = NULL, nodes = nodes)

    return(clause)
  }
  else                 {
    clause <- position_clause("before", anchor = anchor, nodes = nodes)

    return(clause)
  }
}


#' @rdname position_helpers
#' @export
after <- function(anchor = NULL, nodes = NULL) {
  if (is.null(anchor)) {
    clause <- position_clause("last", anchor = NULL, nodes = nodes)

    return(clause)
  }
  else                 {
    clause <- position_clause("after", anchor = anchor, nodes = nodes)

    return(clause)
  }
}


#' Combine position helpers
#'
#' Combines the results of position helpers into one \code{cf_position}
#' object.
#'
#' @param ... \code{cf_position} objects returned by the position helpers.
#' @return A \code{cf_position} object containing the supplied position instructions
#'   in order.
#' @seealso \code{\link{position_helpers}}
#' @export
c.cf_position <- function(...) {
  positions <- structure(do.call("c", lapply(list(...), unclass)), class = "cf_position")

  return(positions)
}


# Resolve position clauses against an existing node order.
#' Place nodes within an existing order
#'
#' @param position Position helpers or their name/value form.
#' @param new_nodes Node names to place.
#' @param spine The existing node order.
#' @return Character vector of node names in the resolved order.
#' @noRd
resolve_position <- function(position, new_nodes, spine) {
  clauses <- position_clauses(position)

  is_empty <- vapply(clauses, function(cl) {
    matches <- is.null(cl$nodes)

    return(matches)
  }, logical(1))
  if (sum(is_empty) > 1L) {
    stop("position: at most one empty (catch-all) clause is allowed; ",
         "only one position helper may omit its nodes argument.",
         call. = FALSE)
  }

  claimed <- unlist(lapply(clauses, function(cl) {
    nodes <- cl$nodes

    return(nodes)
  }), use.names = FALSE)
  dups <- unique(claimed[duplicated(claimed)])
  if (length(dups)) {
    stop("position: each new node may appear in at most one clause; '",
         paste(dups, collapse = "', '"), "' named more than once.", call. = FALSE)
  }
  not_new <- setdiff(claimed, new_nodes)
  if (length(not_new)) {
    stop("position: clause names a node that is not among the new nodes: '",
         paste(not_new, collapse = "', '"), "'.", call. = FALSE)
  }

  remainder <- new_nodes[!(new_nodes %in% claimed)]
  if (any(is_empty)) {
    clauses[[which(is_empty)[1]]]$nodes <- remainder
  } else if (length(remainder)) {
    # Without a catch-all, unclaimed nodes are appended.
    clauses <- c(clauses, list(list(kind = "last", anchor = NULL, nodes = remainder)))
  }

  order <- spine
  after_counts <- list()
  for (cl in clauses) {
    if (cl$kind == "before") {
      if (!(cl$anchor %in% order)) {
        stop("position: before() anchor '", cl$anchor,
             "' is not an existing node in this role.", call. = FALSE)
      }
      order <- append(order, cl$nodes, after = match(cl$anchor, order) - 1L)
    } else if (cl$kind == "after") {
      if (!(cl$anchor %in% order)) {
        stop("position: after() anchor '", cl$anchor,
             "' is not an existing node in this role.", call. = FALSE)
      }
      prior <- after_counts[[cl$anchor]]
      if( is.null(prior) ) {
        prior <- 0L
      }
      order <- append(order, cl$nodes,
                      after = match(cl$anchor, order) + prior)
      after_counts[[cl$anchor]] <- prior + length(cl$nodes)
    }
  }

  firsts <- unlist(lapply(clauses, function(cl){
    if( cl$kind == "first" ){
      nodes <- cl$nodes

      return(nodes)
    }else{
      nodes <- NULL
      return(invisible(nodes))
    }
  }), use.names = FALSE)
  lasts  <- unlist(lapply(clauses, function(cl){
    if( cl$kind == "last" ){
      nodes <- cl$nodes

      return(nodes)
    }else{
      nodes <- NULL
      return(invisible(nodes))
    }
  }), use.names = FALSE)
  positions <- c(firsts, order, lasts)

  return(positions)
}
