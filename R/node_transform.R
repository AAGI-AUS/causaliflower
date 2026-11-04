# ==============================================================================
# node_transform.R


# Node transformation functions and their edge projection helpers
# ==============================================================================


#' Collapse nodes in a dagitty object
#'
#' Removes the named nodes and projects the graph onto the remaining nodes.
#' Directed paths whose intermediate nodes are all removed become directed
#' edges. Omitted common causes can produce bidirected edges. Other nodes
#' and the existing edges between them are kept.
#'
#' @param dag A dagitty object of type `dag`, with an acyclic directed part and
#'   only directed or bidirected edges between different nodes. Adjusted and
#'   selected statuses are not supported because this function does not perform
#'   conditional projection.
#' @param nodes Character vector of node names to remove. Repeated names are
#'   ignored. Names not present in the graph are ignored with a warning.
#' @param print_edges Print newly induced edges in the console. Defaults to
#'   `FALSE`.
#' @returns A dagitty object projected onto the remaining nodes.
#' @seealso [remove_nodes()], [add_nodes()], [get_identification()]
#'
#' @details
#' This operation marginalises the named variables. It does not condition on
#' or intervene on them. Replacement edges depend on the paths through the
#' removed nodes. A directed and a bidirected edge can occur between the same
#' pair of remaining nodes.
#'
#' Only nodes named in `nodes` are removed. Latent status does not select a
#' node automatically. Statuses attached to removed nodes are removed with
#' them. No replacement treatment or outcome is chosen. Remaining nodes keep
#' their coordinates, and surviving edges keep their control points. New
#' edges are straight.
#'
#' A new directed edge may represent a path through omitted variables. It does
#' not establish a direct mechanism in the original graph. The projected graph
#' does not retain the omitted variables' identities, their mechanisms, or
#' every constraint of the original model with latent variables.
#'
#' @examples
#' dag <- dagitty::dagitty("dag { T [exposure] Y [outcome] T -> M -> Y }")
#' dag <- collapse_nodes(dag, "M")
#'
#' dag <- dagitty::dagitty("dag { U [latent] U -> A; U -> B }")
#' dag <- collapse_nodes(dag, "U")
#'
#' @export
collapse_nodes <- function(dag, nodes, print_edges = FALSE){
  existing_edges <- collapse_validate_graph(dag)
  if( !is.character(nodes) ){
    stop("collapse_nodes(): 'nodes' must be a character vector.", call. = FALSE)
  }
  if( length(nodes) == 0L ){
    stop("collapse_nodes(): 'nodes' must contain at least one node name.",
         call. = FALSE)
  }
  if( anyNA(nodes) || any(!nzchar(nodes)) ){
    stop("collapse_nodes(): 'nodes' must contain non-missing, non-empty names.",
         call. = FALSE)
  }
  if( !is.logical(print_edges) || length(print_edges) != 1L ||
      is.na(print_edges) ){
    stop("collapse_nodes(): 'print_edges' must be TRUE or FALSE.", call. = FALSE)
  }
  nodes <- unique(as.character(nodes))
  missing_nodes <- setdiff(nodes, names(dag))
  if( length(missing_nodes) > 0L ){
    warning("collapse_nodes(): nodes not present in the graph were ignored: ",
            paste(missing_nodes, collapse = ", "), ".", call. = FALSE)
  }
  nodes <- intersect(nodes, names(dag))
  if( length(nodes) == 0L ){
    return(dag)
  }

  projected <- collapse_project_edges(names(dag), existing_edges, nodes)
  new_edges <- projected[!edge_key(projected) %in% edge_key(existing_edges), ,
                         drop = FALSE]
  dag <- native_remove_nodes(dag, nodes, "collapse_nodes")
  dag <- native_append_graph(dag, character(), new_edges,
                             caller = "collapse_nodes")
  if( is.null(dag) ){
    stop("collapse_nodes(): the projected edges could not be added while preserving the graph.",
         call. = FALSE)
  }
  if( isTRUE(print_edges) && nrow(new_edges) > 0L ){
    printed <- print_edges_helper(new_edges)
    cat(paste("c(", paste(unlist(printed), collapse = ",\n\n"), ")",
              sep = "\n", collapse = ""))
    message("\nPrinted new edges - copy and paste to use as a vector object.\n")
  }

  return(dag)
}


#' Check the graph supplied for marginalisation
#'
#' @param dag A dagitty object.
#' @return The graph's directed and bidirected edge table.
#' @noRd
collapse_validate_graph <- function(dag){
  graph_type <- dag_graph_type(dag)
  if( !identical(graph_type, "dag") ){
    stop("collapse_nodes(): 'dag' must have native graph type 'dag'.",
         call. = FALSE)
  }
  graph <- dagitty::dagitty(as.character(dag))
  edges <- dag_edges(graph, "collapse_nodes()")
  if( any(!edges$e %in% c("->", "<->")) ){
    stop("collapse_nodes(): only directed and bidirected edges are supported.",
         call. = FALSE)
  }
  if( any(edges$v == edges$w) ){
    stop("collapse_nodes(): self-loops are not supported.", call. = FALSE)
  }
  if( !isTRUE(dag_is_acyclic(graph)) ){
    stop("collapse_nodes(): the directed part of 'dag' must be acyclic.",
         call. = FALSE)
  }
  if( length(dagitty::adjustedNodes(graph)) > 0L ||
      length(selected_nodes(graph)) > 0L ){
    stop("collapse_nodes(): adjusted and selected statuses are not supported; this function does not perform conditional projection.",
         call. = FALSE)
  }

  return(edges)
}


#' Project directed and bidirected edges through omitted vertices
#'
#' @param vertices All graph node names.
#' @param edges Directed and bidirected edge table.
#' @param nodes Existing node names to omit.
#' @return An edge table on the remaining vertices.
#' @noRd
collapse_project_edges <- function(vertices, edges, nodes){
  vertices <- sort(as.character(vertices))
  directed <- matrix(FALSE, length(vertices), length(vertices),
                     dimnames = list(vertices, vertices))
  bidirected <- directed
  endpoints <- cbind(match(edges$v, vertices), match(edges$w, vertices))
  arrows <- edges$e == "->"
  if( any(arrows) ){
    directed[endpoints[arrows, , drop = FALSE]] <- TRUE
  }
  if( any(!arrows) ){
    bidirected[endpoints[!arrows, , drop = FALSE]] <- TRUE
    bidirected <- bidirected | t(bidirected)
  }

  for( node in vertices[vertices %in% nodes] ){
    parents <- which(directed[, node])
    children <- which(directed[node, ])
    spouses <- which(bidirected[node, ])
    if( length(parents) > 0L && length(children) > 0L ){
      directed[parents, children] <- TRUE
    }
    if( length(children) > 1L ){
      bidirected[children, children] <- TRUE
    }
    if( length(spouses) > 0L && length(children) > 0L ){
      bidirected[spouses, children] <- TRUE
      bidirected[children, spouses] <- TRUE
    }
    diag(bidirected) <- FALSE
    directed[node, ] <- FALSE
    directed[, node] <- FALSE
    bidirected[node, ] <- FALSE
    bidirected[, node] <- FALSE
  }

  arrows <- which(directed, arr.ind = TRUE)
  arcs <- which(bidirected & upper.tri(bidirected), arr.ind = TRUE)
  projected <- rbind(
    data.frame(v = vertices[arrows[, 1L]], e = rep("->", nrow(arrows)),
               w = vertices[arrows[, 2L]], stringsAsFactors = FALSE),
    data.frame(v = vertices[arcs[, 1L]], e = rep("<->", nrow(arcs)),
               w = vertices[arcs[, 2L]], stringsAsFactors = FALSE)
  )
  projected <- projected[order(edge_key(projected)), , drop = FALSE]
  rownames(projected) <- NULL

  return(projected)
}


## Merging and recombination ##
#' Merge nodes in a dagitty object
#'
#' Combines named nodes into one node and redirects their edges to that node.
#' Edges within the group are removed. Repeated edges of the same type become
#' one edge. Directed and bidirected edges between the same pair are kept.
#'
#' @param dag An acyclic dagitty DAG with directed or bidirected edges between
#'   different nodes.
#' @param nodes Character vector of at least two existing node names. The
#'   first node supplies the merged node's coordinates and native metadata,
#'   even when `name` reuses another node's name.
#' @param name Name for the merged node. It may reuse a name in `nodes`, but
#'   must not name another node in the graph.
#' @param print_edges Print redirected edges in the console. Defaults to FALSE.
#' @returns A dagitty object with the named nodes combined into one node.
#' @seealso [recombine_nodes()], [collapse_nodes()], [remove_nodes()]
#'
#' @details
#' Use [recombine_nodes()] to merge groups with matching connections
#' automatically, or preview those groups before choosing which to merge.
#' Theoretical support for combining the variables must be assessed separately.
#' Merging changes the variable set and can create paths that were absent from
#' the original graph. A merge that creates a directed cycle is rejected.
#'
#' The selected nodes must have the same treatment, outcome, latent, adjusted
#' and selected statuses. The merged node retains these statuses and the first
#' node's coordinates, labels and notes. Other nodes and their metadata are
#' unchanged. No layout is regenerated.
#'
#' Edges retained from the first node keep their native attributes and control
#' points. Other redirected edges are new, straight edges. If several edges
#' become the same edge, the first node's edge takes precedence. Metadata owned
#' by the other removed nodes and edges is discarded.
#'
#' @examples
#' dag <- dagitty::dagitty("dag { T -> M1 -> Y; T -> M2 -> Y }")
#' dag <- merge_nodes(dag, c("M1", "M2"), "M")
#'
#' @export
merge_nodes <- function(dag, nodes, name, print_edges = FALSE){
  caller <- "merge_nodes"
  edges <- merge_graph_edges(dag, caller)
  nodes <- merge_node_names(nodes, names(dag), caller)
  if( length(nodes) < 2L ){
    stop("merge_nodes(): 'nodes' must contain at least two different names.",
         call. = FALSE)
  }
  if( !is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(name) ){
    stop("merge_nodes(): 'name' must be one non-missing, non-empty name.",
         call. = FALSE)
  }
  if( name %in% setdiff(names(dag), nodes) ){
    stop("merge_nodes(): 'name' already belongs to a node outside the group.",
         call. = FALSE)
  }
  quote_dag_names(name, caller)
  if( !is.logical(print_edges) || length(print_edges) != 1L ||
      is.na(print_edges) ){
    stop("merge_nodes(): 'print_edges' must be TRUE or FALSE.", call. = FALSE)
  }
  statuses <- merge_node_statuses(dag)
  for( status in statuses ){
    flags <- nodes %in% status
    if( any(flags) && !all(flags) ){
      stop("merge_nodes(): nodes in the group must have matching native statuses.",
           call. = FALSE)
    }
  }

  merged_edges <- edges
  merged_edges$v[merged_edges$v %in% nodes] <- name
  merged_edges$w[merged_edges$w %in% nodes] <- name
  changed <- (edges$v %in% nodes & edges$v != name) |
    (edges$w %in% nodes & edges$w != name)
  reported <- merged_edges[changed & merged_edges$v != merged_edges$w, ,
                           drop = FALSE]
  reported <- reported[!duplicated(edge_key(reported)), , drop = FALSE]
  merged_edges <- merged_edges[merged_edges$v != merged_edges$w, , drop = FALSE]
  merged_edges <- merged_edges[!duplicated(edge_key(merged_edges)), , drop = FALSE]
  remaining_nodes <- c(setdiff(names(dag), nodes), name)
  first <- nodes[1L]
  dag <- native_remove_nodes(dag, nodes[-1L], caller)
  dag <- copy_rename_native_nodes(dag, stats::setNames(name, first), caller)
  retained_edges <- dag_edges(dag, caller)
  new_edges <- merged_edges[!edge_key(merged_edges) %in% edge_key(retained_edges), ,
                            drop = FALSE]
  merged <- native_append_graph(dag, character(), new_edges, caller = caller)
  if( is.null(merged) ){
    stop("merge_nodes(): the remaining graph could not be preserved.",
         call. = FALSE)
  }
  if( !setequal(names(merged), remaining_nodes) ||
      !setequal(edge_key(dag_edges(merged, caller)), edge_key(merged_edges)) ){
    stop("merge_nodes(): the requested nodes and edges could not be preserved.",
         call. = FALSE)
  }
  if( !isTRUE(dag_is_acyclic(merged)) ){
    stop("merge_nodes(): merging these nodes would create a directed cycle.",
         call. = FALSE)
  }
  if( isTRUE(print_edges) && nrow(reported) > 0L ){
    printed <- print_edges_helper(reported)
    cat("merge_nodes(): redirected edges:\n")
    print(printed)
  }

  return(merged)
}


#' Recombine nodes with matching connections
#'
#' Finds and merges groups of nodes with identical directed parents and
#' children, matching bidirected neighbours outside the group, and the same
#' native statuses. Within a group, every pair must be bidirected or every
#' pair must be unconnected.
#'
#' @param dag An acyclic dagitty DAG with directed or bidirected edges between
#'   different nodes.
#' @param nodes Optional character vector restricting the candidate nodes.
#'   Defaults to NULL, which considers all eligible nodes.
#' @param preview Return the proposed groups instead of the recombined graph.
#'   Defaults to FALSE. The same validation is performed in either mode.
#' @returns A dagitty object with each matching group combined into one node.
#'   With `preview = TRUE`, a named list of character vectors containing the
#'   original members of each group. List names give the proposed node names. The first member supplies the node's metadata. If no group matches,
#'   returns the original graph, or an empty list when previewing.
#' @seealso [merge_nodes()]
#'
#' @details
#' Treatment, outcome, adjusted and selected nodes are excluded, as are
#' isolated nodes. Latent and observed nodes are considered in separate
#' groups. Neighbours outside `nodes` are included when comparing connections.
#' If `nodes` is supplied, only its members are eligible to be merged.
#'
#' Groups are found once in the original graph. Recombination does not repeat
#' if the resulting graph contains new matching groups. Within each group,
#' node names are sorted and joined with underscores to name the
#' combined node. If that name is already used, suffixes `_2`, `_3`, and so on
#' are tried until the name is unused. Every original node name is reserved,
#' along with names assigned to preceding groups.
#'
#' Each group is combined with [merge_nodes()]. The first member supplies
#' coordinates and native metadata, including attributes and control points
#' on surviving edges. Other redirected edges are straight. No layout is
#' regenerated and no edge report is printed.
#'
#' The function finds groups for the recombination step of ESC-DAGs by
#' comparing their connections. You still need to decide whether the variables
#' should be combined. Use `preview = TRUE` to inspect the groups. To select groups
#' or choose different names, pass their members to [merge_nodes()].
#'
#' @examples
#' dag <- dagitty::dagitty("dag { T -> M1 -> Y; T -> M2 -> Y }")
#' groups <- recombine_nodes(dag, preview = TRUE)
#' merged <- merge_nodes(dag, groups[[1]], name = "M")
#' dag <- recombine_nodes(dag)
#'
#' @export
recombine_nodes <- function(dag, nodes = NULL, preview = FALSE){
  caller <- "recombine_nodes"
  if( !is.logical(preview) || length(preview) != 1L || is.na(preview) ){
    stop("recombine_nodes(): 'preview' must be TRUE or FALSE.", call. = FALSE)
  }
  edges <- merge_graph_edges(dag, caller)
  groups <- merge_groups(dag, nodes, edges, caller)
  group_names <- character(length(groups))
  used_names <- names(dag)
  for( i in seq_along(groups) ){
    base_name <- paste(groups[[i]], collapse = "_")
    name <- base_name
    suffix <- 2L
    while( name %in% used_names ){
      name <- paste(base_name, suffix, sep = "_")
      suffix <- suffix + 1L
    }
    quote_dag_names(name, caller)
    group_names[i] <- name
    used_names <- c(used_names, name)
  }
  if( length(groups) > 0L ){
    names(groups) <- group_names
  }
  merged <- dag
  for( i in seq_along(groups) ){
    merged <- merge_nodes(merged, groups[[i]], name = group_names[i])
  }
  if( preview ){
    return(groups)
  }

  return(merged)
}


#' Find groups with matching connections for recombination
#'
#' @param dag The validated dagitty DAG.
#' @param nodes Candidate node names, or NULL to consider all eligible nodes.
#' @param edges The graph's directed and bidirected edge table.
#' @param caller Function name for errors.
#' @returns A list of character vectors containing matching groups, with
#'   members and groups ordered by node name.
#' @noRd
merge_groups <- function(dag, nodes, edges, caller){
  all_nodes <- as.character(names(dag))
  if( is.null(nodes) ){
    nodes <- all_nodes
  }else{
    nodes <- merge_node_names(nodes, all_nodes, caller)
  }
  statuses <- merge_node_statuses(dag)
  protected <- unique(c(statuses$exposure, statuses$outcome,
                        statuses$adjusted, statuses$selected))
  connected <- unique(c(edges$v, edges$w))
  nodes <- as.character(setdiff(intersect(nodes, connected), protected))
  nodes <- sort(nodes, method = "radix")
  if( length(nodes) < 2L ){
    groups <- list()
    return(groups)
  }

  directed <- edges[edges$e == "->", , drop = FALSE]
  bidirected <- edges[edges$e == "<->", , drop = FALSE]
  parents <- split(match(directed$v, all_nodes), directed$w)
  children <- split(match(directed$w, all_nodes), directed$v)
  spouses <- split(match(c(bidirected$v, bidirected$w), all_nodes),
                   c(bidirected$w, bidirected$v))
  groups <- list()
  for( include_self in c(FALSE, TRUE) ){
    signatures <- vapply(nodes, function(node){
      neighbours <- spouses[[node]]
      if( include_self ){
        neighbours <- c(neighbours, match(node, all_nodes))
      }
      signature <- paste(
        paste(sort(unique(parents[[node]])), collapse = ","),
        paste(sort(unique(children[[node]])), collapse = ","),
        paste(sort(unique(neighbours)), collapse = ","),
        node %in% statuses$latent, sep = "|"
      )

      return(signature)
    }, character(1L), USE.NAMES = FALSE)
    matches <- unname(split(nodes, match(signatures, unique(signatures))))
    groups <- c(groups, matches[lengths(matches) > 1L])
  }
  first <- vapply(groups, function(group){
    index <- match(group[1L], nodes)

    return(index)
  }, integer(1L))
  groups <- groups[order(first)]

  return(groups)
}


#' Validate graphs used for node merging
#'
#' @param dag A dagitty object.
#' @param caller Function name for errors.
#' @returns A table of directed and bidirected edges.
#' @noRd
merge_graph_edges <- function(dag, caller){
  if( dag_graph_type(dag) != "dag" ){
    stop(caller, "(): 'dag' must use the dag graph type.", call. = FALSE)
  }
  edges <- dag_edges(dag, caller)
  if( any(!edges$e %in% c("->", "<->")) ){
    stop(caller, "(): only directed and bidirected edges are supported.",
         call. = FALSE)
  }
  if( any(edges$v == edges$w) ){
    stop(caller, "(): self-loops are not supported.", call. = FALSE)
  }
  if( !isTRUE(dag_is_acyclic(dag)) ){
    stop(caller, "(): 'dag' must be acyclic.", call. = FALSE)
  }

  return(edges)
}


#' Validate existing names for node merging
#'
#' @param nodes Character vector of names.
#' @param existing Available node names.
#' @param caller Function name for errors.
#' @returns Unique node names in the supplied order.
#' @noRd
merge_node_names <- function(nodes, existing, caller){
  if( !is.character(nodes) || length(nodes) == 0L ||
      anyNA(nodes) || any(!nzchar(nodes)) ){
    stop(caller, "(): 'nodes' must contain non-missing, non-empty node names.",
         call. = FALSE)
  }
  absent <- setdiff(nodes, existing)
  if( length(absent) > 0L ){
    stop(caller, "(): nodes not present in the graph: ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  nodes <- unique(as.character(nodes))

  return(nodes)
}


#' Read native statuses for node merging
#'
#' @param dag A dagitty object.
#' @returns A named list of node names for each native status.
#' @noRd
merge_node_statuses <- function(dag){
  parsed <- dagitty::dagitty(as.character(dag))
  statuses <- list(
    exposure = dagitty::exposures(parsed),
    outcome = dagitty::outcomes(parsed),
    latent = dagitty::latents(parsed),
    adjusted = dagitty::adjustedNodes(parsed),
    selected = selected_nodes(parsed)
  )

  return(statuses)
}
