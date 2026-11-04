# ==============================================================================
# draw_funcs.r


# Edge table construction for graph_funcs.r and graph_planner.R
# ==============================================================================


#' Connect every ancestor to every descendant
#'
#' @param ancestors Source node names; missing values are ignored.
#' @param descendants Target node names; missing values are ignored.
#' @param e Edge symbol to use for each pair.
#' @return A data.table with `ancestor`, `edge`, and `descendant` columns.
#' @noRd
cross_edges <- function(ancestors, descendants, e) {
  ancestors   <- ancestors[ !is.na(ancestors) ]
  descendants <- descendants[ !is.na(descendants) ]
  if (length(ancestors) == 0L || length(descendants) == 0L) {
    edges <- data.table::data.table(ancestor = character(0), edge = character(0), descendant = character(0))
    return(edges)
  }
  g <- expand.grid(ancestor = ancestors, descendant = descendants,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  edges <- data.table::data.table(ancestor = as.character(g$ancestor), edge = e, descendant = as.character(g$descendant))

  return(edges)
}


#' Helper for connect_new_nodes()
#'
#' @param order Resolved sequence of existing and new nodes in the graph.
#' @param new_nodes Names of the new nodes.
#' @param type `"chain"` joins adjacent positions; other types join every
#'   earlier and later pair.
#' @param e Edge symbol supplied by the caller.
#' @return A data.table of edges touching at least one new node, with
#'   `ancestor`, `edge`, and `descendant` columns.
#' @noRd
connect_in_order <- function(order, new_nodes, type, e) {
  n <- length(order)
  empty <- data.table::data.table(ancestor = character(0), edge = character(0), descendant = character(0))
  if (n < 2L) {
    return(empty)
  }
  if (type == "chain") {
    i <- seq_len(n - 1L)
    a <- order[i]
    b <- order[i + 1L]
  } else {
    g <- which(upper.tri(matrix(0L, n, n)), arr.ind = TRUE)  # row < col -> earlier, later
    a <- order[g[, "row"]]
    b <- order[g[, "col"]]
  }
  keep <- (a %in% new_nodes) | (b %in% new_nodes)             # only edges involving a new node
  a <- a[keep]
  b <- b[keep]
  if (length(a) == 0L) {
    return(empty)
  }
  edges <- data.table::data.table(ancestor = a, edge = e, descendant = b)

  return(edges)
}


#' Draw observed node edges
#'
#' draw_observed_edges() connects observed nodes to their parents and children.
#'
#' @importFrom data.table as.data.table
#' @param observed Vector of variables without roles, also not unobserved (latent).
#' @param existing_dag An existing dagitty object may be supplied.
#' @param e Edge type, such as "->" for a directed edge or "<->" for a bidirected edge.
#' @returns Data frame of edges.
#' @noRd
draw_observed_edges <- function(observed,
                                existing_dag = NA,
                                e
                                ){

  if( all( complete.cases(observed) ) ) {

    ## connect observed to parents and children ##
    observed_df <- data.frame(NULL)

    # connect children
    existing_edges <- dagitty::edges(existing_dag)
    observed_children <- children_edge_order(existing_edges, observed) # children in dagitty's edge-insertion order (feeds cross_edges row order)

    if( length(observed_children) > 0 ){

      observed_df <- cross_edges(observed, observed_children, e)

    }

    # connect parents
    observed_parents <- parents_edge_order(existing_edges, observed)

    if( length(observed_parents) > 0 ){

      observed_unlist <- cross_edges(observed_parents, observed, e)

      observed_df <- rbind(observed_df, observed_unlist)

    }

    return(observed_df)
  }

  observed_df <- data.frame(NULL)

  return( observed_df )
}


#' Draw instrumental variable edges
#'
#' draw_iv_edges() generates edges from instrumental variable inputs.
#'
#' @importFrom data.table as.data.table
#' @param type Graph construction type: `"full"`, `"chain"`, `"first"`,
#'   `"last"`, `"ordered"`, or `"saturated"`.
#' @param instrumental_variables Inputted list or vector of instrumental variables.
#' @param treatments Vector of treatments.
#' @param e Edge type, such as "->" for a directed edge or "<->" for a bidirected edge.
#' @returns A data frame of instrumental variable edges.
#' @noRd
draw_iv_edges <- function(type,
                          instrumental_variables,
                          treatments,
                          e
                          ){
  empty <- cross_edges(character(), character(), e)
  instrumental_flat <- unlist(instrumental_variables)
  if( length(instrumental_flat) == 0L || all(is.na(instrumental_flat)) ){
    return(empty)
  }
  nested <- length(instrumental_flat) > length(instrumental_variables)
  if( !nested && length(treatments) == 0L ){
    return(empty)
  }

  if( all( complete.cases( unlist(instrumental_variables) ) ) ){

    instrumental_list <- c()

    if( length(unlist(instrumental_variables)) > length(instrumental_variables) ){

      instrumental_list <- suppressWarnings( lapply( seq_along(instrumental_variables), function(x){

        target_edges <- lapply(seq_along( instrumental_variables[[x]][[2]] ), function(y){

          instrument_edges <- list( c( ancestor = instrumental_variables[[x]][[1]], edge = e, descendant = instrumental_variables[[x]][[2]][[y]] ) )

          return(instrument_edges)

        })

        return(invisible(target_edges))

      } ) )

      instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
      if( length(instrumental_list) == 0L ) {
        return(empty)
      }
      instrumental_df <- as.data.table( do.call( rbind, unlist(instrumental_list, recursive = FALSE) ) )

    }else{

      instrumental_vec <- as.vector( unlist( lapply( instrumental_variables, function(x) if( identical( x, character(0) ) ) {
        instrument <- NA_character_

        return(instrument)
      } else {

        return(x)
      } ) ) )

      if( any( type == "full" | type == "saturated" ) ){

        # connect treatments
        instrumental_list <- suppressWarnings( lapply(seq_along(instrumental_vec), function(x){

          target_edges <- lapply(seq_along(treatments), function(y){

            instrument_edges <- list( c( ancestor = instrumental_vec[x], edge = e, descendant = treatments[y]) )

            return(instrument_edges)

          })

          return(invisible(target_edges))

        }) )

        instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
        if( length(instrumental_list) == 0L ) {
          return(empty)
        }
        instrumental_df <- as.data.table( do.call( rbind, unlist(instrumental_list, recursive = FALSE) ) )

      }else{

        len_trt <- length(treatments)

        # connect treatments
        instrumental_list <- suppressWarnings( lapply(seq_along(instrumental_vec), function(x){

          if( x <= len_trt ){

            instrument_edges <- list( c( ancestor = instrumental_vec[x], edge = e, descendant = treatments[x]) )

            return(invisible(instrument_edges))

          }else{

            instrument_edges <- list( c( ancestor = instrumental_vec[x], edge = e, descendant = treatments[1]) )

            return(invisible(instrument_edges))

          }

        } ) )

      }

      instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
      if( length(instrumental_list) == 0L ) {
        return(empty)
      }
      instrumental_df <- as.data.table( do.call( rbind, instrumental_list ) )

    }

    return( instrumental_df )
  }

  return(empty)
}


#' Helper for add_nodes()
#'
#' connect_new_nodes() connects new nodes using explicit parents and children,
#' a position in the graph, or the instrument role.
#'
#' @importFrom data.table as.data.table
#' @param dag An existing dagitty object.
#' @param new_nodes A vector of new nodes.
#' @param ancestors Parent node names for the new nodes.
#' @param descendants Child node names for the new nodes.
#' @param node_roles Optional `"instrument"` role. Other roles are planned by
#'   role_operation_edges().
#' @param type Graph type passed by add_nodes(), including the deprecated
#'   `"first"` and `"last"` values.
#' @param position Optional position in the graph's topological order.
#' @return An edge table with columns v, e and w.
#' @noRd
connect_new_nodes <- function(dag,
                              new_nodes,
                              ancestors,
                              descendants,
                              node_roles,
                              type,
                              position = NULL
                              ){
  e <- "->"

  if(type == "full"){

    e <- "<->"
  }

  if( length(node_roles) > 0 ){

    output_list <- add_nodes_helper(dag = dag,
                                    nodes = new_nodes,
                                    node_role = node_roles,
                                    type = type,
                                    position = position)
    edges <- output_list$new_edges

    colnames(edges) <- c("v", "e", "w")

    return( edges )

  }

  ## without a role, position sets direction in the whole graph's order ##
  ## explicit parents or children take precedence; full edges stay bidirected ##
  if( !is.null(position) && is.null(ancestors) && is.null(descendants) ){
    spine <- names( sort( unlist( dagitty::topologicalOrdering(dag) ) ) )
    ord   <- resolve_position(position, new_nodes, spine)
    edges <- connect_in_order(ord, new_nodes, "saturated", e)
    colnames(edges) <- c("v", "e", "w")
    return( edges )
  }

  node_list <- c()
  edges <- c()

  ## get node names
  node_names <- names( dag )

  ## validate connection targets: unknown names must not silently create nodes or trigger saturation
  unknown_nodes <- c(ancestors, descendants)
  unknown_nodes <- unknown_nodes[ ! unknown_nodes %in% c(node_names, new_nodes) ]

  if( length(unknown_nodes) > 0 ){

    stop(paste0("Parent or child nodes not present in the dag: ",
                paste(unknown_nodes, collapse = ", "),
                ". Please check inputs and try again."))

  }

  if( is.null(ancestors) & is.null(descendants) ){

    descendants <- node_names # default: connect new nodes to every existing node

  }else if( length(ancestors) == 0 & length(descendants) == 0 ){

    stop("'parents' and 'children' are both empty. Supply at least one connection, or leave both NULL to connect new nodes to all existing nodes. Please check inputs and try again.")

  }

  if( length( descendants ) > 0){
    ## connect new_nodes to descendants
    node_list <- suppressWarnings( lapply(seq_along(new_nodes), function(x){

      child_edges <- lapply(seq_along(descendants), function(y){

        child_edges <- list( c( v = new_nodes[x], e = e, w = descendants[y]) )

        return(child_edges)

      })

      return(invisible(child_edges))

    }) )

    node_list <- Filter(Negate(anyNA), unlist(node_list, recursive = FALSE))
    edges <- as.data.table( do.call( rbind, unlist(node_list, recursive = FALSE) ) )
  }

  if( length( ancestors ) > 0){
    ## connect new_nodes to ancestors
    node_list <- suppressWarnings( lapply(seq_along(ancestors), function(x){

      child_edges <- lapply(seq_along(new_nodes), function(y){

        child_edges <- list( c( v = ancestors[x], e = e, w = new_nodes[y]) )

        return(child_edges)

      })

      return(invisible(child_edges))

    }) )

    node_list <- Filter(Negate(anyNA), unlist(node_list, recursive = FALSE))
    node_unlist <- as.data.table( do.call( rbind, unlist(node_list, recursive = FALSE) ) )

    edges <- rbind(edges, node_unlist)
  }

  return( edges )
}
