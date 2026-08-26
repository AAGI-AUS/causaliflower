#' Cross-product edge table (internal)
#'
#' Builds an ancestor-edge-descendant table for every pair, dropping missing
#' endpoints and returning a typed zero-row table when either side is empty.
#' @noRd
.cross_edges <- function(ancestors, descendants, e) {
  ancestors   <- ancestors[ !is.na(ancestors) ]
  descendants <- descendants[ !is.na(descendants) ]
  if (length(ancestors) == 0L || length(descendants) == 0L) {
    return(data.table::data.table(ancestor = character(0), edge = character(0), descendant = character(0)))
  }
  g <- expand.grid(ancestor = ancestors, descendant = descendants,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  data.table::data.table(ancestor = as.character(g$ancestor), edge = e, descendant = as.character(g$descendant))
}


#' Same-role edges along a resolved order (internal)
#'
#' Given the full resolved order of a role's nodes (existing + new interleaved)
#' and the set of new node names, returns the directed same-role edges that
#' touch at least one new node, per `type`. Existing<->existing edges are never
#' emitted because they are already in the graph.
#'   - type "ordered": consecutive pairs only (a chain along the order);
#'   - type "saturated"/"full": every earlier -> later pair.
#' The edge symbol `e` ("->", or "<->" for full) is supplied by the caller so it
#' matches the role builder.
#' @noRd
.connect_in_order <- function(order, new_nodes, type, e) {
  n <- length(order)
  empty <- data.table::data.table(ancestor = character(0), edge = character(0), descendant = character(0))
  if (n < 2L) return(empty)
  if (type == "ordered") {
    i <- seq_len(n - 1L); a <- order[i]; b <- order[i + 1L]
  } else {
    g <- which(upper.tri(matrix(0L, n, n)), arr.ind = TRUE)  # row < col -> earlier, later
    a <- order[g[, "row"]]; b <- order[g[, "col"]]
  }
  keep <- (a %in% new_nodes) | (b %in% new_nodes)             # only edges involving a new node
  a <- a[keep]; b <- b[keep]
  if (length(a) == 0L) return(empty)
  data.table::data.table(ancestor = a, edge = e, descendant = b)
}


#' Place new same-role nodes by directing edges along a node order (internal)
#'
#' Inserts new nodes into the role's existing order and draws same-role edges
#' according to \code{type}. Cross-role edges from the role builder are retained.
#'
#' Returns a list(new_edges, temporal_reference_node).
#' @noRd
.place_role_edges <- function(new_edges, new_nodes, role_vec, nodes_ordered,
                              type, e, position, temporal_reference_node) {
  role_vec <- role_vec[!is.na(role_vec)]
  full_spine <- names(nodes_ordered)[ names(nodes_ordered) %in% role_vec ]
  # connect_nodes() may supply names already present in the role order.
  if (type %in% c("first", "last")) {
    if (!is.null(position)) {
      stop("Use either legacy type = \"first\"/\"last\" or position, not both.",
           call. = FALSE)
    }
    if (length(full_spine) == 0L) {
      return(list(new_edges = new_edges,
                  temporal_reference_node = temporal_reference_node))
    }
    if (type == "first") {
      same <- .cross_edges(new_nodes, full_spine[1L], e)
      temporal_reference_node <- full_spine[1L]
    } else {
      same <- .cross_edges(full_spine[length(full_spine)], new_nodes, e)
      temporal_reference_node <- full_spine[length(full_spine)]
    }
    return(list(new_edges = rbind(new_edges, same),
                temporal_reference_node = temporal_reference_node))
  }

  spine <- full_spine[!full_spine %in% new_nodes]
  ord  <- if (is.null(position)) c(spine, new_nodes) else .resolve_position(position, new_nodes, spine)
  keep <- !(new_edges$ancestor %in% new_nodes & new_edges$descendant %in% new_nodes)
  same <- .connect_in_order(ord, new_nodes, type, e)
  if (length(spine) > 0L) temporal_reference_node <- spine[length(spine)]
  list(new_edges = rbind(new_edges[keep, ], same),
       temporal_reference_node = temporal_reference_node)
}


#' Draw graph edges
#'
#' draw_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param confounders Vector of variable names, treated as confounders. A list can also be supplied. Order determines the assigned coordinates. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param treatments Treatment variable name.
#' @param outcomes Outcome variable name.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @param latent_variables List or vector of additional or already supplied latent (unobserved) variable names.
#' @param instrumental_variables Inputted instrumental variable names.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param competing_cause_vec Vector of competing cause names. An arrow is drawn connecting competing causes to the outcome, with other arrows also connected depending on type of graph specified.
#' @param collider_vec Vector of collider variables, with both treatment and outcome parents.
#' @param observed Vector of variables without roles, also not unobserved (latent).
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param observed_node_names Vector of node names in the dag.
#' @param existing_dag An existing dagitty object may be supplied.
#' @returns Data frame of edges.
#' @noRd
draw_edges <- function(confounders,
                       treatments,
                       outcomes,
                       mediator_vec,
                       latent_vec,
                       latent_variables,
                       instrumental_variables,
                       m_o_confounder_vec,
                       competing_cause_vec,
                       collider_vec,
                       observed,
                       type,
                       observed_node_names,
                       existing_dag = NA){
  confounder_vec <- unlist(confounders)
  e <- "->"
  if(type == "full"){ e <- "<->" }

  ## confounder edges ##
  confounder_df <- draw_confounder_edges(type = type,
                                         confounders = confounders,
                                         confounder_vec = confounder_vec,
                                         treatments = treatments,
                                         outcomes = outcomes,
                                         m_o_confounder_vec = m_o_confounder_vec,
                                         mediator_vec = mediator_vec,
                                         latent_vec = latent_vec,
                                         e = e)
  ## treatment edges ##
  treatment_df <- draw_treatment_edges(type = type,
                                       confounder_vec = confounder_vec,
                                       treatments = treatments,
                                       outcomes = outcomes,
                                       mediator_vec = mediator_vec,
                                       collider_vec = collider_vec,
                                       e = e)
  ## outcome edges ##
  outcome_df <- draw_outcome_edges(type = type,
                                   outcomes = outcomes,
                                   collider_vec = collider_vec,
                                   e = e)
  ## mediator edges ##
  mediator_df <- draw_mediator_edges(type = type,
                                     outcomes = outcomes,
                                     mediator_vec = mediator_vec,
                                     latent_vec = latent_vec,
                                     e = e)
  ## mediator_outcome_confounder edges ##
  moc_df <- draw_moc_edges(type = type,
                           treatments = treatments,
                           outcomes = outcomes,
                           confounder_vec = confounder_vec,
                           m_o_confounder_vec = m_o_confounder_vec, # new nodes as mediator_outcome_confounder,
                           mediator_vec = mediator_vec,
                           latent_vec = latent_vec,
                           e = e)
  ## competing_cause edges ##
  competing_cause_df <- draw_competing_cause_edges(type = type,
                                                         outcomes = outcomes,
                                                         competing_cause_vec,
                                                         e = e)
  ## connect observed to ancestors and descendants ##
  observed_df <- draw_observed_edges(observed,
                                     existing_dag,
                                     e = e)
  ## instrumental_variables edges ##
  instrumental_df <- draw_iv_edges(type = type,
                                   instrumental_variables,
                                   treatments = treatments,
                                   e = e)
  ## latent_variables edges ##
  latent_df <- draw_latent_edges(observed_node_names,
                                 latent_variables,
                                 type,
                                 outcomes,
                                 treatments,
                                 confounder_vec,
                                 m_o_confounder_vec,
                                 mediator_vec,
                                 e)
  ## row bind edges ##
  edges_df <- rbind(treatment_df,
                    outcome_df,
                    confounder_df,
                    moc_df,
                    mediator_df,
                    instrumental_df,
                    latent_df,
                    competing_cause_df,
                    observed_df,
                    fill=TRUE)  # row bind all edge data frames

  edges_df <- unique(edges_df) # remove duplicate edges
  edges_df <- edges_df[ complete.cases(edges_df), ] # remove NAs
  edges_df <- edges_df[edges_df$ancestor != edges_df$descendant, ] # remove edges with identical ancestor and descendant node names
  names(edges_df) <- c("v", "e", "w")

  return(edges_df)
}


#' Draw confounder edges
#'
#' draw_confounder_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param confounders Vector or list of variable names, treated as confounders.
#' @param confounder_vec Vector of variable names, treated as confounders.
#' @param treatments Treatment variable name.
#' @param outcomes Outcome variable name.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns Data frame of edges.
#' @noRd
draw_confounder_edges <- function(type,
                                  confounders,
                                  confounder_vec,
                                  treatments,
                                  outcomes,
                                  m_o_confounder_vec,
                                  mediator_vec,
                                  latent_vec,
                                  e
                                  ){

  constant <- 0
  if( length(confounders) > 1 ){
    constant <- 1
  }

  if( all( complete.cases(confounder_vec) ) ){

    ## confounder edges ##
    confounder_list <- c()

    ## connect treatments to confounders
    confounder_df <- .cross_edges(confounder_vec, treatments, e)

    ## draw edges between confounders
    if( type == "ordered" ){ # only connects consecutive nodes to each other

      confounder_list <- suppressWarnings( lapply(1:(length(confounders) - constant), function(x){

        lapply(seq_along(confounders[[x]]), function(y){

          confounder_list[[x]] <- sapply(seq_along(confounders[[x+constant]]), function(z){

            list( c( ancestor = confounders[[x]][y], edge = "->", descendant = confounders[[x+constant]][z]) )

          })

        })

      }) )

      confounder_unlist <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
      confounder_unlist <- data.table::as.data.table( do.call( rbind, unlist(confounder_unlist, recursive = FALSE) ) )

      confounder_df <- rbind(confounder_df, confounder_unlist)

      }else if( type == "full" | type == "saturated" ){ # connect all node positions to each other

        len_confounders <- length(confounders)

        confounder_list <- suppressWarnings( lapply(1:(length(confounders) - constant), function(x){

          nodes_after_x <- unlist( confounders[(x+constant):len_confounders]) # nodes occurring after current position (x)

          lapply(seq_along(confounders[[x]]), function(y){

            confounder_list[[x]] <- sapply(seq_along(nodes_after_x), function(z){

              list( c( ancestor = confounders[[x]][y], edge = e, descendant = nodes_after_x[z]) )

            })

          })

        }) )

        confounder_unlist <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_unlist <- data.table::as.data.table( do.call( rbind, unlist(confounder_unlist, recursive = FALSE) ) )

        confounder_df <- rbind(confounder_df, confounder_unlist)

        if( type == "full" ){

          confounder_df <- rbind(confounder_df, .cross_edges(confounder_vec, mediator_vec, e))

          # connect mediator-outcome confounders
          confounder_df <- rbind(confounder_df, .cross_edges(confounder_vec, m_o_confounder_vec, e))

          # connect latent variables
          confounder_df <- rbind(confounder_df, .cross_edges(confounder_vec, latent_vec, e))
        }

      }

    ## connect outcomes to confounders
    outcomes <- unlist(outcomes)

    confounder_df <- rbind(confounder_df, .cross_edges(confounder_vec, outcomes, e))

    ## final checks
    confounder_df <- unique(confounder_df) # remove duplicate edges
    confounder_df <- confounder_df[confounder_df$ancestor != confounder_df$descendant, ] # remove edges with identical ancestor and descendant node names

    return(confounder_df)
  }

  return( data.frame(NULL) )
}


#' Draw treatment edges
#'
#' draw_treatment_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param treatments Treatment variable name.
#' @param confounder_vec Vector of variable names, treated as confounders. A list can also be supplied. Order determines the assigned coordinates. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param collider_vec Vector of collider variables, with both treatment and outcome parents.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns Data frame of edges.
#' @noRd
draw_treatment_edges <- function(type,
                                 confounder_vec,
                                 treatments,
                                 outcomes,
                                 mediator_vec,
                                 collider_vec,
                                 e
                                 ){
  ## treatment edges ##
  treatment_list <- c()

  if( length(outcomes) > 1 ){ # connect all node positions to each other

    treatment_list <- suppressWarnings( lapply(seq_along(treatments), function(x){

      lapply(seq_along(outcomes), function(y){

        treatment_list[[x]] <- sapply(seq_along(outcomes[[y]]), function(z){

          list( c( ancestor = treatments[x], edge = e, descendant = outcomes[[y]][z]) )

        })

      })

    }) )

  }else{

    outcomes <- unlist(outcomes) # in case outcomes are enclosed within a list

    treatment_list <- suppressWarnings( lapply(seq_along(treatments), function(x){

      treatment_list[x] <- lapply(seq_along(outcomes), function(y){

        list( c( ancestor = treatments[x], edge = e, descendant = outcomes[y]) )

      })

    }) )

  }

  treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
  treatment_df <- data.table::as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

  # connect mediators
  treatment_unlist <- .cross_edges(treatments, mediator_vec, e)

  treatment_df <- rbind(treatment_df, treatment_unlist)

  # connect colliders
  treatment_unlist <- .cross_edges(treatments, collider_vec, e)

  treatment_df <- rbind(treatment_df, treatment_unlist)

  # connect all treatments (fully connected or saturated graph type)
  if( length(treatments) > 1 ){

    if( type == "full" & length(treatments) > 1 ){

      treatment_unlist <- .cross_edges(treatments, treatments, e)

      treatment_df <- rbind(treatment_df, treatment_unlist)

    }else if( type == "saturated" ){

      treatment_occurrance <- as.numeric(order(match(treatments, treatments)))

      treatment_list <- suppressWarnings(lapply(seq_along(treatments), function(x){

        treatment_list[x] <- lapply(seq_along(treatments), function(y){

          list( c( ancestor = treatments[x], edge = e, descendant = treatments[y],
                   ancestor_order = treatment_occurrance[x], descendant_order = treatment_occurrance[y] ) )

        })

      }))

      treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
      treatment_order_df <- as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

      treatment_order_df <- treatment_order_df[,c(1:3)][!as.integer(treatment_order_df$ancestor_order) > as.integer(treatment_order_df$descendant_order), ] # remove rows where temporal logic is not followed (orders cast from character)
      treatment_df <- rbind(treatment_df, treatment_order_df)

    }

    ## final checks
    treatment_df <- unique(treatment_df) # remove duplicate edges
    treatment_df <- treatment_df[treatment_df$ancestor != treatment_df$descendant, ] # remove edges with identical ancestor and descendant node names

  }

  return(treatment_df)
}


#' Draw outcome edges
#'
#' draw_outcome_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Vector or list of outcomes.
#' @param collider_vec Vector of collider variables, with both treatment and outcome parents.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns Data frame of edges.
#' @noRd
draw_outcome_edges <- function(type,
                               outcomes,
                               collider_vec,
                               e
                               ){
  constant <- 0
  if( length(outcomes) > 1 ){
    constant <- 1
  }

  ## outcome edges ##
  outcome_list <- c()
  outcome_vec <- unlist(outcomes)
  outcome_df <- data.frame(NULL)

  ## draw edges between outcomes ##
  if( type == "ordered" ){ # only connects consecutive nodes to each other

    # connect all outcomes (fully connected or saturated graph type)
    if( length(outcomes) > 1 ){
      outcome_list <- suppressWarnings( lapply(1:(length(outcomes) - constant), function(x){

        lapply(seq_along(outcomes[[x]]), function(y){

          outcome_list[[x]] <- sapply(seq_along(outcomes[[x+constant]]), function(z){

            list( c( ancestor = outcomes[[x]][y], edge = "->", descendant = outcomes[[x+constant]][z]) )

          })

        })

      }) )

      outcome_list <- Filter(Negate(anyNA), unlist(outcome_list, recursive = FALSE))
      outcome_unlist <- data.table::as.data.table( do.call( rbind, unlist(outcome_list, recursive = FALSE) ) )

      outcome_df <- rbind(outcome_df, outcome_unlist)

    }

    }else if( (type == "full" | type == "saturated") & length(outcomes) > 1 ){ # connect all node positions to each other
      len_outcomes <- length(outcomes)

      outcome_list <- suppressWarnings( lapply(1:(length(outcomes) - constant), function(x){

        nodes_after_x <- unlist( outcomes[(x+constant):len_outcomes]) # nodes occurring after current position (x)

        lapply(seq_along(outcomes[[x]]), function(y){

          outcome_list[[x]] <- sapply(seq_along(nodes_after_x), function(z){

            list( c( ancestor = outcomes[[x]][y], edge = e, descendant = nodes_after_x[z]) )

          })

        })

      }) )

      outcome_list <- Filter(Negate(anyNA), unlist(outcome_list, recursive = FALSE))
      outcome_unlist <- data.table::as.data.table( do.call( rbind, unlist(outcome_list, recursive = FALSE) ) )

      outcome_df <- rbind(outcome_df, outcome_unlist)
    }

  # connect colliders
  if( all( complete.cases(collider_vec) ) ){

    outcomes <- unlist(outcomes)

    outcome_unlist <- .cross_edges(outcomes, collider_vec, e)

    outcome_df <- rbind(outcome_df, outcome_unlist)
  }

  ## final checks ##
  outcome_df <- unique(outcome_df) # remove duplicate edges
  outcome_df <- outcome_df[outcome_df$ancestor != outcome_df$descendant, ] # remove edges with identical ancestor and descendant node names

  return(outcome_df)
}


#' Draw mediator edges
#'
#' draw_mediator_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns Data frame of edges.
#' @noRd
draw_mediator_edges <- function(type,
                                outcomes,
                                mediator_vec,
                                latent_vec,
                                e
                                ){

  if( all( complete.cases(mediator_vec) ) ){
    ## mediator edges ##
    mediator_list <- c()

    ## draw edges between treatment(s) and connect outcome
    if( length(outcomes) > 1 ){ # connect all node positions to each other

      mediator_list <- suppressWarnings( lapply(seq_along(mediator_vec), function(x){

        lapply(seq_along(outcomes), function(y){

          mediator_list[[x]] <- sapply(seq_along(outcomes[[y]]), function(z){

            list( c( ancestor = mediator_vec[x], edge = e, descendant = outcomes[[y]][z]) )

          })

        })

      }) )

    }else{

      outcomes <- unlist(outcomes) # in case outcomes are enclosed within a list

      mediator_list <- suppressWarnings( lapply(seq_along(mediator_vec), function(x){

        mediator_list[x] <- lapply(seq_along(outcomes), function(y){

          list( c( ancestor = mediator_vec[x], edge = e, descendant = outcomes[y]) )

        })

      }) )

    }

    mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
    mediator_df <- data.table::as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

    # connect all mediators if the inputted dag type is "full"
    if(type == "full"){

      mediator_unlist <- .cross_edges(mediator_vec, mediator_vec, e)

      mediator_df <- rbind(mediator_df, mediator_unlist)


      mediator_unlist <- .cross_edges(mediator_vec, latent_vec, e)

      mediator_df <- rbind(mediator_df, mediator_unlist)

    }

    if( type == "ordered" ){

      mediator_occurrance <- as.numeric( order( match( mediator_vec, mediator_vec ) ) )

      mediator_list <- suppressWarnings(lapply(seq_along(mediator_vec), function(x){

        mediator_list[x] <- lapply(seq_along(mediator_vec), function(y){

          list( c( ancestor = mediator_vec[x], edge = e, descendant = mediator_vec[y],
                   ancestor_order = mediator_occurrance[x], descendant_order = mediator_occurrance[y] ) )

        })

      }))

      mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
      mediator_order_df <- as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

      mediator_order_df <- mediator_order_df[,c(1:3)][!as.integer(mediator_order_df$ancestor_order) > as.integer(mediator_order_df$descendant_order), ] # remove rows where temporal logic is not followed (orders cast from character)
      mediator_df <- rbind(mediator_df, mediator_order_df)
    }


    mediator_df <- unique(mediator_df) # remove duplicate edges
    mediator_df <- mediator_df[mediator_df$ancestor != mediator_df$descendant, ] # remove edges with identical ancestor and descendant node names

    return(mediator_df)
  }

  mediator_df <- data.frame(NULL)

  return( mediator_df )
}

#' Draw mediator-outcome confounder edges
#'
#' draw_moc_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns Data frame of edges.
#' @noRd
draw_moc_edges <- function(type,
                           treatments,
                           outcomes,
                           confounder_vec,
                           m_o_confounder_vec,
                           mediator_vec,
                           latent_vec,
                           e
                           ){

  if( all( complete.cases(m_o_confounder_vec) ) ){

    ## mediator_outcome_confounder edges ##
    moc_list <- c()

    ## draw edges between treatment(s) and connect outcome
    if( length(outcomes) > 1 ){ # connect all node positions to each other

      moc_list <- suppressWarnings( lapply(seq_along(m_o_confounder_vec), function(x){

        lapply(seq_along(outcomes), function(y){

          moc_list[[x]] <- sapply(seq_along(outcomes[[y]]), function(z){

            list( c( ancestor = m_o_confounder_vec[x], edge = e, descendant = outcomes[[y]][z]) )

          })

        })

      }) )

    }else{

      outcomes <- unlist(outcomes) # in case outcomes are enclosed within a list

      moc_list <- suppressWarnings( lapply(seq_along(m_o_confounder_vec), function(x){

        moc_list[x] <- lapply(seq_along(outcomes), function(y){

          list( c( ancestor = m_o_confounder_vec[x], edge = e, descendant = outcomes[y]) )

        })

      }) )

    }

    moc_list <- Filter(Negate(anyNA), unlist(moc_list, recursive = FALSE))
    moc_df <- data.table::as.data.table( do.call( rbind, unlist(moc_list, recursive = FALSE) ) )

    # connect treatments if the inputted dag type is "full"
    if(type == "full"){

      moc_unlist <- .cross_edges(treatments, m_o_confounder_vec, e)

      moc_df <- rbind(moc_df, moc_unlist)

      # connect all latents if the inputted dag type is "full"
      moc_unlist <- .cross_edges(m_o_confounder_vec, latent_vec, e)

      moc_df <- rbind(moc_df, moc_unlist)

    }

    # connect mediators
    moc_unlist <- .cross_edges(m_o_confounder_vec, mediator_vec, e)

    moc_df <- rbind(moc_df, moc_unlist)

    return(moc_df)
  }

  moc_df <- data.frame(NULL)

  return( moc_df )
}


#' Draw graph edges
#'
#' draw_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param competing_cause_vec Vector of competing cause names. An arrow is drawn connecting competing causes to the outcome, with other arrows also connected depending on type of graph specified.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns Data frame of edges.
#' @noRd
draw_competing_cause_edges <- function(type,
                                          outcomes,
                                          competing_cause_vec,
                                          e
                                          ){

  ## competing_cause edges ##
  if( all( complete.cases(competing_cause_vec) ) ){

    competing_cause_list <- c()

    ## draw edges between treatment(s) and connect outcome
    if( length(outcomes) > 1 ){ # connect all node positions to each other

      competing_cause_list <- suppressWarnings( lapply(seq_along(competing_cause_vec), function(x){

        lapply(seq_along(outcomes), function(y){

          competing_cause_list[[x]] <- sapply(seq_along(outcomes[[y]]), function(z){

            list( c( ancestor = competing_cause_vec[x], edge = e, descendant = outcomes[[y]][z]) )

          })

        })

      }) )

    }else{

      outcomes <- unlist(outcomes) # in case outcomes are enclosed within a list

      competing_cause_list <- suppressWarnings( lapply(seq_along(competing_cause_vec), function(x){

        competing_cause_list[x] <- lapply(seq_along(outcomes), function(y){

          list( c( ancestor = competing_cause_vec[x], edge = e, descendant = outcomes[y]) )

        })

      }) )

    }

    competing_cause_list <- Filter(Negate(anyNA), unlist(competing_cause_list, recursive = FALSE))
    competing_cause_df <- data.table::as.data.table( do.call( rbind, unlist(competing_cause_list, recursive = FALSE) ) )

    return(competing_cause_df)
  }

  competing_cause_df <- data.frame(NULL)

  return( competing_cause_df )
}


#' Draw graph edges
#'
#' draw_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param observed Vector of variables without roles, also not unobserved (latent).
#' @param existing_dag An existing dagitty object may be supplied.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns Data frame of edges.
#' @noRd
draw_observed_edges <- function(observed,
                                existing_dag = NA,
                                e
                                ){

  if( all( complete.cases(observed) ) ) {

    ## connect observed to ancestors and descendants ##
    observed_list <- list()
    observed_df <- data.frame(NULL)

    # connect descendants
    existing_edges <- dagitty::edges(existing_dag)
    observed_descendants <- .children_edge_order(existing_edges, observed) # children in dagitty's edge-insertion order (feeds .cross_edges row order)

    if( length(observed_descendants) > 0 ){

      observed_df <- .cross_edges(observed, observed_descendants, e)

    }

    # connect ancestors
    observed_ancestors <- .parents_edge_order(existing_edges, observed)

    if( length(observed_ancestors) > 0 ){

      observed_unlist <- .cross_edges(observed_ancestors, observed, e)

      observed_df <- rbind(observed_df, observed_unlist)

    }

    return(observed_df)
  }

  observed_df <- data.frame(NULL)

  return( observed_df )
}


#' Draw instrumental variable edges
#'
#' draw_iv_edges() is a helper function for buildGraph().
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param instrumental_variables Inputted list or vector of instrumental variables.
#' @param treatments Vector of treatments.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns A data frame of instrumental variable edges.
#' @noRd
draw_iv_edges <- function(type,
                          instrumental_variables,
                          treatments,
                          e
                          ){
  empty <- .cross_edges(character(), character(), e)
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

        instrumental_list[x] <- lapply(seq_along( instrumental_variables[[x]][[2]] ), function(y){

          list( c( ancestor = instrumental_variables[[x]][[1]], edge = e, descendant = instrumental_variables[[x]][[2]][[y]] ) )

        } )

      } ) )

      instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
      if( length(instrumental_list) == 0L ) return(empty)
      instrumental_df <- as.data.table( do.call( rbind, unlist(instrumental_list, recursive = FALSE) ) )

    }else{

      instrumental_vec <- as.vector( unlist( lapply( instrumental_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

      if( any( type == "full" | type == "saturated" ) ){

      # connect treatments
      instrumental_list <- suppressWarnings( lapply(seq_along(instrumental_vec), function(x){

        instrumental_list[x] <- lapply(seq_along(treatments), function(y){

          list( c( ancestor = instrumental_vec[x], edge = e, descendant = treatments[y]) )

        })

      }) )

      instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
      if( length(instrumental_list) == 0L ) return(empty)
      instrumental_df <- as.data.table( do.call( rbind, unlist(instrumental_list, recursive = FALSE) ) )

      }else{

        len_trt <- length(treatments)

        # connect treatments
        instrumental_list <- suppressWarnings( lapply(seq_along(instrumental_vec), function(x){

          if( x <= len_trt ){

            instrumental_list[x] <- list( c( ancestor = instrumental_vec[x], edge = e, descendant = treatments[x]) )

            }else{

            instrumental_list[x] <- list( c( ancestor = instrumental_vec[x], edge = e, descendant = treatments[1]) )

            }

          } ) )

        }

        instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
        if( length(instrumental_list) == 0L ) return(empty)
        instrumental_df <- as.data.table( do.call( rbind, instrumental_list ) )

        }

    return( instrumental_df )
  }

  empty
}


#' Draw latent variable edges
#'
#' draw_latent_edges() is a helper function for buildGraph().
#'
#' @importFrom data.table as.data.table
#' @param latent_variables Inputted list or vector of latent variables.
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param treatments Treatment variable name.
#' @param confounder_vec Vector of variable names, treated as confounders. A list can also be supplied. Order determines the assigned coordinates. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param e Edge type e.g. "->" (directed arrow) "<->" (nodes connected in both directions).
#' @returns A data frame of latent variable edges.
#' @noRd
draw_latent_edges <- function(observed_node_names,
                              latent_variables,
                              type,
                              outcomes,
                              treatments,
                              confounder_vec,
                              m_o_confounder_vec,
                              mediator_vec,
                              e
                              ){

  if( all( complete.cases( unlist(latent_variables) ) ) ){

    ## latent edges ##
    latent_list <- c()

    if( length(unlist(latent_variables)) > length(latent_variables) ){

      latent_list <- suppressWarnings( lapply( seq_along(latent_variables), function(x){

        if( length(unlist(latent_variables[x])) > length(latent_variables[x]) ){

          latent_list[x] <- lapply(seq_along( latent_variables[[x]][[2]] ), function(y){

            list( c( ancestor = latent_variables[[x]][[1]], edge = e, descendant = latent_variables[[x]][[2]][[y]] ) )

          } )

        }

      } ) )

      latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
      latent_df <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

      return(latent_df)

    }else if( length(latent_variables) > 0 ){

      if( length(outcomes) > 1 ){ # connect all node positions to each other

        latent_list <- suppressWarnings( lapply(seq_along(latent_variables), function(x){

          lapply(seq_along(outcomes), function(y){

            latent_list[[x]] <- sapply(seq_along(outcomes[[y]]), function(z){

              list( c( ancestor = latent_variables[x], edge = e, descendant = outcomes[[y]][z]) )

            })

          })

        }) )

      }else{

        outcomes <- unlist(outcomes) # in case outcomes are enclosed within a list

        latent_list <- suppressWarnings( lapply(seq_along(latent_variables), function(x){

          latent_list[x] <- lapply(seq_along(outcomes), function(y){

            list( c( ancestor = latent_variables[x], edge = e, descendant = outcomes[y]) )

          })

        }) )

      }

      latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
      latent_df <- data.table::as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

      # connect treatment
      latent_unlist <- .cross_edges(latent_variables, treatments, e)

      latent_df <- rbind(latent_df, latent_unlist)


      # connect all latents to confounders, other variables depending on fully connected or saturated graph type
      if( any( type == "full" | type == "saturated" ) ){
        latent_unlist <- .cross_edges(latent_variables, confounder_vec, e)

        latent_df <- rbind(latent_df, latent_unlist)


        # connect confounders to mediators, confounders to mediator-outcome confounders, and confounders to latents if the inputted dag type is "full"
        if( type == "full" ){

          latent_unlist <- .cross_edges(latent_variables, mediator_vec, e)

          latent_df <- rbind(latent_df, latent_unlist)


          # connect mediator-outcome confounders
          latent_unlist <- .cross_edges(latent_variables, m_o_confounder_vec, e)

          latent_df <- rbind(latent_df, latent_unlist)


          # connect latent variables
          latent_unlist <- .cross_edges(latent_variables, latent_variables, e)

          latent_df <- rbind(latent_df, latent_unlist)

        }

      }else if( type == "ordered" ){

        latent_occurrance <- as.numeric(order(match(latent_variables, latent_variables)))

        latent_list <- suppressWarnings(lapply(seq_along(latent_variables), function(x){

          latent_list[x] <- lapply(seq_along(latent_variables), function(y){

            list( c( ancestor = latent_variables[x], edge = e, descendant = latent_variables[y], ancestor_order = latent_occurrance[x], descendant_order = latent_occurrance[y] ) )

          })

        }))

        latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
        latent_order_df <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

        latent_order_df <- latent_order_df[,c(1:3)][!as.integer(latent_order_df$ancestor_order) > as.integer(latent_order_df$descendant_order), ] # remove rows where temporal logic is not followed (orders cast from character)
        latent_df <- rbind(latent_df, latent_order_df)
      }

      latent_df <- unique(latent_df) # remove duplicate edges
      latent_df <- latent_df[latent_df$ancestor != latent_df$descendant, ] # remove edges with identical ancestor and descendant node names

      latent_df <- latent_df[ !unlist(latent_df[,"ancestor"]) %in% observed_node_names ] # remove latent edges containing ancestor nodes with other roles

      return(latent_df)

    }

  }

  latent_df <- data.frame(NULL)

  return( latent_df )
}




#' Fully connect new nodes to others
#'
#' connect_new_nodes() draws edges between new and existing nodes, in both directions.
#'
#' @importFrom data.table as.data.table
#' @param dag An existing dagitty object.
#' @param new_nodes A vector of new nodes.
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

  ## node_role = NULL with a position: still connect the new nodes to ALL existing
  ## nodes, but direct each connection by the new node's resolved slot in the whole-
  ## graph topological order (the same-role spine generalised to every node). Connect-
  ## to-all is preserved; the slot only sets direction (for type = "full" the edge is
  ## "<->", so direction is moot and position has no visible effect). Skipped when
  ## explicit ancestors/descendants are supplied -- those take precedence.
  if( !is.null(position) && is.null(ancestors) && is.null(descendants) ){
    spine <- names( sort( unlist( dagitty::topologicalOrdering(dag) ) ) )
    ord   <- .resolve_position(position, new_nodes, spine)
    edges <- .connect_in_order(ord, new_nodes, "saturated", e)
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

    stop(paste0("Ancestor or descendant nodes not present in the dag: ",
                paste(unknown_nodes, collapse = ", "),
                ". Please check inputs and try again."))

  }

  if( is.null(ancestors) & is.null(descendants) ){

    descendants <- node_names # default: connect new nodes to every existing node

  }else if( length(ancestors) == 0 & length(descendants) == 0 ){

    stop("'ancestors' and 'descendants' are both empty. Supply at least one connection, or leave both NULL to connect new nodes to all existing nodes. Please check inputs and try again.")

  }

  if( length( descendants ) > 0){
    ## connect new_nodes to descendants
    node_list <- suppressWarnings( lapply(seq_along(new_nodes), function(x){

      node_list[x] <- lapply(seq_along(descendants), function(y){

        list( c( v = new_nodes[x], e = e, w = descendants[y]) )

      })

    }) )

    node_list <- Filter(Negate(anyNA), unlist(node_list, recursive = FALSE))
    edges <- as.data.table( do.call( rbind, unlist(node_list, recursive = FALSE) ) )
  }

  if( length( ancestors ) > 0){
    ## connect new_nodes to ancestors
    node_list <- suppressWarnings( lapply(seq_along(ancestors), function(x){

      node_list[x] <- lapply(seq_along(new_nodes), function(y){

        list( c( v = ancestors[x], e = e, w = new_nodes[y]) )

      })

    }) )

    node_list <- Filter(Negate(anyNA), unlist(node_list, recursive = FALSE))
    node_unlist <- as.data.table( do.call( rbind, unlist(node_list, recursive = FALSE) ) )

    edges <- rbind(edges, node_unlist)
  }

  return( edges )
}
