#' Get latent variable names from list
#'
#' get_latent_vec() is a helper function for buildGraph().
#'
#' @importFrom data.table data.table
#' @param latent_variables Inputted list or vector of latent variables.
#' @returns A vector of latent variable names.
#' @noRd
get_latent_vec <- function(latent_variables){
  if( length(unlist(latent_variables)) > length(latent_variables) ){

    latent_vec <- unlist( do.call(cbind, latent_variables)[1,] )

  }else{

    latent_vec <- as.vector( unlist( lapply( latent_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

  }

  return(latent_vec)
}


#' Directed reachability over an edge table
#'
#' directed_reachable() walks the directed ("->") edges of a v/e/w edge table,
#' returning every node with a directed path to (direction = "to") or from
#' (direction = "from") any node in 'nodes'. Nodes in 'avoiding' block traversal.
#' Bidirected ("<->") and undirected ("--") edges do not transmit ancestry.
#' The seed nodes themselves are excluded from the result.
#'
#' @param edges A data frame or data table of edges with columns v, e, w.
#' @param nodes Character vector of seed node names.
#' @param direction "to" for ancestors of nodes, "from" for descendants.
#' @param avoiding Optional character vector of node names paths may not pass through.
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
    # "--" (pdag/CPDAG) edges have unknown orientation within the Markov equivalence class (Meek 1995), ancestry goes in both directions
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

  return(reached[ !reached %in% nodes ])
}


#' Extract node roles
#'
#' extract_node_roles() is a helper function for get_edges().
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom dagitty edges coordinates dagitty
#' @param dag A dagitty object.
#' @returns Data table of edges containing roles for each ancestor and descendant nodes.
#' @noRd
extract_unique_node_roles <- function(dag, .cache = NULL){
  edges <- .cached_edges(dag, .cache)

  if( nrow( as.data.frame(edges) ) == 0 ){
    edges <- data.frame( v = character(0), e = character(0), w = character(0) ) # edge-less graph: edges() returns no v/e/w columns
  }else{
    edges <- edges[, c("v", "e", "w")]
  }

  # treatment
  treatments <- treatments(dag)

  # outcome
  outcomes <- .outcomes(dag)

  # latent variables
  latent_vars <- unobserved(dag)

  # directed reachability sets (BFS over the edge table; "<->" edges do not transmit ancestry)
  ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to")
  ancestors_of_outcome     <- directed_reachable(edges, outcomes, direction = "to")
  descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from")
  ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments)

  # confounders can be common causes, or ancestors of treatment with a directed path to outcome that avoids treatment (VanderWeele & Shpitser 2013)
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- ancestors_of_treatment[ ancestors_of_treatment %in% ancestors_of_outcome_avoiding_treatment &
                                           !ancestors_of_treatment %in% outcomes &
                                           !ancestors_of_treatment %in% latent_vars ]

  # mediators are nodes on a directed path from treatment to outcome
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  mediators <- descendants_of_treatment[ descendants_of_treatment %in% ancestors_of_outcome &
                                           !descendants_of_treatment %in% outcomes &
                                           !descendants_of_treatment %in% confounders &
                                           !descendants_of_treatment %in% latent_vars ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% c(treatments, confounders) & # remove treatment and confounder nodes
                !moc %in% treatment_parents & # double check by removing parents of treatment
                !moc %in% latent_vars &
                !moc %in% mediators &
                !moc %in% outcomes]

  # proxy
  latent_children <- dagitty::children(dag, latent_vars)

  proxy_b <- suppressWarnings( treatment_parents[ treatment_parents %in% latent_children &
                                  !treatment_parents %in% latent_vars &
                                  !treatment_parents %in% treatments ] )# proxy_b (type b confounder)

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c (type a confounder)
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  confounders <- c(confounders, proxy_b, proxy_c) # add to confounders

  # competing cause (computed after the proxy append, so proxy confounders are excluded)
  competing_cause <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% latent_vars &
                                           !outcome_parents %in% moc &
                                           !outcome_parents %in% outcomes ]

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children &
                                   !outcome_children %in% latent_vars &
                                   !outcome_children %in% outcomes &
                                   !outcome_children %in% mediators &
                                   !outcome_children %in% moc &
                                   !outcome_children %in% confounders ]

  # instrumental variables (van der Zander, Textor & Liskiewicz 2015, via dagitty)
  instrumental_vars <- extract_instrumental_variables(dag = dag,
                                                      treatments = treatments,
                                                      outcomes = outcomes)

  # a single role instrument is (a) not in any prior-computed roles (b) never a descendant of treatment (dagitty surfaces conditional instruments that can coincide with a confounder or a post-treatment effect, neither of which is an instrument here)
  instrumental_vars <- instrumental_vars[ !instrumental_vars %in% treatments &
                                            !instrumental_vars %in% outcomes &
                                            !instrumental_vars %in% confounders &
                                            !instrumental_vars %in% mediators &
                                            !instrumental_vars %in% moc &
                                            !instrumental_vars %in% competing_cause &
                                            !instrumental_vars %in% colliders &
                                            !instrumental_vars %in% latent_vars &
                                            !instrumental_vars %in% descendants_of_treatment ]

  # filter latent variables
  latent_vars <- latent_vars[ !latent_vars %in% treatments &
                              !latent_vars %in% outcomes &
                              !latent_vars %in% colliders ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% competing_cause &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]

  # Undetermined (CPDAG "--" edge has unknown orientation) if their relationship to the treatment or outcome exists only by "--" edge(s) (role is not invariant across the class (Perkovic et al. 2017, "Interpreting and using CPDAGs"))
  # In these cases, and no role above is assigned, it is reported as "undetermined" (computed roles are kept IFF definite and possible reachability profiles coincide: "->" plus "--" in either direction (see directed_reachable(possible = TRUE))).
  pos_ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to",   possible = TRUE)
  pos_ancestors_of_outcome     <- directed_reachable(edges, outcomes,   direction = "to",   possible = TRUE)
  pos_descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from", possible = TRUE)
  pos_ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments, possible = TRUE)

  undetermined <- unique( c( setdiff(pos_ancestors_of_treatment,   ancestors_of_treatment),
                             setdiff(pos_ancestors_of_outcome,     ancestors_of_outcome),
                             setdiff(pos_descendants_of_treatment, descendants_of_treatment),
                             setdiff(pos_ancestors_of_outcome_avoiding_treatment,
                                     ancestors_of_outcome_avoiding_treatment) ) )
  undetermined <- undetermined[ !undetermined %in% c(treatments, outcomes, latent_vars) ]

  # a definite role is assigned only when it holds regardless of "--" orientation
  confounders        <- confounders[        !confounders        %in% undetermined ]
  mediators          <- mediators[          !mediators          %in% undetermined ]
  moc                <- moc[                !moc                %in% undetermined ]
  competing_cause <- competing_cause[ !competing_cause %in% undetermined ]
  colliders          <- colliders[          !colliders          %in% undetermined ]
  instrumental_vars  <- instrumental_vars[  !instrumental_vars  %in% undetermined ]
  observed           <- observed[           !observed           %in% undetermined ]


  if( nrow(edges) == 0L ){ # no rows return here
    return(edges)
  }

  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, "v"] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, "v"] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, "v"] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, "v"] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, "v"] %in% mediators &
                             !edges[, "v"] %in% moc ]  <- "mediator"

  edges$ancestor_iv[edges[, "v"] %in% instrumental_vars] <- "instrument"

  edges$ancestor_competing_cause[edges[, "v"] %in% competing_cause &
                                      !edges[, "v"] %in% proxy_c ] <- "competing_cause"

  edges$ancestor_collider[edges[, "v"] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"

  edges$ancestor_undetermined[edges[, "v"] %in% undetermined] <- "undetermined"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, "w"] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, "w"] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, "w"] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, "w"] %in% moc &
                         !edges[, "w"] %in% mediators]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, "w"] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, "w"] %in% instrumental_vars] <- "instrument"

  edges$descendant_competing_cause[edges[, "w"] %in% competing_cause &
                                        !edges[, "w"] %in% proxy_c ] <- "competing_cause"

  edges$descendant_collider[edges[, "w"] %in% colliders] <- "collider"

  edges$descendant_latent[edges$w %in% latent_vars] <- "latent"

  edges$descendant_observed[edges$w %in% observed] <- "observed"

  edges$descendant_undetermined[edges[, "w"] %in% undetermined] <- "undetermined"

  return(edges)
}


#' Extract node roles
#'
#' extract_node_roles() is a helper function for get_edges().
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom dagitty edges coordinates dagitty
#' @param dag A dagitty object.
#' @returns Data table of edges containing roles for each ancestor and descendant nodes.
#' @noRd
extract_node_roles <- function(dag, .cache = NULL){
  edges <- .cached_edges(dag, .cache)

  if( nrow( as.data.frame(edges) ) == 0 ){
    edges <- data.frame( v = character(0), e = character(0), w = character(0) ) # edge-less graph: edges() returns no v/e/w columns
  }else{
    edges <- edges[, c("v", "e", "w")]
  }

  # treatment
  treatments <- treatments(dag)

  # outcome
  outcomes <- .outcomes(dag)

  # latent variables
  latent_vars <- unobserved(dag)

  # directed reachability sets (BFS over the edge table, "<->" edges do not transmit ancestry)
  ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to")
  ancestors_of_outcome     <- directed_reachable(edges, outcomes, direction = "to")
  descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from")
  ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments)

  # confounders can be common causes, or ancestors of treatment with a directed path to outcome that avoids treatment (VanderWeele & Shpitser 2013)
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- ancestors_of_treatment[ ancestors_of_treatment %in% ancestors_of_outcome_avoiding_treatment ]

  # mediators are nodes on a directed path from treatment to outcome
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  mediators <- descendants_of_treatment[ descendants_of_treatment %in% ancestors_of_outcome ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% confounders & # remove treatment and confounder nodes
                !moc %in% treatment_parents ] # double check by removing parents of treatment

  # proxy

  latent_children <- dagitty::children(dag, latent_vars)

  proxy_b <- suppressWarnings( treatment_parents[ treatment_parents %in% latent_children &
                                                    !treatment_parents %in% latent_vars &
                                                    !treatment_parents %in% treatments ] )# proxy_b (type b confounder)

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c (type c confounder)
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  confounders <- c(confounders, proxy_b, proxy_c) # add to confounders

  # competing cause (computed after the proxy append, so proxy confounders are excluded)
  competing_cause <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% outcomes ]

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children ]

  # instrumental variables (van der Zander, Textor & Liskiewicz 2015, via dagitty)
  instrumental_vars <- extract_instrumental_variables(dag = dag,
                                                      treatments = treatments,
                                                      outcomes = outcomes)

  # filter latent variables
  latent_vars <- latent_vars[  !latent_vars %in% treatments &
                                 !latent_vars %in% outcomes ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% competing_cause &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]


  # Undetermined (CPDAG "--" edge has unknown orientation) if their relationship to the treatment or outcome exists only by "--" edge(s) (role is not invariant across the class (Perkovic et al. 2017, "Interpreting and using CPDAGs"))
  # In these cases, and no role above is assigned, it is reported as "undetermined" (computed roles are kept IFF definite and possible reachability profiles coincide: "->" plus "--" in either direction (see directed_reachable(possible = TRUE))).
  pos_ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to",   possible = TRUE)
  pos_ancestors_of_outcome     <- directed_reachable(edges, outcomes,   direction = "to",   possible = TRUE)
  pos_descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from", possible = TRUE)
  pos_ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments, possible = TRUE)

  undetermined <- unique( c( setdiff(pos_ancestors_of_treatment,   ancestors_of_treatment),
                             setdiff(pos_ancestors_of_outcome,     ancestors_of_outcome),
                             setdiff(pos_descendants_of_treatment, descendants_of_treatment),
                             setdiff(pos_ancestors_of_outcome_avoiding_treatment,
                                     ancestors_of_outcome_avoiding_treatment) ) )
  undetermined <- undetermined[ !undetermined %in% c(treatments, outcomes, latent_vars) ]

  # a definite role is assigned only when it holds regardless of "--" orientation
  confounders        <- confounders[        !confounders        %in% undetermined ]
  mediators          <- mediators[          !mediators          %in% undetermined ]
  moc                <- moc[                !moc                %in% undetermined ]
  competing_cause <- competing_cause[ !competing_cause %in% undetermined ]
  colliders          <- colliders[          !colliders          %in% undetermined ]
  instrumental_vars  <- instrumental_vars[  !instrumental_vars  %in% undetermined ]
  observed           <- observed[           !observed           %in% undetermined ]

  if( nrow(edges) == 0L ){ # no rows return here
    return(edges)
  }

  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, "v"] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, "v"] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, "v"] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, "v"] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, "v"] %in% mediators ]  <- "mediator"

  edges$ancestor_iv[edges[, "v"] %in% instrumental_vars] <- "instrument"

  edges$ancestor_competing_cause[edges[, "v"] %in% competing_cause &
                                      !edges[, "v"] %in% proxy_c ] <- "competing_cause"

  edges$ancestor_collider[edges[, "v"] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"

  edges$ancestor_undetermined[edges[, "v"] %in% undetermined] <- "undetermined"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, "w"] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, "w"] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, "w"] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, "w"] %in% moc ]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, "w"] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, "w"] %in% instrumental_vars] <- "instrument"

  edges$descendant_competing_cause[edges[, "w"] %in% competing_cause &
                                        !edges[, "w"] %in% proxy_c ] <- "competing_cause"

  edges$descendant_collider[edges[, "w"] %in% colliders] <- "collider"

  edges$descendant_latent[edges$w %in% latent_vars] <- "latent"

  edges$descendant_observed[edges$w %in% observed] <- "observed"

  edges$descendant_undetermined[edges[, "w"] %in% undetermined] <- "undetermined"

  return(edges)
}


#' Edges to longer format
#'
#' @importFrom data.table as.data.table data.table
#' @param edges Data table of edges, created using extract_node_roles().
#' @returns Data table of edges in longer format.
#' @noRd
edges_longer <- function(edges){
  if( nrow(edges) == 0 ){
    # return empty dataframe if no edges
    return( data.frame( ancestor = character(0),
                        edge = character(0),
                        descendant = character(0),
                        ancestor_role = character(0),
                        descendant_role = character(0) ) )
  }

  edges$id <- 1:nrow(edges)

  ancestor_cols   <- grep("^ancestor_",   names(edges))
  descendant_cols <- grep("^descendant_", names(edges))
  edges_ancestors <- edges[, c(1:3, ancestor_cols)]
  edges_descendants <- edges[, c(1:3, descendant_cols)]

  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:ncol(edges_ancestors)), idvar = "id",
                                      v.names = "ancestor_role", direction = "long")[,c("v", "e", "w", "ancestor_role", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), ]

  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:ncol(edges_descendants)), idvar = "id",
                                        v.names = "descendant_role", direction = "long")[,c("v", "e", "w", "descendant_role", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), ]

  # pair ancestor and descendant roles by original edge row (id), rather than by position,
  edges <- merge(edges_ancestors, edges_descendants, by = c("id", "v", "e", "w"), all = TRUE)
  edges <- edges[order(edges$id), c("v", "e", "w", "ancestor_role", "descendant_role")]

  names(edges) <- c("ancestor", "edge", "descendant", "ancestor_role", "descendant_role")
  rownames(edges) <- NULL

  return(edges)
}


#' Connect new nodes
#'
#' add_nodes_helper() is a helper function for add_nodes() that outputs a list containing temporal_reference_node, existing_node_names, and a data frame of new node edges.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty topologicalOrdering
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param nodes A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param node_role Role assigned to new nodes, from any of the following: c("confounder", "treatment", "outcome", "mediator", "mediator_outcome_confounder", "instrument", "competing_cause", "collider", "latent", "observed").
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, with directed arrows between confounders in input (temporal) order, forming a complete temporal DAG; the bidirected "<->" saturation of the ESC-DAGs Mapping Stage (Ferguson et al. 2020) corresponds to type = "full". Placement is controlled separately by the \code{position} argument; \code{type} sets structure only (\code{"full"}, \code{"saturated"} or \code{"ordered"}).
#' @param position Optional placement within the relevant temporal order,
#'   expressed with \code{first()}, \code{last()},
#'   \code{before()}, and \code{after()} (see
#'   \code{\link{position_helpers}}).
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings uses dagitty::topologicalOrdering() and selects the first of the inputted node_role (e.g., first confounder) as the temporal point of reference.
#' @returns output_list containing temporal_reference_node, existing_node_names, and a data frame of new_edges.
#' @noRd
add_nodes_helper <- function(dag,
                             nodes,
                             node_role,
                             type = "full",
                             position = NULL,
                             temporal_reference_node = NA
                             ){
  e <- "->"

  if(type == "full"){

    e <- "<->"
  }

  ## get initial dag roles ##
  dag_roles <- .get_roles(dag)

  nodes_ordered <- sort( unlist( dagitty::topologicalOrdering(dag) ) ) # ggdag estimated temporal order of new nodes

  dag_roles <- lapply(dag_roles, function(x) {
    x[order(match(x, names(nodes_ordered)))]
  })

  dag_roles <- lapply(dag_roles, function(x) {
    if( identical( x[complete.cases(x)], logical(0) ) ) NA_character_ else x[complete.cases(x)]
  })

  outcomes  <- dag_roles$outcome[complete.cases(dag_roles$outcome)]
  treatments <- dag_roles$treatment[complete.cases(dag_roles$treatment)]
  confounder_vec <- dag_roles$confounder[complete.cases(dag_roles$confounder)]
  m_o_confounder_vec <- dag_roles$mediator_outcome_confounder[complete.cases(dag_roles$mediator_outcome_confounder)]
  mediator_vec <- dag_roles$mediator[complete.cases(dag_roles$mediator)]
  instrumental_variables <- dag_roles$instrument[complete.cases(dag_roles$instrument)]
  competing_cause_vec <- dag_roles$competing_cause[complete.cases(dag_roles$competing_cause)]
  latent_vec <- dag_roles$latent[complete.cases(dag_roles$latent)]
  collider_vec <- dag_roles$collider[complete.cases(dag_roles$collider)]
  observed <- dag_roles$observed[complete.cases(dag_roles$observed)]

  ## variable names ##
  observed_node_names <- unique( c(confounder_vec, m_o_confounder_vec, mediator_vec, competing_cause_vec, collider_vec, instrumental_variables, observed) )
  observed_node_names <- Filter(Negate(anyNA), observed_node_names)
  observed_node_names <- observed_node_names[ !observed_node_names %in% latent_vec ]

  if( length(node_role) == 1 && !node_role %in% .VALID_NODE_ROLES ){
    stop("Invalid node_role: \"", node_role, "\". Must be one of: ",
         paste(.VALID_NODE_ROLES, collapse = ", "), ".", call. = FALSE)
  }

  if( !is.null(position) && length(node_role) == 1 &&
      !node_role %in% c("confounder", "treatment", "mediator") ){
    stop("The position argument is currently supported for node_role ",
         "\"confounder\", \"treatment\" and \"mediator\" only; it would be ",
         "ignored for \"", node_role, "\". Remove position, or use a ",
         "supported role.", call. = FALSE)
  }

  if( length( node_role) > 1 | length( node_role) == 0 ){

    stop("add_nodes() currently only supports single node_role character inputs.")

  }else if( node_role %in% "confounder" ){
    ## confounder edges ##
    new_edges <- draw_confounder_edges(type = type,
                                       confounders = nodes, # new nodes as confounder
                                       confounder_vec = nodes,
                                       treatments = treatments,
                                       outcomes = outcomes,
                                       m_o_confounder_vec = m_o_confounder_vec,
                                       mediator_vec = mediator_vec,
                                       latent_vec = latent_vec,
                                       e = e)

    .placed <- .place_role_edges(new_edges, nodes, confounder_vec, nodes_ordered,
                                 type, e, position, temporal_reference_node)
    new_edges               <- .placed$new_edges
    temporal_reference_node <- .placed$temporal_reference_node

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% confounder_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "treatment" ){
    ## treatment edges ##
    new_edges <- draw_treatment_edges(type = type,
                                         confounder_vec = confounder_vec,
                                         treatments = nodes, # new nodes as treatment
                                         outcomes = outcomes,
                                         mediator_vec = mediator_vec,
                                         collider_vec = collider_vec,
                                         e = e)

    .placed <- .place_role_edges(new_edges, nodes, treatments, nodes_ordered,
                                 type, e, position, temporal_reference_node)
    new_edges               <- .placed$new_edges
    temporal_reference_node <- .placed$temporal_reference_node

    treatments <- c(treatments, nodes)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% treatments ] ) # existing node names in temporal order

  }else if( node_role %in% "outcome" ){
    ## outcome edges ##
    new_edges <- draw_outcome_edges(type = type,
                                     outcomes = nodes, # new nodes as outcome,
                                     collider_vec,
                                     e = e)

    outcomes  <- c(outcomes, nodes)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% outcomes ] ) # existing node names in temporal order

  }else if( node_role %in% "mediator" ){
    ## mediator edges ##
    new_edges <- draw_mediator_edges(type = type,
                                     outcomes = outcomes,
                                     mediator_vec = nodes, # new nodes as mediator
                                     latent_vec = latent_vec,
                                     e = e)

    .placed <- .place_role_edges(new_edges, nodes, mediator_vec, nodes_ordered,
                                 type, e, position, temporal_reference_node)
    new_edges               <- .placed$new_edges
    temporal_reference_node <- .placed$temporal_reference_node

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% mediator_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "mediator_outcome_confounder" ){

    ## mediator_outcome_confounder edges ##
  new_edges <- draw_moc_edges(type = type,
                             treatments = treatments,
                             outcomes = outcomes,
                             confounder_vec = confounder_vec,
                             m_o_confounder_vec = nodes, # new nodes as mediator_outcome_confounder,
                             mediator_vec = mediator_vec,
                             latent_vec = latent_vec,
                             e = e)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% m_o_confounder_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "instrument" ){
    ## instrumental_variables edges ##
    new_edges <- draw_iv_edges(type = type,
                               instrumental_variables = nodes, # new nodes as instrumental
                               treatments = treatments,
                               e = e)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% instrumental_variables ] ) # existing node names in temporal order

  }else if( node_role %in% "competing_cause" ){
    ## competing_cause edges ##
    new_edges <- draw_competing_cause_edges(outcomes = outcomes,
                                               competing_cause_vec = nodes,
                                               e = e) # new nodes as competing_cause

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% competing_cause_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "collider" ){
    ## outcome edges ##
    new_edges <- draw_outcome_edges(type = type,
                                    outcomes = outcomes,
                                    collider_vec = nodes,
                                    e = e) # new nodes as collider

    # connect colliders
    treatment_list <- list()

    treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

      treatment_list[x] <- lapply(1:length(nodes), function(y){
        list( c( ancestor = treatments[x], edge = e, descendant = nodes[y]) )
      })

    }) )

    treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
    treatment_unlist <- data.table::as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

    new_edges <- rbind(new_edges, treatment_unlist)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% outcomes ] ) # existing node names in temporal order

  }else if( node_role %in% "latent" ){
    ## latent_variables edges ##
    new_edges <- draw_latent_edges(observed_node_names = observed_node_names,
                                   latent_variables = nodes,
                                   type = type,
                                   outcomes = outcomes,
                                   treatments = treatments,
                                   confounder_vec = confounder_vec,
                                   m_o_confounder_vec = m_o_confounder_vec,
                                   mediator_vec = mediator_vec,
                                   e = e) # new nodes as latent

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% latent_vec ] ) # existing node names in temporal order

    latent_vec <- c(latent_vec, nodes)

  }else if( node_role %in% "observed" ){
    if( all(nodes %in% names(dag)) ){
      new_edges <- draw_observed_edges(observed = nodes,
                                       existing_dag = dag,
                                       e = e)
    }else{
      new_edges <- data.table::data.table(
        ancestor = character(), edge = character(), descendant = character()
      )
    }

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% observed ] ) # existing node names in temporal order
  }

  output_list <- list(temporal_reference_node = temporal_reference_node,
                      existing_node_names = existing_node_names,
                      new_edges = new_edges,
                      treatments = treatments,
                      outcomes = outcomes,
                      latent_vec = latent_vec)


  return(output_list)
}


#' Fully connect new nodes to others
#'
#' connect_nodes_helper() is a helper function for connect_nodes() that connects new and existing nodes, drawing edges in both directions.
#'
#' @importFrom data.table as.data.table
#' @param dag An existing dagitty object.
#' @param nodes A vector of new nodes.
#' @noRd
connect_nodes_helper <- function(dag, nodes, dag_node_names, type){
  ## get node names ##
  node_names <- dag_node_names

  ## get initial dag roles ##
  dag_roles <- .get_roles(dag)

  outcomes  <- dag_roles$outcome
  treatments <- dag_roles$treatment
  confounders <- dag_roles$confounder
  m_o_confounder_vec <- dag_roles$mediator_outcome_confounder
  mediator_vec <- dag_roles$mediator
  instrumental_variables <- dag_roles$instrument
  competing_cause_vec <- dag_roles$competing_cause
  latent_vec <- dag_roles$latent
  latent_variables <- latent_vec
  collider_vec <- dag_roles$collider
  observed <- dag_roles$observed

  ## get variable names ##
  observed_node_names <- unique( as.vector( c(confounders, m_o_confounder_vec, mediator_vec, competing_cause_vec, collider_vec, instrumental_variables) ) )
  observed_node_names <- Filter(Negate(anyNA), observed_node_names)

  edges <- draw_edges(confounders,
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
                      existing_dag = dag)

  colnames(edges) <- c("v", "e", "w")

  ## keep candidate edges involving a selected node ##
  edges <- edges[ unlist(edges[, "v"]) %in% nodes | unlist(edges[, "w"]) %in% nodes, ]

  return( edges )
}


#' Group edges by unique node name for console display
#'
#' @importFrom data.table as.data.table
#' @param new_edges Data frame / data table of edges.
#' @returns List of edges, grouped by each unique node.
#' @noRd
print_edges_helper <- function(new_edges){
  if( nrow(new_edges) == 0L ){
    return(list())
  }
  # collapse new_edges to a vector and output grouped by nodes
  unique_ancestors <- unique( new_edges[[1]] )
  num_unique_ancestors <- length(unique_ancestors)
  new_edges_list <- list()

  # new_edges grouped by each unique node in a list
  new_edges_list <- suppressWarnings( lapply(seq_len(num_unique_ancestors), function(x){
    new_edges[ unlist(new_edges[,1]) %in% unlist(unique_ancestors)[x], ]
  }) )

  # nodes containing edges are collapsed, outputted in the console to allow easy assessing by copy and paste into .r file
  new_edges_list <- lapply(seq_len(num_unique_ancestors), function(x){
    new_edges_list <- noquote(
      paste0( paste0( "'", sapply(seq_len(nrow(new_edges_list[[x]])), function(y){
        new_edges_list[x] <- noquote( paste( new_edges_list[[x]][y,], collapse=" "  ) )
        }), "'", collapse=", " ), sep = "")
      )
  })

  return(new_edges_list)
}



.empty_edge_frame <- function(){
  data.frame(v = character(), e = character(), w = character(),
             stringsAsFactors = FALSE)
}

.clean_dag_names <- function(x){
  x <- as.character(unlist(x, use.names = FALSE))
  unique(x[!is.na(x) & nzchar(x)])
}

.quote_dag_names <- function(x, caller = "graph construction"){
  x <- as.character(unlist(x, use.names = FALSE))
  if( length(x) == 0L ) return(character())
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
  paste0('"', x, '"')
}

.strip_dag_quotes <- function(x, caller){
  x <- trimws(x)
  starts <- startsWith(x, '"')
  ends <- endsWith(x, '"')
  if( xor(starts, ends) ){
    stop(caller, ": unmatched quote in node name '", x, "'.", call. = FALSE)
  }
  if( starts && ends ) x <- substr(x, 2L, nchar(x) - 1L)
  if( !nzchar(x) ){
    stop(caller, ": node names cannot be empty.", call. = FALSE)
  }
  if( grepl('"', x, fixed = TRUE) ){
    stop(caller, ": embedded quotes are not supported in node names.", call. = FALSE)
  }
  x
}

.normalise_edge_frame <- function(edges, caller = "edge input"){
  if( is.null(edges) ) return(.empty_edge_frame())
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  if( nrow(edges) == 0L && ncol(edges) == 0L ){
    return(.empty_edge_frame())
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
  edges
}

.dag_edges <- function(dag, caller = "graph"){
  .dag_graph_type(dag)
  .normalise_edge_frame(dagitty::edges(dag), caller)
}

.dag_graph_type <- function(dag){
  model <- as.character(dag)
  graph_type <- sub("^\\s*([[:alpha:]]+).*", "\\1", model)
  if( !graph_type %in% c("dag", "pdag", "mag") ){
    stop("causaliflower does not currently support dagitty graph type '",
         graph_type, "'.", call. = FALSE)
  }
  graph_type
}

.selected_nodes <- function(dag){
  lines <- strsplit(as.character(dag), "\n", fixed = TRUE)[[1L]]
  attribute_lines <- lines[grepl("\\[", lines)]
  attribute_text <- sub("^.*\\[", "", attribute_lines)
  attribute_text <- sub("\\].*$", "", attribute_text)
  selected_lines <- attribute_lines[
    grepl("(^|[,[:space:]])selected([,[:space:]]|$)", attribute_text)
  ]
  if( length(selected_lines) == 0L ){
    return(character())
  }
  node_text <- trimws(sub("\\s*\\[.*$", "", selected_lines))
  unique(vapply(node_text, .strip_dag_quotes, character(1L),
                caller = "graph metadata"))
}

.split_graph_statements <- function(text, caller){
  chars <- strsplit(paste(text, collapse = "\n"), "", fixed = TRUE)[[1L]]
  if( length(chars) == 0L ) return(character())

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
  c(statements, paste0(current, collapse = ""))
}

.find_graph_operators <- function(line, caller){
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
        startsWith(remaining, candidate)
      }, logical(1L))]
      if( length(operator) > 0L ){
        operator <- operator[1L]
        positions <- c(positions, i)
        lengths <- c(lengths, nchar(operator))
        operators <- c(operators, operator)
        i <- i + nchar(operator)
        next
      }
      if( chars[i] %in% c("<", ">") ) invalid_marker <- TRUE
    }
    i <- i + 1L
  }
  if( in_quote ){
    stop(caller, ": unmatched quote in graph text.", call. = FALSE)
  }
  list(position = positions, length = lengths, operator = operators,
       invalid_marker = invalid_marker)
}

.parse_graph_text <- function(text, caller, allow_isolated = FALSE,
                              allow_empty = FALSE){
  if( !is.character(text) ){
    stop(caller, ": edge input must be a character vector.", call. = FALSE)
  }
  if( anyNA(text) ){
    stop(caller, ": edge strings cannot contain NA.", call. = FALSE)
  }

  lines <- .split_graph_statements(text, caller)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if( length(lines) == 0L ){
    if( allow_empty ){
      return(list(edges = .empty_edge_frame(), isolated = character()))
    }
    stop(caller, ": no edges were supplied.", call. = FALSE)
  }

  edge_rows <- list()
  isolated <- character()
  for( line in lines ){
    matches <- .find_graph_operators(line, caller)
    if( isTRUE(matches$invalid_marker) ){ # a stray "<" or ">" outside any recognised operator
      stop(caller, ": could not parse '", line,
           "'. Expected 'node op node'.", call. = FALSE)
    }
    if( length(matches$position) == 0L ){
      if( !allow_isolated ){
        stop(caller, ": could not parse '", line,
             "'. Expected 'node op node'.", call. = FALSE)
      }
      isolated <- c(isolated, .strip_dag_quotes(line, caller))
      next
    }

    operators <- matches$operator
    starts <- c(1L, matches$position + matches$length)
    ends <- c(matches$position - 1L, nchar(line))
    nodes <- vapply(seq_along(starts), function(i){
      .strip_dag_quotes(substr(line, starts[i], ends[i]), caller)
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
    .empty_edge_frame()
  }
  list(edges = edges, isolated = unique(isolated))
}

.parse_edge_strings <- function(edge_strings, caller = "edge input"){
  .parse_graph_text(edge_strings, caller, allow_isolated = FALSE,
                    allow_empty = FALSE)$edges
}

.edge_key <- function(edges){
  edges <- .normalise_edge_frame(edges)
  symmetric <- edges$e %in% c("<->", "--")
  left <- edges$v
  right <- edges$w
  if( any(symmetric) ){
    left[symmetric] <- pmin(edges$v[symmetric], edges$w[symmetric])
    right[symmetric] <- pmax(edges$v[symmetric], edges$w[symmetric])
  }
  paste(left, edges$e, right, sep = "\r")
}

.children_edge_order <- function(edges, nodes){
  edges <- .normalise_edge_frame(edges)
  directed <- edges[edges$e == "->", , drop = FALSE]
  Reduce(function(found, node){
    union(found, directed$w[directed$v == node])
  }, nodes, character())
}

.parents_edge_order <- function(edges, nodes){
  edges <- .normalise_edge_frame(edges)
  directed <- edges[edges$e == "->", , drop = FALSE]
  Reduce(function(found, node){
    union(found, directed$v[directed$w == node])
  }, nodes, character())
}

#' Rebuild dag
#'
#' rebuild_dag() rebuilds a dag using a dagitty object and data frame of edges input.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges coordinates dagitty isAcyclic
#' @param dag A dagitty object.
#' @param edges A vector of edges.
#' @param treatments Optional character vector of exposure names to declare; defaults to the declarations in dag.
#' @param outcomes Optional character vector of outcome names to declare; defaults to the declarations in dag.
#' @param latent_vec Optional character vector of latent names to declare; defaults to the declarations in dag.
#' @param extra_nodes Optional character vector of additional node names to declare beyond names(dag) and the edge endpoints.
#' @returns A dagitty object
#' @noRd
rebuild_dag <- function(dag,
                        edges,
                        treatments = NULL,
                        outcomes = NULL,
                        latent_vec = NULL,
                        extra_nodes = NULL,
                        adjusted_vec = NULL,
                        selected_vec = NULL
                        ){
  if( is.null(latent_vec) ) latent_vec <- unobserved(dag)
  if( is.null(treatments) ) treatments <- treatments(dag)
  if( is.null(outcomes) ) outcomes <- .outcomes(dag)
  if( is.null(adjusted_vec) ) adjusted_vec <- dagitty::adjustedNodes(dag)
  if( is.null(selected_vec) ) selected_vec <- .selected_nodes(dag)
  graph_type <- .dag_graph_type(dag)

  edges <- .normalise_edge_frame(edges, "rebuild_dag()")
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
    dagitty::coordinates(rebuilt) <- list(
      x = coordinates$x[positioned],
      y = coordinates$y[positioned]
    )
  }
  rebuilt
}


#' construct dag
#'
#' construct_graph() constructs a daggity object using edges input.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty
#' @param edges A data frame of edges.
#' @param node_names Vector of existing node names, used as a reference for the new graph nodes, e.g., c("Z1", "Z2", "Z3").
#' @param new_dag A dagitty object. Must include exposure and outcome nodes.
#' @param treatments Treatment/exposures in the supplied dag.
#' @param outcomes Outcomes in the supplied dag.
#' @param latent_vec Character or vector of additional or already supplied latent (unobserved) variable names, e.g. "U" or c("U1", "U2", "M1").
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
  edges <- .normalise_edge_frame(edges, "construct_graph()")
  edges <- edges[complete.cases(edges), , drop = FALSE]

  treatments <- .clean_dag_names(treatments)
  outcomes <- .clean_dag_names(outcomes)
  latent_vec <- .clean_dag_names(latent_vec)
  adjusted_vec <- .clean_dag_names(adjusted_vec)
  selected_vec <- .clean_dag_names(selected_vec)
  all_nodes <- unique(c(.clean_dag_names(node_names), edges$v, edges$w,
                        treatments, outcomes, latent_vec, adjusted_vec,
                        selected_vec))
  declarations <- vapply(all_nodes, function(node){
    attributes <- c(
      if( node %in% treatments ) "exposure",
      if( node %in% outcomes ) "outcome",
      if( node %in% latent_vec ) "latent",
      if( node %in% adjusted_vec ) "adjusted",
      if( node %in% selected_vec ) "selected"
    )
    quoted <- .quote_dag_names(node)
    if( length(attributes) ){
      paste0(quoted, " [", paste(attributes, collapse = ","), "]")
    }else{
      quoted
    }
  }, character(1L), USE.NAMES = FALSE)
  edge_lines <- if( nrow(edges) ){
    paste(.quote_dag_names(edges$v), edges$e, .quote_dag_names(edges$w))
  }else{
    character()
  }

  dagitty::dagitty(paste(c(paste(graph_type, "{"), declarations,
                           edge_lines, "}"),
                         collapse = " "))
}

#' string to dag
#'
#' string_to_dag() converts edge text to a dagitty object.
#'
#' Input may contain newline- or semicolon-separated edges, edge chains, and
#' isolated node declarations. Supported operators are \code{->},
#' \code{<-}, \code{<->}, and \code{--}.
#' @param edges_string Character vector containing edge or node declarations.
#' @returns A dagitty object.
#' @export
string_to_dag <- function(edges_string){
  parsed <- .parse_graph_text(edges_string, "string_to_dag()",
                              allow_isolated = TRUE, allow_empty = TRUE)
  node_names <- unique(c(parsed$edges$v, parsed$edges$w, parsed$isolated))
  construct_graph(parsed$edges, node_names, NULL, NULL, NULL)
}










#' Allocate a structural-read cache
#'
#' Returns an empty environment used to memoize structural reads during one
#' public operation. A cache belongs to one graph.
#'
#' @return An empty environment.
#' @noRd
.new_cache <- function() new.env(parent = emptyenv())


#' Memoized fetch from a structural-read cache
#'
#' Returns a cached value when present; otherwise computes and stores it. A NULL
#' cache calls \code{compute()} directly.
#'
#' @param .cache An environment from .new_cache(), or NULL to bypass the cache.
#' @param key Character key naming the fact (e.g. "edges").
#' @param compute Zero-argument function returning the value when it is not cached.
#' @return The cached or freshly computed value.
#' @noRd
.cache_get <- function(.cache, key, compute) {
  if( is.null(.cache) ){
    return( compute() )
  }
  if( exists(key, envir = .cache, inherits = FALSE) ){
    return( get(key, envir = .cache, inherits = FALSE) )
  }
  val <- compute()
  assign(key, val, envir = .cache)
  val
}


#' dagitty::edges(), memoized per dag
#'
#' Returns the cached edge table, computing it on the first request. Callers
#' must copy the result before any in-place modification.
#'
#' @importFrom dagitty edges
#' @param dag dagitty object.
#' @param .cache An environment from .new_cache(), or NULL to read directly.
#' @return The data frame returned by dagitty::edges(dag).
#' @noRd
.cached_edges <- function(dag, .cache = NULL) {
  .cache_get(.cache, "edges", function() .dag_edges(dag))
}


#' Directed-edge adjacency lookup
#'
#' edges_adjacency() builds per-node parent and child lookups from a single
#' dagitty::edges() call, replacing repeated per-node V8 round-trips. Only
#' directed ("->") edges define parent/child relations; bidirected and
#' undirected edges are not parent relations, per the Pearl/dagitty convention.
#'
#' @importFrom data.table as.data.table
#' @param dag dagitty object
#' @return Named list with elements parents and children, each a named list over names(dag).
#' @noRd
edges_adjacency <- function(dag){
  node_names <- names(dag)
  num_node_names <- length(node_names)

  if( num_node_names == 0 ){

    return( list( parents = list(), children = list() ) )
  }

  edges <- data.table::as.data.table( dagitty::edges(dag) )

  if( nrow(edges) == 0 ){

    parents <- rep( list( character(0) ), num_node_names )
    children <- rep( list( character(0) ), num_node_names )

  }else{

    directed_v <- as.character( edges$v[ edges$e == "->" ] )
    directed_w <- as.character( edges$w[ edges$e == "->" ] )

    parents <- lapply(1:num_node_names, function(x){
      sort( unique( directed_v[ directed_w == node_names[x] ] ) )
    })

    children <- lapply(1:num_node_names, function(x){
      sort( unique( directed_w[ directed_v == node_names[x] ] ) )
    })

  }

  names(parents) <- node_names
  names(children) <- node_names

  return( list( parents = parents, children = children ) )
}


#' Longest-path topological rank of a DAG
#'
#' Assigns each node an integer rank: roots (no parents) get rank 0, and every
#' other node is one greater than the maximum rank of its parents. The result is
#' a layered left-to-right ordering in which no directed edge points to an
#' equal-or-lower rank -- the x ordering for the coordinate engine (rank -> x).
#' Consumes a parent list as returned by edges_adjacency()$parents; the rank is
#' invariant to the order of each node's parents.
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
    return( integer( 0L ) )
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

  ## Kahn's algorithm: relax rank along a topological order to obtain the longest
  ## path. A node is enqueued only once all of its parents have been processed, so
  ## by the time it is dequeued its rank is final.
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
    stop( "The graph contains a directed cycle, so a topological rank is undefined; a directed acyclic graph is required." )
  }

  return( rank )
}


#' Space points along one axis to a minimum gap, preserving a given order
#'
#' Returns the nearest (least-squares) sequence that is monotone decreasing with
#' consecutive values at least gap apart, WITHOUT re-sorting -- the externally
#' decided order is kept and only the coordinates are nudged apart. Isotonic
#' regression (PAVA, via stats::isoreg) on the negated sequence with the
#' mandatory gap offsets removed.
#'
#' @param values Numeric vector in the intended order (first = top); names are preserved.
#' @param gap Minimum separation between consecutive values.
#' @return Numeric vector the length of values, monotone decreasing, gap-separated.
#' @noRd
.space_ordered <- function( values, gap ){
  n <- length( values )
  if ( n <= 1L ){
    return( values )
  }
  off <- ( seq_len( n ) - 1L ) * gap
  fit <- stats::isoreg( ( -values ) - off )$yf   # non-decreasing fit to -values, in the given order
  out <- -( fit + off )
  names( out ) <- names( values )
  return( out )
}


#' Arrow-shaped within-rank y spacing with crossing reduction
#'
#' Gives the layout an arrow silhouette: the root bulk stacked wide at x-min,
#' every pathway sloping inward as it flows right, converging on the outcome at
#' the tip. Four stages: (1) within-rank top-to-bottom order decided by
#' barycentric sweeps (a down pass orders a rank by the mean position of its
#' parents, an up pass by that of its children), which reduces edge crossings;
#' (2) each root takes a lane between -1 and 1 in that order, every other node
#' the mean of its parents' lanes; (3) each rank's half-height starts from a
#' straight taper anchored on the widest rank, floored at the space its own node
#' count needs, made non-increasing from root to tip, and the roots lifted a
#' further fixed factor; (4) a node sits at lane * half-height blended by shape
#' with an even fill of the rank's band, and .space_ordered() then holds the
#' order and enforces the minimum gap. x is untouched.
#'
#' @param coordinates List with named numeric x and y. The input y is used only
#'   to seed the within-rank order before the sweeps.
#' @param parents Named list of parent vectors, as produced by edges_adjacency().
#' @param y_step Minimum vertical gap between adjacent nodes within a rank.
#' @param shape Within-rank placement blend from 0 to 1: 0 hugs parents' lanes,
#'   1 fills each rank's band evenly.
#' @param sweeps Number of barycentric ordering passes (alternating down/up).
#' @param node_rank Optional precomputed topological_rank(parents) result.
#' @return The coordinates list with y re-spaced.
#' @noRd
respace_y_taper <- function( coordinates,
                             parents,
                             y_step,
                             shape = 0.5,
                             sweeps = 4,
                             node_rank = NULL
                             ){

  shape       <- max( 0, min( 1, shape ) )   # a blend: outside 0..1 would extrapolate past the band
  if( is.null( node_rank ) ){
    node_rank <- topological_rank( parents )
  }
  y           <- coordinates$y
  node_names  <- names( y )
  ranks       <- node_rank[ node_names ]
  if( anyNA( ranks ) ){
    stop( "respace_y_taper(): every node in the coordinates needs a rank." )
  }
  r_max       <- max( ranks )
  rank_levels <- sort( unique( ranks ) )

  ## directed adjacency restricted to the laid-out nodes (parents in, children derived).
  parents_in <- lapply( parents[ node_names ], function( ps ) ps[ ps %in% node_names ] )
  names( parents_in ) <- node_names
  children <- vector( "list", length( node_names ) )
  names( children ) <- node_names
  for ( nd in node_names ){
    for ( p in parents_in[[ nd ]] ){
      children[[ p ]] <- c( children[[ p ]], nd )
    }
  }

  ## 1. within-rank order (top -> bottom): seed by the incoming y, then refine by
  ##    barycentric sweeps -- a down pass orders each rank by the mean normalised
  ##    position of its parents, an up pass by that of its children.
  layer <- lapply( rank_levels, function( r ){
    column <- node_names[ ranks == r ]
    column[ order( -y[ column ] ) ]
  } )
  names( layer ) <- as.character( rank_levels )

  norm_pos <- function( column ){
    k <- length( column )
    if ( k == 1L ) 0.5 else ( seq_len( k ) - 1L ) / ( k - 1L )
  }

  reorder_layer <- function( layer, neighbours, rank_keys ){
    pos <- rep( 0.5, length( node_names ) )
    names( pos ) <- node_names
    for ( key in names( layer ) ){
      pos[ layer[[ key ]] ] <- norm_pos( layer[[ key ]] )
    }
    for ( r in rank_keys ){
      column <- layer[[ as.character( r ) ]]
      if ( length( column ) >= 2L ){
        key <- vapply( column, function( nd ){
          v <- neighbours[[ nd ]]
          if ( length( v ) == 0L ) NA_real_ else mean( pos[ v ] )
        }, numeric( 1 ) )
        key[ is.na( key ) ] <- pos[ column ][ is.na( key ) ]   # no neighbours: hold position
        column <- column[ order( key, seq_along( column ) ) ]
        layer[[ as.character( r ) ]] <- column
      }
      pos[ column ] <- norm_pos( column )   # next rank in this sweep sees the updated positions
    }
    return( layer )
  }

  down_ranks <- rank_levels[ rank_levels > min( rank_levels ) ]
  up_ranks   <- rev( rank_levels[ rank_levels < max( rank_levels ) ] )
  for ( s in seq_len( sweeps ) ){
    if ( s %% 2L == 1L ){
      layer <- reorder_layer( layer, parents_in, down_ranks )
    } else {
      layer <- reorder_layer( layer, children, up_ranks )
    }
  }

  ## 2. barycentric lanes in [-1, 1]: each root takes a lane by its
  ##    (crossing-reduced) position, every other node the mean of its parents'
  ##    lanes. Placing at lane * half-height (stage 4) ties y to the taper band.
  lane <- rep( NA_real_, length( node_names ) )
  names( lane ) <- node_names

  for ( r in rank_levels ){
    column <- layer[[ as.character( r ) ]]
    k <- length( column )
    if ( r == 0L ){
      lane[ column ] <- if ( k == 1L ) 0 else 1 - ( seq_len( k ) - 1L ) / ( k - 1L ) * 2
    } else {
      for ( nd in column ){
        node_parents <- parents_in[[ nd ]]
        node_parents <- node_parents[ node_parents %in% node_names ]
        lane[ nd ] <- if ( length( node_parents ) > 0L ) mean( lane[ node_parents ] ) else 0
      }
    }
  }

  ## 3. taper envelope: the roots the widest, converging to the outcome at the
  ##    tip, and no rank compressed below the gap its own node count needs.
  root_fill <- 1.2   # roots sit this factor wider than the widest inner rank
  counts    <- vapply( rank_levels, function( r ) sum( ranks == r, na.rm = TRUE ), integer( 1L ) )
  half_need <- ( counts - 1L ) / 2 * y_step
  h_max     <- max( half_need )

  if ( r_max == 0L ){
    half_height <- half_need
  } else {
    half_height <- pmax( h_max * ( 1 - rank_levels / r_max ), half_need )   # straight taper, never cramped
    for ( i in rev( seq_len( length( half_height ) - 1L ) ) ){              # make non-increasing root -> tip
      if ( half_height[ i ] < half_height[ i + 1L ] ){
        half_height[ i ] <- half_height[ i + 1L ]
      }
    }
    half_height[ 1L ] <- half_height[ 1L ] * root_fill
  }
  names( half_height ) <- as.character( rank_levels )

  ## 4. place each rank in its crossing-reduced order: the lane scaled by the
  ##    rank's half-height, blended by shape with an even fill of the band, then
  ##    hold the order and the minimum gap.
  for ( r in rank_levels ){
    column <- layer[[ as.character( r ) ]]
    k <- length( column )
    band <- half_height[[ as.character( r ) ]]

    even <- if ( k == 1L ) 0 else band - ( seq_len( k ) - 1L ) / ( k - 1L ) * 2 * band
    proximity <- unname( lane[ column ] ) * band

    column_y <- ( 1 - shape ) * proximity + shape * even
    names( column_y ) <- column
    y[ column ] <- unname( .space_ordered( column_y, y_step ) )
  }

  ## 5. root-dominance guard on the realised positions: if a crowded inner rank
  ##    reaches past the roots, push the roots out so the ancestors still read as
  ##    the widest part of the arrow. Scaling roots outward only widens their
  ##    gaps, so the overlap floor is preserved.
  root_names  <- layer[[ as.character( rank_levels[ 1L ] ) ]]
  inner_names <- node_names[ !node_names %in% root_names ]
  if ( length( inner_names ) > 0L && length( root_names ) > 1L ){
    inner_reach <- max( abs( y[ inner_names ] ) )
    root_reach  <- max( abs( y[ root_names ] ) )
    if ( root_reach > 0 && root_reach < inner_reach * root_fill ){
      y[ root_names ] <- y[ root_names ] * ( inner_reach * root_fill / root_reach )
    }
  }

  coordinates$y <- y
  return( coordinates )
}


#' Feedback arc set of a directed edge list (DFS back-edges)
#'
#' Returns the indices of the directed edges that must be removed to break every
#' directed cycle -- found by a deterministic depth-first search: an edge u -> v
#' is a back-edge when v is still on the current DFS stack as u is explored.
#' Removing exactly the DFS back-edges always yields an acyclic graph. The search
#' is iterative and visits nodes in names order and each node's out-edges in
#' edge-list order, so the arc set is stable across runs. Self-loops count as
#' back-edges and are removed too.
#'
#' @param node_names Character vector of every node name (the DFS start set).
#' @param from,to Equal-length character vectors of directed-edge tails and heads.
#' @return Integer vector of positions in (from, to) that are back-edges.
#' @noRd
.feedback_arcs <- function( node_names, from, to ){
  n <- length( node_names )
  if( n == 0L || length( from ) == 0L ){ return( integer(0) ) }
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


#' Layout parent list, with cycles isolated
#'
#' Builds the directed parent list used by the coordinate engine. When the
#' directed part of the graph contains a cycle, the cycle-closing edges are set
#' aside so the remainder can be ranked, and the user is warned; the edges
#' themselves are never removed from the graph -- only the layout ignores them.
#'
#' @param dag A dagitty object.
#' @return Named list of parent vectors, acyclic by construction.
#' @noRd
.layout_parents <- function(dag){

  node_names <- names(dag)
  edges <- data.table::as.data.table(.dag_edges(dag, "renew_coords()"))
  directed <- edges[ edges$e == "->", ]
  from <- as.character( directed$v )
  to   <- as.character( directed$w )

  back <- .feedback_arcs( node_names, from, to )
  if( length(back) > 0L ){
    warning("Coordinates were generated using an acyclic subset of the supplied (cyclic) graph.",
            call. = FALSE)
    from <- from[ -back ]
    to   <- to[ -back ]
  }

  parents <- lapply(node_names, function(nm){
    sort( unique( from[ to == nm ] ) )
  })
  names(parents) <- node_names

  return( parents )
}













#' Label DAG nodes helper function
#'
#' @param names_list A list of node names.
#' @return A list of shortened node names.
#' @noRd
get_label_helper <- function(names_list){

  # mostly hard coded (for now) but it should give decent results for common agricultural inputs e.g. location = loc, temperature = temp
  initials_list <- lapply( 1:length(names_list), function(x){
    if(length(names_list[[x]]) == 1){ # single word node labels

      if( nchar(names_list[[x]]) %% 2 == 0 ){ # if even char length
        if(nchar(names_list[[x]]) > 7){
          substr(names_list[[x]], 1, 3)
        }else if (nchar(names_list[[x]]) == 4){
          substr(names_list[[x]], 1, 4)
        }else{
          substr(names_list[[x]], 1, 2)
        }
      }else{ # if odd char length
        if(nchar(names_list[[x]]) > 6){
          substr(names_list[[x]], 1, 4)
        }else if (nchar(names_list[[x]]) == 5){
          substr(names_list[[x]], 1, 2)
        }else{
          substr(names_list[[x]], 1, 1)
        }
      }

    }else if(length(names_list[[x]]) == 2){ # two word labels

      if(nchar(names_list[[x]][[1]]) == 1){
        paste( c(names_list[[x]][[1]], substr(names_list[[x]][[2]], 1, 2) ), collapse = "_" )
      }else if(nchar(names_list[[x]][[1]]) < 5){
        if(nchar(names_list[[x]][[2]]) > 6){
          paste( c(names_list[[x]][[1]], substr(names_list[[x]][[2]], 1, 1) ), collapse = "_" )
        }else{
          paste( c(names_list[[x]][[1]], substr(names_list[[x]][[2]], 1, 4) ), collapse = "_" )
        }
      }else if(nchar(names_list[[x]][[2]]) < 5){
        if(nchar(names_list[[x]][[1]]) > 6){
          paste( c(substr(names_list[[x]][[1]], 1, 1), names_list[[x]][[2]] ), collapse = "_" )
        }else{
          paste( c(substr(names_list[[x]][[1]], 1, 4), names_list[[x]][[2]] ), collapse = "_" )
        }
      }else{
        paste( substr(names_list[[x]], 1,3), collapse = "_" )
      }

    }else{ # three or greater word labels

      if(any(nchar(names_list[[x]]) > 6)){

        if(all(nchar(names_list[[x]]) > 11)){
          paste( substr(names_list[[x]], 1, 1), collapse = "" )
        }else if(any(nchar(names_list[[x]]) > 8)){
          paste( substr(names_list[[x]], 1, 3), collapse = "_" )
        }else{
          paste( substr(names_list[[x]], 1, 2), collapse = "_" )
        }

      }else if(all(nchar(names_list[[x]]) > 4)){
        paste( substr(names_list[[x]], 1,2), collapse = "" )
      }else{
        paste( substr(names_list[[x]], 1, 1), collapse = "" )
      }

    }
  })

  return(initials_list)
}


#' extract_instrumental_variables() is a helper to instrumental_variables()
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
                        error = function(e) list())
        unlist( lapply(ivs, function(iv) iv$I) )
      })

      per_treatment <- unique( unlist(per_treatment) )
      if( usable_instruments ){
        per_treatment <- per_treatment[ !per_treatment %in% dagitty::descendants(dag, treatments[x]) ] # drop post-treatment candidates
      }
      per_treatment

    }) )
  }
  instrumental_vars <- as.character(unique(unlist(instrumental_vars)))
  if( is.null(instrumental_vars) ) instrumental_vars <- character()

  return(instrumental_vars)
}


#' Instrumental variables per treatment and outcome pair
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
                      error = function(e) list())
      lapply(ivs, function(iv){
        c( treatment        = treatments[x],
           outcome          = outcomes[y],
           instrument       = iv$I,
           conditioning_set = paste0(iv$Z, collapse = ", ") )
      })
    })

    per_treatment <- unlist(per_outcome, recursive = FALSE)
    if( usable_instruments ){
      post_treatment <- dagitty::descendants(dag, treatments[x]) # drop post-treatment candidates
      per_treatment  <- Filter(function(iv) !iv[["instrument"]] %in% post_treatment, per_treatment)
    }
    per_treatment

  }) )

  instrumental_variables <- Filter(Negate(anyNA), unlist(instrumental_variables, recursive = FALSE))

  if( length(instrumental_variables) == 0 ){
    return( data.table::data.table( treatment        = character(0),
                                    outcome          = character(0),
                                    instrument       = character(0),
                                    conditioning_set = character(0) ) )
  }

  instrumental_variables <- unique( data.table::as.data.table( do.call( rbind, instrumental_variables ) ) )
  instrumental_variables <- data.table::setorder( instrumental_variables, treatment, outcome, instrument )

  return(instrumental_variables)
}


#' node names in path from treatment to outcome
#'
#' nodes_between_treatment_and_outcome() is a ggdag::dag_paths() wrapper intended for use with multiple treatments and outcomes.
#'
#' @importFrom ggdag dag_paths
#' @importFrom data.table as.data.table data.table
#' @param dag A dagitty object
#' @param treatments A vector of treatment node names.
#' @param outcomes A vector of outcome node names.
#' @param output_list TRUE or FALSE to output a list (default FALSE returns a vector).
#' @param directed When TRUE (default), restricts to directed paths from treatment to outcome; when FALSE, considers all paths.
#' @returns Vector or list of  nodes in the path from treatment to outcome.
#' @noRd
nodes_between_treatment_and_outcome <- function(dag,
                                                treatments,
                                                outcomes,
                                                output_list = FALSE,
                                                directed = TRUE
                                                ){
  nodes <- c()
  treatments <- unlist(treatments)
  outcomes <- unlist(outcomes)

  paths_trt_to_y <- list()

  if( length(treatments) > 0 & length(outcomes) > 0 ){

    paths_trt_to_y <- lapply(1:length(treatments), function(x){

      paths_trt_to_y[[x]] <- lapply(1:length(outcomes), function(y){

        tryCatch({
          paths_trt_y <- NA
          paths_trt_y <- ggdag::dag_paths(dag,
                                          from = treatments[x],
                                          to = outcomes[y],
                                          directed = directed,
                                          paths_only = TRUE)[["data"]]
          if(length(paths_trt_y) > 0 ){
            paths_trt_y <- as.vector(unlist(unique( paths_trt_y[ complete.cases(paths_trt_y["direction"]), "name" ])))
            paths_trt_y <- paths_trt_y[ !paths_trt_y %in% treatments]
          }
        }, error = function(e){
          paths_trt_y <- NA
        })
      })
      names(paths_trt_to_y[[x]]) <- outcomes
      paths_trt_to_y[[x]]
    })

    names(paths_trt_to_y) <- treatments

    if(output_list == TRUE){

      return(paths_trt_to_y)
    }

    paths_trt_to_y <- as.vector(unique(unlist(paths_trt_to_y)))

    paths_trt_to_y <- paths_trt_to_y[ complete.cases(paths_trt_to_y) ]
  }

  return(paths_trt_to_y)
}


#' Position helpers for placing nodes
#'
#' These helpers specify where selected nodes sit in a temporal order when
#' passed to the \code{position} argument of \code{\link{add_nodes}} or
#' \code{\link{connect_nodes}}. Combine clauses with \code{c()}.
#'
#' A helper with \code{nodes} applies only to those names. One helper may omit
#' \code{nodes} to act as a catch-all for the remainder. Unclaimed nodes are
#' appended when no catch-all is supplied.
#'
#' Note that attaching causaliflower masks \code{first()} and \code{last()}
#' from packages loaded earlier (e.g. data.table, dplyr); call those with
#' \code{::} when both are attached.
#'
#' @param nodes Optional character vector of node names the clause applies
#'   to. Omit (the default) to make the clause the catch-all for all otherwise
#'   unplaced new nodes.
#' @param anchor For \code{after()} or \code{before()}, the
#'   existing node to place relative to. Omitting it selects the corresponding
#'   absolute end.
#'
#' @return An object of class \code{cf_position}: a list of one placement
#'   clause. \code{c()}-ing several helpers concatenates their clauses into a
#'   single \code{cf_position} (see \code{\link{c.cf_position}}).
#'
#' @examples
#' # all new nodes before everything in their role
#' first()
#' # Z3 first; Z4 and Z5 last; all remaining nodes immediately after Z2
#' c(first("Z3"), last(c("Z4", "Z5")),
#'   after("Z2"))
#'
#' @name position_helpers
NULL

# Internal clause constructor. Wraps a single clause in a length-1 cf_position
# (a list of clauses) so that a lone helper and a c(...) of helpers share one
# representation: the resolver always receives a list of clauses.
.cf_clause <- function(kind, anchor, nodes) {
  if (!is.null(nodes))  nodes  <- as.character(nodes)
  if (!is.null(anchor)) anchor <- as.character(anchor)[1]
  structure(list(list(kind = kind, anchor = anchor, nodes = nodes)),
            class = "cf_position")
}

#' @rdname position_helpers
#' @export
first <- function(nodes = NULL) {
  .cf_clause("first", anchor = NULL, nodes = nodes)
}

#' @rdname position_helpers
#' @export
last <- function(nodes = NULL) {
  .cf_clause("last", anchor = NULL, nodes = nodes)
}

#' @rdname position_helpers
#' @export
before <- function(anchor = NULL, nodes = NULL) {
  if (is.null(anchor)) .cf_clause("first", anchor = NULL, nodes = nodes)
  else                 .cf_clause("before", anchor = anchor, nodes = nodes)
}

#' @rdname position_helpers
#' @export
after <- function(anchor = NULL, nodes = NULL) {
  if (is.null(anchor)) .cf_clause("last", anchor = NULL, nodes = nodes)
  else                 .cf_clause("after", anchor = anchor, nodes = nodes)
}

#' Combine position clauses
#'
#' Concatenates several position helpers into one \code{cf_position}. Defining
#' \code{c()} for this class keeps the clauses intact.
#'
#' @param ... \code{cf_position} objects returned by the position helpers.
#' @return A single \code{cf_position} holding all the supplied clauses, in
#'   order.
#' @seealso \code{\link{position_helpers}}
#' @export
c.cf_position <- function(...) {
  structure(do.call("c", lapply(list(...), unclass)), class = "cf_position")
}

# Resolve position clauses against an existing node order.
.resolve_position <- function(position, new_nodes, spine) {
  clauses <- unclass(position)

  is_empty <- vapply(clauses, function(cl) is.null(cl$nodes), logical(1))
  if (sum(is_empty) > 1L) {
    stop("position: at most one empty (catch-all) clause is allowed; ",
         "only one position helper may omit its nodes argument.",
         call. = FALSE)
  }

  claimed <- unlist(lapply(clauses, function(cl) cl$nodes), use.names = FALSE)
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
      if( is.null(prior) ) prior <- 0L
      order <- append(order, cl$nodes,
                      after = match(cl$anchor, order) + prior)
      after_counts[[cl$anchor]] <- prior + length(cl$nodes)
    }
  }

  firsts <- unlist(lapply(clauses, function(cl) if (cl$kind == "first") cl$nodes), use.names = FALSE)
  lasts  <- unlist(lapply(clauses, function(cl) if (cl$kind == "last")  cl$nodes), use.names = FALSE)
  c(firsts, order, lasts)
}
