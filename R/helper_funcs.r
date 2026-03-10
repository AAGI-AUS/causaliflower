#' Get latent variable names from list
#'
#' get_latent_vec() is a helper function for buildGraph().
#'
#' @importFrom data.table data.table
#' @param latent_variables Inputted list or vector of latent variables.
#' @returns A vector of latent variable names.
#' @noRd
get_latent_vec <- function(latent_variables){
  .datatable.aware <- TRUE
  if( length(unlist(latent_variables)) > length(latent_variables) ){

    latent_vec <- unlist( do.call(cbind, latent_variables)[1,] )

  }else{

    latent_vec <- as.vector( unlist( lapply( latent_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

  }

  return(latent_vec)

}


#' extract_instrumental_variables() is a helper to instrumental_variables()
#'
#' @importFrom dagitty parents children
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @returns Vector of instrumental variable names.
#' @noRd
extract_instrumental_variables <- function(dag, treatments, outcomes, latent_vars, colliders, outcome_parents, nodes_trt_to_y){
  .datatable.aware <- TRUE
  ## dag node names
  node_names <- names(dag)

  outcome_parents_children <- dagitty::children(dag, outcome_parents)

  outcome_latent_parents <- outcome_parents[ outcome_parents %in% latent_vars ]

  outcome_latent_parents_children <- dagitty::children(dag, outcome_latent_parents)

  outcome_latent_parents_parents <- dagitty::parents(dag, outcome_latent_parents)

  nodes_trt_to_y_parents <- dagitty::parents(dag, nodes_trt_to_y)

  parents_of_nodes_trt_to_y_parents <- dagitty::parents(dag, nodes_trt_to_y_parents)

  nodes_trt_to_y_children <- dagitty::children(dag, nodes_trt_to_y)

  outcome_latent_parents <- dagitty::parents(dag, latent_vars)

  associated_with_treatment <- list()

  associated_with_treatment <- lapply(1:length(treatments), function(x){

    associated_with_treatment[[x]] <- lapply(1:length(node_names), function(n){

      parents_associated_with_treatment <- unique(dagitty::paths(dag, treatments[x], node_names[n], directed = TRUE)$open == TRUE)
      children_associated_with_treatment <- unique(dagitty::paths(dag, node_names[n], treatments[x], directed = TRUE)$open == TRUE)

      associated_with_treatment <- c(parents_associated_with_treatment, children_associated_with_treatment)

      children_associated_with_outcome <- length(unique(dagitty::paths(dag, outcomes[x], node_names[n], directed = TRUE)$open == TRUE))>0

      associated_with_treatment <- length(associated_with_treatment)>0 & children_associated_with_outcome==FALSE

    })

    names(associated_with_treatment[[x]]) <- node_names
    associated_with_treatment[[x]] <- Filter(function(m) m==TRUE,
                                             associated_with_treatment[[x]] )

    associated_with_treatment[[x]] <- names(associated_with_treatment[[x]])

  })

  associated_with_treatment <- unique(unlist(associated_with_treatment))

  instrumental_vars <- associated_with_treatment[  !associated_with_treatment %in% latent_vars &
                                                     !associated_with_treatment %in% colliders  &
                                                     !associated_with_treatment %in% outcomes  &
                                                     !associated_with_treatment %in% treatments &
                                                     !associated_with_treatment %in% outcome_parents &
                                                     !associated_with_treatment %in% outcome_parents_children &
                                                     !associated_with_treatment %in% outcome_latent_parents_children &
                                                     !associated_with_treatment %in% outcome_latent_parents_parents &
                                                     !associated_with_treatment %in% nodes_trt_to_y_children &
                                                     !associated_with_treatment %in% nodes_trt_to_y_parents &
                                                     !associated_with_treatment %in% parents_of_nodes_trt_to_y_parents &
                                                     !associated_with_treatment %in% nodes_trt_to_y ]

  return(instrumental_vars)
}


#' Extract node roles
#'
#' extract_node_roles() is a helper function for get_edges().
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty
#' @param dag A dagitty object.
#' @returns Data table of edges containing roles for each ancestor and descendant nodes.
#' @noRd
extract_unique_node_roles <- function(dag){
  .datatable.aware <- TRUE

  dag <- pdag_to_dag(dag)
  edges <- dagitty::edges(dag)[, c("v", "e", "w")]

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) &
                                     !treatment_parents %in% treatments &
                                     !treatment_parents %in% latent_vars ]

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  nodes_trt_to_y <- nodes_between_treatment_and_outcome(dag = dag,
                                                            treatments = treatments,
                                                            outcomes = outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) &
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% outcomes &
                                !outcome_parents %in% confounders &
                                !outcome_parents %in% latent_vars ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% c(treatments, confounders) & # remove treatment and confounder nodes
                !moc %in% treatment_parents & # double check by removing parents of treatment
                !moc %in% latent_vars &
                !moc %in% mediators &
                !moc %in% outcomes]

  # competing exposure
  competing_exposure <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% latent_vars &
                                           !outcome_parents %in% moc &
                                           #!outcome_parents %in% nodes_trt_to_y  &
                                           !outcome_parents %in% outcomes ]

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

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children &
                                   !outcome_children %in% latent_vars &
                                   !outcome_children %in% outcomes &
                                   !outcome_children %in% mediators &
                                   !outcome_children %in% moc &
                                   !outcome_children %in% confounders ]

  # instrumental variables
  instrumental_vars <- extract_instrumental_variables(dag = dag,
                                                      treatments = treatments,
                                                      outcomes = outcomes,
                                                      latent_vars = latent_vars,
                                                      colliders = colliders,
                                                      outcome_parents = outcome_parents,
                                                      nodes_trt_to_y = nodes_trt_to_y)

  # filter latent variables
  latent_vars <- latent_vars[ !latent_vars %in% treatments &
                              !latent_vars %in% outcomes &
                              !latent_vars %in% colliders ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% competing_exposure &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]

  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, "v"] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, "v"] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, "v"] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, "v"] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, "v"] %in% mediators &
                             !edges[, "v"] %in% moc ]  <- "mediator"

  edges$ancestor_iv[edges[, "v"] %in% instrumental_vars] <- "instrumental"

  edges$ancestor_competing_exposure[edges[, "v"] %in% competing_exposure &
                                      !edges[, "v"] %in% proxy_c ] <- "competing_exposure"

  edges$ancestor_collider[edges[, "v"] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, "w"] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, "w"] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, "w"] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, "w"] %in% moc &
                         !edges[, "w"] %in% mediators]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, "w"] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, "w"] %in% instrumental_vars] <- "instrumental"

  edges$descendant_competing_exposure[edges[, "w"] %in% competing_exposure &
                                        !edges[, "w"] %in% proxy_c ] <- "competing_exposure"

  edges$descendant_collider[edges[, "w"] %in% colliders] <- "collider"

  edges$descendant_latent[edges$w %in% latent_vars] <- "latent"

  edges$descendant_observed[edges$w %in% observed] <- "observed"

  return(edges)

}


#' Extract node roles
#'
#' extract_node_roles() is a helper function for get_edges().
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty
#' @param dag A dagitty object.
#' @returns Data table of edges containing roles for each ancestor and descendant nodes.
#' @noRd
extract_node_roles <- function(dag){
  .datatable.aware <- TRUE

  dag <- pdag_to_dag(dag)
  edges <- dagitty::edges(dag)[, c("v", "e", "w")]

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) ]

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  nodes_trt_to_y <- nodes_between_treatment_and_outcome(dag = dag,
                                                            treatments = treatments,
                                                            outcomes = outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% confounders & # remove treatment and confounder nodes
                !moc %in% treatment_parents ] # double check by removing parents of treatment

  # competing exposure
  competing_exposure <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% outcomes ]

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

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children ]

  # instrumental variables
  instrumental_vars <- extract_instrumental_variables(dag = dag,
                                                      treatments = treatments,
                                                      outcomes = outcomes,
                                                      latent_vars = latent_vars,
                                                      colliders = colliders,
                                                      outcome_parents = outcome_parents,
                                                      nodes_trt_to_y = nodes_trt_to_y)

  # filter latent variables
  latent_vars <- latent_vars[  !latent_vars %in% treatments &
                                 !latent_vars %in% outcomes ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% competing_exposure &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]


  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, "v"] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, "v"] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, "v"] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, "v"] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, "v"] %in% mediators ]  <- "mediator"

  edges$ancestor_iv[edges[, "v"] %in% instrumental_vars] <- "instrumental"

  edges$ancestor_competing_exposure[edges[, "v"] %in% competing_exposure &
                                      !edges[, "v"] %in% proxy_c ] <- "competing_exposure"

  edges$ancestor_collider[edges[, "v"] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, "w"] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, "w"] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, "w"] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, "w"] %in% moc ]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, "w"] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, "w"] %in% instrumental_vars] <- "instrumental"

  edges$descendant_competing_exposure[edges[, "w"] %in% competing_exposure &
                                        !edges[, "w"] %in% proxy_c ] <- "competing_exposure"

  edges$descendant_collider[edges[, "w"] %in% colliders] <- "collider"

  edges$descendant_latent[edges$w %in% latent_vars] <- "latent"

  edges$descendant_observed[edges$w %in% observed] <- "observed"

  return(edges)

}


#' Edges to longer format
#'
#' @importFrom data.table as.data.table data.table
#' @param edges Data table of edges, created using extract_node_roles().
#' @returns Data table of edges in longer format.
#' @noRd
edges_longer <- function(edges){
  .datatable.aware <- TRUE

  edges$id <- 1:nrow(edges)

  edges_ancestors <- edges[,1:13]
  edges_descendants <- edges[,c(1:3,14:23)]

  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:13), idvar = "id",
                                      v.names = "ancestor_role", direction = "long")[,c("v", "e", "w", "ancestor_role", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), ][,c("v", "e", "w", "ancestor_role")]

  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:13), idvar = "id",
                                        v.names = "descendant_role", direction = "long")[,c("v", "e", "w", "descendant_role", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), ][,c("v", "e", "w", "descendant_role")]

  if( nrow(edges_ancestors) != nrow(edges_descendants) ){ # finds missing latent edges

    edges_ancestors <- find_missing_edges(edges_ancestors, edges_descendants)

  }

  edges <- do.call(cbind, list( edges_ancestors, data.table::as.data.table(edges_descendants[,4]) ))
  names(edges) <- c("ancestor", "edge", "descendant", "ancestor_role", "descendant_role")
  rownames(edges) <- NULL

  return(edges)

}


#' Long format node roles
#'
#' @importFrom data.table as.data.table is.data.table
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @noRd
roles_longer <- function(dag){
  .datatable.aware <- TRUE

  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  ## ancestor node edges to list ##
  edges_ancestors <- edges_wide[,1:13]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:13), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "w", "role")] )
  names(edges_ancestors)[1:2] <- c("ancestor", "descendant")

  edges_ancestors$node_id <- "ancestor"

  ## descendant node edges to list ##
  edges_descendants <- edges_wide[,c(1:3,14:23)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:13), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "w", "role")] )
  names(edges_descendants)[1:2] <- c("ancestor", "descendant")

  edges_descendants$node_id <- "descendant"

  edges_long <- rbind(edges_ancestors, edges_descendants)

  return(edges_long)

}


#' Find missing latent edges
#'
#' find_missing_edges() is a helper function for get_edges().
#'
#' @importFrom data.table as.data.table data.table setDT setkey
#' @param edges_ancestors Data table of edges containing ancestor node roles.
#' @param edges_descendants Data table of edges containing descendant node roles.
#' @returns A data table containing ancestor edges.
#' @noRd
find_missing_edges <- function(edges_ancestors, edges_descendants){
  .datatable.aware <- TRUE

  # this was written when I had an issue with latent variables being ignored, unsure if it has been fixed.
  # ideally this would be replaced with R core functions, or removed if no longer necessary

  missing_ancestors <- data.table::setDT(edges_descendants)[!edges_ancestors, on = "v"]

  data.table::setkey(edges_descendants,v)
  common_ancestors <- unique(edges_descendants[edges_ancestors,which=TRUE,allow.cartesian=TRUE])
  common_ancestors <- unique(edges_descendants[common_ancestors])

  #combined_ancestors <- rbind(edges_ancestors[,1:3], common_ancestors[,1:3])

  missing_latent_edges <- data.table::setDT(edges_ancestors)[!common_ancestors, on = c("v", "e", "w")]

  if( nrow(missing_latent_edges) > 0 ){

    missing_latent_rows <- lapply( 1:nrow(missing_latent_edges), function(x){
      missing_row <- missing_latent_edges[x, ]
      missing_row[,"role_descendant"][missing_row[, 1] %in% latent_vars] <- "latent"
      missing_row <- data.table::as.data.table(missing_row)

    })

    missing_latent_rows <- do.call(rbind, missing_latent_rows)
    names(missing_latent_rows) <- c("v", "e", "w", "role_ancestor")
    missing_latent_rows <- missing_latent_rows[ missing_latent_rows$role_ancestor == "latent", ]
    edges_ancestors <- rbind(edges_ancestors, missing_latent_rows)
  }

  return(edges_ancestors)

}


#' Connect new nodes
#'
#' add_nodes_helper() is a helper function for add_nodes() that outputs a list containing temporal_reference_node, existing_node_names, and a data frame of new node edges.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty topologicalOrdering
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param nodes A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param node_role Role assigned to new nodes, from any of the following: c("confounder", "treatment", "outcome", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed").
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'first' or 'last', inputted nodes are ordered first or last if a confounder, mediator, or treatment node_role is selected.
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings uses dagitty::topologicalOrdering() and selects the first of the inputted node_role (e.g., first confounder) as the temporal point of reference. If type = 'last', the last node is used.
#' @returns output_list containing temporal_reference_node, existing_node_names, and a data frame of new_edges.
#' @noRd
add_nodes_helper <- function(dag,
                             nodes,
                             node_role,
                             type = "full",
                             temporal_reference_node = NA
                             ){
  .datatable.aware <- TRUE

  e <- "->"

  if(type == "full"){

    e <- "<->"
  }

  ## get initial dag roles
  dag_roles <- get_roles(dag)

  outcomes  <- dag_roles$outcome
  treatments <- dag_roles$treatment
  confounder_vec <- dag_roles$confounder
  m_o_confounder_vec <- dag_roles$mediator_outcome_confounder
  mediator_vec <- dag_roles$mediator
  instrumental_variables <- dag_roles$instrumental
  competing_exposure_vec <- dag_roles$competing_exposure
  latent_vec <- dag_roles$latent
  collider_vec <- dag_roles$collider
  observed <- dag_roles$observed

  ##variable names
  observed_node_names <- unique( c(confounder_vec, m_o_confounder_vec, mediator_vec, competing_exposure_vec, collider_vec, instrumental_variables, observed) )
  observed_node_names <- Filter(Negate(anyNA), observed_node_names)
  observed_node_names <- observed_node_names[ !observed_node_names %in% latent_vec ]

  nodes_ordered <- sort( unlist( dagitty::topologicalOrdering(dag) ) ) # ggdag estimated temporal order of new nodes

  if( length( node_role) > 1 | length( node_role) == 0 ){

    stop("add_nodes() currently only supports single node_role character inputs.")

  }else if( node_role %in% "confounder" ){
    confounder_occurrance <- as.numeric(order(match(nodes, nodes)))

    ## confounder edges ##
    new_edges <- draw_confounder_edges(type = type,
                                           confounders = confounders,
                                           confounder_vec = nodes, # new nodes as confounder
                                           treatments = treatments,
                                           outcomes = outcomes,
                                           m_o_confounder_vec = m_o_confounder_vec,
                                           mediator_vec = mediator_vec,
                                           latent_vec = latent_vec,
                                           e = e)


    # connect all confounders (fully connected or saturated graph type)
    if( type == "full" | type == "saturated" ){

      confounder_list <- list()

      confounder_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

          list( c( ancestor = nodes[x], edge = e, descendant = confounder_vec[y]) )

        })

      }) )

      confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
      confounder_unlist <- data.table::as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, confounder_unlist)


      confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

        confounder_list[x] <- lapply(1:length(nodes), function(y){

          list( c( ancestor = confounder_vec[x], edge = e, descendant = nodes[y]) )

        })

      }) )

      confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
      confounder_unlist <- data.table::as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, confounder_unlist)



    }else if( type == "first"){

      first_var <- names(nodes_ordered)[ names(nodes_ordered) %in% confounder_vec ][1]

      temporal_reference_node <- first_var

      first_list <- list()

      first_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        first_list[x] <- list( c( ancestor = nodes[x], edge = e, descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), first_list)
      first_unlist <- data.table::as.data.table( do.call( rbind, unlist(first_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, first_unlist)


    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% confounder_vec ][length(confounder_vec) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = e, descendant = nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), last_list)
      last_unlist <- data.table::as.data.table( do.call( rbind, unlist(last_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, last_unlist)


    }

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

    # connect all treatments (fully connected or saturated graph type)
    if( type == "full" ){

      treatment_list <- list()

      treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

        treatment_list[x] <- lapply(1:length(nodes), function(y){

          list( c( ancestor = treatments[x], edge = e, descendant = nodes[y]) )

        })

      }) )

      treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
      treatment_unlist <- data.table::as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, treatment_unlist)



      treatment_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        treatment_list[x] <- lapply(1:length(treatments), function(y){

          list( c( ancestor = nodes[x], edge = e, descendant = treatments[y]) )

        })

      }) )

      treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
      treatment_unlist <- data.table::as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, treatment_unlist)


    }else if( type == "first"){

      first_var <- names(nodes_ordered)[ names(nodes_ordered) %in% treatments ][1]

      temporal_reference_node <- first_var

      first_list <- list()

      first_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        first_list[x] <- list( c( ancestor = nodes[x], edge = e, descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), first_list)
      first_unlist <- data.table::as.data.table( do.call( rbind, unlist(first_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, first_unlist)


    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% treatments ][length(treatments) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = e, descendant = nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), last_list)
      last_unlist <- data.table::as.data.table( do.call( rbind, unlist(last_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, last_unlist)


    }

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

    if(type == "full"){

      mediator_list <- list()

      mediator_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        mediator_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = nodes[x], edge = e, descendant = mediator_vec[y]) )

        })

      }) )

      mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
      mediator_unlist <- data.table::as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, mediator_unlist)


      mediator_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        mediator_list[x] <- lapply(1:length(nodes), function(y){

          list( c( ancestor = mediator_vec[x], edge = e, descendant = nodes[y]) )

        })

      }) )

      mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
      mediator_unlist <- data.table::as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, mediator_unlist)


    }else if( type == "first"){

      first_var <- names(nodes_ordered)[ names(nodes_ordered) %in% mediator_vec ][1]

      temporal_reference_node <- first_var

      first_list <- list()

      first_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        first_list[x] <- list( c( ancestor = nodes[x], edge = e, descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), first_list)
      first_unlist <- data.table::as.data.table( do.call( rbind, unlist(first_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, first_unlist)


    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% mediator_vec ][length(mediator_vec) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = e, descendant = nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), last_list)
      last_unlist <- data.table::as.data.table( do.call( rbind, unlist(last_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, last_unlist)


    }

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

  }else if( node_role %in% "instrumental" ){
    ## instrumental_variables edges ##
    new_edges <- draw_iv_edges(type = type,
                               instrumental_variables = nodes, # new nodes as instrumental
                               treatments = treatments,
                               e = e)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% instrumental_variables ] ) # existing node names in temporal order

  }else if( node_role %in% "competing_exposure" ){
    ## competing_exposure edges ##
    new_edges <- draw_competing_exposure_edges(outcomes = outcomes,
                                               competing_exposure_vec = nodes,
                                               e = e) # new nodes as competing_exposure

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% competing_exposure_vec ] ) # existing node names in temporal order

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
    ## connect observed to ancestors and descendants ##
    new_edges <- draw_observed_edges(observed = nodes, # new nodes as observed
                                     existing_dag = dag,
                                     e = e)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% observed ] ) # existing node names in temporal order

  }else{

    stop("Invalid node_role input. Check documentation for valid currently support inputs and try again.")

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
#' saturate_nodes_helper() is a helper function for saturate_nodes() that connects new and existing nodes, drawing edges in both directions.
#'
#' @importFrom data.table as.data.table
#' @param dag An existing dagitty object.
#' @param nodes A vector of new nodes.
#' @noRd
saturate_nodes_helper <- function(dag, nodes, dag_node_names, type){
  .datatable.aware <- TRUE

  ## get node names
  node_names <- dag_node_names

  ## get initial dag roles
  dag_roles <- get_roles(dag)

  outcomes  <- dag_roles$outcome
  treatments <- dag_roles$treatment
  confounder_vec <- dag_roles$confounder
  m_o_confounder_vec <- dag_roles$mediator_outcome_confounder
  mediator_vec <- dag_roles$mediator
  instrumental_variables <- dag_roles$instrumental
  competing_exposure_vec <- dag_roles$competing_exposure
  latent_vec <- dag_roles$latent
  latent_variables <- latent_vec
  collider_vec <- dag_roles$collider
  observed <- dag_roles$observed

  confounder_occurrance <- as.numeric(order(match(nodes, nodes)))

  ## get variable names ##
  observed_node_names <- unique( as.vector( c(confounder_vec, m_o_confounder_vec, mediator_vec, competing_exposure_vec, collider_vec, instrumental_variables) ) )
  observed_node_names <- Filter(Negate(anyNA), observed_node_names)

  edges <- draw_edges(observed_node_names = observed_node_names,
                      type = type,
                      outcomes = outcomes,
                      treatments = treatments,
                      confounder_vec = confounder_vec,
                      m_o_confounder_vec = m_o_confounder_vec,
                      mediator_vec = mediator_vec,
                      instrumental_variables = instrumental_variables,
                      competing_exposure_vec = competing_exposure_vec,
                      latent_vec = latent_vec,
                      latent_variables = latent_variables,
                      collider_vec = collider_vec,
                      observed = observed,
                      confounder_occurrance = confounder_occurrance,
                      existing_dag = dag)

  colnames(edges) <- c("v", "e", "w")

  return( edges )

}


#' @importFrom data.table as.data.table
#' @param edges Data frame / data table of edges.
#' @returns List of edges, grouped by each unique node.
#' @noRd
print_edges_helper <- function(new_edges){
  .datatable.aware <- TRUE

  ## collapse new_edges to a vector and output grouped by nodes
  unique_ancestors <- unique( new_edges[,1] )
  num_unique_ancestors <- nrow(unique_ancestors)


  new_edges_list <- list()

  # new_edges grouped by each unique node in a list
  new_edges_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    new_edges[ unlist(new_edges[,1]) %in% unlist(unique_ancestors)[x], ]


  }) )

  # nodes containing edges are collapsed, outputted in the console to allow easy assessing by copy and paste into .r file
  new_edges_list <- lapply(1:num_unique_ancestors, function(x){

    new_edges_list <- noquote(
      paste0( paste0( "'", sapply(1:nrow(new_edges_list[[x]]), function(y){

        new_edges_list[x] <- noquote( paste( new_edges_list[[x]][y,], collapse=" "  ) )

      }), "'", collapse=", " ), sep = "") )

  })


  return(new_edges_list)
}


#' Adds nodes to dagitty objects
#'
#' copy_nodes_helper() is a helper function for copy_nodes(). It adds new nodes to a dagitty object by referencing existing nodes.
#'
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param existing_nodes Vector of existing node names, used as a reference for the new graph nodes, e.g., c("Z1", "Z2", "Z3").
#' @param existing_node_type A suffix added to each of the reference node names, e.g. "pre_treatment", or "t0".
#' @param new_node_type A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings use treatment as the temporal point of reference for adding new (post-treatment) nodes, but if other node types are specified it can be useful to specify the existing node name.
#' @param num_repeats Number of additional copies of nodes, such as time points. Each repeat number is included at the end of new node names (new_new_t1, new_node_t2, etc.).
#' @param coords_spec Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures.
#' @returns A dagitty object.
#' @noRd
copy_nodes_helper <- function(dag,
                              edges,
                              node_names,
                              treatments,
                              outcomes,
                              latent_vec,
                              coordinates,
                              existing_nodes,
                              existing_node_type,
                              new_node_type,
                              temporal_reference_node,
                              num_repeats,
                              coords_spec
){

  # create new node names
  new_node_names <- create_names_copy_nodes_helper(existing_nodes, new_node_type, num_repeats)
  new_node_names <- as.data.frame(new_node_names) # convert to data frame (called here instead of in the function to keep a consistent output, will possibly move it in later)

  # update reference node names based on inputted existing_node_type
  existing_node_names <- paste0(existing_nodes, "_", existing_node_type)

  # replace reference nodes in the dag with their new names (coordinates too)
  num_ref_nodes <- length(existing_nodes)

  exclude_names <- c(treatments, outcomes, latent_vec)

  node_names <- node_names[!node_names %in% exclude_names]

  for(i in 1:num_ref_nodes){ # this will be rewritten without a for-loop

    latent_vec <- sapply(   latent_vec, function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    treatments <- sapply(   treatments, function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    outcomes <- sapply(   outcomes, function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    node_names <- sapply(   node_names, function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    edges[,"v"] <- sapply(   as.vector(edges[,"v"]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    edges[,"w"] <- sapply( as.vector(edges[,"w"]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    names(coordinates[["x"]]) <- sapply(   names(coordinates[["x"]]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    names(coordinates[["y"]]) <- sapply(   names(coordinates[["y"]]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

  }

  dag <- construct_graph(edges, node_names, treatments, outcomes, latent_vec)

  # draw edges for all new nodes (includes connecting parent and children nodes, reference nodes to new nodes, and consecutive new nodes)
  new_node_df <- draw_edges_for_copy_nodes_helper(dag, new_node_names, existing_node_names, existing_node_type, new_node_type, temporal_reference_node, num_repeats)


  edges <- rbind(edges, new_node_df)
  edges[] <- lapply(edges, as.character)

  dag <- construct_graph(edges, node_names, treatments, outcomes, latent_vec)


  if(all(!is.na(unlist(coordinates)))){

    dagitty::coordinates(dag) <- coordinates

  }


  new_node_names <- as.vector(unlist(new_node_names))

  num_nodes <- length(length(new_node_names))
  time_limit <- num_nodes + num_nodes*coords_spec[1]

  setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)

  on.exit( {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  } )


  tryCatch({
    new_coordinates <- renew_coords(dag = dag,
                                    new_node_names = new_node_names,
                                    coordinates = coordinates,
                                    coords_spec = coords_spec[ complete.cases(coords_spec) ] )

    dagitty::coordinates(dag) <- new_coordinates

  }, warning = function(w){
    message(paste("Warning:", w, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords_helper(dag)
    return(dag)

  }, error = function(e){
    message(paste("Error:", e, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords_helper(dag)
    return(dag)

  }, finally = return(dag)
  )

}

#' merge dags
#'
#' copy_nodes_helper2() is a helper function for copy_nodes().
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty
#' @param edges A data frame of edges.
#' @returns A dagitty object
#' @noRd
copy_nodes_helper2 <- function(dag,
                              edges,
                              node_names,
                              new_dag,
                              treatments,
                              outcomes,
                              latent_vec,
                              coordinates,
                              coords_spec
                              ){
  .datatable.aware <- TRUE
  # get edges, node names from new daggity object (inputted as existing_nodes)
  new_edges <- pdag_to_dag_edges(new_dag)

  new_node_names <- names(new_dag) # extract new dag node names
  existing_node_names <-  new_node_names[ new_node_names %in% node_names ] # saves duplicate node names
  new_node_names <- new_node_names[ !new_node_names %in% node_names ] # remove duplicate node names

  node_names <- c(node_names, new_node_names) # combine both dag node names

  edges <- rbind(edges, new_edges) # combine both dag edges

  dag <- rebuild_dag(dag, edges)


  num_nodes <- length(length(new_node_names))
  time_limit <- num_nodes + num_nodes*coords_spec[1]

  setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)

  on.exit( {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  } )

  tryCatch({
    new_coordinates <- renew_coords(dag = dag,
                                new_node_names = new_node_names,
                                coordinates = coordinates,
                                coords_spec = coords_spec[ complete.cases(coords_spec) ] )
    dagitty::coordinates(dag) <- new_coordinates

  }, warning = function(w){
    message(paste("Warning:", w, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords_helper(dag)
    return(dag)

  }, error = function(e){
    message(paste("Error:", e, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords_helper(dag)
    return(dag)

  }, finally = return(dag)
  )

}

#' Rebuild dag
#'
#' rebuild_dag() rebuilds a dag using a dagitty object and data frame of edges input.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty isAcyclic
#' @param dag A dagitty object.
#' @param edges A vector of edges.
#' @returns A dagitty object
#' @noRd
rebuild_dag <- function(dag,
                        edges
){
  .datatable.aware <- TRUE

  latent_vec <- dagitty::latents(dag)
  treatments <- dagitty::exposures(dag)
  outcomes <- dagitty::outcomes(dag)

  coordinates <- dagitty::coordinates(dag)

  exclude_names <- c(treatments, outcomes, latent_vec)

  node_names <- unique(edges[,"v"])
  node_names <- node_names[,1][! unlist( node_names[,1] ) %in% exclude_names]

  # construct dag
  dag <- construct_graph(edges, unlist(node_names), treatments, outcomes, latent_vec)

  if(all(!is.na(unlist(coordinates)))){
    dagitty::coordinates(dag) <- coordinates
  }

  return(dag)

}


#' construct dag
#'
#' construct_graph() constructs a daggity object using edges input.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty
#' @param edges A data frame of edges.
#' @returns A dagitty object
#' @noRd
construct_graph <- function(edges,
                            node_names,
                            treatments,
                            outcomes,
                            latent_vec
                            ){
  .datatable.aware <- TRUE

  # handle simultaneous edges
  names(edges) <- c("v", "e", "w")
  edges <- handle_simultaneous_edges(edges)

  outcome_vec <- unlist(outcomes)
  node_name_and_coords_vec <- c()

  if( length(latent_vec) > 0 ){ ########### simplify this section ##############

    if( all(complete.cases(latent_vec)) ) {
      node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                    paste(outcome_vec, " [outcome] ", sep=""),
                                    paste(latent_vec, " [latent] ", sep=""),
                                    paste(node_names, collapse=" "))
    }else{
      node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                    paste(outcome_vec, " [outcome] ", sep=""),
                                    paste(node_names, collapse=" "))
    }

  }else if( length(outcome_vec) + length(treatments) == 0 ){

    node_name_and_coords_vec <-  paste(node_names, collapse=" ")


  }else if( length(treatments) == 0 ){

    node_name_and_coords_vec <- c(paste(outcome_vec, " [outcome] ", sep=""),
                                  paste(node_names, collapse=" "))

  }else if( length(outcome_vec) == 0 ){

    node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                  paste(node_names, collapse=" "))

  }else{
    node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                  paste(outcome_vec, " [outcome] ", sep=""),
                                  paste(node_names, collapse=" "))
  }

  num_edges <- nrow(edges)

  if( num_edges != 0 ){

    edges <- suppressWarnings( sapply(1:num_edges, function(x){

      edges <- paste( edges[x,], collapse=" ")

    }) )

  }

  dag <- paste("dag {", paste(node_name_and_coords_vec, collapse=""), paste(edges, collapse=" "), "}", sep = " ")

  dag <- dagitty::dagitty(dag)

  return(dag)

}


#' Adds coordinates to nodes in a merged dagitty object
#'
#' merged_node_coords_helper() is a helper function for renew_coords().
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_node_name_vec Inputted vector of new node names added to the graph.
#' @param num_nodes Number of new nodes.
#' @param coordinates_x Existing x-coordinates.
#' @param coordinates_y Existing y-coordinates.
#' @param treatments Treatment/exposures in the supplied dag.
#' @param outcomes Outcomes in the supplied dag.
#' @param post_treatment Indicates whether the coordinates of new nodes should be informed by treatment coordinates.
#' @param post_outcome Similarly, outcome coordinates can be used to inform new node coordinates.
#' @param num_vars Total number of variables in the dag.
#' @param lambda Parameter used to generate coordinates. Adjust node placement with lambda; a higher value increases volatility and results in more extreme DAG structures.
#' @param threshold Parameter used to generate coordinates. Threshold controls the closeness of nodes.
#' @return Named list of new coordinates.
#' @noRd
merged_node_coords_helper <- function(dag,
                                      new_node_name_vec,
                                      num_nodes,
                                      coordinates_x,
                                      coordinates_y,
                                      treatments = NA,
                                      outcomes,
                                      post_treatment = FALSE,
                                      post_outcome = FALSE,
                                      num_vars,
                                      lambda,
                                      threshold
){
  .datatable.aware <- TRUE

  existing_node_names <- names(coordinates_x)
  existing_node_names_not_outcome <- existing_node_names[ !existing_node_names %in% outcomes ]

  nodes_parents <- lapply(1:num_nodes, function(x){   # get new node name parents
    nodes_parents <- dagitty::parents(dag, new_node_name_vec[x])
  })

  nodes_children <- lapply( 1:num_nodes, function(x){  # get new node name children
    nodes_children <- dagitty::children(dag, new_node_name_vec[x])
  })

  if( post_treatment == FALSE & all( complete.cases(treatments) ) ){

    existing_node_names_not_outcome <- existing_node_names_not_outcome[ !existing_node_names_not_outcome %in% treatments ]

  }

  if( post_outcome == FALSE ){

    nodes_parents <- lapply(1:num_nodes, function(x){ # parents of new nodes found in both dags prior to merge
      nodes_parents <- existing_node_names_not_outcome[ existing_node_names_not_outcome %in% nodes_parents[[x]] ]
    })

    nodes_children <- lapply(1:num_nodes, function(x){ # children of new nodes found in both dags prior to merge
      nodes_children <- existing_node_names_not_outcome[ existing_node_names_not_outcome %in% nodes_children[[x]] ]
    })

  }

  ## separate x and y coordinates, node names, num nodes
  new_node_name_vec_y <- new_node_name_vec
  num_nodes_y <- num_nodes
  nodes_parents_y <- nodes_parents
  nodes_children_y <- nodes_children

  new_node_name_vec_x <- new_node_name_vec
  num_nodes_x <- num_nodes
  nodes_parents_x <- nodes_parents
  nodes_children_x <- nodes_children

  existing_coordinates_x <- coordinates_x[ !names(coordinates_x) %in% new_node_name_vec ]
  existing_coordinates_y <- coordinates_y[ !names(coordinates_y) %in% new_node_name_vec ]

  if( length(unname(existing_coordinates_y)) == 0 | length(unname(existing_coordinates_x)) == 0 ){
    # if after removing new nodes (as per above) no existing coordinates remain to use as a reference:
    existing_coordinates_y <- 1
    existing_coordinates_x <- 1
    # find root node to assign first coordinates [0,0]
    root_node <- find_root_node(new_node_name_vec,
                                nodes_parents,
                                nodes_children)
    if( length(root_node) !=0 ){

      names(existing_coordinates_y) <- root_node[1]
      names(existing_coordinates_x) <- root_node[1]

    }else{

      names(existing_coordinates_y) <- new_node_name_vec_x[1]
      names(existing_coordinates_x) <- new_node_name_vec_x[1]

    }

    # y-coords
    nodes_children_y <- nodes_children_y[ !new_node_name_vec_y %in% names(existing_coordinates_y) ]
    nodes_parents_y <- nodes_parents_y[ !new_node_name_vec_y %in% names(existing_coordinates_y) ]

    new_node_name_vec_y <- new_node_name_vec_y[ !new_node_name_vec_y %in% names(existing_coordinates_y)]
    num_nodes_y <- length(new_node_name_vec_y)

    # x-coords
    nodes_children_x <- nodes_children_x[ !new_node_name_vec_x %in% names(existing_coordinates_x) ]
    nodes_parents_x <- nodes_parents_x[ !new_node_name_vec_x %in% names(existing_coordinates_x) ]

    new_node_name_vec_x <- new_node_name_vec_x[ !new_node_name_vec_x %in% names(existing_coordinates_x) ]
    num_nodes_x <- length(new_node_name_vec_x)

  }

  quality_check <- FALSE
  iteration <- 1

  while(quality_check == FALSE){

    iteration <- iteration + ( num_nodes*lambda )

    ## get difference between min-max coords
    diff_y_coords <-  abs( min( existing_coordinates_y ) - max( existing_coordinates_y ) )
    diff_x_coords <-  abs( min( existing_coordinates_x ) - max( existing_coordinates_x ) )


    if( num_nodes_y > 0 ){
      ## y coordinates ##

      new_y_coords <- sapply(1:num_nodes_y, function(x){

        if( length( existing_coordinates_y[ names(existing_coordinates_y) %in% nodes_parents_y[[x]] ] ) > 0 ){


          new_y_coords <- max( existing_coordinates_y[ names(existing_coordinates_y) %in% nodes_parents_y[[x]] ]
          ) + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda)/2 )

        }else if( length( existing_coordinates_y[ names(existing_coordinates_y) %in% nodes_children_y[[x]] ] ) > 0 ) {

          new_y_coords <- min( existing_coordinates_y[ names(existing_coordinates_y) %in% nodes_children_y[[x]] ]
          ) - x*iteration - ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) - runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda )/2 )

        }else if( length( existing_coordinates_y[ names(existing_coordinates_y) %in% unlist(nodes_parents_y) ] ) > 0 ){

          new_y_coords <- min( existing_coordinates_y[ names(existing_coordinates_y) %in% unlist(nodes_parents_y) ]
          ) + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda )/2 )

        }else if( length( existing_coordinates_y[ names(existing_coordinates_y) %in% unlist(nodes_children_y) ] ) > 0 ){

          new_y_coords <- min( existing_coordinates_y[ names(existing_coordinates_y) %in% unlist(nodes_children_y) ]
          ) - x*iteration - ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) - runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda )/2 )

        }else{

          new_y_coords <- x + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                                       min = iteration*lambda,
                                                                                                       max = iteration + ( num_nodes*lambda )/2 )
        }

      })

      names(new_y_coords) <- new_node_name_vec_y

    }

    if( num_nodes_x > 0){
      ## x coordinates ##
      new_x_coords <- sapply(1:num_nodes_x, function(x){

        if( length( existing_coordinates_x[ names(existing_coordinates_x) %in% nodes_children_x[[x]] ] ) > 0 ) {

          new_x_coords <- min( existing_coordinates_x[ names(existing_coordinates_x) %in% nodes_children_x[[x]] ]
          ) - x*iteration - (diff_x_coords/num_vars)*num_nodes + diff_x_coords*x - runif(n = 1,
                                                                                                           min = iteration*lambda,
                                                                                                           max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( existing_coordinates_x[ names(existing_coordinates_x) %in% nodes_parents_x[[x]] ] ) > 0 ){

          new_x_coords <- min( existing_coordinates_x[ names(existing_coordinates_x) %in% nodes_parents_x[[x]] ]
          ) + x*iteration + (diff_x_coords/num_vars)*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( existing_coordinates_x[ names(existing_coordinates_x) %in% unlist(nodes_children_x) ] ) > 0 ){

          new_x_coords <- min( existing_coordinates_x[ names(existing_coordinates_x) %in% unlist(nodes_children_x) ]
          ) - x*iteration - (diff_x_coords/num_vars)*num_nodes + diff_x_coords*x - runif(n = 1,
                                                                                                           min = iteration*lambda,
                                                                                                           max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( existing_coordinates_x[ names(existing_coordinates_x) %in% unlist(nodes_parents_x) ] ) > 0 ){

          new_x_coords <- max( existing_coordinates_x[ names(existing_coordinates_x) %in% unlist(nodes_parents_x) ]
          ) + x*iteration + diff_x_coords*x + runif(n = 1,
                                                                    min = iteration*lambda,
                                                                    max = iteration + ( num_nodes*lambda*10 )/2 )
        }else{

          new_x_coords <- x + x*iteration + diff_x_coords*x - runif(n = 1,
                                                                                    min = iteration*lambda,
                                                                                    max = iteration + ( num_nodes*lambda*10 )/2 )
        }

      })

      names(new_x_coords) <- new_node_name_vec_x
    }

    new_coordinates <- quality_check_merged_node_coords_helper(new_node_name_vec_x = new_node_name_vec_x,
                                                               new_node_name_vec_y = new_node_name_vec_y,
                                                               num_nodes_x = num_nodes_x,
                                                               num_nodes_y = num_nodes_y,
                                                               new_x_coords = new_x_coords,
                                                               new_y_coords = new_y_coords,
                                                               existing_coordinates_x = existing_coordinates_x,
                                                               existing_coordinates_y = existing_coordinates_y,
                                                               threshold = threshold)

    ## initialise new variables ##
    # y-coords
    existing_coordinates_y <-  new_coordinates$y

    nodes_children_y <- nodes_children_y[ !new_node_name_vec_y %in% names(existing_coordinates_y) ]
    nodes_parents_y <- nodes_parents_y[ !new_node_name_vec_y %in% names(existing_coordinates_y) ]

    new_node_name_vec_y <- new_node_name_vec_y[ !new_node_name_vec_y %in% names(existing_coordinates_y)]
    num_nodes_y <- length(new_node_name_vec_y)

    # x-coords
    existing_coordinates_x <- new_coordinates$x

    nodes_children_x <- nodes_children_x[ !new_node_name_vec_x %in% names(existing_coordinates_x) ]
    nodes_parents_x <- nodes_parents_x[ !new_node_name_vec_x %in% names(existing_coordinates_x) ]

    new_node_name_vec_x <- new_node_name_vec_x[ !new_node_name_vec_x %in% names(existing_coordinates_x) ]
    num_nodes_x <- length(new_node_name_vec_x)

    if( num_nodes_y + num_nodes_x == 0 ){

      quality_check <- TRUE

    }

  }

  coordinates <- list(x = existing_coordinates_x[!duplicated(names(existing_coordinates_x))], y = existing_coordinates_y[!duplicated(existing_coordinates_y)] )

  return(coordinates)

}


#' Creates latent node coordinates for add_latent_coordinates()
#'
#' latent_new_coordinates_helper() is a helper function for renew_coords().
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @param dag dagitty object
#' @param latent_variables vector of latent variable nodes.
#' @param existing_coordinates list of coordinates from a dagitty object.
#' @param lambda coordinates tuning parameter.
#' @return dagitty objecty with coordinates.
#' @noRd
latent_new_coordinates_helper <- function(dag, latent_variables, coordinates, lambda=2){
  #latent_variables <- get_latent_vec(latent_variables) # used for debugging
  # coordinates <- existing_coords# used for debugging
  if( all( complete.cases( unlist(latent_variables) ) ) ){

    num_latents <- length(latent_variables)

    num_vars <- length(names(dag))

    latent_descendants <- lapply(1:num_latents, function(x){

      latent_descendants <- dagitty::children(dag, latent_variables[x])

    })

    latent_parents <- lapply(1:num_latents, function(x){

      latent_parents <- dagitty::parents(dag, latent_variables[x])

    })

    latent_descendants <- latent_descendants[[1]][ ! latent_descendants[[1]] %in% latent_parents[[1]]]


  }else{

    if( all( complete.cases(latent_variables) ) ){

      warning("NA's present in latent variables (ignore if latent variables were not supplied)")

    }

    return(coordinates)

  }

  quality_check <- FALSE

  iteration <- 0

  lambda <- lambda/num_vars

  while(quality_check == FALSE){

    iteration <- iteration + 1*lambda

    ## y coordinates ##
    new_y_coords <- sapply(1:num_latents, function(x){

      if( length( latent_descendants[[x]] ) > 0 ) {

        new_y_coords <- min( coordinates$y[ names(coordinates$y) %in% latent_descendants[[x]] ]
        ) - runif(n = 1,
                  min = ( num_latents),
                  max = iteration + (num_latents + x) )*lambda - 1

      }else if( length( unlist(latent_descendants) > 0) ){

        new_y_coords <- min( coordinates$y[ names(coordinates$y) %in% unlist(latent_descendants) ]
        ) - runif(n = 1,
                  min = ( num_latents ),
                  max = iteration + (num_latents + x) )*lambda - 1

      }else if( length( unlist(latent_parents) ) > 0 ){

        new_y_coords <- min( coordinates$y[ names(coordinates$y) %in% unlist(latent_parents) ]
        ) + runif(n = 1,
                  min = ( num_latents ),
                  max = iteration + (num_latents + x) )*lambda + 1

      }else{

        new_y_coords <- x - runif(n = 1,
                                  min = ( num_latents/num_vars ),
                                  max = iteration + (num_latents + x) )*lambda

      }

    })

    ## x coordinates ##
    new_x_coords <- sapply(1:num_latents, function(x){

      if( length( latent_descendants[[x]] ) > 0 ) {

        new_node_x_coords <- min( coordinates$x[ names(coordinates$x) %in% latent_descendants[[x]] ]
        ) - runif(n = 1,
                  min = ( num_latents/num_vars ),
                  max = iteration + (num_latents + x)/2 ) - 1

      }else if( length( unlist(latent_descendants) > 0) ){

        new_x_coords <- min( coordinates$x[ names(coordinates$x) %in% unlist(latent_descendants) ]
        ) + runif(n = 1,
                  min = ( num_latents/num_vars ),
                  max = iteration + (num_latents + x) ) + 1

      }else{

        new_x_coords <- x - runif(n = 1,
                                  min = ( num_latents/num_vars ),
                                  max = iteration + (num_latents + x) ) - 1

      }

    })

    new_coordinates <- list(x = new_x_coords, y = new_y_coords)

    quality_check <- quality_check_new_coordinates_helper(dag,
                                                          grouped_nodes = latent_variables,
                                                          num_nodes = num_latents,
                                                          new_coordinates = new_coordinates,
                                                          existing_coords = coordinates)

  }

  names(new_coordinates$x) <- latent_variables
  names(new_coordinates$y) <- latent_variables

  coordinates <- list(x = c(coordinates$x, new_coordinates$x), y = c(coordinates$y, new_coordinates$y))

  return(coordinates)

}


#' Creates coordinates for obserbed nodes.
#'
#' observed_new_coordinates_helper()
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @param dag dagitty object
#' @param latent_variables vector of latent variable nodes.
#' @param existing_coordinates list of coordinates from a dagitty object.
#' @return dagitty objecty with coordinates.
#' @noRd
observed_new_coordinates_helper <- function(dag, observed, existing_coords){

  if( !all(is.na(observed)) ){
    num_nodes <- length(observed)

    num_vars <- length(names(dag))

    nodes_descendants <- lapply(1:num_nodes, function(x){

      nodes_descendants <- dagitty::children(dag, observed[x])

    })

    if( length(unlist(nodes_descendants)) < (num_nodes*2) ){

      nodes_parents <- lapply(1:num_nodes, function(x){

        nodes_parents <- dagitty::parents(dag, observed[x])

      })

      nodes_descendants <- lapply(1:num_nodes, function(x){

        nodes_descendants <- c(unlist(nodes_parents[x]), unlist(nodes_descendants[x]))

      })

    }


    quality_check <- FALSE

    iteration <- 0

    while(quality_check == FALSE){

      iteration <- iteration + 1

      new_group_coords <- lapply(1:num_nodes, function(x){

        nodes_descendants_coords_x <- ( mean(
          replace( existing_coords$x[ names(existing_coords$x) %in% nodes_descendants[[x]] ],
                   is.na( existing_coords$x[ names(existing_coords$x) %in% nodes_descendants[[x]] ] ),
                   0) )

          - runif(n = 1, min = ( 0.1*iteration )*( num_vars/num_nodes ),
                  max = ( ( 0.2*iteration )*( num_vars/num_nodes ) ) ) ** (2*(0.1*iteration))

          - (2*( 1 / x )) - ( (0.2*iteration)*( num_vars ) )
        )

        nodes_descendants_coords_y <- ( mean(
          replace( existing_coords$y[ names(existing_coords$y) %in% nodes_descendants[[x]] ],
                   is.na(existing_coords$y[ names(existing_coords$y) %in% nodes_descendants[[x]] ] ),
                   0) )

          - runif(n = 1, min = ( 0.01*iteration )*( num_vars/num_nodes ),
                  max = ( ( 0.1*iteration )*( num_vars/num_nodes ) ) ) ** (2*(0.1*iteration))

          - (2*( 1 / x )) - ( (0.2*iteration)*( num_vars ) )
        )

        new_group_coords <- list(x = unlist(nodes_descendants_coords_x), y = unlist(nodes_descendants_coords_y))

      })

      # bind rows into a data frame (columns are still classed as lists)
      new_coordinates <- as.data.frame( do.call( rbind, new_group_coords ) )


      ## unlist and replace NaN ##

      # x coordinates
      unlisted_coords_x <- unlist( new_coordinates$x )

      unlisted_coords_x <- sapply(  unlisted_coords_x, function(x){
        replace( x, x == "NaN" , num_vars + runif(1, min=1, max = (num_vars/2) ) )
      } )

      # y coordinates
      unlisted_coords_y <- unlist( new_coordinates$y )

      unlisted_coords_y <- sapply(  unlisted_coords_y, function(x){
        replace( x, x == "NaN" , num_vars + runif(1, min=1, max = (num_vars/2) ) )
      } )



      ## internal check against other generated latent coords ##

      if( !all( unlisted_coords_x >= 0 ) ){

        unlisted_coords_x <- lapply(1:length(unlisted_coords_x), function(x){

          unlisted_coords_x[x] <- sqrt( unlist( unlisted_coords_x[x] )**2 )/2

        })

        unlisted_coords_x <- unlist(unlisted_coords_x)
      }

      if( !all( unlisted_coords_y >= 0 ) ){

        unlisted_coords_y <- lapply(1:length(unlisted_coords_y), function(x){

          unlisted_coords_y[x] <- sqrt( unlist( unlisted_coords_y[x] )**2 )/2

        })

        unlisted_coords_y <- unlist(unlisted_coords_y)

      }

      new_coordinates <- list(x = unlisted_coords_x, y = unlisted_coords_y)

      quality_check <- quality_check_new_coordinates_helper(dag, observed, num_nodes, new_coordinates, existing_coords)

    }
    names(new_coordinates$x) <- observed
    names(new_coordinates$y) <- observed

    coords_list <- list(x = c( na.omit(existing_coords$x), new_coordinates$x), y = c( na.omit(existing_coords$y), new_coordinates$y))

  }else{

    coords_list <- existing_coords
  }


  return(coords_list)

}




#' Returns highest x or y coordinates for a group of nodes
#'
#' @importFrom dagitty parents
#' @param dag A dagittyy object.
#' @param group_node_names Node group of interest.
#' @param num_group_nodes Number of nodes in group of interest.
#' @param coordinates_x_or_y Named vector of either x or y coordinates.
#' @param coord_names Names corrresponding to coordinates_x_or_y input.
#' @return Maximum parent coordinates for a group of nodes.
#' @noRd
parent_node_max_coords <- function(dag, group_node_names, num_group_nodes, coordinates_x_or_y, coord_names){

  group_parents <- lapply(1:num_group_nodes, function(x){   # get new node name parents
    group_parents <- dagitty::parents(dag, group_node_names[x])
  })

  if( length( unlist(group_parents) ) > 0 ){

    group_parent_max_coord <- max( coordinates_x_or_y[ coord_names # highest x or y coordinates for a group of nodes
                                                       %in% unlist(group_parents) ], na.rm = TRUE )

  }else{

    group_parent_max_coord <- 0

  }

  return(group_parent_max_coord)

}



#' dataframe output from a dagitty object
#'
#'Generates a table similar to calling ggdag::tidy_dagitty on a ggdag::dagify object
#'The benefit of this function is that it automatically identifies exposure, outcome, confounder, observed and latent variables inputted from dagitty.net, whereas ggdag::tidy_dagitty only does this for ggdag::dagify objects.
#'Output can be used with ggdag to create better looking DAGs from dagitty.net code.
#'
#' @importFrom ggdag tidy_dagitty
#' @param dag dagitty object
#' @param labels vector of labels for nodes in dagitty object
#' @return dag_df DAG as a dataframe for use with ggdag to create better looking DAGs
#' @noRd
tidy_ggdagitty <- function(dag, labels = NULL){
  # Cleaning the dags and turning it into a data frame.
  dag_df <- data.frame(ggdag::tidy_dagitty(dag))

  # flip y axis for ggplot
  dag_df$y <- dag_df$y*-1
  dag_df$yend <- dag_df$yend*-1

  if( is.null(labels) ){

    return(dag_df)

  }

  dag_df <- add_labels(dag, dag_df, labels)

  return(dag_df)
}


#' add labels to a dag dataframe
#'
#'Generates a table similar to calling ggdag::tidy_dagitty on a ggdag::dagify() object
#'The benefit of this function is that it automatically identifies exposure, outcome, confounder, observed and latent variables inputted from dagitty.net, whereas ggdag::tidy_dagitty only does this for ggdag::dagify objects.
#'Output can be used with ggdag to create better looking DAGs from dagitty.net code.
#'
#' @param dag dagitty object
#' @param dag_df a dag object converted to data frame using the ggdag::tidy_dagitty() function, or similar.
#' @param labels vector of labels for nodes in dagify object
#' @return dagify DAG as a dataframe for use with ggdag to create better looking DAGs
#' @noRd
add_labels <- function(dag, dag_df, labels){
  .datatable.aware <- TRUE

  # Labeling variables

  dag_df$label <- sapply( seq_along( dag_df$name ),
                          function(x) if( dag_df$name[x] %in% attr(labels, "names")) attr(labels[ dag_df$name[x] ], "names") )

  node_roles <- get_roles(dag)

  dag_df$role <- sapply( seq_along(dag_df$name), function(x){

    as.vector( unlist( sapply( seq_along(node_roles),
                               function(n) if( dag_df$name[x] %in% node_roles[[n]] ) names(node_roles[n]) ) ) )
  } )

  dag_df$role <- lapply(dag_df$role, function(x) if( is.null(x)) NA else x)

  dag_df$role <- unlist(dag_df$role)

  return(dag_df)
}


#' Generate new dag coordinates
#'
#' @param dag Dagitty object, with at least a treatment and outcome set. Nodes are automatically placed in the following categories: treatment, outcome, confounder, mediator, latent variable, mediator-outcome-confounder, or instrumental variable.
#' @param coords_spec Adjusts node placement, a higher value increases volatility and results in more extreme DAG structures.
#' @param confounders Vector of confounders in order of occurrence.
#' @return dagitty object with coordinates.
#' @examples
#' dag <- add_coords(dag, coords_spec = 3) # update dagitty object node coordinates
#'
#' plot_dagitty(dag) # check coordinates
#'
#' @export
add_coords_helper <- function(dag,
                       coords_spec = 0.1){

  edges <- get_roles(dag)

  confounders <-  edges$confounder

  mediators <-  edges$mediator

  instrumental_variables <- edges$instrumental

  mediator_outcome_confounders <- edges$mediator_outcome_confounder

  competing_exposures <- edges$competing_exposure

  colliders <- edges$collider

  latent_variables <- edges$latent

  observed <- edges$observed

  treatments <- edges$treatment

  outcomes <- edges$outcome

  dag <- add_coordinates(dag,
                         treatments = treatments,
                         outcomes = outcomes,
                         confounders = confounders,
                         mediators = mediators,
                         instrumental_variables = instrumental_variables,
                         mediator_outcome_confounders = mediator_outcome_confounders,
                         competing_exposures = competing_exposures,
                         latent_variables = latent_variables,
                         observed = observed,
                         colliders = colliders,
                         na.omit(coords_spec))

  return(dag)


}


#' convert pdag edges to dag
#'
#' pdag_to_dag_edges()
#'
#' @importFrom data.table as.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object, must use dag{} instead of pdag{}.
#' @returns A data frame of edges.
#' @noRd
pdag_to_dag_edges <- function(dag){
  .datatable.aware <- TRUE
  edges <- data.table::as.data.table(dagitty::edges(dag))[, c("v", "e", "w")]

  edges <- pdag_to_dag_edges_helper(edges)

  return(edges)

}

#' convert pdag edges to dag edges
#'
#' pdag_to_dag_edges_helper()
#'
#' @importFrom data.table as.data.table
#' @importFrom dagitty edges
#' @param edges A dagitty object, must use dag{} instead of pdag{}.
#' @returns A data frame of edges.
#' @noRd
pdag_to_dag_edges_helper <- function(edges){
  .datatable.aware <- TRUE
  # check if any bi-directional edges are "--" instead of "<->" (pdag)
  if( any( edges$e == "--" | edges$e == "<->" ) ){

    pdag_edges <- edges[ ( edges$e == "--" | edges$e == "<->"), ]
    dag_edges <- edges[ !( edges$e == "--" | edges$e == "<->" ), ]

    pdag_edges <- data.frame( v = c(pdag_edges$v, pdag_edges$w),
                              e = "->",
                              w = c(pdag_edges$w, pdag_edges$v) )

    edges <- rbind(dag_edges, pdag_edges)

  }

  return(edges)

}


#' convert pdag to dag
#'
#' pdag_to_dag()
#'
#' @importFrom data.table as.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object, must use dag{} instead of pdag{}.
#' @returns A data frame of edges.
#' @noRd
pdag_to_dag <- function(dag){
  .datatable.aware <- TRUE
  edges <- data.table::as.data.table(dagitty::edges(dag))[, c("v", "e", "w")]

  # check if any bi-directional edges are "--" instead of "<->" (pdag)
  if( any( edges$e == "--" | edges$e == "<->" ) ){

    pdag_edges <- edges[ ( edges$e == "--" | edges$e == "<->"), ]
    dag_edges <- edges[ !( edges$e == "--" | edges$e == "<->" ), ]

    pdag_edges <- data.frame( v = c(pdag_edges$v, pdag_edges$w),
                              e = "->",
                              w = c(pdag_edges$w, pdag_edges$v) )

    edges <- rbind(dag_edges, pdag_edges)

    dag <- rebuild_dag(dag, edges)

    return(dag)

  }

  return(dag)

}

#' detect simultaneous edges to combine arrows where possible
#'
#' handle_simultaneous_edges()
#'
#' @importFrom data.table as.data.table
#' @param edges Dataframe of graph edges.
#' @returns A data frame of edges.
#' @noRd
handle_simultaneous_edges <- function(edges){
  .datatable.aware <- TRUE

  edges_vw <- paste0(edges$v, edges$w)
  edges_wv <- paste0(edges$w, edges$v)
  check_simultaneous <- edges_vw %in% edges_wv

  # check if any bi-directional edges are "--" instead of "<->" (pdag)
  if( length(edges$e[check_simultaneous]) > 0 & all(edges$e[check_simultaneous] == "->") ){

    edges_bidirected <- edges
    edges_bidirected$e <- "<->"

    directed_edges <- edges[ !check_simultaneous, ]
    simultaneous_edges <- edges_bidirected[ check_simultaneous,]

    # for all [x,v] if [y,v] equals [x,w] and if [y,w] == [x,v] delete [x,w]
    x <- 1
    for(x in 1:(nrow(simultaneous_edges)-1)){
      y <- 2
      for(y in 2:nrow(simultaneous_edges)){
        if( simultaneous_edges[y,"v"] == simultaneous_edges[x,"w"] &
            simultaneous_edges[y,"w"] == simultaneous_edges[x,"v"]
        ){
          simultaneous_edges[x,"v"] <- "<@RM@>"
          simultaneous_edges[x,]
          x <- x + 1
          y <- x + 1
        }else{
          y <- y + 1
        }
      }
      x <- x + 1
    }

    simultaneous_edges <- simultaneous_edges[ !simultaneous_edges$v == "<@RM@>", ]

    edges <- rbind(directed_edges, simultaneous_edges)

    return(edges)

  }

  return(edges)

}



#' node names in path from treatment to outcome
#'
#' @importFrom ggdag dag_paths
#' @importFrom data.table as.data.table data.table
#' @param dag A dagitty object
#' @param output_list TRUE or FALSE to output a list (default TRUE returns a list)
#' @returns Vector of nodes in path from treatment to outcome
#' @noRd
find_where_need_proxy <- function(dag, output_list = TRUE){
  treatments <- dagitty::exposures(dag)
  outcomes <- dagitty::outcomes(dag)
  latent_variables <- dagitty::latents(dag)

  mediators_list <- nodes_between_treatment_and_outcome(dag,
                                                            treatments = treatments,
                                                            outcomes = outcomes,
                                                            output_list = output_list)
  # outcome_parents
  #outcome_parents <- dagitty::parents(dag, outcomes)
  latents_need_proxies <- lapply(1:length(treatments), function(x){

    latents_need_proxies <- lapply(1:length(outcomes), function(y){

      mediators_list[[x]][[y]] <- mediators_list[[x]][[y]][
        mediators_list[[x]][[y]] %in% dagitty::parents(dag, outcomes[y]) ]

      mediators_list[[x]][[y]] <- mediators_list[[x]][[y]][
        mediators_list[[x]][[y]] %in% latent_variables  ]

    })
    names(latents_need_proxies) <- outcomes
    latents_need_proxies
  })

  names(latents_need_proxies) <- treatments


  return(latents_need_proxies)
}

#' Label DAG nodes helper function
#'
#' @importFrom ggdag dag_label
#' @param names_list A list of node names.
#' @return A list of shortened node names.
#' @noRd
get_label_helper <- function(names_list){

  ## mostly hard coded (for now) but it should give decent results for common agricultural inputs e.g. location = loc, temperature = temp
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


#' node names in path from treatment to outcome
#'
#' @importFrom dagitty paths
#' @importFrom data.table as.data.table data.table
#' @param dag A dagitty object
#' @param treatments A vector of treatment node names
#' @param outcomes A vector of outcome node names
#' @returns A nested list of potential intsrumental variables within outcomes, and treatments.
#' @noRd
instrumental_variables_helper <- function(dag, treatments, outcomes){

  ## dag node names
  node_names <- names(dag)
  ## latent variables
  latent_vars <- dagitty::latents(dag)
 ## colliders
  colliders <- colliders(dag)
  # outcome_parents
  outcome_parents <- lapply(1:length(outcomes), function(y){

    outcome_parents <- dagitty::parents(dag, outcomes[y])
  })
  outcome_parents_children <- lapply(1:length(outcomes), function(x){

    dagitty::children(dag, outcome_parents[[x]])

  })
  # outcome_latent_parents_children
  outcome_latent_parents <- lapply(1:length(outcomes), function(y){

    outcome_latent_parents <- outcome_parents[[y]][ outcome_parents[[y]] %in% latent_vars ]
  })
  outcome_latent_parents_children <- lapply(1:length(outcomes), function(x){

    dagitty::children(dag, outcome_latent_parents[[x]])

  })
  outcome_latent_parents_parents <- lapply(1:length(outcomes), function(x){

    dagitty::parents(dag, outcome_latent_parents[[x]])

  })


  nodes_trt_to_outcome <- nodes_between_treatment_and_outcome(dag, treatments = treatments, outcomes = outcomes, output_list = TRUE)

  instrumental_variables <- list()
  associated_with_treatment <- list()

  instrumental_variables <- lapply(1:length(treatments), function(x){

    associated_with_treatment[[x]] <- lapply(1:length(node_names), function(n){

      parents_associated_with_treatment <- unique(dagitty::paths(dag, treatments[x], node_names[n], directed = TRUE)$open == TRUE)
      children_associated_with_treatment <- unique(dagitty::paths(dag, node_names[n], treatments[x], directed = TRUE)$open == TRUE)

      associated_with_treatment <- c(parents_associated_with_treatment, children_associated_with_treatment)

      children_associated_with_outcome <- length(unique(dagitty::paths(dag, outcomes[x], node_names[n], directed = TRUE)$open == TRUE))>0

      associated_with_treatment <- length(associated_with_treatment)>0 & children_associated_with_outcome==FALSE

    })

    names(associated_with_treatment[[x]]) <- node_names
    associated_with_treatment[[x]] <- Filter(function(m) m==TRUE,
                                             associated_with_treatment[[x]] )

    associated_with_treatment[[x]] <- names(associated_with_treatment[[x]])

    associated_with_treatment[[x]] <- associated_with_treatment[[x]][ !associated_with_treatment[[x]] %in% latent_vars &
                                                                        !associated_with_treatment[[x]] %in% colliders  &
                                                                        !associated_with_treatment[[x]] %in% outcomes  &
                                                                        !associated_with_treatment[[x]] %in% treatments ]

    instrumental_variables[[x]] <- lapply(1:length(outcomes), function(y){

      nodes_trt_to_outcome_parents <- dagitty::parents(dag, nodes_trt_to_outcome[[x]][[y]])
      parents_of_nodes_trt_to_outcome_parents <- dagitty::parents(dag, nodes_trt_to_outcome_parents)
      nodes_trt_to_outcome_children <- dagitty::children(dag, nodes_trt_to_outcome[[x]][[y]])

      instrumental_variables[[y]] <- associated_with_treatment[[x]][  !associated_with_treatment[[x]] %in% outcome_parents[[y]] &
                                                                        !associated_with_treatment[[x]] %in% outcome_parents_children[[y]] &
                                                                        !associated_with_treatment[[x]] %in% outcome_latent_parents_children[[y]] &
                                                                        !associated_with_treatment[[x]] %in% outcome_latent_parents_parents[[y]] &
                                                                        !associated_with_treatment[[x]] %in% nodes_trt_to_outcome_children &
                                                                        !associated_with_treatment[[x]] %in% nodes_trt_to_outcome_parents &
                                                                        !associated_with_treatment[[x]] %in% parents_of_nodes_trt_to_outcome_parents &
                                                                        !associated_with_treatment[[x]] %in% nodes_trt_to_outcome[[x]][[y]] ]
    })

  })
  names(instrumental_variables) <- treatments
  instrumental_variables

  i <- 1
  for(i in 1:length(treatments)){ ## only way I could get this to work
    names(instrumental_variables[[i]]) <- outcomes ## I wish it wasn't
    i <- i + 1 ## I really wish it didn't have to be this way
  }

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
#' @returns Vector or list of  nodes in the path from treatment to outcome.
#' @export
nodes_between_treatment_and_outcome <- function(dag,
                                                treatments,
                                                outcomes,
                                                output_list = FALSE
                                                ){
  .datatable.aware <- TRUE
  nodes <- c()
  treatments <- unlist(treatments)
  outcomes <- unlist(outcomes)

  if( length(treatments) > 0 & length(outcomes) > 0 ){

    paths_trt_to_y <- lapply(1:length(treatments), function(x){

      paths_trt_to_y <- lapply(1:length(outcomes), function(y){

        tryCatch({
          paths_trt_to_y <- ggdag::dag_paths(dag,
                                             from = treatments[x],
                                             to = outcomes[y],
                                             directed = TRUE,
                                             paths_only = TRUE)[["data"]]
          paths_trt_to_y <- as.vector(unlist(unique( paths_trt_to_y[ complete.cases(paths_trt_to_y["direction"]), "name" ])))
          paths_trt_to_y <- paths_trt_to_y[ !paths_trt_to_y %in% treatments]
        }, warning = function(w){
          paths_trt_to_y[[x]][y] <- NA
        }, error = function(e){
          paths_trt_to_y[[x]][y] <- NA
        }
        )


      })

      names(paths_trt_to_y) <- outcomes
      paths_trt_to_y
    })

    names(paths_trt_to_y) <- treatments


    if(output_list == TRUE){

      return(paths_trt_to_y)

    }

    paths_trt_to_y <- as.vector(unlist(unique(paths_trt_to_y)))

    paths_trt_to_y <- paths_trt_to_y[ complete.cases(paths_trt_to_y) ]

  }

  return(paths_trt_to_y)

}

