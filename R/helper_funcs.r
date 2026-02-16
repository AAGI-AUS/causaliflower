#' node names in path from treatment to outcome
#'
#' @importFrom ggdag dag_paths
#' @importFrom data.table as.data.table data.table
#' @param dag A dagitty object
#' @param treatments A vector of treatment node names
#' @param outcomes A vector of outcome node names
#' @returns Vector of nodes in path from treatment to outcome
#' @noRd
get_nodes_between_treatment_and_outcome <- function(dag, treatments, outcomes){
  .datatable.aware <- TRUE
  nodes <- c()

  if( length(treatments) > 0 & length(outcomes) > 0 ){

    paths_trt_to_y <- lapply(1:length(treatments), function(x){

      paths_trt_to_y <- lapply(1:length(outcomes), function(y){

        paths_trt_to_y <- na.omit(ggdag::dag_paths(dag,
                                                   from = treatments[x],
                                                   to = outcomes[y],
                                                   directed = TRUE,
                                                   paths_only = TRUE)[["data"]])

        #paths_trt_to_y <- paths_trt_to_y[!grepl(treatments[x], paths_trt_to_y$to),]
      })

    })

    paths_trt_to_y <- data.table::as.data.table( do.call( rbind, unlist(paths_trt_to_y, recursive = FALSE) ) )

    nodes <- as.vector(unique(paths_trt_to_y$name))

    nodes <- unique( c(nodes, as.vector(na.omit(unique(paths_trt_to_y$to)))) )

    nodes <- nodes[ !nodes %in% treatments]

  }

  return(nodes)

}

#' Get latent variable names from list
#'
#' get_latent_vec() is a helper function for buildGraph().
#'
#' @importFrom data.table data.table
#' @param instrumental_variables Inputted list or vector of instrumental variables.
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
extract_instrumental_variables <- function(dag, treatments, outcomes, latent_vars, colliders, treatment_parents, treatment_children, mediator_parents, latent_children, outcome_children, outcome_parents){
  .datatable.aware <- TRUE
  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  nodes_trt_to_y_parents <- dagitty::parents(dag, nodes_trt_to_y)

  # instrumental variables
  latent_parents <- dagitty::parents(dag, latent_vars)

  collider_children <- dagitty::children(dag, colliders)

  exclude_as_instruments <- treatment_parents[

    treatment_parents %in% treatment_children | # not allowed as an instrument
      treatment_parents %in% mediator_parents | # not allowed as an instrument
      treatment_parents %in% latent_parents | # not allowed as an instrument
      treatment_parents %in% latent_children | # not allowed as an instrument
      treatment_parents %in% outcome_children | # cyclic relationship
      treatment_parents %in% outcome_parents | # cyclic relationship
      treatment_parents %in% collider_children | # cyclic relationship
      treatment_parents %in% latent_vars | # exclude latent variables
      treatment_parents %in% nodes_trt_to_y_parents # in path treatment -> ... -> outcome

  ]

  instrumental_vars <- as.data.frame( treatment_parents[!treatment_parents %in% exclude_as_instruments] )
  # instrumental_vars <- rbind(instrumental_vars, instrumental_vars)
  if(length(unlist(instrumental_vars)) > 0){

    instrumental_vars$role <- "instrumental"

  }

  if(length(unname(as.vector(instrumental_vars[,1]))) > 0){


    #instrumental_children <- dagitty::children(dag, unname(as.vector(instrumental_vars[,1]))[1])

    # length(instrumental_children[instrumental_children %in% treatments]) > 1


    instrumental_children_multiple_treatments <- lapply( 1:length(unname(as.vector(instrumental_vars[,1]))), function(x){

      instrumental_children <- list( node = list( as.character(unname(instrumental_vars[,1]))[x], dagitty::children(dag, unname(as.vector(instrumental_vars[,1]))[x]) ) )

    })


    instrumental_children_multiple_treatments <- lapply( 1:length(unname(as.vector(instrumental_vars[,1]))), function(x){

      instrumental_children_multiple_treatments[[x]][["node"]][[1]][ length( instrumental_children_multiple_treatments[[x]][["node"]][[2]][ instrumental_children_multiple_treatments[[x]][["node"]][[2]]  %in% treatments ] ) > 1 ]

    })

    instrumental_children_multiple_treatments <- unlist(instrumental_children_multiple_treatments)

  }else{

    instrumental_children_multiple_treatments <- NULL

  }


  if( length( instrumental_children_multiple_treatments ) > 0  ){

    #warning("Potential instrumental variable(s) connected to multiple treatments. These will be ignored, if intended to be used as instruemnts please check edges and try again.
    #     \nIf inputting multiple instruments/treatments to buildGraph(), supply a nested list to 'instrumental_variables' indicating the descendant treatment,
    #     e.g., list( list('instrument 1', list('treatment 1')), list('instrument 2', list('treatment 2')) ) ")


    instrumental_vars$role[ instrumental_vars[,1] %in% instrumental_children_multiple_treatments] <- "observed"

  }

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

  edges <- pdag_to_dag_edges(dag)

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

  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

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
                                  !treatment_parents %in% treatments ] )# proxy_b

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  proxy <- c(proxy_b, proxy_c)

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children &
                                   !outcome_children %in% latent_vars &
                                   !outcome_children %in% outcomes &
                                   !outcome_children %in% mediators &
                                   !outcome_children %in% moc &
                                   !outcome_children %in% confounders ]

  # instrumental variables
  instrumental_vars <- extract_instrumental_variables(dag,
                                                      treatments,
                                                      outcomes,
                                                      latent_vars,
                                                      colliders,
                                                      treatment_parents,
                                                      treatment_children,
                                                      mediator_parents,
                                                      latent_children,
                                                      outcome_children,
                                                      outcome_parents)


  instrumental_vars <- instrumental_vars[,1][instrumental_vars$role ==  "instrumental"]


  # filter latent variables
  latent_vars <- latent_vars[ !latent_vars %in% treatments &
                              !latent_vars %in% outcomes &
                              !latent_vars %in% colliders ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% proxy &
                         !all_vars %in% competing_exposure &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]


  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, v] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, v] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, v] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, v] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, v] %in% mediators &
                             !edges[, v] %in% moc ]  <- "mediator"

  edges$ancestor_iv[edges[, v] %in% instrumental_vars] <- "instrumental"

  edges$ancestor_competing_exposure[edges[, v] %in% competing_exposure &
                                      !edges[, v] %in% proxy_c ] <- "competing_exposure"

  edges$ancestor_proxy[ ( edges[, v] %in% proxy_b &
                            !edges[, v] %in% confounders ) | # "proxy_b"
                          ( edges[, v] %in% proxy_c &
                              !edges[, v] %in% confounders &
                              !edges[, v] %in% mediators ) ] <- "proxy"

  edges$ancestor_collider[edges[, v] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, w] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, w] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, w] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, w] %in% moc &
                         !edges[, w] %in% mediators]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, w] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, w] %in% instrumental_vars] <- "instrumental"

  edges$descendant_competing_exposure[edges[, w] %in% competing_exposure &
                                        !edges[, w] %in% proxy_c ] <- "competing_exposure"

  edges$descendant_proxy[ ( edges[, w] %in% proxy_b &
                              !edges[, w] %in% confounders ) | # "proxy_b"
                            ( edges[, w] %in% proxy_c &
                                !edges[, w] %in% confounders &
                                !edges[, w] %in% mediators ) ] <- "proxy"

  edges$descendant_collider[edges[, w] %in% colliders] <- "collider"

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

  edges <- pdag_to_dag_edges(dag)

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

  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

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
                                                    !treatment_parents %in% treatments ] )# proxy_b

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  proxy <- c(proxy_b, proxy_c)

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children ]

  # instrumental variables
  instrumental_vars <- extract_instrumental_variables(dag,
                                                      treatments,
                                                      outcomes,
                                                      latent_vars,
                                                      colliders,
                                                      treatment_parents,
                                                      treatment_children,
                                                      mediator_parents,
                                                      latent_children,
                                                      outcome_children,
                                                      outcome_parents)


  instrumental_vars <- instrumental_vars[,1][instrumental_vars$role ==  "instrumental"]


  # filter latent variables
  latent_vars <- latent_vars[  !latent_vars %in% treatments &
                                 !latent_vars %in% outcomes ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% proxy &
                         !all_vars %in% competing_exposure &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]


  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, v] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, v] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, v] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, v] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, v] %in% mediators ]  <- "mediator"

  edges$ancestor_iv[edges[, v] %in% instrumental_vars] <- "instrumental"

  edges$ancestor_competing_exposure[edges[, v] %in% competing_exposure &
                                      !edges[, v] %in% proxy_c ] <- "competing_exposure"

  edges$ancestor_proxy[ ( edges[, v] %in% proxy_b &
                            !edges[, v] %in% confounders ) | # "proxy_b"
                          ( edges[, v] %in% proxy_c &
                              !edges[, v] %in% confounders &
                              !edges[, v] %in% mediators ) ] <- "proxy"

  edges$ancestor_collider[edges[, v] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, w] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, w] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, w] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, w] %in% moc ]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, w] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, w] %in% instrumental_vars] <- "instrumental"

  edges$descendant_competing_exposure[edges[, w] %in% competing_exposure &
                                        !edges[, w] %in% proxy_c ] <- "competing_exposure"

  edges$descendant_proxy[ ( edges[, w] %in% proxy_b &
                              !edges[, w] %in% confounders ) | # "proxy_b"
                            ( edges[, w] %in% proxy_c &
                                !edges[, w] %in% confounders &
                                !edges[, w] %in% mediators ) ] <- "proxy"

  edges$descendant_collider[edges[, w] %in% colliders] <- "collider"

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

  edges_ancestors <- edges[,1:14]
  edges_descendants <- edges[,c(1:3,15:25)]

  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:14), idvar = "id",
                                      v.names = "role_ancestor", direction = "long")[,c("v", "e", "w", "role_ancestor", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), 1:4]

  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:14), idvar = "id",
                                        v.names = "role_descendant", direction = "long")[,c("v", "e", "w", "role_descendant", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), 1:4]


  if( nrow(edges_ancestors) != nrow(edges_descendants) ){ # finds missing latent edges

    edges_ancestors <- find_missing_edges(edges_ancestors, edges_descendants)

  }


  edges <- do.call(cbind, lapply(list(edges_ancestors, edges_descendants[,4]), data.table::setDT))

  names(edges)[1:3] <- c("ancestor", "edge", "descendant")

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
  edges_ancestors <- edges_wide[,1:14]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:14), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), c(1,3,4)]
  names(edges_ancestors)[1:2] <- c("ancestor", "descendant")

  edges_ancestors$role_node_id <- "ancestor"

  ## descendant node edges to list ##
  edges_descendants <- edges_wide[,c(1:3,15:25)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:14), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), c(1,3,4)]
  names(edges_descendants)[1:2] <- c("ancestor", "descendant")

  edges_descendants$role_node_id <- "descendant"

  edges_long <- rbind(edges_ancestors, edges_descendants)

  ## pivot long format

  edges_long <- reshape(edges_long,
                        idvar = c("ancestor", "descendant", "role"),
                        timevar = "role_node_id",
                        direction = "wide",
                        v.names = "role")



  return(edges_long)

}


#' Extract ancestor node edges
#'
#' ancestor_edges() is a dagitty wrapper function that returns a list of node names and edges in a dagitty object.
#'
#' @importFrom dagitty edges
#' @param dag A dagitty object.
#' @returns Named list of edges.
#' @noRd
ancestor_edges <- function(dag){
  .datatable.aware <- TRUE

  edges <- pdag_to_dag_edges(dag)

  ## group by nodes
  unique_ancestors <- unique( edges[,"v"] )
  num_unique_ancestors <- nrow(unique_ancestors)


  edges_list <- list()

  # edges_to_assess grouped by each unique node in a list
  edges_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges[ unlist(edges[,"v"]) %in% unlist(unique_ancestors[x]), ]


  }) )

  names(edges_list) <- unlist(unique_ancestors)

  return(edges_list)

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
                             type,
                             temporal_reference_node = NA
){
  .datatable.aware <- TRUE

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

  nodes_ordered <- sort( unlist( dagitty::topologicalOrdering(dag) ) ) # ggdag estimated temporal order of new nodes

  if( length( node_role) > 1 | length( node_role) == 0){

    stop("add_nodes() currently only supports single node_role character inputs.")

  }else if( node_role %in% "confounder" ){
    confounder_occurrance <- as.numeric(order(match(nodes, nodes)))

    ## confounder edges ##
    new_edges <- draw_confounder_edges(type,
                                       outcomes,
                                       treatments,
                                       confounder_vec = nodes, # new nodes as confounder
                                       m_o_confounder_vec,
                                       mediator_vec,
                                       latent_vec,
                                       confounder_occurrance)

    # connect all confounders (fully connected or saturated graph type)
    if( type == "full" | type == "saturated" ){

      confounder_list <- list()

      confounder_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

          list( c( ancestor = nodes[x], edge = "->", descendant = confounder_vec[y]) )

        })

      }) )

      confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
      confounder_unlist <- data.table::as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, confounder_unlist)


      confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

        confounder_list[x] <- lapply(1:length(nodes), function(y){

          list( c( ancestor = confounder_vec[x], edge = "->", descendant = nodes[y]) )

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

        first_list[x] <- list( c( ancestor = nodes[x], edge = "->", descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), first_list)
      first_unlist <- data.table::as.data.table( do.call( rbind, unlist(first_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, first_unlist)


    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% confounder_vec ][length(confounder_vec) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = "->", descendant = nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), last_list)
      last_unlist <- data.table::as.data.table( do.call( rbind, unlist(last_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, last_unlist)


    }

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% confounder_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "treatment" ){
    ## treatment edges ##
    new_edges <- draw_treatment_edges(type,
                                      outcomes,
                                      treatments = nodes, # new nodes as treatment
                                      confounder_vec,
                                      mediator_vec,
                                      collider_vec)

    # connect all treatments (fully connected or saturated graph type)
    if( type == "full" ){

      treatment_list <- list()

      treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

        treatment_list[x] <- lapply(1:length(nodes), function(y){

          list( c( ancestor = treatments[x], edge = "->", descendant = nodes[y]) )

        })

      }) )

      treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
      treatment_unlist <- data.table::as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, treatment_unlist)



      treatment_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        treatment_list[x] <- lapply(1:length(treatments), function(y){

          list( c( ancestor = nodes[x], edge = "->", descendant = treatments[y]) )

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

        first_list[x] <- list( c( ancestor = nodes[x], edge = "->", descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), first_list)
      first_unlist <- data.table::as.data.table( do.call( rbind, unlist(first_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, first_unlist)


    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% treatments ][length(treatments) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = "->", descendant = nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), last_list)
      last_unlist <- data.table::as.data.table( do.call( rbind, unlist(last_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, last_unlist)


    }


    treatments <- c(treatments, nodes)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% treatments ] ) # existing node names in temporal order

  }else if( node_role %in% "outcome" ){
    ## outcome edges ##
    new_edges <- draw_outcome_edges(type, outcomes = nodes, # new nodes as outcome
                                    collider_vec)

    outcomes  <- c(outcome, nodes)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% outcomes ] ) # existing node names in temporal order

  }else if( node_role %in% "mediator" ){
    ## mediator edges ##
    new_edges <- draw_mediator_edges(type,
                                     outcomes,
                                     mediator_vec = nodes, # new nodes as mediator
                                     latent_vec)

    if(type == "full"){

      mediator_list <- list()

      mediator_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        mediator_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = nodes[x], edge = "->", descendant = mediator_vec[y]) )

        })

      }) )

      mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
      mediator_unlist <- data.table::as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, mediator_unlist)


      mediator_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        mediator_list[x] <- lapply(1:length(nodes), function(y){

          list( c( ancestor = mediator_vec[x], edge = "->", descendant = nodes[y]) )

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

        first_list[x] <- list( c( ancestor = nodes[x], edge = "->", descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), first_list)
      first_unlist <- data.table::as.data.table( do.call( rbind, unlist(first_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, first_unlist)


    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% mediator_vec ][length(mediator_vec) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = "->", descendant = nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), last_list)
      last_unlist <- data.table::as.data.table( do.call( rbind, unlist(last_list, recursive = FALSE) ) )

      new_edges <- rbind(new_edges, last_unlist)


    }

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% mediator_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "mediator_outcome_confounder" ){
    ## mediator_outcome_confounder edges ##
    new_edges <- draw_moc_edges(type,
                                treatments,
                                outcomes,
                                confounder_vec,
                                m_o_confounder_vec = nodes, # new nodes as mediator_outcome_confounder
                                mediator_vec,
                                latent_vec)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% m_o_confounder_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "instrumental" ){
    ## instrumental_variables edges ##
    new_edges <- draw_iv_edges(instrumental_variables = nodes, # new nodes as instrumental
                               treatments)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% instrumental_variables ] ) # existing node names in temporal order

  }else if( node_role %in% "competing_exposure" ){
    ## competing_exposure edges ##
    new_edges <- draw_competing_exposure_edges(outcomes,
                                               competing_exposure_vec = nodes) # new nodes as competing_exposure

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% competing_exposure_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "collider" ){
    ## outcome edges ##
    new_edges <- draw_outcome_edges(type,
                                    outcomes,
                                    collider_vec = nodes) # new nodes as collider

    # connect colliders
    treatment_list <- list()

    treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

      treatment_list[x] <- lapply(1:length(nodes), function(y){

        list( c( ancestor = treatments[x], edge = "->", descendant = nodes[y]) )

      })

    }) )

    treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
    treatment_unlist <- data.table::as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

    new_edges <- rbind(new_edges, treatment_unlist)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% outcomes ] ) # existing node names in temporal order

  }else if( node_role %in% "latent" ){
    ## latent_variables edges ##
    new_edges <- draw_latent_edges(latent_variables = nodes) # new nodes as latent

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% latent_vec ] ) # existing node names in temporal order

    latent_vec <- c(latent_vec, nodes)

  }else if( node_role %in% "observed" ){
    ## connect observed to ancestors and descendants ##
    new_edges <- draw_observed_edges(observed = nodes, # new nodes as observed
                                     existing_dag = dag)

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



  tryCatch({
    new_coordinates <- renew_coords(dag = dag,
                                    new_node_names = new_node_names,
                                    coordinates = coordinates,
                                    coords_spec = coords_spec[ complete.cases(coords_spec) ] )

    dagitty::coordinates(dag) <- new_coordinates

  }, warning = function(w){
    message(paste("Warning:", w, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords(dag, coords_spec = coords_spec[ complete.cases(coords_spec) ] )
    return(dag)

  }, error = function(e){
    message(paste("Error:", e, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords(dag, coords_spec = coords_spec[ complete.cases(coords_spec) ] )
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

  node_name_and_coords_vec <- c()

  if( length(latent_vec) > 0 ){ ########### simplify this section ##############

    if( all(complete.cases(latent_vec)) ) {

      node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                    paste(outcomes, " [outcome] ", sep=""),
                                    paste(latent_vec, " [latent] ", sep=""),
                                    paste(node_names, collapse=" "))

    }else{
      node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                    paste(outcomes, " [outcome] ", sep=""),
                                    paste(node_names, collapse=" "))
    }

  }else{

    node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                  paste(outcomes, " [outcome] ", sep=""),
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

  if(all(!is.na(unlist(coordinates)))){

    dagitty::coordinates(dag) <- coordinates

  }


  tryCatch({
    new_coordinates <- renew_coords(dag = dag,
                                new_node_names = new_node_names,
                                coordinates = coordinates,
                                coords_spec = coords_spec[ complete.cases(coords_spec) ] )
    dagitty::coordinates(dag) <- new_coordinates

  }, warning = function(w){
    message(paste("Warning:", w, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords(dag, coords_spec = coords_spec[ complete.cases(coords_spec) ] )
    return(dag)

  }, error = function(e){
    message(paste("Error:", e, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords(dag, coords_spec = coords_spec[ complete.cases(coords_spec) ] )
    return(dag)

  }, finally = return(dag)
  )

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

  edges <- data.table::as.data.table(dagitty::edges(dag))[, c("v", "e", "w")]

  # check if any bi-directional edges are "--" instead of "<->" (pdag)
  if( any( edges$e == "--" | edges$e == "<->" ) ){

    pdag_edges <- edges[ ( edges$e == "--" | edges$e == "<->"), ]
    dag_edges <- edges[ !( edges$e == "--" | edges$e == "<->" ), ]

    pdag_edges <- data.frame( v = c(pdag_edges$v, pdag_edges$w),
                              e = "->",
                              w = c(pdag_edges$w, pdag_edges$v) )

    edges <- rbind(dag_edges, pdag_edges)

    return(edges)

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

  node_name_and_coords_vec <- c()

  if( length(latent_vec) > 0 ){ ########### simplify this section ##############

    if( all(complete.cases(latent_vec)) ) {
      node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                    paste(outcomes, " [outcome] ", sep=""),
                                    paste(latent_vec, " [latent] ", sep=""),
                                    paste(node_names, collapse=" "))
    }else{
      node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                    paste(outcomes, " [outcome] ", sep=""),
                                    paste(node_names, collapse=" "))
    }

  }else if( length(outcomes) + length(treatments) == 0 ){

    node_name_and_coords_vec <-  paste(node_names, collapse=" ")


  }else if( length(treatments) == 0 ){

    node_name_and_coords_vec <- c(paste(outcomes, " [outcome] ", sep=""),
                                  paste(node_names, collapse=" "))

  }else if( length(outcomes) == 0 ){

    node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                  paste(node_names, collapse=" "))

  }else{
    node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                  paste(outcomes, " [outcome] ", sep=""),
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

    nodes_parents <- lapply(1:num_nodes, function(x){ # get new node name parents in existing nodes (found in both dags prior to merge)
      nodes_parents <- existing_node_names_not_outcome[ existing_node_names_not_outcome %in% nodes_parents[[x]] ]
    })

    nodes_children <- lapply(1:num_nodes, function(x){ # get new node name parents in existing nodes (found in both dags prior to merge)
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


  quality_check <- FALSE
  iteration <- 1

  time_limit <- num_nodes + num_nodes/lambda

  setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)

  on.exit( {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    } )

  while(quality_check == FALSE){

    iteration <- iteration + ( num_nodes*lambda )

    ## get difference between min-max coords
    diff_y_coords <-  abs( min( existing_coordinates_y ) - max( existing_coordinates_y ) )
    diff_x_coords <-  abs( min( existing_coordinates_x ) - max( existing_coordinates_x ) )


    if( num_nodes_y > 0){
      ## y coordinates ##

      new_y_coords <- sapply(1:num_nodes_y, function(x){

        if( length( nodes_parents_y[[x]] ) > 0 ){


          new_y_coords <- max( existing_coordinates_y[ names(existing_coordinates_y) %in% nodes_parents_y[[x]] ]
          ) + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda)/2 )

        }else if( length( nodes_children_y[[x]] ) > 0 ) {

          new_y_coords <- min( existing_coordinates_y[ names(existing_coordinates_y) %in% nodes_children_y[[x]] ]
          ) - x*iteration - ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) - runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda )/2 )

        }else if( length( unlist(nodes_parents_y) ) > 0 ){

          new_y_coords <- min( existing_coordinates_y[ names(existing_coordinates_y) %in% unlist(nodes_parents_y) ]
          ) + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda )/2 )

        }else if( length( unlist(nodes_children_y) ) > 0 ){

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

        if( length( nodes_children_x[[x]] ) > 0 ) {

          new_x_coords <- min( existing_coordinates_x[ names(existing_coordinates_x) %in% nodes_children_x[[x]] ]
          ) - x*iteration - ( diff_x_coords/num_vars )*num_nodes + diff_x_coords*(x/diff_x_coords) - runif(n = 1,
                                                                                                           min = iteration*lambda,
                                                                                                           max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( nodes_parents_x[[x]] ) > 0 ){

          new_x_coords <- min( existing_coordinates_x[ names(existing_coordinates_x) %in% nodes_parents_x[[x]] ]
          ) + x*iteration + ( diff_x_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( unlist(nodes_children_x) ) > 0 ){

          new_x_coords <- min( existing_coordinates_x[ names(existing_coordinates_x) %in% unlist(nodes_children_x) ]
          ) - x*iteration - ( diff_x_coords/num_vars )*num_nodes + diff_x_coords*(x/diff_x_coords) - runif(n = 1,
                                                                                                           min = iteration*lambda,
                                                                                                           max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( unlist(nodes_parents_x) ) > 0 ){

          new_x_coords <- max( existing_coordinates_x[ names(existing_coordinates_x) %in% unlist(nodes_parents_x) ]
          ) + x*iteration + diff_x_coords*(x/diff_x_coords) + runif(n = 1,
                                                                    min = iteration*lambda,
                                                                    max = iteration + ( num_nodes*lambda*10 )/2 )
        }else{

          new_x_coords <- x + x*iteration + diff_x_coords*(x/diff_x_coords) - runif(n = 1,
                                                                                    min = iteration*lambda,
                                                                                    max = iteration + ( num_nodes*lambda*10 )/2 )
        }
        new_x_coords

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

        nodes_descendants_coords_x <- (
          mean( tidyr::replace_na(is.na(existing_coords$x[ names(existing_coords$x) %in% nodes_descendants[[x]] ]), 0))

          - runif(n = 1, min = ( 0.1*iteration )*( num_vars/num_nodes ),
                  max = ( ( 0.2*iteration )*( num_vars/num_nodes ) ) ) ** (2*(0.1*iteration))

          - (2*( 1 / x )) - ( (0.2*iteration)*( num_vars ) )
        )

        nodes_descendants_coords_y <- (
          mean( tidyr::replace_na(is.na(existing_coords$y[ names(existing_coords$y) %in% nodes_descendants[[x]] ] ), 0))

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


