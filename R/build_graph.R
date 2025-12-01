#' Build a dagitty object
#'
#' build_graph() produces a saturated graph by default from the supplied treatments, outcome, confounder, and  mediators (optional) based on the type of graph specified. Optional latent confounder or medfiator variables can be specified.
#'
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when filter
#' @importFrom dagitty dagitty is.dagitty
#' @param type The type of graph generated. Defaults to 'full', producing a fully connected graph with confounders connected in both directions (bi-directional), and to mediators in one direction (uni-directional). If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param variables Dagitty object, or vector of variable names, e.g. "Z" or c("Z1", "Z2", "Z3"). If variable names are inputted the order determines the assigned coordinates. A list can also be supplied. Variables inputted are treated as confounders. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param treatments Treatment variable name, e.g. "X". Must be specified, unless included in the named vector 'variables'.
#' @param outcomes Outcome variable name, e.g. "Y". Must be specified, unless included in the named vector 'variables'.
#' @param mediators Character or vector of mediator variable names, e.g. "M" or c("M1", "M2", "M3").
#' @param latent_variables Character or vector of additional or already supplied latent (unobserved) variable names, e.g. "U" or c("U1", "U2", "M1").
#' @param instrumental_variables Vector of instrumental variable names, e.g. "IV"
#' @param mediator_outcome_confounders Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures. Setting 'lambda_max' generates a DAG at each lambda value between lambda and lambda_max (only used if iterations is supplied). Iterations controls number of repeats for each lambda value (returns the first lambda value if NA).
#' @returns A dagitty object, fully connected (saturated) graph.
#' @examples
#'
#'
#' # There are three main ways to supply variables to build_graph().
#'
#' # Option 1: Inputting a named vector of variables inputted (can either include treatment and outcome or as separate inputs)
#' {
#'   variables <- c(confounders = c("Z1", "Z2", "Z3"),
#'                  treatments = "X",
#'                  outcomes = "Y",
#'                  mediators = "M",
#'                  instrumental_variables = "IV")
#' }
#'
#' # Three graph types can also be generated. First, we build an 'ordered' graph where the supplied vector order is used to determine the temporal order of confounder nodes.
#'   type <- "ordered"
#'
#'   dag <- build_graph(type = type,
#'                     variables = variables)
#'
#' # Plotting the graph using ggdagitty() assigns coordinates based on the variable roles (confounder, treatments, outcomes, mediator, etc.).
#' ggdagitty(dag)
#'
#' #' # Option 2: The 'confounders' input can be used, while 'variables' is ignored.
#' #           Separate inputs are required for treatments, outcomes, and other nodes.
#' {
#'   confounders <- c("Z1", "Z2", "Z3")
#'   treatments <- "X"
#'   outcomes <- "Y"
#'   instrumental_variables <- "IV"
#'   mediator = "M"
#' }
#'
#' A "saturated" graph type bidirectionally connects each of the confounders, in addition to directed arrows to the outcomes and treatments nodes.
#' type <- "saturated"
#'
#' dag <- build_graph(type = type,
#'                   treatments = treatments,
#'                   outcomes = outcomes,
#'                   confounders = confounders,
#'                   mediators = mediators,
#'                   instrumental_variable = instrumental_variable)
#'
#' ggdagitty(dag)
#'
#' # Option 3: The 'variables' input can be used to connect all variables to each other, treating them as confounders, besides treatments and outcomes and any other separate inputs.
#' #           With this configuration, the actual 'confounders' input can be left blank. I decided to keep this as an option for users who may not have much 'exposure' to causal graphs.
#' {
#'   variables <- c("Z1", "Z2", "Z3")
#'   treatments <- "X"
#'   outcomes <- "Y"
#'   instrumental_variables <- "IV"
#'   mediator = "M"
#' }
#'
#' # Using type = "full" generates a fully connected graph: all confounders are connected to each other, in both directions, and also to mediators, treatments and outcomes.
#' type <- "full"
#'
#' dag <- build_graph(type = type,
#'                   variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes,
#'                   mediators = mediators,
#'                   instrumental_variables = instrumental_variables)
#'
#' ggdagitty(dag)
#'
#' @export
build_graph <- function(type = c("full", "saturated", "ordered"),
                       variables,
                       treatments = NA,
                       outcomes = NA,
                       mediators = NA,
                       latent_variables = NA,
                       instrumental_variables = NA,
                       mediator_outcome_confounders = NA,
                       competing_exposures = NA,
                       colliders = NA,
                       coords_spec = c(lambda = 0, lambda_max = NA, iterations = NA)
                       ){

  # check for an existing dag (inputted as confounders)
  if( dagitty::is.dagitty(variables) ){

    existing_dag <- variables
    #existing_dag <- dag
    node_roles <- get_roles(existing_dag)

    treatments <- node_roles$treatment

    outcomes <- node_roles$outcome

    confounder_vec <- variables <- node_roles$confounder
    confounder_occurrance <- as.numeric(order(match(confounder_vec, confounder_vec)))

    m_o_confounder_vec <- mediator_outcome_confounders <- node_roles$mediator_outcome_confounder

    mediator_vec <- mediators <- node_roles$mediator

    instrumental_variables <- node_roles$instrumental

    collider_vec <- colliders <- node_roles$collider

    competing_exposure_vec <- competing_exposures <- node_roles$competing_exposure

    observed <- node_roles$observed

    latent_vec <- node_roles$latent

    if( all(complete.cases(unlist(latent_vec))) ){

      latent_descendants_list <- lapply(1:length(latent_vec), function(x){

        latent_descendants_list <- dagitty::children(existing_dag, latent_vec[x])

      })

      latent_variables_list <- lapply(1:length(latent_vec), function(x){

        latent_variables_list <- latent_vec[x]

      })

      latent_variables <- lapply(1:length(latent_vec), function(x){

        latent_variables_list[[x]] <- c( latent_variables_list[x], list(as.list(unlist(latent_descendants_list[x]))) )

      })

    }


  }else if( all( complete.cases( unlist(treatments) ) ) & all( complete.cases( unlist(outcomes) ) ) ){
    # no existing dag, use other inputs

    existing_dag <- NA

    mediator_vec <- as.vector( unlist( lapply( mediators, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    latent_vec <- get_latent_vec(latent_variables)

    m_o_confounder_vec <- as.vector( unlist( lapply( mediator_outcome_confounders, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    competing_exposure_vec <- as.vector( unlist( lapply( competing_exposures, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    collider_vec <- as.vector( unlist( lapply( colliders, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    # create confounder_vec & confounder_occurance
    if( length(unlist(variables)) > length(variables) ){
      confounders_list <- list()
      confounders_list <- lapply( 1:length(variables), function(x){

        confounders_list[[x]] <- list( as.data.frame( variables[[x]]), rep(x, times = length(variables[[x]]) ) )
        confounders_list <- dplyr::bind_cols(confounders_list[[x]][1], confounders_list[[x]][2])

      } )

      confounders_df <- dplyr::bind_rows( confounders_list )
      confounder_vec <- as.vector( unlist(confounders_df[1]) )
      confounder_occurrance <- as.vector( unlist(confounders_df[2]) )

    }else{

      confounder_vec <- as.vector(variables)
      confounder_occurrance <- as.numeric(order(match(confounder_vec, confounder_vec)))

    }

    # if any mediator-outcome confounders are also inputted as confounders, execution is stopped

    observed <- NA

  }else{

    stop("The 'treatments' and 'outcomes' inputs should be used if a DAG is not provided in the 'variables' input. Please adjust inputs and try again.")

  }


  if( !dagitty::is.dagitty(variables) & all( confounder_vec %in% m_o_confounder_vec ) != FALSE ){

    stop("Confounder names provided in 'variables' have been detected in the 'mediator_outcome_confounder' input. The 'variables' and 'mediator_outcome_confounder'inputs should be mutually exclusive. Please adjust inputs and try again.")

  }

  edges_df <- draw_edges(type,
                         outcomes,
                         treatments,
                         confounder_vec,
                         m_o_confounder_vec,
                         mediator_vec,
                         instrumental_variables,
                         competing_exposure_vec,
                         latent_vec,
                         latent_variables,
                         collider_vec,
                         observed,
                         confounder_occurrance,
                         existing_dag)

  ## get variable names ##
  node_names <- unique( as.vector( c(confounder_vec, m_o_confounder_vec, mediator_vec, unlist(instrumental_variables), competing_exposure_vec, collider_vec) ) )
  node_names <- Filter(Negate(anyNA), node_names)

  exclude_names <- c(treatments, outcomes, latent_vec)
  exclude_names <- exclude_names[ complete.cases(exclude_names) ]
  node_names <- node_names[!node_names %in% exclude_names]


  dag <- construct_graph(edges_df, node_names, treatments, outcomes, latent_vec)

  if( !all(complete.cases(existing_dag)) ){

    dag <- add_coordinates(dag = dag,
                           treatments = treatments,
                           outcomes = outcomes,
                           confounders = variables,
                           mediators = mediators,
                           instrumental_variables = instrumental_variables,
                           mediator_outcome_confounders = mediator_outcome_confounders,
                           competing_exposures = competing_exposures,
                           latent_variables = latent_variables,
                           colliders = colliders,
                           lambda = na.omit(coords_spec))

  }else{

    dagitty::coordinates(dag) <- dagitty::coordinates(existing_dag)

  }


  return(dag)

}


#' Merge two dagitty objects
#'
#' merge_graphs() adds new nodes to the first from the second supplied dagitty object, before generating new coordinates for the merged graph.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty edges exposures outcomes latents coordinates
#' @param dag First dagitty object.
#' @param new_dag Second dagitty object, added to the first.
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings use treatment as the temporal point of reference for adding new (post-treatment) nodes, but if other node types are specified it can be useful to specify the existing node name.
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures.
#' @returns A dagitty object
#' @export
merge_graphs <- function(dag,
                         new_dag,
                         temporal_reference_node = NA,
                         coords_spec = 0.1
){
  .datatable.aware <- TRUE

  edges <- as.data.frame( dagitty::edges(dag) )[, c("v", "e", "w")] # get primary dag edges
  node_names <- names(dag) # extract dag node names

  treatments <- dagitty::exposures(dag) # get treatments
  outcomes <- dagitty::outcomes(dag) # get outcomes
  latent_vec <- unlist(dagitty::latents(dag)) # get dag latent variables
  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  new_edges <- as.data.frame( dagitty::edges(new_dag) )[, c("v", "e", "w")]  # get new dag edges
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

  if( all(!is.na(unlist(coordinates))) ){

    dagitty::coordinates(dag) <- coordinates

  }

  nodes_ordered <- sort( unlist( dagitty::topologicalOrdering(new_dag) ) ) # ggdag estimated temporal order of new nodes
  temporal_reference_node <- names( nodes_ordered[1] ) # first node is selected as the temporal reference node

  existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% existing_node_names ] ) # existing node names in temporal order
  new_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% new_node_names ] ) # new node names in temporal order

  new_coordinates <- add_merged_node_coordinates(dag = dag,
                                                 existing_node_names = existing_node_names,
                                                 new_node_names = new_node_names,
                                                 temporal_reference_node = temporal_reference_node,
                                                 nodes_ordered = nodes_ordered,
                                                 coordinates = coordinates,
                                                 coords_spec = coords_spec[ complete.cases(coords_spec) ] )


  dagitty::coordinates(dag) <- new_coordinates


  return(dag)

}


#' Add nodes to a dagitty object
#'
#' add_nodes() allows multiple nodes to be added to a dagitty object based on existing nodes or another dagitty object.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param existing_nodes Vector of existing node names, used as a reference for the new graph nodes, e.g., c("Z1", "Z2", "Z3").
#' @param existing_node_type A suffix added to each of the reference node names, e.g. "pre_treatment", or "t0".
#' @param new_node_type A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings use treatment as the temporal point of reference for adding new (post-treatment) nodes, but if other node types are specified it can be useful to specify the existing node name.
#' @param num_repeats Number of additional copies of nodes, such as time points. Each repeat number is included at the end of new node names (new_new_t1, new_node_t2, etc.).
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures.
#' @returns A dagitty object.
#' @export
add_nodes <- function(dag,
                      existing_nodes,
                      existing_node_type = "pre_treatment",
                      new_node_type = "post_treatment",
                      temporal_reference_node = NA,
                      num_repeats = NA,
                      coords_spec = 1
                     ){

  edges <- as.data.frame( dagitty::edges(dag) )[, c("v", "e", "w")] # get dag edges
  node_names <- names(dag) # extract dag node names

  treatments <- dagitty::exposures(dag) # get treatments
  outcomes <- dagitty::outcomes(dag) # get outcomes
  latent_vec <- unlist(dagitty::latents(dag)) # get dag latent variables
  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  if( dagitty::is.dagitty(existing_nodes) ){   # check for an existing dag (inputted as existing_nodes)

    dag <- merge_dagitty(dag,
                         edges,
                         node_names,
                         existing_nodes,
                         treatments,
                         outcomes,
                         latent_vec,
                         coordinates,
                         coords_spec
                         )

  }else if( length( node_names[ node_names %in% existing_nodes] ) != length( existing_nodes ) ){  # checks that all specified reference nodes are in the dag

    stop("Please check reference node input and try again.")

  }else{ # create new nodes by referencing existing nodes

    dag <- create_new_node_graph(dag,
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
                                 )

  }


  return(dag)
}



#' Add nodes to a dagitty object
#'
#' create_new_node_graph() is a helper function for add_nodes(). It adds new nodes to a dagitty object by referencing existing nodes.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param existing_nodes Vector of existing node names, used as a reference for the new graph nodes, e.g., c("Z1", "Z2", "Z3").
#' @param existing_node_type A suffix added to each of the reference node names, e.g. "pre_treatment", or "t0".
#' @param new_node_type A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings use treatment as the temporal point of reference for adding new (post-treatment) nodes, but if other node types are specified it can be useful to specify the existing node name.
#' @param num_repeats Number of additional copies of nodes, such as time points. Each repeat number is included at the end of new node names (new_new_t1, new_node_t2, etc.).
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures.
#' @returns A dagitty object.
#' @noRd
create_new_node_graph <- function(dag,
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
  new_node_names <- create_new_node_names(existing_nodes, new_node_type, num_repeats)
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

    edges["v"] <- sapply(   as.vector(edges["v"]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    edges["w"] <- sapply( as.vector(edges["w"]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    names(coordinates[["x"]]) <- sapply(   names(coordinates[["x"]]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

    names(coordinates[["y"]]) <- sapply(   names(coordinates[["y"]]), function(x)
      replace( x, x %in% existing_nodes[i], existing_node_names[i] ) )

  }

  dag <- construct_graph(edges, node_names, treatments, outcomes, latent_vec)

  # draw edges for all new nodes (includes connecting parent and children nodes, reference nodes to new nodes, and consecutive new nodes)
  new_node_df <- draw_new_node_edges(dag, new_node_names, existing_node_names, existing_node_type, new_node_type, temporal_reference_node, num_repeats)


  edges <- rbind(edges, new_node_df)
  edges[] <- lapply(edges, as.character)

  dag <- construct_graph(edges, node_names, treatments, outcomes, latent_vec)


  if(all(!is.na(unlist(coordinates)))){

    dagitty::coordinates(dag) <- coordinates

  }


  new_node_names <- as.vector(unlist(new_node_names))

  coordinates <- new_node_coordinates(dag,
                                      existing_node_names,
                                      new_node_names,
                                      existing_node_type,
                                      num_repeats,
                                      num_ref_nodes,
                                      temporal_reference_node,
                                      coordinates,
                                      coords_spec = coords_spec[ complete.cases(coords_spec) ] )


  dagitty::coordinates(dag) <- coordinates

  return(dag)

}

#' construct dag
#'
#' merge_dagitty() is a helper function for add_nodes() and merge_graph a daggity object using edges input.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty
#' @param edges A data frame of edges.
#' @returns A dagitty object
#' @noRd
merge_dagitty <- function(dag,
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

  new_edges <- as.data.frame( # get edges, node names from new daggity object (inputted as existing_nodes)
    dagitty::edges(new_dag) )[, c("v", "e", "w")]

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

  nodes_ordered <- sort( unlist( dagitty::topologicalOrdering(new_dag) ) ) # ggdag estimated temporal order of new nodes
  temporal_reference_node <- names( nodes_ordered[1] ) # first node is selected as the temporal reference node

  existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% existing_node_names ] ) # existing node names in temporal order
  new_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% new_node_names ] ) # new node names in temporal order

  new_coordinates <- add_merged_node_coordinates(dag = dag,
                                      existing_node_names = existing_node_names,
                                      new_node_names = new_node_names,
                                      temporal_reference_node = temporal_reference_node,
                                      nodes_ordered = nodes_ordered,
                                      coordinates = coordinates,
                                      coords_spec = coords_spec[ complete.cases(coords_spec) ] )


  dagitty::coordinates(dag) <- new_coordinates


  return(dag)

}

#' Rebuild dag
#'
#' rebuild_dag() rebuilds a dag using a dagitty object and data frame of edges input.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty isAcyclic
#' @param dag A dagitty object.
#' @param edges A vector of edges.
#' @returns A dagitty object
#' @noRd
rebuild_dag <- function(dag, edges){
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
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty
#' @param edges A data frame of edges.
#' @returns A dagitty object
#' @noRd
construct_graph <- function(edges, node_names, treatments, outcomes, latent_vec){
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

#' create causal dag with hidden nodes
#'
#' When used on a dag containing competing exposures, identified here as variables with only a causal path affecting outcome, causalify() adds a latent node affecting treatment and the competing exposure so that it is included in minimal sufficient ajustment sets when dagitty::adjustmentSets() is called.
#'
#' @importFrom data.table data.table fcase
#' @param dag dagitty object
#' @return Model object of either
#' @noRd
causalify <- function(dag){

  edges <- get_edges(dag)

  node_names <- unique(edges[c("ancestor", "role_ancestor")])

  treatments <- node_names[,"ancestor"][node_names[, "role_ancestor"] %in% "treatment"]
  outcomes <- dagitty::outcomes(dag)
  latents <- dagitty::latents(dag)

  proxy_var_name <- edges$ancestor[edges$role_ancestor == "competing_exposure"]

  new_edges <- lapply( 1:length(proxy_var_name), function(x){

    new_edges <- c(ancestor = as.character(paste("U",proxy_var_name[x],sep="_")),
                                                  edge = as.character("->"),
                                                  descendant = as.character(proxy_var_name[x]),
                                                  role_ancestor = as.character("latent"),
                                                  role_descendant =as.character("proxy_c"))

  })


  new_outcome_edges <- lapply( 1:length(proxy_var_name), function(x){
    new_outcome_edges <- c(ancestor = as.character(paste("U",proxy_var_name[x],sep="_")),
                           edge = as.character("->"),
                           descendant = as.character(treatments),
                           role_ancestor = as.character("latent"),
                           role_descendant =as.character("treatment"))
  })
  edges <- dplyr::bind_rows(edges, new_edges, new_outcome_edges)


  edges[,"role_ancestor"][edges[, "ancestor"] %in% proxy_var_name] <- "proxy"
  edges[,"role_descendant"][edges[, "descendant"] %in% proxy_var_name] <- "proxy"

  node_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "latent")
  edges_list <- lapply( seq_along(1:9),
                        function(x){

                          edges_list <- edges %>% filter( role_ancestor == node_roles[x] )
                          edges_list[,c(1:3)]

                        } )


  node_names <- unique(edges[c("ancestor", "role_ancestor")])
  node_names <- node_names[,1]

  exclude_names <-c(treatments, outcomes, latents)
  node_names <- node_names[!node_names %in% exclude_names]

  latents <- c(unlist(latents), as.vector(unique(edges[,"ancestor"][edges[, "role_ancestor"] %in% "latent"])))

  coordinates <- dagitty::coordinates(dag)
  node_name_and_coords_vec <- c()

  if(length(latents) > 0) {
    node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                  paste(outcomes, " [outcome] ", sep=""),
                                  paste(latents, " [latent] ", sep=""),
                                  paste(node_names, collapse=" "))
  }else{
    node_name_and_coords_vec <- c(paste(treatments, " [exposure] ", sep=""),
                                  paste(outcomes, " [outcome] ", sep=""),
                                  paste(node_names, collapse=" "))
  }

  edges_unlist <- dplyr::bind_rows(edges_list)

  edges_unlist <- lapply(1:nrow(edges_unlist), function(x){
    edges_unlist <- paste(edges_unlist[x,], collapse=" ")
  })
  edges_vec <- paste(unlist(edges_unlist), collapse=" ")

  dag <- paste("dag {", paste(node_name_and_coords_vec, collapse=""), edges_vec, "}", sep = " ")

  dag <- dagitty::dagitty(dag)

  # sets coordinates for latent nodes
  dag <- add_latent_coordinates(dag, latent_variables = latents)


  return(dag)
}

