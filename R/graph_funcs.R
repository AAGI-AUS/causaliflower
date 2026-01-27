#' Build dagitty objects
#'
#' build_graph() produces a dagitty graph object from inputted parameters (e.g. treatments, outcome, confounders).
#'
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom dagitty is.dagitty children coordinates
#' @param type Type of connected graph (e.g. "full", "saturated", "ordered", "none"). Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param variables Dagitty object, or vector of variable names, e.g. "Z" or c("Z1", "Z2", "Z3"). If variable names are inputted the order determines the assigned coordinates. A list can also be supplied. Variables inputted are treated as confounders. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param treatments Treatment variable name, e.g. "X". Must be specified.
#' @param outcomes Outcome variable name, e.g. "Y". Must be specified.
#' @param mediators Character or vector of mediator variable names, e.g. "M" or c("M1", "M2", "M3").
#' @param latent_variables Character or vector of additional or already supplied latent (unobserved) variable names, e.g. "U" or c("U1", "U2", "M1").
#' @param instrumental_variables Vector of instrumental variable names, e.g. "IV"
#' @param mediator_outcome_confounders Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param competing_exposures Vector of competing exposure names. An arrow is drawn connecting competing exposures to the outcome, with other arrows also connected depending on type of graph specified.
#' @param colliders Vector of collider variables, with both treatment and outcome parents.
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures. Setting 'lambda_max' generates a DAG at each lambda value between lambda and lambda_max (only used if iterations is supplied). Iterations controls number of repeats for each lambda value (returns the first lambda value if NA).
#' @returns A dagitty object, fully connected (saturated) graph.
#' @examples
#' ## initial variables (see above for full list of possible inputs)
#' {
#'   variables <- c("Z1", "Z2", "Z3") # these are treated as confounders
#'   treatments <- "X"
#'   outcomes <- "Y"
#' }
#' # Three types of graphs can be generated using build_graph()
#'
#' ## Option 1: 'ordered' graph, where the supplied vector order is used to determine the temporal order of confounder nodes.
#'
#' type <- "ordered"
#'
#' dag <- build_graph(variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes,
#'                   type = type)
#'
#' plot_dagitty(dag) # Plot the graph with plot_dagitty() (assigns coordinates based on the variable roles, e.g., confounders, treatments, outcomes).
#'
#'
#' ## Option 2: 'saturated' graph connecting each of the confounders (inputted as variables)
#'
#' type <- "saturated"
#'
#' dag <- build_graph(variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes,
#'                   type = type)
#'
#' plot_dagitty(dag)
#'
#'
#' ## Option 3: "full" creates a fully connected graph (bidirectional arrows between confounders, and from confounders to all other nodes except treatment and outcome.)
#'
#' type <- "full"
#'
#' dag <- build_graph(variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes,
#'                   type = type)
#'
#' plot_dagitty(dag)
#'
#' @export
build_graph <- function(variables,
                       treatments = NA,
                       outcomes = NA,
                       mediators = NA,
                       latent_variables = NA,
                       instrumental_variables = NA,
                       mediator_outcome_confounders = NA,
                       competing_exposures = NA,
                       colliders = NA,
                       type = "full",
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

  ## get variable names ##
  observed_node_names <- unique( as.vector( c(confounder_vec, m_o_confounder_vec, mediator_vec, competing_exposure_vec, collider_vec, instrumental_variables) ) )
  observed_node_names <- Filter(Negate(anyNA), observed_node_names)

  edges_df <- draw_edges(observed_node_names,
                         type,
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

  ## remove treatments, outcomes and latents from node names ##
  exclude_names <- c(treatments, outcomes, latent_vec)
  exclude_names <- exclude_names[ complete.cases(exclude_names) ]
  node_names <- observed_node_names[!observed_node_names %in% exclude_names]


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


#' Assess dagitty object edges
#'
#' assess_edges() provides ways to assess connected edges based on causal criteria and/or user inputs.
#'
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param edges_to_keep A vector of directed arrows to be kept between (non-treatment and non-outcome) variables, e.g. c("Z1 -> Z2", "Z2 -> Z3"), c("y", "n", "y"), or c(TRUE, FALSE, TRUE).
#' @param assess_causal_criteria Defaults to FALSE. If TRUE, the user is guided through a sequence that assesses each pair of connected nodes using causal criteria. Based on the Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs) from Ferguson et al. (2020).
#' @returns A list or vector of edges.
#' @examples
#' ## initial dag
#' {
#'   variables <- c("Z1", "Z2", "Z3") # these are treated as confounders
#'   treatments <- "X"
#'   outcomes <- "Y"
#'   type <- "ordered"
#'
#'   dag <- build_graph(type = type,
#'                      variables = variables,
#'                      treatments = treatments,
#'                      outcomes = outcomes)
#' }
#'
#' ## Using a saturated graph as an example, created from an existing dag
#' saturated_graph <- build_graph(type = c('full', 'saturated'), # choose a type
#'                                variables = dag) # existing dagitty object inputted
#'
#' ## Option 1: existing dag as edges to keep
#' edges <- assess_edges(saturated_graph, edges_to_keep = dag)
#'
#' ## Option 2: guided causal criteria sequence (rule out edges already included in existing dag)
#' edges_to_keep <- assess_edges(saturated_graph, edges_to_keep = dag,
#'                               assess_causal_criteria = TRUE)
#'
#' @export
assess_edges <- function(dag, edges_to_keep = NA, assess_causal_criteria = FALSE){
  .datatable.aware <- TRUE
  #edges_to_keep <- initial_dag
  #dag <- saturated_graph
  edges_to_assess  <- data.table::as.data.table(dagitty::edges(dag))[, c("v", "e", "w")]


  if( dagitty::is.dagitty(edges_to_keep) ){

    edges_to_keep <- data.table::as.data.table(dagitty::edges(edges_to_keep))[, c("v", "e", "w")]

    edges_to_assess <-  edges_to_assess[!edges_to_keep, on = c("v", "e", "w")]

  }else if( all( complete.cases( unlist(edges_to_keep) ) ) ){

    if(is.vector( edges_to_keep ) ){

      edges_to_keep <- data.table::as.data.table( do.call( rbind, strsplit(edges_to_keep, " ") ) )

      colnames(edges_to_keep) <- c("v", "e", "w")

    }

    if( is.data.frame(edges_to_keep) | data.table::is.data.table(edges_to_keep) ){

      edges_to_assess <-  edges_to_assess[!edges_to_keep, on = c("v", "e", "w")]

    }else{

      stop("Edges_to_keep must be a data frame, data table, vector, or dagitty object. Please check inputs and try again.")

    }

  }



  if( assess_causal_criteria == TRUE ){

    check_skip_sequence <- FALSE

    #cat("\nThere are", nrow(edges_to_assess), "directed arrows to be assessed:", "\n", "\n", sep=" ")
    #print(edges_to_assess, quote=FALSE)
    #cat("\nAssess the posited causal relationships? (ESC-DAGs causal criteria and counterfactual thought experiment sequence)", "\n")

    edges_to_keep <- causal_criteria_sequence(edges_to_assess, check_skip_sequence, edges_to_keep)

    return(edges_to_keep)

  }

  if( nrow(edges_to_assess) != 0 ){


    if( all( complete.cases( unlist(edges_to_keep) ) ) ){

      ## collapse edges_to_keep to a vector
      num_edges <- nrow(edges_to_keep)

      edges_to_keep <- suppressWarnings( sapply(1:num_edges, function(x){

        edges_to_keep <- paste( edges_to_keep[x,], collapse=" ")

      }) )

    }

    # group edges_to_assess by unique node names
    edges_assess_list <- print_edges_helper(edges_to_assess)

    cat( paste("c(", paste( unlist(edges_assess_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

  }else{

    stop("There are no edges to assess. Please check the supplied dagitty object and try again.")


  }



  if( !all( complete.cases(edges_to_keep) ) & assess_causal_criteria == FALSE ){

    edges_to_assess <- suppressWarnings( sapply(1:nrow(edges_to_assess), function(x){

      edges_to_assess <- paste( edges_to_assess[x,], collapse=" ")

    }) )

    message("\nOutputted edges to assess.
            \n\nPrinted edges to assess - copy and paste in a .R file to use as a vector object.")

    return(edges_to_assess)

  }

  edges_list <- list(
    edges_to_keep = edges_to_keep,
    edges_to_assess = edges_to_assess
  )

  message("\nOutputted edges list (edges_to_keep & edges_to_assess).
          \n\nPrinted edges to assess - copy and paste in a .R file to use as a vector object.")

  return(edges_list)

}


#' Remove dagitty object edges
#'
#' keep_edges() removes edges based on user inputs.
#'
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty isAcyclic
#' @param dag A saturated graph dagitty object. Exposure and outcome must be indicated, and optionally can include assigned coordinates.
#' @param edges_to_keep A vector of directed arrows to be kept between (non-treatment and non-outcome) variables, e.g. c("Z1 -> Z2", "Z2 -> Z3"), c("y", "n", "y"), or c(TRUE, FALSE, TRUE).
#' @returns A dagitty object, with directed arrows removed based on edges_to_keep.
#' @examples
#' ## initial dag
#' {
#'   variables <- c("Z1", "Z2", "Z3") # these are treated as confounders
#'   treatments <- "X"
#'   outcomes <- "Y"
#'   type <- "ordered"
#'
#'   dag <- build_graph(type = type,
#'                      variables = variables,
#'                      treatments = treatments,
#'                      outcomes = outcomes)
#' }
#'
#' ## Using a saturated graph as an example, created from an existing dag
#' saturated_graph <- build_graph(type = c('full', 'saturated'), # choose a type
#'                                variables = dag) # existing dagitty object inputted
#'
#' ## guided causal criteria sequence (rule out edges already included in existing dag)
#' edges_to_keep <- assess_edges(saturated_graph, edges_to_keep = dag,
#'                               assess_causal_criteria = TRUE)
#'
#' dag <- keep_edges(saturated_graph, edges_to_keep)
#'
#' @export
keep_edges <- function(dag, edges_to_keep = NA){
  .datatable.aware <- TRUE

  edges_all  <- as.data.table(dagitty::edges(dag))
  edges <- edges_all[, c("v", "e", "w")]


  if( dagitty::is.dagitty(edges_to_keep)){

    edges_to_keep <- data.table::as.data.table(dagitty::edges(edges_to_keep))[, c("v", "e", "w")]

  }else if( is.vector( edges_to_keep ) ){

    edges_to_keep <- data.table::as.data.table( do.call( rbind, strsplit(edges_to_keep, " ") ) )

    colnames(edges_to_keep) <- c("v", "e", "w")

  }

  if( all(is.null(edges_to_keep)) | all(is.na(edges_to_keep)) ){

    stop("No edges to keep were supplied.")

  }


  if( all( !is.na(edges_to_keep) )){

    if( ( is.data.frame(edges_to_keep) | data.table::is.data.table(edges_to_keep) ) ){

      edges <-  edges[ edges_to_keep, on = c("v", "e", "w")]


    }else{

      stop("Edges_to_keep must be a data frame, data table, vector, or dagitty object. Please check inputs and try again.")

    }
  }

  dag <- rebuild_dag(dag, edges)


  if(dagitty::isAcyclic(dag) == FALSE){
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.")
  }

  return(dag)

}


#' Add nodes to a dagitty object
#'
#' Add nodes to a dagitty object, connecting edges based on the 'type' of graph selected, and generate new node coordinates using existing nodes.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_nodes Inputted vector of node names to be added to the graph.
#' @param node_role Determines which nodes will be connected, defaults to "confounder". Options include the following: c("confounder", "treatment", "outcome", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed").
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'first' or 'last', inputted new_nodes are ordered first or last if a confounder, mediator, or treatment node_role is selected.
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings uses dagitty::topologicalOrdering() and selects the first of the inputted node_role (e.g., first confounder) as the temporal point of reference. If type = 'last', the last node is used.
#' @param print_edges Print new edges in console, defaults to TRUE.
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures.
#' @returns A dagitty object.
#' @examples
#' dag <- add_nodes(dag, new_nodes, node_role)
#'
#' @export
add_nodes <- function(dag,
                      new_nodes,
                      node_role = NULL,
                      type = "full",
                      temporal_reference_node = NA,
                      print_edges = TRUE,
                      coords_spec = c(lambda = 0.1, threshold = 0.5)
                      ){
  .datatable.aware <- TRUE

  ## get initial dag edges
  edges <- as.data.frame( dagitty::edges(dag) )[, c("v", "e", "w")]
  names(edges) <- c("ancestor", "edge", "descendant") # change column names

  dag_node_names <- names(dag) # extract dag node names

  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  new_node_names <- new_nodes[ !new_nodes %in% dag_node_names ] # remove duplicate node names


  ## call helper function to draw new node edges etc.

  if( length(node_role) == 0 ){

    new_edges <- connect_all_nodes_to_new(dag, new_nodes)

    existing_node_names <- dag_node_names[ !dag_node_names %in% new_nodes ] # remove excluded_names from existing node names
    treatments <- dagitty::exposures(dag) # treatment
    outcomes <- dagitty::outcomes(dag)  # outcome
    latent_vec <- dagitty::latents(dag)  # latent variables

  }else{ # node role inputted

    output_list <- add_nodes_helper(dag, new_nodes, node_role, type, temporal_reference_node)

    temporal_reference_node <- output_list$temporal_reference_node
    existing_node_names <- output_list$existing_node_names
    new_edges <- output_list$new_edges
    treatments <- output_list$treatments
    outcomes <- output_list$outcomes
    latent_vec <- output_list$latent_vec

  }

  ## pre-process before merging
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs

  ## merge new and existing dag edges
  edges <- merge(edges, new_edges,
                 by = c("ancestor", "edge", "descendant"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges

  edges <- edges[edges$ancestor != edges$descendant, ] # remove new_edges with identical ancestor and descendant node names

  ## create node names and coordinates vector for dagitty object
  exclude_names <- c(treatments, outcomes, latent_vec) # create excluded_names

  node_names <- dag_node_names[ !dag_node_names %in% exclude_names ] # remove excluded_names from existing node names

  node_names <- c(node_names, new_node_names) # combine both existing and new node names

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

  }else{

    dag <- add_coords(dag, coords_spec = coords_spec)

  }

  if( all( !complete.cases(existing_node_names) ) ){

    existing_node_names <- names( coordinates$y[ order(coordinates$y)] ) # existing dag node names in ascending y-coords order

  }

  if( any( !complete.cases(temporal_reference_node) ) ){

    temporal_reference_node <- existing_node_names[1] # first node is selected as the temporal reference node

  }

  new_coordinates <- renew_coords(dag = dag,
                                 new_node_names = new_nodes,
                                 coordinates = coordinates,
                                 coords_spec = coords_spec[ complete.cases(coords_spec) ] )

  dagitty::coordinates(dag) <- new_coordinates

  if( print_edges == TRUE){

    new_edges_list <- print_edges_helper(new_edges)

    cat( paste("c(", paste( unlist(new_edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    message("\nPrinted new edges to assess. - copy and paste in a .R file to use as a vector object.\n")


  }


  return(dag)
}


#' Connect nodes in a dagitty object
#'
#' Add nodes to a dagitty object, connecting edges based on the 'type' of graph selected, and generate new node coordinates using existing nodes.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges latents exposures outcomes
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param nodes Inputted vector of node names to be added to the graph.
#' @param node_role Determines which nodes will be connected, default fully connects new node to all existing nodes. Other options include the following: c("confounder", "treatment", "outcome", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed").
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'first' or 'last', inputted nodes are ordered first or last if a confounder, mediator, or treatment node_role is selected.
#' @param print_edges Print new edges in console, defaults to TRUE.
#' @returns A dagitty object.
#' @examples
#' dag <- saturate_nodes(dag, nodes)
#'
#' @export
saturate_nodes <- function(dag,
                           nodes = NULL,
                           node_role = NULL,
                           type = "full",
                           print_edges = TRUE
                           ){
  .datatable.aware <- TRUE

  ## get initial dag edges
  edges <- as.data.frame( dagitty::edges(dag) )[, c("v", "e", "w")]
  names(edges) <- c("ancestor", "edge", "descendant") # change column names

  dag_node_names <- names(dag) # extract dag node names

  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  if( length(nodes) == 0 ){

    nodes <- dag_node_names

  }

  if( length(node_role) == 0 ){

    new_edges <- saturate_nodes_helper(dag, nodes, dag_node_names, type)


    # treatment
    treatments <- dagitty::exposures(dag)

    # outcome
    outcomes <- dagitty::outcomes(dag)

    # latent variables
    latent_vec <- dagitty::latents(dag)


  }else{

    ## call helper function to draw new node edges etc.
    output_list <- add_nodes_helper(dag, nodes, node_role, type)

    new_edges <- output_list$new_edges
    treatments <- output_list$treatments
    outcomes <- output_list$outcomes
    latent_vec <- output_list$latent_vec

  }

  ## pre-process before merging
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs

  ## merge new and existing dag edges
  edges <- merge(edges, new_edges,
                 by = c("ancestor", "edge", "descendant"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges

  edges <- edges[edges$ancestor != edges$descendant, ] # remove new_edges with identical ancestor and descendant node names

  ## create node names and coordinates vector for dagitty object
  exclude_names <- c(treatments, outcomes, latent_vec) # create excluded_names

  node_names <- dag_node_names[ !dag_node_names %in% exclude_names ] # remove excluded_names from existing node names

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


  if( print_edges == TRUE){

    new_edges_list <- print_edges_helper(new_edges)

    cat( paste("c(", paste( unlist(new_edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    message("\nPrinted new edges to assess. - copy and paste in a .R file to use as a vector object.\n")


  }


  return(dag)
}


#' Replicate nodes in a dagitty object
#'
#' Create single or multiple copy nodes based on existing nodes, or another dagitty object.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param existing_nodes Vector of existing node names, used as a reference for the new graph nodes, e.g., c("Z1", "Z2", "Z3").
#' @param existing_node_type A suffix added to each of the reference node names, e.g. "pre_treatment", or "t0".
#' @param new_node_type A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings use treatment as the temporal point of reference for adding new (post-treatment) nodes, but if other node types are specified it can be useful to specify the existing node name.
#' @param num_repeats Number of additional copies of nodes, such as time points. Each repeat number is included at the end of new node names (new_new_t1, new_node_t2, etc.).
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures.
#' @returns A dagitty object.
#' @examples
#' dag <- copy_nodes(dag, existing_nodes)
#'
#' @export
copy_nodes <- function(dag,
                       existing_nodes,
                       existing_node_type = "pre_treatment",
                       new_node_type = "post_treatment",
                       temporal_reference_node = NA,
                       num_repeats = NA,
                       coords_spec = c(lambda = 0.1, threshold = 0.5)
){
  ###### NOTE: this was written before add_nodes() and saturate_nodes(), as such needs to be updated (following a similar logic as the other node-adding funcs)
  .datatable.aware <- TRUE
  edges <- as.data.frame( dagitty::edges(dag) )[, c("v", "e", "w")] # get dag edges
  node_names <- names(dag) # extract dag node names

  treatments <- dagitty::exposures(dag) # get treatments
  outcomes <- dagitty::outcomes(dag) # get outcomes
  latent_vec <- unlist(dagitty::latents(dag)) # get dag latent variables
  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  if( dagitty::is.dagitty(existing_nodes) ){   # check for an existing dag (inputted as existing_nodes)

    dag <- copy_nodes_helper(dag,
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

  }else{ # copy existing nodes to a new graph

    dag <- copy_nodes_helper(dag,
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


#' Merge two dagitty objects
#'
#' join_graphs() adds new nodes to the first from the second supplied dagitty object, before generating new coordinates for the merged graph.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty edges exposures outcomes latents coordinates
#' @param dag First dagitty object.
#' @param new_dag Second dagitty object, added to the first.
#' @param coords_spec Parameters used for generating coordinates. Adjust node placement with lambda; a higher value increases volatility and results in more extreme DAG structures. Threshold controls the closeness of nodes.
#' @returns A dagitty object
#' @examples
#' new_graph <- join_graphs(dag, new_dag)
#'
#' @export
join_graphs <- function(dag,
                        new_dag,
                        coords_spec = c(lambda = 0.1, threshold = 0.5)
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
  existing_node_coords_y <- coordinates$y[ names(coordinates$y)  %in% existing_node_names] # saves duplicate node names
  existing_node_names <- names( existing_node_coords_y[ order(existing_node_coords_y) ] ) # existing dag node names in ascending y-coords order

  new_node_names <- new_node_names[ !new_node_names %in% node_names ] # remove duplicate node names


  new_dag_treatments <- dagitty::exposures(new_dag)
  treatments <- c(treatments,
                  new_dag_treatments[ !new_dag_treatments %in% node_names] ) # get treatments

  new_dag_outcomes <- dagitty::outcomes(new_dag) # get outcomes
  outcomes <- c(outcomes,
                new_dag_outcomes[ !new_dag_outcomes %in% node_names] ) # get treatments

  new_dag_latent_vec <- unlist( dagitty::latents(new_dag) ) # get dag latent variables
  latent_vec <- c(latent_vec,
                  new_dag_latent_vec[ !new_dag_latent_vec %in% node_names] ) # get treatments


  coordinates_new_dag <- dagitty::coordinates(new_dag) # extract new dag coordinates
  new_node_coords_y <- coordinates_new_dag$y[ names(coordinates_new_dag$y)  %in% new_node_names] # saves duplicate node names
  new_node_names <- names( new_node_coords_y[ order(new_node_coords_y) ] ) # existing dag node names in ascending y-coords order

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

  coordinates <- renew_coords(dag = dag,
                              new_node_names = new_node_names,
                              coordinates = coordinates,
                              coords_spec = coords_spec[ complete.cases(coords_spec) ] )


  dagitty::coordinates(dag) <- coordinates


  return(dag)

}


#' Assess edges using causal criteria
#'
#' This function uses a written procedure outlined in Ferguson et al. (2019); Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs).
#'
#' @param edges vector of edges whose relationships are to be assessed
#' @param num_edges number of edges to be assessed
#' @param check_skip_sequence TRUE or FALSE depending on prior inputs
#' @noRd
causal_criteria_sequence <- function(edges, check_skip_sequence, edges_to_keep){
  #edges <- edges_to_assess # used for debugging assess_edges()
  #check_skip_sequence <- FALSE # used for debugging assess_edges()
  #edges_to_keep <- NA # used for debugging assess_edges()
  if(check_skip_sequence == FALSE) {

    num_edges <- nrow(edges)

    cat("\nThere are ", num_edges, " directed arrows to be assessed.","\n", "\n", sep="")
    print(edges, quote=FALSE)
    cat("\nAssess the posited causal relationships using causal criteria? (ESC-DAGs protocol)", "\n")

    removed_arrows <- c()
    arrow_count <- 1
    num_arrow_to_remove <- 0

    check_ans <- FALSE

    while(check_ans == FALSE){

      choice <- readline("(y/n/?info): ")

      if(choice == "y"){

        #cat("\n", "\n===================================")
        #cat("\nESC-DAGs (Ferguson et al., 2020)")
        #cat("\nDAGs from background knowledge")
        #cat("\nCode written by Aidan J Moller (2025)")
        #cat("\n===================================", "\n")

        check_ans <- TRUE

      }else if(choice == "n"){

        message("\nSkipped sequence.")

        #edges <- noquote(paste("c('", paste(edges, collapse="', '"), "')", sep = ""))

        # edges_to_keep <- noquote(paste("c('", paste(edges_to_keep, collapse="', '"), "')", sep = ""))

        message("\nOutputting causal criteria assessed edges.")

        edges_list <- list(
          edges_to_keep = edges_to_keep,
          edges_to_assess = edges
        )

        return(edges_list)

      }else if(choice == "?info"){

        cat("\nEach directed edge in the IG is assessed for three causal criteria: temporality; face-validity; and recourse to theory. They are primarily informed by the classic Bradford Hill viewpoints,24 and are compatible with the ‘inference to the best explanation’ approach advocated by Krieger and Davey Smith.1 If a relationship is determined to possess each criterion, a counterfactual thought experiment derived from the POF is used to further explicate the reviewers’ assumptions.25 The translation process thus combines ‘classic’ and ‘modern’ causal thinking and understands DAGs as ‘conceptual tools’1 for exploring causation, rather than substitutes for careful causal thinking.", "\n")
        cat("\nThe ESC-DAGs causal criteria operate sequentially, with each criterion designed to elaborate over the previous. If any criterion on the edge is not present, the edge can be deleted. The exception is the recourse to theory criterion—absence of theory in the study or according to the reviewer does not equate to absence of effect.26 The counterfactual thought experiment is performed after assessing all criteria. All retained directed edges are entered into the directed edge index. However, each edge should be tested in both directions (i.e. with the head and tail of the arrow swapped). If the posited and reverse edges are both retained, then the relationship should be noted as bi-directional in the directed edge index. Reviewers can also note low confidence in particular directed edges.", "\n")

      }else{

        cat("\n", "\nPlease type a valid answer.", "\n")

      }
    }

    for(arrow_count in 1:num_edges){

      criterion_num <- 1
      check_arrow_removed <- FALSE

      check_ans <- FALSE

      arrow <- noquote( paste(edges[arrow_count], collapse=" ") )
      cat("\n", "\nFor the directed arrow '", arrow, "' consider each of the following:", sep="")

      cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
      cat("\n[1/4] Temporality: does the variable to the left of the arrow precede the variable on the right?", "\n")
      cat("\nFor help, enter ?info")

      while(check_ans == FALSE){

        choice <- readline("(y/n/?info): ")

        if(choice == "y"){

          criterion_num <- criterion_num + 1

          check_ans <- TRUE

        }else if(choice == "n"){

          num_arrow_to_remove <- num_arrow_to_remove + 1
          removed_arrows[num_arrow_to_remove] <- arrow_count

          check_arrow_removed <- TRUE

          check_ans <- TRUE

          message("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

        }else if(choice == "?info"){

          cat("\n", "\nCausal criterion 1—temporality:")
          cat("\nOf the Bradford Hill criteria, temporality is the only one not requiring extensive qualification or not yet disproven. (Thomas et al., 2013; DOI: https://doi.org/10.1146/annurev-publhealth-031811-124606) It states that effect cannot precede cause. For example, in Figure 1(A) (Ferguson et al., 2020; DOI: https://doi.org/10.1093/ije/dyz150), adolescent substance use cannot precede historical parental alcohol use, so the relationship would not be temporal. Unless the directed edge is not temporal, we proceed to causal criterion 2.", "\n")
          cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")

        }else{

          cat("\n", "\nPlease type a valid answer.", "\n")

        }
      }

      check_ans <- FALSE

      if(check_arrow_removed == FALSE){

        cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
        cat("\n[2/4] Face-validity: is the posited relationship plausible?", "\n")
        cat("\nFor help, enter ?info")

        while(check_ans == FALSE){

          choice <- readline("(y/n/?info): ")

          if(choice == "y"){

            criterion_num <- criterion_num + 1

            check_ans <- TRUE

          }else if(choice == "n"){

            num_arrow_to_remove <- num_arrow_to_remove + 1
            removed_arrows[num_arrow_to_remove] <- arrow_count

            check_arrow_removed <- TRUE

            check_ans <- TRUE

            cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

          }else if(choice == "?info"){

            cat("\n", "\nCausal criterion 2—face-validity:")
            cat("\nFace-validity is related to the Bradford Hill criterion of (biologic) plausibility. Nested within the wider causal criteria scheme, the face-validity criterion is a rapid means of using reviewer background knowledge to identify implausible relationships, given the temporality established in criterion 1. For example, in Figure 1(A) it is plausible that directed edges originate from sex, but implausible that historical parental alcohol use could influence adolescent sex assignment despite temporal ordering.", "\n")
            cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")


          }else{

            cat("\n", "\nPlease type a valid answer.", "\n")

          }
        }
      }

      check_ans <- FALSE

      if(check_arrow_removed == FALSE){

        cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
        cat("\n[3/4] Recourse to theory—is the posited relationship supported by theory?", "\n")
        cat("\nFor help, enter ?info")

        while(check_ans == FALSE){

          choice <- readline("(y/n/?info): ")

          if(choice == "y"){

            criterion_num <- criterion_num + 1
            check_ans <- TRUE

          }else if(choice == "n"){

            #num_arrow_to_remove <- num_arrow_to_remove + 1
            #removed_arrows[num_arrow_to_remove] <- arrow_count

            #check_arrow_removed <- TRUE

            check_ans <- TRUE

            #cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

            cat("\n", "\nCausal relationship '", arrow, "' assessed; answer recorded.", "\n", sep="")

          }else if(choice == "?info"){

            cat("\n", "\nCausal criterion 3—recourse to theory:")
            cat("\nThe recourse to theory criterion considers background and expert knowledge more overtly. It subsumes the temporality and face-validity criteria and continues to cement a platform for the counterfactual thought experiment. Where the face-validity criterion is concerned with the researcher’s own knowledge, the step assesses whether there is formal theoretical support for the relationship. The decision log for this criterion requires the reviewer to state briefly what theory applies (if any) with space for a reference. As noted above, lack of theory does not equate to lack of effect. As such the purpose of this criterion is not so much falsification as preparation for the next step.", "\n")
            cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")

          }else{

            cat("\n", "\nPlease type a valid answer.", "\n")

          }
        }
      }

      check_ans <- FALSE

      if(check_arrow_removed == FALSE){

        cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
        cat("\n[4/4] Counterfactual thought experiment: is the posited relationship supported by a systematic thought experiment informed by the potential outcomes framework?")
        cat("\n", "\nFor help, enter ?info")

        while(check_ans == FALSE){

          choice <- readline("(y/n/?info): ")

          if(choice == "y"){

            criterion_num <- criterion_num + 1

            check_ans <- TRUE

          }else if(choice == "n"){

            num_arrow_to_remove <- num_arrow_to_remove + 1
            removed_arrows[num_arrow_to_remove] <- arrow_count

            check_arrow_removed <- TRUE

            check_ans <- TRUE

            cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

          }else if(choice == "?info"){

            cat("\n", "\nCounterfactual thought experiment")
            cat("\nFundamentally, potential outcomes compare the outcome that would have occurred if all of the sample had been exposed, with the outcome that would have occurred if all of the sample had not been exposed.3,4,25 The counterfactual thought experiment employs this heuristic in a formulaic and transparent way, comparing two or more ‘counterfactual exposures’ and considering whether their potential outcomes would be different, given the causal criteria. The original study’s measurement of variables should be emulated. See Ferguson et al. (2020) for more details.", "\n")
            cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")

          }else{

            cat("\n", "\nPlease type a valid answer.", "\n")

          }
        }
      }

      arrow_count <- arrow_count + 1
    }
  }

  if(num_arrow_to_remove > 0){

    edges <- edges[-removed_arrows, ]

    edges_to_keep <- rbind(edges_to_keep, edges)

    return(edges_to_keep)

  }else{

    message("\nNo arrows were removed.", "\n")

    edges_list <- list(
      edges_to_keep = edges_to_keep,
      edges_to_assess = edges
    )

    return(edges_list)
  }
}


