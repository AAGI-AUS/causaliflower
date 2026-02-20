#' Build dagitty objects
#'
#' build_graph() produces a dagitty graph object from inputted parameters (e.g. treatments, outcome, confounders).
#'
#' @importFrom dagitty is.dagitty children coordinates
#' @param variables Dagitty object, or vector of variable names, e.g. "Z" or c("Z1", "Z2", "Z3"). If variable names are inputted the order determines the assigned coordinates. A list can also be supplied. Variables inputted are treated as confounders. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param treatments Treatment variable name, e.g. "X". Must be specified.
#' @param outcomes Outcome variable name, e.g. "Y". Must be specified.
#' @param mediators Character or vector of mediator variable names, e.g. "M" or c("M1", "M2", "M3").
#' @param latent_variables Character or vector of additional or already supplied latent (unobserved) variable names, e.g. "U" or c("U1", "U2", "M1").
#' @param instrumental_variables Vector of instrumental variable names, e.g. "IV"
#' @param mediator_outcome_confounders Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param competing_exposures Vector of competing exposure names. An arrow is drawn connecting competing exposures to the outcome, with other arrows also connected depending on type of graph specified.
#' @param colliders Vector of collider variables, with both treatment and outcome parents.
#' @param type Type of connected graph (e.g. "full", "saturated", "ordered", "none"). Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
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
build_graph <- function(variables = NA,
                       treatments = NA,
                       outcomes = NA,
                       mediators = NA,
                       latent_variables = NA,
                       instrumental_variables = NA,
                       mediator_outcome_confounders = NA,
                       competing_exposures = NA,
                       colliders = NA,
                       type = "ordered"
                       ){

  # check for an existing dag (inputted as confounders)
  if( dagitty::is.dagitty(variables) ){

    existing_dag <- variables
    #existing_dag <- dag
    node_roles <- get_roles(existing_dag)

    confounders <- variables <- node_roles$confounder

    treatments <- node_roles$treatment

    outcomes <- node_roles$outcome

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

  }else if( all( complete.cases( unlist(treatments) ) ) & all( complete.cases( unlist(outcomes) ) ) ){ # no existing dag, use other inputs

    confounders <- variables

    existing_dag <- NA

    mediator_vec <- as.vector( unlist( lapply( mediators, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    latent_vec <- get_latent_vec(latent_variables)

    m_o_confounder_vec <- as.vector( unlist( lapply( mediator_outcome_confounders, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    competing_exposure_vec <- as.vector( unlist( lapply( competing_exposures, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    collider_vec <- as.vector( unlist( lapply( colliders, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    observed <- NA

    if( all( complete.cases( unlist(confounders) ) ) & any( unlist(confounders) %in% m_o_confounder_vec ) ){ # if any mediator-outcome confounders are also inputted as confounders, execution is stopped

      stop("Node names inputted in the 'variables' parameter detected in 'mediator_outcome_confounder'. These inputs should be mutually exclusive. Please adjust inputs and try again.")
    }

  }else{

    stop("The 'treatments' and 'outcomes' inputs should be used if a DAG is not provided in the 'variables' input. Please adjust inputs and try again.")
  }

  ## get variable names ##
  observed_node_names <- unique( as.vector( c(unlist(confounders), m_o_confounder_vec, mediator_vec, competing_exposure_vec, collider_vec, instrumental_variables) ) )
  observed_node_names <- Filter(Negate(anyNA), observed_node_names)

  edges <- draw_edges(type = type,
                      confounders = confounders,
                      treatments = treatments,
                      outcomes = outcomes,
                      mediator_vec = mediator_vec,
                      latent_vec = latent_vec,
                      latent_variables = latent_variables,
                      instrumental_variables = instrumental_variables,
                      m_o_confounder_vec = m_o_confounder_vec,
                      competing_exposure_vec = competing_exposure_vec,
                      collider_vec = collider_vec,
                      observed = observed,
                      observed_node_names = observed_node_names,
                      existing_dag = NA)

  ## remove treatments, outcomes and latents from node names ##
  exclude_names <- c(treatments, unlist(outcomes), latent_vec)
  exclude_names <- exclude_names[ complete.cases(exclude_names) ]
  node_names <- observed_node_names[!observed_node_names %in% exclude_names]

  dag <- construct_graph(edges,
                         node_names,
                         treatments,
                         outcomes,
                         latent_vec)

  return(dag)

}

#' Assess dagitty object edges
#'
#' assess_edges() provides ways to assess connected edges based on causal criteria and/or user inputs.
#'
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param edges_to_assess Defaults to "all" edges, "bidirectional" includes only the edges where nodes are connected in both directions.
#' @param edges_to_keep A vector of directed arrows to be kept between (non-treatment and non-outcome) variables, e.g. c("Z1 -> Z2", "Z2 -> Z3"), c("y", "n", "y"), or c(TRUE, FALSE, TRUE).
#' @param assess_causal_criteria Defaults to FALSE. If TRUE, the user is guided through a sequence that assesses each pair of connected nodes using causal criteria. Based on the Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs) from Ferguson et al. (2020).
#' @param causal_criteria = Default "ESCDAGs" (Ferguson et al., 2020) causal criteria considers temporality, face-validity, and recourse to theory. Other causal criteria can be supplied as a data.frame, see summary(ESCDAGs) for column names.
#' @param causal_criteria_answers = NULL,
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
assess_edges <- function(dag,
                         edges_to_assess = "all",
                         edges_to_keep = NA,
                         assess_causal_criteria = FALSE,
                         causal_criteria = ESCDAGs,
                         causal_criteria_answers = NULL
                         ){
  .datatable.aware <- TRUE
  # debugging: edges_to_keep <- initial_dag # dag <- saturated_graph

  # get dag edges
  edges <- pdag_to_dag_edges(dag) # also checks for bi-directional edges "--" (pdag) or "<->" and changes to "->", adding rows for each edge

  if( dagitty::is.dagitty(edges_to_keep) ){

    edges_to_keep <- pdag_to_dag_edges(edges_to_keep)

    edges <-  edges[!edges_to_keep, on = c("v", "e", "w")]

  }else if( all( complete.cases( unlist(edges_to_keep) ) ) ){

    if(is.vector( edges_to_keep ) ){

      edges_to_keep <- data.table::as.data.table( do.call( rbind, strsplit(edges_to_keep, " ") ) )

      colnames(edges_to_keep) <- c("v", "e", "w")

    }

    if( is.data.frame(edges_to_keep) | data.table::is.data.table(edges_to_keep) ){

      edges <-  edges[!edges_to_keep, on = c("v", "e", "w")]

    }else{

      stop("'edges_to_keep' must be a data frame, data table, vector, or dagitty object. Please check inputs and try again.")

    }

  }


  if( assess_causal_criteria == TRUE ){

    if( edges_to_assess == "bidirectional" ){ # if use inputs edges_to_assess = "bidirectional"

      # identify all bidirectional edges
      edges_to_assess <- suppressWarnings( lapply( 1:nrow(edges), function(x){

        lapply( 1:nrow(edges), function(y){

          if( ( identical( edges$v[x], edges$w[y] ) & identical( edges$w[x], edges$v[y] ) ) ){

            edges_to_assess[x] <- edges[x, ]

          }

        } )

      } ) )


      edges_to_assess <- do.call(rbind, unlist(edges_to_assess, recursive = FALSE) ) # expand list elements as dataframe rows

      # remove edges_to_assess from dag edges
      edges <- data.table::setDT(edges)[!edges_to_assess, on = c("v", "e", "w")]

      if( all( complete.cases( edges_to_keep ) ) ){ # check if edges_to_keep specified

        # combine unidirectional edges and edges_to_keep
        edges_to_keep <- rbind(edges, edges_to_keep)

      }else{ # otherwise unidirectional edges become edges_to_keep

        edges_to_keep <- edges
      }

    }else{ # edges_to_assess = all

      edges_to_assess <- edges
    }

    edges <- causal_criteria_sequence(edges = edges_to_assess,
                                      check_skip_sequence = FALSE,
                                      edges_to_keep = edges_to_keep,
                                      causal_criteria = causal_criteria,
                                      causal_criteria_answers = causal_criteria_answers)

    return(edges)

  }

  if( nrow(edges) != 0 ){


    if( all( complete.cases( unlist(edges_to_keep) ) ) ){

      ## collapse edges_to_keep to a vector
      num_edges <- nrow(edges_to_keep)

      edges_to_keep <- suppressWarnings( sapply(1:num_edges, function(x){

        edges_to_keep <- paste( edges_to_keep[x,], collapse=" ")

      }) )

    }

    # group edges by unique node names
    edges_list <- print_edges_helper(edges)

    cat( paste("c(", paste( unlist(edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

  }else{

    stop("There are no edges to assess. Please check the supplied dagitty object and try again.")


  }



  if( !all( complete.cases(edges_to_keep) ) & assess_causal_criteria == FALSE ){

    edges <- suppressWarnings( sapply(1:nrow(edges), function(x){

      edges <- paste( edges[x,], collapse=" ")

    }) )

    message("\nOutputted and printed edges to assess - copy and paste in a .R file to use as a vector object.")

    return(edges)

  }

  edges_list <- list(
    edges_to_assess = edges,
    edges_to_keep = edges_to_keep

  )

  message("\nOutputted edges list contains 'edges_to_assess' and 'edges_to_keep'.
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
keep_edges <- function(dag,
                       edges_to_keep = NA
                       ){
  .datatable.aware <- TRUE

  edges <- pdag_to_dag_edges(dag)

  if( dagitty::is.dagitty(edges_to_keep)){

    edges_to_keep <- pdag_to_dag_edges(edges_to_keep)

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
#' @param coords_spec Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures. Threshold controls the closeness of nodes.
#' @returns A dagitty object.
#' @examples
#' dag <- add_nodes(dag, new_nodes, node_role)
#'
#' @export
add_nodes <- function(dag,
                      new_nodes,
                      ancestors = NULL,
                      descendants = NULL,
                      node_role = NULL,
                      type = NULL,
                      print_edges = TRUE,
                      coords_spec = c(lambda = 0.1, threshold = 0.5)
                      ){
  .datatable.aware <- TRUE

  ## get initial dag edges
  edges <- pdag_to_dag_edges(dag)

  names(edges) <- c("v", "e", "w") # change column names

  dag_node_names <- names(dag) # extract dag node names

  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  new_node_names <- new_nodes[ !new_nodes %in% dag_node_names ] # remove duplicate node names

  new_edges <- connect_new_nodes(dag = dag,
                                 new_nodes = new_node_names,
                                 ancestors = ancestors,
                                 descendants = descendants,
                                 node_role = node_role,
                                 type = type)

  ## pre-process before merging
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs

  ## merge new and existing dag edges
  edges <- merge(edges, new_edges,
                 by = c("v", "e", "w"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges

  edges <- edges[edges$v != edges$w, ] # remove new_edges with identical ancestor and descendant node names

  dag <- rebuild_dag(dag, edges)

  if( print_edges == TRUE){

    new_edges_list <- print_edges_helper(new_edges)

    cat( paste("c(", paste( unlist(new_edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    message("\nPrinted new edges - copy and paste to use as a vector object.\n")

  }

  num_nodes <- length(new_node_names)
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
    dag <- add_coords_helper(dag, coords_spec = coords_spec[1][ complete.cases(coords_spec) ] )
    return(dag)

  }, error = function(e){
    message(paste("Error:", e, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords_helper(dag, coords_spec = coords_spec[1][ complete.cases(coords_spec) ] )
    return(dag)

  }, finally = return(dag)

  )


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
  edges <- pdag_to_dag_edges(dag)
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

    colnames(new_edges) <- c("v", "e", "w")
  }

  ## pre-process before merging
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs

  ## merge new and existing dag edges
  edges <- merge(edges, new_edges,
                 by = c("v", "e", "w"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges

  edges <- edges[edges$v != edges$w, ] # remove new_edges with identical ancestor and descendant node names

  dag <- rebuild_dag(dag, edges)

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
#' @param coords_spec Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures. Threshold controls the closeness of nodes.
#' @returns A dagitty object.
#' @examples
#' dag <- copy_nodes(dag, existing_nodes)
#'
#' @noRd
copy_nodes <- function(dag,
                       existing_nodes,
                       existing_node_type = "pre_treatment",
                       new_node_type = "post_treatment",
                       temporal_reference_node = NA,
                       num_repeats = NA,
                       coords_spec = c(lambda = 0.1, threshold = 0.5)
                       ){
  ###### NOTE: this was written before add_nodes() and saturate_nodes()
  ###### needs to be updated before adding @export
  .datatable.aware <- TRUE
  edges <- pdag_to_dag_edges(dag) # get dag edges
  node_names <- names(dag) # extract dag node names

  treatments <- dagitty::exposures(dag) # get treatments
  outcomes <- dagitty::outcomes(dag) # get outcomes
  latent_vec <- unlist(dagitty::latents(dag)) # get dag latent variables
  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  if( dagitty::is.dagitty(existing_nodes) ){   # check for an existing dag (inputted as existing_nodes)

    dag <- copy_nodes_helper2(dag,
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
#' @param coords_spec Adjust node placement with lambda; a higher value increases volatility and results in more extreme DAG structures. Threshold controls the closeness of nodes.
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

  edges <- pdag_to_dag_edges(dag) # get primary dag edges
  node_names <- names(dag) # extract dag node names

  treatments <- dagitty::exposures(dag) # get treatments
  outcomes <- dagitty::outcomes(dag) # get outcomes
  latent_vec <- unlist(dagitty::latents(dag)) # get dag latent variables
  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  new_edges <- pdag_to_dag_edges(new_dag)  # get new dag edges
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

  dag <- rebuild_dag(dag, edges)


  num_nodes <- length(length(new_node_names))
  time_limit <- num_nodes + num_nodes*coords_spec[1]

  setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)

  on.exit( {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  } )

  tryCatch({
    coordinates <- renew_coords(dag = dag,
                                new_node_names = new_node_names,
                                coordinates = coordinates,
                                coords_spec = coords_spec[ complete.cases(coords_spec) ] )

    dagitty::coordinates(dag) <- coordinates

    }, warning = function(w){
      message(paste("Warning:", w, "\n Using alternative function to generate dag coordinates."))
      dag <- add_coords_helper(dag, coords_spec = coords_spec[ complete.cases(coords_spec) ] )
      return(dag)

    }, error = function(e){
      message(paste("Error:", e, "\n Using alternative function to generate dag coordinates."))
      dag <- add_coords_helper(dag, coords_spec = coords_spec[ complete.cases(coords_spec) ] )
      return(dag)

    }, finally = return(dag)
  )

}


#' Assess graph edges using causal criteria
#'
#' @param edges vector of edges whose relationships are to be assessed
#' @param num_edges number of edges to be assessed
#' @param check_skip_sequence TRUE or FALSE depending on prior inputs
#' @param causal_criteria Set of causal criteria to be used. Can be user-specified, defaults to 'ESCDAGs'.
#' @noRd
causal_criteria_sequence <- function(edges,
                                     check_skip_sequence,
                                     edges_to_keep,
                                     causal_criteria,
                                     causal_criteria_answers
                                     ){
  # debugging: edges <- edges_to_assess # check_skip_sequence <- FALSE # edges_to_keep <- NA
  if(check_skip_sequence == FALSE) {

    num_edges <- nrow(edges)
    cat( "\nThere are ", num_edges, " directed arrows to be assessed.","\n", "\n", sep="" )
    print(edges, quote=FALSE)
    cat( "\nAssess the posited causal relationships using causal criteria?", "\n" )

    check_ans <- FALSE

    while(check_ans == FALSE){

      choice <- readline("(y/n/?info): ")

      if(choice == "y"){

        check_ans <- TRUE

      }else if(choice == "n"){

        message( "Skipped sequence." )

        edges_list <- list(edges = edges_to_keep,
                           edges_to_assess = edges)
        return(edges_list)

      }else if(choice == "?info"){

        cat( "This is a guided sequence for assessing graph edges. By default it uses the causal criteria in ESC-DAGs protocol (Ferguson et al., 2020).
             \nA data frame containing user-specified criteria can be supplied in assess_edges().")
      }else{

        message( "Please type a valid answer." )
      }

    }

    removed_arrows <- c()
    arrow_count <- 1
    num_arrow_to_remove <- 0

    num_criteria <- nrow(causal_criteria)
    criteria_answers <- as.data.frame( matrix(nrow = 0, ncol = num_criteria) )
    colnames(criteria_answers) <- causal_criteria[, "name" ]

    causal_criteria$required <- as.factor(causal_criteria$required)

    for(arrow_count in 1:num_edges){

        criterion_num <- 1
        check_ans <- FALSE

        arrow <- noquote( paste( edges[ arrow_count ], collapse=" " ) )
        cat( "For the directed arrow '", arrow, "' consider each of the following:", "\n", "\n", sep="")

        while( check_ans == FALSE ){

          if( !criterion_num > num_criteria ){

            cat( paste0("'", arrow, "' (", arrow_count, "/", num_edges, ")"),
                 "\n", paste0("\n[", criterion_num, "/", num_criteria, "]"),
                 paste0( causal_criteria[ criterion_num, "name" ], ":"),
                 causal_criteria[ criterion_num, "question" ], "\n",
                 "\nFor help, enter \'?info\'")

            choice <- readline("(y/n/?info): ")

            if( choice == "y" ){

              criteria_answers[ arrow_count, criterion_num] <- "y" # "yes (causal)"

              criterion_num <- criterion_num + 1

            }else if( choice == "n" ){

              if( causal_criteria[ criterion_num, "required" ] == "yes"){

                num_arrow_to_remove <- num_arrow_to_remove + 1
                removed_arrows[ num_arrow_to_remove ] <- arrow_count

                criteria_answers[ arrow_count, criterion_num ] <- "n" # "no - remove edge"

                criterion_num <- num_criteria + 1

                message( "Causal relationship '", arrow, "' assessed; edge removed.", sep="" )

                check_ans <- TRUE

              }else{

                criteria_answers[ arrow_count, criterion_num ] <- "n"

                criterion_num <- criterion_num + 1

                message( "Answer recorded. Causal relationship '", arrow, "' assessed.", sep="" )
              }

            }else if( choice == "?info" ){

              cat( "\n", "\nCausal criterion", criterion_num, "—", causal_criteria[ criterion_num, "name" ],
                   "\n", "\n", paste0("Definition", ":"),
                   "\n", causal_criteria[ criterion_num, "description" ],
                   "\n", "\n", causal_criteria[ criterion_num, "source" ],
                   "\n", "\n" )
            }else{

              cat( "\n", "\nPlease type a valid answer.", "\n" )
            }

          }else{

            check_ans <- TRUE
          }

        }

        arrow_count <- arrow_count + 1
      }



    if(num_arrow_to_remove > 0){

      criteria_answers <- cbind(edges, criteria_answers)

      edges <- edges[ -removed_arrows, ]

      if( all( complete.cases( edges_to_keep ) ) ){

        edges <- rbind(edges_to_keep, edges)

      }

      edges <- list(edges = edges,
                    causal_criteria_answers = criteria_answers)

      message("\nOutputted list containing 'edges' and 'causal_criteria_answers'.")

      return(edges)

    }else{

      message( "No arrows were removed." )

      edges_list <- list(edges = rbind(edges_to_keep, edges),
                         causal_criteria_answers = criteria_answers)

      message("\nOutputted list containing 'edges' and 'causal_criteria_answers'.")

      return(edges_list)

    }

  }

  return(edges_list)
}


