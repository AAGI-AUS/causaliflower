#' Build dagitty objects
#'
#' build_graph() produces a saturated graph by default from the supplied treatments, outcome, confounder, and  mediators (optional) based on the type of graph specified. Optional latent confounder or medfiator variables can be specified.
#'
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom dagitty is.dagitty children coordinates
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
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
#' dag <- build_graph(type = type,
#'                   variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes)
#'
#' ggdagitty(dag) # Plot the graph with ggdagitty() (assigns coordinates based on the variable roles, e.g., confounders, treatments, outcomes).
#'
#'
#' ## Option 2: 'saturated' graph connecting each of the confounders (inputted as variables)
#'
#' type <- "saturated"
#'
#' dag <- build_graph(type = type,
#'                   variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes)
#'
#' ggdagitty(dag)
#'
#'
#' ## Option 3: "full" creates a fully connected graph (bidirectional arrows between confounders, and from confounders to all other nodes except treatment and outcome.)
#'
#' type <- "full"
#'
#' dag <- build_graph(type = type,
#'                   variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes)
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
#' @examples
#' merged_dag <- merge_graphs(dag, new_dag)
#'
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
                                                 coordinates = coordinates,
                                                 coords_spec = coords_spec[ complete.cases(coords_spec) ] )


  dagitty::coordinates(dag) <- new_coordinates


  return(dag)

}


#' Duplicate existing dagitty object nodes
#'
#' Create single or multiple copy nodes based on existing nodes, or another dagitty object.
#'
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
#' Add nodes to a dagitty object, connecting edges based on the 'type' of graph selected, and generate new node coordinates using existing nodes.
#'
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_nodes A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param node_role Role assigned to new nodes, from any of the following: c("confounder", "treatment", "outcome", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed").
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'first' or 'last', inputted new_nodes are ordered first or last if a confounder, mediator, or treatment node_role is selected.
#' @param temporal_reference_node Supply an alternative reference, or simply leave blank. Default settings uses dagitty::topologicalOrdering() and selects the first of the inputted node_role (e.g., first confounder) as the temporal point of reference. If type = 'last', the last node is used.
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures.
#' @returns A dagitty object.
#' @examples
#' dag <- add_nodes(dag,
#'                  new_nodes = new_node_vec,,
#'                  node_role = "confounder",
#'                  type = "last")
#'
#' @export
add_nodes <- function(dag,
                      new_nodes,
                      node_role = c("confounder", "treatment", "outcome", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed"),
                      type = c("full", "saturated", "first", "last"),
                      temporal_reference_node = NA,
                      coords_spec = 1
){
  .datatable.aware <- TRUE

  ## get initial dag edges
  edges <- as.data.frame( dagitty::edges(dag) )[, c("v", "e", "w")]
  names(edges) <- c("ancestor", "edge", "descendant") # change column names

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

  dag_node_names <- names(dag) # extract dag node names

  treatments <- dagitty::exposures(dag) # get treatments
  outcomes <- dagitty::outcomes(dag) # get outcomes
  latent_vec <- unlist(dagitty::latents(dag)) # get dag latent variables

  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  nodes_ordered <- sort( unlist( dagitty::topologicalOrdering(dag) ) ) # ggdag estimated temporal order of new nodes

  new_node_names <- new_nodes[ !new_nodes %in% dag_node_names ] # remove duplicate node names

  if( length( node_role) > 1 | length( node_role) == 0){

    stop("add_nodes() currently only supports single node_role character inputs.")

  }else if( node_role %in% "confounder" ){
    confounder_occurrance <- as.numeric(order(match(new_nodes, new_nodes)))

    ## confounder edges ##
    new_edges <- draw_confounder_edges(type,
                                           outcomes,
                                           treatments,
                                           confounder_vec = new_nodes, # new nodes as confounder
                                           m_o_confounder_vec,
                                           mediator_vec,
                                           latent_vec,
                                           confounder_occurrance)

    # connect all confounders (fully connected or saturated graph type)
    if( type == "full" | type == "saturated" ){

      confounder_list <- list()

      confounder_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

          list( c( ancestor = new_nodes[x], edge = "->", descendant = confounder_vec[y]) )

        })

      }) )

      confounder_list <- Filter(Negate(anyNA), unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE))
      new_edges <- dplyr::bind_rows(new_edges, confounder_list)


      confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

        confounder_list[x] <- lapply(1:length(new_nodes), function(y){

          list( c( ancestor = confounder_vec[x], edge = "->", descendant = new_nodes[y]) )

        })

      }) )

      confounder_list <- Filter(Negate(anyNA), unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE))
      new_edges <- dplyr::bind_rows(new_edges, confounder_list)


    }else if( type == "first"){

      first_var <- names(nodes_ordered)[ names(nodes_ordered) %in% confounder_vec ][1]

      temporal_reference_node <- first_var

      first_list <- list()

      first_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        first_list[x] <- list( c( ancestor = new_nodes[x], edge = "->", descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), unlist(unlist(first_list, recursive = FALSE), recursive = FALSE))

      new_edges <- dplyr::bind_rows(new_edges, first_list)

    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% confounder_vec ][length(confounder_vec) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = "->", descendant = new_nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), unlist(unlist(last_list, recursive = FALSE), recursive = FALSE))

      new_edges <- dplyr::bind_rows(new_edges, last_list)

    }

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% confounder_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "treatment" ){
    ## treatment edges ##
    new_edges <- draw_treatment_edges(type,
                                         outcomes,
                                         treatments = new_nodes, # new nodes as treatment
                                         confounder_vec,
                                         mediator_vec,
                                         collider_vec)

    # connect all treatments (fully connected or saturated graph type)
    if( type == "full" ){

      treatment_list <- list()

      treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

        treatment_list[x] <- lapply(1:length(new_nodes), function(y){

          list( c( ancestor = treatments[x], edge = "->", descendant = new_nodes[y]) )

        })

      }) )

      treatment_list <- Filter(Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE))
      new_edges <- dplyr::bind_rows(new_edges, treatment_list)


      treatment_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        treatment_list[x] <- lapply(1:length(treatments), function(y){

          list( c( ancestor = new_nodes[x], edge = "->", descendant = treatments[y]) )

        })

      }) )

      treatment_list <- Filter(Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE))
      new_edges <- dplyr::bind_rows(new_edges, treatment_list)

    }else if( type == "first"){

      first_var <- names(nodes_ordered)[ names(nodes_ordered) %in% treatments ][1]

      temporal_reference_node <- first_var

      first_list <- list()

      first_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        first_list[x] <- list( c( ancestor = new_nodes[x], edge = "->", descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), unlist(unlist(first_list, recursive = FALSE), recursive = FALSE))

      new_edges <- dplyr::bind_rows(new_edges, first_list)

    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% treatments ][length(treatments) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = "->", descendant = new_nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), unlist(unlist(last_list, recursive = FALSE), recursive = FALSE))

      new_edges <- dplyr::bind_rows(new_edges, last_list)

    }


    treatments <- c(treatments, new_nodes)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% treatments ] ) # existing node names in temporal order

  }else if( node_role %in% "outcome" ){
    ## outcome edges ##
    new_edges <- draw_outcome_edges(type, outcomes = new_nodes, # new nodes as outcome
                                    collider_vec)

    outcomes  <- c(outcome, new_nodes)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% outcomes ] ) # existing node names in temporal order

  }else if( node_role %in% "mediator" ){
    ## mediator edges ##
    new_edges <- draw_mediator_edges(type,
                                       outcomes,
                                       mediator_vec = new_nodes, # new nodes as mediator
                                       latent_vec)

    if(type == "full"){

      mediator_list <- list()

      mediator_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        mediator_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = new_nodes[x], edge = "->", descendant = mediator_vec[y]) )

        })

      }) )

      mediator_list <- Filter( Negate(anyNA), unlist(unlist(mediator_list, recursive = FALSE), recursive = FALSE) )
      new_edges <- dplyr::bind_rows(new_edges, mediator_list)

      mediator_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        mediator_list[x] <- lapply(1:length(new_nodes), function(y){

          list( c( ancestor = mediator_vec[x], edge = "->", descendant = new_nodes[y]) )

        })

      }) )

      mediator_list <- Filter( Negate(anyNA), unlist(unlist(mediator_list, recursive = FALSE), recursive = FALSE) )
      new_edges <- dplyr::bind_rows(new_edges, mediator_list)

    }else if( type == "first"){

      first_var <- names(nodes_ordered)[ names(nodes_ordered) %in% mediator_vec ][1]

      temporal_reference_node <- first_var

      first_list <- list()

      first_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        first_list[x] <- list( c( ancestor = new_nodes[x], edge = "->", descendant = first_var) )


      }) )

      first_list <- Filter(Negate(anyNA), unlist(unlist(first_list, recursive = FALSE), recursive = FALSE))

      new_edges <- dplyr::bind_rows(new_edges, first_list)

    }else if( type == "last" | type == "ordered"){

      last_var <- names(nodes_ordered)[ names(nodes_ordered) %in% mediator_vec ][length(mediator_vec) ]

      temporal_reference_node <- last_var

      last_list <- list()

      last_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

        last_list[x] <- list( c( ancestor = last_var, edge = "->", descendant = new_nodes[x]) )

      }) )

      last_list <- Filter(Negate(anyNA), unlist(unlist(last_list, recursive = FALSE), recursive = FALSE))

      new_edges <- dplyr::bind_rows(new_edges, last_list)

    }

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% mediator_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "mediator_outcome_confounder" ){
    ## mediator_outcome_confounder edges ##
    new_edges <- draw_moc_edges(type,
                             outcomes,
                             confounder_vec,
                             m_o_confounder_vec = new_nodes, # new nodes as mediator_outcome_confounder
                             mediator_vec,
                             latent_vec)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% m_o_confounder_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "instrumental" ){
    ## instrumental_variables edges ##
    new_edges <- draw_iv_edges(instrumental_variables = new_nodes, # new nodes as instrumental
                               treatments)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% instrumental_variables ] ) # existing node names in temporal order

  }else if( node_role %in% "competing_exposure" ){
    ## competing_exposure edges ##
    new_edges <- draw_competing_exposure_edges(outcomes,
                                               competing_exposure_vec = new_nodes) # new nodes as competing_exposure

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% competing_exposure_vec ] ) # existing node names in temporal order

  }else if( node_role %in% "collider" ){
    ## outcome edges ##
    new_edges <- draw_outcome_edges(type,
                                    outcomes,
                                    collider_vec = new_nodes) # new nodes as collider

    # connect colliders
    treatment_list <- list()

    treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

      treatment_list[x] <- lapply(1:length(new_nodes), function(y){

        list( c( ancestor = treatments[x], edge = "->", descendant = new_nodes[y]) )

      })

    }) )

    treatment_list <- Filter( Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE) )
    new_edges <- dplyr::bind_rows(new_edges, treatment_list)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% outcomes ] ) # existing node names in temporal order

  }else if( node_role %in% "latent" ){
    ## latent_variables edges ##
    new_edges <- draw_latent_edges(latent_variables = new_nodes) # new nodes as latent

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% latent_vec ] ) # existing node names in temporal order

    latent_vec <- c(latent_vec, new_nodes)

  }else if( node_role %in% "observed" ){
    ## connect observed to ancestors and descendants ##
    new_edges <- draw_observed_edges(observed = new_nodes, # new nodes as observed
                                       existing_dag)

    existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% observed ] ) # existing node names in temporal order

  }else{

    stop("Invalid node_role input. Check documentation for valid currently support inputs and try again.")

  }

  new_edges <- new_edges[ complete.cases(new_edges), ]

  edges <- merge(edges, new_edges,
                 by = c("ancestor", "edge", "descendant"),
                 all = TRUE) # combine both dag edges

  exclude_names <- c(treatments, outcomes, latent_vec)

  node_names <- dag_node_names[ !dag_node_names %in% exclude_names ]

  node_names <- c(node_names, new_node_names) # combine both dag node names

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

  if( any( !complete.cases(temporal_reference_node) ) ){

    temporal_reference_node <- existing_node_names[1] # first node is selected as the temporal reference node

  }

  #new_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% new_node_names ] ) # new node names in temporal order

  if( length(existing_node_names) == 0 ){

    existing_node_names <- confounder_vec

  }

  #new_node_names <- new_nodes[ !new_nodes %in% dag_node_names ] # remove duplicate node names

  new_node_names <- new_nodes

  new_coordinates <- add_merged_node_coordinates(dag = dag,
                                                 existing_node_names = existing_node_names,
                                                 new_node_names = new_nodes,
                                                 temporal_reference_node = temporal_reference_node,
                                                 coordinates = coordinates,
                                                 coords_spec = coords_spec[ complete.cases(coords_spec) ] )


  dagitty::coordinates(dag) <- new_coordinates


  return(dag)
}



#' Adds nodes to dagitty objects
#'
#' create_new_node_graph() is a helper function for add_nodes(). It adds new nodes to a dagitty object by referencing existing nodes.
#'
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


#' merge dags
#'
#' merge_dagitty() is a helper function for add_nodes() and merge_graph a daggity object using edges input.
#'
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
                                      coordinates = coordinates,
                                      coords_spec = coords_spec[ complete.cases(coords_spec) ] )


  dagitty::coordinates(dag) <- new_coordinates


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

