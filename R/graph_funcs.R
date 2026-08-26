.VALID_GRAPH_TYPES <- c("full", "saturated", "ordered")
.VALID_NODE_ROLES <- c("outcome", "treatment", "confounder", "mediator",
                       "mediator_outcome_confounder", "instrument",
                       "competing_cause", "collider", "latent", "observed")
.ALL_NODE_ROLES <- c(.VALID_NODE_ROLES, "undetermined")


#' Build dagitty objects
#'
#' build_graph() produces a dagitty graph object from inputted parameters (e.g. treatments, outcome, confounders). A node listed in more than one role input receives each listed role's arrows.
#'
#' @importFrom dagitty is.dagitty children coordinates
#' @param variables Dagitty object, or vector of variable names, e.g. "Z" or c("Z1", "Z2", "Z3"). If variable names are inputted the order determines the assigned coordinates. A list can also be supplied. Variables inputted are treated as confounders. If type = "ordered", confounders located in the same list will be assigned similar coordinates. When a dagitty object is supplied, a new DAG is constructed from its inferred node roles; use graph-editing functions when the original graph type and edge semantics must be retained.
#' @param treatments Treatment variable name, e.g. "X". Must be specified.
#' @param outcomes Outcome variable name, e.g. "Y". Must be specified.
#' @param mediators Character or vector of mediator variable names, e.g. "M" or c("M1", "M2", "M3").
#' @param latent_variables Character or vector of additional or already supplied latent (unobserved) variable names, e.g. "U" or c("U1", "U2", "M1").
#' @param instrumental_variables Vector of instrumental variable names, e.g. "IV"
#' @param mediator_outcome_confounders Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param competing_causes Vector of competing cause names. An arrow is drawn from each competing cause to the outcome, with other arrows depending on \code{type}.
#' @param colliders Vector of collider variables, with both treatment and outcome parents.
#' @param type Type of connected graph (e.g. "full", "saturated", "ordered"). Defaults to 'full' (fully connected graph) with bi-directional arrows drawn between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)) and arrows from confounders to mediators. If type = 'saturated', a similar saturated graph is produced except confounders are not connected to mediators or each other bidirectionally. When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennant et al. (2021).
#' @param confounders Optional alias for \code{variables} when supplying confounder names.
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
#' ## Option 1: an 'ordered' graph uses the supplied vector order to determine
#' ## the temporal order of confounder nodes.
#'
#' type <- "ordered"
#'
#' set.seed(1)
#' dag <- build_graph(variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes,
#'                   type = type)
#'
#' # Plot the graph and assign coordinates based on the variable roles.
#' plot_dagitty(dag)
#'
#'
#' ## Option 2: a 'saturated' graph connects each supplied confounder.
#'
#' type <- "saturated"
#'
#' set.seed(1)
#' dag <- build_graph(variables = variables,
#'                   treatments = treatments,
#'                   outcomes = outcomes,
#'                   type = type)
#'
#' plot_dagitty(dag)
#'
#'
#' ## Option 3: 'full' connects confounders bidirectionally and connects them
#' ## to all other nodes except treatment and outcome.
#'
#' type <- "full"
#'
#' set.seed(1)
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
                       competing_causes = NA,
                       colliders = NA,
                       type = "full",
                       confounders = NULL
                       ){

  if( !is.null(confounders) ){
    if( !missing(variables) && !(length(variables) == 1L && is.atomic(variables) && is.na(variables)) ){
      stop("Supply confounder names using either 'variables' or 'confounders', not both.", call. = FALSE)
    }
    variables <- confounders
  }

  if( !type %in% .VALID_GRAPH_TYPES ){
    stop("Invalid type: \"", type, "\". Must be one of: ",
         paste(.VALID_GRAPH_TYPES, collapse = ", "), ".", call. = FALSE)
  }

  if( dagitty::is.dagitty(variables) ){ # check for an existing dag (inputted as confounders)

    existing_dag <- variables
    existing_graph_type <- .dag_graph_type(existing_dag)

    if( existing_graph_type != "dag" ){
      warning("build_graph() constructs a new DAG from inferred node roles; ",
              "the input's mixed graph type and edge semantics are not retained. ",
              "Use keep_edges(), add_nodes(), connect_nodes(), or join_graphs() ",
              "to edit a pdag or mag while retaining its type.", call. = FALSE)
    }

    if( type == "full" ){
      warning("build_graph(type = \"full\") rebuilds an existing graph from ",
              "inferred node roles and does not preserve its directed edges. ",
              "Use connect_nodes(dag, type = \"full\") to retain existing edges.",
              call. = FALSE)
    }

    node_roles <- .get_roles(existing_dag)

    confounders <- variables <- node_roles$confounder

    treatments <- node_roles$treatment

    outcomes <- node_roles$outcome

    m_o_confounder_vec <- mediator_outcome_confounders <- node_roles$mediator_outcome_confounder

    mediator_vec <- mediators <- node_roles$mediator

    instrumental_variables <- node_roles$instrument

    collider_vec <- colliders <- node_roles$collider

    competing_cause_vec <- competing_causes <- node_roles$competing_cause

    observed <- unique(c(node_roles$observed, node_roles$undetermined))

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

    competing_cause_vec <- as.vector( unlist( lapply( competing_causes, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    collider_vec <- as.vector( unlist( lapply( colliders, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

    observed <- NA

    if( all( complete.cases( unlist(confounders) ) ) & any( unlist(confounders) %in% m_o_confounder_vec ) ){ # if any mediator-outcome confounders are also inputted as confounders, execution is stopped

      stop("Node names inputted in the 'variables' parameter detected in 'mediator_outcome_confounder'. These inputs should be mutually exclusive. Please adjust inputs and try again.")
    }

  }else{

    stop("The 'treatments' and 'outcomes' inputs should be used if a DAG is not provided in the 'variables' input. Please adjust inputs and try again.")
  }

  ## get variable names ##
  observed_node_names <- unique( as.vector( c(unlist(confounders), m_o_confounder_vec, mediator_vec, competing_cause_vec, collider_vec, instrumental_variables, unlist(observed)) ) )
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
                      competing_cause_vec = competing_cause_vec,
                      collider_vec = collider_vec,
                      observed = observed,
                      observed_node_names = observed_node_names,
                      existing_dag = existing_dag)

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


#' Add nodes to a dagitty object
#'
#' Add nodes to a dagitty object, connecting edges based on the 'type' of graph selected, and generate new node coordinates using existing nodes.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_nodes Inputted vector of node names to be added to the graph.
#' @param ancestors Optional character vector of existing nodes to set as ancestors (parents) of the new nodes. Defaults to NULL.
#' @param descendants Optional character vector of existing nodes to set as descendants (children) of the new nodes. Defaults to NULL.
#' @param node_role Optional role for the new nodes. When \code{NULL}, each new
#'   node is connected according to \code{ancestors}/\code{descendants}, or to
#'   every existing node when neither is supplied. A \code{node_role} cannot be
#'   combined with explicit \code{ancestors}/\code{descendants}.
#' @param type Graph structure: \code{"full"}, \code{"saturated"}, or
#'   \code{"ordered"}. Defaults to \code{"saturated"}. The legacy values
#'   \code{"first"} and \code{"last"} remain available for compatibility;
#'   \code{position} provides a separate, newer ordering API.
#' @param position Optional placement within the relevant temporal order,
#'   expressed with \code{first()}, \code{last()},
#'   \code{before()}, and \code{after()} (see
#'   \code{\link{position_helpers}}). Supported for node_role
#'   \code{"confounder"}, \code{"treatment"}, and \code{"mediator"}.
#' @param print_edges Print new edges in console, defaults to TRUE.
#' @returns A dagitty object.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' dag <- add_nodes(dag, new_nodes = "W", node_role = "confounder")
#'
#' @export
add_nodes <- function(dag,
                      new_nodes,
                      ancestors = NULL,
                      descendants = NULL,
                      node_role = NULL,
                      type = "saturated",
                      print_edges = TRUE,
                      position = NULL
                      ){
  .datatable.aware <- TRUE

  if( type %in% c("first", "last") ){
    warning("type = \"", type, "\" is retained for compatibility; ",
            "'position' provides the newer ordering API.", call. = FALSE)
  }
  if( !type %in% c(.VALID_GRAPH_TYPES, "first", "last") ){
    stop("Invalid type: \"", type, "\". Must be one of: ",
         paste(.VALID_GRAPH_TYPES, collapse = ", "), ".", call. = FALSE)
  }
  if( length(node_role) > 0L && ( length(ancestors) > 0L || length(descendants) > 0L ) ){
    stop("add_nodes(): supply either a node_role or explicit ancestors/descendants, not both. ",
         "The role template determines the new edges, so the explicit targets would be ignored. ",
         "Please adjust inputs and try again.", call. = FALSE)
  }

  if( length(new_nodes) == 0L ){
    stop("add_nodes(): supply at least one node name in new_nodes.",
         call. = FALSE)
  }
  if( anyNA(new_nodes) || any(!nzchar(as.character(new_nodes))) ){
    stop("add_nodes(): new_nodes must contain non-missing, non-empty names.",
         call. = FALSE)
  }

  ## get initial dag edges ##
  edges <- data.table::as.data.table(.dag_edges(dag, "add_nodes()"))

  names(edges) <- c("v", "e", "w") # change column names

  dag_node_names <- names(dag) # extract dag node names

  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  already_present <- unique(new_nodes[new_nodes %in% dag_node_names])
  if( length(already_present) == length(unique(new_nodes)) ){
    stop("add_nodes(): all new_nodes are already in the graph: ",
         paste(already_present, collapse = ", "),
         ". Use connect_nodes() to wire existing nodes.", call. = FALSE)
  }
  if( length(already_present) > 0L ){
    warning("add_nodes(): skipping names already in the graph: ",
            paste(already_present, collapse = ", "), ".", call. = FALSE)
  }

  new_node_names <- unique(new_nodes[ !new_nodes %in% dag_node_names ])

  new_edges <- connect_new_nodes(dag = dag,
                                 new_nodes = new_node_names,
                                 ancestors = ancestors,
                                 descendants = descendants,
                                 node_roles = node_role,
                                 type = type,
                                 position = position)

  ## pre-process before merging ##
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs

  ## merge new and existing dag edges ##
  edges <- merge(edges, new_edges,
                 by = c("v", "e", "w"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges
  edges <- edges[edges$v != edges$w, ] # remove new_edges with identical ancestor and descendant node names
  updated_treatments <- treatments(dag)
  updated_outcomes <- .outcomes(dag)
  updated_latents <- unobserved(dag)
  if( identical(node_role, "treatment") ){
    updated_treatments <- union(updated_treatments, new_node_names)
  }else if( identical(node_role, "outcome") ){
    updated_outcomes <- union(updated_outcomes, new_node_names)
  }else if( identical(node_role, "latent") ){
    updated_latents <- union(updated_latents, new_node_names)
  }
  dag <- rebuild_dag(
    dag, edges,
    treatments = updated_treatments,
    outcomes = updated_outcomes,
    latent_vec = updated_latents,
    extra_nodes = new_node_names
  )

  if( print_edges == TRUE){
    new_edges_list <- print_edges_helper(new_edges)
    cat( paste("c(", paste( unlist(new_edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    message("\nPrinted new edges - copy and paste to use as a vector object.\n")
  }

  new_coordinates <- renew_coords(dag = dag,
                                  new_node_names = new_node_names,
                                  coordinates = coordinates)
  dagitty::coordinates(dag) <- new_coordinates

  return(dag)
}


#' Connect existing nodes in a dagitty object
#'
#' Connects existing nodes according to the selected graph type and optional
#' causal role while preserving graph metadata and coordinates.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param nodes Existing node names to connect. Defaults to all graph nodes.
#' @param node_role Optional role controlling how the supplied nodes are
#'   connected. The default (\code{NULL}) connects them to all existing nodes.
#' @param type Graph structure: \code{"full"}, \code{"saturated"}, or
#'   \code{"ordered"}. Defaults to \code{"full"}. The legacy values
#'   \code{"first"} and \code{"last"} remain available for compatibility;
#'   \code{position} provides a separate, newer ordering API.
#' @param position Optional placement of the selected nodes within their
#'   \code{node_role} temporal order, expressed with \code{first()},
#'   \code{last()}, \code{before()}, and
#'   \code{after()} (see \code{\link{position_helpers}}). When
#'   supplied, new edges follow the resolved order; \code{NULL} retains the
#'   default behaviour. Requires a \code{node_role}.
#' @param print_edges Print new edges in console, defaults to TRUE.
#' @returns A dagitty object.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' dag <- connect_nodes(dag, nodes = "Z1")
#'
#' @export
connect_nodes <- function(dag,
                           nodes = NULL,
                           node_role = NULL,
                           type = "full",
                           print_edges = TRUE,
                           position = NULL
                           ){
  .datatable.aware <- TRUE

  if( type %in% c("first", "last") ){
    warning("type = \"", type, "\" is retained for compatibility; ",
            "'position' provides the newer ordering API.", call. = FALSE)
  }
  if( !type %in% c(.VALID_GRAPH_TYPES, "first", "last") ){
    stop("Invalid type: \"", type, "\". Must be one of: ",
         paste(.VALID_GRAPH_TYPES, collapse = ", "), ".", call. = FALSE)
  }
  ## get initial dag edges ##
  edges <- data.table::as.data.table(.dag_edges(dag, "connect_nodes()"))
  dag_node_names <- names(dag) # extract dag node names

  if( length(nodes) > 0L ){
    missing_nodes <- unique(nodes[!nodes %in% dag_node_names])
    if( length(missing_nodes) > 0L ){
      stop("connect_nodes(): nodes not in the graph: ",
           paste(missing_nodes, collapse = ", "),
           ". Use add_nodes() to introduce new nodes.", call. = FALSE)
    }
  }

  if( length(nodes) == 0 ){

    nodes <- dag_node_names
  }

  if( length(node_role) == 0 ){
    if( !is.null(position) ){
      stop("connect_nodes(): 'position' requires a 'node_role'. ",
           "Supply a node_role, or remove position.", call. = FALSE)
    }
    new_edges <- connect_nodes_helper(dag, nodes, dag_node_names, type)
  }else{
    ## call helper function to draw new node edges etc. ##
    output_list <- add_nodes_helper(dag, nodes, node_role, type, position = position)

    new_edges <- output_list$new_edges

    colnames(new_edges) <- c("v", "e", "w")
  }

  ## pre-process before merging ##
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs


  ## merge new and existing dag edges ##
  edges <- merge(edges, new_edges,
                 by = c("v", "e", "w"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges

  edges <- edges[edges$v != edges$w, ] # remove new_edges with identical ancestor and descendant node names

  dag <- rebuild_dag(dag, edges)

  if( isFALSE(dagitty::isAcyclic(dag)) ){
    warning("The connected graph contains a directed cycle.", call. = FALSE)
  }

  if( print_edges == TRUE){

    new_edges_list <- print_edges_helper(new_edges)

    cat( paste("c(", paste( unlist(new_edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    message("\nPrinted new edges to assess. - copy and paste in a .R file to use as a vector object.\n")
  }

  return(dag)
}




#' Merge two dagitty objects
#'
#' join_graphs() adds new nodes to the first from the second supplied dagitty object, before generating new coordinates for the merged graph.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty edges coordinates
#' @param dag First dagitty object.
#' @param new_dag Second dagitty object, added to the first.
#' @returns A dagitty object
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y", type = "saturated"
#' )
#' dag2 <- build_graph(
#'   variables = c("Z1", "Z3"), treatments = "X", outcomes = "Y", type = "saturated"
#' )
#' new_graph <- join_graphs(dag, dag2)
#'
#' @export
join_graphs <- function(dag,
                        new_dag
                        ){
  .datatable.aware <- TRUE

  edges <- data.table::as.data.table(.dag_edges(dag, "join_graphs()"))
  node_names <- names(dag) # extract dag node names

  treatments <- treatments(dag) # get treatments
  outcomes <- .outcomes(dag) # get outcomes
  latent_vec <- unlist(unobserved(dag)) # get dag latent variables
  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  new_edges <- data.table::as.data.table(.dag_edges(new_dag, "join_graphs()"))
  new_node_names <- names(new_dag) # extract new dag node names

  existing_node_names <-  new_node_names[ new_node_names %in% node_names ] # saves duplicate node names
  shared_node_names <- existing_node_names # shared nodes, before coordinate filtering
  existing_node_coords_y <- coordinates$y[ names(coordinates$y)  %in% existing_node_names] # saves duplicate node names
  existing_node_names <- names( existing_node_coords_y[ order(existing_node_coords_y) ] ) # existing dag node names in ascending y-coords order

  new_node_names <- new_node_names[ !new_node_names %in% node_names ] # remove duplicate node names


  new_dag_treatments <- treatments(new_dag)
  treatments <- c(treatments,
                  new_dag_treatments[ !new_dag_treatments %in% node_names] ) # get treatments

  new_dag_outcomes <- .outcomes(new_dag) # get outcomes
  outcomes <- c(outcomes,
                new_dag_outcomes[ !new_dag_outcomes %in% node_names] ) # get treatments

  new_dag_latent_vec <- unlist( unobserved(new_dag) ) # get dag latent variables
  latent_vec <- c(latent_vec,
                  new_dag_latent_vec[ !new_dag_latent_vec %in% node_names] ) # get treatments
  adjusted_vec <- union(dagitty::adjustedNodes(dag),
                        dagitty::adjustedNodes(new_dag))
  selected_vec <- union(.selected_nodes(dag), .selected_nodes(new_dag))

  ## warn when the two graphs disagree on a shared node's declarations ##
  if( length(shared_node_names) > 0 ){
    conflicted <- Filter(function(nm){
      ( nm %in% treatments ) != ( nm %in% new_dag_treatments ) |
      ( nm %in% outcomes ) != ( nm %in% new_dag_outcomes ) |
      ( nm %in% latent_vec ) != ( nm %in% new_dag_latent_vec )
    }, shared_node_names)
    if( length(conflicted) > 0 ){
      warning("join_graphs(): the two graphs disagree on the exposure/outcome/latent status of shared node(s): ",
              paste(conflicted, collapse = ", "),
              ". The first graph's declarations were kept. Please check that this is intended.",
              call. = FALSE)
    }
  }


  coordinates_new_dag <- dagitty::coordinates(new_dag) # extract new dag coordinates
  new_node_coords_y <- coordinates_new_dag$y[ names(coordinates_new_dag$y) %in% new_node_names] # saves duplicate node names
  if( length(new_node_coords_y) > 0 ){
    new_node_names <- names( new_node_coords_y[ order(new_node_coords_y) ] ) # new dag node names in ascending y-coords order
  }

  edges <- unique(rbind(edges, new_edges)) # combine both dag edges, dropping duplicates

  ## rebuild dag ##
  dag <- rebuild_dag(dag, edges,
                     treatments = treatments,
                     outcomes = outcomes,
                     latent_vec = latent_vec,
                     extra_nodes = names(new_dag),
                     adjusted_vec = adjusted_vec,
                     selected_vec = selected_vec)

  if( isFALSE(dagitty::isAcyclic(dag)) ){
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.",
            call. = FALSE)
  }

  if( length(new_node_names) == 0L ){
    return(dag)
  }

  coordinates <- renew_coords(dag = dag,
                              new_node_names = new_node_names,
                              coordinates = coordinates)
  dagitty::coordinates(dag) <- coordinates

  return(dag)
}


#' Assess dagitty object edges
#'
#' assess_edges() provides ways to assess connected edges based on causal criteria and/or user inputs.
#'
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param edges_to_assess Defaults to "all" edges, "bidirectional" includes only bidirected ("<->") and undirected ("--") edges.
#' @param edges_to_keep Edges to be kept, e.g. c("Z1 -> Z2", "Z2 -> Z3"). Accepts a character vector of edge strings, a v/e/w edge table, or a dagitty object. Only listed edges are kept, including any into or out of the treatments and outcomes.
#' @param assess_causal_criteria Defaults to FALSE. If TRUE, the user is guided through a sequence that assesses each pair of connected nodes using causal criteria. Based on the Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs) from Ferguson et al. (2020).
#' @param causal_criteria Defaults to "ESCDAGs" (Ferguson et al., 2020) causal criteria considers temporality, face-validity, and recourse to theory. Other causal criteria can be supplied as a data.frame, see summary(ESCDAGs) for column names.
#' @param save_answers Set to \code{TRUE} to return the decision log with the
#'   retained edges, or supply a previous log (v, e, w followed by one column per
#'   criterion) to replay it without prompting.
#' @returns A list or vector of edges.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2", "Z3"), treatments = "X", outcomes = "Y", type = "ordered"
#' )
#'
#' ## a saturated graph built from the existing (sparser) dag
#' saturated_graph <- build_graph(type = "saturated", variables = dag)
#'
#' ## Option 1: keep the edges already present in the existing dag
#' edges <- assess_edges(saturated_graph, edges_to_keep = dag)
#'
#' \dontrun{
#' ## Option 2: interactively rule edges in/out via the guided causal-criteria sequence
#' edges_to_keep <- assess_edges(saturated_graph, edges_to_keep = dag,
#'                               assess_causal_criteria = TRUE)
#' }
#'
#' @export
assess_edges <- function(dag,
                         edges_to_assess = "all",
                         edges_to_keep = NA,
                         assess_causal_criteria = FALSE,
                         save_answers = NULL,
                         causal_criteria = ESCDAGs
                         ){
  .datatable.aware <- TRUE

  allowed_modes <- c("all", "bidirectional", "bidirected",
                     "bi-directional", "bi-directed")
  if( !is.character(edges_to_assess) || length(edges_to_assess) != 1L ||
      is.na(edges_to_assess) || !edges_to_assess %in% allowed_modes ){
    stop("'edges_to_assess' must be one of: ",
         paste(allowed_modes, collapse = ", "), ".", call. = FALSE)
  }

  # get dag edges
  edges <- data.table::as.data.table(.dag_edges(dag, "assess_edges()"))

  if( dagitty::is.dagitty(edges_to_keep) ){

    edges_to_keep <- data.table::as.data.table(
      .dag_edges(edges_to_keep, "assess_edges()")
    )
    edges <-  edges[!edges_to_keep, on = c("v", "e", "w")]

  }else if( all( complete.cases( unlist(edges_to_keep) ) ) ){

    if( is.character(edges_to_keep) ){
      edges_to_keep <- .parse_edge_strings(edges_to_keep, "assess_edges()")
    }

    if( is.data.frame(edges_to_keep) || data.table::is.data.table(edges_to_keep) ){
      edges_to_keep <- data.table::as.data.table(
        .normalise_edge_frame(edges_to_keep, "assess_edges()")
      )
      edges <-  edges[!edges_to_keep, on = c("v", "e", "w")]
    }else{
      stop("'edges_to_keep' must be a data frame, data table, vector, or dagitty object. Please check inputs and try again.")
    }

  }

  if( edges_to_assess == "bidirectional" | edges_to_assess == "bidirected" |
      edges_to_assess == "bi-directional" | edges_to_assess == "bi-directed" ){ # if use inputs edges_to_assess = "bidirectional"

    edges_to_assess <- edges[ unlist(edges[, "e"]) %in% c("<->", "--"), ] # identify all bidirected ("<->") and undirected ("--") edges directly by edge symbol.

    # remove edges_to_assess from dag edges
    edges <- data.table::setDT(edges)[!edges_to_assess, on = c("v", "e", "w")]

    if( all( complete.cases( edges_to_keep ) ) ){ # check if edges_to_keep specified
      edges_to_keep <- rbind(edges, edges_to_keep) # combine unidirectional edges and edges_to_keep
    }else{ # otherwise unidirectional edges become edges_to_keep
      edges_to_keep <- edges
    }
    edges <- edges_to_assess # the detected bidirectional edges are what gets assessed
  }

  if( assess_causal_criteria == TRUE ){

    edges <- causal_criteria_sequence(edges = edges,
                                      check_skip_sequence = FALSE,
                                      edges_to_keep = edges_to_keep,
                                      causal_criteria = causal_criteria,
                                      save_answers = save_answers)

    return(edges)
  }

  if( nrow(edges) != 0 ){
    if( length( unlist(edges_to_keep) ) > 0 ){
      if( all( complete.cases( unlist(edges_to_keep) ) ) ){

        num_edges <- nrow(edges_to_keep)

        edges_to_keep <- suppressWarnings( sapply(1:num_edges, function(x){ # collapse edges_to_keep to a vector
          edges_to_keep <- paste( edges_to_keep[x,], collapse=" ")
        }) )

      }
    }
    edges_list <- print_edges_helper(edges) # group edges by unique node names
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
  edges_list <- list( edges_to_assess = edges,
                      edges_to_keep = edges_to_keep )
  message("\nOutputting list containing 'edges_to_assess' and 'edges_to_keep'.
          \n\nPrinted edges to assess - copy and paste in a .R file to use as a vector object.")

  return(edges_list)
}


#' Remove dagitty object edges
#'
#' keep_edges() removes edges based on user inputs.
#'
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom dagitty edges coordinates dagitty isAcyclic
#' @param dag A saturated graph dagitty object. Exposure and outcome must be indicated, and optionally can include assigned coordinates.
#' @param edges_to_keep Edges to be kept, e.g. c("Z1 -> Z2", "Z2 -> Z3"). Accepts a character vector of edge strings, a v/e/w edge table, or a dagitty object. Only listed edges are kept, including any into or out of the treatments and outcomes.
#' @returns A dagitty object, with directed arrows removed based on edges_to_keep.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2", "Z3"), treatments = "X", outcomes = "Y", type = "ordered"
#' )
#'
#' ## a saturated graph built from the existing (sparser) dag
#' saturated_graph <- build_graph(type = "saturated", variables = dag)
#'
#' ## keep only the edges present in the original dag
#' kept <- keep_edges(dag, saturated_graph)
#'
#' \dontrun{
#' ## interactive: curate edges via the guided causal-criteria sequence, then keep them
#' edges_to_keep <- assess_edges(saturated_graph, edges_to_keep = dag,
#'                               assess_causal_criteria = TRUE)
#' kept <- keep_edges(edges_to_keep, saturated_graph)
#' }
#'
#' @export
keep_edges <- function(edges_to_keep,
                       dag
                       ){
  edges <- .dag_edges(dag, "keep_edges()")

  if( dagitty::is.dagitty(edges_to_keep)){
    edges_to_keep <- .dag_edges(edges_to_keep, "keep_edges()")
  }else if( is.character(edges_to_keep) ){
    edges_to_keep <- .parse_edge_strings(edges_to_keep, "keep_edges()")
  }else if( !(is.data.frame(edges_to_keep) ||
               data.table::is.data.table(edges_to_keep)) ){
    stop("'edges_to_keep' must be a character vector, edge table, or dagitty object.",
         call. = FALSE)
  }

  edges_to_keep <- .normalise_edge_frame(edges_to_keep, "keep_edges()")
  if( nrow(edges_to_keep) == 0L ){
    stop("keep_edges(): no edges were supplied.", call. = FALSE)
  }
  if( anyNA(edges_to_keep) ){
    warning("Rows of 'edges_to_keep' containing NA were ignored.", call. = FALSE)
    edges_to_keep <- edges_to_keep[complete.cases(edges_to_keep), ,
                                   drop = FALSE]
  }
  if( nrow(edges_to_keep) == 0L ){
    stop("keep_edges(): no complete edges were supplied.", call. = FALSE)
  }

  requested_keys <- .edge_key(edges_to_keep)
  existing_keys <- .edge_key(edges)
  missing_edges <- edges_to_keep[!requested_keys %in% existing_keys, ,
                                 drop = FALSE]
  if( nrow(missing_edges) > 0L ){
    warning("Edges not present in the dag were ignored: ",
            paste(missing_edges$v, missing_edges$e, missing_edges$w,
                  collapse = ", "),
            ".", call. = FALSE)
  }

  retained <- edges[existing_keys %in% requested_keys, , drop = FALSE]
  dag <- rebuild_dag(dag, retained)

  if( isFALSE(dagitty::isAcyclic(dag)) ){
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.")
  }

  dag
}


.assessment_is_interactive <- function(){
  interactive()
}

.assessment_readline <- function(prompt){
  readline(prompt)
}


#' Assess graph edges using causal criteria
#' @importFrom data.table data.table
#' @param edges vector of edges whose relationships are to be assessed
#' @param num_edges number of edges to be assessed
#' @param check_skip_sequence TRUE or FALSE depending on prior inputs
#' @param causal_criteria Set of causal criteria to be used. Can be user-specified, defaults to 'ESCDAGs'.
#' @noRd
causal_criteria_sequence <- function(edges,
                                     check_skip_sequence,
                                     edges_to_keep,
                                     causal_criteria,
                                     save_answers
                                     ){
  .datatable.aware <- TRUE
  edges_list <- edges   # ensure the final return() is always bound; degenerate skip path returns edges unchanged

  if( length(save_answers) >= 3 ){

    if( nrow(save_answers) == nrow(edges) ){

      len_answers <- length(save_answers)
      if( all( as.data.frame(save_answers)[,1:3] == as.data.frame(edges) ) & len_answers == (nrow(causal_criteria) + 3) ){
        answers <- data.table::as.data.table(save_answers)[, 4:len_answers] # coerce first so a single criterion keeps its column name
        cols_required <- unlist(causal_criteria[,"name"])[ unlist(causal_criteria[, "required" ]) == "yes" ]
        answers <- answers[, ..cols_required]

        keep_edge <- vapply(seq_len(nrow(save_answers)), function(x){
          isTRUE(all((answers == "y")[x, ]))
        }, logical(1L))
        edges <- edges[keep_edge, ]

        if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
            all(complete.cases(edges_to_keep)) ){
          edges <- rbind(edges_to_keep[, 1:3], edges)
        }

        message("\nOutputting edges. Causal criteria decision log printed.")
        print(answers)

        return(edges)
      }

    }

    stop("Length of answers does not match the causal criteria provided. Please check inputs and try again.")
  }

  if( !.assessment_is_interactive() ){
    stop("Interactive causal assessment requires an interactive R session. ",
         "Supply 'save_answers' to replay a completed assessment.",
         call. = FALSE)
  }

  if(check_skip_sequence == FALSE) {

    num_edges <- nrow(edges)
    if( num_edges < 1 ){
      message( "There are no directed arrows to assess." )
      if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
          all(complete.cases(edges_to_keep)) ){
        edges <- rbind(edges_to_keep[, 1:3], edges)
      }

      return(edges)
    }
    set_noun <- if( any( unlist( edges[, "e"] ) %in% c("<->", "--") ) ) "edges" else "directed arrows"
    cat( "\nThere are ", num_edges, " ", set_noun, " to be assessed.","\n", "\n", sep="" )
    print(edges, quote=FALSE)
    cat( "\nAssess the posited causal relationships using causal criteria?", "\n" )

    check_ans <- FALSE

    while(check_ans == FALSE){

      choice <- .assessment_readline("(y/n/?info): ")
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
        arrow_noun <- if( unlist( edges[ arrow_count, "e" ] ) %in% c("<->", "--") ) "edge" else "directed arrow"
        cat( "For the ", arrow_noun, " '", arrow, "' consider each of the following:", "\n", "\n", sep="")

        while( check_ans == FALSE ){

          if( !criterion_num > num_criteria ){

            cat( paste0("'", arrow, "' (", arrow_count, "/", num_edges, ")"),
                 "\n", paste0("\n[", criterion_num, "/", num_criteria, "]"),
                 paste0( causal_criteria[ criterion_num, "name" ], ":"),
                 causal_criteria[ criterion_num, "question" ], "\n",
                 "\nFor help, enter \'?info\'")
            choice <- .assessment_readline("(y/n/?info): ")
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
              cat( "\n", "\nCausal criterion", criterion_num, "\u2014", causal_criteria[ criterion_num, "name" ],
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

      }

    if(num_arrow_to_remove > 0){

      criteria_answers <- cbind(edges, criteria_answers)
      edges <- edges[ -removed_arrows, ]

      if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
          all(complete.cases(edges_to_keep)) ){
        edges <- rbind(edges_to_keep[, 1:3], edges)
      }

      if( length(save_answers) == 0 ){
        message("\nOutputting edges. Causal criteria decision log printed.")
        print(criteria_answers)

        return(edges)
      }
      edges <- list(edges = edges,
                    answers = criteria_answers)

      message("\nOutputting list containing 'edges' and 'answers'.")

      return(edges)
    }else{

      message( "No arrows were removed." )

      criteria_answers <- cbind(edges, criteria_answers)

      if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
          all(complete.cases(edges_to_keep)) ){
        edges <- rbind(edges_to_keep[, 1:3], edges)
      }

      if( length(save_answers) == 0 ){
        print(criteria_answers)
        return(edges)
      }

      edges_list <- list(edges = edges,
                         answers = criteria_answers)

      message("\nOutputting list containing 'edges' and 'answers'.")
      return(edges_list)

    }

  }

  return(edges_list)
}
