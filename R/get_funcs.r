
#' dagitty node edges and roles
#'
#' get_edges() filters a dagitty object and returns a data frame with edges for specified node roles.
#'
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param selected_nodes Nodes to return edges. Defaults to NULL, or can be a character or vector combination of any of the following: c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent", "observed"))
#' @returns A data frame, data table, or list of edges for the roles specified in selected_nodes.
#' @examples
#' edges <- get_edges(dag)
#'
#' @export
get_edges <- function(dag,
                      selected_nodes = c("outcome",
                                         "treatment",
                                         "confounder",
                                         "mediator",
                                         "mediator_outcome_confounder",
                                         "instrumental",
                                         "proxy",
                                         "competing_exposure",
                                         "collider",
                                         "latent",
                                         "observed")
                      ){
  .datatable.aware <- TRUE

  all_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent", "observed")
  if( identical(selected_nodes, all_roles) ){
    # seems unnecessary but putting this if-statement first prevents unnecessary checks, assuming the selected_nodes input is predominantly left blank
  }else if( any(grepl("!", selected_nodes)) ){
    # if any input includes "!", removes edges with matching roles in selected_nodes (for both parent and children nodes)
    all_nodes <- all_roles
    cleaned_nodes <- gsub("!", "", selected_nodes)
    selected_nodes <- all_nodes[!all_nodes %in% cleaned_nodes]
  }else{
    stop("Please check the selected_nodes input and try again.")
  }

  edges <- extract_unique_node_roles(dag) # add ancestor and descendant nodes (calls a function from later in this file)
  edges <- edges_longer(edges)
  edges <- edges[ edges[, "ancestor_role"] %in% selected_nodes, ]

  # handle simultaneous edges
  names(edges) <- c("v", "e", "w", "ancestor_role", "descendant_role")
  edges <- handle_simultaneous_edges(edges)
  names(edges) <- c("ancestor", "edge", "descendant", "ancestor_role", "descendant_role")

  return(edges)

}


#' dagitty nodes grouped by role
#'
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param multiple_roles Defaults to FALSE (one role per node). If set to TRUE, multiple roles can be returned for some nodes (e.g. latent mediator variable).
#' @return Nested list of nodes and node relationships
#' @examples
#' get_roles(dag)
#'
#' @export
get_roles <- function(dag, multiple_roles = FALSE){
  .datatable.aware <- TRUE

  if(multiple_roles == FALSE){

    edges_wide <- extract_unique_node_roles(dag)

  }else{

    edges_wide <- extract_node_roles(dag)

  }

  node_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent", "observed")
  num_roles <- length(node_roles)

  ## ancestor node edges to list ##
  ancestor_roles <- edges_wide[,1:14]
  ancestor_roles <- na.omit( reshape(ancestor_roles, varying = list(4:14), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  ancestor_roles <- ancestor_roles[order(ancestor_roles$id), c(1,4)]
  names(ancestor_roles)[1] <- "node"

  ## group by nodes
  unique_ancestors <- unique( ancestor_roles[,"node"] ) # vector of unique node names
  num_unique_ancestors <- nrow(unique_ancestors) # count of unique node names

  ## descendant node edges to list ##
  descendant_roles <- edges_wide[,c(1:3,15:25)]
  descendant_roles <- na.omit( reshape(descendant_roles, varying = list(4:14), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  descendant_roles <- descendant_roles[order(descendant_roles$id), c(3,4)]
  names(descendant_roles)[1] <- "node"

  ## find missing edges
  outcomes <- unique( descendant_roles[ descendant_roles$role == "outcome", ] )
  colliders <- unique( descendant_roles[ descendant_roles$role == "collider", ] )
  observed <- unique( descendant_roles[ descendant_roles$role == "observed", ] )
  latent <- unique( descendant_roles[ descendant_roles$role == "latent", ] )

  missing_outcomes <- outcomes[! unlist(outcomes[,1]) %in% unlist(unique_ancestors), ]
  missing_colliders <- colliders[! unlist(colliders[,1]) %in% unlist(unique_ancestors), ]
  missing_observed <- observed[! unlist(observed[,1]) %in% unlist(unique_ancestors), ]
  missing_latent <- latent[! unlist(latent[,1]) %in% unlist(unique_ancestors), ]

  all_roles <- rbind(ancestor_roles, missing_outcomes, missing_colliders, missing_observed, missing_latent)

  # edges grouped by each unique node in a list
  roles_list <- c()

  roles_list <- suppressWarnings( lapply(1:num_roles, function(x){

    roles_list <-  unique(unlist( all_roles[
      unlist(all_roles[,"role"]) %in% unlist(node_roles[x]) , ][, 1] ))

  }) )

  names(roles_list) <- node_roles # assign node roles as list element names

  roles_list <- lapply( roles_list, function(x) if( identical(x, character(0)) ) NA else x )


  return(roles_list)

}


#' dagitty node names and roles
#'
#' get_nodes() is a dagitty wrapper that returns a list of node names and their roles.
#'
#' @importFrom data.table as.data.table is.data.table
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @examples
#' get_nodes(dag)
#'
#' @export
get_nodes <- function(dag){
  .datatable.aware <- TRUE

  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  ## ancestor node edges to list ##
  edges_ancestors <- edges_wide[,1:14]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:14), idvar = "id",
                                      v.names = "role_ancestor", direction = "long")[,c("v", "e", "w", "role_ancestor", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), 1:4]
  names(edges_ancestors)[1:3] <- c("ancestor", "edge", "descendant")

  ## group by nodes
  unique_ancestors <- unique( unlist( edges_ancestors[,"ancestor"] ) )# vector of unique node names
  num_unique_ancestors <- length(unique_ancestors) # count of unique node names

  all_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent", "observed")

  # edges grouped by each unique node in a list
  edges_ancestor_list <- c()

  edges_ancestor_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges_ancestor_list <-  unique(unlist( edges_ancestors[
      unlist(edges_ancestors[,"ancestor"]) %in% unlist(unique_ancestors[x]) , ][, c(4)] ))

  }) )

  names(edges_ancestor_list) <- unlist(unique_ancestors) # assign ancestor node role to name of each element in edges list


  ## descendant node edges to list ##
  edges_descendants <- edges_wide[,c(1:3,15:25)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:14), idvar = "id",
                                        v.names = "role_descendant", direction = "long")[,c("v", "e", "w", "role_descendant", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), c(3,4)]

  ## find missing edges
  outcomes <- unique( edges_descendants[ edges_descendants$role_descendant == "outcome", ] )
  colliders <- unique( edges_descendants[ edges_descendants$role_descendant == "collider", ] )
  observed <- unique( edges_descendants[ edges_descendants$role_descendant == "observed", ] )
  latent <- unique( edges_descendants[ edges_descendants$role_descendant == "latent", ] )

  missing_outcomes <- outcomes[! unlist(outcomes[,1]) %in% unlist(unique_ancestors), ]
  missing_colliders <- colliders[! unlist(colliders[,1]) %in% unlist(unique_ancestors), ]
  missing_observed <- observed[! unlist(observed[,1]) %in% unlist(unique_ancestors), ]
  missing_latent <- latent[! unlist(latent[,1]) %in% unlist(unique_ancestors), ]

  missing_descendant_edges <- rbind(missing_outcomes, missing_colliders, missing_observed, missing_latent)

  unique_descendants <- unique(unlist(missing_descendant_edges[,1]))
  num_unique_descendants <- length(unique_descendants)

  # descendant edges grouped by each unique node in a list
  edges_descendant_list <- c()

  if( num_unique_descendants > 0 ){

    edges_descendant_list <- suppressWarnings( lapply(1:num_unique_descendants, function(x){

      edges_descendant_list <- unique(unlist(
        edges_descendants[ unlist(edges_descendants[,"w"]) %in% unlist(unique_descendants[x]), ][, c(2)] ))


    }) )

    names(edges_descendant_list) <- unlist(unique_descendants) # assign descendant node role to name of each element in edges list

  }



  edges_list <- c(edges_ancestor_list, edges_descendant_list) # combine ancestor and descendant node role lists


  return(edges_list)

}


#' dagitty node names, ancestor/descendant roles
#'
#' node_structure() is a dagitty wrapper function that returns a list of node names extrated from a dagitty object, including the name and role of their ancestor and descendant nodes.
#'
#' @importFrom data.table as.data.table is.data.table
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @examples
#' get_structure(dag)
#'
#' @export
get_structure <- function(dag){
  .datatable.aware <- TRUE

  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  unique_nodes <- unname(unlist(  unique(edges_wide[,1]) ))
  num_unique_nodes <- length(unique_nodes)

  ## ancestor node edges to list ##
  edges_ancestors <- edges_wide[,1:14]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:14), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "w", "role", "id")] )
  names(edges_ancestors)[1:2] <- c("ancestor", "descendant")
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), ][,c(1:3)]

  unique_ancestors <- unname(unlist( unique(edges_ancestors[,1]) ))
  num_unique_ancestors <- length(unique_nodes)

  ## descendant node edges to list ##
  edges_descendants <- edges_wide[,c(1:3,15:25)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:14), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "w", "role", "id")] )
  names(edges_descendants)[1:2] <- c("ancestor", "descendant")
  edges_descendants <- edges_descendants[order(edges_descendants$id), ][,c(1:3)]

  unique_descendant <- unname(unlist(  unique(edges_descendants[,1]) ))
  num_unique_descendant <- length(unique_nodes)

  # ancestor edges grouped by each unique node in a list
  edges_ancestor_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges_ancestors[ unlist(edges_ancestors[,"ancestor"]) %in% unlist(unique_ancestors[x]), ]

  }) )

  # descendant edges grouped by each unique node in a list
  edges_descendant_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges_descendants[ unlist(edges_descendants[,"descendant"]) %in% unlist(unique_ancestors[x]), ]

  }) )


  ## combine ancestor and descendant nodes ##
  edges_list <- Map(rbind, edges_descendant_list, edges_ancestor_list)

  names(edges_list) <- unlist(unique_nodes)

  ancestors_descendants_label_vec <- c("ancestor", "descendant")
  ## next create list with each node as an element
  edges_structure_list <- c()

  edges_structure_list <- suppressWarnings( lapply( 1:num_unique_ancestors , function(x){

    edges_structure_list <- lapply(1:2, function(n){ # edges grouped by unique node name

      edges_structure_list[x][[n]] <- edges_list[[x]][ edges_list[[x]][[n]] %in% names(edges_list[x]), ]

    })

    names(edges_structure_list) <- ancestors_descendants_label_vec

    edges_structure_list

  } ) )

  names(edges_structure_list) <- unlist(unique_ancestors) # assign ancestor node role to name of each element in edges list

  return(edges_structure_list)

}


#' Label DAG nodes
#' Extract node initials to use as labels in plot_dagitty, or dag plots using ggdag and ggplot2. This is essentially a wrapper function for ggdag::dag_label().
#' @importFrom ggdag dag_label
#' @param dag A dagittyy object.
#' @return Labelled vector of node names from the dagitty object
#' @export
get_labels <- function(dag){

  names <- as.data.frame(unique(suppressWarnings(ggdag::dag_label(dag)$data["name"])))
  names <- as.vector(names$name)
  labels <- gsub("(\\b[A-Z])[^A-Z]+", "\\1", names)
  names(labels) <- names

  return(labels)
}

