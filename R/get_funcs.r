
#' dagitty node edges and roles
#'
#' get_edges() filters a dagitty object and returns a data frame with edges for specified node roles.
#'
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param selected_nodes Nodes to return edges. Defaults to NULL, or can be a character or vector combination of any of the following: c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed"))
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
                                         "competing_exposure",
                                         "collider",
                                         "latent",
                                         "observed")
                      ){
  .datatable.aware <- TRUE

  all_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed")
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

  node_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed")
  num_roles <- length(node_roles)

  ## ancestor node edges to list ##
  ancestor_roles <- edges_wide[,1:13]
  ancestor_roles <- na.omit( reshape(ancestor_roles, varying = list(4:13), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  ancestor_roles <- ancestor_roles[order(ancestor_roles$id), c(1,4)]
  names(ancestor_roles)[1] <- "node"

  ## group by nodes
  unique_ancestors <- unique( ancestor_roles[,"node"] ) # vector of unique node names
  num_unique_ancestors <- nrow(unique_ancestors) # count of unique node names

  ## descendant node edges to list ##
  descendant_roles <- edges_wide[,c(1:3,14:23)]
  descendant_roles <- na.omit( reshape(descendant_roles, varying = list(4:13), idvar = "id",
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
  edges_ancestors <- edges_wide[,1:13]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:13), idvar = "id",
                                      v.names = "role_ancestor", direction = "long")[,c("v", "e", "w", "role_ancestor", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), 1:4]
  names(edges_ancestors)[1:3] <- c("ancestor", "edge", "descendant")

  ## group by nodes
  unique_ancestors <- unique( unlist( edges_ancestors[,"ancestor"] ) )# vector of unique node names
  num_unique_ancestors <- length(unique_ancestors) # count of unique node names

  all_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed")

  # edges grouped by each unique node in a list
  edges_ancestor_list <- c()

  edges_ancestor_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges_ancestor_list <-  unique(unlist( edges_ancestors[
      unlist(edges_ancestors[,"ancestor"]) %in% unlist(unique_ancestors[x]) , ][, c(4)] ))

  }) )

  names(edges_ancestor_list) <- unlist(unique_ancestors) # assign ancestor node role to name of each element in edges list


  ## descendant node edges to list ##
  edges_descendants <- edges_wide[,c(1:3,14:23)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:13), idvar = "id",
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
#' @param dag A dagitty object.
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
  edges_ancestors <- edges_wide[,1:13]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:13), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "w", "role", "id")] )
  names(edges_ancestors)[1:2] <- c("ancestor", "descendant")
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), ][,c(1:3)]

  unique_ancestors <- unname(unlist( unique(edges_ancestors[,1]) ))
  num_unique_ancestors <- length(unique_nodes)

  ## descendant node edges to list ##
  edges_descendants <- edges_wide[,c(1:3,14:23)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:13), idvar = "id",
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
#' @param label_type Defaults to using the node names of the inputted dag, or 'initials' if specified. Labels are only assigned by plot_dagitty() if nodes are not already labelled, using a 'ggdag::dag_label()' wrapper..
#' @return Labelled vector of node names from the dagitty object
#' @export
get_labels <- function(dag, label_type){

  if(label_type == "name"){

    labels <- as.data.frame(unique(suppressWarnings(ggdag::dag_label(dag)$data["name"])))
    labels <- as.vector(labels$name)
    names(labels) <- labels

    return(labels)
  }

  names <- as.vector(unlist(unique(suppressWarnings(ggdag::dag_label(dag)$data["name"]))))

  if( length(names) > length( unique(names)) ){
    stop("Identical node names found. Please ensure the dagitty object contains unique names for each node.")
  }

  names_list <- strsplit(names, split = "_")

  if( label_type == "short" ){

    initials_list <- unlist( get_label_helper(names_list) )

  }else if( label_type == "initials" ){

    initials_list <- lapply( names_list, function(x) paste( substr(x, 1, 1), collapse = "" ) )

  }else{

    stop("Invalid label_type input. Please use 'initials' or the default 'name'.")

  }

  if( length(unique(unlist(names))) > length(unique(unlist(initials_list))) ){

    if( label_type == "initials" ){

      replace_initials_list <- unlist( get_label_helper(names_list) )
    }

    new_initials <- lapply(1:length(initials_list), function(x) {

      new_initials <- lapply(1:length(initials_list), function(y) {

        if( initials_list[[x]] == initials_list[[y]] & x != y ){

          if( label_type == "initials" ){

            initials_list[x] <- replace_initials_list[x]

          }else{

            initials_list[x] <- names[x]
          }

        }else{

          initials_list[x] <- NA
        }

        Filter(Negate(is.na), initials_list[[x]])

        unlist(initials_list[[x]])
      })

      new_initials <- Filter(Negate(is.na), unlist(new_initials, recursive = FALSE))
    })

    initials_list <- Map(c, new_initials, initials_list)
    initials_list <- lapply(initials_list, "[[", 1)
  }

  labels <- unlist(initials_list)

  names(labels) <- names

  return(labels)
}


#' Compare dagitty object edges
#'
#' @importFrom data.table setDT
#' @param dag1 A dagitty object.
#' @param dag2 A dagitty object.
#' @param compare_all Return mutually exclusive edges in both dag1 and dag2 (compare_all = TRUE) e.g. dag1 edges not in dag2 AND dag2 edges not in dag1, or only dag1 mutually exclusive edges (compare_all = FALSE), e.g. dag1 edges not in dag2. Defaults to compare_all = TRUE.
#' @param include_roles Determines whether to include node roles in addition to node names and edge direction. Defaults to include_roles = FALSE, e.g. only considers edges where descendant/ancestor node names and edge direction are different.
#' @return Data table of edges.
#' @examples
#' get_diff_edges(dag)
#'
#' @export
get_diff_edges <- function(dag1, dag2, compare_all = TRUE,  include_roles = FALSE){
  .datatable.aware <- TRUE

  if( include_roles == FALSE){
    edges <- get_edges(dag1)[,c("ancestor", "edge", "descendant")]
    edges2 <- get_edges(dag2)[,c("ancestor", "edge", "descendant")]
  }else{
    edges <- get_edges(dag1)
    edges2 <- get_edges(dag2)
  }

  col_names <- names(edges)

  edges_comb <- data.table::setDT(edges)[!edges2, on = col_names]

  if(compare_all == TRUE){
    edges2 <- data.table::setDT(edges2)[!edges, on = col_names]

    edges_comb <- merge(edges_comb, edges2,
                   by = col_names,
                   all = TRUE)
  }

  return(edges_comb)
}


#' Compare dagitty object node roles
#'
#' @param dag1 A dagitty object.
#' @param dag2 A dagitty object.
#' @param compare_all Compare and return mutually all exclusive edges  (compare_all = TRUE) e.g. dag1 edges not in dag2 AND dag2 edges not in dag1, or set compare_all = FALSE to return only dag1 mutually exclusive edges, e.g. dag1 edges not in dag2. Defaults to compare_all = TRUE.
#' @param nodes_in_common Determines whether to include only the nodes found in both DAGs. Defaults to nodes_in_common = FALSE (returns all nodes).
#' @return A list of nodes.
#' @examples
#' get_diff_roles(dag)
#'
#' @export
get_diff_roles <- function(dag1,
                           dag2,
                           compare_all = TRUE,
                           nodes_in_common = FALSE
                           ){
  .datatable.aware <- TRUE

  roles <- get_roles(dag1)
  roles_y <- get_roles(dag2)

  diff_list <- list()

  diff_list <- lapply(1:length(roles), function(x){

    diff_list[[x]] <- roles[[x]][ !roles[[x]] %in% roles_y[[x]] ]

    if( compare_all == TRUE){
      diff_list[[x]] <- c( diff_list[[x]], roles_y[[x]][ !roles_y[[x]] %in% roles[[x]] ] )
    }
    if( nodes_in_common == TRUE){
      dag_names <- names(dag1)
      dag2_names <- names(dag2)
      diff_list[[x]] <- diff_list[[x]][ diff_list[[x]] %in% dag_names &
                                          diff_list[[x]] %in% dag2_names ]
    }
    diff_list[[x]]
  })

  names(diff_list) <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "competing_exposure", "collider", "latent", "observed")

  diff_list <- Filter(function(x) length(x) > 0, diff_list)

  return(diff_list)
}
