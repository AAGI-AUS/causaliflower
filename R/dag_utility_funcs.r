

#' Mediator names in a dagitty object
#'
#' mediators() is a dagitty wrapper that identifies nodes along paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom dagitty exposures outcomes parents
#' @param dag A dagitty object.
#' @returns A vector of mediators, or dataframe of edges.
#' @examples
#' mediators(dag)
#'
#' @export
mediators <- function(dag, get_edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # mediators
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, treatments)]

  if(get_edges == TRUE){
    cat(paste("There are", length(mediators), "mediators in the supplied graph: ", sep = " ", collapse = " "))
    cat(paste("\n", mediators, collapse = "\n"))

    edges <- get_edges(dag, "mediator")

    return(edges)
  }

  return(mediators)
}


#' Confounder names in a dagitty object
#'
#' confounders() is a dagitty wrapper that identifies common causes of treatment and outcome in a directed acyclic graph.
#'
#' @importFrom dagitty exposures outcomes parents
#' @param dag A dagitty object.
#' @returns  A vector of confounder names, or edges in a data table.
#' @examples
#' confounders(dag)
#'
#' @export
confounders <- function(dag, get_edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes)]

  if(get_edges == TRUE){
    cat(paste("There are", length(confounders), "confounders in the supplied graph: ", sep = " ", collapse = " "))
    cat(paste("\n", confounders, collapse = "\n"))

    edges <- get_edges(dag, "confounder")

    return(edges)
  }

  return(confounders)
}

#' Mediator-outcome confounder names in a dagitty object
#'
#' moc() is a dagitty wrapper that identifies nodes along paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom dagitty exposures outcomes latents parents
#' @param dag A dagitty object.
#' @returns A vector of mediator-outcome confounder names, or edges in a data table.
#' @examples
#' mediator_outcome_confounders(dag)
#'
#' @export
mediator_outcome_confounders <- function(dag, get_edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # mediators
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, treatments)]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators)
  moc <- mediator_parents[!mediator_parents %in% c(mediators, treatments)]


  if(get_edges == TRUE){
    cat(paste("There are", length(moc), "mediator-outcome confounders in the supplied graph: ", sep = " ", collapse = " "))
    cat(paste("\n", moc, collapse = "\n"))

    edges <- get_edges(dag, "mediator-outcome-confounder")

    return(edges)
  }

  return(moc)
}


#' Variable names in a dagitty object
#'
#' variables() returns a named vector of variables corresponding to nodes in a dag, providing an easy way to extract node roles e.g. "treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", or "instrumental".
#'
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty
#' @param dag A dagitty object.
#' @param all_info Defaults to TRUE. Set all_info = FALSE to return only variable names.
#' @returns A data table of variable names and their roles, or a vector of variable names.
#' @examples
#' variables(dag)
#'
#' @export
variables <- function(dag, all_info = TRUE){

  variables_df <- unique(get_edges(dag)[,c("ancestor", "role_ancestor")])

  if(all_info == TRUE){

    colnames(variables_df) <- c("name", "role")

    return(variables_df)

  }

  variables <- variables_df$ancestor

  return(variables)
}


#' Competing exposure node names in a dagitty object
#'
#' competing_exposure() is a dagitty wrapper that identifies nodes in a directed acyclic graph connected to outcome, other than indicated exposures.
#'
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty
#' @param dag A dagitty object.
#' @returns A vector of competing exposure names, or edges in a data table.
#' @examples
#' competing_exposures(dag)
#'
#' @export
competing_exposures <- function(dag, get_edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent
  latent_vars <- dagitty::latents(dag)

  # mediators
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, treatments)]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators)
  moc <- mediator_parents[!mediator_parents %in% c(mediators, treatments)]


  # confounder
  confounders <- confounders(dag)

  # competing exposure
  competing_exposure <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% latent_vars &
                                           !outcome_parents %in% moc ]

  if(get_edges == TRUE){

    edges <- get_edges(dag, "competing_exposure")

    return(edges)
  }

  return(competing_exposure)
}


#' Proxy node names in a dagitty object
#'
#' proxies() is a dagitty wrapper that identifies nodes that are proxy variables for indicated unmeasured confounding (existing along a path between a latent variable and outcome) in a directed acyclic graph.
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @returns A vector of proxy variable names, or edges in a data table.
#' @examples
#' proxies(dag)
#'
#' @export
proxies <- function(dag, get_edges = FALSE){
  .datatable.aware <- TRUE
  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) &
                                     !treatment_parents %in% treatments &
                                     !treatment_parents %in% latent_vars]

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  mediators <- outcome_parents[outcome_parents %in% treatment_children &
                                 !outcome_parents %in% treatments &
                                 !outcome_parents %in% confounders &
                                 !outcome_parents %in% latent_vars]

  # latent children
  latent_children <- dagitty::children(dag, latent_vars)

  latent_children <- dagitty::children(dag, latent_vars)
  proxy_b <- treatment_parents[ treatment_parents %in% latent_children &
                                  !treatment_parents %in% latent_vars &
                                  !treatment_parents %in% treatments ] # proxy_b

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars ]

  proxy <- unique( c(proxy_b, proxy_c) )



  if(get_edges == TRUE){

    edges <- data.table::as.data.table(dagitty::edges(dag))[,1:3]

    edges <- edges[, c("v", "e", "w")]

    edges <- edges[ unlist(edges[,"v"]) %in% proxy, ]

    return(edges)
  }

  return(proxy)
}


#' Collider names in a dagitty object
#'
#' colliders() is a dagitty wrapper that identifies nodes in a directed acyclic graph connected to outcome, other than indicated exposures.
#'
#' @importFrom dagitty exposures outcomes children
#' @param dag A dagitty object.
#' @returns A vector of competing exposure names, or edges in a data table.
#' @examples
#' colliders(dag)
#'
#' @export
colliders <- function(dag, get_edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)
  treatment_children <- dagitty::children(dag, treatments)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # collider
  outcome_children <- dagitty::children(dag, outcomes)
  colliders <- outcome_children[outcome_children %in% treatment_children]


  if(get_edges == TRUE){

    edges <- get_edges(dag, "collider")

    return(edges)
  }

  return(colliders)
}

#' Extracts instrumental variable names from a dagitty object
#'
#' getInstrumentalVariables() is a dagitty wrapper function capable of identifying  multiple instrumental variables in multi-treatment and multi-outcome settings.
#'
#' @importFrom dagitty exposures outcomes latents parents children
#' @param dag A dagitty object.
#' @returns Vector of instrumental variable names.
#' @examples
#' instrumental_variables(dag)
#'
#' @export
instrumental_variables <- function(dag){


  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # treatment_parents
  treatment_parents <- dagitty::parents(dag, treatments)

  # outcome_parents
  outcome_parents <- dagitty::parents(dag, outcomes)

  #treatment_children
  treatment_children <- dagitty::children(dag, treatments)

  # mediator_parents
  mediators <- outcome_parents[outcome_parents %in% treatment_children]
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables

  # latent_children
  latent_children <- dagitty::children(dag, latent_vars)

  # outcome_children
  outcome_children <- dagitty::children(dag, outcomes)

  # collider
  colliders <- outcome_children[outcome_children %in% treatment_children]

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
  instrumental_vars <- unlist(instrumental_vars[,1])

  return(instrumental_vars)

}

#' Minimally sufficient adjustment sets
#'
#' minimal_sets() is a dagitty::adjustmentSets() wrapper for obtaining minimal adjustment sets, returning the smallest 5 (if available) by default.
#'
#' @importFrom dagitty adjustmentSets
#' @param dag A dagitty object.
#' @param treatment Vector of treatment(s).
#' @param outcome Vector of outcome(s).
#' @param effect Defaults to total effect, options available as per dagitty::adjustmentSets()
#' @param smallest Number of sets to return, defaults to the smallest five minimally sufficient adjustment sets.
#' @param decreasing Defaults to FALSE (showing smallest). Optionally can be set to filter the largest minimally sufficient adjustment sets.
#' @returns Named list of minimally sufficient adjustment sets.
#' @examples
#' minimal_sets(dag) # defaults to the total effect and 5 smallest sets
#'
#' minimal_sets(dag, effect = "direct") # direct effect
#'
#' minimal_sets(dag, effect = "direct", smallest = 1) # return only the smallest set (direct effect)
#'
#' @export
minimal_sets <- function(dag, treatment = NULL, outcome = NULL, effect = "total", smallest = 5, decreasing = FALSE){

  adjustment_sets <- dagitty::adjustmentSets(dag,
                                             exposure = treatment,
                                             outcome = outcome,
                                             type = "minimal",
                                             effect = effect)

  adj_set_length <- sapply(adjustment_sets, length)
  adjustment_sets <- adjustment_sets[ order(adj_set_length, decreasing = decreasing)]

  adjustment_sets_list <- lapply(1:length(adj_set_length), function(x){

    adjustment_sets[[x]] <- adjustment_sets[[x]]

  })

  names(adjustment_sets_list) <- adj_set_length

  adjustment_sets <- adjustment_sets_list[1:smallest]

  adjustment_sets <- Filter(Negate(is.null), adjustment_sets)

  return(adjustment_sets)

}


#' Remove edges using markov logic
#'
#' markov_graph()
#'
#' @importFrom dagitty exposures outcomes latents parents children
#' @param dag A dagitty object.
#' @returns A dagitty object.
#' @examples
#' markov_dag <- markov_graph(dag) # still in development, may produce some unexpected results
#'
#' @export
markov_graph <- function(dag){


  edges <- data.table::as.data.table(dagitty::edges(dag))[,1:3]
  edges <- edges[, c("v", "e", "w")]


  latent_vec <- dagitty::latents(dag)
  treatments <- dagitty::exposures(dag)
  outcomes <- dagitty::outcomes(dag)


  ## group by nodes
  unique_ancestors <- unique( edges[,"v"] )
  num_unique_ancestors <- nrow(unique_ancestors)


  edges_list <- list()

  # edges_to_assess grouped by each unique node in a list
  edges_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges[ unlist(edges[,"v"]) %in% unlist(unique_ancestors[x]), ]


  }) )


  ## unique children of each ancestor (used to remove edges)
  unique_ancestors_children <- lapply( 1:num_unique_ancestors, function(x){  # get new node name children

    unique_ancestors_children <- dagitty::children(dag, unique_ancestors[x])

  })


  ## remove edges list ##
  remove_edges_list <- list()

  remove_edges_list <- suppressWarnings( sapply( 1:num_unique_ancestors, function(x){

    unique_ancestors_children_vec <- unique_ancestors_children[[x]]

    edges_list[[x]] <- unique_ancestors_children_vec[
      unique_ancestors_children[[x]] %in% dagitty::children( dag, edges_list[[x]][[3]]) &
        !unique_ancestors_children[[x]] %in% treatments &
        !unique_ancestors_children[[x]] %in% outcomes &
        !unique_ancestors_children[[x]] %in% latent_vec
      ]

    } ) )

  ## new edges list ##
  new_edges_list <- list()
  new_edges_list <- suppressWarnings( sapply( 1:num_unique_ancestors, function(x){

    edges_unlisted <- unlist( edges_list[[x]][[3]] )
    new_edges_list[[x]] <- edges_unlisted[ !edges_unlisted %in% remove_edges_list[[x]] ]

  } ) )

  names(new_edges_list) <- unlist(unique_ancestors)

  ## edges output ##
  edges_out <- c()

  edges_out <- lapply( 1:num_unique_ancestors, function(x){

    lapply( 1:length( new_edges_list[[x]] ), function(n){

      edges_out <- data.table::data.table( v = names(new_edges_list[x]), e = "->", w = new_edges_list[[x]][[n]] )

    } )

  } )

  edges_out <- unlist(edges_out, recursive = FALSE)
  edges_out <- do.call( rbind, edges_out )

  dag_out <- rebuild_dag(dag, edges_out)

  return(dag_out)

}


#' dagitty node names and roles
#'
#' node_roles() is a dagitty wrapper that returns a list of node names and their roles.
#'
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @examples
#' node_roles(dag)
#'
#' @export
node_roles <- function(dag){
  .datatable.aware <- TRUE

  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  ## ancestor node edges to list ##
  edges_ancestors <- edges_wide[,1:14]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:14), idvar = "id",
                                      v.names = "role_ancestor", direction = "long")[,c("v", "e", "w", "role_ancestor", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), 1:4]
  names(edges_ancestors)[1:3] <- c("ancestor", "edge", "descendant")

  ## group by nodes
  unique_ancestors <- unique( edges_ancestors[,"ancestor"] ) # vector of unique node names
  num_unique_ancestors <- nrow(unique_ancestors) # count of unique node names

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

  unique_descendants <- unique(missing_descendant_edges[,1])
  num_unique_descendants <- nrow(unique_descendants)


  # descendant edges grouped by each unique node in a list
  edges_descendant_list <- c()

  edges_descendant_list <- suppressWarnings( lapply(1:num_unique_descendants, function(x){

    edges_descendant_list <- unique(unlist(
      edges_descendants[ unlist(edges_descendants[,"w"]) %in% unlist(unique_descendants[x]), ][, c(2)] ))


  }) )

  names(edges_descendant_list) <- unlist(unique_descendants) # assign descendant node role to name of each element in edges list


  edges_list <- c(edges_ancestor_list, edges_descendant_list) # combine ancestor and descendant node role lists


  return(edges_list)

}


#' dagitty node names, ancestor/descendant roles
#'
#' node_structure() is a dagitty wrapper function that returns a list of node names extrated from a dagitty object, including the name and role of their ancestor and descendant nodes.
#'
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @examples
#' node_structure(dag)
#'
#' @export
node_structure <- function(dag){
  .datatable.aware <- TRUE

  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  unique_nodes <- unname(unlist(  unique(edges_wide[,1]) ))
  num_unique_nodes <- length(unique_nodes)


  ## ancestor node edges to list ##
  edges_ancestors <- edges_wide[,1:14]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:14), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), c(1,3,4)]
  names(edges_ancestors)[1:2] <- c("ancestor", "descendant")

  edges_ancestors$role_node_id <- "ancestor"

  unique_ancestors <- unname(unlist(  unique(edges_ancestors[,1]) ))
  num_unique_ancestors <- length(unique_nodes)


  ## descendant node edges to list ##
  edges_descendants <- edges_wide[,c(1:3,15:25)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:14), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), c(1,3,4)]
  names(edges_descendants)[1:2] <- c("ancestor", "descendant")

  edges_descendants$role_node_id <- "descendant"

  unique_descendant <- unname(unlist(  unique(edges_descendants[,1]) ))
  num_unique_descendant <- length(unique_nodes)


  ## combine ancestor and descendant nodes ##
  edges_long <- rbind(edges_ancestors, edges_descendants)

  ancestor_descendant_string_vec <- c("ancestor", "descendant")
  ancestors_descendants_label_vec <- c("ancestors", "descendants")


  # ancestor edges grouped by each unique node in a list
  edges_ancestor_list <- c()

  edges_ancestor_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges_ancestor_list <- edges_long[ unlist(edges_long[,"ancestor"]) %in% unlist(unique_ancestors[x]) &
                  unlist(edges_long[,"role_node_id"]) == "descendant", ]

   #edges_ancestor_list[[x]]$ancestor_role <- unique(unlist( edges_long[,3][ unlist(edges_long[,1]) %in% unlist(unique_ancestors[x]) ] ))


  }) )


  # descendant edges grouped by each unique node in a list
  edges_descendant_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

    edges_long[ unlist(edges_long[,"descendant"]) %in% unlist(unique_ancestors[x]) &
                  unlist(edges_long[,"role_node_id"]) == "ancestor", ]


  }) )

  edges_list <- Map(rbind, edges_descendant_list, edges_ancestor_list)

  names(edges_list) <- unlist(unique_nodes)


  ## next create list with each node as an element
  edges_structure_list <- c()

  edges_structure_list <- lapply( 1:num_unique_ancestors , function(x){

    edges_structure_list <- lapply(1:2, function(n){ # edges grouped by unique node name


      edges_structure_list[[n]] <- edges_list[[x]][ edges_list[[x]][[4]] %in% unlist(ancestor_descendant_string_vec[n]), ]

    })

    names(edges_structure_list) <- ancestors_descendants_label_vec

    edges_structure_list

  } )

  names(edges_structure_list) <- unlist(unique_ancestors) # assign ancestor node role to name of each element in edges list


  return(edges_structure_list)

}


#' Long format node roles
#'
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

  edges <- data.table::as.data.table(dagitty::edges(dag))[,1:3]
  edges <- edges[, c("v", "e", "w")]

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


