

#' Mediator names in a dagitty object
#'
#' mediators() is a dagitty wrapper that identifies nodes along paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of mediators, or dataframe of edges.
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

    edges <- getEdges(dag, "mediator")

    return(edges)
  }

  return(mediators)
}


#' Confounder names in a dagitty object
#'
#' confounders() is a dagitty wrapper that identifies common causes of treatment and outcome in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns  A vector of confounder names, or edges in a data table.
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

    edges <- getEdges(dag, "confounder")

    return(edges)
  }

  return(confounders)
}

#' Mediator-outcome confounder names in a dagitty object
#'
#' moc() is a dagitty wrapper that identifies nodes along paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of mediator-outcome confounder names, or edges in a data table.
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

    edges <- getEdges(dag, "mediator-outcome-confounder")

    return(edges)
  }

  return(moc)
}


#' Variable names in a dagitty object
#'
#' variables() returns a named vector of variables corresponding to nodes in a dag, providing an easy way to extract node roles e.g. "treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", or "instrumental".
#'
#' @param dag A dagitty object.
#' @param all_info Defaults to TRUE. Set all_info = FALSE to return only variable names.
#' @returns A data table of variable names and their roles, or a vector of variable names.
#' @export
variables <- function(dag, all_info = TRUE){

  variables_df <- unique(getEdges(dag)[,c("ancestor", "role_ancestor")])

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
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of competing exposure names, or edges in a data table.
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

    edges <- getEdges(dag, "competing_exposure")

    return(edges)
  }

  return(competing_exposure)
}


#' Proxy node names in a dagitty object
#'
#' proxies() is a dagitty wrapper that identifies nodes that are proxy variables for indicated unmeasured confounding (existing along a path between a latent variable and outcome) in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @returns A vector of proxy variable names, or edges in a data table.
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
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of competing exposure names, or edges in a data table.
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

    edges <- getEdges(dag, "collider")

    return(edges)
  }

  return(colliders)
}

#' Extracts instrumental variable names from a dagitty object
#'
#' getInstrumentalVariables() is a dagitty wrapper function capable of identifying  multiple instrumental variables in multi-treatment and multi-outcome settings.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns Vector of instrumental variable names.
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


