

#' Mediator names in a dagitty object
#'
#' mediators() is a dagitty wrapper that identifies nodes along paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of mediators, or dataframe of edges.
#' @export
mediators <- function(dag, edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # mediators
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, treatments)]

  if(edges == TRUE){
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
confounders <- function(dag, edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes)]

  if(edges == TRUE){
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
mediator_outcome_confounders <- function(dag, edges = FALSE){

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


  if(edges == TRUE){
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
#' @returns A named vector of nodes in a dag.
#' @export
variables <- function(dag, all_info = TRUE){

  variables_df <- unique(getEdges(dag)[c("ancestor", "role_ancestor")])

  if(all_info == TRUE){

    variables <- dplyr::rename(variables_df, variable_names = "ancestor", variable_role = "role_ancestor")

    return(variables)

  }

  variables <- variables_df$ancestor

  names(variables) <- gsub("[0-9]", "", variables_df$role_ancestor)



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
competing_exposures <- function(dag, edges = FALSE){

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

  if(edges == TRUE){

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
#' @param dag A dagitty object.
#' @returns A vector of proxy variable names, or edges in a data table.
#' @export
proxies <- function(dag, edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)

  # latent children
  latent_descendants <- dagitty::children(dag, latent_vars)

  proxy_b <- treatment_parents[ treatment_parents %in% latent_descendants] # proxy_b
  proxy_c <- outcome_parents[ outcome_parents %in% latent_descendants & # proxy_c
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators ]
  proxy <- c(proxy_b, proxy_c)

  if(edges == TRUE){

    edges <- getEdges(dag, "proxy")

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
colliders <- function(dag, edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)
  treatment_children <- dagitty::children(dag, treatments)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # collider
  outcome_children <- dagitty::children(dag, outcomes)
  colliders <- outcome_children[outcome_children %in% treatment_children]


  if(edges == TRUE){

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

  return(instrumental_vars)

}


