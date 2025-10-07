
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
#' @returns A vector of confounders, or dataframe of edges.
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
#' @returns A vector of mediator-outcome confounders, or dataframe of edges.
#' @export
moc <- function(dag, edges = FALSE){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # moc
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, treatments)]

  # mediator-outcome confounders
  mediators_parents <- dagitty::parents(dag, mediators)
  moc <- mediators_parents[!mediators_parents %in% c(mediators, treatments)]


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

  variables_df <- unique(getEdges(dag)[c("v", "role_v", "latent_v")])

  if(all_info == TRUE){

    variables <- dplyr::rename(variables_df, variable_names = "v", variable_role = "role_v", latent_indicator = "latent_v")

    return(variables)

  }

  variables <- variables_df$v

  names(variables) <- gsub("[0-9]", "", variables_df$role_v)



  return(variables)
}
