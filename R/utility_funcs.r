
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
mediators <- function(dag){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]


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
confounders <- function(dag){

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
mediator_outcome_confounders <- function(dag){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) ]

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% confounders & # remove treatment and confounder nodes
                !moc %in% treatment_parents ] # double check by removing parents of treatment


  return(moc)
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
competing_exposures <- function(dag){

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) ]

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]

  # competing exposure
  competing_exposure <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% outcomes ]

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
proxies <- function(dag){
  .datatable.aware <- TRUE
  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) ]

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]

  # proxy

  latent_children <- dagitty::children(dag, latent_vars)

  proxy_b <- suppressWarnings( treatment_parents[ treatment_parents %in% latent_children &
                                                    !treatment_parents %in% latent_vars &
                                                    !treatment_parents %in% treatments ] )# proxy_b

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  proxy <- c(proxy_b, proxy_c)

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
colliders <- function(dag){
  # treatment
  treatments <- dagitty::exposures(dag)

  treatment_children <- dagitty::children(dag, treatments)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  outcome_children <- dagitty::children(dag, outcomes)

  # collider
  colliders <- outcome_children[ outcome_children %in% treatment_children ]


  return(colliders)
}


#' Extracts instrumental variable names from a dagitty object
#'
#' getInstrumentalVariables() is a dagitty wrapper function capable of identifying  multiple instrumental variables in multi-treatment and multi-outcome settings.
#'
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

  # treatment
  treatments <- dagitty::exposures(dag)

  treatment_children <- dagitty::children(dag, treatments)

  treatment_parents <- dagitty::parents(dag, treatments)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  outcome_children <- dagitty::children(dag, outcomes)

  outcome_parents <- dagitty::parents(dag, outcomes)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # confounders
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) ]

  # mediators - first parse (includes mediator-outcome confounders)
  nodes_trt_to_y <- get_nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables

  latent_children <- dagitty::children(dag, latent_vars)

  # collider
  colliders <- outcome_children[ outcome_children %in% treatment_children ]

  # instrumental variables
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


  instrumental_vars <- instrumental_vars[,1][instrumental_vars$role ==  "instrumental"]

  return(instrumental_vars)

}


#' Minimal sufficient adjustment sets
#'
#' minimal_sets() is a dagitty::adjustmentSets() wrapper for obtaining minimal adjustment sets, returning the smallest 5 (if available) by default.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty adjustmentSets
#' @param dag A dagitty object.
#' @param treatment Vector of treatment(s).
#' @param outcome Vector of outcome(s).
#' @param effect Defaults to total effect, options available as per dagitty::adjustmentSets()
#' @param num_sets Number of sets to return, defaults to the smallest five minimally sufficient adjustment sets.
#' @param decreasing Defaults to FALSE (shows minimally sufficient). Optionally can be set to filter the largest minimally sufficient adjustment sets.
#' @returns Named list of minimally sufficient adjustment sets.
#' @examples
#' minimal_sets(dag) # defaults to the total effect and 5 smallest sets
#'
#' minimal_sets(dag, effect = "direct") # direct effect
#'
#' minimal_sets(dag, effect = "direct", num_sets = 1) # return only the smallest set (direct effect)
#'
#' @export
minimal_sets <- function(dag,
                         treatment = NULL,
                         outcome = NULL,
                         effect = "total",
                         num_sets = 5
                         ){
  .datatable.aware <- TRUE
  adjustment_sets <- dagitty::adjustmentSets(dag,
                                             exposure = treatment,
                                             outcome = outcome,
                                             type = "minimal",
                                             effect = effect)
  if( length(adjustment_sets) == 0){

    message("No available adjustment sets. Try adjusting parameters, or assess edges using assess_edges(dag, assess_causal_criteria = TRUE).")

    return(invisible())

  }

  adjustment_sets_list <- lapply(1:length(adjustment_sets), function(x){

    adjustment_sets <- adjustment_sets[[x]]

  })

  set_length <- sapply(adjustment_sets_list, length)
  names(adjustment_sets_list) <- set_length

  adjustment_sets_list <- adjustment_sets_list[order(set_length)]

  adjustment_sets <- adjustment_sets_list[1:(num_sets+1)]

  if( length( adjustment_sets[[1]] ) == length( adjustment_sets[[num_sets]] ) &
      length( adjustment_sets[[num_sets]] ) == length(unlist(adjustment_sets[num_sets+1])) ){

    warning("\nA minimally sufficient adjustment set has been excluded from the minimal_sets() output.
            \nIf this is accidental, it is highly recommended to adjust the num_sets parameter input and run again.")

  }else if(  length( adjustment_sets[[num_sets]] ) > 0 & length( adjustment_sets[[num_sets]] ) == length(unlist(adjustment_sets[num_sets+1])) ){

    warning("Same length adjustment set has been not outputted.")

  }

  adjustment_sets <- adjustment_sets_list[1:(num_sets)]

  adjustment_sets <- Filter(Negate(is.null), adjustment_sets)

  return(adjustment_sets)

}
