
#' Treatment variables in a dagitty object
#'
#' treatments() returns a vector of node names assigned "\code{[exposure]}" in a dagitty object (similar to dagitty::exposures())
#'
#' @param dag A dagitty object.
#' @returns Character vector of treatment names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' treatments(dag)
#'
#' @export
treatments <- function(dag){

  return(.get_dagitty_role(dag, "treatment"))
}


#' Unobserved (latent) variables in a dagitty object
#'
#' unobserved() returns a vector of node names assigned "\code{[latent]}" in a dagitty object (similar to dagitty::latents())
#'
#' @param dag A dagitty object.
#' @returns Character vector of unobserved (latent) names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' unobserved(dag)
#'
#' @export
unobserved <- function(dag){

  return(.get_dagitty_role(dag, "unobserved"))
}


#' Observed variables in a dagitty object
#'
#' observed() returns a vector of variables that are not assigned "\code{[latent]}" in a dagitty object.
#' This is the measurement status; the \code{observed} element of \code{get_roles()} instead lists nodes without another role.
#'
#' @param dag A dagitty object.
#' @returns Character vector of observed names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' observed(dag)
#'
#' @export
observed <- function(dag){

  return(setdiff(.get_node_names(dag), unobserved(dag)))
}


#' Outcome variables in a dagitty object
#'
#' outcomes() returns a vector of node names assigned "\code{[outcome]}" in a dagitty object (similar to dagitty::outcomes()).
#' Attaching causaliflower after dagitty masks \code{dagitty::outcomes()}; both read the same declarations, and the
#' replacement form \code{outcomes(dag) <- value} continues to work.
#'
#' @param dag A dagitty object.
#' @returns Character vector of outcome names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' outcomes(dag)
#'
#' @export
outcomes <- function(dag){

  return(.get_dagitty_role(dag, "outcome"))
}

#' @rdname outcomes
#' @param x A dagitty object.
#' @param value Character vector of outcome names to declare.
#' @export
`outcomes<-` <- function(x, value){

  dagitty::outcomes(x) <- value
  return(x)
}

#' Outcome variables in a dagitty object (internal)
#'
#' Internal function for outcomes - returns a vector of node names assigned "\code{[outcome]}" in a dagitty object (similar to dagitty::outcomes())
#'
#' @param dag A dagitty object.
#' @returns Character vector of outcome names; character(0) if none.
#' @noRd
.outcomes <- function(dag){

  return(.get_dagitty_role(dag, "outcome"))
}


#' Node roles in dagitty
#'
#' Returns node names assigned a given dagitty role, read directly from the dagitty model string
#' (see .dag_to_string() / .lowercase_nodes()) instead of via dagitty functions.
#'
#' @param dag A dagitty object.
#' @param role One of "treatment", "outcome", "unobserved".
#' @returns Character vector of node names, or character(0) if none.
#' @noRd
.get_dagitty_role <- function(dag, role){

  tags <- switch(role,
                 treatment   = c("exposure", "e"),
                 outcome    = c("outcome", "o"),
                 unobserved = c("latent", "unobserved", "u"),
                 stop("Unknown role: ", role, ". Expected treatment, outcome, or unobserved."))
  lowercase_node_names <- .lowercase_nodes(.dag_to_string(dag))
  matched_node_names   <- names(lowercase_node_names)[vapply(lowercase_node_names, function(a) any(tags %in% a), logical(1L))]

  if( length(matched_node_names) > 0L ){
    return(as.character(matched_node_names))
  }

  return(character(0L))
}


#' All variable names in a dagitty object (including isolated nodes)
#'
#' Returns every node name in a dagitty object, including isolated (edge-less) nodes. Used by observed() and the get_* family.
#'
#' @param dag A dagitty object.
#' @returns Character vector of node names, or character(0).
#' @noRd
.get_node_names <- function(dag){

  s  <- .dag_to_string(dag)
  s  <- gsub("\\[[^][]*\\]", " ", s)                            # drop [attribute] blocks (and pos="x,y")
  qn <- regmatches(s, gregexpr("\"[^\"]*\"", s, perl = TRUE))[[1]]   # quoted names (now pos-free, may hold spaces)
  qn <- vapply(qn, .noquote, character(1L))
  s  <- gsub("\"[^\"]*\"", " ", s)                              # remove quoted names from the string
  s  <- sub("^\\s*(dag|pdag|mag|pag)\\b", " ", s, perl = TRUE)  # drop the graph-type keyword
  s  <- gsub("[{}]", " ", s)                                    # drop braces
  s  <- gsub("<->|<-|->|--", " ", s)                            # drop edge operators -> leave bare endpoints
  bn <- strsplit(s, "[[:space:];]+", perl = TRUE)[[1]]          # bare node names
  bn <- bn[nzchar(bn)]
  vals <- unique(c(as.character(qn), as.character(bn)))

  if( length(vals) > 0L ){

    return(vals)
  }

  return(character(0L))
}


#' Obtain the dagitty model string (internal)
#'
#' @param dag A dagitty object.
#' @returns Length-1 character: the dagitty model string.
#' @noRd
.dag_to_string <- function(dag){

  cands <- list(tryCatch(as.character(dag), error = function(e) NULL),
                tryCatch(unclass(dag),      error = function(e) NULL),
                tryCatch(format(dag),       error = function(e) NULL))
  for( s in cands ){
    if( is.null(s) ) next
    s <- paste(as.character(s), collapse = "\n")
    if( grepl("\\{|->|--|<->|\\[", s) ){
      graph_type <- sub("^\\s*([[:alpha:]]+).*", "\\1", s)
      if( !graph_type %in% c("dag", "pdag", "mag") ){
        stop("causaliflower does not currently support dagitty graph type '",
             graph_type, "'.", call. = FALSE)
      }
      return(s)
    }
  }

  stop("Could not read a dagitty model string from the supplied object. ",
       "Please check that a dagitty object is supplied and try again.", call. = FALSE)
}

#' Strip one layer of surrounding double quotes (internal)
#'
#' @param x Character vector.
#' @returns x with a single pair of wrapping double quotes removed where present.
#' @noRd
.noquote <- function(x){

  return(sub('^"(.*)"$', "\\1", x))
}

#' Parse node-attribute declarations from a dagitty model string (internal)
#'
#' Returns a named list of lowercase node names. Helper to .get_dagitty_role().
#'
#' @param s character dagitty model string (length of 1)
#' @returns Named list
#' @noRd
.lowercase_nodes <- function(s){

  pat <- '("[^"]*"|[A-Za-z0-9_.]+)\\s*\\[([^][]*)\\]'
  m   <- gregexpr(pat, s, perl = TRUE)[[1]]
  if( m[1] == -1L ) return(list())
  starts <- m; lens <- attr(m, "match.length")
  out <- list()
  for( i in seq_along(starts) ){
    chunk <- substr(s, starts[i], starts[i] + lens[i] - 1L)
    nm   <- .noquote(sub(pat, "\\1", chunk, perl = TRUE))
    att  <- sub(pat, "\\2", chunk, perl = TRUE)
    att  <- gsub('[A-Za-z0-9_.]+\\s*=\\s*"[^"]*"', "", att)     # drop key="value" pairs, e.g. pos="x,y"
    toks <- trimws(strsplit(att, "[,[:space:]]+", perl = TRUE)[[1]]) # tags are case-sensitive, matching dagitty's grammar
    out[[nm]] <- unique( c( out[[nm]], toks[nzchar(toks)] ) ) # merge tokens across repeated declarations of one node
  }

  return(out)
}


#' Mediator names in a dagitty object
#'
#' mediators() identifies nodes along directed paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom dagitty parents
#' @param dag A dagitty object.
#' @returns Character vector of mediator names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' mediators(dag)
#'
#' @export
mediators <- function(dag){

  mediators <- .get_roles(dag)$mediator

  return(as.character(mediators[ complete.cases(mediators) ]))
}


#' Confounder names in a dagitty object
#'
#' confounders() identifies common causes of treatment and outcome in a directed acyclic graph. Ancestral common causes are included, as are proxy variables of indicated unmeasured confounding.
#'
#' @importFrom dagitty parents
#' @param dag A dagitty object.
#' @returns Character vector of confounder names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' confounders(dag)
#'
#' @export
confounders <- function(dag){

  confounders <- .get_roles(dag)$confounder

  return(as.character(confounders[ complete.cases(confounders) ]))
}

#' Mediator-outcome confounder names in a dagitty object
#'
#' mediator_outcome_confounders() identifies nodes that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y) in a directed acyclic graph.
#'
#' @importFrom dagitty parents
#' @param dag A dagitty object.
#' @returns Character vector of mediator-outcome confounder names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' mediator_outcome_confounders(dag)
#'
#' @export
mediator_outcome_confounders <- function(dag){

  mediator_outcome_confounders <- .get_roles(dag)$mediator_outcome_confounder

  return(as.character(mediator_outcome_confounders[ complete.cases(mediator_outcome_confounders) ]))
}


#' Competing cause node names in a dagitty object
#'
#' competing_causes() returns direct causes of the outcome that are not otherwise classified -- i.e. causes of the outcome lying off the treatment->outcome path and unconnected to the treatment. Sometimes called a competing exposure (Tennant et al., 2021).
#'
#' @importFrom dagitty edges coordinates dagitty
#' @param dag A dagitty object.
#' @returns Character vector of competing cause names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' competing_causes(dag)
#'
#' @export
competing_causes <- function(dag){

  competing_causes <- .get_roles(dag)$competing_cause

  return(as.character(competing_causes[ complete.cases(competing_causes) ]))
}

#' Renamed: competing_exposures
#'
#' competing_exposures() was renamed competing_causes(). Calling this name stops
#' with a pointer to the new one; it will be removed in a later release.
#'
#' @param dag A dagitty object.
#' @returns No return value; calling this function stops with a renaming message.
#' @export
competing_exposures <- function(dag){

  stop("competing_exposures() was renamed competing_causes(); call competing_causes() instead.",
       call. = FALSE)
}


#' Proxy node names in a dagitty object
#'
#' proxies() is a dagitty wrapper that identifies nodes that are proxy variables for indicated unmeasured confounding (direct children of a latent variable that are direct parents of treatment or outcome) in a directed acyclic graph.
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @returns Character vector of proxy variable names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' proxies(dag)
#'
#' @export
proxies <- function(dag){

  # treatment
  treatments <- treatments(dag)

  # outcome
  outcomes <- .outcomes(dag)

  # latent variables
  latent_vars <- unobserved(dag)

  edges <- dagitty::edges(dag)

  treatment_parents <- .parents_edge_order(edges, treatments)

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- .parents_edge_order(edges, outcomes)
  treatment_children <- .children_edge_order(edges, treatments)

  nodes_trt_to_y <- nodes_between_treatment_and_outcome(dag = dag,
                                                            treatments = treatments,
                                                            outcomes = outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]

  # proxy
  latent_children <- .children_edge_order(edges, latent_vars)

  proxy_b <- suppressWarnings( treatment_parents[ treatment_parents %in% latent_children &
                                                    !treatment_parents %in% latent_vars &
                                                    !treatment_parents %in% treatments ] )# proxy_b

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  proxy <- unique( c(proxy_b, proxy_c) )

  return(proxy)
}


#' Collider names in a dagitty object
#'
#' colliders() identifies collider nodes that are common children of both treatment and outcome (treatment -> C <- outcome) in a directed acyclic graph; it does not return arbitrary two-parent colliders.
#'
#' @importFrom dagitty children
#' @param dag A dagitty object.
#' @returns Character vector of collider names; character(0) if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' colliders(dag)
#'
#' @export
colliders <- function(dag){

  colliders <- .get_roles(dag)$collider

  return(as.character(colliders[ complete.cases(colliders) ]))
}


#' Get instrumental variable names
#'
#' instruments() is a dagitty wrapper function capable of identifying  multiple instrumental variables in multi-treatment and multi-outcome settings.
#' Instruments are identified for each treatment-outcome pair using `dagitty::instrumentalVariables()`;
#'
#' @importFrom dagitty parents children descendants
#' @param dag A dagitty object.
#' @param usable_instruments When \code{TRUE}, exclude candidates that are
#'   descendants of their treatment. Defaults to \code{FALSE}.
#' @param details When \code{FALSE} (the default), return the nested-list shape
#'   used by earlier releases. When \code{TRUE}, return a data table containing
#'   treatment, outcome, instrument, and conditioning set.
#' @returns A nested list, or a data table when \code{details = TRUE}.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' instruments(dag)
#' instruments(dag, details = TRUE)
#'
#' @export
instruments <- function(dag,
                        usable_instruments = FALSE,
                        details = FALSE
                        ){

  treatments <- treatments(dag)
  outcomes <- .outcomes(dag)
  instrumental_vars <- instrumental_variables_helper(dag = dag,
                                                     treatments = treatments,
                                                     outcomes = outcomes,
                                                     usable_instruments = usable_instruments)
  if( isTRUE(details) ){
    return(instrumental_vars)
  }

  result <- lapply(treatments, function(treatment){
    per_outcome <- lapply(outcomes, function(outcome){
      values <- instrumental_vars$instrument[
        instrumental_vars$treatment == treatment &
          instrumental_vars$outcome == outcome
      ]
      if( length(values) == 0L ) NULL else unique(as.character(values))
    })
    names(per_outcome) <- outcomes
    per_outcome
  })
  names(result) <- treatments
  result
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
#' @returns Named list of minimally sufficient adjustment sets.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
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

  if( !is.numeric(num_sets) || length(num_sets) != 1L || !is.finite(num_sets) ||
      num_sets < 1 || num_sets != floor(num_sets) ){
    stop("num_sets must be a positive integer of at least 1. \n Please increase num_sets and try again.")
  }
  if( !effect %in% c("total", "direct") ){
    stop("effect must be either \"total\" or \"direct\".", call. = FALSE)
  }
  adjustment_sets <- tryCatch(
    dagitty::adjustmentSets(x = dag,
                            exposure = treatment,
                            outcome = outcome,
                            type = "minimal",
                            effect = effect),
    error = function(e){
      # An error here means the dag argument could not be processed
      stop("Failed to compute adjustment sets: ", conditionMessage(e),
           "\n Please check that a valid dagitty object is supplied and try again.",
           call. = FALSE)
    })

  if( length(adjustment_sets) == 0){

    # Indicates why no set is available (e.g., cycles, latent variables)
    has_latents <- tryCatch(length(unobserved(dag)) > 0, error = function(e) FALSE)
    is_cyclic   <- FALSE
    if( "isAcyclic" %in% getNamespaceExports("dagitty") ){
      is_cyclic <- isFALSE(tryCatch(dagitty::isAcyclic(dag), error = function(e) NA))
    }
    why <- character(0)
    if( isTRUE(is_cyclic) ){ why <- c(why, "the graph is not acyclic (contains a cycle)") }
    if( isTRUE(has_latents) ){ why <- c(why, "latent variables are present") }
    reason <- if( length(why) > 0 ){ paste0("\n This can occur because ", paste(why, collapse = " and "), ".") }else{ "" }
    message("No available adjustment sets for the ", effect, " effect.", reason,
            "\n Try adjusting parameters, or assess edges using assess_edges(dag, assess_causal_criteria = TRUE).")

    return(invisible())
  }
  adjustment_sets_list <- lapply(seq_along(adjustment_sets), function(x){
    adjustment_sets <- adjustment_sets[[x]]
  })

  set_length <- sapply(adjustment_sets_list, length)
  names(adjustment_sets_list) <- set_length
  adjustment_sets_list <- adjustment_sets_list[order(set_length)]
  num_available <- length(adjustment_sets_list)

  # Warning produced when a set is excluded due to a same-length tie
  if( num_available > num_sets &&
      length( adjustment_sets_list[[num_sets]] ) == length( adjustment_sets_list[[num_sets + 1]] ) ){
    warning("A same length adjustment set has been excluded. If this is unintentional, please adjust the num_sets parameter input and run again.")
  }
  adjustment_sets <- adjustment_sets_list[ 1:min(num_sets, num_available) ]

  return(adjustment_sets)
}
