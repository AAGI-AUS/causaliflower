# ==============================================================================
# utility_funcs.r


# Node status and role functions with wrappers for dagitty's identification checks
# ==============================================================================


#' Treatment variables in a dagitty object
#'
#' `treatments()` returns node names marked `[exposure]` in a dagitty object,
#' as does `dagitty::exposures()`.
#'
#' @param dag A dagitty object.
#' @returns A character vector of treatment names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' treatments(dag)
#'
#' @export
treatments <- function(dag){

  nodes <- get_dagitty_role(dag, "treatment")

  return(nodes)
}


#' Unobserved (latent) variables in a dagitty object
#'
#' `unobserved()` returns node names marked `[latent]` in a dagitty object,
#' as does `dagitty::latents()`.
#'
#' @param dag A dagitty object.
#' @returns A character vector of unobserved (latent) names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' unobserved(dag)
#'
#' @export
unobserved <- function(dag){

  nodes <- get_dagitty_role(dag, "unobserved")

  return(nodes)
}


#' Observed variables in a dagitty object
#'
#' `observed()` returns node names without a `[latent]` tag in a dagitty
#' object. This is their measurement status. The `observed` element of
#' `get_roles()` instead lists nodes without another role.
#'
#' @param dag A dagitty object.
#' @returns A character vector of observed names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' observed(dag)
#'
#' @export
observed <- function(dag){

  nodes <- setdiff(get_node_names(dag), unobserved(dag))

  return(nodes)
}


#' Outcome variables in a dagitty object
#'
#' `outcomes()` returns node names marked `[outcome]` in a dagitty object,
#' as does `dagitty::outcomes()`. Attaching causaliflower after dagitty masks
#' that function. The replacement form `outcomes(dag) <- value` still works.
#'
#' @param dag A dagitty object.
#' @returns `outcomes()` returns a character vector of outcome names, or
#'   `character(0)` if none. The replacement function returns the updated
#'   dagitty object.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' outcomes(dag)
#'
#' @export
outcomes <- function(dag){

  nodes <- get_dagitty_role(dag, "outcome")

  return(nodes)
}


#' @rdname outcomes
#' @param x A dagitty object.
#' @param value Character vector of outcome names to declare.
#' @export
`outcomes<-` <- function(x, value){

  dagitty::outcomes(x) <- value
  return(x)
}


#' Declare treatment variables in a dagitty object
#'
#' `set_treatment()` adds `[exposure]` to the named nodes and changes
#' nothing else: no edge, coordinate, curve, label or other node status is
#' altered.
#'
#' @param dag A dagitty object.
#' @param nodes Character vector of node names already in the graph.
#' @returns A dagitty object.
#' @examples
#' set_treatment(dagitty::dagitty("dag { T Y Z T -> Y }"), "T")
#'
#' @export
set_treatment <- function(dag, nodes){

  nodes <- check_status_nodes(dag, nodes, "set_treatment()")

  dag_out <- set_node_status(dag, nodes, "exposure", "set_treatment()")

  return(dag_out)
}


#' Declare outcome variables in a dagitty object
#'
#' `set_outcome()` adds `[outcome]` to the named nodes and changes
#' nothing else: no edge, coordinate, curve, label or other node status is
#' altered.
#'
#' @param dag A dagitty object.
#' @param nodes Character vector of node names already in the graph.
#' @returns A dagitty object.
#' @examples
#' set_outcome(dagitty::dagitty("dag { X Y Z X -> Y }"), "Y")
#'
#' @export
set_outcome <- function(dag, nodes){

  nodes <- check_status_nodes(dag, nodes, "set_outcome()")

  dag_out <- set_node_status(dag, nodes, "outcome", "set_outcome()")

  return(dag_out)
}


#' Declare unobserved (latent) variables in a dagitty object
#'
#' `set_latent()` adds `[latent]` to the named nodes and changes nothing
#' else: no edge, coordinate, curve, label or other node status is altered.
#'
#' @param dag A dagitty object.
#' @param nodes Character vector of node names already in the graph.
#' @returns A dagitty object.
#' @examples
#' set_latent(dagitty::dagitty("dag { X Y Z X -> Y }"), "Z")
#'
#' @export
set_latent <- function(dag, nodes){

  nodes <- check_status_nodes(dag, nodes, "set_latent()")

  dag_out <- set_node_status(dag, nodes, "latent", "set_latent()")

  return(dag_out)
}


#' Check node status setter inputs
#'
#' Checks that `dag` is a dagitty object and that `nodes` names existing nodes,
#' then returns those names.
#'
#' @param dag A dagitty object.
#' @param nodes Character vector of node names already in the graph.
#' @param caller Name of the calling function, used in the messages.
#' @returns Character vector of node names.
#' @noRd
check_status_nodes <- function(dag, nodes, caller){

  if( !inherits(dag, "dagitty") ){
    stop(caller, ": 'dag' must be a dagitty object.", call. = FALSE)
  }
  if( !is.character(nodes) || length(nodes) == 0L ){
    stop(caller, ": 'nodes' must name nodes in the graph: ",
         paste(deparse(nodes), collapse = " "), ".", call. = FALSE)
  }
  missing_nodes <- unique(nodes[!nodes %in% names(dag)])
  if( length(missing_nodes) > 0L ){
    stop(caller, ": 'nodes' must name nodes in the graph: ",
         paste(missing_nodes, collapse = ", "), ".", call. = FALSE)
  }

  return(nodes)
}


#' Write a node status into the native graph text
#'
#' Inserts the dagitty keyword into each node's declaration, keeping every
#' other byte of the graph text, then checks that only the requested status
#' changed.
#'
#' @param dag A dagitty object.
#' @param nodes Character vector of node names already in the graph.
#' @param keyword One of "exposure", "outcome", "latent".
#' @param caller Name of the calling function, used in the messages.
#' @returns A dagitty object.
#' @noRd
set_node_status <- function(dag, nodes, keyword, caller){

  status <- switch(keyword,
                   exposure = treatments(dag),
                   outcome  = get_outcome_names(dag),
                   latent   = unobserved(dag))
  pending <- setdiff(unique(nodes), status)
  if( length(pending) == 0L ){
    return(dag)
  }

  failure <- function(){
    stop(caller, ": the graph could not be updated without changing its ",
         "existing metadata.", call. = FALSE)
  }
  text <- as.character(dag)
  index <- tryCatch(native_layout_index(text, caller),
                    error = function(error) {
                      out <- NULL

                      return(out)
                    })
  if( is.null(index) ) {
    failure()
  }

  edits <- list()
  written <- character()
  for( block in index$blocks ){
    if( block$kind != "node" || is.null(block$node) ) {
      next
    }
    if( !block$node %in% pending || block$node %in% written ) {
      next
    }
    at <- index$tokens$end[block$first] + 1L
    edits[[length(edits) + 1L]] <- list(
      start = at, end = at - 1L,
      text = paste0(keyword, if( length(block$fields) > 0L ) { "," } else { "" })
    )
    written <- c(written, block$node)
  }
  for( declaration in index$declarations ){
    if( !declaration$node %in% pending || declaration$node %in% written ) {
      next
    }
    at <- index$tokens$end[declaration$token]
    edits[[length(edits) + 1L]] <- list(
      start = at + 1L, end = at, text = paste0(" [", keyword, "]")
    )
    written <- c(written, declaration$node)
  }
  undeclared <- setdiff(pending, written)
  if( length(undeclared) > 0L ){
    lines <- paste0(native_layout_quote(undeclared), " [", keyword, "]")
    edits[[length(edits) + 1L]] <- list(
      start = index$closing, end = index$closing - 1L,
      text = paste0("\n", paste(lines, collapse = "\n"), "\n")
    )
  }
  result <- native_layout_value(dag, native_layout_edits(text, edits))

  restored <- as.character(result)
  for( i in order(vapply(edits, `[[`, numeric(1L), "start")) ){
    inserted <- edits[[i]]$text
    at <- edits[[i]]$start
    if( !identical(substr(restored, at, at + nchar(inserted) - 1L), inserted) ){
      failure()
    }
    restored <- paste0(substr(restored, 1L, at - 1L),
                       substring(restored, at + nchar(inserted)))
  }
  if( !identical(restored, text) ) {
    failure()
  }
  if( !status_write_intact(dag, result, pending, keyword, caller) ) {
    failure()
  }

  return(result)
}


#' Check that a node status write changed nothing else
#'
#' Compares the node set, statuses, edges, coordinates and stored curve control
#' points of the rewritten graph with the supplied one.
#'
#' @param dag The supplied dagitty object.
#' @param result The rewritten dagitty object.
#' @param nodes Character vector of node names given the status.
#' @param keyword One of "exposure", "outcome", "latent".
#' @param caller Name of the calling function, used in the messages.
#' @returns TRUE if only the requested status changed, otherwise FALSE.
#' @noRd
status_write_intact <- function(dag, result, nodes, keyword, caller){

  same <- function(x, y){
    matches <- identical(sort(unique(as.character(x))), sort(unique(as.character(y))))

    return(matches)
  }
  keys <- function(x) {
    keys <- edge_control_keys(x$v, x$e, x$w)

    return(keys)
  }

  intact <- tryCatch({
    expected <- list(exposure = as.character(dagitty::exposures(dag)),
                     outcome  = as.character(dagitty::outcomes(dag)),
                     latent   = as.character(dagitty::latents(dag)))
    expected[[keyword]] <- union(expected[[keyword]], nodes)
    graph_nodes <- names(dag)
    before <- dagitty::coordinates(dag)
    after <- dagitty::coordinates(result)
    controls <- edge_controls(dag, caller)
    actual <- edge_controls(result, caller)
    found <- match(keys(controls), keys(actual))
    edges <- dag_edges(dag, caller)
    rewritten <- dag_edges(result, caller)

    same(dagitty::exposures(result), expected$exposure) &&
      same(dagitty::outcomes(result), expected$outcome) &&
      same(dagitty::latents(result), expected$latent) &&
      setequal(names(result), graph_nodes) &&
      length(names(result)) == length(graph_nodes) &&
      setequal(keys(rewritten), keys(edges)) &&
      length(keys(rewritten)) == length(keys(edges)) &&
      identical(after$x[graph_nodes], before$x[graph_nodes]) &&
      identical(after$y[graph_nodes], before$y[graph_nodes]) &&
      nrow(actual) == nrow(controls) && !anyNA(found) &&
      identical(as.numeric(actual$x[found]), as.numeric(controls$x)) &&
      identical(as.numeric(actual$y[found]), as.numeric(controls$y)) &&
      identical(dagitty::adjustedNodes(result), dagitty::adjustedNodes(dag)) &&
      identical(selected_nodes(result), selected_nodes(dag))
  }, error = function(error) {
    matches <- FALSE

    return(matches)
  })

  matches <- isTRUE(intact)

  return(matches)
}


#' Mediator names in a dagitty object
#'
#' `mediators()` identifies nodes along directed paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom dagitty parents
#' @param dag A dagitty object.
#' @returns A character vector of mediator names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' mediators(dag)
#'
#' @export
mediators <- function(dag){

  mediators <- classify_node_roles(dag)$mediator

  nodes <- as.character(mediators[ complete.cases(mediators) ])

  return(nodes)
}


#' Confounder names in a dagitty object
#'
#' `confounders()` identifies common causes of treatment and outcome in a directed acyclic graph. Ancestral common causes are included, as are proxy variables of indicated unmeasured confounding.
#'
#' @importFrom dagitty parents
#' @param dag A dagitty object.
#' @returns A character vector of confounder names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' confounders(dag)
#'
#' @export
confounders <- function(dag){

  confounders <- classify_node_roles(dag)$confounder

  nodes <- as.character(confounders[ complete.cases(confounders) ])

  return(nodes)
}


#' Mediator and outcome confounder names in a dagitty object
#'
#' `mediator_outcome_confounders()` identifies nodes that instead of being common causes of treatment and outcome (T <- X -> Y) are a common cause of mediators and outcome (M <- X -> Y) in a directed acyclic graph.
#'
#' @importFrom dagitty parents
#' @param dag A dagitty object.
#' @returns A character vector of mediator and outcome confounder names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' mediator_outcome_confounders(dag)
#'
#' @export
mediator_outcome_confounders <- function(dag){

  mediator_outcome_confounders <- classify_node_roles(dag)$mediator_outcome_confounder

  nodes <- as.character(mediator_outcome_confounders[ complete.cases(mediator_outcome_confounders) ])

  return(nodes)
}


#' Competing cause node names in a dagitty object
#'
#' `competing_causes()` returns direct causes of the outcome that are not otherwise classified -- i.e. causes of the outcome lying off the treatment->outcome path and unconnected to the treatment. Sometimes called a competing exposure (Tennant et al., 2021).
#'
#' @importFrom dagitty edges coordinates dagitty
#' @param dag A dagitty object.
#' @returns A character vector of competing cause names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' competing_causes(dag)
#'
#' @export
competing_causes <- function(dag){

  competing_causes <- classify_node_roles(dag)$competing_cause

  nodes <- as.character(competing_causes[ complete.cases(competing_causes) ])

  return(nodes)
}


#' Deprecated: competing_exposures
#'
#' `competing_exposures()` is deprecated as of 0.0.1. Use
#' `competing_causes()` instead. Calling `competing_exposures()` warns and
#' returns the result of `competing_causes()`. The deprecated name will be
#' removed in a later release.
#'
#' @param dag A dagitty object.
#' @returns Character vector of competing cause names, as `competing_causes()`.
#' @examples
#' dag <- build_graph(
#'   treatments = "T", outcomes = "Y", competing_causes = "C",
#'   type = "saturated"
#' )
#' competing_exposures(dag)
#' @export
competing_exposures <- function(dag){
  warning(warningCondition(
    paste0("competing_exposures(): 'competing_exposures()' is deprecated as of ",
           "0.0.1. Use 'competing_causes()' instead."),
    class = "deprecatedWarning"
  ))

  nodes <- competing_causes(dag)

  return(nodes)
}


#' Proxy node names in a dagitty object
#'
#' `proxies()` is a dagitty wrapper that identifies nodes that are proxy variables for indicated unmeasured confounding (direct children of a latent variable that are direct parents of treatment or outcome) in a directed acyclic graph.
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @returns A character vector of proxy variable names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' proxies(dag)
#'
#' @export
proxies <- function(dag){

  # treatment
  treatments <- treatments(dag)

  # outcome
  outcomes <- get_outcome_names(dag)

  # latent variables
  latent_vars <- unobserved(dag)

  edges <- dagitty::edges(dag)

  treatment_parents <- parents_edge_order(edges, treatments)

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- parents_edge_order(edges, outcomes)
  treatment_children <- children_edge_order(edges, treatments)

  nodes_trt_to_y <- nodes_between_treatment_and_outcome(dag = dag,
                                                            treatments = treatments,
                                                            outcomes = outcomes)

  mediators <- outcome_parents[ ( outcome_parents %in% treatment_children | outcome_parents %in% nodes_trt_to_y ) ]

  # proxy
  latent_children <- children_edge_order(edges, latent_vars)

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
#' `colliders()` identifies nodes that are common children of both treatment
#' and outcome (`treatment -> C <- outcome`) in a directed acyclic graph. It
#' does not return arbitrary colliders with two parents.
#'
#' @importFrom dagitty children
#' @param dag A dagitty object.
#' @returns A character vector of collider names, or `character(0)` if none.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' colliders(dag)
#'
#' @export
colliders <- function(dag){

  colliders <- classify_node_roles(dag)$collider

  nodes <- as.character(colliders[ complete.cases(colliders) ])

  return(nodes)
}


#' Get instrumental variable names
#'
#' Finds candidate instruments for each treatment and outcome pair using
#' `dagitty::instrumentalVariables()`.
#'
#' @importFrom dagitty parents children descendants
#' @param dag A dagitty object.
#' @param usable_instruments When \code{TRUE}, exclude candidates that are
#'   descendants of their treatment. Defaults to \code{FALSE}.
#' @param details With `FALSE` (the default), the function returns a nested
#'   list grouped by treatment and outcome. With `TRUE`, it returns a data
#'   table containing treatment, outcome, instrument and conditioning set.
#' @returns A nested list, or a data table when \code{details = TRUE}.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
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
  outcomes <- get_outcome_names(dag)
  instrumental_vars <- instrumental_variables_helper(dag = dag,
                                                     treatments = treatments,
                                                     outcomes = outcomes,
                                                     usable_instruments = usable_instruments)
  if( isTRUE(details) ){
    return(instrumental_vars)
  }

  result <- lapply(treatments, function(treatment){
    per_outcome <- lapply(outcomes, function(outcome){
      instrument_names <- instrumental_vars$instrument[
        instrumental_vars$treatment == treatment &
          instrumental_vars$outcome == outcome
      ]
      if( length(instrument_names) == 0L ) {
        out <- NULL

        return(out)
      } else {
        unique_instruments <- unique(as.character(instrument_names))

        return(unique_instruments)
      }
    })
    names(per_outcome) <- outcomes

    return(per_outcome)
  })
  names(result) <- treatments
  return(result)
}


#' Minimal adjustment and single door sets
#'
#' `minimal_sets()` returns sets from `dagitty::adjustmentSets()`, ordered by
#' size. At most `num_sets` are returned.
#'
#' With `effect = "total"`, these are minimal sufficient adjustment sets.
#' With `effect = "direct"`, dagitty's single door criterion is used. These
#' sets identify a direct path coefficient in a linear structural equation
#' model when the corresponding arrow from treatment to outcome exists. The wrapper
#' does not check for that arrow. `get_identification(effect = "direct")`
#' reports the single door check separately from natural direct effects.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty adjustmentSets
#' @param dag A dagitty object.
#' @param treatment Vector of treatment node names.
#' @param outcome Vector of outcome node names.
#' @param effect Either `"total"` (the default) for adjustment sets for the
#'   total effect, or `"direct"` for single door sets.
#' @param num_sets Maximum number of sets to return, ordered by size.
#'   Defaults to five.
#' @returns A list of sets, with set sizes as names, or invisible `NULL` with a
#'   message when no set is available.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' minimal_sets(dag) # defaults to the total effect and 5 smallest sets
#'
#' minimal_sets(dag, effect = "direct") # single door sets
#'
#' minimal_sets(dag, effect = "direct", num_sets = 1) # one single door set
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
      stop("Failed to compute adjustment sets: ", conditionMessage(e),
           "\n Please check that a valid dagitty object is supplied and try again.",
           call. = FALSE)
    })

  if( length(adjustment_sets) == 0){

    # Indicates why no set is available (e.g., cycles, latent variables)
    has_latents <- tryCatch(length(unobserved(dag)) > 0, error = function(e) {
      matches <- FALSE

      return(matches)
    })
    is_cyclic <- isFALSE(tryCatch(dag_is_acyclic(dag), error = function(e) {
      acyclic <- NA

      return(acyclic)
    }))
    why <- character(0)
    if( isTRUE(is_cyclic) ){ why <- c(why, "the graph is not acyclic (contains a cycle)") }
    if( isTRUE(has_latents) ){ why <- c(why, "latent variables are present") }
    reason <- if( length(why) > 0 ){ paste0("\n This can occur because ", paste(why, collapse = " and "), ".") }else{ "" }
    message("No available adjustment sets for the ", effect, " effect.", reason,
            "\n Try adjusting parameters, or assess edges using assess_edges(dag, assess_causal_criteria = TRUE).")

    adjustment_result <- NULL
    return(invisible(adjustment_result))
  }
  adjustment_sets_list <- lapply(seq_along(adjustment_sets), function(x){
    adjustment_sets <- adjustment_sets[[x]]

    return(invisible(adjustment_sets))
  })

  set_length <- sapply(adjustment_sets_list, length)
  names(adjustment_sets_list) <- set_length
  adjustment_sets_list <- adjustment_sets_list[order(set_length)]
  num_available <- length(adjustment_sets_list)

  if( num_available > num_sets &&
      length( adjustment_sets_list[[num_sets]] ) == length( adjustment_sets_list[[num_sets + 1]] ) ){
    warning("A same length adjustment set has been excluded. If this is unintentional, please adjust the num_sets parameter input and run again.")
  }
  adjustment_sets <- adjustment_sets_list[ 1:min(num_sets, num_available) ]

  return(adjustment_sets)
}


#' node names in path from treatment to outcome
#'
#' nodes_between_treatment_and_outcome() is a ggdag::dag_paths() wrapper intended for use with multiple treatments and outcomes.
#'
#' @importFrom ggdag dag_paths
#' @importFrom data.table as.data.table data.table
#' @param dag A dagitty object
#' @param treatments A vector of treatment node names.
#' @param outcomes A vector of outcome node names.
#' @param output_list Whether to keep results grouped by treatment and outcome.
#'   Defaults to FALSE.
#' @param directed When TRUE (default), restricts to directed paths from treatment to outcome; when FALSE, considers all paths.
#' @returns Unique node names, or a nested list when output_list is TRUE.
#'   Returns an empty list if treatments or outcomes is empty.
#' @noRd
nodes_between_treatment_and_outcome <- function(dag,
                                                treatments,
                                                outcomes,
                                                output_list = FALSE,
                                                directed = TRUE
                                                ){
  nodes <- c()
  treatments <- unlist(treatments)
  outcomes <- unlist(outcomes)

  paths_trt_to_y <- list()

  if( length(treatments) > 0 & length(outcomes) > 0 ){

    paths_trt_to_y <- lapply(1:length(treatments), function(x){

      outcome_paths <- lapply(1:length(outcomes), function(y){

        path_nodes <- tryCatch({
          paths_trt_y <- NA
          paths_trt_y <- ggdag::dag_paths(dag,
                                          from = treatments[x],
                                          to = outcomes[y],
                                          directed = directed,
                                          paths_only = TRUE)[["data"]]
          if(length(paths_trt_y) > 0 ){
            paths_trt_y <- as.vector(unlist(unique( paths_trt_y[ complete.cases(paths_trt_y["direction"]), "name" ])))
            paths_trt_y <- paths_trt_y[ !paths_trt_y %in% treatments]
          }
        }, error = function(e){
          paths_trt_y <- NA

          return(invisible(paths_trt_y))
        })

        return(path_nodes)
      })
      names(outcome_paths) <- outcomes

      return(outcome_paths)
    })

    names(paths_trt_to_y) <- treatments

    if(output_list == TRUE){

      return(paths_trt_to_y)
    }

    paths_trt_to_y <- as.vector(unique(unlist(paths_trt_to_y)))

    paths_trt_to_y <- paths_trt_to_y[ complete.cases(paths_trt_to_y) ]
  }

  return(paths_trt_to_y)
}
