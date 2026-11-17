# ==============================================================================
# get_dag_utils.R


# Inspect saved queries and checked alternatives
# ==============================================================================


## Inspection ##
#' Inspect queries and constructed graphs
#'
#' Read the saved records of a constructed DAG without changing it.
#'
#' @param dag A constructed DAG or a dry-run plan.
#' @return
#' - `original_query()`: The supplied query graph and its attributes.
#' - `temporal_dag()`: The temporal graph represented by the current display.
#' - `dag_options()`: The construction arguments, variables kept for each effect,
#'   size (`level$variables`) and number of drawn nodes (`level$nodes`). Size uses
#'   the counting rule in [construct_dag()]. The record also contains the time
#'   summary, identification results and unmeasured declarations.
#'
#' For `dag_options()`, `sizes` lists known sizes and `alternatives` lists DAGs
#' at the current size. `sizes_complete` and `alternatives_complete` say whether
#' these lists are complete within the construction family. `sizes_steps`
#' records the configured ranked search budget. `alternatives_limit` explains
#' why any candidates from the construction rules are unavailable.
#' @export
original_query <- function(dag){
  if( inherits(dag, "cf_dag_plan") ){
    return(original_query(dag$graph))
  }
  query <- .gd_provenance(dag)$original_query

  return(query)
}


#' @rdname original_query
#' @export
dag_options <- function(dag){
  if( inherits(dag, "cf_dag_plan") ){
    return(dag$options)
  }
  construction <- .gd_provenance(dag)$construction

  return(construction)
}


#' @rdname original_query
#' @export
temporal_dag <- function(dag){
  if( inherits(dag, "cf_dag_plan") ){
    return(temporal_dag(dag$graph))
  }
  temporal <- .gd_state(dag)$temporal

  return(temporal)
}


#' Inspect the resolved call
#'
#' Recover the call needed to rebuild the selected DAG.
#'
#' @param dag A constructed DAG or a dry-run plan.
#' @return A call that rebuilds the selected DAG from its saved source and
#'   construction settings. For a `get_dag()` result,
#'   `input_record$resolved_call` stores the original resolved query call and
#'   `input_record$responses` stores the interactive answers.
#' @export
dag_call <- function(dag){
  actual <- if( inherits(dag, "cf_dag_plan") ) dag$graph else dag
  if( inherits(actual, "cf_dag") ){
    state <- .gd_state(actual)
    if( identical(state$construction$metric, "variables-v1") && !is.null(state$construction$candidate) ){
      arguments <- list(graph = state$original_query,
        arguments = state$construction$arguments, candidate = state$construction$candidate,
        dry_run = inherits(dag, "cf_dag_plan"))
      input_record <- .cd_input_record(dag)
      if( !is.null(input_record) ){
        arguments$input_record <- input_record
      }
      return(as.call(c(list(call(":::", as.name("causaliflower"), as.name(".cd_replay"))), arguments)))
    }
    if( !identical(state$construction$metric, "variables-v1") ){
      migrated <- select_dag(dag, quiet = TRUE)
      replay <- dag_call(migrated)
      replay$dry_run <- inherits(dag, "cf_dag_plan")
      return(replay)
    }
  }
  record <- attr(dag, "get_dag_record", exact = TRUE)
  if( !is.null(record) ){
    return(record$resolved_call)
  }
  if( inherits(dag, "cf_dag_plan") ){
    return(.gd_saved_call(dag$arguments))
  }
  state <- .gd_provenance(dag)
  if( is.null(state$input_record$resolved_call) ){
    return(.gd_saved_call(c(list(graph = state$original_query), state$construction$arguments)))
  }
  resolved_call <- state$input_record$resolved_call

  return(resolved_call)
}


#' Inspect adjustment checks
#'
#' Read the saved checks and whether they still describe the supplied graph.
#'
#' @param dag A constructed DAG.
#' @return A list of saved checks for each effect in a DAG built by
#'   `construct_dag()`, recording:
#'   - Whether the retained sets are minimal valid sets of `temporal_dag(dag)`.
#'   - Whether the sets listed on that temporal DAG are valid on the source.
#'     This listing contains at most 50 sets per effect.
#'   - Which declared unmeasured occasions belong to adjustment members or
#'     undecided variables retained for that effect.
#'   - The recorded availability under those declarations for the displayed DAG
#'     (`availability`) and `temporal_dag(dag)` (`temporal_availability`).
#'
#'   The first two checks treat declared unmeasured variables as observed.
#'   Read availability alongside the search `state` and check `status`.
#'   `"available now"` alone does not establish that an adjustment set was found.
#'   The list's `current` attribute is `TRUE` if the graph still matches the saved
#'   structural signature, and `FALSE` otherwise. A mismatch produces a warning.
#' @export
dag_checks <- function(dag){
  state <- .gd_provenance(dag)
  checks <- state$checks
  checked <- state$signature
  current <- !is.null(checked) && identical(checked, .gd_signature(dag))
  attr(checks, "current") <- current
  if( !current ){
    warning("These checks describe the DAG before it was edited. They do not verify the edited DAG.", call. = FALSE)
  }
  return(checks)
}


## Printing methods ##
#' Print a constructed DAG
#'
#' @param x Constructed DAG to print as dagitty graph text.
#' @param ... Additional print arguments, currently unused.
#' @return `x`, invisibly.
#' @noRd
#' @export
print.cf_dag <- function(x, ...){
  cat(as.character(x), "\n", sep = "")
  return(invisible(x))
}


#' Print a construction preview
#'
#' @param x Dry-run plan with graph, options, checks and structural changes.
#' @param ... Additional print arguments, currently unused.
#' @return `x`, invisibly, after printing size, effect checks and vertex counts.
#' @noRd
#' @export
print.cf_dag_plan <- function(x, ...){
  cat("Proposed DAG construction (dry run)\n")
  .cd_print_level(x$options)
  cat("Final checks: ", paste(vapply(names(x$checks), function(name){
    return(paste0(name, " ", x$checks[[name]]$status))
  }, character(1)), collapse = ", "), "\n", sep = "")
  cat("Graph vertices: ", length(names(x$graph)), ", removed from reference: ", length(x$changes$removed_nodes), "\n", sep = "")
  return(invisible(x))
}


.cd_construction_version <- "variables-v1-root-variable-v1-explicit-helpers-v1-static-groups-v1-full-source-v1-member-sets-v1-diagnostic-choices-v1"


#' Rebuild a saved candidate and verify its identity
#'
#' @param graph Original query saved with the candidate.
#' @param arguments Saved construction arguments, including size and time representation.
#' @param candidate Saved source binding and temporal graph signature identifying
#'   the selected alternative.
#' @param dry_run `TRUE` to return a construction plan instead of the graph.
#' @param input_record Optional resolved wrapper call and interactive responses.
#' @return The rebuilt DAG or plan with its input record restored. Unsupported
#'   versions or changes to the saved source, settings or candidate signature
#'   raise an error.
#' @noRd
.cd_replay <- function(graph, arguments, candidate, dry_run = FALSE, input_record = NULL){
  if( !identical(candidate$binding$version, .cd_construction_version) ){
    stop("This saved selection uses an unsupported construction version. Rebuild the DAG.", call. = FALSE)
  }
  result <- .cd_sized(graph, arguments$latents, arguments$keep_nodes, dry_run, TRUE,
    arguments$size, unroll = isTRUE(arguments$unroll), expected = candidate,
    budget = candidate$binding$ranked_steps)
  return(.cd_with_input_record(result, input_record))
}


#' Read saved guided query inputs from a graph or plan
#'
#' @param dag Constructed DAG or dry-run plan.
#' @return The wrapper input record from its attribute or validated provenance,
#'   or `NULL` when no input record was saved.
#' @noRd
.cd_input_record <- function(dag){
  record <- attr(dag, "get_dag_record", exact = TRUE)
  if( inherits(dag, "cf_dag_plan") ){
    dag <- dag$graph
  }
  if( is.null(record) && inherits(dag, "cf_dag") ){
    record <- .gd_state(dag)$input_record
  }
  return(record)
}


#' Preserve guided query inputs through a rebuild
#'
#' @param dag Constructed DAG or dry-run plan receiving the record.
#' @param record Resolved wrapper call and responses, or `NULL` to leave the input unchanged.
#' @return The DAG with updated provenance, or the plan with a `get_dag_record`
#'   attribute. Returns `dag` unchanged when `record` is `NULL`.
#' @noRd
.cd_with_input_record <- function(dag, record){
  if( is.null(record) ){
    return(dag)
  }
  if( inherits(dag, "cf_dag") ){
    state <- .gd_state(dag)
    state$input_record <- record
    dag <- .gd_attach(dag, state)
  }else{
    attr(dag, "get_dag_record") <- record
  }
  return(dag)
}


#' Assemble a construction call from retained public arguments
#'
#' @param arguments Named list of saved construction inputs.
#' @return A `construct_dag()` call containing only arguments accepted by its
#'   current definition.
#' @noRd
.gd_saved_call <- function(arguments){
  arguments <- arguments[intersect(names(formals(construct_dag)), names(arguments))]
  saved_call <- as.call(c(list(as.name("construct_dag")), arguments))

  return(saved_call)
}
