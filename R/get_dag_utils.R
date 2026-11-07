################################################################################
# get_dag_utils.R


# Inspect saved queries, checked alternatives and traced DAG edits
# Original filename section markers below record source history

#####
# Inspection and traced editing
# Historical source: get_dag_utils.R
#####


#' Inspect queries and constructed graphs
#'
#' These functions accept the corresponding recovered graph or collection and
#' read saved records without changing the graph.
#'
#' @param dag The corresponding recovered graph or collection.
#' @param details With `TRUE`, include before and after graphs, typed edge changes
#'   and detailed records.
#' @return
#' - `original_query()`: Return the supplied query graph and its attributes.
#' - `temporal_dag()`: Return the temporal graph represented by the current display,
#'   or a named list for multiple views.
#' - `dag_options()`: Return the saved construction options, including separate
#'   direct and total candidate tables. Their validity refers to the original
#'   selection before subsequent editing commands.
#' - `dag_history()`: Return the transformation history as a data frame.
#'   `details = TRUE` includes before and after graphs, typed edge changes and
#'   detailed records. Multiple views have separate histories.
#' @export
original_query <- function(dag) { # Return the saved query without changing the recovered graph
  if (inherits(dag, "cf_dag_options")) {
    return(dag$graph)
  }
  if (inherits(dag, "cf_dag_plan")) {
    return(original_query(dag$graph))
  }
  return(.gd_provenance(dag)$original_query)
}


#' @rdname original_query
#' @export
dag_options <- function(dag) { # Return construction options and their recorded search assumptions
  if (inherits(dag, "cf_dag_options")) {
    return(dag)
  }
  if (inherits(dag, "cf_dag_plan")) {
    return(dag$options)
  }
  state <- .gd_provenance(dag)
  if (!is.null(state$construction)) {
    return(state$construction)
  }
  options <- state$build$candidates
  attr(options, "search") <- state$build$search
  attr(options, "availability") <- state$build$availability
  return(options)
}


#' Inspect and display quantity sets
#'
#' Use the saved construction search without rerunning the source query or source
#' adjustment search. Direct and total effects remain separate. Options with the
#' same source quantities are combined, and known proper supersets are omitted
#' from the new catalogue. The original search record is retained. An incomplete
#' source search cannot establish that this catalogue is complete.
#'
#' @param dag An unedited constructed graph, or a previous set view of the same source.
#' @param effect `"direct"` or `"total"`.
#' @param index A positive row index in that effect's catalogue.
#' @param simplify With `TRUE`, try summarising repeated time occasions. If the
#'   checks fail, retain the temporal graph and record why. `FALSE` retains exact
#'   temporal names.
#'
#' @details
#' Construction retains occurrences of the selected quantities before attempting
#' a time summary. `temporal_dag()` preserves that temporal graph. A summary is
#' accepted only when the selected source memberships map to minimal valid sets
#' in the fully checked display. This checks the selected effect and sets. It
#' does not establish general causal equivalence, validate the union as an
#' adjustment set, or guarantee preservation of other effects and alternatives.
#' Changes to structure or declarations invalidate switching. Coordinates may change.
#'
#' @return `dag_set_catalogue()` returns direct and total tables, completeness flags,
#' original source sets, omitted supersets and source/method identifiers. Tables
#' retain stable identifiers, row indexes, source variables, concepts, exact sets
#' and original option indexes.
#'
#' `dag_set_view()` returns a graph retaining the chosen quantities, requested
#' nodes, endpoints and required treatment parents. Its saved `set_view` record
#' contains the selection, mapped sets and summary outcome. `dag_checks()` checks
#' only the selected effect for that view. The original construction record is retained.
#' @export
dag_set_catalogue <- function(dag) { # Group the saved source alternatives into canonical families
  state <- .gd_set_source(dag)
  search <- state$build$search
  source_order <- .cd_record(state$original_query)$source_order
  families <- c("direct", "total") # Keep the effect families separate
  grouped <- stats::setNames(lapply(families, function(effect) { # Canonicalise the options within this effect family
    return(.gd_set_family(search[[effect]]$table, effect, source_order))
  }), families)
  complete <- stats::setNames(vapply(search[families], function(family_search) { # Read whether this source family was completely enumerated
    return(isTRUE(family_search$complete))
  }, logical(1)), families)
  return(list(direct = grouped$direct$table, total = grouped$total$table,
    complete = complete, incomplete = !complete,
    source_sets = lapply(search[families], `[[`, "sets"),
    supersets = lapply(grouped, `[[`, "supersets"),
    source_method = state$construction$search_method,
    enumeration = if (identical(state$construction$search_method, "source-witness-v1")) {
      "witness"
    } else if (all(complete)) {
      "complete"
    } else if (identical(state$construction$search_method, "quantity-ranked-v1")) {
      "ranked"
    } else {
      "incomplete"
    },
    method = "canonical-source-members-v1",
    source_signature = state$construction$refinement$source_signature,
    scope = "Source member sets minimal under set inclusion within each returned family. Incomplete searches may omit alternatives or smaller subsets."))
}


#' @rdname dag_set_catalogue
#' @export
dag_set_view <- function(dag, effect, index, simplify = TRUE) { # Build one checked view from a saved canonical option
  .gd_flag(simplify, "simplify")
  effect <- match.arg(effect, c("direct", "total"))
  return(.gd_sets_build(dag, stats::setNames(list(index), effect), simplify, single = TRUE))
}


#' Display several adjustment alternatives together
#'
#' Retains the union of the chosen quantities as graph context. That union is
#' not proposed as an adjustment set. Each alternative is checked separately.
#'
#' @param dag An unedited constructed graph with a saved quantity set catalogue.
#' @param selections A nonempty named list of `direct`, `total`, or both, containing
#'   unique positive catalogue row indexes. An empty integer vector is allowed
#'   only for a selected effect with no available rows. An omitted effect is not checked.
#' @param simplify Try a checked time summary with `TRUE`, or retain exact temporal
#'   occurrences with `FALSE`.
#'
#' @details
#' The source search is not rerun. Each selected effect is checked on the temporal
#' and displayed graphs without a limit on the number of adjustment sets returned
#' by those validation searches. A summary must preserve every selected
#' alternative. Otherwise the view retains time. No adjustment set, incomplete
#' search and a valid empty set remain distinct results. These checks do not
#' establish general causal equivalence.
#'
#' @return A graph with separate selections, exact sets, mappings, completeness
#' and summary records. `dag_checks()` reports each selected effect separately.
#' Original queries, options and option indexes are unchanged.
#' @export
dag_sets_view <- function(dag, selections, simplify = TRUE) { # Build one context retaining separately checked source alternatives
  .gd_flag(simplify, "simplify")
  return(.gd_sets_build(dag, selections, simplify, single = FALSE))
}


#' @rdname original_query
#' @export
temporal_dag <- function(dag) { # Return the saved temporal graph for each recovered view
  if (inherits(dag, "cf_dag_plan")) {
    return(temporal_dag(dag$graph))
  }
  if (inherits(dag, "cf_dags")) {
    return(lapply(dag, temporal_dag))
  }
  return(.gd_state(dag)$temporal)
}


#' @rdname original_query
#' @export
dag_history <- function(dag, details = FALSE) { # Describe the recorded transformations or return their full receipts
  if (inherits(dag, "cf_dag_plan")) {
    return(dag_history(dag$graph, details))
  }
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("details must be TRUE or FALSE.", call. = FALSE)
  }
  if (inherits(dag, "cf_dags")) {
    return(lapply(dag, dag_history, details = details))
  }
  history <- .gd_provenance(dag)$history
  if (details) {
    return(history)
  }
  return(data.frame(step = seq_along(history),
    operation = vapply(history, `[[`, character(1), "operation"),
    added_nodes = vapply(history, function(history_entry) { # Count the nodes added by this recorded operation
      return(length(history_entry$delta$added_nodes))
    }, integer(1)),
    removed_nodes = vapply(history, function(history_entry) { # Count the nodes removed by this recorded operation
      return(length(history_entry$delta$removed_nodes))
    }, integer(1)),
    added_edges = vapply(history, function(history_entry) { # Count the edges added by this recorded operation
      return(nrow(history_entry$delta$added_edges))
    }, integer(1)),
    removed_edges = vapply(history, function(history_entry) { # Count the edges removed by this recorded operation
      return(nrow(history_entry$delta$removed_edges))
    }, integer(1)),
    status_changes = vapply(history, function(history_entry) { # Count the role declarations changed by this operation
      return(length(history_entry$delta$statuses))
    }, integer(1)),
    stringsAsFactors = FALSE))
}



#' Inspect the resolved call
#'
#' These functions accept the corresponding recovered graph or collection and
#' read saved records without changing the graph.
#'
#' @param dag The corresponding recovered graph or collection.
#' @return Return the resolved call, including interactive choices.
#' @export
dag_call <- function(dag) { # Return the resolved construction call and flag obsolete receipts
  construction <- if (inherits(dag, "cf_dag_options")) {
    dag
  } else if (inherits(dag, "cf_dag_plan")) {
    dag$options
  } else {
    attr(dag, "cf_provenance", exact = TRUE)$construction
  }
  if (!is.null(construction) && !isTRUE(construction$search_method %in% c("complete-source-sets-v2", "source-witness-v1", "quantity-ranked-v1"))) {
    warning("This saved call predates the current search method. Its option indexes may now select different sets. Rebuild from original_query(dag) and inspect the new catalogue before replaying it.", call. = FALSE)
  }
  if (identical(construction$search_method, "quantity-ranked-v1") &&
    length(construction$arguments$latents) && is.null(construction$unmeasured)) {
    warning("This saved call predates the rule that declared unmeasured variables do not change the search. Its option indexes may now select different sets. Rebuild from original_query(dag) and inspect the new catalogue before replaying it.", call. = FALSE)
  }
  record <- attr(dag, "get_dag_record", exact = TRUE)
  if (!is.null(record)) {
    return(record$resolved_call)
  }
  if (inherits(dag, "cf_dag_options") || inherits(dag, "cf_dag_plan")) {
    return(as.call(c(list(as.name("construct_dag")), dag$arguments)))
  }
  return(.gd_provenance(dag)$input_record$resolved_call)
}



#' Inspect adjustment checks
#'
#' These functions accept the corresponding recovered graph or collection and
#' read saved records without changing the graph.
#'
#' @param dag The corresponding recovered graph or collection.
#' @return Return the reference and displayed adjustment sets and the
#'   checks of their preservation.
#' @export
dag_checks <- function(dag) { # Return saved adjustment checks with their current validity flag
  state <- .gd_provenance(dag)
  checks <- state$checks
  current <- !is.null(state$checks_signature) &&
    identical(state$checks_signature, .gd_signature(dag))
  attr(checks, "current") <- current
  if (!current) {
    warning("These checks describe the selected graph before later edits. They do not verify the edited graph.", call. = FALSE)
  }
  return(checks)
}


#####
# get_dag_selection.R
#####

#' Select another constructed graph
#'
#' `select_dag()` rebuilds a graph from another choice in its saved catalogue.
#' `direct_select` and `total_select` accept a group and option index, `"auto"`
#' or `NULL`. Groups increase by the number of distinct source quantities.
#' Omitted selectors retain their previous values. `selection` remains an alias
#' for `direct_select`.
#'
#' @param dag A constructed graph with its saved catalogue.
#' @param direct_select,total_select A group and option index, `"auto"` or `NULL`.
#'   Groups increase by the number of distinct source quantities.
#'   Omitted selectors retain their previous values.
#' @param selection An alias for `direct_select`.
#' @return A graph rebuilt from another choice in its saved catalogue.
#' @export
select_dag <- function(dag, selection = NULL, direct_select = NULL, total_select = NULL) { # Reconstruct a saved selection without replacing its source catalogue
  state <- .gd_state(dag)
  if (!is.null(state$construction)) {
    if (identical(state$construction$search_method, "source-witness-v1")) {
      stop("This receipt contains source witnesses, not the complete source catalogue. Use construct_dag(original_query(dag), search = 'complete') to obtain complete option indexes. dag_set_view() can use the saved witness rows.", call. = FALSE)
    }
    if (!isTRUE(state$construction$search_method %in% c("complete-source-sets-v2", "quantity-ranked-v1"))) {
      stop("This saved catalogue predates the current search method. Rebuild from original_query(dag), ",
        "inspect dag_options(), then choose an index from that new catalogue.", call. = FALSE)
    }
    # A ranked catalogue saved before construction ran as observed excluded the
    # options that use a declared unmeasured variable (character() declares none).
    if (identical(state$construction$search_method, "quantity-ranked-v1") &&
      length(state$construction$arguments$latents) && is.null(state$construction$unmeasured)) {
      stop("This ranked receipt was saved before declared unmeasured variables stopped changing the search. ",
        "Its option indexes refer to a different catalogue. Rebuild from original_query(dag) and choose from the new catalogue.", call. = FALSE)
    }
    if (!identical(state$selection_base_signature, .gd_signature(dag))) {
      stop("Commands changed this selected graph. Undo those edits before switching options.", call. = FALSE)
    }
    if (!missing(selection) && !missing(direct_select)) {
      stop("Supply selection or direct_select, not both.", call. = FALSE)
    }
    arguments <- state$construction$arguments
    arguments$dry_run <- FALSE
    if (!missing(selection)) {
      arguments["direct_select"] <- list(selection)
    }
    if (!missing(direct_select)) {
      arguments["direct_select"] <- list(direct_select)
    }
    if (!missing(total_select)) {
      arguments["total_select"] <- list(total_select)
    }
    # Receipts saved before search modes existed used the complete route.
    if (is.null(arguments$search)) {
      arguments$search <- "complete"
    }
    selected_result <- if (identical(state$construction$search_method, "quantity-ranked-v1")) {
      # Indexes refer to the saved ranked catalogue, so reconstruct under its step
      # budget, passed as the steps argument. Receipts saved before steps was recorded
      # keep that budget in their search identity (without it, the option applies).
      if (is.null(arguments$steps)) {
        arguments$steps <- state$construction$search_identity$ranked_steps
      }
      do.call(.cd_reconstruct, c(list(previous = dag), arguments))
    } else {
      do.call(construct_dag, arguments)
    }
    input_history <- c(state$input_history, list(state$input_record))
    if (!inherits(selected_result, "cf_dag")) {
      selected_result$input_history <- input_history
      return(selected_result)
    }
    next_state <- .gd_state(selected_result)
    next_state$input_history <- input_history
    construction_history <- next_state$history
    selected_result <- .gd_attach(selected_result, next_state)
    return(.gd_record(dag, selected_result, "select_dag",
      arguments[c("direct_select", "total_select")],
      list(construction_history = construction_history, selected = next_state$selected)))
  }
  if (is.null(state$catalogue)) {
    stop("This object predates grouped selection. Rerun get_dag.", call. = FALSE)
  }
  if (!identical(state$selection_base_signature, .gd_signature(dag))) {
    stop("Commands changed this selected graph. Undo those edits before switching options.", call. = FALSE)
  }
  record <- state$input_record
  record$responses <- record$responses[0, , drop = FALSE]
  record$values$selection <- selection
  record$resolved_call$selection <- selection
  return(.gd_select_result(state$catalogue, selection, record))
}

#####
# get_dag_utils.R
#####

#' Restore named quantities
#'
#' `restore_nodes()` restores named quantities from the original temporal query.
#' A quantity key or alias restores all its occurrences. An exact temporal name
#' restores that occurrence. Explicit requests can admit `yes:user`. Exclusions
#' recorded as `no:reason` require a metadata correction. `dry_run = TRUE`
#' previews changes without applying them.
#'
#' @param dag A constructed graph with its saved original query.
#' @param node_name Named quantities to restore. A quantity key or alias restores
#'   all its occurrences. An exact temporal name restores that occurrence.
#' @param dry_run With `TRUE`, preview changes without applying them.
#' @return A graph. Preview modes return the proposed changes invisibly.
#' @export
restore_nodes <- function(dag, node_name, dry_run = FALSE) { # Restore requested quantities from the saved temporal query
  .gd_flag(dry_run, "dry_run")
  state <- .gd_state(dag)
  requested <- .gd_names(node_name, "node_name")
  if (!length(requested)) {
    stop("Supply node_name.", call. = FALSE)
  }
  nodes <- .gd_match_reference(state$original_query, requested, state$controls)
  added <- setdiff(nodes, state$retained)
  if (!length(added)) {
    if (dry_run) {
      return(.gd_preview(dag, dag, "restore_nodes", list(requested = requested)))
    }
    return(dag)
  }
  restored_graph <- .gd_restore_build(dag, state, added)
  details <- list(restored = added, source = "original_query", summary = .gd_state(restored_graph)$summary)
  if (dry_run) {
    return(.gd_preview(dag, restored_graph, "restore_nodes", details))
  }
  return(.gd_record(dag, restored_graph, "restore_nodes", list(node_name = requested), details))
}


#' Restore omitted common causes
#'
#' `add_confounders()` searches the saved query for omitted common causes with
#' separate directed branches reaching a selected node and another displayed
#' quantity. Paths stop at displayed nodes. A chain through one displayed node
#' alone is not treated as confounding, and two occasions of one quantity do not
#' count as two variables. This restores supplied graph structure, rather than
#' inferring a new causal relationship or recommending adjustment.
#'
#' @param dag A constructed graph with its saved original query.
#' @param role A displayed role, defaulting to `"outcome"`.
#' @param node_name A quantity to select instead.
#' @param unmeasured With `TRUE`, restrict results to common causes absent from
#'   supplied columns. `FALSE` also includes supplied common causes. Restored
#'   quantities absent from columns are marked latent, including when the
#'   initial search assumed full observation.
#' @param include Restrict restoration to named matches.
#' @param dry_run Preview the changes.
#' @return A graph. Preview modes return the proposed changes invisibly.
#' @export
add_confounders <- function(dag, role = NULL, node_name = NULL,
  unmeasured = TRUE, dry_run = FALSE, include = NULL) { # Restore omitted causes with separate branches to displayed quantities
  .gd_flag(unmeasured, "unmeasured")
  .gd_flag(dry_run, "dry_run")
  state <- .gd_state(dag)
  reference <- state$analysis
  targets <- .gd_target_nodes(dag, state, role, node_name)
  ids <- .gd_node_ids(reference)
  current_concepts <- unique(unname(ids[state$retained]))
  frontiers <- .gd_visible_frontiers(reference, state$retained)
  witnesses <- data.frame(cause = character(), target = character(), other = character())
  for (cause in names(frontiers)) {
    # A prior occasion of an already displayed quantity is explored explicitly
    # through restore_nodes/temporal_dag, not minted as a new confounder.
    if (ids[[cause]] %in% current_concepts) {
      next
    }
    if (unmeasured && ids[[cause]] %in% state$build$observed_ids) {
      next
    }
    reached <- frontiers[[cause]]
    for (target in intersect(reached, targets)) {
      others <- reached[ids[reached] != ids[[target]]]
      if (length(others)) {
        witnesses <- rbind(witnesses,
          data.frame(cause = cause, target = target, other = others))
      }
    }
  }
  witnesses <- unique(witnesses)
  nodes <- unique(witnesses$cause)
  if (!is.null(include)) {
    include <- .gd_names(include, "include")
    if (!length(include)) {
      stop("include must name at least one common cause.", call. = FALSE)
    }
    matches <- lapply(include, function(name) { # Match a requested cause within the eligible omitted causes
      return(intersect(.gd_match_reference(reference, name, state$controls,
        .cd_record(state$original_query)), nodes))
    })
    absent <- include[!lengths(matches)]
    if (length(absent)) {
      stop("include contains quantities that do not match the selected common causes: ",
        paste(absent, collapse = ", "), call. = FALSE)
    }
    nodes <- unique(unlist(matches, use.names = FALSE))
    witnesses <- witnesses[witnesses$cause %in% nodes, , drop = FALSE]
  }
  details <- list(witnesses = witnesses, unmeasured = unmeasured,
    availability = "Restored nodes are unmeasured unless supplied through columns.",
    source = "full recovered temporal reference")
  if (!length(nodes)) {
    if (dry_run) {
      return(.gd_preview(dag, dag, "add_confounders", details))
    }
    message("No omitted ", if (unmeasured) {
      "unmeasured "
    } else {
      ""
    },
      "common causes matched this selection.")
    return(dag)
  }
  expanded_graph <- .gd_restore_build(dag, state, nodes, confounders = TRUE)
  details$restored <- nodes
  details$summary <- .gd_state(expanded_graph)$summary
  if (dry_run) {
    return(.gd_preview(dag, expanded_graph, "add_confounders", details))
  }
  return(.gd_record(dag, expanded_graph, "add_confounders",
    list(role = role, node_name = node_name, unmeasured = unmeasured, include = include), details))
}


#' Undo the most recent editing command
#'
#' `undo_dag()` restores the graph and temporal reference from before the most
#' recent editing command and records the reversal. `dry_run` previews it.
#' Initial query transformations can be inspected through the history and
#' temporal graph, but are not individually undone by this command.
#'
#' @param dag A constructed graph with its recorded editing commands.
#' @param dry_run Preview the reversal.
#' @return A graph. Preview modes return the proposed changes invisibly.
#' @export
undo_dag <- function(dag, dry_run = FALSE) { # Restore the preceding graph and record the reversal
  .gd_flag(dry_run, "dry_run")
  state <- .gd_state(dag)
  if (is.null(state$undo) || !length(state$undo)) {
    stop("There is no editing command to undo.", call. = FALSE)
  }
  previous <- state$undo[[length(state$undo)]]
  next_state <- previous$state
  next_state$history <- state$history
  next_state$undo <- head(state$undo, -1L)
  restored_graph <- .gd_attach(previous$graph, next_state)
  if (dry_run) {
    return(.gd_preview(dag, restored_graph, "undo_dag", list()))
  }
  return(.gd_record(dag, restored_graph, "undo_dag", list(), list(), undo = FALSE))
}


#####
# Printing methods
# Historical source: get_dag_utils.R
#####


#' @export
print.cf_dag <- function(x, ...) { # Print the recovered DAG and return it invisibly
  cat(as.character(x), "\n", sep = "")
  return(invisible(x))
}


#' @export
print.cf_dags <- function(x, ...) { # Print each named recovered DAG and return the collection invisibly
  for (view_index in seq_along(x)) {
    cat("# ", names(x)[view_index], "\n", sep = "")
    print(x[[view_index]], ...)
    if (view_index < length(x)) {
      cat("\n")
    }
  }
  return(invisible(x))
}


#' @export
print.cf_dag_options <- function(x, ...) { # Print construction choices and their search certificates
  cat("Construction options\n")
  if (identical(x$search_method, "quantity-ranked-v1")) {
    for (effect in c("direct", "total")) {
      ranked <- x$search[[effect]]$ranked
      if (!is.null(ranked)) {
        cat(effect, " options ranked by distinct source members: ",
          if (identical(ranked$certificate, "complete")) {
            "complete"
          } else if (identical(ranked$certificate, "validation")) {
            "not certified (an option failed validation and was dropped)"
          } else if (!isTRUE(ranked$complete_up_to >= 0)) {
            "not certified (the step budget stopped the search before any member count was completed)"
          } else {
            paste0("certified complete up to ", ranked$complete_up_to, " members")
          }, ".\n", sep = "")
      }
    }
  }
  if (identical(x$search_method, "source-witness-v1")) {
    cat("Deterministic source witnesses. Minimum quantity cost and the full source catalogue were not searched. Final checks cover the constructed graph.\n")
  }
  cat("Counts refer to distinct source nodes: time repeats count once. Recorded groups are expanded.\n")
  if (length(x$unmeasured$declared)) {
    cat("Declared unmeasured variables: ", paste(x$unmeasured$declared, collapse = ", "),
      " (marked on the finished DAG after construction). They never change the options or the selection.\n", sep = "")
  }
  if (length(x$refinement$passes)) {
    before <- x$refinement$passes[[1L]]$selected_before
    after <- x$refinement$selected_counts
    improved <- names(before)[is.finite(before) & is.finite(after) & after < before]
    if (length(improved)) {
      cat("Verification against the original query accepted smaller options: ",
        paste(paste0(improved, " ", before[improved], " to ", after[improved]), collapse = ", "),
        " distinct nodes. Search remains bounded.\n", sep = "")
    }
  }
  for (effect in c("direct", "total")) {
    cat(effect, ": ", x$status[[effect]], "\n", sep = "")
    .gd_print_options(x[[effect]], x$choices[[effect]]$index)
    if (length(x$choices[[effect]]$ties) > 1L) {
      cat("  ", length(x$choices[[effect]]$ties), " tied choices. The first is selected.\n", sep = "")
    }
  }
  return(invisible(x))
}


#' @export
print.cf_dag_plan <- function(x, ...) { # Print a proposed construction and its final checks
  cat("Proposed DAG construction (dry run)\n")
  print(x$options)
  cat("Final checks: ", paste(vapply(names(x$checks), function(name) { # Describe the saved status of this effect family check
    return(paste0(name, " ", x$checks[[name]]$status))
  }, character(1)), collapse = ", "), "\n", sep = "")
  cat("Graph vertices: ", length(names(x$graph)), ", removed from reference: ", length(x$changes$removed_nodes), "\n", sep = "")
  return(invisible(x))
}



#####
# dag_provenance.R
#####

#' @export
print.cf_dag_change <- function(x, ...) { # Print the proposed structural and role changes
  cat(x$operation, " (dry run)\n", sep = "")
  show_names <- function(title, values) { # Print a nonempty collection of changed node names
    return(if (length(values)) {
      cat(title, paste(values, collapse = ", "), "\n", sep = "")
    })
  }
  show_edges <- function(title, edges) { # Print the typed edge changes when any are present
    return(if (nrow(edges)) {
      cat(title, "\n", sep = "")
      cat(paste0("  ", edges$v, " ", edges$e, " ", edges$w), sep = "\n")
      cat("\n")
    })
  }
  show_names("Add nodes: ", x$delta$added_nodes)
  show_names("Remove nodes: ", x$delta$removed_nodes)
  show_edges("Add edges:", x$delta$added_edges)
  show_edges("Remove edges:", x$delta$removed_edges)
  for (name in names(x$delta$statuses)) {
    show_names(paste0("Set ", name, ": "), x$delta$statuses[[name]]$added)
    show_names(paste0("Clear ", name, ": "), x$delta$statuses[[name]]$removed)
  }
  if (!length(x$delta$added_nodes) && !length(x$delta$removed_nodes) &&
    !nrow(x$delta$added_edges) && !nrow(x$delta$removed_edges) && !length(x$delta$statuses)) {
    cat("No structural or status changes.\n")
  }
  return(invisible(x))
}



#####
# dag_candidates.R
#####

#' @export
print.cf_dag_views <- function(x, ...) { # Print the recovered views and their availability limits
  cat("DAG views from query_graph()\n", x$availability$statement, "\n", sep = "")
  cat(x$recovery$coverage, "\n", sep = "")
  if (!nrow(x$candidates)) {
    cat("No adjustment set view found in this analysis graph. Use focus for a named process view.\n")
  } else {
    print(x$candidates[, intersect(c("id", "n_concepts", "n_resolved", "full_graph_valid", "status", "display_nodes", "label"), names(x$candidates)), drop = FALSE], row.names = FALSE)
  }
  if (x$search$mode == "sets" && identical(x$search$scope, "paths")) {
    cat("Process options based on paths: full_graph_valid checks their raw sets in the larger analysis graph.\n")
  }
  if (identical(x$search$complete, FALSE)) {
    cat("Search limit reached. The table may omit alternatives.\n")
  }
  cat("Use $dags[[i]] for display. $selection, $analysis and $sets retain the temporal checks.\n")
  return(invisible(x))
}


#####
# Traced editing internals
# Historical source: get_dag_utils.R
#####


.gd_set_source <- function(dag) { # Validate the saved source and catalogue before constructing a view
  state <- .gd_state(dag)
  construction <- state$construction
  if (is.null(construction) ||
    !(identical(construction$search_method, "source-verified-refinement-v1") ||
    identical(construction$search_method, "complete-source-sets-v2") ||
    identical(construction$search_method, "quantity-ranked-v1") ||
    identical(construction$search_method, "source-witness-v1"))) {
    stop("This saved source uses an unsupported search method. Reconstruct it before requesting set views.", call. = FALSE)
  }
  if (!identical(state$selection_base_signature, .gd_signature(dag))) {
    stop("Commands changed this selected graph. Undo those edits before requesting set views.", call. = FALSE)
  }
  # A receipt constructed as observed searched the analysis without its declared marks.
  source <- if (identical(construction$unmeasured$semantics, "observed-as-if-v1")) {
    .cd_unmarked(state$analysis)
  } else {
    state$analysis
  }
  if (!identical(state$original_query, construction$arguments$graph) ||
    is.null(construction$refinement$source_signature) ||
    !identical(construction$refinement$source_signature, .gd_signature(source)) ||
    !identical(state$build$search, construction$search)) {
    stop("The saved source or search receipt no longer matches this catalogue. Reconstruct it first.", call. = FALSE)
  }
  for (effect in c("direct", "total")) {
    item <- state$build$search[[effect]]
    table <- item$table
    if (!is.list(item) || !is.logical(item$complete) || length(item$complete) != 1L ||
      is.na(item$complete) || !is.data.frame(table) ||
      !all(c("resolved", "concepts", "raw_sets", "id", "group", "option") %in% names(table))) {
      stop("The saved ", effect, " candidate receipt is incomplete. Reconstruct it first.", call. = FALSE)
    }
  }
  return(state)
}


.gd_sets_build <- function(dag, selections, simplify, single) { # Build and check the selected canonical source alternatives
  catalogue <- dag_set_catalogue(dag)
  if (!is.list(selections) || is.data.frame(selections) || !length(selections) ||
    is.null(names(selections)) || anyNA(names(selections)) ||
    anyDuplicated(names(selections)) ||
    any(!names(selections) %in% c("total", "direct"))) {
    stop("selections must be a nonempty named list with distinct total/direct families.", call. = FALSE)
  }
  effects <- intersect(c("total", "direct"), names(selections))
  selections <- selections[effects]
  for (effect in effects) {
    selection_indices <- selections[[effect]]
    valid <- is.numeric(selection_indices) && !is.complex(selection_indices) && !anyNA(selection_indices) && all(is.finite(selection_indices)) &&
      all(selection_indices >= 1 & selection_indices == floor(selection_indices) & selection_indices <= nrow(catalogue[[effect]])) &&
      !anyDuplicated(selection_indices)
    if (single) {
      valid <- valid && length(selection_indices) == 1L
    } else {
      valid <- valid && (length(selection_indices) > 0L || nrow(catalogue[[effect]]) == 0L)
    }
    if (!valid) {
      stop(if (single) {
        "index"
      } else {
        paste0("selections$", effect)
      },
        " must identify ", if (single) {
          "an available row"
        } else {
          "unique available rows (or integer() for an unavailable family)"
        },
        " in the ", effect, " set catalogue.", call. = FALSE)
    }
    selections[[effect]] <- sort(as.integer(selection_indices))
  }
  state <- .gd_state(dag)
  uses_witnesses <- identical(state$construction$search_method, "source-witness-v1")
  arguments <- state$construction$arguments
  # Receipts constructed as observed build their views the same way and mark the
  # declared unmeasured variables afterwards. Other receipts are unchanged.
  as_observed <- identical(state$construction$unmeasured$semantics, "observed-as-if-v1")
  prepared <- .cd_prepare(state$original_query, if (as_observed) {
    NULL
  } else {
    arguments$latents
  }, arguments$keep_nodes)
  latent_nodes <- if (as_observed) {
    .cd_match_nodes(state$original_query, arguments$latents,
      .cd_record(state$original_query))
  } else {
    character()
  }
  if (!identical(.gd_signature(prepared$selection), catalogue$source_signature) ||
    !identical(attr(prepared$analysis, "cf_node_map"), attr(state$analysis, "cf_node_map"))) {
    stop("The saved source identities or preparation no longer match this catalogue. Reconstruct it first.", call. = FALSE)
  }
  selected_rows <- stats::setNames(lapply(effects, function(effect) { # Read the selected rows for this effect family
    return(catalogue[[effect]][selections[[effect]], , drop = FALSE])
  }), effects)
  source_sets_by_effect <- lapply(selected_rows, function(table) { # Retain distinct exact source memberships from these rows
    sets <- unlist(table$raw_sets, recursive = FALSE)
    if (is.null(sets)) {
      return(list())
    }
    return(sets[!duplicated(vapply(sets, .gd_set_key, character(1)))])
  })
  concepts <- as.character(unique(unlist(lapply(selected_rows, function(table) { # Read the canonical concepts retained by these rows
    return(unlist(table$concepts, use.names = FALSE))
  }), use.names = FALSE)))
  selected_nodes <- as.character(unique(unlist(source_sets_by_effect, recursive = TRUE, use.names = FALSE)))
  if (!all(selected_nodes %in% names(prepared$analysis))) {
    stop("Saved adjustment memberships are absent from the prepared reference.", call. = FALSE)
  }
  quantity_nodes <- names(prepared$analysis)[.gd_node_ids(prepared$analysis) %in% concepts]
  retained <- union(quantity_nodes, prepared$required)
  identification <- .cd_identification_retention(prepared$analysis,
    state$construction$choices$total)
  identification$added_nodes <- setdiff(identification$parents, retained)
  retained <- union(retained, identification$parents)
  omitted_nodes <- setdiff(names(prepared$analysis), retained)
  graph <- if (length(omitted_nodes)) {
    collapse_nodes(prepared$analysis, omitted_nodes)
  } else {
    prepared$analysis
  }
  graph <- .cd_working_graph(graph, prepared$map[prepared$map$node %in% names(graph), , drop = FALSE])
  if (isTRUE(identification$applied)) {
    edges <- dag_edges(graph)
    treatment <- dagitty::exposures(prepared$analysis)
    required <- union(treatment, union(dagitty::outcomes(prepared$analysis), identification$parents))
    if (!all(required %in% names(graph)) || length(dagitty::latents(graph)) ||
      any(edges$e == "<->" & (edges$v %in% treatment | edges$w %in% treatment))) {
      stop("Retention of exact treatment parents failed its structural invariant.", call. = FALSE)
    }
    identification$status <- "preserved by retaining exact treatment parents"
  }
  if (as_observed) {
    identification <- .cd_unmeasured_parents(identification, latent_nodes)
  }
  temporal <- graph
  required_nodes <- union(union(selected_nodes, prepared$required), identification$parents)
  summary <- if (simplify) {
    .gd_summary(temporal, keep = required_nodes,
      preserve_occasions = isTRUE(identification$applied),
      source_check = .cd_summary_source_check(prepared$selection, effects))
  } else {
    list(applied = FALSE, reason = "Time summarisation was not requested.",
      graph = temporal, temporal = temporal)
  }
  graph <- summary$graph
  source_status <- stats::setNames(vapply(effects, function(effect) { # Describe the source family's completeness and availability
    return(if (!isTRUE(catalogue$complete[[effect]])) {
      "incomplete"
    } else if (!nrow(catalogue[[effect]])) {
      "no_sets"
    } else {
      "available"
    })
  }, character(1)), effects)
  check_graph <- function(graph) { # Check the selected source memberships on the supplied graph
    return(stats::setNames(lapply(effects, function(effect) { # Check the selected memberships within this effect family
      final_search <- .gd_adjustment_sets(graph, effect, Inf)
      check <- .gd_check_family(prepared$selection, graph, source_sets_by_effect[[effect]], final_search$sets,
        effect, final_search$complete, source_sets_by_effect[[effect]])
      check$reference_sets <- source_sets_by_effect[[effect]]
      check$required <- source_sets_by_effect[[effect]]
      check$final_sets <- final_search$sets
      check$reference_complete <- catalogue$complete[[effect]]
      check$final_complete <- final_search$complete
      check$scope <- if (single) {
        "Selected canonical option's mapped source memberships only. This does not cover all source alternatives or the other effect. Summary checks do not establish general causal equivalence."
      } else {
        "Each included canonical option's mapped source memberships separately. The retained union is not an adjustment set. Unselected alternatives and disabled effects are not certified."
      }
      if (!single) {
        check$source_status <- source_status[[effect]]
        check$selected_indices <- selections[[effect]]
        if (!length(source_sets_by_effect[[effect]]) && check$status != "changed") {
          check$status <- if (!isTRUE(final_search$complete) || !isTRUE(catalogue$complete[[effect]])) {
            "not fully checked"
          } else {
            "no eligible reference set"
          }
        }
      }
      if (uses_witnesses) {
        check <- .cd_output_source_check(check, prepared$selection, graph, effect)
      }
      return(check)
    }), effects))
  }
  checks <- check_graph(graph)
  temporal_checks <- if (if (uses_witnesses) {
    identical(graph, temporal)
  } else {
    identical(.gd_signature(graph), .gd_signature(temporal))
  }) {
    checks
  } else {
    check_graph(temporal)
  }
  accepted <- vapply(checks, function(check) { # Check whether this effect family retained its selected memberships
    return(check$status %in% c("preserved", "no eligible reference set"))
  }, logical(1))
  if (isTRUE(summary$applied) && !all(accepted)) {
    rejected_status <- if (single) {
      checks[[1L]]$status
    } else {
      paste(paste0(names(checks)[!accepted], ": ", vapply(checks[!accepted], `[[`, character(1), "status")), collapse = ", ")
    }
    attempted <- list(summary = summary, checks = if (single) {
      checks[[1L]]
    } else {
      checks
    })
    summary <- list(applied = FALSE,
      reason = paste0("Summary checks for the selected family were ", rejected_status,
        ". Exact temporal occurrences were retained."),
      graph = temporal, temporal = temporal, attempted = attempted)
    graph <- temporal
    checks <- temporal_checks
  }
  if (uses_witnesses && any(!vapply(checks, function(check) { # Require complete original source validation for this witness family
    return(check$status %in% c("preserved", "no eligible reference set") &&
      isTRUE(check$output_reference_complete))
  }, logical(1)))) {
    stop("The witness view could not certify exact temporal memberships and all Output alternatives on the original source. No complete source search was attempted.", call. = FALSE)
  }
  if (as_observed) {
    marked_graph <- .cd_mark_unmeasured(summary, latent_nodes)
    summary <- marked_graph$summary
    graph <- summary$graph
    projected_latent <- intersect(names(temporal), latent_nodes)
    if (length(projected_latent)) {
      temporal <- .cd_set_latents(temporal, projected_latent)
    }
    checks <- .cd_annotate_unmeasured(checks, marked_graph$display)
    temporal_checks <- .cd_annotate_unmeasured(temporal_checks, projected_latent)
  }
  .gd_segment_tie_notice(summary)
  summary$requested <- simplify
  graph_size <- function(graph) { # Count the vertices and edges of this graph
    return(c(nodes = length(names(graph)), edges = nrow(dag_edges(graph))))
  }
  summary$counts <- list(temporal = graph_size(temporal),
    summary_temporal = graph_size(summary$temporal), displayed = graph_size(graph))
  mapped_sets_by_effect <- lapply(checks, function(check) { # Read the mapped memberships from the effect family checks
    return(lapply(check$records, `[[`, "mapped"))
  })
  if (single) {
    effect <- effects[[1L]]
    selected_row <- selected_rows[[effect]]
    selection_indices <- selections[[effect]]
    view <- list(id = selected_row$id[[1L]], effect = effect, index = selection_indices,
      variables = selected_row$variables[[1L]], concepts = selected_row$concepts[[1L]],
      raw_sets = source_sets_by_effect[[effect]], mapped_sets = mapped_sets_by_effect[[effect]], source_options = selected_row$source_options[[1L]],
      source_complete = catalogue$complete[[effect]], source_method = catalogue$source_method,
      catalogue_method = catalogue$method, identification = identification,
      summary = summary[c("requested", "applied", "reason", "counts")], scope = checks[[effect]]$scope)
  } else {
    view <- list(mode = "combined", selections = selections, alternatives = selected_rows,
      variables = as.character(sort(unique(unlist(lapply(selected_rows, function(table) { # Read the retained source variables from these selected rows
        return(unlist(table$variables, use.names = FALSE))
      }), use.names = FALSE)))),
      concepts = concepts, raw_sets = source_sets_by_effect, mapped_sets = mapped_sets_by_effect,
      source_complete = catalogue$complete[effects], source_status = source_status,
      source_method = catalogue$source_method, catalogue_method = catalogue$method,
      identification = identification, summary = summary[c("requested", "applied", "reason", "counts")],
      scope = "Retained context for separately checked source alternatives. This does not certify a union adjustment set or general causal equivalence.")
  }
  state$temporal <- temporal
  state$retained <- names(temporal)
  state$summary <- summary
  state$controls$simplify <- isTRUE(summary$applied)
  state$focus <- if (single) {
    selected_rows[[1L]]$concepts[[1L]]
  } else {
    concepts
  }
  state$selected <- list(direct = NULL, total = NULL)
  for (effect in effects) {
    state$selected[[effect]] <- if (single) {
      list(id = view$id, index = selections[[effect]], raw_sets = source_sets_by_effect[[effect]])
    } else {
      list(ids = selected_rows[[effect]]$id, indices = selections[[effect]], raw_sets = source_sets_by_effect[[effect]])
    }
  }
  state$set_view <- view
  state$checks <- checks
  state$temporal_checks <- temporal_checks
  state$checks_signature <- .gd_signature(graph)
  state$selection_base_signature <- .gd_signature(graph)
  view_graph <- .gd_attach(graph, state)
  operation <- if (single) {
    "dag_set_view"
  } else {
    "dag_sets_view"
  }
  operation_args <- if (single) {
    list(effect = effects[[1L]], index = selections[[1L]], simplify = simplify)
  } else {
    list(selections = selections, simplify = simplify)
  }
  view_graph <- .gd_record(dag, view_graph, operation, operation_args, view)
  if (any(!vapply(checks, function(check) { # Check whether this final family retained its selected memberships
    return(check$status %in% c("preserved", "no eligible reference set"))
  }, logical(1)))) {
    warning("Selected set view memberships changed or were not fully verified. Inspect dag_checks().", call. = FALSE)
  }
  return(view_graph)
}


.gd_set_family <- function(table, effect, source_order = NULL) { # Combine equal source member options and record proper supersets
  family_members <- unique(c(unlist(table$resolved, use.names = FALSE),
    unlist(table$concepts, use.names = FALSE), unlist(table$unmeasured, use.names = FALSE)))
  if( !.query_order_applies(family_members, source_order) ){
    source_order <- NULL
  }
  family_table <- data.frame(id = character(), index = integer(), label = character(),
    n_variables = integer(), stringsAsFactors = FALSE)
  family_table$variables <- list()
  family_table$resolved <- list()
  family_table$concepts <- list()
  family_table$raw_sets <- list()
  family_table$source_options <- list()
  # Ranked tables list the declared unmeasured members each option needs.
  has_unmeasured <- "unmeasured" %in% names(table)
  if (has_unmeasured) {
    family_table$unmeasured <- list()
  }
  membership_keys <- vapply(table$resolved, .gd_set_key, character(1))
  for (key in unique(membership_keys)) {
    option_rows <- which(membership_keys == key)
    variables <- .query_sort(unique(table$resolved[[option_rows[1L]]]), source_order)
    raw_sets <- unlist(table$raw_sets[option_rows], recursive = FALSE)
    raw_sets <- raw_sets[!duplicated(vapply(raw_sets, .gd_set_key, character(1)))]
    family_row <- data.frame(id = "", index = 0L,
      label = if (length(variables)) {
        paste(variables, collapse = " + ")
      } else {
        "No adjustment nodes"
      },
      n_variables = length(variables), stringsAsFactors = FALSE)
    family_row$variables <- family_row$resolved <- list(variables)
    family_row$concepts <- list(.query_sort(unique(unlist(table$concepts[option_rows], use.names = FALSE)), source_order))
    family_row$raw_sets <- list(raw_sets)
    family_row$source_options <- list(table[option_rows, c("id", "group", "option"), drop = FALSE])
    if (has_unmeasured) {
      family_row$unmeasured <- list(.query_sort(unique(as.character(unlist(table$unmeasured[option_rows],
        use.names = FALSE))), source_order))
    }
    family_table <- rbind(family_table, family_row)
  }
  if (nrow(family_table)) {
    label_order <- family_table$label
    if( !is.null(source_order) ){
      label_order <- vapply(family_table$variables, .query_set_key, character(1), source_order = source_order)
    }
    family_table <- family_table[order(family_table$n_variables, label_order, method = "radix"), , drop = FALSE]
    redundant <- vapply(seq_len(nrow(family_table)), function(candidate_index) { # Check whether this option contains a known smaller member set
      return(any(vapply(family_table$variables, function(other_variables) { # Compare another option with the current candidate membership
        return(length(other_variables) < family_table$n_variables[candidate_index] && all(other_variables %in% family_table$variables[[candidate_index]]))
      }, logical(1))))
    }, logical(1))
  } else {
    redundant <- logical()
  }
  supersets <- family_table[redundant, , drop = FALSE]
  family_table <- family_table[!redundant, , drop = FALSE]
  family_table$index <- seq_len(nrow(family_table))
  family_table$id <- if (nrow(family_table)) {
    paste0(effect, "-set-", family_table$index)
  } else {
    character()
  }
  rownames(family_table) <- rownames(supersets) <- NULL
  return(list(table = family_table, supersets = supersets))
}


# Directed reachability that stops at any currently displayed temporal vertex.
# This identifies the first visible descendants on each branch, avoiding serial
# chains U -> A -> B being described as a fork between A and B.
.gd_visible_frontiers <- function(reference, retained) { # Find the first displayed descendants along each omitted node branch
  edges <- dag_edges(reference)
  edges <- edges[edges$e == "->", , drop = FALSE]
  children <- split(edges$w, edges$v)
  frontier_cache <- new.env(parent = emptyenv())
  frontier <- function(node) { # Find and cache the first retained descendants of this vertex
    if (node %in% retained) {
      return(node)
    }
    if (exists(node, envir = frontier_cache, inherits = FALSE)) {
      return(get(node, frontier_cache))
    }
    next_nodes <- children[[node]]
    visible_descendants <- unique(unlist(lapply(next_nodes, frontier), use.names = FALSE))
    if (is.null(visible_descendants)) {
      visible_descendants <- character()
    }
    assign(node, visible_descendants, envir = frontier_cache)
    return(visible_descendants)
  }
  return(stats::setNames(lapply(setdiff(names(reference), retained), frontier),
    setdiff(names(reference), retained)))
}


.gd_target_nodes <- function(dag, state, role, node_name) { # Resolve the selected displayed roles to retained temporal nodes
  if (!is.null(role) && !is.null(node_name)) {
    stop("Choose role or node_name, not both.", call. = FALSE)
  }
  if (is.null(role) && is.null(node_name)) {
    role <- "outcome"
  }
  if (!is.null(role)) {
    role <- .gd_names(role, "role")
    allowed <- c("treatment", "outcome", "mediator", "confounder", "competing_cause",
      "instrument", "collider", "latent", "observed")
    if (!length(role) || length(setdiff(role, allowed))) {
      stop("Unknown or empty role. Use ", paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    members <- lapply(role, function(selected_role) { # Read the vertices carrying this requested causal role
      return(switch(selected_role,
        treatment = dagitty::exposures(dag), outcome = dagitty::outcomes(dag),
        latent = dagitty::latents(dag), observed = setdiff(names(dag), dagitty::latents(dag)),
        get_roles(dag, multiple_roles = TRUE)[[selected_role]]))
    })
    targets <- unique(unlist(members, use.names = FALSE))
    targets <- targets[!is.na(targets)]
  } else {
    node_name <- .gd_names(node_name, "node_name")
    if (!length(node_name)) {
      stop("Supply node_name.", call. = FALSE)
    }
    targets <- .gd_match_reference(dag, node_name, state$controls,
      .cd_record(state$original_query))
  }
  if (!length(targets)) {
    return(character())
  }
  # Map display concepts back to their retained temporal vertices.
  node_ids <- .gd_node_ids(dag)
  exact_nodes <- targets[targets != .strip_slice(targets)]
  quantity_ids <- unique(unname(node_ids[setdiff(targets, exact_nodes)]))
  temporal_ids <- .gd_node_ids(state$analysis)
  return(intersect(state$retained, union(exact_nodes, names(temporal_ids)[temporal_ids %in% quantity_ids])))
}


# Editing a recovered graph always starts from its saved temporal reference.
# No new causal facts are inferred from a display edge.

.gd_match_reference <- function(reference, requested, controls, record = NULL) { # Resolve requested names against the saved query and its aliases
  if (!is.null(attr(reference, "query_record", exact = TRUE))) {
    return(.cd_match_nodes(reference, requested))
  }
  requested <- .gd_names(requested, "node_name")
  quantity_ids <- .gd_node_ids(reference)
  resolved_nodes <- character()
  for (name in requested) {
    if (name %in% names(reference)) {
      resolved_nodes <- union(resolved_nodes, name)
      next
    }
    matching <- names(quantity_ids)[quantity_ids == name | .strip_slice(names(quantity_ids)) == name]
    if (!length(matching)) {
      canonical <- if (!is.null(record)) {
        resolved <- .resolve_columns(name, .build_alias_index(record$aliases),
          controls$disambiguate, character())
        if (length(resolved$ambiguous)) {
          stop("Ambiguous node name: ", name, call. = FALSE)
        }
        unname(resolved$resolved)
      } else {
        .gd_resolve_focus(name, controls$water_model,
          controls$disambiguate, strict = FALSE)
      }
      matching <- names(quantity_ids)[quantity_ids %in% canonical]
      node_map <- attr(reference, "cf_node_map", exact = TRUE)
      if (is.data.frame(node_map)) {
        matching <- union(matching,
          node_map$node[node_map$member_id %in% canonical | node_map$member_id == name])
      }
    }
    if (!length(matching)) {
      stop("Quantity absent from the saved query: ", name, call. = FALSE)
    }
    resolved_nodes <- union(resolved_nodes, matching)
  }
  return(resolved_nodes)
}


.gd_restore_reference <- function(state, requested) { # Recover requested quantities using the saved policy and measurements
  if (!is.null(state$construction)) {
    # Reuse the constructor's saved policy and explicit measurement assumptions.
    # Restore commands never query a possibly changed package node table.
    recovery_receipt <- state$build$recovery
    prepared <- .cd_prepare(state$original_query,
      state$construction$arguments$latents, union(recovery_receipt$requested_ids, requested))
    reference <- prepared$analysis
    known_nodes <- intersect(names(reference), names(state$analysis))
    latent <- union(setdiff(dagitty::latents(reference), known_nodes),
      intersect(dagitty::latents(state$analysis), known_nodes))
    reference <- .cd_set_latents(reference, latent)
    state$analysis <- .gd_canonical(reference)
    state$build$recovery$requested_ids <- unique(unname(
      .gd_node_ids(prepared$reference)[prepared$required]))
    state$build$recovery$ledger <- prepared$ledger
    return(state)
  }
  # Explicit requests can recover yes:user concepts, but not no:reason concepts.
  ids <- unique(unname(.gd_node_ids(state$original_query)[requested]))
  recovery_receipt <- state$build$recovery
  recovered <- .gd_recover(state$original_query,
    union(recovery_receipt$requested_ids, ids), state$controls$water_model,
    state$controls$disambiguate, apply_policy = recovery_receipt$policy_applied)
  reference <- recovered$graph
  known_nodes <- intersect(names(reference), names(state$analysis))
  availability_mode <- state$build$availability$mode
  latent <- if (availability_mode == "all") {
    character()
  } else if (availability_mode == "columns") {
    union(intersect(dagitty::latents(state$analysis), known_nodes),
      setdiff(names(reference), names(state$analysis)))
  } else {
    dagitty::latents(reference)
  }
  # Keep statuses from the active reference for quantities already present.
  latent <- union(setdiff(latent, known_nodes), intersect(dagitty::latents(state$analysis), known_nodes))
  reference <- .cd_set_latents(reference, latent)
  state$analysis <- .gd_canonical(reference)
  state$build$recovery <- recovered
  return(state)
}


.gd_restore_build <- function(dag, state, nodes, confounders = FALSE) { # Rebuild a restored projection and retain its checked temporal source
  state <- .gd_restore_reference(state, nodes)
  # Author's measurement rule: a restored quantity is unmeasured unless its
  # name was supplied in columns, irrespective of the earlier planning
  # scenario used to enumerate process views.
  if (confounders) {
    measured_ids <- state$build$observed_ids
    source_ids <- .gd_node_ids(state$analysis)
    unmeasured <- nodes[!source_ids[nodes] %in% measured_ids]
    measured <- setdiff(nodes, unmeasured)
    latent <- union(setdiff(dagitty::latents(state$analysis), measured), unmeasured)
    state$analysis <- .cd_set_latents(state$analysis, latent)
  }
  retained <- union(state$retained, nodes)
  missing_nodes <- setdiff(retained, names(state$analysis))
  if (length(missing_nodes)) {
    stop("These requested quantities were folded by the recovery policy: ",
      paste(missing_nodes, collapse = ", "), call. = FALSE)
  }
  omitted_nodes <- setdiff(names(state$analysis), retained)
  # A receipt constructed as observed folds the unmarked analysis, as its
  # construction did, and marks the result afterwards.
  as_observed <- identical(state$construction$unmeasured$semantics, "observed-as-if-v1")
  analysis_graph <- if (as_observed) {
    .cd_unmarked(state$analysis)
  } else {
    state$analysis
  }
  temporal <- if (length(omitted_nodes)) {
    collapse_nodes(analysis_graph, omitted_nodes)
  } else {
    analysis_graph
  }
  temporal <- .gd_canonical(temporal)
  # Every adjustment set the edited display would report is validated again on
  # the selection reference (the prepared analysis graph with its latent
  # statuses, unmarked for receipts constructed as observed), as construction
  # and views do.
  summary <- if (state$controls$simplify) {
    .gd_summary(temporal, keep = retained,
      source_check = .cd_summary_source_check(.gd_canonical(analysis_graph)))
  } else {
    NULL
  }
  if (!is.null(summary)) {
    .gd_segment_tie_notice(summary)
  }
  graph <- if (is.null(summary)) {
    temporal
  } else {
    summary$graph
  }
  if (as_observed) {
    latent_nodes <- dagitty::latents(state$analysis)
    marked_graph <- .cd_mark_unmeasured(if (is.null(summary)) {
      list(graph = graph, temporal = temporal)
    } else {
      summary
    }, latent_nodes)
    graph <- marked_graph$summary$graph
    if (!is.null(summary)) {
      summary <- marked_graph$summary
    }
    projected_latent <- intersect(names(temporal), latent_nodes)
    if (length(projected_latent)) {
      temporal <- .cd_set_latents(temporal, projected_latent)
    }
  }
  state$temporal <- temporal
  state$retained <- names(temporal)
  state$summary <- summary
  state$focus <- unique(unname(.gd_node_ids(graph)))
  state$edit_count <- if (is.null(state$edit_count)) {
    1L
  } else {
    state$edit_count + 1L
  }
  return(.gd_attach(graph, state))
}
