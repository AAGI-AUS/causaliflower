################################################################################
# get_dag.R


# Guided DAG construction, temporal summaries and their saved provenance
# Original filename section markers below record source history


#####
# Guided query and construction
# Historical source: get_dag.R
#####


#####
# get_dag.R
#####
#' Construct a DAG through a guided query
#'
#' Calls `query_graph()` and then `construct_dag()`. Collects missing query inputs,
#' resolves deferred choices interactively and saves the resolved call and answers.
#' Use the underlying functions to work with an existing or deliberately
#' restricted query.
#'
#' @param crop,treatments,outcomes Crop and query endpoints.
#' @param randomised_treatment Required experimental design assumption, `TRUE` or `FALSE`.
#' @param latents Quantities or exact occasions declared unmeasured. Their treatment
#'   depends on the search mode described in `construct_dag()`. `NULL` treats all
#'   quantities as observed. The original query's tags remain saved.
#' @param label Naming option passed to `query_graph()`. `NULL` retains quantity keys.
#' @param query_graph Further query inputs: `columns`, `time_unroll` and `params`.
#' @param construct_dag Construction inputs: `direct_select`, `total_select`,
#'   `search` and `steps`. Selectors default to `"auto"`. `NULL` defers a choice.
#'   `keep_nodes` is available when supplied explicitly.
#' @param dry_run Return a proposed construction for inspection.
#' @param quiet Suppress routine progress when `TRUE`. Saved checks remain available.
#'
#' @return The graph, deferred choices or a construction preview, with its saved
#' query and construction records.
#' @export
get_dag <- function(crop, treatments, outcomes, randomised_treatment,
  latents = NULL, label = NULL, query_graph = list(),
  construct_dag = list(), dry_run = FALSE, quiet = NULL) { # Recover the guided query and construct its selected DAG
  return(.gd_get_dag(crop, treatments, outcomes, randomised_treatment, latents, label,
    query_graph, construct_dag, dry_run, quiet, original_call = match.call()))
}


# Private continuation boundary. The public query contract and the new query's
# receipts are unchanged: only a certified source enumeration may be reused.
.gd_get_dag <- function(crop, treatments, outcomes, randomised_treatment,
  latents = NULL, label = NULL, query_graph = list(),
  construct_dag = list(), dry_run = FALSE, quiet = NULL,
  previous = NULL, original_call = NULL) { # Collect query inputs and preserve the construction call and responses
  if (is.null(original_call)) {
    original_call <- match.call()
    original_call$previous <- original_call$original_call <- NULL
    original_call[[1L]] <- as.name("get_dag")
  }
  supplied <- intersect(names(formals(get_dag)), names(original_call))
  required <- c("crop", "treatments", "outcomes", "randomised_treatment")
  missing_inputs <- setdiff(required, supplied)
  if (length(missing_inputs) && !.gd_interactive()) {
    stop("get_dag requires: ", paste(missing_inputs, collapse = ", "),
      ". Supply these inputs or call interactively.", call. = FALSE)
  }
  responses <- new.env(parent = emptyenv())
  responses$rows <- data.frame(field = character(), prompt = character(), response = character())
  values <- list()
  if ("crop" %in% supplied) {
    values$crop <- .gd_query_value("crop", crop)
  } else {
    values$crop <- .gd_prompt_value("crop", responses)
  }
  if ("treatments" %in% supplied) {
    values$treatments <- .gd_query_value("treatments", treatments)
  } else {
    values$treatments <- .gd_prompt_value("treatments", responses)
  }
  if ("outcomes" %in% supplied) {
    values$outcomes <- .gd_query_value("outcomes", outcomes)
  } else {
    values$outcomes <- .gd_prompt_value("outcomes", responses)
  }
  if ("randomised_treatment" %in% supplied) {
    values$randomised_treatment <- .gd_query_value("randomised_treatment", randomised_treatment)
  } else {
    values$randomised_treatment <- .gd_prompt_value("randomised_treatment", responses)
  }
  .gd_flag(dry_run, "dry_run")
  .gd_flag(quiet, "quiet", nullable = TRUE)
  progress <- .gd_progress(quiet, "Query")
  query_engine <- get("query_graph", envir = asNamespace("causaliflower"))
  engine <- .gd_named_inputs(query_graph,
    setdiff(names(formals(query_engine)), c(required, "label", "quiet")), "query_graph")
  choices <- .gd_named_inputs(construct_dag,
    c("direct_select", "total_select", "search", "steps"), "construct_dag")
  selectors <- list(direct_select = "auto", total_select = "auto")
  for (key in intersect(names(selectors), names(choices))) {
    selectors[key] <- choices[key]
  }
  if (!"params" %in% names(engine) || is.null(engine$params)) {
    engine$params <- query_params()
  } else {
    engine$params <- do.call(query_params,
      .gd_named_inputs(engine$params, names(formals(query_params)), "query_graph$params"))
  }
  if (!is.infinite(engine$params$query_depth) || !isTRUE(engine$params$show_hidden_nodes)) {
    stop("get_dag uses the complete query: query_depth = Inf and show_hidden_nodes = TRUE. ",
      "For a deliberately restricted graph, call query_graph() then construct_dag().", call. = FALSE)
  }
  .gd_notify(progress, "Building the source query")
  graph <- do.call("query_graph", c(values, list(label = label, quiet = quiet), engine),
    envir = asNamespace("causaliflower"))
  .gd_notify(progress, paste0("Source query ready (", length(names(graph)), " nodes)"))
  construct_arguments <- c(list(graph = graph, latents = latents), selectors,
    list(dry_run = dry_run, quiet = quiet))
  if ("search" %in% names(choices)) {
    construct_arguments$search <- choices$search
  }
  if ("steps" %in% names(choices)) {
    construct_arguments$steps <- choices$steps
  }
  constructor <- if (is.null(previous)) {
    "construct_dag"
  } else {
    ".cd_reconstruct"
  }
  if (!is.null(previous)) {
    construct_arguments$previous <- previous
  }
  result <- do.call(constructor, construct_arguments, envir = asNamespace("causaliflower"))
  if (.gd_interactive()) {
    for (attempt in seq_len(3L)) {
      if (!inherits(result, "cf_dag_options")) {
        break
      }
      next_selectors <- .gd_ui_options(result, selectors, responses)
      if (identical(next_selectors, selectors)) {
        break
      }
      selectors <- next_selectors
      for (key in names(selectors)) {
        construct_arguments[key] <- selectors[key]
      }
      result <- do.call(constructor, construct_arguments, envir = asNamespace("causaliflower"))
    }
  }
  effective_choices <- selectors
  completed <- if (inherits(result, "cf_dag")) {
    .gd_state(result)$selected
  } else if (inherits(result, "cf_dag_plan")) {
    result$selected
  } else {
    NULL
  }
  method <- if (inherits(result, "cf_dag")) {
    .gd_state(result)$construction$search_method
  } else if (inherits(result, "cf_dag_plan")) {
    result$options$search_method
  } else if (inherits(result, "cf_dag_options")) {
    result$search_method
  }
  witness <- identical(method, "source-witness-v1")
  if (!witness && !is.null(completed)) {
    for (family in names(completed)) {
      if (!is.null(completed[[family]])) {
        effective_choices[paste0(family, "_select")] <- completed[family]
      }
    }
  }
  # The saved call names the route that actually produced this result.
  modes <- c("source-witness-v1" = "witness", "quantity-ranked-v1" = "auto",
    "complete-source-sets-v2" = "complete")
  if (isTRUE(method %in% names(modes))) {
    effective_choices$search <- unname(modes[method])
  } else if ("search" %in% names(choices)) {
    effective_choices$search <- choices$search
  }
  # A ranked result records the step budget it used, so a replay reproduces its catalogue.
  budget <- if (inherits(result, "cf_dag")) {
    .gd_state(result)$construction$arguments$steps
  } else if (inherits(result, "cf_dag_plan")) {
    result$options$arguments$steps
  } else if (inherits(result, "cf_dag_options")) {
    result$arguments$steps
  }
  if (!is.null(budget)) {
    effective_choices$steps <- budget
  } else if ("steps" %in% names(choices)) {
    effective_choices$steps <- choices$steps
  }
  query_record <- attr(graph, "query_record", exact = TRUE)
  if (!is.null(query_record$arguments)) {
    engine <- query_record$arguments[
      intersect(c("columns", "time_unroll", "params"), names(query_record$arguments))]
  }
  resolved <- c(values, list(latents = latents, label = label, query_graph = engine,
    construct_dag = effective_choices, dry_run = dry_run, quiet = quiet))
  record <- list(original_call = original_call,
    resolved_call = as.call(c(list(as.name("get_dag")), resolved)),
    values = resolved, supplied = supplied, responses = responses$rows,
    requested_selectors = choices, requested_query_graph = query_graph,
    query_record = query_record)
  if (inherits(result, "cf_dag")) {
    state <- .gd_state(result)
    state$input_record <- record
    result <- .gd_attach(result, state)
  } else {
    attr(result, "get_dag_record") <- record
  }
  if (nrow(responses$rows)) {
    message("Your responses and resolved call were saved to the output object. Use dag_call().")
  }
  return(result)
}


# Progress is transient console output, never part of the graph or its receipts.
# A complete separator count is unknown until enumeration finishes, so no ETA is
# inferred. Throttled loop updates inspect the existing count and elapsed time only.
.gd_progress <- function(quiet, operation = "Construction") { # Create a throttled progress reporter without changing saved results
  if (isTRUE(quiet)) {
    return(NULL)
  }
  started <- proc.time()[["elapsed"]]
  last <- started
  return(function(stage, force = TRUE) { # Report elapsed time when an update is due
    now <- proc.time()[["elapsed"]]
    if (!force && now - last < 2) {
      return(invisible(NULL))
    }
    message(operation, " [", format(round(now - started, 1), nsmall = 1), "s]: ", stage)
    flush.console()
    last <<- now
    return(invisible(NULL))
  })
}


.gd_notify <- function(progress, stage, force = TRUE) { # Send a progress update when a reporter is present
  if (!is.null(progress)) {
    progress(stage, force = force)
  }
  return(invisible(NULL))
}


#####
# construct_dag.R
#####
#' Construct a compact DAG from a supplied query
#'
#' Uses the supplied `query_graph()` result as its reference graph.
#'
#' @param graph The supplied query result.
#' @param latents Quantities or exact occasions declared unmeasured. With
#'   `search = "auto"`, they do not change the ranked search or automatic
#'   selection. Construction treats the graph as fully observed, then marks
#'   the declared quantities on the finished graph. With `"complete"` or
#'   `"witness"`, declarations are applied before searching. `NULL` treats
#'   all quantities as observed.
#' @param direct_select,total_select `"auto"`, `NULL` to list choices, or a
#'   group and option index.
#' @param keep_nodes Additional quantities or exact occasions to retain.
#' @param dry_run Return the proposed graph, construction plan and checks.
#' @param quiet `NULL` prints concise progress. `TRUE` suppresses it. `FALSE`
#'   also prints detailed construction records.
#' @param search `"auto"` ranks adjustment sets that are minimal under set
#'   inclusion by their number of distinct source quantities. Repeated time
#'   occasions count once. `"complete"` enumerates all source adjustment sets
#'   minimal under set inclusion, retaining exact occasions, and can be slow.
#'   `"witness"` retains one qualifying set minimal
#'   under set inclusion per effect and requires automatic selectors.
#' @param steps Step budget for the ranked search (`search = "auto"`) for each
#'   effect. It is not accepted for the other search modes. `NULL` uses
#'   `getOption("causaliflower.ranked_steps", 20000)`. The saved call retains this budget.
#'
#' @details
#' The options table's `unmeasured` field describes declared occasions used by
#' each option. `dag_checks()` describes the marked nodes used in the displayed
#' graph. A displayed node is marked if any occasion it represents is declared
#' unmeasured, so these records have different scopes.
#'
#' When no total adjustment option exists on a fully observed directed reference,
#' exact treatment parents are also retained to preserve joint identification
#' through the g-formula. This retains temporal occurrences and does not certify
#' a graph combined across time. Declared unmeasured treatment parents are reported.
#'
#' @return A `cf_dag`, `cf_dag_options` or `cf_dag_plan`, with the source and checks saved.
#' @export
construct_dag <- function(graph, latents = NULL, direct_select = "auto", total_select = "auto",
  keep_nodes = NULL, dry_run = FALSE, quiet = NULL,
  search = c("auto", "complete", "witness"), steps = NULL) { # Construct the requested DAG and retain its reference and checks
  return(.cd_construct(graph, latents, direct_select, total_select, keep_nodes, dry_run,
    quiet, original_call = match.call(), search_mode = match.arg(search), steps = steps))
}


.cd_reconstruct <- function(graph, previous, latents = NULL, direct_select = "auto",
  total_select = "auto", keep_nodes = NULL,
  dry_run = FALSE, quiet = NULL, search = c("auto", "complete", "witness"),
  steps = NULL) { # Reconstruct from a previous result while preserving the public call
  original_call <- match.call()
  original_call$previous <- NULL
  original_call[[1L]] <- as.name("construct_dag")
  return(.cd_construct(graph, latents, direct_select, total_select, keep_nodes, dry_run,
    quiet, previous = previous, original_call = original_call, search_mode = match.arg(search),
    steps = steps))
}


.cd_construct <- function(graph, latents, direct_select, total_select, keep_nodes,
  dry_run, quiet, previous = NULL, original_call, search_mode = "complete",
  steps = NULL) { # Prepare, search, select and verify the supplied query graph
  .gd_flag(dry_run, "dry_run")
  if (!is.null(quiet)) {
    .gd_flag(quiet, "quiet")
  }
  if (!is.null(steps) && !identical(search_mode, "auto")) {
    stop("steps applies to the ranked search (search = 'auto') only.", call. = FALSE)
  }
  # One step budget serves the ranked search, its identity and the saved arguments.
  budget <- if (identical(search_mode, "auto")) {
    .cd_ranked_steps(steps)
  }
  max_sets <- Inf
  progress <- .gd_progress(quiet, "Construction")
  .gd_notify(progress, "Preparing the reference graph")
  prepared <- .cd_prepare(graph, latents, keep_nodes)
  .gd_notify(progress, paste0("Reference ready (", length(names(prepared$selection)), " nodes)"))
  witness <- identical(search_mode, "witness")
  if (witness && !(identical(direct_select, "auto") && identical(total_select, "auto"))) {
    stop("search = 'witness' requires automatic direct and total selectors.", call. = FALSE)
  }
  ranked <- identical(search_mode, "auto")
  # The ranked route constructs exactly as for a fully observed graph (the same
  # cached preparation as latents = NULL). Declared unmeasured variables are
  # marked on the finished DAG. Witness and complete search the declared graph.
  observed <- if (ranked && length(prepared$latent_nodes)) {
    .cd_prepare(graph, NULL, keep_nodes)
  } else {
    prepared
  }
  search_identity <- .cd_search_identity(observed)
  if (witness && !is.null(search_identity)) {
    search_identity$method <- "source-witness-v1"
    search_identity$witness_adapter <- paste(deparse(body(.cd_witness_search)), collapse = "\n")
  }
  if (ranked && !is.null(search_identity)) {
    search_identity$method <- "quantity-ranked-v1"
    search_identity$ranked_adapter <- lapply(list(.cd_ranked_search, .cd_ranked_enumerate,
      .cd_ranked_path, .cd_ranked_total_region, .cd_ranked_direct_region, .cd_ranked_graph,
      .cd_ranked_reach, .cd_ranked_without, .cd_ranked_moral, .cd_ranked_separator,
      .cd_ranked_labels), function(adapter) { # Record the implementation used by this ranked adapter
        return(deparse(body(adapter)))
      })
    search_identity$ranked_steps <- budget
  }
  if (witness) {
    search <- .cd_reuse_witness(previous, prepared, search_identity, progress)
    if (is.null(search)) {
      search <- stats::setNames(lapply(c("direct", "total"), function(effect) { # Find one source witness for this effect family
        return(.cd_witness_search(prepared, effect, progress))
      }), c("direct", "total"))
    }
    refined <- list(search = search, refinement = list(method = "source-witness-v1",
      source_signature = .gd_signature(prepared$selection),
      source_sets = lapply(search, `[[`, "sets"), global_minimum_certified = FALSE,
      passes = list(), stop_reason = "One source witness per effect. Refinement and full source enumeration were not requested."))
  } else if (ranked) {
    search <- .cd_reuse_ranked(previous, observed, search_identity, progress)
    if (is.null(search)) {
      search <- stats::setNames(lapply(c("direct", "total"), function(effect) { # Rank adjustment options for this effect family
        return(.cd_ranked_search(observed, effect, progress, steps = budget))
      }), c("direct", "total"))
    }
    # Each option lists the declared unmeasured nodes (exact occasions or groups)
    # its realised set uses. It is never removed.
    declared_members <- unique(prepared$map$member_id[prepared$map$node %in% prepared$latent_nodes])
    for (effect in names(search)) {
      search[[effect]]$table$unmeasured <- lapply(search[[effect]]$table$raw_sets, function(sets) { # List declared unmeasured nodes used by these realised sets
        return(sort(unique(intersect(as.character(unlist(sets, use.names = FALSE)), prepared$latent_nodes)), method = "radix"))
      })
    }
    refined <- list(search = search, refinement = list(method = "quantity-ranked-v1",
      source_signature = .gd_signature(observed$selection),
      source_sets = lapply(search, `[[`, "sets"),
      global_minimum_certified = all(vapply(search, function(x) { # Check whether this effect search is complete
        return(isTRUE(x$complete))
      }, logical(1))),
      passes = list(), stop_reason = "Ranked by distinct source members. Each family records how far its ranking is certified."))
  } else {
    search <- .cd_reuse_search(previous, prepared, search_identity, progress)
    if (is.null(search)) {
      search <- stats::setNames(lapply(c("direct", "total"), function(effect) { # Enumerate adjustment options for this effect family
        return(.cd_search(prepared, effect, max_sets, progress))
      }), c("direct", "total"))
    }
    refined <- .cd_refine_search(prepared, search, max_sets)
  }
  search <- refined$search
  direct <- .cd_select(search$direct$table, direct_select, "direct")
  retained_ids <- unique(observed$map$member_id[observed$map$node %in% observed$required])
  if (length(direct$index)) {
    retained_ids <- union(retained_ids, search$direct$table$resolved[[direct$index]])
  }
  total <- .cd_select(search$total$table, total_select, "total", retained_ids,
    pending = identical(direct$status, "deferred"))
  choices <- list(direct = direct, total = total)
  if (ranked) {
    for (effect in names(choices)) {
      if (choices[[effect]]$status == "no_sets" && !isTRUE(search[[effect]]$complete)) {
        choices[[effect]]$reason <- "No option was found within the ranked search's step budget. This is not evidence that no set exists."
      }
    }
  }
  if (witness) {
    for (effect in names(choices)) {
      if (choices[[effect]]$status == "selected") {
        choices[[effect]]$reason <- "Deterministic source witness minimal under set inclusion. Minimum quantity cost and other alternatives were not searched."
      }
    }
  }
  arguments <- list(graph = graph, latents = latents, direct_select = direct_select,
    total_select = total_select, keep_nodes = keep_nodes, dry_run = dry_run, quiet = quiet,
    search = if (witness) {
      "witness"
    } else if (ranked) {
      "auto"
    } else {
      "complete"
    })
  if (ranked) {
    arguments$steps <- budget
  }
  options <- structure(list(graph = graph, direct = search$direct$table, total = search$total$table,
    selected = lapply(choices, `[[`, "selection"), status = vapply(choices, `[[`, character(1), "status"),
    choices = choices, search = search, refinement = refined$refinement,
    search_method = if (witness) {
      "source-witness-v1"
    } else if (ranked) {
      "quantity-ranked-v1"
    } else {
      "complete-source-sets-v2"
    },
    enumeration = if (witness) {
      "witness"
    } else if (ranked) {
      "ranked"
    } else {
      "complete"
    },
    search_identity = search_identity,
    arguments = arguments, limitations = prepared$limitations,
    scope = "The supplied emitted graph, with explicit observation and saved constraints on recovery eligibility.",
    original_call = original_call), class = "cf_dag_options")
  if (ranked) {
    options$unmeasured <- list(semantics = "observed-as-if-v1", declared = latents,
      nodes = prepared$latent_nodes, members = declared_members, display_nodes = character(),
      partial = list(), selected = list(direct = character(), total = character()),
      scope = "Declared unmeasured variables never change the search, the catalogue or the selection. They are marked on the finished DAG, and identification and the checks report which members they affect.")
  }
  for (effect in names(search)) {
    if (ranked && !isTRUE(search[[effect]]$complete)) {
      certificate <- search[[effect]]$ranked
      if (identical(certificate$certificate, "validation")) {
        warning(effect, " options are not certified: ", certificate$dropped,
          " option(s) failed validation on the selection reference and were dropped.", call. = FALSE)
      } else if (!isTRUE(certificate$complete_up_to >= 0)) {
        warning(effect, " options are ranked but not certified: the step budget stopped the search before any member count was completed.", call. = FALSE)
      } else {
        warning(effect, " options are ranked and certified complete up to ",
          certificate$complete_up_to, " distinct source members. Larger options may exist",
          if (isTRUE(search[[effect]]$omitted_equal_smallest)) {
            ". Alternatives with the same member count as the smallest shown may be missing"
          } else {
            ""
          },
          " (step budget ", certificate$step_budget, ").", call. = FALSE)
      }
    } else if (!witness && !isTRUE(search[[effect]]$complete)) {
      search_warning <- paste0(effect, " initial search stopped after ", max_sets,
        " sets with exact occasions, before all sets were examined. ",
        "Smaller sets or tied alternatives may remain after removing repeated time names.")
      if (length(search[[effect]]$refinement_sets)) {
        search_warning <- paste0(search_warning, " Additional smaller sets were admitted after verification on the fixed selection reference.")
      }
      if (isTRUE(search[[effect]]$omitted_equal_smallest)) {
        search_warning <- paste0(search_warning, " The next set has the same count of distinct nodes as the smallest group shown.")
      }
      warning(search_warning, call. = FALSE)
    }
  }
  if (length(prepared$limitations)) {
    warning(paste(prepared$limitations, collapse = " "), call. = FALSE)
  }
  if (any(options$status == "deferred")) {
    .gd_notify(progress, "Adjustment searches complete. Waiting for an option selection")
    if (isFALSE(quiet)) {
      print(options)
    }
    return(options)
  }
  chosen_sets <- lapply(names(search), function(effect) { # Read the selected memberships for this effect family
    return(.cd_selected_sets(search[[effect]], choices[[effect]]))
  })
  .gd_notify(progress, "Constructing the selected graph")
  names(chosen_sets) <- names(search)
  selected_nodes <- unique(unlist(chosen_sets, recursive = TRUE, use.names = FALSE))
  selected_ids <- unique(unname(.gd_node_ids(observed$reference)[selected_nodes]))
  # Recover selected quantities across occasions, preserving every exact chosen membership.
  quantity_nodes <- names(observed$analysis)[.gd_node_ids(observed$analysis) %in% selected_ids]
  retained <- union(quantity_nodes, prepared$required)
  identification <- .cd_identification_retention(observed$analysis, choices$total)
  identification$added_nodes <- setdiff(identification$parents, retained)
  retained <- union(retained, identification$parents)
  omit <- setdiff(names(observed$analysis), retained)
  temporal <- if (length(omit)) {
    collapse_nodes(observed$analysis, omit)
  } else {
    observed$analysis
  }
  temporal <- .cd_working_graph(temporal, prepared$map[prepared$map$node %in% names(temporal), , drop = FALSE],
    source_order = attr(observed$analysis, "cf_source_order", exact = TRUE))
  pins <- union(union(selected_nodes, prepared$required), identification$parents)
  summary <- .gd_summary(temporal, keep = pins,
    preserve_occasions = isTRUE(identification$applied),
    source_check = .cd_summary_source_check(observed$selection))
  # On the ranked route this graph is unmarked until .cd_mark_unmeasured below.
  final <- summary$graph
  if (isTRUE(identification$applied)) {
    edges <- dag_edges(final)
    treatment <- dagitty::exposures(observed$analysis)
    required <- union(treatment, union(dagitty::outcomes(observed$analysis), identification$parents))
    verified <- all(required %in% names(final)) &&
    !length(dagitty::latents(final)) &&
    !any(edges$e == "<->" & (edges$v %in% treatment | edges$w %in% treatment))
    if (!verified) {
      stop("Retention of exact treatment parents failed its structural invariant.", call. = FALSE)
    }
    identification$status <- "preserved by retaining exact treatment parents"
    # Later restore_nodes() operations must not silently undo this time policy.
    prepared$controls$simplify <- FALSE
  }
  if (ranked) {
    identification <- .cd_unmeasured_parents(identification, prepared$latent_nodes)
  }
  options$identification <- identification
  .gd_notify(progress, "Checking adjustment alternatives on the constructed graph")
  checks <- .cd_checks(observed, final, search, choices, max_sets, progress, verify_all = witness)
  # Temporal fallback returns this exact same graph, so reuse its completed checks.
  temporal_checks <- if (identical(final, summary$temporal)) {
    checks
  } else {
    .cd_checks(observed, summary$temporal, search, choices, max_sets, progress, verify_all = witness)
  }
  if (witness) {
    accepted <- function(x) { # Check that source memberships and output alternatives remain certified
      return(all(vapply(x, function(check) { # Check that source memberships and output alternatives remain certified
        return(check$status %in% c("preserved", "no eligible reference set") &&
          isTRUE(check$output_reference_complete))
      }, logical(1))))
    }
    if (!accepted(checks) && !identical(final, temporal)) {
      attempted <- list(summary = summary, checks = checks)
      summary <- .gd_summary(temporal, keep = pins, preserve_occasions = TRUE)
      summary$reason <- "Compact alternatives could not be certified with exact source temporal memberships. Exact temporal occurrences were retained."
      summary$attempted <- attempted
      final <- summary$graph
      checks <- .cd_checks(prepared, final, search, choices, max_sets, progress, verify_all = TRUE)
      temporal_checks <- checks
      prepared$controls$simplify <- FALSE
    }
    if (!accepted(checks) || !accepted(temporal_checks)) {
      stop("Witness construction could not certify its selected memberships and all Output alternatives on the unchanged source. No complete source search was attempted. Inspect an explicit search = 'complete' construction separately.", call. = FALSE)
    }
  }
  # Ranked: mark the declared unmeasured variables on the finished construction.
  # The checks above ran on it unmarked and now say which members they need.
  shown <- if (ranked) {
    .cd_mark_unmeasured(summary, prepared$latent_nodes)
  }
  if (ranked) {
    checks <- .cd_annotate_unmeasured(checks, shown$display)
    temporal_checks <- .cd_annotate_unmeasured(temporal_checks, shown$temporal)
  }
  .gd_segment_tie_notice(summary)
  requested_ledger <- data.frame(node = prepared$requested,
    in_direct = prepared$requested %in% unique(unlist(chosen_sets$direct)),
    in_total = prepared$requested %in% unique(unlist(chosen_sets$total)),
    endpoint = prepared$requested %in% c(dagitty::exposures(prepared$reference), dagitty::outcomes(prepared$reference)),
    stringsAsFactors = FALSE)
  request_notice <- data.frame(input = character(), origin = character(), reason = character(), stringsAsFactors = FALSE)
  notice_inputs <- list(columns = names(prepared$record$requests$columns$resolved), keep_nodes = keep_nodes)
  for (origin in names(notice_inputs)) {
    for (input in notice_inputs[[origin]]) {
      nodes <- .cd_match_nodes(graph, input, prepared$record)
      rows <- match(nodes, requested_ledger$node)
      rows <- rows[!is.na(rows)]
      if (!length(rows)) {
        next
      }
      represented <- requested_ledger$in_direct[rows] | requested_ledger$in_total[rows] | requested_ledger$endpoint[rows]
      if (!all(represented)) {
        next
      }
      reasons <- c(if (any(requested_ledger$in_direct[rows])) {
        "selected direct membership"
      },
        if (any(requested_ledger$in_total[rows])) {
          "selected total membership"
        },
        if (any(requested_ledger$endpoint[rows])) {
          "query endpoint"
        })
      rules <- prepared$ledger$recovered_node[match(nodes, prepared$ledger$node)]
      if (any(grepl("yes:user", rules, fixed = TRUE))) {
        reasons <- c(reasons, "eligibility enabled by the request")
      }
      request_notice <- rbind(request_notice, data.frame(input = input, origin = origin,
        reason = paste(reasons, collapse = ", "), stringsAsFactors = FALSE))
    }
  }
  recovery <- list(graph = prepared$analysis, ledger = prepared$ledger,
    requested_ids = unique(unname(.gd_node_ids(prepared$reference)[prepared$required])),
    policy_applied = prepared$policy_applied, coverage = "Saved query recovery policy and explicit requests.")
  availability <- list(mode = "explicit", latent_nodes = prepared$latent_nodes,
    ledger = data.frame(node = names(prepared$reference), query_latent = names(prepared$reference) %in% dagitty::latents(graph),
    analysis_latent = names(prepared$reference) %in% prepared$latent_nodes))
  if (ranked) {
    availability$display_latent <- shown$display
    availability$partial <- shown$partial
    options$unmeasured$display_nodes <- shown$display
    options$unmeasured$partial <- shown$partial
    options$unmeasured$selected <- lapply(stats::setNames(nm = names(choices)), function(effect) { # Read the unmeasured nodes used by the selected option
      return(if (length(choices[[effect]]$index)) {
        search[[effect]]$table$unmeasured[[choices[[effect]]$index]]
      } else {
        character()
      })
    })
  }
  resolved <- arguments
  for (effect in names(choices)) {
    if (!witness && identical(choices[[effect]]$status, "selected")) {
      resolved[[paste0(effect, "_select")]] <- choices[[effect]]$selection
    }
  }
  resolved$dry_run <- FALSE
  options$resolved_arguments <- resolved
  input_record <- list(original_call = original_call, resolved_call = as.call(c(list(as.name("construct_dag")), resolved)),
    values = resolved, responses = data.frame(field = character(), prompt = character(), response = character()))
  # The ranked route sets nothing before its search: its history starts from the
  # supplied query, records the observed construction and ends by marking.
  history <- c(if (ranked) {
    list(.gd_transition(graph, observed$analysis, "recover_nodes", list(requested = recovery$requested_ids), recovery[setdiff(names(recovery), "graph")]))
  } else {
    list(.gd_transition(graph, prepared$reference, "set_latents", list(latents = latents), availability),
      .gd_transition(prepared$reference, prepared$analysis, "recover_nodes", list(requested = recovery$requested_ids), recovery[setdiff(names(recovery), "graph")]))
  },
    list(.gd_transition(observed$analysis, observed$analysis, "refine_adjustment_search",
    details = refined$refinement),
    .gd_transition(observed$analysis, temporal, "select_direct_and_total", options$selected,
    list(memberships = chosen_sets, requested = requested_ledger,
    identification_retention = identification)),
    .gd_transition(temporal, summary$temporal, "trim_time_routes", details = list(ledger = summary$trim_ledger)),
    .gd_transition(summary$temporal, final, "summarise_time", details = summary[setdiff(names(summary), c("graph", "temporal"))])))
  if (ranked) {
    history <- c(history, list(.gd_transition(final, shown$summary$graph, "mark_unmeasured", list(latents = latents),
      list(nodes = shown$display, temporal_nodes = shown$temporal, partial = shown$partial,
      availability = lapply(checks, `[[`, "availability")))))
    summary <- shown$summary
    final <- summary$graph
  }
  state <- list(original_query = graph, analysis = prepared$analysis, temporal = summary$temporal,
    retained = names(summary$temporal), summary = summary, focus = selected_ids, controls = prepared$controls,
    build = list(candidates = list(direct = options$direct, total = options$total), search = search,
    recovery = recovery[setdiff(names(recovery), "graph")], availability = availability, observed_ids = prepared$columns),
    observed_ids = prepared$columns, history = history, undo = list(), construction = options,
    selected = options$selected, checks = checks, temporal_checks = temporal_checks,
    checks_signature = .gd_signature(final), selection_base_signature = .gd_signature(final),
    input_record = input_record, requested = requested_ledger, request_notice = request_notice)
  result <- .gd_attach(final, state)
  .gd_notify(progress, paste0("Complete (", length(names(final)), " nodes)"))
  # A ranked family stopped by its step budget is not evidence of absence.
  unfinished <- names(options$status)[options$status == "no_sets" &
    vapply(names(options$status), function(effect_name) { # Identify a family whose step budget ended before finding an option
      return(ranked && !isTRUE(search[[effect_name]]$complete))
    }, logical(1))]
  absent <- setdiff(names(options$status)[options$status == "no_sets"], unfinished)
  if (length(unfinished)) {
    warning("No ", paste(unfinished, collapse = "/"),
      " adjustment option was found within the ranked search's step budget. This is not evidence that no set exists.", call. = FALSE)
  }
  if (length(absent)) {
    warning("No eligible ", paste(absent, collapse = "/"),
      " adjustment set was found. Returning a diagnostic graph. This is not a valid empty set.", call. = FALSE)
  }
  if (any(vapply(checks, function(x) { # Identify a changed or incompletely verified adjustment family
    return(x$status %in% c("changed", "not fully checked"))
  }, logical(1)))) {
    warning("Selected adjustment memberships changed or were not fully verified in the final graph. Inspect dag_checks() and temporal_dag().", call. = FALSE)
  }
  if (ranked) {
    for (effect in names(options$unmeasured$selected)) {
      if (length(options$unmeasured$selected[[effect]])) {
        warning("The selected ", effect, " adjustment set relies on declared unmeasured variable(s) ",
          paste(options$unmeasured$selected[[effect]], collapse = ", "), ". They are marked on the DAG.", call. = FALSE)
      }
    }
  }
  if (isFALSE(quiet)) {
    cat("Causal query: ", paste(c(paste0(dagitty::exposures(graph), " (treatment)"),
      paste0(dagitty::outcomes(graph), " (outcome)")), collapse = ", "), "\n", sep = "")
    print(options)
    if (nrow(request_notice)) {
      for (i in seq_len(nrow(request_notice))) {
        message("Requested ", request_notice$input[i], " already retained: ", request_notice$reason[i], ".")
      }
    }
  }
  if (dry_run) {
    plan <- structure(list(graph = result, options = options, selected = options$selected, status = options$status,
      changes = .gd_delta(graph, result), checks = checks, temporal_checks = temporal_checks,
      requested = requested_ledger, request_notice = request_notice, arguments = arguments), class = "cf_dag_plan")
    if (isFALSE(quiet)) {
      print(plan)
    }
    return(plan)
  }
  return(result)
}


#####
# Input, construction and provenance internals
# Historical source: get_dag.R
#####


#####
# get_dag_inputs.R
#####
# Guided input collection never evaluates console answers as R code.
.gd_interactive <- function() { # Read whether input collection can use the current console
  return(interactive())
}


.gd_readline <- function(prompt) { # Read one console response without evaluating it as code
  return(readline(prompt))
}


.gd_named_inputs <- function(x, allowed, name) { # Validate a named list against the supported input names
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
    any(!nzchar(names(x))) || anyDuplicated(names(x))))) {
    stop(name, " must be a named list with unique names.", call. = FALSE)
  }
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    stop("Unsupported ", name, " input: ",
      paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  return(x)
}


.gd_query_value <- function(name, value) { # Validate one required causal query input
  if (name == "randomised_treatment") {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("randomised_treatment must be TRUE or FALSE.", call. = FALSE)
    }
  } else if (!is.character(value) || !length(value) || anyNA(value) ||
    any(!nzchar(trimws(value))) || (name == "crop" && length(value) != 1L)) {
    stop(name, " must contain ", if (name == "crop") {
      "one nonempty crop name."
    } else {
      "nonempty node names."
    }, call. = FALSE)
  }
  return(value)
}


.gd_answer <- function(name, answer) { # Parse a console response into a validated query value or selector
  if (!is.character(answer) || length(answer) != 1L || is.na(answer) ||
    !nzchar(trimws(answer))) {
    stop("An explicit answer is required.", call. = FALSE)
  }
  answer <- trimws(answer)
  if (name == "randomised_treatment") {
    if (!toupper(answer) %in% c("TRUE", "FALSE")) {
      stop("Enter TRUE or FALSE.", call. = FALSE)
    }
    return(toupper(answer) == "TRUE")
  }
  if (name %in% c("direct_select", "total_select")) {
    if (tolower(answer) == "auto") {
      return("auto")
    }
    answer <- sub("^c\\((.*)\\)$", "\\1", answer)
    answer <- sub("^\\[(.*)\\]$", "\\1", answer)
    values <- suppressWarnings(as.numeric(trimws(strsplit(answer, ",", fixed = TRUE)[[1L]])))
    if (length(values) != 2L || anyNA(values) || any(!is.finite(values)) ||
      any(values < 1 | values != floor(values))) {
      stop("Enter group, option (for example 1, 2), or auto.", call. = FALSE)
    }
    return(as.integer(values))
  }
  values <- if (name == "crop") {
    answer
  } else {
    trimws(strsplit(answer, ",", fixed = TRUE)[[1L]])
  }
  return(.gd_query_value(name, values))
}


.gd_prompt_value <- function(name, responses, validate = identity) { # Collect and record responses until one passes validation
  prompt <- switch(name,
    crop = "Crop: ", treatments = "Treatment nodes (comma-separated): ",
    outcomes = "Outcome nodes (comma-separated): ",
    randomised_treatment = "Was the treatment randomised? Enter TRUE or FALSE: ",
    paste0(if (name == "direct_select") {
      "Direct"
    } else {
      "Total"
    },
    " set: enter group, option, or auto: "))
  return(repeat {
    answer <- .gd_readline(prompt)
    if (!is.character(answer) || length(answer) != 1L || is.na(answer)) {
      stop("Input collection ended without an answer for ", name, ".", call. = FALSE)
    }
    responses$rows <- rbind(responses$rows,
      data.frame(field = name, prompt = prompt, response = answer))
    value <- tryCatch(validate(.gd_answer(name, answer)), error = identity)
    if (!inherits(value, "error")) {
      return(value)
    }
    message(conditionMessage(value))
  })
}


.gd_flag <- function(x, name, nullable = FALSE) { # Validate a scalar logical flag with optional NULL support
  if (nullable && is.null(x)) {
    return(invisible(x))
  }
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE", if (nullable) {
      " or NULL"
    } else {
      ""
    }, ".", call. = FALSE)
  }
  return(invisible(x))
}


.gd_ui_options <- function(options, selectors, responses) { # Resolve one deferred selector and return to construction
  for (family in c("direct", "total")) {
    key <- paste0(family, "_select")
    if (!is.null(selectors[[key]])) {
      next
    }
    table <- options[[family]]
    if (is.null(table) || !nrow(table)) {
      selectors[[key]] <- "auto"
    } else if (nrow(table) == 1L) {
      selectors[[key]] <- as.integer(c(table$group[1L], table$option[1L]))
    } else {
      cat(if (family == "direct") {
        "Direct"
      } else {
        "Total"
      }, " adjustment choices:\n", sep = "")
      .gd_print_options(table)
      valid <- function(value) { # Check that the selected group and option exist
        if (identical(value, "auto")) {
          return(value)
        }
        if (!any(table$group == value[1L] & table$option == value[2L])) {
          stop("That group and option are absent from the displayed catalogue.", call. = FALSE)
        }
        return(value)
      }
      selectors[[key]] <- .gd_prompt_value(key, responses, valid)
    }
    # Return to construction after choosing direct: total compatibility depends
    # on that actual choice, not on an undisclosed provisional selection.
    return(selectors)
  }
  return(selectors)
}


# Constructor helpers deliberately operate on the supplied emitted graph.
.cd_record <- function(graph) { # Read the saved query receipt without partial attribute matching
  return(attr(graph, "query_record", exact = TRUE))
}


.cd_set_latents <- function(graph, nodes) { # Set latent nodes while preserving identity metadata
  labels <- attr(graph, "input_labels", exact = TRUE)
  map <- attr(graph, "cf_node_map", exact = TRUE)
  source_order <- attr(graph, "cf_source_order", exact = TRUE)
  if (is.null(source_order)) {
    source_order <- .cd_record(graph)$source_order
  }
  dagitty::latents(graph) <- nodes
  attr(graph, "input_labels") <- labels
  attr(graph, "cf_node_map") <- map
  attr(graph, "cf_source_order") <- source_order
  return(graph)
}


# The same graph with every latent status cleared. Attributes are kept.
.cd_unmarked <- function(graph) { # Clear latent status while retaining the graph and its attributes
  if (!length(dagitty::latents(graph))) {
    return(graph)
  }
  return(.cd_set_latents(graph, character()))
}


# Declared unmeasured variables are marked on a finished construction that ran
# as for a fully observed graph. A display node is unmeasured when any occasion
# it stands for is declared: adjusting for it would need that occasion
# observed. A display node that also stands for observed occasions is recorded
# in `partial` (its unmeasured and observed occasions) and warned about.
.cd_mark_unmeasured <- function(summary, latent_nodes) { # Mark declared unmeasured occasions on the finished temporal and display graphs
  map <- summary$node_map
  nodes <- if (is.data.frame(map)) {
    as.character(map$node)
  } else {
    names(summary$graph)
  }
  display <- if (is.data.frame(map) && "display" %in% names(map)) {
    as.character(map$display)
  } else {
    nodes
  }
  display_latent <- unique(display[nodes %in% latent_nodes])
  temporal_latent <- intersect(names(summary$temporal), latent_nodes)
  partial <- list()
  for (node in display_latent) {
    occasions <- unique(nodes[display == node])
    declared <- intersect(occasions, latent_nodes)
    if (length(declared) < length(occasions)) {
      partial[[node]] <- list(unmeasured = declared, observed = setdiff(occasions, declared))
    }
  }
  if (length(partial)) {
    warning(paste(vapply(names(partial), function(node) { # Explain a display node that represents mixed observation status
      return(paste0(node,
        " stands for both observed and unmeasured occasions (unmeasured: ",
        paste(partial[[node]]$unmeasured, collapse = ", "), ", observed: ",
        paste(partial[[node]]$observed, collapse = ", "), ")"))
    }, character(1)), collapse = ". "),
      if (length(partial) == 1L) {
        ". It is shown unmeasured."
      } else {
        ". They are shown unmeasured."
      }, call. = FALSE)
  }
  if (length(display_latent)) {
    summary$graph <- .cd_set_latents(summary$graph, display_latent)
  }
  if (length(temporal_latent)) {
    summary$temporal <- .cd_set_latents(summary$temporal, temporal_latent)
  }
  return(list(summary = summary, display = display_latent, temporal = temporal_latent, partial = partial))
}


# Structural checks keep their status vocabulary. Each record and family also
# says which of its nodes are marked unmeasured on the checked graph. A family
# with no chosen set represented on that graph (no set, a budget stop, or an
# unrepresented choice) has no selected set to make available.
.cd_annotate_unmeasured <- function(checks, latent) { # Annotate selected memberships with the unmeasured nodes they require
  for (effect in names(checks)) {
    check <- checks[[effect]]
    check$records <- lapply(check$records, function(record) { # Record the unmeasured nodes used by one checked membership
      record$unmeasured <- if (is.null(record$mapped)) {
        character()
      } else {
        intersect(record$mapped, latent)
      }
      return(record)
    })
    chosen <- Filter(function(record) { # Select records belonging to the chosen option
      return(isTRUE(record$chosen))
    }, check$records)
    represented <- Filter(function(record) { # Keep chosen records represented on the checked graph
      return(!is.null(record$mapped))
    }, chosen)
    needs <- sort(unique(as.character(unlist(lapply(chosen, `[[`, "unmeasured"), use.names = FALSE))), method = "radix")
    check$unmeasured <- needs
    check$availability <- if (!length(represented)) {
      "no selected set"
    } else if (length(needs)) {
      paste0("needs unmeasured: ", paste(needs, collapse = ", "))
    } else {
      "available now"
    }
    checks[[effect]] <- check
  }
  return(checks)
}


# Exact treatment parents retained for joint total identification but declared
# unmeasured leave that safeguard uncertified. The occasions stay retained.
.cd_unmeasured_parents <- function(identification, latent_nodes) { # Record unmeasured treatment parents without removing their occasions
  unmeasured <- intersect(identification$parents, latent_nodes)
  if (!isTRUE(identification$applied) || !length(unmeasured)) {
    return(identification)
  }
  identification$unmeasured_parents <- unmeasured
  identification$status <- paste0("not certified: treatment parent(s) ", paste(unmeasured, collapse = ", "),
    " are declared unmeasured. Exact treatment parents and temporal occurrences are retained to show what must be observed")
  return(identification)
}


.cd_node_map <- function(graph, record = .cd_record(graph)) { # Read or reconstruct the supplied graph's canonical node map
  map <- record$node_map
  if (is.null(map)) {
    map <- attr(graph, "cf_node_map", exact = TRUE)
  }
  if (is.data.frame(map) && all(c("node", "quantity", "member_id", "occurrence") %in% names(map))) {
    map <- map[map$node %in% names(graph), , drop = FALSE]
    if (!setequal(unique(map$node), names(graph))) {
      stop("The query receipt's node map does not match the supplied graph.", call. = FALSE)
    }
    return(map)
  }
  ids <- .gd_node_ids(graph)
  return(data.frame(node = names(graph), quantity = unname(ids), member_id = unname(ids),
    occurrence = names(graph), stringsAsFactors = FALSE))
}


.cd_working_graph <- function(graph, map, source_order = NULL) { # Attach the validated canonical identities to a working graph
  result <- .gd_canonical(graph)
  if (!is.null(source_order)) {
    attr(result, "cf_source_order") <- source_order
  }
  # Complete reverse labels, including human labels, are analysis metadata.
  pairs <- unique(map[c("node", "quantity")])
  if (anyDuplicated(pairs$node)) {
    stop("A displayed node maps to multiple canonical quantity IDs.", call. = FALSE)
  }
  attr(result, "input_labels") <- stats::setNames(pairs$quantity, .strip_slice(pairs$node))
  attr(result, "cf_node_map") <- map
  return(result)
}


.cd_match_nodes <- function(graph, requested, record = .cd_record(graph)) { # Resolve requested quantities and exact occasions against the saved receipt
  requested <- .gd_names(requested, "node names")
  map <- .cd_node_map(graph, record)
  result <- character()
  for (value in requested) {
    exact_time <- value != .strip_slice(value)
    hit <- if (exact_time) {
      unique(map$node[map$node == value | map$occurrence == value])
    } else {
      unique(map$node[map$quantity == value | map$member_id == value |
        .strip_slice(map$node) == value])
    }
    if (!length(hit) && value %in% names(graph)) {
      hit <- value
    }
    # Recorded input aliases are safe to replay without consulting changed tables.
    if (!length(hit)) {
      requests <- record$requests
      resolved <- unlist(lapply(requests, function(x) { # Read the canonical identities resolved for one input family
        return(x$resolved)
      }), use.names = FALSE)
      aliases <- unlist(lapply(requests, function(x) { # Read the recorded names from one input family
        return(names(x$resolved))
      }), use.names = FALSE)
      canonical <- unique(resolved[aliases == value])
      hit <- unique(map$node[map$member_id %in% canonical | map$quantity %in% canonical])
    }
    if (!length(hit) && length(record$aliases)) {
      resolved <- .resolve_columns(value, .build_alias_index(record$aliases),
        record$arguments$params$disambiguate, character())
      if (length(resolved$ambiguous)) {
        stop("Ambiguous node name: ", value, call. = FALSE)
      }
      canonical <- unname(resolved$resolved)
      hit <- unique(map$node[map$member_id %in% canonical | map$quantity %in% canonical])
    }
    if (exact_time && length(hit) && any(vapply(hit, function(node) { # Check whether a displayed node merges multiple exact occasions
      return(length(unique(map$occurrence[map$node == node])) > 1L)
    }, logical(1)))) {
      stop("The supplied graph merges multiple occasions for ", value, ". A request for an exact occasion cannot distinguish them.", call. = FALSE)
    }
    if (!length(hit)) {
      stop("Quantity absent from the supplied graph: ", value, call. = FALSE)
    }
    result <- union(result, hit)
  }
  return(result)
}


.cd_prepare <- function(graph, latents, keep_nodes) { # Reuse exact preparation when acceleration and the graph contract allow it
  if (!.cf_acceleration_enabled() || !.cf_standard_graph(graph)) {
    return(.cd_prepare_uncached(graph, latents, keep_nodes))
  }
  return(.cf_cache$call("preparation", list(graph, latents, keep_nodes),
    function() { # Prepare the supplied graph after a cache miss
      return(.cd_prepare_uncached(graph, latents, keep_nodes))
    }))
}


.cd_prepare_uncached <- function(graph, latents, keep_nodes) { # Validate the supplied graph and prepare its policy and selection references
  if (!inherits(graph, "dagitty")) {
    stop("graph must be a dagitty query_graph result.", call. = FALSE)
  }
  if (!length(dagitty::exposures(graph)) || !length(dagitty::outcomes(graph))) {
    stop("The supplied graph must tag treatments and outcomes.", call. = FALSE)
  }
  record <- .cd_record(graph)
  map <- .cd_node_map(graph, record)
  reference <- .cd_working_graph(graph, map)
  if (!dag_is_acyclic(reference)) {
    stop("The supplied graph has directed cycles. Provide an acyclic query configuration.", call. = FALSE)
  }
  if (length(intersect(dagitty::exposures(reference), dagitty::outcomes(reference)))) {
    stop("Treatment and outcome declarations overlap. Construction requires disjoint endpoint sets.", call. = FALSE)
  }
  latent_nodes <- .cd_match_nodes(graph, latents, record)
  if (length(intersect(latent_nodes, c(dagitty::exposures(graph), dagitty::outcomes(graph))))) {
    stop("Treatment/outcome nodes cannot be declared latent for these adjustment searches.", call. = FALSE)
  }
  reference <- .cd_set_latents(reference, latent_nodes)
  endpoint <- unique(c(dagitty::exposures(reference), dagitty::outcomes(reference)))
  columns <- unique(unname(record$requests$columns$resolved))
  column_nodes <- unique(map$node[map$member_id %in% columns | map$quantity %in% columns])
  missing_columns <- setdiff(columns, unique(c(map$member_id, map$quantity)))
  extra_nodes <- .cd_match_nodes(graph, keep_nodes, record)
  requested <- union(column_nodes, extra_nodes)
  required <- union(endpoint, requested)
  policy <- record$policy
  crop <- record$arguments$crop
  if (is.null(crop)) {
    crop <- attr(graph, "crop", exact = TRUE)
  }
  applied <- (isTRUE(record$policy_applied) || identical(crop, "barley")) && is.data.frame(policy) &&
  all(c("node_id", "recovered_node") %in% names(policy))
  member_rule <- rep("yes", nrow(map))
  if (applied) {
    policy_rows <- match(map$member_id, policy$node_id)
    member_rule[!is.na(policy_rows)] <- as.character(policy$recovered_node[policy_rows[!is.na(policy_rows)]])
  }
  if (anyNA(member_rule) || any(!member_rule %in% c("yes", "yes:user") & !startsWith(member_rule, "no:"))) {
    stop("The query receipt contains invalid recovery policy values.", call. = FALSE)
  }
  explicit_members <- unique(map$member_id[map$node %in% required])
  allowed_member <- member_rule == "yes" | (member_rule == "yes:user" & map$member_id %in% explicit_members)
  eligible <- vapply(names(reference), function(node) { # Check whether every source member of this node is eligible
    return(all(allowed_member[map$node == node]))
  }, logical(1))
  prohibited <- intersect(required, names(reference)[!eligible])
  if (length(prohibited)) {
    stop("Requested nodes conflict with the saved recovery policy: ", paste(prohibited, collapse = ", "), call. = FALSE)
  }
  rules <- vapply(names(reference), function(node) { # Collect the recovery rules represented by this node
    return(paste(unique(member_rule[map$node == node]), collapse = "; "))
  }, character(1))
  ledger <- data.frame(node = names(reference), node_id = unname(.gd_node_ids(reference)),
    metadata_row = names(reference) %in% map$node, recovered_node = rules,
    explicitly_requested = names(reference) %in% required, policy_applied = applied,
    eligible = eligible, action = ifelse(eligible, "retain", "fold"), stringsAsFactors = FALSE)
  excluded <- ledger$node[!eligible]
  analysis <- .cd_fold_policy(reference, excluded)
  analysis <- .cd_working_graph(analysis, map[map$node %in% names(analysis), , drop = FALSE])
  # Recovery exclusions are metadata simplifications, not latent variables.
  # Search and certificates use that graph folded through directed paths. The supplied graph
  # and reference before folding remain available for provenance and edit history.
  selection <- .gd_canonical(analysis)
  source_order <- attr(selection, "cf_source_order", exact = TRUE)
  if (!is.null(source_order) && !.query_order_applies(c(unname(.gd_node_ids(selection)),
    map$quantity, map$member_id), source_order)) {
    source_order <- NULL
  }
  registry <- list()
  for (id in unique(map$quantity)) {
    members <- .query_sort(unique(map$member_id[map$quantity == id]), source_order)
    if (!identical(members, id)) {
      registry[[id]] <- list(members = members,
        reason = "Exact source members recorded by query_graph.", source = "query_record$node_map")
    }
  }
  depth <- record$arguments$params$query_depth
  if (is.null(depth)) {
    depth <- attr(graph, "query_depth", exact = TRUE)
  }
  limitations <- character()
  if (length(excluded)) {
    limitations <- c(limitations,
      "Recovery policy exclusions were folded through directed paths, not treated as latent variables. Adjustment checks use the folded graph. The supplied query is preserved separately.")
  }
  if (is.null(record)) {
    limitations <- c(limitations, "No query receipt: original requests and policy provenance are unavailable. Quantity counts use stored node identities.")
  }
  if (length(depth) && is.finite(depth)) {
    limitations <- c(limitations,
      paste0("Checks use the supplied query_depth=", depth, " graph only. Excluded context was not substituted from a full snapshot."))
  }
  if (length(missing_columns)) {
    limitations <- c(limitations,
      paste0("Requested columns absent from the supplied graph: ", paste(missing_columns, collapse = ", "), "."))
  }
  unresolved <- record$requests$columns$unresolved
  if (is.list(unresolved)) {
    unresolved_names <- unique(c(unresolved$unmatched, names(unresolved$ambiguous)))
    if (length(unresolved_names)) {
      limitations <- c(limitations,
        paste0("Column requests unresolved by query_graph: ", paste(unresolved_names, collapse = ", "), "."))
    }
  }
  return(list(graph = graph, record = record, map = map, reference = reference,
    selection = selection, analysis = analysis, ledger = ledger, registry = registry,
    requested = requested, required = required, columns = columns, column_nodes = column_nodes,
    extra_nodes = extra_nodes, missing_columns = missing_columns, latent_nodes = latent_nodes,
    limitations = limitations, policy_applied = applied,
    controls = list(crop = crop, water_model = if (!is.null(record$arguments$params$water_model)) {
      record$arguments$params$water_model
    } else {
      attr(graph, "water_model")
    },
    disambiguate = record$arguments$params$disambiguate, simplify = TRUE, keep_common_causes = FALSE)))
}


.cd_fold_policy <- function(graph, excluded) { # Fold nodes excluded by the policy through directed paths only
  if (!length(excluded)) {
    return(graph)
  }
  edges <- dag_edges(graph)
  incident <- edges$v %in% excluded | edges$w %in% excluded
  if (any(incident & edges$e != "->")) {
    stop("Recovery policy directed folding does not support edges that are not directed and are incident to excluded nodes.",
      call. = FALSE)
  }
  arrows <- edges[edges$e == "->", , drop = FALSE]
  folded <- .collapse_hidden(names(graph),
    data.table::data.table(from = arrows$v, to = arrows$w),
    stats::setNames(names(graph) %in% excluded, names(graph)))
  replacements <- data.frame(v = folded$edges$from, e = rep("->", nrow(folded$edges)),
    w = folded$edges$to, stringsAsFactors = FALSE)
  replacements <- replacements[!edge_key(replacements) %in% edge_key(edges), , drop = FALSE]
  result <- native_remove_nodes(graph, excluded, ".cd_fold_policy")
  result <- native_append_graph(result, character(), replacements, caller = ".cd_fold_policy")
  if (is.null(result)) {
    stop("Recovery policy directed folding could not preserve the remaining graph.", call. = FALSE)
  }
  source_order <- attr(graph, "cf_source_order", exact = TRUE)
  if (is.null(source_order)) {
    source_order <- .cd_record(graph)$source_order
  }
  attr(result, "cf_source_order") <- source_order
  return(result)
}


.cd_candidates <- function(raw, prepared, progress = NULL) { # Group adjustment sets by their exact source quantities
  source_order <- attr(prepared$selection, "cf_source_order", exact = TRUE)
  if (!is.null(source_order)) {
    family_ids <- .strip_slice(unlist(raw, use.names = FALSE))
    labels <- attr(prepared$selection, "input_labels")
    labelled <- family_ids %in% names(labels)
    family_ids[labelled] <- unname(labels[family_ids[labelled]])
    if (!.query_order_applies(c(unname(.gd_node_ids(prepared$selection)), family_ids,
      prepared$map$quantity, prepared$map$member_id), source_order)) {
      source_order <- NULL
    }
  }
  table <- .gd_candidates(raw, prepared$selection, list(), progress = progress, source_order = source_order)
  if (!nrow(table)) {
    return(table)
  }
  # Saved groups expand once to source quantities, rather than recursively
  # interpreting a group as its own member when a real row shares its name.
  table$resolved <- lapply(table$concepts, function(ids) { # Expand these recorded quantities to their source members
    return(.query_sort(unique(prepared$map$member_id[prepared$map$quantity %in% ids]), source_order))
  })
  table$n_resolved <- lengths(table$resolved)
  catalogue_keys <- table$label
  if (!is.null(source_order)) {
    catalogue_keys <- vapply(table$concepts, .query_set_key, character(1), source_order = source_order)
  }
  table <- table[order(table$n_resolved, catalogue_keys, method = "radix"), , drop = FALSE]
  table$group <- match(table$n_resolved, unique(table$n_resolved))
  table$option <- as.integer(ave(seq_len(nrow(table)), table$group, FUN = seq_along))
  table$group_size <- table$n_resolved
  table$id <- paste0("set_", seq_len(nrow(table)))
  rownames(table) <- NULL
  return(table)
}


.cd_witness_search <- function(prepared, effect, progress = NULL) { # Find and validate one deterministic source witness
  .gd_notify(progress, paste0("Finding one ", effect, " source witness"))
  condition_free <- TRUE
  raw <- withCallingHandlers(.gd_adjustment_engine(prepared$selection, effect, 1L, progress),
    condition = function(c) { # Record that the witness search signalled a condition
      return(condition_free <<- FALSE)
    })
  sets <- unname(lapply(raw, function(x) { # Represent one returned witness as a character vector
    return(as.character(unlist(x, use.names = FALSE)))
  }))
  if (length(sets) > 1L || !all(.gd_validate_sets(prepared$selection, sets, effect))) {
    stop("The ", effect, " source witness failed validation. No complete search was attempted.", call. = FALSE)
  }
  complete <- !length(sets) || !length(sets[[1L]])
  return(list(sets = sets, source_sets = sets, complete = complete, exhaustive_requested = FALSE, probe = NULL, limit = 1L,
    condition_free = condition_free, table = .cd_candidates(sets, prepared, progress),
    probe_quantity_count = NA_integer_, omitted_equal_smallest = FALSE,
    scope = "One deterministic source witness minimal under set inclusion. Other source alternatives and minimum quantity cost were not searched.",
    optimality = if (!length(sets)) {
      "no eligible set"
    } else if (!length(sets[[1L]])) {
      "valid empty set is the unique set minimal under set inclusion"
    } else {
      "not assessed"
    }))
}


.cd_reuse_witness <- function(previous, prepared, identity, progress = NULL) { # Reuse only an unchanged and fully validated source witness receipt
  if (is.null(previous) || is.null(identity) || !isTRUE(getOption("causaliflower.accelerate", TRUE))) {
    return(NULL)
  }
  result <- tryCatch({
    state <- .gd_set_source(previous)
    if (!identical(state$construction$search_method, "source-witness-v1") ||
      !identical(state$construction$search_identity, identity)) {
      return(NULL)
    }
    arguments <- state$construction$arguments
    original <- .cd_prepare(state$original_query, arguments$latents, arguments$keep_nodes)
    original_identity <- .cd_search_identity(original)
    if (is.null(original_identity)) {
      return(NULL)
    }
    original_identity$method <- "source-witness-v1"
    original_identity$witness_adapter <- paste(deparse(body(.cd_witness_search)), collapse = "\n")
    if (!identical(original_identity, identity) ||
      !identical(.gd_signature(original$analysis), .gd_signature(state$analysis)) ||
      !identical(attr(original$analysis, "cf_node_map"), attr(state$analysis, "cf_node_map"))) {
      return(NULL)
    }
    search <- state$build$search
    for (effect in c("direct", "total")) {
      item <- search[[effect]]
      sets <- item$sets
      if (!is.list(sets) || length(sets) > 1L ||
        !all(vapply(sets, function(x) { # Check that a saved witness has valid unique node names
          return(is.character(x) && !anyNA(x) && !anyDuplicated(x))
        }, logical(1))) ||
        !identical(item$source_sets, sets) || !identical(item$limit, 1L) ||
        !is.null(item$probe) || !isTRUE(item$condition_free) ||
        !identical(item$exhaustive_requested, FALSE) ||
        !identical(item$complete, !length(sets) || !length(sets[[1L]])) ||
        !identical(item$table, .cd_candidates(sets, original)) ||
        !all(.gd_validate_sets(prepared$selection, sets, effect))) {
        return(NULL)
      }
      item$table <- .cd_candidates(sets, prepared)
      search[[effect]] <- item
    }
    search
  }, error = function(e) { # Decline reuse after a signalled condition
    return(NULL)
  }, warning = function(w) { # Decline reuse after a signalled condition
    return(NULL)
  }, message = function(m) { # Decline reuse after a signalled condition
    return(NULL)
  })
  if (!is.null(result)) {
    .gd_notify(progress, "Reusing exact source witnesses. Rebuilding construction and Output checks")
  }
  return(result)
}


.cd_search <- function(prepared, effect, max_sets, progress = NULL) { # Enumerate adjustment alternatives and record search conditions
  .gd_notify(progress, paste0("Searching ", effect, " adjustment alternatives"))
  condition_free <- TRUE
  result <- withCallingHandlers(
    if (is.null(progress)) {
      .gd_adjustment_sets(prepared$selection, effect, max_sets)
    } else {
      .gd_adjustment_sets(prepared$selection, effect, max_sets, progress = progress)
    },
    condition = function(c) { # Record that the adjustment search signalled a condition
      return(condition_free <<- FALSE)
    })
  result$condition_free <- condition_free
  .gd_notify(progress, paste0(effect, " search complete (", length(result$sets), " alternatives). Grouping source quantities"))
  return(.cd_group_search(result, prepared, effect, progress))
}


.cd_group_search <- function(result, prepared, effect, progress = NULL) { # Group search results and describe any unfinished boundary
  table <- .cd_candidates(result$sets, prepared, progress)
  .gd_notify(progress, paste0(effect, " catalogue ready (", nrow(table), " alternatives by source quantity)"))
  probe_count <- if (!isTRUE(result$complete)) {
    .cd_candidates(list(result$probe), prepared)$n_resolved
  } else {
    NA_integer_
  }
  result$probe_quantity_count <- probe_count
  result$omitted_equal_smallest <- !isTRUE(result$complete) && nrow(table) > 0L &&
  !is.na(probe_count) && !is.na(table$n_resolved[1L]) && probe_count == table$n_resolved[1L]
  result$table <- table
  return(result)
}


# Retention requests are deliberately absent. Everything that defines the
# search itself remains exact, including graph text/order, temporal identities,
# declarations, policy decisions, query scope and separator implementation.
# Fresh preparation precedes this comparison, so columns cannot silently admit
# a node folded by the policy or enlarge a restricted query while retaining old sets.
.cd_search_identity <- function(prepared) { # Read an exact search identity or decline unsafe reuse
  return(tryCatch(.cd_search_identity_value(prepared), error = function(e) { # Decline an identity after a signalled condition
    return(NULL)
  },
    warning = function(w) { # Decline an identity after a signalled condition
      return(NULL)
    }, message = function(m) { # Decline an identity after a signalled condition
      return(NULL)
    }))
}


.cd_search_identity_value <- function(prepared) { # Record every input and implementation that defines the search
  record <- prepared$record
  arguments <- record$arguments
  arguments[c("columns", "quiet")] <- NULL
  policy <- prepared$ledger[setdiff(names(prepared$ledger), "explicitly_requested")]
  script <- system.file("js", "dagitty-alg.js", package = "dagitty")
  if (!nzchar(script) || !file.exists(script)) {
    return(NULL)
  }
  identity <- list(version = 1L, method = "complete-source-sets-v2",
    engine = list(dagitty = as.character(utils::packageVersion("dagitty")),
    native = readBin(script, "raw", n = file.info(script)$size),
    adapter = deparse(body(.gd_adjustment_sets)),
    implementation = deparse(body(.gd_adjustment_engine)),
    separator = "Installed dagitty GraphAnalyzer.listMinimalSeparators", indexed = FALSE),
    # Alias/case handling has already run in fresh preparation. The separator
    # sees UTF-8 graph text and uses JavaScript ordering, independent of CTYPE.
    # Existing saved signatures/candidate receipts still use R collation.
    collation = Sys.getlocale("LC_COLLATE"),
    selection = charToRaw(enc2utf8(as.character(prepared$selection))),
    signature = .gd_signature(prepared$selection),
    reference = charToRaw(enc2utf8(as.character(prepared$reference))),
    map = prepared$map, analysis_map = attr(prepared$analysis, "cf_node_map", exact = TRUE),
    policy = policy, policy_applied = prepared$policy_applied,
    source_policy = record$source_policy, query_policy = record$policy,
    source_node_map = record$source_node_map, scope = record$scope,
    arguments = arguments, controls = prepared$controls)
  source_order <- attr(prepared$selection, "cf_source_order", exact = TRUE)
  if (!is.null(source_order)) {
    identity$source_order <- source_order
    identity$ordering_adapter <- list(
      order_key = deparse(body(.query_order_key)), sort = deparse(body(.query_sort)),
      applies = deparse(body(.query_order_applies)), set_key = deparse(body(.query_set_key)))
  }
  return(identity)
}


.cd_reuse_search <- function(previous, prepared, identity, progress = NULL) { # Reuse a certified complete search only for an identical prepared reference
  if (is.null(previous) || is.null(identity) ||
    !isTRUE(getOption("causaliflower.accelerate", TRUE))) {
    return(NULL)
  }
  # Malformed, old, manually changed or incomplete receipts take the ordinary
  # constructor path. They never become evidence for a complete new catalogue.
  state <- tryCatch(.gd_set_source(previous), error = function(e) { # Decline reuse when previous provenance cannot be read
    return(NULL)
  })
  if (is.null(state) || !identical(state$construction$search_identity, identity)) {
    return(NULL)
  }
  arguments <- state$construction$arguments
  original_preparation <- tryCatch(.cd_prepare(state$original_query, arguments$latents, arguments$keep_nodes),
    error = function(e) { # Decline reuse when the prior reference cannot be prepared cleanly
      return(NULL)
    }, warning = function(w) { # Decline reuse when the prior reference cannot be prepared cleanly
      return(NULL)
    }, message = function(m) { # Decline reuse when the prior reference cannot be prepared cleanly
      return(NULL)
    })
  if (is.null(original_preparation) || !identical(.cd_search_identity(original_preparation), identity) ||
    !identical(.gd_signature(original_preparation$analysis), .gd_signature(state$analysis)) ||
    !identical(attr(original_preparation$analysis, "cf_node_map", exact = TRUE),
    attr(state$analysis, "cf_node_map", exact = TRUE))) {
    return(NULL)
  }
  search <- state$build$search
  for (effect in c("direct", "total")) {
    item <- search[[effect]]
    if (!isTRUE(item$complete) || !identical(item$limit, Inf) ||
      !is.null(item$probe) || !isTRUE(item$condition_free) ||
      !is.list(item$sets) || !identical(item$source_sets, item$sets) ||
      length(item$refinement_sets) ||
      !all(vapply(item$sets, function(z) { # Check the saved set's names and membership in the selection reference
        return(is.character(z) && !anyNA(z) &&
          !anyDuplicated(z) && all(z %in% names(prepared$selection)))
      }, logical(1)))) {
      return(NULL)
    }
    if (!identical(item$table, .cd_candidates(item$sets, original_preparation))) {
      return(NULL)
    }
  }
  .gd_notify(progress, "Reusing complete adjustment searches for the identical prepared reference")
  return(stats::setNames(lapply(c("direct", "total"), function(effect) { # Rebuild the catalogue from this family's saved source sets
    # Rebuild grouping, automatic choices, projection and checks from new inputs.
    raw <- search[[effect]][c("sets", "complete", "probe", "limit", "condition_free")]
    return(.cd_group_search(raw, prepared, effect, progress))
  }), c("direct", "total")))
}


# Search refinement verified against the source ---------------------------------------------

# Propose smaller sets from temporal projections, but certify every admitted set
# against the unchanged selection reference. This is a bounded search heuristic.
.cd_refine_search <- function(prepared, search, max_sets) { # Certify smaller proposals from temporal projection on the unchanged source
  families <- c("direct", "total")
  source_signature <- .gd_signature(prepared$selection)
  source_sets <- lapply(search, function(x) { # Preserve this family's original source sets
    return(x$sets)
  })
  incomplete <- families[!vapply(search[families], function(x) { # Check whether this source family was fully enumerated
    return(isTRUE(x$complete))
  }, logical(1))]
  for (effect in families) {
    search[[effect]]$source_sets <- source_sets[[effect]]
  }
  receipt <- list(method = "source-verified temporal-projection refinement",
    source_signature = source_signature, source_sets = source_sets,
    proposal_limit = if (is.finite(max_sets)) {
      as.integer(max_sets)
    } else {
      Inf
    }, pass_limit = 3L,
    certification_limit = 1L, certification_probe = 1L,
    objective = "Sum of selected counts of distinct source nodes, counting the two families separately.",
    global_minimum_certified = FALSE, passes = list(),
    stop_reason = "source search complete")
  if (!length(incomplete)) {
    return(list(search = search, refinement = receipt))
  }
  choose <- function(search) { # Choose compatible direct and total options and count their source members
    direct <- .cd_select(search$direct$table, "auto", "direct")
    retained <- unique(prepared$map$member_id[prepared$map$node %in% prepared$required])
    if (length(direct$index)) {
      retained <- union(retained, search$direct$table$resolved[[direct$index]])
    }
    total <- .cd_select(search$total$table, "auto", "total", retained)
    choices <- list(direct = direct, total = total)
    counts <- vapply(families, function(effect) { # Read the selected option's count of distinct source members
      selected_index <- choices[[effect]]$index
      return(if (length(selected_index)) {
        search[[effect]]$table$n_resolved[selected_index]
      } else {
        Inf
      })
    }, numeric(1))
    if (length(total$index)) {
      retained <- union(retained, search$total$table$resolved[[total$index]])
    }
    return(list(choices = choices, counts = counts, objective = sum(counts), retained_union_count = length(retained)))
  }
  key <- function(membership) { # Build a stable identity for one adjustment set
    return(paste(sort(unique(membership), method = "radix"), collapse = "\r"))
  }
  selected <- choose(search)
  if (any(!is.finite(selected$counts))) {
    receipt$stop_reason <- "no automatic source choice available"
    return(list(search = search, refinement = receipt))
  }
  for (pass in seq_len(receipt$pass_limit)) {
    chosen <- lapply(families, function(effect) { # Read the chosen memberships for this effect family
      return(.cd_selected_sets(search[[effect]], selected$choices[[effect]]))
    })
    nodes <- unique(unlist(chosen, recursive = TRUE, use.names = FALSE))
    ids <- unique(unname(.gd_node_ids(prepared$reference)[nodes]))
    retained <- union(names(prepared$analysis)[.gd_node_ids(prepared$analysis) %in% ids], prepared$required)
    omit <- setdiff(names(prepared$analysis), retained)
    projection <- if (length(omit)) {
      collapse_nodes(prepared$analysis, omit)
    } else {
      prepared$analysis
    }
    projection <- .cd_working_graph(projection, prepared$map[prepared$map$node %in% names(projection), , drop = FALSE],
      source_order = attr(prepared$analysis, "cf_source_order", exact = TRUE))
    step <- list(pass = pass, source_signature = source_signature,
      projection_signature = .gd_signature(projection),
      projection_nodes = names(projection), selected_before = selected$counts,
      objective_before = selected$objective, union_count_before = selected$retained_union_count,
      families = list())
    updated <- search
    for (effect in incomplete) {
      proposals <- .gd_adjustment_sets(projection, effect, max_sets)
      table <- .cd_candidates(proposals$sets, prepared)
      ledger <- lapply(proposals$sets, function(membership) { # Start the assessment record for one proposed membership
        return(list(nodes = membership, status = "not assessed",
          reason = "Not in an assessed smaller group."))
      })
      record <- list(proposal_complete = proposals$complete,
        proposal_count = length(proposals$sets), proposal_probe = proposals$probe,
        proposal_table = table, candidates = ledger,
        accepted = list(), accepted_group = NULL)
      keys <- vapply(proposals$sets, key, character(1))
      groups <- unique(table$group[table$n_resolved < selected$counts[[effect]]])
      for (group in groups) {
        rows <- which(table$group == group)
        raw <- unlist(table$raw_sets[rows], recursive = FALSE)
        raw <- raw[!duplicated(vapply(raw, key, character(1)))]
        accepted <- list()
        for (membership in raw) {
          certificate <- .cd_certify_refinement(prepared$selection, membership, effect)
          index <- which(keys == key(membership))
          for (i in index) {
            record$candidates[[i]] <- c(list(nodes = membership), certificate)
          }
          if (identical(certificate$status, "accepted")) {
            accepted[[length(accepted) + 1L]] <- membership
          }
        }
        if (length(accepted)) {
          pooled <- c(updated[[effect]]$sets, accepted)
          pooled <- pooled[!duplicated(vapply(pooled, key, character(1)))]
          updated[[effect]]$sets <- pooled
          updated[[effect]]$table <- .cd_candidates(pooled, prepared)
          updated[[effect]]$refinement_sets <- c(updated[[effect]]$refinement_sets, accepted)
          record$accepted <- accepted
          record$accepted_group <- group
          break
        }
      }
      step$families[[effect]] <- record
    }
    after <- choose(updated)
    step$selected_after <- after$counts
    step$objective_after <- after$objective
    step$union_count_after <- after$retained_union_count
    step$strict_improvement <- after$objective < selected$objective
    receipt$passes[[pass]] <- step
    if (!step$strict_improvement) {
      receipt$stop_reason <- "no smaller group verified on the source within this pass"
      break
    }
    search <- updated
    selected <- after
    receipt$stop_reason <- "internal pass limit reached"
  }
  stopifnot(identical(source_signature, .gd_signature(prepared$selection)))
  for (effect in families) {
    # Completeness describes original enumeration, never the proposal catalogue.
    stopifnot(identical(search[[effect]]$source_sets, source_sets[[effect]]))
    probe_count <- search[[effect]]$probe_quantity_count
    table <- search[[effect]]$table
    search[[effect]]$omitted_equal_smallest <- !isTRUE(search[[effect]]$complete) &&
    nrow(table) > 0L && length(probe_count) == 1L && !is.na(probe_count) &&
    !is.na(table$n_resolved[1L]) && probe_count == table$n_resolved[1L]
  }
  receipt$selected_counts <- selected$counts
  return(list(search = search, refinement = receipt))
}


.cd_certify_refinement <- function(reference, nodes, effect) { # Certify eligibility, sufficiency and inclusion minimality on the source
  eligible <- all(nodes %in% names(reference)) &&
  !any(nodes %in% c(dagitty::exposures(reference), dagitty::outcomes(reference), dagitty::latents(reference)))
  if (!eligible) {
    return(list(status = "rejected", eligible = FALSE, sufficient = FALSE,
      minimal = FALSE, reason = "Nodes are missing, are endpoints, or are ineligible on the source."))
  }
  sufficient <- isTRUE(.gd_validate_sets(reference, list(nodes), effect)[[1L]])
  if (!sufficient) {
    return(list(status = "rejected", eligible = TRUE, sufficient = FALSE,
      minimal = FALSE, reason = "Does not satisfy the original direct/total criterion."))
  }
  # With only this candidate selectable, a returned proper subset disproves its
  # inclusion minimality. An identical minimal subset certifies it. No edges change.
  constrained <- .cd_set_latents(reference, setdiff(names(reference),
    c(nodes, dagitty::exposures(reference), dagitty::outcomes(reference))))
  certificate <- .gd_adjustment_sets(constrained, effect, 1L)
  witness <- if (length(certificate$sets)) {
    certificate$sets[[1L]]
  } else {
    NULL
  }
  minimal <- length(certificate$sets) > 0L && setequal(witness, nodes)
  return(list(status = if (minimal) {
    "accepted"
  } else {
    "rejected"
  }, eligible = TRUE,
    sufficient = TRUE, minimal = minimal, witness = witness,
    certificate_complete = certificate$complete,
    reason = if (minimal) {
      "Exact minimal subset returned on original edges with only candidate nodes eligible."
    } else {
      "Source search restricted to candidates returned a proper subset or no certificate."
    }))
}


.cd_select <- function(table, selector, family, retained = character(), pending = FALSE) { # Select a requested option or a stable compatible automatic choice
  if (!is.null(selector) && !identical(selector, "auto") &&
    (!is.numeric(selector) || length(selector) != 2L || anyNA(selector) ||
    any(!is.finite(selector)) || any(selector < 1 | selector != floor(selector)))) {
    stop(family, "_select must be 'auto', NULL, or c(group, option).", call. = FALSE)
  }
  if (!nrow(table)) {
    return(list(status = "no_sets", selection = NULL, index = integer(), ties = integer(), reason = "No eligible set found."))
  }
  if (is.null(selector) || (pending && identical(selector, "auto"))) {
    return(list(status = "deferred", selection = NULL, index = integer(), ties = integer(),
      reason = if (pending && identical(selector, "auto")) {
        "Pending direct selection."
      } else {
        "Selection deferred."
      }))
  }
  if (identical(selector, "auto")) {
    if (is.na(table$n_resolved[1L])) {
      stop("Automatic selection requires resolved quantity counts.", call. = FALSE)
    }
    candidates <- which(table$group == table$group[1L])
    if (family == "total") {
      overlap <- vapply(table$resolved[candidates], function(x) { # Count this option's overlap with already retained quantities
        return(length(intersect(x, retained)))
      }, integer(1))
      candidates <- candidates[overlap == max(overlap)]
    }
    index <- candidates[1L]
  } else {
    index <- which(table$group == selector[1L] & table$option == selector[2L])
    if (length(index) != 1L) {
      stop("That ", family, " [group, option] is absent from its catalogue.", call. = FALSE)
    }
    candidates <- index
  }
  return(list(status = "selected", selection = as.integer(c(table$group[index], table$option[index])),
    index = index, ties = candidates, reason = if (identical(selector, "auto")) {
      "Automatic stable choice."
    } else {
      "Explicit choice."
    }))
}


.cd_selected_sets <- function(search, choice) { # Read the exact temporal memberships of the selected option
  if (!length(choice$index)) {
    return(list())
  }
  return(search$table$raw_sets[[choice$index]])
}


.cd_identification_retention <- function(graph, total_choice) { # Retain exact treatment parents when ordinary total adjustment is unavailable
  result <- list(applied = FALSE, parents = character(), added_nodes = character(),
    status = "not required: a total adjustment option was selected",
    scope = "Joint total effect on the supplied working reference. This does not cover all adjustment alternatives or other effects.")
  if (!identical(total_choice$status, "no_sets")) {
    return(result)
  }
  result$status <- "not certified: the reference is not a fully observed directed DAG"
  if (dag_graph_type(graph) != "dag" ||
    length(dagitty::latents(graph)) || length(dagitty::adjustedNodes(graph)) ||
    length(selected_nodes(graph)) || any(dag_edges(graph)$e != "->")) {
    return(result)
  }
  # From the observed DAG g-formula, every intervention denominator depends
  # only on X and Pa(X). Keeping these exact vertices allows the other observed
  # variables to be marginalised without losing the joint total functional.
  # Pa(X) need not be an ordinary adjustment set (it may follow another X).
  treatment <- dagitty::exposures(graph)
  result$parents <- sort(unique(unlist(lapply(treatment, function(x) { # Read the exact parents of one treatment occurrence
    return(dagitty::parents(graph, x))
  }), use.names = FALSE)))
  result$applied <- TRUE
  result$status <- "pending structural verification"
  result$reason <- if (isTRUE(grepl("step budget", total_choice$reason, fixed = TRUE))) {
    "No total adjustment option was found within the ranked search's step budget (absence is not certified). Exact treatment parents and temporal occurrences are retained for joint total identification."
  } else {
    "No total adjustment option. Retain exact treatment parents and temporal occurrences for joint total identification."
  }
  return(result)
}


.cd_checks <- function(prepared, final, search, choices, max_sets, progress = NULL, verify_all = FALSE) { # Check selected memberships and final alternatives against the source
  family_checks <- lapply(c("direct", "total"), function(effect) { # Check one effect family on the constructed graph
    base <- search[[effect]]
    required <- .cd_selected_sets(base, choices[[effect]])
    .gd_notify(progress, paste0("Checking ", effect, " adjustment alternatives"))
    final_search <- if (is.null(progress)) {
      .gd_adjustment_sets(final, effect, max_sets)
    } else {
      .gd_adjustment_sets(final, effect, max_sets, progress = progress)
    }
    .gd_notify(progress, paste0("Verifying ", length(final_search$sets), " final ", effect, " alternatives against the reference"))
    check <- .gd_check_family(prepared$selection, final, base$sets, final_search$sets,
      effect, final_search$complete, required)
    check$catalogue_status <- check$status
    selected <- Filter(function(x) { # Select records belonging to the chosen membership
      return(isTRUE(x$chosen))
    }, check$records)
    introduced_bad <- any(vapply(check$introduced, function(x) { # Detect an introduced alternative that is invalid on the reference
      return(x$status != "valid on reference; minimality not certified")
    }, logical(1)))
    # A missing source family cannot excuse an invalid introduced set, and a
    # preserved selected set does not complete the search for new alternatives.
    check$status <- if (introduced_bad) {
      "changed"
    } else if (!length(required)) {
      "no eligible reference set"
    } else if (isTRUE(final_search$complete) && all(vapply(selected, function(x) { # Check whether the selected membership is preserved
      return(x$status == "preserved")
    }, logical(1)))) {
      "preserved"
    } else if (all(vapply(selected, function(x) { # Check whether the selected membership has a valid incomplete certificate
      return(x$status %in% c("preserved", "minimality not verified"))
    }, logical(1)))) {
      "not fully checked"
    } else {
      "changed"
    }
    check$reference_sets <- base$sets
    check$final_sets <- final_search$sets
    check$reference_complete <- base$complete
    check$final_complete <- final_search$complete
    check$required <- required
    if (verify_all) {
      check <- .cd_output_source_check(check, prepared$selection, final, effect)
    }
    .gd_notify(progress, paste0(effect, " check complete (", check$status, ")"))
    return(check)
  })
  return(stats::setNames(family_checks, c("direct", "total")))
}


# A repeated or grouped occurrence cannot be recovered by matching its label.
# Every Output alternative, including a mapped selected witness, is checked.
.cd_lift_exact <- function(set, from, to) { # Lift a set only when exact source identities and occasions agree
  mapped <- .gd_map_set(set, from, to)
  if (is.null(mapped) || length(mapped) != length(unique(set))) {
    return(NULL)
  }
  identity <- function(graph, node) { # Read the canonical source rows represented by one node
    map <- attr(graph, "cf_node_map", exact = TRUE)
    columns <- c("quantity", "member_id", "occurrence")
    if (is.data.frame(map) && all(c("node", columns) %in% names(map))) {
      value <- unique(map[map$node == node, columns, drop = FALSE])
      if (nrow(value)) {
        value[] <- lapply(value, as.character)
        value <- value[do.call(order, c(value, list(method = "radix"))), , drop = FALSE]
        rownames(value) <- NULL
        return(value)
      }
    }
    id <- unname(.gd_node_ids(graph)[node])
    return(data.frame(quantity = id, member_id = id, occurrence = node, stringsAsFactors = FALSE))
  }
  for (node in unique(set)) {
    target <- .gd_map_set(node, from, to)
    if (length(target) != 1L ||
      !identical(identity(from, node), identity(to, target))) {
      return(NULL)
    }
  }
  return(mapped)
}


# Reads every adjustment set of a display merged across time as all the source
# occurrences its nodes stand for (the summary node map's display column, so a
# time segment reads as its own occurrences only) and validates it on the
# unchanged selection reference. Returns NULL when all pass, else a reason.
.cd_summary_source_check <- function(reference, effects = c("direct", "total")) { # Validate adjustment alternatives from a display merged across time on the source
  return(function(display, mapping) { # Validate every displayed alternative after expanding its source occurrences
    for (effect in effects) {
      sets <- .gd_adjustment_sets(display, effect, Inf)$sets
      if (!length(sets)) {
        next
      }
      expanded <- lapply(sets, function(set) { # Expand merged nodes to their exact source occurrences
        return(as.character(unique(unlist(lapply(set, function(node) { # Expand merged nodes to their exact source occurrences
          occurrences <- mapping$node[mapping$display == node]
          return(if (length(occurrences)) {
            occurrences
          } else {
            node
          })
        }), use.names = FALSE))))
      })
      if (!all(.gd_validate_sets(reference, expanded, effect))) {
        return(paste0("A temporally merged ", effect, " adjustment set is not valid on the source when its merged nodes are read as all represented occurrences. The temporal display was retained."))
      }
    }
    return(NULL)
  })
}


.cd_output_source_check <- function(check, reference, final, effect) { # Certify final alternatives using exact temporal interpretations
  lifted <- lapply(check$final_sets, .cd_lift_exact, from = final, to = reference)
  liftable <- !vapply(lifted, is.null, logical(1))
  valid <- rep(FALSE, length(lifted))
  if (any(liftable)) {
    valid[liftable] <- .gd_validate_sets(reference, lifted[liftable], effect)
  }
  final_valid <- .gd_validate_sets(final, check$final_sets, effect)
  check$output_reference_sets <- lifted
  check$output_reference_valid <- valid
  check$output_reference_status <- ifelse(!liftable, "exact temporal interpretation unavailable",
    ifelse(valid, "valid on reference; source minimality not assessed", "not valid on reference"))
  check$output_valid <- final_valid
  check$output_reference_complete <- isTRUE(check$final_complete) && all(valid) && all(final_valid)
  check$scope <- "Selected exact source witnesses and every alternative on this constructed graph. Other source alternatives and global minimum quantity cost are not certified."
  if (!isTRUE(check$output_reference_complete)) {
    check$status <- "changed"
  }
  return(check)
}



#####
# get_dag_transforms.R
#####
# Internal graph transformations used by get_dag(). Analysis remains temporal.
# Summaries of repeated processes are display abstractions with separate ledgers.

.gd_names <- function(x, argument) { # Validate and deduplicate supplied node names
  if (is.null(x)) {
    return(character())
  }
  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    stop(argument, " must contain non-missing, non-empty names.", call. = FALSE)
  }
  return(unique(unname(x)))
}


.gd_canonical <- function(dag) { # Rebuild canonical graph syntax and retain applicable identity metadata
  dag_graph_type(dag)
  result <- dagitty::dagitty(as.character(dag))
  labels <- attr(dag, "input_labels")
  if (length(labels)) {
    attr(result, "input_labels") <- labels[names(labels) %in% .strip_slice(names(result))]
  }
  map <- attr(dag, "cf_node_map", exact = TRUE)
  if (is.data.frame(map)) {
    attr(result, "cf_node_map") <- map[map$node %in% names(result), , drop = FALSE]
  }
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  recorded_order <- .cd_record(dag)$source_order
  if (!is.null(source_order) && !is.null(recorded_order) && !identical(source_order, recorded_order)) {
    stop("The graph and its query record have different source ordering metadata.", call. = FALSE)
  }
  if (is.null(source_order)) {
    source_order <- recorded_order
  }
  if (!is.null(source_order)) {
    attr(result, "cf_source_order") <- source_order
  }
  return(result)
}


.gd_node_ids <- function(dag) { # Read canonical node identities with optional exact metadata reuse
  if (!.cf_acceleration_enabled() || !.cf_standard_graph(dag)) {
    return(.gd_node_ids_uncached(dag))
  }
  return(.gd_metadata(dag)$ids)
}


.gd_metadata <- function(dag) { # Reuse deterministic node metadata for the supplied graph
  return(.cf_cache$call("metadata", list(dag), function() { # Read graph metadata after a cache miss
    metadata <- list(ids = .gd_node_ids_uncached(dag), nodes = names(dag))
    source_order <- attr(dag, "cf_source_order", exact = TRUE)
    if (is.null(source_order)) {
      source_order <- .cd_record(dag)$source_order
    }
    if (!is.null(source_order)) {
      metadata$source_order <- source_order
    }
    return(metadata)
  }))
}


.gd_node_ids_uncached <- function(dag) { # Resolve graph vertices to canonical node identities
  ids <- .strip_slice(names(dag))
  labels <- attr(dag, "input_labels")
  matched <- ids %in% names(labels)
  ids[matched] <- unname(labels[ids[matched]])
  map <- attr(dag, "cf_node_map", exact = TRUE)
  if (is.data.frame(map)) {
    rows <- match(names(dag), map$node)
    ids[!is.na(rows)] <- map$quantity[rows[!is.na(rows)]]
  }
  return(stats::setNames(ids, names(dag)))
}


.gd_new_edges <- function(before, after) { # Return edges introduced by a graph transformation
  previous_edges <- dag_edges(before)
  current_edges <- dag_edges(after)
  return(current_edges[!edge_key(current_edges) %in% edge_key(previous_edges), , drop = FALSE])
}


.gd_recover <- function(dag, requested, water_model = "SoilWater",
  disambiguate = NULL, apply_policy = TRUE) { # Apply recovery eligibility and record every retained or folded quantity
  requested <- if (is.null(requested)) {
    character()
  } else {
    requested
  }
  .gd_names(requested, "requested")
  if (!is.logical(apply_policy) || length(apply_policy) != 1L || is.na(apply_policy)) {
    stop("apply_policy must be TRUE or FALSE.", call. = FALSE)
  }
  water_model <- match.arg(water_model, c("SoilWater", "SWIM"))
  nodes <- getExportedValue("causaliflower", "apsim_dag_nodes")
  aliases <- getExportedValue("causaliflower", "node_aliases")
  if (!is.null(.query_source_order(nodes)) && identical(water_model, "SWIM")) {
    stop("water_model = \"SWIM\" was requested, but the bundled slice has no SWIM family rows. Use water_model = \"SoilWater\".", call. = FALSE)
  }
  if (!"recovered_node" %in% names(nodes)) {
    stop("The package node table lacks recovered_node metadata.", call. = FALSE)
  }
  rule <- as.character(nodes$recovered_node)
  if (anyNA(rule) || any(!rule %in% c("yes", "yes:user") &
    !grepl("^no:[[:space:]]*[^[:space:]].*$", rule))) {
    stop("recovered_node must be yes, yes:user, or no: followed by a reason.",
      call. = FALSE)
  }
  inactive <- if (water_model == "SoilWater") {
    "SWIM."
  } else {
    "SoilWater."
  }
  aliases <- aliases[!startsWith(names(aliases), inactive)]
  resolved <- .resolve_columns(requested, .build_alias_index(aliases),
    disambiguate, nodes$node_id[nodes$hidden_node == "TRUE"])
  if (length(resolved$ambiguous)) {
    stop("Ambiguous explicit input: ",
      paste(names(resolved$ambiguous), collapse = ", "),
      ". Use canonical node IDs or disambiguate.", call. = FALSE)
  }
  requested_ids <- unique(c(unname(resolved$resolved),
    unname(requested)[unname(requested) %in% nodes$node_id]))
  ids <- .gd_node_ids(dag)
  endpoint <- names(dag) %in% c(dagitty::exposures(dag), dagitty::outcomes(dag))
  requested_ids <- union(requested_ids, unname(ids[endpoint]))
  prohibited <- intersect(requested_ids, nodes$node_id[startsWith(rule, "no:")])
  if (apply_policy && length(prohibited)) {
    stop("Explicit request conflicts with recovery policy: ",
      paste(paste0(prohibited, " (", rule[match(prohibited, nodes$node_id)], ")"),
      collapse = ", "), call. = FALSE)
  }

  row <- match(unname(ids), nodes$node_id)
  policy <- rep("yes", length(row))
  policy[!is.na(row)] <- rule[row[!is.na(row)]]
  explicit <- unname(ids) %in% requested_ids
  retain <- !apply_policy | policy == "yes" | (policy == "yes:user" & explicit)
  ledger <- data.frame(node = names(dag), node_id = unname(ids),
    metadata_row = !is.na(row), recovered_node = policy,
    explicitly_requested = explicit, policy_applied = rep(apply_policy, length(row)),
    action = ifelse(retain, "retain", "fold"),
    stringsAsFactors = FALSE)
  reference <- .gd_canonical(dag)
  omit <- ledger$node[!retain]
  graph <- if (length(omit)) {
    collapse_nodes(reference, omit)
  } else {
    reference
  }
  attr(graph, "cf_source_order") <- attr(reference, "cf_source_order", exact = TRUE)
  return(list(graph = .gd_canonical(graph), ledger = ledger,
    introduced_edges = .gd_new_edges(reference, graph),
    requested_ids = requested_ids,
    policy_applied = apply_policy,
    coverage = if (apply_policy) {
      "Manual recovered_node policy applied."
    } else {
      "Manual recovered_node policy not applied."
    },
    unresolved = list(ambiguous = resolved$ambiguous,
    unmatched = resolved$unmatched)))
}


.gd_common_causes <- function(dag, groups) { # Find nodes ancestral to at least two selected groups
  ancestors <- lapply(groups, function(group) { # Collect the proper ancestors of this selected group
    return(unique(unlist(lapply(group, function(node) { # Read the proper ancestors of one group member
      return(setdiff(dagitty::ancestors(dag, node), node))
    }), use.names = FALSE)))
  })
  counts <- vapply(names(dag), function(node) { # Count the selected groups for which this node is ancestral
    return(sum(vapply(ancestors, function(group) { # Check whether this ancestor set contains the node
      return(node %in% group)
    }, logical(1))))
  }, integer(1))
  return(names(counts)[counts >= 2L])
}


.gd_focus <- function(dag, ids, keep_common_causes = FALSE) { # Focus on requested quantities while recording folded nodes and new paths
  ids <- .gd_names(ids, "ids")
  if (!is.logical(keep_common_causes) || length(keep_common_causes) != 1L ||
    is.na(keep_common_causes)) {
    stop("keep_common_causes must be TRUE or FALSE.", call. = FALSE)
  }
  node_ids <- .gd_node_ids(dag)
  missing <- setdiff(ids, unname(node_ids))
  if (length(missing)) {
    stop("Focus quantities absent from this graph: ", paste(missing, collapse = ", "),
      call. = FALSE)
  }
  focus <- names(dag)[unname(node_ids) %in% ids]
  endpoints <- as.character(c(dagitty::exposures(dag), dagitty::outcomes(dag)))
  causes <- if (keep_common_causes) {
    .gd_common_causes(dag, split(focus, node_ids[focus]))
  } else {
    character()
  }
  keep <- unique(c(endpoints, focus, causes))
  ledger <- data.frame(node = names(dag), node_id = unname(node_ids),
    focus = names(dag) %in% focus, common_cause = names(dag) %in% causes,
    endpoint = names(dag) %in% endpoints,
    action = ifelse(names(dag) %in% keep, "retain", "fold"),
    stringsAsFactors = FALSE)
  reference <- .gd_canonical(dag)
  omit <- setdiff(names(dag), keep)
  graph <- if (length(omit)) {
    collapse_nodes(reference, omit)
  } else {
    reference
  }
  attr(graph, "cf_source_order") <- attr(reference, "cf_source_order", exact = TRUE)
  return(list(graph = .gd_canonical(graph), ledger = ledger,
    introduced_edges = .gd_new_edges(reference, graph), common_causes = causes))
}


.gd_time <- function(nodes) { # Read numeric time occasions from vertex names
  result <- rep(NA_real_, length(nodes))
  timed <- grepl("(?:_|@)t[0-9]+$", nodes, perl = TRUE)
  result[timed] <- as.numeric(sub("^.*(?:_|@)t([0-9]+)$", "\\1", nodes[timed],
    perl = TRUE))
  return(result)
}


.gd_edge_ledger <- function(dag) { # Record edge endpoints, time occasions and initial display roles
  edges <- dag_edges(dag)
  edges$from <- .strip_slice(edges$v)
  edges$to <- .strip_slice(edges$w)
  edges$time_from <- .gd_time(edges$v)
  edges$time_to <- .gd_time(edges$w)
  edges$display_role <- rep("shown", nrow(edges))
  return(edges)
}


# Proper endpoint paths exclude another intervention from their interior.
# Strip time only in the comparison keys. Traverse exact directed vertices.
.gd_endpoint_paths <- function(dag) { # Find proper directed endpoint paths without crossing another intervention
  treatments <- as.character(dagitty::exposures(dag))
  outcomes <- as.character(dagitty::outcomes(dag))
  edges <- dag_edges(dag)
  children <- split(edges$w[edges$e == "->"], edges$v[edges$e == "->"])
  paths <- logical()
  for (treatment in treatments) {
    blocked <- setdiff(treatments, treatment)
    reached <- character()
    pending <- treatment
    while (length(pending)) {
      node <- pending[[1L]]
      pending <- pending[-1L]
      if (node %in% reached || node %in% blocked) {
        next
      }
      reached <- c(reached, node)
      pending <- union(pending, setdiff(children[[node]], c(reached, blocked)))
    }
    keys <- paste(.strip_slice(treatment), .strip_slice(outcomes), sep = "\r")
    paths[keys] <- outcomes %in% setdiff(reached, treatment)
  }
  return(paths[order(names(paths), method = "radix")])
}


# Quantities on a common directed cycle of the time fold: the strongly
# connected components with two or more quantities. Members and groups are
# sorted by name, so the result never depends on input order.
.gd_feedback_groups <- function(from, to, source_order = NULL, node_ids = NULL) { # Find quantity groups linked by mutual directed reachability
  quantities <- unique(c(from, to))
  if (!.query_order_applies(quantities, source_order, unname(node_ids[quantities]))) {
    source_order <- NULL
  }
  quantities <- .query_sort(quantities, source_order, unname(node_ids[quantities]))
  children <- split(to, factor(from, levels = quantities))
  reach <- lapply(stats::setNames(nm = quantities), function(start) { # Collect the quantities reachable from one starting quantity
    seen <- character()
    frontier <- start
    while (length(frontier)) {
      frontier <- setdiff(unique(unlist(children[frontier], use.names = FALSE)), seen)
      seen <- c(seen, frontier)
    }
    return(seen)
  })
  groups <- lapply(quantities, function(quantity) { # Collect quantities that can return to this quantity
    returning <- Filter(function(other) { # Check for a return path from the other quantity
      return(quantity %in% reach[[other]])
    }, reach[[quantity]])
    members <- unique(c(quantity, returning))
    return(.query_sort(members, source_order, unname(node_ids[members])))
  })
  groups <- unique(groups[lengths(groups) > 1L])
  first_members <- vapply(groups, `[[`, character(1), 1L)
  return(groups[order(.query_order_key(first_members, source_order, unname(node_ids[first_members])), method = "radix")])
}


# Cut points for k contiguous segments over sorted occasions: the last occasion
# of every segment but the final one, listed in lexicographic order.
.gd_cut_candidates <- function(occasions, k) { # List candidate cut positions for contiguous temporal segments
  if (k == 1L) {
    return(list(numeric()))
  }
  closers <- occasions[-length(occasions)]
  if (length(closers) < k - 1L) {
    return(list())
  }
  index <- utils::combn(length(closers), k - 1L)
  return(lapply(seq_len(ncol(index)), function(combination) { # Read the occasion values of one candidate cut combination
    return(closers[index[, combination]])
  }))
}


# Display node of every temporal node (mapping rows). A cut quantity's
# occurrences fall into contiguous segments: a single segment keeps the
# quantity name, a segment of one occurrence keeps its exact temporal name and
# a longer one is named <quantity>_t<first>to<last>. An untimed occurrence of a
# cut quantity is never placed by guesswork. It stays a segment of its own.
.gd_segment_display <- function(mapping, cuts) { # Map temporal nodes to contiguous display segments
  display <- mapping$quantity
  time <- .gd_time(mapping$node)
  occasion_text <- function(nodes) { # Read the numeric occasion text without changing its representation
    return(sub("^.*(?:_|@)t([0-9]+)$", "\\1", nodes, perl = TRUE))
  }
  for (quantity in names(cuts)) {
    rows <- which(mapping$quantity == quantity)
    segment <- vapply(time[rows], function(occasion) { # Place this occasion in the corresponding contiguous segment
      return(if (is.na(occasion)) {
        0L
      } else {
        1L + sum(occasion > cuts[[quantity]])
      })
    }, integer(1))
    if (length(unique(segment)) < 2L) {
      next
    }
    for (segment_id in unique(segment)) {
      members <- rows[segment == segment_id]
      display[members] <- if (length(members) == 1L) {
        mapping$node[members]
      } else {
        paste0(quantity, "_t", occasion_text(mapping$node[members][which.min(time[members])]),
          "to", occasion_text(mapping$node[members][which.max(time[members])]))
      }
    }
  }
  return(display)
}


# Typed display edges once each temporal node is read as its display node.
# Edges from repeated passes stay hidden, an edge inside one display node vanishes,
# and an edge between two segments of one quantity is a real displayed edge.
.gd_display_edges <- function(links, display) { # Map typed temporal edges to distinct display edges
  kept <- links$display_role %in% c("shown", "within_quantity")
  edges <- data.frame(v = unname(display[links$v[kept]]), e = links$e[kept],
    w = unname(display[links$w[kept]]), stringsAsFactors = FALSE)
  edges <- edges[edges$v != edges$w, , drop = FALSE]
  symmetric <- edges$e == "<->"
  left <- pmin(edges$v[symmetric], edges$w[symmetric])
  right <- pmax(edges$v[symmetric], edges$w[symmetric])
  edges$v[symmetric] <- left
  edges$w[symmetric] <- right
  edges <- edges[!duplicated(edge_key(edges)), , drop = FALSE]
  rownames(edges) <- NULL
  return(edges)
}


# Kahn's algorithm over the directed rows of an edge table (as dag_is_acyclic).
.gd_edges_acyclic <- function(edges) { # Check the directed display edges by topological elimination
  directed <- edges[edges$e == "->", , drop = FALSE]
  nodes <- unique(c(directed$v, directed$w))
  indegree <- tabulate(match(directed$w, nodes), length(nodes))
  children <- split(match(directed$w, nodes), factor(directed$v, levels = nodes))
  frontier <- which(indegree == 0L)
  removed <- 0L
  while (length(frontier)) {
    node <- frontier[[1L]]
    frontier <- frontier[-1L]
    removed <- removed + 1L
    for (child in children[[node]]) {
      indegree[[child]] <- indegree[[child]] - 1L
      if (indegree[[child]] == 0L) {
        frontier <- c(frontier, child)
      }
    }
  }
  return(removed == length(nodes))
}


# Cut positions for every feedback group, chosen together. A group's fewest
# segment count depends only on its own edges, so each group first collects the
# cuts at its smallest k that leave its own part of the display acyclic (every
# member of the group shares the cuts). The groups' cuts are then combined: the
# combination with the fewest display edges wins, then the earliest cuts (all
# cut positions in ascending order, compared position by position), then the
# group order by quantity name: the first group by name takes its earliest cut,
# then the next group, and so on. `alternatives` counts the combinations tied
# on every criterion but the last, so a tie settled by name is never silent.
.gd_segment_choice <- function(links, mapping, groups) { # Choose acyclic segment cuts by the recorded rules for breaking ties
  time <- .gd_time(mapping$node)
  options <- lapply(groups, function(group) { # Find this feedback group's smallest acyclic segment alternatives
    occasions <- sort(unique(time[mapping$quantity %in% group & !is.na(time)]))
    for (k in seq_len(max(1L, length(occasions)))) {
      acyclic <- Filter(function(candidate) { # Check whether these cuts make the group's display acyclic
        display <- stats::setNames(.gd_segment_display(mapping,
          stats::setNames(rep(list(candidate), length(group)), group)), mapping$node)
        edges <- .gd_display_edges(links, display)
        members <- unique(display[mapping$quantity %in% group])
        return(.gd_edges_acyclic(edges[edges$v %in% members & edges$w %in% members, , drop = FALSE]))
      }, .gd_cut_candidates(occasions, k))
      if (length(acyclic)) {
        return(acyclic)
      }
    }
    # Only untimed occurrences or two occurrences at one occasion (A_t0 and
    # A@t0) defeat every cut. Keep each occasion apart and let the display's
    # acyclicity check refuse the summary.
    return(list(occasions[-length(occasions)]))
  })
  if (!length(options)) {
    return(list(cuts = list(), alternatives = 1L))
  }
  grid <- as.matrix(expand.grid(lapply(options, seq_along), KEEP.OUT.ATTRS = FALSE))
  cuts_of <- function(row) { # Collect the cut choice for every quantity in each group
    cuts <- list()
    for (i in seq_along(groups)) {
      cuts[groups[[i]]] <- rep(list(options[[i]][[row[[i]]]]), length(groups[[i]]))
    }
    return(cuts)
  }
  edges <- integer(nrow(grid))
  positions <- NULL
  for (combination in seq_len(nrow(grid))) {
    row <- grid[combination, ]
    display <- stats::setNames(.gd_segment_display(mapping, cuts_of(row)), mapping$node)
    edges[[combination]] <- nrow(.gd_display_edges(links, display))
    positions <- rbind(positions, sort(unlist(lapply(seq_along(groups), function(i) { # Read this group's selected cut positions
      return(options[[i]][[row[[i]]]])
    }))))
  }
  ranking <- do.call(order, c(list(edges), lapply(seq_len(ncol(positions)), function(column) { # Read this column of cut positions for lexicographic ordering
    return(positions[, column])
  }),
    lapply(seq_len(ncol(grid)), function(column) { # Read this column of group choices for the final stable tie break
      return(grid[, column])
    })))
  best <- ranking[[1L]]
  same_positions <- if (ncol(positions)) {
    rowSums(positions != matrix(positions[best, ], nrow(positions),
      ncol(positions), byrow = TRUE)) == 0
  } else {
    rep(TRUE, nrow(grid))
  }
  return(list(cuts = cuts_of(grid[best, ]), alternatives = sum(edges == edges[[best]] & same_positions)))
}


.gd_segment_cuts <- function(links, mapping, groups) { # Read the cuts chosen for the feedback groups
  return(.gd_segment_choice(links, mapping, groups)$cuts)
}


# A tie between feedback groups settled by quantity name order is reported by
# the caller that finally shows the segmented display. A caller that discards
# an applied summary (family checks, witness certification) never reports it.
.gd_segment_tie_notice <- function(summary) { # Report applied segment choices decided by quantity name order
  if (!isTRUE(summary$applied)) {
    return(invisible(summary))
  }
  parts <- c(
    if (isTRUE(summary$segment_alternatives > 1L)) {
      paste0(summary$segment_alternatives, " equally ranked cut choices")
    },
    if (isTRUE(summary$segment_edge_ties > 0L)) {
      paste0(summary$segment_edge_ties, " edge(s) from repeated stages left omitted only because an edge of the same lag and occasion came first by name")
    })
  if (length(parts)) {
    warning("Time segments: ", paste(parts, collapse = ", "), ". Quantity name order decided ",
      "(summary segment_alternatives / segment_edge_ties). A tie between feedback quantities is structural information.",
      call. = FALSE)
  }
  return(invisible(summary))
}


# Exact source rows of every temporal node a display node stands for, keyed by
# the display node: one row per represented occurrence and recorded member.
.gd_display_map <- function(trimmed, display) { # Map each display node back to every represented source row
  columns <- c("node", "quantity", "member_id", "occurrence")
  map <- attr(trimmed, "cf_node_map", exact = TRUE)
  rows <- if (is.data.frame(map) && all(columns %in% names(map))) {
    map[map$node %in% names(trimmed), columns, drop = FALSE]
  } else {
    NULL
  }
  missing <- setdiff(names(trimmed), rows$node)
  if (length(missing)) {
    ids <- unname(.gd_node_ids(trimmed)[missing])
    rows <- rbind(rows, data.frame(node = missing, quantity = ids, member_id = ids,
      occurrence = missing, stringsAsFactors = FALSE))
  }
  rows$node <- unname(display[rows$node])
  rownames(rows) <- NULL
  return(rows)
}


# One row per time segment of every quantity shown as more than one node.
.gd_segment_table <- function(mapping, source_order = NULL) { # Describe quantities displayed as multiple temporal segments
  if (!.query_order_applies(mapping$node, source_order, mapping$node_id)) {
    source_order <- NULL
  }
  pieces <- tapply(mapping$display, mapping$quantity, function(x) { # Count distinct display segments for one quantity
    return(length(unique(x)))
  })
  cut <- names(pieces)[pieces > 1L]
  if (!length(cut)) {
    return(NULL)
  }
  rows <- mapping[mapping$quantity %in% cut, , drop = FALSE]
  time <- .gd_time(rows$node)
  keys <- unique(rows[c("quantity", "display")])
  table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) { # Describe one segment and its ordered temporal members
    members <- which(rows$quantity == keys$quantity[i] & rows$display == keys$display[i])
    member_keys <- .query_order_key(rows$node[members], source_order, rows$node_id[members])
    members <- members[order(time[members], member_keys, method = "radix")]
    occasions <- time[members][!is.na(time[members])]
    return(data.frame(quantity = keys$quantity[i], segment = keys$display[i],
      first = if (length(occasions)) {
        min(occasions)
      } else {
        NA_real_
      },
      last = if (length(occasions)) {
        max(occasions)
      } else {
        NA_real_
      },
      members = paste(rows$node[members], collapse = ", "), stringsAsFactors = FALSE))
  }))
  quantity_ids <- rows$node_id[match(table$quantity, rows$quantity)]
  quantity_keys <- .query_order_key(table$quantity, source_order, quantity_ids)
  segment_keys <- .query_order_key(table$segment, source_order, quantity_ids)
  table <- table[order(quantity_keys, table$first, segment_keys, method = "radix"), , drop = FALSE]
  rownames(table) <- NULL
  return(table)
}


.gd_summary <- function(dag, keep = character(), preserve_occasions = FALSE, source_check = NULL) { # Summarise temporal routes only when all structural safeguards pass
  keep <- .gd_names(keep, "keep")
  all_nodes <- names(dag)
  all_ids <- unname(.gd_node_ids(dag))
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  if (is.null(source_order)) {
    source_order <- .cd_record(dag)$source_order
  }
  recorded_order <- source_order
  if (!.query_order_applies(all_nodes, source_order, all_ids)) {
    source_order <- NULL
  }
  missing <- setdiff(keep, c(all_nodes, all_ids))
  if (length(missing)) {
    stop("Requested summary quantities absent from this graph: ",
      paste(missing, collapse = ", "), call. = FALSE)
  }
  fallback <- function(reason) { # Retain the exact temporal graph and record why summarisation was refused
    ledger <- .gd_edge_ledger(dag)
    ledger$display_role <- rep("temporal_retained", nrow(ledger))
    return(list(graph = dag, temporal = dag,
      trim_ledger = data.frame(node = all_nodes, node_id = all_ids,
      action = rep("retain", length(all_nodes)), stringsAsFactors = FALSE),
      node_map = data.frame(node = all_nodes, quantity = .strip_slice(all_nodes),
      node_id = all_ids,
      stringsAsFactors = FALSE),
      edge_ledger = ledger, applied = FALSE, reason = reason))
  }
  if (isTRUE(preserve_occasions)) {
    return(fallback("Exact temporal occurrences retained for the joint total identification safeguard. Time merging is not certified."))
  }
  if (dag_graph_type(dag) != "dag") {
    return(fallback("Time summarisation requires a DAG with directed and bidirected edges."))
  }
  if (!isTRUE(dag_is_acyclic(dag))) {
    return(fallback("The input has a directed cycle. Its temporal display was retained."))
  }
  if (any(!dag_edges(dag)$e %in% c("->", "<->"))) {
    return(fallback("Unsupported edge types prevent time summarisation."))
  }
  if (length(dagitty::adjustedNodes(dag)) || length(selected_nodes(dag))) {
    return(fallback("Adjusted or selected statuses require a conditional projection. The temporal display was retained."))
  }
  if (any(grepl("_(pre|post)_treatment$", all_nodes))) {
    return(fallback("Numeric time occasions are required to recognise complete repeated passes. The pre/post display was retained."))
  }
  treatment <- as.character(dagitty::exposures(dag))
  outcome <- as.character(dagitty::outcomes(dag))
  if (!length(treatment) || !length(outcome)) {
    return(fallback("Treatment and outcome statuses are required for time summarisation based on routes."))
  }
  downstream <- dagitty::descendants(dag, treatment)
  upstream <- dagitty::ancestors(dag, outcome)
  routes <- unique(c(treatment, outcome, intersect(downstream, upstream)))
  exact <- intersect(keep, all_nodes)
  concept_pins <- setdiff(keep, exact)
  off_route <- setdiff(concept_pins, all_ids[all_nodes %in% routes])
  pinned <- union(exact, all_nodes[all_ids %in% off_route])
  retained <- union(routes, pinned)
  omitted <- setdiff(all_nodes, retained)
  reference <- .gd_canonical(dag)
  trimmed <- if (length(omitted)) {
    collapse_nodes(reference, omitted)
  } else {
    reference
  }
  attr(trimmed, "cf_source_order") <- attr(reference, "cf_source_order", exact = TRUE)
  mapping <- data.frame(node = names(trimmed), quantity = .strip_slice(names(trimmed)),
    node_id = unname(.gd_node_ids(trimmed)),
    stringsAsFactors = FALSE)
  # An intervention on an exact occasion cannot be replaced by an intervention on a
  # merged quantity. Equal adjustment results alone cannot certify that change.
  endpoint_ids <- unique(.strip_slice(c(treatment, outcome)))
  repeated_endpoints <- intersect(endpoint_ids,
    names(which(table(mapping$quantity) > 1L)))
  if (length(repeated_endpoints)) {
    return(fallback(paste0("Time merging would combine distinct occasions of an endpoint: ",
      paste(repeated_endpoints, collapse = ", "),
      ". The temporal display was retained.")))
  }
  latent <- names(trimmed) %in% dagitty::latents(trimmed)
  mixed <- tapply(latent, mapping$quantity, function(x) { # Detect mixed observation status within one quantity
    return(length(unique(x)) > 1L)
  })
  if (any(mixed)) {
    return(fallback(paste0("Mixed observation status across time for ",
      paste(names(mixed)[mixed], collapse = ", "), ". The temporal display was retained.")))
  }
  treatment_ids <- unique(.strip_slice(treatment))
  outcome_ids <- unique(.strip_slice(outcome))
  conflict <- intersect(treatment_ids, outcome_ids)
  if (length(conflict)) {
    return(fallback(paste0("Time merging would combine treatment and outcome: ",
      paste(conflict, collapse = ", "), ". The temporal display was retained.")))
  }

  links <- .gd_edge_ledger(trimmed)
  directed <- links[links$e == "->", , drop = FALSE]
  process_key <- function(from, time_from, to, time_to) { # Identify a directed process edge and its endpoint occasions
    return(paste(from, time_from, to, time_to, sep = "\r"))
  }
  keys <- process_key(directed$from, directed$time_from,
    directed$to, directed$time_to)
  first <- process_key(links$to, links$time_from, links$from, links$time_from)
  later <- process_key(links$to, links$time_to, links$from, links$time_to)
  repeated <- links$e == "->" & links$from != links$to &
  !is.na(links$time_from) & !is.na(links$time_to) &
  links$time_to > links$time_from & first %in% keys & later %in% keys
  repeated[is.na(repeated)] <- FALSE
  links$display_role <- ifelse(links$from == links$to, "within_quantity",
    ifelse(repeated, "between_repetitions", "shown"))
  # Fold every occasion of each quantity. Quantities on a common directed cycle
  # of that fold are cut into the fewest contiguous time segments that restore
  # acyclicity. A treatment or outcome in such a group has exactly one
  # occurrence (distinct endpoint occasions were refused above), so it is never
  # cut: its partners are segmented around it. Analysis stays on the temporal
  # graph.
  folded <- links[links$display_role == "shown" & links$e == "->", , drop = FALSE]
  quantity_ids <- stats::setNames(mapping$node_id, mapping$quantity)
  groups <- .gd_feedback_groups(folded$from, folded$to, source_order, quantity_ids)
  choice <- .gd_segment_choice(links, mapping, groups)
  mapping$display <- .gd_segment_display(mapping, choice$cuts)
  owners <- tapply(mapping$quantity, mapping$display, function(x) { # Count the quantities claiming one display label
    return(length(unique(x)))
  })
  if (any(owners > 1L)) {
    return(fallback("A time segment name coincides with another quantity. The temporal display was retained."))
  }
  segment_of <- stats::setNames(mapping$display, mapping$node)
  links$display_role[links$display_role == "within_quantity" &
    segment_of[links$v] != segment_of[links$w]] <- "shown"
  edges <- .gd_display_edges(links, segment_of)
  # Edges from repeated stages (lagged copies closing a loop within an occasion) stay
  # omitted only where showing them would close a cycle of the display. Between
  # time segments they are real edges and are shown (ledger role
  # "between_segments"). Edges are tried by tier: shorter lags first, then
  # earlier source occasions, then edge key. An edge left omitted only because
  # an edge of the same tier came first by name is a tie decided by name. It
  # is counted (segment_edge_ties) and reported by the caller that shows the
  # display.
  repeated_rows <- which(links$display_role == "between_repetitions")
  edge_ties <- 0L
  if (length(repeated_rows)) {
    tier <- paste(links$time_to[repeated_rows] - links$time_from[repeated_rows], links$time_from[repeated_rows])
    node_keys <- stats::setNames(.query_order_key(names(trimmed), source_order,
      unname(.gd_node_ids(trimmed))), names(trimmed))
    ordering <- order(links$time_to[repeated_rows] - links$time_from[repeated_rows],
      links$time_from[repeated_rows], node_keys[links$v[repeated_rows]], node_keys[links$w[repeated_rows]], method = "radix")
    repeated_rows <- repeated_rows[ordering]
    tier <- tier[ordering]
    for (level in unique(tier)) {
      before <- edges
      for (i in repeated_rows[tier == level]) {
        candidate <- data.frame(v = unname(segment_of[links$v[i]]), e = links$e[i],
          w = unname(segment_of[links$w[i]]), stringsAsFactors = FALSE)
        if (edge_key(candidate) %in% edge_key(edges)) {
          links$display_role[i] <- "between_segments"
        } else if (.gd_edges_acyclic(rbind(edges, candidate))) {
          edges <- rbind(edges, candidate)
          links$display_role[i] <- "between_segments"
        } else if (.gd_edges_acyclic(rbind(before, candidate))) {
          edge_ties <- edge_ties + 1L
        }
      }
    }
    rownames(edges) <- NULL
  }
  # A display node is unobserved only when every occurrence it stands for is.
  hidden <- tapply(latent, mapping$display, all)
  display <- construct_graph(edges, unique(mapping$display), treatment_ids,
    outcome_ids, names(hidden)[hidden])
  attr(display, "input_labels") <- attr(trimmed, "input_labels")
  attr(display, "cf_node_map") <- .gd_display_map(trimmed, segment_of)
  attr(display, "cf_source_order") <- recorded_order
  if (!isTRUE(dag_is_acyclic(display))) {
    return(fallback("Unresolved directed cycles remain after compression of repeated passes. The temporal display was retained."))
  }
  if (!identical(.gd_endpoint_paths(dag), .gd_endpoint_paths(display))) {
    return(fallback("Time merging changes a proper directed path from treatment to outcome. The temporal display was retained."))
  }
  # Dropping edges from repeated stages must not hide or invent a route role: a
  # display node (a quantity or a time segment) lies on a directed
  # route from treatment to outcome exactly when one of the occurrences it stands
  # for does in the temporal graph. Endpoints are never segmented.
  display_routes <- unique(c(treatment_ids, outcome_ids,
    intersect(dagitty::descendants(display, treatment_ids), dagitty::ancestors(display, outcome_ids))))
  source_routes <- unique(mapping$display[mapping$node %in% routes])
  if (!setequal(source_routes, display_routes)) {
    return(fallback("Time merging changes which quantities lie on a directed route from treatment to outcome. The temporal display was retained."))
  }
  # A merged node stands for every occurrence it absorbed. Each adjustment set
  # the display would report must remain valid on the source when read that way.
  if (!is.null(source_check)) {
    refusal <- source_check(.gd_canonical(display), mapping)
    if (!is.null(refusal)) {
      return(fallback(refusal))
    }
  }
  return(list(graph = .gd_canonical(display), temporal = .gd_canonical(trimmed),
    trim_ledger = data.frame(node = all_nodes, node_id = all_ids,
    route = all_nodes %in% routes, pinned = all_nodes %in% pinned,
    action = ifelse(all_nodes %in% retained, "retain", "fold"),
    stringsAsFactors = FALSE),
    node_map = mapping, edge_ledger = links, segments = .gd_segment_table(mapping, source_order),
    segment_alternatives = choice$alternatives, segment_edge_ties = edge_ties,
    trim_introduced_edges = .gd_new_edges(reference, trimmed),
    applied = TRUE, reason = NULL))
}



# Selection is grouped by distinct quantities. Graphs retain temporal evidence.
.gd_set_key <- function(x) { # Build a stable identity for a set of source quantities
  return(paste(sort(unique(x), method = "radix"), collapse = "\r"))
}


.gd_map_set <- function(set, from, to) { # Map a membership between graphs without inventing temporal equivalence
  source <- .gd_node_ids(from)
  target <- .gd_node_ids(to)
  cached_target <- .cf_acceleration_enabled() && .cf_standard_graph(to)
  target_nodes <- if (cached_target) {
    .gd_metadata(to)$nodes
  } else {
    NULL
  }
  result <- character()
  for (node in set) {
    if (node %in% if (cached_target) {
      target_nodes
    } else {
      names(to)
    }) {
      source_map <- attr(from, "cf_node_map", exact = TRUE)
      target_map <- attr(to, "cf_node_map", exact = TRUE)
      if (is.data.frame(source_map) && is.data.frame(target_map)) {
        source_times <- unique(source_map$occurrence[source_map$node == node])
        target_times <- unique(target_map$occurrence[target_map$node == node])
        if (length(source_times) == 1L && length(target_times) == 1L &&
          !identical(source_times, target_times)) {
          return(NULL)
        }
      }
      result <- union(result, node)
    } else {
      id <- unname(source[node])
      matches <- names(target)[target == id]
      # Time segments of one quantity share its id. Only the segment holding
      # this exact occurrence can represent it.
      occurrence_map <- attr(to, "cf_node_map", exact = TRUE)
      if (length(matches) > 1L && is.data.frame(occurrence_map)) {
        matches <- intersect(matches, occurrence_map$node[occurrence_map$occurrence == node])
      }
      if (length(matches) != 1L || is.na(id)) {
        return(NULL)
      }
      # A quantity summary can combine occasions. It cannot rename t2 as t3.
      source_suffix <- substring(node, nchar(.strip_slice(node)) + 1L)
      target_suffix <- substring(matches, nchar(.strip_slice(matches)) + 1L)
      source_map <- attr(from, "cf_node_map", exact = TRUE)
      target_map <- attr(to, "cf_node_map", exact = TRUE)
      if (is.data.frame(source_map) && is.data.frame(target_map)) {
        source_times <- unique(source_map$occurrence[source_map$node == node])
        target_times <- unique(target_map$occurrence[target_map$node == matches])
        if (length(source_times) == 1L && length(target_times) == 1L &&
          !identical(source_times, target_times)) {
          return(NULL)
        }
      }
      if (nzchar(source_suffix) && nzchar(target_suffix) && source_suffix != target_suffix) {
        return(NULL)
      }
      result <- union(result, matches)
    }
  }
  return(result)
}


.gd_check_family <- function(reference, final, base_sets, final_sets, effect,
  final_complete, required = list()) { # Compare source memberships and introduced alternatives on the final graph
  required_keys <- vapply(required, .gd_set_key, character(1))
  all_sets <- c(base_sets, required)
  all_keys <- vapply(all_sets, .gd_set_key, character(1))
  unique_sets <- !duplicated(all_keys)
  all_sets <- all_sets[unique_sets]
  chosen_sets <- all_keys[unique_sets] %in% required_keys
  mapped_sets <- lapply(all_sets, .gd_map_set, from = reference, to = final)
  represented <- !vapply(mapped_sets, is.null, logical(1))
  validation <- rep(FALSE, length(all_sets))
  if (any(represented)) {
    validation[represented] <- .gd_validate_sets(final, mapped_sets[represented], effect)
  }
  keys <- vapply(final_sets, .gd_set_key, character(1))
  exact_match <- rep(FALSE, length(mapped_sets))
  exact_match[represented] <- vapply(mapped_sets[represented], .gd_set_key, character(1)) %in% keys
  check_one <- function(i) { # Classify one source membership and its representation on the final graph
    raw <- all_sets[[i]]
    chosen <- chosen_sets[[i]]
    mapped <- mapped_sets[[i]]
    if (is.null(mapped)) {
      return(list(raw = raw, mapped = NULL, chosen = chosen,
        status = if (chosen) {
          "not represented"
        } else {
          "alternative omitted"
        }))
    }
    valid <- validation[[i]]
    subset_seen <- valid && !exact_match[[i]] && !isTRUE(final_complete) &&
    any(vapply(final_sets, function(candidate_set) { # Check whether an enumerated set is a proper subset of this membership
      return(length(candidate_set) < length(mapped) && all(candidate_set %in% mapped))
    }, logical(1)))
    status <- if (!valid) {
      "not valid"
    } else if (exact_match[[i]]) {
      "preserved"
    } else if (subset_seen || isTRUE(final_complete)) {
      "valid but not minimal"
    } else {
      "minimality not verified"
    }
    return(list(raw = raw, mapped = mapped, chosen = chosen, status = status))
  }
  records <- lapply(seq_along(all_sets), check_one)
  base_mapped_keys <- vapply(base_sets, function(source_set) { # Read the identity in the final graph of one source membership
    mapped <- .gd_map_set(source_set, reference, final)
    return(if (is.null(mapped)) {
      NA_character_
    } else {
      .gd_set_key(mapped)
    })
  }, character(1))
  new_sets <- final_sets[!keys %in% base_mapped_keys]
  lifted_sets <- lapply(new_sets, .gd_map_set, from = final, to = reference)
  liftable <- !vapply(lifted_sets, is.null, logical(1))
  lift_valid <- rep(FALSE, length(new_sets))
  if (any(liftable)) {
    lift_valid[liftable] <- .gd_validate_sets(reference, lifted_sets[liftable], effect)
  }
  introduced <- lapply(seq_along(new_sets), function(i) { # Classify a newly introduced final alternative on the source reference
    introduced_set <- new_sets[[i]]
    lifted <- lifted_sets[[i]]
    return(list(set = introduced_set, lifted = lifted, status = if (is.null(lifted)) {
      "temporal interpretation requires review"
    } else if
      (isTRUE(lift_valid[[i]])) {
        "valid on reference; minimality not certified"
      } else {
        "not valid on reference"
      }))
  })
  introduced <- Filter(Negate(is.null), introduced)
  changed <- any(vapply(records, function(x) { # Detect changed or omitted source memberships
    return(x$status %in% c("not valid", "valid but not minimal", "not represented", "alternative omitted"))
  }, logical(1)))
  introduced_bad <- any(vapply(introduced, function(x) { # Detect an introduced alternative without a valid reference interpretation
    return(x$status != "valid on reference; minimality not certified")
  }, logical(1)))
  uncertain <- any(vapply(records, function(x) { # Detect a membership whose minimality remains unverified
    return(x$status == "minimality not verified")
  }, logical(1)))
  return(list(effect = effect, records = records, introduced = introduced,
    status = if (changed || introduced_bad) {
      "changed"
    } else if
    (!isTRUE(final_complete) || uncertain) {
      "not fully checked"
    } else {
      "preserved"
    }))
}


.gd_final_checks <- function(built, final, chosen) { # Check both effect families for a selected graph view
  reference <- built$selection
  limit <- built$search$limit
  result <- lapply(c("direct", "total"), function(effect) { # Compare one final effect family with its reference search
    base <- if (built$search$mode == "sets" && identical(effect, built$search$effect)) {
      list(sets = built$sets, complete = built$search$complete)
    } else {
      .gd_adjustment_sets(reference, effect, limit)
    }
    final_search <- .gd_adjustment_sets(.gd_plain(final), effect, limit)
    required <- if (built$search$mode == "sets" && identical(effect, built$search$effect)) {
      chosen
    } else {
      list()
    }
    check <- .gd_check_family(reference, final, base$sets, final_search$sets, effect,
      final_search$complete, required)
    if (!isTRUE(base$complete) && identical(check$status, "preserved")) {
      check$status <- "not fully checked"
    }
    check$reference_sets <- base$sets
    check$final_sets <- final_search$sets
    check$reference_complete <- base$complete
    check$final_complete <- final_search$complete
    return(check)
  })
  return(stats::setNames(result, c("direct", "total")))
}


.gd_fallback_view <- function(built) { # Build a diagnostic endpoint view when no adjustment set exists
  focused <- .gd_focus(built$analysis, character(), keep_common_causes = FALSE)
  summary <- if (built$controls$simplify) {
    .gd_summary(focused$graph,
      source_check = .cd_summary_source_check(built$selection))
  } else {
    NULL
  }
  if (!is.null(summary)) {
    .gd_segment_tie_notice(summary)
  }
  graph <- if (is.null(summary)) {
    focused$graph
  } else {
    summary$graph
  }
  candidate <- .gd_candidates(list(character()), built$analysis, built$components, focus = TRUE)
  candidate$id <- "no_set"
  candidate$status <- "No adjustment set found"
  candidate$raw_sets <- list(list())
  candidate$full_graph_valid <- NA
  candidate$raw_full_graph_valid <- list(logical())
  candidate$display_nodes <- length(names(graph))
  candidate$temporal_nodes <- length(names(focused$graph))
  built$candidates <- candidate
  built$views <- list(no_set = list(graph = graph, temporal = focused$graph,
    focus = focused$ledger, introduced_edges = focused$introduced_edges, summary = summary))
  built$dags <- list(no_set = graph)
  built$fallback <- TRUE
  return(built)
}


.gd_print_options <- function(candidates, selected = integer()) { # Print grouped options and their differences from the selected membership
  if (!nrow(candidates)) {
    return(invisible(candidates))
  }
  for (group in unique(candidates$group)) {
    rows <- which(candidates$group == group)
    count <- candidates$n_resolved[rows[1L]]
    cat("Group ", group, " (", if (is.na(count)) {
      "node count unresolved"
    } else {
      paste(count, if (count == 1L) {
        "distinct node"
      } else {
        "distinct nodes"
      })
    }, ")\n", sep = "")
    for (i in rows) {
      tag <- paste0("  [", candidates$group[i], ",", candidates$option[i], "]")
      if (i %in% selected) {
        tag <- paste0(tag, " selected")
      }
      values <- candidates$concepts[[i]]
      if (length(selected) == 1L && i != selected) {
        prior <- candidates$concepts[[selected]]
        added_members <- setdiff(values, prior)
        removed_members <- setdiff(prior, values)
        description <- paste(c(if (length(removed_members)) {
          paste0("remove {", paste(removed_members, collapse = ", "), "}")
        },
          if (length(added_members)) {
            paste0("add {", paste(added_members, collapse = ", "), "}")
          }), collapse = ", ")
        if (!nzchar(description)) {
          description <- "same source nodes, different recorded temporal memberships"
        }
      } else {
        description <- paste0("{", paste(values, collapse = ", "), "}")
      }
      needs <- if ("unmeasured" %in% names(candidates)) {
        candidates$unmeasured[[i]]
      } else {
        character()
      }
      if (length(needs)) {
        description <- paste0(description, " needs unmeasured {", paste(needs, collapse = ", "), "}")
      }
      cat(tag, ": ", description, "\n", sep = "")
    }
  }
  return(invisible(candidates))
}


.gd_select_result <- function(built, selection, input_record) { # Choose a saved graph view and attach its checks and input record
  fallback <- !nrow(built$candidates)
  if (fallback) {
    built <- .gd_fallback_view(built)
  }
  table <- built$candidates
  # Internal callers predating grouped selection use the same explicit accounting.
  if (!all(c("group", "option") %in% names(table))) {
    sizes <- unique(table$n_resolved)
    table$group <- match(table$n_resolved, sizes)
    table$option <- ave(seq_len(nrow(table)), table$group, FUN = seq_along)
    built$candidates <- table
  }
  if (identical(selection, "ask")) {
    if (!.gd_interactive()) {
      stop("selection='ask' needs an interactive session. Supply c(group, option).", call. = FALSE)
    }
    .gd_print_options(table)
    answer <- .gd_readline("Choose group, option (for example 1, 2): ")
    selection <- .gd_answer("direct_select", if (identical(trimws(answer), "smallest")) {
      "auto"
    } else {
      answer
    })
    if (identical(selection, "auto")) {
      selection <- "smallest"
    }
    input_record$responses <- rbind(input_record$responses,
      data.frame(field = "selection", prompt = "Choose group, option", response = answer))
  }
  if (identical(selection, "smallest")) {
    index <- 1L
  } else {
    if (!is.numeric(selection) || length(selection) != 2L || anyNA(selection) ||
      any(!is.finite(selection)) || any(selection < 1 | selection != floor(selection))) {
      stop("selection must be c(group, option), 'smallest', or 'ask'.", call. = FALSE)
    }
    index <- which(table$group == selection[1L] & table$option == selection[2L])
    if (length(index) != 1L) {
      stop("That [group, option] is absent from the saved options.", call. = FALSE)
    }
  }
  id <- table$id[index]
  graph <- .gd_attach(built$views[[id]]$graph, .gd_initial_state(built, id))
  state <- .gd_state(graph)
  pair <- as.integer(c(table$group[index], table$option[index]))
  input_record$requested_selection <- input_record$values$selection
  input_record$values$selection <- pair
  input_record$resolved_call$selection <- pair
  state$input_record <- input_record
  state$selected <- pair
  state$catalogue <- built
  state$checks <- .gd_final_checks(built, graph, table$raw_sets[[index]])
  state$checks_signature <- .gd_signature(graph)
  state$selection_base_signature <- .gd_signature(graph)
  graph <- .gd_attach(graph, state)
  call_args <- input_record$values
  treatments <- if (!is.null(call_args$treatments)) {
    call_args$treatments
  } else {
    dagitty::exposures(built$query)
  }
  outcomes <- if (!is.null(call_args$outcomes)) {
    call_args$outcomes
  } else {
    dagitty::outcomes(built$query)
  }
  cat("Causal query: ", paste(c(paste0(treatments, " (treatment)"),
    paste0(outcomes, " (outcome)")), collapse = ", "), "\n", sep = "")
  if (fallback) {
    warning("No adjustment set found. Returning the endpoint graph. This is not an empty adjustment set.", call. = FALSE)
  } else {
    .gd_print_options(table, index)
  }
  if (anyNA(table$n_resolved)) {
    warning("Some quantity counts remain unresolved. See dag_options().", call. = FALSE)
  }
  if (identical(built$search$complete, FALSE)) {
    warning("Options are the sets found within the supplied search limit. The search limit may omit alternatives or smaller sets.", call. = FALSE)
  }
  statuses <- vapply(state$checks, `[[`, character(1), "status")
  cat("Adjustment checks: ", paste(paste(names(statuses), statuses), collapse = ", "), ".\n", sep = "")
  if (any(statuses != "preserved")) {
    warning("Adjustment sets changed or were not fully verified after recovery. Inspect dag_checks().", call. = FALSE)
  }
  if (nrow(input_record$responses)) {
    message("Your responses and the resolved call were saved to the output object. Use dag_call().")
  }
  return(graph)
}



#####
# get_dag_search.R
#####
# Full construction searches use Inf. Finite private searches are needed only
# for a certificate witness and retained tests of the separator contract.
# The public dagitty 0.3-4 direct branch ignores max.results, so its disposable
# context receives the same limit in the underlying separator routine.
.gd_adjustment_sets <- function(dag, effect, max_sets = Inf, progress = NULL) { # Enumerate adjustment sets and record any finite search boundary
  effect <- match.arg(effect, c("total", "direct"))
  if (!is.numeric(max_sets) || length(max_sets) != 1L || is.na(max_sets) ||
    max_sets < 1 || (is.finite(max_sets) &&
    (max_sets != floor(max_sets) || max_sets >= .Machine$integer.max))) {
    stop("The private search limit must be Inf or a positive whole number below .Machine$integer.max.",
      call. = FALSE)
  }
  limit <- if (is.finite(max_sets)) {
    as.integer(max_sets) + 1L
  } else {
    Inf
  }
  dag <- dagitty::as.dagitty(dag)
  if (!length(dagitty::exposures(dag)) || !length(dagitty::outcomes(dag))) {
    stop("Both treatment and outcome tags are required to enumerate adjustment sets.",
      call. = FALSE)
  }
  raw <- .gd_adjustment_engine(dag, effect, limit, progress)
  sets <- unname(lapply(raw, function(x) { # Represent one returned adjustment set as a character vector
    return(as.character(unlist(x, use.names = FALSE)))
  }))
  complete <- length(sets) <= max_sets
  probe <- if (!complete) {
    sets[[as.integer(max_sets) + 1L]]
  } else {
    NULL
  }
  if (!complete) {
    sets <- sets[seq_len(as.integer(max_sets))]
  }
  # A capped prefix is not guaranteed to contain the globally smallest sets.
  # An empty set remains list(character(0)). No sets remains list().
  return(list(sets = sets, complete = complete, probe = probe,
    limit = if (is.finite(max_sets)) {
      as.integer(max_sets)
    } else {
      Inf
    }))
}


# Validate exact sets with time indices. This does not claim inclusion minimality.
# In particular, direct validation is NOT dagitty::isAdjustmentSet(), which
# checks total effect adjustment. No separator family is enumerated here.
.gd_validate_sets <- function(dag, sets, effect) { # Validate exact memberships with the verified total or direct criterion
  effect <- match.arg(effect, c("total", "direct"))
  dag <- dagitty::as.dagitty(dag)
  if (!length(dagitty::exposures(dag)) || !length(dagitty::outcomes(dag))) {
    stop("Both treatment and outcome tags are required to validate adjustment sets.",
      call. = FALSE)
  }
  if (!is.list(sets) || any(!vapply(sets, function(z) { # Check that this membership contains character names without missing values
    return(is.character(z) && !anyNA(z))
  }, logical(1)))) {
    stop("sets must be a list of character vectors without missing values.",
      call. = FALSE)
  }
  if (!length(sets)) {
    return(logical())
  }
  endpoints <- c(dagitty::exposures(dag), dagitty::outcomes(dag))
  eligible <- vapply(sets, function(z) { # Check source membership and endpoint and latent exclusions
    return(all(z %in% names(dag)) && !any(z %in% endpoints) && !any(z %in% dagitty::latents(dag)))
  }, logical(1))
  if (effect == "total") {
    return(unname(vapply(seq_along(sets), function(i) { # Validate this eligible set for total effect adjustment
      return(eligible[[i]] && isTRUE(dagitty::isAdjustmentSet(dag, sets[[i]])))
    }, logical(1))))
  }
  if (utils::packageVersion("dagitty") != "0.3.4") {
    stop("Direct set validation currently supports dagitty 0.3-4 only. ",
      "The installed version is ", as.character(utils::packageVersion("dagitty")),
      ".", call. = FALSE)
  }
  script <- system.file("js", "dagitty-alg.js", package = "dagitty")
  if (!nzchar(script)) {
    stop("The installed dagitty JavaScript engine is unavailable.", call. = FALSE)
  }
  context <- V8::v8()
  context$source(script)
  context$eval(paste0(
    "var gd_supported = typeof DAGitty !== 'undefined' && ",
    "typeof DAGitty.GraphParser.parseGuess === 'function' && ",
    "typeof DAGitty.GraphTransformer.indirectGraph === 'function' && ",
    "typeof DAGitty.GraphAnalyzer.listMsasDirectEffect === 'function' && ",
    "typeof DAGitty.GraphAnalyzer.listMinimalSeparators === 'function' && ",
    "typeof DAGitty.GraphAnalyzer.dConnected === 'function' && ",
    "DAGitty.GraphAnalyzer.listMinimalSeparators.length === 4;"
    ))
  if (!isTRUE(context$get("gd_supported"))) {
    stop("The installed dagitty engine lacks the verified API for validating direct sets.",
      call. = FALSE)
  }
  context$assign("gd_graph_text", as.character(dag))
  context$assign("gd_sets", lapply(sets, as.list))
  context$eval(paste(
    "var gd_valid = (function() {",
    "  var g = DAGitty.GraphParser.parseGuess(gd_graph_text);",
    "  var analyzer = DAGitty.GraphAnalyzer;",
    "  var constraints = null;",
    # Ask the original direct routine for its required/forbidden variables.
    # This hook captures the contract and returns without enumerating any sets.
    "  analyzer.listMinimalSeparators = function(graph, must, mustNot, max) {",
    "    constraints = {",
    "      must: must.map(function(v) { return v.id; }),",
    "      mustNot: mustNot.map(function(v) { return v.id; })",
    "    };",
    "    return [];",
    "  };",
    "  analyzer.listMsasDirectEffect(g, [], []);",
    "  if (constraints === null) return gd_sets.map(function() { return false; });",
    # Include every candidate in the conditioning operation. Restricting the
    # check to the original endpoint ancestor graph could miss a collider
    # opened by an additional variable outside that ancestor graph.
    "  var indirect = DAGitty.GraphTransformer.indirectGraph(g);",
    "  return gd_sets.map(function(z) {",
    "    if (!constraints.must.every(function(id) { return z.indexOf(id) >= 0; }) ||",
    "        constraints.mustNot.some(function(id) { return z.indexOf(id) >= 0; }) ||",
    "        z.some(function(id) { return !g.getVertex(id); })) return false;",
    "    return !analyzer.dConnected(indirect, indirect.getSources(),",
    "      indirect.getTargets(), z);",
    "  });",
    "})();",
    sep = "\n"
    ))
  valid <- as.logical(unlist(context$get("gd_valid", simplifyVector = FALSE),
    use.names = FALSE))
  return(unname(eligible & valid))
}



# Graph code is the public object. Reconstruction and edit records are attributes.
# Stored graphs never contain this provenance again.
.gd_plain <- function(graph) { # Remove construction provenance from a stored graph value
  attr(graph, "cf_provenance") <- NULL
  class(graph) <- setdiff(class(graph), "cf_dag")
  return(graph)
}


.gd_signature <- function(graph) { # Record graph structure and statuses for provenance checks
  edges <- dag_edges(graph)
  symmetric <- edges$e %in% c("<->", "--")
  if (any(symmetric)) {
    left_endpoint <- pmin(edges$v[symmetric], edges$w[symmetric])
    right_endpoint <- pmax(edges$v[symmetric], edges$w[symmetric])
    edges$v[symmetric] <- left_endpoint
    edges$w[symmetric] <- right_endpoint
  }
  labels <- attr(graph, "input_labels")
  if (length(labels)) {
    labels <- labels[order(names(labels))]
  }
  return(list(type = dag_graph_type(graph), nodes = sort(names(graph)),
    edges = sort(unique(edge_key(edges))),
    exposures = sort(as.character(dagitty::exposures(graph))),
    outcomes = sort(as.character(dagitty::outcomes(graph))),
    latent = sort(as.character(dagitty::latents(graph))),
    adjusted = sort(as.character(dagitty::adjustedNodes(graph))),
    selected = sort(as.character(selected_nodes(graph))), labels = labels))
}


.gd_provenance <- function(dag) { # Read and validate the graph's construction provenance
  state <- attr(dag, "cf_provenance", exact = TRUE)
  if (is.null(state) || !is.list(state) || is.null(state$original_query)) {
    stop("This graph has no construction provenance. Start from get_dag() or construct_dag() to use traced commands.", call. = FALSE)
  }
  return(state)
}


.gd_state <- function(dag) { # Refuse traced operations when graph structure no longer matches its receipt
  if (!inherits(dag, "dagitty")) {
    stop("Choose one DAG before applying this command.", call. = FALSE)
  }
  state <- .gd_provenance(dag)
  if (!identical(.gd_signature(dag), state$signature)) {
    stop("The DAG was changed outside the traced commands. Its saved reference no longer matches. Reconstruct it before restoring nodes. Changes to coordinates alone are allowed.", call. = FALSE)
  }
  return(state)
}


.gd_attach <- function(graph, state) { # Attach construction state and a fresh structural signature
  if (!inherits(graph, "dagitty")) {
    stop("graph must inherit dagitty.", call. = FALSE)
  }
  graph <- .gd_plain(graph)
  state$checkpoint <- graph
  state$signature <- .gd_signature(graph)
  attr(graph, "cf_provenance") <- state
  class(graph) <- unique(c("cf_dag", class(graph)))
  return(graph)
}


.gd_delta <- function(before, after) { # Describe changes to nodes, edges and causal statuses
  before_sig <- .gd_signature(before)
  after_sig <- .gd_signature(after)
  statuses <- c("exposures", "outcomes", "latent", "adjusted", "selected")
  status_changes <- lapply(statuses, function(name) { # Record the added and removed nodes for one causal status
    return(list(
      added = setdiff(after_sig[[name]], before_sig[[name]]),
      removed = setdiff(before_sig[[name]], after_sig[[name]])))
  })
  names(status_changes) <- statuses
  status_changes <- status_changes[vapply(status_changes,
    function(x) { # Keep only statuses with an actual membership change
      return(length(x$added) + length(x$removed) > 0L)
    }, logical(1))]
  return(list(added_nodes = setdiff(names(after), names(before)),
    removed_nodes = setdiff(names(before), names(after)),
    added_edges = .gd_new_edges(before, after),
    removed_edges = .gd_new_edges(after, before),
    statuses = status_changes))
}


.gd_transition <- function(before, after, operation, arguments = list(), details = list()) { # Record a transformation with its before and after graphs
  return(list(operation = operation, arguments = arguments, delta = .gd_delta(before, after),
    before = .gd_canonical(.gd_plain(before)),
    after = .gd_canonical(.gd_plain(after)), details = details))
}


.gd_record <- function(before, after, operation, arguments = list(), details = list(), undo = TRUE) { # Append a traced transformation while retaining the original query
  previous_state <- .gd_state(before)
  state <- attr(after, "cf_provenance", exact = TRUE)
  if (is.null(state)) {
    state <- previous_state
  }
  if (isTRUE(undo)) {
    snapshot <- previous_state[setdiff(names(previous_state), c("undo", "history"))]
    state$undo <- c(previous_state$undo, list(list(graph = .gd_plain(before), state = snapshot)))
  }
  # A command may rebuild state, but it cannot replace the original query.
  state$original_query <- previous_state$original_query
  state$history <- c(previous_state$history,
    list(.gd_transition(before, after, operation, arguments, details)))
  return(.gd_attach(after, state))
}


.gd_preview <- function(before, after, operation, details = list()) { # Print and return an uncommitted preview of graph changes
  .gd_state(before)
  result <- structure(list(operation = operation, delta = .gd_delta(before, after),
    before = .gd_plain(before), after = after, details = details), class = "cf_dag_change")
  print(result)
  return(invisible(result))
}


.gd_initial_state <- function(built, id) { # Build provenance and history for one constructed graph view
  view <- built$views[[id]]
  summary <- view$summary
  temporal <- if (is.null(summary)) {
    view$temporal
  } else {
    summary$temporal
  }
  recovery <- built$recovery
  build <- built[setdiff(names(built), c("query", "analysis", "dags", "views", "recovery"))]
  build$recovery <- recovery[setdiff(names(recovery), "graph")]
  history <- list(
    .gd_transition(built$query, recovery$graph, "recover_nodes",
    list(requested = recovery$requested_ids), recovery[setdiff(names(recovery), "graph")]),
    .gd_transition(recovery$graph, built$analysis, "set_availability",
    list(mode = built$availability$mode), built$availability),
    .gd_transition(built$analysis, view$temporal, "focus_graph",
    list(focus = built$candidates$concepts[[match(id, built$candidates$id)]]),
    list(ledger = view$focus, introduced_edges = view$introduced_edges)))
  if (!is.null(summary)) {
    history <- c(history, list(
      .gd_transition(view$temporal, summary$temporal, "trim_time_routes", details = list(
      ledger = summary$trim_ledger, introduced_edges = summary$trim_introduced_edges)),
      .gd_transition(summary$temporal, summary$graph, "summarise_time", details = list(
      applied = summary$applied, reason = summary$reason, node_map = summary$node_map,
      edge_ledger = summary$edge_ledger))))
  }
  return(list(original_query = built$query, build = build, view_id = id,
    analysis = built$analysis, temporal = temporal, retained = names(temporal),
    summary = summary, view = view[c("focus", "introduced_edges")],
    focus = built$candidates$concepts[[match(id, built$candidates$id)]],
    controls = built$controls, observed_ids = built$observed_ids,
    history = history, undo = list()))
}


.gd_graph_result <- function(built) { # Attach provenance to each constructed graph view
  graphs <- lapply(names(built$views), function(id) { # Attach the saved initial state to this graph view
    return(.gd_attach(built$views[[id]]$graph, .gd_initial_state(built, id)))
  })
  names(graphs) <- names(built$views)
  if (length(graphs) == 1L) {
    return(graphs[[1L]])
  }
  class(graphs) <- c("cf_dags", "list")
  attr(graphs, "cf_provenance") <- list(original_query = built$query,
    build = built[setdiff(names(built), c("query", "dags", "views"))])
  return(graphs)
}



# Internal identity and counts of distinct quantities used by DAG construction.
.gd_resolve_focus <- function(focus, water_model, disambiguate, strict = TRUE) { # Resolve focus requests against the active alias vocabulary
  aliases <- getExportedValue("causaliflower", "node_aliases")
  nodes <- getExportedValue("causaliflower", "apsim_dag_nodes")
  if (!is.null(.query_source_order(nodes)) && identical(water_model, "SWIM")) {
    stop("water_model = \"SWIM\" was requested, but the bundled slice has no SWIM family rows. Use water_model = \"SoilWater\".", call. = FALSE)
  }
  inactive <- if (water_model == "SoilWater") {
    "SWIM."
  } else {
    "SoilWater."
  }
  aliases <- aliases[!startsWith(names(aliases), inactive)]
  found <- .resolve_columns(focus, .build_alias_index(aliases), disambiguate,
    hidden_ids = nodes$node_id[nodes$hidden_node == "TRUE"])
  if (length(found$ambiguous) || (strict && length(found$unmatched))) {
    stop("Unresolved focus: ", paste(c(names(found$ambiguous), found$unmatched), collapse = ", "), call. = FALSE)
  }
  return(unique(unname(found$resolved)))
}


.gd_components <- function(x) { # Validate declared component members and their supporting records
  if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))) {
    stop("components must be a named list with unique canonical concept IDs.", call. = FALSE)
  }
  for (name in names(x)) {
    item <- x[[name]]
    if (!is.list(item) || !all(c("members", "reason", "source") %in% names(item)) ||
      !is.character(item$members) || anyNA(item$members) || any(!nzchar(item$members)) ||
      !is.character(item$reason) || length(item$reason) != 1L || is.na(item$reason) || !nzchar(trimws(item$reason)) ||
      !is.character(item$source) || length(item$source) != 1L || is.na(item$source) || !nzchar(trimws(item$source))) {
      stop("Each components entry needs character members, a reason and a source: ", name, call. = FALSE)
    }
  }
  return(x)
}


.gd_resolve_components <- function(ids, registry, source_order = NULL) { # Expand component memberships while recording cycles and unresolved definitions
  trace <- list()
  visit <- function(id, path = character()) { # Resolve one component and record the traversal path
    item <- registry[[id]]
    if (is.null(item)) {
      return(list(members = id, unresolved = character()))
    }
    repeated <- id %in% path
    status <- if (repeated) {
      "cycle"
    } else if (!length(item$members)) {
      "unresolved"
    } else {
      "expanded"
    }
    trace[[length(trace) + 1L]] <<- data.frame(concept = id, path = paste(c(path, id), collapse = " -> "),
      status = status, reason = item$reason, source = item$source)
    if (status != "expanded") {
      return(list(members = character(), unresolved = paste0("unresolved:", id)))
    }
    children <- lapply(unique(item$members), visit, path = c(path, id))
    return(list(members = unique(unlist(lapply(children, `[[`, "members"), use.names = FALSE)),
      unresolved = unique(unlist(lapply(children, `[[`, "unresolved"), use.names = FALSE))))
  }
  results <- lapply(ids, visit)
  return(list(members = .query_sort(as.character(unique(unlist(lapply(results, `[[`, "members"), use.names = FALSE))), source_order),
    unresolved = sort(as.character(unique(unlist(lapply(results, `[[`, "unresolved"), use.names = FALSE))), method = "radix"),
    trace = if (length(trace)) {
      do.call(rbind, trace)
    } else {
      data.frame()
    }))
}


.gd_candidates <- function(raw, dag, registry, focus = FALSE, progress = NULL, source_order = NULL) { # Group source sets by their distinct quantities and resolved components
  mapping <- attr(dag, "input_labels")
  if (missing(source_order)) {
    source_order <- attr(dag, "cf_source_order", exact = TRUE)
    if (is.null(source_order)) {
      source_order <- .cd_record(dag)$source_order
    }
  }
  base_id <- function(x) { # Read canonical base identities using the graph's saved labels
    x <- .strip_slice(x)
    has <- x %in% names(mapping)
    x[has] <- unname(mapping[x[has]])
    return(x)
  }
  if (!is.null(source_order)) {
    family_ids <- unique(base_id(unlist(raw, use.names = FALSE)))
    resolved_ids <- .gd_resolve_components(family_ids, registry)$members
    if (!.query_order_applies(c(unname(.gd_node_ids(dag)), family_ids, resolved_ids), source_order)) {
      source_order <- NULL
    }
  }
  sets <- lapply(raw, function(x) { # Sort the distinct canonical quantities in one source set
    return(.query_sort(unique(base_id(x)), source_order))
  })
  keys <- vapply(sets, function(x) { # Build the grouping key for one canonical quantity set
    return(paste(x, collapse = "\r"))
  }, character(1))
  groups <- unname(split(seq_along(keys), match(keys, unique(keys))))
  result <- data.frame(id = character(), group = integer(), option = integer(), group_size = integer(),
    n_concepts = integer(), n_resolved = integer(),
    status = character(), label = character(), stringsAsFactors = FALSE)
  result$concepts <- list()
  result$raw_sets <- list()
  result$resolved <- list()
  result$unresolved <- list()
  result$resolution <- list()
  rows <- vector("list", length(groups))
  for (index in seq_along(groups)) {
    source_indices <- groups[[index]]
    ids <- sets[[source_indices[1L]]]
    resolution <- .gd_resolve_components(ids, registry, source_order)
    status <- if (length(resolution$unresolved)) {
      "needs definition"
    } else if (nrow(resolution$trace)) {
      "resolved count only"
    } else {
      "candidate"
    }
    row <- data.frame(id = "", group = 0L, option = 0L, group_size = NA_integer_, n_concepts = length(ids),
      n_resolved = if (length(resolution$unresolved)) {
        NA_integer_
      } else {
        length(resolution$members)
      },
      status = if (focus && status == "candidate") {
        "requested focus"
      } else {
        status
      },
      label = if (length(ids)) {
        paste(ids, collapse = " + ")
      } else {
        "No adjustment nodes"
      }, stringsAsFactors = FALSE)
    row$concepts <- list(ids)
    row$raw_sets <- list(raw[source_indices])
    row$resolved <- list(resolution$members)
    row$unresolved <- list(resolution$unresolved)
    row$resolution <- list(resolution$trace)
    rows[[index]] <- row
    if (index %% 256L == 0L) {
      .gd_notify(progress,
        paste0("Grouped ", index, " of ", length(groups), " alternatives by source quantity"), force = FALSE)
    }
  }
  if (length(rows)) {
    result <- do.call(rbind, rows)
  }
  if (nrow(result)) {
    catalogue_keys <- result$label
    if (!is.null(source_order)) {
      catalogue_keys <- vapply(result$concepts, .query_set_key, character(1), source_order = source_order)
    }
    result <- result[order(is.na(result$n_resolved), result$n_resolved, catalogue_keys,
      method = "radix"), , drop = FALSE]
    result$group <- match(result$n_resolved, unique(result$n_resolved))
    result$option <- as.integer(ave(seq_len(nrow(result)), result$group, FUN = seq_along))
    result$group_size <- result$n_resolved
    result$id <- paste0(if (focus) {
      "focus_"
    } else {
      "set_"
    }, seq_len(nrow(result)))
    rownames(result) <- NULL
  }
  return(result)
}
