# ==============================================================================
# dagby_query.R


# Query graphs from the crop structural causal model
# ==============================================================================


#' Query the source graph for a recipe
#'
#' Builds the source graph described by a [dagby()] recipe from the bundled crop
#' structural causal model derived from APSIM. Inspect or reuse this graph before
#' passing it to [dagby_construct()] to construct a DAG of the requested size.
#'
#' @param recipe A `dagby_recipe` created by [dagby()].
#' @param quiet Suppress the build summary when `TRUE`.
#'
#' @details
#' Set the crop, treatments, outcomes and query controls in [dagby()]. The
#' recipe's construction settings (`size`, `keep_nodes`, `latents` and `unroll`)
#' are carried forward for [dagby_construct()] and do not change this source query.
#'
#' Requested endpoints are retained, are not recombined, and receive exposure or
#' outcome tags. With `randomised_treatment = TRUE`, incoming directed arrows are
#' removed from every treatment occasion in ordinary queries. Queries anchored
#' to growth stages remove incoming arrows only from their designated treatment
#' occasions. Use this only when the treatment was randomised. Otherwise it can
#' remove genuine confounding paths and invalidate the adjustment sets. Removed
#' arrows and the assumption are saved.
#'
#' The recipe's `label = "id"` keeps quantity identifiers, with suffixes for
#' distinct time occasions. `"input"` uses resolved input names, including
#' supplied display names. `"readable"` uses display labels with the same input
#' override, falling back to quantity keys when labels are unavailable or conflict.
#' The legacy `"APSIM"` mode retains dotted
#' keys where available and otherwise follows `"input"`. The label to quantity
#' map is saved.
#'
#' The `time_unroll` setting in [dagby()] controls the temporal representation.
#' When omitted, it uses `"all"` if a treatment and lagged edges are present,
#' otherwise `"none"`. Explicit `NULL` means `"none"`. Temporal modes:
#'
#' - `"none"`: Keep one occurrence of each quantity and represent lagged arrows
#'   without separate time occasions. This can produce cycles. A warning reports
#'   collapsed lags and cycles. A cyclic result cannot be used as a causal DAG.
#' - `"all"`: Keep separate occasions for every included quantity within the
#'   finite temporal window. Requires treatments.
#' - `"selective"`: Use the same temporal window, retaining separate occasions
#'   for quantities with occurrences both affected and unaffected by treatment.
#'   Additional quantities retain separate occasions if combining them leaves
#'   a cycle. Requires treatments.
#' - `"pre_post"`: Initially group occurrences according to whether they are
#'   descendants of treatment, using `_pre_treatment` and `_post_treatment`
#'   suffixes. If grouping leaves cycles, the function removes intermediate
#'   occasions while retaining directed paths through them, or restores separate
#'   `_t{s}` occasions where needed. Requires treatments.
#' - `"collapse"`: Construct the `"all"` representation, then combine each
#'   quantity's occasions. The result may contain cycles and is not used for adjustment.
#' - Logical `FALSE` and `TRUE` retain their earlier meanings of `"none"` and
#'   `"all"`. Treatment and outcome each remain anchored to one occasion.
#'
#' The temporal representation does not itself establish identification. Check
#' adjustment sets for the resulting graph. The finite window uses recorded lags
#' between treatment and outcome, with preceding history where required. Paths
#' outside that window are not recovered. Growth stage queries instead carry
#' each lagged arrow forward one developmental phase. This is a discretisation
#' choice, not a claim that every process takes a whole phase.
#'
#' A hidden quantity is mechanistic detail, not necessarily unobserved. In
#' [dagby_query_params()], `show_hidden_nodes = FALSE` removes eligible hidden
#' quantities while retaining directed paths through them. This operation does
#' not introduce bidirected arrows. Explicitly requested quantities are retained.
#' Recovery policy exclusions use the same directed operation regardless of this option.
#'
#' The recipe's `columns` resolve management inputs and identify observed
#' quantities. Included quantities not matched by those columns are marked latent,
#' subject to the endpoint and recovery rules. Without `columns`, the default
#' `observe_apsim_outputs = TRUE` treats ordinary included quantities as obtainable.
#' This is a working assumption, not verified availability for every APSIM
#' output. With `FALSE`, included hidden quantities are marked latent.
#' Supplied columns take precedence over these automatic
#' observation tags. Query endpoints retain their observed status.
#'
#' The recovery policy distinguishes `yes`, quantities admitted only on explicit
#' request (`yes:user`), and exclusions (`no:reason`). Explicit endpoints or
#' columns can admit `yes:user`. Requesting an excluded quantity raises an error.
#' The current interpretation policy was reviewed for barley. Recorded provenance
#' outside APSIM applies across crops. Observation assumptions, conflicts, output
#' mappings and removed graph structure are saved with the query.
#'
#' Quantities left without any arrows are dropped unless requested as
#' treatments, outcomes or columns. Other occasions of an outcome in a growth stage query
#' remain distinct causal states. They are not automatically outcomes,
#' interventions or observed measurements.
#'
#' The bundled data support only SoilWater and individual readable quantity keys.
#' SWIM, coarser dotted grouping and supernode overrides raise errors. Supplied
#' columns and endpoint names must resolve to the bundled keys or aliases.
#'
#' @return A dagitty object with its recipe, query record, resolved quantities,
#' edges, roles, effective settings and source ordering metadata. The records
#' describe the supplied recipe and query, including any deliberately restricted
#' scope. Construction does not silently replace it with a fuller query. Use the
#' returned `edges` for the current graph. Saved source structures may describe
#' the graph before excluded quantities were removed.
#' @examples
#' \dontrun{
#' recipe <- dagby(crop = "wheat", treatments = "N_fertiliser_applied",
#'                 outcomes = "grain_wt", randomised_treatment = FALSE)
#' graph <- dagby_query(recipe)
#' dag <- dagby_construct(graph)
#'
#' # Anchor treatment and outcome to developmental phases
#' recipe <- dagby(crop = "barley", treatments = "N_fertiliser_applied",
#'                 outcomes = "grain_protein", randomised_treatment = FALSE,
#'                 params = dagby_query_params(treatment_growth_stage = "flowering",
#'                                             outcome_growth_stage = "harvest_ripe"))
#' graph <- dagby_query(recipe)
#' }
#' @export
dagby_query <- function(recipe, quiet = FALSE){
  recipe <- .dagby_validate_recipe(recipe)
  .gd_flag(quiet, "quiet")
  graph <- do.call(".dagby_query_source", c(recipe$query, list(quiet = quiet)),
                   envir = asNamespace("causaliflower"))
  query_call <- as.call(list(quote(causaliflower::dagby_query), .dagby_recipe_call(recipe), quiet = quiet))
  record <- attr(graph, "query_record", exact = TRUE)
  record$call <- query_call
  attr(graph, "query_record") <- record
  attr(graph, "dagby_record") <- list(recipe = recipe, resolved_call = query_call)

  return(graph)
}


#' Build the source graph from supplied query arguments
#'
#' @inheritParams dagby
#' @param label NULL keeps source model identifiers. Other supported modes are
#'   `"input"`, `"readable"` and `"APSIM"`.
#' @param quiet Suppress the build summary when `TRUE`.
#' @return The source graph with its scientific query record and effective settings.
#' @noRd
.dagby_query_source <- function(crop,
                        columns              = NULL,
                        treatments           = NULL,
                        outcomes             = NULL,
                        randomised_treatment = FALSE,
                        time_unroll          = NULL,
                        label                = NULL,
                        quiet                = FALSE,
                        params               = dagby_query_params()
                        ){

  .datatable.aware <- TRUE
  query_call <- match.call(expand.dots = FALSE)
  query_input <- list(crop = crop, columns = columns, treatments = treatments,
                      outcomes = outcomes, randomised_treatment = randomised_treatment,
                      time_unroll = time_unroll, label = label, quiet = quiet, params = params)
  # Capture omission before resolving time_unroll. Omission selects the
  # default frame later, whereas explicit NULL means no unrolling.
  time_unroll_defaulted <- missing(time_unroll)
  # Bind package data explicitly so manual and wrapper calls cannot read a
  # table with the same name from the caller's workspace.
  nodes <- getExportedValue("causaliflower", "apsim_dag_nodes")
  source_order <- .query_source_order(nodes)
  aliases <- getExportedValue("causaliflower", "node_aliases")
  management_defaults <- getExportedValue("causaliflower", "management_defaults")
  crop_phases <- getExportedValue("causaliflower", "crop_phases")
  phase_of <- getExportedValue("causaliflower", "phase_of")

  query_depth <- params$query_depth
  water_model <- params$water_model
  apsim_nodes_detail <- params$apsim_nodes_detail
  supernode_overrides <- params$supernode_overrides
  force_acyclic <- params$force_acyclic
  show_hidden_nodes <- params$show_hidden_nodes
  observe_apsim_outputs <- if( is.null(params$observe_apsim_outputs) ){ TRUE }else{ params$observe_apsim_outputs }
  node_recombination <- params$node_recombination
  disambiguate <- params$disambiguate
  treatment_growth_stage <- params$treatment_growth_stage
  outcome_growth_stage <- params$outcome_growth_stage
  phase_window <- params$phase_window
  if( !is.null(source_order) ){
    if( identical(water_model, "SWIM") ){
      stop("water_model = \"SWIM\" was requested, but the bundled slice has no SWIM family rows. Use water_model = \"SoilWater\".", call. = FALSE)
    }
    if( !identical(apsim_nodes_detail, "variable") || length(supernode_overrides) > 0L ){
      stop("The bundled slice has readable variable keys without dotted groups. Use apsim_nodes_detail = \"variable\" without supernode_overrides.", call. = FALSE)
    }
  }
  # The default NULL preserves readable node keys. Input labels and readable labels
  # remain explicit choices. The legacy APSIM mode preserves dotted identifiers when present.
  if( !is.null(label) ){
    label <- match.arg(label, c("APSIM", "input", "readable"))
  }
  crops <- sort(unique(unlist(strsplit(nodes$crop_type[ nodes$crop_type != "agnostic" ], ";", fixed = TRUE))))
  if( !crop %in% crops ){
    stop(paste0("crop must be one of: ", paste(crops, collapse = ", "), ". Please supply a supported crop and try again."))
  }

  # Hidden node IDs cannot break an alias tie as a uniquely observed column
  hidden_ids <- nodes$node_id[ nodes$hidden_node == "TRUE" ]

  # Resolve aliases only within the selected water model so an inactive model
  # cannot claim an input name
  inactive_prefix <- if( identical(water_model, "SoilWater") ){ "SWIM." }else{ "SoilWater." }
  aliases_active <- aliases[ !startsWith(names(aliases), inactive_prefix) ]
  alias_index <- .build_alias_index(aliases_active)
  # Normalise explicit time modes now. Leave omitted time_unroll as NULL until
  # the retained lagged graph is known. The mode definitions belong to
  # dagby_query() help.
  if( is.logical(time_unroll) ){
    time_unroll <- if( isTRUE(time_unroll) ){ "all" }else{ "none" }
  }else if( time_unroll_defaulted ){
    time_unroll <- NULL                                    # sentinel: resolve to all / none below
  }else if( is.null(time_unroll) ){
    time_unroll <- "none"                                  # explicit NULL -> none
  }else{
    time_unroll <- match.arg(time_unroll, c("none", "all", "selective", "collapse", "pre_post"))
  }
  if( !is.null(time_unroll) && time_unroll %in% c("selective", "all", "collapse", "pre_post") ){
    if( length(treatments) == 0L ){
      stop(sprintf("time_unroll = \"%s\" needs at least one treatment: the unrolling window is sized from the treatment to outcome lag distance, and the split into occasions before and after treatment is defined relative to the treatment. Please pass treatments = and try again.", time_unroll))
    }
    # show_hidden_nodes = TRUE is supported: the slicer skips its internal .collapse_hidden and
    # hidden occasions are carried through to the emitted graph, marked [latent] downstream.
  }

  # Phase anchored (developmental time) mode: triggered whenever treatment_growth_stage or
  # outcome_growth_stage is supplied (both now via dagby_query_params()). It runs the absolute phase selective
  # frame (.selective_phase_slice, with B0 at the sowing origin) instead of the frame anchored
  # to treatment. When both growth stage params are NULL, dagby_query uses time_unroll.
  # dagby_query_params() validates phase_window, which is used only here.
  phase_mode <- !is.null(treatment_growth_stage) || !is.null(outcome_growth_stage)
  if( phase_mode ){
    if( length(treatments) == 0L ){
      stop("A query anchored to developmental phases needs at least one treatment: treatment_growth_stage / outcome_growth_stage anchor the treatment and outcome to developmental phases. Please pass treatments = and try again.")
    }
    # show_hidden_nodes = TRUE is supported: the phase slicer skips its internal .collapse_hidden
    # and hidden occasions are carried through to the emitted graph, marked [latent] downstream.
    if( !time_unroll_defaulted && !identical(time_unroll, "selective") ){
      stop("A query anchored to developmental phases uses the selective frame in developmental time. Leave time_unroll unset (the default) or set it to \"selective\", not \"none\" / \"all\" / \"collapse\" / \"pre_post\". Please adjust time_unroll and try again.")
    }
    # the phase frame is selective in developmental time: resolve the sentinel so downstream
    # references and the emitted attribute read "selective".
    time_unroll <- "selective"
    if( !exists("crop_phases") ){
      stop("A query anchored to developmental phases needs the crop_phases table, which is not loaded. Build the node layer (its mobilise step writes crop_phases.rda) and load it, then try again.")
    }
    if( !(crop %in% crop_phases$crop) ){
      stop(sprintf("crop '%s' has no crop_phases partition, so growth stages cannot be resolved (partitioned crops: %s). Please supply a crop with a phase partition and try again.",
                   crop, paste(sort(unique(crop_phases$crop)), collapse = ", ")))
    }
  }
  has_columns <- length(columns) > 0L
  col_matched <- character(0L)
  col_res <- list(resolved = character(0L), ambiguous = list(), unmatched = character(0L))
  if( has_columns ){
    col_res <- .resolve_columns(columns, alias_index, disambiguate, hidden_ids)
    col_matched <- unique(unname(col_res$resolved))
    if( length(col_res$ambiguous) > 0L ){
      warning(paste0("These columns are ambiguous and were left unresolved: ",
                     paste(names(col_res$ambiguous), collapse = ", "),
                     ". Rename them to a node_id or set params = dagby_query_params(disambiguate = ...) in dagby() to resolve them."))
    }
    if( !is.null(source_order) && length(col_res$unmatched) > 0L ){
      stop("These columns match no node in the bundled slice: ",
           paste(col_res$unmatched, collapse = ", "),
           ". Use a readable node key or an unambiguous alias.", call. = FALSE)
    }
  }

  # treatment / outcome node ids (resolved by alias or supplied directly as node_ids)
  treat_ids <- intersect(unique(c(treatments, .match_columns(treatments, alias_index))), nodes$node_id)
  out_ids <- intersect(unique(c(outcomes,   .match_columns(outcomes,   alias_index))), nodes$node_id)
  if( !is.null(source_order) ){
    unresolved_endpoints <- character()
    for( endpoint in as.character(unlist(c(treatments, outcomes))) ){
      if( !endpoint %in% nodes$node_id && !length(.match_columns(endpoint, alias_index)) ){
        unresolved_endpoints <- c(unresolved_endpoints, endpoint)
      }
    }
    if( length(unresolved_endpoints) > 0L ){
      stop("These treatments or outcomes resolve to no single node in the bundled slice: ",
           paste(unique(unresolved_endpoints), collapse = ", "),
           ". Use a readable node key or an unambiguous alias.", call. = FALSE)
    }
  }

  # Resolve display labels separately from structural node selection.
  # Treatment and outcome labels take precedence over column labels, with the
  # first input retained within each source. A node grouped at a coarser
  # detail keeps its group label.
  input_pairs <- character(0L)                              # node_id -> input string
  if( !is.null(label) ){
    .add_first <- function( pairs, resolved ){ # Keep the first resolved input label for each node
      # resolved: NAMED input -> node_id. Accumulate node_id -> input, first writer wins.
      for( inp in names(resolved) ){
        nid <- unname(resolved[[ inp ]])
        if( !( nid %in% names(pairs) ) ){
          pairs[[ nid ]] <- inp
        }
      }
      return(pairs)
    }
    col_pairs <- if( has_columns ){ .add_first(character(0L), col_res$resolved) }else{ character(0L) }
    to_pairs <- .add_first(character(0L), .resolve_columns(treatments, alias_index, disambiguate, hidden_ids)$resolved)
    to_pairs <- .add_first(to_pairs,      .resolve_columns(outcomes,   alias_index, disambiguate, hidden_ids)$resolved)
    input_pairs <- col_pairs
    for( nid in names(to_pairs) ){                          # treatments / outcomes override columns
      input_pairs[[ nid ]] <- to_pairs[[ nid ]]
    }
    # label = "input" shows the supplied labels: the names in columns = c(display = token),
    # or the token itself when the element is unnamed, and likewise for treatments / outcomes.
    # These labels need not be registered aliases. Keep a custom name even when it cannot be
    # resolved back through the alias map. .sanitise_label
    # strips any character that would corrupt the emitted dagitty text, and the uniqueness guard in
    # .apply_input_labels reverts any label collision back to the node_id.
  }

  resolved <- .resolve_edges(nodes, crop = crop)
  if( !.is_acyclic(nodes$node_id, resolved$contemporaneous) ){
    stop("The contemporaneous structural graph must be acyclic. Please check function_of / composition_of in the node table and try again.")
  }

  present <- .resolve_management_presence(nodes, management_defaults$core, c(treat_ids, out_ids), col_matched, alias_index)
  kept_ids <- .filter_nodes(nodes, crop, water_model, present)

  # Group dotted IDs at the requested detail level, then redirect their edges.
  # Overrides choose a different level for named first-segment groups. Grouping
  # can create cycles, so the contracted graph is checked below.
  group_map <- .apsim_truncate(nodes$node_id, apsim_nodes_detail, overrides = supernode_overrides)
  g_kept <- unique(unname(group_map[ kept_ids ]))
  g_contemp <- .contract_edges(resolved$contemporaneous[ from %in% kept_ids & to %in% kept_ids ], group_map)
  g_lagged <- .contract_edges(resolved$lagged[ from %in% kept_ids & to %in% kept_ids ], group_map)

  # Resolve treatments / outcomes to the underlying node_id SET (stable across any later
  # refinement), so a supernode named as treatment / outcome maps to whatever group label(s)
  # it ends up as. treat_ids / out_ids are the direct node_id matches. A supernode name
  # additionally pulls in every node_id its current label covers (strict: only resolves at
  # the detail level where that group label exists).
  first_segment <- function( ids ){ # Extract the leading segment of each node identifier
    return(sapply(strsplit(ids, ".", fixed = TRUE), `[`, 1L))
  }
  group_labels <- unique(unname(group_map))
  treat_named <- intersect(unique(c(treatments, .match_columns(treatments, alias_index))), group_labels)
  out_named <- intersect(unique(c(outcomes,   .match_columns(outcomes,   alias_index))), group_labels)
  treat_id_set <- unique(c(treat_ids, names(group_map)[ unname(group_map) %in% treat_named ]))
  out_id_set <- unique(c(out_ids,   names(group_map)[ unname(group_map) %in% out_named ]))

  # Recovery eligibility is structural, independent of observation tags.
  # Requests permit yes:user, never no:reason. Keep the original row policy
  # in the receipt, including rows which will disappear from the returned DAG.
  output_policy <- .query_output_policy(nodes, crop)
  requested_ids <- unique(c(treat_id_set, out_id_set, col_matched))
  prohibited <- output_policy$node_id[startsWith(output_policy$recovered_node, "no:")]
  forbidden_requests <- intersect(requested_ids, prohibited)
  if( length(forbidden_requests) ){
    reasons <- output_policy$recovered_node[match(forbidden_requests, output_policy$node_id)]
    stop(paste0("Cannot recover requested node(s): ",
                paste(paste0(forbidden_requests, " (", reasons, ")"), collapse = ", ")), call. = FALSE)
  }
  excluded_ids <- union(prohibited, setdiff(
    output_policy$node_id[output_policy$recovered_node == "yes:user"], requested_ids))

  # A bare supernode name supplied at a finer than supernode level resolves to nothing
  # (strict). Flag that case so the user is not left wondering why.
  all_supers <- unique(first_segment(nodes$node_id))
  unmatched_super <- setdiff(intersect(c(treatments, outcomes), all_supers),
                             unique(c(treat_named, out_named, unname(group_map[ c(treat_ids, out_ids) ]))))
  if( length(unmatched_super) > 0L ){
    message(paste0("Treatment / outcome name(s) ", paste(unmatched_super, collapse = ", "),
                   " are supernodes, but apsim_nodes_detail = \"", apsim_nodes_detail,
                   "\" is finer than \"supernode\" so they are not single nodes. Name a node id, ",
                   "set apsim_nodes_detail = \"supernode\", or pass supernode_overrides to recover them."))
  }

  # supernodes we prefer not to split when forcing acyclicity: those named as treatment /
  # outcome, and the present management levers (all keyed by first segment).
  pinned_supers <- unique(first_segment(c(treat_id_set, out_id_set, present)))

  if( !.is_acyclic(g_kept, g_contemp) ){
    if( isTRUE(force_acyclic) ){
      auto_overrides <- .refine_to_acyclic(nodes$node_id, kept_ids, resolved$contemporaneous,
                                                apsim_nodes_detail, supernode_overrides, pinned_supers)
      supernode_overrides <- auto_overrides
      group_map <- .apsim_truncate(nodes$node_id, apsim_nodes_detail, overrides = supernode_overrides)
      g_kept <- unique(unname(group_map[ kept_ids ]))
      g_contemp <- .contract_edges(resolved$contemporaneous[ from %in% kept_ids & to %in% kept_ids ], group_map)
      g_lagged <- .contract_edges(resolved$lagged[ from %in% kept_ids & to %in% kept_ids ], group_map)
      if( isFALSE(quiet) ){
        split_txt <- if( length(auto_overrides) > 0L ){ paste(names(auto_overrides), "->", unname(auto_overrides), collapse = ", ") }else{ "none" }
        message(paste0("force_acyclic split supernode(s) to remove the cycle: ", split_txt, "."))
      }
    }else{
      cyc <- merge(g_contemp, g_contemp, by.x = c("from", "to"), by.y = c("to", "from"))
      cyc <- unique(cyc[ from < to, list(from, to) ])
      eg <- if( nrow(cyc) > 0L ){ paste0(" (e.g. \"", cyc$from[ 1 ], "\" <-> \"", cyc$to[ 1 ], "\")") }else{ "" }
      warning(paste0("apsim_nodes_detail = \"", apsim_nodes_detail, "\" contracts the graph into a cycle between supernodes", eg,
                     ", so the result is not a DAG. It is still returned. To obtain an acyclic graph you can: set force_acyclic = TRUE ",
                     "to refine the groups to an acyclic representation, or select a finer apsim_nodes_detail (\"node\" or \"subnode\")."))
    }
  }

  # A coarse label must not combine an excluded interpretation with allowed
  # quantities (for example spike and Spike.*). Refine only those mixed groups.
  group_members <- split(kept_ids, unname(group_map[kept_ids]))
  mixed <- vapply(group_members, function(ids) { # Find groups containing both excluded and admitted quantities
    return(any(ids %in% excluded_ids) && any(!ids %in% excluded_ids))
  }, logical(1L))
  refined_ids <- unlist(group_members[mixed], use.names = FALSE)
  if( length(refined_ids) ){
    group_map[refined_ids] <- refined_ids
    g_kept <- unique(unname(group_map[kept_ids]))
    g_contemp <- .contract_edges(resolved$contemporaneous[from %in% kept_ids & to %in% kept_ids], group_map)
    g_lagged <- .contract_edges(resolved$lagged[from %in% kept_ids & to %in% kept_ids], group_map)
  }
  g_excluded <- unique(unname(group_map[intersect(excluded_ids, kept_ids)]))
  g_requested <- unique(unname(group_map[intersect(requested_ids, kept_ids)]))

  # grouped names of the treatment / outcome nodes (after refinement), used to pin and tag
  # them below. Mapping the stable set of ids through the final group_map means a named
  # supernode that force_acyclic had to split is recovered as ALL of its split pieces, each of
  # which is then pinned and tagged: the multiple treatment / multiple outcome representation
  # required once an aggregate cannot be a single acyclic node.
  treat_grouped <- unique(unname(group_map[ intersect(treat_id_set, names(group_map)) ]))
  out_grouped <- unique(unname(group_map[ intersect(out_id_set,   names(group_map)) ]))

  # a grouped node is hidden only when every node it groups is hidden (hidden_ids computed
  # above). It counts as matched (has data) when any node it groups matched a supplied column.
  g_hidden <- names(which(tapply(nodes$node_id %in% hidden_ids, unname(group_map[ nodes$node_id ]), all)))
  g_matched <- unique(unname(group_map[ col_matched ]))
  output_mappings <- .query_apsim_outputs(crop)
  output_ids <- if( isTRUE(observe_apsim_outputs) ){
    output_policy$node_id[output_policy$recovered_node == "yes"]
  }else{ character(0L) }
  # A grouped node is eligible only if every represented member is admitted.
  g_outputs <- names(which(tapply(nodes$node_id %in% output_ids,
                                 unname(group_map[nodes$node_id]), all)))

  # resolve the defaulted sentinel now that g_lagged is available: if there are lagged edges in the
  # kept subgraph and treatments are specified, the default is "all": an uncompressed, acyclic
  # frame anchored to treatment within the adaptive finite window. Otherwise "none".
  # (An explicit NULL was already resolved to "none" above. Only the unset default reaches here as NULL.)
  if( !phase_mode && is.null(time_unroll) ){
    if( nrow(g_lagged) > 0L && length(treat_grouped) > 0L ){
      time_unroll <- "all"
      if( isFALSE(quiet) ){
        message("Lagged edges detected in kept subgraph. Defaulting to time_unroll = \"all\" (uncompressed occasions in the query's finite time window). Pass another mode explicitly for a compressed or static view.")
      }
    }else{
      time_unroll <- "none"
    }
  }

  # Path dependent stocks (persistent_ids) carry a temporal succession edge in selective mode
  # (the earlier occasion causes the later one). Read from a path_dependent column if the node
  # table carries one. Absent that column, no succession edges are emitted (graceful fallback).
  persistent_ids <- if( "path_dependent" %in% names(nodes) ){
    pd_nodes <- nodes$node_id[ toupper(as.character(nodes$path_dependent)) == "TRUE" ]
    unique(c(pd_nodes, unname(group_map[ intersect(pd_nodes, names(group_map)) ])))
  }else{
    character(0L)
  }

  if( phase_mode ){
    # developmental phase selective frame. Resolve each treatment / outcome to its anchored phase
    # index(es). Every treatment must be anchored, an unanchored outcome defaults to the terminal
    # (measurement) phase. B0 flows through via the contracted init edges.
    g_initial <- .contract_edges(resolved$initial[ from %in% kept_ids & to %in% kept_ids ], group_map)
    n_ph <- max(crop_phases$phase_index[ crop_phases$crop == crop ]) + 1L
    ts_map <- .stages_to_phase_indices(.normalize_stage_spec(treatment_growth_stage, treat_grouped), crop, crop_phases)
    os_map <- .stages_to_phase_indices(.normalize_stage_spec(outcome_growth_stage,   out_grouped),   crop, crop_phases)
    unanchored_t <- setdiff(treat_grouped, names(ts_map))
    if( length(unanchored_t) > 0L ){
      stop(sprintf("Every treatment must be given a growth stage in a query anchored to developmental phases. Missing: %s. Pass treatment_growth_stage as a single stage for all, or a named list per treatment, and try again.",
                   paste(unanchored_t, collapse = ", ")))
    }
    for( o in setdiff(out_grouped, names(os_map)) ){ os_map[[ o ]] <- n_ph - 1L }
    pa <- if( exists("phase_of") && is.data.frame(phase_of) && nrow(phase_of) > 0L ){
      .resolve_phase_active(phase_of, crop, crop_phases)
    }else{
      NULL
    }
    sliced <- .selective_phase_slice(g_kept, g_contemp, g_lagged, g_initial, g_hidden,
                                     ts_map, os_map, n_ph, persistent_ids = persistent_ids, window = phase_window,
                                     phase_active = pa, show_hidden_nodes = TRUE)
  }else if( time_unroll %in% c("selective", "all", "collapse", "pre_post") ){
    # Same slicer, same window, same scope. "all" is the UNCOMPRESSED frame: every included base is
    # split at every slice rather than only the straddlers, so no lagged edge lacks an occasion to
    # land on. "selective" compresses bases that do not straddle treatment back to a single bare node. "collapse" runs the
    # SAME uncompressed frame as "all" (so every lagged edge has an occasion to land on) and defers the
    # static projection to a final step in the dagby_query tail. It is uncompressed here for the same
    # reason "all" is. "pre_post" is the straddler split of "selective", but each straddler collapses to
    # two occasions (<id>_pre_treatment / <id>_post_treatment by descendant status) instead of one per
    # slice. label_mode = "pre_post" selects that inside the slicer, with the promote loop keeping it
    # acyclic.
    sliced <- .selective_slice(g_kept, g_contemp, g_lagged, g_hidden, treat_grouped, out_grouped,
                               persistent_ids,
                               split_all  = time_unroll %in% c("all", "collapse"),
                               label_mode = if( identical(time_unroll, "pre_post") ){ "pre_post" }else{ "slice" },
                               show_hidden_nodes = TRUE)
  }else{
    sliced <- .contemporaneous_slice(g_kept, g_contemp, g_lagged)
    # Combining time occasions can turn feedback into a cycle. Self-loops are
    # omitted and duplicate arrows combined, so this view does not preserve
    # every temporal edge. Report the result before any adjustment analysis.
    if( isTRUE(sliced$n_collapsed > 0L) ){
      warning(sprintf(paste0("time_unroll = \"none\" does not unroll time, so %d lagged edge%s been ",
                             "COLLAPSED onto the single slice. Self-loops are omitted and duplicate ",
                             "arrows are combined. The result %s. Collapsing a feedback loop (e.g. the ",
                             "fertiliser decision rule reads the soil N it changes) onto one slice makes ",
                             "it a cycle. Use time_unroll = \"selective\" (or \"all\") to retain time ",
                             "occasions, then check adjustment validity on the resulting graph."),
                     sliced$n_collapsed,
                     if( sliced$n_collapsed == 1L ){ " has" }else{ "s have" },
                     if( isTRUE(sliced$acyclic) ){ "happens to be acyclic here" }else{ "IS CYCLIC" }),
              call. = FALSE)
    }
  }

  hidden_vertex <- stats::setNames(.strip_slice(sliced$vertices) %in% g_hidden, sliced$vertices)

  # Fold after persistence and treatment surgery, so paths through temporal
  # occurrences are retained too. Slicers above deliberately defer this step.
  recovered_vertices <- sliced$vertices
  recovered_edges <- sliced$edges
  # Phase endpoints are selected occasions, not every occurrence of their base.
  # Histories of an outcome remain ordinary states with their own identities.
  phase_treatments <- if( phase_mode ){ sliced$treatments }else{ character(0L) }
  phase_outcomes <- if( phase_mode ){ sliced$outcomes }else{ character(0L) }
  requested_occurrences <- if( phase_mode ){
    union(c(phase_treatments, phase_outcomes),
          recovered_vertices[.strip_slice(recovered_vertices) %in% g_matched])
  }else{ recovered_vertices[.strip_slice(recovered_vertices) %in% g_requested] }

  # Endpoints, requested columns and present management levers survive scoping
  # and recombination. The final isolated-node removal also preserves explicit
  # endpoints and columns.
  pin_ids <- unique(c(treat_grouped, out_grouped, g_requested, unname(group_map[ present ])))
  pinned <- recovered_vertices[ .strip_slice(recovered_vertices) %in% pin_ids ]
  if( phase_mode ){
    pinned <- union(requested_occurrences,
                    recovered_vertices[.strip_slice(recovered_vertices) %in% unname(group_map[present])])
  }

  # Temporal succession edges (selective mode): the earlier occasion of a path dependent stock
  # causally precedes the later one (forward Euler identity: stock(t+1) = stock(t) + net_flux(t)).
  # Add these before scoping so the ancestral scope retains the baseline occasion of a stock
  # on the path X_pre -> X_post -> ... -> outcome. A stock outside that path whose later occasion
  # reaches nothing (X_t1) is still pruned because it is not an ancestor, rather than kept as an orphan sink.
  # Direction is always earlier -> later, so these cannot create a cycle.
  if( (phase_mode || time_unroll %in% c("selective", "all", "collapse", "pre_post")) && !is.null(sliced$persistence) && nrow(sliced$persistence) > 0L ){
    p <- sliced$persistence
    hit <- p[ from %in% recovered_vertices & to %in% recovered_vertices ]
    if( nrow(hit) > 0L ){
      recovered_edges <- unique(rbind(recovered_edges, hit))
    }
  }

  # Randomisation removes incoming edges to the selected treatment occasions
  # before recording or scoping the graph. Both operations must see the graph
  # under that design assumption. Cycles that do not pass through treatment
  # are unaffected.
  randomised_cut <- character(0L)
  if( isTRUE(randomised_treatment) && length(treat_grouped) > 0L && nrow(recovered_edges) > 0L ){
    into_treat <- if( phase_mode ){ recovered_edges$to %in% phase_treatments }else{
      .strip_slice(recovered_edges$to) %in% treat_grouped
    }
    if( any(into_treat) ){
      randomised_cut <- sort(unique(.strip_slice(recovered_edges$from[ into_treat ])))
      recovered_edges <- recovered_edges[ !into_treat ]
      if( isFALSE(quiet) ){
        message(sprintf("randomised_treatment = TRUE: removed %d edge%s into the treatment (do() surgery, parents cut: %s). This is only valid if the treatment is randomised.",
                        sum(into_treat), if( sum(into_treat) == 1L ){ "" }else{ "s" }, paste(randomised_cut, collapse = ", ")))
      }
    }
  }

  source_node_map <- .query_identity_rows(recovered_vertices, group_map, kept_ids)
  excluded_occurrences <- recovered_vertices[.strip_slice(recovered_vertices) %in% g_excluded]
  hidden_occurrences <- if( isTRUE(show_hidden_nodes) ){ character(0L) }else{
    setdiff(recovered_vertices[.strip_slice(recovered_vertices) %in% g_hidden], requested_occurrences)
  }
  fold_nodes <- union(excluded_occurrences, hidden_occurrences)
  recovery_folding <- .query_fold_nodes(recovered_vertices, recovered_edges, fold_nodes)
  recovery_folding$policy_removed <- excluded_occurrences
  recovery_folding$hidden_removed <- setdiff(hidden_occurrences, excluded_occurrences)
  recovery_folding$refined_members <- refined_ids
  recovered_vertices <- recovery_folding$vertices
  recovered_edges <- recovery_folding$edges
  pinned <- intersect(pinned, recovered_vertices)

  full_vertices <- recovered_vertices
  full_edges <- recovered_edges
  # query_depth starts at 1: 1 = the adjustment scope (.scope_adjust), 2 = the ancestral
  # subgraph. Each further step adds a generation of parents and children (.scope_to_query with
  # graph_order = query_depth - 2), Inf = the full graph (no scoping). Scoping needs a query.
  if( is.finite(query_depth) && length(treat_grouped) > 0L && length(out_grouped) > 0L ){
    treat_v_pre <- if( phase_mode ){ intersect(recovered_vertices, phase_treatments) }else{
      recovered_vertices[.strip_slice(recovered_vertices) %in% treat_grouped]
    }
    out_v_pre <- if( phase_mode ){ intersect(recovered_vertices, phase_outcomes) }else{
      recovered_vertices[.strip_slice(recovered_vertices) %in% out_grouped]
    }
    scoped <- if( query_depth == 1L ){
      .scope_adjust(recovered_vertices, recovered_edges, treat_v_pre, out_v_pre, pinned, quiet)
    }else{
      .scope_to_query(recovered_vertices, recovered_edges, treat_v_pre, out_v_pre,
                      pinned, query_depth - 2L, quiet)
    }
    recovered_vertices <- scoped$vertices
    recovered_edges <- scoped$edges
    pinned <- intersect(pinned, recovered_vertices)
  }

  latent_set <- .assign_latent(recovered_vertices, hidden_vertex, has_columns, g_matched, show_hidden_nodes,
                                always_observed = if( phase_mode ){ character(0L) }else{ unique(c(treat_grouped, out_grouped)) },
                                apsim_observed_ids = union(g_outputs, if( phase_mode ){ g_matched }else{ g_requested }))
  if( phase_mode ){ latent_set <- setdiff(latent_set, c(phase_treatments, phase_outcomes)) }

  # Record identity before graph grouping or display relabelling. One emitted
  # vertex may represent several source rows. Retain every member explicitly.
  identity_rows <- .query_identity_rows(recovered_vertices, group_map, kept_ids)
  recombined <- list()
  if( isTRUE(node_recombination) ){
    recomb <- .esc_recombine(recovered_vertices, recovered_edges, latent_set, pinned)
    recovered_vertices <- recomb$vertices
    recovered_edges <- recomb$edges
    latent_set <- recomb$latent
    recombined <- recomb$supernodes
    for( combined in names(recombined) ){
      take <- identity_rows$node %in% recombined[[ combined ]]
      identity_rows$node[ take ] <- combined
      identity_rows$quantity[ take ] <- .strip_slice(combined)
      identity_rows$occurrence[ take ] <- combined
    }
  }

  # The collapse mode projects the completed temporal query onto variable
  # names after scoping and recombination. Self-loops disappear and parallel
  # edges are deduplicated. Feedback between variables can reappear as a
  # directed cycle, so this view is not necessarily a DAG.
  if( identical(time_unroll, "collapse") ){
    identity_rows$node <- .strip_slice(identity_rows$node)
    ef <- .strip_slice(recovered_edges$from)
    et <- .strip_slice(recovered_edges$to)
    keep <- ef != et                                       # drop collapsed self-loops (succession/persistence)
    recovered_edges <- unique(data.table::data.table(from = ef[ keep ], to = et[ keep ]))
    recovered_vertices <- unique(.strip_slice(recovered_vertices))
    latent_set <- intersect(unique(.strip_slice(latent_set)), recovered_vertices)
    pinned <- intersect(unique(.strip_slice(pinned)),     recovered_vertices)
  }

  # Final step: drop orphaned nodes (no incoming or outgoing edge) from the output. This is a
  # documented internal decision, no parameter. Management levers are admitted and pinned into
  # scope so their causal edges are recovered. Explicit endpoints and columns
  # are exempt, so their requested representation is not silently discarded.
  connected <- unique(c(recovered_edges$from, recovered_edges$to))
  endpoint_pins <- intersect(recovered_vertices, requested_occurrences)
  if( identical(time_unroll, "collapse") ){
    endpoint_pins <- recovered_vertices[.strip_slice(recovered_vertices) %in% g_requested]
  }
  recovered_vertices <- intersect(recovered_vertices, union(connected, endpoint_pins))
  latent_set <- intersect(latent_set, recovered_vertices)
  identity_rows <- identity_rows[ identity_rows$node %in% recovered_vertices, , drop = FALSE ]

  # Remove a temporal suffix when only one occasion of that variable remains
  # and its bare name is free. Keep the recorded source identity through the
  # relabelling.
  suffixed <- recovered_vertices[ recovered_vertices != .strip_slice(recovered_vertices) ]
  occ_count <- table(.strip_slice(suffixed))
  lone_base <- names(occ_count)[ occ_count == 1L ]
  relabel <- stats::setNames(recovered_vertices, recovered_vertices)
  for( v in suffixed ){
    bare <- .strip_slice(v)
    if( ( bare %in% lone_base ) && !( bare %in% recovered_vertices ) ){
      relabel[[ v ]] <- bare
    }
  }
  if( any(relabel != names(relabel)) ){
    identity_rows$node <- unname(relabel[ identity_rows$node ])
    recovered_vertices <- unname(relabel[ recovered_vertices ])
    recovered_edges <- unique(data.table::data.table(from = unname(relabel[ recovered_edges$from ]),
                                                       to   = unname(relabel[ recovered_edges$to ])))
    latent_set <- unique(unname(relabel[ latent_set ]))
    pinned <- unname(relabel[ pinned ])
  }

  treatments_v <- if( phase_mode ){ intersect(recovered_vertices, unname(relabel[phase_treatments])) }else{
    recovered_vertices[.strip_slice(recovered_vertices) %in% treat_grouped]
  }
  outcomes_v <- if( phase_mode ){ intersect(recovered_vertices, unname(relabel[phase_outcomes])) }else{
    recovered_vertices[.strip_slice(recovered_vertices) %in% out_grouped]
  }

  # These legacy rolled attributes describe the source process graph. A folded
  # path can cross multiple time steps, so filtering source arrows to retained
  # endpoints would silently omit its replacement. Preserve the source view
  # explicitly whenever folding occurred. The returned DAG is in `edges`.
  .pg_bases <- if( length(fold_nodes) ){ unique(.strip_slice(recovery_folding$source_vertices)) }else{ unique(.strip_slice(recovered_vertices)) }
  .pg_contemp <- g_contemp[ from %in% .pg_bases & to %in% .pg_bases ]
  .pg_lagged <- g_lagged[  from %in% .pg_bases & to %in% .pg_bases ]
  .pg_stocks <- intersect(persistent_ids, .pg_bases)
  .pg_treatments <- unique(.strip_slice(treatments_v))
  .pg_outcomes <- unique(.strip_slice(outcomes_v))
  .pg_latent <- unique(.strip_slice(latent_set))

  # label= feature: relabel the recovered graph so a resolved input alias becomes the node's display
  # label (see .apply_input_labels). Applied to vertices / edges / latent / exposure / outcome sets
  # BEFORE .emit_dagitty, so the emitted text, the dagitty object, and the *_v attributes all carry
  # the display labels. input_labels (display label -> node_id) is used to restore node identities.
  # The NULL setting skips relabelling entirely, leaving the emitted text unchanged.
  input_labels <- character(0L)
  if( !is.null(label) ){
    original_vertices <- recovered_vertices
    labelled_graph <- .apply_input_labels(recovered_vertices, recovered_edges, latent_set,
                                             treatments_v, outcomes_v, nodes, input_pairs, label)
    recovered_vertices <- labelled_graph$V
    recovered_edges <- labelled_graph$E
    latent_set <- labelled_graph$lat
    treatments_v <- labelled_graph$X
    outcomes_v <- labelled_graph$Y
    input_labels <- labelled_graph$input_labels
    display_map <- stats::setNames(recovered_vertices, original_vertices)
    identity_rows$node <- unname(display_map[ identity_rows$node ])
  }

  dag_text <- .emit_dagitty(recovered_vertices, recovered_edges, latent_set, treatments_v, outcomes_v)

  dag <- dag_text
  if( requireNamespace("dagitty", quietly = TRUE) ){
    dag <- dagitty::dagitty(dag_text)
  }
  attr(dag, "dag_text") <- dag_text
  attr(dag, "crop") <- crop
  attr(dag, "water_model") <- water_model
  attr(dag, "latent") <- latent_set
  attr(dag, "recovered_vertices") <- recovered_vertices
  attr(dag, "edges") <- recovered_edges
  attr(dag, "treatments") <- treatments_v
  attr(dag, "outcomes") <- outcomes_v
  attr(dag, "show_hidden_nodes") <- isTRUE(show_hidden_nodes)
  attr(dag, "observe_apsim_outputs") <- isTRUE(observe_apsim_outputs)
  attr(dag, "apsim_output_mappings") <- output_mappings[
    output_mappings$node_id %in% identity_rows$member_id, , drop = FALSE]
  present_policy <- output_policy[output_policy$node_id %in% source_node_map$member_id, , drop = FALSE]
  attr(dag, "apsim_output_assumption") <- list(
    applied = isTRUE(observe_apsim_outputs) && !has_columns,
    basis = "Ordinary recovery policy nodes presumed obtainable. Recorded exceptions excluded. Not an exhaustive output audit.",
    policy = present_policy,
    conflicts = present_policy[present_policy$metadata_conflict, , drop = FALSE])
  attr(dag, "time_unroll") <- if( phase_mode ){ "selective" }else{ time_unroll }
  attr(dag, "phase_window") <- if( phase_mode ){ phase_window }else{ NA_character_ }
  attr(dag, "apsim_nodes_detail") <- apsim_nodes_detail
  attr(dag, "supernode_overrides") <- supernode_overrides
  attr(dag, "force_acyclic") <- isTRUE(force_acyclic)
  attr(dag, "node_recombination") <- isTRUE(node_recombination)
  attr(dag, "randomised_treatment") <- isTRUE(randomised_treatment)
  attr(dag, "randomised_parents_cut") <- randomised_cut
  attr(dag, "query_depth") <- query_depth
  # Save the unscoped graph separately, with latent status assigned over its
  # full vertex set
  attr(dag, "full_vertices") <- full_vertices
  attr(dag, "full_edges") <- full_edges
  attr(dag, "full_latent") <- .assign_latent(full_vertices, hidden_vertex, has_columns, g_matched,
                                                    show_hidden_nodes,
                                                    always_observed = if( phase_mode ){ character(0L) }else{ unique(c(treat_grouped, out_grouped)) },
                                                    apsim_observed_ids = union(g_outputs, if( phase_mode ){ g_matched }else{ g_requested }))
  if( phase_mode ){
    attr(dag, "full_latent") <- setdiff(attr(dag, "full_latent"), c(phase_treatments, phase_outcomes))
  }
  attr(dag, "recombined_nodes") <- recombined
  attr(dag, "input_labels") <- input_labels
  attr(dag, "contemporaneous_edges") <- .pg_contemp
  attr(dag, "lagged_edges") <- .pg_lagged
  attr(dag, "persistent_stocks") <- .pg_stocks
  attr(dag, "rolled_treatments") <- .pg_treatments
  attr(dag, "rolled_outcomes") <- .pg_outcomes
  attr(dag, "rolled_latent") <- .pg_latent
  attr(dag, "rolled_structure") <- if( length(fold_nodes) ){ "source_before_folding" }else{ "retained_source" }
  attr(dag, "recovery_folding") <- recovery_folding
  # Record the temporal window for consumers of the source process graph. If
  # distance calculation fails, save the existing fallback of two slices with
  # reachable = FALSE.
  attr(dag, "unroll_depth") <- tryCatch(
    .query_unroll_depth(.pg_bases, .pg_contemp, .pg_lagged, .pg_treatments, .pg_outcomes),
    error = function( e ){ # Retain the recorded fallback when depth calculation fails
      return(list(reachable = FALSE, fwd = 1L, conf_depth = 0L,
                                treat_slice = 0L, out_slice = 1L, K = 2L))
    }
  )
  attr(dag, "management_options") <- get_management_options(crop = crop, nodes = nodes, aliases = aliases, defaults = management_defaults)
  # Versioned, complete construction receipt. This captures the query
  # actually returned, not a different full graph or a reconstruction from tags.
  request_record <- function(raw, resolution, ids = unname(resolution$resolved)){ # Record supplied names and their resolved node identities
    ids <- unique(ids)
    return(list(raw = raw, resolved = resolution$resolved, resolved_ids = ids,
         present = intersect(ids, identity_rows$member_id),
         absent = setdiff(ids, identity_rows$member_id),
         unresolved = resolution[c("ambiguous", "unmatched")]))
  }
  requests <- list(
    columns = request_record(columns, col_res),
    treatments = request_record(treatments,
      .resolve_columns(treatments, alias_index, disambiguate, hidden_ids), treat_id_set),
    outcomes = request_record(outcomes,
      .resolve_columns(outcomes, alias_index, disambiguate, hidden_ids), out_id_set))
  receipt_ids <- unique(c(identity_rows$member_id, source_node_map$member_id,
                         unlist(lapply(requests, `[[`, "resolved_ids"), use.names = FALSE)))
  policy_columns <- intersect(c("node_id", "concept_id", "node_name_readable", "recovered_node",
                                "hidden_node", "crop_type", "provenance", "tier", "structural_category",
                                "level", "path_dependent"), names(nodes))
  policy <- as.data.frame(as.data.frame(nodes)[nodes$node_id %in% receipt_ids, policy_columns, drop = FALSE],
                          stringsAsFactors = FALSE)
  source_policy <- policy
  policy$recovered_node <- output_policy$recovered_node[match(policy$node_id, output_policy$node_id)]
  effective <- query_input
  effective$time_unroll <- time_unroll
  effective$params["supernode_overrides"] <- list(supernode_overrides)
  effective$params$observe_apsim_outputs <- isTRUE(observe_apsim_outputs)
  attr(dag, "query_record") <- list(
    version = 1L, call = query_call, input = query_input, arguments = effective,
    time_unroll_defaulted = time_unroll_defaulted,
    requests = requests, node_map = unique(identity_rows), policy = policy,
    source_node_map = source_node_map, recovery_folding = recovery_folding,
    source_policy = source_policy,
    policy_scope = "Barley recovery review plus explicit non-APSIM provenance across crops", policy_applied = TRUE,
    aliases = aliases_active[names(aliases_active) %in% receipt_ids],
    scope = list(query_depth = query_depth, time_unroll = time_unroll,
                 phase_window = if( phase_mode ){ phase_window }else{ NULL },
                 apsim_nodes_detail = apsim_nodes_detail,
                 node_recombination = isTRUE(node_recombination),
                 full_query_defaults = is.infinite(query_depth) && !phase_mode &&
                   time_unroll %in% c("all", "none") && !isTRUE(force_acyclic) &&
                   identical(apsim_nodes_detail, "variable") && isTRUE(show_hidden_nodes) &&
                   !isTRUE(node_recombination)))
  if( !is.null(source_order) ){
    attr(dag, "query_record")$source_order <- source_order
  }
  attr(dag, "cf_source_order") <- source_order

  if( isFALSE(quiet) ){
    optional <- sum(attr(dag, "management_options")$status == "optional")
    basis <- if( phase_mode ){
      paste0("anchored to phases (", phase_window, " window)")
    }else{
      switch(time_unroll,
             "all"       = paste0("uncompressed adaptive depth K=", sliced$K),
             "collapse"  = paste0("collapsed static projection (from uncompressed K=", sliced$K, " unroll)"),
             "selective" = paste0("selective adaptive depth K=", sliced$K),
             "pre_post"  = paste0("pre/post relative to treatment (from adaptive depth K=", sliced$K, " unroll)"),
             "contemporaneous")
    }
    detail <- if( length(recombined) > 0L ){ paste0(apsim_nodes_detail, ", ", length(recombined), " recombined") }else{ apsim_nodes_detail }
    message(paste0("Recovered DAG for crop '", crop, "' [", basis, ", apsim detail: ", detail, "]: ",
                   length(recovered_vertices), " nodes, ", nrow(recovered_edges), " edges, ", length(latent_set), " latent. ",
                   length(present), " management levers present, ", optional, " more available. See get_management_options()."))
  }
  return(dag)
}


#' Options for a crop query
#'
#' Creates the specialist options supplied to [dagby()] through its `params`
#' argument. Defaults can be used without supplying any arguments.
#'
#' @param query_depth `Inf` retains the configured graph within the chosen crop and
#'   temporal window. `2L` retains ancestors of treatment and outcome. Each
#'   further integer adds a generation of parents and children. `1L` retains
#'   directed treatment to outcome paths and the endpoints' direct parents.
#'   That restricted scope can omit confounding paths, so adjustment sets need
#'   checking against the fuller graph. Values must be at least one. Without
#'   both treatment and outcome, this option does not restrict the graph.
#' @param water_model The bundled data require `"SoilWater"`. Requesting `"SWIM"`
#'   with these data raises an error.
#' @param apsim_nodes_detail The bundled data require `"variable"`, retaining each
#'   readable quantity key. Coarser grouping settings raise an error.
#' @param supernode_overrides Must be `NULL` for the bundled data.
#' @param force_acyclic Refines cyclic contracted groups when that grouping is
#'   supported. Defaults to `FALSE`. The bundled data do not support those
#'   grouping settings. This is not a general option for repairing a cyclic graph.
#' @param show_hidden_nodes Retain eligible hidden mechanistic quantities when
#'   `TRUE`, the default. `FALSE` retains directed paths through removed hidden
#'   quantities. Recovery exclusions are removed regardless. This controls graph
#'   inclusion separately from observation status.
#' @param node_recombination Defaults to `FALSE`. With `TRUE`, merge quantities
#'   with identical parents and identical children, preserving acyclicity.
#' @param disambiguate Named mappings from ambiguous column names to quantity keys.
#' @param treatment_growth_stage A growth stage applied to all treatments, or named
#'   stages for individual treatments. Supplying a treatment or outcome growth
#'   stage enables developmental phase queries. Every treatment must then have
#'   an anchor. The crop must have a phase partition.
#' @param outcome_growth_stage The corresponding outcome stages. An unspecified
#'   outcome in a developmental phase query uses the terminal measurement phase.
#' @param phase_window `"origin"` retains the sowing origin. `"compact"` restricts
#'   the window to the anchored span. Used only for growth stage queries.
#' @param observe_apsim_outputs Defaults to `TRUE`. Without columns, included
#'   ordinary hidden quantities are treated as observed under the working
#'   availability assumption. `FALSE` marks them latent. Explicit columns take
#'   precedence and endpoints remain observed. Existing output mappings provide
#'   additional evidence. Recovery exclusions are removed regardless of this setting.
#'
#' @return A validated options list. Without growth stage anchors, temporal
#' representation follows `time_unroll`.
#' @export
dagby_query_params <- function( query_depth            = Inf,
                          water_model            = c("SoilWater", "SWIM"),
                          apsim_nodes_detail     = c("variable", "subnode", "node", "supernode"),
                          supernode_overrides    = NULL,
                          force_acyclic          = FALSE,
                          show_hidden_nodes      = TRUE,
                          node_recombination     = FALSE,
                          disambiguate           = NULL,
                          treatment_growth_stage = NULL,
                          outcome_growth_stage   = NULL,
                          phase_window           = c("origin", "compact"),
                          observe_apsim_outputs  = TRUE
                          ){
  # Validate and collect the expert query options
  water_model <- match.arg(water_model)
  apsim_nodes_detail <- match.arg(apsim_nodes_detail)
  phase_window <- match.arg(phase_window)
  if( !is.logical(observe_apsim_outputs) || length(observe_apsim_outputs) != 1L || is.na(observe_apsim_outputs) ){
    stop("observe_apsim_outputs must be TRUE or FALSE.")
  }
  if( !is.numeric(query_depth) || length(query_depth) != 1L || is.na(query_depth) || query_depth < 1 ){
    stop("query_depth must be a single number >= 1 (1 = the adjustment scope of causal paths + parents of treatment and outcome, 2 = the full ancestral subgraph, higher adds successive generations of context, Inf = the full graph). Please supply a valid value and try again.")
  }
  query_depth <- if( is.infinite(query_depth) ){ Inf }else{ as.integer(query_depth) }
  return(list( query_depth            = query_depth,
               water_model            = water_model,
               apsim_nodes_detail     = apsim_nodes_detail,
               supernode_overrides    = supernode_overrides,
               force_acyclic          = isTRUE(force_acyclic),
               show_hidden_nodes      = isTRUE(show_hidden_nodes),
               node_recombination     = isTRUE(node_recombination),
               disambiguate           = disambiguate,
               treatment_growth_stage = treatment_growth_stage,
               outcome_growth_stage   = outcome_growth_stage,
               phase_window           = phase_window,
               observe_apsim_outputs  = observe_apsim_outputs ))
}


## Structural model conversion ##
#' Convert a configured structural node table to a DAG
#'
#' Converts a supplied node table without requiring treatment or outcome endpoints.
#' Every quantity in the chosen crop and water configuration is retained,
#' including management, hidden and isolated quantities. Identifiers are preserved.
#' This is a structural conversion, not a claim that the table is a scientifically
#' complete structural causal model.
#'
#' @param scm A nonempty data frame with unique `node_id` values and `function_of`.
#'   Optional columns include `concept_id`, `composition_of`, `init_of`,
#'   `crop_type`, `tier`, `hidden_node` and `path_dependent`.
#' @param crop Required when crop specific rows or parent references occur.
#' @param water_model SoilWater or SWIM. Select one when both families occur.
#'   A sole family is inferred. The unused family's rows and edges are excluded.
#' @param time_slices Positive finite window length. Required for lagged edges,
#'   initial parents or persistent stocks, and must exceed the largest retained
#'   lag. One slice has no persistence arrows between slices.
#'
#' @details
#' Parent references use quantity keys or `c:concept_id`, separated by semicolons,
#' with optional crop qualifiers and recorded lags. Composition references add
#' arrows from parts to the whole at the same occasion. Unresolved or ambiguous
#' references raise errors. Crop qualifiers are applied before duplicate arrows
#' are removed.
#'
#' Within an explicit window, all selected quantities are retained at every
#' slice. A lag of `k` connects slice `s` to `s + k` when both are within the
#' window. Initial parents enter slice zero. Persistent stocks receive arrows
#' between consecutive slices. Without those temporal features, omitting the
#' window returns one graph without time suffixes. Cycles are rejected.
#' Alternative crop and water configurations are not combined into one model.
#'
#' @return A dagitty object with source quantities, resolved edges, latent
#' status, configuration and separate contemporaneous, lagged and initial edge
#' records. It has no exposure or outcome tags. Optional attributes with `NULL`
#' values are absent.
#' @examples
#' table <- data.frame(node_id = c("weather", "growth", "record"),
#'                     function_of = c("NA", "weather", "NA"))
#' scm_to_dag(table)
#' table$function_of[2] <- "weather@t-2"
#' scm_to_dag(table, time_slices = 3L)
#' @export
scm_to_dag <- function( scm, crop = NULL, water_model = NULL, time_slices = NULL ){
  .datatable.aware <- TRUE
  if( !is.data.frame(scm) || nrow(scm) == 0L ||
       !all(c("node_id", "function_of") %in% names(scm)) ){
    stop("scm must be a nonempty data.frame with node_id and function_of columns.")
  }
  nodes <- as.data.frame(scm, stringsAsFactors = FALSE)
  defaults <- list(concept_id = "", composition_of = "", crop_type = "agnostic",
                   tier = "", hidden_node = "FALSE")
  for( column_name in names(defaults) ){
    if( !column_name %in% names(nodes) ){ nodes[[ column_name ]] <- defaults[[ column_name ]] }
  }
  text_cols <- intersect(c("node_id", "concept_id", "function_of", "composition_of",
                           "init_of", "crop_type", "tier", "hidden_node"), names(nodes))
  for( column_name in text_cols ){ nodes[[ column_name ]] <- as.character(nodes[[ column_name ]]) }
  if( anyNA(nodes$node_id) || any(!nzchar(nodes$node_id)) || anyDuplicated(nodes$node_id) ){
    stop("scm$node_id must contain unique, nonempty ids.")
  }
  if( any(grepl('["\\\\\n\r]', nodes$node_id)) ){
    stop("scm$node_id cannot contain quotes, backslashes or line breaks.")
  }
  for( column_name in c("concept_id", "tier", "hidden_node") ){
    nodes[[ column_name ]][ is.na(nodes[[ column_name ]]) ] <- defaults[[ column_name ]]
  }
  nodes$crop_type[ is.na(nodes$crop_type) | !nzchar(nodes$crop_type) ] <- "agnostic"
  cells <- unlist(nodes[ intersect(c("function_of", "composition_of", "init_of"), names(nodes)) ],
                  use.names = FALSE)
  parsed <- lapply(cells, .parse_parents)
  qualifiers <- unique(unlist(lapply(parsed, function(parsed_parents) { # Collect the crop qualifiers present in this parent cell
    return(parsed_parents$crops[ !is.na(parsed_parents$crops) ])
  })))
  crops <- unique(trimws(unlist(strsplit(nodes$crop_type, ";", fixed = TRUE))))
  crops <- setdiff(crops, "agnostic")
  qualified_crops <- trimws(unlist(strsplit(qualifiers, "+", fixed = TRUE)))
  choices <- sort(unique(c(crops, qualified_crops)))
  if( is.null(crop) && length(choices) > 0L ){
    stop("Select crop explicitly when rows or parent qualifiers name particular crops: ",
         paste(choices, collapse = ", "), ".")
  }
  if( !is.null(crop) && (!is.character(crop) || length(crop) != 1L || is.na(crop) ||
                          !nzchar(crop) || (length(choices) > 0L && !crop %in% choices)) ){
    stop("crop must select one crop represented by the table or its parent qualifiers.")
  }
  families <- c("SoilWater", "SWIM")
  families <- families[ vapply(families, function(x) { # Check whether this water model family occurs in the table
    return(any(startsWith(nodes$node_id, paste0(x, "."))))
  }, logical(1L)) ]
  if( is.null(water_model) && length(families) > 1L ){
    stop("Select water_model = \"SoilWater\" or \"SWIM\". Alternative models cannot be combined.")
  }
  if( !is.null(water_model) ){
    water_model <- match.arg(water_model, c("SoilWater", "SWIM"))
  }else if( length(families) == 1L ){
    water_model <- families
  }
  kept <- .filter_nodes(nodes, if( is.null(crop) ){ "agnostic" }else{ crop },
                        water_model, nodes$node_id)
  if( length(kept) == 0L ){ stop("The selected configuration contains no nodes.") }
  # Resolve using the complete vocabulary, then retain only the selected configuration.
  # This permits legitimate references to alternative water families to be filtered out.
  active <- nodes[ nodes$node_id %in% kept, , drop = FALSE ]
  concept <- sub("^c:", "", active$concept_id)
  ambiguous <- unique(concept[ duplicated(concept) & nzchar(concept) ])
  for( cell in unlist(active[ intersect(c("function_of", "composition_of", "init_of"), names(active)) ], use.names = FALSE) ){
    parsed_parents <- .parse_parents(cell)
    if( nrow(parsed_parents) == 0L ){ next }
    parsed_parents <- parsed_parents[ vapply(parsed_parents$crops, .qualifier_admits, logical(1L), crop = crop), , drop = FALSE ]
    collisions <- intersect(parsed_parents$token, intersect(nodes$node_id, concept))
    conflicting <- vapply(collisions, function(token) { # Detect a token that identifies different nodes and concepts
      return(any(active$node_id[ concept == token ] != token))
    }, logical(1L))
    if( any(conflicting) ){
      stop("A parent identifier has a concept/node_id collision and can resolve to multiple nodes. Make the conflicting identifiers distinct.")
    }
    if( any(parsed_parents$token %in% ambiguous) ){
      stop("A parent concept resolves to multiple nodes. Use explicit node_id references.")
    }
  }
  for( column_name in intersect(c("composition_of", "init_of"), names(active)) ){
    if( any(vapply(active[[ column_name ]], function(cell) { # Detect lag annotations in this parent cell
      return(any(.parse_parents(cell)$lag > 0L))
    }, logical(1L))) ){
      stop("Lag annotations are supported only in function_of. Revise ", column_name, ".")
    }
  }
  # Prefer the selected configuration when alternative families share a concept.
  # Inactive rows remain vocabulary entries so explicit references to an inactive
  # node can be resolved and excluded with that configuration, rather than guessed.
  nodes <- rbind(active, nodes[ !nodes$node_id %in% kept, , drop = FALSE ])
  # Inactive rows are vocabulary entries only: their parent cells need not be resolved.
  inactive <- !nodes$node_id %in% kept
  for( column_name in intersect(c("function_of", "composition_of", "init_of"), names(nodes)) ){
    nodes[[ column_name ]][ inactive ] <- ""
  }
  resolved <- withCallingHandlers(.resolve_edges(nodes, crop, preserve_lags = TRUE),
                                  warning = function(w) { # Reject unresolved structural references
                                    stop(conditionMessage(w), call. = FALSE)
                                  })
  selected <- lapply(resolved, function(e) { # Retain edges within the selected configuration
    return(e[ from %in% kept & to %in% kept ])
  })
  if( !.is_acyclic(kept, selected$contemporaneous) ){
    stop("The selected contemporaneous graph contains a cycle. Revise the table or select a compatible configuration.")
  }
  stocks <- if( "path_dependent" %in% names(active) ){
    active$node_id[ !is.na(active$path_dependent) & toupper(as.character(active$path_dependent)) == "TRUE" ]
  }else{ character(0L) }
  dynamic <- nrow(selected$lagged) > 0L || nrow(selected$initial) > 0L || length(stocks) > 0L
  max_lag <- if( nrow(selected$lagged) > 0L ){ max(selected$lagged$lag) }else{ 0L }
  if( is.null(time_slices) && dynamic ){
    stop("Supply time_slices explicitly for a table with lagged, initial or persistent structure.")
  }
  if( !is.null(time_slices) && (!is.numeric(time_slices) || length(time_slices) != 1L ||
                                 is.na(time_slices) || !is.finite(time_slices) || time_slices < 1 ||
                                 time_slices != floor(time_slices) || time_slices > .Machine$integer.max || time_slices <= max_lag) ){
    stop("time_slices must be a positive integer greater than the largest retained lag (", max_lag, ").")
  }
  if( is.null(time_slices) ){
    vertices <- kept
    edges <- selected$contemporaneous
  }else{
    time_slices <- as.integer(time_slices)
    lags <- selected$lagged
    if( length(stocks) > 0L ){
      if( !"lag" %in% names(lags) ){ lags$lag <- integer(0L) }
      lags <- unique(data.table::rbindlist(list(lags, data.table::data.table(from = stocks, to = stocks, lag = 1L))))
    }
    sliced <- .unroll_k_slices(kept, selected$contemporaneous, lags, time_slices)
    vertices <- sliced$vertices
    edges <- sliced$edges
    if( nrow(selected$initial) > 0L ){
      edges <- unique(rbind(edges, data.table::data.table(from = paste0(selected$initial$from, "@t0"),
                                                  to = paste0(selected$initial$to, "@t0"))))
    }
  }
  if( !.is_acyclic(vertices, edges) ){
    stop("The configured temporal graph contains a cycle, including initial parents. Revise the table.")
  }
  hidden <- active$node_id[ toupper(active$hidden_node) == "TRUE" ]
  latent <- if( is.null(time_slices) ){
    intersect(vertices, hidden)
  }else if( length(hidden) == 0L ){
    character(0L)
  }else{
    intersect(vertices, as.vector(vapply(seq_len(time_slices) - 1L,
                                  function(s) { # Name the hidden nodes at this time slice
                                    return(paste0(hidden, "@t", s))
                                  }, character(length(hidden)))))
  }
  code <- .emit_dagitty(vertices, edges, latent, character(0L), character(0L))
  dag <- if( requireNamespace("dagitty", quietly = TRUE) ){ dagitty::dagitty(code) }else{ code }
  attr(dag, "dag_text") <- code
  attr(dag, "recovered_vertices") <- vertices
  attr(dag, "edges") <- edges
  attr(dag, "latent") <- latent
  attr(dag, "source_nodes") <- active
  attr(dag, "crop") <- crop
  attr(dag, "water_model") <- water_model
  attr(dag, "time_slices") <- time_slices
  attr(dag, "contemporaneous_edges") <- selected$contemporaneous
  attr(dag, "lagged_edges") <- selected$lagged
  attr(dag, "initial_edges") <- selected$initial
  return(dag)
}


#' Use display labels for a query graph
#'
#' Returns a separate graph using `node_name_readable` labels while retaining edges,
#' roles, query records and endpoint identities. The input is unchanged. Supply
#' a `dagby_query()` result with its saved vertices and edges. A bare dagitty
#' object or graph string is insufficient.
#'
#' @param dag A `dagby_query()` result with its saved vertices and edges.
#' @param nodes Table with `node_id` and `node_name_readable` columns.
#'   Defaults to `apsim_dag_nodes`.
#' @details
#' Labels are matched to quantity keys with time suffixes removed, then the
#' original suffix is restored. Empty or conflicting labels fall back to the
#' quantity key so distinct nodes remain distinct. A display label identical
#' to the key leaves the name unchanged.
#' @return A relabelled dagitty object retaining the query's edges, roles,
#'   records and endpoint identities.
#' @export
readable_dag <- function( dag, nodes = apsim_dag_nodes ){
  .datatable.aware <- TRUE
  vertices <- attr(dag, "recovered_vertices")
  edges <- attr(dag, "edges")
  latent <- attr(dag, "latent")
  exposures <- attr(dag, "treatments")
  outcomes <- attr(dag, "outcomes")
  if( is.null(vertices) || is.null(edges) ){
    stop("readable_dag() needs a DAG produced by dagby_query() (it reads the recovered_vertices and edges attributes). Please pass a dagby_query() result and try again.")
  }
  if( !("node_name_readable" %in% names(nodes)) ){
    stop("The node table must carry a 'node_name_readable' column to build display labels. Please pass nodes = apsim_dag_nodes and try again.")
  }

  # Keep the exact suffix so relabelling preserves the occurrence names
  base <- .strip_slice(vertices)
  suffix <- substring(vertices, nchar(base) + 1L)

  # Resolve missing and conflicting labels before restoring temporal suffixes
  human_labels <- stats::setNames(as.character(nodes$node_name_readable), nodes$node_id)
  unique_bases <- unique(base)
  candidate_labels <- unname(human_labels[ unique_bases ])
  missing_labels <- is.na(candidate_labels) | !nzchar(candidate_labels)
  candidate_labels[ missing_labels ] <- unique_bases[ missing_labels ]
  candidate_labels <- stats::setNames(candidate_labels, unique_bases)
  candidate_labels <- .query_unique_labels(candidate_labels, base, suffix)

  relabel <- stats::setNames(paste0(unname(candidate_labels[ base ]), suffix), vertices)
  new_vertices <- unname(relabel[ vertices ])
  new_edges <- data.table::data.table(from = unname(relabel[ edges$from ]), to = unname(relabel[ edges$to ]))
  new_latent <- unname(relabel[ intersect(latent, vertices) ])
  new_exposures <- unname(relabel[ intersect(exposures, vertices) ])
  new_outcomes <- unname(relabel[ intersect(outcomes, vertices) ])

  graph_text <- .emit_dagitty(new_vertices, new_edges, new_latent, new_exposures, new_outcomes)
  result <- if( requireNamespace("dagitty", quietly = TRUE) ){ dagitty::dagitty(graph_text) }else{ graph_text }
  for( key in setdiff(names(attributes(dag)), c("class", "names")) ){
    attr(result, key) <- attr(dag, key)
  }
  attr(result, "dag_text") <- graph_text
  attr(result, "recovered_vertices") <- new_vertices
  attr(result, "edges") <- new_edges
  attr(result, "latent") <- new_latent
  attr(result, "treatments") <- new_exposures
  attr(result, "outcomes") <- new_outcomes
  record <- attr(dag, "query_record")
  if( !is.null(record) ){
    record$node_map$node <- unname(relabel[record$node_map$node])
    record$display_transform <- c(record$display_transform, list(list(
      operation = "readable_dag", mapping = relabel)))
    attr(result, "query_record") <- record
    canonical <- unique(record$node_map[c("node", "quantity")])
    attr(result, "input_labels") <- stats::setNames(canonical$quantity, .strip_slice(canonical$node))
  }
  return(result)
}


#' List available management practices
#'
#' Returns each management quantity's label, aliases and status: `core` for
#' default inclusion or `optional` for explicit selection. `crop` restricts the
#' result to applicable quantities. These defaults select quantities for
#' inclusion in a query. They do not specify management values.
#'
#' @param crop Restrict the result to applicable quantities.
#' @param nodes Defaults to `apsim_dag_nodes`.
#' @param aliases Defaults to `node_aliases`.
#' @param defaults Defaults to `management_defaults`.
#' @return A data.table with `node_id`, `label`, `aliases`, `status` and
#'   `crop_type`. Rows marked `optional` are the available explicit choices.
#'   If no management quantity qualifies, returns these columns with no rows.
#' @examples
#' \dontrun{
#' get_management_options(crop = "wheat")
#' }
#' @export
get_management_options <- function( crop     = NULL,
                                    nodes    = apsim_dag_nodes,
                                    aliases  = node_aliases,
                                    defaults = management_defaults
                                    ){

  .datatable.aware <- TRUE
  management <- nodes[ nodes$tier == "management", ]
  if( nrow(management) == 0L ){
    return(data.table::data.table(node_id = character(), label = character(),
                                  aliases = character(), status = character(),
                                  crop_type = character()))
  }
  rows <- lapply(seq_len(nrow(management)), function( i ){ # Describe one management node and its searchable aliases
    id <- management$node_id[ i ]
    management_aliases <- aliases[[ id ]]
    if( is.null(management_aliases) ){
      management_aliases <- unique(c(id, management$node_name_readable[ i ]))
    }
    return(data.table::data.table(node_id   = id,
                           label     = management$node_name_readable[ i ],
                           aliases   = paste(management_aliases, collapse = "; "),
                           status    = if( id %in% defaults$core ){ "core" }else{ "optional" },
                           crop_type = management$crop_type[ i ]))
  })
  options <- data.table::rbindlist(rows)
  if( length(crop) > 0L ){
    keep <- options$crop_type == "agnostic" |
            sapply(strsplit(options$crop_type, ";", fixed = TRUE), function( z ){ # Check whether the listed crops include the requested crop
              return(crop %in% trimws(z))
            })
    options <- options[ keep ]
  }
  return(options[])
}


## Query engine internals ##
#' Validate carried source ordering against the current quantity identities
#'
#' @param nodes Prepared node data frame with `node_id`, `concept_id` and an
#'   optional `cf_source_order` attribute.
#' @return The original ordering list if its version, source ID, quantity
#'   identities and positive, unique integer ranks pass validation against
#'   the current table. Otherwise, `NULL`.
#' @noRd
.query_source_order <- function(nodes){
  source_order <- attr(nodes, "cf_source_order", exact = TRUE)
  if( !is.list(source_order) || !identical(source_order$version, 1L) ||
      !identical(source_order$source, "apsim-slice-2026-11-05") ){
    return(NULL)
  }
  identities <- source_order$nodes
  if( !is.data.frame(identities) || !all(c("node_id", "concept_id") %in% names(identities)) ||
      nrow(identities) != nrow(nodes) || anyDuplicated(identities$node_id) ||
      !setequal(identities$node_id, nodes$node_id) ){
    return(NULL)
  }
  matching <- match(nodes$node_id, identities$node_id)
  if( !identical(nodes$concept_id, identities$concept_id[matching]) ){
    return(NULL)
  }
  rank_fields <- c("base_order", "at_order", "post_order", "pre_order", "time_order") # Distinct occurrence intervals
  ranks <- source_order$ranks
  if( !is.data.frame(ranks) || !all(c("node_id", rank_fields) %in% names(ranks)) ||
      nrow(ranks) != nrow(nodes) || anyDuplicated(ranks$node_id) ||
      !setequal(ranks$node_id, nodes$node_id) ){
    return(NULL)
  }
  values <- unlist(ranks[rank_fields], use.names = FALSE)
  if( !is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 1 | values != floor(values)) || anyDuplicated(values) ){
    return(NULL)
  }
  return(source_order)
}


#' Parse one structural parent cell
#'
#' A trailing `@t-<k>` marks a lag of k steps. A leading `c:` marks a concept
#' reference. A crop qualifier in brackets appears before any lag annotation.
#'
#' @param cell One semicolon-separated string, missing value or empty vector.
#' @return A data frame with character `token` and integer `lag` columns. A
#'   nonempty result also has `crops`, with `NA` for an unqualified reference.
#'   Concept prefixes and trailing lag annotations are removed from tokens.
#' @noRd
.parse_parents <- function( cell ){
  if( length(cell) == 0L || is.na(cell) || !nzchar(cell) ){
    return(data.frame(token = character(0L), lag = integer(0L), stringsAsFactors = FALSE))
  }
  tokens <- trimws(strsplit(cell, ";", fixed = TRUE)[[ 1L ]])
  tokens <- tokens[ nzchar(tokens) & !is.na(tokens) & tokens != "NA" ]  # "NA" = marker for a confirmed root
  if( length(tokens) == 0L ){
    return(data.frame(token = character(0L), lag = integer(0L), stringsAsFactors = FALSE))
  }
  lags <- rep(0L, length(tokens))
  has_lag <- grepl("@t-[0-9]+$", tokens)
  if( any(has_lag) ){
    lags[ has_lag ] <- as.integer(sub("^.*@t-([0-9]+)$", "\\1", tokens[ has_lag ]))
    tokens[ has_lag ] <- sub("@t-[0-9]+$", "", tokens[ has_lag ])
  }
  # The crop qualifier is now at the end of the token because the lag has been
  # removed. NA means every crop admitted by the child. Otherwise, crop names
  # are separated by plus signs.
  crops <- rep(NA_character_, length(tokens))
  has_qualifier <- grepl("\\[[^]]*\\]$", tokens, perl = TRUE)
  if( any(has_qualifier) ){
    qualifiers <- regmatches(tokens[ has_qualifier ], regexpr("\\[[^]]*\\]$", tokens[ has_qualifier ], perl = TRUE))
    crops[ has_qualifier ] <- substr(qualifiers, 2L, nchar(qualifiers) - 1L)
    tokens[ has_qualifier ] <- sub("\\[[^]]*\\]$", "", tokens[ has_qualifier ], perl = TRUE)
  }
  tokens <- trimws(sub("^c:", "", tokens))
  return(data.frame(token = tokens, lag = lags, crops = crops, stringsAsFactors = FALSE))
}


#' Test a parent reference's crop qualifier
#'
#' @param crops One plus-separated crop qualifier, or `NA` for no restriction.
#' @param crop One selected crop, or `NULL` to admit every qualifier.
#' @return `TRUE` if no crop is selected, the reference has no crop qualifier,
#'   or the selected crop is included in the qualifier. Otherwise, `FALSE`.
#' @noRd
.qualifier_admits <- function( crops, crop ){
  if( length(crop) == 0L || is.na(crops) ){
    return(TRUE)
  }
  return( crop %in% trimws(strsplit(crops, "+", fixed = TRUE)[[ 1L ]]) )
}


# Resolve function_of and composition_of references using node IDs or concept
# keys. Only positive lags from function_of enter the lagged edge table.
# Preserve their magnitude when preserve_lags is TRUE. Filter crop qualifiers
# before deduplication so a path cannot cross an edge absent from the selected
# crop. Unresolved references are warned about and excluded.
#' Resolve structural references to directed quantity edges
#'
#' @param nodes Prepared node table containing `node_id`, `concept_id`,
#'   `function_of`, `composition_of` and optionally `init_of`.
#' @param crop Selected crop, or `NULL` to admit every edge qualifier.
#' @param preserve_lags Retain recorded lag magnitudes when `TRUE`.
#' @return A list of `contemporaneous`, `lagged` and `initial` data.tables with
#'   `from` and `to` columns. A resolved lagged table also includes `lag` when
#'   requested. Unresolved references are warned about and excluded.
#' @noRd
.resolve_edges <- function( nodes, crop = NULL, preserve_lags = FALSE ){
  .datatable.aware <- TRUE
  by_id <- nodes$node_id
  concept_stripped <- sub("^c:", "", nodes$concept_id)
  first_seen <- !duplicated(concept_stripped) & nzchar(concept_stripped)
  by_concept <- stats::setNames(nodes$node_id[ first_seen ], concept_stripped[ first_seen ])
  has_init <- "init_of" %in% names(nodes)

  all_tokens <- unique(unlist(lapply(seq_len(nrow(nodes)), function( i ){ # Collect the admitted parent tokens for one node
    cells <- c(nodes$function_of[ i ], nodes$composition_of[ i ],
               if( has_init ){ nodes$init_of[ i ] }else{ NULL })
    return(unlist(lapply(cells, function( cell ){ # Keep parent tokens admitted by the crop qualifier
      parsed_parents <- .parse_parents(cell)
      return(if( nrow(parsed_parents) > 0L ){
        parsed_parents$token[ vapply(parsed_parents$crops, .qualifier_admits, logical(1L), crop = crop) ]
      }else{ character(0L) })
    })))
  })))
  unresolved <- all_tokens[ !(all_tokens %in% by_id) &
                            !(all_tokens %in% names(by_concept)) ]
  if( length(unresolved) > 0L ){
    warning(paste0("These parent references resolve to no node and were excluded: ",
                   paste(sort(unresolved), collapse = ", "),
                   ". Check function_of / composition_of in the node table and try again."))
  }

  resolve_token <- function( tok ){ # Resolve one parent token through the node and concept indexes
    if( tok %in% by_id ){
      return(tok)
    }
    hit <- by_concept[ tok ]
    if( is.na(hit) ){
      return(NA_character_)
    }
    return(unname(hit))
  }

  # Keep the initial network separate from the transition edges so its parent
  # relationships can be applied only at the origin slice.
  initial_edges <- data.table::data.table(from = character(0L), to = character(0L))
  if( has_init ){
    init_records <- unlist(
      lapply(seq_len(nrow(nodes)), function( i ){ # Resolve the initial parents of one node
        child <- nodes$node_id[ i ]
        parsed <- .parse_parents(nodes$init_of[ i ])
        if( nrow(parsed) > 0L ){
          parsed <- parsed[ vapply(parsed$crops, .qualifier_admits, logical(1L), crop = crop), , drop = FALSE ]
        }
        if( nrow(parsed) == 0L ){
          return(NULL)
        }
        return(lapply(seq_len(nrow(parsed)), function( r ){ # Pair an initial parent with its child
          return(c(from = resolve_token(parsed$token[ r ]), to = child))
        }))
      }),
      recursive = FALSE
    )
    init_records <- Filter(Negate(anyNA), init_records)
    if( length(init_records) > 0L ){
      initial_edges <- unique(data.table::as.data.table(do.call(rbind, init_records)))
    }
  }

  records <- suppressWarnings(
    unlist(

      lapply(seq_len(nrow(nodes)), function( i ){ # Resolve the structural parents of one node
        child <- nodes$node_id[ i ]
        cells <- c(function_of = nodes$function_of[ i ], composition_of = nodes$composition_of[ i ])
        # cells is always length 2 (function_of, composition_of), so 1:length is safe
        per_node <- lapply(1:length(cells), function( k ){ # Resolve parents from one structural column
          parsed <- .parse_parents(cells[[ k ]])
          if( nrow(parsed) > 0L ){
            parsed <- parsed[ vapply(parsed$crops, .qualifier_admits, logical(1L), crop = crop), , drop = FALSE ]
          }
          if( nrow(parsed) == 0L ){
            return(NULL)
          }
          return(lapply(1:nrow(parsed), function( r ){ # Record a resolved parent edge with its lag and source column
            return(c(from = resolve_token(parsed$token[ r ]),
              to   = child,
              lag  = as.character(parsed$lag[ r ]),
              kind = names(cells)[ k ]))
          }))
        })
        return(unlist(per_node, recursive = FALSE))
      }),
      recursive = FALSE
    )
  )

  records <- Filter(Negate(anyNA), records)
  if( length(records) == 0L ){
    empty <- data.table::data.table(from = character(0L), to = character(0L))
    return(list(contemporaneous = empty, lagged = data.table::copy(empty), initial = initial_edges))
  }
  edges <- data.table::as.data.table(do.call(rbind, records))
  edges[, lag := as.integer(lag) ]
  edges[, is_lagged := lag > 0L & kind == "function_of" ]
  contemporaneous <- unique(edges[ is_lagged == FALSE, list(from, to) ])
  lagged_edges <- if( isTRUE(preserve_lags) ){
    unique(edges[ is_lagged == TRUE, list(from, to, lag) ])
  }else{
    unique(edges[ is_lagged == TRUE, list(from, to) ])
  }
  return(list(contemporaneous = contemporaneous, lagged = lagged_edges, initial = initial_edges))
}


# Kahn topological sort over the contemporaneous (lag 0) edges: TRUE iff every
# node is ordered, i.e. the lag-0 structural graph is acyclic.
#' Check a directed graph by topological traversal
#'
#' @param node_ids Unique names of all vertices, including isolated vertices.
#' @param contemporaneous Directed edge table with `from` and `to` in `node_ids`.
#' @return `TRUE` if the directed graph is acyclic, otherwise `FALSE`.
#' @noRd
.is_acyclic <- function( node_ids, contemporaneous ){
  children_by_parent <- split(contemporaneous$to, contemporaneous$from)
  indegree <- stats::setNames(rep(0L, length(node_ids)), node_ids)
  parent_counts <- table(contemporaneous$to)
  indegree[ names(parent_counts) ] <- as.integer(parent_counts)
  queue <- node_ids[ indegree[ node_ids ] == 0L ]
  visited <- 0L
  while( length(queue) > 0L ){
    parent <- queue[ 1L ]
    queue <- queue[ -1L ]
    visited <- visited + 1L
    children <- children_by_parent[[ parent ]]
    if( !is.null(children) ){
      for( child in children ){
        indegree[ child ] <- indegree[ child ] - 1L
        if( indegree[ child ] == 0L ){
          queue <- c(queue, child)
        }
      }
    }
    children_by_parent[[ parent ]] <- NULL
  }
  return(visited == length(node_ids))
}


#' Retain only unambiguous column-to-quantity matches
#'
#' @param columns Supplied column names, or `NULL`.
#' @param alias_index Alias lookup returned by `.build_alias_index()`.
#' @param hidden_ids Hidden quantity IDs used by the resolver's tie rule.
#' @return Unique resolved node IDs, omitting ambiguous and unmatched columns.
#'   Returns `character(0)` for `NULL` or empty columns and `NULL` if supplied
#'   columns resolve to no IDs.
#' @noRd
.match_columns <- function( columns, alias_index, hidden_ids = character(0L) ){
  if( length(columns) == 0L ){
    return(character(0L))
  }
  resolved_columns <- lapply(columns, function( col ){ # Resolve one supplied column through the alias index
    return(.resolve_one(col, alias_index, hidden_ids))
  })
  unambiguous_columns <- Filter(function( v ){ # Retain only an unambiguous column match
    return(length(v) == 1L)
  }, resolved_columns)
  return(unique(unlist(unambiguous_columns, use.names = FALSE)))
}


#' Select management quantities admitted to the query
#'
#' @param nodes Prepared node table with `node_id` and `tier`.
#' @param management_practices Default or explicitly selected management names.
#' @param forced Endpoint names or IDs that require inclusion.
#' @param col_matched Quantity IDs resolved from supplied data columns.
#' @param alias_index Alias lookup returned by `.build_alias_index()`.
#' @return Unique IDs in the management tier admitted by any supplied source.
#' @noRd
.resolve_management_presence <- function( nodes, management_practices, forced, col_matched, alias_index ){
  management_ids <- nodes$node_id[ nodes$tier == "management" ]
  candidates <- unique(c(management_practices,
                         .match_columns(management_practices, alias_index),
                         forced,
                         .match_columns(forced, alias_index),
                         col_matched))
  return(intersect(candidates, management_ids))
}


# Node deletion: drop nodes incompatible with the crop, the unused soil water family,
# and management nodes that are not present. Returns the surviving node_ids.
#' Filter the node configuration without changing the source table
#'
#' @param nodes Prepared table with `node_id`, `crop_type` and `tier`.
#' @param crop One crop name matched to semicolon-separated crop memberships.
#' @param water_model Selected `SoilWater` or `SWIM` family.
#' @param management_present Management IDs admitted to the query.
#' @return Surviving node IDs in their source order.
#' @noRd
.filter_nodes <- function( nodes, crop, water_model, management_present ){
  keep <- rep(TRUE, nrow(nodes))
  crop_ok <- nodes$crop_type == "agnostic" |
             sapply(strsplit(nodes$crop_type, ";", fixed = TRUE), function( z ){ # Check whether the listed crops include the requested crop
               return(crop %in% trimws(z))
             })
  keep <- keep & crop_ok
  if( identical(water_model, "SoilWater") ){
    keep <- keep & !startsWith(nodes$node_id, "SWIM.")
  }else if( identical(water_model, "SWIM") ){
    keep <- keep & !startsWith(nodes$node_id, "SoilWater.")
  }
  is_management <- nodes$tier == "management"
  keep <- keep & ( !is_management | nodes$node_id %in% management_present )
  return(nodes$node_id[ keep ])
}


#' Expand quantities over a finite time window
#'
#' Positive lags always move forward in time. The result is acyclic if the
#' contemporaneous graph is acyclic.
#'
#' @param kept_ids Unique quantity IDs retained at every slice.
#' @param contemporaneous Directed `from`/`to` data.table copied within slices.
#' @param lagged Directed `from`/`to` data.table with optional positive `lag`.
#'   Without that column each edge advances one slice.
#' @param K Positive integer number of slices. Indices run from zero to `K - 1`.
#' @return A list with `vertices` named using `@t` suffixes and a directed
#'   `edges` data.table. Edges outside the finite window are omitted.
#' @noRd
.unroll_k_slices <- function( kept_ids, contemporaneous, lagged, K ){
  .datatable.aware <- TRUE
  K <- as.integer(K)
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
  slices <- seq_len(K) - 1L                                   # 0 .. K-1
  vertices <- as.vector(vapply(slices, function( slice ){ # Name the retained nodes at this time slice
    return(paste0(kept_ids, "@t", slice))
  },
                               character(length(kept_ids))))
  within <- if( nrow(contemporaneous_edges) > 0L ){
    data.table::rbindlist(lapply(slices, function( slice ){ # Copy contemporaneous edges into this time slice
      return(data.table::data.table(from = paste0(contemporaneous_edges$from, "@t", slice), to = paste0(contemporaneous_edges$to, "@t", slice)))
    }))
  }else{
    # Construct the empty table explicitly because paste0() would turn empty
    # endpoint vectors into a spurious edge between suffix-only names.
    data.table::data.table(from = character(0L), to = character(0L))
  }
  across <- if( K >= 2L && nrow(lagged_edges) > 0L ){
    edge_lags <- if( "lag" %in% names(lagged_edges) ){ lagged_edges$lag }else{ rep(1L, nrow(lagged_edges)) }
    data.table::rbindlist(lapply(slices[ -length(slices) ], function( slice ){ # Connect lagged edges that end within the finite window
      take <- which(slice + edge_lags < K)
      if( length(take) == 0L ){ return(NULL) }
      return(data.table::data.table(from = paste0(lagged_edges$from[ take ], "@t", slice),
                            to = paste0(lagged_edges$to[ take ], "@t", slice + edge_lags[ take ])))
    }))
  }else{
    data.table::data.table(from = character(0L), to = character(0L))
  }
  edges <- unique(data.table::rbindlist(list(within, across)))
  return(list(vertices = vertices, edges = edges))
}


# Resolve PHASE_OF active_from / active_to tokens to integer phase index windows for a single crop.
# Returns a named list node_id -> c(from_idx, to_idx), suitable for passing as phase_active to
# .unroll_phase_slices. "terminal" is resolved to max(phase_index) for the crop. Nodes whose
# crop_type is not agnostic and does not include the crop are silently skipped (they will already be
# removed by .filter_nodes). Tokens are resolved via the canonical_stage vocabulary exactly as
# .stages_to_phase_indices does (canonical stage, APSIM phase name, or bare integer).
#' Resolve crop-specific active phase intervals
#'
#' @param phase_of Table with `node_id`, `crop_type`, `active_from`, `active_to`.
#' @param crop Selected crop with an available phase partition.
#' @param crop_phases Table of crop, phase index, phase name and canonical stage.
#' @return A list named by node ID, each value containing inclusive first and
#'   last phase indices. Inapplicable or unresolved intervals are omitted.
#'   The phase unroller treats an absent entry as always active.
#' @noRd
.resolve_phase_active <- function( phase_of, crop, crop_phases ) {
  crop_partition <- crop_phases[ crop_phases$crop == crop, ]
  last_phase <- max(crop_partition$phase_index)
  has_canonical_stage <- nzchar(crop_partition$canonical_stage)
  canonical_to_index <- stats::setNames(crop_partition$phase_index[ has_canonical_stage ], crop_partition$canonical_stage[ has_canonical_stage ])
  name_to_index <- stats::setNames(crop_partition$phase_index, crop_partition$phase_name)
  resolve_tok <- function( tok ) { # Resolve a phase token to its integer index
    if( identical(tok, "terminal") ){ return( last_phase ) }
    if( grepl("^[0-9]+$", tok) ){ return( as.integer(tok) ) }
    if( tok %in% names(canonical_to_index) ){ return( unname(canonical_to_index[ tok ]) ) }
    if( tok %in% names(name_to_index) ){ return( unname(name_to_index[ tok ]) ) }
    return( NA_integer_ )
  }
  active_intervals <- list()
  for( i in seq_len(nrow(phase_of)) ){
    node_id <- phase_of$node_id[ i ]
    crop_type <- phase_of$crop_type[ i ]
    if( nzchar(crop_type) && crop_type != "agnostic" ){
      crops <- trimws(strsplit(crop_type, ";", fixed = TRUE)[[ 1L ]])
      if( !(crop %in% crops) ){ next }
    }
    first_phase <- resolve_tok(phase_of$active_from[ i ])
    last_active_phase <- resolve_tok(phase_of$active_to[ i ])
    if( is.na(first_phase) || is.na(last_active_phase) ){ next }    # unresolvable token for this crop -> always active
    active_intervals[[ node_id ]] <- c(first_phase, last_active_phase)
  }
  return( active_intervals )
}


#' Expand the configured graph over developmental phases
#'
#' Phase zero represents sowing. The initial network supplies stock parents at
#' that origin, alongside contemporaneous edges. Later phases repeat the
#' contemporaneous structure. The result is acyclic when the contemporaneous
#' graph and its union with the initial network are both acyclic. This helper
#' does not check either condition.
#'
#' @param kept_ids Unique quantity IDs admitted to the phase frame.
#' @param contemporaneous,lagged,initial Directed `from`/`to` data.tables. Initial
#'   edges occur only at phase zero. Every lagged edge advances one phase.
#' @param n_phases Positive integer number of phases, including the sowing origin.
#' @param phase_active Optional list of inclusive active intervals by node ID.
#'   Nodes absent from this list remain active in every phase.
#' @return A list with `vertices` and `edges`, whose `edge_type` distinguishes
#'   contemporaneous, lagged and initial arrows. With masking, the list also
#'   carries an `active_per_slice` count attribute.
#' @noRd
.unroll_phase_slices <- function( kept_ids, contemporaneous, lagged, initial, n_phases,
                                  phase_active = NULL ){
  .datatable.aware <- TRUE
  n_phases <- as.integer(n_phases)
  if( length(n_phases) != 1L || is.na(n_phases) || n_phases < 1L ){
    stop("n_phases must be a single positive whole number of developmental phase slices. Please supply a phase count of at least 1 and try again.")
  }
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
  initial_edges <- initial[ from %in% kept_ids & to %in% kept_ids ]
  slices <- seq_len(n_phases) - 1L                            # 0 .. n_phases-1 (0 = sowing origin)

  # Active node sets for each slice: if phase_active is NULL (no masking), every kept_id is active in every
  # slice. Otherwise, compute the active subset per slice.
  if( is.null(phase_active) ){
    active_per_slice <- lapply(slices, function( slice ) { # Keep every node when no phase mask is supplied
      return(kept_ids)
    })
  }else{
    active_per_slice <- lapply(slices, function( slice ) { # Select the nodes active at this developmental phase
      mask <- vapply(kept_ids, function( node_id ) { # Check this node's active phase interval
        active_interval <- phase_active[[ node_id ]]
        if( is.null(active_interval) ){ return( TRUE ) }           # not in phase_of -> always active
        return( active_interval[1L] <= slice && slice <= active_interval[2L] )
      }, logical(1L))
      return(kept_ids[ mask ])
    })
  }

  vertices <- unlist(lapply(slices, function( slice ) { # Name the active nodes at this developmental phase
    return(paste0(active_per_slice[[ slice + 1L ]], "@t", slice))
  }),
                     use.names = FALSE)

  # within slice contemporaneous edges, per slice filtered to active nodes
  within <- if( nrow(contemporaneous_edges) > 0L ){
    data.table::rbindlist(lapply(slices, function( slice ) { # Copy contemporaneous edges between active nodes at this phase
      active_ids <- active_per_slice[[ slice + 1L ]]
      active_edges <- contemporaneous_edges[ from %in% active_ids & to %in% active_ids ]
      if( nrow(active_edges) == 0L ){ return( data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L)) ) }
      return(data.table::data.table(from = paste0(active_edges$from, "@t", slice), to = paste0(active_edges$to, "@t", slice),
                             edge_type = "contemporaneous"))
    }))
  }else{
    data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
  }

  # Carry every lagged relationship into the next phase, including stock
  # accumulation, management feedback and faster physiological feedback.
  # This retains forward causal paths without putting feedback within a slice.
  # The one-phase lag is a discretisation choice, not an estimate of duration.
  carry <- lagged_edges
  across <- if( n_phases >= 2L && nrow(carry) > 0L ){
    data.table::rbindlist(lapply(slices[ -length(slices) ], function( slice ) { # Carry lagged edges between adjacent active phases
      active_sources <- active_per_slice[[ slice + 1L ]]
      active_targets <- active_per_slice[[ slice + 2L ]]
      active_edges <- carry[ from %in% active_sources & to %in% active_targets ]
      if( nrow(active_edges) == 0L ){ return( data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L)) ) }
      return(data.table::data.table(from = paste0(active_edges$from, "@t", slice), to = paste0(active_edges$to, "@t", slice + 1L),
                             edge_type = "lagged"))
    }))
  }else{
    data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
  }

  # initial network B0: materialised at the origin slice only (parent@t0 -> stock@t0), filtered to
  # nodes active in slice 0
  origin_nodes <- active_per_slice[[ 1L ]]
  origin_network <- if( nrow(initial_edges) > 0L ){
    origin_edges <- initial_edges[ from %in% origin_nodes & to %in% origin_nodes ]
    if( nrow(origin_edges) == 0L ){
      data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
    }else{
      data.table::data.table(from = paste0(origin_edges$from, "@t0"), to = paste0(origin_edges$to, "@t0"),
                             edge_type = "initial")
    }
  }else{
    data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
  }

  edges <- unique(data.table::rbindlist(list(within, across, origin_network)))
  phase_graph <- list(vertices = vertices, edges = edges)
  if( !is.null(phase_active) ){
    attr(phase_graph, "active_per_slice") <- vapply(active_per_slice, length, integer(1L))
  }
  return( phase_graph )
}


# A single contemporaneous slice has no @t suffix. Lagged edges become ordinary
# edges. Self-loops are omitted and duplicate arrows are combined.
# Feedback can therefore produce cycles: fertiliser and irrigation decisions
# read the soil state they change, and soil nitrogen and microbial pools also
# contain feedback. The caller is warned and can unroll time instead.
#' Combine temporal arrows into a single occurrence per quantity
#'
#' @param kept_ids Quantity IDs retained as unsuffixed vertices.
#' @param contemporaneous Directed `from`/`to` data.table.
#' @param lagged Optional directed `from`/`to` data.table.
#' @return A list with `vertices`, `edges` without duplicates or self-loops,
#'   an `acyclic` flag and `n_collapsed`, the number of distinct retained lagged
#'   endpoint pairs before self-loops and duplicate arrows are removed.
#' @noRd
.contemporaneous_slice <- function( kept_ids, contemporaneous, lagged = NULL ){
  .datatable.aware <- TRUE
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  edges <- data.table::data.table(from = contemporaneous_edges$from, to = contemporaneous_edges$to)

  n_collapsed <- 0L
  if( !is.null(lagged) && nrow(lagged) > 0L ){
    lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
    if( nrow(lagged_edges) > 0L ){
      n_collapsed <- nrow(unique(data.table::data.table(from = lagged_edges$from, to = lagged_edges$to)))
      edges <- rbind(edges, data.table::data.table(from = lagged_edges$from, to = lagged_edges$to))
    }
  }
  edges <- unique(edges)
  # Self-loops can arise when temporal occurrences collapse to one node
  edges <- edges[ from != to ]

  return(list(vertices    = kept_ids,
              edges       = edges,
              n_collapsed = n_collapsed,
              acyclic     = .is_acyclic(kept_ids, edges)))
}


#' Choose a query window from endpoint lag distances
#'
#' Use the greatest shortest lag distance among reachable endpoint pairs so
#' the window includes a path for each pair. Add one slice of history when a
#' treatment has a lagged parent, preserving the state before the decision.
#'
#' @param kept_ids Quantity IDs in the configured graph.
#' @param contemporaneous,lagged Directed `from`/`to` data.tables. A lagged edge
#'   costs one slice, irrespective of any separate lag magnitude column.
#' @param treat_ids,outcome_ids Resolved treatment and outcome quantity IDs.
#' @return A list with `reachable`, forward depth `fwd`, prior history
#'   `conf_depth`, endpoint anchors `treat_slice` and `out_slice`, and window
#'   length `K`. `reachable` is `TRUE` only when every retained endpoint pair
#'   is connected. If none is connected, forward depth is one. If either
#'   endpoint set is absent, returns a two-slice window with `reachable = FALSE`
#'   and an unreachable-distance sentinel in `fwd`.
#' @noRd
.query_unroll_depth <- function( kept_ids, contemporaneous, lagged, treat_ids, outcome_ids ){
  .datatable.aware <- TRUE
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
  unreachable_distance <- .Machine$integer.max %/% 2L
  treatment_seeds <- intersect(treat_ids, kept_ids)
  selected_outcomes <- intersect(outcome_ids, kept_ids)
  if( length(treatment_seeds) == 0L || length(selected_outcomes) == 0L ){
    return(list(reachable = FALSE, fwd = unreachable_distance, conf_depth = 0L, treat_slice = 0L, out_slice = 1L, K = 2L))
  }

  .lag_dist <- function( src ){ # Find the shortest lag distance from one treatment
    distances <- stats::setNames(rep(unreachable_distance, length(kept_ids)), kept_ids)
    distances[ src ] <- 0L
    repeat {
      changed <- FALSE
      if( nrow(contemporaneous_edges) > 0L ){
        better <- which(distances[ contemporaneous_edges$from ] < distances[ contemporaneous_edges$to ])
        for( i in better ){ distances[[ contemporaneous_edges$to[i] ]] <- distances[[ contemporaneous_edges$from[i] ]] }
        changed <- changed || length(better) > 0L
      }
      if( nrow(lagged_edges) > 0L ){
        better <- which(distances[ lagged_edges$from ] + 1L < distances[ lagged_edges$to ])
        for( i in better ){ distances[[ lagged_edges$to[i] ]] <- distances[[ lagged_edges$from[i] ]] + 1L }
        changed <- changed || length(better) > 0L
      }
      if( !changed ){
        break
      }
    }
    return(distances)
  }

  # relax from each treatment separately, then take the maximum over every (treatment, outcome)
  # pair. A pair with no path in the base graph stays disconnected for any K,
  # so it is excluded from the max and reported via reachable = FALSE.
  pair_distances <- integer(0L)
  has_unreachable_pair <- FALSE
  for( s in treatment_seeds ){
    distances <- .lag_dist(s)
    for( o in selected_outcomes ){
      pair_distance <- distances[[ o ]]
      if( pair_distance >= unreachable_distance ){
        has_unreachable_pair <- TRUE
      }else{
        pair_distances <- c(pair_distances, as.integer(pair_distance))
      }
    }
  }
  if( length(pair_distances) == 0L ){
    # Keep a finite window when no pair is connected. The FALSE flag records
    # that this window does not establish a treatment-to-outcome path.
    forward_depth <- 1L
    reachable <- FALSE
  }else{
    forward_depth <- max(pair_distances)
    reachable <- !has_unreachable_pair
  }
  history_depth <- as.integer(any(lagged_edges$to %in% treat_ids))
  treatment_slice <- history_depth
  outcome_slice <- treatment_slice + forward_depth
  return(list(reachable = reachable, fwd = as.integer(forward_depth), conf_depth = history_depth,
              treat_slice = treatment_slice, out_slice = outcome_slice, K = outcome_slice + 1L))
}


#' Project an endpoint-anchored temporal window
#'
#' Determine treatment descendants before removing hidden nodes so paths through
#' hidden occasions still contribute. By default, keep occasions separate for a
#' quantity with both descendant and nondescendant occasions. Each treatment and
#' outcome is anchored once (VanderWeele & Robins 2007).
#'
#' If projection creates a cycle, split further quantities when possible.
#' Splitting stops when no further split is available, so the result is not
#' guaranteed to be acyclic.
#'
#' @param kept_ids Quantity IDs retained in the configured graph.
#' @param contemporaneous,lagged Directed `from`/`to` data.tables.
#' @param hidden_ids Hidden quantity IDs eligible for directed removal.
#' @param treat_ids,outcome_ids Resolved endpoint IDs, each anchored once.
#' @param persistent_ids Quantity IDs whose retained occasions need carryover.
#' @param split_all Keep all eligible quantities split across the finite window.
#' @param label_mode Use separate `slice` names or initial `pre_post` grouping.
#' @param show_hidden_nodes Keep hidden quantities instead of removing their
#'   interiors while retaining directed paths.
#' @return A list containing projected `vertices` and `edges`, the quantity IDs
#'   in `split` and `promoted`, endpoint slice anchors, `reachable`, `K` and
#'   `persistence` arrows.
#' @noRd
.selective_slice <- function( kept_ids, contemporaneous, lagged, hidden_ids, treat_ids, outcome_ids,
                              persistent_ids = character(0L), split_all = FALSE,
                              label_mode = c("slice", "pre_post"),
                              show_hidden_nodes = FALSE ){
  .datatable.aware <- TRUE
  label_mode <- match.arg(label_mode)
  never_split <- union(treat_ids, outcome_ids)

  depth <- .query_unroll_depth(kept_ids, contemporaneous, lagged, treat_ids, outcome_ids)
  K <- depth$K
  treat_slice <- depth$treat_slice
  out_slice <- depth$out_slice

  # full unroll of K slices: the working graph, never emitted
  unrolled <- .unroll_k_slices(kept_ids, contemporaneous, lagged, K)
  hidden_vertex <- stats::setNames(.strip_slice(unrolled$vertices) %in% hidden_ids, unrolled$vertices)

  # forward reachability from each treatment's anchor (treat_slice) vertex over the unrolled edges
  # (through hidden interiors): each vertex's treatment descendant status.
  children_by_parent <- split(unrolled$edges$to, unrolled$edges$from)
  reach_forward <- function( sources ){ # Collect descendants reached from the selected treatment occasions
    seen <- character(0L)
    stack <- unlist(children_by_parent[ sources ], use.names = FALSE)
    while( length(stack) > 0L ){
      vertex <- stack[ length(stack) ]
      stack <- stack[ -length(stack) ]
      if( vertex %in% seen ){
        next
      }
      seen <- c(seen, vertex)
      children <- children_by_parent[[ vertex ]]
      if( !is.null(children) ){
        stack <- c(stack, children)
      }
    }
    return(seen)
  }
  treatment_sources <- intersect(paste0(treat_ids, "@t", treat_slice), unrolled$vertices)
  treatment_descendants <- union(treatment_sources, reach_forward(treatment_sources))
  descendant <- stats::setNames(unrolled$vertices %in% treatment_descendants, unrolled$vertices)

  # A base is SPLIT into occasions for each slice <id>_t{s} when it straddles the intervention: some
  # occasion is not a descendant of the treatment (before treatment, off the effect) and some occasion
  # is a descendant (the effect has reached it). A base that is entirely before treatment (a static / exogenous
  # confounder, one value throughout) or never straddles collapses to a single bare node. Treatment
  # and outcome never split. Each is anchored to a single slice below.
  base_ids <- unique(.strip_slice(unrolled$vertices))
  straddles <- function( b ){ # Check whether a base node straddles the treatment
    occasions <- intersect(paste0(b, "@t", seq_len(K) - 1L), unrolled$vertices)
    return(any(!descendant[ occasions ]) && any(descendant[ occasions ]))
  }
  straddlers <- setdiff(base_ids[ vapply(base_ids, straddles, logical(1L)) ], never_split)

  # Detect straddlers before hiding nodes so the test still sees paths through
  # hidden occasions. show_hidden_nodes = TRUE preserves those occasions for
  # later tagging.
  collapsed <- if( isTRUE(show_hidden_nodes) ){
    list(vertices = unrolled$vertices, edges = unrolled$edges)
  }else{
    .collapse_hidden(unrolled$vertices, unrolled$edges, hidden_vertex)
  }
  vertex_bases <- .strip_slice(collapsed$vertices)
  vertex_slices <- as.integer(sub("^.*@t([0-9]+)$", "\\1", collapsed$vertices))

  # Anchor the treatment to treat_slice and the outcome to out_slice: an intervention is applied
  # once, and the outcome is a terminal measurement read once, at the slice the effect arrives
  # (VanderWeele & Robins 2007). Their occasions at every other slice, and the edges on them, are
  # dropped before contracting, leaving a single bare treatment / outcome node.
  anchor_drop <- c(collapsed$vertices[ vertex_bases %in% treat_ids   & vertex_slices != treat_slice ],
                   collapsed$vertices[ vertex_bases %in% outcome_ids & vertex_slices != out_slice ])
  keep_v0 <- setdiff(collapsed$vertices, anchor_drop)

  # Project split variables to one node per occasion. Other variables use a
  # single node. If projection creates a cycle, split the affected variables
  # further and try again.
  #
  # The all mode starts with every included variable split. The pre_post mode
  # instead groups a straddler's occasions by whether they descend from
  # treatment. Its cycle handling can remove interior occasions or restore
  # individual occasion labels.
  post_treatment <- descendant[ collapsed$vertices ]
  promoted <- character(0L)                          # bases that fall back to per slice labelling (a fold could not break the cycle)
  folded_drop <- character(0L)                          # interior occasions marginalised out by the pre_post pre/post fold
  base_split_set <- if( isTRUE(split_all) ){ setdiff(base_ids, never_split) }else{ straddlers }
  slice_of <- stats::setNames(vertex_slices, collapsed$vertices)
  repeat {
    split_set <- union(base_split_set, promoted)
    keep_now <- setdiff(keep_v0, folded_drop)
    is_split <- vertex_bases %in% split_set
    proj_label <- vertex_bases
    if( identical(label_mode, "pre_post") ){
      # straddlers collapse to _pre_treatment / _post_treatment by descendant status. A base whose
      # fold could not break a residual cycle falls back to per slice (promoted). Interior occasions
      # dropped by the fold (folded_drop) are excluded from keep_now and marginalised out below.
      use_slice <- is_split & (vertex_bases %in% promoted)
      use_pp <- is_split & !(vertex_bases %in% promoted)
      proj_label[ use_slice ] <- paste0(vertex_bases[ use_slice ], "_t", vertex_slices[ use_slice ])
      proj_label[ use_pp ] <- paste0(vertex_bases[ use_pp ],
                                        ifelse(post_treatment[ use_pp ], "_post_treatment", "_pre_treatment"))
    }else{
      proj_label[ is_split ] <- paste0(vertex_bases[ is_split ], "_t", vertex_slices[ is_split ])
    }
    projection_map <- stats::setNames(proj_label, collapsed$vertices)[ keep_now ]
    kept_edges <- collapsed$edges[ collapsed$edges$from %in% keep_v0 & collapsed$edges$to %in% keep_v0 ]
    kept_edges <- .marginalise(kept_edges, folded_drop)
    kept_edges <- kept_edges[ kept_edges$from %in% keep_now & kept_edges$to %in% keep_now ]
    projected_edges <- .contract_edges(kept_edges, projection_map)
    projected_vertices <- unique(unname(projection_map))
    cyclic <- Filter(function( s ){ # Select components containing a directed cycle
      return(length(s) > 1L)
    }, .scc_of(projected_vertices, projected_edges))
    if( length(cyclic) == 0L ){
      break
    }
    cyc_bases <- setdiff(unique(.strip_slice(unlist(cyclic))), never_split)
    if( identical(label_mode, "pre_post") ){
      # FOLD the cyclic bases to pre/post rather than promoting them to per slice. Drop each cyclic
      # base's interior occasions (keep the earliest -> _pre_treatment and the latest -> _post_treatment).
      # Marginalising the interiors turns the cross slice feedback that had folded to a within group
      # cycle back into a forward _pre_treatment -> _post_treatment edge, which is acyclic. A base with
      # < 3 kept occasions has no interior to drop, so it cannot be folded. Fall back to per slice for
      # it and restore its interiors. The loop checks the next projection and
      # stops if it cannot make further progress, even if a cycle remains.
      progress <- FALSE
      for( b in setdiff(cyc_bases, promoted) ){
        occasions <- intersect(collapsed$vertices[ vertex_bases == b ], keep_now)
        occasions <- occasions[ order(slice_of[ occasions ]) ]
        if( length(occasions) >= 3L ){
          folded_drop <- union(folded_drop, occasions[ -c(1L, length(occasions)) ])
          progress <- TRUE
        }else{
          promoted <- union(promoted, b)
          folded_drop <- setdiff(folded_drop, collapsed$vertices[ vertex_bases == b ])
          progress <- TRUE
        }
      }
      if( !progress ){
        break
      }
    }else{
      new_promote <- setdiff(cyc_bases, split_set)
      if( length(new_promote) == 0L ){
        break
      }
      promoted <- union(promoted, new_promote)
    }
  }
  keep_now <- setdiff(keep_v0, folded_drop)

  # temporal persistence edges (returned separately, applied after any query_depth scoping). The
  # earlier occasion of a path dependent STATE variable causally precedes the next (forward Euler
  # identity: stock(t+1) = stock(t) + net_flux(t)). A flux or instantaneous ratio gets none. Emitted
  # between consecutive surviving occasions, always earlier -> later, so they cannot create a cycle.
  # Persistence edges are emitted in the SAME labelling the projection used, so they attach to real
  # vertices in either mode: each stock's surviving occasions are taken in slice order, mapped through
  # the final projection (pmap), deduplicated where consecutive occasions share a label (two slices in
  # the same pre/post bucket), and linked earlier -> later. In slice mode every occasion has its own
  # <id>_t{s}, giving an <id>_t{s} -> <id>_t{s+1} chain. In pre_post mode a straddler stock
  # carries a single <id>_pre_treatment -> <id>_post_treatment. Always forward, so it cannot add a cycle.
  persist_from <- character(0L)
  persist_to <- character(0L)
  for( base in intersect(split_set, persistent_ids) ){
    occasion_vertices <- intersect(collapsed$vertices[ vertex_bases == base ], keep_now)
    if( length(occasion_vertices) >= 2L ){
      occasion_labels <- unname(projection_map[ occasion_vertices[ order(slice_of[ occasion_vertices ]) ] ])            # projected labels, slice order
      occasion_labels <- occasion_labels[ c(TRUE, occasion_labels[-1L] != occasion_labels[-length(occasion_labels)]) ]           # drop consecutive duplicates
      if( length(occasion_labels) >= 2L ){
        persist_from <- c(persist_from, occasion_labels[ -length(occasion_labels) ])
        persist_to <- c(persist_to,   occasion_labels[ -1L ])
      }
    }
  }
  persistence <- unique(data.table::data.table(from = persist_from, to = persist_to))

  return(list(vertices = projected_vertices, edges = projected_edges, split = split_set, promoted = promoted,
              treat_slice = treat_slice, out_slice = out_slice, reachable = depth$reachable,
              K = K, persistence = persistence))
}


#' Expand growth-stage requests into endpoint-specific tokens
#'
#' @param spec `NULL`, stage tokens to apply to every endpoint, or a list
#'   named by endpoint IDs. Endpoints absent from a supplied list are omitted.
#' @param node_ids Resolved endpoint IDs allowed as list names.
#' @return A named list of trimmed character stage vectors, or an empty named
#'   list for `NULL` or an empty specification. Unknown endpoint names raise
#'   an error.
#' @noRd
.normalize_stage_spec <- function( spec, node_ids ){
  if( is.null(spec) || length(spec) == 0L ){
    return(stats::setNames(list(), character(0L)))
  }
  if( is.list(spec) ){
    nm <- names(spec)
    if( is.null(nm) || any(!nzchar(nm)) ){
      stop("A growth stage list must be named by treatment/outcome node id, e.g. list(nitrogen_applied = \"flowering\"). Please name every element and try again.")
    }
    unknown <- setdiff(nm, node_ids)
    if( length(unknown) > 0L ){
      stop(sprintf("Growth stage list names must be among the resolved treatment/outcome nodes (%s). Unrecognised: %s. Please correct the names and try again.",
                   paste(node_ids, collapse = ", "), paste(unknown, collapse = ", ")))
    }
    return(lapply(spec, function( v ) { # Trim the stage tokens supplied for one endpoint
      return(trimws(as.character(v)))
    }))
  }
  stages <- trimws(as.character(spec))
  return(stats::setNames(rep(list(stages), length(node_ids)), node_ids))
}


# Resolve a per node canonical stage map to per node integer phase indices for `crop`, via crop_phases.
# Each token is matched against the crop's canonical_stage vocabulary (the interpretable, public
# handle). A bare integer token is accepted as a literal phase index, and an APSIM phase name resolves
# too (both supported in the backend). Returns a named list node_id -> sorted unique integer indices.
# Errors name the valid canonical stages for the crop.
#' Resolve endpoint stage tokens in a selected crop partition
#'
#' @param stage_map Named list of stage tokens returned by `.normalize_stage_spec()`.
#' @param crop Selected crop with an available phase partition.
#' @param crop_phases Table of crop, phase index, phase name and canonical stage.
#' @return A list with the same names and sorted unique integer phase indices.
#'   Unknown tokens and indices outside the crop partition raise errors.
#' @noRd
.stages_to_phase_indices <- function( stage_map, crop, crop_phases ){
  crop_partition <- crop_phases[ crop_phases$crop == crop, ]
  if( nrow(crop_partition) == 0L ){
    stop(sprintf("crop '%s' has no crop_phases partition, so growth stages cannot be resolved. Please supply a crop with a phase partition and try again.", crop))
  }
  has_canonical_stage <- nzchar(crop_partition$canonical_stage)
  canonical_to_index <- stats::setNames(crop_partition$phase_index[ has_canonical_stage ], crop_partition$canonical_stage[ has_canonical_stage ])
  name_to_index <- stats::setNames(crop_partition$phase_index, crop_partition$phase_name)
  valid_stages <- sort(unique(crop_partition$canonical_stage[ has_canonical_stage ]))
  resolve_one <- function( tok ){ # Resolve and validate one growth stage token
    if( grepl("^[0-9]+$", tok) ){
      phase_index <- as.integer(tok)
      if( !(phase_index %in% crop_partition$phase_index) ){
        stop(sprintf("phase index %d is out of range for %s (valid 0..%d). Please supply a valid phase and try again.",
                     phase_index, crop, max(crop_partition$phase_index)))
      }
      return(phase_index)
    }
    if( tok %in% names(canonical_to_index) ){ return(unname(canonical_to_index[ tok ])) }
    if( tok %in% names(name_to_index) ){ return(unname(name_to_index[ tok ])) }
    stop(sprintf("growth stage '%s' is not a canonical stage for %s. Valid canonical stages: %s. Please use one of these (or a phase index) and try again.",
                 tok, crop, paste(valid_stages, collapse = ", ")))
  }
  return(lapply(stage_map, function( toks ) { # Resolve the distinct phases requested for one endpoint
    return(sort(unique(vapply(toks, resolve_one, integer(1L)))))
  } ))
}


# Phase anchored selective slice: the developmental phase sibling of .selective_slice. Instead of a
# floating window of K slices anchored to treatment it unrolls the crop's phase partition (via
# .unroll_phase_slices, so the initial network B0 is in the frame with the sowing origin at slice 0),
# and anchors each treatment and outcome at its selected developmental phases.
# A base is an endpoint only at those phases. Outcome histories remain distinct states.
# Splitting, hidden node collapse, cycle promotion and persistence then use these phase occasions.
# `window` chooses the retained span: "origin" keeps slices 0..hi (initial conditions through B0).
# "compact" keeps (earliest treatment phase - 1)..hi (one phase of confounder history before treatment,
# no origin/B0). Returns the visible projected graph plus the split/promoted sets, the window, and
# persistence edges. This is the same shape .selective_slice returns, so the dagby_query tail consumes it as is.
#' Project developmental phases while preserving endpoint occasions
#'
#' @param kept_ids Quantity IDs retained in the configured graph.
#' @param contemporaneous,lagged,initial Directed `from`/`to` data.tables.
#' @param hidden_ids Hidden quantity IDs eligible for directed removal.
#' @param treat_phase_map,out_phase_map Lists of selected phase indices by
#'   treatment and outcome ID. Only the selected occasions are endpoint tags.
#' @param n_phases Positive integer size of the crop's phase partition.
#' @param persistent_ids Quantity IDs whose retained occasions need carryover.
#' @param window `origin` retains sowing. `compact` starts one phase before
#'   the earliest treatment, bounded below by zero.
#' @param phase_active Optional inclusive active intervals by quantity ID.
#' @param show_hidden_nodes Retain hidden quantities when `TRUE`.
#' @return A list with projected `vertices` and `edges`, the quantity IDs in
#'   `split` and `promoted`, retained `window`, `n_phases`, `reachable = TRUE`,
#'   `persistence` arrows, and `treatments` and `outcomes` mapped to their
#'   retained occasions.
#' @noRd
.selective_phase_slice <- function( kept_ids, contemporaneous, lagged, initial, hidden_ids,
                                    treat_phase_map, out_phase_map, n_phases,
                                    persistent_ids = character(0L), window = c("origin", "compact"),
                                    phase_active = NULL, show_hidden_nodes = FALSE ){
  .datatable.aware <- TRUE
  window <- match.arg(window)
  treat_ids <- names(treat_phase_map)
  outcome_ids <- names(out_phase_map)
  # Single anchors used only as treatments can remain bare. An outcome's other occasions
  # are still causal states: keep all of them split so they cannot be merged
  # into the measured occasion. This also separates same base treatment/outcome
  # anchors. Final single occasion relabelling may remove an unnecessary suffix.
  single_anchor <- c(treat_ids[ lengths(treat_phase_map) == 1L ], outcome_ids[ lengths(out_phase_map) == 1L ])
  multi_anchor <- c(treat_ids[ lengths(treat_phase_map) >  1L ], outcome_ids[ lengths(out_phase_map) >  1L ])
  never_split <- setdiff(single_anchor, outcome_ids)

  # slice window [lo..hi]: hi = latest anchored phase, lo = origin (0) under "origin", else one phase of
  # history before treatment (earliest treatment phase - 1, floored at 0) under "compact".
  all_anchor <- c(unlist(treat_phase_map, use.names = FALSE), unlist(out_phase_map, use.names = FALSE))
  last_slice <- max(all_anchor)
  first_slice <- if( identical(window, "origin") ){ 0L }else{ max(0L, min(unlist(treat_phase_map, use.names = FALSE)) - 1L) }

  # full phase unroll (working graph, never emitted). edge_type dropped: the recovered query graph is a
  # plain DAG, and B0 init edges flow through as ordinary from->to edges into the origin slice stock.
  full <- .unroll_phase_slices(kept_ids, contemporaneous, lagged, initial, n_phases,
                                     phase_active = phase_active)
  all_vertex_slices <- as.integer(sub("^.*@t([0-9]+)$", "\\1", full$vertices))
  window_vertices <- full$vertices[ all_vertex_slices >= first_slice & all_vertex_slices <= last_slice ]
  window_edges <- full$edges[ full$edges$from %in% window_vertices & full$edges$to %in% window_vertices, list(from, to) ]

  hidden_vertex <- stats::setNames(.strip_slice(window_vertices) %in% hidden_ids, window_vertices)

  # descendants of treatment over the windowed edges, from every treatment's anchored phase occasions
  children_by_parent <- split(window_edges$to, window_edges$from)
  reach_forward <- function( sources ){ # Collect descendants reached from the selected treatment occasions
    seen <- character(0L)
    stack <- unlist(children_by_parent[ sources ], use.names = FALSE)
    while( length(stack) > 0L ){
      vertex <- stack[ length(stack) ]
      stack <- stack[ -length(stack) ]
      if( vertex %in% seen ){
        next
      }
      seen <- c(seen, vertex)
      children <- children_by_parent[[ vertex ]]
      if( !is.null(children) ){
        stack <- c(stack, children)
      }
    }
    return(seen)
  }
  treatment_sources <- character(0L)
  for( t in treat_ids ){ treatment_sources <- c(treatment_sources, paste0(t, "@t", treat_phase_map[[ t ]])) }
  treatment_sources <- intersect(treatment_sources, window_vertices)
  treatment_descendants <- union(treatment_sources, reach_forward(treatment_sources))
  descendant <- stats::setNames(window_vertices %in% treatment_descendants, window_vertices)

  # straddlers: a base with occasions both outside and within treatment's descendants in the window
  base_ids <- unique(.strip_slice(window_vertices))
  window_slices <- first_slice:last_slice
  straddles <- function( b ){ # Check whether a base node straddles the treatment
    node_occasions <- intersect(paste0(b, "@t", window_slices), window_vertices)
    return(any(!descendant[ node_occasions ]) && any(descendant[ node_occasions ]))
  }
  straddlers <- setdiff(base_ids[ vapply(base_ids, straddles, logical(1L)) ], never_split)

  # collapse hidden on the windowed (ordered in time) graph (skip when showing hidden nodes)
  collapsed <- if( isTRUE(show_hidden_nodes) ){
    list(vertices = window_vertices, edges = window_edges)
  }else{
    .collapse_hidden(window_vertices, window_edges, hidden_vertex)
  }
  vertex_bases <- .strip_slice(collapsed$vertices)
  vertex_slices <- as.integer(sub("^.*@t([0-9]+)$", "\\1", collapsed$vertices))

  # Preserve the existing event selection rule for treatment only bases. Never
  # delete an outcome's other occasions: they may mediate its response or be
  # common causes. Selecting an endpoint does not erase its state history.
  drop_anchor <- character(0L)
  for( t in setdiff(treat_ids, outcome_ids) ){
    drop_anchor <- c(drop_anchor, collapsed$vertices[ vertex_bases == t & !(vertex_slices %in% treat_phase_map[[ t ]]) ])
  }
  keep_v0 <- setdiff(collapsed$vertices, drop_anchor)

  # Split straddlers, outcome histories and endpoints selected at several phases.
  # If projection creates a cycle, split further quantities where possible.
  # Stop when acyclic or no further split is available, even if a cycle remains.
  promoted <- character(0L)
  split_set <- union(union(straddlers, multi_anchor), outcome_ids)
  repeat {
    proj_label <- vertex_bases
    is_split <- vertex_bases %in% split_set
    proj_label[ is_split ] <- paste0(vertex_bases[ is_split ], "_t", vertex_slices[ is_split ])
    projection_map <- stats::setNames(proj_label, collapsed$vertices)[ keep_v0 ]
    kept_edges <- collapsed$edges[ collapsed$edges$from %in% keep_v0 & collapsed$edges$to %in% keep_v0 ]
    projected_edges <- .contract_edges(kept_edges, projection_map)
    projected_vertices <- unique(unname(projection_map))
    cyclic <- Filter(function( s ){ # Select components containing a directed cycle
      return(length(s) > 1L)
    }, .scc_of(projected_vertices, projected_edges))
    if( length(cyclic) == 0L ){
      break
    }
    promote <- setdiff(unique(.strip_slice(unlist(cyclic))), union(never_split, split_set))
    if( length(promote) == 0L ){
      break
    }
    promoted <- union(promoted, promote)
    split_set <- union(union(union(straddlers, multi_anchor), outcome_ids), promoted)
  }

  # persistence for split path dependent stocks (earlier phase occasion -> next), forward only
  persist_from <- character(0L)
  persist_to <- character(0L)
  for( base in intersect(split_set, persistent_ids) ){
    occasion_slices <- sort(unique(vertex_slices[ vertex_bases == base & (collapsed$vertices %in% keep_v0) ]))
    if( length(occasion_slices) >= 2L ){
      for( i in seq_len(length(occasion_slices) - 1L) ){
        persist_from <- c(persist_from, paste0(base, "_t", occasion_slices[i]))
        persist_to <- c(persist_to,   paste0(base, "_t", occasion_slices[i + 1L]))
      }
    }
  }
  persistence <- data.table::data.table(from = persist_from, to = persist_to)

  endpoint_vertices <- function(phase_map){ # Map selected endpoint phases to the retained graph vertices
    occasions <- unlist(lapply(names(phase_map), function(id) { # Name the selected occasions of one endpoint
      return(paste0(id, "@t", phase_map[[id]]))
    }), use.names = FALSE)
    return(unique(unname(projection_map[intersect(occasions, names(projection_map))])))
  }
  return(list(vertices = projected_vertices, edges = projected_edges, split = split_set, promoted = promoted,
              window = c(first_slice, last_slice), n_phases = n_phases, reachable = TRUE, persistence = persistence,
              treatments = endpoint_vertices(treat_phase_map), outcomes = endpoint_vertices(out_phase_map)))
}


# Directed collapse: remove hidden vertices, and for each visible vertex add an
# edge to every visible vertex reachable from it through hidden interiors
# (any chain length). No bidirected edges are produced. Isolated visible
# vertices are retained (declared by the emitter) but gain no edges.
#' Record a directed fold through removed vertices
#'
#' @param vertices Unique graph vertex names.
#' @param edges Directed `from`/`to` data.table on those vertices.
#' @param remove Vertex names to remove, intersected with the graph.
#' @return A list of retained `vertices` and `edges`, copied source structures,
#'   removed nodes, removed and added arrows, and an operation description.
#'   The operation reconnects directed paths. It does not perform latent projection.
#' @noRd
.query_fold_nodes <- function(vertices, edges, remove){
  remove <- intersect(vertices, remove)
  result <- if( length(remove) ){
    .collapse_hidden(vertices, edges, stats::setNames(vertices %in% remove, vertices))
  }else{ list(vertices = vertices, edges = data.table::copy(edges)) }
  edge_key <- function(e) { # Build an unambiguous key for each directed edge
    return(paste(e$from, e$to, sep = "\r"))
  }
  result$source_vertices <- vertices
  result$source_edges <- data.table::copy(edges)
  result$removed_nodes <- remove
  result$removed_edges <- edges[!(edge_key(edges) %in% edge_key(result$edges)), ]
  result$added_edges <- result$edges[!(edge_key(result$edges) %in% edge_key(edges)), ]
  result$operation <- "Directed paths through removed interiors. No latent projection"
  return(result)
}


#' Replace directed paths through hidden interiors with arrows
#'
#' @param vertices Unique graph vertex names, including isolated vertices.
#' @param edges Directed `from`/`to` data.table on those vertices.
#' @param hidden_vertex Nonmissing logical vector named by every vertex.
#' @return A list with visible `vertices` and directed `edges`. No bidirected
#'   arrows are introduced. Self-loops and duplicate arrows are removed.
#' @noRd
.collapse_hidden <- function( vertices, edges, hidden_vertex ){
  .datatable.aware <- TRUE
  children_by_parent <- split(edges$to, edges$from)
  visible <- vertices[ !hidden_vertex[ vertices ] ]
  records <- lapply(visible, function( source ){ # Find visible targets reachable through hidden interiors
    seen <- character(0L)
    stack <- children_by_parent[[ source ]]
    targets <- character(0L)
    while( length(stack) > 0L ){
      vertex <- stack[ length(stack) ]
      stack <- stack[ -length(stack) ]
      if( vertex %in% seen ){
        next
      }
      seen <- c(seen, vertex)
      if( hidden_vertex[[ vertex ]] ){
        children <- children_by_parent[[ vertex ]]
        if( !is.null(children) ){
          stack <- c(stack, children)
        }
      }else{
        targets <- c(targets, vertex)
      }
    }
    targets <- unique(targets)
    if( length(targets) == 0L ){
      return(NULL)
    }
    return(data.table::data.table(from = source, to = targets))
  })
  kept <- Filter(Negate(is.null), records)
  if( length(kept) == 0L ){
    edges_out <- data.table::data.table(from = character(0L), to = character(0L))
  }else{
    edges_out <- unique(data.table::rbindlist(kept))
    edges_out <- edges_out[ from != to ]
  }
  return(list(vertices = visible, edges = edges_out))
}


# Truncate dotted node ids to the requested APSIM detail level: "supernode" keeps
# the segment before the first dot, "node" the first two segments, "subnode" the
# first three. Ids with fewer segments (or no dot) pass through unchanged. overrides is
# an optional named vector (supernode -> level) setting a per supernode level that
# replaces the global one for that first segment group (supernode_overrides).
# Returns a node_id -> grouped_id map (a character vector named by node_id).
#' Map dotted identifiers to the requested detail level
#'
#' @param ids Nonempty quantity identifiers.
#' @param level One of `supernode`, `node`, `subnode` or `variable`.
#' @param overrides Optional detail levels named by the first dotted segment.
#' @return Character group names indexed by original IDs. IDs with fewer
#'   segments than the requested depth remain unchanged.
#' @noRd
.apsim_truncate <- function( ids, level, overrides = NULL ){
  base_depth <- switch(level, "supernode" = 1L, "node" = 2L, "subnode" = 3L, "variable" = Inf)
  depth_of <- c(supernode = 1L, node = 2L, subnode = 3L, variable = Inf)
  return(sapply(ids, function( id ){ # Apply the global or overridden detail level to one identifier
    parts <- strsplit(id, ".", fixed = TRUE)[[ 1L ]]
    depth <- base_depth
    if( length(overrides) > 0L && parts[ 1L ] %in% names(overrides) ){
      depth <- unname(depth_of[ overrides[[ parts[ 1L ] ]] ])
    }
    if( length(parts) <= depth ){
      return(id)
    }
    return(paste(parts[ 1:depth ], collapse = "."))
  }))
}


#' Contract directed edges through a complete grouping map
#'
#' @param edges Directed `from`/`to` data.table.
#' @param group_map Group names indexed by every endpoint in `edges`.
#' @return Unique inter-group directed edges, with self-loops removed. An empty
#'   input edge table is returned unchanged.
#' @noRd
.contract_edges <- function( edges, group_map ){
  .datatable.aware <- TRUE
  if( nrow(edges) == 0L ){
    return(edges)
  }
  contracted_edges <- data.table::data.table(from = unname(group_map[ edges$from ]),
                                to   = unname(group_map[ edges$to ]))
  contracted_edges <- contracted_edges[ from != to ]
  return(unique(contracted_edges))
}


#' Reconnect directed parents and children while removing vertices
#'
#' Acyclicity of the result alone does not establish causal equivalence.
#'
#' @param edges Directed `from`/`to` data.table.
#' @param drop Vertex names to remove in the supplied order.
#' @return Remaining and reconnected directed arrows. This adds no bidirected
#'   arrows and does not generally preserve d-separation.
#' @noRd
.marginalise <- function( edges, drop ){
  .datatable.aware <- TRUE
  if( nrow(edges) == 0L || length(drop) == 0L ){
    return(edges)
  }
  remaining_edges <- data.table::as.data.table(edges)
  for( removed_vertex in drop ){
    parents <- unique(remaining_edges$from[ remaining_edges$to == removed_vertex ])
    children <- unique(remaining_edges$to[ remaining_edges$from == removed_vertex ])
    remaining_edges <- remaining_edges[ from != removed_vertex & to != removed_vertex ]
    if( length(parents) > 0L && length(children) > 0L ){
      reconnected_edges <- data.table::CJ(from = parents, to = children)
      reconnected_edges <- reconnected_edges[ from != to ]
      remaining_edges <- unique(data.table::rbindlist(list(remaining_edges, reconnected_edges)))
    }
  }
  return(remaining_edges)
}


# Two vertices belong to the same component if each can reach the other.
# After .contract_edges() removes self-loops, a component with at least two
# vertices identifies a cycle that grouping must resolve.
#' Partition a directed graph into strongly connected components
#'
#' @param labels Vertex names, including isolated vertices.
#' @param edges Directed `from`/`to` edge table. Extra endpoints participate in
#'   traversal, but do not extend the requested labels.
#' @return A list of character component vectors containing only `labels`,
#'   including singleton components.
#' @noRd
.scc_of <- function( labels, edges ){
  .datatable.aware <- TRUE
  if( nrow(edges) == 0L ){
    return(as.list(labels))
  }

  vertices <- unique(c(labels, as.character(edges$from), as.character(edges$to)))
  vertex_count <- length(vertices)
  vertex_index <- stats::setNames(seq_len(vertex_count), vertices)
  from <- as.integer(vertex_index[ as.character(edges$from) ])
  to <- as.integer(vertex_index[ as.character(edges$to) ])
  children_by_parent <- split(to, factor(from, levels = seq_len(vertex_count)))

  index <- integer(vertex_count)
  lowlink <- integer(vertex_count)
  component <- integer(vertex_count)
  on_stack <- logical(vertex_count)

  stack <- integer(vertex_count)
  stack_top <- 0L
  vertex_stack <- integer(vertex_count + 1L)
  child_stack <- integer(vertex_count + 1L)
  counter <- 0L
  component_count <- 0L

  for( root in seq_len(vertex_count) ){
    if( index[ root ] != 0L ){
      next
    }
    counter <- counter + 1L
    index[ root ] <- counter
    lowlink[ root ] <- counter
    stack_top <- stack_top + 1L
    stack[ stack_top ] <- root
    on_stack[ root ] <- TRUE
    depth <- 1L
    vertex_stack[ 1L ] <- root
    child_stack[ 1L ] <- 1L

    while( depth > 0L ){
      vertex <- vertex_stack[ depth ]
      children <- children_by_parent[[ vertex ]]
      child_index <- child_stack[ depth ]
      if( child_index <= length(children) ){
        child_stack[ depth ] <- child_index + 1L
        child <- children[ child_index ]
        if( index[ child ] == 0L ){
          counter <- counter + 1L
          index[ child ] <- counter
          lowlink[ child ] <- counter
          stack_top <- stack_top + 1L
          stack[ stack_top ] <- child
          on_stack[ child ] <- TRUE
          depth <- depth + 1L
          vertex_stack[ depth ] <- child
          child_stack[ depth ] <- 1L
        }else if( on_stack[ child ] ){
          if( index[ child ] < lowlink[ vertex ] ){
            lowlink[ vertex ] <- index[ child ]
          }
        }
      }else{
        if( lowlink[ vertex ] == index[ vertex ] ){
          component_count <- component_count + 1L
          repeat {
            child <- stack[ stack_top ]
            stack_top <- stack_top - 1L
            on_stack[ child ] <- FALSE
            component[ child ] <- component_count
            if( child == vertex ){
              break
            }
          }
        }
        depth <- depth - 1L
        if( depth > 0L ){
          parent <- vertex_stack[ depth ]
          if( lowlink[ vertex ] < lowlink[ parent ] ){
            lowlink[ parent ] <- lowlink[ vertex ]
          }
        }
      }
    }
  }

  label_components <- component[ as.integer(vertex_index[ labels ]) ]
  components <- vector("list", 0L)
  emitted <- logical(component_count)
  for( k in seq_along(labels) ){
    component_id <- label_components[ k ]
    if( emitted[ component_id ] ){
      next
    }
    emitted[ component_id ] <- TRUE
    components[[ length(components) + 1L ]] <- labels[ label_components == component_id ]
  }
  return(components)
}


#' Refine cyclic dotted groups within the available detail levels
#'
#' @param node_ids All source quantity IDs used to define grouping depths.
#' @param kept_ids IDs included in the graph being checked.
#' @param contemp Directed `from`/`to` data.table before grouping.
#' @param level Requested global dotted detail level.
#' @param base_overrides Optional initial detail levels by first dotted segment.
#' @param pinned_supers First-segment groups to avoid refining when an effective
#'   unpinned split is available.
#' @return Detail overrides relative to `level`, or `NULL` when no group differs
#'   from that level. Refinement stops when the graph is acyclic or no finer
#'   group is available. The caller must still check acyclicity. This greedy
#'   search does not guarantee the coarsest acyclic partition.
#' @noRd
.refine_to_acyclic <- function( node_ids, kept_ids, contemp, level, base_overrides, pinned_supers ){
  .datatable.aware <- TRUE
  depth_of <- c(supernode = 1L, node = 2L, subnode = 3L, variable = 4L)
  first_segment <- function( ids ){ # Extract the leading segment of each node identifier
    return(sapply(strsplit(ids, ".", fixed = TRUE), `[`, 1L))
  }
  segment_count <- function( ids ){ # Count the dotted segments in each node identifier
    return(lengths(strsplit(ids, ".", fixed = TRUE)))
  }

  supers <- unique(first_segment(node_ids))
  first_segments <- first_segment(node_ids)
  segment_counts <- segment_count(node_ids)
  maximum_depths <- stats::setNames(sapply(supers, function( s ){ # Find the finest available depth of this supernode
    return(min(4L, max(segment_counts[ first_segments == s ])))
  }), supers)

  current_depths <- stats::setNames(rep(unname(depth_of[ level ]), length(supers)), supers)
  if( length(base_overrides) > 0L ){
    overridden_groups <- intersect(names(base_overrides), supers)
    current_depths[ overridden_groups ] <- unname(depth_of[ base_overrides[ overridden_groups ] ])
  }
  to_overrides <- function( candidate_depths ){ # Keep only depths that differ from the global level
    base <- unname(depth_of[ level ])
    keep <- candidate_depths[ candidate_depths != base ]
    if( length(keep) == 0L ){
      return(NULL)
    }
    return(stats::setNames(names(depth_of)[ match(keep, depth_of) ], names(keep)))
  }
  # the cyclic state at a given per supernode depth: number of vertices left in a cycle
  # (0 = acyclic), with the contracted graph and the cyclic vertex set.
  cyclic_state <- function( candidate_depths ){ # Measure the cyclic vertices under a candidate grouping
    group_map <- .apsim_truncate(node_ids, level, overrides = to_overrides(candidate_depths))
    grouped_vertices <- unique(unname(group_map[ kept_ids ]))
    grouped_edges <- .contract_edges(contemp[ from %in% kept_ids & to %in% kept_ids ], group_map)
    if( .is_acyclic(grouped_vertices, grouped_edges) ){
      return(list(n = 0L, gmap = group_map, cyclic = character(0L)))
    }
    scc <- .scc_of(grouped_vertices, grouped_edges)
    return(list(n = length(unlist(scc[ sapply(scc, length) >= 2L ])),
                gmap = group_map, cyclic = unlist(scc[ sapply(scc, length) >= 2L ])))
  }

  repeat {
    cyclic_graph <- cyclic_state(current_depths)
    if( cyclic_graph$n == 0L ){
      break
    }
    cyc_supers <- unique(first_segment(names(cyclic_graph$gmap)[ unname(cyclic_graph$gmap) %in% cyclic_graph$cyclic ]))
    refinable <- cyc_supers[ current_depths[ cyc_supers ] < maximum_depths[ cyc_supers ] ]
    if( length(refinable) == 0L ){
      break  # unreachable while the finest graph is acyclic. Guards against a spin
    }
    # Prefer splits that reduce the cyclic set, first unpinned and then
    # pinned. If none helps immediately, refine an unpinned group when
    # possible, otherwise any refinable group.
    helps <- sapply(refinable, function( s ){ # Check whether refining this supernode reduces the cyclic set
      trial <- current_depths
      trial[ s ] <- trial[ s ] + 1L
      return(cyclic_state(trial)$n < cyclic_graph$n)
    })
    effective_groups <- refinable[ helps ]
    effective_unpinned <- setdiff(effective_groups, pinned_supers)
    nonpinned <- setdiff(refinable, pinned_supers)
    pool <- if( length(effective_unpinned) > 0L ){
              effective_unpinned
            }else if( length(effective_groups) > 0L ){
              effective_groups
            }else if( length(nonpinned) > 0L ){
              nonpinned
            }else{
              refinable
            }
    pick <- sort(pool)[ 1L ]
    current_depths[ pick ] <- current_depths[ pick ] + 1L
  }
  return(to_overrides(current_depths))
}


#' Name a group by its common dotted prefix
#'
#' @param ids Nonempty character vector of nonempty node identifiers.
#' @return The common leading dotted segments, or the first ID if none are shared.
#'   The caller must resolve collisions between different groups.
#' @noRd
.dotted_lcp <- function( ids ){
  parts <- strsplit(ids, ".", fixed = TRUE)
  n <- min(lengths(parts))
  prefix <- character(0L)
  for( i in 1:n ){                              # n >= 1L: node ids are nonempty
    segment <- unique(sapply(parts, function( p ){ # Read this segment from one split identifier
      return(p[ i ])
    }))
    if( length(segment) == 1L ){
      prefix <- c(prefix, segment)
    }else{
      break
    }
  }
  if( length(prefix) == 0L ){
    # no shared leading segment: fall back to the first id. Two disjoint groups can hit this
    # and collide, so .esc_recombine's #k disambiguation must run downstream to keep names unique
    return(ids[ 1 ])
  }
  return(paste(prefix, collapse = "."))
}


# ESC-DAGs recombination (Ferguson et al. 2019, Integration 2): merge vertices that
# share an identical set of parents AND an identical set of children into one node.
# This preserves acyclicity. Two such vertices can never lie on a directed path
# through other vertices (that would need a child of the set to reach a parent of the
# set, impossible in a DAG), so the contracted graph stays a DAG. Fully isolated
# vertices and any pinned vertices (named treatments/outcomes) are never merged. Each
# merged node is labelled by the longest common dotted prefix of its members (slice
# suffix reattached, made unique with a #k suffix). Latent carries over only when
# every member of a merged node was latent.
#' Recombine unpinned quantities with identical directed neighbours
#'
#' @param vertices Unique vertex names of the input DAG.
#' @param edges Directed `from`/`to` data.table on those vertices.
#' @param latent_set Vertex names currently tagged latent.
#' @param pinned Vertex names that must remain separate.
#' @return A list with `vertices`, `edges`, `latent` and a `supernodes` member
#'   list. A merged node is latent only when all members were latent. Isolated
#'   vertices remain separate.
#' @noRd
.esc_recombine <- function( vertices, edges, latent_set, pinned = character(0L) ){
  parents_by_child <- split(edges$from, factor(edges$to,   levels = vertices))
  children_by_parent <- split(edges$to,   factor(edges$from, levels = vertices))
  signatures <- sapply(vertices, function( v ){ # Encode this vertex's parent and child sets without collisions
    parents <- if( is.null(parents_by_child[[ v ]]) ){ character(0L) }else{ sort(unique(parents_by_child[[ v ]])) }
    children <- if( is.null(children_by_parent[[ v ]]) ){ character(0L) }else{ sort(unique(children_by_parent[[ v ]])) }
    # "\u0001"/"\u0002" cannot occur in a node id, and prefixing with the length keeps the parent
    # and child blocks unambiguous, so structurally distinct vertices never collide here
    return(paste0(length(parents), "\u0001", paste(parents, collapse = "\u0001"),
                  "\u0002", length(children), "\u0001", paste(children, collapse = "\u0001")))
  })
  empty_sig <- paste0(0L, "\u0001", "", "\u0002", 0L, "\u0001", "")
  groups <- lapply(split(vertices, signatures), function( group ){ # Exclude pinned vertices from a possible merge group
    return(setdiff(group, pinned))
  })
  classes <- Filter(function( group ){ # Keep connected groups with multiple members
    return(length(group) >= 2L && signatures[[ group[ 1 ] ]] != empty_sig)
  }, groups)

  map <- stats::setNames(vertices, vertices)
  members <- list()
  taken <- setdiff(vertices, unlist(classes))
  for( group in classes ){
    bases <- .strip_slice(group)
    suffix <- unique(regmatches(group, regexpr("(_pre_treatment|_post_treatment|_t[0-9]+|@t[0-9]+)$", group)))
    suffix <- if( length(suffix) == 1L ){ suffix }else{ "" }
    name <- paste0(.dotted_lcp(bases), suffix)
    base <- name
    k <- 1L
    # .dotted_lcp falls back to its first member when a group shares no dotted prefix, so
    # two disjoint groups can propose the same name. This #k disambiguation is needed
    # for label uniqueness and must not be removed.
    while( name %in% taken ){
      k <- k + 1L
      name <- paste0(base, "#", k)
    }
    taken <- c(taken, name)
    map[ group ] <- name
    members[[ name ]] <- group
  }
  new_vertices <- unique(unname(map[ vertices ]))
  new_edges <- .contract_edges(edges, map)
  latent_by_group <- tapply(vertices %in% latent_set, unname(map[ vertices ]), all)
  new_latent <- names(which(latent_by_group))
  return(list(vertices = new_vertices, edges = new_edges, latent = new_latent, supernodes = members))
}


# Broader query depths use .scope_to_query()
#' Retain the restricted depth-one query scope
#'
#' @param recovered_vertices Unique vertex names in the recovered graph.
#' @param recovered_edges Directed `from`/`to` table on those vertices.
#' @param treat_v,out_v Treatment and outcome occurrence names.
#' @param pinned Vertex names that must remain included.
#' @param quiet Suppress the retained-node summary when `TRUE`.
#' @return A list of scoped `vertices` and induced `edges`. Directed endpoint
#'   paths, endpoint parents and ancestral occasions of those quantities are
#'   retained. This restricted scope does not certify adjustment validity.
#' @noRd
.scope_adjust <- function( recovered_vertices, recovered_edges, treat_v, out_v, pinned, quiet ){
  from <- as.character(recovered_edges$from)
  to <- as.character(recovered_edges$to)
  parents_by_child <- split(from, to)
  children_by_parent <- split(to, from)

  # Intersect forward reach from treatment with backward reach from outcome to
  # find directed path nodes
  forward_reach <- character(0L)
  stack <- treat_v
  while( length(stack) > 0L ){
    vertex <- stack[[ length(stack) ]]
    stack <- stack[ -length(stack) ]
    if( vertex %in% forward_reach ){ next }
    forward_reach <- c(forward_reach, vertex)
    stack <- c(stack, children_by_parent[[ vertex ]])
  }
  backward_reach <- character(0L)
  stack <- out_v
  while( length(stack) > 0L ){
    vertex <- stack[[ length(stack) ]]
    stack <- stack[ -length(stack) ]
    if( vertex %in% backward_reach ){ next }
    backward_reach <- c(backward_reach, vertex)
    stack <- c(stack, parents_by_child[[ vertex ]])
  }
  path_nodes <- intersect(forward_reach, backward_reach)                       # on at least one directed X -> Y path


  treatment_parents <- unique(unlist(parents_by_child[ treat_v ], use.names = FALSE))
  outcome_parents <- unique(unlist(parents_by_child[ out_v ],   use.names = FALSE))

  base_scope <- unique(c(path_nodes, treatment_parents, outcome_parents, treat_v, out_v))
  base_scope <- intersect(base_scope, recovered_vertices)


  endpoint_ancestors <- character(0L)
  stack <- unique(c(treat_v, out_v))
  while( length(stack) > 0L ){
    vertex <- stack[[ length(stack) ]]
    stack <- stack[ -length(stack) ]
    if( vertex %in% endpoint_ancestors ){ next }
    endpoint_ancestors <- c(endpoint_ancestors, vertex)
    stack <- c(stack, parents_by_child[[ vertex ]])
  }

  # Retain ancestral occasions of every selected variable so omitted
  # intermediate occasions cannot break a retained temporal path. Keep pinned
  # nodes as well.
  bases <- unique(.strip_slice(base_scope))
  current <- recovered_vertices[ .strip_slice(recovered_vertices) %in% bases ]
  current <- intersect(current, endpoint_ancestors)
  current <- unique(c(current, intersect(pinned, recovered_vertices)))

  keep <- recovered_edges$from %in% current & recovered_edges$to %in% current
  if( isFALSE(quiet) ){
    message(length(current), " of ", length(recovered_vertices),
            " nodes retained (adjust scope: causal paths + parents of treatment/outcome).")
  }
  return(list(vertices = current, edges = recovered_edges[ keep, , drop = FALSE ]))
}


# Start from the treatment and outcome ancestors. Each extra generation adds
# their parents and children. Pinned nodes remain in the result.
#' Extend endpoint ancestry by additional parent-child generations
#'
#' @param recovered_vertices Unique vertex names in the recovered graph.
#' @param recovered_edges Directed `from`/`to` table on those vertices.
#' @param treat_v,out_v Treatment and outcome occurrence names.
#' @param pinned Vertex names that must remain included.
#' @param graph_order Nonnegative extra-generation count, or `Inf` to closure.
#' @param quiet Suppress the retained-node summary when `TRUE`.
#' @return A list of retained `vertices` and induced `edges` after ancestry,
#'   requested expansion and pinned-node inclusion.
#' @noRd
.scope_to_query <- function( recovered_vertices, recovered_edges, treat_v, out_v, pinned, graph_order, quiet ){
  parents_by_child <- split(as.character(recovered_edges$from), as.character(recovered_edges$to))
  children_by_parent <- split(as.character(recovered_edges$to),   as.character(recovered_edges$from))


  current <- unique(c(treat_v, out_v))
  repeat {
    new_ancestors <- setdiff(unique(unlist(parents_by_child[ current ], use.names = FALSE)), current)
    new_ancestors <- intersect(new_ancestors, recovered_vertices)
    if( length(new_ancestors) == 0L ){
      break
    }
    current <- c(current, new_ancestors)
  }

  # Expand: each generation adds parents AND children of current boundary
  generations <- 0L
  if( graph_order > 0L ){
    limit <- if( is.infinite(graph_order) ){ length(recovered_vertices) }else{ as.integer(graph_order) }
    for( i in seq_len(limit) ){
      new_parents <- setdiff(unique(unlist(parents_by_child[ current ], use.names = FALSE)), current)
      new_children <- setdiff(unique(unlist(children_by_parent[ current ], use.names = FALSE)), current)
      new_nodes <- intersect(unique(c(new_parents, new_children)), recovered_vertices)
      if( length(new_nodes) == 0L ){
        break
      }
      current <- c(current, new_nodes)
      generations <- generations + 1L
    }
  }


  current <- unique(c(current, intersect(pinned, recovered_vertices)))


  keep <- recovered_edges$from %in% current & recovered_edges$to %in% current

  if( isFALSE(quiet) ){
    message(length(current), " of ", length(recovered_vertices), " nodes retained",
            if( generations > 0L ){ paste0(" (", generations, " generation(s) expanded)") }else{ "" },
            ".")
  }

  return(list(vertices = current, edges = recovered_edges[ keep, ]))
}


# Reviewed public output paths, qualified by crop. Empty/missing metadata never
# certifies availability. This register is independent of recovery policy and
# the older apsim_obtainable module locators. No node layer changes are implied.
#' Read verified public APSIM outputs for one crop
#'
#' @param crop One crop name in the bundled output-availability register.
#' @return A data frame of register rows for the crop with status
#'   `source_verified_public` and a nonempty source expression. If none match,
#'   the result has no rows and supplies no evidence of public output
#'   availability. A missing register raises an error.
#' @noRd
.query_apsim_outputs <- function(crop){
  path <- system.file("extdata", "apsim-output-availability.csv", package = "causaliflower")
  if( !nzchar(path) ){ stop("The APSIM output availability register is missing. Reload the updated package.") }
  mappings <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  return(mappings[mappings$crop == crop & !is.na(mappings$source_expression) &
             nzchar(trimws(mappings$source_expression)) &
             mappings$status == "source_verified_public", , drop = FALSE])
}


#' Resolve recovery rules and record conflicting output locators
#'
#' @param nodes Prepared node table with `node_id` and optional `recovered_node`,
#'   `provenance` and `apsim_obtainable` columns.
#' @param crop One crop name. Recovery rules in `nodes` apply to barley, while
#'   explicit non-APSIM provenance applies across crops.
#' @return A data frame with `node_id`, effective `recovered_node`, `reason`,
#'   `apsim_locator` and `metadata_conflict`. The source table is unchanged.
#' @noRd
.query_output_policy <- function(nodes, crop){
  ids <- nodes$node_id
  rule <- rep("yes", length(ids))
  if( identical(crop, "barley") && "recovered_node" %in% names(nodes) ){
    rule <- as.character(nodes$recovered_node)
  }
  if( anyNA(rule) || any(!rule %in% c("yes", "yes:user") & !grepl("^no:.+", rule)) ){
    stop("Invalid recovered_node policy in the node table.")
  }
  non_apsim <- if( "provenance" %in% names(nodes) ){ nodes$provenance %in% "not_apsim" }else{ rep(FALSE, length(ids)) }
  rule[non_apsim & !startsWith(rule, "no:")] <- "yes:user"
  locator <- if( "apsim_obtainable" %in% names(nodes) ){ as.character(nodes$apsim_obtainable) }else{ rep("", length(ids)) }
  return(data.frame(node_id = ids, recovered_node = rule,
    reason = ifelse(non_apsim, "Provenance outside APSIM. The user must supply its interpretation/data",
      ifelse(rule == "yes:user", "Explicit quantity definition/request required",
        ifelse(startsWith(rule, "no:"), substring(rule, 4L), "Presumed obtainable. Not individually certified"))),
    apsim_locator = locator,
    metadata_conflict = rule != "yes" & !is.na(locator) & nzchar(trimws(locator)),
    stringsAsFactors = FALSE))
}


#' Assign observation status after query inclusion is resolved
#'
#' @param recovered_vertices Retained vertex names, possibly with time suffixes.
#' @param hidden_vertex Nonmissing logical hidden-node flags by vertex name.
#' @param has_columns `TRUE` if explicit data column names were supplied.
#' @param matched_ids Quantity IDs observed in those columns.
#' @param show_hidden_nodes `TRUE` if hidden quantities remain in the graph.
#' @param always_observed Quantity IDs exempt from latent tagging.
#' @param apsim_observed_ids Quantity IDs presumed obtainable without columns.
#' @param unavailable_ids Quantity IDs marked latent without columns when
#'   hidden quantities are retained.
#' @return A character vector of vertex names to tag latent, empty if none
#'   qualify. Explicit columns override automatic availability rules.
#'   `always_observed` exemptions are applied last.
#' @noRd
.assign_latent <- function( recovered_vertices, hidden_vertex, has_columns, matched_ids, show_hidden_nodes,
                            always_observed = character(0L), apsim_observed_ids = character(0L),
                            unavailable_ids = character(0L) ){
  if( isFALSE(has_columns) ){
    if( isFALSE(show_hidden_nodes) ){
      latent <- character(0L)
    }else{
      latent <- recovered_vertices[ hidden_vertex[ recovered_vertices ] ]
      latent <- latent[!(.strip_slice(latent) %in% apsim_observed_ids)]
      latent <- union(latent, recovered_vertices[.strip_slice(recovered_vertices) %in% unavailable_ids])
    }
  }else{
    base_ids <- .strip_slice(recovered_vertices)
    latent <- recovered_vertices[ !(base_ids %in% matched_ids) ]
  }
  if( length(always_observed) > 0L ){
    latent <- setdiff(latent, recovered_vertices[ .strip_slice(recovered_vertices) %in% always_observed ])
  }
  return(latent)
}


#' Emit a directed graph with explicit vertex and role declarations
#'
#' @param vertices Unique vertex names suitable for quoted dagitty identifiers.
#' @param edges Directed `from`/`to` table on those vertices.
#' @param latent_set,exposures,outcomes Vertex names carrying the corresponding roles.
#' @return One dagitty graph string. All vertices are declared, including
#'   isolated ones. Role names absent from the graph are ignored.
#' @noRd
.emit_dagitty <- function( vertices, edges, latent_set, exposures = character(0L), outcomes = character(0L) ){
  quote_v <- function( v ){ # Quote a vertex identifier for dagitty syntax
    return(paste0("\"", v, "\""))
  }
  tags <- stats::setNames(vector("list", length(vertices)), vertices)
  for( v in intersect(latent_set, vertices) ){
    tags[[ v ]] <- c(tags[[ v ]], "latent")
  }
  for( v in intersect(exposures, vertices) ){
    tags[[ v ]] <- c(tags[[ v ]], "exposure")
  }
  for( v in intersect(outcomes, vertices) ){
    tags[[ v ]] <- c(tags[[ v ]], "outcome")
  }
  node_lines <- sapply(vertices, function( v ){ # Attach the assigned causal roles to one declared vertex
    vertex_tags <- tags[[ v ]]
    if( is.null(vertex_tags) ){
      return(quote_v(v))
    }
    return(paste0(quote_v(v), " [", paste(vertex_tags, collapse = ","), "]"))
  })
  edge_lines <- character(0L)
  if( nrow(edges) > 0L ){
    edge_lines <- paste0(quote_v(edges$from), " -> ", quote_v(edges$to))
  }
  body <- paste(c(node_lines, edge_lines), collapse = "\n")
  return(paste0("dag {\n", body, "\n}\n"))
}


#' Remove quotes and line breaks from a requested display label
#'
#' @param x Character label supplied as an input-name override.
#' @return Trimmed text without double quotes, carriage returns or newlines.
#'   An empty result tells the caller to retain its previous label.
#' @noRd
.sanitise_label <- function( x ){
  return(trimws(gsub("[\"\r\n]", "", x)))
}


#' Record canonical quantity membership before display renaming
#'
#' @param vertices Current occurrence names, with optional time suffixes.
#' @param group_map Grouped quantity IDs indexed by source quantity IDs.
#' @param kept_ids Source quantity IDs admitted to the query.
#' @return A data frame with `node`, `quantity`, `member_id`, `occurrence` and
#'   `member_occurrence`, one row per vertex/member pair. An empty graph returns
#'   the same five character columns. This ledger does not alter the graph.
#' @noRd
.query_identity_rows <- function(vertices, group_map, kept_ids){
  grouped_members <- split(kept_ids, unname(group_map[ kept_ids ]))
  rows <- lapply(vertices, function(vertex){ # Describe the canonical members of one emitted vertex
    quantity <- .strip_slice(vertex)
    members <- grouped_members[[ quantity ]]
    if( length(members) == 0L ){ members <- quantity }
    suffix <- substring(vertex, nchar(quantity) + 1L)
    return(data.frame(node = vertex, quantity = quantity, member_id = members,
               occurrence = vertex, member_occurrence = paste0(members, suffix),
               stringsAsFactors = FALSE))
  })
  if( length(rows) == 0L ){
    return(data.frame(node = character(), quantity = character(), member_id = character(),
                      occurrence = character(), member_occurrence = character()))
  }
  return(do.call(rbind, rows))
}


# Apply input aliases as node labels to a recovered graph (the label= feature). For each distinct base
# (the suffix stripped vertex, which equals a node_id) a base label is chosen by mode:
#   "input"          -> node_id                            (an input override may rename any node)
#   "APSIM"          -> node_id                            (an input override renames only
#                                                           undotted keys)
#   "readable" -> node_name_readable, falling back to node_id where the human label is empty
# The legacy APSIM label mode uses has_dot to protect dotted keys from input overrides.
# An input override from input_pairs (node_id -> input string) replaces the base label wherever
# override_allowed = !( label == "APSIM" && has_dot ), after sanitising (guard b). An empty
# sanitised label is ignored. A label claimed by more than one base reverts those bases to the
# node_id until both base labels and suffixed vertex names are unique. The temporal
# suffix is then reattached and every vertex / edge / latent / exposure / outcome remapped. Returns
# the relabelled sets plus input_labels (displayed label -> node_id) for the nodes an input renamed.
#' Apply query display labels while preserving vertex identities
#'
#' @param V Unique canonical vertex names, with optional temporal suffixes.
#' @param E Directed `from`/`to` data.table on `V`.
#' @param lat,X,Y Latent, treatment and outcome vertex names respectively.
#' @param nodes Node table with `node_id` and `node_name_readable`.
#' @param input_pairs Supplied input strings named by their resolved quantity IDs.
#' @param label One of `input`, `APSIM` or `readable`.
#' @return A list with renamed `V`, `E`, `lat`, `X` and `Y`, plus `input_labels`
#'   mapping retained input overrides back to their canonical quantity IDs.
#'   Conflicting labels fall back to IDs until final vertex names are unique.
#' @noRd
.apply_input_labels <- function( V, E, lat, X, Y, nodes, input_pairs, label ){
  .datatable.aware <- TRUE
  base <- .strip_slice(V)
  suffix <- substring(V, nchar(base) + 1L)
  unique_bases <- unique(base)
  has_dot <- grepl(".", unique_bases, fixed = TRUE)

  # label for each distinct base by mode
  if( identical(label, "readable") ){
    human_labels <- stats::setNames(as.character(nodes$node_name_readable), nodes$node_id)
    candidate_labels <- unname(human_labels[ unique_bases ])
    missing_labels <- is.na(candidate_labels) | !nzchar(candidate_labels)
    candidate_labels[ missing_labels ] <- unique_bases[ missing_labels ]
  }else{
    candidate_labels <- unique_bases                                             # "APSIM" and "input": base = node_id
  }
  candidate_labels <- stats::setNames(candidate_labels, unique_bases)

  # input override, per base, where allowed. Record which bases an input actually renamed
  override_allowed <- stats::setNames(if( identical(label, "APSIM") ){ !has_dot }else{ rep(TRUE, length(unique_bases)) }, unique_bases)
  renamed <- character(0L)                                 # node_id -> displayed label
  for( node_id in intersect(names(input_pairs), unique_bases) ){
    if( isTRUE(override_allowed[[ node_id ]]) ){
      input_label <- .sanitise_label(input_pairs[[ node_id ]])
      if( nzchar(input_label) ){
        candidate_labels[[ node_id ]] <- input_label
        renamed[[ node_id ]] <- input_label
      }
    }
  }

  candidate_labels <- .query_unique_labels(candidate_labels, base, suffix)
  # a base that the guard reverted no longer holds its input label
  if( length(renamed) > 0L ){
    still <- vapply(names(renamed), function( node_id ){ # Check whether the uniqueness guard retained an input label
      return(identical(unname(candidate_labels[[ node_id ]]), unname(renamed[[ node_id ]])))
    }, logical(1L))
    renamed <- renamed[ still ]
  }

  relabel <- stats::setNames(paste0(unname(candidate_labels[ base ]), suffix), V)
  new_vertices <- unname(relabel[ V ])
  new_edges <- if( nrow(E) > 0L ){
    data.table::data.table(from = unname(relabel[ E$from ]), to = unname(relabel[ E$to ]))
  }else{
    E
  }
  new_latent <- unname(relabel[ intersect(lat, V) ])
  new_exposures <- unname(relabel[ intersect(X, V) ])
  new_outcomes <- unname(relabel[ intersect(Y, V) ])

  input_labels <- if( length(renamed) > 0L ){ stats::setNames(names(renamed), unname(renamed)) }else{ character(0L) }
  return(list(V = new_vertices, E = new_edges, lat = new_latent, X = new_exposures, Y = new_outcomes, input_labels = input_labels))
}


#' Keep display labels injective before and after restoring temporal suffixes
#'
#' @param candidate_labels Nonempty display labels named by unique quantity IDs.
#' @param base Quantity ID for each unique original vertex, indexing the labels.
#' @param suffix Temporal suffix removed from each original vertex. Reattaching
#'   these suffixes to `base` must recover unique original vertex names.
#' @return The named labels, with colliding quantities reverted to their IDs.
#'   Each iteration reverts at least one additional label, so the process ends
#'   after at most one reversion per quantity. Unaffected labels are retained.
#' @noRd
.query_unique_labels <- function(candidate_labels, base, suffix){
  repeat{
    duplicate_labels <- duplicated(candidate_labels) | duplicated(candidate_labels, fromLast = TRUE)
    vertex_labels <- paste0(unname(candidate_labels[base]), suffix)
    duplicate_vertices <- duplicated(vertex_labels) | duplicated(vertex_labels, fromLast = TRUE)
    conflicting_bases <- union(names(candidate_labels)[duplicate_labels], base[duplicate_vertices])
    if( length(conflicting_bases) == 0L ){
      return(candidate_labels)
    }
    candidate_labels[conflicting_bases] <- conflicting_bases
  }
}
