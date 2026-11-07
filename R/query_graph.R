################################################################################
# query_graph.R


# Public graph queries and their temporal projection algorithms
# Original filename section markers below record source history


#####
# Query and source graph entry points
# Historical source: dag_by_query.R
#####
#' Build a graph from a crop query
#'
#' Constructs a causal graph for one crop from the bundled structural causal model
#' derived from APSIM. The query selects the crop, management inputs and a finite
#' temporal window. Defaults retain eligible hidden quantities and separate time
#' occasions, without recombination or an assumption of randomisation.
#'
#' @param crop Crop name. The bundled data support barley, wheat and maize.
#' @param columns Optional data column names. These resolve management inputs and
#'   identify observed quantities. Unmatched quantities are marked latent,
#'   subject to the endpoint and recovery rules below.
#' @param treatments,outcomes Quantity names or aliases. Requested endpoints are
#'   retained, are not recombined, and receive exposure or outcome tags.
#' @param randomised_treatment Defaults to `FALSE`. With `TRUE`, incoming directed
#'   arrows are removed from every treatment occasion. Use this only when the
#'   treatment was randomised. Otherwise it can remove genuine confounding paths
#'   and invalidate the adjustment sets. Removed arrows and the assumption are saved.
#' @param time_unroll Representation of lagged relationships. When omitted, uses
#'   `"all"` if a treatment and lagged edges are present, otherwise `"none"`.
#'   Explicit `NULL` means `"none"`.
#' @param label `NULL` retains readable quantity keys. `"input"` uses a resolved
#'   input alias as the label. `"readable"` uses the display label, falling
#'   back to the quantity key when unavailable, with the same input override.
#'   The legacy `"APSIM"` mode retains dotted keys where
#'   available and otherwise follows `"input"`. The label to quantity map is saved.
#' @param quiet Suppress the build summary when `TRUE`.
#' @param params Further options supplied by `query_params()`.
#'
#' @details
#' Temporal modes:
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
#' A hidden quantity is mechanistic detail, not necessarily unobserved.
#' `show_hidden_nodes = FALSE` removes eligible hidden quantities while retaining
#' directed paths through them. This operation does not introduce bidirected
#' arrows. Explicitly requested quantities are retained. Recovery policy
#' exclusions are removed through the same directed operation regardless of this option.
#'
#' Without `columns`, the default `observe_apsim_outputs = TRUE` treats ordinary
#' included quantities as obtainable. This is a working assumption, not verified
#' availability for every APSIM output. With `FALSE`, included hidden quantities
#' are marked latent. Supplied columns take precedence over these automatic
#' observation tags. Query endpoints retain their observed status.
#'
#' The recovery policy distinguishes `yes`, quantities admitted only on explicit
#' request (`yes:user`), and exclusions (`no:reason`). Explicit endpoints or
#' columns can admit `yes:user`. Requesting an excluded quantity raises an error.
#' The current interpretation policy was reviewed for barley. Recorded provenance
#' outside APSIM applies across crops. Observation assumptions, conflicts, output
#' mappings and removed graph structure are saved with the query.
#'
#' Quantities left without any arrows are dropped unless they are requested
#' treatments or outcomes. Other occasions of an outcome in a growth stage query
#' remain distinct causal states. They are not automatically outcomes,
#' interventions or observed measurements.
#'
#' The bundled data support only SoilWater and individual readable quantity keys.
#' SWIM, coarser dotted grouping and supernode overrides raise errors. Supplied
#' columns and endpoint names must resolve to the bundled keys or aliases.
#'
#' @return A dagitty object with its query record, resolved quantities, edges,
#' roles, effective settings and source ordering metadata. The record describes
#' the supplied query, including any deliberately restricted scope. Construction
#' does not silently replace it with a fuller query. Use the returned `edges`
#' for the current graph. Saved source structures may describe the graph before
#' excluded quantities were removed.
#' @examples
#' \dontrun{
#' dag <- query_graph(crop = "wheat")
#' dag <- query_graph(crop = "barley", columns = c("plant_population", "soil_no3", "grain_yield"))
#' # broad default: full configured graph in an uncompressed finite time frame
#' dag <- query_graph(crop = "wheat", treatments = "N_fertiliser_applied",
#'                    outcomes = "grain_wt")
#' # complete ancestry relevant to the query
#' dag <- query_graph(crop = "wheat", treatments = "N_fertiliser_applied",
#'                    outcomes = "grain_wt", params = query_params(query_depth = 2L))
#' # previous compact defaults, explicitly requested (validate sets on the full graph)
#' dag <- query_graph(crop = "wheat", treatments = "N_fertiliser_applied",
#'                    outcomes = "grain_wt", time_unroll = "pre_post",
#'                    params = query_params(query_depth = 1L, show_hidden_nodes = FALSE))
#' # experimental trial: randomised treatment has no parents (do() graph surgery)
#' dag <- query_graph(crop = "wheat", treatments = "N_fertiliser_applied",
#'                    outcomes = "grain_wt", randomised_treatment = TRUE)
#' # developmental phase anchoring via query_params()
#' dag <- query_graph(crop = "barley", treatments = "N_fertiliser_applied",
#'                    outcomes = "grain_protein",
#'                    params = query_params(treatment_growth_stage = "flowering",
#'                                          outcome_growth_stage   = "harvest_ripe"))
#' }
#' @export
query_graph <- function( crop,
                                columns              = NULL,
                                treatments           = NULL,
                                outcomes             = NULL,
                                randomised_treatment = FALSE,
                                time_unroll          = NULL,
                                label                = NULL,
                                quiet                = FALSE,
                                params               = query_params()
                                ){ # Recover the configured graph and record the requested causal query

  .datatable.aware <- TRUE
  query_call <- match.call(expand.dots = FALSE)
  query_input <- list(crop = crop, columns = columns, treatments = treatments,
                      outcomes = outcomes, randomised_treatment = randomised_treatment,
                      time_unroll = time_unroll, label = label, quiet = quiet, params = params)
  ## Whether time_unroll was left at its default (vs supplied, even as NULL): the default resolves to
  ## "all" (uncompressed occasions) once the kept subgraph is known, whereas an
  ## explicit NULL means "none". Captured before time_unroll is touched.
  time_unroll_defaulted <- missing(time_unroll)
  ## Bind package data explicitly so manual and wrapper calls cannot read a
  ## table with the same name from the caller's workspace.
  nodes <- getExportedValue("causaliflower", "apsim_dag_nodes")
  source_order <- .query_source_order(nodes)
  aliases <- getExportedValue("causaliflower", "node_aliases")
  management_defaults <- getExportedValue("causaliflower", "management_defaults")
  crop_phases <- getExportedValue("causaliflower", "crop_phases")
  phase_of <- getExportedValue("causaliflower", "phase_of")
  ## Read query_params() options into the local names used below
  ## (water_model / apsim_nodes_detail / phase_window are already validated by match.arg inside query_params.)
  query_depth <- params$query_depth
  water_model <- params$water_model
  apsim_nodes_detail <- params$apsim_nodes_detail
  supernode_overrides <- params$supernode_overrides
  force_acyclic <- params$force_acyclic
  show_hidden_nodes <- params$show_hidden_nodes
  observe_apsim_outputs <- if (is.null(params$observe_apsim_outputs)) { TRUE } else { params$observe_apsim_outputs }
  node_recombination <- params$node_recombination
  disambiguate <- params$disambiguate
  treatment_growth_stage <- params$treatment_growth_stage
  outcome_growth_stage <- params$outcome_growth_stage
  phase_window <- params$phase_window
  ## Guard: query_params() passed positionally lands on `columns`. Catch this common mistake early.
  if ( is.list(columns) && "query_depth" %in% names(columns) ){
    stop("It looks like query_params() was passed as a positional argument and matched to 'columns'. Use the named argument instead: params = query_params(...). Please fix the call and try again.")
  }
  if( !is.null(source_order) ){
    if( identical(water_model, "SWIM") ){
      stop("water_model = \"SWIM\" was requested, but the bundled slice has no SWIM family rows. Use water_model = \"SoilWater\".", call. = FALSE)
    }
    if( !identical(apsim_nodes_detail, "variable") || length(supernode_overrides) > 0L ){
      stop("The bundled slice has readable variable keys without dotted groups. Use apsim_nodes_detail = \"variable\" without supernode_overrides.", call. = FALSE)
    }
  }
  ## The default NULL preserves readable node keys. Input labels and readable labels
  ## remain explicit choices. The legacy APSIM mode preserves dotted identifiers when present.
  if ( !is.null(label) ){
    label <- match.arg(label, c("APSIM", "input", "readable"))
  }
  crops <- sort(unique(unlist(strsplit(nodes$crop_type[ nodes$crop_type != "agnostic" ], ";", fixed = TRUE))))
  if( !crop %in% crops ){
    stop(paste0("crop must be one of: ", paste(crops, collapse = ", "), ". Please supply a supported crop and try again."))
  }

  ## Hidden node ids, used as an ambiguity rule for breaking ties in name resolution (a hidden node is never a
  ## uniquely observed data column). Computed here so it is available at column resolution time.
  hidden_ids <- nodes$node_id[ nodes$hidden_node == "TRUE" ]

  ## Resolution for the selected water_model: the two soil water models are mutually exclusive (.filter_nodes
  ## drops the unused family), so a supplied name must never resolve to a node from an inactive model.
  ## Restrict the alias vocabulary to the active model before building the resolution index. The
  ## full `aliases` object is still used elsewhere (e.g. get_management_options).
  inactive_prefix <- if( identical(water_model, "SoilWater") ){ "SWIM." }else{ "SoilWater." }
  aliases_active <- aliases[ !startsWith(names(aliases), inactive_prefix) ]
  alias_index <- .build_alias_index(aliases_active)
  ## time_unroll controls how lagged (@t-1) relationships are represented. The DEFAULT (argument unset)
  ## resolves to "all" when a treatment anchors the split and the kept subgraph carries lagged edges
  ## (otherwise "none"): an uncompressed finite time frame. Passing time_unroll = NULL
  ## EXPLICITLY means "none" (do not unroll). The legacy logical (FALSE -> "none", TRUE -> "all") is still
  ## accepted, as are the named modes "none", "all", "selective", "collapse", "pre_post":
  ##   * "pre_post" splits each straddler into <id>_pre_treatment / <id>_post_treatment (residual cycle
  ##      bases fall back to <id>_t{s}). Acyclic. The historical compact view.
  ##   * "selective" splits straddlers at every slice they occupy (<id>_t{s}). Acyclic, ready for adjustment.
  ##   * "all" is the uncompressed frame (every included base split at every slice). Acyclic.
  ##   * "collapse" runs "all" then projects every occasion back onto its static base. May be cyclic.
  ##   * "none" does not unroll. Every lagged edge collapses onto one slice. Generally cyclic.
  ## The defaulted case is resolved to a concrete mode below, once g_lagged is known. Here we only
  ## normalise the explicit forms and validate them. The NULL sentinel marks "still defaulted".
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
    ## show_hidden_nodes = TRUE is supported: the slicer skips its internal .collapse_hidden and
    ## hidden occasions are carried through to the emitted graph, marked [latent] downstream.
  }

  ## Phase anchored (developmental time) mode: triggered whenever treatment_growth_stage or
  ## outcome_growth_stage is supplied (both now via query_params()). It runs the absolute phase selective
  ## frame (.selective_phase_slice, with B0 at the sowing origin) instead of the frame anchored
  ## to treatment. When both growth stage params are NULL, query_graph uses time_unroll.
  ## query_params() validates phase_window, which is used only here.
  phase_mode <- !is.null(treatment_growth_stage) || !is.null(outcome_growth_stage)
  if ( phase_mode ){
    if ( length(treatments) == 0L ){
      stop("A query anchored to developmental phases needs at least one treatment: treatment_growth_stage / outcome_growth_stage anchor the treatment and outcome to developmental phases. Please pass treatments = and try again.")
    }
    ## show_hidden_nodes = TRUE is supported: the phase slicer skips its internal .collapse_hidden
    ## and hidden occasions are carried through to the emitted graph, marked [latent] downstream.
    if ( !time_unroll_defaulted && !identical(time_unroll, "selective") ){
      stop("A query anchored to developmental phases uses the selective frame in developmental time. Leave time_unroll unset (the default) or set it to \"selective\", not \"none\" / \"all\" / \"collapse\" / \"pre_post\". Please adjust time_unroll and try again.")
    }
    ## the phase frame is selective in developmental time: resolve the sentinel so downstream
    ## references and the emitted attribute read "selective".
    time_unroll <- "selective"
    if ( !exists("crop_phases") ){
      stop("A query anchored to developmental phases needs the crop_phases table, which is not loaded. Build the node layer (its mobilise step writes crop_phases.rda) and load it, then try again.")
    }
    if ( !(crop %in% crop_phases$crop) ){
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
                     ". Rename them to a node_id or pass disambiguate = to resolve them."))
    }
    if( !is.null(source_order) && length(col_res$unmatched) > 0L ){
      stop("These columns match no node in the bundled slice: ",
           paste(col_res$unmatched, collapse = ", "),
           ". Use a readable node key or an unambiguous alias.", call. = FALSE)
    }
  }

  ## treatment / outcome node ids (resolved by alias or supplied directly as node_ids)
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

  ## input -> node_id pairings for the label= feature: a resolved input alias becomes the node's
  ## display label. Built here (skipped entirely under the NULL setting). The pairing draws on the
  ## SAME resolver as the structural resolution but through SEPARATE .resolve_columns() calls kept
  ## purely for the input -> node_id pairing, so the structural resolution above (col_matched,
  ## treat_ids, out_ids) is untouched. Precedence (proposal handover 2.4): treatments / outcomes beat
  ## columns. First occurrence wins within a source (two inputs -> one node keeps the first input).
  ## input_pairs is keyed by node_id at the underlying (subnode) level. When apsim_nodes_detail is
  ## coarser the resolved node_id may be folded into a group label the base no longer equals, so that
  ## node keeps its group label - labelling is most meaningful at the "variable" default.
  input_pairs <- character(0L)                              # node_id -> input string
  if ( !is.null(label) ){
    .add_first <- function( pairs, resolved ){ # Keep the first resolved input label for each node
      ## resolved: NAMED input -> node_id. Accumulate node_id -> input, first writer wins.
      for ( inp in names(resolved) ){
        nid <- unname(resolved[[ inp ]])
        if ( !( nid %in% names(pairs) ) ){
          pairs[[ nid ]] <- inp
        }
      }
      return(pairs)
    }
    col_pairs <- if ( has_columns ){ .add_first(character(0L), col_res$resolved) }else{ character(0L) }
    to_pairs <- .add_first(character(0L), .resolve_columns(treatments, alias_index, disambiguate, hidden_ids)$resolved)
    to_pairs <- .add_first(to_pairs,      .resolve_columns(outcomes,   alias_index, disambiguate, hidden_ids)$resolved)
    input_pairs <- col_pairs
    for ( nid in names(to_pairs) ){                          # treatments / outcomes override columns
      input_pairs[[ nid ]] <- to_pairs[[ nid ]]
    }
    ## label = "input" shows the supplied labels: the names in columns = c(display = token),
    ## or the token itself when the element is unnamed, and likewise for treatments / outcomes.
    ## These labels need not be registered aliases. Keep a custom name even when it cannot be
    ## resolved back through the alias map. .sanitise_label
    ## strips any character that would corrupt the emitted dagitty text, and the uniqueness guard in
    ## .apply_input_labels reverts any label collision back to the node_id.
  }

  resolved <- .resolve_edges(nodes, crop = crop)
  if( !.is_acyclic(nodes$node_id, resolved$contemporaneous) ){
    stop("The contemporaneous structural graph must be acyclic. Please check function_of / composition_of in the node table and try again.")
  }

  present <- .resolve_management_presence(nodes, management_defaults$core, c(treat_ids, out_ids), col_matched, alias_index)
  kept_ids <- .filter_nodes(nodes, crop, water_model, present)

  ## APSIM detail grouping: "variable" (default) keeps complete ids. "subnode" keeps three
  ## dotted segments. "node" and "supernode" contract dotted ids to the first two / one
  ## segment(s) and rewire edges through the grouping with the same directed,
  ## structure preserving logic as the hidden collapse. supernode_overrides (a named
  ## supernode -> level vector) sets a finer level for specific first segment groups on top
  ## of the global level, as set by force_acyclic.
  group_map <- .apsim_truncate(nodes$node_id, apsim_nodes_detail, overrides = supernode_overrides)
  g_kept <- unique(unname(group_map[ kept_ids ]))
  g_contemp <- .contract_edges(resolved$contemporaneous[ from %in% kept_ids & to %in% kept_ids ], group_map)
  g_lagged <- .contract_edges(resolved$lagged[ from %in% kept_ids & to %in% kept_ids ], group_map)

  ## Resolve treatments / outcomes to the underlying node_id SET (stable across any later
  ## refinement), so a supernode named as treatment / outcome maps to whatever group label(s)
  ## it ends up as. treat_ids / out_ids are the direct node_id matches. A supernode name
  ## additionally pulls in every node_id its current label covers (strict: only resolves at
  ## the detail level where that group label exists).
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
  if (length(forbidden_requests)) {
    reasons <- output_policy$recovered_node[match(forbidden_requests, output_policy$node_id)]
    stop(paste0("Cannot recover requested node(s): ",
                paste(paste0(forbidden_requests, " (", reasons, ")"), collapse = ", ")), call. = FALSE)
  }
  excluded_ids <- union(prohibited, setdiff(
    output_policy$node_id[output_policy$recovered_node == "yes:user"], requested_ids))

  ## A bare supernode name supplied at a finer than supernode level resolves to nothing
  ## (strict). Flag that case so the user is not left wondering why.
  all_supers <- unique(first_segment(nodes$node_id))
  unmatched_super <- setdiff(intersect(c(treatments, outcomes), all_supers),
                             unique(c(treat_named, out_named, unname(group_map[ c(treat_ids, out_ids) ]))))
  if( length(unmatched_super) > 0L ){
    message(paste0("Treatment / outcome name(s) ", paste(unmatched_super, collapse = ", "),
                   " are supernodes, but apsim_nodes_detail = \"", apsim_nodes_detail,
                   "\" is finer than \"supernode\" so they are not single nodes. Name a node id, ",
                   "set apsim_nodes_detail = \"supernode\", or pass supernode_overrides to recover them."))
  }

  ## supernodes we prefer not to split when forcing acyclicity: those named as treatment /
  ## outcome, and the present management levers (all keyed by first segment).
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
                     "to reduce automatically to the coarsest acyclic representation, or select a finer apsim_nodes_detail (\"node\" or \"subnode\")."))
    }
  }

  # A coarse label must not combine an excluded interpretation with allowed
  # quantities (for example spike and Spike.*). Refine only those mixed groups.
  group_members <- split(kept_ids, unname(group_map[kept_ids]))
  mixed <- vapply(group_members, function(ids) { # Find groups containing both hidden and visible nodes
    return(any(ids %in% excluded_ids) && any(!ids %in% excluded_ids))
  }, logical(1L))
  refined_ids <- unlist(group_members[mixed], use.names = FALSE)
  if (length(refined_ids)) {
    group_map[refined_ids] <- refined_ids
    g_kept <- unique(unname(group_map[kept_ids]))
    g_contemp <- .contract_edges(resolved$contemporaneous[from %in% kept_ids & to %in% kept_ids], group_map)
    g_lagged <- .contract_edges(resolved$lagged[from %in% kept_ids & to %in% kept_ids], group_map)
  }
  g_excluded <- unique(unname(group_map[intersect(excluded_ids, kept_ids)]))
  g_requested <- unique(unname(group_map[intersect(requested_ids, kept_ids)]))

  ## grouped names of the treatment / outcome nodes (after refinement), used to pin and tag
  ## them below. Mapping the stable set of ids through the final group_map means a named
  ## supernode that force_acyclic had to split is recovered as ALL of its split pieces, each of
  ## which is then pinned and tagged - the multiple treatment / multiple outcome representation
  ## required once an aggregate cannot be a single acyclic node.
  treat_grouped <- unique(unname(group_map[ intersect(treat_id_set, names(group_map)) ]))
  out_grouped <- unique(unname(group_map[ intersect(out_id_set,   names(group_map)) ]))

  ## a grouped node is hidden only when every node it groups is hidden (hidden_ids computed
  ## above). It counts as matched (has data) when any node it groups matched a supplied column.
  g_hidden <- names(which(tapply(nodes$node_id %in% hidden_ids, unname(group_map[ nodes$node_id ]), all)))
  g_matched <- unique(unname(group_map[ col_matched ]))
  output_mappings <- .query_apsim_outputs(crop)
  output_ids <- if (isTRUE(observe_apsim_outputs)) {
    output_policy$node_id[output_policy$recovered_node == "yes"]
  } else { character(0L) }
  # A grouped node is eligible only if every represented member is admitted.
  g_outputs <- names(which(tapply(nodes$node_id %in% output_ids,
                                 unname(group_map[nodes$node_id]), all)))

  ## resolve the defaulted sentinel now that g_lagged is available: if there are lagged edges in the
  ## kept subgraph and treatments are specified, the default is "all" -- an uncompressed, acyclic,
  ## frame anchored to treatment within the adaptive finite window. Otherwise "none".
  ## (An explicit NULL was already resolved to "none" above. Only the unset default reaches here as NULL.)
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

  ## Path dependent stocks (persistent_ids) carry a temporal succession edge in selective mode
  ## (the earlier occasion causes the later one). Read from a path_dependent column if the node
  ## table carries one. Absent that column, no succession edges are emitted (graceful fallback).
  persistent_ids <- if ( "path_dependent" %in% names(nodes) ){
    pd_nodes <- nodes$node_id[ toupper(as.character(nodes$path_dependent)) == "TRUE" ]
    unique(c(pd_nodes, unname(group_map[ intersect(pd_nodes, names(group_map)) ])))
  }else{
    character(0L)
  }

  if( phase_mode ){
    ## developmental phase selective frame. Resolve each treatment / outcome to its anchored phase
    ## index(es). Every treatment must be anchored, an unanchored outcome defaults to the terminal
    ## (measurement) phase. B0 flows through via the contracted init edges.
    g_initial <- .contract_edges(resolved$initial[ from %in% kept_ids & to %in% kept_ids ], group_map)
    n_ph <- max(crop_phases$phase_index[ crop_phases$crop == crop ]) + 1L
    ts_map <- .stages_to_phase_indices(.normalize_stage_spec(treatment_growth_stage, treat_grouped), crop, crop_phases)
    os_map <- .stages_to_phase_indices(.normalize_stage_spec(outcome_growth_stage,   out_grouped),   crop, crop_phases)
    unanchored_t <- setdiff(treat_grouped, names(ts_map))
    if ( length(unanchored_t) > 0L ){
      stop(sprintf("Every treatment must be given a growth stage in a query anchored to developmental phases. Missing: %s. Pass treatment_growth_stage as a single stage for all, or a named list per treatment, and try again.",
                   paste(unanchored_t, collapse = ", ")))
    }
    for ( o in setdiff(out_grouped, names(os_map)) ){ os_map[[ o ]] <- n_ph - 1L }
    stock_ids <- intersect(persistent_ids, g_kept)
    mgmt_ids <- if ( "temporal_class" %in% names(nodes) ){
      intersect(nodes$node_id[ nodes$temporal_class == "per_event" ], g_kept)
    } else {
      character(0L)
    }
    pa <- if ( exists("phase_of") && is.data.frame(phase_of) && nrow(phase_of) > 0L ) {
      .resolve_phase_active(phase_of, crop, crop_phases)
    } else {
      NULL
    }
    sliced <- .selective_phase_slice(g_kept, g_contemp, g_lagged, g_initial, g_hidden,
                                     ts_map, os_map, n_ph, stock_ids, mgmt_ids, persistent_ids, phase_window,
                                     phase_active = pa, show_hidden_nodes = TRUE)
  }else if( time_unroll %in% c("selective", "all", "collapse", "pre_post") ){
    ## Same slicer, same window, same scope. "all" is the UNCOMPRESSED frame: every included base is
    ## split at every slice rather than only the straddlers, so no lagged edge lacks an occasion to
    ## land on. "selective" compresses bases that do not straddle treatment back to a single bare node. "collapse" runs the
    ## SAME uncompressed frame as "all" (so every lagged edge has an occasion to land on) and defers the
    ## static projection to a final step in the query_graph tail. It is uncompressed here for the same
    ## reason "all" is. "pre_post" is the straddler split of "selective", but each straddler collapses to
    ## two occasions (<id>_pre_treatment / <id>_post_treatment by descendant status) instead of one per
    ## slice. label_mode = "pre_post" selects that inside the slicer, with the promote loop keeping it
    ## acyclic.
    sliced <- .selective_slice(g_kept, g_contemp, g_lagged, g_hidden, treat_grouped, out_grouped,
                               persistent_ids,
                               split_all  = time_unroll %in% c("all", "collapse"),
                               label_mode = if ( identical(time_unroll, "pre_post") ){ "pre_post" }else{ "slice" },
                               show_hidden_nodes = TRUE)
  }else{
    sliced <- .contemporaneous_slice(g_kept, g_contemp, g_lagged)
    ## "none" does not unroll, so lagged edges have nowhere to go but onto the single slice. Every one
    ## of them is kept -- and collapsing a decision rule feedback (fertiliser reads the soil N it
    ## changes) or a biogeochemical loop onto one slice necessarily produces a cycle. Say so loudly:
    ## the result is a faithful picture of the lag topology but is NOT a DAG ready for adjustment.
    if( isTRUE(sliced$n_collapsed > 0L) ){
      warning(sprintf(paste0("time_unroll = \"none\" does not unroll time, so %d lagged edge%s been ",
                             "COLLAPSED onto the single slice and %s emitted as ordinary edge%s. No edge ",
                             "was dropped. The result %s -- collapsing a feedback loop (e.g. the ",
                             "fertiliser decision rule reads the soil N it changes) onto one slice makes ",
                             "it a cycle. Use time_unroll = \"selective\" (or \"all\") for an acyclic, ",
                             "DAG ready for adjustment."),
                     sliced$n_collapsed,
                     if( sliced$n_collapsed == 1L ){ " has" }else{ "s have" },
                     if( sliced$n_collapsed == 1L ){ "is" }else{ "are" },
                     if( sliced$n_collapsed == 1L ){ "" }else{ "s" },
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
  phase_treatments <- if (phase_mode) { sliced$treatments } else { character(0L) }
  phase_outcomes <- if (phase_mode) { sliced$outcomes } else { character(0L) }
  requested_occurrences <- if (phase_mode) {
    union(c(phase_treatments, phase_outcomes),
          recovered_vertices[.strip_slice(recovered_vertices) %in% g_matched])
  } else { recovered_vertices[.strip_slice(recovered_vertices) %in% g_requested] }

  ## Pinned nodes: treatments, outcomes, and present management levers. These are never
  ## dropped by scoping or recombination. Management levers are additionally subject to the
  ## final removal of isolated nodes below (only treatments / outcomes are exempt from that).
  pin_ids <- unique(c(treat_grouped, out_grouped, g_requested, unname(group_map[ present ])))
  pinned <- recovered_vertices[ .strip_slice(recovered_vertices) %in% pin_ids ]
  if (phase_mode) {
    pinned <- union(requested_occurrences,
                    recovered_vertices[.strip_slice(recovered_vertices) %in% unname(group_map[present])])
  }

  ## Temporal succession edges (selective mode): the earlier occasion of a path dependent stock
  ## causally precedes the later one (forward Euler identity: stock(t+1) = stock(t) + net_flux(t)).
  ## Add these before scoping so the ancestral scope retains the baseline occasion of a stock
  ## on the path X_pre -> X_post -> ... -> outcome. A stock outside that path whose later occasion
  ## reaches nothing (X_t1) is still pruned because it is not an ancestor, rather than kept as an orphan sink.
  ## Direction is always earlier -> later, so these cannot create a cycle.
  if ( (phase_mode || time_unroll %in% c("selective", "all", "collapse", "pre_post")) && !is.null(sliced$persistence) && nrow(sliced$persistence) > 0L ){
    p <- sliced$persistence
    hit <- p[ from %in% recovered_vertices & to %in% recovered_vertices ]
    if ( nrow(hit) > 0L ){
      recovered_edges <- unique(rbind(recovered_edges, hit))
    }
  }

  ## Scope the graph to the causal query subgraph when treatments and outcomes
  ## are given and query_depth is finite. Skipped when no query is specified.
  ## The unrestricted recovered graph is retained (full_vertices / full_edges) so get_valid_sets() can
  ## validate again candidate adjustment sets against the complete ancestral structure -- the "adjust"
  ## scope masks some backdoors, so a set that looks valid in the scoped graph may not be, and
  ## validation must use the unscoped graph.
  ## randomised_treatment = TRUE: the treatment is set by the experiment's randomisation, not by the
  ## modelled decision rule, so it has NO causal parents. This is do() graph surgery -- every
  ## directed edge INTO any treatment occasion is removed -- and it is ONLY sound when the treatment is
  ## randomised (as N rate / seed rate / variety are across the trial plots). It is off by
  ## default. If a treatment was not randomised, declaring it randomised removes its real
  ## confounders and yields adjustment sets that are too small.
  ##
  ## It runs here, before the full graph snapshot and before scoping, because do() graph surgery defines the
  ## graph identification is performed on: the scope (.scope_adjust / .scope_to_query) must select the
  ## relevant subgraph of the MUTILATED graph, and the full_vertices/full_edges snapshot (which
  ## get_valid_sets validates candidate sets against) must be the mutilated graph too, or a randomised
  ## query would be scoped and validated against a graph that still carries the cut confounder edges.
  ## Cutting the incoming edges also removes any decision rule feedback into the treatment (e.g. the fertiliser
  ## rule reading soil N), so the treatment itself is never caught in a cycle in the scoped/emitted graph.
  ## This removes only cycles that pass through the treatment. Feedback loops that bypass it (e.g.
  ## the soil water / mineralisation blob that time_unroll = "none" collapses onto one slice) are
  ## untouched -- those are broken only by unrolling time (selective / pre_post / all). Any parent left
  ## edgeless by the cut is removed by the orphan sweep below (a pure cause of the treatment drops out, while a
  ## confounder that also causes the outcome stays, now with only its outcome edge).
  randomised_cut <- character(0L)
  if( isTRUE(randomised_treatment) && length(treat_grouped) > 0L && nrow(recovered_edges) > 0L ){
    into_treat <- if (phase_mode) { recovered_edges$to %in% phase_treatments } else {
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
  hidden_occurrences <- if (isTRUE(show_hidden_nodes)) { character(0L) } else {
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
  ## query_depth starts at 1: 1 = the adjustment scope (.scope_adjust), 2 = the ancestral
  ## subgraph. Each further step adds a generation of parents and children (.scope_to_query with
  ## graph_order = query_depth - 2), Inf = the full graph (no scoping). Scoping needs a query.
  if ( is.finite(query_depth) && length(treat_grouped) > 0L && length(out_grouped) > 0L ){
    treat_v_pre <- if (phase_mode) { intersect(recovered_vertices, phase_treatments) } else {
      recovered_vertices[.strip_slice(recovered_vertices) %in% treat_grouped]
    }
    out_v_pre <- if (phase_mode) { intersect(recovered_vertices, phase_outcomes) } else {
      recovered_vertices[.strip_slice(recovered_vertices) %in% out_grouped]
    }
    scoped <- if ( query_depth == 1L ){
      .scope_adjust(recovered_vertices, recovered_edges, treat_v_pre, out_v_pre, pinned, quiet)
    } else {
      .scope_to_query(recovered_vertices, recovered_edges, treat_v_pre, out_v_pre,
                      pinned, query_depth - 2L, quiet)
    }
    recovered_vertices <- scoped$vertices
    recovered_edges <- scoped$edges
    pinned <- intersect(pinned, recovered_vertices)
  }

  latent_set <- .assign_latent(recovered_vertices, hidden_vertex, has_columns, g_matched, show_hidden_nodes,
                                always_observed = if (phase_mode) { character(0L) } else { unique(c(treat_grouped, out_grouped)) },
                                apsim_observed_ids = union(g_outputs, if (phase_mode) { g_matched } else { g_requested }))
  if (phase_mode) { latent_set <- setdiff(latent_set, c(phase_treatments, phase_outcomes)) }

  ## Record identity before graph grouping or display relabelling. One emitted
  ## vertex may represent several source rows. Retain every member explicitly.
  identity_rows <- .query_identity_rows(recovered_vertices, group_map, kept_ids)
  recombined <- list()
  if( isTRUE(node_recombination) ){
    recomb <- .esc_recombine(recovered_vertices, recovered_edges, latent_set, pinned)
    recovered_vertices <- recomb$vertices
    recovered_edges <- recomb$edges
    latent_set <- recomb$latent
    recombined <- recomb$supernodes
    for ( combined in names(recombined) ){
      take <- identity_rows$node %in% recombined[[ combined ]]
      identity_rows$node[ take ] <- combined
      identity_rows$quantity[ take ] <- .strip_slice(combined)
      identity_rows$occurrence[ take ] <- combined
    }
  }

  ## time_unroll = "collapse": the static projection. The "all" frame above recovered the lagged
  ## structure as per occasion vertices on an acyclic graph, and scoping / persistence / recombination
  ## have all run on that acyclic graph. As the FINAL graph transformation step, project every occasion back
  ## onto its static base node: strip the temporal suffix (_pre_treatment / _post_treatment / _t{s} /
  ## @t{s}) from every vertex and endpoint via .strip_slice, redirect every edge onto base ids, DROP the
  ## resulting self-loops -- these are exactly the collapsed succession / persistence edges (X_t ->
  ## X_{t+1}) and any same slice mediator that mapped onto itself, meaningless once time is removed --
  ## and dedupe parallel edges. A lagged feedback (X_t -> Y contemporaneously, Y -> X_{t+1} across the
  ## lag) reappears as a directed 2-cycle X -> Y -> X: the emitted graph MAY therefore be cyclic,
  ## which represents the static feedback topology (dagitty
  ## parses a cyclic dag{} string, and the downstream plotting layer warns rather than errors). This
  ## runs after the succession edges are added (so they collapse to dropped self-loops) and before the
  ## orphan sweep below (a collapsed base rarely orphans, but the sweep still runs). The single occasion
  ## relabelling below has no effect after this: no vertex carries a suffix any more.
  if ( identical(time_unroll, "collapse") ){
    identity_rows$node <- .strip_slice(identity_rows$node)
    ef <- .strip_slice(recovered_edges$from)
    et <- .strip_slice(recovered_edges$to)
    keep <- ef != et                                       # drop collapsed self-loops (succession/persistence)
    recovered_edges <- unique(data.table::data.table(from = ef[ keep ], to = et[ keep ]))
    recovered_vertices <- unique(.strip_slice(recovered_vertices))
    latent_set <- intersect(unique(.strip_slice(latent_set)), recovered_vertices)
    pinned <- intersect(unique(.strip_slice(pinned)),     recovered_vertices)
  }

  ## Final step: drop orphaned nodes (no incoming or outgoing edge) from the output - a
  ## documented internal decision, no parameter. Management levers are admitted and pinned into
  ## scope so their causal edges are recovered. Explicit endpoints and columns
  ## are exempt, so their requested representation is not silently discarded.
  connected <- unique(c(recovered_edges$from, recovered_edges$to))
  endpoint_pins <- intersect(recovered_vertices, requested_occurrences)
  if (identical(time_unroll, "collapse")) {
    endpoint_pins <- recovered_vertices[.strip_slice(recovered_vertices) %in% g_requested]
  }
  recovered_vertices <- intersect(recovered_vertices, union(connected, endpoint_pins))
  latent_set <- intersect(latent_set, recovered_vertices)
  identity_rows <- identity_rows[ identity_rows$node %in% recovered_vertices, , drop = FALSE ]

  ## Remove temporal suffixes from lone occasions: a temporal suffix only distinguishes a node from its sibling, so a
  ## base that survives with just ONE of its _t{s} occasions carries a suffix that marks nothing -
  ## strip it to the bare node_id. A base keeping >= 2 occasions stays suffixed. Applied to the final
  ## graph, after scoping / persistence / isolated node removal. Identification is unaffected:
  ## .strip_slice maps a bare id to itself, and the bare form of a lone suffixed node is never already
  ## present (a base is projected as either bare or suffixed, not both).
  suffixed <- recovered_vertices[ recovered_vertices != .strip_slice(recovered_vertices) ]
  occ_count <- table(.strip_slice(suffixed))
  lone_base <- names(occ_count)[ occ_count == 1L ]
  relabel <- stats::setNames(recovered_vertices, recovered_vertices)
  for ( v in suffixed ){
    bare <- .strip_slice(v)
    if ( ( bare %in% lone_base ) && !( bare %in% recovered_vertices ) ){
      relabel[[ v ]] <- bare
    }
  }
  if ( any(relabel != names(relabel)) ){
    identity_rows$node <- unname(relabel[ identity_rows$node ])
    recovered_vertices <- unname(relabel[ recovered_vertices ])
    recovered_edges <- unique(data.table::data.table(from = unname(relabel[ recovered_edges$from ]),
                                                       to   = unname(relabel[ recovered_edges$to ])))
    latent_set <- unique(unname(relabel[ latent_set ]))
    pinned <- unname(relabel[ pinned ])
  }

  treatments_v <- if (phase_mode) { intersect(recovered_vertices, unname(relabel[phase_treatments])) } else {
    recovered_vertices[.strip_slice(recovered_vertices) %in% treat_grouped]
  }
  outcomes_v <- if (phase_mode) { intersect(recovered_vertices, unname(relabel[phase_outcomes])) } else {
    recovered_vertices[.strip_slice(recovered_vertices) %in% out_grouped]
  }

  # These legacy rolled attributes describe the source process graph. A folded
  # path can cross multiple time steps, so filtering source arrows to retained
  # endpoints would silently omit its replacement. Preserve the source view
  # explicitly whenever folding occurred. The returned DAG is in `edges`.
  .pg_bases <- if (length(fold_nodes)) { unique(.strip_slice(recovery_folding$source_vertices)) } else { unique(.strip_slice(recovered_vertices)) }
  .pg_contemp <- g_contemp[ from %in% .pg_bases & to %in% .pg_bases ]
  .pg_lagged <- g_lagged[  from %in% .pg_bases & to %in% .pg_bases ]
  .pg_stocks <- intersect(persistent_ids, .pg_bases)
  .pg_treatments <- unique(.strip_slice(treatments_v))
  .pg_outcomes <- unique(.strip_slice(outcomes_v))
  .pg_latent <- unique(.strip_slice(latent_set))

  ## label= feature: relabel the recovered graph so a resolved input alias becomes the node's display
  ## label (see .apply_input_labels). Applied to vertices / edges / latent / exposure / outcome sets
  ## BEFORE .emit_dagitty, so the emitted text, the dagitty object, and the *_v attributes all carry
  ## the display labels. input_labels (display label -> node_id) is used to restore node identities.
  ## The NULL setting skips relabelling entirely, leaving the emitted text unchanged.
  input_labels <- character(0L)
  if ( !is.null(label) ){
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
  attr(dag, "time_unroll") <- if ( phase_mode ) { "selective" } else { time_unroll }
  attr(dag, "phase_window") <- if ( phase_mode ) { phase_window } else { NA_character_ }
  attr(dag, "apsim_nodes_detail") <- apsim_nodes_detail
  attr(dag, "supernode_overrides") <- supernode_overrides
  attr(dag, "force_acyclic") <- isTRUE(force_acyclic)
  attr(dag, "node_recombination") <- isTRUE(node_recombination)
  attr(dag, "randomised_treatment") <- isTRUE(randomised_treatment)
  attr(dag, "randomised_parents_cut") <- randomised_cut
  attr(dag, "query_depth") <- query_depth
  ## The full unrestricted ancestral graph, for get_valid_sets() to validate against. Latents are
  ## recomputed over the full vertex set (a node observed in the scoped graph is observed here too.
  ## the columns/hidden semantics are identical, only the vertex set is larger).
  attr(dag, "full_vertices") <- full_vertices
  attr(dag, "full_edges") <- full_edges
  attr(dag, "full_latent") <- .assign_latent(full_vertices, hidden_vertex, has_columns, g_matched,
                                                    show_hidden_nodes,
                                                    always_observed = if (phase_mode) { character(0L) } else { unique(c(treat_grouped, out_grouped)) },
                                                    apsim_observed_ids = union(g_outputs, if (phase_mode) { g_matched } else { g_requested }))
  if (phase_mode) {
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
  attr(dag, "rolled_structure") <- if (length(fold_nodes)) { "source_before_folding" } else { "retained_source" }
  attr(dag, "recovery_folding") <- recovery_folding
  ## Forward causal distance (treatment -> outcome) on the rolled base graph, for the process graph
  ## view's K control: the app reads $K as the minimum unroll window (the slider floor) and $fwd for
  ## the "reaches the outcome in N lag steps" label. Additive and read only, exactly like the
  ## process graph attributes above. The fallback keeps a sensible minimal two slice window if the
  ## depth computation errors, so the app can always render.
  attr(dag, "unroll_depth") <- tryCatch(
    .query_unroll_depth(.pg_bases, .pg_contemp, .pg_lagged, .pg_treatments, .pg_outcomes),
    error = function( e ){ # Retain the recorded fallback when depth calculation fails
      return(list(reachable = FALSE, fwd = 1L, conf_depth = 0L,
                                treat_slice = 0L, out_slice = 1L, K = 2L))
    }
  )
  attr(dag, "management_options") <- get_management_options(crop = crop, nodes = nodes, aliases = aliases, defaults = management_defaults)
  ## Versioned, complete construction receipt. This captures the query
  ## actually returned, not a different full graph or a reconstruction from tags.
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
                 phase_window = if (phase_mode) { phase_window } else { NULL },
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
  class(dag) <- unique(c("dag_by_query", class(dag)))

  if( isFALSE(quiet) ){
    optional <- sum(attr(dag, "management_options")$status == "optional")
    basis <- if ( phase_mode ){
      paste0("anchored to phases (", phase_window, " window)")
    } else {
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
                   length(present), " management levers present, ", optional, " more available - see get_management_options()."))
  }
  return(dag)
}


#' Options for a crop query
#'
#' Creates the options list accepted by `query_graph()`. Defaults can be used
#' without supplying any arguments.
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
query_params <- function( query_depth            = Inf,
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
                          ){ # Validate and collect the expert query options
  water_model <- match.arg(water_model)
  apsim_nodes_detail <- match.arg(apsim_nodes_detail)
  phase_window <- match.arg(phase_window)
  if (!is.logical(observe_apsim_outputs) || length(observe_apsim_outputs) != 1L || is.na(observe_apsim_outputs)) {
    stop("observe_apsim_outputs must be TRUE or FALSE.")
  }
  if ( !is.numeric(query_depth) || length(query_depth) != 1L || is.na(query_depth) || query_depth < 1 ){
    stop("query_depth must be a single number >= 1 (1 = the adjustment scope of causal paths + parents of treatment and outcome, 2 = the full ancestral subgraph, higher adds successive generations of context, Inf = the full graph). Please supply a valid value and try again.")
  }
  query_depth <- if ( is.infinite(query_depth) ){ Inf }else{ as.integer(query_depth) }
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



#####
# scm_to_dag.R
#####
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
scm_to_dag <- function( scm, crop = NULL, water_model = NULL, time_slices = NULL ){ # Convert the selected structural configuration into a finite DAG
  .datatable.aware <- TRUE
  if ( !is.data.frame(scm) || nrow(scm) == 0L ||
       !all(c("node_id", "function_of") %in% names(scm)) ){
    stop("scm must be a nonempty data.frame with node_id and function_of columns.")
  }
  nodes <- as.data.frame(scm, stringsAsFactors = FALSE)
  defaults <- list(concept_id = "", composition_of = "", crop_type = "agnostic",
                   tier = "", hidden_node = "FALSE")
  for ( column_name in names(defaults) ){
    if ( !column_name %in% names(nodes) ){ nodes[[ column_name ]] <- defaults[[ column_name ]] }
  }
  text_cols <- intersect(c("node_id", "concept_id", "function_of", "composition_of",
                           "init_of", "crop_type", "tier", "hidden_node"), names(nodes))
  for ( column_name in text_cols ){ nodes[[ column_name ]] <- as.character(nodes[[ column_name ]]) }
  if ( anyNA(nodes$node_id) || any(!nzchar(nodes$node_id)) || anyDuplicated(nodes$node_id) ){
    stop("scm$node_id must contain unique, nonempty ids.")
  }
  if ( any(grepl('["\\\\\n\r]', nodes$node_id)) ){
    stop("scm$node_id cannot contain quotes, backslashes or line breaks.")
  }
  for ( column_name in c("concept_id", "tier", "hidden_node") ){
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
  if ( is.null(crop) && length(choices) > 0L ){
    stop("Select crop explicitly when rows or parent qualifiers name particular crops: ",
         paste(choices, collapse = ", "), ".")
  }
  if ( !is.null(crop) && (!is.character(crop) || length(crop) != 1L || is.na(crop) ||
                          !nzchar(crop) || (length(choices) > 0L && !crop %in% choices)) ){
    stop("crop must select one crop represented by the table or its parent qualifiers.")
  }
  families <- c("SoilWater", "SWIM")
  families <- families[ vapply(families, function(x) { # Check whether this water model family occurs in the table
    return(any(startsWith(nodes$node_id, paste0(x, "."))))
  }, logical(1L)) ]
  if ( is.null(water_model) && length(families) > 1L ){
    stop("Select water_model = \"SoilWater\" or \"SWIM\". Alternative models cannot be combined.")
  }
  if ( !is.null(water_model) ){
    water_model <- match.arg(water_model, c("SoilWater", "SWIM"))
  }else if ( length(families) == 1L ){
    water_model <- families
  }
  kept <- .filter_nodes(nodes, if ( is.null(crop) ){ "agnostic" }else{ crop },
                        water_model, nodes$node_id)
  if ( length(kept) == 0L ){ stop("The selected configuration contains no nodes.") }
  ## Resolve using the complete vocabulary, then retain only the selected configuration.
  ## This permits legitimate references to alternative water families to be filtered out.
  active <- nodes[ nodes$node_id %in% kept, , drop = FALSE ]
  concept <- sub("^c:", "", active$concept_id)
  ambiguous <- unique(concept[ duplicated(concept) & nzchar(concept) ])
  for ( cell in unlist(active[ intersect(c("function_of", "composition_of", "init_of"), names(active)) ], use.names = FALSE) ){
    parsed_parents <- .parse_parents(cell)
    if ( nrow(parsed_parents) == 0L ){ next }
    parsed_parents <- parsed_parents[ vapply(parsed_parents$crops, .qualifier_admits, logical(1L), crop = crop), , drop = FALSE ]
    collisions <- intersect(parsed_parents$token, intersect(nodes$node_id, concept))
    conflicting <- vapply(collisions, function(token) { # Detect a token that identifies different nodes and concepts
      return(any(active$node_id[ concept == token ] != token))
    }, logical(1L))
    if ( any(conflicting) ){
      stop("A parent identifier has a concept/node_id collision and can resolve to multiple nodes. Make the conflicting identifiers distinct.")
    }
    if ( any(parsed_parents$token %in% ambiguous) ){
      stop("A parent concept resolves to multiple nodes. Use explicit node_id references.")
    }
  }
  for ( column_name in intersect(c("composition_of", "init_of"), names(active)) ){
    if ( any(vapply(active[[ column_name ]], function(cell) { # Detect lag annotations in this parent cell
      return(any(.parse_parents(cell)$lag > 0L))
    }, logical(1L))) ){
      stop("Lag annotations are supported only in function_of. Revise ", column_name, ".")
    }
  }
  ## Prefer the selected configuration when alternative families share a concept.
  ## Inactive rows remain vocabulary entries so explicit references to an inactive
  ## node can be resolved and excluded with that configuration, rather than guessed.
  nodes <- rbind(active, nodes[ !nodes$node_id %in% kept, , drop = FALSE ])
  ## Inactive rows are vocabulary entries only: their parent cells need not be resolved.
  inactive <- !nodes$node_id %in% kept
  for ( column_name in intersect(c("function_of", "composition_of", "init_of"), names(nodes)) ){
    nodes[[ column_name ]][ inactive ] <- ""
  }
  resolved <- withCallingHandlers(.resolve_edges(nodes, crop, preserve_lags = TRUE),
                                  warning = function(w) { # Reject unresolved structural references
                                    stop(conditionMessage(w), call. = FALSE)
                                  })
  selected <- lapply(resolved, function(e) { # Retain edges within the selected configuration
    return(e[ from %in% kept & to %in% kept ])
  })
  if ( !.is_acyclic(kept, selected$contemporaneous) ){
    stop("The selected contemporaneous graph contains a cycle. Revise the table or select a compatible configuration.")
  }
  stocks <- if ( "path_dependent" %in% names(active) ){
    active$node_id[ !is.na(active$path_dependent) & toupper(as.character(active$path_dependent)) == "TRUE" ]
  }else{ character(0L) }
  dynamic <- nrow(selected$lagged) > 0L || nrow(selected$initial) > 0L || length(stocks) > 0L
  max_lag <- if ( nrow(selected$lagged) > 0L ){ max(selected$lagged$lag) }else{ 0L }
  if ( is.null(time_slices) && dynamic ){
    stop("Supply time_slices explicitly for a table with lagged, initial or persistent structure.")
  }
  if ( !is.null(time_slices) && (!is.numeric(time_slices) || length(time_slices) != 1L ||
                                 is.na(time_slices) || !is.finite(time_slices) || time_slices < 1 ||
                                 time_slices != floor(time_slices) || time_slices > .Machine$integer.max || time_slices <= max_lag) ){
    stop("time_slices must be a positive integer greater than the largest retained lag (", max_lag, ").")
  }
  if ( is.null(time_slices) ){
    vertices <- kept
    edges <- selected$contemporaneous
  }else{
    time_slices <- as.integer(time_slices)
    lags <- selected$lagged
    if ( length(stocks) > 0L ){
      if ( !"lag" %in% names(lags) ){ lags$lag <- integer(0L) }
      lags <- unique(data.table::rbindlist(list(lags, data.table::data.table(from = stocks, to = stocks, lag = 1L))))
    }
    sliced <- .unroll_k_slices(kept, selected$contemporaneous, lags, time_slices)
    vertices <- sliced$vertices
    edges <- sliced$edges
    if ( nrow(selected$initial) > 0L ){
      edges <- unique(rbind(edges, data.table::data.table(from = paste0(selected$initial$from, "@t0"),
                                                  to = paste0(selected$initial$to, "@t0"))))
    }
  }
  if ( !.is_acyclic(vertices, edges) ){
    stop("The configured temporal graph contains a cycle, including initial parents. Revise the table.")
  }
  hidden <- active$node_id[ toupper(active$hidden_node) == "TRUE" ]
  latent <- if ( is.null(time_slices) ){
    intersect(vertices, hidden)
  }else if ( length(hidden) == 0L ){
    character(0L)
  }else{
    intersect(vertices, as.vector(vapply(seq_len(time_slices) - 1L,
                                  function(s) { # Name the hidden nodes at this time slice
                                    return(paste0(hidden, "@t", s))
                                  }, character(length(hidden)))))
  }
  code <- .emit_dagitty(vertices, edges, latent, character(0L), character(0L))
  dag <- if ( requireNamespace("dagitty", quietly = TRUE) ){ dagitty::dagitty(code) }else{ code }
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
#' a `query_graph()` result with its saved vertices and edges. A bare dagitty
#' object or graph string is insufficient. `nodes` defaults to `apsim_dag_nodes`
#' and must contain the display label column.
#'
#' @param dag A `query_graph()` result with its saved vertices and edges.
#' @param nodes Defaults to `apsim_dag_nodes` and must contain the display label column.
#' @details
#' Labels are matched to quantity keys with time suffixes removed, then the
#' `_t{s}` suffix is restored. Empty or conflicting labels fall back to the
#' quantity key so distinct nodes remain distinct. A display label identical
#' to the key leaves the name unchanged.
#' @return A dagitty object.
#' @export
readable_dag <- function( dag, nodes = apsim_dag_nodes ){ # Relabel a query graph while retaining its structure and receipt
  .datatable.aware <- TRUE
  vertices <- attr(dag, "recovered_vertices")
  edges <- attr(dag, "edges")
  latent <- attr(dag, "latent")
  exposures <- attr(dag, "treatments")
  outcomes <- attr(dag, "outcomes")
  if( is.null(vertices) || is.null(edges) ){
    stop("readable_dag() needs a DAG produced by query_graph() (it reads the recovered_vertices and edges attributes). Please pass a query_graph() result and try again.")
  }
  if( !("node_name_readable" %in% names(nodes)) ){
    stop("The node table must carry a 'node_name_readable' column to build display labels. Please pass nodes = apsim_dag_nodes and try again.")
  }

  ## suffix preserving split: base is the suffix stripped id, sfx is exactly the removed suffix
  base <- .strip_slice(vertices)
  suffix <- substring(vertices, nchar(base) + 1L)

  ## candidate label for each distinct base: node_name_readable, falling back to the node_id where the
  ## human label is missing. Then any label claimed by more than one base reverts to the node_id
  ## so the relabelling stays a bijection (a second pass covers the rare case where a human label
  ## coincides with another base's node_id).
  human_labels <- stats::setNames(as.character(nodes$node_name_readable), nodes$node_id)
  unique_bases <- unique(base)
  candidate_labels <- unname(human_labels[ unique_bases ])
  missing_labels <- is.na(candidate_labels) | !nzchar(candidate_labels)
  candidate_labels[ missing_labels ] <- unique_bases[ missing_labels ]
  candidate_labels <- stats::setNames(candidate_labels, unique_bases)
  for( pass in 1:2 ){
    duplicate_labels <- candidate_labels[ duplicated(candidate_labels) | duplicated(candidate_labels, fromLast = TRUE) ]
    if( length(duplicate_labels) == 0L ){
      break
    }
    clash <- candidate_labels %in% duplicate_labels
    candidate_labels[ clash ] <- names(candidate_labels)[ clash ]
  }

  relabel <- stats::setNames(paste0(unname(candidate_labels[ base ]), suffix), vertices)
  new_vertices <- unname(relabel[ vertices ])
  new_edges <- data.table::data.table(from = unname(relabel[ edges$from ]), to = unname(relabel[ edges$to ]))
  new_latent <- unname(relabel[ intersect(latent, vertices) ])
  new_exposures <- unname(relabel[ intersect(exposures, vertices) ])
  new_outcomes <- unname(relabel[ intersect(outcomes, vertices) ])

  graph_text <- .emit_dagitty(new_vertices, new_edges, new_latent, new_exposures, new_outcomes)
  result <- if (requireNamespace("dagitty", quietly = TRUE)) { dagitty::dagitty(graph_text) } else { graph_text }
  for ( key in setdiff(names(attributes(dag)), c("class", "names")) ){
    attr(result, key) <- attr(dag, key)
  }
  attr(result, "dag_text") <- graph_text
  attr(result, "recovered_vertices") <- new_vertices
  attr(result, "edges") <- new_edges
  attr(result, "latent") <- new_latent
  attr(result, "treatments") <- new_exposures
  attr(result, "outcomes") <- new_outcomes
  record <- attr(dag, "query_record")
  if (!is.null(record)){
    record$node_map$node <- unname(relabel[record$node_map$node])
    record$display_transform <- c(record$display_transform, list(list(
      operation = "readable_dag", mapping = relabel)))
    attr(result, "query_record") <- record
    canonical <- unique(record$node_map[c("node", "quantity")])
    attr(result, "input_labels") <- stats::setNames(canonical$quantity, .strip_slice(canonical$node))
  }
  if (inherits(dag, "dag_by_query")) { class(result) <- unique(c("dag_by_query", class(result))) }
  return(result)
}


#' List available management practices
#'
#' Returns each management quantity's label, aliases and status: `core` for
#' default inclusion or `optional` for explicit selection. `crop` restricts the
#' result to applicable quantities. The `nodes`, `aliases` and `defaults`
#' arguments use `apsim_dag_nodes`, `node_aliases` and `management_defaults`
#' unless supplied. Defaults establish management presence only.
#'
#' @param crop Restrict the result to applicable quantities.
#' @param nodes Defaults to `apsim_dag_nodes`.
#' @param aliases Defaults to `node_aliases`.
#' @param defaults Defaults to `management_defaults`.
#' @return a data.table with `node_id`, `label`, `aliases`, `status` and
#' `crop_type`. Rows marked `optional` are the available explicit choices.
#' @examples
#' \dontrun{
#' get_management_options(crop = "wheat")
#' }
#' @export
get_management_options <- function( crop     = NULL,
                                    nodes    = apsim_dag_nodes,
                                    aliases  = node_aliases,
                                    defaults = management_defaults
                                    ){ # List the available management nodes and their defaults

  .datatable.aware <- TRUE
  management <- nodes[ nodes$tier == "management", ]
  ## the management tier is nonempty by construction, so 1:nrow is safe
  rows <- lapply(1:nrow(management), function( i ){ # Describe one management node and its searchable aliases
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


#####
# Query engine internals
# Historical source: dag_by_query.R
#####


.query_source_order <- function(nodes){ # Retain slice ordering only when its complete node identities still match
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


## R/dag_by_query.R
## DAG-by-query layer for causaliflower: build a dagitty syntax causal DAG for a
## crop from the frozen APSIM-derived node table, then query it for covariate
## adjustment sets. Ported from the build_dag.py prototype with two corrections
## carried from the build plan (the soil_thermal_regime concept resolver fix and the
## function_of union composition_of structural edge set), and Step 5 implemented
## as a directed collapse of hidden nodes rather than a latent projection.
##
## Hidden nodes (hidden_node == "TRUE") are mechanistic detail. With
## show_hidden_nodes = FALSE they are removed from the recovered DAG but the directed
## reachability they induce is retained: a visible A gains an edge to every
## visible node reachable from it through hidden interiors. No bidirected
## edges are produced. The default show_hidden_nodes = TRUE retains them instead.
## Observation tags use supplied columns or the recorded availability assumption
## (see the matrix in query_graph()). Inclusion alone does not establish availability.

## Parse one function_of / composition_of cell into parent tokens and lags. A
## trailing @t-<k> marks a lag of k steps. A leading c: marks a concept reference.
.parse_parents <- function( cell ){ # Parse parent identifiers, crop qualifiers and time lags
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
  ## Crop qualifier, stripped AFTER the lag -- see the grammar note in build_node_layer.R.
  ## NA_character_ means unqualified, i.e. "every crop the child admits". A "+"-joined string
  ## names the crops the edge is restricted to.
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


## TRUE when a crop qualifier admits `crop`. An unqualified token (NA) admits every crop, and a
## NULL crop (no crop filtering, e.g. process_graph's whole model unroll) admits every edge.
.qualifier_admits <- function( crops, crop ){ # Check whether a parent qualifier admits the requested crop
  if( length(crop) == 0L || is.na(crops) ){
    return(TRUE)
  }
  return( crop %in% trimws(strsplit(crops, "+", fixed = TRUE)[[ 1L ]]) )
}


## Resolve structural edges as function_of UNION composition_of (part -> whole,
## lag 0). Tokens resolve to a node_id directly, or via a concept index keyed on
## the c:-stripped concept_id. The concept lookup supports tables whose concept
## keys differ from their node keys, including the bundled numbered concepts.
## function_of @t-1 tokens become lagged edges. Everything else is lag 0.
## Unresolved tokens carry NA and are dropped by the assembly filter.
## `crop` filters parent tokens by their crop qualifiers during resolution. This matters on
## two counts: it precedes the unique(from, to) dedup below, which would otherwise discard the crop
## information and force a union rule for a parent asserted both qualified and unqualified. It
## precedes every downstream consumer, so hidden node collapse can never recover reachability through
## an edge that does not exist for the queried crop. crop = NULL keeps every edge (the union), which
## is what process_graph.R's whole model unroll still calls with -- logged as an open item.
.resolve_edges <- function( nodes, crop = NULL, preserve_lags = FALSE ){ # Resolve parent references into contemporaneous, lagged and initial edges
  .datatable.aware <- TRUE
  by_id <- nodes$node_id
  concept_stripped <- sub("^c:", "", nodes$concept_id)
  first_seen <- !duplicated(concept_stripped) & nzchar(concept_stripped)
  by_concept <- stats::setNames(nodes$node_id[ first_seen ], concept_stripped[ first_seen ])
  has_init <- "init_of" %in% names(nodes)

  all_tokens <- unique(unlist(lapply(seq_len(nrow(nodes)), function( i ){ # Collect the admitted parent tokens for one node
    cells <- c(nodes$function_of[ i ], nodes$composition_of[ i ],
               if ( has_init ) { nodes$init_of[ i ] } else { NULL })
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

  ## init_of parse: the initial network B0. init_of cells carry c: tokens with no lag annotation. Each
  ## resolves to a lag-0 edge parent -> stock (the parents of a stock's origin slice value). Kept as a
  ## SEPARATE edge set so the daily contemporaneous / lagged transition graph is untouched - only the
  ## phase unroll consumes it, materialising these edges at the origin slice alone.
  initial_edges <- data.table::data.table(from = character(0L), to = character(0L))
  if ( has_init ){
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
    if ( length(init_records) > 0L ){
      initial_edges <- unique(data.table::as.data.table(do.call(rbind, init_records)))
    }
  }

  records <- suppressWarnings(
    unlist(
      ## 1:nrow safe by construction: nodes is the frozen node table, and an empty
      ## table is caught by the length(records) == 0L guard further down
      lapply(seq_len(nrow(nodes)), function( i ){ # Resolve the structural parents of one node
        child <- nodes$node_id[ i ]
        cells <- c(function_of = nodes$function_of[ i ], composition_of = nodes$composition_of[ i ])
        ## cells is always length 2 (function_of, composition_of), so 1:length is safe
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
  lagged_edges <- if ( isTRUE(preserve_lags) ){
    unique(edges[ is_lagged == TRUE, list(from, to, lag) ])
  }else{
    unique(edges[ is_lagged == TRUE, list(from, to) ])
  }
  return(list(contemporaneous = contemporaneous, lagged = lagged_edges, initial = initial_edges))
}


## Kahn topological sort over the contemporaneous (lag 0) edges: TRUE iff every
## node is ordered, i.e. the lag-0 structural graph is acyclic.
.is_acyclic <- function( node_ids, contemporaneous ){ # Check acyclicity by consuming vertices in topological order
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


## Node_ids UNAMBIGUOUSLY resolved from the supplied columns: a column counts only when it
## resolves to exactly one node (after the .resolve_one rules for breaking ties). Ambiguous columns (several
## surviving candidates) and unmatched columns are dropped here so they never silently bind to the
## wrong node. .resolve_columns() reports them for the user. NULL columns -> none.
.match_columns <- function( columns, alias_index, hidden_ids = character(0L) ){ # Keep supplied columns that resolve to exactly one node
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


## Management nodes present in the recovered DAG: those named in management_practices
## (which now defaults to the core set, so the core levers are present unless the user
## overrides it), plus any management node forced present by being a treatment/outcome
## or a supplied data column. Names resolve by alias or as direct node_ids.
.resolve_management_presence <- function( nodes, management_practices, forced, col_matched, alias_index ){ # Resolve management nodes admitted by defaults and explicit requests
  management_ids <- nodes$node_id[ nodes$tier == "management" ]
  candidates <- unique(c(management_practices,
                         .match_columns(management_practices, alias_index),
                         forced,
                         .match_columns(forced, alias_index),
                         col_matched))
  return(intersect(candidates, management_ids))
}


## Node deletion: drop nodes incompatible with the crop, the unused soil water family,
## and management nodes that are not present. Returns the surviving node_ids.
.filter_nodes <- function( nodes, crop, water_model, management_present ){ # Keep nodes admitted by the crop, water model and management choices
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


## Unroll K slices: each kept node becomes id@t0 .. id@t(K-1). Contemporaneous edges are copied
## within every slice. Each lagged edge advances by its recorded lag, defaulting to one when
## no lag column is supplied, and is retained only when its target is inside the window.
## Acyclic by construction for any K, since every lagged edge strictly increases the slice index and
## the contemporaneous subgraph is itself a DAG.
.unroll_k_slices <- function( kept_ids, contemporaneous, lagged, K ){ # Build the finite time slices and connect retained lagged edges
  .datatable.aware <- TRUE
  K <- as.integer(K)
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
  slices <- seq_len(K) - 1L                                   # 0 .. K-1
  vertices <- as.vector(vapply(slices, function( slice ){ # Name the retained nodes at this time slice
    return(paste0(kept_ids, "@t", slice))
  },
                               character(length(kept_ids))))
  within <- if ( nrow(contemporaneous_edges) > 0L ){
    data.table::rbindlist(lapply(slices, function( slice ){ # Copy contemporaneous edges into this time slice
      return(data.table::data.table(from = paste0(contemporaneous_edges$from, "@t", slice), to = paste0(contemporaneous_edges$to, "@t", slice)))
    }))
  } else {
    ## no contemporaneous edges among the kept nodes: emit none. Guards against paste0() recycling an
    ## empty base vector to "" (which would fabricate a phantom "@ts" self-loop). query_graph
    ## always supplies at least one contemporaneous edge, so this branch is only reached by callers
    ## such as the process graph unroll of a pure dynamics graph.
    data.table::data.table(from = character(0L), to = character(0L))
  }
  across <- if ( K >= 2L && nrow(lagged_edges) > 0L ){
    edge_lags <- if ( "lag" %in% names(lagged_edges) ){ lagged_edges$lag }else{ rep(1L, nrow(lagged_edges)) }
    data.table::rbindlist(lapply(slices[ -length(slices) ], function( slice ){ # Connect lagged edges that end within the finite window
      take <- which(slice + edge_lags < K)
      if ( length(take) == 0L ){ return(NULL) }
      return(data.table::data.table(from = paste0(lagged_edges$from[ take ], "@t", slice),
                            to = paste0(lagged_edges$to[ take ], "@t", slice + edge_lags[ take ])))
    }))
  } else {
    data.table::data.table(from = character(0L), to = character(0L))
  }
  edges <- unique(data.table::rbindlist(list(within, across)))
  return(list(vertices = vertices, edges = edges))
}


## Resolve PHASE_OF active_from / active_to tokens to integer phase index windows for a single crop.
## Returns a named list node_id -> c(from_idx, to_idx), suitable for passing as phase_active to
## .unroll_phase_slices. "terminal" is resolved to max(phase_index) for the crop. Nodes whose
## crop_type is not agnostic and does not include the crop are silently skipped (they will already be
## removed by .filter_nodes). Tokens are resolved via the canonical_stage vocabulary exactly as
## .stages_to_phase_indices does (canonical stage, APSIM phase name, or bare integer).
.resolve_phase_active <- function( phase_of, crop, crop_phases ) { # Resolve each node's active developmental phase interval
  crop_partition <- crop_phases[ crop_phases$crop == crop, ]
  last_phase <- max(crop_partition$phase_index)
  has_canonical_stage <- nzchar(crop_partition$canonical_stage)
  canonical_to_index <- stats::setNames(crop_partition$phase_index[ has_canonical_stage ], crop_partition$canonical_stage[ has_canonical_stage ])
  name_to_index <- stats::setNames(crop_partition$phase_index, crop_partition$phase_name)
  resolve_tok <- function( tok ) { # Resolve a phase token to its integer index
    if ( identical(tok, "terminal") ) { return( last_phase ) }
    if ( grepl("^[0-9]+$", tok) )    { return( as.integer(tok) ) }
    if ( tok %in% names(canonical_to_index) ) { return( unname(canonical_to_index[ tok ]) ) }
    if ( tok %in% names(name_to_index) )  { return( unname(name_to_index[ tok ]) ) }
    return( NA_integer_ )
  }
  active_intervals <- list()
  for ( i in seq_len(nrow(phase_of)) ) {
    node_id <- phase_of$node_id[ i ]
    crop_type <- phase_of$crop_type[ i ]
    if ( nzchar(crop_type) && crop_type != "agnostic" ) {
      crops <- trimws(strsplit(crop_type, ";", fixed = TRUE)[[ 1L ]])
      if ( !(crop %in% crops) ) { next }
    }
    first_phase <- resolve_tok(phase_of$active_from[ i ])
    last_active_phase <- resolve_tok(phase_of$active_to[ i ])
    if ( is.na(first_phase) || is.na(last_active_phase) ) { next }    # unresolvable token for this crop -> always active
    active_intervals[[ node_id ]] <- c(first_phase, last_active_phase)
  }
  return( active_intervals )
}


## Phase indexed unroll for the process graph DBN (absolute developmental time, origin = sowing). Unlike
## .unroll_k_slices, whose frame is anchored to treatment and has no sowing origin, this
## anchors slice 0 at sowing and materialises the initial network B0 there: for each stock, its init_of
## parents feed the stock's origin slice occasion. Slices 1..n_phases-1 are the crop's remaining
## developmental phases (count from crop_phases). Within a slice the contemporaneous (lag-0) structure is
## copied. Every lagged edge is carried forward across one phase boundary k -> k+1.
## This includes stock carryover, management feedback and faster physiological feedback.
## The phase lag is a discretisation choice, not an estimate of each process's duration.
## Acyclic by construction: each slice is a DAG,
## every cross slice edge strictly increases the phase index, and init edges live only at slice 0. Returns
## vertices and a typed edge table (edge_type in contemporaneous / lagged / initial).
.unroll_phase_slices <- function( kept_ids, contemporaneous, lagged, initial, n_phases,
                                  stock_ids, mgmt_ids = character(0L),
                                  phase_active = NULL ){ # Build developmental phase slices with initial and cross phase edges
  ## NOTE: stock_ids / mgmt_ids are retained in the signature (callers still compute and pass them, and
  ## .selective_phase_slice uses stock_ids for persistence) but are no longer consulted here - all lagged
  ## edges now carry forward regardless of target class (see the carry comment below).
  ## phase_active (Stage C masking): named list node_id -> c(from_idx, to_idx). When supplied, a node
  ## is included in slice s only if from_idx <= s <= to_idx. Nodes absent from phase_active are treated
  ## as always active (safe default for exogenous weather, static parameters, etc.). Masking is applied
  ## per slice to vertices and edges: an edge whose source or target is masked in a slice is dropped
  ## from that slice automatically (the from/to %in% active_ids filter handles it).
  .datatable.aware <- TRUE
  n_phases <- as.integer(n_phases)
  if ( length(n_phases) != 1L || is.na(n_phases) || n_phases < 1L ){
    stop("n_phases must be a single positive whole number of developmental phase slices. Please supply a phase count of at least 1 and try again.")
  }
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
  initial_edges <- initial[ from %in% kept_ids & to %in% kept_ids ]
  slices <- seq_len(n_phases) - 1L                            # 0 .. n_phases-1 (0 = sowing origin)

  ## Active node sets for each slice: if phase_active is NULL (no masking), every kept_id is active in every
  ## slice. Otherwise, compute the active subset per slice.
  if ( is.null(phase_active) ) {
    active_per_slice <- lapply(slices, function( slice ) { # Keep every node when no phase mask is supplied
      return(kept_ids)
    })
  } else {
    active_per_slice <- lapply(slices, function( slice ) { # Select the nodes active at this developmental phase
      mask <- vapply(kept_ids, function( node_id ) { # Check this node's active phase interval
        active_interval <- phase_active[[ node_id ]]
        if ( is.null(active_interval) ) { return( TRUE ) }           # not in phase_of -> always active
        return( active_interval[1L] <= slice && slice <= active_interval[2L] )
      }, logical(1L))
      return(kept_ids[ mask ])
    })
  }

  vertices <- unlist(lapply(slices, function( slice ) { # Name the active nodes at this developmental phase
    return(paste0(active_per_slice[[ slice + 1L ]], "@t", slice))
  }),
                     use.names = FALSE)

  ## within slice contemporaneous edges, per slice filtered to active nodes
  within <- if ( nrow(contemporaneous_edges) > 0L ) {
    data.table::rbindlist(lapply(slices, function( slice ) { # Copy contemporaneous edges between active nodes at this phase
      active_ids <- active_per_slice[[ slice + 1L ]]
      active_edges <- contemporaneous_edges[ from %in% active_ids & to %in% active_ids ]
      if ( nrow(active_edges) == 0L ) { return( data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L)) ) }
      return(data.table::data.table(from = paste0(active_edges$from, "@t", slice), to = paste0(active_edges$to, "@t", slice),
                             edge_type = "contemporaneous"))
    }))
  } else {
    data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
  }

  ## across phase carry: every lagged relationship becomes a forward phase boundary edge (from@ts ->
  ## to@t{s+1}) for each adjacent pair. This includes edges into a stock (accumulation: stock(t+1) =
  ## stock(t) + net_flux(t)), into a per_event management lever (observational management feedback, cut
  ## only under do()), AND the fast within phase physiological feedbacks (e.g. live_leaf_n ->
  ## n_translocation, rue_actual -> cgr_radn). The last group cannot be represented within a phase without
  ## recreating a within phase cycle, and dropping it severs real forward causal chains (N -> ... -> grain
  ## runs through these edges), so at phase resolution the acyclic approximation is a single phase
  ## forward lag - the same treatment the day frame gives every lagged edge. Because every carried edge is
  ## strictly forward (t -> t+1) and never enters the within slice contemporaneous graph, the unroll stays
  ## a DAG regardless. NOTE (for agronomy check): the one phase lag on the fast feedbacks is a
  ## discretisation choice, not a claim that the coupling literally takes a full phase. Verify the
  ## direction and presence of each carried edge against agronomy sources, not its exact lag length.
  carry <- lagged_edges
  across <- if ( n_phases >= 2L && nrow(carry) > 0L ) {
    data.table::rbindlist(lapply(slices[ -length(slices) ], function( slice ) { # Carry lagged edges between adjacent active phases
      active_sources <- active_per_slice[[ slice + 1L ]]
      active_targets <- active_per_slice[[ slice + 2L ]]
      active_edges <- carry[ from %in% active_sources & to %in% active_targets ]
      if ( nrow(active_edges) == 0L ) { return( data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L)) ) }
      return(data.table::data.table(from = paste0(active_edges$from, "@t", slice), to = paste0(active_edges$to, "@t", slice + 1L),
                             edge_type = "lagged"))
    }))
  } else {
    data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
  }

  ## initial network B0: materialised at the origin slice only (parent@t0 -> stock@t0), filtered to
  ## nodes active in slice 0
  origin_nodes <- active_per_slice[[ 1L ]]
  origin_network <- if ( nrow(initial_edges) > 0L ) {
    origin_edges <- initial_edges[ from %in% origin_nodes & to %in% origin_nodes ]
    if ( nrow(origin_edges) == 0L ) {
      data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
    } else {
      data.table::data.table(from = paste0(origin_edges$from, "@t0"), to = paste0(origin_edges$to, "@t0"),
                             edge_type = "initial")
    }
  } else {
    data.table::data.table(from = character(0L), to = character(0L), edge_type = character(0L))
  }

  edges <- unique(data.table::rbindlist(list(within, across, origin_network)))
  phase_graph <- list(vertices = vertices, edges = edges)
  if ( !is.null(phase_active) ) {
    attr(phase_graph, "active_per_slice") <- vapply(active_per_slice, length, integer(1L))
  }
  return( phase_graph )
}


## A single contemporaneous slice has no @t suffix. Every lagged edge becomes an
## ordinary edge alongside the contemporaneous edges. No edge is dropped.
## Feedback can therefore produce cycles: fertiliser and irrigation decisions
## read the soil state they change, and soil nitrogen and microbial pools also
## contain feedback. The caller is warned and can unroll time instead.
.contemporaneous_slice <- function( kept_ids, contemporaneous, lagged = NULL ){ # Collapse lagged edges onto a single contemporaneous view
  .datatable.aware <- TRUE
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  edges <- data.table::data.table(from = contemporaneous_edges$from, to = contemporaneous_edges$to)

  n_collapsed <- 0L
  if ( !is.null(lagged) && nrow(lagged) > 0L ){
    lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
    if ( nrow(lagged_edges) > 0L ){
      n_collapsed <- nrow(unique(data.table::data.table(from = lagged_edges$from, to = lagged_edges$to)))
      edges <- rbind(edges, data.table::data.table(from = lagged_edges$from, to = lagged_edges$to))
    }
  }
  edges <- unique(edges)
  ## Defensive only: .resolve_edges()$lagged carries no self-loops (persistence
  ## is derived separately, from path_dependent), so nothing is removed here today. It guards the case
  ## where a future X@t-1 -> X lagged edge would collapse to an unrepresentable X -> X on one slice.
  edges <- edges[ from != to ]

  return(list(vertices    = kept_ids,
              edges       = edges,
              n_collapsed = n_collapsed,
              acyclic     = .is_acyclic(kept_ids, edges)))
}


## Adaptive unrolling depth for one query. The forward reach of a treatment's effect is the shortest
## lag distance to an outcome (contemporaneous edges cost 0 slices, lagged edges cost 1). The window
## also carries one slice of history before treatment when the treatment has a lagged parent (the prior
## state that informs the decision, i.e. a confounder occasion). The treatment is anchored at
## treat_slice = conf_depth, the outcome arrives at out_slice = treat_slice + fwd_dist, and the window
## spans t0 .. t(K-1) with K = out_slice + 1. Distances are computed before unrolling.
## Bellman-Ford converges because contemporaneous edges are a DAG and lagged edges are strictly
## positive, so there is no zero/negative cycle.
##
## fwd_dist uses the slowest treatment and outcome pair: max over treatments of max over outcomes.
## Taking the nearest pair can truncate other causal paths at the window boundary.
## For example, adding grain_protein must not remove a slower path to grain_yield.
## Taking the maximum leaves enough slices for every reachable pair.
## Single treatment and single outcome queries are unaffected (max of one element is that element).
.query_unroll_depth <- function( kept_ids, contemporaneous, lagged, treat_ids, outcome_ids ){ # Size the temporal window for all treatment and outcome pairs
  .datatable.aware <- TRUE
  contemporaneous_edges <- contemporaneous[ from %in% kept_ids & to %in% kept_ids ]
  lagged_edges <- lagged[ from %in% kept_ids & to %in% kept_ids ]
  unreachable_distance <- .Machine$integer.max %/% 2L
  treatment_seeds <- intersect(treat_ids, kept_ids)
  selected_outcomes <- intersect(outcome_ids, kept_ids)
  if ( length(treatment_seeds) == 0L || length(selected_outcomes) == 0L ){
    return(list(reachable = FALSE, fwd = unreachable_distance, conf_depth = 0L, treat_slice = 0L, out_slice = 1L, K = 2L))
  }

  ## shortest lag distance from one treatment to every kept node
  .lag_dist <- function( src ){ # Find the shortest lag distance from one treatment
    distances <- stats::setNames(rep(unreachable_distance, length(kept_ids)), kept_ids)
    distances[ src ] <- 0L
    repeat {
      changed <- FALSE
      if ( nrow(contemporaneous_edges) > 0L ){
        better <- which(distances[ contemporaneous_edges$from ] < distances[ contemporaneous_edges$to ])
        for ( i in better ){ distances[[ contemporaneous_edges$to[i] ]] <- distances[[ contemporaneous_edges$from[i] ]] }
        changed <- changed || length(better) > 0L
      }
      if ( nrow(lagged_edges) > 0L ){
        better <- which(distances[ lagged_edges$from ] + 1L < distances[ lagged_edges$to ])
        for ( i in better ){ distances[[ lagged_edges$to[i] ]] <- distances[[ lagged_edges$from[i] ]] + 1L }
        changed <- changed || length(better) > 0L
      }
      if ( !changed ){
        break
      }
    }
    return(distances)
  }

  ## relax from each treatment separately, then take the maximum over every (treatment, outcome)
  ## pair. A pair with no path in the base graph stays disconnected for any K,
  ## so it is excluded from the max and reported via reachable = FALSE.
  pair_distances <- integer(0L)
  has_unreachable_pair <- FALSE
  for ( s in treatment_seeds ){
    distances <- .lag_dist(s)
    for ( o in selected_outcomes ){
      pair_distance <- distances[[ o ]]
      if ( pair_distance >= unreachable_distance ){
        has_unreachable_pair <- TRUE
      } else {
        pair_distances <- c(pair_distances, as.integer(pair_distance))
      }
    }
  }
  if ( length(pair_distances) == 0L ){
    ## no outcome is reachable from any treatment: fall back to a single post slice so the
    ## query still returns a graph (the reachability harness assertion is what flags this case).
    forward_depth <- 1L
    reachable <- FALSE
  } else {
    forward_depth <- max(pair_distances)
    reachable <- !has_unreachable_pair
  }
  history_depth <- as.integer(any(lagged_edges$to %in% treat_ids))
  treatment_slice <- history_depth
  outcome_slice <- treatment_slice + forward_depth
  return(list(reachable = reachable, fwd = as.integer(forward_depth), conf_depth = history_depth,
              treat_slice = treatment_slice, out_slice = outcome_slice, K = outcome_slice + 1L))
}


## Selective time unroll (adaptive depth). Build the full unroll of K slices internally (K, treat_slice,
## out_slice from .query_unroll_depth), collapse hidden nodes on it so cross slice reachability
## through hidden interiors stays ordered in time, then project back to a graph in which only the
## bases that straddle the intervention keep occasions for each slice. Treatment descendant status is
## read off the unrolled graph BEFORE the hidden collapse. A base with occasions both outside
## and within the treatment's descendants straddles treatment and is split into <id>_t{s} (absolute slice s).
## A confounder entirely before treatment, or a base that never straddles, collapses to a single bare <id>. Treatment
## and outcome are anchored to one slice each (treat_slice / out_slice) and stay bare - an
## intervention is applied once, an outcome is read once (VanderWeele & Robins 2007). Any base left
## in a directed cycle after projection is promoted (also split) and the projection recomputed. The
## fully split unroll is acyclic, so the loop terminates. Returns the visible projected graph
## plus the split / promoted base sets, the slice anchors, reachability, and persistence edges.
.selective_slice <- function( kept_ids, contemporaneous, lagged, hidden_ids, treat_ids, outcome_ids,
                              persistent_ids = character(0L), split_all = FALSE,
                              label_mode = c("slice", "pre_post"),
                              show_hidden_nodes = FALSE ){ # Project the temporal frame while retaining required occasions
  .datatable.aware <- TRUE
  label_mode <- match.arg(label_mode)
  never_split <- union(treat_ids, outcome_ids)

  ## Adaptive window depth: the treatment is anchored at treat_slice, the outcome arrives at
  ## out_slice, the window spans t0 .. t(K-1). Interventions have no instantaneous effect (their
  ## effect edges are lagged in the base graph), so the effect first reaches the outcome fwd_dist
  ## lag hops after application. One slice of history before treatment is carried when the treatment has
  ## a lagged parent (the prior state that informs the decision - a confounder occasion).
  depth <- .query_unroll_depth(kept_ids, contemporaneous, lagged, treat_ids, outcome_ids)
  K <- depth$K
  treat_slice <- depth$treat_slice
  out_slice <- depth$out_slice

  ## full unroll of K slices: the working graph, never emitted
  unrolled <- .unroll_k_slices(kept_ids, contemporaneous, lagged, K)
  hidden_vertex <- stats::setNames(.strip_slice(unrolled$vertices) %in% hidden_ids, unrolled$vertices)

  ## forward reachability from each treatment's anchor (treat_slice) vertex over the unrolled edges
  ## (through hidden interiors): each vertex's treatment descendant status.
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

  ## A base is SPLIT into occasions for each slice <id>_t{s} when it straddles the intervention: some
  ## occasion is not a descendant of the treatment (before treatment, off the effect) and some occasion
  ## is a descendant (the effect has reached it). A base that is entirely before treatment (a static / exogenous
  ## confounder, one value throughout) or never straddles collapses to a single bare node. Treatment
  ## and outcome never split - each is anchored to a single slice below.
  base_ids <- unique(.strip_slice(unrolled$vertices))
  straddles <- function( b ){ # Check whether a base node straddles the treatment
    occasions <- intersect(paste0(b, "@t", seq_len(K) - 1L), unrolled$vertices)
    return(any(!descendant[ occasions ]) && any(descendant[ occasions ]))
  }
  straddlers <- setdiff(base_ids[ vapply(base_ids, straddles, logical(1L)) ], never_split)

  ## collapse hidden on the unrolled (ordered in time) graph: visible, still acyclic.
  ## When show_hidden_nodes = TRUE, skip the collapse so hidden occasions stay in the graph
  ## and are marked [latent] downstream. The reachability from treatment already
  ## traversed through hidden interiors (line 654), so straddler detection is unaffected.
  collapsed <- if ( isTRUE(show_hidden_nodes) ){
    list(vertices = unrolled$vertices, edges = unrolled$edges)
  } else {
    .collapse_hidden(unrolled$vertices, unrolled$edges, hidden_vertex)
  }
  vertex_bases <- .strip_slice(collapsed$vertices)
  vertex_slices <- as.integer(sub("^.*@t([0-9]+)$", "\\1", collapsed$vertices))

  ## Anchor the treatment to treat_slice and the outcome to out_slice: an intervention is applied
  ## once, and the outcome is a terminal measurement read once, at the slice the effect arrives
  ## (VanderWeele & Robins 2007). Their occasions at every other slice, and the edges on them, are
  ## dropped before contracting, leaving a single bare treatment / outcome node.
  anchor_drop <- c(collapsed$vertices[ vertex_bases %in% treat_ids   & vertex_slices != treat_slice ],
                   collapsed$vertices[ vertex_bases %in% outcome_ids & vertex_slices != out_slice ])
  keep_v0 <- setdiff(collapsed$vertices, anchor_drop)

  ## Project to a single graph. Straddlers and mediators promoted after feedback retain <id>_t{s},
  ## everything else is single/bare. Any base left in a directed cycle is promoted and the projection
  ## recomputed. The fully split unroll is acyclic, so this terminates.
  ##
  ## split_all = TRUE is time_unroll = "all": the UNCOMPRESSED frame. Same window K, same scope, same
  ## anchoring, same hidden collapse -- the ONLY difference is that every included base is split at
  ## every slice instead of only the straddlers. Compression is what makes a lagged edge lack an
  ## occasion to land on. With nothing compressed, every lagged edge between two included occasions
  ## materialises. (Edges can still be absent because an OCCASION is out of scope -- e.g. a descendant
  ## of the treatment that is not an ancestor of the outcome is pruned by the ancestral restriction,
  ## and its edges go with it. That is scoping, not deletion, and it is correct: keeping such a node
  ## would invite adjusting for a variable measured after treatment.) The promotion loop has no effect here: with
  ## every base split, the projection IS the unrolled graph, which is acyclic by construction.
  ##
  ## label_mode = "pre_post" is time_unroll = "pre_post": the same straddler split, but each straddler
  ## collapses to just TWO occasions -- <id>_pre_treatment and <id>_post_treatment -- bucketed by whether
  ## the treatment's effect has reached that occasion (descendant status), rather than one <id>_t{s} per
  ## slice. A 2-bucket collapse is NOT acyclic in general: a feedback that recurs entirely within the
  ## period after treatment (X@t1 -> Y@t1 contemporaneous, Y@t1 -> X@t2 lagged, both descendants) folds to
  ## X_post <-> Y_post. Any base caught in a residual cycle
  ## is promoted back to per slice labelling (which is acyclic on the unrolled frame), so the emitted
  ## graph is always a DAG -- mostly pre/post, with the rare base with feedback within a group kept at <id>_t{s}.
  ## descendant status per collapsed vertex (a subset of the unrolled vertices, so every id is present).
  post_treatment <- descendant[ collapsed$vertices ]
  promoted <- character(0L)                          # bases that fall back to per slice labelling (a fold could not break the cycle)
  folded_drop <- character(0L)                          # interior occasions marginalised out by the pre_post pre/post fold
  base_split_set <- if ( isTRUE(split_all) ){ setdiff(base_ids, never_split) } else { straddlers }
  slice_of <- stats::setNames(vertex_slices, collapsed$vertices)
  repeat {
    split_set <- union(base_split_set, promoted)
    keep_now <- setdiff(keep_v0, folded_drop)
    is_split <- vertex_bases %in% split_set
    proj_label <- vertex_bases
    if ( identical(label_mode, "pre_post") ){
      ## straddlers collapse to _pre_treatment / _post_treatment by descendant status. A base whose
      ## fold could not break a residual cycle falls back to per slice (promoted). Interior occasions
      ## dropped by the fold (folded_drop) are excluded from keep_now and marginalised out below.
      use_slice <- is_split & (vertex_bases %in% promoted)
      use_pp <- is_split & !(vertex_bases %in% promoted)
      proj_label[ use_slice ] <- paste0(vertex_bases[ use_slice ], "_t", vertex_slices[ use_slice ])
      proj_label[ use_pp ] <- paste0(vertex_bases[ use_pp ],
                                        ifelse(post_treatment[ use_pp ], "_post_treatment", "_pre_treatment"))
    } else {
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
    if ( identical(label_mode, "pre_post") ){
      ## FOLD the cyclic bases to pre/post rather than promoting them to per slice. Drop each cyclic
      ## base's interior occasions (keep the earliest -> _pre_treatment and the latest -> _post_treatment).
      ## Marginalising the interiors turns the cross slice feedback that had folded to a within group
      ## cycle back into a forward _pre_treatment -> _post_treatment edge, which is acyclic. A base with
      ## < 3 kept occasions has no interior to drop, so it cannot be folded -- fall back to per slice for
      ## it (correctness before compactness) and restore its interiors. The loop projects again and
      ## checks again, so the emitted graph is always a DAG.
      progress <- FALSE
      for ( b in setdiff(cyc_bases, promoted) ){
        occasions <- intersect(collapsed$vertices[ vertex_bases == b ], keep_now)
        occasions <- occasions[ order(slice_of[ occasions ]) ]
        if ( length(occasions) >= 3L ){
          folded_drop <- union(folded_drop, occasions[ -c(1L, length(occasions)) ])
          progress <- TRUE
        } else {
          promoted <- union(promoted, b)
          folded_drop <- setdiff(folded_drop, collapsed$vertices[ vertex_bases == b ])
          progress <- TRUE
        }
      }
      if( !progress ){
        break
      }
    } else {
      new_promote <- setdiff(cyc_bases, split_set)
      if( length(new_promote) == 0L ){
        break
      }
      promoted <- union(promoted, new_promote)
    }
  }
  keep_now <- setdiff(keep_v0, folded_drop)

  ## temporal persistence edges (returned separately, applied after any query_depth scoping). The
  ## earlier occasion of a path dependent STATE variable causally precedes the next (forward Euler
  ## identity: stock(t+1) = stock(t) + net_flux(t)). A flux or instantaneous ratio gets none. Emitted
  ## between consecutive surviving occasions, always earlier -> later, so they cannot create a cycle.
  ## Persistence edges are emitted in the SAME labelling the projection used, so they attach to real
  ## vertices in either mode: each stock's surviving occasions are taken in slice order, mapped through
  ## the final projection (pmap), deduplicated where consecutive occasions share a label (two slices in
  ## the same pre/post bucket), and linked earlier -> later. In slice mode every occasion has its own
  ## <id>_t{s}, so this is the old <id>_t{s} -> <id>_t{s+1} chain. In pre_post mode a straddler stock
  ## carries a single <id>_pre_treatment -> <id>_post_treatment. Always forward, so it cannot add a cycle.
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


## Normalise a growth stage spec into a canonical stage map for each node. `spec` may be NULL (no anchoring),
## a bare character scalar/vector of canonical stage tokens, applied to EVERY id in `node_ids`, or a named
## list mapping node_id -> one or more canonical stage tokens (ids absent from the list are left out).
## Returns a named list node_id -> character vector of stage tokens (empty when spec is NULL).
.normalize_stage_spec <- function( spec, node_ids ){ # Normalise growth stage requests into a map for each endpoint
  if ( is.null(spec) || length(spec) == 0L ){
    return(stats::setNames(list(), character(0L)))
  }
  if ( is.list(spec) ){
    nm <- names(spec)
    if ( is.null(nm) || any(!nzchar(nm)) ){
      stop("A growth stage list must be named by treatment/outcome node id, e.g. list(nitrogen_applied = \"flowering\"). Please name every element and try again.")
    }
    unknown <- setdiff(nm, node_ids)
    if ( length(unknown) > 0L ){
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


## Resolve a per node canonical stage map to per node integer phase indices for `crop`, via crop_phases.
## Each token is matched against the crop's canonical_stage vocabulary (the interpretable, public
## handle). A bare integer token is accepted as a literal phase index, and an APSIM phase name resolves
## too (both supported in the backend). Returns a named list node_id -> sorted unique integer indices.
## Errors name the valid canonical stages for the crop.
.stages_to_phase_indices <- function( stage_map, crop, crop_phases ){ # Resolve growth stage tokens to sorted phase indices
  crop_partition <- crop_phases[ crop_phases$crop == crop, ]
  if ( nrow(crop_partition) == 0L ){
    stop(sprintf("crop '%s' has no crop_phases partition, so growth stages cannot be resolved. Please supply a crop with a phase partition and try again.", crop))
  }
  has_canonical_stage <- nzchar(crop_partition$canonical_stage)
  canonical_to_index <- stats::setNames(crop_partition$phase_index[ has_canonical_stage ], crop_partition$canonical_stage[ has_canonical_stage ])
  name_to_index <- stats::setNames(crop_partition$phase_index, crop_partition$phase_name)
  valid_stages <- sort(unique(crop_partition$canonical_stage[ has_canonical_stage ]))
  resolve_one <- function( tok ){ # Resolve and validate one growth stage token
    if ( grepl("^[0-9]+$", tok) ){
      phase_index <- as.integer(tok)
      if ( !(phase_index %in% crop_partition$phase_index) ){
        stop(sprintf("phase index %d is out of range for %s (valid 0..%d). Please supply a valid phase and try again.",
                     phase_index, crop, max(crop_partition$phase_index)))
      }
      return(phase_index)
    }
    if ( tok %in% names(canonical_to_index) ){ return(unname(canonical_to_index[ tok ])) }
    if ( tok %in% names(name_to_index)  ){ return(unname(name_to_index[ tok ])) }
    stop(sprintf("growth stage '%s' is not a canonical stage for %s. Valid canonical stages: %s. Please use one of these (or a phase index) and try again.",
                 tok, crop, paste(valid_stages, collapse = ", ")))
  }
  return(lapply(stage_map, function( toks ) { # Resolve the distinct phases requested for one endpoint
    return(sort(unique(vapply(toks, resolve_one, integer(1L)))))
  } ))
}


## Phase anchored selective slice: the developmental phase sibling of .selective_slice. Instead of a
## floating window of K slices anchored to treatment it unrolls the crop's phase partition (via
## .unroll_phase_slices, so the initial network B0 is in the frame with the sowing origin at slice 0),
## and anchors each treatment and outcome at its selected developmental phases.
## A base is an endpoint only at those phases. Outcome histories remain distinct states.
## Splitting, hidden node collapse, cycle promotion and persistence then use these phase occasions.
## `window` chooses the retained span: "origin" keeps slices 0..hi (initial conditions through B0).
## "compact" keeps (earliest treatment phase - 1)..hi (one phase of confounder history before treatment,
## no origin/B0). Returns the visible projected graph plus the split/promoted sets, the window, and
## persistence edges - the same shape .selective_slice returns, so the query_graph tail consumes it as is.
.selective_phase_slice <- function( kept_ids, contemporaneous, lagged, initial, hidden_ids,
                                    treat_phase_map, out_phase_map, n_phases, stock_ids, mgmt_ids,
                                    persistent_ids = character(0L), window = c("origin", "compact"),
                                    phase_active = NULL, show_hidden_nodes = FALSE ){ # Project the phase frame while preserving selected endpoint histories
  .datatable.aware <- TRUE
  window <- match.arg(window)
  treat_ids <- names(treat_phase_map)
  outcome_ids <- names(out_phase_map)
  ## Single anchors used only as treatments can remain bare. An outcome's other occasions
  ## are still causal states: keep all of them split so they cannot be merged
  ## into the measured occasion. This also separates same base treatment/outcome
  ## anchors. Final single occasion relabelling may remove an unnecessary suffix.
  single_anchor <- c(treat_ids[ lengths(treat_phase_map) == 1L ], outcome_ids[ lengths(out_phase_map) == 1L ])
  multi_anchor <- c(treat_ids[ lengths(treat_phase_map) >  1L ], outcome_ids[ lengths(out_phase_map) >  1L ])
  never_split <- setdiff(single_anchor, outcome_ids)

  ## slice window [lo..hi]: hi = latest anchored phase, lo = origin (0) under "origin", else one phase of
  ## history before treatment (earliest treatment phase - 1, floored at 0) under "compact".
  all_anchor <- c(unlist(treat_phase_map, use.names = FALSE), unlist(out_phase_map, use.names = FALSE))
  last_slice <- max(all_anchor)
  first_slice <- if ( identical(window, "origin") ){ 0L } else { max(0L, min(unlist(treat_phase_map, use.names = FALSE)) - 1L) }

  ## full phase unroll (working graph, never emitted). edge_type dropped - the recovered query graph is a
  ## plain DAG, and B0 init edges flow through as ordinary from->to edges into the origin slice stock.
  full <- .unroll_phase_slices(kept_ids, contemporaneous, lagged, initial, n_phases, stock_ids, mgmt_ids,
                                     phase_active = phase_active)
  all_vertex_slices <- as.integer(sub("^.*@t([0-9]+)$", "\\1", full$vertices))
  window_vertices <- full$vertices[ all_vertex_slices >= first_slice & all_vertex_slices <= last_slice ]
  window_edges <- full$edges[ full$edges$from %in% window_vertices & full$edges$to %in% window_vertices, list(from, to) ]

  hidden_vertex <- stats::setNames(.strip_slice(window_vertices) %in% hidden_ids, window_vertices)

  ## descendants of treatment over the windowed edges, from every treatment's anchored phase occasions
  children_by_parent <- split(window_edges$to, window_edges$from)
  reach_forward <- function( sources ){ # Collect descendants reached from the selected treatment occasions
    seen <- character(0L)
    stack <- unlist(children_by_parent[ sources ], use.names = FALSE)
    while ( length(stack) > 0L ){
      vertex <- stack[ length(stack) ]
      stack <- stack[ -length(stack) ]
      if ( vertex %in% seen ){
        next
      }
      seen <- c(seen, vertex)
      children <- children_by_parent[[ vertex ]]
      if ( !is.null(children) ){
        stack <- c(stack, children)
      }
    }
    return(seen)
  }
  treatment_sources <- character(0L)
  for ( t in treat_ids ){ treatment_sources <- c(treatment_sources, paste0(t, "@t", treat_phase_map[[ t ]])) }
  treatment_sources <- intersect(treatment_sources, window_vertices)
  treatment_descendants <- union(treatment_sources, reach_forward(treatment_sources))
  descendant <- stats::setNames(window_vertices %in% treatment_descendants, window_vertices)

  ## straddlers: a base with occasions both outside and within treatment's descendants in the window
  base_ids <- unique(.strip_slice(window_vertices))
  window_slices <- first_slice:last_slice
  straddles <- function( b ){ # Check whether a base node straddles the treatment
    node_occasions <- intersect(paste0(b, "@t", window_slices), window_vertices)
    return(any(!descendant[ node_occasions ]) && any(descendant[ node_occasions ]))
  }
  straddlers <- setdiff(base_ids[ vapply(base_ids, straddles, logical(1L)) ], never_split)

  ## collapse hidden on the windowed (ordered in time) graph (skip when showing hidden nodes)
  collapsed <- if ( isTRUE(show_hidden_nodes) ){
    list(vertices = window_vertices, edges = window_edges)
  } else {
    .collapse_hidden(window_vertices, window_edges, hidden_vertex)
  }
  vertex_bases <- .strip_slice(collapsed$vertices)
  vertex_slices <- as.integer(sub("^.*@t([0-9]+)$", "\\1", collapsed$vertices))

  ## Preserve the existing event selection rule for treatment only bases. Never
  ## delete an outcome's other occasions: they may mediate its response or be
  ## common causes. Selecting an endpoint does not erase its state history.
  drop_anchor <- character(0L)
  for ( t in setdiff(treat_ids, outcome_ids) ){
    drop_anchor <- c(drop_anchor, collapsed$vertices[ vertex_bases == t & !(vertex_slices %in% treat_phase_map[[ t ]]) ])
  }
  keep_v0 <- setdiff(collapsed$vertices, drop_anchor)

  ## Project straddlers and mediators promoted after feedback to <id>_t{phase}, leaving everything else bare. Any base
  ## left in a directed cycle after projection is promoted and the projection recomputed. The fully split
  ## unroll is acyclic, so this terminates.
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
    if ( length(cyclic) == 0L ){
      break
    }
    promote <- setdiff(unique(.strip_slice(unlist(cyclic))), union(never_split, split_set))
    if ( length(promote) == 0L ){
      break
    }
    promoted <- union(promoted, promote)
    split_set <- union(union(union(straddlers, multi_anchor), outcome_ids), promoted)
  }

  ## persistence for split path dependent stocks (earlier phase occasion -> next), forward only
  persist_from <- character(0L)
  persist_to <- character(0L)
  for ( base in intersect(split_set, persistent_ids) ){
    occasion_slices <- sort(unique(vertex_slices[ vertex_bases == base & (collapsed$vertices %in% keep_v0) ]))
    if ( length(occasion_slices) >= 2L ){
      for ( i in seq_len(length(occasion_slices) - 1L) ){
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


## Directed collapse: remove hidden vertices, and for each visible vertex add an
## edge to every visible vertex reachable from it through hidden interiors
## (any chain length). No bidirected edges are produced. Isolated visible
## vertices are retained (declared by the emitter) but gain no edges.
.query_fold_nodes <- function(vertices, edges, remove){ # Fold removed nodes and record the changed directed paths
  remove <- intersect(vertices, remove)
  result <- if (length(remove)) {
    .collapse_hidden(vertices, edges, stats::setNames(vertices %in% remove, vertices))
  } else { list(vertices = vertices, edges = data.table::copy(edges)) }
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


.collapse_hidden <- function( vertices, edges, hidden_vertex ){ # Connect visible nodes through hidden path interiors
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


## Truncate dotted node ids to the requested APSIM detail level: "supernode" keeps
## the segment before the first dot, "node" the first two segments, "subnode" the
## first three. Ids with fewer segments (or no dot) pass through unchanged. overrides is
## an optional named vector (supernode -> level) setting a per supernode level that
## replaces the global one for that first segment group (supernode_overrides).
## Returns a node_id -> grouped_id map (a character vector named by node_id).
.apsim_truncate <- function( ids, level, overrides = NULL ){ # Group node identifiers at the requested dotted detail level
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


## Contract an edge table through a node grouping: rewrite both endpoints to their
## group, drop edges that fall inside a single group (would be self-loops), and
## deduplicate. Directed contraction - it preserves reachability between groups and
## never produces a bidirected edge. Isolated groups simply gain no edges.
.contract_edges <- function( edges, group_map ){ # Contract edges through a node grouping and remove self-loops
  .datatable.aware <- TRUE
  if( nrow(edges) == 0L ){
    return(edges)
  }
  contracted_edges <- data.table::data.table(from = unname(group_map[ edges$from ]),
                                to   = unname(group_map[ edges$to ]))
  contracted_edges <- contracted_edges[ from != to ]
  return(unique(contracted_edges))
}


## Directed marginalisation ("latent projection", directed part only): remove every vertex in
## `drop` and reconnect each removed vertex's parents to its children, so a path pa -> m -> ch
## survives as pa -> ch. Applied iteratively, so a chain of removed vertices (m1 -> m2, both in
## `drop`) transitively reconnects pa(m1) -> ch(m2). The result is independent of order (confluent).
## Self-loops from the reconnection are dropped. NOTE this emits ONLY directed edges: it is faithful
## (preserves d-separation) exactly when no removed vertex is a common cause of two retained vertices
## that lack a retained common ancestor. In the treatment anchored unroll that condition holds for the
## interior occasions folded by pre_post -- every interior occasion descends from the retained boundary
## (_pre_treatment) occasion via the persistence / lag chain, so the confounding is represented at that
## retained boundary rather than being lost. The pre_post fold asserts acyclicity of the result and
## falls back to per slice labelling for any base where the fold would not break the cycle.
.marginalise <- function( edges, drop ){ # Reconnect parents and children through removed vertices
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


## Strongly connected components of a contracted graph, via mutual reachability (two
## labels share an SCC iff each reaches the other). The contracted graph carries no
## self-loops (.contract_edges drops them), so a nontrivial SCC has length >= 2 and is
## exactly a directed cycle of groups. Returns a list of label vectors (singletons too).
.scc_of <- function( labels, edges ){ # Find strongly connected components using an explicit traversal stack
  .datatable.aware <- TRUE
  if ( nrow(edges) == 0L ) {
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

  for ( root in seq_len(vertex_count) ) {
    if ( index[ root ] != 0L ) {
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

    while ( depth > 0L ) {
      vertex <- vertex_stack[ depth ]
      children <- children_by_parent[[ vertex ]]
      child_index <- child_stack[ depth ]
      if ( child_index <= length(children) ) {
        child_stack[ depth ] <- child_index + 1L
        child <- children[ child_index ]
        if ( index[ child ] == 0L ) {
          counter <- counter + 1L
          index[ child ] <- counter
          lowlink[ child ] <- counter
          stack_top <- stack_top + 1L
          stack[ stack_top ] <- child
          on_stack[ child ] <- TRUE
          depth <- depth + 1L
          vertex_stack[ depth ] <- child
          child_stack[ depth ] <- 1L
        } else if ( on_stack[ child ] ) {
          if ( index[ child ] < lowlink[ vertex ] ) {
            lowlink[ vertex ] <- index[ child ]
          }
        }
      } else {
        if ( lowlink[ vertex ] == index[ vertex ] ) {
          component_count <- component_count + 1L
          repeat {
            child <- stack[ stack_top ]
            stack_top <- stack_top - 1L
            on_stack[ child ] <- FALSE
            component[ child ] <- component_count
            if ( child == vertex ) {
              break
            }
          }
        }
        depth <- depth - 1L
        if ( depth > 0L ) {
          parent <- vertex_stack[ depth ]
          if ( lowlink[ vertex ] < lowlink[ parent ] ) {
            lowlink[ parent ] <- lowlink[ vertex ]
          }
        }
      }
    }
  }

  label_components <- component[ as.integer(vertex_index[ labels ]) ]
  components <- vector("list", 0L)
  emitted <- logical(component_count)
  for ( k in seq_along(labels) ) {
    component_id <- label_components[ k ]
    if ( emitted[ component_id ] ) {
      next
    }
    emitted[ component_id ] <- TRUE
    components[[ length(components) + 1L ]] <- labels[ label_components == component_id ]
  }
  return(components)
}


## Greedy refinement to the coarsest acyclic mixed granularity partition (force_acyclic).
## Starting from the contracted (cyclic) graph, repeatedly refine one supernode in a
## cycle by one detail level - preferring unpinned supernodes, splitting a pinned one
## only when a cycle has no other refinable member - and check again, until acyclic. The
## variable level is the original (acyclic) graph, so this always terminates. Returns a
## supernode_overrides vector (the per supernode levels that differ from the global one).
## This is a maximal reduction heuristic, not a provably coarsest partition.
.refine_to_acyclic <- function( node_ids, kept_ids, contemp, level, base_overrides, pinned_supers ){ # Refine grouped nodes until their contracted graph is acyclic
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
  ## the cyclic state at a given per supernode depth: number of vertices left in a cycle
  ## (0 = acyclic), with the contracted graph and the cyclic vertex set.
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
      break  ## unreachable while the finest graph is acyclic. Guards against a spin
    }
    ## only commit a split that actually shrinks the cyclic set, so an unpinned supernode
    ## that happens to sit in the cycle but does not break it is never wasted. Preference:
    ## an effective unpinned split, then an effective pinned one (a named supernode must
    ## split when it is the only thing breaking the cycle), then make progress otherwise.
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


## Longest common dotted segment prefix of a set of ids. The first id if there is no
## shared leading segment. Used to label a recombined supernode.
.dotted_lcp <- function( ids ){ # Find the longest common dotted prefix of node identifiers
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
    ## no shared leading segment: fall back to the first id. Two disjoint groups can hit this
    ## and collide, so .esc_recombine's #k disambiguation must run downstream to keep names unique
    return(ids[ 1 ])
  }
  return(paste(prefix, collapse = "."))
}


## ESC-DAGs recombination (Ferguson et al. 2019, Integration 2): merge vertices that
## share an identical set of parents AND an identical set of children into one node.
## This preserves acyclicity - two such vertices can never lie on a directed path
## through other vertices (that would need a child of the set to reach a parent of the
## set, impossible in a DAG) - so the contracted graph stays a DAG. Fully isolated
## vertices and any pinned vertices (named treatments/outcomes) are never merged. Each
## merged node is labelled by the longest common dotted prefix of its members (slice
## suffix reattached, made unique with a #k suffix). Latent carries over only when
## every member of a merged node was latent.
.esc_recombine <- function( vertices, edges, latent_set, pinned = character(0L) ){ # Merge unpinned vertices with identical parents and children
  parents_by_child <- split(edges$from, factor(edges$to,   levels = vertices))
  children_by_parent <- split(edges$to,   factor(edges$from, levels = vertices))
  signatures <- sapply(vertices, function( v ){ # Encode this vertex's parent and child sets without collisions
    parents <- if( is.null(parents_by_child[[ v ]]) ){ character(0L) }else{ sort(unique(parents_by_child[[ v ]])) }
    children <- if( is.null(children_by_parent[[ v ]]) ){ character(0L) }else{ sort(unique(children_by_parent[[ v ]])) }
    ## "\u0001"/"\u0002" cannot occur in a node id, and prefixing with the length keeps the parent
    ## and child blocks unambiguous, so structurally distinct vertices never collide here
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
    ## .dotted_lcp falls back to its first member when a group shares no dotted prefix, so
    ## two disjoint groups can propose the same name - this #k disambiguation is needed
    ## for label uniqueness and must not be removed
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


## Scope the recovered graph to nodes structurally relevant to the causal query.
## query_depth = 1 (historical compact option): the minimal readable scope that still answers the adjustment
## question. Keeps only the nodes on a directed X -> Y causal path, plus the DIRECT parents of X and
## of Y, and then FILLS every time slice occasion of each included base. The filling of occasions is not
## optional: a temporally unrolled graph with occasion gaps hides backdoor paths (e.g. NO3@t0 reaches the
## outcome through NO3@t1 -> NH4@t1 -> ... which, if @t1 is absent, makes an open backdoor invisible
## and dagitty then reports the empty set as sufficient -- wrong). Parents of X and Y are exactly the
## confounders (common parents) and the immediate mediators, so this scope recovers the structure a
## reader needs while cutting the ancestral tail that only inflates the node count.
##
## This scope is a fixed neighbourhood with query_depth = 1. query_depth >= 2
## switches to the ancestral "explore" scope, widened by generations. It can produce some minimal sets that are NOT valid in
## the full graph -- the scoped subgraph masks some backdoors -- but at least one valid set is always
## present (verified across the Fert/grain_n x GP/GY queries). Callers who need valid sets pipe
## the result through get_valid_sets(), which validates again against the full ancestral graph.
.scope_adjust <- function( recovered_vertices, recovered_edges, treat_v, out_v, pinned, quiet ){ # Retain the compact query scope and its ancestral occasions
  from <- as.character(recovered_edges$from)
  to <- as.character(recovered_edges$to)
  parents_by_child <- split(from, to)
  children_by_parent <- split(to, from)

  ## directed X -> Y paths: forward DFS from each treatment occasion, keeping only vertices that lie
  ## on a path reaching some outcome occasion. Two passes find nodes reachable forward from X
  ## and backward from Y. A node is on a causal path if and only if it is in both.
  forward_reach <- character(0L)
  stack <- treat_v
  while ( length(stack) > 0L ){
    vertex <- stack[[ length(stack) ]]
    stack <- stack[ -length(stack) ]
    if ( vertex %in% forward_reach ){ next }
    forward_reach <- c(forward_reach, vertex)
    stack <- c(stack, children_by_parent[[ vertex ]])
  }
  backward_reach <- character(0L)
  stack <- out_v
  while ( length(stack) > 0L ){
    vertex <- stack[[ length(stack) ]]
    stack <- stack[ -length(stack) ]
    if ( vertex %in% backward_reach ){ next }
    backward_reach <- c(backward_reach, vertex)
    stack <- c(stack, parents_by_child[[ vertex ]])
  }
  path_nodes <- intersect(forward_reach, backward_reach)                       # on at least one directed X -> Y path

  ## direct parents of X and of Y (confounders + immediate causes of the outcome)
  treatment_parents <- unique(unlist(parents_by_child[ treat_v ], use.names = FALSE))
  outcome_parents <- unique(unlist(parents_by_child[ out_v ],   use.names = FALSE))

  base_scope <- unique(c(path_nodes, treatment_parents, outcome_parents, treat_v, out_v))
  base_scope <- intersect(base_scope, recovered_vertices)

  ## ancestors of X or Y (inclusive), the set relevant to identification. A backward reach from the
  ## treatment and outcome occasions along parents.
  endpoint_ancestors <- character(0L)
  stack <- unique(c(treat_v, out_v))
  while ( length(stack) > 0L ){
    vertex <- stack[[ length(stack) ]]
    stack <- stack[ -length(stack) ]
    if ( vertex %in% endpoint_ancestors ){ next }
    endpoint_ancestors <- c(endpoint_ancestors, vertex)
    stack <- c(stack, parents_by_child[[ vertex ]])
  }

  ## fill every occasion of each included base, then keep only occasions that are themselves
  ## ancestors of X or Y (plus pinned). This preserves backdoors through intermediate occasions --
  ## those ARE ancestors of the outcome -- while dropping leaf occasions that are not ancestors, which are
  ## childless colliders that can never enter an adjustment set (dagitty already ignores them. They
  ## only clutter the recovered graph).
  bases <- unique(.strip_slice(base_scope))
  current <- recovered_vertices[ .strip_slice(recovered_vertices) %in% bases ]
  current <- intersect(current, endpoint_ancestors)
  current <- unique(c(current, intersect(pinned, recovered_vertices)))

  keep <- recovered_edges$from %in% current & recovered_edges$to %in% current
  if ( isFALSE(quiet) ){
    message(length(current), " of ", length(recovered_vertices),
            " nodes retained (adjust scope: causal paths + parents of treatment/outcome).")
  }
  return(list(vertices = current, edges = recovered_edges[ keep, , drop = FALSE ]))
}


## Base (graph_order = 0): ancestors of treatments and outcomes â€” the standard
## ancestral subgraph for covariate adjustment. Each additional generation adds
## the parents and children of the current boundary, widening the view.
## Pinned nodes (treatments, outcomes, present management levers) always survive.
.scope_to_query <- function( recovered_vertices, recovered_edges, treat_v, out_v, pinned, graph_order, quiet ){ # Retain endpoint ancestors and the requested extra generations
  parents_by_child <- split(as.character(recovered_edges$from), as.character(recovered_edges$to))
  children_by_parent <- split(as.character(recovered_edges$to),   as.character(recovered_edges$from))

  ## Base: An(X âˆª Y)
  current <- unique(c(treat_v, out_v))
  repeat {
    new_ancestors <- setdiff(unique(unlist(parents_by_child[ current ], use.names = FALSE)), current)
    new_ancestors <- intersect(new_ancestors, recovered_vertices)
    if ( length(new_ancestors) == 0L ){
      break
    }
    current <- c(current, new_ancestors)
  }

  ## Expand: each generation adds parents AND children of current boundary
  generations <- 0L
  if ( graph_order > 0L ){
    limit <- if ( is.infinite(graph_order) ){ length(recovered_vertices) } else { as.integer(graph_order) }
    for ( i in seq_len(limit) ){
      new_parents <- setdiff(unique(unlist(parents_by_child[ current ], use.names = FALSE)), current)
      new_children <- setdiff(unique(unlist(children_by_parent[ current ], use.names = FALSE)), current)
      new_nodes <- intersect(unique(c(new_parents, new_children)), recovered_vertices)
      if ( length(new_nodes) == 0L ){
        break
      }
      current <- c(current, new_nodes)
      generations <- generations + 1L
    }
  }

  ## Pinned nodes always survive
  current <- unique(c(current, intersect(pinned, recovered_vertices)))

  ## Trim edges to retained vertices
  keep <- recovered_edges$from %in% current & recovered_edges$to %in% current

  if ( isFALSE(quiet) ){
    message(length(current), " of ", length(recovered_vertices), " nodes retained",
            if ( generations > 0L ){ paste0(" (", generations, " generation(s) expanded)") } else { "" },
            ".")
  }

  return(list(vertices = current, edges = recovered_edges[ keep, ]))
}


## Reviewed public output paths, qualified by crop. Empty/missing metadata never
## certifies availability. This register is independent of recovery policy and
## the older apsim_obtainable module locators. No node layer changes are implied.
.query_apsim_outputs <- function(crop){ # Read the verified public output mappings for the requested crop
  path <- system.file("extdata", "apsim-output-availability.csv", package = "causaliflower")
  if (!nzchar(path)) { stop("The APSIM output availability register is missing. Reload the updated package.") }
  mappings <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  return(mappings[mappings$crop == crop & !is.na(mappings$source_expression) &
             nzchar(trimws(mappings$source_expression)) &
             mappings$status == "source_verified_public", , drop = FALSE])
}


## Reuse deliberate recovery exceptions rather than require 1,000 output mappings.
## Preserve original metadata and expose conflicting broad module locators.
.query_output_policy <- function(nodes, crop){ # Resolve recovery assumptions and conflicting output metadata
  ids <- nodes$node_id
  rule <- rep("yes", length(ids))
  if (identical(crop, "barley") && "recovered_node" %in% names(nodes)) {
    rule <- as.character(nodes$recovered_node)
  }
  if (anyNA(rule) || any(!rule %in% c("yes", "yes:user") & !grepl("^no:.+", rule))) {
    stop("Invalid recovered_node policy in the node table.")
  }
  non_apsim <- if ("provenance" %in% names(nodes)) { nodes$provenance %in% "not_apsim" } else { rep(FALSE, length(ids)) }
  rule[non_apsim & !startsWith(rule, "no:")] <- "yes:user"
  locator <- if ("apsim_obtainable" %in% names(nodes)) { as.character(nodes$apsim_obtainable) } else { rep("", length(ids)) }
  return(data.frame(node_id = ids, recovered_node = rule,
    reason = ifelse(non_apsim, "Provenance outside APSIM. The user must supply its interpretation/data",
      ifelse(rule == "yes:user", "Explicit quantity definition/request required",
        ifelse(startsWith(rule, "no:"), substring(rule, 4L), "Presumed obtainable. Not individually certified"))),
    apsim_locator = locator,
    metadata_conflict = rule != "yes" & !is.na(locator) & nzchar(trimws(locator)),
    stringsAsFactors = FALSE))
}


## Explicit columns take precedence over the policy based availability assumption.
## matched_ids is the resolved set of node_ids present in the columns.
##
## always_observed is EXEMPT from all of the above. A declared treatment or outcome is observed BY
## CONSTRUCTION: you know the dose you applied and you measured the thing you are explaining. Marking
## either [latent] makes the query self-contradictory ("estimate the effect of X on Y, given that I
## never observed X"), and dagitty then correctly reports NO adjustment set for it -- silently, and
## for a reason that has nothing to do with the causal structure. Before this, supplying
##   columns = c("grain protein")   with   treatments = "N_fertiliser_applied"
## emitted  "nitrogen_applied" [latent,exposure]  -- a latent exposure -- and every
## adjustment query on it returned zero sets no matter how well identified the effect actually was.
.assign_latent <- function( recovered_vertices, hidden_vertex, has_columns, matched_ids, show_hidden_nodes,
                            always_observed = character(0L), apsim_observed_ids = character(0L),
                            unavailable_ids = character(0L) ){ # Assign latent status from explicit columns or availability assumptions
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


## Emit dagitty syntax. Every vertex is declared (so isolated nodes survive) and
## quoted (node_ids carry dots and an @-slice suffix). Tags [latent] / [exposure]
## / [outcome] are attached where assigned.
.emit_dagitty <- function( vertices, edges, latent_set, exposures = character(0L), outcomes = character(0L) ){ # Emit every vertex, edge and causal role in dagitty syntax
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


## Directed reachability from src to dst over an edge table.
.reachable <- function( edges, src, dst ){ # Check directed reachability between the supplied endpoints
  if( nrow(edges) == 0L ){
    return(FALSE)
  }
  children_by_parent <- split(edges$to, edges$from)
  seen <- src
  stack <- unlist(children_by_parent[ src ], use.names = FALSE)
  while( length(stack) > 0L ){
    vertex <- stack[ length(stack) ]
    stack <- stack[ -length(stack) ]
    if( vertex %in% dst ){
      return(TRUE)
    }
    if( !(vertex %in% seen) ){
      seen <- c(seen, vertex)
      children <- children_by_parent[[ vertex ]]
      if( !is.null(children) ){
        stack <- c(stack, children)
      }
    }
  }
  return(any(dst %in% seen))
}


## Normalise either backend's adjustment set return into one common
## shape: a plain, unnamed list of character vectors, one per minimal adjustment
## set (the empty set is character(0)), and an empty list when no adjustment
## set exists. dagitty::adjustmentSets() returns a classed list of character
## vectors (empty list when none). minimal_sets() returns a named list and
## invisible(NULL) when none. Consumers then see one shape from either backend.
.normalise_adjustment_sets <- function( sets ){ # Normalise either adjustment backend to an unnamed list of sets
  if( is.null(sets) ){
    return(list())
  }
  sets <- unclass(sets)
  if( !is.list(sets) ){
    sets <- list(sets)
  }
  normalised_sets <- lapply(sets, function( s ){ # Represent one adjustment set as a character vector
    return(as.character(s))
  })
  names(normalised_sets) <- NULL
  return(normalised_sets)
}


## Adjustment set backend dispatcher. Default is dagitty::adjustmentSets. Set
## options(causaliflower.adjustment_backend = "minimal_sets") to route through
## causaliflower::minimal_sets instead (the alternative). Both backends' returns
## pass through .normalise_adjustment_sets so callers see one shape either way.
##
## WHY THIS IS BUDGETED. type = "minimal" ENUMERATES EVERY minimal adjustment set, and on a
## temporally unrolled DAG that enumeration does not terminate in practice. Unrolling splits a base into
## occasions for each slice, and every occasion of a base is an INTERCHANGEABLE blocker for the same
## backdoor path, so the number of distinct minimal sets multiplies across split bases. On the barley
## worked example (Fert -> grain_protein, selective, no columns) that is 40 split bases -- 10 with 3
## occasions, 30 with 2 -- for a crude bound of ~6.3e13 candidate sets. Measured: type = "minimal"
## had not returned after 4 minutes. type = "canonical" returned in 0.20 s. This is INHERENT to
## unrolling. Enumeration cost does not establish whether the graph is scientifically correct.
##
## Run "minimal" in a separate Rscript process with a time limit. V8 does not honour
## R's setTimeLimit, so the child process must be stopped when time runs out.
## Request a canonical adjustment set where supported. A canonical set is valid if any valid
## adjustment set exists. Identification by another strategy does not guarantee adjustment.
## Supplying `columns` states measurement availability and can change the search. In the
## historical worked example with the treatment's confounders declared, the query returns
##   1 minimal set: { soil_nh4_t0, soil_no3_t0 }
## in 0.20 s -- the before treatment occasion of the soil N that drives the fertiliser decision rule.
## Run an R expression in a subprocess with a time limit. Unlike parallel::mcparallel (which
## forks), this spawns a fresh Rscript process, so the V8 isolate is created from scratch instead of
## inheriting a broken copy from the parent's fork. The expression must evaluate to an R object, which
## is returned via a temp RDS file. Returns NULL on timeout, crash, or missing result.
.run_budgeted <- function( expr_lines, budget ){ # Run an adjustment calculation in a subprocess with a time budget
  result_file <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(result_file, script)), add = TRUE)
  writeLines(c(expr_lines, sprintf("saveRDS(.RESULT., %s)", deparse(result_file))), script)
  ## Use coreutils `timeout` for the time limit. Rscript runs in its own process group,
  ## which the timeout command stops along with V8.
  env_args <- character(0L)
  library_path <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
  if ( nzchar(library_path) ) { env_args <- c("env", paste0("LD_LIBRARY_PATH=", library_path)) }
  status <- system2("timeout",
                    args = c("--signal=KILL", as.character(ceiling(budget)),
                             env_args, "Rscript", script),
                    stdout = FALSE, stderr = FALSE)
  if ( status == 0L && file.exists(result_file) ){
    return(readRDS(result_file))
  }
  return(NULL)
}


.adjustment_sets_budgeted <- function( dag_text, exposure, outcome, effect, type ){ # Bound minimal set enumeration and retain the existing fallback
  budget <- getOption("causaliflower.adjustment_timeout", default = 30)
  direct <- function(){ # Run the requested adjustment calculation in the current process
    graph <- dagitty::dagitty(dag_text)
    return(dagitty::adjustmentSets(graph, exposure = exposure, outcome = outcome,
                                   effect = effect, type = type))
  }
  ## Only "minimal" needs an enumeration budget. Other types run in the current process.
  if ( !identical(type, "minimal") ){
    return(direct())
  }
  ## Spawn a subprocess with a time limit because V8's isolate is not safe to fork
  minimal_result <- .run_budgeted(c(
    "suppressMessages(library(dagitty))",
    sprintf("g <- dagitty(%s)", deparse(dag_text)),
    sprintf(".RESULT. <- adjustmentSets(g, exposure = %s, outcome = %s, effect = %s, type = \"minimal\")",
            deparse(exposure), deparse(outcome), deparse(effect))
  ), budget)
  if ( !is.null(minimal_result) ){
    return(minimal_result)
  }
  ## Timed out or crashed. Fall back.
  graph <- dagitty::dagitty(dag_text)
  ## Canonical adjustment sets exist only for TOTAL effects -- dagitty supports type = "minimal"
  ## exclusively for direct effects (adjustmentSets(effect = "direct", type = "canonical") throws
  ## "Only minimal adjustment sets are supported for direct effects!"). So the canonical fallback is
  ## available for total effects only. For a direct effect whose minimal enumeration times out there
  ## is no polynomial alternative, so we return an empty result with an actionable warning rather
  ## than crash.
  if ( identical(effect, "direct") ){
    warning(sprintf(paste0("Enumerating minimal adjustment sets for the DIRECT effect did not finish ",
                           "within %g s, and dagitty offers no canonical fallback for direct effects. ",
                           "No adjustment set is returned. Narrow the graph (query_depth = 1) or ",
                           "declare observed variables with columns = to shrink the search, or raise the ",
                           "budget with options(causaliflower.adjustment_timeout = <seconds>)."),
                   budget),
            call. = FALSE)
    return(structure(list(), class = "dagitty.sets"))
  }
  canonical_result <- dagitty::adjustmentSets(graph, exposure = exposure, outcome = outcome,
                                 effect = effect, type = "canonical")
  warning(sprintf(paste0("Enumerating every MINIMAL adjustment set did not finish within %g s, so the ",
                         "CANONICAL adjustment set is returned instead (%s). Time unrolling splits a ",
                         "variable into occasions for each slice, and each occasion blocks the same backdoor ",
                         "path, so the number of minimal sets grows combinatorially. To get small, ",
                         "usable sets back, pass columns = <the variables you actually have>: that ",
                         "restricts the search to observed variables, which makes it both tractable and ",
                         "meaningful. Raise the budget with ",
                         "options(causaliflower.adjustment_timeout = <seconds>)."),
                 budget,
                 if ( length(canonical_result) ) { sprintf("%d variables", length(canonical_result[[ 1L ]])) } else { "none exists" }),
          call. = FALSE)
  return(canonical_result)
}


.adjustment_backend <- function( dag, exposure, outcome, effect = "total", type = "minimal" ){ # Dispatch to the selected adjustment backend and normalise its result
  backend <- getOption("causaliflower.adjustment_backend", default = "dagitty")
  if( !requireNamespace("dagitty", quietly = TRUE) ){
    stop("Computing adjustment sets needs the dagitty package. Please install dagitty and try again.")
  }
  ## Both backends operate on a real dagitty object, but the carrier passed in here holds only
  ## the graph text (no structure). Rebuild the object once and dispatch on it. The minimal_sets
  ## branch previously passed the structureless carrier straight to minimal_sets(), which made
  ## dagitty::adjustmentSets() fail internally and silently report "no sets".
  dag_text <- attr(dag, "dag_text")
  if( identical(backend, "minimal_sets") && exists("minimal_sets", mode = "function") ){
    ## minimal_sets() is a thin wrapper that calls adjustmentSets(type = "minimal") and then keeps the
    ## n smallest, so it has the same enumeration cost and needs the same budget. Its formals
    ## are `treatments` / `outcomes` (plural). Pass them exactly, not on R's partial matching.
    raw <- .with_budget_wrapper(dag_text, exposure, outcome, effect)
    return(.normalise_adjustment_sets(raw))
  }
  raw <- .adjustment_sets_budgeted(dag_text, exposure, outcome, effect, type)
  return(.normalise_adjustment_sets(raw))
}


## Apply the dagitty backend's time limit and canonical fallback to minimal_sets()
.with_budget_wrapper <- function( dag_text, exposure, outcome, effect ){ # Apply the same time budget to the minimal_sets wrapper
  budget <- getOption("causaliflower.adjustment_timeout", default = 30)
  ## minimal_sets() is called in a subprocess (not a fork) with a time limit, exactly like
  ## .adjustment_sets_budgeted, because V8's isolate is not safe to fork.
  minimal_result <- .run_budgeted(c(
    "suppressMessages(library(dagitty))",
    sprintf("g <- dagitty(%s)", deparse(dag_text)),
    sprintf(".RESULT. <- minimal_sets(g, treatments = %s, outcomes = %s, effect = %s)",
            deparse(exposure), deparse(outcome), deparse(effect))
  ), budget)
  if ( !is.null(minimal_result) ){
    return(minimal_result)
  }
  return(.adjustment_sets_budgeted(dag_text, exposure, outcome, effect, type = "canonical"))
}


## Sanitise a supplied label so it cannot corrupt the emitted dagitty text. .emit_dagitty()
## quotes every vertex with double quotes, so a literal double quote or a newline in a label would
## break the model string. Strip both (plus a carriage return) and trim. An empty result signals the
## caller to fall back to the node_id. (label= feature, guard b)
.sanitise_label <- function( x ){ # Remove characters that would corrupt a quoted graph label
  return(trimws(gsub("[\"\r\n]", "", x)))
}


## Canonical quantity membership of every emitted vertex, before naming changes.
## This ledger is metadata only: it never changes the graph's nodes or arrows.
.query_identity_rows <- function(vertices, group_map, kept_ids){ # Record each emitted vertex's canonical quantity and members
  grouped_members <- split(kept_ids, unname(group_map[ kept_ids ]))
  rows <- lapply(vertices, function(vertex){ # Describe the canonical members of one emitted vertex
    quantity <- .strip_slice(vertex)
    members <- grouped_members[[ quantity ]]
    if ( length(members) == 0L ){ members <- quantity }
    suffix <- substring(vertex, nchar(quantity) + 1L)
    return(data.frame(node = vertex, quantity = quantity, member_id = members,
               occurrence = vertex, member_occurrence = paste0(members, suffix),
               stringsAsFactors = FALSE))
  })
  if ( length(rows) == 0L ){
    return(data.frame(node = character(), quantity = character(), member_id = character(),
                      occurrence = character(), member_occurrence = character()))
  }
  return(do.call(rbind, rows))
}


## Apply input aliases as node labels to a recovered graph (the label= feature). For each distinct base
## (the suffix stripped vertex, which equals a node_id) a base label is chosen by mode:
##   "input"          -> node_id                            (an input override may rename any node)
##   "APSIM"          -> node_id                            (an input override renames only
##                                                           undotted keys)
##   "readable" -> node_name_readable, falling back to node_id where the human label is empty
## The legacy APSIM label mode uses has_dot to protect dotted keys from input overrides.
## An input override from input_pairs (node_id -> input string) replaces the base label wherever
## override_allowed = !( label == "APSIM" && has_dot ), after sanitising (guard b). An empty
## sanitised label is ignored. A label claimed by more than one base reverts those bases to the
## node_id (two passes, mirroring readable_dag) so the relabelling stays a bijection. The temporal
## suffix is then reattached and every vertex / edge / latent / exposure / outcome remapped. Returns
## the relabelled sets plus input_labels (displayed label -> node_id) for the nodes an input renamed.
.apply_input_labels <- function( V, E, lat, X, Y, nodes, input_pairs, label ){ # Apply display labels while preserving a unique mapping to node identities
  .datatable.aware <- TRUE
  base <- .strip_slice(V)
  suffix <- substring(V, nchar(base) + 1L)
  unique_bases <- unique(base)
  has_dot <- grepl(".", unique_bases, fixed = TRUE)

  ## label for each distinct base by mode
  if ( identical(label, "readable") ){
    human_labels <- stats::setNames(as.character(nodes$node_name_readable), nodes$node_id)
    candidate_labels <- unname(human_labels[ unique_bases ])
    missing_labels <- is.na(candidate_labels) | !nzchar(candidate_labels)
    candidate_labels[ missing_labels ] <- unique_bases[ missing_labels ]
  }else{
    candidate_labels <- unique_bases                                             # "APSIM" and "input": base = node_id
  }
  candidate_labels <- stats::setNames(candidate_labels, unique_bases)

  ## input override, per base, where allowed. Record which bases an input actually renamed
  override_allowed <- stats::setNames(if ( identical(label, "APSIM") ){ !has_dot }else{ rep(TRUE, length(unique_bases)) }, unique_bases)
  renamed <- character(0L)                                 # node_id -> displayed label
  for ( node_id in intersect(names(input_pairs), unique_bases) ){
    if ( isTRUE(override_allowed[[ node_id ]]) ){
      input_label <- .sanitise_label(input_pairs[[ node_id ]])
      if ( nzchar(input_label) ){
        candidate_labels[[ node_id ]] <- input_label
        renamed[[ node_id ]] <- input_label
      }
    }
  }

  ## uniqueness guard: any label claimed by more than one base reverts to the node_id (two passes,
  ## covering the rare case where a chosen label coincides with another base's node_id)
  for ( pass in 1:2 ){
    duplicate_labels <- candidate_labels[ duplicated(candidate_labels) | duplicated(candidate_labels, fromLast = TRUE) ]
    if ( length(duplicate_labels) == 0L ){
      break
    }
    clash <- candidate_labels %in% duplicate_labels
    candidate_labels[ clash ] <- names(candidate_labels)[ clash ]
  }
  ## a base that the guard reverted no longer holds its input label
  if ( length(renamed) > 0L ){
    still <- vapply(names(renamed), function( node_id ){ # Check whether the uniqueness guard retained an input label
      return(identical(unname(candidate_labels[[ node_id ]]), unname(renamed[[ node_id ]])))
    }, logical(1L))
    renamed <- renamed[ still ]
  }

  relabel <- stats::setNames(paste0(unname(candidate_labels[ base ]), suffix), V)
  new_vertices <- unname(relabel[ V ])
  new_edges <- if ( nrow(E) > 0L ){
    data.table::data.table(from = unname(relabel[ E$from ]), to = unname(relabel[ E$to ]))
  }else{
    E
  }
  new_latent <- unname(relabel[ intersect(lat, V) ])
  new_exposures <- unname(relabel[ intersect(X, V) ])
  new_outcomes <- unname(relabel[ intersect(Y, V) ])

  input_labels <- if ( length(renamed) > 0L ){ stats::setNames(names(renamed), unname(renamed)) }else{ character(0L) }
  return(list(V = new_vertices, E = new_edges, lat = new_latent, X = new_exposures, Y = new_outcomes, input_labels = input_labels))
}
