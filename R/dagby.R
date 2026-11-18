# ==============================================================================
# dagby.R


# Query recipes and construction
# ==============================================================================


#' Describe a DAG query
#'
#' Create a recipe for a crop question without running the query. Save it as an
#' object, then use `recipe |> dagby_query() |> dagby_construct()` to build a DAG.
#' The recipe holds both the query inputs and preferences for construction.
#'
#' @param crop Crop whose process-based model supplies the graph, such as
#'   `"barley"` or `"wheat"`.
#' @param treatments,outcomes Variables whose causal relationship you want to
#'   examine. Supply model identifiers or recognised aliases. Named vectors such
#'   as `c(N = "nitrogen_applied")` also supply display names for `label = "input"`
#'   or `"readable"`. Explicit `NULL` omits that endpoint for source exploration.
#' @param randomised_treatment Whether the treatment was randomised. Supply
#'   `TRUE` or `FALSE` explicitly. `TRUE` removes incoming arrows to treatment
#'   occasions according to the temporal representation in [dagby_query()].
#' @param size Preferred number of source variables in the constructed DAG.
#'   `"smallest"` (default) requests the smallest available DAG and `"largest"`
#'   the largest. A number requests the closest available size, choosing the
#'   larger on a tie. `0` and `Inf` are equivalent to `"smallest"` and `"largest"`.
#'   See [dagby_construct()] for counting and retention rules.
#' @param keep_nodes Additional variables or exact time occasions to retain
#'   during construction, even if they are not needed by its selection rules.
#' @param latents Variables or exact occasions to mark as unmeasured, as in
#'   [dagby_construct()]. They are marked after construction and do not affect
#'   which variables are kept. `NULL` treats all variables as observed.
#' @param columns Observed variables in your data, used to tag the source query.
#'   `NULL` uses the automatic observation assumptions in [dagby_query()]. These
#'   source tags remain saved. Use `latents` to declare unmeasured variables on
#'   the constructed DAG.
#' @param time_unroll Temporal representation of the source query. If omitted,
#'   lagged relationships retain separate time occasions when treatments are
#'   supplied. Explicit `NULL` or `"none"` combines time occasions in the query
#'   and can introduce cycles.
#'   `"all"` retains separate occasions. Specialist modes `"selective"`,
#'   `"collapse"` and `"pre_post"` are described in [dagby_query()].
#' @param label Names to display on queried nodes. `"id"` (default) keeps model
#'   identifiers such as `nitrogen_applied`. `"readable"` uses descriptive model
#'   names, with your supplied names taking precedence. `"input"` uses the names
#'   or aliases supplied in `columns`, `treatments` and `outcomes`. `"APSIM"` keeps
#'   dotted source identifiers where present and otherwise follows `"input"`.
#'   The source identities remain saved.
#' @param unroll Whether construction should retain separate time occasions.
#'   `FALSE` requests the existing checked static view. `TRUE` retains occasions.
#'   This controls construction's output, independently of the source query's
#'   `time_unroll` setting.
#' @param params Specialist query settings created with [dagby_query_params()],
#'   including query depth, growth-stage endpoints and alias disambiguation.
#'
#' @return A `dagby_recipe` containing query inputs and construction preferences.
#'   Creating or printing it does not query the model or construct a DAG.
#'   Source names and graph-dependent assumptions are checked when the relevant
#'   stage runs.
#' @examples
#' recipe <- dagby(
#'   crop = "barley",
#'   treatments = "nitrogen_applied",
#'   outcomes = "grain_protein",
#'   randomised_treatment = FALSE
#' )
#' recipe
#' \dontrun{
#' dag <- recipe |> dagby_query() |> dagby_construct()
#' }
#' @export
dagby <- function(crop,
                  treatments,
                  outcomes,
                  randomised_treatment,
                  size = "smallest",
                  keep_nodes = NULL,
                  latents = NULL,
                  columns = NULL,
                  time_unroll = NULL,
                  label = "id",
                  unroll = FALSE,
                  params = dagby_query_params()
                  ){
  required <- c("crop", "treatments", "outcomes", "randomised_treatment")
  absent <- setdiff(required, names(match.call()))
  if( length(absent) ){
    stop("dagby() requires: ", paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  if( !is.character(label) || length(label) != 1L || is.na(label) ||
    !label %in% c("id", "readable", "input", "APSIM") ){
    stop('label must be "id", "readable", "input" or "APSIM".', call. = FALSE)
  }
  query <- list(crop = crop, columns = columns, treatments = treatments,
    outcomes = outcomes, randomised_treatment = randomised_treatment,
    label = if( label == "id" ) NULL else label, params = params)
  if( !missing(time_unroll) ){
    query["time_unroll"] <- list(time_unroll)
  }
  recipe <- structure(list(query = query, construction = list(size = size,
    latents = latents, keep_nodes = keep_nodes, unroll = unroll)), class = "dagby_recipe")

  return(.dagby_validate_recipe(recipe))
}


#' Construct a DAG from a query
#'
#' Builds a DAG of the requested size from the supplied graph. It retains
#' variables used by selected adjustment sets and construction rules, reducing
#' the graph while checking the retained adjustment sets against the source.
#'
#' @param graph A `dagby_query()` result, a dagitty DAG or a constructed DAG.
#'   A query result supplies its recipe's construction preferences for omitted
#'   arguments. Explicit arguments override those preferences.
#'   A constructed DAG is rebuilt from its saved query, with its saved
#'   `latents` and `keep_nodes` unless they are supplied.
#' @param size The requested number of source variables, counted as described
#'   below. `"smallest"` (default) and `"largest"` request the limits. A number
#'   requests the closest available size, choosing the larger on a tie.
#'   `0` and `Inf` also request the smallest and largest. An omitted size inherits
#'   the recipe's preference when `graph` is a query result.
#' @param latents Variables or exact occasions to mark as unmeasured. These
#'   markings do not change which nodes are kept. `NULL` treats all variables as
#'   observed. When omitted for a plain DAG, its latent tags are used.
#' @param keep_nodes Additional quantities or exact occasions to retain.
#' @param dry_run Compute a preview containing the proposed graph, construction
#'   plan and checks. This runs construction, unlike creating a recipe.
#' @param quiet `NULL` prints progress and how to change the size. `TRUE`
#'   suppresses them. `FALSE` also prints the construction record.
#' @param unroll `FALSE` requests a checked static view. `TRUE` keeps separate
#'   time occasions. If omitted for a constructed DAG, its current view is kept.
#'
#' @details
#' Size counts each source variable once, even when it appears at several
#' occasions or in several display segments. Variables combined in one displayed
#' node still count separately. For a plain DAG without recorded identities,
#' each exact node name counts separately.
#'
#' Smaller candidates keep variables used by ranked minimal adjustment sets for
#' the direct and total effects. Ranking certificates bound the sizes offered.
#' The largest candidates also retain variables left undecided by the step budget.
#'
#' Larger candidates add paths from treatment to outcome, other treatment
#' descendants and other parents of treatment, outcome and nodes on those paths.
#' The largest also includes relevant upstream variables with at least one
#' occasion that has no parent from a different variable. A parent at another
#' occasion does not disqualify the variable. The query and retention rules still
#' limit which occasions are kept.
#'
#' Every DAG keeps the endpoints and requested nodes. When no total adjustment
#' set is found, exact treatment parents are also kept where the identification
#' safeguard applies. When no direct set is found, the first node of each indirect
#' path is kept. Reaching the step budget does not establish that no set exists.
#'
#' A static view must pass the source checks. If it fails, construction explains
#' why. Use `unroll = TRUE` to request the DAG with separate time occasions.
#' `dagby_prune()` and `dagby_grow()` change the size. `dagby_select()` chooses an
#' alternative at the current or requested size, listed by `dagby_options()`.
#'
#' `options(causaliflower.ranked_steps = n)` sets the ranked search budget for
#' smaller sizes. The default is 20000 steps.
#'
#' @return A constructed `cf_dag`, or a `cf_dag_plan` when `dry_run = TRUE`.
#'   Both retain the query and checks.
#' @export
dagby_construct <- function(graph,
                          size = "smallest",
                          latents = NULL,
                          keep_nodes = NULL,
                          dry_run = FALSE,
                          quiet = NULL,
                          unroll = FALSE
                          ){
  input_record <- .cd_input_record(graph)
  if( inherits(graph, "cf_dag_plan") ){
    graph <- graph$graph
  }
  if( !inherits(graph, "cf_dag") && !is.null(input_record$recipe) ){
    preferences <- .dagby_validate_recipe(input_record$recipe)$construction
    if( missing(size) ){
      size <- preferences$size
    }
    if( missing(latents) ){
      latents <- preferences$latents
    }
    if( missing(keep_nodes) ){
      keep_nodes <- preferences$keep_nodes
    }
    if( missing(unroll) ){
      unroll <- preferences$unroll
    }
  }
  if( inherits(graph, "cf_dag") ){
    saved <- dagby_options(graph)
    if( missing(latents) ){
      latents <- saved$arguments$latents
    }
    if( missing(keep_nodes) ){
      keep_nodes <- saved$arguments$keep_nodes
    }
    if( missing(unroll) ){
      unroll <- .cd_output_mode(graph)
    }
    graph <- dagby_source(graph)
  }else if( missing(latents) && inherits(graph, "dagitty") &&
    !inherits(graph, "dag_by_query") && is.null(attr(graph, "query_record", exact = TRUE)) ){
    latents <- as.character(dagitty::latents(graph))
    if( !length(latents) ){
      latents <- NULL
    }
  }
  size <- .dagby_size(size)
  dag <- .cd_sized(graph, latents, keep_nodes, dry_run, quiet, size, unroll = unroll)

  return(.cd_with_input_record(dag, input_record))
}


## Recipe validation ##
#' Validate a recipe without querying its source model
#'
#' @param recipe A `dagby_recipe` containing query inputs and construction preferences.
#' @return The recipe with validated query parameters. Invalid inputs raise an error.
#' @noRd
.dagby_validate_recipe <- function(recipe){
  if( !inherits(recipe, "dagby_recipe") || !is.list(recipe) ||
    !setequal(names(recipe), c("query", "construction")) ){
    stop("Supply a recipe created by dagby().", call. = FALSE)
  }
  .gd_named_inputs(recipe, c("query", "construction"), "recipe")
  required <- c("crop", "columns", "treatments", "outcomes",
    "randomised_treatment", "label", "params")
  query <- .gd_named_inputs(recipe$query, c(required, "time_unroll"), "recipe query")
  preferences <- .gd_named_inputs(recipe$construction,
    c("size", "latents", "keep_nodes", "unroll"), "construction preferences")
  if( !all(required %in% names(query)) ||
    !all(c("size", "latents", "keep_nodes", "unroll") %in% names(preferences)) ){
    stop("The recipe is incomplete. Create it with dagby().", call. = FALSE)
  }
  for( name in c("crop", "treatments", "outcomes", "randomised_treatment") ){
    .gd_query_value(name, query[[name]])
  }
  for( name in c("columns", "latents", "keep_nodes") ){
    value <- if( name == "columns" ) query$columns else preferences[[name]]
    if( !is.null(value) && (!is.character(value) || anyNA(value) ||
      any(!nzchar(trimws(value)))) ){
      stop(name, " must contain nonempty node names, or be NULL.", call. = FALSE)
    }
  }
  if( !is.null(query$label) && (!is.character(query$label) ||
    length(query$label) != 1L || is.na(query$label) ||
    !query$label %in% c("readable", "input", "APSIM")) ){
    stop("The recipe has an invalid label setting. Create it with dagby().", call. = FALSE)
  }
  time <- query$time_unroll
  if( !is.null(time) ){
    if( is.logical(time) ){
      .gd_flag(time, "time_unroll")
    }else if( !is.character(time) || length(time) != 1L || is.na(time) ||
      !time %in% c("none", "all", "selective", "collapse", "pre_post") ){
      stop('time_unroll must be "none", "all", "selective", "collapse", "pre_post", TRUE, FALSE or NULL.',
        call. = FALSE)
    }
  }
  params <- .gd_named_inputs(query$params, names(formals(dagby_query_params)), "params")
  recipe$query$params <- do.call(dagby_query_params, params)
  .dagby_size(preferences$size)
  .gd_flag(preferences$unroll, "unroll")

  return(recipe)
}


#' Translate a public size preference for the construction engine
#'
#' @param size "smallest", "largest" or a nonnegative numeric size.
#' @return A numeric size, with zero and infinity representing the limits.
#' @noRd
.dagby_size <- function(size){
  if( is.character(size) ){
    if( length(size) != 1L || is.na(size) || !size %in% c("smallest", "largest") ){
      stop('size must be "smallest", "largest" or a number of at least 0.', call. = FALSE)
    }
    return(if( size == "smallest" ) 0 else Inf)
  }
  return(.cd_size_value(size))
}


#' Create a throttled progress reporter
#'
#' Reports elapsed time without estimating when separator enumeration will finish.
#' Progress messages are not saved in graph records.
#' @param quiet `TRUE` to suppress routine progress.
#' @param operation Label preceding the elapsed time and stage message.
#' @return `NULL` if `quiet` is `TRUE`. Otherwise, a function accepting `stage`
#'   and `force`. Forced calls report immediately. Other calls report at most
#'   once every two seconds.
#' @noRd
.gd_progress <- function(quiet, operation = "Construction"){
  if( isTRUE(quiet) ){
    return(NULL)
  }
  started <- proc.time()[["elapsed"]]
  last <- started
  progress <- function(stage, force = TRUE){
    now <- proc.time()[["elapsed"]]
    if( !force && now - last < 2 ){
      return(invisible(NULL))
    }
    message(operation, " [", format(round(now - started, 1), nsmall = 1), "s]: ", stage)
    flush.console()
    last <<- now
    return(invisible(NULL))
  }

  return(progress)
}


#' Report a construction stage when progress is enabled
#'
#' @param progress A function from `.gd_progress()`, or `NULL`.
#' @param stage Text describing the current operation.
#' @param force Whether to bypass the reporter's time throttle.
#' @return `NULL`, invisibly.
#' @noRd
.gd_notify <- function(progress, stage, force = TRUE){
  if( !is.null(progress) ){
    progress(stage, force = force)
  }
  return(invisible(NULL))
}


#' Validate named settings
#'
#' @param x A list, empty or with unique nonempty names.
#' @param allowed Argument names accepted by the underlying function.
#' @param name Input name used in errors.
#' @return `x`, unchanged. Invalid names or unsupported inputs raise an error.
#' @noRd
.gd_named_inputs <- function(x, allowed, name){
  if( !is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
    any(!nzchar(names(x))) || anyDuplicated(names(x)))) ){
    stop(name, " must be a named list with unique names.", call. = FALSE)
  }
  unknown <- setdiff(names(x), allowed)
  if( length(unknown) ){
    stop("Unsupported ", name, " input: ",
      paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  return(x)
}


#' Validate a required query input
#'
#' @param name One of `crop`, `treatments`, `outcomes` or `randomised_treatment`.
#' @param value The supplied character names or scalar logical design assumption.
#' @return `value`, unchanged. Missing or invalid values raise an error.
#' @noRd
.gd_query_value <- function(name, value){
  if( name %in% c("treatments", "outcomes") && is.null(value) ){
    return(NULL)
  }
  if( name == "randomised_treatment" ){
    if( !is.logical(value) || length(value) != 1L || is.na(value) ){
      stop("randomised_treatment must be TRUE or FALSE.", call. = FALSE)
    }
  }else if( !is.character(value) || !length(value) || anyNA(value) ||
    any(!nzchar(trimws(value))) || (name == "crop" && length(value) != 1L) ){
    stop(name, " must contain ", if( name == "crop" ){
      "one nonempty crop name."
    }else{
      "nonempty node names."
    }, call. = FALSE)
  }
  return(value)
}


#' Validate a scalar logical option
#'
#' @param x Value to validate.
#' @param name Option name used in errors.
#' @param nullable `TRUE` to accept `NULL` as well as a scalar logical value.
#' @return `x`, invisibly. Invalid values raise an error.
#' @noRd
.gd_flag <- function(x, name, nullable = FALSE){
  if( nullable && is.null(x) ){
    return(invisible(x))
  }
  if( !is.logical(x) || length(x) != 1L || is.na(x) ){
    stop(name, " must be TRUE or FALSE", if( nullable ){
      " or NULL"
    }else{
      ""
    }, ".", call. = FALSE)
  }
  return(invisible(x))
}


#' Read the supplied graph's query record
#'
#' @param graph The emitted query graph used for construction.
#' @return Its exact `query_record` attribute, or `NULL` when absent.
#' @noRd
.cd_record <- function(graph){
  return(attr(graph, "query_record", exact = TRUE))
}


#' Replace latent tags while retaining source identities
#'
#' @param graph A dagitty graph.
#' @param nodes Exact graph node names to mark as latent.
#' @return The graph with latent tags replaced and its input labels, node map
#'   and recorded source order preserved.
#' @noRd
.cd_set_latents <- function(graph, nodes){
  labels <- attr(graph, "input_labels", exact = TRUE)
  map <- attr(graph, "cf_node_map", exact = TRUE)
  source_order <- attr(graph, "cf_source_order", exact = TRUE)
  if( is.null(source_order) ){
    source_order <- .cd_record(graph)$source_order
  }
  dagitty::latents(graph) <- nodes
  attr(graph, "input_labels") <- labels
  attr(graph, "cf_node_map") <- map
  attr(graph, "cf_source_order") <- source_order
  return(graph)
}


#' Mark declared unmeasured occasions after construction
#'
#' A display node is unmeasured if any occasion it represents is declared.
#' Mixed observed and unmeasured membership produces a warning.
#' @param summary A `.gd_summary()` result with display and temporal graphs
#'   and a node map relating their exact occasions.
#' @param latent_nodes Exact source occasions declared unmeasured.
#' @return A list with the marked `summary`, display and temporal latent names,
#'   and `partial` membership split into observed and unmeasured occasions.
#' @noRd
.cd_mark_unmeasured <- function(summary, latent_nodes){
  map <- summary$node_map
  nodes <- if( is.data.frame(map) ){
    as.character(map$node)
  }else{
    names(summary$graph)
  }
  display <- if( is.data.frame(map) && "display" %in% names(map) ){
    as.character(map$display)
  }else{
    nodes
  }
  display_latent <- unique(display[nodes %in% latent_nodes])
  temporal_latent <- intersect(names(summary$temporal), latent_nodes)
  partial <- list()
  for( node in display_latent ){
    occasions <- unique(nodes[display == node])
    declared <- intersect(occasions, latent_nodes)
    if( length(declared) < length(occasions) ){
      partial[[node]] <- list(unmeasured = declared, observed = setdiff(occasions, declared))
    }
  }
  if( length(partial) ){
    warning(paste(vapply(names(partial), function(node){
      return(paste0(node,
        " stands for both observed and unmeasured occasions (unmeasured: ",
        paste(partial[[node]]$unmeasured, collapse = ", "), ", observed: ",
        paste(partial[[node]]$observed, collapse = ", "), ")"))
    }, character(1)), collapse = ". "),
      if( length(partial) == 1L ){
        ". It is shown unmeasured."
      }else{
        ". They are shown unmeasured."
      }, call. = FALSE)
  }
  if( length(display_latent) ){
    summary$graph <- .cd_set_latents(summary$graph, display_latent)
  }
  if( length(temporal_latent) ){
    summary$temporal <- .cd_set_latents(summary$temporal, temporal_latent)
  }
  shown <- list(summary = summary, display = display_latent, temporal = temporal_latent, partial = partial)

  return(shown)
}


#' Record unmeasured parents in the identification safeguard
#'
#' @param identification Retention record from `.cd_identification_retention()`.
#' @param latent_nodes Exact source occasions declared unmeasured.
#' @return The record, with an uncertified status and `unmeasured_parents` when
#'   an applied safeguard retains declared unmeasured parents. Retention is unchanged.
#' @noRd
.cd_unmeasured_parents <- function(identification, latent_nodes){
  unmeasured <- intersect(identification$parents, latent_nodes)
  if( !isTRUE(identification$applied) || !length(unmeasured) ){
    return(identification)
  }
  identification$unmeasured_parents <- unmeasured
  identification$status <- paste0("not certified: treatment parent(s) ", paste(unmeasured, collapse = ", "),
    " are declared unmeasured. Exact treatment parents and temporal occurrences are retained to show what must be observed")
  return(identification)
}


#' Read and validate occurrence identities
#'
#' @param graph A query graph whose current nodes must all be represented.
#' @param record Its saved query record, checked before the graph attribute.
#' @return A data frame with `node`, `quantity`, `member_id` and `occurrence`
#'   columns. Without a saved map, each exact name supplies all four identities.
#'   An invalid map raises an error.
#' @noRd
.cd_node_map <- function(graph, record = .cd_record(graph)){
  map <- record$node_map
  if( is.null(map) ){
    map <- attr(graph, "cf_node_map", exact = TRUE)
  }
  if( is.data.frame(map) && all(c("node", "quantity", "member_id", "occurrence") %in% names(map)) ){
    map <- map[map$node %in% names(graph), , drop = FALSE]
    if( !setequal(unique(map$node), names(graph)) || anyNA(map[c("node", "quantity", "member_id", "occurrence")]) ||
      any(!nzchar(as.character(map$member_id))) || any(!nzchar(as.character(map$quantity))) ){
      stop("The query receipt's node map does not match the supplied graph.", call. = FALSE)
    }
    return(map)
  }
  if( !is.null(map) ){
    stop("The query's source identity map is incomplete or invalid.", call. = FALSE)
  }
  map <- data.frame(node = names(graph), quantity = names(graph), member_id = names(graph),
    occurrence = names(graph), stringsAsFactors = FALSE)

  return(map)
}


#' Map graph nodes to source variables for size counting
#'
#' @param dag A graph with a node map attribute or saved query record.
#' @return A data frame of unique node and variable pairs, sorted in radix order
#'   when a map is present.
#'   Several variables may belong to one node. Without a map, each exact node
#'   name is its own variable and graph node order is retained. An incomplete
#'   or invalid map raises an error.
#' @noRd
.cd_variable_map <- function(dag){
  map <- attr(dag, "cf_node_map", exact = TRUE)
  if( is.null(map) ){
    map <- .cd_record(dag)$node_map
  }
  if( is.null(map) ){
    return(data.frame(node = names(dag), variable = names(dag), stringsAsFactors = FALSE))
  }
  if( !is.data.frame(map) || !all(c("node", "member_id") %in% names(map)) ){
    stop("The graph's source identity map is incomplete or invalid.", call. = FALSE)
  }
  map <- map[map$node %in% names(dag), c("node", "member_id"), drop = FALSE]
  if( !setequal(map$node, names(dag)) || anyNA(map) || any(!nzchar(as.character(map$member_id))) ){
    stop("The source identity map has missing nodes or invalid values.", call. = FALSE)
  }
  mapped <- unique(data.frame(node = as.character(map$node), variable = as.character(map$member_id),
    stringsAsFactors = FALSE))
  mapped <- mapped[order(mapped$node, mapped$variable, method = "radix"), , drop = FALSE]
  rownames(mapped) <- NULL
  return(mapped)
}


#' Count distinct source variables represented by nodes
#'
#' @param dag Graph providing the variable map.
#' @param nodes Exact node names whose variables are counted.
#' @return The number of distinct variables. Unknown nodes raise an error.
#' @noRd
.cd_variable_count <- function(dag, nodes = names(dag)){
  map <- .cd_variable_map(dag)
  if( !all(nodes %in% map$node) ){
    stop("A retained node has no recorded variable identity.", call. = FALSE)
  }
  return(length(unique(map$variable[map$node %in% nodes])))
}


#' Attach source identities to a canonical working graph
#'
#' @param graph Graph to canonicalise without altering its node identities.
#' @param map Occurrence map with one canonical quantity per displayed node.
#' @param source_order Optional recorded ordering for source quantities.
#' @return The canonical graph with reverse labels, node map and source order.
#'   A node assigned to several canonical quantities raises an error.
#' @noRd
.cd_working_graph <- function(graph, map, source_order = NULL){
  result <- .gd_canonical(graph)
  if( !is.null(source_order) ){
    attr(result, "cf_source_order") <- source_order
  }
  # Complete reverse labels, including human labels, are analysis metadata.
  pairs <- unique(map[c("node", "quantity")])
  if( anyDuplicated(pairs$node) ){
    stop("A displayed node maps to multiple canonical quantity IDs.", call. = FALSE)
  }
  attr(result, "input_labels") <- stats::setNames(pairs$quantity, .strip_slice(pairs$node))
  attr(result, "cf_node_map") <- map
  return(result)
}


#' Resolve retention or latent requests against saved identities
#'
#' @param graph Graph containing the requested quantities or occasions.
#' @param requested Character node, quantity, member or occurrence names.
#'   `NULL` means no requests. Saved aliases are also accepted.
#' @param record Saved query record supplying maps, aliases and disambiguation.
#' @return The ordered union of matching exact node names, or `character()` for
#'   an empty request. Absent or ambiguous names raise an error, as does an exact
#'   occasion already merged with others.
#' @noRd
.cd_match_nodes <- function(graph, requested, record = .cd_record(graph)){
  requested <- .gd_names(requested, "node names")
  map <- .cd_node_map(graph, record)
  resolved_nodes <- character()
  for( requested_name in requested ){
    exact_time <- requested_name %in% map$occurrence && !requested_name %in% c(map$quantity, map$member_id)
    matched_nodes <- unique(map$node[map$quantity == requested_name | map$member_id == requested_name |
      map$node == requested_name | map$occurrence == requested_name])
    if( !length(matched_nodes) && requested_name %in% names(graph) ){
      matched_nodes <- requested_name
    }
    # Recorded input aliases are safe to replay without consulting changed tables.
    if( !length(matched_nodes) ){
      requests <- record$requests
      resolved <- unlist(lapply(requests, function(request){
        return(request$resolved)
      }), use.names = FALSE)
      aliases <- unlist(lapply(requests, function(request){
        return(names(request$resolved))
      }), use.names = FALSE)
      canonical <- unique(resolved[aliases == requested_name])
      matched_nodes <- unique(map$node[map$member_id %in% canonical | map$quantity %in% canonical])
    }
    if( !length(matched_nodes) && length(record$aliases) ){
      resolved <- .resolve_columns(requested_name, .build_alias_index(record$aliases),
        record$arguments$params$disambiguate, character())
      if( length(resolved$ambiguous) ){
        stop("Ambiguous node name: ", requested_name, call. = FALSE)
      }
      canonical <- unname(resolved$resolved)
      matched_nodes <- unique(map$node[map$member_id %in% canonical | map$quantity %in% canonical])
    }
    if( exact_time && length(matched_nodes) && any(vapply(matched_nodes, function(node){
      return(length(unique(map$occurrence[map$node == node])) > 1L)
    }, logical(1))) ){
      stop("The supplied graph merges multiple occasions for ", requested_name, ". A request for an exact occasion cannot distinguish them.", call. = FALSE)
    }
    if( !length(matched_nodes) ){
      stop("Quantity absent from the supplied graph: ", requested_name, call. = FALSE)
    }
    resolved_nodes <- union(resolved_nodes, matched_nodes)
  }
  return(resolved_nodes)
}


#' Prepare construction inputs with exact cache reuse
#'
#' @inheritParams .cd_prepare_uncached
#' @return The preparation record from `.cd_prepare_uncached()`. Only standard
#'   graphs use the cache, keyed by the complete graph and supplied requests.
#' @noRd
.cd_prepare <- function(graph, latents, keep_nodes){
  if( !.cf_acceleration_enabled() || !.cf_standard_graph(graph) ){
    return(.cd_prepare_uncached(graph, latents, keep_nodes))
  }
  prepared <- .cf_cache$call("preparation", list(graph, latents, keep_nodes),
    function(){
      return(.cd_prepare_uncached(graph, latents, keep_nodes))
    })

  return(prepared)
}


#' Prepare the supplied graph for adjustment search
#'
#' @param graph An acyclic dagitty graph with disjoint treatment and outcome tags.
#' @param latents Variables or exact occasions to declare unmeasured.
#' @param keep_nodes Additional variables or occasions to retain.
#' @return A list containing the source graph and record, identity map, reference,
#'   analysis and selection graphs, policy ledger, required and latent nodes,
#'   limitations and query controls. Policy exclusions are folded through directed
#'   paths. Requests conflicting with that policy raise an error.
#' @noRd
.cd_prepare_uncached <- function(graph, latents, keep_nodes){
  if( !inherits(graph, "dagitty") ){
    stop("graph must be a dagitty object.", call. = FALSE)
  }
  if( !length(dagitty::exposures(graph)) || !length(dagitty::outcomes(graph)) ){
    stop("The supplied graph must tag treatments and outcomes.", call. = FALSE)
  }
  record <- .cd_record(graph)
  map <- .cd_node_map(graph, record)
  reference <- .cd_working_graph(graph, map)
  if( !dag_is_acyclic(reference) ){
    stop("The supplied graph has directed cycles. Provide an acyclic query configuration.", call. = FALSE)
  }
  if( length(intersect(dagitty::exposures(reference), dagitty::outcomes(reference))) ){
    stop("Treatment and outcome declarations overlap. Construction requires disjoint endpoint sets.", call. = FALSE)
  }
  latent_nodes <- .cd_match_nodes(graph, latents, record)
  if( length(intersect(latent_nodes, c(dagitty::exposures(graph), dagitty::outcomes(graph)))) ){
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
  if( is.null(crop) ){
    crop <- attr(graph, "crop", exact = TRUE)
  }
  applied <- (isTRUE(record$policy_applied) || identical(crop, "barley")) && is.data.frame(policy) &&
  all(c("node_id", "recovered_node") %in% names(policy))
  member_rule <- rep("yes", nrow(map))
  if( applied ){
    policy_rows <- match(map$member_id, policy$node_id)
    member_rule[!is.na(policy_rows)] <- as.character(policy$recovered_node[policy_rows[!is.na(policy_rows)]])
  }
  if( anyNA(member_rule) || any(!member_rule %in% c("yes", "yes:user") & !startsWith(member_rule, "no:")) ){
    stop("The query receipt contains invalid recovery policy values.", call. = FALSE)
  }
  explicit_members <- unique(map$member_id[map$node %in% required])
  allowed_member <- member_rule == "yes" | (member_rule == "yes:user" & map$member_id %in% explicit_members)
  eligible <- vapply(names(reference), function(node){
    return(all(allowed_member[map$node == node]))
  }, logical(1))
  prohibited <- intersect(required, names(reference)[!eligible])
  if( length(prohibited) ){
    stop("Requested nodes conflict with the saved recovery policy: ", paste(prohibited, collapse = ", "), call. = FALSE)
  }
  rules <- vapply(names(reference), function(node){
    return(paste(unique(member_rule[map$node == node]), collapse = "; "))
  }, character(1))
  ledger <- data.frame(node = names(reference), node_id = unname(.gd_node_ids(reference)),
    metadata_row = names(reference) %in% map$node, recovered_node = rules,
    explicitly_requested = names(reference) %in% required, policy_applied = applied,
    eligible = eligible, action = ifelse(eligible, "retain", "fold"), stringsAsFactors = FALSE)
  excluded <- ledger$node[!eligible]
  analysis <- .cd_fold_policy(reference, excluded)
  analysis <- .cd_working_graph(analysis, map[map$node %in% names(analysis), , drop = FALSE])
  # Recovery policy folds excluded nodes through directed paths instead of
  # marking them latent. Search and certificates use the folded graph. The
  # supplied graph and the reference before folding are retained for
  # provenance.
  selection <- .gd_canonical(analysis)
  depth <- record$arguments$params$query_depth
  if( is.null(depth) ){
    depth <- attr(graph, "query_depth", exact = TRUE)
  }
  limitations <- character()
  if( length(excluded) ){
    limitations <- c(limitations,
      "Recovery policy exclusions were folded through directed paths, not treated as latent variables. Adjustment checks use the folded graph. The supplied query is preserved separately.")
  }
  if( is.null(record) ){
    limitations <- c(limitations, "No saved query record. The original requests and policy information are unavailable. Size uses recorded variable identities or exact node names.")
  }
  if( length(depth) && is.finite(depth) ){
    limitations <- c(limitations,
      paste0("Checks use the supplied query_depth=", depth, " graph only. Excluded context was not substituted from a full snapshot."))
  }
  if( length(missing_columns) ){
    limitations <- c(limitations,
      paste0("Requested columns absent from the supplied graph: ", paste(missing_columns, collapse = ", "), "."))
  }
  unresolved <- record$requests$columns$unresolved
  if( is.list(unresolved) ){
    unresolved_names <- unique(c(unresolved$unmatched, names(unresolved$ambiguous)))
    if( length(unresolved_names) ){
      limitations <- c(limitations,
        paste0("Column requests unresolved by dagby_query: ", paste(unresolved_names, collapse = ", "), "."))
    }
  }
  prepared <- list(graph = graph, record = record, map = map, reference = reference,
    selection = selection, analysis = analysis, ledger = ledger,
    required = required, latent_nodes = latent_nodes,
    limitations = limitations, policy_applied = applied,
    controls = list(crop = crop, water_model = if( !is.null(record$arguments$params$water_model) ){
      record$arguments$params$water_model
    }else{
      attr(graph, "water_model")
    },
    disambiguate = record$arguments$params$disambiguate, simplify = TRUE, keep_common_causes = FALSE))

  return(prepared)
}


#' Fold recovery-policy exclusions through directed paths
#'
#' @param graph Reference graph before policy exclusions.
#' @param excluded Exact node names prohibited by the saved policy.
#' @return A graph without those nodes, preserving directed paths through them
#'   and recorded source order. Other incident edge types raise an error.
#' @noRd
.cd_fold_policy <- function(graph, excluded){
  if( !length(excluded) ){
    return(graph)
  }
  edges <- dag_edges(graph)
  incident <- edges$v %in% excluded | edges$w %in% excluded
  if( any(incident & edges$e != "->") ){
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
  if( is.null(result) ){
    stop("Recovery policy directed folding could not preserve the remaining graph.", call. = FALSE)
  }
  source_order <- attr(graph, "cf_source_order", exact = TRUE)
  if( is.null(source_order) ){
    source_order <- .cd_record(graph)$source_order
  }
  attr(result, "cf_source_order") <- source_order
  return(result)
}


#' Catalogue adjustment sets using the prepared source order
#'
#' @param raw List of adjustment sets expressed as exact node names.
#' @param prepared Preparation record from `.cd_prepare()`.
#' @return The `.gd_candidates()` table. Recorded order is used only when it
#'   covers the entire graph, identity map and set family.
#' @noRd
.cd_candidates <- function(raw, prepared){
  source_order <- attr(prepared$selection, "cf_source_order", exact = TRUE)
  if( !is.null(source_order) ){
    family_ids <- unname(.gd_node_ids(prepared$selection)[unlist(raw, use.names = FALSE)])
    if( !.query_order_applies(c(unname(.gd_node_ids(prepared$selection)), family_ids,
      prepared$map$quantity, prepared$map$member_id), source_order) ){
      source_order <- NULL
    }
  }
  table <- .gd_candidates(raw, prepared$selection, source_order = source_order)
  return(table)
}


#' Obtain an exact search identity without affecting construction
#'
#' @param prepared Preparation record from `.cd_prepare()`.
#' @return The identity from `.cd_search_identity_value()`, or `NULL` on any
#'   error, warning or message so that uncertain identity disables reuse.
#' @noRd
.cd_search_identity <- function(prepared){
  return(tryCatch(.cd_search_identity_value(prepared), error = function(e){
    return(NULL)
  },
    warning = function(w){
      return(NULL)
    }, message = function(m){
      return(NULL)
    }))
}


#' Bind search reuse to source state and implementation
#'
#' The identity includes graph text and order, occurrence identities, declarations,
#' policy, query scope, collation and separator implementation. It omits retention
#' requests themselves. Their effects on the prepared graph and recovery policy
#' remain included.
#' @param prepared Preparation record from `.cd_prepare()`.
#' @return An identity list including engine code and source records, or `NULL`
#'   when the installed dagitty script cannot be read.
#' @noRd
.cd_search_identity_value <- function(prepared){
  record <- prepared$record
  arguments <- record$arguments
  arguments[c("columns", "quiet")] <- NULL
  policy <- prepared$ledger[setdiff(names(prepared$ledger), "explicitly_requested")]
  script <- system.file("js", "dagitty-alg.js", package = "dagitty")
  if( !nzchar(script) || !file.exists(script) ){
    return(NULL)
  }
  identity <- list(version = 2L, method = "complete-source-sets-v2",
    construction = .cd_construction_version,
    engine = list(dagitty = as.character(utils::packageVersion("dagitty")),
    native = readBin(script, "raw", n = file.info(script)$size),
    adapter = deparse(body(.gd_adjustment_sets)),
    implementation = deparse(body(.gd_adjustment_engine)),
    separator = "Installed dagitty GraphAnalyzer.listMinimalSeparators", indexed = FALSE),
    # Aliases and case were resolved during preparation. The separator
    # sees UTF-8 graph text and uses JavaScript ordering, independent of CTYPE.
    # Saved signatures and candidate records still use R collation.
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
  if( !is.null(source_order) ){
    identity$source_order <- source_order
    identity$ordering_adapter <- list(
      order_key = deparse(body(.query_order_key)), sort = deparse(body(.query_sort)),
      applies = deparse(body(.query_order_applies)), set_key = deparse(body(.query_set_key)))
  }
  return(identity)
}


#' Retain exact treatment parents when no total adjustment option is found
#'
#' @param graph Working reference graph before construction.
#' @param total_choice Total-effect choice with search status and reason.
#' @return A safeguard record with `applied`, `status` and exact treatment parents.
#'   The safeguard applies only to fully observed directed DAGs. The caller must
#'   verify that these parents were retained before certifying joint total
#'   identification.
#' @noRd
.cd_identification_retention <- function(graph, total_choice){
  result <- list(applied = FALSE, parents = character(), added_nodes = character(),
    status = "not required: a total adjustment option was selected",
    scope = "Joint total effect on the supplied working reference. This does not cover all adjustment alternatives or other effects.")
  if( !identical(total_choice$status, "no_sets") ){
    return(result)
  }
  result$status <- "not certified: the reference is not a fully observed directed DAG"
  if( dag_graph_type(graph) != "dag" ||
    length(dagitty::latents(graph)) || length(dagitty::adjustedNodes(graph)) ||
    length(selected_nodes(graph)) || any(dag_edges(graph)$e != "->") ){
    return(result)
  }
  # From the observed DAG g-formula, every intervention denominator depends
  # only on X and Pa(X). Keeping these exact vertices allows the other observed
  # variables to be marginalised without losing the joint total functional.
  # Pa(X) need not be an ordinary adjustment set (it may follow another X).
  treatment <- dagitty::exposures(graph)
  result$parents <- sort(unique(unlist(lapply(treatment, function(treatment_node){
    return(dagitty::parents(graph, treatment_node))
  }), use.names = FALSE)))
  result$applied <- TRUE
  result$status <- "pending structural verification"
  result$reason <- if( isTRUE(grepl("step budget", total_choice$reason, fixed = TRUE)) ){
    "No total adjustment option was found within the ranked search's step budget (absence is not certified). Exact treatment parents and temporal occurrences are retained for joint total identification."
  }else{
    "No total adjustment option. Retain exact treatment parents and temporal occurrences for joint total identification."
  }
  return(result)
}


## Graph transformations and time summaries ##
# Internal graph transformations used by dagby_construct(). Analysis remains temporal.
# Summaries of repeated processes are display abstractions with separate ledgers.

#' Validate and deduplicate node requests
#'
#' @param x Character names, or `NULL` for an empty request.
#' @param argument Input name used in errors.
#' @return Unique unnamed character values in supplied order, or `character()`
#'   for `NULL`. Missing or empty names raise an error.
#' @noRd
.gd_names <- function(x, argument){
  if( is.null(x) ){
    return(character())
  }
  if( !is.character(x) || anyNA(x) || any(!nzchar(x)) ){
    stop(argument, " must contain non-missing, non-empty names.", call. = FALSE)
  }
  node_names <- unique(unname(x))

  return(node_names)
}


#' Reparse graph text while preserving analysis metadata
#'
#' @param dag A graph accepted by `dag_graph_type()` and dagitty.
#' @return A dagitty graph with labels and map rows restricted to its nodes.
#'   Source ordering is retained. Conflicting graph and record orderings raise
#'   an error.
#' @noRd
.gd_canonical <- function(dag){
  dag_graph_type(dag)
  result <- dagitty::dagitty(as.character(dag))
  labels <- attr(dag, "input_labels")
  if( length(labels) ){
    attr(result, "input_labels") <- labels[names(labels) %in% .strip_slice(names(result))]
  }
  map <- attr(dag, "cf_node_map", exact = TRUE)
  if( is.data.frame(map) ){
    attr(result, "cf_node_map") <- map[map$node %in% names(result), , drop = FALSE]
  }
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  recorded_order <- .cd_record(dag)$source_order
  if( !is.null(source_order) && !is.null(recorded_order) && !identical(source_order, recorded_order) ){
    stop("The graph and its query record have different source ordering metadata.", call. = FALSE)
  }
  if( is.null(source_order) ){
    source_order <- recorded_order
  }
  if( !is.null(source_order) ){
    attr(result, "cf_source_order") <- source_order
  }
  return(result)
}


#' Read canonical quantity identities with optional cache reuse
#'
#' @inheritParams .gd_node_ids_uncached
#' @return The named quantity vector from `.gd_node_ids_uncached()`.
#' @noRd
.gd_node_ids <- function(dag){
  if( !.cf_acceleration_enabled() || !.cf_standard_graph(dag) ){
    return(.gd_node_ids_uncached(dag))
  }
  ids <- .gd_metadata(dag)$ids

  return(ids)
}


#' Cache graph identities and source ordering
#'
#' @param dag A standard graph whose complete value forms the cache key.
#' @return A list with named quantity `ids`, exact `nodes` and recorded
#'   `source_order` when available.
#' @noRd
.gd_metadata <- function(dag){
  return(.cf_cache$call("metadata", list(dag), function(){
    metadata <- list(ids = .gd_node_ids_uncached(dag), nodes = names(dag))
    source_order <- attr(dag, "cf_source_order", exact = TRUE)
    if( is.null(source_order) ){
      source_order <- .cd_record(dag)$source_order
    }
    if( !is.null(source_order) ){
      metadata$source_order <- source_order
    }
    return(metadata)
  }))
}


#' Resolve exact nodes to canonical quantities
#'
#' @param dag Graph with optional reverse input labels and occurrence map.
#' @return A character vector named by exact graph nodes. Map quantities take
#'   precedence over reverse labels, then names with temporal suffixes removed.
#' @noRd
.gd_node_ids_uncached <- function(dag){
  ids <- .strip_slice(names(dag))
  labels <- attr(dag, "input_labels")
  matched <- ids %in% names(labels)
  ids[matched] <- unname(labels[ids[matched]])
  map <- attr(dag, "cf_node_map", exact = TRUE)
  if( is.data.frame(map) ){
    rows <- match(names(dag), map$node)
    ids[!is.na(rows)] <- map$quantity[rows[!is.na(rows)]]
  }
  node_ids <- stats::setNames(ids, names(dag))

  return(node_ids)
}


#' Find typed edges introduced by a graph change
#'
#' @param before,after Graphs before and after the change.
#' @return Rows of `dag_edges(after)` whose typed edge keys are absent before.
#' @noRd
.gd_new_edges <- function(before, after){
  previous_edges <- dag_edges(before)
  current_edges <- dag_edges(after)
  new_edges <- current_edges[!edge_key(current_edges) %in% edge_key(previous_edges), , drop = FALSE]

  return(new_edges)
}


#' Read numeric occasions from exact node suffixes
#'
#' @param nodes Character node names ending in `_tN` or `@tN` when timed.
#' @return A numeric vector in input order, with `NA` for unrecognised suffixes.
#' @noRd
.gd_time <- function(nodes){
  occasions <- rep(NA_real_, length(nodes))
  timed <- grepl("(?:_|@)t[0-9]+$", nodes, perl = TRUE)
  occasions[timed] <- as.numeric(sub("^.*(?:_|@)t([0-9]+)$", "\\1", nodes[timed],
    perl = TRUE))
  return(occasions)
}


#' Group temporal nodes by represented source variables
#'
#' @param dag Graph with optional source variable identities.
#' @return A data frame with `node`, `quantity` and `identity` columns. Identity
#'   encodes complete variable membership. Quantity supplies a display label that
#'   distinguishes groups. Without a map, canonical quantity identities and
#'   node labels without temporal suffixes are used.
#' @noRd
.gd_summary_groups <- function(dag){
  nodes <- names(dag)
  map <- attr(dag, "cf_node_map", exact = TRUE)
  if( is.null(map) ){
    map <- .cd_record(dag)$node_map
  }
  encode <- function(values){
    return(paste0(nchar(values, type = "bytes"), ":", values, collapse = ""))
  }
  if( is.null(map) ){
    return(data.frame(node = nodes, quantity = .strip_slice(nodes),
      identity = vapply(unname(.gd_node_ids(dag)), encode, character(1)),
      stringsAsFactors = FALSE))
  }
  variables <- .cd_variable_map(dag)
  identity <- stats::setNames(vapply(nodes, function(node){
    return(encode(sort(variables$variable[variables$node == node], method = "radix")))
  }, character(1)), nodes)
  groups <- split(nodes, identity)
  representative <- vapply(groups, function(group){
    return(sort(group, method = "radix")[[1L]])
  }, character(1))
  label <- vapply(names(groups), function(key){
    group <- groups[[key]]
    members <- unique(variables$variable[variables$node %in% group])
    if( length(group) == 1L && identical(members, group) ){
      return(group)
    }
    preferred <- sort(unique(.strip_slice(group)), method = "radix")[[1L]]
    if( preferred %in% nodes && identity[[preferred]] != key ){
      return(representative[[key]])
    }
    return(preferred)
  }, character(1))
  names(label) <- names(groups)
  collision <- duplicated(label) | duplicated(label, fromLast = TRUE)
  label[collision] <- representative[collision]
  return(data.frame(node = nodes, quantity = unname(label[identity]),
    identity = unname(identity), stringsAsFactors = FALSE))
}


#' Record temporal edges and their proposed display groups
#'
#' @param dag Temporal graph before summarisation.
#' @param groups Node grouping from `.gd_summary_groups()`.
#' @return Typed edge rows with the quantity and numeric occasion at each end,
#'   and an initial `display_role` of `"shown"`.
#' @noRd
.gd_edge_ledger <- function(dag, groups = .gd_summary_groups(dag)){
  edges <- dag_edges(dag)
  edges$from <- groups$quantity[match(edges$v, groups$node)]
  edges$to <- groups$quantity[match(edges$w, groups$node)]
  edges$time_from <- .gd_time(edges$v)
  edges$time_to <- .gd_time(edges$w)
  edges$display_role <- rep("shown", nrow(edges))
  return(edges)
}


#' Record proper directed paths between endpoint identities
#'
#' @param dag Graph with treatment and outcome tags and source variable identities.
#'   Endpoint identities must be unique, as checked by `.gd_summary()`.
#' @return A logical vector named by treatment and outcome identity pairs in radix
#'   order. Each value is `TRUE` if a directed path connects that pair without
#'   passing through another treatment, and `FALSE` otherwise.
#' @noRd
.gd_endpoint_paths <- function(dag){
  treatments <- as.character(dagitty::exposures(dag))
  outcomes <- as.character(dagitty::outcomes(dag))
  edges <- dag_edges(dag)
  children <- split(edges$w[edges$e == "->"], edges$v[edges$e == "->"])
  groups <- .gd_summary_groups(dag)
  identity <- stats::setNames(groups$identity, groups$node)
  paths <- logical()
  for( treatment in treatments ){
    blocked <- setdiff(treatments, treatment)
    reached <- character()
    pending <- treatment
    while( length(pending) ){
      node <- pending[[1L]]
      pending <- pending[-1L]
      if( node %in% reached || node %in% blocked ){
        next
      }
      reached <- c(reached, node)
      pending <- union(pending, setdiff(children[[node]], c(reached, blocked)))
    }
    keys <- paste0(nchar(identity[[treatment]], type = "bytes"), ":",
      identity[[treatment]], identity[outcomes])
    paths[keys] <- outcomes %in% setdiff(reached, treatment)
  }
  paths <- paths[order(names(paths), method = "radix")]

  return(paths)
}


#' Find feedback groups in a folded quantity graph
#'
#' @param from,to Quantity names at the ends of directed folded edges.
#' @param source_order Optional recorded source order covering all quantities.
#' @param node_ids Optional canonical quantity IDs named by folded quantity.
#' @return Strongly connected groups with at least two quantities. Members and
#'   groups follow recorded order when applicable, otherwise radix name order.
#' @noRd
.gd_feedback_groups <- function(from, to, source_order = NULL, node_ids = NULL){
  quantities <- unique(c(from, to))
  if( !.query_order_applies(quantities, source_order, unname(node_ids[quantities])) ){
    source_order <- NULL
  }
  quantities <- .query_sort(quantities, source_order, unname(node_ids[quantities]))
  children <- split(to, factor(from, levels = quantities))
  reach <- lapply(stats::setNames(nm = quantities), function(start){
    seen <- character()
    frontier <- start
    while( length(frontier) ){
      frontier <- setdiff(unique(unlist(children[frontier], use.names = FALSE)), seen)
      seen <- c(seen, frontier)
    }
    return(seen)
  })
  groups <- lapply(quantities, function(quantity){
    returning <- Filter(function(other){
      return(quantity %in% reach[[other]])
    }, reach[[quantity]])
    members <- unique(c(quantity, returning))
    group <- .query_sort(members, source_order, unname(node_ids[members]))

    return(group)
  })
  groups <- unique(groups[lengths(groups) > 1L])
  first_members <- vapply(groups, `[[`, character(1), 1L)
  groups <- groups[order(.query_order_key(first_members, source_order, unname(node_ids[first_members])), method = "radix")]

  return(groups)
}


#' Enumerate contiguous time cuts
#'
#' @param occasions Sorted numeric occasions.
#' @param k Number of desired segments.
#' @return A list of cut vectors in lexicographic order. Cuts end every segment
#'   except the last. One segment gives `list(numeric())`. For more than one
#'   segment, requesting more segments than occasions gives `list()`.
#' @noRd
.gd_cut_candidates <- function(occasions, k){
  if( k == 1L ){
    return(list(numeric()))
  }
  closers <- occasions[-length(occasions)]
  if( length(closers) < k - 1L ){
    return(list())
  }
  index <- utils::combn(length(closers), k - 1L)
  candidates <- lapply(seq_len(ncol(index)), function(combination){
    return(closers[index[, combination]])
  })

  return(candidates)
}


#' Name display nodes for contiguous time segments
#'
#' @param mapping Data frame relating exact temporal `node` to display `quantity`.
#' @param cuts Named list of cut occasions for each segmented quantity.
#' @return Display names in mapping row order. An uncut quantity keeps its name.
#'   After cutting, a segment with one temporal node keeps its exact name and a
#'   segment with several temporal nodes uses a time range.
#'   Untimed occurrences of a cut quantity remain in a separate segment.
#' @noRd
.gd_segment_display <- function(mapping, cuts){
  display <- mapping$quantity
  time <- .gd_time(mapping$node)
  occasion_text <- function(nodes){
    return(sub("^.*(?:_|@)t([0-9]+)$", "\\1", nodes, perl = TRUE))
  }
  for( quantity in names(cuts) ){
    rows <- which(mapping$quantity == quantity)
    segment <- vapply(time[rows], function(occasion){
      return(if( is.na(occasion) ){
        0L
      }else{
        1L + sum(occasion > cuts[[quantity]])
      })
    }, integer(1))
    if( length(unique(segment)) < 2L ){
      next
    }
    for( segment_id in unique(segment) ){
      members <- rows[segment == segment_id]
      display[members] <- if( length(members) == 1L ){
        mapping$node[members]
      }else{
        paste0(quantity, "_t", occasion_text(mapping$node[members][which.min(time[members])]),
          "to", occasion_text(mapping$node[members][which.max(time[members])]))
      }
    }
  }
  return(display)
}


#' Map retained temporal edges onto display nodes
#'
#' @param links Edge ledger with `v`, `e`, `w` and `display_role` columns.
#' @param display Character display names keyed by temporal node.
#' @return Unique typed edges selected from rows whose `display_role` is `"shown"`
#'   or `"within_quantity"`, excluding self edges. Bidirected endpoints are
#'   canonicalised. Other roles stay hidden.
#' @noRd
.gd_display_edges <- function(links, display){
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


#' Check directed edges for cycles
#'
#' @param edges Typed edge data frame with `v`, `e` and `w` columns.
#' @return `TRUE` if the directed edges are acyclic, and `FALSE` otherwise.
#'   Uses Kahn's algorithm and ignores other edge types.
#' @noRd
.gd_edges_acyclic <- function(edges){
  directed <- edges[edges$e == "->", , drop = FALSE]
  nodes <- unique(c(directed$v, directed$w))
  indegree <- tabulate(match(directed$w, nodes), length(nodes))
  children <- split(match(directed$w, nodes), factor(directed$v, levels = nodes))
  frontier <- which(indegree == 0L)
  removed <- 0L
  while( length(frontier) ){
    node <- frontier[[1L]]
    frontier <- frontier[-1L]
    removed <- removed + 1L
    for( child in children[[node]] ){
      indegree[[child]] <- indegree[[child]] - 1L
      if( indegree[[child]] == 0L ){
        frontier <- c(frontier, child)
      }
    }
  }
  acyclic <- removed == length(nodes)

  return(acyclic)
}


#' Choose contiguous cuts for all feedback groups
#'
#' Each group first finds the fewest segments that make its own edges acyclic.
#' Its quantities share the cuts. Combinations across groups are then compared
#' by display edge count, followed by sorted cut occasions. Remaining ties use
#' group order, taking the earliest cuts for the first group, then the next.
#' `alternatives` records these remaining ties so the caller can report them.
#' @param links Temporal edge ledger with display roles.
#' @param mapping Temporal node-to-quantity mapping.
#' @param groups Ordered feedback groups from `.gd_feedback_groups()`.
#' @return A list of per-quantity `cuts` and the number of `alternatives` tied
#'   before group order decides. Fewest segments, fewest edges and earliest
#'   cut occasions take precedence over group order.
#' @noRd
.gd_segment_choice <- function(links, mapping, groups){
  time <- .gd_time(mapping$node)
  options <- lapply(groups, function(group){
    occasions <- sort(unique(time[mapping$quantity %in% group & !is.na(time)]))
    for( k in seq_len(max(1L, length(occasions))) ){
      acyclic <- Filter(function(candidate){
        display <- stats::setNames(.gd_segment_display(mapping,
          stats::setNames(rep(list(candidate), length(group)), group)), mapping$node)
        edges <- .gd_display_edges(links, display)
        members <- unique(display[mapping$quantity %in% group])
        group_acyclic <- .gd_edges_acyclic(edges[edges$v %in% members & edges$w %in% members, , drop = FALSE])

        return(group_acyclic)
      }, .gd_cut_candidates(occasions, k))
      if( length(acyclic) ){
        return(acyclic)
      }
    }
    # Only untimed occurrences or two occurrences at one occasion (A_t0 and
    # A@t0) defeat every cut. Keep each occasion apart and let the display's
    # acyclicity check refuse the summary.
    candidates <- list(occasions[-length(occasions)])

    return(candidates)
  })
  if( !length(options) ){
    return(list(cuts = list(), alternatives = 1L))
  }
  grid <- as.matrix(expand.grid(lapply(options, seq_along), KEEP.OUT.ATTRS = FALSE))
  cuts_of <- function(row){
    cuts <- list()
    for( i in seq_along(groups) ){
      cuts[groups[[i]]] <- rep(list(options[[i]][[row[[i]]]]), length(groups[[i]]))
    }
    return(cuts)
  }
  edges <- integer(nrow(grid))
  positions <- NULL
  for( combination in seq_len(nrow(grid)) ){
    row <- grid[combination, ]
    display <- stats::setNames(.gd_segment_display(mapping, cuts_of(row)), mapping$node)
    edges[[combination]] <- nrow(.gd_display_edges(links, display))
    positions <- rbind(positions, sort(unlist(lapply(seq_along(groups), function(i){
      return(options[[i]][[row[[i]]]])
    }))))
  }
  ranking <- do.call(order, c(list(edges), lapply(seq_len(ncol(positions)), function(column){
    return(positions[, column])
  }),
    lapply(seq_len(ncol(grid)), function(column){
      return(grid[, column])
    })))
  best <- ranking[[1L]]
  same_positions <- if( ncol(positions) ){
    rowSums(positions != matrix(positions[best, ], nrow(positions),
      ncol(positions), byrow = TRUE)) == 0
  }else{
    rep(TRUE, nrow(grid))
  }
  choice <- list(cuts = cuts_of(grid[best, ]), alternatives = sum(edges == edges[[best]] & same_positions))

  return(choice)
}


#' Report ties only for the selected segmented display
#'
#' @param summary Applied summary carrying segment and repeated-edge tie counts.
#' @return `summary`, invisibly, after warning when quantity order broke a tie.
#'   Callers discarding a trial summary should not call this reporter.
#' @noRd
.gd_segment_tie_notice <- function(summary){
  if( !isTRUE(summary$applied) ){
    return(invisible(summary))
  }
  parts <- c(
    if( isTRUE(summary$segment_alternatives > 1L) ){
      paste0(summary$segment_alternatives, " equally ranked cut choices")
    },
    if( isTRUE(summary$segment_edge_ties > 0L) ){
      paste0(summary$segment_edge_ties, " edge(s) from repeated stages left omitted only because an edge of the same lag and occasion came first in quantity order")
    })
  if( length(parts) ){
    warning("Time segments: ", paste(parts, collapse = ", "), ". Quantity order decided ",
      "(summary segment_alternatives / segment_edge_ties). A tie between feedback quantities is structural information.",
      call. = FALSE)
  }
  return(invisible(summary))
}


#' Map display nodes to exact source occurrences and members
#'
#' @param trimmed Temporal graph retained by construction.
#' @param display Display names keyed by exact temporal nodes.
#' @return A data frame with `node`, `quantity`, `member_id` and `occurrence`
#'   columns, with one row per represented occurrence and member. `node` names
#'   the display node. Unmapped temporal nodes receive their canonical quantity
#'   identity and exact occurrence name.
#' @noRd
.gd_display_map <- function(trimmed, display){
  columns <- c("node", "quantity", "member_id", "occurrence")
  map <- attr(trimmed, "cf_node_map", exact = TRUE)
  rows <- if( is.data.frame(map) && all(columns %in% names(map)) ){
    map[map$node %in% names(trimmed), columns, drop = FALSE]
  }else{
    NULL
  }
  missing <- setdiff(names(trimmed), rows$node)
  if( length(missing) ){
    ids <- unname(.gd_node_ids(trimmed)[missing])
    rows <- rbind(rows, data.frame(node = missing, quantity = ids, member_id = ids,
      occurrence = missing, stringsAsFactors = FALSE))
  }
  rows$node <- unname(display[rows$node])
  rownames(rows) <- NULL
  return(rows)
}


#' Describe quantities displayed in several time segments
#'
#' @param mapping Data frame with temporal node, quantity, node_id and display.
#' @param source_order Optional recorded quantity order.
#' @return A data frame with quantity, segment, first and last occasions and exact
#'   members, or `NULL` when no quantity has several display nodes.
#' @noRd
.gd_segment_table <- function(mapping, source_order = NULL){
  if( !.query_order_applies(mapping$node, source_order, mapping$node_id) ){
    source_order <- NULL
  }
  pieces <- tapply(mapping$display, mapping$quantity, function(display_nodes){
    return(length(unique(display_nodes)))
  })
  cut <- names(pieces)[pieces > 1L]
  if( !length(cut) ){
    return(NULL)
  }
  rows <- mapping[mapping$quantity %in% cut, , drop = FALSE]
  time <- .gd_time(rows$node)
  keys <- unique(rows[c("quantity", "display")])
  table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i){
    members <- which(rows$quantity == keys$quantity[i] & rows$display == keys$display[i])
    member_keys <- .query_order_key(rows$node[members], source_order, rows$node_id[members])
    members <- members[order(time[members], member_keys, method = "radix")]
    occasions <- time[members][!is.na(time[members])]
    segment_row <- data.frame(quantity = keys$quantity[i], segment = keys$display[i],
      first = if( length(occasions) ){
        min(occasions)
      }else{
        NA_real_
      },
      last = if( length(occasions) ){
        max(occasions)
      }else{
        NA_real_
      },
      members = paste(rows$node[members], collapse = ", "), stringsAsFactors = FALSE)

    return(segment_row)
  }))
  quantity_ids <- rows$node_id[match(table$quantity, rows$quantity)]
  quantity_keys <- .query_order_key(table$quantity, source_order, quantity_ids)
  segment_keys <- .query_order_key(table$segment, source_order, quantity_ids)
  table <- table[order(quantity_keys, table$first, segment_keys, method = "radix"), , drop = FALSE]
  rownames(table) <- NULL
  return(table)
}


#' Build a checked static display of a temporal graph
#'
#' @param dag Constructed temporal graph to summarise.
#' @param keep Exact nodes or canonical quantities whose routes must remain.
#' @param preserve_occasions Whether identification requires the temporal view.
#' @param source_check Optional callback accepting a proposed display and its
#'   temporal-to-display mapping. It returns `NULL` if source checks pass, or
#'   text explaining why the proposed display failed them.
#' @return A list with display `graph`, retained `temporal` graph, node map,
#'   trimming and edge ledgers and an `applied` flag. If a summary is refused,
#'   `graph` keeps the temporal graph and `reason` explains the refusal.
#'   Successful summaries also record segment choices.
#' @noRd
.gd_summary <- function(dag, keep = character(), preserve_occasions = FALSE, source_check = NULL){
  keep <- .gd_names(keep, "keep")
  all_nodes <- names(dag)
  all_ids <- unname(.gd_node_ids(dag))
  variable_groups <- .gd_summary_groups(dag)
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  if( is.null(source_order) ){
    source_order <- .cd_record(dag)$source_order
  }
  recorded_order <- source_order
  if( !.query_order_applies(all_nodes, source_order, all_ids) ){
    source_order <- NULL
  }
  missing <- setdiff(keep, c(all_nodes, all_ids))
  if( length(missing) ){
    stop("Requested summary quantities absent from this graph: ",
      paste(missing, collapse = ", "), call. = FALSE)
  }
  fallback <- function(reason){
    ledger <- .gd_edge_ledger(dag, variable_groups)
    ledger$display_role <- rep("temporal_retained", nrow(ledger))
    temporal_summary <- list(graph = dag, temporal = dag,
      trim_ledger = data.frame(node = all_nodes, node_id = all_ids,
      action = rep("retain", length(all_nodes)), stringsAsFactors = FALSE),
      node_map = data.frame(node = all_nodes, quantity = variable_groups$quantity,
      node_id = all_ids,
      stringsAsFactors = FALSE),
      edge_ledger = ledger, applied = FALSE, reason = reason)

    return(temporal_summary)
  }
  if( isTRUE(preserve_occasions) ){
    return(fallback("Exact temporal occurrences retained for the joint total identification safeguard. Time merging is not certified."))
  }
  if( dag_graph_type(dag) != "dag" ){
    return(fallback("Time summarisation requires a DAG with directed and bidirected edges."))
  }
  if( !isTRUE(dag_is_acyclic(dag)) ){
    return(fallback("The input has a directed cycle. Its temporal display was retained."))
  }
  if( any(!dag_edges(dag)$e %in% c("->", "<->")) ){
    return(fallback("Unsupported edge types prevent time summarisation."))
  }
  if( length(dagitty::adjustedNodes(dag)) || length(selected_nodes(dag)) ){
    return(fallback("Adjusted or selected statuses require a conditional projection. The temporal display was retained."))
  }
  if( any(grepl("_(pre|post)_treatment$", all_nodes)) ){
    return(fallback("Numeric time occasions are required to recognise complete repeated passes. The pre/post display was retained."))
  }
  treatment <- as.character(dagitty::exposures(dag))
  outcome <- as.character(dagitty::outcomes(dag))
  if( !length(treatment) || !length(outcome) ){
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
  if( is.null(attr(reference, "cf_node_map", exact = TRUE)) &&
    !is.null(.cd_record(dag)$node_map) ){
    attr(reference, "cf_node_map") <- .cd_node_map(dag)
  }
  trimmed <- if( length(omitted) ){
    collapse_nodes(reference, omitted)
  }else{
    reference
  }
  attr(trimmed, "cf_source_order") <- attr(reference, "cf_source_order", exact = TRUE)
  mapping <- data.frame(node = names(trimmed),
    quantity = variable_groups$quantity[match(names(trimmed), variable_groups$node)],
    node_id = unname(.gd_node_ids(trimmed)),
    stringsAsFactors = FALSE)
  # An intervention on an exact occasion cannot be replaced by an intervention on a
  # merged quantity. Equal adjustment results alone cannot certify that change.
  endpoint_ids <- unique(variable_groups$quantity[match(c(treatment, outcome), variable_groups$node)])
  repeated_endpoints <- intersect(endpoint_ids,
    names(which(table(mapping$quantity) > 1L)))
  if( length(repeated_endpoints) ){
    return(fallback(paste0("Time merging would combine distinct occasions of an endpoint: ",
      paste(repeated_endpoints, collapse = ", "),
      ". The temporal display was retained.")))
  }
  latent <- names(trimmed) %in% dagitty::latents(trimmed)
  mixed <- tapply(latent, mapping$quantity, function(statuses){
    return(length(unique(statuses)) > 1L)
  })
  if( any(mixed) ){
    return(fallback(paste0("Mixed observation status across time for ",
      paste(names(mixed)[mixed], collapse = ", "), ". The temporal display was retained.")))
  }
  treatment_ids <- unique(variable_groups$quantity[match(treatment, variable_groups$node)])
  outcome_ids <- unique(variable_groups$quantity[match(outcome, variable_groups$node)])
  conflict <- intersect(treatment_ids, outcome_ids)
  if( length(conflict) ){
    return(fallback(paste0("Time merging would combine treatment and outcome: ",
      paste(conflict, collapse = ", "), ". The temporal display was retained.")))
  }

  links <- .gd_edge_ledger(trimmed, variable_groups)
  directed <- links[links$e == "->", , drop = FALSE]
  process_key <- function(from, time_from, to, time_to){
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
  owners <- tapply(mapping$quantity, mapping$display, function(quantities){
    return(length(unique(quantities)))
  })
  if( any(owners > 1L) ){
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
  # earlier source occasions, then quantity order (recorded source order when
  # available, otherwise name order). An edge left omitted only because an
  # edge of the same tier came first is counted as a tie (segment_edge_ties)
  # and reported by the caller that shows the display.
  repeated_rows <- which(links$display_role == "between_repetitions")
  edge_ties <- 0L
  if( length(repeated_rows) ){
    tier <- paste(links$time_to[repeated_rows] - links$time_from[repeated_rows], links$time_from[repeated_rows])
    node_keys <- stats::setNames(.query_order_key(names(trimmed), source_order,
      unname(.gd_node_ids(trimmed))), names(trimmed))
    ordering <- order(links$time_to[repeated_rows] - links$time_from[repeated_rows],
      links$time_from[repeated_rows], node_keys[links$v[repeated_rows]], node_keys[links$w[repeated_rows]], method = "radix")
    repeated_rows <- repeated_rows[ordering]
    tier <- tier[ordering]
    for( level in unique(tier) ){
      before <- edges
      for( i in repeated_rows[tier == level] ){
        candidate <- data.frame(v = unname(segment_of[links$v[i]]), e = links$e[i],
          w = unname(segment_of[links$w[i]]), stringsAsFactors = FALSE)
        if( edge_key(candidate) %in% edge_key(edges) ){
          links$display_role[i] <- "between_segments"
        }else if( .gd_edges_acyclic(rbind(edges, candidate)) ){
          edges <- rbind(edges, candidate)
          links$display_role[i] <- "between_segments"
        }else if( .gd_edges_acyclic(rbind(before, candidate)) ){
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
  if( !isTRUE(dag_is_acyclic(display)) ){
    return(fallback("Unresolved directed cycles remain after compression of repeated passes. The temporal display was retained."))
  }
  if( !identical(.gd_endpoint_paths(dag), .gd_endpoint_paths(display)) ){
    return(fallback("Time merging changes a proper directed path from treatment to outcome. The temporal display was retained."))
  }
  # Dropping edges from repeated stages must not hide or invent a route role: a
  # display node (a quantity or a time segment) lies on a directed
  # route from treatment to outcome exactly when one of the occurrences it stands
  # for does in the temporal graph. Endpoints are never segmented.
  display_routes <- unique(c(treatment_ids, outcome_ids,
    intersect(dagitty::descendants(display, treatment_ids), dagitty::ancestors(display, outcome_ids))))
  source_routes <- unique(mapping$display[mapping$node %in% routes])
  if( !setequal(source_routes, display_routes) ){
    return(fallback("Time merging changes which quantities lie on a directed route from treatment to outcome. The temporal display was retained."))
  }
  # A merged node stands for every occurrence it absorbed. Each adjustment set
  # the display would report must remain valid on the source when read that way.
  if( !is.null(source_check) ){
    refusal <- source_check(.gd_canonical(display), mapping)
    if( !is.null(refusal) ){
      return(fallback(refusal))
    }
  }
  summary <- list(graph = .gd_canonical(display), temporal = .gd_canonical(trimmed),
    trim_ledger = data.frame(node = all_nodes, node_id = all_ids,
    route = all_nodes %in% routes, pinned = all_nodes %in% pinned,
    action = ifelse(all_nodes %in% retained, "retain", "fold"),
    stringsAsFactors = FALSE),
    node_map = mapping, edge_ledger = links, segments = .gd_segment_table(mapping, source_order),
    segment_alternatives = choice$alternatives, segment_edge_ties = edge_ties,
    trim_introduced_edges = .gd_new_edges(reference, trimmed),
    applied = TRUE, reason = NULL)

  return(summary)
}


## Adjustment search and validation ##
#' Enumerate a bounded prefix of exact minimal adjustment sets
#'
#' @param dag Graph with treatment and outcome tags.
#' @param effect `"total"` or `"direct"`.
#' @param max_sets Maximum retained sets, a positive whole number or `Inf`.
#' @return A list with `sets`, `complete`, `probe` and `limit`. One extra set is
#'   fetched to detect truncation. A capped prefix need not contain the globally
#'   smallest sets. `sets` is `list()` when none exist and `list(character())`
#'   when the empty set is valid.
#' @noRd
.gd_adjustment_sets <- function(dag, effect, max_sets = Inf){
  effect <- match.arg(effect, c("total", "direct"))
  if( !is.numeric(max_sets) || length(max_sets) != 1L || is.na(max_sets) ||
    max_sets < 1 || (is.finite(max_sets) &&
    (max_sets != floor(max_sets) || max_sets >= .Machine$integer.max)) ){
    stop("The private search limit must be Inf or a positive whole number below .Machine$integer.max.",
      call. = FALSE)
  }
  limit <- if( is.finite(max_sets) ){
    as.integer(max_sets) + 1L
  }else{
    Inf
  }
  dag <- dagitty::as.dagitty(dag)
  if( !length(dagitty::exposures(dag)) || !length(dagitty::outcomes(dag)) ){
    stop("Both treatment and outcome tags are required to enumerate adjustment sets.",
      call. = FALSE)
  }
  raw <- .gd_adjustment_engine(dag, effect, limit)
  sets <- unname(lapply(raw, function(x){
    return(as.character(unlist(x, use.names = FALSE)))
  }))
  complete <- length(sets) <= max_sets
  probe <- if( !complete ){
    sets[[as.integer(max_sets) + 1L]]
  }else{
    NULL
  }
  if( !complete ){
    sets <- sets[seq_len(as.integer(max_sets))]
  }
  listed <- list(sets = sets, complete = complete, probe = probe,
    limit = if( is.finite(max_sets) ){
      as.integer(max_sets)
    }else{
      Inf
    })

  return(listed)
}


#' Validate exact adjustment sets for one effect
#'
#' @param dag Graph with treatment and outcome tags.
#' @param sets List of character vectors containing exact temporal node names.
#' @param effect `"total"` or `"direct"`, each with its own adjustment criterion.
#' @return A logical vector with one value per set: `TRUE` if valid for this effect,
#'   and `FALSE` otherwise. This check does not establish minimality.
#' @noRd
.gd_validate_sets <- function(dag, sets, effect){
  effect <- match.arg(effect, c("total", "direct"))
  valid <- .gd_validator(dag)(sets, effect)

  return(valid)
}


#' Create an adjustment validator with reusable parsed graphs
#'
#' @param dag Graph with treatment and outcome tags, retained by the validator.
#' @return A function accepting `sets` (exact-node character vectors in a list)
#'   and `effect` (`"total"` or `"direct"`). It returns `TRUE` for each valid set
#'   and `FALSE` otherwise. Parsed V8 graph state is reused for each effect.
#'   Validity does not establish minimality.
#' @noRd
.gd_validator <- function(dag){
  force(dag)
  context <- NULL
  total <- NULL
  parsed <- FALSE
  direct <- FALSE
  validate <- function(sets, effect){
    effect <- match.arg(effect, c("total", "direct"))
    dag <- dagitty::as.dagitty(dag)
    if( !length(dagitty::exposures(dag)) || !length(dagitty::outcomes(dag)) ){
      stop("Both treatment and outcome tags are required to validate adjustment sets.",
        call. = FALSE)
    }
    if( !is.list(sets) || any(!vapply(sets, function(z){
      return(is.character(z) && !anyNA(z))
    }, logical(1))) ){
      stop("sets must be a list of character vectors without missing values.",
        call. = FALSE)
    }
    if( !length(sets) ){
      return(logical())
    }
    endpoints <- c(dagitty::exposures(dag), dagitty::outcomes(dag))
    latent <- dagitty::latents(dag)
    eligible <- unname(vapply(sets, function(z){
      return(all(z %in% names(dag)) && !any(z %in% endpoints) && !any(z %in% latent))
    }, logical(1)))
    script <- system.file("js", "dagitty-alg.js", package = "dagitty")
    if( effect == "total" ){
      if( !any(eligible) ){
        return(eligible)
      }
      if( is.null(total) ){
        total <<- utils::packageVersion("dagitty") == "0.3.4" && nzchar(script) &&
          dagitty::graphType(dag) %in% c("dag", "mag", "pdag", "pag")
      }
      if( !total ){
        return(unname(vapply(seq_along(sets), function(i){
          return(eligible[[i]] && isTRUE(dagitty::isAdjustmentSet(dag, sets[[i]])))
        }, logical(1))))
      }
      if( is.null(context) ){
        context <<- V8::v8()
        context$source(script)
      }
      if( !parsed ){

        restated <- dag
        dagitty::exposures(restated) <- dagitty::exposures(dag)
        dagitty::outcomes(restated) <- dagitty::outcomes(dag)
        context$assign("gd_total_text", as.character(restated))
        context$eval("var gd_total = DAGitty.GraphParser.parseGuess(gd_total_text);")
        parsed <<- TRUE
      }
      context$assign("gd_sets", lapply(unname(sets[eligible]), as.list))
      context$eval(paste0("var gd_valid = gd_sets.map(function(z) { ",
        "return DAGitty.GraphAnalyzer.isAdjustmentSet(gd_total, z); });"))
      eligible[eligible] <- vapply(context$get("gd_valid", simplifyVector = FALSE), isTRUE, logical(1))
      return(eligible)
    }
    if( utils::packageVersion("dagitty") != "0.3.4" ){
      stop("Direct set validation currently supports dagitty 0.3-4 only. ",
        "The installed version is ", as.character(utils::packageVersion("dagitty")),
        ".", call. = FALSE)
    }
    if( !nzchar(script) ){
      stop("The installed dagitty JavaScript engine is unavailable.", call. = FALSE)
    }
    if( is.null(context) ){
      context <<- V8::v8()
      context$source(script)
    }
    if( !direct ){
      context$eval(paste0(
        "var gd_supported = typeof DAGitty !== 'undefined' && ",
        "typeof DAGitty.GraphParser.parseGuess === 'function' && ",
        "typeof DAGitty.GraphTransformer.indirectGraph === 'function' && ",
        "typeof DAGitty.GraphAnalyzer.listMsasDirectEffect === 'function' && ",
        "typeof DAGitty.GraphAnalyzer.listMinimalSeparators === 'function' && ",
        "typeof DAGitty.GraphAnalyzer.dConnected === 'function' && ",
        "DAGitty.GraphAnalyzer.listMinimalSeparators.length === 4;"
        ))
      if( !isTRUE(context$get("gd_supported")) ){
        stop("The installed dagitty engine lacks the verified API for validating direct sets.",
          call. = FALSE)
      }
      context$assign("gd_graph_text", as.character(dag))
      context$eval(paste(
        "var gd_direct = (function() {",
        "  var g = DAGitty.GraphParser.parseGuess(gd_graph_text);",
        "  var analyzer = DAGitty.GraphAnalyzer;",
        "  var separators = analyzer.listMinimalSeparators;",
        "  var constraints = null;",
        # Read the required and forbidden variables from Dagitty's direct
        # effect routine without enumerating sets
        "  analyzer.listMinimalSeparators = function(graph, must, mustNot, max) {",
        "    constraints = {",
        "      must: must.map(function(v) { return v.id; }),",
        "      mustNot: mustNot.map(function(v) { return v.id; })",
        "    };",
        "    return [];",
        "  };",
        "  try {",
        "    analyzer.listMsasDirectEffect(g, [], []);",
        "  } finally {",
        "    analyzer.listMinimalSeparators = separators;",
        "  }",
        # Condition on the full candidate set. Checking only endpoint
        # ancestors could miss a collider opened by an additional variable
        # outside those ancestors.
        "  return {g: g, constraints: constraints,",
        "    indirect: constraints === null ? null : DAGitty.GraphTransformer.indirectGraph(g)};",
        "})();",
        sep = "\n"
        ))
      direct <<- TRUE
    }
    context$assign("gd_sets", lapply(sets, as.list))
    context$eval(paste(
      "var gd_valid = (function() {",
      "  var g = gd_direct.g;",
      "  var analyzer = DAGitty.GraphAnalyzer;",
      "  var constraints = gd_direct.constraints;",
      "  if (constraints === null) return gd_sets.map(function() { return false; });",
      "  var indirect = gd_direct.indirect;",
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
    accepted <- eligible & valid

    return(accepted)
  }

  return(validate)
}


#' Remove nested construction provenance
#'
#' @param graph Graph to store inside a construction record.
#' @return The graph without `cf_provenance` or the `cf_dag` class.
#' @noRd
.gd_plain <- function(graph){
  attr(graph, "cf_provenance") <- NULL
  class(graph) <- setdiff(class(graph), "cf_dag")
  return(graph)
}


#' Summarise structural graph identity
#'
#' @param graph Graph whose structure and declarations are compared.
#' @return A list of graph type, sorted nodes, canonical typed edges, statuses
#'   and input labels. Coordinates and construction records are excluded.
#' @noRd
.gd_signature <- function(graph){
  edges <- dag_edges(graph)
  symmetric <- edges$e %in% c("<->", "--")
  if( any(symmetric) ){
    left_endpoint <- pmin(edges$v[symmetric], edges$w[symmetric])
    right_endpoint <- pmax(edges$v[symmetric], edges$w[symmetric])
    edges$v[symmetric] <- left_endpoint
    edges$w[symmetric] <- right_endpoint
  }
  labels <- attr(graph, "input_labels")
  if( length(labels) ){
    labels <- labels[order(names(labels))]
  }
  signature <- list(type = dag_graph_type(graph), nodes = sort(names(graph)),
    edges = sort(unique(edge_key(edges))),
    exposures = sort(as.character(dagitty::exposures(graph))),
    outcomes = sort(as.character(dagitty::outcomes(graph))),
    latent = sort(as.character(dagitty::latents(graph))),
    adjusted = sort(as.character(dagitty::adjustedNodes(graph))),
    selected = sort(as.character(selected_nodes(graph))), labels = labels)

  return(signature)
}


#' Require a saved construction record
#'
#' @param dag Graph expected to have construction provenance.
#' @return The provenance list. An absent or invalid record raises an error.
#' @noRd
.gd_provenance <- function(dag){
  state <- attr(dag, "cf_provenance", exact = TRUE)
  if( is.null(state) || !is.list(state) || is.null(state$original_query) ){
    stop("This DAG has no saved construction record. Run dagby_construct() on a query result or DAG first.", call. = FALSE)
  }
  return(state)
}


#' Validate a constructed graph against its saved structure
#'
#' @param dag A dagitty graph with construction provenance.
#' @return Its provenance after checking the structural signature. Structural
#'   edits raise an error. Changes to coordinates alone are allowed.
#' @noRd
.gd_state <- function(dag){
  if( !inherits(dag, "dagitty") ){
    stop("Supply a DAG.", call. = FALSE)
  }
  state <- .gd_provenance(dag)
  if( !identical(.gd_signature(dag), state$signature) ){
    stop("The DAG has changed and no longer matches its saved reference. Changes to coordinates alone are allowed.", call. = FALSE)
  }
  return(state)
}


#' Attach construction provenance to a graph
#'
#' @param graph A dagitty graph to tag as constructed.
#' @param state Construction record to save, with its signature refreshed.
#' @return The graph with new provenance and `cf_dag` class. Provenance from any
#'   earlier construction is removed first.
#' @noRd
.gd_attach <- function(graph, state){
  if( !inherits(graph, "dagitty") ){
    stop("graph must inherit dagitty.", call. = FALSE)
  }
  graph <- .gd_plain(graph)
  state$signature <- .gd_signature(graph)
  attr(graph, "cf_provenance") <- state
  class(graph) <- unique(c("cf_dag", class(graph)))
  return(graph)
}


#' Record structural changes between graphs
#'
#' @param before,after Graphs before and after a transformation.
#' @return A list of added and removed nodes and typed edges, plus changed status
#'   memberships. Coordinates and unchanged status groups are omitted.
#' @noRd
.gd_delta <- function(before, after){
  before_sig <- .gd_signature(before)
  after_sig <- .gd_signature(after)
  statuses <- c("exposures", "outcomes", "latent", "adjusted", "selected")
  status_changes <- lapply(statuses, function(name){
    return(list(
      added = setdiff(after_sig[[name]], before_sig[[name]]),
      removed = setdiff(before_sig[[name]], after_sig[[name]])))
  })
  names(status_changes) <- statuses
  status_changes <- status_changes[vapply(status_changes,
    function(changes){
      return(length(changes$added) + length(changes$removed) > 0L)
    }, logical(1))]
  delta <- list(added_nodes = setdiff(names(after), names(before)),
    removed_nodes = setdiff(names(before), names(after)),
    added_edges = .gd_new_edges(before, after),
    removed_edges = .gd_new_edges(after, before),
    statuses = status_changes)

  return(delta)
}


#' Group exact adjustment sets by source-variable membership
#'
#' @param raw List of exact-node adjustment sets, including an empty set if valid.
#' @param dag Graph providing canonical quantity and member identities.
#' @param source_order Optional recorded order covering all represented members.
#' @return A candidate data frame with IDs, size groups, labels and list columns
#'   for concepts, resolved members and their exact raw sets. Rows rank by member
#'   count, then recorded or label order. No sets gives an empty table.
#' @noRd
.gd_candidates <- function(raw, dag, source_order = NULL){
  mapping <- .cd_node_map(dag)
  ids_for <- function(nodes, column){
    return(unique(as.character(mapping[[column]][mapping$node %in% nodes])))
  }
  if( !is.null(source_order) ){
    if( !.query_order_applies(c(mapping$quantity, mapping$member_id), source_order) ){
      source_order <- NULL
    }
  }
  concepts <- lapply(raw, function(adjustment_set){
    return(.query_sort(ids_for(adjustment_set, "quantity"), source_order))
  })
  members <- lapply(raw, function(adjustment_set){
    return(.query_sort(ids_for(adjustment_set, "member_id"), source_order))
  })
  keys <- vapply(members, function(member_ids){
    return(paste0(nchar(member_ids, type = "bytes"), ":", member_ids, collapse = ""))
  }, character(1))
  groups <- unname(split(seq_along(keys), match(keys, unique(keys))))
  catalogue <- data.frame(id = character(), group = integer(), option = integer(), group_size = integer(),
    n_concepts = integer(), n_resolved = integer(),
    label = character(), stringsAsFactors = FALSE)
  catalogue$concepts <- list()
  catalogue$raw_sets <- list()
  catalogue$resolved <- list()
  rows <- vector("list", length(groups))
  for( index in seq_along(groups) ){
    source_indices <- groups[[index]]
    ids <- concepts[[source_indices[1L]]]
    resolved <- members[[source_indices[1L]]]
    row <- data.frame(id = "", group = 0L, option = 0L, group_size = NA_integer_, n_concepts = length(ids),
      n_resolved = length(resolved),
      label = if( length(ids) ){
        paste(ids, collapse = " + ")
      }else{
        "No adjustment nodes"
      }, stringsAsFactors = FALSE)
    row$concepts <- list(ids)
    row$raw_sets <- list(raw[source_indices])
    row$resolved <- list(resolved)
    rows[[index]] <- row
  }
  if( length(rows) ){
    catalogue <- do.call(rbind, rows)
  }
  if( nrow(catalogue) ){
    catalogue_keys <- catalogue$label
    if( !is.null(source_order) ){
      catalogue_keys <- vapply(catalogue$resolved, .query_set_key, character(1), source_order = source_order)
    }
    member_keys <- vapply(catalogue$resolved, function(member_ids){
      return(paste0(nchar(member_ids, type = "bytes"), ":", member_ids, collapse = ""))
    }, character(1))
    catalogue <- catalogue[order(catalogue$n_resolved, catalogue_keys, member_keys, method = "radix"), , drop = FALSE]
    catalogue$group <- match(catalogue$n_resolved, unique(catalogue$n_resolved))
    catalogue$option <- as.integer(ave(seq_len(nrow(catalogue)), catalogue$group, FUN = seq_along))
    catalogue$group_size <- catalogue$n_resolved
    catalogue$id <- paste0("set_", seq_len(nrow(catalogue)))
    rownames(catalogue) <- NULL
  }
  return(catalogue)
}
