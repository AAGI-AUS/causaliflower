# ==============================================================================
# dag_representations.R


# DAG sizes of a query
# ==============================================================================


## Changing the size ##
#' Change the size of a constructed DAG
#'
#' `prune_dag()` returns a smaller DAG and `grow_dag()` a larger one. Both rebuild
#' from the saved query.
#'
#' @param dag A DAG returned by `get_dag()`, `construct_dag()`, `prune_dag()`
#'   or `grow_dag()`, or a dry-run plan. A dagitty DAG is also accepted: its
#'   latent nodes are treated as unmeasured.
#' @param size The number of nodes by which to reduce or increase the size.
#'   The nearest available size in that direction is returned. A tie gives the
#'   larger change. `Inf` gives the smallest DAG for `prune_dag()` or the largest
#'   for `grow_dag()`. Nodes are counted as described in [construct_dag()].
#' @param quiet `NULL` prints progress and how to change the size. `TRUE`
#'   suppresses them. `FALSE` also prints the construction record.
#' @param unroll `TRUE` keeps separate time occasions. `FALSE` requests a checked
#'   static view. `NULL` keeps the current view.
#'
#' @details
#' The sizes are those of `construct_dag()`.
#'
#' @return A dagitty DAG. If the size does not change, the DAG is returned
#'   unchanged unless a different time view is requested.
#' @examples
#' \dontrun{
#' query <- query_graph(crop = "barley", treatments = "cultivar",
#'                      outcomes = "grain_yield", randomised_treatment = TRUE)
#' largest <- construct_dag(query)
#' smallest <- prune_dag(largest)
#' larger <- grow_dag(smallest, 5)
#' }
#' @export
prune_dag <- function(dag, size = Inf, quiet = NULL, unroll = NULL){
  return(.cd_resize(dag, size, quiet, "prune", unroll))
}


#' @rdname prune_dag
#' @export
grow_dag <- function(dag, size = Inf, quiet = NULL, unroll = NULL){
  return(.cd_resize(dag, size, quiet, "grow", unroll))
}


#' Select a DAG at a given size
#'
#' Choose another DAG from the same query at the current or requested size.
#' @param dag A constructed DAG or construction plan.
#' @param size The requested number of nodes. `NULL` keeps the current size.
#'   Counting and nearest-size rules follow [construct_dag()].
#' @param select The alternative number listed by `dag_options()`. `NULL` keeps
#'   the current candidate if it is available at the returned size.
#' @param unroll `TRUE` keeps separate time occasions. `FALSE` requests a checked
#'   static view. `NULL` keeps the current view.
#' @param quiet Suppress routine guidance when `TRUE`.
#' @details
#' Alternatives are numbered separately at each size in this order:
#' 1. The candidate using all construction rules.
#' 2. The candidate containing the union of minimal adjustment set members.
#' 3. Ranked candidates in increasing order of certified set size.
#' 4. Candidates from the remaining construction rules in rule order.
#'
#' Identical temporal graphs are listed once. Saved selections require the
#' same query and construction settings.
#' @return A constructed DAG.
#' @export
select_dag <- function(dag, size = NULL, select = NULL, unroll = NULL, quiet = NULL){
  input_record <- .cd_input_record(dag)
  if( inherits(dag, "cf_dag_plan") ){
    dag <- dag$graph
  }
  if( !inherits(dag, "cf_dag") ){
    stop("Call construct_dag() before selecting an alternative.", call. = FALSE)
  }
  state <- .gd_state(dag)
  current <- .cd_variable_count(state$temporal)
  if( is.null(size) ){
    size <- current
  }
  if( is.null(unroll) ){
    unroll <- .cd_output_mode(dag)
  }
  args <- state$construction$arguments
  result <- .cd_sized(state$original_query, args$latents, args$keep_nodes, FALSE, quiet, size,
    request = "select", current = current, unchanged = dag, select = select, unroll = unroll,
    prefer_current = is.null(select))
  return(.cd_with_input_record(result, input_record))
}


#' Rebuild a graph at a larger or smaller size
#'
#' @inheritParams prune_dag
#' @param request `"prune"` or `"grow"`.
#' @return A constructed DAG in the requested direction, or the unchanged DAG
#'   when none is available. A changed time view can still trigger rebuilding.
#' @noRd
.cd_resize <- function(dag, size, quiet, request, unroll = NULL){
  input_record <- .cd_input_record(dag)
  size <- .cd_size_value(size)
  if( !is.null(quiet) ){
    .gd_flag(quiet, "quiet")
  }
  if( inherits(dag, "cf_dag_plan") ){
    dag <- dag$graph
  }
  if( !inherits(dag, "dagitty") ){
    stop("Supply a DAG.", call. = FALSE)
  }
  .gd_flag(unroll, "unroll", nullable = TRUE)
  if( is.null(attr(dag, "cf_provenance", exact = TRUE)) ){
    warning(request, "_dag() is intended for working with DAGs created by querying the process-based structural causal model. It may have unexpected results when used with knowledge-based DAGs.", call. = FALSE)
    graph <- dag
    latents <- as.character(dagitty::latents(dag))
    if( !length(latents) ){
      latents <- NULL
    }
    keep_nodes <- NULL
    if( is.null(unroll) ){
      unroll <- FALSE
    }
  }else{
    state <- .gd_state(dag)
    arguments <- state$construction$arguments
    graph <- state$original_query
    latents <- arguments$latents
    keep_nodes <- arguments$keep_nodes
    if( is.null(unroll) ){
      unroll <- .cd_output_mode(dag)
    }
  }
  if( !size ){
    if( inherits(dag, "cf_dag") && !identical(unroll, .cd_output_mode(dag)) ){
      return(.cd_with_input_record(select_dag(dag, unroll = unroll, quiet = quiet), input_record))
    }
    return(.cd_with_input_record(dag, input_record))
  }
  current <- if( inherits(dag, "cf_dag") ) .cd_variable_count(state$temporal) else .cd_variable_count(dag)
  resized <- .cd_sized(graph, latents, keep_nodes, FALSE, quiet, size, request, current, NULL, dag,
    unroll = unroll)

  return(.cd_with_input_record(resized, input_record))
}


#' Explain an unavailable size change
#'
#' @param request `"prune"` or `"grow"`.
#' @return A character string explaining that no size in that direction is available.
#' @noRd
.cd_unchanged <- function(request){
  return(paste0("No ", if( request == "prune" ){
    "smaller"
  }else{
    "larger"
  }, " DAG is available. The DAG was returned unchanged."))
}


## Sizes on offer ##
#' Validate a size request
#'
#' @param size One nonnegative numeric value, including `Inf`.
#' @return `size` as a numeric scalar. Missing, negative or nonscalar values raise
#'   an error.
#' @noRd
.cd_size_value <- function(size){
  if( !is.numeric(size) || length(size) != 1L || is.na(size) || size < 0 ){
    stop("size must be one number of at least 0, or Inf.", call. = FALSE)
  }
  size_value <- as.numeric(size)

  return(size_value)
}


#' Select adjustment families at a ranked size threshold
#'
#' @param definitions Ranked levels and effect families from `.cd_level_definitions()`.
#' @param index Position of the chosen level.
#' @return A list named by effect, with each entry recording state, members, exact
#'   sets and an empty `undecided` vector. Certified nonempty families also record
#'   their set-size cap.
#' @noRd
.cd_level_found <- function(definitions, index){
  level <- definitions$levels[[index]]
  return(lapply(definitions$effects, function(effect){
    set_size <- max(level$cap, effect$smallest, na.rm = TRUE)
    kept <- effect$sizes <= set_size
    family <- list(state = effect$state,
      members = as.character(unique(unlist(effect$members[kept], use.names = FALSE))),
      undecided = character(), sets = unlist(effect$sets[kept], recursive = FALSE))
    if( effect$state == "sets" && effect$certified ){
      family$set_size <- set_size
    }
    return(family)
  }))
}


#' Describe sizes offered by ranked and rule candidates
#'
#' @param query Observed preparation and reuse key from `.cd_size_query()`.
#' @param rules Optional rule candidates from `.cd_rule_rows()`.
#' @param definitions Optional ranked levels from `.cd_level_definitions()`.
#' @return A list of candidate rows. Each row records kind, index, cap, rule,
#'   certification, retained nodes and distinct variable count. Candidates using
#'   all construction rules are flagged by `full`.
#' @noRd
.cd_choice_rows <- function(query, rules = NULL, definitions = NULL){
  variable_map <- NULL
  count_variables <- function(nodes){
    if( is.null(variable_map) ){
      variable_map <<- .cd_variable_map(query$observed$analysis)
    }
    if( !all(nodes %in% variable_map$node) ){
      stop("A retained node has no recorded variable identity.", call. = FALSE)
    }
    return(length(unique(variable_map$variable[variable_map$node %in% nodes])))
  }
  rows <- list()
  if( !is.null(definitions) ){
    for( index in seq_along(definitions$levels) ){
      found <- .cd_level_found(definitions, index)
      pins <- .cd_level_pins(query$observed, found)$nodes
      rows[[length(rows) + 1L]] <- list(kind = "ranked", index = index,
        cap = definitions$levels[[index]]$cap, rule = 0L,
        certified = definitions$levels[[index]]$certified, pins = pins,
        variables = count_variables(pins))
    }
  }
  if( !is.null(rules) ){
    pins <- .cd_level_pins(query$observed, rules$found)$nodes
    for( index in seq_along(rules$rows) ){
      nodes <- union(pins, rules$rows[[index]]$extra)
      rows[[length(rows) + 1L]] <- list(kind = "rule", index = index,
        cap = Inf, rule = rules$rows[[index]]$rule, certified = rules$certified,
        pins = nodes, variables = count_variables(nodes),
        full = isTRUE(rules$rows[[index]]$top))
    }
  }
  return(rows)
}


#' Choose the nearest eligible variable count
#'
#' @param rows Candidate rows with a `variables` count.
#' @param target Desired absolute variable count.
#' @param request Size, select, prune or grow operation.
#' @param current Current count, used to restrict directional changes.
#' @return A list with selected `count` and any tied alternative count, or `NULL`
#'   when none is eligible. Directional ties favour the larger change.
#' @noRd
.cd_choice_count <- function(rows, target, request, current){
  counts <- sort(unique(vapply(rows, function(row) row$variables, integer(1))))
  if( request == "prune" ){
    counts <- counts[counts < current]
  }else if( request == "grow" ){
    counts <- counts[counts > current]
  }
  if( !length(counts) ){
    return(NULL)
  }
  if( target >= max(counts) || target <= min(counts) ){
    return(list(count = if( target >= max(counts) ) max(counts) else min(counts), tie = numeric()))
  }
  distance <- abs(counts - target)
  nearest <- counts[distance == min(distance)]
  selected <- if( request == "prune" ) min(nearest) else max(nearest)
  return(list(count = selected, tie = setdiff(nearest, selected)))
}


#' Validate a ranked prefix before using its size bound
#'
#' @param query Observed preparation and exact reuse key.
#' @param definitions Certified ranked prefix from `.cd_level_definitions()`.
#' @param map Map of nodes to variables in the observed analysis graph.
#' @return A list with candidate `rows` and a `lower` bound on variable count,
#'   or `NULL` if an exact set fails minimality. Reuse requires unchanged helper
#'   implementations as well as unchanged inputs.
#' @noRd
.cd_prefix_proof <- function(query, definitions, map){
  key <- query$key
  implementation <- list()
  namespace <- asNamespace("causaliflower")
  pending <- c(".cd_prefix_proof", ".cd_choice_rows", ".cd_level_found", ".cd_maximal_minimal_sets", ".gd_validator")
  while( !is.null(key) && length(pending) ){
    name <- pending[[1L]]
    pending <- pending[-1L]
    if( name %in% names(implementation) ){
      next
    }
    fun <- get(name, envir = namespace, inherits = FALSE)
    if( !identical(environment(fun), namespace) ){
      key <- NULL
      break
    }
    implementation[[name]] <- list(formals = deparse(formals(fun)), body = deparse(body(fun)))
    references <- all.names(body(fun), functions = TRUE, unique = TRUE)
    references <- references[vapply(references, function(reference){
      return(is.function(get0(reference, envir = namespace, inherits = FALSE)))
    }, logical(1))]
    pending <- union(pending, setdiff(references, names(implementation)))
  }
  .cd_kept(key, list("prefix-proof", 200L, definitions, implementation,
    query$observed$selection, query$observed$analysis, query$observed$required, query$observed$map, map), function(){
    rows <- .cd_choice_rows(query, definitions = definitions)
    found <- .cd_level_found(definitions, length(definitions$levels))
    validate <- .gd_validator(query$observed$selection)
    for( effect in names(found) ){
      if( !all(.cd_maximal_minimal_sets(query$observed$selection,
        found[[effect]]$sets, effect, validate)) ){
        return(NULL)
      }
    }
    members <- unique(unlist(lapply(found, `[[`, "members"), use.names = FALSE))
    pins <- union(query$observed$required, map$node[map$variable %in% members])
    list(rows = rows, lower = length(unique(map$variable[map$node %in% pins])))
  })
}


#' Try a short ranked search with a proven nearest size
#'
#' @param query Observed preparation and reuse key.
#' @param budget Configured ranked search budget. The trial uses 200 steps.
#' @param target Desired absolute variable count.
#' @param request Size, select, prune or grow operation.
#' @param current Current count for directional requests.
#' @param quick Whether the first ranked size supplies the target.
#' @param retain `TRUE` to return a usable prefix even if it does not yet prove
#'   which candidate is nearest to the requested size.
#' @return A list with prefix `definitions`, candidate `rows`, a `closed` flag
#'   indicating that later candidates cannot change the selected size, and a
#'   `complete` flag. Returns `NULL` if the short search cannot justify reuse
#'   for this request.
#' @noRd
.cd_choice_prefix <- function(query, budget, target, request, current, quick, retain = FALSE){
  if( is.null(budget) || budget <= 200L || is.infinite(target) && target > 0 ){
    return(NULL)
  }
  captured <- .cd_captured(suppressMessages(.cd_level_definitions(query, 200L)))
  definitions <- captured$value
  if( length(captured$warnings) || !definitions$certified ||
    any(!vapply(definitions$effects, function(effect){
      return(effect$state %in% c("sets", "empty_set") && effect$certified &&
        !identical(effect$ranked$certificate, "validation"))
    }, logical(1))) ){
    return(NULL)
  }
  map <- .cd_variable_map(query$observed$analysis)
  if( anyDuplicated(map$node) ){
    return(NULL)
  }
  proof <- .cd_prefix_proof(query, definitions, map)
  if( is.null(proof) ){
    return(NULL)
  }
  rows <- proof$rows
  if( quick ){
    target <- rows[[1L]]$variables
  }
  definitions$steps <- budget
  prefix <- list(definitions = definitions, rows = rows, closed = FALSE, complete = FALSE)
  chosen <- .cd_choice_count(rows, target, request, current)
  if( is.null(chosen) || length(chosen$tie) ){
    return(if( retain ) prefix else NULL)
  }
  counts <- vapply(rows, function(row) row$variables, integer(1))
  larger <- counts[counts > chosen$count]
  if( !length(larger) ){
    return(if( retain ) prefix else NULL)
  }
  # Later ranked levels and rule levels retain these minimal set members
  if( proof$lower < min(larger) ){
    return(if( retain ) prefix else NULL)
  }
  prefix$closed <- TRUE
  return(prefix)
}


#' Bound variable counts of omitted ranked candidates
#'
#' @param query Observed preparation and reuse key.
#' @param definitions Ranked prefix and effect certificates.
#' @return One lower bound per effect, adding required variables outside the
#'   removable vocabulary. `Inf` means no omitted candidate for that effect.
#' @noRd
.cd_prefix_bounds <- function(query, definitions){
  map <- .cd_variable_map(query$observed$analysis)
  required <- unique(map$variable[map$node %in% query$observed$required])
  cap <- definitions$levels[[length(definitions$levels)]]$cap
  return(vapply(names(definitions$effects), function(effect){
    family <- definitions$effects[[effect]]
    omitted <- family$sizes[family$sizes > max(cap, family$smallest)]
    size <- min(c(Inf, omitted, if( !family$complete ) family$certified_to + 1))
    if( is.infinite(size) ){
      return(Inf)
    }
    region <- if( effect == "total" ) .cd_ranked_total_region(query$observed$selection) else{
      .cd_ranked_direct_region(query$observed$selection)
    }
    members <- .cd_ranked_labels(query$observed, region$graph$nodes)
    members[region$graph$synthetic] <- list(character())
    vocabulary <- unique(unlist(members[region$removable], use.names = FALSE))
    return(size + length(setdiff(required, vocabulary)))
  }, numeric(1)))
}


#' Close a ranked prefix using rule candidates and omission bounds
#'
#' @param query Observed preparation and reuse key.
#' @param prefix Short ranked prefix from `.cd_choice_prefix()`.
#' @param rules Rule candidates from `.cd_rule_rows()`.
#' @param target Desired absolute variable count.
#' @param request Size, select, prune or grow operation.
#' @param current Current count for directional requests.
#' @param quick Whether the first ranked size supplies the target.
#' @return The prefix with combined candidate rows, `closed = TRUE`, an omission
#'   bound `lower` and a `complete` flag. Returns `NULL` if omitted candidates
#'   could change the selected size.
#' @noRd
.cd_choice_extend <- function(query, prefix, rules, target, request, current, quick){
  rows <- .cd_choice_rows(query, rules, prefix$definitions)
  if( quick ){
    target <- rows[[1L]]$variables
  }
  chosen <- .cd_choice_count(rows, target, request, current)
  lower <- min(.cd_prefix_bounds(query, prefix$definitions))
  if( is.null(chosen) || length(chosen$tie) || chosen$count >= lower ){
    return(NULL)
  }
  if( is.finite(lower) ){
    possible <- c(rows, list(list(variables = as.integer(lower))))
    if( !identical(chosen, .cd_choice_count(possible, target, request, current)) ){
      return(NULL)
    }
  }
  larger <- vapply(rows, function(row) row$variables, integer(1))
  larger <- larger[larger > chosen$count]
  if( is.finite(lower) && target > chosen$count && (!length(larger) || min(larger) > lower) ){
    return(NULL)
  }
  prefix$rows <- rows
  prefix$closed <- TRUE
  prefix$complete <- is.infinite(lower)
  prefix$lower <- lower
  return(prefix)
}


#' Encode retained membership without delimiter collisions
#'
#' @param nodes Character member names whose order and duplication are irrelevant.
#' @return One character key with a version prefix and members in radix order.
#'   Each member is preceded by its byte length to prevent delimiter collisions.
#' @noRd
.cd_choice_key <- function(nodes){
  nodes <- sort(unique(nodes), method = "radix")
  return(paste0("variables-v1:", paste0(nchar(nodes, type = "bytes"), ":", nodes, collapse = "")))
}


#' Bind a candidate to its temporal graph and identities
#'
#' @param dag An unchanged constructed DAG.
#' @return Compressed raw serialisation of the temporal structural signature
#'   and complete occurrence map, used to deduplicate candidates.
#' @noRd
.cd_choice_signature <- function(dag){
  state <- .gd_state(dag)
  return(memCompress(serialize(list(graph = .gd_signature(state$temporal),
    map = .cd_node_map(state$temporal)), NULL, version = 3), "gzip"))
}


#' Recover the selected time representation
#'
#' @param dag An unchanged constructed DAG.
#' @return `TRUE` if the saved `unroll` argument is `TRUE`, and `FALSE` if it is
#'   present with any other value. Without a saved argument, returns `TRUE` if
#'   no summary was applied and the temporal graph has repeated variable
#'   membership, and `FALSE` otherwise.
#' @noRd
.cd_output_mode <- function(dag){
  state <- .gd_state(dag)
  mode <- state$construction$arguments$unroll
  if( !is.null(mode) ){
    return(isTRUE(mode))
  }
  map <- .cd_variable_map(state$temporal)
  return(!isTRUE(state$construction$summary$applied) && any(table(map$variable) > 1L))
}


#' Bind saved candidate selection to source and construction settings
#'
#' @param graph Original supplied query graph.
#' @param query Observed preparation, exact reuse key and optional saved signatures.
#' @param latents,keep_nodes Original unmeasured and retention requests.
#' @param budget Configured ranked-search step budget.
#' @return A list containing the construction version, compressed source records
#'   and the ranked and membership search budgets. Saved signatures are reused
#'   only when graph text, class, labels and locale match.
#' @noRd
.cd_choice_binding <- function(graph, query, latents, keep_nodes, budget){
  signature <- function(name){
    dag <- query$observed[[name]]
    saved <- query$signatures[[name]]
    if( !is.null(saved$value) && identical(query$signatures$locale, .cf_cache_locale()) &&
      identical(saved$class, class(dag)) && identical(saved$text, as.character(dag)) &&
      identical(saved$labels, attr(dag, "input_labels")) ){
      return(saved$value)
    }
    return(.gd_signature(dag))
  }
  record <- .cd_record(graph)
  source_order <- attr(query$observed$analysis, "cf_source_order", exact = TRUE)
  if( is.null(source_order) ){
    source_order <- record$source_order
  }
  query_depth <- record$arguments$params$query_depth
  if( is.null(query_depth) ){
    query_depth <- attr(graph, "query_depth", exact = TRUE)
  }
  source <- list(
    original = .gd_signature(graph), original_map = .cd_node_map(graph),
    original_record = record[intersect(c("version", "input", "requests", "source_node_map",
      "recovery_folding", "scope", "aliases", "policy", "policy_scope", "policy_applied",
      "time_unroll_defaulted"), names(record))],
    source = signature("analysis"), map = .cd_variable_map(query$observed$analysis),
    selection = signature("selection"), policy = query$observed$ledger,
    source_policy = query$observed$record$source_policy, arguments = query$observed$record$arguments,
    source_order = source_order, query_depth = query_depth,
    controls = query$observed$controls, keep_nodes = keep_nodes, latents = latents)
  return(list(version = .cd_construction_version,
    source = memCompress(serialize(source, NULL, version = 3), "gzip"),
    ranked_steps = budget, member_steps = .cd_maximal_steps()))
}


#' Construct and select a source-bound DAG at the nearest size
#'
#' @inheritParams construct_dag
#' @param request Size, select, prune or grow operation.
#' @param size Requested distinct variable count for `"size"` or `"select"`, or change
#'   for `"prune"` or `"grow"`. `0` and `Inf` follow the operation's boundary rules.
#' @param current Current variable count for directional requests.
#' @param known Reserved previous-size input, currently unused.
#' @param unchanged Input DAG to return when no directional change is available.
#' @param select Alternative number at the selected size, or `NULL` for the default.
#' @param prefer_current Whether to retain the current candidate when available.
#' @param expected Saved candidate and binding required for exact replay.
#' @param budget Ranked search step budget, or `NULL` for the package default.
#' @param quick Whether the first ranked size supplies the initial target.
#' @param requested_size Original size request recorded for guidance.
#' @return A constructed DAG, or a dry-run plan when `dry_run = TRUE`, with
#'   candidate, size, alternatives and completeness records. Returns `unchanged`
#'   if no directional change is available. An invalid replay binding or an
#'   unavailable checked static representation raises an error.
#' @noRd
.cd_sized <- function(graph, latents, keep_nodes, dry_run, quiet, size, request = "size",
                      current = NA_integer_, known = NULL, unchanged = NULL,
                      select = NULL, unroll = FALSE, prefer_current = FALSE,
                      expected = NULL, budget = NULL, quick = FALSE, requested_size = size){
  .gd_flag(dry_run, "dry_run")
  .gd_flag(quiet, "quiet", nullable = TRUE)
  .gd_flag(unroll, "unroll")
  size <- .cd_size_value(size)
  invalid_select <- !is.null(select) && (!is.numeric(select) || length(select) != 1L || is.na(select) ||
    !is.finite(select) || select < 1 || select != floor(select))
  query <- .cd_size_query(graph, keep_nodes)
  if( length(dagitty::adjustedNodes(query$observed$selection)) || length(selected_nodes(query$observed$selection)) ){
    stop("The DAG has adjusted or selected nodes. Remove those tags first.", call. = FALSE)
  }
  if( is.null(budget) && is.null(expected) ){
    budget <- tryCatch(.cd_ranked_steps(), error = function(e) NULL)
  }
  target <- switch(request, size = size, select = size, prune = current - size, grow = current + size)
  prefix <- .cd_choice_prefix(query, budget, target, request, current,
    quick && is.na(requested_size), retain = TRUE)
  refused <- NULL
  rules <- if( is.null(prefix) || !prefix$closed ) tryCatch(.cd_rule_rows(query, quiet), cd_several_members = function(e){
    refused <<- e
    return(NULL)
  })
  if( is.infinite(target) && target > 0 && !is.null(refused) ){
    stop(refused)
  }
  if( !is.null(prefix) && !prefix$closed ){
    prefix <- if( !is.null(rules) ) .cd_choice_extend(query, prefix, rules, target,
      request, current, quick && is.na(requested_size))
  }
  rows <- if( is.null(prefix) ) .cd_choice_rows(query, rules) else prefix$rows
  rule_counts <- vapply(Filter(function(row) row$kind == "rule", rows),
    function(row) row$variables, integer(1))
  # Ranked candidates may also retain treatment parents and the first nodes of
  # indirect paths
  ranked_bound <- if( !is.null(rules) && length(rows) ){
    diagnostic <- lapply(rules$found, function(family){
      family$state <- "not_found"
      return(family)
    })
    .cd_variable_count(query$observed$analysis, .cd_level_pins(query$observed, diagnostic)$nodes)
  }else{
    Inf
  }
  definitions <- if( !is.null(prefix) ) prefix$definitions
  need_list <- is.null(prefix) && (!length(rows) || target <= ranked_bound || quick)
  read_list <- function(){
    if( is.null(budget) ){
      budget <<- .cd_ranked_steps()
    }
    definitions <<- .cd_level_definitions(query, budget)
    rows <<- .cd_choice_rows(query, rules, definitions)
  }
  if( need_list ){
    read_list()
  }
  if( quick && is.na(requested_size) ){
    target <- rows[[1L]]$variables
  }
  chosen <- .cd_choice_count(rows, target, request, current)
  if( is.null(definitions) && (is.null(chosen) || chosen$count <= ranked_bound ||
    is.finite(target) && abs(chosen$count - target) >= target - ranked_bound) ){
    read_list()
    chosen <- .cd_choice_count(rows, target, request, current)
  }
  if( is.null(chosen) ){
    if( !is.null(unchanged) && inherits(unchanged, "cf_dag") &&
      !identical(unroll, .cd_output_mode(unchanged)) ){
      return(select_dag(unchanged, unroll = unroll, quiet = quiet))
    }
    if( !isTRUE(quiet) ){
      message(.cd_unchanged(request))
    }
    return(unchanged)
  }
  binding <- .cd_choice_binding(graph, query, latents, keep_nodes, budget)
  if( !is.null(expected) && !identical(expected$binding, binding) ){
    stop("The saved DAG's source or construction settings no longer match. Rebuild the DAG.", call. = FALSE)
  }
  group <- Filter(function(row) row$variables == chosen$count, rows)
  full <- vapply(group, function(row) isTRUE(row$full), logical(1))
  ranked <- vapply(group, function(row) identical(row$kind, "ranked"), logical(1))
  union <- vapply(group, function(row) row$kind == "rule" && row$rule == 0L, logical(1))
  group <- group[order(!full, !union, !ranked, vapply(group, function(row) row$cap, numeric(1)),
    vapply(group, function(row) row$rule, integer(1)), method = "radix")]
  build <- function(row, output, silent, plan = FALSE){
    settings <- list(unroll = output)
    if( row$kind == "rule" ){
      return(.cd_maximal(graph, latents, keep_nodes, plan, silent, settings, query,
        rules, rules$rows[[row$index]]))
    }
    return(.cd_level(graph, latents, keep_nodes, query, definitions, row$index,
      plan, silent, settings, smallest = row$index == 1L))
  }
  signatures <- list()
  kept <- list()
  available <- logical()
  quick_row <- NULL
  quick_choice <- NULL
  for( row in group ){
    captured <- .cd_captured(suppressMessages(build(row, TRUE, TRUE)))
    signature <- .cd_choice_signature(captured$value)
    duplicate <- which(vapply(signatures, identical, logical(1), y = signature))
    if( row$kind == "ranked" && row$index == 1L ){
      quick_row <- row
      quick_choice <- if( length(duplicate) ) duplicate[[1L]] else length(kept) + 1L
    }
    if( length(duplicate) ){
      next
    }
    if( .cd_variable_count(temporal_dag(captured$value)) != row$variables ){
      stop("The constructed DAG has a different size from the one predicted.", call. = FALSE)
    }
    signatures[[length(signatures) + 1L]] <- signature
    kept[[length(kept) + 1L]] <- row
    available <- c(available, isTRUE(dag_options(captured$value)$level$static_available))
  }
  group <- kept
  if( invalid_select ){
    stop("select must be one positive whole number. Available alternatives at this size: ",
      paste(seq_along(group), collapse = ", "), ".", call. = FALSE)
  }
  selected <- 1L
  quick <- quick || !any(full) && !is.null(definitions) &&
    chosen$count == min(vapply(rows, function(row) row$variables, integer(1)))
  if( quick && !is.null(quick_choice) ){
    selected <- quick_choice
  }
  reference <- if( !is.null(expected) ) expected$signature else if( prefer_current && !is.null(unchanged) ){
    .cd_choice_signature(unchanged)
  }
  if( !is.null(reference) ){
    matches <- which(vapply(signatures, identical, logical(1), y = reference))
    if( length(matches) ){
      selected <- matches[[1L]]
    }else if( !is.null(expected) || prefer_current && chosen$count == current ){
      stop("The saved DAG does not exactly match an available alternative. Rebuild the DAG.", call. = FALSE)
    }
  }
  if( !is.null(select) ){
    if( select > length(group) ){
      stop("Available alternatives at this size: ", paste(seq_along(group), collapse = ", "), ".", call. = FALSE)
    }
    selected <- as.integer(select)
  }
  recipe <- if( quick && identical(selected, quick_choice) ) quick_row else group[[selected]]
  result <- tryCatch(build(recipe, unroll, quiet, dry_run), cd_static_unavailable = function(condition){
    retry <- switch(request,
      size = paste0("construct_dag(graph, size = ", chosen$count, ", unroll = TRUE)"),
      select = paste0("select_dag(dag, size = ", chosen$count, ", select = ", selected, ", unroll = TRUE)"),
      paste0(request, "_dag(dag, size = ", size, ", unroll = TRUE)"))
    stop(errorCondition(paste0(conditionMessage(condition), " Use ", retry, "."),
      class = "cd_static_unavailable"))
  })
  dag <- if( inherits(result, "cf_dag_plan") ) result$graph else result
  state <- .gd_state(dag)
  counts <- sort(unique(vapply(rows, function(row) row$variables, integer(1))))
  state$construction$metric <- "variables-v1"
  state$construction$level$variables <- chosen$count
  state$construction$level$smallest <- !is.null(definitions) && chosen$count == min(counts)
  state$construction$level$top <- !is.null(rules) && chosen$count == max(rule_counts)
  state$construction$level$select <- selected
  state$construction$arguments$size <- chosen$count
  state$construction$arguments$unroll <- unroll
  state$construction$candidate <- list(id = .cd_choice_key(group[[selected]]$pins), binding = binding,
    signature = signatures[[selected]])
  ids <- vapply(group, function(row) .cd_choice_key(row$pins), character(1))
  state$construction$alternatives <- data.frame(select = seq_along(group), variables = chosen$count,
    selected = seq_along(group) == selected, id = ids,
    description = vapply(group, function(row) if( row$kind == "rule" ){
      c("Adjustment members", "First and last route steps", "Paths from treatment to outcome",
        "First descendants outside paths from treatment to outcome", "Other treatment descendants", "Outside parents", "Relevant upstream variables")[[row$rule + 1L]]
    }else{
      paste("Adjustment DAG", row$index)
    }, character(1)), static_available = available, stringsAsFactors = FALSE)
  state$construction$alternatives_complete <- !is.null(rules) || !is.null(prefix)
  state$construction$alternatives_limit <- if( !is.null(refused) ) conditionMessage(refused)
  state$construction$sizes <- data.frame(variables = counts)
  state$construction$sizes_complete <- if( !is.null(prefix) ) prefix$complete else{
    !is.null(definitions) && !is.null(rules)
  }
  state$construction$sizes_dropped <- NULL
  state$construction$sizes_steps <- if( !is.null(definitions) ) budget else NULL
  dag <- .gd_attach(dag, state)
  if( inherits(result, "cf_dag_plan") ){
    result$graph <- dag
    result$options <- state$construction
    result$arguments <- c(list(graph = graph), state$construction$arguments,
      list(dry_run = TRUE, quiet = quiet))
  }else{
    result <- dag
  }
  if( length(chosen$tie) ){
    other <- chosen$tie[[1L]]
    warning("Sizes ", min(other, chosen$count), " and ", max(other, chosen$count),
      " are equally close to ", target, " nodes. Returned ", chosen$count, ". Use ",
      if( other < chosen$count ) "prune_dag(dag, " else "grow_dag(dag, ",
      abs(other - chosen$count), ") for ", other, ".", call. = FALSE)
  }
  if( !isTRUE(quiet) ){
    smaller <- counts[counts < chosen$count]
    larger <- counts[counts > chosen$count]
    if( is.null(definitions) && length(smaller) && max(smaller) < ranked_bound ){
      smaller <- integer()
    }
    if( !is.null(prefix$lower) && length(larger) && min(larger) > prefix$lower ){
      larger <- integer()
    }
    message(.cd_size_message(result, request, requested_size, current,
      if( length(smaller) ) max(smaller) else NA_integer_,
      if( length(larger) ) min(larger) else NA_integer_))
  }
  return(result)
}


#' Evaluate a candidate while retaining its warnings
#'
#' @param expr Unevaluated construction expression.
#' @return A list with its `value` and warning condition objects, muffled during
#'   evaluation so callers can report only warnings for the chosen candidate.
#' @noRd
.cd_captured <- function(expr){
  held <- list()
  value <- withCallingHandlers(expr, warning = function(condition){
    held[[length(held) + 1L]] <<- condition
    invokeRestart("muffleWarning")
  })
  captured <- list(value = value, warnings = held)

  return(captured)
}


#' Read a construction record from a DAG or plan
#'
#' @param result An unchanged constructed DAG or dry-run plan.
#' @return The saved construction list after validating its graph signature.
#' @noRd
.cd_construction <- function(result){
  if( inherits(result, "cf_dag_plan") ){
    result <- result$graph
  }
  construction <- .gd_state(result)$construction

  return(construction)
}


#' Describe the returned size and available changes
#'
#' @param result Constructed DAG or dry-run plan.
#' @param request Size, select, prune or grow operation.
#' @param size Original requested count or directional change.
#' @param current Previous variable count for directional requests.
#' @param smaller,larger Known neighbouring counts, or `NA` when unknown.
#' @return Guidance text including size limits, alternatives and any rule limitation.
#' @noRd
.cd_size_message <- function(result,
                             request = "size",
                             size = NA_real_,
                             current = NA_integer_,
                             smaller = NA_integer_,
                             larger = NA_integer_
                             ){
  level <- .cd_construction(result)$level
  nodes <- level$variables
  text <- paste0("DAG with ", nodes, " nodes", if( isTRUE(level$smallest) && isTRUE(level$top) ){
    ". This is the only size"
  }else if( isTRUE(level$smallest) ){
    ". This is the smallest size"
  }else if( isTRUE(level$top) ){
    ". This is the largest size"
  }, ".")
  changed <- abs(nodes - current)
  if( request %in% c("prune", "grow") && is.finite(size) && changed != size ){
    text <- paste0(text, " ", format(changed, scientific = 15), if( changed == 1 ){
      " node was "
    }else{
      " nodes were "
    }, if( request == "prune" ){
      "removed"
    }else{
      "added"
    }, ", not ", format(size, scientific = 15), ".")
  }
  if( request %in% c("size", "select") && is.finite(size) && size != nodes ){
    text <- paste0(text, " Requested ", size, " nodes. Returned ", nodes, ".")
  }
  if( !is.na(smaller) ){
    text <- paste0(text, " Next smaller size: ", smaller, " nodes. Use prune_dag(dag, ", nodes - smaller, ").")
  }else if( !isTRUE(level$smallest) ){
    text <- paste0(text, " Smaller, if any: prune_dag(dag, 1).")
  }
  if( !is.na(larger) ){
    text <- paste0(text, " Next larger size: ", larger, " nodes. Use grow_dag(dag, ", larger - nodes, ").")
  }else if( !isTRUE(level$top) ){
    text <- paste0(text, " Larger, if any: grow_dag(dag, 1).")
  }
  alternatives <- .cd_construction(result)$alternatives
  if( !is.null(alternatives) && nrow(alternatives) > 1L ){
    other <- alternatives$select[!alternatives$selected]
    text <- paste0(text, if( isTRUE(.cd_construction(result)$alternatives_complete) ){
      " Alternative "
    }else{
      " Known alternative "
    }, level$select, " of ", nrow(alternatives), ".\nOther alternatives:\n",
      paste0("select_dag(dag, select = ", other, ")", collapse = "\n"), ".")
  }
  if( !is.null(.cd_construction(result)$alternatives_limit) ){
    text <- paste0(text, " Some alternatives from the construction rules are unavailable. ",
      .cd_construction(result)$alternatives_limit)
  }
  return(text)
}


## Level definitions ##
#' Prepare a fully observed query for size discovery
#'
#' @param graph Supplied query graph, which may contain live-reference attributes.
#' @param keep_nodes Additional requested variables or exact occasions.
#' @return A list containing observed preparation and a cache key based on search
#'   identity (`NULL` when reuse is unsafe). It also contains reusable structural
#'   signatures if the locale remained unchanged while the search identity and
#'   structural signatures were computed.
#' @noRd
.cd_size_query <- function(graph, keep_nodes){
  observed <- .cd_prepare(graph, NULL, keep_nodes)
  locale <- .cf_cache_locale()
  identity <- .cd_search_identity(observed)
  query <- list(observed = observed, key = if( !is.null(identity) ){
    list(.cd_construction_version, identity, .gd_signature(observed$analysis), keep_nodes)
  })
  if( !is.null(query$key) && identical(locale, .cf_cache_locale()) ){
    signature <- function(dag, value){
      return(list(value = value, text = as.character(dag), class = class(dag),
        labels = attr(dag, "input_labels")))
    }
    query$signatures <- list(locale = locale,
      analysis = signature(observed$analysis, query$key[[3L]]),
      selection = signature(observed$selection, identity$signature))
  }

  return(query)
}


#' Obtain observed and declared-unmeasured construction references
#'
#' @param graph Supplied query graph.
#' @param latents,keep_nodes Unmeasured and retention requests.
#' @param query Optional observed preparation from `.cd_size_query()`.
#' @return A list containing `prepared`, with the requested unmeasured declarations,
#'   and `observed`, treating all variables as observed. The `key` is retained
#'   when the observed preparation matches the query, and is `NULL` otherwise.
#' @noRd
.cd_size_prepared <- function(graph, latents, keep_nodes, query = NULL){
  observed <- query$observed
  prepared <- if( is.null(latents) && !is.null(observed) ){
    observed
  }else{
    .cd_prepare(graph, latents, keep_nodes)
  }
  if( !length(prepared$latent_nodes) ){
    observed <- prepared
  }else if( is.null(observed) ){
    observed <- .cd_prepare(graph, NULL, keep_nodes)
  }
  references <- list(prepared = prepared, observed = observed, key = if( identical(observed, query$observed) ){
    query$key
  })

  return(references)
}


#' Reuse a size-discovery calculation for an exact query
#'
#' @param query Observed query with a reuse key. A `NULL` key disables caching.
#' @param key Additional identity of this level calculation.
#' @param compute Zero-argument function computing the value on a cache miss.
#' @return The cached or freshly computed value.
#' @noRd
.cd_level_cache <- function(query, key, compute){
  if( is.null(query$key) ){
    return(compute())
  }
  cached <- .cf_cache$call("levels", c(query$key, key), compute)

  return(cached)
}


#' Reuse a construction step that signals only progress
#'
#' @param key Exact query key, or `NULL` to disable caching.
#' @param step Additional identity of the build step.
#' @param compute Zero-argument function computing the step on a cache miss.
#' @return The cached or freshly computed value. Conditions other than construction progress
#'   prevent storage, so reuse cannot hide a warning or other diagnostic.
#' @noRd
.cd_kept <- function(key, step, compute){
  if( is.null(key) ){
    return(compute())
  }
  key <- c(key, step)
  if( .cf_cache$has("builds", key) ){
    return(.cf_cache$call("builds", key, function(){
      return(NULL)
    }))
  }
  plain <- TRUE
  value <- withCallingHandlers(compute(), condition = function(condition){
    if( !inherits(condition, "message") || !startsWith(conditionMessage(condition), "Construction [") ){
      plain <<- FALSE
    }
  })
  if( plain ){
    .cf_cache$call("builds", key, function(){
      return(value)
    })
  }
  return(value)
}


#' Define certified ranked levels for both effects
#'
#' @param query Fully observed preparation and reuse key.
#' @param steps Ranked-search budget for each effect.
#' @return A list containing effect families with exact sets, members and
#'   certificates, along with level caps and overall completeness. Levels selecting
#'   the same families share their highest cap. An uncertified smallest family
#'   permits only one level.
#' @noRd
.cd_level_definitions <- function(query, steps){
  return(.cd_level_cache(query, list("definitions", steps), function(){
    validate <- .gd_validator(query$observed$selection)
    effects <- stats::setNames(lapply(c("direct", "total"), function(effect){
      search <- .cd_ranked_search(query$observed, effect, steps = steps, validate = validate)
      table <- search$table
      sizes <- as.numeric(table$n_resolved)
      state <- if( !length(sizes) ){
        if( search$complete ){
          "no_sets"
        }else{
          "not_found"
        }
      }else if( !sizes[[1L]] ){
        "empty_set"
      }else{
        "sets"
      }
      certified_to <- if( search$complete ){
        Inf
      }else{
        search$ranked$complete_up_to
      }
      smallest <- if( length(sizes) ){
        min(sizes)
      }else{
        NA_real_
      }
      family <- list(state = state, sizes = sizes, members = lapply(table$resolved, as.character),
        sets = lapply(table$raw_sets, function(raw){
          return(lapply(raw, as.character))
        }), smallest = smallest, certified_to = certified_to,
        certified = state %in% c("no_sets", "empty_set") || (state == "sets" && certified_to >= smallest),
        complete = search$complete, ranked = search$ranked,
        omitted_equal_smallest = search$omitted_equal_smallest)

      return(family)
    }), c("direct", "total"))
    present <- Filter(function(effect){
      return(effect$state != "no_sets")
    }, effects)
    limit <- min(c(Inf, vapply(present, `[[`, numeric(1), "certified_to")))
    smallest <- suppressWarnings(min(vapply(effects, `[[`, numeric(1), "smallest"), na.rm = TRUE))
    if( !is.finite(smallest) ){
      smallest <- 0
    }
    certified <- all(vapply(effects, `[[`, logical(1), "certified"))
    thresholds <- smallest
    if( certified ){
      thresholds <- sort(unique(c(smallest, unlist(lapply(effects, `[[`, "sizes")))))
      thresholds <- thresholds[thresholds <= max(limit, smallest)]
    }
    levels <- list()
    previous <- NULL
    for( threshold in thresholds ){
      chosen <- lapply(effects, function(effect){
        return(which(effect$sizes <= max(threshold, effect$smallest, na.rm = TRUE)))
      })
      last <- length(levels)
      if( last && identical(chosen, previous) ){
        levels[[last]]$cap <- threshold
        next
      }
      kept <- unique(unlist(lapply(names(effects), function(effect){
        return(effects[[effect]]$members[chosen[[effect]]])
      }), use.names = FALSE))
      levels[[last + 1L]] <- list(cap = threshold, quantities = as.character(kept), certified = certified)
      previous <- chosen
    }
    definitions <- list(steps = steps, effects = effects, levels = levels, certified = certified, limit = limit,
      complete = certified && all(vapply(effects, `[[`, logical(1), "complete")))

    return(definitions)
  }))
}


#' Report the limit of an incomplete ranked family
#'
#' @param effect `"total"` or `"direct"`.
#' @param family Effect family carrying ranked certificate and omission details.
#' @return `NULL`, invisibly, after warning about failed validation or the certified
#'   member-count limit under the configured step budget.
#' @noRd
.cd_ranked_warning <- function(effect, family){
  certificate <- family$ranked
  if( identical(certificate$certificate, "validation") ){
    warning(effect, " options are not certified: ", certificate$dropped,
      " option(s) failed validation on the selection reference and were dropped.", call. = FALSE)
  }else if( !isTRUE(certificate$complete_up_to >= 0) ){
    warning(effect, " options are ranked but not certified: the step budget stopped the search before any member count was completed.", call. = FALSE)
  }else{
    warning(effect, " options are ranked and certified complete up to ",
      certificate$complete_up_to, " distinct source members. Larger options may exist",
      if( isTRUE(family$omitted_equal_smallest) ){
        ". Alternatives with the same member count as the smallest shown may be missing"
      }else{
        ""
      },
      " (step budget ", certificate$step_budget, ").", call. = FALSE)
  }
  return(invisible(NULL))
}


#' Build a DAG from one ranked level
#'
#' @inheritParams construct_dag
#' @param query Observed preparation and reuse key.
#' @param definitions Ranked levels and effect families.
#' @param index Selected level position.
#' @param sizes Known available sizes to save with the result.
#' @param smallest,top Whether this is the smallest or largest offered level.
#' @return A constructed DAG or dry-run plan from `.cd_level_dag()`.
#' @noRd
.cd_level <- function(graph,
                      latents,
                      keep_nodes,
                      query,
                      definitions,
                      index,
                      dry_run,
                      quiet,
                      sizes = NULL,
                      smallest = FALSE,
                      top = FALSE
                      ){
  level <- definitions$levels[[index]]
  progress <- .gd_progress(quiet, "Construction")
  .gd_notify(progress, "Preparing the reference graph")
  references <- .cd_size_prepared(graph, latents, keep_nodes, query)
  prepared <- references$prepared
  .gd_notify(progress, paste0("Reference ready (", length(names(prepared$selection)), " nodes)"))
  found <- lapply(definitions$effects, function(effect){
    set_size <- max(level$cap, effect$smallest, na.rm = TRUE)
    kept <- effect$sizes <= set_size
    family <- list(state = effect$state, members = as.character(unique(unlist(effect$members[kept], use.names = FALSE))),
      undecided = character(), sets = unlist(effect$sets[kept], recursive = FALSE))
    if( effect$state == "sets" && effect$certified && !top ){
      family$set_size <- set_size
    }
    return(family)
  })
  for( effect in names(definitions$effects) ){
    if( definitions$effects[[effect]]$state == "sets" && !definitions$effects[[effect]]$certified ){
      .cd_ranked_warning(effect, definitions$effects[[effect]])
    }
  }
  dag <- .cd_level_dag(graph, latents, keep_nodes, dry_run, quiet, progress, prepared, references$observed, found,
    references$key, list(method = "quantity-ranked-v1", representation = if( top ){
      "maximal"
    }else{
      "level"
    },
    arguments = list(latents = latents, keep_nodes = keep_nodes, search = "auto", steps = definitions$steps,
      size = NA_real_),
    level = list(cap = if( top ){
      Inf
    }else{
      level$cap
    }, rule = 0L, certified = level$certified, smallest = smallest, top = top), sizes = sizes))

  return(dag)
}


## Step budget ##
#' Read the variable-membership search budget
#'
#' @return `256L`, the step budget for each search for one quantity.
#' @noRd
.cd_maximal_steps <- function(){
  return(256L)
}


## One effect's search region ##
#' Prepare an effect region for variable-membership search
#'
#' @param prepared Fully observed construction preparation.
#' @param effect `"total"` or `"direct"`.
#' @return A list with region, indexed graph, moral adjacency, variable vocabulary,
#'   vertex labels, and masks for removable and ancestral vertices. Adjusted or
#'   selected nodes and removable nodes representing several members raise errors.
#' @noRd
.cd_maximal_setup <- function(prepared, effect){
  if( length(dagitty::adjustedNodes(prepared$selection)) || length(selected_nodes(prepared$selection)) ){
    stop("The DAG has adjusted or selected nodes. Remove those tags first.", call. = FALSE)
  }
  region <- if( effect == "total" ){
    .cd_ranked_total_region(prepared$selection)
  }else{
    .cd_ranked_direct_region(prepared$selection)
  }
  if( is.null(region) ){
    stop("The installed dagitty engine did not report direct effect constraints.", call. = FALSE)
  }
  if( length(region$must) ){
    stop("The DAG has adjusted or selected nodes. Remove those tags first.", call. = FALSE)
  }
  graph <- region$graph
  members <- .cd_ranked_labels(prepared, graph$nodes)
  members[graph$synthetic] <- list(character())
  if( any(lengths(members[region$removable]) != 1L) ){
    stop(errorCondition("A node of this query stands for several source members, so the largest DAG cannot be built. Use a smaller size.",
      class = "cd_several_members"))
  }
  source_order <- prepared$record$source_order
  if( !.query_order_applies(names(prepared$selection), source_order,
    unname(.gd_node_ids(prepared$selection))) ){
    source_order <- NULL
  }
  vocabulary <- .query_sort(unique(as.character(unlist(members[region$removable]))), source_order)
  label_of <- rep(NA_integer_, graph$n)
  label_of[region$removable] <- match(unlist(members[region$removable]), vocabulary)
  removable <- logical(graph$n)
  removable[region$removable] <- TRUE
  inside <- logical(graph$n)
  inside[region$ancestral] <- TRUE
  setup <- list(effect = effect, region = region, graph = graph,
    moral = .cd_ranked_moral(graph, region$ancestral), vocabulary = vocabulary,
    label_of = label_of, removable = removable, inside = inside)

  return(setup)
}


#' Traverse an effect region with whole variables open
#'
#' @param setup Region and labels from `.cd_maximal_setup()`.
#' @param start Starting vertex indices.
#' @param open Logical mask over the variable vocabulary. All occurrences of a
#'   variable are opened together.
#' @return A logical reached mask over graph vertices, restricted to the ancestral region.
#' @noRd
.cd_maximal_reach <- function(setup, start, open){
  passable <- setup$inside
  passable[setup$removable] <- open[setup$label_of[setup$removable]]
  seen <- logical(length(passable))
  seen[start] <- TRUE
  frontier <- start
  while( length(frontier) ){
    neighbours <- unlist(setup$moral[frontier], use.names = FALSE)
    neighbours <- unique(neighbours[!seen[neighbours] & passable[neighbours]])
    seen[neighbours] <- TRUE
    frontier <- neighbours
  }
  return(seen)
}


#' Check endpoint connectivity with whole variables open
#'
#' @inheritParams .cd_maximal_reach
#' @return `TRUE` if an outcome is reachable from the treatments, and `FALSE`
#'   otherwise. The ancestral moral graph is fixed, so opening more variables
#'   cannot disconnect the endpoints.
#' @noRd
.cd_maximal_connects <- function(setup, open){
  return(any(.cd_maximal_reach(setup, setup$region$x, open)[setup$region$y]))
}


## Exact occurrences: the graph between the endpoints ##
#' Eliminate nonadjustable vertices from the separator graph
#'
#' @param setup Effect region from `.cd_maximal_setup()`.
#' @return A list containing an adjacency matrix, neighbour lists and two added
#'   endpoint indices. Nonadjustable components join their adjustable neighbours.
#'   Each endpoint set becomes one vertex, preserving exact minimal separators.
#' @noRd
.cd_maximal_closure <- function(setup){
  count <- setup$graph$n
  adjacency <- matrix(FALSE, count, count)
  for( vertex in which(setup$inside) ){
    adjacency[vertex, setup$moral[[vertex]]] <- TRUE
  }
  adjacency <- adjacency | t(adjacency)
  treatment <- count + 1L
  outcome <- count + 2L
  closed <- matrix(FALSE, count + 2L, count + 2L)
  kept <- which(setup$removable)
  closed[kept, kept] <- adjacency[kept, kept]
  fixed <- setup$inside & !setup$removable
  grouped <- logical(count)
  for( vertex in which(fixed) ){
    if( grouped[vertex] ){
      next
    }
    group <- vertex
    frontier <- vertex
    grouped[vertex] <- TRUE
    while( length(frontier) ){
      neighbours <- which(colSums(adjacency[frontier, , drop = FALSE]) > 0L & fixed & !grouped)
      grouped[neighbours] <- TRUE
      group <- c(group, neighbours)
      frontier <- neighbours
    }
    ends <- which(colSums(adjacency[group, , drop = FALSE]) > 0L & setup$removable)
    if( any(group %in% setup$region$x) ){
      ends <- c(ends, treatment)
    }
    if( any(group %in% setup$region$y) ){
      ends <- c(ends, outcome)
    }
    closed[ends, ends] <- TRUE
  }
  diag(closed) <- FALSE
  closure <- list(adjacency = closed, neighbours = lapply(seq_len(count + 2L), function(vertex){
    return(which(closed[vertex, ]))
  }), treatment = treatment, outcome = outcome)

  return(closure)
}


#' Traverse the separator graph around blocked vertices
#'
#' @param closure Adjacency and endpoint indices from `.cd_maximal_closure()`.
#' @param start Starting vertex indices.
#' @param banned Logical mask of vertices that cannot be entered.
#' @return A logical reached mask including the starting vertices.
#' @noRd
.cd_maximal_closure_reach <- function(closure, start, banned){
  seen <- logical(length(closure$neighbours))
  seen[start] <- TRUE
  frontier <- start
  while( length(frontier) ){
    neighbours <- unlist(closure$neighbours[frontier], use.names = FALSE)
    neighbours <- unique(neighbours[!seen[neighbours] & !banned[neighbours]])
    seen[neighbours] <- TRUE
    frontier <- neighbours
  }
  return(seen)
}


# A vertex lies in an exact minimal set exactly when a chordless path between
# the endpoints passes through it. Such a path enters through one neighbour
# and leaves through another that is not adjacent to the first. The part before
# the vertex avoids the vertex's and the leaving neighbour's other neighbours,
# and the part after it avoids the vertex's and the entering neighbour's. This
# necessary condition excludes a vertex without any search.
#' Screen an occurrence for possible minimal-separator membership
#'
#' @param closure Adjacency and endpoint indices from `.cd_maximal_closure()`.
#' @param vertex Exact occurrence index to screen.
#' @return `FALSE` if neighbour reachability excludes a chordless endpoint path
#'   through the vertex. Otherwise returns `TRUE`, meaning membership still
#'   needs search.
#' @noRd
.cd_maximal_screen <- function(closure, vertex){
  adjacency <- closure$adjacency
  beside <- adjacency[vertex, ]
  around <- beside
  around[vertex] <- TRUE
  neighbours <- which(beside)
  outcome_side <- list()
  for( leaving in neighbours ){
    if( leaving == closure$treatment ){
      next
    }
    entering <- neighbours[!adjacency[leaving, neighbours] & neighbours != leaving & neighbours != closure$outcome]
    banned <- around | adjacency[leaving, ]
    banned[leaving] <- TRUE
    if( banned[closure$treatment] ){
      entering <- intersect(entering, closure$treatment)
    }else if( length(entering) ){
      reached <- .cd_maximal_closure_reach(closure, closure$treatment, banned)
      entering <- entering[colSums(adjacency[reached, entering, drop = FALSE]) > 0L]
    }
    for( first in entering ){
      if( leaving == closure$outcome ){
        return(TRUE)
      }
      key <- as.character(first)
      if( is.null(outcome_side[[key]]) ){
        banned_first <- around | adjacency[first, ]
        banned_first[first] <- TRUE
        outcome_side[[key]] <- if( banned_first[closure$outcome] ){
          logical(nrow(adjacency))
        }else{
          .cd_maximal_closure_reach(closure, closure$outcome, banned_first)
        }
      }
      if( any(outcome_side[[key]] & adjacency[leaving, ]) ){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}


## Two complete searches for one member ##
# Grows a chordless path from one endpoint, through an occurrence of the member,
# to the other endpoint. The other members on the path are opened with all their
# occurrences, and a path is abandoned as soon as they connect the endpoints by
# themselves. A finished path shows that blocking every member not on it
# separates the endpoints together with this member, and not without it.
#' Search for a membership witness along a chordless path
#'
#' @param setup Effect region and variable labels.
#' @param closure Adjacency and endpoint indices from `.cd_maximal_closure()`.
#' @param label Variable-vocabulary index whose membership is tested.
#' @param budget Maximum path-extension steps.
#' @param from_outcome Whether to grow the path from the outcome endpoint.
#' @return A list with `result` (`"member"`, `"excluded"` or `"undecided"`), steps
#'   used and an `open` variable mask for a witness. The mask is `NULL` when no
#'   witness is found. Budget exhaustion never establishes exclusion.
#' @noRd
.cd_maximal_route <- function(setup, closure, label, budget, from_outcome = FALSE){
  adjacency <- closure$adjacency
  start <- if( from_outcome ){
    closure$outcome
  }else{
    closure$treatment
  }
  goal <- if( from_outcome ){
    closure$treatment
  }else{
    closure$outcome
  }
  label_of <- c(setup$label_of, NA_integer_, NA_integer_)
  carriers <- which(label_of == label)
  distance <- function(from){
    steps_away <- rep(Inf, nrow(adjacency))
    steps_away[from] <- 0
    frontier <- from
    while( length(frontier) ){
      neighbours <- unique(unlist(closure$neighbours[frontier], use.names = FALSE))
      neighbours <- neighbours[!is.finite(steps_away[neighbours])]
      steps_away[neighbours] <- steps_away[frontier[[1L]]] + 1
      frontier <- neighbours
    }
    return(steps_away)
  }
  to_goal <- distance(goal)
  to_member <- distance(carriers)
  steps <- 0L
  exhausted <- FALSE
  extend <- function(head, shadow, open, passed){
    steps <<- steps + 1L
    if( steps > budget ){
      exhausted <<- TRUE
      return(NULL)
    }
    banned <- shadow
    banned[head] <- FALSE
    reachable <- .cd_maximal_closure_reach(closure, head, banned)
    if( !reachable[goal] || (!passed && !any(reachable[carriers])) ){
      return(NULL)
    }
    if( adjacency[head, goal] ){
      return(if( passed ){
        open
      }else{
        NULL
      })
    }
    candidates <- which(adjacency[head, ] & !shadow)
    candidates <- candidates[order(if( passed ){
      to_goal[candidates]
    }else{
      to_member[candidates]
    })]
    shadow <- shadow | adjacency[head, ]
    shadow[head] <- TRUE
    for( candidate in candidates ){
      carrier <- label_of[candidate] == label
      trial <- open
      if( !carrier && !open[label_of[candidate]] ){
        trial[label_of[candidate]] <- TRUE
        if( .cd_maximal_connects(setup, trial) ){
          next
        }
      }
      found <- extend(candidate, shadow, trial, passed || carrier)
      if( !is.null(found) || exhausted ){
        return(found)
      }
    }
    return(NULL)
  }
  first_shadow <- logical(nrow(adjacency))
  first_shadow[start] <- TRUE
  open <- extend(start, first_shadow, logical(length(setup$vocabulary)), FALSE)
  route <- list(result = if( !is.null(open) ){
    "member"
  }else if( exhausted ){
    "undecided"
  }else{
    "excluded"
  }, steps = min(steps, budget), open = open)

  return(route)
}


# Blocks members until the endpoints are separated once this member is blocked
# too, never letting the blocked members separate them without it. Each step
# blocks one more member that a remaining connecting path needs.
#' Search for a separator requiring one variable
#'
#' @param setup Effect region and variable labels.
#' @param label Variable-vocabulary index whose membership is tested.
#' @param budget Maximum blocking-search steps.
#' @param known Previously found minimal options, expressed as vocabulary indices.
#' @return A list with `result` (`"member"`, `"excluded"` or `"undecided"`), steps
#'   used and `blocked` indices for a witness. These indices are `NULL` when no
#'   witness is found. The caller adds the tested label to them.
#' @noRd
.cd_maximal_blocking <- function(setup, label, budget, known){
  count <- length(setup$vocabulary)
  weight <- as.integer(setup$removable)
  carriers <- which(setup$label_of == label)
  steps <- 0L
  exhausted <- FALSE
  visited <- new.env(hash = TRUE, parent = emptyenv())
  block <- function(blocked){
    key <- paste0("k", paste(blocked, collapse = ","))
    if( exists(key, envir = visited, inherits = FALSE) ){
      return(NULL)
    }
    assign(key, TRUE, envir = visited)
    for( option in known ){
      if( !label %in% option && all(option %in% blocked) ){
        return(NULL)
      }
    }
    steps <<- steps + 1L
    if( steps > budget ){
      exhausted <<- TRUE
      return(NULL)
    }
    open <- rep(TRUE, count)
    open[blocked] <- FALSE
    reached <- .cd_maximal_reach(setup, setup$region$x, open)
    if( !any(reached[setup$region$y]) || !any(reached[carriers]) ){
      return(NULL)
    }
    bypass <- .cd_ranked_path(setup$moral, setup$region$x, setup$region$y,
      which(setup$removable & setup$label_of %in% c(blocked, label)), weight)
    if( is.null(bypass) ){
      return(blocked)
    }
    for( move in unique(setup$label_of[bypass[setup$removable[bypass]]]) ){
      found <- block(sort(c(blocked, move)))
      if( !is.null(found) || exhausted ){
        return(found)
      }
    }
    return(NULL)
  }
  blocked <- block(integer())
  blocking <- list(result = if( !is.null(blocked) ){
    "member"
  }else if( exhausted ){
    "undecided"
  }else{
    "excluded"
  }, steps = min(steps, budget), blocked = blocked)

  return(blocking)
}


#' Reduce a variable separator while keeping its required witness
#'
#' @param setup Effect region and variable labels.
#' @param blocked Vocabulary indices forming a separator.
#' @param label Required witness index, exempt from removal.
#' @return Sorted indices after removing other variables whose release leaves
#'   the endpoints separated. All occurrences of a variable are released together.
#' @noRd
.cd_maximal_minimise <- function(setup, blocked, label){
  count <- length(setup$vocabulary)
  for( member in rev(setdiff(blocked, label)) ){
    open <- rep(TRUE, count)
    open[setdiff(blocked, member)] <- FALSE
    if( !.cd_maximal_connects(setup, open) ){
      blocked <- setdiff(blocked, member)
    }
  }
  option <- sort(blocked)

  return(option)
}


## Every member of one effect ##
# A `member` belongs to a set minimal under set inclusion, checked on the
# unchanged source. An `excluded` member failed the screen or a completed
# search found no such set. An `undecided` member exhausted every search
# budget, which gives no evidence either way. Each member has at most three
# searches, each limited by `steps`.
#' Classify variables used by minimal adjustment sets
#'
#' @param prepared Fully observed construction preparation.
#' @param effect `"total"` or `"direct"`.
#' @param steps Budget for each of up to three searches per variable.
#' @param progress Optional function reporting construction progress.
#' @param validate Exact-source validator accepting sets and effect.
#' @return An effect record with state, members, excluded and undecided variables,
#'   options expressed as members, validated exact sets, decision reasons and total
#'   steps. The states `"no_sets"` and `"empty_set"` remain distinct. Failed source
#'   validation raises an error.
#' @noRd
.cd_maximal_members <- function(prepared,
                                effect,
                                steps,
                                progress = NULL,
                                validate = .gd_validator(prepared$selection)
                                ){
  setup <- .cd_maximal_setup(prepared, effect)
  vocabulary <- setup$vocabulary
  count <- length(vocabulary)
  result <- list(effect = effect, state = "sets", members = character(), excluded = character(),
    undecided = character(), options = list(), sets = list(), decided_by = character(), steps = 0L)
  if( .cd_maximal_connects(setup, logical(count)) ){
    result$state <- "no_sets"
    result$excluded <- vocabulary
    return(result)
  }
  if( !.cd_maximal_connects(setup, rep(TRUE, count)) ){
    result$state <- "empty_set"
    result$excluded <- vocabulary
    result$options <- list(character())
    result$sets <- list(character())
    return(result)
  }
  closure <- .cd_maximal_closure(setup)
  carriers <- split(which(setup$removable), setup$label_of[setup$removable])
  possible <- logical(setup$graph$n)
  possible[setup$removable] <- vapply(which(setup$removable), .cd_maximal_screen, logical(1), closure = closure)
  status <- rep("pending", count)
  decided_by <- rep("", count)
  options <- list()
  for( label in seq_len(count) ){
    if( status[[label]] != "pending" ){
      next
    }
    .gd_notify(progress, paste0(effect, ": ", label, " of ", count, " quantities (",
      length(options), " sets found)"), force = FALSE)
    if( !any(possible[carriers[[as.character(label)]]]) ){
      status[[label]] <- "excluded"
      decided_by[[label]] <- "screen"
      next
    }
    attempts <- list(
      list(by = "path from treatment", run = function(){
        return(.cd_maximal_route(setup, closure, label, steps))
      }),
      list(by = "path from outcome", run = function(){
        return(.cd_maximal_route(setup, closure, label, steps, from_outcome = TRUE))
      }),
      list(by = "blocking search", run = function(){
        return(.cd_maximal_blocking(setup, label, steps, options))
      }))
    status[[label]] <- "undecided"
    for( attempt in attempts ){
      found <- attempt$run()
      result$steps <- result$steps + found$steps
      if( found$result == "undecided" ){
        next
      }
      decided_by[[label]] <- attempt$by
      if( found$result == "excluded" ){
        status[[label]] <- "excluded"
        break
      }
      blocked <- if( is.null(found$blocked) ){
        which(!found$open)
      }else{
        c(found$blocked, label)
      }
      option <- .cd_maximal_minimise(setup, blocked, label)
      if( any(status[option] == "excluded") ){
        stop("A set found by the search contains a quantity excluded earlier. No graph was built.", call. = FALSE)
      }
      options[[length(options) + 1L]] <- option
      covered <- option[status[option] != "member"]
      decided_by[setdiff(covered, label)] <- paste0("set found for ", vocabulary[[label]])
      status[option] <- "member"
      break
    }
  }
  sets <- lapply(options, function(option){
    allowed <- which(setup$removable & setup$label_of %in% option)
    separator <- .cd_ranked_separator(setup$moral, setup$region$x, setup$region$y, allowed)
    set <- setup$graph$nodes[sort(separator)]

    return(set)
  })
  if( length(sets) && !all(validate(sets, effect)) ){
    stop("A ", effect, " set found by the search failed validation on the source. No graph was built.", call. = FALSE)
  }
  result$members <- vocabulary[status == "member"]
  result$excluded <- vocabulary[status == "excluded"]
  result$undecided <- vocabulary[status == "undecided"]
  result$options <- lapply(options, function(option){
    return(vocabulary[option])
  })
  result$sets <- sets
  result$decided_by <- stats::setNames(decided_by, vocabulary)
  return(result)
}


## Building one size ##
#' Test exact adjustment-set minimality
#'
#' Uses the separator-neighbour criterion of TestMinAdj (van der Zander,
#' Liskiewicz and Textor, 2019) when its region and set guards hold.
#' @param graph Reference graph for the sets.
#' @param sets List of character vectors naming exact graph occurrences.
#' @param effect `"total"` or `"direct"`.
#' @param validate Validator accepting sets and effect for fallback cases.
#' @return A logical vector with one value per set: `TRUE` if it is a minimal
#'   valid set of exact occurrences, and `FALSE` otherwise. Unsupported fast
#'   cases use validity checks and single-occurrence deletions. This does not
#'   test minimality by source variable.
#' @noRd
.cd_maximal_minimal_sets <- function(graph, sets, effect, validate){
  minimal <- logical(length(sets))
  old <- rep(TRUE, length(sets))
  region <- if( !is.list(sets) || !inherits(graph, "dagitty") || !identical(dagitty::graphType(graph), "dag") ){
    NULL
  }else if( identical(effect, "total") ){
    tryCatch(.cd_ranked_total_region(graph), error = function(e) NULL)
  }else if( identical(effect, "direct") ){
    tryCatch(.cd_ranked_direct_region(graph), error = function(e) NULL)
  }else{
    NULL
  }
  if( !is.null(region) && !length(region$must) && length(region$x) && length(region$y) &&
    !length(intersect(region$x, region$y)) ){
    adjustable <- intersect(region$graph$nodes[region$removable], names(graph))
    moral <- .cd_ranked_moral(region$graph, region$ancestral)
    for( index in seq_along(sets) ){
      set <- sets[[index]]
      if( !is.character(set) || anyNA(set) || anyDuplicated(set) || !all(set %in% adjustable) ){
        next
      }
      members <- unname(region$graph$id[set])
      from_x <- logical(region$graph$n)
      from_x[.cd_ranked_reach(region$x, moral, members)] <- TRUE
      old[index] <- FALSE
      if( !any(from_x[region$y]) ){
        from_y <- logical(region$graph$n)
        from_y[.cd_ranked_reach(region$y, moral, members)] <- TRUE
        minimal[index] <- all(vapply(members, function(member){
          return(any(from_x[moral[[member]]]) && any(from_y[moral[[member]]]))
        }, logical(1)))
      }
    }
  }
  old <- which(old)
  if( length(old) ){
    checked <- old[validate(sets[old], effect)]
    minimal[checked] <- vapply(checked, function(index){
      set <- sets[[index]]
      set_minimal <- !length(set) || !any(validate(lapply(seq_along(set), function(drop){
        return(set[-drop])
      }), effect))

      return(set_minimal)
    }, logical(1))
  }
  return(minimal)
}


#' Require every supplied set to remain minimal
#'
#' @inheritParams .cd_maximal_minimal_sets
#' @return `TRUE` if every exact set is a minimal valid set, and `FALSE` otherwise.
#'   An empty list returns `TRUE`.
#' @noRd
.cd_maximal_minimal_on <- function(graph, sets, effect, validate = .gd_validator(graph)){
  if( !length(sets) ){
    return(TRUE)
  }
  minimal <- all(.cd_maximal_minimal_sets(graph, sets, effect, validate))

  return(minimal)
}


#' Check whether a display only renames temporal vertices
#'
#' @param temporal Temporal graph before summarisation.
#' @param display Proposed display graph.
#' @param mapping Data frame mapping temporal node to display name.
#' @return `TRUE` if the mapping is one-to-one and preserves graph type, typed
#'   edges and node roles, and `FALSE` otherwise. A display that passes has the
#'   same adjustment sets as the temporal graph, with nodes renamed.
#' @noRd
.cd_maximal_renamed <- function(temporal, display, mapping){
  if( anyDuplicated(mapping$display) || !setequal(mapping$node, names(temporal)) ){
    return(FALSE)
  }
  shown <- stats::setNames(as.character(mapping$display), mapping$node)
  renamed <- .gd_signature(temporal)
  edges <- dag_edges(temporal)
  edges$v <- unname(shown[edges$v])
  edges$w <- unname(shown[edges$w])
  symmetric <- edges$e %in% c("<->", "--")
  if( any(symmetric) ){
    left_endpoint <- pmin(edges$v[symmetric], edges$w[symmetric])
    right_endpoint <- pmax(edges$v[symmetric], edges$w[symmetric])
    edges$v[symmetric] <- left_endpoint
    edges$w[symmetric] <- right_endpoint
  }
  renamed$edges <- sort(unique(edge_key(edges)))
  for( role in c("nodes", "exposures", "outcomes", "latent", "adjusted", "selected") ){
    renamed[[role]] <- sort(unname(shown[renamed[[role]]]))
  }
  shown_signature <- .gd_signature(display)
  roles <- c("type", "nodes", "edges", "exposures", "outcomes", "latent", "adjusted", "selected")
  same <- identical(renamed[roles], shown_signature[roles])

  return(same)
}


## Sizes by rule ##
#' Find cumulative construction-rule nodes
#'
#' @param dag Fully observed analysis graph with endpoint and variable identities.
#' @return Six cumulative node vectors, each in graph order: immediate route steps,
#'   full routes, routes plus the first other descendants, routes plus all other
#'   descendants, those nodes plus outside parents, and those nodes plus upstream
#'   variables with an occasion having no parent from another variable.
#'   Endpoints are excluded.
#' @noRd
.cd_rule_nodes <- function(dag){
  edges <- dag_edges(dag)
  treatment <- as.character(dagitty::exposures(dag))
  outcome <- as.character(dagitty::outcomes(dag))
  ends <- c(treatment, outcome)
  after <- directed_reachable(edges, treatment, "from")
  route <- setdiff(intersect(after, directed_reachable(edges, outcome)), ends)
  dead <- setdiff(after, c(route, ends))
  first <- intersect(dead, children_edge_order(edges, c(treatment, route)))
  outside <- setdiff(parents_edge_order(edges, c(treatment, route, outcome)), c(route, ends))
  map <- .cd_variable_map(dag)
  ancestors <- setdiff(directed_reachable(edges, c(treatment, route, outcome)), ends)
  qualifying <- vapply(seq_len(nrow(map)), function(i){
    if( !map$node[[i]] %in% ancestors ){
      return(FALSE)
    }
    parents <- edges$v[edges$e == "->" & edges$w == map$node[[i]]]
    return(!any(map$variable[map$node %in% parents] != map$variable[[i]]))
  }, logical(1))
  exogenous <- setdiff(map$node[map$variable %in% map$variable[qualifying]], ends)
  steps <- intersect(route, union(children_edge_order(edges, treatment), parents_edge_order(edges, outcome)))
  rules <- list(steps, route, c(route, first), c(route, dead), c(route, dead, outside),
    c(route, dead, outside, exogenous))
  rule_nodes <- lapply(rules, function(nodes){
    return(intersect(names(dag), nodes))
  })

  return(rule_nodes)
}


#' Combine membership search with cumulative construction rules
#'
#' @param query Observed preparation and reuse key from `.cd_size_query()`.
#' @param quiet Whether to suppress routine construction progress.
#' @return A list containing effect `found` records, a source validator, distinct
#'   rule `rows` and a certification flag. Equal sets of added nodes share their
#'   lowest rule number. Undecided variables remain retained.
#' @noRd
.cd_rule_rows <- function(query, quiet){
  observed <- query$observed
  variable_map <- NULL
  count_variables <- function(nodes){
    if( is.null(variable_map) ){
      variable_map <<- .cd_variable_map(observed$analysis)
    }
    if( !all(nodes %in% variable_map$node) ){
      stop("A retained node has no recorded variable identity.", call. = FALSE)
    }
    return(length(unique(variable_map$variable[variable_map$node %in% nodes])))
  }
  budget <- .cd_maximal_steps()
  progress <- .gd_progress(quiet, "Construction")
  validate <- .gd_validator(observed$selection)
  found <- stats::setNames(lapply(c("direct", "total"), function(effect){
    return(.cd_kept(query$key, list("members", effect, budget), function(){
      .gd_notify(progress, paste0("Deciding the quantities used by ", effect, " adjustment sets"))
      family <- .cd_maximal_members(observed, effect, budget, progress, validate)

      return(family)
    }))
  }), c("direct", "total"))
  pins <- .cd_level_pins(observed, found)$nodes
  rule_nodes <- c(list(character()), .cd_rule_nodes(observed$analysis))
  extra <- lapply(rule_nodes, setdiff, pins)
  label <- match(extra, extra) - 1L
  kinds <- c("route", "route", "dead_ends", "dead_ends", "outside_causes", "exogenous_causes")
  first <- vapply(extra[[7L]], function(node){
    return(min(which(vapply(extra[-1L], `%in%`, logical(1), x = node))))
  }, integer(1))
  rows <- lapply(unique(label), function(rule){
    nodes <- extra[[rule + 1L]]
    kept <- table(factor(kinds[first[nodes]], levels = unique(kinds)))
    row <- list(rule = rule, extra = nodes,
      rule_nodes = rule_nodes[[if( rule == label[[7L]] ) 7L else rule + 1L]],
      variables = count_variables(c(pins, nodes)),
      top = rule == label[[7L]],
      kept = stats::setNames(as.integer(kept), names(kept))[kept > 0L])

    return(row)
  })
  rules <- list(found = found, validate = validate, rows = rows,
    certified = !length(unlist(lapply(found, `[[`, "undecided"))))

  return(rules)
}


#' Build a DAG from a variable-membership rule candidate
#'
#' @inheritParams construct_dag
#' @param sizes Known size information and selected time-view setting.
#' @param query Observed preparation and reuse key.
#' @param rules Membership results and rule rows from `.cd_rule_rows()`.
#' @param row Selected rule row, including extra nodes and largest-size flag.
#' @return A constructed DAG or dry-run plan retaining members and undecided
#'   variables, so the search budget cannot exclude a possible adjustment member.
#' @noRd
.cd_maximal <- function(graph,
                        latents,
                        keep_nodes,
                        dry_run,
                        quiet,
                        sizes,
                        query,
                        rules,
                        row
                        ){
  progress <- .gd_progress(quiet, "Construction")
  .gd_notify(progress, "Preparing the reference graph")
  references <- .cd_size_prepared(graph, latents, keep_nodes, query)
  prepared <- references$prepared
  .gd_notify(progress, paste0("Reference ready (", length(names(prepared$selection)), " nodes)"))
  found <- rules$found
  # The adjustment family has one size when each effect has no set, the empty
  # set, or a found set containing all its members with none undecided
  only <- all(vapply(found, function(family){
    return(family$state %in% c("no_sets", "empty_set") || (!length(family$undecided) &&
      any(vapply(family$options, setequal, logical(1), family$members))))
  }, logical(1)))
  dag <- .cd_level_dag(graph, latents, keep_nodes, dry_run, quiet, progress, prepared, references$observed, found,
    references$key, list(method = "maximal-quantities-v1", representation = "maximal",
    arguments = list(latents = latents, keep_nodes = keep_nodes, search = "exact", steps = .cd_maximal_steps(),
      size = if( row$top ){
        Inf
      }else{
        NA_real_
      }),
    level = list(cap = Inf, rule = row$rule, certified = rules$certified, smallest = only && !row$rule, top = row$top),
    sizes = sizes, extra = row$extra, kept = row$kept, rule_nodes = row$rule_nodes), rules$validate)

  return(dag)
}


#' Retain adjustment members and identification safeguards
#'
#' @param observed Fully observed construction preparation.
#' @param found Effect families with state, members, undecided variables and sets.
#' @return Pinned exact `nodes`, a total-effect identification record and
#'   `mediating` nodes preventing indirect paths from collapsing into direct edges.
#' @noRd
.cd_level_pins <- function(observed, found){
  none <- vapply(found, function(family){
    return(family$state %in% c("no_sets", "not_found"))
  }, logical(1))
  kept <- unique(unlist(lapply(found, function(family){
    return(c(family$members, family$undecided))
  }), use.names = FALSE))
  if( !all(kept %in% observed$map$member_id) ){
    stop("A retained quantity has no recorded source member in this query.", call. = FALSE)
  }
  quantity_nodes <- intersect(names(observed$analysis),
    observed$map$node[observed$map$member_id %in% kept])
  identification <- .cd_identification_retention(observed$analysis,
    list(status = if( none[["total"]] ){
      "no_sets"
    }else{
      "selected"
    }, reason = if( found$total$state == "not_found" ){
      "No total set was found within the step budget."
    }else{
      "No eligible set found."
    }))
  treatment <- as.character(dagitty::exposures(observed$analysis))
  outcome <- as.character(dagitty::outcomes(observed$analysis))
  # Without a direct set no quantity is retained for that effect. Keep the first
  # node of every indirect path, so collapsing cannot turn one into a direct edge.
  mediating <- if( none[["direct"]] ){
    setdiff(intersect(dagitty::children(observed$analysis, treatment),
      dagitty::ancestors(observed$analysis, outcome)), c(treatment, outcome))
  }else{
    character()
  }
  identification$added_nodes <- setdiff(identification$parents,
    union(union(quantity_nodes, observed$required), mediating))
  pinned <- list(nodes = union(union(union(quantity_nodes, observed$required), identification$parents), mediating),
    identification = identification, mediating = mediating)

  return(pinned)
}


#' Collapse and check a fully observed construction candidate
#'
#' @param observed Fully observed construction preparation.
#' @param found Effect families with state, members, undecided variables and sets.
#' @param method Construction method name recorded in check scope.
#' @param progress Optional function reporting construction progress.
#' @param validate Exact-source validator accepting sets and effect.
#' @param extra Exact nodes added by the chosen construction rule.
#' @return A list containing identification and mediating safeguards, checks,
#'   a checked summary, a merging flag and the number of displayed sets checked.
#'   Temporal listings are bounded. A summary that merges occasions requires
#'   every displayed set to pass validation on the source.
#' @noRd
.cd_level_core <- function(observed, found, method, progress, validate, extra){
  listed_limit <- 50L # Sets listed on a constructed graph for the checks
  merged_limit <- 1000L # Sets a time summary that merges occasions may list and still be checked in full
  effects <- c("direct", "total")
  pinned <- .cd_level_pins(observed, found)
  identification <- pinned$identification
  mediating <- pinned$mediating
  pins <- union(pinned$nodes, extra)
  treatment <- as.character(dagitty::exposures(observed$analysis))
  outcome <- as.character(dagitty::outcomes(observed$analysis))
  omitted <- setdiff(names(observed$analysis), pins)
  .gd_notify(progress, paste0("Collapsing ", length(omitted), " of ", length(names(observed$analysis)), " nodes"))
  temporal <- if( length(omitted) ){
    collapse_nodes(observed$analysis, omitted)
  }else{
    observed$analysis
  }
  temporal <- .cd_working_graph(temporal, observed$map[observed$map$node %in% names(temporal), , drop = FALSE],
    source_order = attr(observed$analysis, "cf_source_order", exact = TRUE))
  if( isTRUE(identification$applied) ){
    edges <- dag_edges(temporal)
    required <- union(treatment, union(outcome, identification$parents))
    if( !all(required %in% names(temporal)) || length(dagitty::latents(temporal)) ||
      any(edges$e == "<->" & (edges$v %in% treatment | edges$w %in% treatment)) ){
      stop("Retention of exact treatment parents failed its structural invariant.", call. = FALSE)
    }
    identification$status <- "preserved by retaining exact treatment parents"
  }
  .gd_notify(progress, "Checking the sets found and the first listed sets on the temporal graph")
  on_temporal <- .gd_validator(temporal)
  checks <- stats::setNames(lapply(effects, function(effect){
    family <- found[[effect]]
    listed <- .gd_adjustment_sets(temporal, effect, listed_limit)
    check <- list(effect = effect, state = family$state, sets_found = length(family$sets),
      sets_minimal = .cd_maximal_minimal_on(temporal, family$sets, effect, on_temporal),
      listed = length(listed$sets), listing_complete = listed$complete,
      listed_valid_on_source = !length(listed$sets) ||
      all(validate(listed$sets, effect)),
      no_set_kept = family$state != "no_sets" || !length(listed$sets))
    check$status <- if( family$state == "no_sets" ){
      if( check$no_set_kept ){
        "no eligible reference set"
      }else{
        "changed"
      }
    }else if( !check$sets_minimal || !check$listed_valid_on_source ){
      "changed"
    }else if( !check$listing_complete || family$state == "not_found" ){
      "not fully checked"
    }else{
      "preserved"
    }
    check$scope <- paste0("Asks whether the sets ", if( identical(method, "maximal-quantities-v1") ){
      "found while deciding the quantities"
    }else{
      "kept"
    }, " are minimal valid sets of the temporal graph read as fully observed, and whether the sets listed on it are valid on the source. ",
      if( listed$complete ){
        "Every set of that graph was listed."
      }else{
        paste0("Sets beyond the first ", listed_limit, " were not listed.")
      })
    return(check)
  }), effects)
  displayed <- 0L
  merges <- FALSE
  summary <- .gd_summary(temporal, keep = pins, preserve_occasions = isTRUE(identification$applied),
    source_check = function(display, mapping){
      # A summary that only renames nodes is the temporal graph itself. One that
      # merges occasions is not a projection, so every set it lists is checked.
      merges <<- !.cd_maximal_renamed(temporal, display, mapping)
      displayed <<- 0L
      on_display <- .gd_validator(display)
      for( effect in effects ){
        listed <- .gd_adjustment_sets(display, effect, if( merges ){
          merged_limit
        }else{
          listed_limit
        })
        if( merges && !listed$complete ){
          return(paste0("The time summary merges occasions and lists more than ", merged_limit, " ", effect,
            " adjustment sets, so not every one could be checked on the source. The temporal display was retained."))
        }
        sets <- listed$sets
        if( found[[effect]]$state == "no_sets" && length(sets) ){
          return(paste0("The source has no ", effect, " adjustment set, but the time summary lists one. The temporal display was retained."))
        }
        shown_sets <- lapply(found[[effect]]$sets, function(set){
          return(as.character(unique(mapping$display[match(set, mapping$node)])))
        })
        if( merges && (anyNA(unlist(shown_sets, use.names = FALSE)) ||
          !.cd_maximal_minimal_on(display, shown_sets, effect, on_display)) ){
          return(paste0("A ", effect, " adjustment set found on the source is not a minimal valid set of the time summary. The temporal display was retained."))
        }
        expanded <- lapply(sets, function(set){
          return(as.character(unique(unlist(lapply(set, function(node){
            occurrences <- mapping$node[mapping$display == node]
            represented <- if( length(occurrences) ){
              occurrences
            }else{
              node
            }

            return(represented)
          }), use.names = FALSE))))
        })
        displayed <<- displayed + length(sets)
        if( length(sets) && !all(validate(expanded, effect)) ){
          return(paste0("A temporally merged ", effect, " adjustment set is not valid on the source when its merged nodes are read as all represented occurrences. The temporal display was retained."))
        }
      }
      return(NULL)
    })
  core <- list(identification = identification, mediating = mediating, checks = checks, summary = summary,
    merges = merges, displayed = displayed)

  return(core)
}


#' Identify the reusable fully observed construction step
#'
#' @param observed Preparation supplying the required-node set.
#' @param found Effect families to retain and check.
#' @param method Construction method name.
#' @param extra Exact nodes added by construction rules.
#' @return A list key for use alongside the query's exact search identity.
#' @noRd
.cd_core_step <- function(observed, found, method, extra){
  return(list("core", found, observed$required, method, extra))
}


#' Finish a selected construction and attach its records
#'
#' @inheritParams construct_dag
#' @param progress Optional function reporting construction progress.
#' @param prepared Preparation carrying declared unmeasured occasions and limitations.
#' @param observed Fully observed preparation used for construction and search.
#' @param found Effect families with state, members, undecided variables and sets.
#' @param key Exact query key for reuse of the observed construction core.
#' @param record Candidate method, arguments, level, size, time view and rule information.
#' @param validate Exact-source validator accepting sets and effect.
#' @return A constructed DAG, or a dry-run plan when `dry_run = TRUE`. Construction
#'   treats all variables as observed. Unmeasured tags and availability checks
#'   are added afterwards.
#'   Provenance retains the source, temporal graph, decisions, checks and limitations.
#' @noRd
.cd_level_dag <- function(graph,
                          latents,
                          keep_nodes,
                          dry_run,
                          quiet,
                          progress,
                          prepared,
                          observed,
                          found,
                          key,
                          record,
                          validate = .gd_validator(observed$selection)
                          ){
  effects <- c("direct", "total")
  core <- .cd_kept(key, .cd_core_step(observed, found, record$method, record$extra), function(){
    return(.cd_level_core(observed, found, record$method, progress, validate, record$extra))
  })
  identification <- .cd_unmeasured_parents(core$identification, prepared$latent_nodes)
  identification$added_nodes <- setdiff(identification$added_nodes, record$rule_nodes)
  mediating <- core$mediating
  checks <- core$checks
  summary <- core$summary
  temporal_map <- .cd_variable_map(summary$temporal)
  repeated_variables <- any(table(temporal_map$variable) > 1L)
  plain <- is.null(.cd_record(graph)$node_map) && is.null(attr(graph, "cf_node_map", exact = TRUE))
  static_available <- isTRUE(summary$applied) || !repeated_variables
  static_reason <- summary$reason
  unroll <- isTRUE(record$sizes$unroll)
  if( !unroll && !static_available ){
    reason <- sub(" (The|Its) temporal display was retained[.]$", "", static_reason)
    stop(errorCondition(paste0("A checked static view is unavailable. ", reason),
      class = "cd_static_unavailable"))
  }
  if( unroll || plain ){
    temporal <- summary$temporal
    summary <- list(graph = temporal, temporal = temporal,
      node_map = data.frame(node = names(temporal), display = names(temporal), stringsAsFactors = FALSE),
      applied = FALSE, reason = if( unroll ) "Time-unrolled output requested." else "No time summary was applied.")
  }
  .gd_segment_tie_notice(summary)
  shown <- .cd_mark_unmeasured(summary, prepared$latent_nodes)
  summary <- shown$summary
  if( isTRUE(summary$applied) ){
    summary$merges_occasions <- core$merges
    summary$displayed_sets_checked <- core$displayed
  }
  final <- summary$graph
  # The declared unmeasured variables are marked now. Each effect's check lists
  # the declared occasions of the quantities it keeps, and says whether the
  # returned graph and the temporal graph still list a set that uses no
  # declared variable. Neither statement covers quantities these graphs do not
  # keep: the source may have a set without the declared variables that uses one.
  checks <- lapply(checks, function(check){
    family <- found[[check$effect]]
    family_nodes <- observed$map$node[observed$map$member_id %in% c(family$members, family$undecided)]
    check$unmeasured <- sort(intersect(family_nodes, shown$temporal), method = "radix")
    remains <- function(marked){
      return(if( family$state == "no_sets" ){
        "no adjustment set"
      }else if( !length(check$unmeasured) || length(.gd_adjustment_sets(marked, check$effect, 1L)$sets) ){
        "available now"
      }else{
        "needs an unmeasured variable"
      })
    }
    check$availability <- remains(final)
    check$temporal_availability <- if( isTRUE(summary$applied) ){
      remains(summary$temporal)
    }else{
      check$availability
    }
    return(check)
  })
  arguments <- record$arguments
  level <- record$level
  level$nodes <- length(names(final))
  level$variables <- .cd_variable_count(summary$temporal)
  level$static_available <- static_available
  level$static_reason <- static_reason
  arguments$unroll <- unroll
  if( is.na(arguments$size) ){
    arguments$size <- if( level$smallest && !identical(record$method, "maximal-quantities-v1") ){
      0
    }else{
      as.numeric(level$variables)
    }
  }
  quantities <- lapply(found, function(family){
    return(family[intersect(c("state", "members", "undecided", "set_size"), names(family))])
  })
  unmeasured <- list(semantics = "observed-as-if-v1", declared = latents, nodes = prepared$latent_nodes,
    display_nodes = shown$display, temporal_nodes = shown$temporal, partial = shown$partial)
  state <- .cd_representation_state(graph, summary, record$representation, record$method, arguments,
    quantities, checks, identification, unmeasured, prepared$limitations)
  state$construction$direct_path_nodes <- mediating
  if( length(record$kept) ){
    state$construction$kept_by_rule <- record$kept
  }
  state$construction$level <- level
  state$construction$metric <- "variables-v1"
  state$construction$sizes <- data.frame(variables = level$variables)
  result <- .gd_attach(final, state)
  .gd_notify(progress, paste0("Complete (", length(names(final)), " nodes)"))
  if( length(prepared$limitations) ){
    warning(paste(prepared$limitations, collapse = " "), call. = FALSE)
  }
  for( effect in effects ){
    family <- found[[effect]]
    if( length(family$undecided) ){
      single <- length(family$undecided) == 1L
      warning(length(family$undecided), " ", effect, if( single ){
        " quantity was not decided within the step budget ("
      }else{
        " quantities were not decided within the step budget ("
      }, arguments$steps, if( single ){
        ") and was retained. This is not evidence that a set uses it."
      }else{
        ") and were retained. This is not evidence that a set uses them."
      }, call. = FALSE)
    }
    check <- checks[[effect]]
    if( !check$no_set_kept ){
      warning("The source has no ", effect, " adjustment set, but the temporal graph lists one. The graph is not reliable for this effect.", call. = FALSE)
    }else if( identical(check$status, "changed") ){
      warning("The temporal graph changed the ", effect, " adjustment sets. Inspect dag_checks() and temporal_dag().", call. = FALSE)
    }
    if( identical(check$availability, "needs an unmeasured variable") ){
      warning("Every ", effect, " adjustment set of the returned DAG uses a declared unmeasured variable (kept for this effect: ",
        paste(check$unmeasured, collapse = ", "), ").",
        if( identical(check$temporal_availability, "available now") ){
          " temporal_dag() lists a set that uses none."
        }, call. = FALSE)
    }
  }
  absent <- effects[vapply(found, function(family){
    return(family$state == "no_sets")
  }, logical(1))]
  if( length(absent) ){
    warning("No eligible ", paste(absent, collapse = "/"),
      " adjustment set was found. Returning a diagnostic graph. This is not a valid empty set.", call. = FALSE)
  }
  unfound <- effects[vapply(found, function(family){
    return(family$state == "not_found")
  }, logical(1))]
  if( length(unfound) ){
    warning("No ", paste(unfound, collapse = "/"),
      " adjustment set was found within the step budget. This is not evidence that no set exists.", call. = FALSE)
  }
  if( isFALSE(quiet) ){
    cat("Causal query: ", paste(c(paste0(dagitty::exposures(graph), " (treatment)"),
      paste0(dagitty::outcomes(graph), " (outcome)")), collapse = ", "), "\n", sep = "")
    if( !dry_run ){
      .cd_print_level(state$construction)
    }
  }
  if( dry_run ){
    plan <- structure(list(graph = result, options = state$construction,
      status = vapply(found, `[[`, character(1), "state"), changes = .gd_delta(graph, result),
      checks = checks,
      arguments = c(list(graph = graph), arguments, list(dry_run = TRUE, quiet = quiet))), class = "cf_dag_plan")
    if( isFALSE(quiet) ){
      print(plan)
    }
    return(plan)
  }
  return(result)
}


## The record kept with each size ##
#' Assemble provenance for a constructed representation
#'
#' @param graph Original supplied query graph.
#' @param summary Checked summary containing temporal and display graphs.
#' @param representation Level or maximal representation label.
#' @param method Search or construction method name.
#' @param arguments Resolved construction arguments.
#' @param quantities Per-effect state, members, undecided variables and set-size cap.
#' @param checks Per-effect preservation and availability checks.
#' @param identification Exact-parent retention safeguard record.
#' @param unmeasured Declared, temporal, displayed and partial unmeasured membership.
#' @param limitations Query and preparation limitations to retain with the result.
#' @return A provenance list with original query, temporal graph, construction
#'   details and checks. The final graph signature is attached separately.
#' @noRd
.cd_representation_state <- function(graph,
                                     summary,
                                     representation,
                                     method,
                                     arguments,
                                     quantities,
                                     checks,
                                     identification,
                                     unmeasured,
                                     limitations
                                     ){
  return(list(original_query = graph, temporal = summary$temporal,
    construction = list(search_method = method, representation = representation, arguments = arguments,
    quantities = quantities,
    summary = summary[intersect(c("applied", "reason", "segment_alternatives", "segment_edge_ties",
    "merges_occasions", "displayed_sets_checked"), names(summary))],
    identification = identification, unmeasured = unmeasured,
    limitations = limitations, scope = if( identical(method, "maximal-quantities-v1") ){
      "Retains the endpoints, requested nodes, nodes specified by the size's rule and variables used by any minimal adjustment set for the direct or total effect. Minimality uses set inclusion of distinct source members. Variables left undecided by the step budget are retained. Nodes kept for identification are recorded under identification and direct_path_nodes. Checks cover the sets found during search and up to 50 sets from the temporal graph for each effect. dag_checks() reports whether that listing was complete."
    }else{
      "The endpoints, requested nodes and the quantities of every minimal adjustment set up to the set size recorded under level (cap), for each effect at least its smallest sets. Nodes kept for identification are recorded under identification and direct_path_nodes. The checks cover the sets kept and up to 50 sets of the temporal graph for each effect."
    }),
    checks = checks))
}


#' Print the retained adjustment families and rule nodes
#'
#' @param construction Saved construction record with level, quantities and rule counts.
#' @return `construction`, invisibly, after printing its size and effect summaries.
#' @noRd
.cd_print_level <- function(construction){
  largest <- identical(construction$representation, "maximal")
  cat(if( isTRUE(construction$level$top) ){
    "Largest DAG\n"
  }else{
    paste0("DAG with ", construction$level$variables, " nodes\n")
  })
  for( effect in names(construction$quantities) ){
    family <- construction$quantities[[effect]]
    quantities <- paste0(length(family$members), if( length(family$members) == 1L ){
      " quantity"
    }else{
      " quantities"
    })
    cat(effect, ": ", switch(family$state,
      sets = if( largest ){
        paste0(quantities, " used by minimal adjustment sets")
      }else if( is.null(family$set_size) ){
        paste0("the smallest adjustment sets found (", quantities, ")")
      }else{
        paste0("every minimal adjustment set of at most ", family$set_size, if( family$set_size == 1 ){
          " member ("
        }else{
          " members ("
        }, quantities, ")")
      },
      empty_set = "the empty set is the only minimal adjustment set",
      not_found = "no adjustment set found within the step budget",
      "no adjustment set exists"),
      if( length(family$undecided) ){
        paste0(", ", length(family$undecided), " undecided and retained")
      }, "\n", sep = "")
  }
  kept <- construction$kept_by_rule
  if( length(kept) ){
    words <- c(route = "route node", dead_ends = "dead end", outside_causes = "outside cause",
      exogenous_causes = "upstream node")
    cat("also kept: ", paste0(kept, " ", words[names(kept)], ifelse(kept == 1L, "", "s"), collapse = ", "), "\n",
      sep = "")
  }
  return(invisible(construction))
}
