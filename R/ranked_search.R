################################################################################
# ranked_search.R


# Ranked adjustment search by quantity (method "quantity-ranked-v1")
#
# Options are sets minimal under set inclusion, using resolved source members
# (time indices removed, recorded groups expanded), the unit the catalogue counts.
# A member set is feasible when some valid adjustment set uses only its allowed
# occurrences. Candidates lie in the ancestors of the endpoints, so the moral graph
# of that set is fixed and feasibility is one separation test (van der Zander,
# Liskiewicz and Textor 2019). Feasibility is monotone, so a smallest first
# hitting set search over paths connecting the endpoints returns exactly the
# options minimal under set inclusion, ranked by member count, with a certificate.
# Each option is realised as one set of exact temporal occurrences that is minimal
# under set inclusion and checked on the unchanged selection reference.


# Step budget per effect family: an explicit steps argument, else the option
.cd_ranked_steps <- function(steps = NULL) { # Validate the explicit or configured search budget
  origin <- if (is.null(steps)) {
    "Option causaliflower.ranked_steps"
  } else {
    "Argument steps"
  }
  if (is.null(steps)) {
    steps <- getOption("causaliflower.ranked_steps", 20000L)
  }
  if (!is.numeric(steps) || length(steps) != 1L || is.na(steps) || steps < 1 || steps > .Machine$integer.max ||
    steps != floor(steps)) {
    stop(origin, " must be one positive whole number no larger than ", .Machine$integer.max, ".", call. = FALSE)
  }
  return(as.integer(steps))
}


.cd_ranked_graph <- function(dag) { # Index graph vertices and edges for adjustment search
  edges <- as.data.frame(dagitty::edges(dag))
  edges <- if (all(c("v", "w", "e") %in% names(edges))) {
    edges[, c("v", "w", "e")]
  } else {
    data.frame(v = character(), w = character(), e = character())
  }
  edges[] <- lapply(edges, as.character)
  if (any(!edges$e %in% c("->", "<->"))) {
    stop("The ranked adjustment search needs directed or bidirected edges only.", call. = FALSE)
  }
  nodes <- names(dag)
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  if( .query_order_applies(nodes, source_order, unname(.gd_node_ids(dag))) ){
    node_keys <- stats::setNames(.query_order_key(nodes, source_order, unname(.gd_node_ids(dag))), nodes)
    nodes <- nodes[order(node_keys[nodes], method = "radix")]
    reversed <- edges$e == "<->" & node_keys[edges$v] > node_keys[edges$w]
    first_end <- edges$v[reversed]
    edges$v[reversed] <- edges$w[reversed]
    edges$w[reversed] <- first_end
    edges <- edges[order(node_keys[edges$v], node_keys[edges$w], edges$e, method = "radix"), , drop = FALSE]
  }
  if (length(dagitty::adjustedNodes(dag)) || length(selected_nodes(dag))) {
    stop("The ranked adjustment search needs a graph without adjusted or selected nodes. Use search = 'complete' for such graphs.", call. = FALSE)
  }
  latent <- as.character(unlist(dagitty::latents(dag)))
  directed <- edges[edges$e == "->", c("v", "w")]
  bidirected <- edges[edges$e == "<->", c("v", "w")]
  # A bidirected edge is an unobserved common parent of its two ends.
  if (nrow(bidirected)) {
    hidden <- paste0(".bidirected_", seq_len(nrow(bidirected)))
    nodes <- c(nodes, hidden)
    latent <- c(latent, hidden)
    directed <- rbind(directed, data.frame(v = c(hidden, hidden), w = c(bidirected$v, bidirected$w)))
  }
  node_index <- stats::setNames(seq_along(nodes), nodes)
  parent_index <- unname(node_index[directed$v])
  child_index <- unname(node_index[directed$w])
  node_count <- length(nodes)
  return(list(nodes = nodes, id = node_index, n = node_count, from = parent_index, to = child_index,
    parents = split(parent_index, factor(child_index, levels = seq_len(node_count))),
    children = split(child_index, factor(parent_index, levels = seq_len(node_count))),
    latent = unname(node_index[latent]),
    exposures = unname(node_index[as.character(unlist(dagitty::exposures(dag)))]),
    outcomes = unname(node_index[as.character(unlist(dagitty::outcomes(dag)))])))
}


.cd_ranked_reach <- function(start, adjacency, blocked = integer()) { # Find vertices reachable without crossing blocked vertices
  node_count <- length(adjacency)
  seen <- logical(node_count)
  seen[start] <- TRUE
  stopped <- logical(node_count)
  stopped[blocked] <- TRUE
  frontier <- start
  while (length(frontier)) {
    next_nodes <- unlist(adjacency[frontier], use.names = FALSE)
    next_nodes <- unique(next_nodes[!seen[next_nodes] & !stopped[next_nodes]])
    seen[next_nodes] <- TRUE
    frontier <- next_nodes
  }
  return(which(seen))
}


.cd_ranked_without <- function(graph, keep) { # Retain selected directed edges and rebuild adjacency lists
  graph$from <- graph$from[keep]
  graph$to <- graph$to[keep]
  graph$parents <- split(graph$from, factor(graph$to, levels = seq_len(graph$n)))
  graph$children <- split(graph$to, factor(graph$from, levels = seq_len(graph$n)))
  return(graph)
}


# Adjustment criterion: candidates are endpoint ancestors in the proper
# back door graph, excluding descendants of proper causal path nodes.
.cd_ranked_total_region <- function(dag) { # Prepare the ancestral region for total effect adjustment
  graph <- .cd_ranked_graph(dag)
  treatments <- graph$exposures
  outcomes <- graph$outcomes
  into_x_removed <- lapply(graph$children, function(children) { # Remove children that are treatment vertices
    return(children[!children %in% treatments])
  })
  out_of_x_removed <- lapply(graph$parents, function(parents) { # Remove parents that are treatment vertices
    return(parents[!parents %in% treatments])
  })
  causal <- intersect(setdiff(.cd_ranked_reach(treatments, into_x_removed), treatments),
    .cd_ranked_reach(outcomes, out_of_x_removed))
  forbidden <- .cd_ranked_reach(causal, graph$children)
  graph <- .cd_ranked_without(graph, !(graph$from %in% treatments & graph$to %in% causal))
  ancestral <- .cd_ranked_reach(c(treatments, outcomes), graph$parents)
  return(list(graph = graph, ancestral = ancestral, x = treatments, y = outcomes, must = integer(),
    removable = setdiff(ancestral, c(treatments, outcomes, forbidden, graph$latent))))
}


# Linear single door direct effect, exactly as .gd_validate_sets defines it:
# dagitty's own listMsasDirectEffect constraints and indirect graph.
.cd_ranked_direct_region <- function(dag) { # Read direct effect constraints from the installed dagitty engine
  script <- system.file("js", "dagitty-alg.js", package = "dagitty")
  if (!nzchar(script)) {
    stop("The installed dagitty JavaScript engine is unavailable.", call. = FALSE)
  }
  context <- V8::v8()
  context$source(script)
  context$assign("gd_graph_text", as.character(dag))
  context$eval(paste(
    "var gd_region = (function() {",
    "  var g = DAGitty.GraphParser.parseGuess(gd_graph_text);",
    "  var analyzer = DAGitty.GraphAnalyzer;",
    "  var constraints = null;",
    "  analyzer.listMinimalSeparators = function(graph, must, mustNot, max) {",
    "    constraints = {",
    "      must: must.map(function(v) { return v.id; }),",
    "      mustNot: mustNot.map(function(v) { return v.id; })",
    "    };",
    "    return [];",
    "  };",
    "  analyzer.listMsasDirectEffect(g, [], []);",
    "  return { must: constraints ? constraints.must : null,",
    "           mustNot: constraints ? constraints.mustNot : null,",
    "           indirect: DAGitty.GraphTransformer.indirectGraph(g).toString() };",
    "})();",
    sep = "\n"
    ))
  region <- context$get("gd_region", simplifyVector = FALSE)
  if (is.null(region$must)) {
    return(NULL)
  }
  indirect <- dagitty::dagitty(region$indirect)
  dagitty::latents(indirect) <- intersect(as.character(unlist(dagitty::latents(dag))), names(indirect))
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  if( .query_order_applies(names(dag), source_order, unname(.gd_node_ids(dag))) ){
    attr(indirect, "cf_source_order") <- source_order
  }
  attr(indirect, "input_labels") <- attr(dag, "input_labels", exact = TRUE)
  attr(indirect, "cf_node_map") <- attr(dag, "cf_node_map", exact = TRUE)
  graph <- .cd_ranked_graph(indirect)
  treatments <- graph$exposures
  outcomes <- graph$outcomes
  ancestral <- .cd_ranked_reach(c(treatments, outcomes), graph$parents)
  must_not <- unname(graph$id[as.character(unlist(region$mustNot))])
  return(list(graph = graph, ancestral = ancestral, x = treatments, y = outcomes,
    must = unname(graph$id[as.character(unlist(region$must))]),
    removable = setdiff(ancestral, c(treatments, outcomes, must_not, graph$latent))))
}


# Undirected moral graph of the ancestral subgraph
.cd_ranked_moral <- function(graph, ancestral) { # Connect the ancestral graph's neighbours and parents of the same child
  inside <- logical(graph$n)
  inside[ancestral] <- TRUE
  pairs <- lapply(ancestral, function(vertex) { # Connect this vertex to its parents and connect each parent pair
    parents <- graph$parents[[vertex]]
    parents <- parents[inside[parents]]
    if (!length(parents)) {
      return(NULL)
    }
    links <- cbind(vertex, parents)
    if (length(parents) > 1L) {
      links <- rbind(links, t(utils::combn(parents, 2L)))
    }
    return(links)
  })
  pairs <- do.call(rbind, pairs)
  if (is.null(pairs) || !nrow(pairs)) {
    return(vector("list", graph$n))
  }
  ends <- c(pairs[, 1L], pairs[, 2L])
  others <- c(pairs[, 2L], pairs[, 1L])
  return(lapply(split(others, factor(ends, levels = seq_len(graph$n))), unique))
}


# Level by level 0-1 search: a path connecting the endpoints through unblocked nodes
# with the fewest removable nodes. NULL when the endpoints are separated.
.cd_ranked_path <- function(moral, x, y, blocked, weight) { # Find an endpoint path crossing the fewest removable vertices
  node_count <- length(moral)
  parent <- integer(node_count)
  seen <- logical(node_count)
  seen[x] <- TRUE
  stopped <- logical(node_count)
  stopped[blocked] <- TRUE
  target <- logical(node_count)
  target[y] <- TRUE
  expand <- function(level, cost) { # Find the first unvisited neighbour reached at this cost
    to <- unlist(moral[level], use.names = FALSE)
    from <- rep.int(level, lengths(moral[level]))
    keep <- !seen[to] & !stopped[to] & weight[to] == cost
    to <- to[keep]
    from <- from[keep]
    first <- !duplicated(to)
    return(list(to = to[first], from = from[first]))
  }
  trace <- function(node) { # Trace the recorded parents back to a treatment vertex
    path <- node
    while (parent[path[[1L]]] != 0L) {
      path <- c(parent[path[[1L]]], path)
    }
    return(path)
  }
  level <- x
  repeat {
    frontier <- level
    while (length(frontier)) {
      hit <- frontier[target[frontier]]
      if (length(hit)) {
        return(trace(hit[[1L]]))
      }
      step <- expand(frontier, 0L)
      seen[step$to] <- TRUE
      parent[step$to] <- step$from
      level <- c(level, step$to)
      frontier <- step$to
    }
    step <- expand(level, 1L)
    if (!length(step$to)) {
      return(NULL)
    }
    seen[step$to] <- TRUE
    parent[step$to] <- step$from
    level <- step$to
  }
  return(invisible(NULL))
}


# Inclusion minimal separator inside `allowed` (closest to the treatment, then
# to the outcome).
.cd_ranked_separator <- function(moral, x, y, allowed) { # Find a minimal separator among the allowed vertices
  near <- .cd_ranked_reach(x, moral, blocked = allowed)
  first <- allowed[vapply(allowed, function(vertex) { # Keep allowed vertices adjacent to the treatment region
    return(any(moral[[vertex]] %in% near))
  }, logical(1))]
  near <- .cd_ranked_reach(y, moral, blocked = first)
  return(first[vapply(first, function(vertex) { # Keep the first separator's vertices adjacent to the outcome region
    return(any(moral[[vertex]] %in% near))
  }, logical(1))])
}


# Enumerate feasible member sets minimal under set inclusion, starting with the
# smallest and increasing the member count. `labels` lists member ids per node.
.cd_ranked_enumerate <- function(region, labels, step_budget) { # Enumerate minimal member sets within the explicit step budget
  graph <- region$graph
  moral <- .cd_ranked_moral(graph, region$ancestral)
  node_count <- graph$n
  node_labels <- vector("list", node_count)
  node_labels[region$removable] <- labels[region$removable]
  need <- lengths(node_labels)
  weight <- as.integer(need > 0L)
  label_nodes <- split(rep(seq_len(node_count), need), unlist(node_labels))
  base <- if (length(region$must)) {
    sort(unique(unlist(node_labels[region$must])))
  } else {
    integer()
  }
  blocked_by <- function(set) { # Find vertices whose source members are all covered by this set
    covered <- integer(node_count)
    for (label in as.character(set)) {
      member_vertices <- label_nodes[[label]]
      if (length(member_vertices)) {
        covered[member_vertices] <- covered[member_vertices] + 1L
      }
    }
    return(which(need > 0L & covered == need))
  }
  found <- list()
  steps <- 0L
  certificate <- "budget"
  reached <- -1L
  root <- .cd_ranked_path(moral, region$x, region$y, blocked_by(base), weight)
  if (is.null(root)) {
    found <- list(base)
    certificate <- "complete"
    reached <- length(base)
  } else if (!any(need[root] > 0L)) {
    certificate <- "complete"
    reached <- length(unique(unlist(node_labels)))
  } else {
    for (bound in seq.int(length(base), max(length(base), length(label_nodes)))) {
      truncated <- FALSE
      exhausted <- FALSE
      visited <- new.env(hash = TRUE, parent = emptyenv())
      expand <- function(set) { # Search this member set and its bounded extensions
        if (exhausted) {
          return(invisible())
        }
        key <- paste0("k", paste(set, collapse = ","))
        if (exists(key, envir = visited, inherits = FALSE)) {
          return(invisible())
        }
        assign(key, TRUE, envir = visited)
        if (any(vapply(found, function(existing_option) { # Check whether a recorded option is contained in this set
          return(all(existing_option %in% set))
        }, logical(1)))) {
          return(invisible())
        }
        steps <<- steps + 1L
        if (steps > step_budget) {
          exhausted <<- TRUE
          return(invisible())
        }
        path <- .cd_ranked_path(moral, region$x, region$y, blocked_by(set), weight)
        if (is.null(path)) {
          found[[length(found) + 1L]] <<- set
          return(invisible())
        }
        moves <- unique(lapply(path[need[path] > 0L], function(vertex) { # Collect the additional members needed to block this vertex
          return(setdiff(node_labels[[vertex]], set))
        }))
        for (move in moves) {
          candidate <- sort(c(set, move))
          if (length(candidate) > bound) {
            truncated <<- TRUE
            next
          }
          expand(candidate)
        }
        return(invisible(NULL))
      }
      expand(base)
      if (exhausted) {
        break
      }
      reached <- bound
      if (!truncated) {
        certificate <- "complete"
        break
      }
    }
  }
  realised <- lapply(found, function(set) { # Realise this member option as an adjustment set of exact temporal occurrences
    allowed <- region$removable[vapply(labels[region$removable], function(vertex_members) { # Check whether this vertex's members are allowed by the option
      return(all(vertex_members %in% set))
    }, logical(1))]
    separator <- union(.cd_ranked_separator(moral, region$x, region$y, allowed), region$must)
    return(graph$nodes[sort(separator)])
  })
  return(list(options = found, realised = realised, certificate = certificate,
    complete_up_to = reached, steps = min(steps, step_budget)))
}


.cd_ranked_labels <- function(prepared, nodes) { # Resolve each graph vertex to its distinct source members
  inputs <- attr(prepared$selection, "input_labels")
  base <- .strip_slice(nodes)
  mapped <- base %in% names(inputs)
  base[mapped] <- unname(inputs[base[mapped]])
  return(lapply(base, function(concept) { # Use recorded members when this quantity has a mapping
    members <- unique(prepared$map$member_id[prepared$map$quantity %in% concept])
    return(if (length(members)) {
      members
    } else {
      concept
    })
  }))
}


.cd_ranked_search <- function(prepared, effect, progress = NULL, steps = NULL) { # Rank and validate adjustment options by distinct source members
  .gd_notify(progress, paste0("Ranking ", effect, " adjustment options by distinct source members"))
  region <- if (effect == "total") {
    .cd_ranked_total_region(prepared$selection)
  } else {
    .cd_ranked_direct_region(prepared$selection)
  }
  if (is.null(region)) {
    stop("The installed dagitty engine did not report direct effect constraints.", call. = FALSE)
  }
  members <- .cd_ranked_labels(prepared, region$graph$nodes)
  members[startsWith(region$graph$nodes, ".bidirected_")] <- list(character())
  source_order <- prepared$record$source_order
  if( !.query_order_applies(names(prepared$selection), source_order,
                          unname(.gd_node_ids(prepared$selection))) ){
    source_order <- NULL
  }
  vocabulary <- .query_sort(unique(as.character(unlist(members[region$removable]))), source_order)
  labels <- lapply(members, function(vertex_members) { # Index each vertex's members in the common vocabulary
    return(match(vertex_members, vocabulary))
  })
  budget <- .cd_ranked_steps(steps)
  ranked <- .cd_ranked_enumerate(region, labels, budget)
  sets <- ranked$realised
  valid <- if (length(sets)) {
    .gd_validate_sets(prepared$selection, sets, effect)
  } else {
    logical()
  }
  complete <- identical(ranked$certificate, "complete")
  dropped <- 0L
  if (!all(valid)) {
    # A realised option that fails the source check is removed, and nothing is
    # certified: the dropped option may have been the smallest.
    dropped <- sum(!valid)
    warning(effect, " ranked search dropped ", dropped,
      " option(s) that failed validation on the selection reference. The ranking is not certified.", call. = FALSE)
    sets <- sets[valid]
    ranked$options <- ranked$options[valid]
    complete <- FALSE
    ranked$certificate <- "validation"
    ranked$complete_up_to <- -1L
  }
  sizes <- lengths(ranked$options)
  order <- order(sizes, vapply(sets, function(adjustment_set) { # Build the deterministic key for breaking ties for this realised set
    ids <- unname(.gd_node_ids(prepared$selection)[adjustment_set])
    return(paste(.query_order_key(adjustment_set, source_order, ids), collapse = "\r"))
  }, character(1)), method = "radix")
  sets <- unname(sets[order])
  .gd_notify(progress, paste0(effect, " ranking ", if (complete) {
    "complete"
  } else {
    paste0("certified up to ", ranked$complete_up_to, " members")
  }, " (", length(sets), " options). Grouping source quantities"))
  table <- .cd_candidates(sets, prepared, progress)
  return(list(sets = sets, source_sets = sets, complete = complete, exhaustive_requested = FALSE,
    probe = NULL, limit = budget, condition_free = TRUE, table = table,
    probe_quantity_count = NA_integer_,
    omitted_equal_smallest = !complete && length(sizes) > 0L && ranked$complete_up_to < min(sizes),
    ranked = list(method = "quantity-ranked-v1", certificate = if (complete) {
      "complete"
    } else {
      ranked$certificate
    },
    complete_up_to = ranked$complete_up_to, steps = ranked$steps, step_budget = budget,
    dropped = dropped),
    scope = if (complete) {
      "Every option minimal under set inclusion by distinct source members."
    } else if (identical(ranked$certificate, "validation")) {
      "Options that failed validation on the selection reference were dropped. The ranking is not certified."
    } else if (!isTRUE(ranked$complete_up_to >= 0)) {
      "The step budget stopped the search before any member count was completed. Listed options are valid but the ranking is not certified."
    } else {
      paste0("Every option minimal under set inclusion with at most ", ranked$complete_up_to,
        " distinct source members. Larger options may exist.")
    },
    optimality = if (!length(sets)) {
      (if (complete) {
        "no eligible set"
      } else if (identical(ranked$certificate, "validation")) {
        "every option failed validation on the selection reference"
      } else {
        "none found within the step budget"
      })
    } else if (!length(sets[[1L]])) {
      "valid empty set is the unique set minimal under set inclusion"
    } else if (identical(ranked$certificate, "validation")) {
      "smallest member count not certified"
    } else if (!isTRUE(ranked$complete_up_to >= min(sizes))) {
      "smallest member count certified. Alternatives of that size may be missing"
    } else {
      "smallest member count certified"
    }))
}


.cd_reuse_ranked <- function(previous, prepared, identity, progress = NULL) { # Reuse ranked options only when query identity and validation agree
  if (is.null(previous) || is.null(identity) || !isTRUE(getOption("causaliflower.accelerate", TRUE))) {
    return(NULL)
  }
  result <- tryCatch({
    state <- .gd_set_source(previous)
    if (!identical(state$construction$search_method, "quantity-ranked-v1") ||
      !identical(state$construction$search_identity, identity)) {
      return(NULL)
    }
    arguments <- state$construction$arguments
    # The ranked search runs as fully observed: compare the observed preparation
    # of the saved query with the saved analysis without its declared marks.
    original <- .cd_prepare(state$original_query, NULL, arguments$keep_nodes)
    original_identity <- .cd_search_identity(original)
    if (is.null(original_identity)) {
      return(NULL)
    }
    original_identity$method <- identity$method
    original_identity$ranked_adapter <- identity$ranked_adapter
    original_identity$ranked_steps <- identity$ranked_steps
    unmarked <- .cd_unmarked(state$analysis)
    if (!identical(original_identity, identity) ||
      !identical(.gd_signature(original$analysis), .gd_signature(unmarked)) ||
      !identical(attr(original$analysis, "cf_node_map"), attr(state$analysis, "cf_node_map")) ||
      !identical(.gd_signature(prepared$analysis), .gd_signature(unmarked)) ||
      !identical(attr(prepared$analysis, "cf_node_map"), attr(state$analysis, "cf_node_map"))) {
      return(NULL)
    }
    search <- state$build$search
    for (effect in c("direct", "total")) {
      item <- search[[effect]]
      if (!is.list(item$sets) || !identical(item$source_sets, item$sets) ||
        !identical(item$ranked$method, "quantity-ranked-v1") ||
        (length(item$sets) && !all(.gd_validate_sets(prepared$selection, item$sets, effect)))) {
        return(NULL)
      }
      item$table <- .cd_candidates(item$sets, prepared)
      search[[effect]] <- item
    }
    search
  }, error = function(e) { # Decline reuse when checking the previous result raises a condition
    return(NULL)
  }, warning = function(w) { # Decline reuse when checking the previous result raises a condition
    return(NULL)
  }, message = function(m) { # Decline reuse when checking the previous result raises a condition
    return(NULL)
  })
  if (!is.null(result)) {
    .gd_notify(progress, "Reusing exact ranked source options. Rebuilding construction and checks")
  }
  return(result)
}
