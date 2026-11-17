# ==============================================================================
# ranked_search.R


# Ranked adjustment search by quantity (method "quantity-ranked-v1")
# ==============================================================================


# Options are minimal under set inclusion of recorded source members, with
# temporal occurrences counted together and recorded groups expanded. A member
# set is feasible if some valid adjustment set uses only its allowed
# occurrences. Candidates lie in the endpoints' ancestors, whose fixed moral
# graph makes feasibility one separation test (van der Zander, Liskiewicz and
# Textor 2019). The bounded hitting set search visits member sets in
# increasing size and records how far the ranking is complete. Each option is
# realised as one minimal set of exact temporal occurrences and checked on the
# unchanged selection reference.


#' Read the ranked-search budget
#'
#' @return The positive integer step budget for each effect. It is read from
#'   `causaliflower.ranked_steps`, or is 20000 when that option is unset.
#'   An invalid option raises an error.
#' @noRd
.cd_ranked_steps <- function(){
  steps <- getOption("causaliflower.ranked_steps", 20000L)
  if( !is.numeric(steps) || length(steps) != 1L || is.na(steps) || steps < 1 || steps > .Machine$integer.max ||
    steps != floor(steps) ){
    stop("Option causaliflower.ranked_steps must be one positive whole number no larger than ", .Machine$integer.max,
      ".", call. = FALSE)
  }
  budget <- as.integer(steps)

  return(budget)
}


#' Index a directed or bidirected graph for separator search
#'
#' @param dag Graph with treatment and outcome tags and no adjusted or selected nodes.
#' @return A list of node names, indices, parent and child lists, directed edges
#'   and status indices. Bidirected edges become explicit synthetic latent parents.
#'   `synthetic` distinguishes them from real latent nodes without name guessing.
#' @noRd
.cd_ranked_graph <- function(dag){
  edges <- as.data.frame(dagitty::edges(dag))
  edges <- if( all(c("v", "w", "e") %in% names(edges)) ){
    edges[, c("v", "w", "e")]
  }else{
    data.frame(v = character(), w = character(), e = character())
  }
  edges[] <- lapply(edges, as.character)
  if( any(!edges$e %in% c("->", "<->")) ){
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
  if( length(dagitty::adjustedNodes(dag)) || length(selected_nodes(dag)) ){
    stop("The ranked adjustment search needs a graph without adjusted or selected nodes.", call. = FALSE)
  }
  latent <- as.character(unlist(dagitty::latents(dag)))
  directed <- edges[edges$e == "->", c("v", "w")]
  bidirected <- edges[edges$e == "<->", c("v", "w")]
  hidden <- character()
  if( nrow(bidirected) ){
    canonical <- dagitty::canonicalize(dag)
    hidden <- as.character(unlist(canonical$L))
    inserted <- as.data.frame(dagitty::edges(canonical$g))
    inserted <- inserted[inserted$v %in% hidden, c("v", "w")]
    pair_key <- function(ends){
      return(paste(sort(match(ends, nodes)), collapse = ":"))
    }
    hidden_keys <- vapply(split(inserted$w, factor(inserted$v, levels = hidden)), pair_key, character(1))
    edge_keys <- vapply(seq_len(nrow(bidirected)), function(i){
      return(pair_key(c(bidirected$v[i], bidirected$w[i])))
    }, character(1))
    hidden <- hidden[match(edge_keys, hidden_keys)]
    # Preserve source traversal order after Dagitty adds latent parents
    inserted_keys <- paste(match(inserted$v, hidden), match(inserted$w, nodes), sep = ":")
    wanted_keys <- paste(rep(seq_along(hidden), 2L), match(c(bidirected$v, bidirected$w), nodes), sep = ":")
    inserted <- inserted[match(wanted_keys, inserted_keys), , drop = FALSE]
    nodes <- c(nodes, hidden)
    latent <- c(latent, hidden)
    directed <- rbind(directed, inserted)
  }
  node_index <- stats::setNames(seq_along(nodes), nodes)
  parent_index <- unname(node_index[directed$v])
  child_index <- unname(node_index[directed$w])
  node_count <- length(nodes)
  graph <- list(nodes = nodes, id = node_index, n = node_count, from = parent_index, to = child_index,
    parents = split(parent_index, factor(child_index, levels = seq_len(node_count))),
    children = split(child_index, factor(parent_index, levels = seq_len(node_count))),
    latent = unname(node_index[latent]), synthetic = unname(node_index[hidden]),
    exposures = unname(node_index[as.character(unlist(dagitty::exposures(dag)))]),
    outcomes = unname(node_index[as.character(unlist(dagitty::outcomes(dag)))]))

  return(graph)
}


#' Traverse indexed adjacency lists
#'
#' @param start Starting vertex indices.
#' @param adjacency One neighbour-index vector per vertex.
#' @param blocked Vertex indices that cannot be entered.
#' @return Reached indices in vertex order, including the starting vertices.
#' @noRd
.cd_ranked_reach <- function(start, adjacency, blocked = integer()){
  node_count <- length(adjacency)
  seen <- logical(node_count)
  seen[start] <- TRUE
  stopped <- logical(node_count)
  stopped[blocked] <- TRUE
  frontier <- start
  while( length(frontier) ){
    next_nodes <- unlist(adjacency[frontier], use.names = FALSE)
    next_nodes <- unique(next_nodes[!seen[next_nodes] & !stopped[next_nodes]])
    seen[next_nodes] <- TRUE
    frontier <- next_nodes
  }
  reached <- which(seen)

  return(reached)
}


#' Retain selected directed edges in an indexed graph
#'
#' @param graph Indexed graph from `.cd_ranked_graph()`.
#' @param keep Logical mask or indices selecting rows of its directed edges.
#' @return The graph with selected edges and rebuilt parent and child lists.
#' @noRd
.cd_ranked_without <- function(graph, keep){
  graph$from <- graph$from[keep]
  graph$to <- graph$to[keep]
  graph$parents <- split(graph$from, factor(graph$to, levels = seq_len(graph$n)))
  graph$children <- split(graph$to, factor(graph$from, levels = seq_len(graph$n)))
  return(graph)
}


#' Derive the total-effect adjustment region
#'
#' @param dag Graph accepted by `.cd_ranked_graph()`.
#' @return A list containing the indexed proper back-door graph, ancestral vertices, endpoint
#'   indices, required vertices and removable candidates. Descendants of proper
#'   causal-path nodes, endpoints and latent vertices cannot be removed.
#' @noRd
.cd_ranked_total_region <- function(dag){
  graph <- .cd_ranked_graph(dag)
  treatments <- graph$exposures
  outcomes <- graph$outcomes
  into_x_removed <- lapply(graph$children, function(children){
    return(children[!children %in% treatments])
  })
  out_of_x_removed <- lapply(graph$parents, function(parents){
    return(parents[!parents %in% treatments])
  })
  causal <- intersect(setdiff(.cd_ranked_reach(treatments, into_x_removed), treatments),
    .cd_ranked_reach(outcomes, out_of_x_removed))
  forbidden <- .cd_ranked_reach(causal, graph$children)
  graph <- .cd_ranked_without(graph, !(graph$from %in% treatments & graph$to %in% causal))
  ancestral <- .cd_ranked_reach(c(treatments, outcomes), graph$parents)
  region <- list(graph = graph, ancestral = ancestral, x = treatments, y = outcomes, must = integer(),
    removable = setdiff(ancestral, c(treatments, outcomes, forbidden, graph$latent)))

  return(region)
}


#' Derive the installed engine's direct-effect adjustment region
#'
#' @param dag Graph supplying endpoints, latent tags and source identities.
#' @return A list containing the indexed indirect graph, ancestral vertices,
#'   endpoints, required vertices and removable candidates from dagitty's linear
#'   single-door routine. Returns `NULL` if the engine reports no constraints.
#'   An unavailable engine script raises an error.
#' @noRd
.cd_ranked_direct_region <- function(dag){
  script <- system.file("js", "dagitty-alg.js", package = "dagitty")
  if( !nzchar(script) ){
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
  if( is.null(region$must) ){
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
  search_region <- list(graph = graph, ancestral = ancestral, x = treatments, y = outcomes,
    must = unname(graph$id[as.character(unlist(region$must))]),
    removable = setdiff(ancestral, c(treatments, outcomes, must_not, graph$latent)))

  return(search_region)
}


#' Moralise an indexed ancestral subgraph
#'
#' @param graph Indexed directed graph.
#' @param ancestral Indices retained in the ancestral subgraph.
#' @return One undirected neighbour list per graph vertex, joining each child's
#'   parents. Vertices outside the region have empty lists.
#' @noRd
.cd_ranked_moral <- function(graph, ancestral){
  inside <- logical(graph$n)
  inside[ancestral] <- TRUE
  pairs <- lapply(ancestral, function(vertex){
    parents <- graph$parents[[vertex]]
    parents <- parents[inside[parents]]
    if( !length(parents) ){
      return(NULL)
    }
    links <- cbind(vertex, parents)
    if( length(parents) > 1L ){
      links <- rbind(links, t(utils::combn(parents, 2L)))
    }
    return(links)
  })
  pairs <- do.call(rbind, pairs)
  if( is.null(pairs) || !nrow(pairs) ){
    return(vector("list", graph$n))
  }
  ends <- c(pairs[, 1L], pairs[, 2L])
  others <- c(pairs[, 2L], pairs[, 1L])
  moral <- lapply(split(others, factor(ends, levels = seq_len(graph$n))), unique)

  return(moral)
}


#' Find a connecting path with the fewest removable vertices
#'
#' @param moral Undirected neighbour lists for the search region.
#' @param x,y Treatment and outcome vertex indices.
#' @param blocked Vertex indices that cannot be entered.
#' @param weight A cost of 0 or 1 per vertex, with 1 marking removable candidates.
#' @return A path of vertex indices from a treatment to an outcome, or `NULL`
#'   when the blocked vertices separate the endpoints.
#' @noRd
.cd_ranked_path <- function(moral, x, y, blocked, weight){
  node_count <- length(moral)
  parent <- integer(node_count)
  seen <- logical(node_count)
  seen[x] <- TRUE
  stopped <- logical(node_count)
  stopped[blocked] <- TRUE
  target <- logical(node_count)
  target[y] <- TRUE
  expand <- function(level, cost){
    to <- unlist(moral[level], use.names = FALSE)
    from <- rep.int(level, lengths(moral[level]))
    keep <- !seen[to] & !stopped[to] & weight[to] == cost
    to <- to[keep]
    from <- from[keep]
    first <- !duplicated(to)
    step <- list(to = to[first], from = from[first])

    return(step)
  }
  trace <- function(node){
    path <- node
    while( parent[path[[1L]]] != 0L ){
      path <- c(parent[path[[1L]]], path)
    }
    return(path)
  }
  level <- x
  repeat{
    frontier <- level
    while( length(frontier) ){
      hit <- frontier[target[frontier]]
      if( length(hit) ){
        return(trace(hit[[1L]]))
      }
      step <- expand(frontier, 0L)
      seen[step$to] <- TRUE
      parent[step$to] <- step$from
      level <- c(level, step$to)
      frontier <- step$to
    }
    step <- expand(level, 1L)
    if( !length(step$to) ){
      return(NULL)
    }
    seen[step$to] <- TRUE
    parent[step$to] <- step$from
    level <- step$to
  }
  return(invisible(NULL))
}


#' Reduce an allowed separator to inclusion minimality
#'
#' @param moral Undirected neighbour lists for the search region.
#' @param x,y Treatment and outcome vertex indices.
#' @param allowed Vertex indices known to separate the endpoints when blocked.
#' @return An inclusion-minimal separator within `allowed`, taking the treatment
#'   boundary first and then the outcome boundary.
#' @noRd
.cd_ranked_separator <- function(moral, x, y, allowed){
  near <- .cd_ranked_reach(x, moral, blocked = allowed)
  first <- allowed[vapply(allowed, function(vertex){
    return(any(moral[[vertex]] %in% near))
  }, logical(1))]
  near <- .cd_ranked_reach(y, moral, blocked = first)
  separator <- first[vapply(first, function(vertex){
    return(any(moral[[vertex]] %in% near))
  }, logical(1))]

  return(separator)
}


#' Enumerate minimal member sets in increasing member count
#'
#' @param region Effect region with graph, endpoints, removable and required indices.
#' @param labels List of member-vocabulary indices represented by each graph node.
#' @param step_budget Maximum search expansions across all size bounds.
#' @return A list containing member options, separators of exact nodes, the
#'   certificate, `complete_up_to` and steps used. `complete_up_to` is the largest
#'   member-count bound completed, or `-1` if none was completed. Budget exhaustion
#'   retains partial results without claiming completeness beyond that bound.
#' @noRd
.cd_ranked_enumerate <- function(region, labels, step_budget){
  graph <- region$graph
  moral <- .cd_ranked_moral(graph, region$ancestral)
  node_count <- graph$n
  node_labels <- vector("list", node_count)
  node_labels[region$removable] <- labels[region$removable]
  need <- lengths(node_labels)
  weight <- as.integer(need > 0L)
  label_nodes <- split(rep(seq_len(node_count), need), unlist(node_labels))
  base <- if( length(region$must) ){
    sort(unique(unlist(node_labels[region$must])))
  }else{
    integer()
  }
  plain_members <- is.integer(base) && is.null(attributes(base)) &&
    all(vapply(node_labels[region$removable], function(members){
      return(is.integer(members) && is.null(attributes(members)) &&
        !anyNA(members) && !anyDuplicated(members) && all(members > 0L))
    }, logical(1)))
  blocked_by <- function(set){
    covered <- integer(node_count)
    for( label in as.character(set) ){
      member_vertices <- label_nodes[[label]]
      if( length(member_vertices) ){
        covered[member_vertices] <- covered[member_vertices] + 1L
      }
    }
    blocked <- which(need > 0L & covered == need)

    return(blocked)
  }
  found <- list()
  expansions <- new.env(hash = TRUE, parent = emptyenv())
  cached_entries <- 0L
  cached_bytes <- 0
  cache_limit <- 8 * 1024^2
  steps <- 0L
  certificate <- "budget"
  reached <- -1L
  root <- .cd_ranked_path(moral, region$x, region$y, blocked_by(base), weight)
  if( is.null(root) ){
    found <- list(base)
    certificate <- "complete"
    reached <- length(base)
  }else if( !any(need[root] > 0L) ){
    certificate <- "complete"
    reached <- length(unique(unlist(node_labels)))
  }else{
    for( bound in seq.int(length(base), max(length(base), length(label_nodes))) ){
      truncated <- FALSE
      exhausted <- FALSE
      visited <- new.env(hash = TRUE, parent = emptyenv())
      expand <- function(set){
        if( exhausted ){
          return(invisible())
        }
        key <- paste0("k", paste(set, collapse = ","))
        if( exists(key, envir = visited, inherits = FALSE) ){
          return(invisible())
        }
        assign(key, TRUE, envir = visited)
        if( any(vapply(found, function(existing_option){
          return(all(existing_option %in% set))
        }, logical(1))) ){
          return(invisible())
        }
        steps <<- steps + 1L
        if( steps > step_budget ){
          exhausted <<- TRUE
          return(invisible())
        }
        expansion <- if( exists(key, envir = expansions, inherits = FALSE) ){
          get(key, envir = expansions, inherits = FALSE)
        }else{
          path <- .cd_ranked_path(moral, region$x, region$y, blocked_by(set), weight)
          moves <- if( !is.null(path) ) unique(lapply(path[need[path] > 0L], function(vertex){
            members <- node_labels[[vertex]]
            return(if( plain_members ) members[!members %in% set] else setdiff(members, set))
          })) else list()
          value <- list(separated = is.null(path), candidates = lapply(moves, function(move){
            return(if( !plain_members ){
              sort(c(set, move))
            }else if( !length(move) ){
              set
            }else if( length(move) == 1L ){
              c(set[set < move], move, set[set > move])
            }else{
              sort.int(c(set, move), method = "quick")
            })
          }))
          if( cached_entries < 4096L && cached_bytes < cache_limit ){
            bytes <- as.numeric(utils::object.size(list(key, value)))
            if( bytes <= cache_limit - cached_bytes ){
              assign(key, value, envir = expansions)
              cached_entries <<- cached_entries + 1L
              cached_bytes <<- cached_bytes + bytes
            }
          }
          value
        }
        if( expansion$separated ){
          found[[length(found) + 1L]] <<- set
          return(invisible())
        }
        for( candidate in expansion$candidates ){
          if( length(candidate) > bound ){
            truncated <<- TRUE
            next
          }
          expand(candidate)
        }
        return(invisible(NULL))
      }
      expand(base)
      if( exhausted ){
        break
      }
      reached <- bound
      if( !truncated ){
        certificate <- "complete"
        break
      }
    }
  }
  realised <- lapply(found, function(set){
    allowed <- region$removable[vapply(labels[region$removable], function(vertex_members){
      return(all(vertex_members %in% set))
    }, logical(1))]
    separator <- union(.cd_ranked_separator(moral, region$x, region$y, allowed), region$must)
    separator_nodes <- graph$nodes[sort(separator)]

    return(separator_nodes)
  })
  ranked <- list(options = found, realised = realised, certificate = certificate,
    complete_up_to = reached, steps = min(steps, step_budget))

  return(ranked)
}


#' Resolve search vertices to recorded source members
#'
#' @param prepared Construction preparation with a map of occurrences to members.
#' @param nodes Exact vertex names in search order.
#' @return One character member vector per vertex. An unmapped name stands for
#'   itself. Callers separately clear the explicitly indexed synthetic vertices.
#' @noRd
.cd_ranked_labels <- function(prepared, nodes){
  labels <- lapply(nodes, function(node){
    members <- unique(prepared$map$member_id[prepared$map$node %in% node])
    concept_members <- if( length(members) ){
      members
    }else{
      node
    }

    return(concept_members)
  })

  return(labels)
}


#' Search and validate ranked adjustment options for one effect
#'
#' @param prepared Construction preparation supplying selection graph and identities.
#' @param effect `"total"` or `"direct"`.
#' @param steps Ranked-search step budget.
#' @param validate Exact-source validator accepting sets and effect.
#' @return A list containing exact sets, a candidate table, completeness and the
#'   ranking certificate. `omitted_equal_smallest` reports whether additional
#'   options of the smallest size may exist. Options that fail source validation
#'   are dropped with a warning, and the remaining ranking is not certified.
#' @noRd
.cd_ranked_search <- function(prepared, effect, steps = .cd_ranked_steps(), validate = .gd_validator(prepared$selection)){
  region <- if( effect == "total" ){
    .cd_ranked_total_region(prepared$selection)
  }else{
    .cd_ranked_direct_region(prepared$selection)
  }
  if( is.null(region) ){
    stop("The installed dagitty engine did not report direct effect constraints.", call. = FALSE)
  }
  members <- .cd_ranked_labels(prepared, region$graph$nodes)
  members[region$graph$synthetic] <- list(character())
  source_order <- prepared$record$source_order
  if( !.query_order_applies(names(prepared$selection), source_order,
                          unname(.gd_node_ids(prepared$selection))) ){
    source_order <- NULL
  }
  vocabulary <- .query_sort(unique(as.character(unlist(members[region$removable]))), source_order)
  labels <- lapply(members, function(vertex_members){
    return(match(vertex_members, vocabulary))
  })
  budget <- steps
  ranked <- .cd_ranked_enumerate(region, labels, budget)
  sets <- ranked$realised
  valid <- if( length(sets) ){
    validate(sets, effect)
  }else{
    logical()
  }
  complete <- identical(ranked$certificate, "complete")
  dropped <- 0L
  if( !all(valid) ){
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
  order <- order(sizes, vapply(sets, function(adjustment_set){
    ids <- unname(.gd_node_ids(prepared$selection)[adjustment_set])
    key <- paste(.query_order_key(adjustment_set, source_order, ids), collapse = "\r")

    return(key)
  }, character(1)), method = "radix")
  sets <- unname(sets[order])
  table <- .cd_candidates(sets, prepared)
  search <- list(sets = sets, complete = complete, table = table,
    omitted_equal_smallest = !complete && length(sizes) > 0L && ranked$complete_up_to < min(sizes),
    ranked = list(method = "quantity-ranked-v1", certificate = if( complete ){
      "complete"
    }else{
      ranked$certificate
    },
    complete_up_to = ranked$complete_up_to, steps = ranked$steps, step_budget = budget,
    dropped = dropped))

  return(search)
}
