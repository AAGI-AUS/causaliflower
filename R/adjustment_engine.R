# ==============================================================================
# adjustment_engine.R


# Call the installed dagitty separator and preserve its result order
# ==============================================================================


#' Enumerate minimal sets with the installed dagitty separator
#'
#' @param dag Graph with endpoint tags and optional recorded occurrence order.
#' @param effect `"total"` or `"direct"`.
#' @param limit Maximum number of separator results, or `Inf`. The direct-effect adapter
#'   passes this bound to dagitty's underlying separator routine.
#' @return The engine's list of sets in result order, with original node names
#'   restored after any ordering rename. The private adapter requires dagitty
#'   0.3-4 and its verified separator API. Unsupported installations raise errors.
#' @noRd
.gd_adjustment_engine <- function(dag, effect, limit){
  original_names <- NULL
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  if( !is.null(source_order) ){
    nodes <- names(dag)
    ordered_names <- .query_order_key(nodes, source_order, .gd_node_ids(dag))
    if( !identical(ordered_names, nodes) ){
      if( anyDuplicated(ordered_names) ){
        stop("The query's recorded occurrence order does not identify unique nodes.", call. = FALSE)
      }
      original_names <- stats::setNames(nodes, ordered_names)
      dag <- copy_rename_native_nodes(dag, stats::setNames(ordered_names, nodes),
                                     "Construction enumeration")
    }
  }
  restore_members <- function(value){
    if( is.null(original_names) ){
      return(value)
    }
    if( is.list(value) ){
      for( index in seq_along(value) ){
        value[[index]] <- restore_members(value[[index]])
      }
      return(value)
    }
    if( is.character(value) ){
      restored <- unname(original_names[value])
      if( anyNA(restored) ){
        stop("Construction enumeration returned an unrecognised internal node.", call. = FALSE)
      }
      attributes(restored) <- attributes(value)
      return(restored)
    }
    return(value)
  }
  use_private_adapter <- isTRUE(getOption("causaliflower.accelerate", TRUE)) # Retain the reference total effect call and its graph type errors
  if( effect == "total" && (!use_private_adapter ||
      !dagitty::graphType(dag) %in% c("dag", "mag", "pdag", "pag")) ){
    result <- dagitty::adjustmentSets(dag, type = "minimal", effect = "total", max.results = limit)
    return(restore_members(result))
  }
  if( utils::packageVersion("dagitty") != "0.3.4" ){
    stop("Construction enumeration currently supports dagitty 0.3-4 only. ",
      "The installed version is ", as.character(utils::packageVersion("dagitty")),
      ".", call. = FALSE)
  }
  script <- system.file("js", "dagitty-alg.js", package = "dagitty")
  if( !nzchar(script) ){
    stop("The installed dagitty JavaScript engine is unavailable.", call. = FALSE)
  }
  context <- V8::v8()
  context$source(script)
  context$eval(paste0(
    "var gd_supported = typeof DAGitty.GraphParser.parseGuess === 'function' && ",
    "typeof DAGitty.GraphAnalyzer.listMsasDirectEffect === 'function' && ",
    "typeof DAGitty.GraphAnalyzer.listMsasTotalEffect === 'function' && ",
    "DAGitty.GraphAnalyzer.listMinimalSeparators.length === 4;"))
  if( !isTRUE(context$get("gd_supported")) ){
    stop("The installed dagitty engine lacks the verified separator API.", call. = FALSE)
  }
  context$assign("gd_graph_text", as.character(dag))
  if( is.finite(limit) ){
    context$assign("gd_limit", limit)
  }else{
    context$eval("var gd_limit=Infinity;") # V8 serializes R NULL as an object
  }
  # Dagitty's direct effect function has no max.results argument, so the
  # separator receives the limit here
  method <- if( effect == "direct" ){ "listMsasDirectEffect" }else{ "listMsasTotalEffect" }
  member_order <- if( effect == "total" ){ ".sort()" }else{ "" }
  context$eval(paste(
    "var gd_result = (function() {",
    "  var g = DAGitty.GraphParser.parseGuess(gd_graph_text);",
    "  var analyzer = DAGitty.GraphAnalyzer;",
    "  var originalSeparators = analyzer.listMinimalSeparators;",
    "  analyzer.listMinimalSeparators = function(graph, must, mustNot, max) {",
    "    var cap = typeof max === 'number' ? Math.min(max, gd_limit) : gd_limit;",
    "    return originalSeparators.call(this, graph, must, mustNot, cap);",
    "  };",
    paste0("  return analyzer.", method, "(g, [], [], gd_limit).map(function(set) {"),
    paste0("    return set.map(function(v) { return v.id; })", member_order, ";"),
    "  });",
    "})();", sep = "\n"))
  context$eval("var gd_count = gd_result.length;")
  count <- context$get("gd_count")
  chunk_size <- 1000L # Bound each JSON transfer while preserving every result and its order
  if( count <= chunk_size ){
    return(restore_members(context$get("gd_result", simplifyVector = FALSE)))
  }
  result <- vector("list", count)
  for( first in seq.int(0, count - 1, by = chunk_size) ){
    last <- min(first + chunk_size, count)
    context$eval(sprintf("var gd_chunk=gd_result.slice(%.0f,%.0f);void 0;", first, last))
    result[seq.int(first + 1, last)] <- context$get("gd_chunk", simplifyVector = FALSE)
  }
  sets <- restore_members(result)

  return(sets)
}
