# ==============================================================================
# adjustment_engine.R


# Call the installed dagitty separator and preserve its result order
# ==============================================================================


.gd_adjustment_engine <- function(dag, effect, limit, progress = NULL) { # Enumerate through the installed dagitty implementation
  original_names <- NULL
  source_order <- attr(dag, "cf_source_order", exact = TRUE)
  if (!is.null(source_order)) {
    nodes <- names(dag)
    ordered_names <- .query_order_key(nodes, source_order, .gd_node_ids(dag))
    if (!identical(ordered_names, nodes)) {
      if (anyDuplicated(ordered_names)) {
        stop("The query's recorded occurrence order does not identify unique nodes.", call. = FALSE)
      }
      original_names <- stats::setNames(nodes, ordered_names)
      dag <- copy_rename_native_nodes(dag, stats::setNames(ordered_names, nodes),
                                     "Construction enumeration")
    }
  }
  restore_members <- function(value) { # Restore public names without reordering results or dropping their attributes
    if (is.null(original_names)) {
      return(value)
    }
    if (is.list(value)) {
      for (index in seq_along(value)) {
        value[[index]] <- restore_members(value[[index]])
      }
      return(value)
    }
    if (is.character(value)) {
      restored <- unname(original_names[value])
      if (anyNA(restored)) {
        stop("Construction enumeration returned an unrecognised internal node.", call. = FALSE)
      }
      attributes(restored) <- attributes(value)
      return(restored)
    }
    return(value)
  }
  use_private_adapter <- isTRUE(getOption("causaliflower.accelerate", TRUE)) # Retain the reference total effect call and its graph type errors
  if (effect == "total" && ((!use_private_adapter && is.null(progress)) ||
      !dagitty::graphType(dag) %in% c("dag", "mag", "pdag", "pag"))) {
    result <- dagitty::adjustmentSets(dag, type = "minimal", effect = "total", max.results = limit)
    return(restore_members(result))
  }
  if (utils::packageVersion("dagitty") != "0.3.4") {
    stop("Construction enumeration currently supports dagitty 0.3-4 only. ",
      "The installed version is ", as.character(utils::packageVersion("dagitty")),
      ".", call. = FALSE)
  }
  script <- system.file("js", "dagitty-alg.js", package = "dagitty")
  if (!nzchar(script)) {
    stop("The installed dagitty JavaScript engine is unavailable.", call. = FALSE)
  }
  context <- V8::v8()
  context$source(script)
  if (!is.null(progress)) {
    code <- paste(readLines(script, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    separator <- context$eval("DAGitty.GraphAnalyzer.listMinimalSeparators.toString()")
    pattern <- "([[:alnum:]_$]+=function\\([^)]*\\)\\{)if\\(!\\(([[:alnum:]_$]+)\\.length>=[[:alnum:]_$]+\\)\\)" # Find the pinned implementation's capped recursion entry
    locations <- gregexpr(pattern, separator, perl = TRUE)[[1L]]
    if (length(locations) != 1L || locations[[1L]] < 1L) {
      stop("The separator source does not match the verified progress hook.", call. = FALSE)
    }
    matched <- regmatches(separator, regexec(pattern, separator, perl = TRUE))[[1L]]
    anchor <- matched[[1L]]
    opening <- matched[[2L]]
    results_name <- matched[[3L]]
    locations <- gregexpr(anchor, code, fixed = TRUE)[[1L]]
    if (length(locations) != 1L || locations[[1L]] < 1L) {
      stop("The separator source does not match the verified progress hook.", call. = FALSE)
    }
    guard <- substring(anchor, nchar(opening) + 1L)
    replacement <- paste0(opening, "gd_progress_tick(", results_name, ".length);", guard)
    instrumented <- sub(anchor, replacement, code, fixed = TRUE)
    # Reload the installed script to preserve private lexical bindings. The tick
    # observes existing counts without changing traversal, membership or stopping.
    context$eval(instrumented)
    context$assign("gd_effect", effect)
    context$eval(paste(
      "var gd_steps = 0;",
      "var gd_started = Date.now();",
      "var gd_last = gd_started;",
      "function gd_progress_tick(count) {",
      "  gd_steps++;",
      "  if ((gd_steps & 255) === 0 && Date.now() - gd_last >= 2000) {",
      "    console.log(gd_effect + ' search [' + ((Date.now() - gd_started) / 1000).toFixed(1) +",
      "      's]: ' + count + ' alternatives found');",
      "    gd_last = Date.now();",
      "  }",
      "}", sep = "\n"))
  }
  context$eval(paste0(
    "var gd_supported = typeof DAGitty.GraphParser.parseGuess === 'function' && ",
    "typeof DAGitty.GraphAnalyzer.listMsasDirectEffect === 'function' && ",
    "typeof DAGitty.GraphAnalyzer.listMsasTotalEffect === 'function' && ",
    "DAGitty.GraphAnalyzer.listMinimalSeparators.length === 4;"))
  if (!isTRUE(context$get("gd_supported"))) {
    stop("The installed dagitty engine lacks the verified separator API.", call. = FALSE)
  }
  context$assign("gd_graph_text", as.character(dag))
  if (is.finite(limit)) {
    context$assign("gd_limit", limit)
  } else {
    context$eval("var gd_limit=Infinity;") # V8 serializes R NULL as an object
  }
  # The installed direct effect API omits max.results. This wrapper applies the
  # limit used to find a witness. Complete construction uses Inf.
  method <- if (effect == "direct") { "listMsasDirectEffect" } else { "listMsasTotalEffect" }
  member_order <- if (effect == "total") { ".sort()" } else { "" }
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
  if (count <= chunk_size) {
    return(restore_members(context$get("gd_result", simplifyVector = FALSE)))
  }
  result <- vector("list", count)
  for (first in seq.int(0, count - 1, by = chunk_size)) {
    last <- min(first + chunk_size, count)
    context$eval(sprintf("var gd_chunk=gd_result.slice(%.0f,%.0f);void 0;", first, last))
    result[seq.int(first + 1, last)] <- context$get("gd_chunk", simplifyVector = FALSE)
  }
  return(restore_members(result))
}
