# ==============================================================================
# acceleration.R


# Optional native kernels and bounded reuse of deterministic graph work
# ==============================================================================


.cf_acceleration_enabled <- function() { # Read the shared switch for caches and native kernels
  return(isTRUE(getOption("causaliflower.accelerate", TRUE)))
}


.cf_cache_locale <- function() { # Include case conversion and ordering conventions in cache inputs
  return(c(Sys.getlocale("LC_CTYPE"), Sys.getlocale("LC_COLLATE")))
}


.cf_native_symbol <- function(name) { # Find a registered kernel, retaining the R fallback when unavailable
  if (!.cf_acceleration_enabled()) {
    return(NULL)
  }
  symbol <- get0(name, envir = environment(.cf_native_symbol), inherits = FALSE)
  return(symbol)
}


.cf_standard_graph <- function(graph) { # Restrict cached graph dispatch to ordinary package classes
  standard <- identical(class(graph), "dagitty") ||
    identical(class(graph), c("cf_dag", "dagitty")) ||
    identical(class(graph), c("dag_by_query", "dagitty"))
  return(standard)
}


# Serialized inputs include every attribute. Serialized outputs prevent callers
# from changing cached results. Calls that signal conditions or contain live
# references run without caching.
.cf_cache <- local({
  entries <- list()
  bytes <- 0
  limits <- c(metadata = 8L, preparation = 2L) # Maximum retained entries of each kind
  maximum <- 64 * 1024^2 # Maximum accounted cache size in bytes
  overhead <- 1024 # Accounted overhead for each serialized input/output pair
  counts <- c(hits = 0L, misses = 0L, bypasses = 0L)
  inert_classes <- list("dagitty", c("cf_dag", "dagitty"),
    c("dag_by_query", "dagitty"), "data.frame", c("data.table", "data.frame"),
    "factor", c("ordered", "factor")) # Classes whose dispatch does not require external state

  safe <- function(value) { # Reject live references, custom classes and attributes that depend on external state
    if (typeof(value) %in% c("environment", "externalptr", "weakref", "closure",
                            "builtin", "special", "promise") || isS4(value)) {
      return(FALSE)
    }
    classes <- attr(value, "class", exact = TRUE)
    if (!is.null(classes) && !any(vapply(inert_classes, identical, logical(1), classes))) {
      return(FALSE)
    }
    if (is.recursive(value) && length(value)) {
      for (index in seq_along(value)) {
        if (!safe(.subset2(value, index))) {
          return(FALSE)
        }
      }
    }
    value_attributes <- attributes(value)
    return(is.null(value_attributes) || all(vapply(value_attributes, safe, logical(1))))
  }

  encode <- function(value) { # Cache only serialization that completes without conditions
    encoded <- tryCatch(
      withCallingHandlers(serialize(value, NULL, version = 3),
        warning = function(warning) {
          stop(warning)
        },
        message = function(message) {
          stop(message)
        }),
      error = function(error) {
        return(NULL)
      })
    return(encoded)
  }

  drop <- function(index) { # Evict one entry and update accounting without an intervening interrupt
    remaining <- suspendInterrupts({
      bytes <<- bytes - entries[[index]]$bytes
      entries <<- entries[-index]
    })
    return(invisible(remaining))
  }

  call <- function(kind, arguments, compute) { # Reuse exact inputs while preserving output visibility and conditions
    if (!.cf_acceleration_enabled() ||
        !isTRUE(tryCatch(safe(arguments), error = function(error) {
          return(FALSE)
        }))) {
      counts[["bypasses"]] <<- counts[["bypasses"]] + 1L
      return(compute())
    }
    input <- encode(list(locale = .cf_cache_locale(), arguments = arguments))
    if (is.null(input) || length(input) + overhead > maximum) {
      counts[["bypasses"]] <<- counts[["bypasses"]] + 1L
      return(compute())
    }
    for (index in seq_along(entries)) {
      saved <- entries[[index]]
      if (identical(saved$kind, kind) && identical(saved$input, input)) {
        counts[["hits"]] <<- counts[["hits"]] + 1L
        entries <<- c(entries[-index], list(saved))
        value <- unserialize(saved$output)
        if (value$visible) {
          return(value$value)
        }
        return(invisible(value$value))
      }
    }
    counts[["misses"]] <<- counts[["misses"]] + 1L
    signalled <- FALSE
    value <- withCallingHandlers(withVisible(compute()), condition = function(condition) {
      signalled <<- TRUE
      return(invisible(NULL))
    })
    output <- NULL
    if (!signalled && isTRUE(tryCatch(safe(value), error = function(error) {
      return(FALSE)
    }))) {
      output <- encode(value)
    }
    size <- Inf
    if (!is.null(output)) {
      size <- length(input) + length(output) + overhead
    }
    if (!is.finite(size) || size > maximum) {
      counts[["bypasses"]] <<- counts[["bypasses"]] + 1L
    } else {
      same <- which(vapply(entries, function(entry) {
        return(identical(entry$kind, kind))
      }, logical(1)))
      while (length(same) >= limits[[kind]]) {
        drop(same[[1L]])
        same <- which(vapply(entries, function(entry) {
          return(identical(entry$kind, kind))
        }, logical(1)))
      }
      while (length(entries) && bytes + size > maximum) {
        drop(1L)
      }
      suspendInterrupts({
        entries[[length(entries) + 1L]] <<- list(kind = kind, input = input,
          output = output, bytes = size)
        bytes <<- bytes + size
      })
    }
    if (value$visible) {
      return(value$value)
    }
    return(invisible(value$value))
  }

  clear <- function() { # Clear entries, accounting and counters together
    suspendInterrupts({
      entries <<- list()
      bytes <<- 0
      counts[] <<- 0L
    })
    return(invisible(NULL))
  }

  stats <- function() { # Report cache usage and its fixed bounds
    entries_by_kind <- vapply(names(limits), function(kind) {
      return(sum(vapply(entries, function(entry) {
        return(identical(entry$kind, kind))
      }, logical(1))))
    }, integer(1))
    return(list(counts = counts, entries = length(entries), accounted_bytes = bytes,
      maximum_bytes = maximum, entries_by_kind = entries_by_kind, maximum_entries = limits))
  }

  list(call = call, clear = clear, stats = stats)
})
