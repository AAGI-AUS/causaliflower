# ==============================================================================
# acceleration.R


# Optional native kernels and bounded reuse of deterministic graph work
# ==============================================================================


#' Read the optional acceleration switch
#'
#' @return `TRUE` if `causaliflower.accelerate` is exactly `TRUE` or unset, and
#'   `FALSE` otherwise.
#' @noRd
.cf_acceleration_enabled <- function(){
  return(isTRUE(getOption("causaliflower.accelerate", TRUE)))
}


#' Bind cached work to character interpretation and collation
#'
#' @return A character vector containing the current `LC_CTYPE` and `LC_COLLATE`
#'   settings, in that order.
#' @noRd
.cf_cache_locale <- function(){
  return(c(Sys.getlocale("LC_CTYPE"), Sys.getlocale("LC_COLLATE")))
}


#' Find an optional native symbol in the package namespace
#'
#' @param name Registered symbol binding to look up without inheritance.
#' @return The symbol binding, or `NULL` if it is absent or acceleration is disabled.
#' @noRd
.cf_native_symbol <- function(name){
  if( !.cf_acceleration_enabled() ){
    return(NULL)
  }
  symbol <- get0(name, envir = environment(.cf_native_symbol), inherits = FALSE)
  return(symbol)
}


#' Recognise graph classes safe for the standard acceleration path
#'
#' @param graph Object whose exact class vector is checked.
#' @return `TRUE` if the class vector is exactly `"dagitty"` or
#'   `c("cf_dag", "dagitty")`, and `FALSE` otherwise.
#' @noRd
.cf_standard_graph <- function(graph){
  standard <- identical(class(graph), "dagitty") ||
    identical(class(graph), c("cf_dag", "dagitty"))
  return(standard)
}


#' Keep bounded session caches for deterministic graph work
#'
#' Cache keys include the arguments, their attributes and the current locale.
#' Serialised outputs prevent callers from changing stored values. An entry is
#' not stored if computation signals a condition, or if its arguments or result
#' contain live references or unsupported classes. The cache limits both entry
#' counts and accounted bytes.
#' @return A list of functions sharing the cache:
#'   - `call(kind, arguments, compute)` returns a cached value or calls the
#'     zero-argument function `compute`. It preserves return visibility.
#'   - `has(kind, arguments)` returns `TRUE` if that entry is cached, and
#'     `FALSE` otherwise.
#'   - `clear()` removes all entries, resets counters and returns `NULL` invisibly.
#'   - `stats()` returns entry counts, byte totals and configured limits.
#'   The supported kinds are `"metadata"`, `"preparation"`, `"levels"` and `"builds"`.
#' @noRd
.cf_cache <- local({
  entries <- list()
  bytes <- 0
  limits <- c(metadata = 8L, preparation = 2L, levels = 8L, builds = 24L) # Maximum retained entries of each kind
  maximum <- 64 * 1024^2 # Maximum accounted cache size in bytes
  overhead <- 1024 # Accounted overhead for each serialised input and output pair
  counts <- c(hits = 0L, misses = 0L, bypasses = 0L)
  inert_classes <- list("dagitty", c("cf_dag", "dagitty"),
    "data.frame", c("data.table", "data.frame"),
    "factor", c("ordered", "factor")) # Classes whose dispatch does not require external state

  safe <- function(value){
    if( typeof(value) %in% c("environment", "externalptr", "weakref", "closure",
                            "builtin", "special", "promise") || isS4(value) ){
      return(FALSE)
    }
    classes <- attr(value, "class", exact = TRUE)
    if( !is.null(classes) && !any(vapply(inert_classes, identical, logical(1), classes)) ){
      return(FALSE)
    }
    if( is.recursive(value) && length(value) ){
      for( index in seq_along(value) ){
        if( !safe(.subset2(value, index)) ){
          return(FALSE)
        }
      }
    }
    value_attributes <- attributes(value)
    attributes_safe <- is.null(value_attributes) || all(vapply(value_attributes, safe, logical(1)))

    return(attributes_safe)
  }

  encode <- function(value){
    encoded <- tryCatch(
      withCallingHandlers(serialize(value, NULL, version = 3),
        warning = function(warning){
          stop(warning)
        },
        message = function(message){
          stop(message)
        }),
      error = function(error){
        return(NULL)
      })
    return(encoded)
  }

  drop <- function(index){
    remaining <- suspendInterrupts({
      bytes <<- bytes - entries[[index]]$bytes
      entries <<- entries[-index]
    })
    return(invisible(remaining))
  }

  call <- function(kind, arguments, compute){
    if( !.cf_acceleration_enabled() ||
        !isTRUE(tryCatch(safe(arguments), error = function(error){
          return(FALSE)
        })) ){
      counts[["bypasses"]] <<- counts[["bypasses"]] + 1L
      return(compute())
    }
    input <- encode(list(locale = .cf_cache_locale(), arguments = arguments))
    if( is.null(input) || length(input) + overhead > maximum ){
      counts[["bypasses"]] <<- counts[["bypasses"]] + 1L
      return(compute())
    }
    for( index in seq_along(entries) ){
      saved <- entries[[index]]
      if( identical(saved$kind, kind) && identical(saved$input, input) ){
        counts[["hits"]] <<- counts[["hits"]] + 1L
        entries <<- c(entries[-index], list(saved))
        value <- unserialize(saved$output)
        if( value$visible ){
          return(value$value)
        }
        return(invisible(value$value))
      }
    }
    counts[["misses"]] <<- counts[["misses"]] + 1L
    signalled <- FALSE
    value <- withCallingHandlers(withVisible(compute()), condition = function(condition){
      signalled <<- TRUE
      return(invisible(NULL))
    })
    output <- NULL
    if( !signalled && isTRUE(tryCatch(safe(value), error = function(error){
      return(FALSE)
    })) ){
      output <- encode(value)
    }
    size <- Inf
    if( !is.null(output) ){
      size <- length(input) + length(output) + overhead
    }
    if( !is.finite(size) || size > maximum ){
      counts[["bypasses"]] <<- counts[["bypasses"]] + 1L
    }else{
      same <- which(vapply(entries, function(entry){
        return(identical(entry$kind, kind))
      }, logical(1)))
      while( length(same) >= limits[[kind]] ){
        drop(same[[1L]])
        same <- which(vapply(entries, function(entry){
          return(identical(entry$kind, kind))
        }, logical(1)))
      }
      while( length(entries) && bytes + size > maximum ){ # Evict stored builds first. A build never displaces another cache entry type.
        builds <- which(vapply(entries, function(entry){
          return(identical(entry$kind, "builds"))
        }, logical(1)))
        if( length(builds) ){
          drop(builds[[1L]])
        }else if( identical(kind, "builds") ){
          break
        }else{
          drop(1L)
        }
      }
      if( bytes + size > maximum ){
        counts[["bypasses"]] <<- counts[["bypasses"]] + 1L
      }else{
        suspendInterrupts({
          entries[[length(entries) + 1L]] <<- list(kind = kind, input = input,
            output = output, bytes = size)
          bytes <<- bytes + size
        })
      }
    }
    if( value$visible ){
      return(value$value)
    }
    return(invisible(value$value))
  }

  has <- function(kind, arguments){
    if( !.cf_acceleration_enabled() ||
        !isTRUE(tryCatch(safe(arguments), error = function(error){
          return(FALSE)
        })) ){
      return(FALSE)
    }
    input <- encode(list(locale = .cf_cache_locale(), arguments = arguments))
    if( is.null(input) ){
      return(FALSE)
    }
    for( saved in entries ){
      if( identical(saved$kind, kind) && identical(saved$input, input) ){
        return(TRUE)
      }
    }
    return(FALSE)
  }

  clear <- function(){
    suspendInterrupts({
      entries <<- list()
      bytes <<- 0
      counts[] <<- 0L
    })
    return(invisible(NULL))
  }

  stats <- function(){
    entries_by_kind <- vapply(names(limits), function(kind){
      return(sum(vapply(entries, function(entry){
        return(identical(entry$kind, kind))
      }, logical(1))))
    }, integer(1))
    statistics <- list(counts = counts, entries = length(entries), accounted_bytes = bytes,
      maximum_bytes = maximum, entries_by_kind = entries_by_kind, maximum_entries = limits)

    return(statistics)
  }

  list(call = call, has = has, clear = clear, stats = stats)
})
