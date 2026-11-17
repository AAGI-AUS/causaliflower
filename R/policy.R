# ==============================================================================
# policy.R


# Apply the reviewed recovery interpretation without changing prepared source data
# ==============================================================================


#' Describe the reviewed recovery exceptions
#'
#' @return A list with the policy `version`, a `default` admission rule and
#'   exception `rules` named by quantity ID.
#' @noRd
.cf_recovery_policy <- function(){
  policy <- list(
    version = "dagby-recovery-20261019-v1",
    default = "yes",
    rules = c(
      "cgr_actual" = "yes:user",
      "unspecified_organ_live_n" = "no:Class-level live N does not identify an organ instance",
      "n_translocation" = "yes:user",
      "canopy_senescence" = "yes:user",
      "frac_radn_intercepted" = "yes:user",
      "spike" = "no:Whole organ label does not specify a numeric quantity",
      "canopy_conductance" = "yes:user"
    )
  )
  return(policy)
}


#' Validate or attach the prepared table's recovery policy
#'
#' @param nodes Nonempty plain data frame with unique, nonmissing character
#'   `node_id` and `concept_id` values and optional `recovered_node` metadata.
#' @return `nodes` unchanged when its `recovered_node` column is valid. If that
#'   column is absent, returns a copy with the reviewed policy added. Invalid
#'   identities or rules raise errors, as do missing quantities required to
#'   apply the reviewed policy.
#' @noRd
.cf_effective_nodes <- function(nodes){
  if( !identical(class(nodes), "data.frame") || !nrow(nodes) ||
      anyDuplicated(names(nodes)) || !all(c("node_id", "concept_id") %in% names(nodes)) ){
    stop("Prepared node schema is incomplete.", call. = FALSE)
  }
  for( field in c("node_id", "concept_id") ){
    identities <- nodes[[field]]
    if( !is.character(identities) || anyNA(identities) ||
        any(!nzchar(identities)) || anyDuplicated(identities) ){
      stop("Invalid node identity column: ", field, call. = FALSE)
    }
  }
  if( "recovered_node" %in% names(nodes) ){
    recovery <- nodes$recovered_node
    if( !is.character(recovery) || anyNA(recovery) ||
        any(!recovery %in% c("yes", "yes:user") & !grepl("^no:[[:space:]]*[^[:space:]].*$", recovery)) ){
      stop("Invalid recovered_node policy in the prepared node table.", call. = FALSE)
    }
    return(nodes)
  }
  policy <- .cf_recovery_policy()
  missing_nodes <- setdiff(names(policy$rules), nodes$node_id)
  if( length(missing_nodes) ){
    stop("Reviewed policy quantities are missing. Review the policy for this table: ",
         paste(missing_nodes, collapse = ", "), call. = FALSE)
  }
  expected <- rep(policy$default, nrow(nodes))
  expected[match(names(policy$rules), nodes$node_id)] <- unname(policy$rules)
  # Replacing data frame values copies the table, preserving the caller's columns and attributes
  nodes$recovered_node <- expected
  return(nodes)
}


#' Bind the same effective node table for namespace and exported data access
#'
#' @param namespace Package namespace environment receiving the effective table.
#' @param lazydata Lazy-data environment containing `apsim_dag_nodes`.
#' @return `NULL`, invisibly, after assigning the validated effective table in
#'   both environments. Missing prepared data raise an error.
#' @noRd
.cf_bind_recovery_policy <- function(namespace,
                                     lazydata = getNamespaceInfo(namespace, "lazydata")){
  name <- "apsim_dag_nodes"
  if( !is.environment(lazydata) || !exists(name, envir = lazydata, inherits = FALSE) ){
    stop("The prepared package node table is missing.", call. = FALSE)
  }
  nodes <- .cf_effective_nodes(get(name, envir = lazydata, inherits = FALSE))
  # R and pkgload populate lazy data before .onLoad. Namespace callers need the same effective table.
  assign(name, nodes, envir = namespace)
  assign(name, nodes, envir = lazydata)
  return(invisible(NULL))
}
