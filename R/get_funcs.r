# ==============================================================================
# get_funcs.r


# Graph inspection functions and the role classification code they use
# ==============================================================================


#' dagitty node edges and roles
#'
#' `get_edges()` filters a dagitty object and returns a data frame with edges
#' for specified node roles.
#'
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param selected_nodes Node roles whose edges should be returned. Defaults to all roles. Supply
#'   one or more of the following: c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrument", "competing_cause", "collider", "latent", "observed", "undetermined"). Prefix a role with "!" to exclude edges touching that role instead.
#' @returns A data frame of edges for the roles specified in
#'   \code{selected_nodes}.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' edges <- get_edges(dag)
#'
#' @export
get_edges <- function(dag,
                      selected_nodes = c("outcome", "treatment", "confounder",
                                         "mediator", "mediator_outcome_confounder",
                                         "instrument", "competing_cause", "collider",
                                         "latent", "observed", "undetermined")
                      ){
  all_roles <- ALL_NODE_ROLES
  excluded_nodes <- NULL
  if( !identical(selected_nodes, all_roles) ){
    # Separate into exclusions ("!") and inclusions
    is_excluded     <- grepl("^!", selected_nodes)                 # exclusion is a leading-"!" prefix only:
    excluded_tokens <- sub("^!", "", selected_nodes[ is_excluded ]) # a "!" inside a node name is part of the name
    included_tokens <- selected_nodes[ !is_excluded ]

    unknown_roles <- setdiff(c(excluded_tokens, included_tokens), all_roles)
    if( length(unknown_roles) > 0 ){
      stop("The following selected_nodes entries are not valid roles: ", paste(unknown_roles, collapse = ", "), ". Please check the selected_nodes input against the available roles and try again.")
    }

    if( length(excluded_tokens) > 0 && length(included_tokens) > 0 ){
      stop("The selected_nodes input mixes excluded roles (prefixed with \"!\") and included roles in a single call, which is ambiguous. Please supply either all-excluded or all-included role names and try again.")
    }

    if( length(excluded_tokens) > 0 ){
      # exclusion, drop edges touching any excluded role on either side
      excluded_nodes <- excluded_tokens
      selected_nodes <- all_roles[ !all_roles %in% excluded_tokens ]
    }else{
      # positive selection, keep edges touching any included role on either side
      selected_nodes <- included_tokens
    }
  }

  edges <- extract_unique_node_roles(dag) # add ancestor and descendant nodes (calls a function from later in this file)
  edges <- edges_longer(edges)

  if( is.null(excluded_nodes) ){
    edges <- edges[ edges[, "ancestor_role"] %in% selected_nodes | edges[, "descendant_role"] %in% selected_nodes, ]
  }else{
    edges <- edges[ !edges[, "ancestor_role"] %in% excluded_nodes & !edges[, "descendant_role"] %in% excluded_nodes, ]
  }


  return(edges)
}


#' dagitty nodes grouped by role
#'
#' `get_roles()` groups node names by their causal roles in the supplied graph.
#'
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param multiple_roles Defaults to `FALSE` (one role per node). With `TRUE`, a node can be
#'   returned under more than one role, for example as both a latent variable
#'   and a mediator.
#' @return A named list of node names grouped by role. With the default
#'   \code{multiple_roles = FALSE} each node is reported under its primary role
#'   only. \code{get_nodes()} reports every role a node holds. \code{undetermined} is
#'   appended for partially directed graphs.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' get_roles(dag)
#'
#' @export
get_roles <- function(dag, multiple_roles = FALSE){
  roles <- classify_node_roles(dag, multiple_roles = multiple_roles)
  return(roles)
}


#' Classify nodes by role
#'
#' @param dag A dagitty object.
#' @param multiple_roles Whether a node may be returned under more than one
#'   role.
#' @param .cache An environment for cached graph reads, or NULL to create one.
#' @return Named list of node names for each role, with NA for empty roles.
#' @noRd
classify_node_roles <- function(dag, multiple_roles = FALSE, .cache = NULL){
  # Reuse structural reads within a public operation.
  if( is.null(.cache) ) {
    .cache <- graph_cache()
  }

  ## classify edgeless nodes by their declared status ##
  if( nrow( as.data.frame( cached_edges(dag, .cache) ) ) == 0L ){
    node_roles <- ALL_NODE_ROLES
    all_nodes  <- get_node_names(dag)
    role_of    <- ifelse(all_nodes %in% unobserved(dag), "latent",
                  ifelse(all_nodes %in% treatments(dag), "treatment",
                  ifelse(all_nodes %in% get_outcome_names(dag),  "outcome", "observed")))
    roles_list <- lapply(node_roles, function(role){
      nodes <- all_nodes[role_of == role]
      if( length(nodes) == 0L ) {
        nodes <- NA
      }

      return(nodes)
    })
    names(roles_list) <- node_roles
    return(roles_list)
  }

  if(multiple_roles == FALSE){
    edges_wide <- extract_unique_node_roles(dag, .cache = .cache)
  }else{
    edges_wide <- extract_node_roles(dag, .cache = .cache)
  }

  node_roles <- ALL_NODE_ROLES
  num_roles <- length(node_roles)

  ## select role columns by name ##
  edges_wide <- as.data.frame(edges_wide)
  ancestor_cols   <- grep("^ancestor_",   names(edges_wide))
  descendant_cols <- grep("^descendant_", names(edges_wide))

  ## ancestor node edges to list ##
  ancestor_roles <- edges_wide[, c(1:3, ancestor_cols)]
  ancestor_roles <- na.omit( reshape(ancestor_roles, varying = list(4:ncol(ancestor_roles)), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  ancestor_roles <- ancestor_roles[order(ancestor_roles$id), c(1,4)]

  ## descendant node edges to list ##
  descendant_roles <- edges_wide[, c(1:3, descendant_cols)]
  descendant_roles <- na.omit( reshape(descendant_roles, varying = list(4:ncol(descendant_roles)), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  descendant_roles <- descendant_roles[order(descendant_roles$id), c(3,4)]

  ## find missing edges ##
  outcomes <- unique( descendant_roles[ descendant_roles$role == "outcome", ] )
  colliders <- unique( descendant_roles[ descendant_roles$role == "collider", ] )
  observed <- unique( descendant_roles[ descendant_roles$role == "observed", ] )
  latent <- unique( descendant_roles[ descendant_roles$role == "latent", ] )
  instruments <- unique( descendant_roles[ descendant_roles$role == "instrument", ] )
  undetermined_nodes <- unique( descendant_roles[ descendant_roles$role == "undetermined", ] )

  # group by nodes
  names(ancestor_roles)[1] <- "w"
  unique_ancestors <- unique( ancestor_roles[,"w"] ) # vector of unique node names

  missing_outcomes <- outcomes[! unlist(outcomes[,1]) %in% unlist(unique_ancestors), ]
  missing_colliders <- colliders[! unlist(colliders[,1]) %in% unlist(unique_ancestors), ]
  missing_observed <- observed[! unlist(observed[,1]) %in% unlist(unique_ancestors), ]
  missing_latent <- latent[! unlist(latent[,1]) %in% unlist(unique_ancestors), ]
  missing_instruments <- instruments[! unlist(instruments[,1]) %in% unlist(unique_ancestors), ]
  missing_undetermined <- undetermined_nodes[! unlist(undetermined_nodes[,1]) %in% unlist(unique_ancestors), ]

  all_roles <- rbind(ancestor_roles, missing_outcomes, missing_colliders, missing_observed, missing_latent, missing_instruments, missing_undetermined)

  ## include isolated nodes using their declared status ##
  all_nodes <- get_node_names(dag)
  isolated  <- setdiff(all_nodes, as.character(all_roles[, 1]))
  if( length(isolated) > 0L ){
    iso_role <- ifelse(isolated %in% unobserved(dag), "latent",
                ifelse(isolated %in% treatments(dag), "treatment",
                ifelse(isolated %in% get_outcome_names(dag),  "outcome", "observed")))
    iso_rows <- data.frame(isolated, iso_role, stringsAsFactors = FALSE)
    names(iso_rows) <- names(all_roles)
    all_roles[, 1]  <- as.character(all_roles[, 1])
    all_roles <- rbind(all_roles, iso_rows)
  }

  # edges grouped by each unique node in a list
  roles_list <- c()

  roles_list <- suppressWarnings( lapply(1:num_roles, function(x){

    roles_list <-  unique(unlist( all_roles[
      unlist(all_roles[,"role"]) %in% unlist(node_roles[x]) , ][, 1] ))

    return(invisible(roles_list))

  }) )

  names(roles_list) <- node_roles # assign node roles as list element names

  roles_list <- lapply(roles_list, function(nodes){
    if( identical(nodes, character(0)) ) {
      nodes <- NA
    }

    return(nodes)
  })
  return(roles_list)
}


#' dagitty node names and roles
#'
#' `get_nodes()` returns a list of node names and their roles.
#'
#' @importFrom data.table as.data.table is.data.table
#' @param dag A dagitty object.
#' @return A named list giving the roles of each node.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' get_nodes(dag)
#'
#' @export
get_nodes <- function(dag){
  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  ## classify edgeless nodes by their declared status ##
  if( nrow(edges_wide) == 0L ){
    all_nodes <- get_node_names(dag)
    role_of   <- ifelse(all_nodes %in% unobserved(dag), "latent",
                 ifelse(all_nodes %in% treatments(dag), "treatment",
                 ifelse(all_nodes %in% get_outcome_names(dag),  "outcome", "observed")))
    out <- as.list(role_of)
    names(out) <- all_nodes
    return(out)
  }

  ## ancestor node edges to list ##
  edges_wide <- as.data.frame(edges_wide)
  ancestor_cols   <- grep("^ancestor_",   names(edges_wide))
  descendant_cols <- grep("^descendant_", names(edges_wide))
  edges_ancestors <- edges_wide[, c(1:3, ancestor_cols)]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:ncol(edges_ancestors)), idvar = "id",
                                      v.names = "role_ancestor", direction = "long")[,c("v", "e", "w", "role_ancestor", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), 1:4]
  names(edges_ancestors)[1:3] <- c("ancestor", "edge", "descendant")

  ## group by nodes
  unique_ancestors <- unique( unlist( edges_ancestors[,"ancestor"] ) )# vector of unique node names
  num_unique_ancestors <- length(unique_ancestors) # count of unique node names


  # edges grouped by each unique node in a list
  edges_ancestor_list <- c()

  edges_ancestor_list <- suppressWarnings( lapply(seq_len(num_unique_ancestors), function(x){ # seq_len(): 1:0 would index backwards when there are no edges

    edges_ancestor_list <-  unique(unlist( edges_ancestors[
      unlist(edges_ancestors[,"ancestor"]) %in% unlist(unique_ancestors[x]) , ][, c(4)] ))

    return(invisible(edges_ancestor_list))

  }) )

  names(edges_ancestor_list) <- unlist(unique_ancestors) # assign ancestor node role to name of each element in edges list

  ## descendant node edges to list ##
  edges_descendants <- edges_wide[, c(1:3, descendant_cols)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:ncol(edges_descendants)), idvar = "id",
                                        v.names = "role_descendant", direction = "long")[,c("v", "e", "w", "role_descendant", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), c(3,4)]

  ## find missing edges
  outcomes <- unique( edges_descendants[ edges_descendants$role_descendant == "outcome", ] )
  colliders <- unique( edges_descendants[ edges_descendants$role_descendant == "collider", ] )
  observed <- unique( edges_descendants[ edges_descendants$role_descendant == "observed", ] )
  latent <- unique( edges_descendants[ edges_descendants$role_descendant == "latent", ] )
  instruments <- unique( edges_descendants[ edges_descendants$role_descendant == "instrument", ] )
  undetermined <- unique( edges_descendants[ edges_descendants$role_descendant == "undetermined", ] )

  missing_outcomes <- outcomes[! unlist(outcomes[,1]) %in% unlist(unique_ancestors), ]
  missing_colliders <- colliders[! unlist(colliders[,1]) %in% unlist(unique_ancestors), ]
  missing_observed <- observed[! unlist(observed[,1]) %in% unlist(unique_ancestors), ]
  missing_latents <- latent[! unlist(latent[,1]) %in% unlist(unique_ancestors), ]
  missing_instruments <- instruments[! unlist(instruments[,1]) %in% unlist(unique_ancestors), ]
  missing_undetermined <- undetermined[! unlist(undetermined[,1]) %in% unlist(unique_ancestors), ]

  missing_descendant_edges <- rbind(
    missing_outcomes, missing_colliders, missing_observed, missing_latents,
    missing_instruments, missing_undetermined
  )

  unique_descendants <- unique(unlist(missing_descendant_edges[,1]))
  num_unique_descendants <- length(unique_descendants)

  # descendant edges grouped by each unique node in a list
  edges_descendant_list <- c()

  if( num_unique_descendants > 0 ){

    edges_descendant_list <- suppressWarnings( lapply(1:num_unique_descendants, function(x){

      edges_descendant_list <- unique(unlist(
        edges_descendants[ unlist(edges_descendants[,"w"]) %in% unlist(unique_descendants[x]), ][, c(2)] ))

      return(invisible(edges_descendant_list))

    }) )

    names(edges_descendant_list) <- unlist(unique_descendants) # assign descendant node role to name of each element in edges list
  }

  edges_list <- c(edges_ancestor_list, edges_descendant_list) # combine ancestor and descendant node role lists

  ## include isolated nodes using their declared status ##
  isolated <- setdiff(get_node_names(dag), names(edges_list))
  if( length(isolated) > 0L ){
    iso_role <- ifelse(isolated %in% unobserved(dag), "latent",
                ifelse(isolated %in% treatments(dag), "treatment",
                ifelse(isolated %in% get_outcome_names(dag),  "outcome", "observed")))
    iso_list <- as.list(iso_role)
    names(iso_list) <- isolated
    edges_list <- c(edges_list, iso_list)
  }


  return(edges_list)
}


#' Graph structure by focal node
#'
#' `get_structure()` reports each node's role, direct parents and children,
#' nodes joined by undirected or bidirected edges, proper ancestors and
#' descendants, and maximal directed paths. Ancestors are grouped by the number
#' of edges in their shortest directed path to the focal node, and descendants
#' by their shortest directed path from it. Letters
#' beside each ancestor or descendant identify the paths listed under `paths`.
#'
#' @param dag A dagitty object.
#' @param max_paths Maximum paths retained separately for each focal node's
#'   incoming and outgoing collections. A positive whole number or `Inf`.
#' @param max_length Maximum number of directed edges in a retained path. A
#'   positive whole number or `Inf`.
#' @return A nested list indexed by focal node. Each entry contains `role`,
#'   `parents`, `children`, `undirected`, `bidirected`, `ancestors`,
#'   `descendants`, and `paths`. The `truncated` attribute on each `paths`
#'   component records whether `max_paths` or `max_length` limited its incoming
#'   or outgoing paths.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' get_structure(dag)
#'
#' @export
get_structure <- function(dag, max_paths = 100L, max_length = Inf){
  max_paths <- inspection_limit(max_paths, "max_paths", "get_structure()")
  max_length <- inspection_limit(max_length, "max_length", "get_structure()")
  node_names <- inspection_nodes(dag, NULL, "get_structure()")
  edges <- dag_edges(dag, "get_structure()")
  adjacency <- directed_adjacency(edges, node_names)
  node_roles <- get_nodes(dag)
  path_index <- paths_from_edges(
    edges, node_names, node_names, max_paths = max_paths,
    max_length = max_length
  )

  result <- lapply(node_names, function(focal){
    focal_paths <- path_index[[focal]]
    incoming <- attr(focal_paths, "incoming")
    outgoing <- attr(focal_paths, "outgoing")
    paths <- c(incoming, outgoing)
    attr(paths, "direction") <- c(
      stats::setNames(rep("incoming", length(incoming)), names(incoming)),
      stats::setNames(rep("outgoing", length(outgoing)), names(outgoing))
    )
    attr(paths, "truncated") <- attr(focal_paths, "truncated")
    attr(paths, "max_paths") <- max_paths
    attr(paths, "max_length") <- max_length

    node_structure <- list(
      role = node_roles[[focal]],
      parents = adjacency$parents[[focal]],
      children = adjacency$children[[focal]],
      undirected = symmetric_neighbours(
        edges, focal, "--", node_names
      ),
      bidirected = symmetric_neighbours(
        edges, focal, "<->", node_names
      ),
      ancestors = distance_path_groups(
        directed_distances(focal, adjacency$parents), incoming, node_names
      ),
      descendants = distance_path_groups(
        directed_distances(focal, adjacency$children), outgoing, node_names
      ),
      paths = paths
    )

    return(node_structure)
  })
  names(result) <- node_names
  class(result) <- c("causaliflower_structure", "list")
  return(result)
}


#' Label DAG nodes
#'
#' Extract node names, initials, or shortened labels for plotting.
#' @param dag A dagitty object.
#' @param label_type Defaults to \code{"name"} (the graph's node names). Alternatively \code{"initials"} uses the first letter of each word and
#'   \code{"short"} generates shortened unique labels.
#' @return A named character vector of labels, with the original node names as names.
#' @examples
#' dag <- string_to_dag(c(
#'   '"soil moisture" -> "crop growth"',
#'   '"crop growth" -> "grain yield"'
#' ))
#' get_labels(dag, label_type = "initials")
#'
#' @export
get_labels <- function(dag, label_type = "name"){
  if( !label_type %in% c("name", "initials", "short") ){
    stop("Invalid 'label_type'. Use 'name', 'initials', or 'short'.",
         call. = FALSE)
  }

  node_names <- as.character(names(dag))
  if( length(node_names) == 0L ){
    labels_out <- stats::setNames(character(), character())
    return(labels_out)
  }
  if( label_type == "name" ){
    labels_out <- stats::setNames(node_names, node_names)
    return(labels_out)
  }

  words <- strsplit(node_names, "[_[:space:]]+", perl = TRUE)
  labels <- if( label_type == "initials" ){
    vapply(words, function(x) {
      text <- paste0(substr(x, 1L, 1L), collapse = "")

      return(text)
    },
           character(1L))
  }else{
    as.character(unlist(get_label_helper(words), use.names = FALSE))
  }
  labels <- make.unique(labels, sep = "_")
  labels_out <- stats::setNames(labels, node_names)

  return(labels_out)
}


#' Compare dagitty object edges
#'
#' `get_diff_edges()` returns differences between the supplied graphs' edges,
#' optionally including endpoint roles in the comparison.
#'
#' @importFrom data.table setDT
#' @param dag1 A dagitty object.
#' @param dag2 A dagitty object.
#' @param compare_all With `TRUE` (the default), the function returns edges
#'   from either graph that have no match in the other. With `FALSE`, it
#'   returns only edges from `dag1` that have no match in `dag2`.
#' @param include_roles With `FALSE` (the default), edges are compared by
#'   node names and edge types. With `TRUE`, the roles of the two nodes are
#'   also compared.
#' @return Data table of edges.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y", type = "saturated"
#' )
#' dag2 <- build_graph(
#'   variables = c("X1", "X3"), treatments = "T", outcomes = "Y", type = "saturated"
#' )
#' get_diff_edges(dag, dag2)
#'
#' @export
get_diff_edges <- function(dag1, dag2, compare_all = TRUE,  include_roles = FALSE){
  if( include_roles == FALSE){
    edges <- get_edges(dag1)[,c("ancestor", "edge", "descendant")]
    edges2 <- get_edges(dag2)[,c("ancestor", "edge", "descendant")]
  }else{
    edges <- get_edges(dag1)
    edges2 <- get_edges(dag2)
  }

  col_names <- names(edges)

  edges_comb <- data.table::setDT(edges)[!edges2, on = col_names]

  if(compare_all == TRUE){
    edges2 <- data.table::setDT(edges2)[!edges, on = col_names]

    edges_comb <- merge(edges_comb, edges2,
                   by = col_names,
                   all = TRUE)
  }

  return(edges_comb)
}


#' Compare dagitty object node roles
#'
#' `get_diff_roles()` compares the node names assigned to each causal role
#' in the two supplied graphs.
#'
#' @param dag1 A dagitty object.
#' @param dag2 A dagitty object.
#' @param compare_all With `TRUE` (the default), the result includes nodes
#'   assigned to a role in one graph but not the other. With `FALSE`, it
#'   includes only nodes assigned to a role in `dag1` but not that role in
#'   `dag2`.
#' @param nodes_in_common With `TRUE`, the result is restricted to nodes
#'   present in both graphs. The default, `FALSE`, allows nodes found in only
#'   one graph.
#' @return A list of nodes.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y", type = "saturated"
#' )
#' dag2 <- build_graph(
#'   variables = c("X1", "X3"), treatments = "T", outcomes = "Y", type = "saturated"
#' )
#' get_diff_roles(dag, dag2)
#'
#' @export
get_diff_roles <- function(dag1,
                           dag2,
                           compare_all = TRUE,
                           nodes_in_common = FALSE
                           ){
  roles <- get_roles(dag1)
  roles_y <- get_roles(dag2)

  diff_list <- lapply(1:length(roles), function(x){
    nodes <- roles[[x]][ !roles[[x]] %in% roles_y[[x]] ]

    if( compare_all == TRUE){
      nodes <- c( nodes, roles_y[[x]][ !roles_y[[x]] %in% roles[[x]] ] )
    }
    if( nodes_in_common == TRUE){
      dag_names <- names(dag1)
      dag2_names <- names(dag2)
      nodes <- nodes[ nodes %in% dag_names & nodes %in% dag2_names ]
    }
    nodes <- nodes[ complete.cases(nodes) ] # drop NA placeholders from empty role buckets

    return(nodes)
  })

  names(diff_list) <- names(roles)

  diff_list <- Filter(function(x) {
    matches <- length(x) > 0

    return(matches)
  }, diff_list)

  return(diff_list)
}


#' Extract ancestor node edges
#'
#' `get_ancestor_edges()` is deprecated as of 0.0.2. Use
#' `get_node_edges(direction = "any")` to include directed and symmetric
#' edges. The deprecated function keeps its original grouping: each edge is
#' listed once under the first endpoint in the stored edge. `get_node_edges()` instead
#' lists each edge under both endpoints when `direction = "any"`.
#'
#' @param dag A dagitty object.
#' @returns Named list of edges.
#' @examples
#' dag <- string_to_dag(c(
#'   "rainfall -> irrigation",
#'   "rainfall -> yield",
#'   "irrigation -> yield"
#' ))
#' get_ancestor_edges(dag)
#'
#' @export
get_ancestor_edges <- function(dag){
  warning(warningCondition(
    paste0("get_ancestor_edges(): 'get_ancestor_edges()' is deprecated as of ",
           "0.0.2. Use 'get_node_edges(direction = \"any\")' instead."),
    class = "deprecatedWarning"
  ))
  edges <- data.table::as.data.table(dag_edges(dag, "get_ancestor_edges()"))

  ## group by each unique left-hand node, in order of first appearance ##
  unique_ancestors <- unique( edges[["v"]] )
  if( length(unique_ancestors) == 0L ){
    edges_out <- list()
    return(edges_out)
  }
  edges_list <- lapply(unique_ancestors, function(node){
    node_edges <- edges[ edges[["v"]] == node, ]

    return(node_edges)
  })
  names(edges_list) <- unique_ancestors

  return(edges_list)
}


## Role and label helpers ##
#' Extract node roles
#'
#' extract_unique_node_roles() assigns one role to each edge endpoint.
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom dagitty edges coordinates dagitty
#' @param dag A dagitty object.
#' @param .cache Optional environment for cached values from this graph.
#' @returns Edge table with roles for each ancestor and descendant node.
#' @noRd
extract_unique_node_roles <- function(dag, .cache = NULL){
  edges <- cached_edges(dag, .cache)

  if( nrow( as.data.frame(edges) ) == 0 ){
    edges <- data.frame( v = character(0), e = character(0), w = character(0) ) # a graph without edges has no v/e/w columns
  }else{
    edges <- edges[, c("v", "e", "w")]
  }

  # treatment
  treatments <- treatments(dag)

  # outcome
  outcomes <- get_outcome_names(dag)

  # latent variables
  latent_vars <- unobserved(dag)

  # directed reachability sets (BFS over the edge table; "<->" edges do not transmit ancestry)
  ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to")
  ancestors_of_outcome     <- directed_reachable(edges, outcomes, direction = "to")
  descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from")
  ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments)

  # confounders can be common causes, or ancestors of treatment with a directed path to outcome that avoids treatment (VanderWeele & Shpitser 2013)
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- ancestors_of_treatment[ ancestors_of_treatment %in% ancestors_of_outcome_avoiding_treatment &
                                           !ancestors_of_treatment %in% outcomes &
                                           !ancestors_of_treatment %in% latent_vars ]

  # mediators are nodes on a directed path from treatment to outcome
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  mediators <- descendants_of_treatment[ descendants_of_treatment %in% ancestors_of_outcome &
                                           !descendants_of_treatment %in% outcomes &
                                           !descendants_of_treatment %in% confounders &
                                           !descendants_of_treatment %in% latent_vars ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% c(treatments, confounders) & # remove treatment and confounder nodes
                !moc %in% treatment_parents & # double check by removing parents of treatment
                !moc %in% latent_vars &
                !moc %in% mediators &
                !moc %in% outcomes]

  # proxy
  latent_children <- dagitty::children(dag, latent_vars)

  proxy_b <- suppressWarnings( treatment_parents[ treatment_parents %in% latent_children &
                                  !treatment_parents %in% latent_vars &
                                  !treatment_parents %in% treatments ] )# proxy_b (type b confounder)

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c (type a confounder)
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  confounders <- c(confounders, proxy_b, proxy_c) # add to confounders

  # competing cause (computed after the proxy append, so proxy confounders are excluded)
  competing_cause <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% latent_vars &
                                           !outcome_parents %in% moc &
                                           !outcome_parents %in% outcomes ]

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children &
                                   !outcome_children %in% latent_vars &
                                   !outcome_children %in% outcomes &
                                   !outcome_children %in% mediators &
                                   !outcome_children %in% moc &
                                   !outcome_children %in% confounders ]

  # instrumental variables (van der Zander, Textor & Liskiewicz 2015, via dagitty)
  instrumental_vars <- extract_instrumental_variables(dag = dag,
                                                      treatments = treatments,
                                                      outcomes = outcomes)

  # a single role instrument is (a) not in any prior-computed roles (b) never a descendant of treatment (dagitty surfaces conditional instruments that can coincide with a confounder or a post-treatment effect, neither of which is an instrument here)
  instrumental_vars <- instrumental_vars[ !instrumental_vars %in% treatments &
                                            !instrumental_vars %in% outcomes &
                                            !instrumental_vars %in% confounders &
                                            !instrumental_vars %in% mediators &
                                            !instrumental_vars %in% moc &
                                            !instrumental_vars %in% competing_cause &
                                            !instrumental_vars %in% colliders &
                                            !instrumental_vars %in% latent_vars &
                                            !instrumental_vars %in% descendants_of_treatment ]

  # filter latent variables
  latent_vars <- latent_vars[ !latent_vars %in% treatments &
                              !latent_vars %in% outcomes &
                              !latent_vars %in% colliders ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% competing_cause &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]

  # Mark roles undetermined when reachability changes after allowing both
  # orientations of "--" edges (Perkovic et al. 2017, "Interpreting and using
  # CPDAGs").
  pos_ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to",   possible = TRUE)
  pos_ancestors_of_outcome     <- directed_reachable(edges, outcomes,   direction = "to",   possible = TRUE)
  pos_descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from", possible = TRUE)
  pos_ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments, possible = TRUE)

  undetermined <- unique( c( setdiff(pos_ancestors_of_treatment,   ancestors_of_treatment),
                             setdiff(pos_ancestors_of_outcome,     ancestors_of_outcome),
                             setdiff(pos_descendants_of_treatment, descendants_of_treatment),
                             setdiff(pos_ancestors_of_outcome_avoiding_treatment,
                                     ancestors_of_outcome_avoiding_treatment) ) )
  undetermined <- undetermined[ !undetermined %in% c(treatments, outcomes, latent_vars) ]

  # a definite role is assigned only when it holds regardless of "--" orientation
  confounders        <- confounders[        !confounders        %in% undetermined ]
  mediators          <- mediators[          !mediators          %in% undetermined ]
  moc                <- moc[                !moc                %in% undetermined ]
  competing_cause <- competing_cause[ !competing_cause %in% undetermined ]
  colliders          <- colliders[          !colliders          %in% undetermined ]
  instrumental_vars  <- instrumental_vars[  !instrumental_vars  %in% undetermined ]
  observed           <- observed[           !observed           %in% undetermined ]


  if( nrow(edges) == 0L ){ # no rows return here
    return(edges)
  }

  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, "v"] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, "v"] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, "v"] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, "v"] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, "v"] %in% mediators &
                             !edges[, "v"] %in% moc ]  <- "mediator"

  edges$ancestor_iv[edges[, "v"] %in% instrumental_vars] <- "instrument"

  edges$ancestor_competing_cause[edges[, "v"] %in% competing_cause &
                                      !edges[, "v"] %in% proxy_c ] <- "competing_cause"

  edges$ancestor_collider[edges[, "v"] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"

  edges$ancestor_undetermined[edges[, "v"] %in% undetermined] <- "undetermined"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, "w"] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, "w"] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, "w"] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, "w"] %in% moc &
                         !edges[, "w"] %in% mediators]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, "w"] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, "w"] %in% instrumental_vars] <- "instrument"

  edges$descendant_competing_cause[edges[, "w"] %in% competing_cause &
                                        !edges[, "w"] %in% proxy_c ] <- "competing_cause"

  edges$descendant_collider[edges[, "w"] %in% colliders] <- "collider"

  edges$descendant_latent[edges$w %in% latent_vars] <- "latent"

  edges$descendant_observed[edges$w %in% observed] <- "observed"

  edges$descendant_undetermined[edges[, "w"] %in% undetermined] <- "undetermined"

  return(edges)
}


#' Helper for classify_node_roles() and get_nodes()
#'
#' extract_node_roles() assigns roles to both endpoints of each edge.
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom dagitty edges coordinates dagitty
#' @param dag A dagitty object.
#' @param .cache Optional environment for cached values from this graph.
#' @returns Edge table with roles for each ancestor and descendant node.
#' @noRd
extract_node_roles <- function(dag, .cache = NULL){
  edges <- cached_edges(dag, .cache)

  if( nrow( as.data.frame(edges) ) == 0 ){
    edges <- data.frame( v = character(0), e = character(0), w = character(0) ) # a graph without edges has no v/e/w columns
  }else{
    edges <- edges[, c("v", "e", "w")]
  }

  # treatment
  treatments <- treatments(dag)

  # outcome
  outcomes <- get_outcome_names(dag)

  # latent variables
  latent_vars <- unobserved(dag)

  # directed reachability sets (BFS over the edge table, "<->" edges do not transmit ancestry)
  ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to")
  ancestors_of_outcome     <- directed_reachable(edges, outcomes, direction = "to")
  descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from")
  ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments)

  # confounders can be common causes, or ancestors of treatment with a directed path to outcome that avoids treatment (VanderWeele & Shpitser 2013)
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- ancestors_of_treatment[ ancestors_of_treatment %in% ancestors_of_outcome_avoiding_treatment ]

  # mediators are nodes on a directed path from treatment to outcome
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  mediators <- descendants_of_treatment[ descendants_of_treatment %in% ancestors_of_outcome ]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% confounders & # remove treatment and confounder nodes
                !moc %in% treatment_parents ] # double check by removing parents of treatment

  # proxy

  latent_children <- dagitty::children(dag, latent_vars)

  proxy_b <- suppressWarnings( treatment_parents[ treatment_parents %in% latent_children &
                                                    !treatment_parents %in% latent_vars &
                                                    !treatment_parents %in% treatments ] )# proxy_b (type b confounder)

  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c (type c confounder)
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars &
                                !outcome_parents %in% outcomes ]
  confounders <- c(confounders, proxy_b, proxy_c) # add to confounders

  # competing cause (computed after the proxy append, so proxy confounders are excluded)
  competing_cause <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% outcomes ]

  # collider
  outcome_children <- dagitty::children(dag, outcomes)

  colliders <- outcome_children[ outcome_children %in% treatment_children ]

  # instrumental variables (van der Zander, Textor & Liskiewicz 2015, via dagitty)
  instrumental_vars <- extract_instrumental_variables(dag = dag,
                                                      treatments = treatments,
                                                      outcomes = outcomes)

  # filter latent variables
  latent_vars <- latent_vars[  !latent_vars %in% treatments &
                                 !latent_vars %in% outcomes ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% competing_cause &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]


  # Mark roles undetermined when reachability changes after allowing both
  # orientations of "--" edges (Perkovic et al. 2017, "Interpreting and using
  # CPDAGs").
  pos_ancestors_of_treatment   <- directed_reachable(edges, treatments, direction = "to",   possible = TRUE)
  pos_ancestors_of_outcome     <- directed_reachable(edges, outcomes,   direction = "to",   possible = TRUE)
  pos_descendants_of_treatment <- directed_reachable(edges, treatments, direction = "from", possible = TRUE)
  pos_ancestors_of_outcome_avoiding_treatment <- directed_reachable(edges, outcomes, direction = "to", avoiding = treatments, possible = TRUE)

  undetermined <- unique( c( setdiff(pos_ancestors_of_treatment,   ancestors_of_treatment),
                             setdiff(pos_ancestors_of_outcome,     ancestors_of_outcome),
                             setdiff(pos_descendants_of_treatment, descendants_of_treatment),
                             setdiff(pos_ancestors_of_outcome_avoiding_treatment,
                                     ancestors_of_outcome_avoiding_treatment) ) )
  undetermined <- undetermined[ !undetermined %in% c(treatments, outcomes, latent_vars) ]

  # a definite role is assigned only when it holds regardless of "--" orientation
  confounders        <- confounders[        !confounders        %in% undetermined ]
  mediators          <- mediators[          !mediators          %in% undetermined ]
  moc                <- moc[                !moc                %in% undetermined ]
  competing_cause <- competing_cause[ !competing_cause %in% undetermined ]
  colliders          <- colliders[          !colliders          %in% undetermined ]
  instrumental_vars  <- instrumental_vars[  !instrumental_vars  %in% undetermined ]
  observed           <- observed[           !observed           %in% undetermined ]

  if( nrow(edges) == 0L ){ # no rows return here
    return(edges)
  }

  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, "v"] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, "v"] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, "v"] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, "v"] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, "v"] %in% mediators ]  <- "mediator"

  edges$ancestor_iv[edges[, "v"] %in% instrumental_vars] <- "instrument"

  edges$ancestor_competing_cause[edges[, "v"] %in% competing_cause &
                                      !edges[, "v"] %in% proxy_c ] <- "competing_cause"

  edges$ancestor_collider[edges[, "v"] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"

  edges$ancestor_undetermined[edges[, "v"] %in% undetermined] <- "undetermined"


  # assign roles for all v in edges
  edges$descendant_outcome[edges[, "w"] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, "w"] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, "w"] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, "w"] %in% moc ]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, "w"] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, "w"] %in% instrumental_vars] <- "instrument"

  edges$descendant_competing_cause[edges[, "w"] %in% competing_cause &
                                        !edges[, "w"] %in% proxy_c ] <- "competing_cause"

  edges$descendant_collider[edges[, "w"] %in% colliders] <- "collider"

  edges$descendant_latent[edges$w %in% latent_vars] <- "latent"

  edges$descendant_observed[edges$w %in% observed] <- "observed"

  edges$descendant_undetermined[edges[, "w"] %in% undetermined] <- "undetermined"

  return(edges)
}


#' Helper for get_edges()
#'
#' @importFrom data.table as.data.table data.table
#' @param edges Edge data frame from extract_unique_node_roles().
#' @returns A data frame with ancestor, edge, descendant, ancestor_role
#'   and descendant_role columns.
#' @noRd
edges_longer <- function(edges){
  if( nrow(edges) == 0 ){
    edges_out <- data.frame( ancestor = character(0),
                             edge = character(0),
                             descendant = character(0),
                             ancestor_role = character(0),
                             descendant_role = character(0) )
    return(edges_out)
  }

  edges$id <- 1:nrow(edges)

  ancestor_cols   <- grep("^ancestor_",   names(edges))
  descendant_cols <- grep("^descendant_", names(edges))
  edges_ancestors <- edges[, c(1:3, ancestor_cols)]
  edges_descendants <- edges[, c(1:3, descendant_cols)]

  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:ncol(edges_ancestors)), idvar = "id",
                                      v.names = "ancestor_role", direction = "long")[,c("v", "e", "w", "ancestor_role", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), ]

  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:ncol(edges_descendants)), idvar = "id",
                                        v.names = "descendant_role", direction = "long")[,c("v", "e", "w", "descendant_role", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), ]

  # Pair ancestor and descendant roles by the original edge row.
  edges <- merge(edges_ancestors, edges_descendants, by = c("id", "v", "e", "w"), all = TRUE)
  edges <- edges[order(edges$id), c("v", "e", "w", "ancestor_role", "descendant_role")]

  names(edges) <- c("ancestor", "edge", "descendant", "ancestor_role", "descendant_role")
  rownames(edges) <- NULL

  return(edges)
}


#' Helper for get_labels()
#'
#' @param names_list A list of character vectors, one per node name split into words.
#' @return A list of shortened node names.
#' @noRd
get_label_helper <- function(names_list){

  # mostly hard coded (for now) but it should give decent results for common agricultural inputs e.g. location = loc, temperature = temp
  initials_list <- lapply( 1:length(names_list), function(x){
    if(length(names_list[[x]]) == 1){ # single word node labels

      if( nchar(names_list[[x]]) %% 2 == 0 ){ # if even char length
        if(nchar(names_list[[x]]) > 7){
          text <- substr(names_list[[x]], 1, 3)

          return(text)
        }else if (nchar(names_list[[x]]) == 4){
          text <- substr(names_list[[x]], 1, 4)

          return(text)
        }else{
          text <- substr(names_list[[x]], 1, 2)

          return(text)
        }
      }else{ # if odd char length
        if(nchar(names_list[[x]]) > 6){
          text <- substr(names_list[[x]], 1, 4)

          return(text)
        }else if (nchar(names_list[[x]]) == 5){
          text <- substr(names_list[[x]], 1, 2)

          return(text)
        }else{
          text <- substr(names_list[[x]], 1, 1)

          return(text)
        }
      }

    }else if(length(names_list[[x]]) == 2){ # two word labels

      if(nchar(names_list[[x]][[1]]) == 1){
        text <- paste( c(names_list[[x]][[1]], substr(names_list[[x]][[2]], 1, 2) ), collapse = "_" )

        return(text)
      }else if(nchar(names_list[[x]][[1]]) < 5){
        if(nchar(names_list[[x]][[2]]) > 6){
          text <- paste( c(names_list[[x]][[1]], substr(names_list[[x]][[2]], 1, 1) ), collapse = "_" )

          return(text)
        }else{
          text <- paste( c(names_list[[x]][[1]], substr(names_list[[x]][[2]], 1, 4) ), collapse = "_" )

          return(text)
        }
      }else if(nchar(names_list[[x]][[2]]) < 5){
        if(nchar(names_list[[x]][[1]]) > 6){
          text <- paste( c(substr(names_list[[x]][[1]], 1, 1), names_list[[x]][[2]] ), collapse = "_" )

          return(text)
        }else{
          text <- paste( c(substr(names_list[[x]][[1]], 1, 4), names_list[[x]][[2]] ), collapse = "_" )

          return(text)
        }
      }else{
        text <- paste( substr(names_list[[x]], 1,3), collapse = "_" )

        return(text)
      }

    }else{ # three or greater word labels

      if(any(nchar(names_list[[x]]) > 6)){

        if(all(nchar(names_list[[x]]) > 11)){
          text <- paste( substr(names_list[[x]], 1, 1), collapse = "" )

          return(text)
        }else if(any(nchar(names_list[[x]]) > 8)){
          text <- paste( substr(names_list[[x]], 1, 3), collapse = "_" )

          return(text)
        }else{
          text <- paste( substr(names_list[[x]], 1, 2), collapse = "_" )

          return(text)
        }

      }else if(all(nchar(names_list[[x]]) > 4)){
        text <- paste( substr(names_list[[x]], 1,2), collapse = "" )

        return(text)
      }else{
        text <- paste( substr(names_list[[x]], 1, 1), collapse = "" )

        return(text)
      }

    }
  })

  return(initials_list)
}


## Paths and graph structure ##
#' Check the focal nodes for graph inspection
#'
#' @param dag A dagitty object.
#' @param nodes Character vector of node names.
#' @param caller Name of the calling function for messages.
#' @return Unique focal node names, or all graph node names when `nodes` is
#'   `NULL`.
#' @noRd
inspection_nodes <- function(dag, nodes, caller){
  dag_graph_type(dag)
  graph_nodes <- get_node_names(dag)

  if( is.null(nodes) ) {
    return(graph_nodes)
  }
  if( !is.character(nodes) || anyNA(nodes) || any(!nzchar(nodes)) ){
    stop(caller, ": 'nodes' must contain non-missing, non-empty node names.",
         call. = FALSE)
  }

  missing_nodes <- setdiff(nodes, graph_nodes)
  if( length(missing_nodes) > 0L ){
    stop(caller, ": nodes not in the graph: ",
         paste(missing_nodes, collapse = ", "), ".", call. = FALSE)
  }

  nodes_out <- unique(nodes)

  return(nodes_out)
}


#' Check a graph inspection limit
#'
#' @param value Path-count or path-length limit to check.
#' @param argument Argument name for messages.
#' @param caller Name of the calling function for messages.
#' @return A positive numeric whole number or `Inf`.
#' @noRd
inspection_limit <- function(value, argument, caller){
  valid <- is.numeric(value) && !is.logical(value) && length(value) == 1L &&
    !is.na(value) && value > 0 &&
    (is.infinite(value) || (is.finite(value) && value == floor(value)))
  if( !valid ){
    stop(caller, ": '", argument,
         "' must be a positive whole number or Inf.", call. = FALSE)
  }
  if( is.infinite(value) ) {
    limit <- Inf

    return(limit)
  } else {
    limit <- as.numeric(value)

    return(limit)
  }
}


#' Group directed parents and children by node
#'
#' @param edges Edge table with columns `v`, `e`, and `w`.
#' @param node_names Graph node names in their display order.
#' @return A list containing named `parents` and `children` lists.
#' @noRd
directed_adjacency <- function(edges, node_names){
  directed <- edges[edges$e == "->", c("v", "w"), drop = FALSE]

  in_node_order <- function(values){
    values <- unique(as.character(values))
    nodes <- node_names[node_names %in% values]

    return(nodes)
  }

  parents <- lapply(node_names, function(node){
    parents <- in_node_order(directed$v[directed$w == node])

    return(parents)
  })
  children <- lapply(node_names, function(node){
    children <- in_node_order(directed$w[directed$v == node])

    return(children)
  })
  names(parents) <- node_names
  names(children) <- node_names

  adjacency <- list(parents = parents, children = children)

  return(adjacency)
}


#' Enumerate maximal paths from a focal node
#'
#' @param focal Name of the focal node.
#' @param adjacency Named adjacency list for the traversal direction.
#' @param direction `"incoming"` to prepend parents or `"outgoing"` to append
#'   children.
#' @param max_paths Maximum retained paths in each direction, or `Inf`.
#' @param max_length Maximum edges in each retained path, or `Inf`.
#' @return A list containing node-vector `paths` and `truncated` limit flags.
#' @noRd
enumerate_maximal_paths <- function(focal, adjacency, direction,
                                     max_paths, max_length){
  paths <- list()
  path_limited <- FALSE
  length_limited <- FALSE
  stack <- list(list(node = focal, path = focal))

  while( length(stack) > 0L && !path_limited ){
    state <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL

    next_nodes <- adjacency[[state$node]]
    next_nodes <- next_nodes[!next_nodes %in% state$path]
    path_length <- length(state$path) - 1L
    at_length_limit <- is.finite(max_length) && path_length >= max_length

    if( at_length_limit && length(next_nodes) > 0L ){
      length_limited <- TRUE
      next
    }
    terminal <- length(next_nodes) == 0L && path_length > 0L

    if( terminal ){
      if( length(paths) < max_paths ){
        paths[[length(paths) + 1L]] <- state$path
      }else{
        path_limited <- TRUE
      }
      next
    }

    if( length(next_nodes) == 0L ) {
      next
    }

    ## reverse stack insertion to visit paths in graph node order ##
    for( next_node in rev(next_nodes) ){
      next_path <- if( direction == "incoming" ){
        c(next_node, state$path)
      }else{
        c(state$path, next_node)
      }
      stack[[length(stack) + 1L]] <- list(node = next_node,
                                         path = next_path)
    }
  }

  paths_out <- list(
    paths = paths,
    truncated = c(max_paths = path_limited, max_length = length_limited)
  )

  return(paths_out)
}


#' Assign letter keys to numbered paths
#'
#' @param index Positive whole-number path indices.
#' @return A character vector of keys such as `"a"`, `"z"`, and `"aa"`.
#' @noRd
path_key <- function(index){
  keys <- vapply(index, function(value){
    key <- character()
    while( value > 0L ){
      remainder <- (value - 1L) %% 26L
      key <- c(letters[remainder + 1L], key)
      value <- (value - 1L) %/% 26L
    }
    text <- paste0(key, collapse = "")

    return(text)
  }, character(1L))

  return(keys)
}


#' Select paths using a predicate
#'
#' @param paths Named list of node vectors describing paths.
#' @param predicate Function returning one logical value for a node-vector
#'   path.
#' @return The named subset of paths satisfying the predicate.
#' @noRd
subset_paths <- function(paths, predicate){
  if( length(paths) == 0L ) {
    paths_out <- list()
    return(paths_out)
  }
  paths_out <- paths[vapply(paths, predicate, logical(1L))]

  return(paths_out)
}


#' Index a focal node's paths by other nodes
#'
#' @param focal Name of the focal node.
#' @param node_names Graph node names in their display order.
#' @param incoming Named list of paths ending at the focal node.
#' @param outgoing Named list of paths starting at the focal node.
#' @param truncated Incoming and outgoing flags for path-count and path-length
#'   truncation.
#' @return A `causaliflower_path_index` list with incoming, outgoing, and
#'   truncation attributes.
#' @noRd
focal_path_index <- function(focal, node_names, incoming, outgoing,
                              truncated){
  all_paths <- c(incoming, outgoing)
  other_nodes <- node_names[node_names != focal &
    vapply(node_names, function(node){
      matches <- any(vapply(all_paths, function(path) {
        matches <- node %in% path

        return(matches)
      }, logical(1L)))

      return(matches)
    }, logical(1L))]

  index <- lapply(other_nodes, function(other){
    incoming_other <- subset_paths(incoming, function(path) {
      matches <- other %in% path

      return(matches)
    })
    outgoing_other <- subset_paths(outgoing, function(path) {
      matches <- other %in% path

      return(matches)
    })
    any_other <- c(incoming_other, outgoing_other)

    node_paths <- list(
      incoming = incoming_other,
      outgoing = outgoing_other,
      starts_at = subset_paths(any_other,
                                function(path) {
                                  matches <- identical(path[[1L]], other)

                                  return(matches)
                                }),
      ends_at = subset_paths(any_other, function(path){
        matches <- identical(path[[length(path)]], other)

        return(matches)
      }),
      any = any_other
    )

    return(node_paths)
  })
  names(index) <- other_nodes
  class(index) <- c("causaliflower_path_index", "list")
  attr(index, "incoming") <- incoming
  attr(index, "outgoing") <- outgoing
  attr(index, "truncated") <- truncated

  return(index)
}


#' Build the indexed path collection from graph edges
#'
#' @param edges Edge table with columns `v`, `e`, and `w`.
#' @param node_names Graph node names in their display order.
#' @param focal_nodes Names of nodes whose incoming and outgoing paths are
#'   collected.
#' @param max_paths Maximum retained paths in each direction, or `Inf`.
#' @param max_length Maximum edges in each retained path, or `Inf`.
#' @return A `causaliflower_paths` list grouped by focal node, with path-limit
#'   attributes.
#' @noRd
paths_from_edges <- function(edges, node_names, focal_nodes,
                              max_paths, max_length){
  adjacency <- directed_adjacency(edges, node_names)

  result <- lapply(focal_nodes, function(focal){
    incoming <- enumerate_maximal_paths(
      focal, adjacency$parents, "incoming", max_paths, max_length
    )
    outgoing <- enumerate_maximal_paths(
      focal, adjacency$children, "outgoing", max_paths, max_length
    )

    all_paths <- c(incoming$paths, outgoing$paths)
    keys <- path_key(seq_along(all_paths))
    incoming_count <- length(incoming$paths)
    incoming_paths <- if( incoming_count > 0L ){
      stats::setNames(incoming$paths, keys[seq_len(incoming_count)])
    }else{
      list()
    }
    outgoing_paths <- if( length(outgoing$paths) > 0L ){
      positions <- incoming_count + seq_along(outgoing$paths)
      stats::setNames(outgoing$paths, keys[positions])
    }else{
      list()
    }

    path_index <- focal_path_index(
      focal, node_names, incoming_paths, outgoing_paths,
      list(incoming = incoming$truncated, outgoing = outgoing$truncated)
    )

    return(path_index)
  })
  names(result) <- focal_nodes
  class(result) <- c("causaliflower_paths", "list")
  attr(result, "max_paths") <- max_paths
  attr(result, "max_length") <- max_length

  return(result)
}


#' Directed paths grouped by focal node
#'
#' `get_paths()` returns maximal directed paths into and out of each selected
#' focal node. Incoming paths run from an upstream root to the focal node, and
#' outgoing paths run from the focal node to a downstream leaf. Only `->` edges
#' are followed. Paths never repeat a node.
#'
#' @param dag A dagitty object.
#' @param nodes Focal node names. `NULL` uses every node in the graph.
#' @param max_paths Maximum paths retained separately for each focal node's
#'   incoming and outgoing collections. A positive whole number or `Inf`.
#' @param max_length Maximum number of directed edges in a retained path. A
#'   positive whole number or `Inf`.
#' @return A nested list indexed first by focal node and then by another node on
#'   one of its incoming or outgoing paths. Each entry for another node contains
#'   `incoming`, `outgoing`, `starts_at`, `ends_at`, and `any` path lists.
#'   The `truncated` attribute of the focal node's entry records whether
#'   `max_paths` or `max_length` excluded incoming or outgoing paths.
#' @examples
#' dag <- string_to_dag(c("A -> B", "B -> C"))
#' paths <- get_paths(dag)
#' paths$B$A$incoming
#' paths$B$C$outgoing
#' @export
get_paths <- function(dag, nodes = NULL, max_paths = 100L,
                      max_length = Inf){
  max_paths <- inspection_limit(max_paths, "max_paths", "get_paths()")
  max_length <- inspection_limit(max_length, "max_length", "get_paths()")
  focal_nodes <- inspection_nodes(dag, nodes, "get_paths()")
  node_names <- get_node_names(dag)
  edges <- dag_edges(dag, "get_paths()")

  paths <- paths_from_edges(edges, node_names, focal_nodes, max_paths, max_length)

  return(paths)
}


# Exact `$` lookup avoids partial matching of node names.
#' Read a graph inspection entry by its exact name
#'
#' @param x A `causaliflower_paths` object.
#' @param name Exact node name to retrieve.
#' @return The matching list element, or `NULL` when no exact name is present.
#' @noRd
#' @export
`$.causaliflower_paths` <- function(x, name){
  paths <- unclass(x)[[name, exact = TRUE]]

  return(paths)
}


# Exact `$` lookup avoids partial matching of node names.
#' Read a graph inspection entry by its exact name
#'
#' @param x A `causaliflower_path_index` object.
#' @param name Exact node name to retrieve.
#' @return The matching list element, or `NULL` when no exact name is present.
#' @noRd
#' @export
`$.causaliflower_path_index` <- function(x, name){
  paths <- unclass(x)[[name, exact = TRUE]]

  return(paths)
}


#' Print named directed paths
#'
#' @param paths Named list of node vectors describing paths.
#' @return Invisible `NULL`.
#' @noRd
print_path_rows <- function(paths){
  if( length(paths) == 0L ){
    cat("  None\n")
    out <- NULL
    return(invisible(out))
  }
  for( key in names(paths) ){
    cat("  (", key, ") ", paste(paths[[key]], collapse = " -> "), "\n",
        sep = "")
  }
  out <- NULL

  return(invisible(out))
}


#' Print path truncation notices
#'
#' @param truncated Incoming and outgoing flags for path-count and path-length
#'   truncation.
#' @param max_paths Maximum retained paths in each direction, or `Inf`.
#' @param max_length Maximum edges in each retained path, or `Inf`.
#' @param indent Text prepended to each notice.
#' @return Invisible `NULL`.
#' @noRd
print_path_truncation <- function(truncated, max_paths, max_length,
                                   indent = "  "){
  for( direction in c("incoming", "outgoing") ){
    label <- paste0(toupper(substr(direction, 1L, 1L)),
                    substring(direction, 2L))
    if( isTRUE(truncated[[direction]][["max_paths"]]) ){
      cat(indent, label, " paths truncated at max_paths = ", max_paths,
          ".\n", sep = "")
    }
    if( isTRUE(truncated[[direction]][["max_length"]]) ){
      cat(indent, label, " paths truncated at max_length = ", max_length,
          ".\n", sep = "")
    }
  }
  paths <- NULL

  return(invisible(paths))
}


#' Print graph directed paths
#'
#' @param x A `causaliflower_paths` object.
#' @param ... Additional arguments, unused.
#' @return The supplied object, invisibly.
#' @noRd
#' @export
print.causaliflower_paths <- function(x, ...){
  if( length(x) == 0L ){
    cat("No focal nodes.\n")
    return(invisible(x))
  }

  for( focal in names(x) ){
    focal_paths <- x[[focal]]
    cat(focal, "\n", sep = "")
    cat("Incoming\n")
    print_path_rows(attr(focal_paths, "incoming"))
    cat("Outgoing\n")
    print_path_rows(attr(focal_paths, "outgoing"))
    print_path_truncation(
      attr(focal_paths, "truncated"), attr(x, "max_paths"),
      attr(x, "max_length")
    )
    if( focal != tail(names(x), 1L) ) {
      cat("\n")
    }
  }

  return(invisible(x))
}


#' Incident graph edges grouped by focal node
#'
#' `get_node_edges()` groups incident edges under each selected focal node.
#' `direction = "any"` includes directed, bidirected and undirected edges.
#' `"incoming"` and `"outgoing"` include only directed arrows whose head or
#' tail is the focal node. Undirected and bidirected edges are excluded from
#' those directional subsets.
#'
#' @param dag A dagitty object.
#' @param nodes Focal node names. `NULL` uses every node in the graph.
#' @param direction One of `"any"`, `"incoming"`, or `"outgoing"`.
#' @return A named list of `v`, `e`, `w` edge tables, including empty tables for
#'   selected nodes without matching edges.
#' @examples
#' dag <- string_to_dag(c("A -> B", "B -> C"))
#' get_node_edges(dag, "B")
#' get_node_edges(dag, "B", direction = "incoming")
#' @export
get_node_edges <- function(dag, nodes = NULL, direction = "any"){
  if( !is.character(direction) || length(direction) != 1L ||
      is.na(direction) || !direction %in% c("any", "incoming", "outgoing") ){
    stop("get_node_edges(): 'direction' must be one of: any, incoming, outgoing.",
         call. = FALSE)
  }

  focal_nodes <- inspection_nodes(dag, nodes, "get_node_edges()")
  edges <- data.table::as.data.table(dag_edges(dag, "get_node_edges()"))

  result <- lapply(focal_nodes, function(focal){
    keep <- if( direction == "incoming" ){
      edges$e == "->" & edges$w == focal
    }else if( direction == "outgoing" ){
      edges$e == "->" & edges$v == focal
    }else{
      edges$v == focal | edges$w == focal
    }
    node_edges <- edges[keep, c("v", "e", "w"), with = FALSE]

    return(node_edges)
  })
  names(result) <- focal_nodes
  return(result)
}


#' Find minimum directed distances from a focal node
#'
#' @param focal Name of the focal node.
#' @param adjacency Named adjacency list for the traversal direction.
#' @return Named numeric distances to reachable nodes, excluding the focal
#'   node.
#' @noRd
directed_distances <- function(focal, adjacency){
  distance <- stats::setNames(0, focal)
  frontier <- focal
  current_distance <- 0L

  while( length(frontier) > 0L ){
    current_distance <- current_distance + 1L
    candidates <- unique(unlist(adjacency[frontier], use.names = FALSE))
    candidates <- candidates[!candidates %in% names(distance)]
    if( length(candidates) == 0L ) {
      break
    }
    distance[candidates] <- current_distance
    frontier <- candidates
  }

  distance_out <- distance[names(distance) != focal]

  return(distance_out)
}


#' Group reachable nodes and their paths by distance
#'
#' @param distances Named minimum directed distances for reachable nodes.
#' @param paths Named list of node vectors describing paths.
#' @param node_names Graph node names in their display order.
#' @return A list grouped by distance, with named path keys for each node.
#' @noRd
distance_path_groups <- function(distances, paths, node_names){
  if( length(distances) == 0L ) {
    groups_out <- list()
    return(groups_out)
  }

  levels <- sort(unique(as.integer(distances)))
  groups <- lapply(levels, function(distance){
    members <- node_names[node_names %in%
      names(distances)[as.integer(distances) == distance]]
    keys <- lapply(members, function(node){
      path_keys <- names(paths)[vapply(paths, function(path) {
        matches <- node %in% path

        return(matches)
      }, logical(1L))]

      return(path_keys)
    })
    names(keys) <- members

    return(keys)
  })
  names(groups) <- paste0("d", levels)

  return(groups)
}


#' Find neighbours joined by a symmetric edge
#'
#' @param edges Edge table with columns `v`, `e`, and `w`.
#' @param focal Name of the focal node.
#' @param edge_type Symmetric edge type, `"--"` or `"<->"`.
#' @param node_names Graph node names in their display order.
#' @return Neighbour names in graph node order, excluding the focal node.
#' @noRd
symmetric_neighbours <- function(edges, focal, edge_type, node_names){
  selected <- edges[edges$e == edge_type &
                      (edges$v == focal | edges$w == focal), , drop = FALSE]
  neighbours <- ifelse(selected$v == focal, selected$w, selected$v)
  neighbours <- unique(as.character(neighbours))
  nodes <- node_names[node_names %in% neighbours & node_names != focal]

  return(nodes)
}


#' Print nodes grouped by directed distance
#'
#' @param label Heading to print.
#' @param groups Distance groups containing node names and their path keys.
#' @return Invisible `NULL`.
#' @noRd
print_distance_groups <- function(label, groups){
  cat(label, "\n", sep = "")
  if( length(groups) == 0L ){
    cat("  None\n")
    out <- NULL
    return(invisible(out))
  }

  for( group_name in names(groups) ){
    group <- groups[[group_name]]
    entries <- vapply(names(group), function(node){
      keys <- group[[node]]
      key_text <- if( length(keys) > 0L ){
        paste0(" ", paste0("(", keys, ")", collapse = ""))
      }else{
        ""
      }
      text <- paste0("[", node, key_text, "]")

      return(text)
    }, character(1L))
    cat("  ", sub("^d", "", group_name), ". ",
        paste(entries, collapse = " "), "\n", sep = "")
  }
  out <- NULL

  return(invisible(out))
}


#' Print graph node roles and connections
#'
#' @param x A `causaliflower_structure` object.
#' @param ... Additional arguments, unused.
#' @return The supplied object, invisibly.
#' @noRd
#' @export
print.causaliflower_structure <- function(x, ...){
  if( length(x) == 0L ){
    cat("No nodes.\n")
    return(invisible(x))
  }

  print_values <- function(label, values){
    printed <- cat(label, ": ", if( length(values) > 0L ){
               paste(values, collapse = ", ")
             }else{
               "None"
             }, "\n", sep = "")

    return(invisible(printed))
  }

  for( focal in names(x) ){
    node <- x[[focal]]
    cat(focal, "\n", sep = "")
    print_values("Role", node$role)
    print_values("Parents", node$parents)
    print_values("Children", node$children)
    print_values("Undirected", node$undirected)
    print_values("Bidirected", node$bidirected)
    print_distance_groups("Ancestors", node$ancestors)
    print_distance_groups("Descendants", node$descendants)
    cat("Paths\n")
    print_path_rows(node$paths)
    print_path_truncation(
      attr(node$paths, "truncated"), attr(node$paths, "max_paths"),
      attr(node$paths, "max_length")
    )
    if( focal != tail(names(x), 1L) ) {
      cat("\n")
    }
  }

  return(invisible(x))
}


# Exact `$` lookup avoids partial matching of node names.
#' Read a graph inspection entry by its exact name
#'
#' @param x A `causaliflower_structure` object.
#' @param name Exact node name to retrieve.
#' @return The matching list element, or `NULL` when no exact name is present.
#' @noRd
#' @export
`$.causaliflower_structure` <- function(x, name){
  value <- unclass(x)[[name, exact = TRUE]]

  return(value)
}
