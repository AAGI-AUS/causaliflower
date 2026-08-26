

#' dagitty node edges and roles
#'
#' get_edges() filters a dagitty object and returns a data frame with edges for specified node roles.
#'
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param selected_nodes Nodes to return edges. Defaults to all roles, or can be a character or vector combination of any of the following: c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrument", "competing_cause", "collider", "latent", "observed", "undetermined"). Prefix a role with "!" to exclude edges touching that role instead.
#' @returns A data frame of edges for the roles specified in
#'   \code{selected_nodes}.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
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
  all_roles <- .ALL_NODE_ROLES
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
#' @importFrom data.table data.table
#' @param dag A dagitty object.
#' @param multiple_roles Defaults to FALSE (one role per node). If set to TRUE, multiple roles can be returned for some nodes (e.g. latent mediator variable).
#' @return Nested list of nodes and node relationships. With the default
#'   \code{multiple_roles = FALSE} each node is reported under its primary role
#'   only; \code{get_nodes()} reports every role a node holds. \code{undetermined} is
#'   appended for partially directed graphs.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' get_roles(dag)
#'
#' @export
get_roles <- function(dag, multiple_roles = FALSE){
  roles <- .get_roles(dag, multiple_roles = multiple_roles)
  roles
}

.get_roles <- function(dag, multiple_roles = FALSE, .cache = NULL){
  # Reuse structural reads within a public operation.
  if( is.null(.cache) ) .cache <- .new_cache()

  ## edge-less graphs: extract_unique_node_roles()/extract_node_roles() assign roles by
  ## projecting onto the edge table, which cannot operate on a 0-row table (get_roles never
  ## reaches its reshape). With no edges a node can hold no structural role, so classify
  ## every node directly from .get_node_names() and the dagitty-status primitives:
  ## unobserved -> latent, else exposure -> treatment, else outcome -> outcome, else observed.
  if( nrow( as.data.frame( .cached_edges(dag, .cache) ) ) == 0L ){
    node_roles <- .ALL_NODE_ROLES
    all_nodes  <- .get_node_names(dag)
    role_of    <- ifelse(all_nodes %in% unobserved(dag), "latent",
                  ifelse(all_nodes %in% treatments(dag), "treatment",
                  ifelse(all_nodes %in% .outcomes(dag),  "outcome", "observed")))
    roles_list <- lapply(node_roles, function(r){ v <- all_nodes[role_of == r]; if( length(v) > 0L ) v else NA })
    names(roles_list) <- node_roles
    return(roles_list)
  }

  if(multiple_roles == FALSE){
    edges_wide <- extract_unique_node_roles(dag, .cache = .cache)
  }else{
    edges_wide <- extract_node_roles(dag, .cache = .cache)
  }

  node_roles <- .ALL_NODE_ROLES
  num_roles <- length(node_roles)

  ## select role columns by name so the fixed v/e/w prefix and any added role column (e.g.
  ## "undetermined") are handled without relying on hard-coded column positions
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

  ## fold in isolated nodes -- present in the graph (.get_node_names(), which includes nodes with
  ## no edges) but absent from the edge table, since extract_*() projects roles onto edges and
  ## a node in no edge has no row. A disconnected node holds no structural role, so classify it
  ## by status: unobserved -> latent, else exposure -> treatment, else outcome -> outcome, else observed.
  all_nodes <- .get_node_names(dag)
  isolated  <- setdiff(all_nodes, as.character(all_roles[, 1]))
  if( length(isolated) > 0L ){
    iso_role <- ifelse(isolated %in% unobserved(dag), "latent",
                ifelse(isolated %in% treatments(dag), "treatment",
                ifelse(isolated %in% .outcomes(dag),  "outcome", "observed")))
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

  }) )

  names(roles_list) <- node_roles # assign node roles as list element names

  roles_list <- lapply( roles_list, function(x) if( identical(x, character(0)) ) NA else x )
  return(roles_list)
}


#' dagitty node names and roles
#'
#' get_nodes() is a dagitty wrapper that returns a list of node names and their roles.
#'
#' @importFrom data.table as.data.table is.data.table
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' get_nodes(dag)
#'
#' @export
get_nodes <- function(dag){
  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  ## fully edge-less graph: every node is isolated and reshape() cannot take a zero-row table,
  ## so classify all nodes by status here (same idiom as get_roles + the fold-in below)
  if( nrow(edges_wide) == 0L ){
    all_nodes <- .get_node_names(dag)
    role_of   <- ifelse(all_nodes %in% unobserved(dag), "latent",
                 ifelse(all_nodes %in% treatments(dag), "treatment",
                 ifelse(all_nodes %in% .outcomes(dag),  "outcome", "observed")))
    out <- as.list(role_of); names(out) <- all_nodes
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

  edges_ancestor_list <- suppressWarnings( lapply(seq_len(num_unique_ancestors), function(x){ # seq_len(): 1:0 would index backwards on an edge-less dag

    edges_ancestor_list <-  unique(unlist( edges_ancestors[
      unlist(edges_ancestors[,"ancestor"]) %in% unlist(unique_ancestors[x]) , ][, c(4)] ))

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

    }) )

    names(edges_descendant_list) <- unlist(unique_descendants) # assign descendant node role to name of each element in edges list
  }

  edges_list <- c(edges_ancestor_list, edges_descendant_list) # combine ancestor and descendant node role lists

  ## fold in isolated nodes -- present in the graph but in no edge (extract_*() projects roles
  ## onto edges, so a node in no edge has no row); classify by status, as get_roles() does
  isolated <- setdiff(.get_node_names(dag), names(edges_list))
  if( length(isolated) > 0L ){
    iso_role <- ifelse(isolated %in% unobserved(dag), "latent",
                ifelse(isolated %in% treatments(dag), "treatment",
                ifelse(isolated %in% .outcomes(dag),  "outcome", "observed")))
    iso_list <- as.list(iso_role); names(iso_list) <- isolated
    edges_list <- c(edges_list, iso_list)
  }


  return(edges_list)
}


#' dagitty node names, ancestor/descendant roles
#'
#' get_structure() is a dagitty wrapper function that returns a list of node names extracted from a dagitty object, including each node's direct edges and its own role. For each node, the ancestor component holds the edges where the node is the source and the descendant component the edges where it is the target; the role column reports the focal node's role, and the edge operator is not reported. Symmetric edges appear once, on their stored orientation. Isolated nodes are included with empty components.
#'
#' @importFrom data.table as.data.table is.data.table
#' @param dag A dagitty object.
#' @return Nested list of nodes and node relationships
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' get_structure(dag)
#'
#' @export
get_structure <- function(dag){
  edges_wide <- extract_node_roles(dag) # extract node roles from daggity object (wide format)

  if( nrow(edges_wide) == 0L ){
    empty_edges <- data.frame(
      ancestor = character(), descendant = character(), role = character()
    )
    node_names <- .get_node_names(dag)
    structure <- lapply(node_names, function(node){
      list(ancestor = empty_edges, descendant = empty_edges)
    })
    names(structure) <- node_names
    return(structure)
  }

  unique_nodes <- unname(unlist(  unique(c(edges_wide[,1],edges_wide[,3])) ))
  unique_nodes <- c(unique_nodes, setdiff(.get_node_names(dag), unique_nodes)) # include isolated nodes (in the graph but in no edge)
  num_unique_nodes <- length(unique_nodes)

  ## ancestor node edges to list ##
  edges_wide <- as.data.frame(edges_wide)
  ancestor_cols   <- grep("^ancestor_",   names(edges_wide))
  descendant_cols <- grep("^descendant_", names(edges_wide))
  edges_ancestors <- edges_wide[, c(1:3, ancestor_cols)]
  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:ncol(edges_ancestors)), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "w", "role", "id")] )
  names(edges_ancestors)[1:2] <- c("ancestor", "descendant")
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), ][,c(1:3)]

  ## descendant node edges to list ##
  edges_descendants <- edges_wide[, c(1:3, descendant_cols)]
  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:ncol(edges_descendants)), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "w", "role", "id")] )
  names(edges_descendants)[1:2] <- c("ancestor", "descendant")
  edges_descendants <- edges_descendants[order(edges_descendants$id), ][,c(1:3)]

  # ancestor edges grouped by each unique node in a list
  edges_ancestor_list <- suppressWarnings( lapply(1:num_unique_nodes, function(x){

    edges_ancestors[ unlist(edges_ancestors[,"ancestor"]) %in% unlist(unique_nodes[x]), ]

  }) )

  # descendant edges grouped by each unique node in a list
  edges_descendant_list <- suppressWarnings( lapply(1:num_unique_nodes, function(x){

    edges_descendants[ unlist(edges_descendants[,"descendant"]) %in% unlist(unique_nodes[x]), ]

  }) )


  ## combine ancestor and descendant nodes ##
  edges_list <- Map(rbind, edges_descendant_list, edges_ancestor_list)

  names(edges_list) <- unlist(unique_nodes)

  ancestors_descendants_label_vec <- c("ancestor", "descendant")
  ## next create list with each node as an element
  edges_structure_list <- c()

  edges_structure_list <- suppressWarnings( lapply( 1:num_unique_nodes , function(x){

    edges_structure_list <- lapply(1:2, function(n){ # edges grouped by unique node name

      edges_structure_list[x][[n]] <- edges_list[[x]][ edges_list[[x]][[n]] %in% names(edges_list[x]), ]

    })

    names(edges_structure_list) <- ancestors_descendants_label_vec

    edges_structure_list

  } ) )

  names(edges_structure_list) <- unlist(unique_nodes) # assign ancestor node role to name of each element in edges list

  return(edges_structure_list)
}


#' Label DAG nodes
#' Extract node names, initials, or shortened labels for plotting.
#' @param dag A dagitty object.
#' @param label_type Defaults to \code{"name"} (the graph's node names);
#'   alternatively \code{"initials"} uses the first letter of each word and
#'   \code{"short"} generates shortened unique labels.
#' @return Labelled vector of node names from the dagitty object
#' @export
get_labels <- function(dag, label_type = "name"){
  if( !label_type %in% c("name", "initials", "short") ){
    stop("Invalid 'label_type'. Use 'name', 'initials', or 'short'.",
         call. = FALSE)
  }

  node_names <- as.character(names(dag))
  if( length(node_names) == 0L ){
    return(stats::setNames(character(), character()))
  }
  if( label_type == "name" ){
    return(stats::setNames(node_names, node_names))
  }

  words <- strsplit(node_names, "[_[:space:]]+", perl = TRUE)
  labels <- if( label_type == "initials" ){
    vapply(words, function(x) paste0(substr(x, 1L, 1L), collapse = ""),
           character(1L))
  }else{
    as.character(unlist(get_label_helper(words), use.names = FALSE))
  }
  labels <- make.unique(labels, sep = "_")
  stats::setNames(labels, node_names)
}


#' Compare dagitty object edges
#'
#' @importFrom data.table setDT
#' @param dag1 A dagitty object.
#' @param dag2 A dagitty object.
#' @param compare_all Return mutually exclusive edges in both dag1 and dag2 (compare_all = TRUE) e.g. dag1 edges not in dag2 AND dag2 edges not in dag1, or only dag1 mutually exclusive edges (compare_all = FALSE), e.g. dag1 edges not in dag2. Defaults to compare_all = TRUE.
#' @param include_roles Determines whether to include node roles in addition to node names and edge direction. Defaults to include_roles = FALSE, e.g. only considers edges where descendant/ancestor node names and edge direction are different.
#' @return Data table of edges.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y", type = "saturated"
#' )
#' dag2 <- build_graph(
#'   variables = c("Z1", "Z3"), treatments = "X", outcomes = "Y", type = "saturated"
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
#' @param dag1 A dagitty object.
#' @param dag2 A dagitty object.
#' @param compare_all Compare and return mutually all exclusive roles  (compare_all = TRUE) e.g. dag1 roles not in dag2 AND dag2 roles not in dag1, or set compare_all = FALSE to return only dag1 mutually exclusive roles, e.g. dag1 roles not in dag2. Defaults to compare_all = TRUE.
#' @param nodes_in_common Determines whether to include only the nodes found in both DAGs. Defaults to nodes_in_common = FALSE (returns all nodes).
#' @return A list of nodes.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y", type = "saturated"
#' )
#' dag2 <- build_graph(
#'   variables = c("Z1", "Z3"), treatments = "X", outcomes = "Y", type = "saturated"
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

  diff_list <- list()

  diff_list <- lapply(1:length(roles), function(x){

    diff_list[[x]] <- roles[[x]][ !roles[[x]] %in% roles_y[[x]] ]

    if( compare_all == TRUE){
      diff_list[[x]] <- c( diff_list[[x]], roles_y[[x]][ !roles_y[[x]] %in% roles[[x]] ] )
    }
    if( nodes_in_common == TRUE){
      dag_names <- names(dag1)
      dag2_names <- names(dag2)
      diff_list[[x]] <- diff_list[[x]][ diff_list[[x]] %in% dag_names &
                                          diff_list[[x]] %in% dag2_names ]
    }
    diff_list[[x]] <- diff_list[[x]][ complete.cases(diff_list[[x]]) ] # drop NA placeholders from empty role buckets
    diff_list[[x]]
  })

  names(diff_list) <- names(roles)

  diff_list <- Filter(function(x) length(x) > 0, diff_list)

  return(diff_list)
}


#' Extract ancestor node edges
#'
#' get_ancestor_edges() returns a list of edges grouped by node names for an input dagitty object. Symmetric edges ("<->", "--") are listed once, under their stored left-hand endpoint.
#'
#' @importFrom dagitty edges
#' @param dag A dagitty object.
#' @returns Named list of edges.
#' @export
get_ancestor_edges <- function(dag){
  edges <- data.table::as.data.table(.dag_edges(dag, "get_ancestor_edges()"))

  ## group by nodes (extract column as a vector so this holds for a plain data.frame as
  ## well as a data.table; nrow(unique(DT[,"v"])) silently relied on data.table input)
  unique_ancestors <- unique( edges[["v"]] )
  if( length(unique_ancestors) == 0L ){
    return(list())
  }

  # edges_to_assess grouped by each unique node in a list
  edges_list <- suppressWarnings( lapply(seq_along(unique_ancestors), function(x){

    edges[ unlist(edges[,"v"]) %in% unlist(unique_ancestors[x]), ]


  }) )

  names(edges_list) <- unlist(unique_ancestors)

  return(edges_list)
}
