# ==============================================================================
# graph_funcs.r


# Public graph construction, editing and edge assessment functions
# ==============================================================================


## Package constants ##
VALID_GRAPH_TYPES <- c("full", "chain", "ordered", "saturated") # Graph types accepted by node editing


VALID_BUILD_TYPES <- c("full", "chain", "first", "last", "ordered", "saturated") # Graph types accepted by build_graph()


VALID_NODE_ROLES <- c( # Roles accepted by node editing
  "outcome", "treatment", "confounder", "mediator",
  "mediator_outcome_confounder", "instrument", "competing_cause",
  "collider", "latent", "observed"
)


ALL_NODE_ROLES <- c(VALID_NODE_ROLES, "undetermined") # Inspection roles in output order


## Graph construction ##
#' Build a graph from causal roles
#'
#' `build_graph()` creates a dagitty DAG from node names grouped by causal role.
#' A node listed in more than one role input receives the arrows generated for
#' each supplied role.
#'
#' @importFrom dagitty is.dagitty children coordinates
#' @param variables Confounder names, or a dagitty object to rebuild from its
#'   inferred roles. Use `confounders` as a more explicit name for a new graph.
#'   Rebuilding regenerates the edges from those roles. Use `keep_edges()`,
#'   `add_nodes()`, `connect_nodes()`, or `join_graphs()` for changes to an
#'   existing graph.
#' @param treatments Treatment node names. Required for a new graph.
#' @param outcomes Outcome node names. Required for a new graph.
#' @param mediators Mediator node names.
#' @param latents Latent node names. To map one to specific children,
#'   use a form such as `list(list("U", c("M1", "M2")))`.
#' @param instrumental_variables Instrumental variable names, for example
#'   `"IV"`.
#' @param mediator_outcome_confounders Names of common causes of mediators and
#'   outcomes, for example `M <- X -> Y`.
#' @param competing_causes Names of other direct causes of an outcome.
#' @param colliders Collider node names.
#' @param type One of `"full"`, `"chain"`, `"first"`, `"last"`, `"ordered"`, or
#'   `"saturated"`.
#' @param confounders Optional alias for `variables` when supplying confounder
#'   names.
#' @param latent_variables Deprecated as of 0.0.2. Named alias for `latents`.
#'
#' @details
#' For role inputs other than instrumental variables, a vector gives successive
#' temporal positions. Nodes in the same list element share a temporal
#' position. For example, `list("M1", c("M2", "M3"), "M4")` places `M2` and
#' `M3` together, after `M1` and before `M4`.
#'
#' `type = "full"` uses bidirected edges following the ESC-DAGs Mapping Stage in
#' Ferguson et al. (2020). `"chain"` connects adjacent temporal positions within each role except instrumental variables
#' and adds no arrows between roles apart from explicit
#' latent targets and instrument edges. `"last"` adds arrows between roles from
#' the last position of the sending role to the first position of the receiving
#' role. `"first"` sends from the first position instead. `"ordered"` connects
#' every earlier position to every later position within each role except
#' instrumental variables. `"saturated"` also adds every directed arrow allowed by the package's
#' order of causal roles. Use `incoming()` and `outgoing()` to choose where
#' generated arrows enter or leave a role.
#'
#' @returns A dagitty object.
#' @examples
#' dag <- build_graph(
#'   confounders = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = list("M1", c("M2", "M3"), "M4"), type = "last"
#' )
#'
#' routed <- build_graph(
#'   confounders = "X", treatments = "T", outcomes = "Y",
#'   mediators = c(incoming("M1"), "M2", outgoing("M3")),
#'   type = "saturated"
#' )
#'
#' @export
build_graph <- function(variables = NA,
                       treatments = NA,
                       outcomes = NA,
                       mediators = NA,
                       latents = NA,
                       instrumental_variables = NA,
                       mediator_outcome_confounders = NA,
                       competing_causes = NA,
                       colliders = NA,
                       type = "full",
                       confounders = NULL,
                       latent_variables = NULL
                       ){

  if( !is.null(confounders) ){
    if( !missing(variables) && !(length(variables) == 1L && is.atomic(variables) && is.na(variables)) ){
      stop("build_graph(): supply either 'variables' or 'confounders', ",
           "not both.", call. = FALSE)
    }
    variables <- confounders
  }

  if( !missing(latent_variables) ){
    if( !missing(latents) && !(length(latents) == 1L &&
                               is.atomic(latents) && is.na(latents)) ){
      stop("build_graph(): supply either 'latents' or its deprecated alias ",
           "'latent_variables', not both.", call. = FALSE)
    }
    warning("build_graph(): 'latent_variables' is deprecated as of 0.0.2. ",
            "Use 'latents' instead.", call. = FALSE)
    latents <- latent_variables
  }

  if( length(type) != 1L || !type %in% VALID_BUILD_TYPES ){
    stop("Invalid type: \"", paste(type, collapse = ", "), "\". Must be one of: ",
         paste(VALID_BUILD_TYPES, collapse = ", "), ".", call. = FALSE)
  }

  if( dagitty::is.dagitty(variables) ){ # check for an existing dag (inputted as confounders)

    existing_dag <- variables
    existing_graph_type <- dag_graph_type(existing_dag)

    if( existing_graph_type != "dag" ){
      warning("build_graph() constructs a new DAG from inferred node roles; ",
              "the input's mixed graph type and edge semantics are not retained. ",
              "Use keep_edges(), add_nodes(), connect_nodes(), or join_graphs() ",
              "to edit a pdag or mag while retaining its type.", call. = FALSE)
    }

    if( type == "full" ){
      warning("build_graph(type = \"full\") rebuilds an existing graph from ",
              "inferred node roles and does not preserve its directed edges. ",
              "Use connect_nodes(dag, type = \"full\") to retain existing edges.",
              call. = FALSE)
    }

    node_roles <- classify_node_roles(existing_dag)

    confounders <- variables <- node_roles$confounder

    treatments <- node_roles$treatment

    outcomes <- node_roles$outcome

    m_o_confounder_vec <- mediator_outcome_confounders <- node_roles$mediator_outcome_confounder

    mediator_vec <- mediators <- node_roles$mediator

    instrumental_variables <- node_roles$instrument

    collider_vec <- colliders <- node_roles$collider

    competing_cause_vec <- competing_causes <- node_roles$competing_cause

    observed <- unique(c(node_roles$observed, node_roles$undetermined))
    observed <- observed[!is.na(observed)]

    latent_vec <- node_roles$latent

    if( all(complete.cases(unlist(latent_vec))) ){

      latent_descendants_list <- lapply(1:length(latent_vec), function(x){
        latent_descendants_list <- dagitty::children(existing_dag, latent_vec[x])

        return(invisible(latent_descendants_list))
      })

      latent_variables_list <- lapply(1:length(latent_vec), function(x){
        latent_variables_list <- latent_vec[x]

        return(invisible(latent_variables_list))
      })

      latents <- lapply(1:length(latent_vec), function(x){
        if( length(unlist(latent_descendants_list[x])) == 0 ){
          latent_mapping <- latent_variables_list[[x]]

          return(latent_mapping)
        }else{
          latent_mapping <- c( latent_variables_list[x], list(as.list(unlist(latent_descendants_list[x]))) )

          return(invisible(latent_mapping))
        }
      })

    }

  }else if( !absent_role_input(treatments) && !absent_role_input(outcomes) &&
            all( complete.cases( unlist(treatments) ) ) && all( complete.cases( unlist(outcomes) ) ) ){ # no existing dag, use other inputs

    confounders <- variables

    existing_dag <- NA

    mediator_vec <- mediators

    latent_vec <- latents

    m_o_confounder_vec <- mediator_outcome_confounders

    competing_cause_vec <- competing_causes

    collider_vec <- colliders

    observed <- NA

    if( all( complete.cases( unlist(confounders) ) ) &
        any( unlist(confounders, use.names = FALSE) %in%
             unlist(m_o_confounder_vec, use.names = FALSE) ) ){ # if any mediator-outcome confounders are also inputted as confounders, execution is stopped

      stop("Node names inputted in the 'variables' parameter detected in 'mediator_outcome_confounder'. These inputs should be mutually exclusive. Please adjust inputs and try again.")
    }

  }else{

    stop("The 'treatments' and 'outcomes' inputs should be used if a DAG is not provided in the 'variables' input. Please adjust inputs and try again.")
  }

  role_specs <- role_specs(
    confounders = confounders,
    treatments = treatments,
    outcomes = outcomes,
    mediators = mediators,
    mediator_outcome_confounders = mediator_outcome_confounders,
    competing_causes = competing_causes,
    colliders = colliders,
    latent_variables = latents,
    observed = observed,
    caller = "build_graph()"
  )

  edges <- plan_graph_edges(
    role_specs,
    type = type,
    instrumental_variables = instrumental_variables,
    existing_dag = existing_dag
  )

  latent_vec <- role_specs$latent$nodes
  treatment_vec <- role_specs$treatment$nodes
  outcome_vec <- role_specs$outcome$nodes
  observed_node_names <- unique(c(
    unlist(lapply(role_specs[names(role_specs) != "latent"], `[[`, "nodes"),
           use.names = FALSE),
    unlist(instrumental_variables, recursive = TRUE, use.names = FALSE)
  ))
  observed_node_names <- observed_node_names[
    !is.na(observed_node_names) & nzchar(observed_node_names)
  ]
  node_names <- setdiff(observed_node_names,
                        c(treatment_vec, outcome_vec, latent_vec))

  dag <- construct_graph(edges, node_names, treatment_vec, outcome_vec,
                         latent_vec)
  if( isFALSE(dag_is_acyclic(dag)) ){
    warning("build_graph(): the supplied roles and temporal positions create a directed cycle.",
            call. = FALSE)
  }

  return(dag)
}


#' Add nodes to a dagitty object
#'
#' Adds nodes, generates their edges, and places them using the existing graph
#' coordinates. The graph name, bounding box, labels, notes, and stored curves
#' are kept.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param nodes Node names to add. When a `role` other than `"instrument"` is supplied,
#'   use a vector or list to set temporal positions and `incoming()` or
#'   `outgoing()`, or the names `incoming` and `outgoing`, to choose where
#'   arrows enter or leave the role.
#' @param parents Optional existing nodes to use as direct parents.
#' @param children Optional existing nodes to use as direct children.
#' @param role Optional causal role for the new nodes. Only `"treatment"`,
#'   `"outcome"`, and `"latent"` can be combined with explicit `parents` or
#'   `children`. When `parents` or `children` is supplied, those arguments
#'   determine the connections and `role` sets only the node status.
#' @param new_nodes,ancestors,descendants,node_role Deprecated as of 0.0.2.
#'   Named aliases for `nodes`, `parents`, `children`, and `role`.
#' @param type One of `"full"`, `"chain"`, `"ordered"`, or `"saturated"`.
#'   Defaults to `"saturated"`. `"first"` and `"last"` are deprecated as of
#'   0.0.1. Use `position` instead. With explicit `parents` or `children`,
#'   `"full"` creates bidirected edges. The other types create directed edges.
#' @param position Optional placement in the temporal order for a role other than
#'   `"instrument"`, expressed with [first()], [last()], [before()], and [after()], or
#'   with `"first"`, `"last"`, `c(before = "M2")`, and
#'   `c(after = "M2")`. Without a role, placement in the whole graph's order
#'   sets the direction of connections to existing nodes. Explicit `parents`
#'   or `children` take precedence over `position`. With `type = "full"`,
#'   positions do not change the direction of bidirected edges. Position is
#'   not supported for the instrument role.
#' @param print_edges Print new edges in console, defaults to FALSE.
#' @param renew_coords Regenerate the coordinates of every node after adding.
#'   Defaults to `FALSE`, which places only the new nodes and keeps existing
#'   coordinates. If a new node would sit on an existing node, every coordinate
#'   is regenerated with a warning. With `TRUE`, stored curves are carried to
#'   the new coordinates.
#' @returns A dagitty object.
#' @examples
#' dag <- build_graph(
#'   confounders = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", type = "saturated"
#' )
#' dag <- add_nodes(dag, nodes = "W", role = "confounder")
#' dag <- add_nodes(dag, nodes = "U", role = "latent",
#'                  children = c("X1", "X2"))
#'
#' @export
add_nodes <- function(dag,
                      nodes,
                      parents = NULL,
                      children = NULL,
                      role = NULL,
                      type = "saturated",
                      print_edges = FALSE,
                      position = NULL,
                      renew_coords = FALSE,
                      new_nodes = NULL,
                      ancestors = NULL,
                      descendants = NULL,
                      node_role = NULL
                      ){
  .datatable.aware <- TRUE

  role_alias_used <- FALSE
  nodes_missing <- missing(nodes)
  if( !missing(new_nodes) ){
    if( !nodes_missing ){
      stop("add_nodes(): supply either 'nodes' or its deprecated alias ",
           "'new_nodes', not both.", call. = FALSE)
    }
    warning("add_nodes(): 'new_nodes' is deprecated as of 0.0.2. ",
            "Use 'nodes' instead.", call. = FALSE)
    nodes <- new_nodes
    nodes_missing <- FALSE
  }
  if( nodes_missing ){
    stop("add_nodes(): supply at least one node name in 'nodes'.",
         call. = FALSE)
  }
  if( !missing(ancestors) ){
    if( !missing(parents) ){
      stop("add_nodes(): supply either 'parents' or its deprecated alias ",
           "'ancestors', not both.", call. = FALSE)
    }
    warning("add_nodes(): 'ancestors' is deprecated as of 0.0.2. ",
            "Use 'parents' instead.", call. = FALSE)
    parents <- ancestors
  }
  if( !missing(descendants) ){
    if( !missing(children) ){
      stop("add_nodes(): supply either 'children' or its deprecated alias ",
           "'descendants', not both.", call. = FALSE)
    }
    warning("add_nodes(): 'descendants' is deprecated as of 0.0.2. ",
            "Use 'children' instead.", call. = FALSE)
    children <- descendants
  }
  if( !missing(node_role) ){
    if( !missing(role) ){
      stop("add_nodes(): supply either 'role' or its deprecated alias ",
           "'node_role', not both.", call. = FALSE)
    }
    warning("add_nodes(): 'node_role' is deprecated as of 0.0.2. ",
            "Use 'role' instead.", call. = FALSE)
    role <- node_role
    role_alias_used <- TRUE
  }

  if( length(type) != 1L || !type %in% c(VALID_GRAPH_TYPES, "first", "last") ){
    stop("Invalid type: \"", paste(type, collapse = ", "), "\". Must be one of: ",
         paste(VALID_GRAPH_TYPES, collapse = ", "), ".", call. = FALSE)
  }
  if( !is.logical(renew_coords) || length(renew_coords) != 1L ||
      is.na(renew_coords) ){
    stop("add_nodes(): renew_coords must be TRUE or FALSE.", call. = FALSE)
  }
  if( type %in% c("first", "last") ){
    warning("add_nodes(): 'type = \"", type, "\"' is deprecated as of 0.0.1. ",
            "Use 'position = ", type, "()' instead.", call. = FALSE)
    if( !is.null(position) ){
      stop("add_nodes(): supply either 'position' or the deprecated ",
           "'type = \"first\"/\"last\"', not both.", call. = FALSE)
    }
  }
  if( length(role) > 1L ){
    stop("add_nodes(): 'role' must be one role name.", call. = FALSE)
  }
  if( length(role) == 1L && !role %in% VALID_NODE_ROLES ){
    role_name <- if( role_alias_used ) { "node_role" } else { "role" }
    stop("Invalid ", role_name, ": \"", role, "\". Must be one of: ",
         paste(VALID_NODE_ROLES, collapse = ", "), ".", call. = FALSE)
  }
  explicit_edges <- !is.null(parents) || !is.null(children)
  if( length(role) > 0L && explicit_edges ){
    if( !role %in% c("treatment", "outcome", "latent") ){
      stop("add_nodes(): explicit parents/children can only be combined with ",
           "role = \"treatment\", \"outcome\", or \"latent\".", call. = FALSE)
    }
    warning("add_nodes(): parents/children input overrides role-based edges.",
            call. = FALSE)
  }

  supplied_spec <- NULL
  if( length(role) == 1L && !identical(role, "instrument") ){
    supplied_spec <- normalise_role_input(
      nodes, role, caller = "add_nodes() 'nodes'"
    )
    supplied_node_names <- supplied_spec$nodes
  }else{
    supplied_node_names <- as.character(unlist(nodes, recursive = TRUE,
                                               use.names = FALSE))
  }
  if( length(supplied_node_names) == 0L ){
    stop("add_nodes(): supply at least one node name in 'nodes'.",
         call. = FALSE)
  }
  if( anyNA(supplied_node_names) || any(!nzchar(supplied_node_names)) ){
    stop("add_nodes(): 'nodes' must contain non-missing, non-empty names.",
         call. = FALSE)
  }

  ## get initial dag edges ##
  edges <- data.table::as.data.table(dag_edges(dag, "add_nodes()"))

  names(edges) <- c("v", "e", "w") # change column names

  dag_node_names <- names(dag) # extract dag node names

  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  if( !explicit_edges && !is.null(supplied_spec) &&
      nrow(supplied_spec$mapping_edges) > 0L ){
    unknown_targets <- setdiff(unique(supplied_spec$mapping_edges$w),
                               dag_node_names)
    if( length(unknown_targets) > 0L ){
      stop("add_nodes(): latent mapping targets are not in the graph: ",
           paste(unknown_targets, collapse = ", "), ".", call. = FALSE)
    }
  }

  already_present <- unique(supplied_node_names[
    supplied_node_names %in% dag_node_names
  ])
  if( length(already_present) == length(unique(supplied_node_names)) ){
    stop("add_nodes(): all supplied nodes are already in the graph: ",
         paste(already_present, collapse = ", "),
         ". Use connect_nodes() to connect existing nodes.", call. = FALSE)
  }
  if( length(already_present) > 0L ){
    warning("add_nodes(): skipping names already in the graph: ",
            paste(already_present, collapse = ", "), ".", call. = FALSE)
  }

  new_node_names <- unique(supplied_node_names[
    !supplied_node_names %in% dag_node_names
  ])

  if( explicit_edges ){
    new_edges <- connect_new_nodes(dag = dag,
                                   new_nodes = new_node_names,
                                   ancestors = parents,
                                   descendants = children,
                                   node_roles = NULL,
                                   type = type)
  }else if( !is.null(supplied_spec) ){
    supplied_spec <- subset_role_spec(
      supplied_spec, new_node_names, preserve_routing = TRUE
    )
    planned_specs <- specs_after_role_operation(
      dag, supplied_spec, role, type, position = position
    )
    new_edges <- role_operation_edges(
      dag, supplied_spec, role, type, position = position,
      specs = planned_specs
    )
  }else{
    new_edges <- connect_new_nodes(dag = dag,
                                   new_nodes = new_node_names,
                                   ancestors = parents,
                                   descendants = children,
                                   node_roles = role,
                                   type = type,
                                   position = position)
  }

  ## pre-process before merging ##
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs
  new_edges <- new_edges[ new_edges$v != new_edges$w, ] # remove new_edges with identical ancestor and descendant node names

  ## merge new and existing dag edges ##
  edges <- merge(edges, new_edges,
                 by = c("v", "e", "w"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges
  edges <- edges[edges$v != edges$w, ] # remove new_edges with identical ancestor and descendant node names
  updated_treatments <- treatments(dag)
  updated_outcomes <- get_outcome_names(dag)
  updated_latents <- unobserved(dag)
  new_status <- NULL
  if( identical(role, "treatment") ){
    updated_treatments <- union(updated_treatments, new_node_names)
    new_status <- "exposure"
  }else if( identical(role, "outcome") ){
    updated_outcomes <- union(updated_outcomes, new_node_names)
    new_status <- "outcome"
  }else if( identical(role, "latent") ){
    updated_latents <- union(updated_latents, new_node_names)
    new_status <- "latent"
  }

  ## append to the native text; rebuild only when the append cannot be verified ##
  appended <- native_append_graph(
    dag,
    new_nodes = new_node_names,
    new_edges = new_edges,
    statuses = stats::setNames(rep(list(new_status), length(new_node_names)),
                               new_node_names),
    caller = "add_nodes()"
  )
  if( is.null(appended) ){
    appended <- rebuild_dag(
      dag, edges,
      treatments = updated_treatments,
      outcomes = updated_outcomes,
      latent_vec = updated_latents,
      extra_nodes = new_node_names,
      .curve_caller = "add_nodes"
    )
  }
  dag <- appended
  if( isFALSE(dag_is_acyclic(dag)) ){
    warning("add_nodes(): the requested edges create a directed cycle.",
            call. = FALSE)
  }

  if( print_edges == TRUE){
    new_edges_list <- print_edges_helper(new_edges)
    cat( paste("c(", paste( unlist(new_edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    message("\nPrinted new edges - copy and paste to use as a vector object.\n")
  }

  new_coordinates <- NULL
  if( !isTRUE(renew_coords) ){
    new_coordinates <- renew_coords(dag = dag,
                                    new_node_names = new_node_names,
                                    coordinates = coordinates)
    positioned <- names(coordinates$x)[is.finite(coordinates$x) & is.finite(coordinates$y)]
    new_coordinates$x[positioned] <- coordinates$x[positioned]
    new_coordinates$y[positioned] <- coordinates$y[positioned]
    if( copy_coincident_copies(new_coordinates, new_node_names) ){
      warning("add_nodes(): the new nodes could not be placed apart from existing nodes; ",
              "coordinates were regenerated for every node.", call. = FALSE)
      new_coordinates <- NULL
    }
  }
  if( is.null(new_coordinates) ){
    new_coordinates <- renew_coords(dag)
  }
  dag <- set_coords_preserving_curves(dag, new_coordinates, "add_nodes")

  return(dag)
}


#' Connect existing nodes in a dagitty object
#'
#' Adds edges to existing nodes. Graph metadata and coordinates are preserved.
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param nodes Existing nodes to include in new connections. Each new edge
#'   must involve at least one selected node. Defaults to all graph nodes.
#' @param type One of `"full"`, `"chain"`, `"ordered"`, or `"saturated"`.
#'   Defaults to `"full"`. `"first"` and `"last"` are deprecated as of 0.0.1
#'   and are treated as `"ordered"`.
#' @param print_edges Print new edges in console, defaults to FALSE.
#' @param ... Deprecated arguments. `role`, `node_role`, and `position` are
#'   deprecated as of 0.0.2 and are ignored with a warning. Any other argument causes an error.
#' @returns A dagitty object.
#' @examples
#' dag <- build_graph(
#'   confounders = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   mediators = "M", type = "saturated"
#' )
#' dag <- connect_nodes(dag, nodes = "X1")
#'
#' @export
connect_nodes <- function(dag,
                           nodes = NULL,
                           type = "full",
                           print_edges = FALSE,
                           ...
                           ){
  .datatable.aware <- TRUE

  extra <- list(...)
  extra_names <- names(extra)
  if( is.null(extra_names) ){
    extra_names <- rep("", length(extra))
  }
  for( argument in c("role", "node_role", "position") ){
    if( argument %in% extra_names ){
      warning("connect_nodes(): '", argument,
              "' is deprecated as of 0.0.2 and was ignored.", call. = FALSE)
    }
  }
  unused <- which(!extra_names %in% c("role", "node_role", "position"))
  if( length(unused) > 0L ){
    labels <- ifelse(nzchar(extra_names[unused]), extra_names[unused],
                     paste("argument", unused))
    stop("connect_nodes(): unused argument(s): ",
         paste(labels, collapse = ", "), ".", call. = FALSE)
  }

  if( length(type) != 1L || !type %in% c(VALID_GRAPH_TYPES, "first", "last") ){
    stop("Invalid type: \"", paste(type, collapse = ", "), "\". Must be one of: ",
         paste(VALID_GRAPH_TYPES, collapse = ", "), ".", call. = FALSE)
  }
  if( type %in% c("first", "last") ){
    warning("connect_nodes(): 'type = \"", type,
            "\"' is deprecated as of 0.0.1. ",
            "Use 'type = \"ordered\"' instead.", call. = FALSE)
  }
  ## get initial dag edges ##
  edges <- data.table::as.data.table(dag_edges(dag, "connect_nodes()"))
  dag_node_names <- names(dag) # extract dag node names

  if( length(nodes) == 0L ){
    nodes <- dag_node_names
  }
  selected_node_names <- as.character(unlist(nodes, recursive = TRUE,
                                             use.names = FALSE))

  if( length(selected_node_names) == 0L ){
    stop("connect_nodes(): supply at least one node name in 'nodes'.",
         call. = FALSE)
  }

  missing_nodes <- unique(selected_node_names[
    !selected_node_names %in% dag_node_names
  ])
  if( length(missing_nodes) > 0L ){
    stop("connect_nodes(): nodes not in the graph: ",
         paste(missing_nodes, collapse = ", "),
         ". Use add_nodes() to introduce new nodes.", call. = FALSE)
  }

  new_edges <- connect_nodes_helper(
    dag, selected_node_names, dag_node_names, type
  )

  ## pre-process before merging ##
  new_edges <- new_edges[ complete.cases(new_edges), ] # remove NAs


  ## merge new and existing dag edges ##
  edges <- merge(edges, new_edges,
                 by = c("v", "e", "w"),
                 all = TRUE) # combine both dag edges

  edges <- unique(edges) # remove duplicate new_edges

  edges <- edges[edges$v != edges$w, ] # remove new_edges with identical ancestor and descendant node names

  updated_treatments <- treatments(dag)
  updated_outcomes <- get_outcome_names(dag)
  updated_latents <- unobserved(dag)

  ## append to the native text; rebuild only when the append cannot be verified ##
  appended <- native_append_graph(
    dag,
    new_nodes = character(),
    new_edges = new_edges[ new_edges$v != new_edges$w, ],
    caller = "connect_nodes()"
  )
  if( is.null(appended) ){
    appended <- rebuild_dag(
      dag, edges,
      treatments = updated_treatments,
      outcomes = updated_outcomes,
      latent_vec = updated_latents,
      .curve_caller = "connect_nodes"
    )
  }
  dag <- appended

  if( isFALSE(dag_is_acyclic(dag)) ){
    warning("The connected graph contains a directed cycle.", call. = FALSE)
  }

  if( print_edges == TRUE){

    new_edges_list <- print_edges_helper(new_edges)

    cat( paste("c(", paste( unlist(new_edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    message("\nPrinted new edges to assess - copy and paste in a .R file to use as a vector object.\n")
  }

  return(dag)
}


#' Merge two dagitty objects
#'
#' `join_graphs()` adds the nodes and edges of the second dagitty object to the
#' first. By default, new nodes are placed using the first graph's existing
#' coordinates.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty dagitty edges coordinates
#' @param dag First dagitty object.
#' @param new_dag Second dagitty object, added to the first.
#' @param renew_coords Regenerate the coordinates of every node after joining.
#'   Defaults to `FALSE`, which places only the new nodes and keeps existing
#'   coordinates. If a new node would sit on an existing node, every coordinate
#'   is regenerated with a warning. With `TRUE`, stored curves are carried to
#'   the new coordinates.
#'
#' @details
#' Stored edge control points are kept for retained edges. When both graphs
#' contain the same edge, including its edge type, its path comes from `dag`. A straight edge in `dag` stays straight. Curved edges added from `new_dag` move
#' with their endpoints, preserving their bend relative to the edge length.
#' They are not rerouted around nodes or other edges in the merged graph.
#'
#' @returns A dagitty object
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y", type = "saturated"
#' )
#' dag2 <- build_graph(
#'   variables = c("X1", "X3"), treatments = "T", outcomes = "Y", type = "saturated"
#' )
#' new_graph <- join_graphs(dag, dag2)
#'
#' @export
join_graphs <- function(dag,
                        new_dag,
                        renew_coords = FALSE
                        ){
  .datatable.aware <- TRUE

  if( !dagitty::is.dagitty(dag) ){
    stop("join_graphs(): 'dag' must be a dagitty object.", call. = FALSE)
  }
  if( missing(new_dag) || !dagitty::is.dagitty(new_dag) ){
    stop("join_graphs(): 'new_dag' must be a dagitty object.", call. = FALSE)
  }
  if( !is.logical(renew_coords) || length(renew_coords) != 1L ||
      is.na(renew_coords) ){
    stop("join_graphs(): renew_coords must be TRUE or FALSE.", call. = FALSE)
  }

  edges <- data.table::as.data.table(dag_edges(dag, "join_graphs()"))
  first_edge_keys <- edge_control_keys(edges$v, edges$e, edges$w)
  donor_controls <- edge_controls(new_dag, "join_graphs")
  donor_controls <- donor_controls[
    !edge_control_keys(donor_controls$v, donor_controls$e, donor_controls$w) %in% first_edge_keys,
    , drop = FALSE
  ]
  node_names <- names(dag) # extract dag node names

  treatments <- treatments(dag) # get treatments
  outcomes <- get_outcome_names(dag) # get outcomes
  latent_vec <- unlist(unobserved(dag)) # get dag latent variables
  coordinates <- dagitty::coordinates(dag) # extract dag coordinates

  new_edges <- data.table::as.data.table(dag_edges(new_dag, "join_graphs()"))
  new_node_names <- names(new_dag) # extract new dag node names

  existing_node_names <-  new_node_names[ new_node_names %in% node_names ] # saves duplicate node names
  shared_node_names <- existing_node_names # shared nodes, before coordinate filtering
  existing_node_coords_y <- coordinates$y[ names(coordinates$y)  %in% existing_node_names] # saves duplicate node names
  existing_node_names <- names( existing_node_coords_y[ order(existing_node_coords_y) ] ) # existing dag node names in ascending y-coords order

  new_node_names <- new_node_names[ !new_node_names %in% node_names ] # remove duplicate node names


  new_dag_treatments <- treatments(new_dag)
  treatments <- c(treatments,
                  new_dag_treatments[ !new_dag_treatments %in% node_names] ) # get treatments

  new_dag_outcomes <- get_outcome_names(new_dag) # get outcomes
  outcomes <- c(outcomes,
                new_dag_outcomes[ !new_dag_outcomes %in% node_names] ) # get treatments

  new_dag_latent_vec <- unlist( unobserved(new_dag) ) # get dag latent variables
  latent_vec <- c(latent_vec,
                  new_dag_latent_vec[ !new_dag_latent_vec %in% node_names] ) # get treatments
  adjusted_vec <- union(dagitty::adjustedNodes(dag),
                        dagitty::adjustedNodes(new_dag))
  selected_vec <- union(selected_nodes(dag), selected_nodes(new_dag))

  ## warn when the two graphs disagree on a shared node's declarations ##
  if( length(shared_node_names) > 0 ){
    conflicted <- Filter(function(nm){
      matches <- ( nm %in% treatments ) != ( nm %in% new_dag_treatments ) |
      ( nm %in% outcomes ) != ( nm %in% new_dag_outcomes ) |
      ( nm %in% latent_vec ) != ( nm %in% new_dag_latent_vec )

      return(matches)
    }, shared_node_names)
    if( length(conflicted) > 0 ){
      warning("join_graphs(): the two graphs disagree on the exposure/outcome/latent status of shared node(s): ",
              paste(conflicted, collapse = ", "),
              ". The first graph's declarations were kept. Please check that this is intended.",
              call. = FALSE)
    }
  }


  coordinates_new_dag <- dagitty::coordinates(new_dag) # extract new dag coordinates
  new_node_coords_y <- coordinates_new_dag$y[ names(coordinates_new_dag$y) %in% new_node_names] # saves duplicate node names
  if( length(new_node_coords_y) > 0 ){
    new_node_names <- names( new_node_coords_y[ order(new_node_coords_y) ] ) # new dag node names in ascending y-coords order
  }

  edges <- unique(rbind(edges, new_edges)) # combine both dag edges, dropping duplicates

  ## rebuild dag ##
  dag <- rebuild_dag(dag, edges,
                     treatments = treatments,
                     outcomes = outcomes,
                     latent_vec = latent_vec,
                     extra_nodes = names(new_dag),
                     adjusted_vec = adjusted_vec,
                     selected_vec = selected_vec,
                     .curve_caller = "join_graphs")

  if( isFALSE(dag_is_acyclic(dag)) ){
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.",
            call. = FALSE)
  }

  relayout <- isTRUE(renew_coords)
  new_coordinates <- NULL
  if( !relayout && length(new_node_names) > 0L ){
    new_coordinates <- renew_coords(dag = dag,
                                     new_node_names = new_node_names,
                                     coordinates = coordinates)
    positioned <- names(coordinates$x)[is.finite(coordinates$x) & is.finite(coordinates$y)]
    new_coordinates$x[positioned] <- coordinates$x[positioned]
    new_coordinates$y[positioned] <- coordinates$y[positioned]
    if( copy_coincident_copies(new_coordinates, new_node_names) ){
      warning("join_graphs(): the new nodes could not be placed apart from existing nodes; ",
              "coordinates were regenerated for every node.", call. = FALSE)
      relayout <- TRUE
    }
  }
  if( relayout ){
    new_coordinates <- renew_coords(dag)
  }
  if( !is.null(new_coordinates) ){
    dag <- set_coords_preserving_curves(dag, new_coordinates, "join_graphs")
  }
  if( nrow(donor_controls) > 0L ){
    donor_controls <- transfer_edge_controls(donor_controls, coordinates_new_dag,
                                               dagitty::coordinates(dag), "join_graphs")
    dag <- set_edge_controls(dag, donor_controls, caller = "join_graphs")
  }

  return(dag)
}


#' Assess dagitty object edges
#'
#' `assess_edges()` lists graph edges for review and can guide you through an
#' assessment using causal criteria.
#'
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object containing the edges to assess.
#' @param edges_to_assess Which remaining edges to assess after applying
#'   `keep_edges` and `remove_edges`. Use `"all"` (the default) or
#'   `"bidirectional"` to assess only bidirected (`<->`) and undirected (`--`)
#'   edges. `"bidirected"`, `"bi-directional"`, and `"bi-directed"` are
#'   accepted aliases.
#' @param keep_edges Edges to retain without assessment. Supply one or more
#'   edges as character strings (for example, `"X -> Y"`), an edge table such
#'   as the output from [get_edges()], or a dagitty object whose edges will be
#'   used. These edges are excluded from assessment and retained in the
#'   output. Supplied edges not present in `dag` are ignored with a warning.
#' @param remove_edges Edges to discard without assessment. Accepts the same
#'   input forms and edge syntax as `keep_edges`. These edges are excluded
#'   from both the assessment and the output. Supplied edges not present in
#'   `dag` are ignored with a warning.
#' @param edges_to_keep Deprecated as of 0.0.2. Named alias for `keep_edges`. Supplying both is an error.
#' @param assess_causal_criteria Defaults to `FALSE`. With `TRUE`, the user is guided through a sequence
#'   that assesses each pair of connected nodes using causal criteria. Based on the Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs) from Ferguson et al. (2020).
#' @param causal_criteria The default `ESCDAGs` criteria consider temporality,
#'   face validity, and recourse to theory (Ferguson et al., 2020). To use other
#'   criteria, supply a data frame with the same columns as `ESCDAGs`.
#' @param save_answers Set to \code{TRUE} to return the decision log with the
#'   retained edges, or supply a previous log (v, e, w followed by one column per
#'   criterion) to replay it without prompting.
#' @returns With `assess_causal_criteria = FALSE`, a character vector of edges
#'   to assess. If edges were retained, removed, or selected by direction, a
#'   list containing `edges_to_assess` (an edge table) and `edges_to_keep` (a
#'   character vector) is returned instead.
#'
#'   With `assess_causal_criteria = TRUE`, an edge table of retained edges is
#'   returned. Setting `save_answers = TRUE` returns a list containing `edges`
#'   and `answers` (the decision log). Replaying a log returns the edge table.
#'   Declining the guided sequence returns a list containing `edges` (retained
#'   without assessment) and `edges_to_assess` (unassessed), both edge tables.
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2", "X3"), treatments = "T", outcomes = "Y", type = "last"
#' )
#'
#' ## a saturated graph built from the existing (sparser) dag
#' saturated_graph <- build_graph(type = "saturated", variables = dag)
#'
#' ## Option 1: keep the edges already present in the existing dag
#' edges <- assess_edges(saturated_graph, keep_edges = dag)
#'
#' \dontrun{
#' ## Option 2: interactively rule edges in/out via the guided causal-criteria sequence
#' edges_to_keep <- assess_edges(saturated_graph, keep_edges = dag,
#'                               assess_causal_criteria = TRUE)
#' }
#'
#' @export
assess_edges <- function(dag,
                         edges_to_assess = "all",
                         keep_edges = NULL,
                         assess_causal_criteria = FALSE,
                         save_answers = NULL,
                         causal_criteria = ESCDAGs,
                         remove_edges = NULL,
                         edges_to_keep = NULL
                         ){
  .datatable.aware <- TRUE

  keep_missing <- missing(keep_edges)
  if( !missing(edges_to_keep) ){
    if( !keep_missing ){
      stop("assess_edges(): supply either 'keep_edges' or its deprecated alias ",
           "'edges_to_keep', not both.", call. = FALSE)
    }
    warning("assess_edges(): 'edges_to_keep' is deprecated as of 0.0.2. ",
            "Use 'keep_edges' instead.", call. = FALSE)
    keep_edges <- edges_to_keep
    keep_missing <- FALSE
  }
  if( keep_missing ) {
    keep_edges <- NULL
  }

  allowed_modes <- c("all", "bidirectional", "bidirected",
                     "bi-directional", "bi-directed")
  if( !is.character(edges_to_assess) || length(edges_to_assess) != 1L ||
      is.na(edges_to_assess) || !edges_to_assess %in% allowed_modes ){
    stop("'edges_to_assess' must be one of: ",
         paste(allowed_modes, collapse = ", "), ".", call. = FALSE)
  }

  if( !dagitty::is.dagitty(dag) ){
    stop("assess_edges(): 'dag' must be a dagitty object.", call. = FALSE)
  }
  graph_edges <- data.table::as.data.table(dag_edges(dag, "assess_edges()"))
  empty_edges <- graph_edges[0L]

  normalise_selection <- function(value, argument){
    if( is.null(value) || length(value) == 0L ) {
      return(empty_edges)
    }
    if( length(value) == 1L && is.atomic(value) && is.na(value) ){
      return(empty_edges)
    }
    if( dagitty::is.dagitty(value) ){
      selected <- dag_edges(value, "assess_edges()")
    }else if( is.character(value) ){
      incomplete <- is.na(value)
      if( any(incomplete) ){
        warning("assess_edges(): entries in '", argument,
                "' containing missing values were ignored.", call. = FALSE)
        value <- value[!incomplete]
      }
      if( length(value) == 0L ){
        stop("assess_edges(): '", argument,
             "' contains no complete edges.", call. = FALSE)
      }
      selected <- parse_edge_strings(value, "assess_edges()")
    }else if( is.data.frame(value) || data.table::is.data.table(value) ){
      selected <- normalise_edge_frame(value, "assess_edges()")
      incomplete <- !complete.cases(selected)
      if( any(incomplete) ){
        warning("assess_edges(): entries in '", argument,
                "' containing missing values were ignored.", call. = FALSE)
        selected <- selected[!incomplete, , drop = FALSE]
      }
      if( nrow(selected) == 0L && nrow(value) > 0L ){
        stop("assess_edges(): '", argument,
             "' contains no complete edges.", call. = FALSE)
      }
    }else{
      stop("assess_edges(): '", argument,
           "' must be a character vector, an edge table, or a dagitty object.",
           call. = FALSE)
    }
    if( nrow(selected) == 0L ) {
      return(empty_edges)
    }
    selected <- selected[!duplicated(edge_key(selected)), , drop = FALSE]
    selected_edges <- data.table::as.data.table(selected)

    return(selected_edges)
  }

  requested_keep <- normalise_selection(keep_edges, "keep_edges")
  requested_remove <- normalise_selection(remove_edges, "remove_edges")
  keep_keys <- edge_key(requested_keep)
  remove_keys <- edge_key(requested_remove)
  overlap <- intersect(keep_keys, remove_keys)
  if( length(overlap) > 0L ){
    conflicted <- requested_keep[keep_keys %in% overlap]
    stop("assess_edges(): edges cannot be supplied in both 'keep_edges' and ",
         "'remove_edges': ",
         paste(conflicted$v, conflicted$e, conflicted$w, collapse = ", "),
         ".", call. = FALSE)
  }

  graph_keys <- edge_key(graph_edges)
  select_graph_edges <- function(requested, requested_keys, argument){
    missing_selection <- requested[!requested_keys %in% graph_keys]
    if( nrow(missing_selection) > 0L ){
      warning("assess_edges(): edges in '", argument,
              "' not present in 'dag' were ignored: ",
              paste(missing_selection$v, missing_selection$e,
                    missing_selection$w, collapse = ", "),
              ".", call. = FALSE)
    }
    selected_edges <- graph_edges[graph_keys %in% requested_keys]

    return(selected_edges)
  }
  kept_edges <- select_graph_edges(requested_keep, keep_keys, "keep_edges")
  removed_edges <- select_graph_edges(requested_remove, remove_keys,
                                      "remove_edges")
  classified_keys <- c(edge_key(kept_edges), edge_key(removed_edges))
  assessment_edges <- graph_edges[!graph_keys %in% classified_keys]

  bidirectional <- edges_to_assess != "all"
  if( bidirectional ){
    symmetric <- assessment_edges$e %in% c("<->", "--")
    kept_edges <- data.table::rbindlist(
      list(kept_edges, assessment_edges[!symmetric]), use.names = TRUE
    )
    assessment_edges <- assessment_edges[symmetric]
  }
  edges <- assessment_edges

  edge_strings <- function(value){
    if( nrow(value) == 0L ) {
      edge_text <- character()
      return(edge_text)
    }
    edge_text <- paste(value$v, value$e, value$w)

    return(edge_text)
  }

  if( nrow(edges) == 0L && nrow(graph_edges) == 0L ){
    stop("There are no edges to assess. Please check the supplied dagitty object and try again.")
  }

  if( assess_causal_criteria == TRUE ){

    if( nrow(edges) == 0L ){
      message("No edges remain to assess.")
      if( isTRUE(save_answers) ){
        answers <- data.frame(v = character(), e = character(), w = character())
        for( criterion in as.character(causal_criteria$name) ){
          answers[[criterion]] <- character()
        }
        assessment <- list(edges = kept_edges, answers = answers)
        return(assessment)
      }
      return(kept_edges)
    }

    edges <- causal_criteria_sequence(edges = edges,
                                      check_skip_sequence = FALSE,
                                      edges_to_keep = kept_edges,
                                      causal_criteria = causal_criteria,
                                      save_answers = save_answers)

    return(edges)
  }

  if( nrow(edges) != 0 ){
    edges_list <- print_edges_helper(edges) # group edges by unique node names
    cat( paste("c(", paste( unlist(edges_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )
  }else{
    message("No edges remain to assess.")
  }

  if( nrow(kept_edges) == 0L && nrow(removed_edges) == 0L &&
      !bidirectional ){

    edges <- edge_strings(edges)

    message("\nOutputted and printed edges to assess - copy and paste in a .R file to use as a vector object.")

    return(edges)
  }
  edges_list <- list( edges_to_assess = edges,
                      edges_to_keep = edge_strings(kept_edges) )
  if( nrow(edges) > 0L ){
    message("\nOutputting list containing 'edges_to_assess' and 'edges_to_keep'.
            \n\nPrinted edges to assess - copy and paste in a .R file to use as a vector object.")
  }

  return(edges_list)
}


#' Keep selected edges in a dagitty object
#'
#' `keep_edges()` removes every edge in a graph except the supplied edges.
#'
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom dagitty edges coordinates dagitty
#' @param edges Edges to keep. Supply one or more edges as character strings
#'   (for example, `"X -> Y"`), an edge table such as the output from
#'   [get_edges()], or a dagitty object whose edges will be used. Character
#'   input accepts `->`, `<-`, `<->`, and `--`, including chains and multiple
#'   edges separated by semicolons or line breaks. All other edges in `dag`
#'   are removed. Supplied edges not present in `dag` are ignored with a
#'   warning.
#' @param dag A dagitty object.
#' @param edges_to_keep Deprecated as of 0.0.2. Named alias for `edges`. Supplying both is an error.
#' @returns A dagitty object with all unlisted edges removed.
#' @seealso [add_edges()], [remove_edges()], [assess_edges()]
#'
#' @details
#' Existing graph metadata, node statuses, coordinates, labels, notes, R
#' attributes, and stored control points remain unchanged apart from those
#' attached to the edges that are not kept. A node left without edges stays
#' declared in the graph.
#'
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("X1", "X2", "X3"), treatments = "T", outcomes = "Y", type = "last"
#' )
#'
#' ## a saturated graph built from the existing (sparser) dag
#' saturated_graph <- build_graph(type = "saturated", variables = dag)
#'
#' ## keep only the edges present in the original dag
#' kept <- keep_edges(dag, saturated_graph)
#'
#' \dontrun{
#' ## interactive: curate edges via the guided causal-criteria sequence, then keep them
#' edges_to_keep <- assess_edges(saturated_graph, keep_edges = dag,
#'                               assess_causal_criteria = TRUE)
#' kept <- keep_edges(edges_to_keep, saturated_graph)
#' }
#'
#' @export
keep_edges <- function(edges,
                       dag,
                       edges_to_keep = NULL
                       ){
  edges_missing <- missing(edges)
  if( !missing(edges_to_keep) ){
    if( !edges_missing ){
      stop("keep_edges(): supply either 'edges' or its deprecated alias ",
           "'edges_to_keep', not both.", call. = FALSE)
    }
    warning("keep_edges(): 'edges_to_keep' is deprecated as of 0.0.2. ",
            "Use 'edges' instead.", call. = FALSE)
    edges <- edges_to_keep
    edges_missing <- FALSE
  }
  if( edges_missing ){
    stop("keep_edges(): supply at least one edge in 'edges'.", call. = FALSE)
  }
  if( !dagitty::is.dagitty(dag) ){
    stop("keep_edges(): 'dag' must be a dagitty object.", call. = FALSE)
  }

  existing_edges <- dag_edges(dag, "keep_edges()")

  if( dagitty::is.dagitty(edges)){
    edges <- dag_edges(edges, "keep_edges()")
  }else if( is.character(edges) ){
    incomplete <- is.na(edges)
    if( any(incomplete) ){
      warning("keep_edges(): rows in 'edges' containing missing values were ignored.",
              call. = FALSE)
      edges <- edges[!incomplete]
      if( length(edges) == 0L ){
        stop("keep_edges(): no complete edges were supplied.", call. = FALSE)
      }
    }
    edges <- parse_edge_strings(edges, "keep_edges()")
  }else if( !(is.data.frame(edges) || data.table::is.data.table(edges)) ){
    stop("'edges' must be a character vector, edge table, or dagitty object.",
         call. = FALSE)
  }

  edges <- normalise_edge_frame(edges, "keep_edges()")
  if( nrow(edges) == 0L ){
    stop("keep_edges(): no edges were supplied.", call. = FALSE)
  }
  if( anyNA(edges) ){
    warning("keep_edges(): rows in 'edges' containing missing values were ignored.",
            call. = FALSE)
    edges <- edges[complete.cases(edges), , drop = FALSE]
  }
  if( nrow(edges) == 0L ){
    stop("keep_edges(): no complete edges were supplied.", call. = FALSE)
  }

  requested_keys <- edge_key(edges)
  existing_keys <- edge_key(existing_edges)
  missing_edges <- edges[!requested_keys %in% existing_keys, , drop = FALSE]
  if( nrow(missing_edges) > 0L ){
    warning("Edges not present in the dag were ignored: ",
            paste(missing_edges$v, missing_edges$e, missing_edges$w,
                  collapse = ", "),
            ".", call. = FALSE)
  }

  dropped <- existing_edges[!existing_keys %in% requested_keys, , drop = FALSE]
  dag <- native_remove_edges(dag, dropped, "keep_edges")

  if( isFALSE(dag_is_acyclic(dag)) ){
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.")
  }

  return(dag)
}


#' Add edges to a dagitty object
#'
#' `add_edges()` adds specified edges between nodes already present in the
#' supplied graph.
#'
#' @param dag A dagitty object.
#' @param edges Edges to add between nodes already present in `dag`. Supply one
#'   or more edges as character strings (for example, `"X -> Y"`), an edge
#'   table such as the output from [get_edges()], or a dagitty object whose
#'   edges will be used. Character input accepts `->`, `<-`, `<->`, and `--`,
#'   including chains and multiple edges separated by semicolons or line
#'   breaks. Edges already present in `dag` are ignored with a warning.
#' @returns A dagitty object with the specified edges added.
#'
#' @details
#' Edge endpoints must be different nodes already present in `dag`. Use
#' [add_nodes()] to introduce a new node. If `edges` is a dagitty object, only
#' its edge endpoints and types are used.
#'
#' Existing graph metadata, node statuses, coordinates, labels, notes, R
#' attributes, and stored control points remain unchanged. Newly added edges
#' are straight. New edges are rejected if they would introduce a directed
#' cycle into an acyclic graph. A graph that already contains a cycle is
#' edited with a warning.
#'
#' @examples
#' dag <- string_to_dag(c("X -> T", "T -> Y"))
#' dag <- add_edges(dag, "X -> Y")
#' get_edges(dag)
#'
#' @seealso [keep_edges()], [remove_edges()], [assess_edges()], [add_nodes()]
#' @export
add_edges <- function(dag, edges){
  if( !dagitty::is.dagitty(dag) ){
    stop("add_edges(): 'dag' must be a dagitty object.", call. = FALSE)
  }
  if( missing(edges) || is.null(edges) || length(edges) == 0L ){
    stop("add_edges(): supply at least one edge in 'edges'.", call. = FALSE)
  }

  existing_edges <- dag_edges(dag, "add_edges()")
  if( dagitty::is.dagitty(edges) ){
    requested <- dag_edges(edges, "add_edges()")
  }else if( is.character(edges) ){
    incomplete <- is.na(edges)
    if( any(incomplete) ){
      warning("add_edges(): rows in 'edges' containing missing values were ignored.",
              call. = FALSE)
      edges <- edges[!incomplete]
    }
    if( length(edges) == 0L || all(!nzchar(trimws(edges))) ){
      stop("add_edges(): 'edges' contains no complete edges.", call. = FALSE)
    }
    requested <- parse_edge_strings(edges, "add_edges()")
  }else if( is.data.frame(edges) || data.table::is.data.table(edges) ){
    requested <- normalise_edge_frame(edges, "add_edges()")
    incomplete <- !complete.cases(requested)
    if( any(incomplete) ){
      warning("add_edges(): rows in 'edges' containing missing values were ignored.",
              call. = FALSE)
      requested <- requested[!incomplete, , drop = FALSE]
    }
  }else{
    stop("add_edges(): 'edges' must be a character vector, a v/e/w edge table, or a dagitty object.",
         call. = FALSE)
  }

  requested <- normalise_edge_frame(requested, "add_edges()")
  if( nrow(requested) == 0L ){
    stop("add_edges(): 'edges' contains no complete edges.", call. = FALSE)
  }
  requested_keys <- edge_key(requested)
  requested <- requested[!duplicated(requested_keys), , drop = FALSE]
  requested_keys <- edge_key(requested)

  unknown <- setdiff(unique(c(requested$v, requested$w)), names(dag))
  if( length(unknown) > 0L ){
    stop("add_edges(): edge endpoints are not in the graph: ",
         paste(unknown, collapse = ", "),
         ". Use add_nodes() to introduce new nodes.", call. = FALSE)
  }
  self_edges <- requested$v == requested$w
  if( any(self_edges) ){
    stop("add_edges(): edges must connect different nodes: ",
         paste(requested$v[self_edges], requested$e[self_edges],
               requested$w[self_edges], collapse = ", "),
         ".", call. = FALSE)
  }

  existing_keys <- edge_key(existing_edges)
  already_present <- requested_keys %in% existing_keys
  present_edges <- requested[already_present, , drop = FALSE]
  new_edges <- requested[!already_present, , drop = FALSE]
  warn_present <- function(){
    if( nrow(present_edges) > 0L ){
      out <- warning("add_edges(): edges already present in the graph were ignored: ",
                     paste(present_edges$v, present_edges$e, present_edges$w,
                           collapse = ", "),
                     ".", call. = FALSE)

      return(invisible(out))
    }else{
      out <- NULL
      return(invisible(out))
    }
  }
  if( nrow(new_edges) == 0L ){
    warn_present()
    return(dag)
  }

  candidate <- tryCatch({
    text <- as.character(dag)
    closing <- native_layout_index(text, "add_edges()")$closing
    lines <- paste(native_layout_quote(new_edges$v), new_edges$e,
                   native_layout_quote(new_edges$w))
    text <- native_layout_edits(text, list(list(
      start = closing,
      end = closing - 1L,
      text = paste0("\n", paste(lines, collapse = "\n"), "\n")
    )))
    native_layout_value(dag, text)
  }, error = function(error){
    stop("add_edges(): the graph could not be updated without changing its existing metadata.",
         call. = FALSE)
  })

  actual_edges <- tryCatch(dag_edges(candidate, "add_edges()"),
                           error = function(error) {
                             out <- NULL

                             return(out)
                           })
  candidate_nodes <- tryCatch(names(candidate), error = function(error) {
    out <- NULL

    return(out)
  })
  expected_keys <- unique(c(existing_keys, edge_key(new_edges)))
  actual_keys <- if( is.null(actual_edges) ) { character() } else { edge_key(actual_edges) }
  if( is.null(candidate_nodes) ||
      !setequal(candidate_nodes, names(dag)) ||
      length(candidate_nodes) != length(names(dag)) ||
      !setequal(actual_keys, expected_keys) ||
      length(unique(actual_keys)) != length(expected_keys) ){
    stop("add_edges(): dagitty could not represent the requested edges exactly.",
         call. = FALSE)
  }
  if( !isTRUE(dag_is_acyclic(candidate)) ){
    if( isTRUE(dag_is_acyclic(dag)) ){
      stop("add_edges(): the resulting graph would contain a directed cycle.",
           call. = FALSE)
    }
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.",
            call. = FALSE)
  }

  warn_present()
  return(candidate)
}


#' Remove edges from a dagitty object
#'
#' `remove_edges()` removes specified edges and any nodes left without edges
#' by that removal. Nodes that were already isolated are kept.
#'
#' @param dag A dagitty object.
#' @param edges Edges to remove. Supply one or more edges as character strings
#'   (for example, `"X -> Y"`), an edge table such as the output from
#'   [get_edges()], or a dagitty object whose edges will be used. Character
#'   input accepts `->`, `<-`, `<->`, and `--`, including chains and multiple
#'   edges separated by semicolons or line breaks. Supplied edges not present
#'   in `dag` are ignored with a warning.
#' @returns A dagitty object with the specified edges and newly isolated nodes
#'   removed.
#' @seealso [add_edges()], [keep_edges()], [remove_nodes()], [assess_edges()]
#'
#' @details
#' Existing graph metadata, node statuses, coordinates, labels, notes, R
#' attributes, and stored control points remain unchanged apart from those
#' attached to the removed edges and nodes. Nodes removed automatically are
#' named in a warning. A node is kept if it has any remaining directed,
#' bidirected, or undirected edge. Use [keep_edges()] to discard edges while
#' retaining every node, or [remove_nodes()] to remove named nodes directly.
#'
#' @examples
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
#'   type = "saturated"
#' )
#'
#' dag <- dag |>
#'   remove_edges(c("X1 -> X2", "X2 -> T"))
#'
#' @export
remove_edges <- function(dag,
                         edges
                         ){
  if( !dagitty::is.dagitty(dag) ){
    stop("Invalid input detected: only dagitty objects may be supplied as 'dag' input. Please check the input and try again.",
         call. = FALSE)
  }

  existing_edges <- dag_edges(dag, "remove_edges()")

  if( dagitty::is.dagitty(edges) ){
    edges <- dag_edges(edges, "remove_edges()")
  }else if( is.character(edges) ){
    incomplete <- is.na(edges)
    if( any(incomplete) ){
      warning("remove_edges(): rows in 'edges' containing missing values were ignored.",
              call. = FALSE)
      edges <- edges[!incomplete]
    }
    if( length(edges) == 0L || all(!nzchar(trimws(edges))) ){
      stop("Invalid input detected: 'edges' must contain at least one edge. Please check the input and try again.",
           call. = FALSE)
    }
    edges <- parse_edge_strings(edges, "remove_edges()")
  }else if( !(is.data.frame(edges) || data.table::is.data.table(edges)) ){
    stop("Invalid input detected: 'edges' must be a character vector, an edge table, or a dagitty object. Please check the input and try again.",
         call. = FALSE)
  }

  edges <- normalise_edge_frame(edges, "remove_edges()")
  if( nrow(edges) == 0L ){
    stop("Invalid input detected: 'edges' must contain at least one edge. Please check the input and try again.",
         call. = FALSE)
  }
  if( anyNA(edges) ){
    warning("remove_edges(): rows in 'edges' containing missing values were ignored.",
            call. = FALSE)
    edges <- edges[complete.cases(edges), , drop = FALSE]
  }
  if( nrow(edges) == 0L ){
    stop("Invalid input detected: 'edges' contains no complete edges. Please supply at least one complete edge and try again.",
         call. = FALSE)
  }

  requested_keys <- edge_key(edges)
  existing_keys <- edge_key(existing_edges)
  missing_edges <- edges[
    !requested_keys %in% existing_keys & !duplicated(requested_keys), ,
    drop = FALSE
  ]
  if( nrow(missing_edges) > 0L ){
    warning("Edges not present in the supplied dagitty object were ignored: ",
            paste(missing_edges$v, missing_edges$e, missing_edges$w,
                  collapse = ", "),
            ".", call. = FALSE)
  }

  matched <- existing_keys %in% requested_keys
  if( !any(matched) ){
    return(dag)
  }

  dag <- native_remove_edges(dag, existing_edges[matched, , drop = FALSE],
                              "remove_edges")
  remaining_edges <- existing_edges[!matched, , drop = FALSE]
  connected <- unique(c(existing_edges$v, existing_edges$w))
  remaining_connected <- unique(c(remaining_edges$v, remaining_edges$w))
  isolated <- names(dag)[names(dag) %in% setdiff(connected, remaining_connected)]
  if( length(isolated) > 0L ){
    dag <- native_remove_nodes(dag, isolated, "remove_edges")
    warning("remove_edges(): removed nodes with no remaining edges: ",
            paste(isolated, collapse = ", "), ".", call. = FALSE)
  }

  if( isFALSE(dag_is_acyclic(dag)) ){
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.",
            call. = FALSE)
  }

  return(dag)
}


#' Remove nodes from a dagitty object
#'
#' Removes named nodes and all their incident edges. Other nodes are kept,
#' including neighbours left without edges.
#'
#' @param dag A dagitty object.
#' @param nodes Character vector of node names to remove. Repeated names are
#'   ignored. Names not present in the graph are ignored with a warning.
#' @returns A dagitty object without the specified nodes and their edges.
#' @seealso [add_nodes()], [remove_edges()], [keep_edges()]
#'
#' @details
#' Statuses, coordinates, labels, notes, and stored control points attached
#' to removed nodes or edges are removed with them. The remaining graph keeps
#' its metadata and coordinates, and no layout is regenerated. Removing every node
#' returns an empty dagitty object of the same graph type.
#'
#' @examples
#' dag <- build_graph(confounders = "X", treatments = "T", outcomes = "Y")
#' dag <- remove_nodes(dag, "X")
#'
#' @export
remove_nodes <- function(dag, nodes){
  dag_graph_type(dag)
  if( !is.character(nodes) ){
    stop("remove_nodes(): 'nodes' must be a character vector.", call. = FALSE)
  }
  if( length(nodes) == 0L ){
    stop("remove_nodes(): 'nodes' must contain at least one node name.",
         call. = FALSE)
  }
  if( anyNA(nodes) || any(!nzchar(nodes)) ){
    stop("remove_nodes(): 'nodes' must contain non-missing, non-empty names.",
         call. = FALSE)
  }
  nodes <- unique(as.character(nodes))
  missing_nodes <- setdiff(nodes, names(dag))
  if( length(missing_nodes) > 0L ){
    warning("remove_nodes(): nodes not present in the graph were ignored: ",
            paste(missing_nodes, collapse = ", "), ".", call. = FALSE)
  }
  nodes <- intersect(nodes, names(dag))
  if( length(nodes) == 0L ){
    return(dag)
  }

  dag <- native_remove_nodes(dag, nodes, "remove_nodes")
  if( isFALSE(dag_is_acyclic(dag)) ){
    warning("The outputted graph contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.",
            call. = FALSE)
  }

  return(dag)
}


## Edge assessment helpers ##
#' Check whether edge assessment can prompt
#'
#' @return A logical value indicating whether the R session is interactive.
#' @noRd
assessment_is_interactive <- function(){
  interactive_out <- interactive()

  return(interactive_out)
}


#' Read an edge assessment response
#'
#' @param prompt Text displayed before reading the response.
#' @return The response as a character string.
#' @noRd
assessment_readline <- function(prompt){
  answer <- readline(prompt)

  return(answer)
}


#' Assess graph edges using causal criteria
#' @importFrom data.table data.table
#' @param edges Edge table with columns v, e and w to assess.
#' @param check_skip_sequence Whether to skip the guided assessment.
#' @param edges_to_keep Edge table retained without assessment.
#' @param causal_criteria Table of criteria supplied by assess_edges().
#' @param save_answers A decision-log table to replay, a non-empty scalar
#'   to return the new log, or an empty vector to return only edges.
#' @return A retained edge table, or a list containing edges and answers when
#'   recording the assessment. Declining the assessment returns a list
#'   containing edges and edges_to_assess. Skipping returns the supplied edges.
#' @noRd
causal_criteria_sequence <- function(edges,
                                     check_skip_sequence,
                                     edges_to_keep,
                                     causal_criteria,
                                     save_answers
                                     ){
  .datatable.aware <- TRUE
  edges_list <- edges   # ensure the final return() is always bound; degenerate skip path returns edges unchanged

  if( length(save_answers) >= 3 ){

    if( is.data.frame(save_answers) && nrow(save_answers) == nrow(edges) ){

      len_answers <- length(save_answers)
      if( all( as.data.frame(save_answers)[,1:3] == as.data.frame(edges) ) & len_answers == (nrow(causal_criteria) + 3) ){
        answers <- data.table::as.data.table(save_answers)[, 4:len_answers] # coerce first so a single criterion keeps its column name
        cols_required <- unlist(causal_criteria[,"name"])[ unlist(causal_criteria[, "required" ]) == "yes" ]
        if( length(cols_required) == 0L ){
          keep_edge <- rep(TRUE, nrow(save_answers))
        }else{
          answers <- answers[, ..cols_required]
          keep_edge <- vapply(seq_len(nrow(save_answers)), function(x){
            matches <- isTRUE(all((answers == "y")[x, ]))

            return(matches)
          }, logical(1L))
        }
        edges <- edges[keep_edge, ]

        if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
            all(complete.cases(edges_to_keep)) ){
          edges <- rbind(edges_to_keep[, 1:3], edges)
        }

        message("\nOutputting edges. Causal criteria decision log printed.")
        print(answers)

        return(edges)
      }

    }

    stop("Length of answers does not match the causal criteria provided. Please check inputs and try again.")
  }

  if( !assessment_is_interactive() ){
    stop("Interactive causal assessment requires an interactive R session. ",
         "Supply 'save_answers' to replay a completed assessment.",
         call. = FALSE)
  }

  if(check_skip_sequence == FALSE) {

    num_edges <- nrow(edges)
    if( num_edges < 1 ){
      message( "There are no directed arrows to assess." )
      if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
          all(complete.cases(edges_to_keep)) ){
        edges <- rbind(edges_to_keep[, 1:3], edges)
      }

      return(edges)
    }
    set_noun <- if( any( unlist( edges[, "e"] ) %in% c("<->", "--") ) ) { "edges" } else { "directed arrows" }
    cat( "\nThere are ", num_edges, " ", set_noun, " to be assessed.","\n", "\n", sep="" )
    print(edges, quote=FALSE)
    cat( "\nAssess the posited causal relationships using causal criteria?", "\n" )

    check_ans <- FALSE

    while(check_ans == FALSE){

      choice <- assessment_readline("(y/n/?info): ")
      if(choice == "y"){
        check_ans <- TRUE
      }else if(choice == "n"){
        message( "Skipped sequence." )
        edges_list <- list(edges = edges_to_keep,
                           edges_to_assess = edges)

        return(edges_list)
      }else if(choice == "?info"){
        cat( "This is a guided sequence for assessing graph edges. By default it uses the causal criteria in ESC-DAGs protocol (Ferguson et al., 2020).
             \nA data frame containing user-specified criteria can be supplied in assess_edges().")
      }else{
        message( "Please type a valid answer." )
      }

    }

    removed_arrows <- c()
    arrow_count <- 1
    num_arrow_to_remove <- 0
    num_criteria <- nrow(causal_criteria)
    criteria_answers <- as.data.frame( matrix(nrow = 0, ncol = num_criteria) )
    colnames(criteria_answers) <- causal_criteria[, "name" ]
    causal_criteria$required <- as.factor(causal_criteria$required)

    for(arrow_count in 1:num_edges){

        criterion_num <- 1
        check_ans <- FALSE
        arrow <- noquote( paste( edges[ arrow_count ], collapse=" " ) )
        arrow_noun <- if( unlist( edges[ arrow_count, "e" ] ) %in% c("<->", "--") ) { "edge" } else { "directed arrow" }
        cat( "For the ", arrow_noun, " '", arrow, "' consider each of the following:", "\n", "\n", sep="")

        while( check_ans == FALSE ){

          if( !criterion_num > num_criteria ){

            cat( paste0("'", arrow, "' (", arrow_count, "/", num_edges, ")"),
                 "\n", paste0("\n[", criterion_num, "/", num_criteria, "]"),
                 paste0( causal_criteria[ criterion_num, "name" ], ":"),
                 causal_criteria[ criterion_num, "question" ], "\n",
                 "\nFor help, enter \'?info\'")
            choice <- assessment_readline("(y/n/?info): ")
            if( choice == "y" ){
              criteria_answers[ arrow_count, criterion_num] <- "y" # "yes (causal)"
              criterion_num <- criterion_num + 1
            }else if( choice == "n" ){

              if( causal_criteria[ criterion_num, "required" ] == "yes"){
                num_arrow_to_remove <- num_arrow_to_remove + 1
                removed_arrows[ num_arrow_to_remove ] <- arrow_count

                criteria_answers[ arrow_count, criterion_num ] <- "n" # "no - remove edge"
                criterion_num <- num_criteria + 1

                message( "Causal relationship '", arrow, "' assessed; edge removed.", sep="" )

                check_ans <- TRUE
              }else{
                criteria_answers[ arrow_count, criterion_num ] <- "n"
                criterion_num <- criterion_num + 1

                message( "Answer recorded. Causal relationship '", arrow, "' assessed.", sep="" )
              }

            }else if( choice == "?info" ){
              cat( "\n", "\nCausal criterion", criterion_num, "\u2014", causal_criteria[ criterion_num, "name" ],
                   "\n", "\n", paste0("Definition", ":"),
                   "\n", causal_criteria[ criterion_num, "description" ],
                   "\n", "\n", causal_criteria[ criterion_num, "source" ],
                   "\n", "\n" )
            }else{
              cat( "\n", "\nPlease type a valid answer.", "\n" )
            }

          }else{
            check_ans <- TRUE
          }

        }

      }

    if(num_arrow_to_remove > 0){

      criteria_answers <- cbind(edges, criteria_answers)
      edges <- edges[ -removed_arrows, ]

      if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
          all(complete.cases(edges_to_keep)) ){
        edges <- rbind(edges_to_keep[, 1:3], edges)
      }

      if( length(save_answers) == 0 ){
        message("\nOutputting edges. Causal criteria decision log printed.")
        print(criteria_answers)

        return(edges)
      }
      edges <- list(edges = edges,
                    answers = criteria_answers)

      message("\nOutputting list containing 'edges' and 'answers'.")

      return(edges)
    }else{

      message( "No arrows were removed." )

      criteria_answers <- cbind(edges, criteria_answers)

      if( is.data.frame(edges_to_keep) && ncol(edges_to_keep) >= 3L &&
          all(complete.cases(edges_to_keep)) ){
        edges <- rbind(edges_to_keep[, 1:3], edges)
      }

      if( length(save_answers) == 0 ){
        print(criteria_answers)
        return(edges)
      }

      edges_list <- list(edges = edges,
                         answers = criteria_answers)

      message("\nOutputting list containing 'edges' and 'answers'.")
      return(edges_list)

    }

  }

  return(edges_list)
}


## Node connection helpers ##
#' Instrument helper for connect_new_nodes()
#'
#' Called by connect_new_nodes() when add_nodes() adds instruments.
#' Other roles are planned by role_operation_edges().
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty topologicalOrdering
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param nodes New node names.
#' @param node_role Role assigned to the new nodes. Must be "instrument".
#' @param type Graph type passed by add_nodes(), including the deprecated
#'   `"first"` and `"last"` values.
#' @param position Optional placement. Not supported for instruments.
#' @param temporal_reference_node Returned unchanged.
#' @returns output_list containing temporal_reference_node, existing_node_names, and a data frame of new_edges.
#' @noRd
add_nodes_helper <- function(dag,
                             nodes,
                             node_role,
                             type = "full",
                             position = NULL,
                             temporal_reference_node = NA
                             ){
  e <- "->"

  if(type == "full"){

    e <- "<->"
  }

  ## get initial dag roles ##
  dag_roles <- classify_node_roles(dag)

  nodes_ordered <- sort( unlist( dagitty::topologicalOrdering(dag) ) ) # Existing graph nodes in topological order

  dag_roles <- lapply(dag_roles, function(x) {
    ordered_nodes <- x[order(match(x, names(nodes_ordered)))]

    return(ordered_nodes)
  })

  dag_roles <- lapply(dag_roles, function(x) {
    nodes <- x[complete.cases(x)]
    if( identical(nodes, logical(0)) ){
      nodes <- NA_character_
    }

    return(nodes)
  })

  outcomes  <- dag_roles$outcome[complete.cases(dag_roles$outcome)]
  treatments <- dag_roles$treatment[complete.cases(dag_roles$treatment)]
  instrumental_variables <- dag_roles$instrument[complete.cases(dag_roles$instrument)]
  latent_vec <- dag_roles$latent[complete.cases(dag_roles$latent)]

  if( length(node_role) == 1 && !node_role %in% VALID_NODE_ROLES ){
    stop("Invalid node_role: \"", node_role, "\". Must be one of: ",
         paste(VALID_NODE_ROLES, collapse = ", "), ".", call. = FALSE)
  }

  if( !is.null(position) && length(node_role) == 1 &&
      !node_role %in% c("confounder", "treatment", "mediator") ){
    stop("The position argument is not supported for role \"", node_role,
         "\"; it would be ignored. Remove position, or use a supported role.",
         call. = FALSE)
  }

  if( length( node_role) > 1 | length( node_role) == 0 ){

    stop("add_nodes() currently only supports single node_role character inputs.")

  }
  if( !identical(node_role, "instrument") ){
    stop("add_nodes_helper() only connects instrument nodes; add_nodes() plans ",
         "other roles with role_operation_edges().", call. = FALSE)
  }

  ## instrumental_variables edges ##
  new_edges <- draw_iv_edges(type = type,
                             instrumental_variables = nodes, # new nodes as instrumental
                             treatments = treatments,
                             e = e)

  existing_node_names <- names( nodes_ordered[ names(nodes_ordered) %in% instrumental_variables ] ) # existing node names in temporal order

  output_list <- list(temporal_reference_node = temporal_reference_node,
                      existing_node_names = existing_node_names,
                      new_edges = new_edges,
                      treatments = treatments,
                      outcomes = outcomes,
                      latent_vec = latent_vec)


  return(output_list)
}


#' Generate edges involving selected nodes
#'
#' connect_nodes_helper() generates candidate edges involving selected nodes
#' according to the construction type.
#'
#' @importFrom data.table as.data.table
#' @param dag An existing dagitty object.
#' @param nodes Existing node names selected for new connections.
#' @param dag_node_names Existing graph node names; retained but unused.
#' @param type Construction type. `"first"` and `"last"` use `"ordered"`.
#' @return An edge table with columns v, e and w.
#' @noRd
connect_nodes_helper <- function(dag, nodes, dag_node_names, type){
  dag_roles <- classify_node_roles(dag)
  instrumental_variables <- dag_roles$instrument
  instrumental_variables <- instrumental_variables[
    !is.na(instrumental_variables)
  ]
  planner_type <- if( type %in% c("first", "last") ) { "ordered" } else { type }
  edges <- if( identical(type, "full") ){
    full_pair_edges(dag, nodes)
  }else{
    plan_graph_edges(
      existing_role_specs(dag),
      type = planner_type,
      instrumental_variables = instrumental_variables,
      existing_dag = dag,
      include_iv = FALSE
    )
  }

  iv_source <- unlist(edges[, "v"]) %in% instrumental_variables |
    (unlist(edges[, "e"]) %in% c("<->", "--") &
       unlist(edges[, "w"]) %in% instrumental_variables)
  edges <- edges[!iv_source, ]

  ## keep candidate edges involving a selected node ##
  edges <- edges[ unlist(edges[, "v"]) %in% nodes | unlist(edges[, "w"]) %in% nodes, ]

  return( edges )
}
