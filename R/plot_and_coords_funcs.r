
# plot_dagitty()
#' ggdag plot from a dagitty object
#'
#' Generates a ggdag plot from a dagitty object, adding coordinates and node
#' labels when needed. Essentially a wrapper for dagitty and ggdag plotting
#' functions with added features.
#'
#' @importFrom dagitty coordinates
#' @importFrom ggplot2 ggplot scale_shape_manual scale_size_manual scale_colour_manual labs guides guide_legend theme
#' @importFrom ggdag geom_dag_edges geom_dag_text geom_dag_label_repel geom_dag_point theme_dag
#' @importFrom ggraph circle
#' @importFrom grid unit arrow
#' @param dag A dagitty object
#' @param include_legend Node roles to include in legend. This parameter controls the DAG node colour scale and can be overridden by setting include_legend = FALSE. To include all possible node roles use include_legend = c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrument", "competing_cause", "collider", "latent", "observed", "undetermined").
#' @param labels A vector of labels for nodes in dagitty object
#' @param label_type Label style passed to \code{get_labels()} when \code{labels}
#'   is not supplied. Defaults to \code{"name"}; \code{"initials"} and
#'   \code{"short"} are also available.
#' @param label_placement Labels are placed in a text box adjacent each node by default (label_placement = "text_box"), or can be placed inside each node using 'node'. It is recommended to also specify label_type = 'initials' where label_placement = 'node' is used.
#' @param seed Numeric input, sets seed for label placement (passed to ggdag::geom_dag_label_repel() seed parameter).
#' @param x_step Horizontal spacing between consecutive temporal ranks when coordinates are generated. Single positive number.
#' @param y_step Minimum vertical gap between nodes within a rank when coordinates are generated. Single positive number.
#' @param layout_shape Within-rank placement between 0 and 1: 0 places nodes nearest their parents' pathways, 1 spreads each rank evenly.
#' @return ggdag plot
#' @examples
#' set.seed(1)
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' plot_dagitty(dag)
#'
#' plot_dagitty(dag, y_step = 1.5) # more vertical spacing between nodes
#'
#' @export
plot_dagitty <- function(dag,
                         seed = NULL,
                         x_step = 2,
                         y_step = 1,
                         layout_shape = 0.5,
                         labels = NULL,
                         label_type = "name",
                         label_placement = "text_box",
                         include_legend = c("outcome",
                                            "treatment",
                                            "confounder",
                                            "mediator",
                                            "instrument",
                                            "competing_cause",
                                            "collider",
                                            "latent",
                                            "undetermined")
                         ){

  # Share structural reads between the edge check and role extraction.
  .cache <- .new_cache()

  if( nrow( data.table::as.data.table( .cached_edges(dag, .cache) ) ) == 0 ){
    stop("plot_dagitty() requires a graph with at least one edge. Roles and coordinates for edge-less graphs remain available via get_roles() and add_coords().")
  }

  roles <- .get_roles(dag, .cache = .cache)

  graph_type <- .dag_graph_type(dag)
  if( graph_type != "dag" ){
    original_coordinates <- dagitty::coordinates(dag)
    dag <- construct_graph(
      .dag_edges(dag, "plot_dagitty()"), names(dag),
      treatments(dag), .outcomes(dag), unobserved(dag),
      adjusted_vec = dagitty::adjustedNodes(dag),
      selected_vec = .selected_nodes(dag),
      graph_type = "dag"
    )
    coordinate_names <- intersect(names(original_coordinates$x),
                                  names(original_coordinates$y))
    positioned <- coordinate_names[
      !is.na(original_coordinates$x[coordinate_names]) &
        !is.na(original_coordinates$y[coordinate_names])
    ]
    if( length(positioned) > 0L ){
      dagitty::coordinates(dag) <- list(
        x = original_coordinates$x[positioned],
        y = original_coordinates$y[positioned]
      )
    }
  }
  existing_coordinates <- dagitty::coordinates(dag)

  if( any( is.na( existing_coordinates$x ) ) | any( is.na( existing_coordinates$y ) ) ){

    dag <- .add_coords(dag, x_step = x_step, y_step = y_step, layout_shape = layout_shape)
  }

  # aesthetics for ggdag legend
  col <- c("outcome"="deepskyblue",
           "treatment"="darkolivegreen2",
           "confounder"="coral2",
           "mediator"="darkorchid1",
           "mediator_outcome_confounder"="magenta4",
           "instrument"="deeppink1",
           "competing_cause"="darkseagreen4",
           "collider"="darkred",
           "latent"="black",
           "undetermined"="grey60",
           "observed"="#111111")

  col <- sapply(1:length(col), function(x){
    if( !names(col[x]) %in% include_legend ){
      col[x] <- "#111111"
    }else{ col[x] <- col[x] }
    col[x]
  })

  shape <- c("outcome"=19,
             "treatment"=19,
             "confounder"=19,
             "mediator"=19,
             "mediator_outcome_confounder"=19,
             "instrument"=19,
             "collider"=19,
             "competing_cause"=19,
             "latent"=21,
             "undetermined"=19,
             "observed"=19)

  # variable for legend order
  order_col <- c("outcome",
                 "treatment",
                 "confounder",
                 "mediator",
                 "mediator_outcome_confounder",
                 "instrument",
                 "competing_cause",
                 "collider",
                 "latent",
                 "undetermined",
                 "observed")
  order_col <- order_col[ order_col %in% include_legend ]

  if( is.null(labels) ){
    labels <- get_labels(dag, label_type)
  }

  if(!is.null(labels) & length(labels) != length(names(dag))){
    stop("The length of supplied labels does not equal the number of nodes in the graph. Please check labels input and try again.")
  }

  dag_df <- tidy_ggdagitty(dag, labels, roles = roles)
  dag_df_complete_cases <- dag_df[complete.cases(dag_df[, "role"]), ]

  if( nrow(dag_df) > nrow(dag_df_complete_cases) ){
    message("Nodes without arrows removed and not displayed in the plotted graph.")
    dag_df <- dag_df_complete_cases
  }

  dag_df[sapply(dag_df, is.character)] <- lapply( dag_df[sapply(dag_df, is.character)], as.factor )

  if(label_placement == "text_box"){

    label_col <- c("outcome"="#FFFFFF",
                   "treatment"="#FFFFFF",
                   "confounder"="#FFFFFF",
                   "mediator"="#FFFFFF",
                   "mediator_outcome_confounder"="#FFFFFF",
                   "instrument"="#FFFFFF",
                   "competing_cause"="#FFFFFF",
                   "collider"="#FFFFFF",
                   "latent"="#FFFFFF",
                   "undetermined"="#FFFFFF",
                   "observed"="#FFFFFF")

    ggdag <- ggplot2::ggplot(data = dag_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color=role, shape = role, fill = role)) +
      ggdag::geom_dag_edges(arrow_directed = grid::arrow(angle = 25, length = grid::unit(5, "pt"), type = "closed"),
                            arrow_bidirected = grid::arrow(angle = 25, length = grid::unit(5, "pt"), ends = "both",
                                                           type = "closed"),
                            start_cap = ggraph::circle(5.3, 'mm'),
                            end_cap = ggraph::circle(5.8, 'mm')) +
      ggdag::geom_dag_point(size=13) +
      ggdag::geom_dag_label_repel(ggplot2::aes(label = label), colour = "black", alpha = 0.65,
                                  box.padding = grid::unit(1, "lines"),
                                  label.padding = grid::unit(0.2, "lines"),
                                  point.padding = grid::unit(1.9, "lines"),
                                  label.r = grid::unit(0.1, "lines"),
                                  seed = seed,
                                  label.size = 0.1,
                                  show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = col,
                                   name = "Group",
                                   breaks = order_col) +
      ggplot2::scale_shape_manual(values = shape,
                                  name = "Group",
                                  breaks = order_col) +
      ggplot2::scale_fill_manual(values = label_col,
                                 name = "Group",
                                 breaks = order_col
      ) +
      ggdag::theme_dag() +
      ggplot2::labs(x = "", y="") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 8))) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        # title = element_text(size = 16),
        #legend.text = ggplot2::element_text(size = 6),               # Increase the legend text size
      )

    return(ggdag)

  }else if(label_placement == "node"){

    dag_df_subset <- subset(dag_df, !duplicated(dag_df$label))
    dag_df_subset$label_col <- unlist( lapply(1:nrow(dag_df_subset), function(x){
      if(dag_df_subset[x,"role"] == "observed"){
        dag_df_subset[x,"label_col"] <- "latent"
      }else{
        dag_df_subset[x,"label_col"] <- "observed"
      }
      dag_df_subset[x,"label_col"]
    }) )

    ggdag <- ggplot2::ggplot(data = dag_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color=role, shape = role)) +
      ggdag::geom_dag_edges(arrow_directed = grid::arrow(angle = 25, length = grid::unit(5, "pt"), type = "closed"),
                            arrow_bidirected = grid::arrow(angle = 25, length = grid::unit(5, "pt"), ends = "both",
                                                           type = "closed"),
                            start_cap = ggraph::circle(5.1, 'mm'),
                            end_cap = ggraph::circle(5.9, 'mm')) +
      ggdag::geom_dag_point(size=13) +
      ggdag::geom_dag_text(data = dag_df_subset, mapping = ggplot2::aes(label = label, color = label_col),
                           size = 3, show.legend = FALSE)+
      ggplot2::scale_colour_manual(values = col,
                                   name = "Group",
                                   breaks = order_col) +
      ggplot2::scale_shape_manual(values = shape,
                                  name = "Group",
                                  breaks = order_col) +
      ggdag::theme_dag() +
      ggplot2::labs(x = "", y="") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 8))) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        #legend.text = ggplot2::element_text(size = 6),               # Increase the legend text size
      )

    return(ggdag)

  }else{

    stop("The label_placement input is invalid. Please use either 'node' or 'text_box'.")
  }
}


#' Generate new dag coordinates
#'
#' Generates coordinates deterministically: nodes are ranked left to right by
#' longest directed path (temporal order) and spaced vertically within each
#' rank. The same graph always receives the same coordinates. A graph containing
#' a directed cycle is laid out on an acyclic subset of its edges, with a warning.
#'
#' @param dag A dagitty object.
#' @param x_step Horizontal spacing between consecutive temporal ranks. Single positive number.
#' @param y_step Minimum vertical gap between nodes within a rank. Single positive number.
#' @param layout_shape Within-rank placement between 0 and 1: 0 places nodes nearest their parents' pathways, 1 spreads each rank evenly.
#' @return dagitty object with coordinates.
#' @examples
#' dag <- build_graph(
#'   variables = c("Z1", "Z2"), treatments = "X", outcomes = "Y",
#'   mediators = "M", instrumental_variables = "IV", type = "saturated"
#' )
#' dag <- add_coords(dag) # update dagitty object node coordinates
#'
#' plot_dagitty(dag) # check coordinates
#'
#' @export
add_coords <- function(dag,
                       x_step = 2,
                       y_step = 1,
                       layout_shape = 0.5
                       ){
  .add_coords(dag, x_step = x_step, y_step = y_step, layout_shape = layout_shape)
}

.add_coords <- function(dag,
                       x_step = 2,
                       y_step = 1,
                       layout_shape = 0.5
                       ){

  if( !is.numeric(x_step) || length(x_step) != 1L || !is.finite(x_step) || x_step <= 0 ){
    stop("add_coords(): x_step must be a single positive number.", call. = FALSE)
  }
  if( !is.numeric(y_step) || length(y_step) != 1L || !is.finite(y_step) || y_step <= 0 ){
    stop("add_coords(): y_step must be a single positive number.", call. = FALSE)
  }
  if( !is.numeric(layout_shape) || length(layout_shape) != 1L || !is.finite(layout_shape) ||
      layout_shape < 0 || layout_shape > 1 ){
    stop("add_coords(): layout_shape must be a single number between 0 and 1.", call. = FALSE)
  }

  node_names <- names(dag)
  if( length(node_names) == 0L ){
    return(dag)
  }
  if( nrow(.dag_edges(dag, "add_coords()")) == 0L ){
    dagitty::coordinates(dag) <- list(
      x = stats::setNames(seq_along(node_names), node_names),
      y = stats::setNames(rep(0, length(node_names)), node_names)
    )
    return(dag)
  }

  dagitty::coordinates(dag) <- renew_coords(dag,
                                            x_step = x_step,
                                            y_step = y_step,
                                            layout_shape = layout_shape)

  return(dag)
}




#' Create new dag node coordinates
#'
#' Deterministic layout engine. Nodes are ranked left to right by longest
#' directed path (temporal order); vertical placement within each rank reduces
#' edge crossings and keeps a minimum gap. When existing coordinates are
#' supplied only the new nodes are placed, on the existing horizontal scale;
#' existing positions are never moved.
#'
#' @param dag A dagitty object.
#' @param new_node_names Vector of new node names. Existing nodes are used as
#'   reference points when supplied; otherwise all coordinates are regenerated.
#' @param coordinates list of coordinates from a dagitty object.
#' @param x_step Horizontal spacing between consecutive temporal ranks.
#' @param y_step Minimum vertical gap between nodes within a rank.
#' @param layout_shape Within-rank placement between 0 and 1.
#' @return Dagitty coordinates.
#' @noRd
renew_coords <- function(dag,
                         new_node_names = NULL,
                         coordinates = NULL,
                         x_step = 2,
                         y_step = 1,
                         layout_shape = 0.5
                         ){

  edges <- data.table::as.data.table(.dag_edges(dag, "renew_coords()"))
  pdag_edges <- edges[ ( edges$e == "--" | edges$e == "<->"), ]

  dag_node_names <- names(dag)
  parents <- .layout_parents(dag)
  node_rank <- topological_rank(parents)

  new_node_names <- as.vector( unlist(new_node_names) )

  if( is.null(coordinates) || length(unlist(coordinates)) == 0 ){ # if no coordinates exist

    ## full relayout: rank -> x grid, then vertical spacing within each rank ##
    newly_placed <- dag_node_names
    coordinates <- list(
      x = stats::setNames( node_rank[ dag_node_names ] * x_step, dag_node_names ),
      y = stats::setNames( rep(0, length(dag_node_names)), dag_node_names )
    )
    coordinates <- respace_y_taper(coordinates, parents,
                                   y_step = y_step,
                                   shape = layout_shape,
                                   node_rank = node_rank)

  }else{ # partial/complete coordinates: place only the unpositioned nodes

    x <- coordinates$x
    y <- coordinates$y
    placed <- intersect( names(x)[ !is.na(x) ], names(y)[ !is.na(y) ] )
    placed <- placed[ placed %in% dag_node_names ]
    to_place <- unique( c( new_node_names, dag_node_names ) )
    to_place <- to_place[ to_place %in% dag_node_names & !to_place %in% placed ]

    ## match the existing horizontal scale: least-squares map from rank to x ##
    slope <- x_step
    intercept <- 0
    if( length(placed) >= 2L ){
      rank_placed <- node_rank[ placed ]
      if( stats::var(rank_placed) > 0 ){
        fitted_slope <- stats::cov( rank_placed, x[ placed ] ) / stats::var( rank_placed )
        if( is.finite(fitted_slope) && fitted_slope != 0 ){
          slope <- fitted_slope
        }
      }
      intercept <- mean( x[ placed ] ) - slope * mean( node_rank[ placed ] )
    }else if( length(placed) == 1L ){
      intercept <- x[[ placed ]] - slope * node_rank[[ placed ]]
    }

    ## children, for seeding a new node's height from its neighbours ##
    child_map <- vector("list", length(dag_node_names))
    names(child_map) <- dag_node_names
    for( nd in dag_node_names ){
      for( pr in parents[[nd]] ){
        child_map[[pr]] <- c( child_map[[pr]], nd )
      }
    }

    newly_placed <- to_place[ order( node_rank[ to_place ], to_place ) ]
    for( nd in newly_placed ){
      x[ nd ] <- intercept + slope * node_rank[[ nd ]]
      node_parents <- parents[[ nd ]]
      node_parents <- node_parents[ node_parents %in% names(y) ]
      node_parents <- node_parents[ !is.na( y[ node_parents ] ) ]
      node_children <- child_map[[ nd ]]
      node_children <- node_children[ node_children %in% names(y) ]
      node_children <- node_children[ !is.na( y[ node_children ] ) ]
      seed_y <- if( length(node_parents) > 0L ){
        mean( y[ node_parents ] )
      }else if( length(node_children) > 0L ){
        mean( y[ node_children ] )
      }else{
        0
      }
      same_rank <- names(x)[ names(x) %in% dag_node_names & names(x) != nd ]
      same_rank <- same_rank[ node_rank[ same_rank ] == node_rank[[ nd ]] & !is.na( x[ same_rank ] ) ]
      taken <- y[ same_rank ]
      taken <- taken[ !is.na(taken) ]
      while( any( abs( seed_y - taken ) < y_step * 0.999 ) ){ # keep the minimum gap to already placed nodes
        seed_y <- seed_y - y_step
      }
      y[ nd ] <- seed_y
    }
    coordinates <- list( x = x, y = y )
  }

  ## symmetric edges ("--", "<->") sit at the same height where the two ends
  ## occupy different ranks; existing positions are never moved ##
  if( nrow(pdag_edges) > 0 ){
    for( i in seq_len( nrow(pdag_edges) ) ){
      lhs <- as.character( pdag_edges$v[i] )
      rhs <- as.character( pdag_edges$w[i] )
      if( !lhs %in% names(coordinates$y) || !rhs %in% names(coordinates$y) ){ next }
      if( is.na( coordinates$y[[lhs]] ) || is.na( coordinates$y[[rhs]] ) ){ next }
      if( node_rank[[lhs]] == node_rank[[rhs]] ){ next } # same rank: keep the vertical gap
      lhs_new <- lhs %in% newly_placed
      rhs_new <- rhs %in% newly_placed
      if( lhs_new && rhs_new ){
        level <- mean( c( coordinates$y[[lhs]], coordinates$y[[rhs]] ) )
        coordinates$y[lhs] <- level
        coordinates$y[rhs] <- level
      }else if( lhs_new ){
        coordinates$y[lhs] <- coordinates$y[[rhs]]
      }else if( rhs_new ){
        coordinates$y[rhs] <- coordinates$y[[lhs]]
      }
    }
  }

  coordinates$y <- coordinates$y[!duplicated(names(coordinates$y))]
  coordinates$x <- coordinates$x[!duplicated(names(coordinates$x))]

  coordinates <- list(x = round( coordinates$x[order(names(coordinates$x))], 3), y = round( coordinates$y[order(names(coordinates$y))], 3) )

  return(coordinates)
}



#' dataframe output from a dagitty object
#'
#' Generates a table similar to calling ggdag::tidy_dagitty on a ggdag::dagify object
#' The benefit of this function is that it automatically identifies exposure, outcome, confounder, observed and latent variables inputted from dagitty.net, whereas ggdag::tidy_dagitty only does this for ggdag::dagify objects.
#' Output can be used with ggdag to create better looking DAGs from dagitty.net code.
#'
#' @importFrom ggdag tidy_dagitty
#' @param dag dagitty object
#' @param labels vector of labels for nodes in dagitty object
#' @param roles Optional precomputed get_roles(dag) output. Computed internally when NULL.
#' @return dag_df DAG as a dataframe for use with ggdag to create better looking DAGs
#' @noRd
tidy_ggdagitty <- function(dag, labels = NULL, roles = NULL){
  # Cleaning the dags and turning it into a data frame.
  dag_df <- data.frame(ggdag::tidy_dagitty(dag))

  # flip y axis for ggplot
  dag_df$y <- dag_df$y*-1
  dag_df$yend <- dag_df$yend*-1

  if( is.null(labels) ){

    return(dag_df)
  }

  dag_df <- add_labels(dag, dag_df, labels, roles = roles)

  return(dag_df)
}


#' add labels to a dag dataframe
#'
#' Generates a table similar to calling ggdag::tidy_dagitty on a ggdag::dagify() object
#' The benefit of this function is that it automatically identifies exposure, outcome, confounder, observed and latent variables inputted from dagitty.net, whereas ggdag::tidy_dagitty only does this for ggdag::dagify objects.
#' Output can be used with ggdag to create better looking DAGs from dagitty.net code.
#'
#' @param dag dagitty object
#' @param dag_df a dag object converted to data frame using the ggdag::tidy_dagitty() function, or similar.
#' @param labels vector of labels for nodes in dagify object
#' @param roles Optional precomputed get_roles(dag) output. Computed internally when NULL.
#' @return dagify DAG as a dataframe for use with ggdag to create better looking DAGs
#' @noRd
add_labels <- function(dag, dag_df, labels, roles = NULL){

  # Labeling variables
  dag_df$label <- sapply( seq_along( dag_df$name ),
                          function(x) if( dag_df$name[x] %in% attr(labels, "names")) labels[[ dag_df$name[x] ]] )

  if( is.null(roles) ){
    roles <- .get_roles(dag)
  }

  node_roles <- roles
  dag_df$role <- vapply(as.character(dag_df$name), function(node){
    matches <- names(node_roles)[vapply(
      node_roles, function(role_nodes) node %in% role_nodes, logical(1L)
    )]
    if( length(matches) == 0L ) NA_character_ else matches[1L]
  }, character(1L), USE.NAMES = FALSE)

  return(dag_df)
}
