# ==============================================================================
# coords_funcs.r


# Coordinate calculations shared by graph editing and plotting
# ==============================================================================


#' Space ordered points along one axis
#'
#' Uses isotonic regression to find the closest decreasing sequence with the
#' required gap, without changing the supplied order.
#'
#' @param values Numeric vector in the intended order, from top to bottom.
#' @param gap Minimum separation between consecutive values.
#' @return A numeric vector with the same length and names as `values`, with
#'   consecutive values separated by at least `gap`.
#' @noRd
space_ordered <- function( values, gap ){
  n <- length( values )
  if ( n <= 1L ){
    return( values )
  }
  off <- ( seq_len( n ) - 1L ) * gap
  fit <- stats::isoreg( ( -values ) - off )$yf   # non-decreasing fit to -values, in the given order
  out <- -( fit + off )
  names( out ) <- names( values )
  return( out )
}


#' Arrange vertical coordinates in a tapered layout
#'
#' Orders nodes within each rank using their parents' and children's positions
#' to reduce crossings. Vertical bands narrow towards the final rank, allowing
#' enough space for that rank's nodes. Placement blends inherited parent lanes
#' with evenly spaced positions; `y_step` is then enforced. Horizontal
#' coordinates are kept.
#'
#' @param coordinates List with named numeric x and y vectors. The input y
#'   values set the initial order within each rank.
#' @param parents Named list of parent vectors from layout_parents().
#' @param y_step Minimum vertical gap within a rank.
#' @param shape Placement blend from 0 (parent lanes) to 1 (evenly spaced).
#' @param sweeps Number of ordering passes, alternating down and up the graph.
#' @param node_rank Optional precomputed topological_rank(parents) result.
#' @return The coordinates list with adjusted y values.
#' @noRd
respace_y_taper <- function( coordinates,
                             parents,
                             y_step,
                             shape = 0.5,
                             sweeps = 4,
                             node_rank = NULL
                             ){

  shape       <- max( 0, min( 1, shape ) )   # a blend: outside 0..1 would extrapolate past the band
  if( is.null( node_rank ) ){
    node_rank <- topological_rank( parents )
  }
  y           <- coordinates$y
  node_names  <- names( y )
  ranks       <- node_rank[ node_names ]
  if( anyNA( ranks ) ){
    stop( "respace_y_taper(): every node in the coordinates needs a rank." )
  }
  r_max       <- max( ranks )
  rank_levels <- sort( unique( ranks ) )

  ## directed adjacency restricted to the laid-out nodes (parents in, children derived).
  parents_in <- lapply( parents[ node_names ], function( ps ) {
    value <- ps[ ps %in% node_names ]

    return(value)
  } )
  names( parents_in ) <- node_names
  children <- vector( "list", length( node_names ) )
  names( children ) <- node_names
  for ( nd in node_names ){
    for ( p in parents_in[[ nd ]] ){
      children[[ p ]] <- c( children[[ p ]], nd )
    }
  }

  ## order each rank by incoming y, then refine using neighbouring ranks ##
  layer <- lapply( rank_levels, function( r ){
    column <- node_names[ ranks == r ]
    value <- column[ order( -y[ column ] ) ]

    return(value)
  } )
  names( layer ) <- as.character( rank_levels )

  norm_pos <- function( column ){ # Express a column's order on a 0 to 1 scale
    k <- length( column )
    if ( k == 1L ) {
      positions <- 0.5

      return(positions)
    } else {
      positions <- ( seq_len( k ) - 1L ) / ( k - 1L )

      return(positions)
    }
  }

  reorder_layer <- function( layer, neighbours, rank_keys ){ # Order each column by the mean position of its parents or children
    pos <- rep( 0.5, length( node_names ) )
    names( pos ) <- node_names
    for ( key in names( layer ) ){
      pos[ layer[[ key ]] ] <- norm_pos( layer[[ key ]] )
    }
    for ( r in rank_keys ){
      column <- layer[[ as.character( r ) ]]
      if ( length( column ) >= 2L ){
        key <- vapply( column, function( nd ){
          v <- neighbours[[ nd ]]
          if ( length( v ) == 0L ) {
            value <- NA_real_

            return(value)
          } else {
            value <- mean( pos[ v ] )

            return(value)
          }
        }, numeric( 1 ) )
        key[ is.na( key ) ] <- pos[ column ][ is.na( key ) ]   # no neighbours: hold position
        column <- column[ order( key, seq_along( column ) ) ]
        layer[[ as.character( r ) ]] <- column
      }
      pos[ column ] <- norm_pos( column )   # next rank in this sweep sees the updated positions
    }
    return( layer )
  }

  down_ranks <- rank_levels[ rank_levels > min( rank_levels ) ]
  up_ranks   <- rev( rank_levels[ rank_levels < max( rank_levels ) ] )
  for ( s in seq_len( sweeps ) ){
    if ( s %% 2L == 1L ){
      layer <- reorder_layer( layer, parents_in, down_ranks )
    } else {
      layer <- reorder_layer( layer, children, up_ranks )
    }
  }

  ## assign root lanes by order and other lanes by the mean of parent lanes ##
  lane <- rep( NA_real_, length( node_names ) )
  names( lane ) <- node_names

  for ( r in rank_levels ){
    column <- layer[[ as.character( r ) ]]
    k <- length( column )
    if ( r == 0L ){
      lane[ column ] <- if ( k == 1L ) 0 else 1 - ( seq_len( k ) - 1L ) / ( k - 1L ) * 2
    } else {
      for ( nd in column ){
        node_parents <- parents_in[[ nd ]]
        node_parents <- node_parents[ node_parents %in% node_names ]
        lane[ nd ] <- if ( length( node_parents ) > 0L ) mean( lane[ node_parents ] ) else 0
      }
    }
  }

  ## taper towards the outcome while allowing the required gap in each rank ##
  root_fill <- 1.2   # roots sit this factor wider than the widest inner rank
  counts    <- vapply( rank_levels, function( r ) {
    value <- sum( ranks == r, na.rm = TRUE )

    return(value)
  }, integer( 1L ) )
  half_need <- ( counts - 1L ) / 2 * y_step
  h_max     <- max( half_need )

  if ( r_max == 0L ){
    half_height <- half_need
  } else {
    half_height <- pmax( h_max * ( 1 - rank_levels / r_max ), half_need )   # straight taper, never cramped
    for ( i in rev( seq_len( length( half_height ) - 1L ) ) ){              # make non-increasing root -> tip
      if ( half_height[ i ] < half_height[ i + 1L ] ){
        half_height[ i ] <- half_height[ i + 1L ]
      }
    }
    half_height[ 1L ] <- half_height[ 1L ] * root_fill
  }
  names( half_height ) <- as.character( rank_levels )

  ## place each rank within its height limit, keeping node order and gaps ##
  for ( r in rank_levels ){
    column <- layer[[ as.character( r ) ]]
    k <- length( column )
    band <- half_height[[ as.character( r ) ]]

    even <- if ( k == 1L ) 0 else band - ( seq_len( k ) - 1L ) / ( k - 1L ) * 2 * band
    proximity <- unname( lane[ column ] ) * band

    column_y <- ( 1 - shape ) * proximity + shape * even
    names( column_y ) <- column
    y[ column ] <- unname( space_ordered( column_y, y_step ) )
  }

  ## widen the root positions when an inner rank extends beyond them ##
  root_names  <- layer[[ as.character( rank_levels[ 1L ] ) ]]
  inner_names <- node_names[ !node_names %in% root_names ]
  if ( length( inner_names ) > 0L && length( root_names ) > 1L ){
    inner_reach <- max( abs( y[ inner_names ] ) )
    root_reach  <- max( abs( y[ root_names ] ) )
    if ( root_reach > 0 && root_reach < inner_reach * root_fill ){
      y[ root_names ] <- y[ root_names ] * ( inner_reach * root_fill / root_reach )
    }
  }

  coordinates$y <- y
  return( coordinates )
}


#' Build an acyclic parent list for layout
#'
#' Warns when directed cycles are present and omits their back edges from the
#' layout parent list. The graph itself is unchanged.
#'
#' @param dag A dagitty object.
#' @return A named list of parent vectors with no directed cycles.
#' @noRd
layout_parents <- function(dag){

  node_names <- names(dag)
  edges <- data.table::as.data.table(dag_edges(dag, "renew_coords()"))
  directed <- edges[ edges$e == "->", ]
  from <- as.character( directed$v )
  to   <- as.character( directed$w )

  back <- feedback_arcs( node_names, from, to )
  if( length(back) > 0L ){
    warning("Coordinates were generated using an acyclic subset of the supplied (cyclic) graph.",
            call. = FALSE)
    from <- from[ -back ]
    to   <- to[ -back ]
  }

  parents <- lapply(node_names, function(nm){
    value <- sort( unique( from[ to == nm ] ) )

    return(value)
  })
  names(parents) <- node_names

  return( parents )
}


#' Generate new dag coordinates
#'
#' Generates the same coordinates for the same graph and layout settings.
#' A node's horizontal position is based on the longest directed path leading
#' to it, so parents appear to the left of their children. Nodes in the same
#' column are spaced vertically. A graph containing a directed cycle is laid
#' out on an acyclic subset of its edges, with a warning.
#'
#' @param dag A dagitty object.
#' @param x_step Horizontal spacing between node columns. Single positive number.
#' @param y_step Minimum vertical gap between nodes in the same column. Single
#'   positive number.
#' @param layout_shape Vertical placement of nodes in the same column, between
#'   0 and 1: 0 places nodes nearest their parents' paths, 1 spreads the column
#'   evenly.
#'
#' @details
#' Curved edges move with their endpoints, preserving their bend relative to
#' the edge length. This does not call `add_curves()` or check whether the
#' moved paths cross nodes or other edges. To set an exact curve height, call
#' `add_curves()` afterwards with `mode = "absolute"`.
#'
#' @return A dagitty object with coordinates.
#' @examples
#' dag <- build_graph(
#'   variables = c("X1", "X2"), treatments = "T", outcomes = "Y",
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
  dag_out <- layout_coords(dag, x_step = x_step, y_step = y_step, layout_shape = layout_shape)

  return(dag_out)
}


#' Generate graph coordinates and move stored curves with the nodes
#'
#' @param dag A dagitty object.
#' @param x_step Horizontal spacing between node columns.
#' @param y_step Minimum vertical spacing between nodes in a column.
#' @param layout_shape Vertical placement: `0` follows parent paths and `1`
#'   spreads each column evenly.
#' @return A dagitty object with generated node coordinates and transferred
#'   curve control points.
#' @noRd
layout_coords <- function(dag,
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
  if( nrow(dag_edges(dag, "add_coords()")) == 0L ){
    dag <- set_coords_preserving_curves(dag, list(
      x = stats::setNames(1 + (seq_along(node_names) - 1) * x_step, node_names),
      y = stats::setNames(rep(0, length(node_names)), node_names)
    ), "add_coords")
    return(dag)
  }

  new_coordinates <- renew_coords(dag,
                                    x_step = x_step,
                                    y_step = y_step,
                                    layout_shape = layout_shape)
  dag <- set_coords_preserving_curves(dag, new_coordinates, "add_coords")

  return(dag)
}


#' Create new dag node coordinates
#'
#' Without supplied coordinates, nodes are placed in columns according to
#' the longest directed path leading to them. Within each column, vertical
#' placement reduces crossings and keeps a minimum gap. With supplied
#' coordinates, only nodes missing either coordinate are placed on the
#' existing horizontal scale. All returned coordinates are rounded to
#' three decimal places.
#'
#' @param dag A dagitty object.
#' @param new_node_names Optional new node names. Placement also includes any
#'   other graph nodes missing coordinates.
#' @param coordinates Existing coordinates. If NULL or empty, all coordinates
#'   are regenerated.
#' @param x_step Horizontal spacing between node columns.
#' @param y_step Minimum vertical gap between nodes in the same column.
#' @param layout_shape Vertical placement of nodes in the same column, between
#'   0 and 1.
#' @return A list of named x and y vectors rounded to three decimal places.
#' @noRd
renew_coords <- function(dag,
                         new_node_names = NULL,
                         coordinates = NULL,
                         x_step = 2,
                         y_step = 1,
                         layout_shape = 0.5
                         ){

  edges <- data.table::as.data.table(dag_edges(dag, "renew_coords()"))
  pdag_edges <- edges[ ( edges$e == "--" | edges$e == "<->"), ]

  dag_node_names <- names(dag)
  parents <- layout_parents(dag)
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

  ## level symmetric edges across ranks only where the destination is free ##
  ## keep all existing node positions ##
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
        if( column_height_free(coordinates, lhs, level, node_rank, dag_node_names, y_step) &&
            column_height_free(coordinates, rhs, level, node_rank, dag_node_names, y_step) ){
          coordinates$y[lhs] <- level
          coordinates$y[rhs] <- level
        }
      }else if( lhs_new ){
        level <- coordinates$y[[rhs]]
        if( column_height_free(coordinates, lhs, level, node_rank, dag_node_names, y_step) ){
          coordinates$y[lhs] <- level
        }
      }else if( rhs_new ){
        level <- coordinates$y[[lhs]]
        if( column_height_free(coordinates, rhs, level, node_rank, dag_node_names, y_step) ){
          coordinates$y[rhs] <- level
        }
      }
    }
  }

  coordinates$y <- coordinates$y[!duplicated(names(coordinates$y))]
  coordinates$x <- coordinates$x[!duplicated(names(coordinates$x))]

  coordinates <- list(x = round( coordinates$x[order(names(coordinates$x))], 3), y = round( coordinates$y[order(names(coordinates$y))], 3) )

  return(coordinates)
}


#' Check a column height for the minimum vertical gap
#'
#' Tests whether a node can take a height without coming within `y_step` of
#' another positioned node in its column, that is a node of the same rank or
#' at the same horizontal position.
#'
#' @param coordinates list of x and y coordinates.
#' @param node Name of the node to move.
#' @param height Candidate vertical position.
#' @param node_rank Named topological ranks.
#' @param node_names Node names of the graph.
#' @param y_step Minimum vertical gap between nodes in the same column.
#' @return TRUE when the height is free.
#' @noRd
column_height_free <- function(coordinates, node, height, node_rank,
                                node_names, y_step){
  x <- coordinates$x
  y <- coordinates$y
  if( !node %in% names(x) ){ free <- TRUE
  return(free) }
  others <- setdiff( intersect( names(x), names(y) ), node )
  others <- others[ others %in% node_names & !is.na( x[others] ) & !is.na( y[others] ) ]
  others <- others[ node_rank[ others ] == node_rank[[ node ]] |
                    abs( x[ others ] - x[[ node ]] ) < 1e-9 ]

  free <- all( abs( y[ others ] - height ) >= y_step * 0.999 )

  return(free)
}


#' Helper for add_curves() and plot_dagitty()
#'
#' Rebuilds a pdag or mag with graph type dag for `ggdag::tidy_dagitty()`.
#' Nodes, edges, statuses, complete node coordinates and R attributes are
#' retained. Stored edge controls are handled separately by the caller.
#' A dag is returned unchanged.
#'
#' @param dag A dagitty object.
#' @param caller Function label for messages.
#' @return A dagitty object of graph type dag.
#' @noRd
dag_for_layout <- function(dag, caller){
  graph_type <- dag_graph_type(dag)
  if( graph_type == "dag" ){
    return(dag)
  }
  original_attributes <- attributes(dag)
  original_coordinates <- dagitty::coordinates(dag)
  dag <- construct_graph(
    dag_edges(dag, paste0(caller, "()")), names(dag),
    treatments(dag), get_outcome_names(dag), unobserved(dag),
    adjusted_vec = dagitty::adjustedNodes(dag),
    selected_vec = selected_nodes(dag),
    graph_type = "dag"
  )
  attributes(dag) <- original_attributes
  coordinate_names <- intersect(names(original_coordinates$x),
                                names(original_coordinates$y))
  positioned <- coordinate_names[
    !is.na(original_coordinates$x[coordinate_names]) &
      !is.na(original_coordinates$y[coordinate_names])
  ]
  if( length(positioned) > 0L ){
    dag <- set_coords_preserving_curves(dag, list(
      x = original_coordinates$x[positioned],
      y = original_coordinates$y[positioned]
    ), caller = caller)
  }
  return(dag)
}
