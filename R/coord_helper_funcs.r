
#' Add coordinates to nodes in a merged dagitty object
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_node_name_vec Inputted vector of new node names added to the graph.
#' @param num_nodes Number of new nodes.
#' @param nodes_children Descendants of new nodes.
#' @param nodes_parents Ancestors of new nodes.
#' @param coordinates_x Existing x-coordinates.
#' @param coordinates_y Existing y-coordinates.
#' @param treatments Treatment/exposures in the supplied dag.
#' @param outcomes Outcomes in the supplied dag.
#' @param post_treatment Indicates whether the coordinates of new nodes should be informed by treatment coordinates.
#' @param post_outcome Similarly, outcome coordinates can be used to inform new node coordinates.
#' @param num_vars Total number of variables in the dag.
#' @param lambda Parameter used to generate coordinates. Adjust node placement with lambda; a higher value increases volatility and results in more extreme DAG structures.
#' @param lambda Parameter used to generate coordinates. Threshold controls the closeness of nodes.
#' @return Named list of new coordinates.
#' @noRd
merged_node_coords_helper <- function(dag,
                                      new_node_name_vec,
                                      num_nodes,
                                      nodes_children,
                                      nodes_parents,
                                      coordinates_x,
                                      coordinates_y,
                                      treatments = NA,
                                      outcomes,
                                      post_treatment = FALSE,
                                      post_outcome = FALSE,
                                      num_vars,
                                      lambda,
                                      threshold
                                      ){
  .datatable.aware <- TRUE
  existing_node_names <- names(coordinates_x)
  existing_node_names_not_outcome <- existing_node_names[ !existing_node_names %in% outcomes ]

  nodes_parents <- lapply(1:num_nodes, function(x){   # get new node name parents
    nodes_parents <- dagitty::parents(dag, new_node_name_vec[x])
  })

  nodes_children <- lapply( 1:num_nodes, function(x){  # get new node name children
    nodes_children <- dagitty::children(dag, new_node_name_vec[x])
  })

  if( post_treatment == FALSE & all( complete.cases(treatments) ) ){

    existing_node_names_not_outcome <- existing_node_names_not_outcome[ !existing_node_names_not_outcome %in% treatments ]

  }

  if( post_outcome == FALSE ){

    nodes_parents <- lapply(1:num_nodes, function(x){ # get new node name parents in existing nodes (found in both dags prior to merge)
      nodes_parents <- existing_node_names_not_outcome[ existing_node_names_not_outcome %in% nodes_parents[[x]] ]
    })

    nodes_children <- lapply(1:num_nodes, function(x){ # get new node name parents in existing nodes (found in both dags prior to merge)
      nodes_children <- existing_node_names_not_outcome[ existing_node_names_not_outcome %in% nodes_children[[x]] ]
    })

  }

  ## separate x and y coordinates, node names, num nodes
  new_node_name_vec_y <- new_node_name_vec
  num_nodes_y <- num_nodes
  nodes_parents_y <- nodes_parents
  nodes_children_y <- nodes_children
  new_node_name_vec_x <- new_node_name_vec

  num_nodes_x <- num_nodes
  nodes_parents_x <- nodes_parents
  nodes_children_x <- nodes_children
  all_node_names_vec <- c(names(coordinates_x), new_node_name_vec)

  ## get difference between min-max coords
  diff_y_coords <-  abs( min( coordinates_y ) - max( coordinates_y ) )
  diff_x_coords <-  abs( min( coordinates_x ) - max( coordinates_x ) )

  quality_check <- FALSE
  iteration <- 1

  while(quality_check == FALSE){

    iteration <- iteration + ( num_nodes*lambda )

    if( num_nodes_y > 0){
      ## y coordinates ##

      new_y_coords <- sapply(1:num_nodes_y, function(x){

        if( length( nodes_parents_y[[x]] ) > 0 ){


          new_y_coords <- max( coordinates_y[ names(coordinates_y) %in% nodes_parents_y[[x]] ]
          ) + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                       min = iteration*lambda,
                                                                                       max = iteration + ( num_nodes*lambda)/2 )

          }else if( length( nodes_children_y[[x]] ) > 0 ) {

          new_y_coords <- min( coordinates_y[ names(coordinates_y) %in% nodes_children_y[[x]] ]
          ) - x*iteration - ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) - runif(n = 1,
                                                                             min = iteration*lambda,
                                                                             max = iteration + ( num_nodes*lambda )/2 )

          }else if( length( unlist(nodes_parents_y) ) > 0 ){

          new_y_coords <- min( coordinates_y[ names(coordinates_y) %in% unlist(nodes_parents_y) ]
          ) + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                             min = iteration*lambda,
                                                                             max = iteration + ( num_nodes*lambda )/2 )

          }else if( length( unlist(nodes_children_y) ) > 0 ){

          new_y_coords <- min( coordinates_y[ names(coordinates_y) %in% unlist(nodes_children_y) ]
          ) - x*iteration - ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) - runif(n = 1,
                                                                             min = iteration*lambda,
                                                                             max = iteration + ( num_nodes*lambda )/2 )

          }else{

          new_y_coords <- x + x*iteration + ( diff_y_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                                             min = iteration*lambda,
                                                                                             max = iteration + ( num_nodes*lambda )/2 )
        }

      })

      names(new_y_coords) <- new_node_name_vec_y

    }

    if( num_nodes_x > 0){
      ## x coordinates ##
      new_x_coords <- sapply(1:num_nodes_x, function(x){

        if( length( nodes_children_x[[x]] ) > 0 ) {

          new_x_coords <- min( coordinates_x[ names(coordinates_x) %in% nodes_children_x[[x]] ]
          ) - x*iteration - ( diff_x_coords/num_vars )*num_nodes + diff_x_coords*(x/diff_x_coords) - runif(n = 1,
                                                                                               min = iteration*lambda,
                                                                                               max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( nodes_parents_x[[x]] ) > 0 ){

          new_x_coords <- min( coordinates_x[ names(coordinates_x) %in% nodes_parents_x[[x]] ]
          ) + x*iteration + ( diff_x_coords/num_vars )*num_nodes*(x/num_nodes) + runif(n = 1,
                                                                             min = iteration*lambda,
                                                                             max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( unlist(nodes_children_x) ) > 0 ){

          new_x_coords <- min( coordinates_x[ names(coordinates_x) %in% unlist(nodes_children_x) ]
          ) - x*iteration - ( diff_x_coords/num_vars )*num_nodes + diff_x_coords*(x/diff_x_coords) - runif(n = 1,
                                                                                               min = iteration*lambda,
                                                                                               max = iteration + ( num_nodes*lambda*10 )/2 )

        }else if( length( unlist(nodes_parents_x) ) > 0 ){

          new_x_coords <- max( coordinates_x[ names(coordinates_x) %in% unlist(nodes_parents_x) ]
          ) + x*iteration + diff_x_coords*(x/diff_x_coords) + runif(n = 1,
                                                      min = iteration*lambda,
                                                      max = iteration + ( num_nodes*lambda*10 )/2 )
        }else{

          new_x_coords <- x + x*iteration + diff_x_coords*(x/diff_x_coords) - runif(n = 1,
                                                                      min = iteration*lambda,
                                                                      max = iteration + ( num_nodes*lambda*10 )/2 )
        }
        new_x_coords

      })

      names(new_x_coords) <- new_node_name_vec_x
    }

    new_coordinates <- quality_check_coords_2(new_node_name_vec_x = new_node_name_vec_x,
                                              new_node_name_vec_y = new_node_name_vec_y,
                                              num_nodes_x = num_nodes_x,
                                              num_nodes_y = num_nodes_y,
                                              new_x_coords = new_x_coords,
                                              new_y_coords = new_y_coords,
                                              coordinates_x = coordinates_x,
                                              coordinates_y = coordinates_y,
                                              threshold = threshold)

    ## initialise new variables ##
    # y-coords
    coordinates_y <-  new_coordinates$y

    nodes_children_y <- nodes_children_y[ !new_node_name_vec_y %in% names(coordinates_y) ]
    nodes_parents_y <- nodes_parents_y[ !new_node_name_vec_y %in% names(coordinates_y) ]

    new_node_name_vec_y <- new_node_name_vec_y[ !new_node_name_vec_y %in% names(coordinates_y)]
    num_nodes_y <- length(new_node_name_vec_y)

    # x-coords
    coordinates_x <- new_coordinates$x

    nodes_children_x <- nodes_children_x[ !new_node_name_vec_x %in% names(coordinates_x) ]
    nodes_parents_x <- nodes_parents_x[ !new_node_name_vec_x %in% names(coordinates_x) ]

    new_node_name_vec_x <- new_node_name_vec_x[ !new_node_name_vec_x %in% names(coordinates_x) ]
    num_nodes_x <- length(new_node_name_vec_x)

    if( num_nodes_y + num_nodes_x == 0 ){

      quality_check <- TRUE

    }

  }

  coordinates <- list(x = coordinates_x[!duplicated(coordinates_x)], y = coordinates_y[!duplicated(coordinates_y)])

  return(coordinates)

}


#' Checks node coordinates are not overlapping
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @return dagitty objecty with coordinates.
#' @noRd
quality_check_coords_2 <- function(new_node_name_vec_x, new_node_name_vec_y, num_nodes_x, num_nodes_y, new_x_coords, new_y_coords, coordinates_x, coordinates_y, threshold){

  # remove group node names from existing coords vectors
  existing_coordinates_x <- coordinates_x[ !names(coordinates_x) %in% new_node_name_vec_x]
  existing_coordinates_y <- coordinates_y[ !names(coordinates_y) %in% new_node_name_vec_y]

  num_vars_x <- length(existing_coordinates_x)
  num_vars_y <- length(existing_coordinates_y)


  ## internal check against other generated coords ##

  # x-coords
  if( num_nodes_x > 0 ){

    coords_df <- data.table::as.data.table(

      sapply(1:num_nodes_x, function(a){
      sapply(1:num_nodes_x, function(b){

        sqrt( diff( range( c(new_x_coords[a],
                             new_x_coords[b]) ) ) )**2

      })
    }) + diag(nrow = num_nodes_x, ncol = num_nodes_x) )


    coords_df$ID <- 1:num_nodes_x
    coords_df <- coords_df[apply(coords_df >= threshold, 1, all)]

    keep_coord_names_x <- new_node_name_vec_x[ coords_df$ID ]

  }

  # y-coords
  if( num_nodes_y > 0 ){

    coords_df <-  data.table::as.data.table(

      sapply(1:num_nodes_y, function(a){
        sapply(1:num_nodes_y, function(b){
          sqrt( diff( range( c(new_y_coords[a],
                               new_y_coords[b]) ) ) )**2

      })
    }) + diag(nrow = num_nodes_y, ncol = num_nodes_y) )

    coords_df$ID <- 1:num_nodes_y
    coords_df <- coords_df[apply(coords_df >= threshold, 1, all)]

    keep_coord_names_y <- new_node_name_vec_y[ coords_df$ID ]

  }

  ## check against existing coordinates ##
  if( num_nodes_x > 0 ){

    coords_df <- data.table::as.data.table(

      sapply(1:num_vars_x, function(a){
        sapply(1:num_nodes_x, function(b){
          sqrt( diff( range( c(existing_coordinates_x[a],
                               new_x_coords[b]) ) ) )**2

        })
      }) )

    if( num_nodes_x == 1 ){

      new_x_coords <- new_x_coords[ sum(coords_df) / length( existing_coordinates_x ) >= threshold ]

    }else{

      coords_df$ID <- 1:num_nodes_x
      coords_df <- coords_df[apply(coords_df >= threshold, 1, all)]

      new_x_coords <- new_x_coords[ coords_df$ID ]

      new_x_coords <- new_x_coords[ names(new_x_coords) %in% keep_coord_names_x]

    }


  }

  if( num_nodes_y > 0 ){

    coords_df <- data.table::as.data.table(

      sapply(1:num_vars_y, function(a){
        sapply(1:num_nodes_y, function(b){
          sqrt( diff( range( c(existing_coordinates_y[a],
                               new_y_coords[b]) ) ) )**2

        })
      }) )

    if( num_nodes_y == 1 ){

      new_y_coords <- new_y_coords[ sum(coords_df) / length( existing_coordinates_y ) >= threshold ]

    }else{

      coords_df$ID <- 1:num_nodes_y
      coords_df <- coords_df[apply(coords_df >= threshold, 1, all)]

      new_y_coords <- new_y_coords[ coords_df$ID ]

      new_y_coords <- new_y_coords[ names(new_y_coords) %in% keep_coord_names_y]

    }

  }

  new_coordinates <- list( x = c(existing_coordinates_x, new_x_coords),  y = c(existing_coordinates_y, new_y_coords) )

  return(new_coordinates)

}

#' Returns highest x or y coordinates for a group of nodes
#'
#' @importFrom dagitty parents
#' @param dag A dagittyy object.
#' @param group_node_names Node group of interest.
#' @param num_group_nodes Number of nodes in group of interest.
#' @param coordinates_x_or_y Named vector of either x or y coordinates.
#' @param coord_names Names corrresponding to coordinates_x_or_y input.
#' @return Maximum parent coordinates for a group of nodes.
#' @noRd
parent_node_max_coords <- function(dag, group_node_names, num_group_nodes, coordinates_x_or_y, coord_names){

  group_parents <- lapply(1:num_group_nodes, function(x){   # get new node name parents
    group_parents <- dagitty::parents(dag, group_node_names[x])
  })

  group_parent_max_coord <- max( coordinates_x_or_y[ coord_names # highest x or y coordinates for a group of nodes
                                             %in% unlist(group_parents) ], na.rm = TRUE )

  return(group_parent_max_coord)

}

