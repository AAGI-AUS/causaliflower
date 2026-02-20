
# plot_dagitty()
#' ggdag plot from a dagitty object
#'
#'Generates a ggdag plot from a dagitty object, adding coordinates and node labels using initials (capital letters).
#'Essentially a wrapper for dagitty and ggdag plotting functions with added features.
#'
#' @importFrom dagitty coordinates
#' @importFrom ggplot2 ggplot scale_shape_manual scale_size_manual scale_colour_manual labs guides guide_legend theme
#' @importFrom ggdag geom_dag_edges geom_dag_text geom_dag_label_repel geom_dag_point theme_dag
#' @importFrom ggraph circle
#' @param dag A dagitty object
#' @param coords_spec Adjusts node placement, a higher value increases volatility and results in more extreme DAG structures.
#' @param labels A vector of labels for nodes in dagitty object
#' @param label_type Defaults to using the node names of the inputted dag, or 'initials' if specified. Labels are only assigned by plot_dagitty() if nodes are not already labelled, using a 'ggdag::dag_label()' wrapper..
#' @param label_placement Labels are placed in a text box adjacent each node by default, or can be placed inside each node using 'label_node'. It is recommended to also specify label_type = 'initials' where label_placement = 'label_node' is used.
#' @return ggdag plot
#' @examples
#' plot_dagitty(dag)
#'
#' plot_dagitty(dag, coords_spec = 3) # increasing spec for more variation in node placement
#'
#' @export
plot_dagitty <- function(dag,
                         coords_spec =  0.1,
                         labels = NULL,
                         label_type = "name",
                         label_placement = "label_repel"){

  dag <- pdag_to_dag(dag)

  existing_coordinates <- dagitty::coordinates(dag)

  if( any( is.na( existing_coordinates$x ) ) | any( is.na( existing_coordinates$y ) ) ){

    dag <- add_coords(dag, coords_spec = coords_spec)

  }

  # aesthetics for ggdag legend
  col <- c("outcome"="deepskyblue",
           "treatment"="darkolivegreen3",
           "confounder"="coral2",
           "mediator"="darkorchid1",
           "mediator_outcome_confounder"="magenta4",
           "instrumental"="deeppink1",
           "proxy"="darkorange",
           "competing_exposure"="darkseagreen4",
           "collider"="darkred",
           "latent"="grey",
           "observed"="#222222")

  shape <- c("outcome"=19,
             "treatment"=19,
             "confounder"=19,
             "mediator"=19,
             "mediator_outcome_confounder"=19,
             "instrumental"=19,
             "proxy"=19,
             "collider"=19,
             "competing_exposure"=19,
             "latent"=1,
             "observed"=19)

  # variable for legend order
  order_col <- c("outcome",
                 "treatment",
                 "confounder",
                 "mediator",
                 "mediator_outcome_confounder",
                 "instrumental",
                 "proxy",
                 "competing_exposure",
                 "collider",
                 "latent",
                 "observed")

  if( is.null(labels) ){

    if(label_type == "name"){

      labels <- as.data.frame(unique(suppressWarnings(ggdag::dag_label(dag)$data["name"])))
      labels <- as.vector(labels$name)
      names(labels) <- labels

    }else if(label_type == "initials"){


      labels <- get_labels(dag)


    }else{

      stop("Invalid label_type input. Please use 'initials' or the default 'name'.")

    }

  }else if(!is.null(labels) & length(labels) != length(names(dag))){

    stop("The length of supplied labels does not equal the number of nodes in the graph. Please check labels and try again.")

  }


  dag_df <- tidy_ggdagitty(dag, labels)

  dag_df_complete_cases <- dag_df[complete.cases(dag_df[, "role"]), ]

  if( nrow(dag_df) > nrow(dag_df_complete_cases) ){
    message("Nodes without arrows removed and not displayed in the plotted graph.")
    dag_df <- dag_df_complete_cases
  }

  dag_df[sapply(dag_df, is.character)] <- lapply( dag_df[sapply(dag_df, is.character)], as.factor )


  if(label_placement == "label_repel"){


    ggdag <- ggplot2::ggplot(data = dag_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color=role, shape = role, fill = role)) +
      ggdag::geom_dag_edges(start_cap = ggraph::circle(5, 'mm'),
                            end_cap = ggraph::circle(5, 'mm')) +
      ggdag::geom_dag_point(size=10) +
      ggdag::geom_dag_label_repel(ggplot2::aes(label = label), colour = "white", alpha = 0.7, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = col,
                                   name = "Group",
                                   breaks = order_col) +
      ggplot2::scale_shape_manual(values = shape,
                                  name = "Group",
                                  breaks = order_col) +
      ggplot2::scale_fill_manual(values = col,
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

  }else if(label_type == "label_node"){


    ggdag <- ggplot2::ggplot(data = dag_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color=role, shape = role)) +
      ggdag::geom_dag_edges(start_cap = ggraph::circle(5, 'mm'),
                            end_cap = ggraph::circle(5, 'mm')) +
      ggdag::geom_dag_point(size=10) +
      ggdag::geom_dag_text(data = subset(dag_df, !duplicated(dag_df$label)), ggplot2::aes(label = label), colour = "grey10", size = 3) +
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

    stop("The label_placement input is invalid. Please use either 'label_node' or 'label_repel'.")
  }

}


#' Generate new dag coordinates
#'
#' @param dag Dagitty object, with at least a treatment and outcome set. Nodes are automatically placed in the following categories: treatment, outcome, confounder, mediator, latent variable, mediator-outcome-confounder, or instrumental variable.
#' @param coords_spec Adjusts node placement, a higher value increases volatility and results in more extreme DAG structures.
#' @param confounders Vector of confounders in order of occurrence.
#' @return dagitty object with coordinates.
#' @examples
#' dag <- add_coords(dag, coords_spec = 3) # update dagitty object node coordinates
#'
#' plot_dagitty(dag) # check coordinates
#'
#' @export
add_coords <- function(dag,
                       coords_spec = 0.1){

  num_nodes <- length(names(dag))
  time_limit <- num_nodes + num_nodes*coords_spec[1]

  setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)

  on.exit( {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  } )

  tryCatch({
    new_coordinates <- renew_coords(dag = dag,
                                    coords_spec = c( coords_spec[ complete.cases(coords_spec) ], threshold = 0.9 ) )

    dagitty::coordinates(dag) <- new_coordinates
    return(dag)
  }, warning = function(w){
    message(paste("Warning:", w, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords_helper(dag, coords_spec = unname( coords_spec[1][ complete.cases(coords_spec) ] ) )
    return(dag)

  }, error = function(e){
    message(paste("Error:", e, "\n Using alternative function to generate dag coordinates."))
    dag <- add_coords_helper(dag, coords_spec = unname( coords_spec[1][ complete.cases(coords_spec) ] ) )
    return(dag)

  }, finally = return(dag)

  )

}


#' generate coordinates for a dagitty object
#'
#' @importFrom dagitty coordinates parents
#' @importFrom ggdag time_ordered_coords
#' @param dag dagitty object
#' @param confounders vector of confounders nodes in the supplied dag.
#' @param mediators vector of mediator nodes in the supplied dag.
#' @param instrumental_variables vector of instrumental variable nodes in the supplied dag.
#' @param lambda adjusts sensitivity of node placement. Higher lambda introduces more volatility in confounder nodes and treatment along the y-axis.
#' @return dagitty objecty with coordinates.
#' @noRd
add_coordinates <- function(dag, treatments, outcomes, confounders, mediators, instrumental_variables, mediator_outcome_confounders, competing_exposures, latent_variables, colliders,
                            lambda = 1,
                            observed = NA
){

  o <- length(outcomes)

  t <- length(treatments)


  if( all( complete.cases( unlist(confounders) ) ) ){

    num_confounders <- c <- length(confounders)

  }else{

    confounders <- NA

    num_confounders <- c <- 0
  }

  if( all( complete.cases( unlist(instrumental_variables) ) ) ){

    i <- length(instrumental_variables)

  }else{

    instrumental_variables <- NA

    i <- 0
  }

  if( all( complete.cases( unlist(latent_variables) ) ) ){

    l <- length(latent_variables)

  }else{

    latent_variables <- NA

    l <- 0
  }

  if( all( complete.cases( mediators) ) ){

    m <- length(mediators)

  }else{

    m <- 0
  }

  if( all( complete.cases( mediator_outcome_confounders) ) ){

    moc <- length(mediator_outcome_confounders)

  }else{

    moc <- 0
  }

  if( all( complete.cases( competing_exposures) ) ){

    comp <- length(competing_exposures)

  }else{

    comp <- 0
  }

  if( all( complete.cases( colliders) ) ){

    coll <- length(colliders)

  }else{

    coll <- 0
  }


  ## variable order (used for generating coordinates)
  var_order <- c( unlist(confounders, recursive = FALSE), mediator_outcome_confounders, mediators, instrumental_variables, competing_exposures, colliders, treatments, outcomes)

  var_order <- var_order[ complete.cases(var_order) ]

  num_vars <- n <- length(unlist(var_order))


  ## lambda tuning parameter
  lambda <- lambda/(c + (n-i)) #print(paste("Coordinates added with randomization multiplier set to:", lambda/(num_confounders + (n-i))))


  ## get default ordered coordinates
  coords <- ggdag::time_ordered_coords(unlist(var_order), direction = c("x"), auto_sort_direction = c("right"))
  coords$y <- (ggdag::time_ordered_coords(unlist(var_order), direction = c("y"), auto_sort_direction = c("right")))$y


  ## randomize default coordinates
  rand_num_x <- runif(num_vars, min = -lambda, max = lambda)
  rand_num_y <- runif(num_vars, min = -lambda, max = lambda)

  coords$y <- coords$y + rand_num_y
  coords$x <- coords$x + rand_num_x

  ## initial coordinates ##

  ## confounders
  if( num_confounders > 0 ){
    rand_conf_x_y <- sort(rlnorm(num_confounders*2), decreasing = FALSE)*(lambda/(n-i))

    rand_conf_x <- rev(rand_conf_x_y)
    rand_conf_x <- rand_conf_x[1:(num_confounders+1)]
    rand_conf_y <- rand_conf_x_y[1:(num_confounders+1)]

    rand_conf_x <- c(0, seq(c*(lambda/3), (c)*(lambda), length.out = num_confounders))[1:(num_confounders)] + rand_conf_x[1:(num_confounders)]
    rand_conf_y <- c(0, seq(-(c)*(lambda/3), -c*(lambda^2), length.out = num_confounders/2),
                     seq(c*(lambda/3), c*(lambda), length.out = num_confounders/2))[1:(num_confounders)] + (rand_conf_y*(n-i))[1:(num_confounders)]

    conf_num_vars_modifier <- seq(0, num_vars, length.out = c)

    seq_conf <- seq(0, c, length.out = c)

    coords[1:num_confounders,]$x <- coords[1:num_confounders,]$x + conf_num_vars_modifier*0.6 + rand_conf_x[1:num_confounders]

    coords[1:num_confounders,]$y <- coords[1:num_confounders,]$y + conf_num_vars_modifier*0.2 + rand_conf_y[1:num_confounders] - c/num_vars*seq_conf

  }


  ## mediator-outcome-confounders
  if( all( complete.cases( mediator_outcome_confounders) ) ){

    a <- 1 + as.integer(!is.null((lambda/lambda)))
    b <- 1 + as.integer(is.null((lambda/lambda)))
    x <- lambda + (as.integer(is.null((lambda/lambda))) / 2)


    rand_moc_y <- seq( -( (((b*m^2)*(x^2))/(a*(n-i))) + ((b*m)*x)/a - (m/2) ),
                       -( (((m^2)*(x^2))/(n-i)) + (m)*x - (m/2) ),  length.out = moc)

    rand_moc_x <- seq( ( c/(moc + m)),
                       x*(m*x + m ),  length.out = moc)

    coords[(num_confounders+1):(num_confounders+moc),]$y <- coords[(num_confounders+1):(num_confounders+moc),]$y + rand_moc_y
    coords[(num_confounders+1):(num_confounders+moc),]$x <- coords[(num_confounders+1):(num_confounders+moc),]$x + rand_moc_x

  }


  ## treatments
  trt_y_max <- seq( 0 + length(t)/3,
                    1.1 + length(t)/n,  length.out = t )

  trt_x_min <- seq( ( ((num_confounders/num_vars)*(lambda*num_vars)) - (i/n) ) ,
                    ( ((num_confounders/num_vars)*(lambda*num_vars)) + (i/3) ) ,  length.out = t )

  coords[ ( nrow(coords) - t - o + 1 ):( nrow(coords) - o ), ]$x <- coords[ ( nrow(coords) - t - o + 1):(nrow(coords) - o ), ]$x + trt_x_min
  coords[ ( nrow(coords) - t - o + 1 ):( nrow(coords) - o ), ]$y <- trt_y_max


  ## outcomes
  outcome_x_min <- ((m-i)*(lambda)*2) + m - i
  outcome_y_min <- ((m-i)*(lambda)*2) + i

  coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$x <- coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$x + outcome_x_min
  coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$y <- coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$y - outcome_y_min


  ## mediators
  if( all( complete.cases( mediators) ) ){

    a <- 1 + as.integer(!is.na((lambda/lambda)))
    b <- 1 + as.integer(is.na((lambda/lambda)))
    x <- lambda + (as.integer(is.na((lambda/lambda))) / 2)

    rand_med_y <- seq( ( ( trt_y_max[1] +  (((b*m^2)*(x^2))/(a*(n-i))) + ((b*m)*x)/a + (m/2) ) ),
                       ( ( trt_y_max[t] +  (((m^2)*(x^2))/(n-i)) + (m)*x + (m/2) ) ),  length.out = m)

    rand_med_x <- seq( ( coords[ ( nrow(coords) -t-o+1 ), ]$x - m - (n*(m/n)) - (m/(4/b)) ),
                       ( coords[ (nrow(coords) -o ), ]$x - ( x*(m*x + m) + m ) ),  length.out = m)

    coords[(num_confounders+moc+1):(num_confounders+moc+m),]$y <- coords[(num_confounders+moc+1):(num_confounders+moc+m),]$y + rand_med_y
    coords[(num_confounders+moc+1):(num_confounders+moc+m),]$x <- rand_med_x

  }


  ## fine-tune treatment coordinates (based on parent nodes)
  grouped_nodes <- treatments
  len <- t

  group_parent_max_y_coord <- parent_node_max_coords(dag,
                                                     group_node_names = treatments,
                                                     num_group_nodes = t,
                                                     coordinates_x_or_y = coords$y,
                                                     coord_names = coords$name)

  # update treatment y-coords
  coords[ ( nrow(coords) - t - o + 1):(nrow(coords) - o ), ]$y <- coords[ ( nrow(coords) - t - o + 1):(nrow(coords) - o ), ]$y + group_parent_max_y_coord

  # update treatment x-coords
  coords[ ( nrow(coords) - t - o + 1):(nrow(coords) - o ), ]$x <- max( coords$x ) + ( max( coords$x ) -
                                                                                        coords[ ( nrow(coords) - t - o + 1):(nrow(coords) - o ), ]$x )


  ## fine-tune outcome coordinates (based on parent nodes)
  group_parent_max_y_coord <- parent_node_max_coords(dag,
                                                     group_node_names = outcomes,
                                                     num_group_nodes = o,
                                                     coordinates_x_or_y = coords$y,
                                                     coord_names = coords$name)

  # update outcome y-coords
  coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$y <- group_parent_max_y_coord + abs( ( group_parent_max_y_coord - coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$y ) )

  # update outcome x-coords
  coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$x <- max( coords$x ) + abs( ( max( coords$x ) - coords[ ( nrow(coords) - o + 1 ):( nrow(coords) ), ]$x ) )


  ## competing exposures
  if( all( complete.cases( competing_exposures) ) ){

    trt_updat_x_max <- max( coords[ ( nrow(coords) - t - o + 1):(nrow(coords) - o ), ]$x )


    comp_y <- seq( ( coords[nrow(coords),]$y - (comp*(3/comp)) ) , (coords[nrow(coords),]$y - (comp*(2/comp))*(c/n) ), length.out = comp)


    comp_x <- seq( ( trt_updat_x_max + ( max( coords$x ) / num_vars ) + 1 ), ( trt_updat_x_max + ( max( coords$x ) / num_vars ) + comp + 1 ), length.out = comp)

    coords[(num_confounders+moc+m+i+1):(num_confounders+moc+m+i+comp),]$y <- comp_y
    coords[(num_confounders+moc+m+i+1):(num_confounders+moc+m+i+comp),]$x <- comp_x

  }


  ## colliders
  if( all( complete.cases( colliders) ) ){

    coll_y <- seq( ( coords[nrow(coords),]$y + ( max( coords$y ) / num_vars )*lambda + 1 ), (coords[nrow(coords),]$y + ( max( coords$y ) / num_vars )*lambda + coll ), length.out = coll )

    coll_x <- seq( ( coords[nrow(coords),]$x + coll*lambda ), ( coords[nrow(coords),]$x + ( coll*(coll/2)*lambda ) ), length.out = coll)

    coords[(num_confounders+moc+m+i+comp+1):(num_confounders+moc+m+i+comp+coll),]$y <- coll_y
    coords[(num_confounders+moc+m+i+comp+1):(num_confounders+moc+m+i+comp+coll),]$x <- coll_x

  }


  ## instrumental variables
  if( all( complete.cases( unlist(instrumental_variables) ) ) ){

    instr_y <- seq( ( coords[nrow(coords)-t-o+1,]$y - ( 1 + (i/n) )*1.5 ) , ( coords[nrow(coords)-o-1,]$y - i ), length.out = i) # fix these coords

    instr_x <- seq( ( coords[nrow(coords)-t-o+1,]$x + i*(c/n) ) , ( coords[nrow(coords)-o-1,]$x + i ),  length.out = i)

    coords[(num_confounders+moc+m+1):(num_confounders+moc+m+i),]$y <- instr_y
    coords[(num_confounders+moc+m+1):(num_confounders+moc+m+i),]$x <- instr_x

  }


  ## final coordinates ##
  coords_x <- coords$x
  coords_y <- coords$y

  names(coords_x) <- coords$name
  names(coords_y) <- coords$name

  coords_list <- list(
    x = coords_x,
    y = coords_y
  )

  # save existing coords
  existing_coords <- dagitty::coordinates(dag)

  # check if no existing coords and no coords for nodes with observed role
  if( all( is.na( unlist(existing_coords) ) ) & all( !is.na(observed) ) ){
    # generate coords for observed nodes
    coords_list <- observed_new_coordinates_helper(dag, observed, coords_list)

  }

  # generate coords for latent nodes
  coords_list <- latent_new_coordinates_helper(dag, latent_variables, coords_list, lambda)

  dagitty::coordinates(dag) <- coords_list

  return(dag)
}


#' Create new dag node coordinates
#'
#' @importFrom ggdag time_ordered_coords
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_node_names Vector of new node names. When new nodes are specified, existing nodes can be used as a reference point. Otherwise, all new coordinates are generated using ggdag::automatically generates
#' @param coordinates list of coordinates from a dagitty object.
#' @param coords_spec Parameters used for generating coordinates. Adjust node placement with lambda; a higher value increases volatility and results in more extreme DAG structures. Threshold controls the closeness of nodes.
#' @return Dagitty coordinates.
#' @examples
#' coordinates <- renew_coords(dag, coords_spec = c(0.1, 0.5)) # generate new coordinates
#'
#' dagitty::coordinates(dag) <- coordinates # add coordinates to dagitty object
#'
#' @noRd
renew_coords <- function(dag,
                         new_node_names = NULL,
                         coordinates = NA,
                         coords_spec = c(lambda = 0.1, threshold = 0.9)
){

  new_node_names <- as.vector( unlist( new_node_names ) ) # new node names as vector
  num_nodes <- length( new_node_names ) # number of new nodes
  dag_node_names <- names(dag)
  num_vars <- length( dag_node_names ) # total number of nodes (combined dag)

  lambda <- coords_spec[1]/(num_nodes + (num_vars)) # calc lambda value (controls volatility of generated coordinates)
  threshold <- coords_spec[2]

  if( num_nodes < 1 ){

    new_node_names <- dag_node_names

  }

  node_roles <- get_roles(dag)

  confounders <- node_roles$confounder
  mediators <- node_roles$mediator
  mediator_outcome_confounders <- node_roles$mediator_outcome_confounder
  instrumental_variables <- node_roles$instrumental
  proxy_variables <- node_roles$proxy
  competing_exposures <- node_roles$competing_exposure
  colliders <- node_roles$collider
  latent_variables <- node_roles$latent
  observed <- node_roles$observed
  outcomes <- node_roles$outcome
  treatments <- node_roles$treatment

  keep_existing_coords <- TRUE

  existing_coordinates <- coordinates

  if( all( !complete.cases( unlist(coordinates) ) ) ){

    keep_existing_coords <- FALSE

    coordinates <- dagitty::coordinates(dag)

    coordinates <- list( x = coordinates$x[ complete.cases( unlist(coordinates$x) ) ],
                         y = coordinates$y[ complete.cases( unlist(coordinates$y) ) ] )

    if( all( !complete.cases( unlist(coordinates) ) ) ){

      var_order <- c(confounders, instrumental_variables, mediator_outcome_confounders,
                     treatments, observed, competing_exposures,
                     proxy_variables, latent_variables, mediators,
                     outcomes, colliders)

      var_order <- var_order[complete.cases(var_order)]

      ## get default ordered coordinates
      coordinates <- ggdag::time_ordered_coords(unlist(var_order), direction = c("x"), auto_sort_direction = c("right"))

      coordinates_x <- coordinates$x
      names(coordinates_x) <- coordinates$name

      coordinates_y <- (ggdag::time_ordered_coords(unlist(var_order), direction = c("y"), auto_sort_direction = c("right")))$y
      names(coordinates_y) <- coordinates$name

      coordinates <- list(x = coordinates_x, y = coordinates_y)

    }

    existing_coordinates <- coordinates

  }

  ordered_nodes <- unlist( dagitty::topologicalOrdering(dag) )
  ordered_nodes <- order(ordered_nodes)
  new_node_names <- new_node_names[ordered_nodes]

  if( all( complete.cases(confounders) ) ){

    ## confounders

    new_node_name_vec <- new_node_names[ new_node_names %in% confounders ]
    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates (without treatments/outcomes)
      coordinates_x <- coordinates$x[ !names(coordinates$x) %in% c(treatments, outcomes)]
      coordinates_y <- coordinates$y[ !names(coordinates$y) %in% c(treatments, outcomes)]

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 treatments = treatments,
                                                                 outcomes = outcomes,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )

    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% confounders ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% confounders ],
              coordinates$y)
      )

    }


  }

  if( all( complete.cases(treatments) ) ){

    ## treatments
    new_node_name_vec <- new_node_names[ new_node_names %in% treatments ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( treatments, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates (without treatments/outcomes)
      coordinates_x <- coordinates$x[ !names(coordinates$x) %in% c(outcomes)]
      coordinates_y <- coordinates$y[ !names(coordinates$y) %in% c(outcomes)]

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 outcomes = outcomes,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )

    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% treatments ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% treatments ],
              coordinates$y)
      )

    }

  }

  if(  all( complete.cases(instrumental_variables) ) ){

    ## instrumental_variables
    new_node_name_vec <- new_node_names[ new_node_names %in% instrumental_variables ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( instrumental_variables, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates (without treatments/outcomes)
      coordinates_x <- coordinates$x[ !names(coordinates$x) %in% c(outcomes)]
      coordinates_y <- coordinates$y[ !names(coordinates$y) %in% c(outcomes)]

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 outcomes = outcomes,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )


    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% instrumental_variables ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% instrumental_variables ],
              coordinates$y)
      )

    }

  }

  if( all( complete.cases(mediator_outcome_confounders) ) ){

    ## mediator_outcome_confounders
    new_node_name_vec <- new_node_names[ new_node_names %in% mediator_outcome_confounders ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( mediator_outcome_confounders, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates (without outcomes)
      coordinates_x <- coordinates$x[ !names(coordinates$x) %in% c(outcomes)]
      coordinates_y <- coordinates$y[ !names(coordinates$y) %in% c(outcomes)]

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 outcomes = outcomes,
                                                                 post_treatment = TRUE,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )

    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% mediator_outcome_confounders ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% mediator_outcome_confounders ],
              coordinates$y)
      )

    }

  }

  if( all( complete.cases(proxy_variables) ) ){

    ## proxy_variables
    new_node_name_vec <- new_node_names[ new_node_names %in% proxy_variables ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( proxy_variables, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates (without outcomes)
      coordinates_x <- coordinates$x[ !names(coordinates$x) %in% c(outcomes)]
      coordinates_y <- coordinates$y[ !names(coordinates$y) %in% c(outcomes)]

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 outcomes = outcomes,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )


    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% proxy_variables ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% proxy_variables ],
              coordinates$y)
      )

    }

  }

  if( all( complete.cases(mediators) ) ){

    ## mediators
    new_node_name_vec <- new_node_names[ new_node_names %in% mediators ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( mediators, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates (without outcomes)
      coordinates_x <- coordinates$x[ !names(coordinates$x) %in% c(outcomes)]
      coordinates_y <- coordinates$y[ !names(coordinates$y) %in% c(outcomes)]

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 outcomes = outcomes,
                                                                 post_treatment = TRUE,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )

    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% mediators ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% mediators ],
              coordinates$y)
      )

    }

  }

  if( all( complete.cases(competing_exposures) ) ){

    ## competing_exposures
    new_node_name_vec <- new_node_names[ new_node_names %in% competing_exposures ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( competing_exposures, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates$x,
                                                                 coordinates_y = coordinates$y,
                                                                 outcomes = outcomes,
                                                                 post_outcome = FALSE,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )


    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% competing_exposures ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% competing_exposures ],
              coordinates$y)
      )

    }

  }


  if( all( complete.cases(outcomes) ) ){

    ## outcomes
    new_node_name_vec <- new_node_names[ new_node_names %in% outcomes ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( outcomes, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates$x,
                                                                 coordinates_y = coordinates$y,
                                                                 outcomes = outcomes,
                                                                 post_treatment = TRUE,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )

    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% outcomes ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% outcomes ],
              coordinates$y)
      )

    }

  }

  if( all( complete.cases(colliders) ) ){

    ## colliders
    new_node_name_vec <- new_node_names[ new_node_names %in% colliders ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( colliders, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates$x,
                                                                 coordinates_y = coordinates$y,
                                                                 outcomes = outcomes,
                                                                 post_outcome = TRUE,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )

    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% colliders ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% colliders ],
              coordinates$y)
      )

    }

  }


  if( all( complete.cases(observed) ) ){

    ## observed
    new_node_name_vec <- new_node_names[ new_node_names %in% observed ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( observed, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates
      coordinates_x <- coordinates$x
      coordinates_y <- coordinates$y

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 post_outcome = TRUE,
                                                                 outcomes = outcomes,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )

    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% observed ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% observed ],
              coordinates$y)
      )

    }


  if( all( complete.cases(latent_variables) ) ){

    ## latent_variables
    new_node_name_vec <- new_node_names[ new_node_names %in% latent_variables ]

    if(  keep_existing_coords == FALSE ){

      new_node_name_vec <- c( latent_variables, new_node_name_vec )
      new_node_name_vec <- new_node_name_vec[ !duplicated(new_node_name_vec) ]

    }

    num_nodes <- length(new_node_name_vec)

    if( num_nodes > 0 ){

      # new node coordinates
      coordinates_x <- coordinates$x
      coordinates_y <- coordinates$y

      coordinates <- suppressWarnings( merged_node_coords_helper(dag,
                                                                 new_node_name_vec = new_node_name_vec,
                                                                 num_nodes = num_nodes,
                                                                 coordinates_x = coordinates_x,
                                                                 coordinates_y = coordinates_y,
                                                                 post_outcome = TRUE,
                                                                 outcomes = outcomes,
                                                                 num_vars = num_vars,
                                                                 lambda = lambda,
                                                                 threshold = threshold)
      )


    }else{

      coordinates <- list(
        x = c(existing_coordinates$x[ names(existing_coordinates$x) %in% latent_variables ],
              coordinates$x),
        y = c(existing_coordinates$y[ names(existing_coordinates$y) %in% latent_variables ],
              coordinates$y)
      )

    }

  }

  }

  coordinates <- list(x = round( coordinates$x[order(names(coordinates$x))], 3), y = round( coordinates$y[order(names(coordinates$y))], 3) )

  return(coordinates)

}

