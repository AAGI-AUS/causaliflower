
# ggdagitty()
#' ggdag plot from a dagitty object
#'
#'Generates a ggdag plot from a dagitty object, adding coordinates and node labels using initials (capital letters).
#'Essentially a wrapper for dagitty and ggdag plotting functions with added features.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom dagitty coordinates
#' @importFrom ggplot2 ggplot scale_shape_manual scale_size_manual scale_colour_manual labs guides guide_legend theme
#' @importFrom ggdag geom_dag_edges geom_dag_text geom_dag_label_repel geom_dag_point theme_dag
#' @importFrom ggraph circle
#' @param dag A dagitty object
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures. Setting 'lambda_max' generates a DAG at each lambda value between lambda and lambda_max (only used if iterations is supplied). Iterations controls number of repeats for each lambda value (returns the first lambda value if NULL).
#' @param labels A vector of labels for nodes in dagitty object
#' @param label_type Defaults to using the node names of the inputted dag, or 'initials' if specified. Labels are only assigned by ggdagitty() if nodes are not already labelled, using a 'ggdag::dag_label()' wrapper..
#' @param label_placement Labels are placed in a text box adjacent each node by default, or can be placed inside each node using 'label_node'. It is recommended to also specify label_type = 'initials' where label_placement = 'label_node' is used.
#' @return ggdag plot
#' @importFrom magrittr %>%
#' @export
ggdagitty <- function(dag,
                      coords_spec = c(lambda = 0, lambda_max = NA, iterations = NA),
                      labels = NULL, label_type = "name", label_placement = "label_repel"){


  if(any(is.na(dagitty::coordinates(dag)))){

    dag <- getCoords(dag, coords_spec = coords_spec)

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
  order_col <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent")
  if( is.null(labels) ){

    if(label_type == "name"){

      if(is.null(labels)){
        labels <- as.data.frame(unique(suppressWarnings(ggdag::dag_label(dag)$data["name"])))
        labels <- as.vector(labels$name)
        names(labels) <- labels
      }
    }else if(label_type == "initials"){

      if(is.null(labels)){
        labels <- getLabels(dag)
      }

    }else{

      stop("Invalid label_type input. Please use 'initials' or the default 'name'.")

    }

  }else if(!is.null(labels) & length(labels) != length(names(dag))){

    stop("The length of supplied labels does not equal the number of nodes in the graph. Please check labels and try again.")

  }


  dag_df <- tidy_ggdagitty(dag, labels)

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

  #return(ggdag)
}


#' Label DAG nodes
#' Extract node initials to use as labels in ggdagitty, or dag plots using ggdag and ggplot2. This is essentially a wrapper function for ggdag::dag_label().
#' @importFrom ggdag dag_label
#' @param dag A dagittyy object.
#' @return Labelled vector of node names from the dagitty object
#' @export
getLabels <- function(dag){

  names <- as.data.frame(unique(suppressWarnings(ggdag::dag_label(dag)$data["name"])))
  names <- as.vector(names$name)
  labels <- gsub("(\\b[A-Z])[^A-Z]+", "\\1", names)
  names(labels) <- names

  return(labels)
}


#' Add coordinates to a dagitty object
#'
#' @importFrom magrittr %>%
#' @importFrom ggdag time_ordered_coords
#' @importFrom dplyr filter select
#' @param dag Dagitty object, with at least a treatment and outcome set. Nodes are automatically placed in the following categories: treatment, outcome, confounder, mediator, latent variable, mediator-outcome-confounder, or instrumental variable.
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures. Setting 'lambda_max' generates a DAG at each lambda value between lambda and lambda_max (only used if iterations is supplied). Iterations controls number of repeats for each lambda value (returns the first lambda value if NULL).
#' @param confounders Vector of confounders in order of occurrence. If left blank, the order is inferred with dagitty::topologicalOrdering().
#' @return dagitty object with coordinates.
#' @export
getCoords <- function(dag,
                      coords_spec = c(lambda = 0, lambda_max = NA, iterations = NA),
                      confounders = NULL){

  edges <- getFeatureMap(dag)

  if(is.null(confounders)){

    confounders <-  as.vector( unlist( lapply(edges$confounder,
                                              function(x) if(identical(x, character(0))) NA_character_ else x) ) )

      as.vector( unlist(edges$confounder) )

  }

  mediators <-  as.vector( unlist( lapply(edges$mediator,
                                          function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  instrumental_vars <- as.vector( unlist( lapply(edges$instrumental,
                                                 function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  mediator_outcome_confounders <- as.vector( unlist( lapply(edges$mediator_outcome_confounder,
                                                            function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  competing_exposures <-  as.vector( unlist( lapply(edges$competing_exposure,
                                                    function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  colliders <- as.vector( unlist( lapply(edges$collider,
                                         function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  latent_vars <- as.vector( unlist( lapply(edges$latent,
                                           function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  observed <- as.vector( unlist( lapply(edges$observed,
                                        function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  treatments <- as.vector( unlist( lapply(edges$treatment,
                                          function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  outcomes <- as.vector( unlist( lapply(edges$outcome,
                              function(x) if(identical(x, character(0))) NA_character_ else x) ) )

  if(length(na.omit(coords_spec)) == 3){
    message("Generating coordinates for ", coords_spec["iterations"], " DAGs for each lambda value between ", coords_spec["lambda"], " and ", coords_spec["lambda_max"], ".")
    dag <- auto_coords(dag, confounders, mediators, instrumental_vars, mediator_outcome_confounders, coords_spec)

    return(dag)

  }else if(length(na.omit(coords_spec)) == 1){

    tryCatch(
      { # try

        dag <- add_coordinates(dag,
                               treatments = treatments,
                               outcomes = outcomes,
                               confounders = confounders,
                               mediators = mediators,
                               instrumental_vars = instrumental_vars,
                               mediator_outcome_confounders = mediator_outcome_confounders,
                               competing_exposures = competing_exposures,
                               latent_vars = latent_vars,
                               observed = observed,
                               colliders = colliders,
                               na.omit(coords_spec))

      },
      # error = function(cond) {
      #  message("An unexpected error occurred:")
      #  message(conditionMessage(cond))
      #  NA
      # },
      #warning = function(cond) {
      #  message("A warning occurred:")
      #  message(conditionMessage(cond))
      #},
      finally = {
        return(dag)
      }
    )
  }
}


#' Generate multiple versions of DAG coordinates
#'
#' @importFrom magrittr %>%
#' @importFrom gridExtra grid.arrange
#' @importFrom dagitty dagitty
#' @param dag dagitty object
#' @param confounders Vector of confounder variables, e.g. c("Z1", "Z2", "Z3").
#' @param mediators Vector of mediator variables, e.g. c("M1", "M2", "M3").
#' @param instrumental_vars vector of instrumental variable nodes in the supplied dag.
#' @param coords_spec Vector containing 'iterations' for specifying number of training epochs and 'lambda_range' to control the volatility of random coordinates.
#' @return dag objecty with coordinates.
#' @noRd
auto_coords <- function(dag,
                          confounders,
                          mediators = NULL,
                          instrumental_vars = NULL,
                          coords_spec = c(lambda = 0, lambda_max = 10, iterations = NULL)){

  labels <- getLabels(dag)

  lambda_min <- coords_spec[1]
  lambda_max <- coords_spec[2]
  iterations <- coords_spec[3]

  if(lambda_min == 0){
    lambda_range <- lambda_max + 1
    l_increment <- 1
  }else if(is.integer(lambda_min)){
    lambda_range <- lambda_max - lambda_min
    l_increment <- 1
  }else if(is.double(lambda_min)){
    lambda_range <- (lambda_max*10) - (lambda_min*10)
    l_increment <- 0.1
  }else{
    stop("Please enter a valid lambda number.")
  }


  lambda <- lambda_min
  l <- 1
  dagitty_list <- list()
  dag_list <- list()
  for(l in 1:lambda_range){
    epoch <- 1
    dags <- list()
    temp_list <- list()
    for(epoch in 1:iterations){
      dag <- add_coordinates(dag, confounders, mediators, lambda)
      dags[[epoch]] <- dag
      ggdagitty <- ggdagitty(dag, labels)
      temp_list[[epoch]] <- ggdagitty
      epoch <- epoch + 1
    }
    dag_list[[l]] <- dags
    dagitty_list[[l]] <- temp_list
    lambda <- lambda + l_increment
    l <- l + 1
  }

  #lambda <- lambda_min
  l <- 1
  dag_to_save <- list()
  dagitty_to_save <- list()
  for(l in 1:lambda_range){
    temp_list <- dagitty_list[[l]]
    temp_dag <- dag_list[[l]]
    n <- length(temp_list)
    nCol <- floor(sqrt(n))
    do.call("grid.arrange", c(temp_list, ncol=c(nCol)))
    epoch <- as.numeric(readline(cat("Select a graph (1 to ", iterations, "): ")))
    if(is.na(epoch)){
      cat("\nSkipping lambda value", l, "\n")
      l <- l + 1
    }else if(epoch <=iterations && epoch >= 1){
      cat("\nSaving...\n")
      dagitty_to_save[l] <- temp_list[as.numeric(epoch)]
      dag_to_save[l] <- temp_dag[as.numeric(epoch)]
      l <- l + 1
    }else{
      cat("\nSkipping lambda value", l, "\n")
      l <- l + 1
    }
  }

  dagitty_to_save <- dagitty_to_save[!sapply(dagitty_to_save,is.null)]
  dag_to_save <- dag_to_save[!sapply(dag_to_save,is.null)]
  ?string()
  dag_out <- list()
  dagitty_out <- list()
  if(length(dag_to_save) == 1){
    dagitty_to_save
    dag <- dagitty::dagitty(as.string(dag_to_save))
    return(dag)
  }else if(length(dag_to_save) >5){
    n <- length(dag_to_save)
    n_1 <- as.integer(n/2)
    n_2 <- n
    check_two_iterations <- TRUE
  }else{
    check_two_iterations <- FALSE
    print(check_two_iterations)
    print("only one iteration")
  }
  m <- 1
  if(check_two_iterations == TRUE){
    dag_list_2 <- list()
    dagitty_list_2 <- list()
    n <- n_1
    l <- 1
    for(m in 1:2){
      check_ans <- FALSE
      while(check_ans == FALSE){
        temp_list <- dagitty_to_save[l:n]
        n <- length(temp_list)
        nCol <- floor(sqrt(n))
        do.call("grid.arrange", c(temp_list, ncol=c(nCol)))
        choice <- as.numeric(readline(cat("Select a graph (1 to ", n, "): ")))
        if(is.na(choice)){
          cat("\nSkipping selection. \n")

          m <- m + 1

          check_ans <- TRUE

        }else if(choice <=n && choice >= 1){
          cat("\nSaving...\n")

          dag_list_2[m] <- dag_to_save[as.numeric(choice)]
          dagitty_list_2[m] <- dagitty_to_save[as.numeric(choice)]

          m <- m + 1

          check_ans <- TRUE
        }

      }
      l <- n
      n <- n_2

    }

    if(length(dag_list_2) == 1){
      dagitty_list_2
      dag <- dagitty::dagitty(dag_list_2)
      return(dag)
    }else{
      #suppressWarnings(plot(NULL))
      cat("/nChoose an output from the following graphs.")
      temp_list <- dagitty_list_2
      temp_dag <- dag_list_2
      n <- length(temp_list)
      nCol <- floor(sqrt(n))
      do.call("grid.arrange", c(temp_list, ncol=c(nCol)))
      choice <- as.numeric(readline(cat("Select a graph (1 to ", n, "): ")))
      if(is.na(choice)){
        cat("\nSkipping selection. \n")
      }else if(choice <=n && choice >= 1){
        cat("\nOutputting selected graph.\n")
        dag <- dag_to_save[[as.numeric(choice)]]
        dagitty_to_save <- dagitty_to_save[[as.numeric(choice)]]
      }else{
        cat("\nSkipping selection. \n")
      }
      dagitty_to_save
      dag <- dagitty::dagitty(dag)
      return(dag)
    }
    dagitty_to_save
    dag <- dagitty::dagitty(dag)
    return(dag)

  }else if(check_two_iterations == FALSE){
    #suppressWarnings(plot(NULL))
    cat("\nChoose an output from the following graphs.")
    temp_list <- dagitty_to_save
    n <- length(temp_list)
    nCol <- floor(sqrt(n))
    do.call("grid.arrange", c(temp_list, ncol=c(nCol)))
    choice <- as.numeric(readline(cat("Select a graph (1 to ", iterations, "): ")))
    if(is.na(choice)){
      cat("\nSkipping selection. \n")
    }else if(choice <=n && choice >= 1){
      cat("\nOutputting selected graph.\n")
      dag <- dag_to_save[as.numeric(choice)]
      dagitty_to_save <- dagitty_to_save[as.numeric(choice)]
      dagitty_to_save
    }else{
      cat("\nSkipping selection. \n")
    }
    dagitty_to_save
    return(dag)
  }
}


#' Add coordinates to a dagitty object
#'
#' @importFrom dagitty exposures outcomes coordinates latents topologicalOrdering
#' @importFrom ggdag time_ordered_coords
#' @param dag dagitty object
#' @param confounders vector of confounders nodes in the supplied dag.
#' @param mediators vector of mediator nodes in the supplied dag.
#' @param instrumental_vars vector of instrumental variable nodes in the supplied dag.
#' @param lambda adjusts sensitivity of node placement. Higher lambda introduces more volatility in confounder nodes and treatment along the y-axis.
#' @return dagitty objecty with coordinates.
#' @noRd
add_coordinates <- function(dag, treatments, outcomes, confounders, mediators, instrumental_vars, mediator_outcome_confounders, competing_exposures, latent_vars, colliders, lambda, observed = NULL){

  o <- length(outcomes)
  t <- length(treatments)

  num_confounders <- c <- length(confounders)

  if(!all(is.null(instrumental_vars))){

    instrumental_variables <- unique( draw_iv_edges(instrumental_vars, treatments)$ancestor )

    i <- length(instrumental_variables)

  }else{

    instrumental_variables <- NULL

    i <- 0
  }
  if(!all(is.null(latent_vars))){

    latent_variables <- get_latent_vec(latent_vars)

    l <- length(latent_variables)

  }else{

    latent_variables <- NULL

    l <- 0
  }

  if(!all(is.null(mediators))){

    m <- length(mediators)

  }else{
    m <- 0
  }
  if(!all(is.null(mediator_outcome_confounders))){

    moc <- length(mediator_outcome_confounders)

  }else{
    moc <- 0
  }
  if(!all(is.null(competing_exposures))){

    comp <- length(competing_exposures)

  }else{
    comp <- 0
  }
  if(!all(is.null(colliders))){

    coll <- length(colliders)

  }else{
    coll <- 0
  }


  var_order <- na.omit( c(unlist(confounders, recursive = FALSE), mediator_outcome_confounders, mediators, instrumental_variables, competing_exposures, colliders, treatments, outcomes) )


  num_vars <- n <- length(unlist(var_order))

  lambda <- lambda/(c + (n-i))

  #print(paste("Coordinates added with randomization multiplier set to:", lambda/(num_confounders + (n-i))))

  # get default ordered coordinates
  coords <- ggdag::time_ordered_coords(unlist(var_order), direction = c("x"), auto_sort_direction = c("right"))
  coords$y <- (ggdag::time_ordered_coords(unlist(var_order), direction = c("y"), auto_sort_direction = c("right")))$y

  # randomize default coordinates
  rand_num_x <- runif(num_vars, min = -lambda/2, max = lambda/2)
  rand_num_y <- runif(num_vars, min = -lambda/2, max = lambda/2)

  coords$y <- coords$y + rand_num_y
  coords$x <- coords$x + rand_num_x

  # confounders
  rand_conf_x_y <- sort(rlnorm(num_confounders*2), decreasing = FALSE)*(lambda/(n-i))

  rand_conf_x <- rev(rand_conf_x_y)
  rand_conf_x <- rand_conf_x[1:(num_confounders+1)]
  rand_conf_y <- rand_conf_x_y[1:(num_confounders+1)]

  rand_conf_x <- c(0, seq(c*(lambda/3), (c)*(lambda), length.out = num_confounders))[1:(num_confounders)] + rand_conf_x[1:(num_confounders)]
  rand_conf_y <- c(0, seq(-(c)*(lambda/3), -c*(lambda^2), length.out = num_confounders/2),
                   seq(c*(lambda/3), c*(lambda), length.out = num_confounders/2))[1:(num_confounders)] + (rand_conf_y*(n-i))[1:(num_confounders)]

  coords[1:num_confounders,]$x <- coords[1:num_confounders,]$x + rand_conf_x[1:num_confounders]
  coords[1:num_confounders,]$y <- coords[1:num_confounders,]$y + rand_conf_y[1:num_confounders]

  # mediator-outcome-confounders
  if( all( !is.null(mediator_outcome_confounders) ) ){

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

  trt_y_max <- seq( 1.1 - length(t)/3,
                    1.1 + length(t)/n,  length.out = t )

  trt_x_min <- seq( ( ((num_confounders/num_vars)*(lambda*num_vars)) - (i/n) ) ,
                    ( ((num_confounders/num_vars)*(lambda*num_vars)) + (i/3) ) ,  length.out = t )

  coords[ ( nrow(coords) -t):(nrow(coords)-1), ]$x <- ( coords[ ( nrow(coords) -t):(nrow(coords)-1), ]$x + trt_x_min )
  coords[ ( nrow(coords) -t):(nrow(coords)-1), ]$y <- trt_y_max

  outcome_x_min <- ((m-i)*(lambda)*2) + m - i
  outcome_y_min <- ((m-i)*(lambda)*2) + i

  coords[nrow(coords),]$x <- coords[nrow(coords),]$x + outcome_x_min
  coords[nrow(coords),]$y <- coords[nrow(coords),]$y - outcome_y_min

  # mediators
  if(all(!is.null(mediators))){

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

  # instrumental variables
  if(all(!is.null(instrumental_variables))){

    instr_y <- seq( ( trt_y_max - (1 + (i/n)) )[1] , (trt_y_max - i )[t], length.out = i) # fix these coords

    instr_x <- seq( ( ( coords[nrow(coords)-t-o+1,]$x ) - i*(c/n) ) , ( ( coords[nrow(coords)-o-1,]$x ) + i) ,  length.out = i)

    coords[(num_confounders+moc+m+1):(num_confounders+moc+m+i),]$y <- instr_y
    coords[(num_confounders+moc+m+1):(num_confounders+moc+m+i),]$x <- instr_x

  }

  # competing exposures
  if(all(!is.null(competing_exposures))){

    comp_y <- seq( ( coords[nrow(coords),]$y - (comp*(3/comp)) ) , (coords[nrow(coords),]$y - (comp*(2/comp))*(c/n) ), length.out = comp)

    comp_x <- seq( ( coords[nrow(coords),]$x ) - comp , ( coords[nrow(coords),]$x ) - (comp*(c/n) )  , length.out = comp)

    coords[(num_confounders+moc+m+i+1):(num_confounders+moc+m+i+comp),]$y <- comp_y
    coords[(num_confounders+moc+m+i+1):(num_confounders+moc+m+i+comp),]$x <- comp_x

  }

  # colliders
  if(all(!is.null(colliders))){

    coll_y <- seq( ( coords[nrow(coords),]$y + (coll*(coll/3)) ) , (coords[nrow(coords),]$y + (coll*(coll/2))*(c/n) ), length.out = coll)

    coll_x <- seq( ( coords[nrow(coords),]$x ) + coll + 0.5 , ( coords[nrow(coords),]$x ) + (coll*(coll/2) )  , length.out = coll)

    coords[(num_confounders+moc+m+i+comp+1):(num_confounders+moc+m+i+comp+coll),]$y <- coll_y
    coords[(num_confounders+moc+m+i+comp+1):(num_confounders+moc+m+i+comp+coll),]$x <- coll_x

  }
  coords_x <- coords$x
  coords_y <- coords$y
  names(coords_x) <- coords$name
  names(coords_y) <- coords$name

  coords_list <- list(
    x = coords_x,
    y = coords_y
  )


  existing_coords <- dagitty::coordinates(dag)

  coords_list <-  check_existing_coordinates(dag, observed, existing_coords, coords_list)

  coords_list <- generate_latent_coordinates(dag, latent_variables, coords_list)

  dagitty::coordinates(dag) <- coords_list

  return(dag)
}


#' Checks if dag previously had coordinates
#'
#' @param dag dagitty object
#' @param grouped_nodes vector of grouped variable nodes.
#' @param existing_coords existing dag coordinates
#' @param coords_list add_coordinates() generated coordinates
#' @return coords_list object with new group coordinates.
#' @noRd
check_existing_coordinates <-  function(dag, grouped_nodes, existing_coords, coords_list){

  if( all( is.na( unlist(existing_coords) ) ) & all( !is.na(grouped_nodes) ) ){

    coords_list <- generate_grouped_coordinates(dag, grouped_nodes, coords_list)

  }


  return(coords_list)
}


#' Creates latent node coordinates for add_latent_coordinates()
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @param dag dagitty object
#' @param latent_variables vector of latent variable nodes.
#' @param existing_coordinates vector of latent variable nodes.
#' @return dagitty objecty with coordinates.
#' @noRd
generate_grouped_coordinates <- function(dag, grouped_nodes, existing_coords){
  #grouped_nodes <- observed

  if( !all(is.na(grouped_nodes))){
  num_nodes <- length(grouped_nodes)

  num_vars <- length(names(dag))

  nodes_descendants <- lapply(1:num_nodes, function(x){

    nodes_descendants <- dagitty::children(dag, grouped_nodes[x])

  })

  if( length(unlist(nodes_descendants)) < (num_nodes*2) ){

    nodes_parents <- lapply(1:num_nodes, function(x){

      nodes_parents <- dagitty::parents(dag, grouped_nodes[x])

    })

    nodes_descendants <- lapply(1:num_nodes, function(x){

      nodes_descendants <- c(unlist(nodes_parents[x]), unlist(nodes_descendants[x]))

    })

  }


  quality_check <- FALSE

  iteration <- 0

  while(quality_check == FALSE){

    iteration <- iteration + 1

    new_group_coords <- lapply(1:num_nodes, function(x){

      nodes_descendants_coords_x <- (
        mean( tidyr::replace_na(is.na(existing_coords$x[ names(existing_coords$x) %in% nodes_descendants[[x]] ]), 0))

        - runif(n = 1, min = ( 0.1*iteration )*( num_vars/num_nodes ),
                max = ( ( 0.2*iteration )*( num_vars/num_nodes ) ) ) ** (2*(0.1*iteration))

        - (2*( 1 / x )) - ( (0.2*iteration)*( num_vars ) )
      )

      nodes_descendants_coords_y <- (
        mean( tidyr::replace_na(is.na(existing_coords$y[ names(existing_coords$y) %in% nodes_descendants[[x]] ] ), 0))

        - runif(n = 1, min = ( 0.01*iteration )*( num_vars/num_nodes ),
                max = ( ( 0.1*iteration )*( num_vars/num_nodes ) ) ) ** (2*(0.1*iteration))

        - (2*( 1 / x )) - ( (0.2*iteration)*( num_vars ) )
      )

      new_group_coords <- list(x = unlist(nodes_descendants_coords_x), y = unlist(nodes_descendants_coords_y))

    })

    # bind rows into a data frame (columns are still classed as lists)
    new_coordinates <- as.data.frame( do.call( rbind, new_group_coords ) )


    ## unlist and replace NaN ##

    # x coordinates
    unlisted_coords_x <- unlist( new_coordinates$x )

    unlisted_coords_x <- sapply(  unlisted_coords_x, function(x){
      replace( x, x == "NaN" , num_vars + runif(1, min=1, max = (num_vars/2) ) )
    } )

    # y coordinates
    unlisted_coords_y <- unlist( new_coordinates$y )

    unlisted_coords_y <- sapply(  unlisted_coords_y, function(x){
      replace( x, x == "NaN" , num_vars + runif(1, min=1, max = (num_vars/2) ) )
    } )



    ## internal check against other generated latent coords ##

    if( any( unlisted_coords_x) < 0 ){

      unlisted_coords_x <- lapply(1:length(unlisted_coords_x), function(x){

        unlisted_coords_x[x] <- sqrt( unlist( unlisted_coords_x[x] )**2 )/2

      })

      unlisted_coords_x <- unlist(unlisted_coords_x)
    }

    if( any( unlisted_coords_y < 0 ) ){

      unlisted_coords_y <- lapply(1:length(unlisted_coords_y), function(x){

        unlisted_coords_y[x] <- sqrt( unlist( unlisted_coords_y[x] )**2 )/2

      })

      unlisted_coords_y <- unlist(unlisted_coords_y)

    }

    new_coordinates <- list(x = unlisted_coords_x, y = unlisted_coords_y)

    quality_check <- quality_check_coords(dag, grouped_nodes, num_nodes, new_coordinates, existing_coords)

  }
  names(new_coordinates$x) <- grouped_nodes
  names(new_coordinates$y) <- grouped_nodes

  coords_list <- list(x = c( na.omit(existing_coords$x), new_coordinates$x), y = c( na.omit(existing_coords$y), new_coordinates$y))

  }else{

    coords_list <- existing_coords
  }


  return(coords_list)

}


#' Add coordinates to latent nodes
#'
#' @importFrom dagitty coordinates
#' @param dag dagitty object
#' @param latent_variables vector of latent variable nodes.
#' @return dagitty objecty with coordinates.
#' @noRd
add_latent_coordinates <- function(dag, latent_variables){

  coordinates <- dagitty::coordinates(dag)

  coordinates <- generate_latent_coordinates(dag, latent_variables, coordinates)

  dagitty::coordinates(dag) <- coordinates

  return(dag)
}

#' Creates latent node coordinates for add_latent_coordinates()
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @param dag dagitty object
#' @param latent_variables vector of latent variable nodes.
#' @param existing_coordinates vector of latent variable nodes.
#' @return dagitty objecty with coordinates.
#' @noRd
generate_latent_coordinates <- function(dag, latent_variables, existing_coords){
  #existing_coords <- coords_list

  if( !all(is.na(latent_variables))){

    num_latents <- length(latent_variables)

    num_vars <- length(names(dag))

    latent_descendants <- lapply(1:num_latents, function(x){

      latent_descendants <- dagitty::children(dag, latent_variables[x])

    })


    quality_check <- FALSE

    iteration <- 0

    while(quality_check == FALSE){

      iteration <- iteration + 1

      latent_coords <- lapply(1:num_latents, function(x){

        latent_descendants_coords_x <- ( mean( existing_coords$x[ names(existing_coords$x) %in% latent_descendants[[x]] ])
                                         - runif(n=1, min=0.1*(num_vars/num_latents ), max= (0.2*( num_vars/num_latents ) ) ) ** (2*( (0.1*iteration) / x ))
                                         - (2*( (0.1*iteration) / x ))
                                         - ( 0.1*( num_vars ) )
        )

        latent_descendants_coords_y <- ( mean( existing_coords$y[ names(existing_coords$y) %in% latent_descendants[[x]] ] )
                                         - runif(n=1, min=0.01*( num_vars/num_latents ), max= (0.1*( num_vars/num_latents)) ) ** (2*( (0.1*iteration) / x ))
                                         - (2*( (0.1*iteration) / x ))  - ( 0.1*( num_vars ) )
        )

        latent_coords <- list(x = unlist(latent_descendants_coords_x), y = unlist(latent_descendants_coords_y))

      })

      # bind rows into a data frame (columns are still classed as lists)
      new_coordinates <- as.data.frame( do.call( rbind, latent_coords ) )

      ## unlist and replace NaN ##

      # x coordinates
      unlisted_coords_x <- unlist( new_coordinates$x )

      unlisted_coords_x <- sapply(  unlisted_coords_x, function(x){
        replace( x, x == "NaN" , num_vars + runif(1, min=1, max = (num_vars/2) ) )
      } )

      # y coordinates
      unlisted_coords_y <- unlist( new_coordinates$y )

      unlisted_coords_y <- sapply(  unlisted_coords_y, function(x){
        replace( x, x == "NaN" , num_vars + runif(1, min=1, max = (num_vars/2) ) )
      } )


      ## internal check against other generated latent coords ##
      if( any(unlisted_coords_x < 0 ) ){

        unlisted_coords_x <- lapply(1:length(unlisted_coords_x), function(x){

          unlisted_coords_x[x] <- sqrt( unlist( unlisted_coords_x[x] )**2 )/2

        })

        unlisted_coords_x <- unlist(unlisted_coords_x)
      }

      if( any(unlisted_coords_y < 0 ) ){

        unlisted_coords_y <- lapply(1:length(unlisted_coords_y), function(x){

          unlisted_coords_y[x] <- sqrt( unlist( unlisted_coords_y[x] )**2 )/2

        })

        unlisted_coords_y <- unlist(unlisted_coords_y)

      }

      latent_coords <- list(x = unlisted_coords_x, y = unlisted_coords_y)


      quality_check <- quality_check_coords(dag, latent_variables, num_latents, latent_coords, existing_coords)

    }

    names(latent_coords$x) <- latent_variables
    names(latent_coords$y) <- latent_variables

    coordinates <- list(x = c(existing_coords$x, latent_coords$x), y = c(existing_coords$y, latent_coords$y))

  }else{

    coordinates <- existing_coords
  }

  return(coordinates)

}


#' Checks node coordinates are not overlapping
#'
#' @importFrom dagitty exposures outcomes coordinates latents
#' @param dag dagitty object
#' @param latent_variables vector of latent variable nodes.
#' @param existing_coordinates vector of latent variable nodes.
#' @return dagitty objecty with coordinates.
#' @noRd
quality_check_coords <- function(dag, grouped_nodes, num_nodes, new_coordinates, existing_coords){

  num_vars <- length(names(dag))

  ## internal check against other generated latent coords ##
  unlisted_coords_x <- new_coordinates$x
  unlisted_coords_y <- new_coordinates$y

  check_coords <- lapply(1:num_nodes, function(a){

    lapply(1:num_nodes, function(b){

      diff_list <-
        sum(unlist( sqrt( diff( range( c(unlisted_coords_x[a], unlisted_coords_x[b]) ) ) )**2 < 10/num_vars ), na.rm = TRUE)
      + sum(unlist( sqrt( diff( range( c(unlisted_coords_y[a], unlisted_coords_y[b]) ) ) )**2 < 10/num_vars ), na.rm = TRUE)

    })

  })

  if( sum(unlist(check_coords)) <= num_nodes ){

    return(TRUE)

  }


  ## external check against existing coords ##
  existing_coords_x <- existing_coords$x
  existing_coords_y <- existing_coords$y

  # remove group node names from existing coords vectors
  existing_coords_x[ !names(existing_coords_x) %in% grouped_nodes]
  existing_coords_y[ !names(existing_coords_y) %in% grouped_nodes]

  check_coords <- lapply(1:num_nodes, function(a){

    lapply(1:length(existing_coords_x), function(b){

      diff_list <-
        sum(unlist( sqrt( diff( range( c(unlisted_coords_x[a], existing_coords_x[b]) ) ) )**2 < 6/num_vars ), na.rm = TRUE)
      + sum(unlist( sqrt( diff( range( c(unlisted_coords_y[a], existing_coords_y[b]) ) ) )**2 < 6/num_vars ), na.rm = TRUE)

    })

  })

  if( sum(unlist(check_coords)) == 0 ){

    return(TRUE)

  }

  return(FALSE)

}

#' dataframe output from a dagitty object
#'
#'Generates a table similar to calling ggdag::tidy_dagitty on a ggdag::dagify object
#'The benefit of this function is that it automatically identifies exposure, outcome, confounder, observed and latent variables inputted from dagitty.net, whereas ggdag::tidy_dagitty only does this for ggdag::dagify objects.
#'Output can be used with ggdag to create better looking DAGs from dagitty.net code.
#'
#' @importFrom ggdag tidy_dagitty
#' @param dag dagitty object
#' @param labels vector of labels for nodes in dagitty object
#' @return dag_df DAG as a dataframe for use with ggdag to create better looking DAGs
#' @noRd
tidy_ggdagitty <- function(dag, labels = NULL){
  # Cleaning the dags and turning it into a data frame.
  dag_df <- data.frame(ggdag::tidy_dagitty(dag))

  # flip y axis for ggplot
  dag_df$y <- dag_df$y*-1
  dag_df$yend <- dag_df$yend*-1

  if(is.null(labels)){
    return(dag_df)
  }

  dag_df <- add_labels(dag, dag_df, labels)

  return(dag_df)
}


#' add labels to a dag dataframe
#'
#'Generates a table similar to calling ggdag::tidy_dagitty on a ggdag::dagify() object
#'The benefit of this function is that it automatically identifies exposure, outcome, confounder, observed and latent variables inputted from dagitty.net, whereas ggdag::tidy_dagitty only does this for ggdag::dagify objects.
#'Output can be used with ggdag to create better looking DAGs from dagitty.net code.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom dagitty exposures outcomes latents adjustmentSets
#' @param dag dagitty object
#' @param dag_df a dag object converted to data frame using the ggdag::tidy_dagitty() function, or similar.
#' @param labels vector of labels for nodes in dagify object
#' @return dagify DAG as a dataframe for use with ggdag to create better looking DAGs
#' @noRd
add_labels <- function(dag, dag_df, labels){

  # Labeling variables
  dag_df <- dag_df %>%
    dplyr::mutate(label = dplyr::case_when(
      name %in% attr(labels, "names") ~ as.data.frame(labels)[name,],
    ))

  dag_na <- dag_df %>% dplyr::filter(is.na(direction))

  dag_df <- dag_df %>% dplyr::filter(!is.na(direction))

  node_labels <- getEdges(dag)


  #node_labels <- node_labels[,c(1:3)]
  #dag_df_subset <- dag_df[,c(1,4,5)]
  #rownames(dag_df_subset) <- NULL
  #rownames(dag_df_subset) <- NULL
  #node_labels[ !node_labels %in% dag_df_subset]

  dag_df <- dag_df[order(dag_df$name),]


  dag_df <- dag_df %>%
    dplyr::mutate(role = dplyr::case_when(
      name == node_labels$ancestor ~ node_labels$role_ancestor,
    ))

  dag_na <- dag_na %>%
    dplyr::mutate(role = dplyr::case_when(
      name %in% dagitty::outcomes(dag) ~ "outcome",
      name %in% colliders(dag) ~ "collider",
      .default = "observed",
    ))



  dag_df <- rbind(dag_df, dag_na)


  return(dag_df)
}


