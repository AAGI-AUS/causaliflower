
################################################################################
# Functions:
# ggdagitty()
# bigDagitty()
# bigDagify()
################################################################################

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
  col <- c("outcome"="cadetblue2",
           "treatment"="darkolivegreen3",
           "confounder"="coral2",
           "mediator"="cornflowerblue",
           "mediator-outcome confounder"="darkorange",
           "instrumental"="magenta",
          # "collider"="darkred",
           "latent"="grey")

  shape <- c("outcome"=19,
             "treatment"=19,
             "confounder"=19,
             "mediator"=19,
             "mediator-outcome confounder"=19,
             "instrumental"=19,
            # "collider"=19,
             "latent"=1)

  # variable for legend order
  order_col <- c("outcome", "treatment", "confounder", "mediator", "mediator-outcome confounder", "instrumental", #"collider",
                 "latent")

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


  dag_df <- dagittyData(dag, labels)

  if(label_placement == "label_repel"){


    ggdag <- ggplot2::ggplot(data = dag_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color=role, shape = role, fill = role)) +
      ggdag::geom_dag_edges(start_cap = ggraph::circle(5, 'mm'),
                            end_cap = ggraph::circle(5, 'mm')) +
      ggdag::geom_dag_point(size=10) +
      ggdag::geom_dag_label_repel(ggplot2::aes(label = label), colour = "white", alpha = 0.8, show.legend = FALSE) +
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

  dag_df <- dagittyData(dag)

  #dag_na <- dag_df %>% dplyr::filter(is.na(direction))

  edges <- getEdges(dag)

  node_labels <- unique(edges[c("ancestor", "role_ancestor")])

  #node_labels <- rbind(node_labels, dag_na)


  if(is.null(confounders)){

    confounders <-  as.vector(unlist(lapply(as.vector( ( node_labels %>% dplyr::filter(role_ancestor == "confounder") )[,1]), function(x) if(identical(x, character(0))) NA_character_ else x)))
  }

  mediators <-  as.vector(unlist(lapply(as.vector( ( node_labels %>% dplyr::filter(role_ancestor == "mediator") )[,1]), function(x) if(identical(x, character(0))) NA_character_ else x)))
  instrumental_vars <- as.vector(unlist(lapply(as.vector( ( node_labels %>% dplyr::filter(role_ancestor == "instrumental") )[,1]), function(x) if(identical(x, character(0))) NA_character_ else x)))
  latent_vars <- as.vector(unlist(lapply(as.vector( ( node_labels %>% dplyr::filter(role_ancestor == "latent")  )[,1]), function(x) if(identical(x, character(0))) NA_character_ else x)))
  mediator_outcome_confounders <- as.vector(unlist(lapply(as.vector( ( node_labels %>% dplyr::filter(role_ancestor == "mediator-outcome-confounder")  )[,1]), function(x) if(identical(x, character(0))) NA_character_ else x)))


  if(length(na.omit(coords_spec)) == 3){
    message("Generating coordinates for ", coords_spec["iterations"], " DAGs for each lambda value between ", coords_spec["lambda"], " and ", coords_spec["lambda_max"], ".")
    dag <- addAutoCoords(dag, confounders, mediators, instrumental_vars, mediator_outcome_confounders, coords_spec)

    return(dag)

  }else if(length(na.omit(coords_spec)) == 1){

    tryCatch(
      { # try

        dag <- addCoords(dag, confounders, mediators, instrumental_vars, mediator_outcome_confounders, na.omit(coords_spec))

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
addAutoCoords <- function(dag,
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
      dag <- addCoords(dag, confounders, mediators, lambda)
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
addCoords <- function(dag, confounders, mediators, instrumental_vars, mediator_outcome_confounders, lambda){


  if(is.null(confounders)){

    var_order <- sort(unlist(dagitty::topologicalOrdering(dag)))

  }else{
    #confounders_order <- ggdag::time_ordered_coords(confounders, direction = c("x"), auto_sort_direction = c("right"))$name


    var_order <- sort(unlist(dagitty::topologicalOrdering(dag)))
    var_order <- as.vector(names(var_order))
    var_order <- c(confounders, var_order[( length(confounders) +1 ):length(var_order)])


  }


  var_names <- names(var_order)
  confounders_order <- var_names[var_names %in% confounders]
  confounders <- c(confounders[order(match(confounders, confounders_order))])

  treatment <- dagitty::exposures(dag)
  outcome <- dagitty::outcomes(dag)
  latent_vars <- dagitty::latents(dag)

  if(!all(is.null(confounders))){

    c <- length(confounders)

    unlisted_c <- length(unlist(confounders))

  }else{
    c <- 0
    unlisted_c <- 0
  }
  if(!all(is.null(mediators))){

    m <- length(mediators)

  }else{
    m <- 0
  }
  if(!all(is.null(instrumental_vars))){

    i <- length(instrumental_vars)

  }else{
    i <- 0
  }
  if(!all(is.null(mediator_outcome_confounders))){

    moc <- length(mediator_outcome_confounders)

  }else{
    moc <- 0
  }


  if(unlisted_c > c){

    confounders <- na.omit(confounders)
    vars <- confounders

    length <- c
    if(any(!is.null(mediator_outcome_confounders))){
      mediator_outcome_confounders <- na.omit(mediator_outcome_confounders)
      vars[[length + 1]] <- mediator_outcome_confounders
      length <- length + 1
    }
    if(any(!is.null(mediators))){
      mediators <- na.omit(mediators)
      vars[[length + 1]] <- mediators
      length <- length + 1
    }
    if(any(!is.null(instrumental_vars))){
      instrumental_vars <- na.omit(instrumental_vars)
      vars[[length + 1]] <- instrumental_vars
      length <- length + 1
    }
    vars[[length + 1]] <- treatment
    vars[[length + 2]] <- outcome
  }else{
    vars <- confounders
    if(any(!is.null(mediator_outcome_confounders))){
      mediator_outcome_confounders <- na.omit(mediator_outcome_confounders)
      vars <- c(vars, mediator_outcome_confounders)
    }
    if(any(!is.null(mediators))){
      mediators <- na.omit(mediators)
      vars <- c(vars, mediators)
    }
    if(any(!is.null(instrumental_vars))){
      instrumental_vars <- na.omit(instrumental_vars)
      vars <- c(vars, instrumental_vars)
    }
    vars <- c(vars, treatment, outcome)
  }


  num_vars <- n <- length(unlist(vars))
  num_confounders <- length(unlist(confounders))

  lambda <- lambda/(num_confounders + (n-i))

  #print(paste("Coordinates added with randomization multiplier set to:", lambda/(num_confounders + (n-i))))

  # get default ordered coordinates
  coords <- ggdag::time_ordered_coords(vars, direction = c("x"), auto_sort_direction = c("right"))
  coords$y <- (ggdag::time_ordered_coords(unlist(vars), direction = c("y"), auto_sort_direction = c("right")))$y

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
  if(all(!is.null(mediator_outcome_confounders))){

    a <- 1 + as.integer(!is.na((lambda/lambda)))
    b <- 1 + as.integer(is.na((lambda/lambda)))
    x <- lambda + (as.integer(is.na((lambda/lambda))) / 2)


    rand_moc_y <- seq( -( (((b*m^2)*(x^2))/(a*(n-i))) + ((b*m)*x)/a - (m/2) ),
                       -( (((m^2)*(x^2))/(n-i)) + (m)*x - (m/2) ),  length.out = moc)

    rand_moc_x <- seq( ( c/(moc + m)),
                       ( x*(m*x + m) + c/m),  length.out = moc)

    coords[(num_confounders+1):(num_confounders+moc),]$y <- coords[(num_confounders+1):(num_confounders+moc),]$y + rand_moc_y
    coords[(num_confounders+1):(num_confounders+moc),]$x <- coords[(num_confounders+1):(num_confounders+moc),]$x + rand_moc_x

  }

  # mediators
  if(all(!is.null(mediators))){

    a <- 1 + as.integer(!is.na((lambda/lambda)))
    b <- 1 + as.integer(is.na((lambda/lambda)))
    x <- lambda + (as.integer(is.na((lambda/lambda))) / 2)


    rand_med_y <- seq( -( (((b*m^2)*(x^2))/(a*(n-i))) + ((b*m)*x)/a + (m/2) ),
                       -( (((m^2)*(x^2))/(n-i)) + (m)*x + (m/2) ),  length.out = m)

    rand_med_x <- seq( ( (((2*m^2)*(x^2))/(n-i) ) + (m)*x + (m/(4/b)) ),
                       ( x*(m*x + m) + m ),  length.out = m)

    coords[(num_confounders+moc+1):(num_confounders+moc+m),]$y <- coords[(num_confounders+moc+1):(num_confounders+moc+m),]$y + rand_med_y
    coords[(num_confounders+moc+1):(num_confounders+moc+m),]$x <- coords[(num_confounders+moc+1):(num_confounders+moc+m),]$x + rand_med_x

  }

  trt_y_max <- 1.1 + length(instrumental_vars)/2
  trt_x_min <- (((num_confounders/num_vars)*(lambda*num_vars)) - m)

  coords[nrow(coords)-1,]$x <- coords[nrow(coords)-1,]$x + trt_x_min
  coords[nrow(coords)-1,]$y <- trt_y_max

  outcome_x_min <- ((m-i)*(lambda)*2) + m - i
  outcome_y_min <- ((m-i)*(lambda)*2) + i

  coords[nrow(coords),]$x <- coords[nrow(coords),]$x + outcome_x_min
  coords[nrow(coords),]$y <- coords[nrow(coords),]$y - outcome_y_min

  # instrumental variables
  if(!is.null(instrumental_vars)){

    instr_y <- seq( ( trt_y_max - i ) , (trt_y_max + i*(c/n) ), length.out = i)

    instr_x <- seq( ( ( coords[nrow(coords)-1,]$x ) - i) , ( ( coords[nrow(coords)-1,]$x ) + i*(c/n) ) , length.out = i)

    coords[(num_confounders+moc+m+1):(num_confounders+moc+m+i),]$y <- instr_y
    coords[(num_confounders+moc+m+1):(num_confounders+moc+m+i),]$x <- instr_x

  }

  coords_x <- coords$x
  coords_y <- coords$y
  names(coords_x) <- coords$name
  names(coords_y) <- coords$name

  coords_list <- list(
    x = coords_x,
    y = coords_y
  )

  dagitty::coordinates(dag) <- coords_list

  return(dag)
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
dagittyData <- function(dag, labels = NULL){
  # Cleaning the dags and turning it into a data frame.
  dag_df <- data.frame(ggdag::tidy_dagitty(dag))

  # flip y axis for ggplot
  dag_df$y <- dag_df$y*-1
  dag_df$yend <- dag_df$yend*-1

  if(is.null(labels)){
    return(dag_df)
  }

  dag_df <- addLabels(dag, dag_df, labels)

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
addLabels <- function(dag, dag_df, labels){

  # Labeling variables
  dag_df <- dag_df %>%
    dplyr::mutate(label = dplyr::case_when(
      name %in% attr(labels, "names") ~ as.data.frame(labels)[name,],
    ))

  dag_na <- dag_df %>% dplyr::filter(is.na(direction))

  dag_df <- dag_df %>% dplyr::filter(!is.na(direction))

  node_labels <- getEdges(dag)

  dag_df <- dag_df[order(dag_df$name),]


  dag_df <- dag_df %>%
    dplyr::mutate(role = dplyr::case_when(
      name == node_labels$ancestor ~ node_labels$role_ancestor,
    ))

  dag_na <- dag_na %>%
    dplyr::mutate(role = dplyr::case_when(
      name %in% dagitty::outcomes(dag) ~ "outcome",
      #name %in% name ~ "collider",
    ))



  dag_df <- rbind(dag_df, dag_na)


  return(dag_df)
}


