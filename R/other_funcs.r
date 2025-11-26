
#' Generate multiple versions of DAG coordinates
#'
#' @importFrom magrittr %>%
#' @importFrom gridExtra grid.arrange
#' @importFrom dagitty dagitty
#' @param dag dagitty object
#' @param confounders Vector of confounder variables, e.g. c("Z1", "Z2", "Z3").
#' @param mediators Vector of mediator variables, e.g. c("M1", "M2", "M3").
#' @param instrumental_variables vector of instrumental variable nodes in the supplied dag.
#' @param coords_spec Vector containing 'iterations' for specifying number of training epochs and 'lambda_range' to control the volatility of random coordinates.
#' @return dag objecty with coordinates.
#' @noRd
auto_coords <- function(dag,
                        confounders,
                        mediators = NULL,
                        instrumental_variables = NULL,
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
