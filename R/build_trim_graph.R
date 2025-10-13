#' Builds a dagitty object saturated (fully connected) graph
#'
#' buildGraph() produces a saturated graph by default from the supplied treatment, outcome, confounder, and  mediators (optional) based on the type of graph specified. Optional latent confounder or medfiator variables can be specified.
#'
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when filter
#' @importFrom dagitty dagitty
#' @param type The type of graph generated. Defaults to 'full', producing a fully connected graph with confounders connected in both directions (bi-directional), and to mediators in one direction (uni-directional). If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param variables Vector or list of variables to be assigned nodes in a graph. Order of variables determines the assigned coordinates. A named vector can be supplied, containing any of the other input variable roles. A list can also be supplied.
#' @param treatment Treatment variable name, e.g. "X". Must be specified, unless included in the named vector 'variables'.
#' @param outcome Outcome variable name, e.g. "Y". Must be specified, unless included in the named vector 'variables'.
#' @param confounders Character or vector of confounder variable names, e.g. "Z" or c("Z1", "Z2", "Z3"). Order of variables determines the assigned coordinates. A list can also be supplied. If type = "ordered", confounders in the same list will be assigned similar coordinates.
#' @param mediators Character or vector of mediator variable names, e.g. "M" or c("M1", "M2", "M3").
#' @param latent_variables Character or vector of additional or already supplied latent (unobserved) variable names, e.g. "U" or c("U1", "U2", "M1").
#' @param instrumental_variables Vector of instrumental variable names, e.g. "IV"
#' @param mediator_outcome_confounders Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param coords_spec Set of parameters for generating coordinates. Adjust node placement with lambda, a higher value increases volatility and results in more extreme DAG structures. Setting 'lambda_max' generates a DAG at each lambda value between lambda and lambda_max (only used if iterations is supplied). Iterations controls number of repeats for each lambda value (returns the first lambda value if NULL).
#' @returns A dagitty object, fully connected (saturated) graph.
#' @examples
#'
#'
#' # There are three main ways to supply variables to buildGraph().
#'
#' # Option 1: Inputting a named vector of variables inputted (can either include treatment and outcome or as separate inputs)
#' {
#'   variables <- c(confounders = c("Z1", "Z2", "Z3"),
#'                  treatment = "X",
#'                  outcome = "Y",
#'                  mediators = "M",
#'                  instrumental_variables = "IV")
#' }
#'
#' # Three graph types can also be generated. First, we build an 'ordered' graph where the supplied vector order is used to determine the temporal order of confounder nodes.
#'   type <- "ordered"
#'
#'   dag <- buildGraph(type = type,
#'                     variables = variables)
#'
#' # Plotting the graph using ggdagitty() assigns coordinates based on the variable roles (confounder, treatment, outcome, mediator, etc.).
#' ggdagitty(dag)
#'
#' #' # Option 2: The 'confounders' input can be used, while 'variables' is ignored.
#' #           Separate inputs are required for treatment, outcome, and other nodes.
#' {
#'   confounders <- c("Z1", "Z2", "Z3")
#'   treatment <- "X"
#'   outcome <- "Y"
#'   instrumental_variables <- "IV"
#'   mediator = "M"
#' }
#'
#' A "saturated" graph type bidirectionally connects each of the confounders, in addition to directed arrows to the outcome and treatment nodes.
#' type <- "saturated"
#'
#' dag <- buildGraph(type = type,
#'                   treatment = treatment,
#'                   outcome = outcome,
#'                   confounders = confounders,
#'                   mediators = mediators,
#'                   instrumental_variable = instrumental_variable)
#'
#' ggdagitty(dag, labels)
#'
#' # Option 3: The 'variables' input can be used to connect all variables to each other, treating them as confounders, besides treatment and outcome and any other separate inputs.
#' #           With this configuration, the actual 'confounders' input can be left blank. I decided to keep this as an option for users who may not have much 'exposure' to causal graphs.
#' {
#'   variables <- c("Z1", "Z2", "Z3")
#'   treatment <- "X"
#'   outcome <- "Y"
#'   instrumental_variables <- "IV"
#'   mediator = "M"
#' }
#'
#' # Using type = "full" generates a fully connected graph: all confounders are connected to each other, in both directions, and also to mediators, treatment and outcome.
#' type <- "full"
#'
#' dag <- buildGraph(type = type,
#'                   variables = variables,
#'                   treatment = treatment,
#'                   outcome = outcome,
#'                   mediators = mediators,
#'                   instrumental_variables = instrumental_variables)
#'
#' ggdagitty(dag)
#'
#' @export
buildGraph <- function(type = c("full", "saturated", "ordered"),
                       variables = NA,
                       treatment = NA,
                       outcome = NA,
                       confounders = NA,
                       mediators = NA,
                       latent_variables = NA,
                       instrumental_variables = NA,
                       mediator_outcome_confounders = NA,
                       competing_exposure = NA,
                       colliders = NA,
                       coords_spec = c(lambda = 0, lambda_max = NA, iterations = NA)){

  # Option 1: Named list of variables inputted, corresponding to role input options
  #           (treatment and outcome can be included in the named list, or supplied in their own inputs)

  if( all( is.na(confounders) ) & !all( is.na(variables) ) & !all( is.na(names(unlist(variables)))) ){

    variables_df <- extract_roles(variables)


    confounder_df <- variables_df %>% dplyr::filter(role == "confounder")
    confounder_vec <- confounder_df$variables

    if( length(unlist(confounder_vec)) > length(confounder_vec) ){
      confounders_list <- list()
      confounders_list <- lapply( 1:length(confounder_vec), function(x){

        confounders_list[[x]] <- list( as.data.frame( confounder_vec[[x]]), rep(x, times = length(confounder_vec[[x]]) ) )
        confounders_list <- dplyr::bind_cols(confounders_list[[x]][1], confounders_list[[x]][2])

      } )

      confounders_df <- dplyr::bind_rows( confounders_list )

      confounder_occurrance <- as.vector( unlist(confounders_df[2]) )

    }else{

      confounder_occurrance <- as.numeric(order(match(confounder_vec, confounder_vec)))

    }

    outcome_df <- variables_df %>% dplyr::filter(role == "outcome")
    outcome <- as.vector(unlist(outcome_df["variables"]))

    treatment_df <- variables_df %>% dplyr::filter(role == "treatment")
    treatment <- as.vector(unlist(treatment_df["variables"]))

    mediators <- variables_df %>% dplyr::filter(role == "mediator") %>% select(variables)
    latent_variables <- variables_df %>% dplyr::filter(role == "latent") %>% select(variables)
    instrumental_variables <- variables_df %>% dplyr::filter(role == "instrumental") %>% select(variables)
    mediator_outcome_confounders <- variables_df %>% dplyr::filter(role == "mediator_outcome_confounder") %>% select(variables)
    competing_exposure <- variables_df %>% dplyr::filter(role == "competing_exposure") %>% select(variables)
    colliders <- variables_df %>% dplyr::filter(role == "collider") %>% select(variables)

  }

  mediator_vec <- as.vector( unlist( lapply( mediators, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )
  latent_vec <-as.vector( unlist( lapply( latent_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )
  instrumental_vec <- as.vector( unlist( lapply( instrumental_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )
  m_o_confounder_vec <- as.vector( unlist( lapply( mediator_outcome_confounders, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )
  competing_exposure_vec <- as.vector( unlist( lapply( competing_exposure, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )
  collider_vec <- as.vector( unlist( lapply( colliders, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

  # Option 2: Vector of variables or confounders are inputted (treated as confounder nodes / assumed common causes of treatment and outcome),
  #           Separate inputs for treatment, outcome, and other nodes.(e.g. treatment, outcome, mediators) while confounders input is left blank/unused.

  if( any( is.na(confounders) ) & !any( is.na(variables) ) & any( is.na(names(unlist(variables))) ) ){
    # variables input is used while confounders input is not

    confounders <- variables

    variables <- NULL

  }

  if( !any( is.na(confounders) ) ){

    if( length(unlist(confounders)) > length(confounders) ){
      confounders_list <- list()
      confounders_list <- lapply( 1:length(confounders), function(x){

        confounders_list[[x]] <- list( as.data.frame( confounders[[x]]), rep(x, times = length(confounders[[x]]) ) )
        confounders_list <- dplyr::bind_cols(confounders_list[[x]][1], confounders_list[[x]][2])

      } )

      confounders_df <- dplyr::bind_rows( confounders_list )
      confounders_vec <- as.vector( unlist(confounders_df[1]) )
      confounder_occurrance <- as.vector( unlist(confounders_df[2]) )

    }else{

      variables_vec <- as.vector(unlist(confounders))
      confounder_vec <- as.vector(variables_vec[!variables_vec %in% mediator_vec & !variables_vec %in% instrumental_vec & !variables_vec %in% m_o_confounder_vec])
      confounder_occurrance <- as.numeric(order(match(confounder_vec, confounder_vec)))

    }

  }

  if( ( any( is.na(treatment) ) | any(is.na(outcome) ) & any( is.na(names(unlist(variables))) ) ) ){

    stop("Missing treatment and/or outcome inputs, or the supplied variables are not labelled. Both treatment and outcome must be given as inputs or supplied in a named vector for graphs using buildGraph().")

  }

  # if all mediator-outcome confounders are not confounders, execution is stopped
  if( all( confounder_vec %in% m_o_confounder_vec ) != FALSE ){

    stop("Confounders detected in 'mediator_outcome_confounder' input. These roles should be mutually exclusive. Please adjust supplied parameters and try again.")

  }


  ## outcome edges ##
  outcome_list <- c()
  # connect colliders
  outcome_list <- lapply(1:length(outcome), function(x){

    outcome_list[x] <- lapply(1:length(collider_vec), function(y){

      list( c( ancestor = outcome[x], edge = "->", descendant = collider_vec[y]) )

    })

  })

  outcome_list <- Filter( Negate(anyNA), unlist(unlist(outcome_list, recursive = FALSE), recursive = FALSE) )
  outcome_df <- dplyr::bind_rows(outcome_list)


  ## treatment edges ##
  treatment_list <- c()
  # connect outcome
  treatment_list <- lapply(1:length(treatment), function(x){

    treatment_list[x] <- list( c( ancestor = treatment[x], edge = "->", descendant = outcome) )

  })

  treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
  treatment_df <- dplyr::bind_rows(treatment_list)

  # connect mediators
  treatment_list <- suppressWarnings( lapply(1:length(treatment), function(x){

    treatment_list[x] <- lapply(1:length(mediator_vec), function(y){

      list( c( ancestor = treatment[x], edge = "->", descendant = mediator_vec[y]) )

    })

  }) )

  treatment_list <- Filter( Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE) )
  treatment_df <- dplyr::bind_rows(treatment_df, treatment_list)

  # connect colliders
  treatment_list <- suppressWarnings( lapply(1:length(treatment), function(x){

    treatment_list[x] <- lapply(1:length(collider_vec), function(y){

      list( c( ancestor = treatment[x], edge = "->", descendant = collider_vec[y]) )

    })

  }) )

  treatment_list <- Filter( Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE) )
  treatment_df <- dplyr::bind_rows(treatment_df, treatment_list)

  ## confounder edges ##
  confounder_list <- c()

  # connect outcome
  confounder_list <- lapply(1:length(confounder_vec), function(x){

    confounder_list[x] <- list( c( ancestor = confounder_vec[x], edge = "->", descendant = outcome) )

  })

  confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
  confounder_df <- dplyr::bind_rows(confounder_list)

  # connect treatment
  confounder_list <- lapply(1:length(confounder_vec), function(x){

    confounder_list[x] <- list( c( ancestor = confounder_vec[x], edge = "->", descendant = treatment) )

  })

  confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
  confounder_df <- dplyr::bind_rows(confounder_df, confounder_list)

  # connect all confounders
  if( type == "full" | type == "saturated" ){

    confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

      confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

        list( c( ancestor = confounder_vec[x], edge = "->", descendant = confounder_vec[y]) )

      })

    }) )

    confounder_list <- Filter(Negate(anyNA), unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE))
    confounder_df <- dplyr::bind_rows(confounder_df, confounder_list)

    # connect mediators if the inputted dag type is "full"
    if( type == "full" ){

      confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

        confounder_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = confounder_vec[x], edge = "->", descendant = mediator_vec[y]) )

        })

      }) )

      confounder_list <- Filter(Negate(anyNA), unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE))
      confounder_df <- dplyr::bind_rows(confounder_df, confounder_list)

    }

  }else if( type == "ordered" ){

    confounder_list <- suppressWarnings(lapply(1:length(confounder_vec), function(x){

      confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

        list( c( ancestor = confounder_vec[x], edge = "->", descendant = confounder_vec[y], ancestor_order = confounder_occurrance[x], descendant_order = confounder_occurrance[y] ) )

      })

    }))

    confounder_list <- Filter(Negate(anyNA), unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE)) # unnested list and remove NA's
    confounder_order_df <- dplyr::bind_rows(confounder_list) # bind list elements to a single data frame
    confounder_order_df <- confounder_order_df[,c(1:3)][!confounder_order_df$ancestor_order > confounder_order_df$descendant_order, ] # remove rows where temporal logic is not followed
    confounder_df <- rbind(confounder_df, confounder_order_df)
  }

  confounder_df <- unique(confounder_df) # remove duplicate edges
  confounder_df <- confounder_df[confounder_df$ancestor != confounder_df$descendant, ] # remove edges with identical ancestor and descendant node names


  ## mediator edges ##
  mediator_list <- c()

  # connect outcome
  mediator_list <- lapply(1:length(mediator_vec), function(x){

    mediator_list[x] <- list( c( ancestor = mediator_vec[x], edge = "->", descendant = outcome) )

  })

  mediator_list <- Filter( Negate(anyNA), unlist(mediator_list, recursive = FALSE) )
  mediator_df <- dplyr::bind_rows(mediator_list)

  # connect all mediators if the inputted dag type is "full"
  if(type == "full"){

    mediator_list <- suppressWarnings( lapply(1:length(m_o_confounder_vec), function(x){

      mediator_list[x] <- lapply(1:length(mediator_vec), function(y){

        list( c( ancestor = mediator_vec[x], edge = "->", descendant = mediator_vec[y]) )

      })

    }) )

    mediator_list <- Filter( Negate(anyNA), unlist(unlist(mediator_list, recursive = FALSE), recursive = FALSE) )
    mediator_df <- dplyr::bind_rows(mediator_df, mediator_list)

  }


  ## mediator_outcome_confounder edges ##
  moc_list <- c()

  # connect outcome
  moc_list <- lapply(1:length(m_o_confounder_vec), function(x){

    moc_list[x] <- list( c( ancestor = m_o_confounder_vec[x], edge = "->", descendant = outcome) )

  })

  moc_list <- Filter(Negate(anyNA), unlist(moc_list, recursive = FALSE))
  moc_df <- dplyr::bind_rows(moc_list)

  # connect treatment if the inputted dag type is "full"
  if(type == "full"){

    moc_list <- lapply(1:length(m_o_confounder_vec), function(x){

      moc_list[x] <- list( c( ancestor = treatment, edge = "->", descendant = m_o_confounder_vec[x] ) )

    })

    moc_list <- Filter(Negate(anyNA), unlist(moc_list, recursive = FALSE))
    moc_df <- dplyr::bind_rows(moc_df, moc_list)

  }

  # connect mediators
  moc_list <- suppressWarnings( lapply(1:length(m_o_confounder_vec), function(x){

    moc_list[x] <- lapply(1:length(mediator_vec), function(y){

      list( c( ancestor = m_o_confounder_vec[x], edge = "->", descendant = mediator_vec[y]) )

    })

  }) )

  moc_list <- Filter(Negate(anyNA), unlist(unlist(moc_list, recursive = FALSE), recursive = FALSE))
  moc_df <- dplyr::bind_rows(moc_df, moc_list)


  ## instrumental variable edges ##
  instrumental_list <- c()

  # connect outcome
  instrumental_list <- lapply(1:length(instrumental_vec), function(x){

    instrumental_list[x] <- list( c( ancestor = instrumental_vec[x], edge = "->", descendant = treatment) )

  })

  instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
  instrumental_df <- dplyr::bind_rows(instrumental_list)


  ## competing_exposure edges ##
  competing_exposure_list <- c()

  # connect outcome
  competing_exposure_list <- lapply(1:length(competing_exposure_vec), function(x){

    competing_exposure_list[x] <- list( c( ancestor = competing_exposure_vec[x], edge = "->", descendant = outcome) )

  })

  competing_exposure_list <- Filter(Negate(anyNA), unlist(competing_exposure_list, recursive = FALSE))
  competing_exposure_df <- dplyr::bind_rows(competing_exposure_list)


  ## row bind all edge data frames ##
  edges_df <- rbind(treatment_df, outcome_df, confounder_df, moc_df, mediator_df, instrumental_df, competing_exposure_df)

  edges_df <- unique(edges_df) # remove duplicate edges
  edges_df <- edges_df[edges_df$ancestor != edges_df$descendant, ] # remove edges with identical ancestor and descendant node names

  ## get variable names ##
  var_names <- as.vector(variables)

  exclude_names <-c(treatment, outcome, latent_vec)
  var_names <- var_names[!var_names %in% exclude_names]

  node_name_and_coords_vec <- c() # empty vector for variables names / coordinates

  ## add dagitty [treatment]/[outcome]/[latent] identifiers to variable names ##
  if(!is.na(latent_vec)) {
    node_name_and_coords_vec <- c(paste(treatment, " [exposure] ", sep=""),
                                  paste(outcome, " [outcome] ", sep=""),
                                  paste(latent_vec, " [latent] ", sep=""),
                                  paste(var_names, collapse=" "))
  }else{
    node_name_and_coords_vec <- c(paste(treatment, " [exposure] ", sep=""),
                                  paste(outcome, " [outcome] ", sep=""),
                                  paste(var_names, collapse=" "))
  }


    edges_unlist <- lapply(1:nrow(edges_df), function(x){

    edges_unlist <- paste(edges_df[x,], collapse=" ")

  })

  edges_vec <- paste(unlist(edges_unlist), collapse=" ")

  dag_code <- paste("dag {", paste(node_name_and_coords_vec, collapse=""), edges_vec, "}", sep = " ")

  dag <- dagitty::dagitty(dag_code)

  dag <- getCoords(dag, confounders = as.vector(confounder_vec))

  return(dag)

}

#' Trim edges in a graph
#'
#' trimGraph() trims connected edges based on causal criteria and/or user inputs.
#'
#' @importFrom magrittr %>%
#' @param dag A saturated graph dagitty object. Exposure and outcome must be indicated, and optionally can include assigned coordinates.
#' @param edges_to_keep A vector of directed arrows to be kept between (non-treatment and non-outcome) variables, e.g. c("Z1 -> Z2", "Z2 -> Z3"), c("y", "n", "y"), or c(TRUE, FALSE, TRUE).
#' @returns A dagitty object, with directed arrows removed based on edges_to_keep.
#' @examples
#' #######################
#'
#' # a fully connected graph
#' latent_variables <- c("Z3", "M1")
#' mediators <- c("M1", "M2")
#' confounders <- c("Z1", "Z2", "Z3")
#' treatment <- "X"
#' outcome <- "Y"
#' type <-"full"
#' dag <- buildGraph(type = type,
#'                                  treatment = treatment,
#'                                  outcome = outcome,
#'                                  confounders = confounders)
#'
#' # get graph edges
#' edges_vec <- trimGraph(dag, edges_to_keep = FALSE) # enter 'n' to save edges
#'
#' # Option 1: remove all edges between confounders
#' edges_to_keep <- c()
#' dag <- trimGraph(dag, edges_to_keep)
#' ggDagitty(dag, labels)
#'
#' # Option 2: supply vector of kept edges
#' edges_to_keep <- c('Z1 -> Z2', 'Z3 -> Z2')
#' dag <- trimGraph(dag, edges_to_keep)
#' labels <- getLabels(dag)
#' ggDagitty(dag, labels)
#'
#' # Option 3: run ESC-DAGs protocol and assess each edge with causal criteria
#' dag <- trimGraph(dag)
#' labels <- getLabels(dag)
#' ggDagitty(dag, labels)
#'
#' # Option 4: keep edges 'y' or 'n'
#' edges_to_keep <- c("n", "n", "y", "n", "n", "y")
#' dag <- trimGraph(dag, edges_to_keep)
#' labels <- getLabels(dag)
#' ggDagitty(dag, labels)
#'
#' # a saturated graph (confounders are not connected to mediators)
#' treatment <- "X"
#' outcome <- "Y"
#' confounders <- c("Z1", "Z2", "Z3")
#' mediators <- c("M1", "M2")
#' latent_variables <- c("Z3", "M1")
#' type <- "saturated"
#' coords_spec <- 3
#'
#' dag <- causaliflower::buildGraph(type = type,
#'                                  treatment = treatment,
#'                                  outcome = outcome,
#'                                  confounders = confounders)
#'
#' edges_to_keep <- c('Z1 -> Z2', 'Z3 -> Z2')
#' dag <- trimGraph(dag, edges_to_keep)
#'
#' dag <- getCoords(dag,
#'                     confounders,
#'                     mediators = mediators,
#'                     coords_spec = coords_spec)
#'
#' labels <- getLabels(dag)
#' ggDagitty(dag, labels)
#'
#'
#' # an ordered graph (confounders are connected in order of the supplied vector)
#' # adding coordinates and two mediator variables
#' confounders <- c("Z1", "Z2", "Z3")
#' treatment <- "X"
#' outcome <- "Y"
#' latents <- c("Z3", "M1") # latent confounder and mediator variables specified
#' mediators <- c("M1", "M2")
#' type <- "ordered"
#'
#' dag <- buildGraph(type = type,
#'                  treatment = treatment,
#'                  outcome = outcome,
#'                  confounders = confounders,
#'                  mediators = mediators,
#'                  latent_variables = latent_variables)
#'
#' dag <- trimGraph(dag)
#' labels <- getLabels(dag)
#' ggDagitty(dag, labels)
#'
#' @export
trimGraph <- function(dag, edges_to_keep = NA){

  treatment_or_outcome_edges <- getEdges(dag, c("treatment", "outcome"))[,1:3]

  treatment_or_outcome_edges_vec <- c()
  num_edges <- nrow(treatment_or_outcome_edges)
  arrow_count <- 1

  for(arrow_count in 1:num_edges){

    treatment_or_outcome_edges_vec[arrow_count] <- paste(treatment_or_outcome_edges[arrow_count,], collapse=" ")
    arrow_count <-arrow_count + 1

  }



  edges <- getEdges(dag, c("!treatment", "!outcome"))[,1:3]

  edges_vec <- c()
  num_edges <- nrow(edges)
  arrow_count <- 1

  for(arrow_count in 1:num_edges){

    edges_vec[arrow_count] <- paste(edges[arrow_count,], collapse=" ")
    arrow_count <-arrow_count + 1

  }

  latent_vec <- dagitty::latents(dag)

  check_ans <- FALSE
  if(any(is.null(edges_to_keep))){

    cat("\nEmpty 'edges_to_keep' vector supplied.", num_edges, "directed arrows will be removed:", "\n", sep=" ")
    print(edges, quote=FALSE)

  }else if(any(is.na(edges_to_keep))){

    check_skip_sequence <- FALSE

    cat("\nThere are", num_edges, "directed arrows to be assessed:", "\n", "\n", sep=" ")
    print(edges, quote=FALSE)
    cat("\nAssess the posited causal relationships? (ESC-DAGs causal criteria and counterfactual thought experiment sequence)", "\n")

    check_ans <- TRUE

  }else if(any(edges_to_keep == FALSE)){

    message("\nSkipped sequence.")

    edges_vec <- noquote(paste("c('", paste(edges_vec, collapse="', '"), "')", sep = ""))

    message("\nOutputting vector of edges.")

    return(edges_vec)

  }else{

    check_skip_sequence <- TRUE

    ans <- 1
    edge_name <- c()

    if(length(edges_to_keep) == num_edges){
      while(check_ans == FALSE){
        for(ans in 1:length(edges_to_keep)){
          if(edges_to_keep[ans] == "y" | edges_to_keep[ans] == "Y" | edges_to_keep[ans] == "yes" | edges_to_keep[ans] == TRUE){

            edge_name[ans] <- edges_vec[ans]
            ans <- ans + 1

          }else if(edges_to_keep[ans] == "n" | edges_to_keep[ans] == "N" | edges_to_keep[ans] == "no" | edges_to_keep[ans] == FALSE){

            ans <- ans + 1

          }else{

            stop("\nEntries in 'edges_to_keep' are not the correct syntax e.g. 'x -> y', or do not match 'y', 'Y', 'yes', TRUE or 'n', 'N', 'no', FALSE for kept/discarded arrows.
           \nA supplied vector must match the style of dagitty edges, or provide valid entries in 'edges_to_keep' for the number of edges to be assessed.")

            return(edges_vec)

          }
        }

        edges_vec <- na.omit(edge_name)

        check_ans <- TRUE
      }
    }
  }

  while(check_ans == FALSE){
    check_skip_sequence <- TRUE

    edges_vec <- edges_to_keep[edges_to_keep %in% edges_vec]

    if(length(edges_to_keep) != length(edges_vec)){

      stop("\nEntries in 'edges_to_keep' are not the correct syntax e.g. 'x -> y', or do not match 'y', 'Y', 'yes', TRUE or 'n', 'N', 'no', FALSE for kept/discarded arrows.
           \nA supplied vector must match the style of dagitty edges, or provide valid entries in 'edges_to_keep' for the number of edges to be assessed.")

      return(edges_vec)

    }

    check_ans <- TRUE
  }

  edges_vec <- ESC_DAGs_sequence(edges_vec, treatment_or_outcome_edges_vec, num_edges, check_skip_sequence)

  if( length(edges_vec[!treatment_or_outcome_edges_vec %in% edges_vec]) != 0){

    message("\nOutputting edges to be assessed.")

    return(edges_vec)
  }


  node_names <- names(dag)
  exclude_names <-   c(treatment, outcome, latent_vec)
  var_names <- node_names[!node_names %in% exclude_names]

  coordinates <- dagitty::coordinates(dag)
  node_name_and_coords_vec <- c()
  if(length(latent_vec) > 0) {
    node_name_and_coords_vec <- c(paste(dagitty::exposures(dag), " [exposure] ", sep=""),
                                  paste(dagitty::outcomes(dag), " [outcome] ", sep=""),
                                  paste(dagitty::latents(dag), " [latent] ", sep=""),
                                  paste(var_names, collapse=" "))
  }else{
    node_name_and_coords_vec <- c(paste(dagitty::exposures(dag), " [exposure] ", sep=""),
                                  paste(dagitty::outcomes(dag), " [outcome] ", sep=""),
                                  paste(var_names, collapse=" "))
  }

  dag <- paste("dag {", paste(node_name_and_coords_vec, collapse=""), paste(edges_vec, collapse=" "), "}", sep = " ")
  dag <- dagitty::dagitty(dag)

  if(all(!is.na(unlist(coordinates)))){
    dagitty::coordinates(dag) <- coordinates
  }

  if(dagitty::isAcyclic(dag) == FALSE){
    warning("The output is not a directed acyclic graph (DAG). Relationships may need to be further assessed.")
  }

  return(dag)

}

#' ESC-DAGs causal criteria for removinng edges
#' @param edges_vec vector of edges whose relationships are to be assessed
#' @param treatment_or_outcome_edges_vec vector of all other edges supplied by the dag
#' @param num_edges number of edges to be assessed
#' @param check_skip_sequence TRUE or FALSE depending on prior inputs
#' @noRd
ESC_DAGs_sequence <- function(edges_vec, treatment_or_outcome_edges_vec, num_edges, check_skip_sequence){

  if(check_skip_sequence == TRUE){

    cat("\nSkipped ESC-DAGs protocol.", "\n")

    edges_vec <- c(treatment_or_outcome_edges_vec, edges_vec)

    return(edges_vec)

  }else if(check_skip_sequence == FALSE) {

    removed_arrows <- c()
    arrow_count <- 1
    num_arrow_to_remove <- 0

    check_ans <- FALSE

    while(check_ans == FALSE){

      choice <- readline("(y/n/?info): ")

      if(choice == "y"){

        #cat("\n", "\n===================================")
        #cat("\nESC-DAGs (Ferguson et al., 2020)")
        #cat("\nDAGs from background knowledge")
        #cat("\nCode written by AJ Moller (2025)")
        #cat("\n===================================", "\n")

        check_ans <- TRUE

      }else if(choice == "n"){

        cat("\nSee you next time!")

        return(edges_vec)

      }else if(choice == "?info"){

        cat("\nEach directed edge in the IG is assessed for three causal criteria: temporality; face-validity; and recourse to theory. They are primarily informed by the classic Bradford Hill viewpoints,24 and are compatible with the ‘inference to the best explanation’ approach advocated by Krieger and Davey Smith.1 If a relationship is determined to possess each criterion, a counterfactual thought experiment derived from the POF is used to further explicate the reviewers’ assumptions.25 The translation process thus combines ‘classic’ and ‘modern’ causal thinking and understands DAGs as ‘conceptual tools’1 for exploring causation, rather than substitutes for careful causal thinking.", "\n")
        cat("\nThe ESC-DAGs causal criteria operate sequentially, with each criterion designed to elaborate over the previous. If any criterion on the edge is not present, the edge can be deleted. The exception is the recourse to theory criterion—absence of theory in the study or according to the reviewer does not equate to absence of effect.26 The counterfactual thought experiment is performed after assessing all criteria. All retained directed edges are entered into the directed edge index. However, each edge should be tested in both directions (i.e. with the head and tail of the arrow swapped). If the posited and reverse edges are both retained, then the relationship should be noted as bi-directional in the directed edge index. Reviewers can also note low confidence in particular directed edges.", "\n")

      }else{

        cat("\n", "\nPlease type a valid answer.", "\n")

      }
    }

    for(arrow_count in 1:num_edges){

      criterion_num <- 1
      check_arrow_removed <- FALSE

      check_ans <- FALSE

      arrow <- edges_vec[arrow_count]
      cat("\n", "\nFor the directed arrow '", arrow, "' consider each of the following:", sep="")

      cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
      cat("\n[1/4] Temporality: does the variable to the left of the arrow precede the variable on the right?", "\n")
      cat("\nFor help, enter ?info")

      while(check_ans == FALSE){

        choice <- readline("(y/n/?info): ")

        if(choice == "y"){

          criterion_num <- criterion_num + 1

          check_ans <- TRUE

        }else if(choice == "n"){

          num_arrow_to_remove <- num_arrow_to_remove + 1
          removed_arrows[num_arrow_to_remove] <- arrow_count

          check_arrow_removed <- TRUE

          check_ans <- TRUE

          cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

        }else if(choice == "?info"){

          cat("\n", "\nCausal criterion 1—temporality:")
          cat("\nOf the Bradford Hill criteria, temporality is the only one not requiring extensive qualification or not yet disproven. (Thomas et al., 2013; DOI: https://doi.org/10.1146/annurev-publhealth-031811-124606) It states that effect cannot precede cause. For example, in Figure 1(A) (Ferguson et al., 2020; DOI: https://doi.org/10.1093/ije/dyz150), adolescent substance use cannot precede historical parental alcohol use, so the relationship would not be temporal. Unless the directed edge is not temporal, we proceed to causal criterion 2.", "\n")
          cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")

        }else{

          cat("\n", "\nPlease type a valid answer.", "\n")

        }
      }

      check_ans <- FALSE

      if(check_arrow_removed == FALSE){

        cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
        cat("\n[2/4] Face-validity: is the posited relationship plausible?", "\n")
        cat("\nFor help, enter ?info")

        while(check_ans == FALSE){

          choice <- readline("(y/n/?info): ")

          if(choice == "y"){

            criterion_num <- criterion_num + 1

            check_ans <- TRUE

          }else if(choice == "n"){

            num_arrow_to_remove <- num_arrow_to_remove + 1
            removed_arrows[num_arrow_to_remove] <- arrow_count

            check_arrow_removed <- TRUE

            check_ans <- TRUE

            cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

          }else if(choice == "?info"){

            cat("\n", "\nCausal criterion 2—face-validity:")
            cat("\nFace-validity is related to the Bradford Hill criterion of (biologic) plausibility. Nested within the wider causal criteria scheme, the face-validity criterion is a rapid means of using reviewer background knowledge to identify implausible relationships, given the temporality established in criterion 1. For example, in Figure 1(A) it is plausible that directed edges originate from sex, but implausible that historical parental alcohol use could influence adolescent sex assignment despite temporal ordering.", "\n")
            cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")


          }else{

            cat("\n", "\nPlease type a valid answer.", "\n")

          }
        }
      }

      check_ans <- FALSE

      if(check_arrow_removed == FALSE){

        cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
        cat("\n[3/4] Recourse to theory—is the posited relationship supported by theory?", "\n")
        cat("\nFor help, enter ?info")

        while(check_ans == FALSE){

          choice <- readline("(y/n/?info): ")

          if(choice == "y"){

            criterion_num <- criterion_num + 1
            check_ans <- TRUE

          }else if(choice == "n"){

            num_arrow_to_remove <- num_arrow_to_remove + 1
            removed_arrows[num_arrow_to_remove] <- arrow_count

            check_arrow_removed <- TRUE

            check_ans <- TRUE

            cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

          }else if(choice == "?info"){

            cat("\n", "\nCausal criterion 3—recourse to theory:")
            cat("\nThe recourse to theory criterion considers background and expert knowledge more overtly. It subsumes the temporality and face-validity criteria and continues to cement a platform for the counterfactual thought experiment. Where the face-validity criterion is concerned with the researcher’s own knowledge, the step assesses whether there is formal theoretical support for the relationship. The decision log for this criterion requires the reviewer to state briefly what theory applies (if any) with space for a reference. As noted above, lack of theory does not equate to lack of effect. As such the purpose of this criterion is not so much falsification as preparation for the next step.", "\n")
            cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")

          }else{

            cat("\n", "\nPlease type a valid answer.", "\n")

          }
        }
      }

      check_ans <- FALSE

      if(check_arrow_removed == FALSE){

        cat("\n", "\n'", arrow, "' (", arrow_count, "/", num_edges, ")",  sep="")
        cat("\n[4/4] Counterfactual thought experiment: is the posited relationship supported by a systematic thought experiment informed by the potential outcomes framework?")
        cat("\n", "\nFor help, enter ?info")

        while(check_ans == FALSE){

          choice <- readline("(y/n/?info): ")

          if(choice == "y"){

            criterion_num <- criterion_num + 1

            check_ans <- TRUE

          }else if(choice == "n"){

            num_arrow_to_remove <- num_arrow_to_remove + 1
            removed_arrows[num_arrow_to_remove] <- arrow_count

            check_arrow_removed <- TRUE

            check_ans <- TRUE

            cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

          }else if(choice == "?info"){

            cat("\n", "\nCounterfactual thought experiment")
            cat("\nFundamentally, potential outcomes compare the outcome that would have occurred if all of the sample had been exposed, with the outcome that would have occurred if all of the sample had not been exposed.3,4,25 The counterfactual thought experiment employs this heuristic in a formulaic and transparent way, comparing two or more ‘counterfactual exposures’ and considering whether their potential outcomes would be different, given the causal criteria. The original study’s measurement of variables should be emulated. See Ferguson et al. (2020) for more details.", "\n")
            cat("\nSource: Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', DOI: https://doi.org/10.1093/ije/dyz150)")

          }else{

            cat("\n", "\nPlease type a valid answer.", "\n")

          }
        }
      }

      arrow_count <- arrow_count + 1
    }
  }
  if(num_arrow_to_remove > 0){

    edges_vec <- edges_vec[-removed_arrows]

    edges_vec <- c(treatment_or_outcome_edges_vec, edges_vec)

    return(edges_vec)

  }else{

    cat("\nNo arrows were removed.", "\n")

    edges_vec <- c(treatment_or_outcome_edges_vec, edges_vec)

    return(edges_vec)
  }
}

#' causalify a dagitty object
#'
#' causalify()
#'
#' @importFrom data.table data.table fcase
#' @param dag dagitty object
#' @return Model object of either
#' @export
causalify <- function(dag){

  edges <- getEdges(dag)

  var_names <- unique(edges[c("ancestor", "role_ancestor")])

  treatment <- var_names[,"ancestor"][var_names[, "role_ancestor"] %in% "treatment"]
  outcome <- dagitty::outcomes(dag)
  latents <- dagitty::latents(dag)

  proxy_var_name <- edges$ancestor[edges$role_ancestor == "competing_exposure"]

  new_edges <- lapply( 1:length(proxy_var_name), function(x){

    new_edges <- c(ancestor = as.character(paste("U",proxy_var_name[x],sep="_")),
                                                  edge = as.character("->"),
                                                  descendant = as.character(proxy_var_name[x]),
                                                  role_ancestor = as.character("latent"),
                                                  role_descendant =as.character("proxy_c"))

  })


  new_outcome_edges <- lapply( 1:length(proxy_var_name), function(x){
    new_outcome_edges <- c(ancestor = as.character(paste("U",proxy_var_name[x],sep="_")),
                           edge = as.character("->"),
                           descendant = as.character(treatment),
                           role_ancestor = as.character("latent"),
                           role_descendant =as.character("treatment"))
  })
  edges <- dplyr::bind_rows(edges, new_edges, new_outcome_edges)


  edges[,"role_ancestor"][edges[, "ancestor"] %in% proxy_var_name] <- "proxy"
  edges[,"role_descendant"][edges[, "descendant"] %in% proxy_var_name] <- "proxy"

  node_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "latent")
  edges_list <- lapply( seq_along(1:9),
                        function(x){

                          edges_list <- edges %>% filter( role_ancestor == node_roles[x] )
                          edges_list[,c(1:3)]

                        } )


  var_names <- unique(edges[c("ancestor", "role_ancestor")])
  var_names <- var_names[,1]

  exclude_names <-c(treatment, outcome, latents)
  var_names <- var_names[!var_names %in% exclude_names]

  latents <- c(unlist(latents), as.vector(unique(edges[,"ancestor"][edges[, "role_ancestor"] %in% "latent"])))

  coordinates <- dagitty::coordinates(dag)
  node_name_and_coords_vec <- c()

  if(length(latents) > 0) {
    node_name_and_coords_vec <- c(paste(treatment, " [exposure] ", sep=""),
                                  paste(outcome, " [outcome] ", sep=""),
                                  paste(latents, " [latent] ", sep=""),
                                  paste(var_names, collapse=" "))
  }else{
    node_name_and_coords_vec <- c(paste(treatment, " [exposure] ", sep=""),
                                  paste(outcome, " [outcome] ", sep=""),
                                  paste(var_names, collapse=" "))
  }

  edges_unlist <- dplyr::bind_rows(edges_list)

  edges_unlist <- lapply(1:nrow(edges_unlist), function(x){
    edges_unlist <- paste(edges_unlist[x,], collapse=" ")
  })
  edges_vec <- paste(unlist(edges_unlist), collapse=" ")

  dag <- paste("dag {", paste(node_name_and_coords_vec, collapse=""), edges_vec, "}", sep = " ")

  dag <- dagitty::dagitty(dag)

  # sets coordinates for latent nodes
  dag <- add_latent_coordinates(dag, latents = latents, existing_coords = coordinates)


  return(dag)
}

