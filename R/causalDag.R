#' Builds a dagitty object saturated (fully connected) graph
#'
#' buildGraph() produces a saturated graph by default from the supplied treatment, outcome, confounder, and  mediators (optional) based on the type of graph specified. Optional latent confounder or medfiator variables can be specified.
#'
#'
#' @importFrom magrittr %>%
#' @param type The type of graph generated. Defaults to 'full', producing a fully connected graph with confounders connected in both directions (bi-directional), and to mediators in one direction (uni-directional). If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param variables Vector of variables to be assigned nodes in a graph.  A named vector can be supplied, containing any of the other input variable types
#' @param treatment Treatment variable, e.g. "X". Must be specified, unless included in the named vector 'variables'.
#' @param outcome Outcome variable, e.g. "Y". Must be specified, unless included in the named vector 'variables'.
#' @param confounders Vector of confounder variables, e.g. c("Z1", "Z2", "Z3").
#' @param mediators Vector of mediator variables, e.g. c("M1", "M2", "M3").
#' @param latent_variables Vector of already supplied variable names to be labelled latent (unobserved), e.g. c("L1", "L2", "M1").
#' @param instrumental_variables Vector of instrumental variables, e.g. "IV"
#' @param mediator_outcome_confounders Vector of already supplied confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list of vectors can instead be supplied with separate vectors for each confounder name, and their connected mediator(s).
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
#' # Option 3: The 'variables' input can be used to connect all variables to each other, with treatment, outcome and other nodes inputted separately.
#' #           Inputted variables are treated as confounders, while the actual 'confounders' input can be left blank.
#' #           We decided to keep this option, rather than remove it, to provide easier access to users who may not have much 'exposure' to causal graphs and terminology.
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
                       variables = NULL,
                       treatment = NULL,
                       outcome = NULL,
                       confounders = NULL,
                       mediators = NULL,
                       latent_variables = NULL,
                       instrumental_variables = NULL,
                       mediator_outcome_confounders = NA,
                       coords_spec = c(lambda = 0, lambda_max = NA, iterations = NA)){


  # Option 1: named vector of variables inputted (can either include treatment and outcome or as separate inputs)
  if(is.null(confounders) & !is.null(variables) & !is.null(names(unlist(variables)))){

    variables_df <- data.frame(variables = as.vector(unlist(variables)),
                               name = names(unlist(variables)),
                               stringsAsFactors = FALSE)

    variables_df <- variables_df %>%
      dplyr::mutate(role = dplyr::case_when(
        grepl("mediator-outcome", name) ~ "mediator-outcome-confounder",
        grepl("mediator_outcome", name) ~ "mediator-outcome-confounder",
        grepl("mediator outcome", name) ~ "mediator-outcome-confounder",
        grepl("outcome", name) ~ "outcome",
        grepl("treatment", name) ~ "treatment",
        grepl("confounder", name) ~ "confounder",
        grepl("instrumental", name) ~ "instrumental",
        grepl("collider", name) ~ "collider",
        grepl("latent", name) ~ "latent",
        grepl("mediator", name) ~ "mediator",
      )) %>%
      dplyr::mutate(name = dplyr::case_when(
        role == "outcome" ~ gsub("[^0-9.-]", "", name),
        role == "outcomes" ~ gsub("[^0-9.-]", "", name),
        grepl("treatment", role) ~ gsub("[^0-9.-]", "", name),
        role == "confounder" ~ gsub("[^0-9.-]", "", name),
        role == "confounders" ~ gsub("[^0-9.-]", "", name),
        grepl("instrumental", role) ~ gsub("[^0-9.-]", "", name),
        grepl("collider", role) ~ gsub("[^0-9.-]", "", name),
        grepl("latent", role) ~ gsub("[^0-9.-]", "", name),
        role == "mediator" ~ gsub("[^0-9.-]", "", name),
        role == "mediators" ~ gsub("[^0-9.-]", "", name),
        grepl("mediator-outcome", role) ~ gsub("[^0-9.-]", "", name),
      )) %>% dplyr::rename("order" = name)

    confounder_df <- variables_df %>% dplyr::filter(role == "confounder")
    confounder_vec <-  confounder_df$variables

    outcome_df <- variables_df %>% dplyr::filter(role == "outcome")
    outcome <- as.vector(outcome_df["variables"])

    treatment_df <- variables_df %>% dplyr::filter(role == "treatment")
    treatment <- as.vector(treatment_df["variables"])

    mediator_df <- variables_df %>% dplyr::filter(role == "mediator")
    latent_df <- variables_df %>% dplyr::filter(role == "latent")
    instrumental_df <- variables_df %>% dplyr::filter(role == "instrumental")
    m_o_confounder_df <- variables_df %>% dplyr::filter(role == "mediator-outcome-confounder")

    if(nrow(mediator_df) == 0){
      mediator_vec <- NULL
    }else{
      mediator_vec <- as.vector(mediator_df$variables)
    }

    if(nrow(latent_df) == 0){
      latent_vec <- NULL
    }else{
      latent_vec <- as.vector(latent_df$variables)
    }

    if(nrow(instrumental_df) == 0){
      instrumental_vec <- NULL
    }else{
      instrumental_vec <- as.vector(instrumental_df$variables)
    }

    if(nrow(m_o_confounder_df) == 0){
      m_o_confounder_vec <- as.vector(NA)
    }else{
      m_o_confounder_vec <- as.vector(m_o_confounder_df$variables)
    }
  }else{
    mediator_vec <- as.vector(mediators)
    latent_vec <- as.vector(latent_variables)
    instrumental_vec <- as.vector(instrumental_variables)
    m_o_confounder_vec <- as.vector(mediator_outcome_confounders)

    # Option 2: variables inputted (assumed common causes of treatment and outcome),
    #           Separate inputs for treatment, outcome and other nodes while blank confounders input
    if(is.null(confounders) & !is.null(variables)){

      variables_vec <- as.vector(unlist(variables))
      confounder_vec <- variables_vec[!variables_vec %in% mediator_vec & !variables_vec %in% instrumental_vec & !variables_vec %in% m_o_confounder_vec]

      confounder_df <- data.frame(variables = confounder_vec,
                                  role = "confounder",
                                  stringsAsFactors = FALSE)

    }else if(!is.null(confounders)){

      # Option 3: confounders inputted, variables input is ignored.
      #           Separate inputs for treatment, outcome, and other nodes.
      confounder_vec <- as.vector(unlist(confounders))

      confounder_df <- data.frame(variables = confounder_vec,
                                  role = "confounder",
                                  stringsAsFactors = FALSE)
    }


    # mediator variables data frame
    mediator_df <- data.frame(variables = mediator_vec,
                              stringsAsFactors = FALSE)

  }



  if( ( is.null(treatment) | is.null(outcome) ) & is.null(names(unlist(variables)))){

    stop("Missing treatment and/or outcome inputs, or the supplied variables are not labelled. Both treatment and outcome must be given as inputs or supplied in a named vector for graphs using buildGraph().")

  }

  saturated_vec <- c()
  treatment_vec <- c()
  outcome_vec <- c()

  confounder_df$moc <- FALSE

  m_o_confounder <- 1
  confounder <- 1
  arrow_count <- 1
  treatment_count <- 1
  outcome_count <- 1

  for(m_o_confounder in 1:length(m_o_confounder_vec)){

    if(all(!is.na(m_o_confounder_vec))){

      # if all mediator-outcome confounders are not confounders, edges are added to the vector
      if( length(confounder_vec[!confounder_vec %in% m_o_confounder_vec[m_o_confounder]]) == length(confounder_vec) ){

        outcome_vec[outcome_count] <- paste(m_o_confounder_vec[m_o_confounder], "->", outcome)
        outcome_count <- outcome_count + 1

        if(type == "full"){
          treatment_vec[treatment_count] <- paste(treatment, "->", m_o_confounder_vec[m_o_confounder])
          treatment_count <- treatment_count + 1


        }


        if(all(!is.null(mediator_vec))){

          mediator_count <- 1

          for(mediator_count in 1:length(mediator_vec)){

            saturated_vec[arrow_count] <- paste(m_o_confounder_vec[m_o_confounder], "->", mediator_vec[mediator_count])
            arrow_count <- arrow_count + 1
            mediator_count <- mediator_count + 1
          }
        }

      }

    }


    confounder <- 1

    for(confounder in 1:length(confounder_vec)){

      # if mediator-outcome confounders are supplied, each confounder is checked and assigned
      if(all(!is.na(m_o_confounder_vec))){
        if(confounder_vec[confounder] ==  m_o_confounder_vec[m_o_confounder]){
          confounder_df[confounder,"moc"] <- TRUE
        }else{
          treatment_vec[treatment_count] <- paste(confounder_vec[confounder], "->", treatment)
          outcome_vec[outcome_count] <- paste(confounder_vec[confounder], "->", outcome)

          treatment_count <- treatment_count + 1
          outcome_count <- outcome_count + 1
        }
      }else{
        treatment_vec[treatment_count] <- paste(confounder_vec[confounder], "->", treatment)
        outcome_vec[outcome_count] <- paste(confounder_vec[confounder], "->", outcome)

        treatment_count <- treatment_count + 1
        outcome_count <- outcome_count + 1
      }

      saturated <- 1


      if(confounder_df[confounder,"moc"] == FALSE){
        if(all(!is.null(mediator_vec))){

          mediator_count <- 1

          for(mediator_count in 1:length(mediator_vec)){

            outcome_vec[outcome_count] <- paste(mediator_vec[mediator_count], "->", outcome)

            saturated_vec[arrow_count] <- paste(treatment, "->", mediator_vec[mediator_count])

            outcome_count <- outcome_count + 1
            arrow_count <- arrow_count + 1
            saturated <- 1
            sat_mediator <- 1

            for(saturated in 1:length(confounder_vec)){

              if(!is.na(confounder_vec[saturated])){

                if(type == "full" | type == "saturated"){

                  if(confounder_df[saturated,"variables"] != confounder_df[confounder,"variables"]){

                    saturated_vec[arrow_count] <- paste(confounder_df[confounder,"variables"], "->", confounder_df[saturated,1])

                    arrow_count <- arrow_count + 1
                  }

                  if(type == "full"){

                    saturated_vec[arrow_count] <- paste(confounder_vec[saturated], "->", mediator_vec[mediator_count])

                    arrow_count <- arrow_count + 1
                  }

                  saturated <- saturated + 1

                }else if(type == "ordered"){

                  order_occurrance <- as.numeric(order(match(confounder_df$variables, confounder_df$variables)))

                  if(confounder_vec[saturated] != confounder_vec[confounder] & order_occurrance[saturated] > order_occurrance[confounder]){

                    saturated_vec[arrow_count] <- paste(confounder_vec[confounder], "->", confounder_vec[saturated])

                    saturated <- saturated + 1
                    arrow_count <- arrow_count + 1
                  }
                }
              }

            }

            mediator_count <- mediator_count + 1
          }

        }else if(all(is.null(mediator_vec))){

          for(saturated in 1:nrow(confounder_df)){

            if(type == "full" | type == "saturated"){

              if(confounder_df[saturated,"variables"] != confounder_df[confounder,"variables"]){

                saturated_vec[arrow_count] <- paste(confounder_df[confounder,"variables"], "->", confounder_df[saturated,"variables"])

                saturated <- saturated + 1
                arrow_count <- arrow_count + 1
              }
            }else if(type == "ordered"){

              order_occurrance <- as.numeric(order(match(confounder_df$variables, confounder_df$variables)))

              if(confounder_vec[saturated] != confounder_vec[confounder] & order_occurrance[saturated] > order_occurrance[confounder]){

                saturated_vec[arrow_count] <- paste(confounder_vec[confounder], "->", confounder_vec[saturated])

                saturated <- saturated + 1
                arrow_count <- arrow_count + 1
              }
            }
          }
        }

      }else if(confounder_df[confounder,"moc"] == TRUE){
        if(all(!is.null(mediator_vec))){

          mediator_count <- 1

          for(mediator_count in 1:length(mediator_vec)){

            outcome_vec[outcome_count] <- paste(mediator_vec[mediator_count], "->", outcome)
            saturated_vec[arrow_count] <- paste(treatment, "->", mediator_vec[mediator_count])

            outcome_count <- outcome_count + 1
            arrow_count <- arrow_count + 1
            saturated <- 1

            if(type == "full"){
              for(saturated in 1:length(confounder_vec)){

                saturated_vec[arrow_count] <- paste(confounder_vec[saturated], "->", mediator_vec[mediator_count])

                saturated <- saturated + 1
                arrow_count <- arrow_count + 1
              }
            }

            mediator_count <- mediator_count + 1
          }
        }
      }

      confounder <- confounder + 1
    }

    m_o_confounder <- m_o_confounder +1
  }

  instrumental_treatment_vec <- c()

  if(all(!is.null(instrumental_vec))){

    instrumental <- 1

    for(instrumental in 1:length(instrumental_vec)){
      instrumental_treatment_vec[instrumental] <- paste(instrumental_vec[instrumental], "->", treatment)

      instrumental <- instrumental + 1
    }

  }



  if(all(!is.null(latent_vec))){

    latents <- 1

    for(latents in 1:length(latent_vec)){
      if(length(latent_vec[latent_vec[latents] %in% mediator_vec]) > 0){

        mediator_count <- 1

        for(mediator_count in 1:length(mediator_vec)){
          if(latent_vec[latents] == mediator_vec[mediator_count]){
            mediator_vec[mediator_count] <- NA
          }
          mediator_count <- mediator_count + 1
        }

        mediator_vec <- na.omit(mediator_vec)

      }
      if(length(latent_vec[latent_vec[latents] %in% confounder_vec]) > 0){

        confounder <- 1

        for(confounder in 1:length(confounder_vec)){
          if(latent_vec[latents] == confounder_vec[confounder]){
            confounder_vec[confounders] <- NA
          }
          confounder <- confounder + 1
        }

        confounder_vec <- na.omit(confounder_vec)
      }

      latent_vec[latents] <- paste(latent_vec[latents], "[latent]")
      latents <- latents + 1
    }
  }


  treatment_to_outcome <- paste(treatment, "->", outcome)

  treatment <- paste(treatment, "[exposure]")
  outcome <- paste(outcome, "[outcome]")

  outcome_vec <- na.omit(outcome_vec)


  treatment_vec <- paste(treatment_vec, collapse=" ")
  outcome_vec <- paste(outcome_vec, collapse=" ")
  saturated_vec <- paste(saturated_vec, collapse=" ")
  confounder_vec <- paste(confounder_vec, collapse=" ")
  mediator_vec <- paste(mediator_vec, collapse=" ")
  instrumental_treatment_vec <- paste(instrumental_treatment_vec, collapse=" ")
  latent_vec <- paste(latent_vec, collapse=" ")

  dag_code <- paste("dag {", treatment, outcome, latent_vec, confounder_vec, mediator_vec, instrumental_vec, treatment_to_outcome, treatment_vec, instrumental_treatment_vec, outcome_vec, saturated_vec, "}", sep = " ")

  dag <- dagitty::dagitty(dag_code)

  dag <- getCoords(dag)

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


#' Gets edges between nodes and their roles
#'
#' getEdges() filters a dagitty object and returns a data frame with edges for specified node roles.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @param selected_nodes Nodes to return edges. Defaults to NULL, or can be a character or vector combination of any of the following: "treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", or "instrumental".
#' @returns A dataframe containing the specified node edges.
#' @export
getEdges <- function(dag, selected_nodes = NULL){

  edges <- (dagitty::edges(dag))[1:3]
  edges <- edges %>% dplyr::relocate(e, .before = w)

  edges$role_v <- NA
  edges$role_w <- NA

  edges$latent_v <- FALSE
  edges$latent_w <- FALSE

  # treatment
  exposure <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, exposure)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes)]

  # instrumental variables
  instrumental_vars <- as.vector(unlist(dagitty::instrumentalVariables(dag)))

  # labelling instrumental variable nodes in edges
  if(length(instrumental_vars) > 0){

    instrumental <- 1

    for(instrumental in 1:length(instrumental_vars)){

      iv_edges <- edges %>% dplyr::filter(v == instrumental_vars[instrumental])

      edges <- edges %>%
        mutate(role_v = case_when(v == iv_edges$v ~ "instrumental")
        )

      instrumental <- instrumental + 1
    }

  }


  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, exposure)]

  # mediator-outcome confounders
  mediators_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  mediator_outcome_confounders <- mediators_parents[mediators_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  mediator_outcome_confounders <- mediators_parents[!mediators_parents %in% c(exposure, confounders)] # remove treatment and confounder nodes
  mediator_outcome_confounders <- mediator_outcome_confounders[!mediator_outcome_confounders %in% treatment_parents] # double check by removing parents of treatment

  # mediators - second parse (removing mediator-outcome confounders)
  mediators <- mediators[!mediators %in% mediator_outcome_confounders]


  # labelling each node in edges
  edges$role_v[is.na(edges$role_v) & edges$v %in% exposure] <- "treatment"
  edges$role_w[is.na(edges$role_w) & edges$w %in% exposure] <- "treatment"

  edges$role_v[is.na(edges$role_v) & edges$v %in% outcomes] <- "outcome"
  edges$role_w[is.na(edges$role_w) & edges$w %in% outcomes] <- "outcome"

  edges$role_v[is.na(edges$role_v) & edges$v %in% mediator_outcome_confounders] <- "mediator-outcome-confounder"
  edges$role_w[is.na(edges$role_w) & edges$w %in% mediator_outcome_confounders] <- "mediator-outcome-confounder"

  edges$role_v[is.na(edges$role_v) & edges$v %in% mediators] <- "mediator"
  edges$role_w[is.na(edges$role_w) & edges$w %in% mediators] <- "mediator"

  edges$role_v[is.na(edges$role_v) & edges$v %in% confounders] <- "confounder"
  edges$role_w[is.na(edges$role_w) & edges$w %in% confounders] <- "confounder"

  edges$latent_v[edges$v %in% latent_vars] <- TRUE
  edges$latent_w[edges$w %in% latent_vars] <- TRUE


  selected_remove <- NULL

  if(any(grepl("!", selected_nodes))){
    # if any input includes "!", removes edges with matching roles in selected_nodes (for both parent and children nodes)
    all_nodes <- c("treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", "instrumental")
    cleaned_nodes <- gsub("!", "", selected_nodes)

    selected_remove <- cleaned_nodes
    selected_nodes <- all_nodes[!all_nodes %in% cleaned_nodes]

  }else if(any(grepl(">", selected_nodes))){
    # if any input includes ">", the function removes all edges without roles that match the selected_nodes input
    all_nodes <- c("treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", "instrumental")
    cleaned_nodes <- gsub(">", "", selected_nodes)

    selected_remove <- all_nodes[!all_nodes %in% cleaned_nodes]
    selected_nodes <- cleaned_nodes

  }

  selected <- 1

  if(!is.null(selected_remove)){
    # if selected_nodes input does not contain "!" or ">", then for each selected_nodes input the edges containing a matching parent or children node are kept.
    for(selected in 1:length(selected_remove)){

      edges <- edges %>% dplyr::filter(role_v != selected_remove[selected] & role_w != selected_remove[selected])

      selected <- selected + 1

    }



  }else if(all(!is.null(selected_nodes))){
    # if selected_nodes is not inputted, the function returns all edges with node labels in the 'role_v' and 'role_w' columns.
    edges_out <- data.frame(v = character(),
                            e = character(),
                            w = character(),
                            role_v = character(),
                            role_w = character(),
                            stringsAsFactors = FALSE)

    for(selected in 1:length(selected_nodes)){

      selected_edges <- edges %>% dplyr::filter(role_v == selected_nodes[selected] | role_w == selected_nodes[selected])

      selected_edges <- suppressMessages(dplyr::anti_join(selected_edges, edges_out))

      edges_out <- rbind(edges_out, selected_edges)

      selected <- selected + 1

    }

    edges <- edges_out

  }



  edges <- edges %>% dplyr::relocate(role_v, .before = role_w)



  return(edges)
}

#' Mediator names in a dagitty object
#'
#' mediators() is a dagitty wrapper that identifies nodes along paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of mediators, or dataframe of edges.
#' @export
mediators <- function(dag, edges = FALSE){

  # treatment
  exposure <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # mediators
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, exposure)]

  if(edges == TRUE){
    cat(paste("There are", length(mediators), "mediators in the supplied graph: ", sep = " ", collapse = " "))
    cat(paste("\n", mediators, collapse = "\n"))

    edges <- getEdges(dag, "mediator")

    return(edges)
  }

  return(mediators)
}


#' Confounder names in a dagitty object
#'
#' confounders() is a dagitty wrapper that identifies common causes of treatment and outcome in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of confounders, or dataframe of edges.
#' @export
confounders <- function(dag, edges = FALSE){

  # treatment
  exposure <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, exposure)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes)]

  if(edges == TRUE){
    cat(paste("There are", length(confounders), "confounders in the supplied graph: ", sep = " ", collapse = " "))
    cat(paste("\n", confounders, collapse = "\n"))

    edges <- getEdges(dag, "confounder")

    return(edges)
  }

  return(confounders)
}

#' Mediator-outcome confounder names in a dagitty object
#'
#' moc() is a dagitty wrapper that identifies nodes along paths between treatment and outcome in a directed acyclic graph.
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns A vector of mediator-outcome confounders, or dataframe of edges.
#' @export
moc <- function(dag, edges = FALSE){

  # treatment
  exposure <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # moc
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, exposure)]

  # mediator-outcome confounders
  mediators_parents <- dagitty::parents(dag, mediators)
  moc <- mediators_parents[!mediators_parents %in% c(mediators, exposure)]


  if(edges == TRUE){
    cat(paste("There are", length(moc), "mediator-outcome confounders in the supplied graph: ", sep = " ", collapse = " "))
    cat(paste("\n", moc, collapse = "\n"))

    edges <- getEdges(dag, "mediator-outcome-confounder")

    return(edges)
  }

  return(moc)
}

#' Variable names in a dagitty object
#'
#' variables() returns a named vector of variables corresponding to nodes in a dag, providing an easy way to extract node roles e.g. "treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", or "instrumental".
#'
#' @param dag A dagitty object.
#' @returns A named vector of nodes in a dag.
#' @export
variables <- function(dag, all_info = TRUE){

  variables_df <- unique(getEdges(dag)[c("v", "role_v", "latent_v")])

  if(all_info == TRUE){

    variables <- dplyr::rename(variables_df, variable_names = "v", variable_role = "role_v", latent_indicator = "latent_v")

    return(variables)

  }

  variables <- variables_df$v

  #names(variables) <- variables_df$role_v

  return(variables)
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



#' Treatment names from a dagitty object
#'
#' treatments() identifies treatment nodes in a directed acyclic graph dagitty object.
#'
#' @param dag A dagitty object.
#' @returns Name of treatment node(s), as a character vector.
#' @noRd
treatments <- function(dag){


  dag_vec <- as.vector(dag)
  dag_vec <- unlist(strsplit(dag_vec, "exposure"))

  treatments <- c()

  treatment_count <- 1
  for(treatment in 1:( length(dag_vec) -1 )){

    treatments[treatment_count] <-  rev(sub(".*[(.*?)\\\b.*", "\\1", rev(dag_vec[treatment_count])))

  }

  treatments <- substr(treatments, 1, nchar(treatments)-1)

  return(treatments)
}


#' Outcome names from a dagitty object
#'
#' outcomes() identifies treatment nodes in a directed acyclic graph dagitty object.
#'
#' @param dag A dagitty object.
#' @returns Name of outcomes node(s), as a character vector.
#' @noRd
outcomes <- function(dag){


  dag_vec <- as.vector(dag)
  dag_vec <- unlist(strsplit(dag_vec, "outcome"))

  outcomes <- c()

  outcome_count <- 1
  for(outcome in 1:( length(dag_vec) -1 )){

    outcomes[outcome_count] <-  rev(sub(".*[(.*?)\\\b.*", "\\1", rev(dag_vec[outcome_count])))

  }

  outcomes <- substr(outcomes, 1, nchar(outcomes)-1)

  return(outcomes)
}



