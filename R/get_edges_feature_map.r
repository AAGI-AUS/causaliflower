#' Convert dagitty object to list
#'
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @export
getFeatureMap <- function(dag, include_descendant_nodes = FALSE){
  .datatable.aware <- TRUE

  edges <- getEdges(dag, output_format = "default")
  feature_map_list <- getEdges(dag, output_structure = "list")

  names(feature_map_list) <-  c("outcome",
                                "treatment",
                                "confounder",
                                "mediator",
                                "mediator_outcome_confounder",
                                "instrumental")

  feature_map_new <- c()

  feature_map_new <- lapply( seq_along( names(feature_map_list) ), function(role){

    feature_map_list[[role]] <- list(ancestor = unique(feature_map_list[[role]]$ancestor))

    feature_map_list <- lapply(feature_map_list[[role]], function(x) if(identical(x, character(0))) NA_character_ else x)


  })

  names(feature_map_new) <-  c("outcome",
                               "treatment",
                               "confounder",
                               "mediator",
                               "mediator_outcome_confounder",
                               "instrumental")

  outcomes <- unique( edges[,"descendant"][ edges[,"role_descendant"] == "outcome" ] )

  if(is.null(nrow(feature_map_new[1]))){

    feature_map_new[[1]] <- list(ancestor = outcomes)

  }


  if(include_descendant_nodes == TRUE){

    nestedFeatureMap <- function(feature_map_list){

      num_roles <- length(names(feature_map_list))

      role <- 1

      for(role in 1:num_roles){

        num_cols <- length(as.data.frame(feature_map_list[role]))

        unique_features <- unique(feature_map_list[[role]][[1]])
        num_features <- length(unique_features)
        #num_edges <- nrow(as.data.frame(feature_map_list[role]))

        feature_map_new <- list()

        feature <- 1

        if(num_cols > 2) {


          for(feature in 1:num_features){

            if(all(!is.na(feature_map_list[[role]]))){

              ancestor <- unique_features[feature]
              descendants <- feature_map_list[[role]][[3]][  feature_map_list[[role]][[1]] %in% ancestor]

              feature_map_new[feature] <- list( features = list( ancestor = ancestor,
                                                                 descendant = descendants ) )

            }

            feature <- feature + 1

          }

          feature_map_list[[role]] <- feature_map_new

        }

        role <- role + 1

      }

      return(feature_map_list)

    }

    nestedFeatureMap(feature_map_list)

    feature_map_new <- nestedFeatureMap(feature_map_list)

  }


  return(feature_map_new)

}

#' Gets edges between nodes and their roles
#'
#' getEdges() filters a dagitty object and returns a data frame with edges for specified node roles.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @param selected_nodes Nodes to return edges. Defaults to NULL, or can be a character or vector combination of any of the following: c("treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", "instrumental")
#' @param output_format Options for outputted data format: "default" returns a data frame with separate columns for ancestor and descendant roles, "wide" returns separate columns for each unique role for both ancestor and descendant nodes ("list" returns the same wide data format except as a list), and "long" assigns separate rows to ancestor and descendant nodes.
#' @param output_structure Outputted data can be a "data.table", "data.frame", or "list".
#' @returns A data frame, data table, or list of edges for the roles specified in selected_nodes.
#' @export
getEdges <- function(dag, selected_nodes = c("outcome", "treatment", "confounder", "mediator", "mediator-outcome-confounder", "instrumental"), variables = NULL, output_format = "default", output_structure = "data.table"){
  .datatable.aware <- TRUE

  edges <- data.table::as.data.table(dagitty::edges(dag))[,1:3]
  edges <- edges %>% dplyr::relocate(e, .before = w)

  if( identical(selected_nodes, c("outcome", "treatment", "confounder", "mediator", "mediator-outcome-confounder", "instrumental")) ){
    # seems unnecessary but putting this if-statement first prevents unnecessary checks, assuming the selected_nodes input is predominantly left blank
    selected_remove <- NULL

  }else if(any(grepl("!", selected_nodes))){
    # if any input includes "!", removes edges with matching roles in selected_nodes (for both parent and children nodes)
    all_nodes <- c("outcome", "treatment", "confounder", "mediator", "mediator-outcome-confounder", "instrumental")
    cleaned_nodes <- gsub("!", "", selected_nodes)

    selected_remove <- cleaned_nodes
    selected_nodes <- all_nodes[!all_nodes %in% cleaned_nodes]

  }else if(any(grepl(">", selected_nodes))){
    # if any input includes ">", the function removes all edges without roles that match the selected_nodes input
    all_nodes <- c("outcome", "treatment", "confounder", "mediator", "mediator-outcome-confounder", "instrumental")
    cleaned_nodes <- gsub(">", "", selected_nodes)

    selected_remove <- all_nodes[!all_nodes %in% cleaned_nodes]
    selected_nodes <- cleaned_nodes
    is.character(selected_nodes)
  }else if(is.character(selected_nodes) & length(selected_nodes) < 7){
    # this ensures the order of variable names in selected_nodes is consistent
    all_nodes <- c("outcome", "treatment", "confounder", "mediator", "mediator-outcome-confounder", "instrumental")

    selected_remove <- NULL
    selected_nodes <- all_nodes[all_nodes %in% selected_nodes]
  }else{

    stop("Please check the selected_nodes input and try again.")

  }

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes)]

  # instrumental variables
  instrumental_vars <- as.vector(unlist(dagitty::instrumentalVariables(dag)))

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  mediators <- outcome_parents[outcome_parents %in% dagitty::children(dag, treatments)]

  # mediator-outcome confounders
  mediators_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediators_parents[mediators_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- mediators_parents[!mediators_parents %in% c(treatments, confounders)] # remove treatment and confounder nodes
  moc <- moc[!moc %in% treatment_parents] # double check by removing parents of treatment

  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, v] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, v] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, v] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, v] %in% moc]  <- "mediator-outcome-confounder"

  edges$ancestor_mediator[edges[, v] %in% mediators & !edges[, v] %in% moc]  <- "mediator"

  edges$ancestor_iv[edges[, v] %in% instrumental_vars] <- "instrumental"

  edges$ancestor_latent[edges$v %in% latent_vars] <- NA

  # assign roles for all v in edges
  edges$descendant_outcome[edges[, w] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, w] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, w] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, w] %in% moc & !edges[, w] %in% mediators]  <- "mediator-outcome-confounder"

  edges$descendant_mediator[edges[, w] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, w] %in% instrumental_vars] <- "instrumental"

  edges$descendant_latent[edges$w %in% latent_vars] <- NA


  if(output_format == "default"){

    edges_ancestors <- edges[,1:10]
    edges_descendants <- edges[,c(1:3,11:17)]

    edges_ancestors <- na.omit(tidyr::pivot_longer(edges_ancestors, 4:10))[,c(1:3,5)]

    edges_descendants <- na.omit(tidyr::pivot_longer(edges_descendants, 4:10))[,5]

    edges <- cbind(edges_ancestors, edges_descendants)

    colnames(edges) <- c("ancestor",
                         "edge" ,
                         "descendant",
                         "role_ancestor",
                         "role_descendant")

    edges_list <- lapply( seq_along(selected_nodes),
                          function(x){

                            edges_list <- edges %>% filter(role_ancestor == selected_nodes[x])

                          } )

    if(output_structure == "list"){

      return(edges_list)
    }

    edges <-  Reduce(function(x,y) merge(x,y,all=T),edges_list)

    return(edges)

  }


  if(output_format == "wide"){

    edges_list <- lapply( seq_along(selected_nodes),
                          function(x){

                            edges_list <- edges %>% filter(  ancestor_treatment == selected_nodes[x] |
                                                               ancestor_outcome == selected_nodes[x] |
                                                               ancestor_confounder == selected_nodes[x] |
                                                               ancestor_moc == selected_nodes[x] |
                                                               ancestor_mediator == selected_nodes[x] |
                                                               ancestor_iv == selected_nodes[x]  )

                            edges_list <- edges_list[,which(unlist(lapply(edges_list, function(y)!all(is.na(y))))),with=F]

                          } )

    if(output_structure == "list"){

      return(edges_list)
    }

    edges_list <- edges_list[sapply(edges_list, function(x) dim(x)[1]) > 0]

    edges_wide <- ( dplyr::bind_rows(edges_list, id = edges_list) )[1:nrow(edges),]

    names(edges_wide)[names(edges_wide) == 'v'] <- 'ancestor'
    names(edges_wide)[names(edges_wide) == 'e'] <- 'edge'
    names(edges_wide)[names(edges_wide) == 'w'] <- 'descendant'

    return(edges_wide)


  }


  if(output_format == "long"){

    # output based on output format
    edges_long <- na.omit(tidyr::pivot_longer(edges, 4:17))

    colnames(edges_long) <- c("ancestor",
                              "edge" ,
                              "descendant",
                              "ancestor_descendant",
                              "role")

    edges_list <- lapply( seq_along(selected_nodes),
                          function(x){

                            edges_list <- edges_long %>% filter(role == selected_nodes[x])

                          }
    )

    edges_list <- lapply(edges_list, function(x) cbind(x))

    if(output_structure == "list"){

      return(edges_list)
    }

    edges_long <- Reduce(function(x,y) merge(x,y,all=T),edges_list)

    edges_long$ancestor_descendant <- sub(paste0("_", ".*"), "", edges_long$ancestor_descendant)

    return(edges_long)

  }

  stop("Invalid 'output_format' - check input and try again.")

}

