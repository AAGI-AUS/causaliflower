
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
getEdges <- function(dag, selected_nodes = c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent"), variables = NULL, output_format = "default", output_structure = "data.table"){
  .datatable.aware <- TRUE

  all_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent")

  edges_dagitty <- data.table::as.data.table(dagitty::edges(dag))[,1:3]
  edges <- edges_dagitty %>% dplyr::relocate(e, .before = w)

  if( identical(selected_nodes, all_roles) ){
    # seems unnecessary but putting this if-statement first prevents unnecessary checks, assuming the selected_nodes input is predominantly left blank
    selected_remove <- NULL

  }else if(any(grepl("!", selected_nodes))){
    # if any input includes "!", removes edges with matching roles in selected_nodes (for both parent and children nodes)
    all_nodes <- all_roles
    cleaned_nodes <- gsub("!", "", selected_nodes)

    selected_remove <- cleaned_nodes
    selected_nodes <- all_nodes[!all_nodes %in% cleaned_nodes]

  }else if(any(grepl(">", selected_nodes))){
    # if any input includes ">", the function removes all edges without roles that match the selected_nodes input
    all_nodes <- all_roles
    cleaned_nodes <- gsub(">", "", selected_nodes)

    selected_remove <- all_nodes[!all_nodes %in% cleaned_nodes]
    selected_nodes <- cleaned_nodes
    is.character(selected_nodes)
  }else if(is.character(selected_nodes) & length(selected_nodes) < 7){
    # this ensures the order of variable names in selected_nodes is consistent
    all_nodes <- all_roles

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
  treatment_children <- dagitty::children(dag, treatments)
  mediators <- outcome_parents[outcome_parents %in% treatment_children]

  # mediator-outcome confounders
  mediators_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediators_parents[mediators_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- mediators_parents[!mediators_parents %in% c(treatments, confounders)] # remove treatment and confounder nodes
  moc <- moc[!moc %in% treatment_parents] # double check by removing parents of treatment

  # competing exposure
  competing_exposure <- outcome_parents[ !outcome_parents %in% mediators &
                                         !outcome_parents %in% treatments &
                                         !outcome_parents %in% confounders &
                                         !outcome_parents %in% latent_vars &
                                         !outcome_parents %in% moc ]

  # proxy

  latent_descendants <- dagitty::children(dag, latent_vars)
  proxy_b <- treatment_parents[ treatment_parents %in% latent_descendants] # proxy_b
  proxy_c <- outcome_parents[ outcome_parents %in% latent_descendants & # proxy_c
                                  !outcome_parents %in% treatments &
                              !outcome_parents %in% mediators ]
  proxy <- c(proxy_b, proxy_c)

  # collider
  outcome_children <- dagitty::children(dag, outcomes)
  colliders <- outcome_children[outcome_children %in% treatment_children]


  # assign roles for all v in edges
  edges$ancestor_outcome[edges[, v] %in% outcomes]  <- "outcome"

  edges$ancestor_treatment[edges[, v] %in% treatments]  <- "treatment"

  edges$ancestor_confounder[edges[, v] %in% confounders]  <- "confounder"

  edges$ancestor_moc[edges[, v] %in% moc]  <- "mediator_outcome_confounder"

  edges$ancestor_mediator[ edges[, v] %in% mediators &
                             !edges[, v] %in% moc ]  <- "mediator"

  edges$ancestor_iv[edges[, v] %in% instrumental_vars] <- "instrumental"

  edges$ancestor_competing_exposure[edges[, v] %in% competing_exposure &
                                      !edges[, v] %in% proxy_c ] <- "competing_exposure"

  edges$ancestor_proxy[ ( edges[, v] %in% proxy_b &
                         !edges[, v] %in% confounders ) | # "proxy_b"
                        ( edges[, v] %in% proxy_c &
                         !edges[, v] %in% confounders &
                         !edges[, v] %in% mediators ) ] <- "proxy"

  edges$ancestor_collider[edges[, v] %in% colliders] <- "collider"

  edges$ancestor_latent[edges$v %in% latent_vars] <- NA

  # assign roles for all v in edges
  edges$descendant_outcome[edges[, w] %in% outcomes]  <- "outcome"

  edges$descendant_treatment[edges[, w] %in% treatments]  <- "treatment"

  edges$descendant_confounder[edges[, w] %in% confounders]  <- "confounder"

  edges$descendant_moc[edges[, w] %in% moc &
                         !edges[, w] %in% mediators]  <- "mediator_outcome_confounder"

  edges$descendant_mediator[edges[, w] %in% mediators]  <- "mediator"

  edges$descendant_iv[edges[, w] %in% instrumental_vars] <- "instrumental"

  edges$descendant_competing_exposure[edges[, w] %in% competing_exposure &
                                        !edges[, w] %in% proxy_c ] <- "competing_exposure"

  edges$descendant_proxy[ ( edges[, w] %in% proxy_b &
                           !edges[, w] %in% confounders ) | # "proxy_b"
                          ( edges[, w] %in% proxy_c &
                           !edges[, w] %in% confounders &
                           !edges[, w] %in% mediators ) ] <- "proxy"

  edges$descendant_collider[edges[, w] %in% colliders] <- "collider"

  edges$descendant_latent[edges$w %in% latent_vars] <- NA

  if(output_format == "default"){

    edges_ancestors <- edges[,1:13]
    edges_descendants <- edges[,c(1:3,14:23)]

    edges_ancestors <- stats::na.omit(tidyr::pivot_longer(edges_ancestors, 4:13))#[,c(1:3,5)]
    edges_descendants <- stats::na.omit(tidyr::pivot_longer(edges_descendants, 4:13))#[,c(1:3,5)]

    if(nrow(edges_ancestors) != nrow(edges_descendants)){


      missing_ancestors <- dplyr::anti_join(edges_descendants, edges_ancestors, by = c("v" = "v"))

      missing_latent_rows <- NULL


      missing_latent_rows <- lapply( 1:nrow(missing_ancestors), function(x){
        missing_row <- missing_ancestors[x, ]
        missing_row[,"value"][missing_row[, 1] %in% latent_vars] <- "latent"
        missing_row <- as.data.frame(missing_row)

        })

      edges_ancestors <- dplyr::bind_rows(edges_ancestors, missing_latent_rows)

    }

    edges <- cbind(edges_ancestors[,c(1:3,5)], edges_descendants[,5])

    colnames(edges) <- c("ancestor",
                         "edge" ,
                         "descendant",
                         "role_ancestor",
                         "role_descendant")
    #edges %>% filter(role_ancestor == selected_nodes[9])
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

                            edges_list <- edges %>% filter( ancestor_treatment == selected_nodes[x] |
                                                               ancestor_outcome == selected_nodes[x] |
                                                               ancestor_confounder == selected_nodes[x] |
                                                               ancestor_moc == selected_nodes[x] |
                                                               ancestor_mediator == selected_nodes[x] |
                                                               ancestor_iv == selected_nodes[x]  |
                                                               ancestor_proxy == selected_nodes[x]  |
                                                               ancestor_competing_exposure == selected_nodes[x]  |
                                                               ancestor_latent == selected_nodes[x] )

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
    edges_long <- na.omit(tidyr::pivot_longer(edges, 4:23))

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

#' Convert dagitty object to list
#'
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @noRd
getFeatureMap <- function(dag, include_descendant_nodes = FALSE){
  .datatable.aware <- TRUE

  edges <- getEdges(dag, output_format = "default")
  feature_map_list <- getEdges(dag, output_structure = "list")

  feature_map_new <- c()

  feature_map_new <- lapply( 1:( length(feature_map_list) ), function(role){

    feature_map_list[[role]] <- list(ancestor = unique(feature_map_list[[role]]$ancestor))

    feature_map_list <- lapply(feature_map_list[[role]], function(x) if(identical(x, character(0))) NA_character_ else x)


  })

  names(feature_map_new) <- c("outcome",
                              "treatment",
                              "confounder",
                              "mediator",
                              "mediator_outcome_confounder",
                              "instrumental",
                              "proxy",
                              "competing_exposure",
                              "collider",
                              "latent")

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
