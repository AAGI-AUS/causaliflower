
#' Draw instrumental variable edges
#'
#' draw_edges_instrumental_variables() is a helper function for buildGraph().
#'
#' @param instrumental_variables Inputted list or vector of instrumental variables.
#' @returns A vector of latent variable names.
#' @noRd
get_latent_vec <- function(latent_variables){

  if( length(unlist(latent_variables)) > length(latent_variables) ){
    latent_variables_list <- list()
    latent_variables_list <- lapply( 1:length(latent_variables), function(x){

      latent_variables_list[[x]] <- as.data.frame( latent_variables[[x]])[1]

    } )


    latent_variables_df <- dplyr::bind_cols( latent_variables_list )
    latent_vec <- as.vector( unlist(latent_variables_df) )

  }else{

    latent_vec <- as.vector( unlist( lapply( latent_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )


  }

  return(latent_vec)

}

#' extract_instrumental_variables() is a helper to instrumental_variables()
#'
#' @importFrom magrittr %>%
#' @param dag A dagitty object.
#' @returns Vector of instrumental variable names.
#' @noRd
extract_instrumental_variables <- function(dag, treatments, outcomes, latent_vars, colliders, treatment_parents, treatment_children, mediator_parents, latent_children, outcome_children, outcome_parents){

  nodes_trt_to_y <- nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  nodes_trt_to_y_parents <- dagitty::parents(dag, nodes_trt_to_y)

  # instrumental variables
  latent_parents <- dagitty::parents(dag, latent_vars)

  collider_children <- dagitty::children(dag, colliders)

  exclude_as_instruments <- treatment_parents[

    treatment_parents %in% treatment_children | # not allowed as an instrument
      treatment_parents %in% mediator_parents | # not allowed as an instrument
      treatment_parents %in% latent_parents | # not allowed as an instrument
      treatment_parents %in% latent_children | # not allowed as an instrument
      treatment_parents %in% outcome_children | # cyclic relationship
      treatment_parents %in% outcome_parents | # cyclic relationship
      treatment_parents %in% collider_children | # cyclic relationship
      treatment_parents %in% latent_vars | # exclude latent variables
      treatment_parents %in% nodes_trt_to_y_parents # in path treatment -> ... -> outcome

  ]

  instrumental_vars <- as.data.frame( treatment_parents[!treatment_parents %in% exclude_as_instruments] )
  # instrumental_vars <- rbind(instrumental_vars, instrumental_vars)
  if(length(unlist(instrumental_vars)) > 0){

    instrumental_vars$role <- "instrumental"

  }

  if(length(unname(as.vector(instrumental_vars[,1]))) > 0){


    #instrumental_children <- dagitty::children(dag, unname(as.vector(instrumental_vars[,1]))[1])

    # length(instrumental_children[instrumental_children %in% treatments]) > 1


    instrumental_children_multiple_treatments <- lapply( 1:length(unname(as.vector(instrumental_vars[,1]))), function(x){

      instrumental_children <- list( node = list( as.character(unname(instrumental_vars[,1]))[x], dagitty::children(dag, unname(as.vector(instrumental_vars[,1]))[x]) ) )

    })


    instrumental_children_multiple_treatments <- lapply( 1:length(unname(as.vector(instrumental_vars[,1]))), function(x){

      instrumental_children_multiple_treatments[[x]][["node"]][[1]][ length( instrumental_children_multiple_treatments[[x]][["node"]][[2]][ instrumental_children_multiple_treatments[[x]][["node"]][[2]]  %in% treatments ] ) > 1 ]

    })

    instrumental_children_multiple_treatments <- unlist(instrumental_children_multiple_treatments)

  }else{

    instrumental_children_multiple_treatments <- NULL

  }


  if( length( instrumental_children_multiple_treatments ) > 0  ){

    #warning("Potential instrumental variable(s) connected to multiple treatments. These will be ignored, if intended to be used as instruemnts please check edges and try again.
    #     \nIf inputting multiple instruments/treatments to buildGraph(), supply a nested list to 'instrumental_variables' indicating the descendant treatment,
    #     e.g., list( list('instrument 1', list('treatment 1')), list('instrument 2', list('treatment 2')) ) ")


    instrumental_vars$role[ instrumental_vars[,1] %in% instrumental_children_multiple_treatments] <- "observed"

  }

  return(instrumental_vars)

}

#' Extract node roles
#'
#' extract_node_roles() is a helper function for getEdges().
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @returns Data table of edges containing roles for each ancestor and descendant nodes.
#' @noRd
extract_node_roles <- function(dag){

  edges_dagitty <- data.table::as.data.table(dagitty::edges(dag))[,1:3]

  edges <- edges_dagitty[, c("v", "e", "w")]

  # treatment
  treatments <- dagitty::exposures(dag)

  # outcome
  outcomes <- dagitty::outcomes(dag)

  # latent variables
  latent_vars <- dagitty::latents(dag)

  # confounders
  treatment_parents <- dagitty::parents(dag, treatments)
  confounders <- treatment_parents[treatment_parents %in% dagitty::parents(dag, outcomes) &
                                     !treatment_parents %in% treatments &
                                     !treatment_parents %in% latent_vars]

  # mediators - first parse (includes mediator-outcome confounders)
  outcome_parents <- dagitty::parents(dag, outcomes)
  treatment_children <- dagitty::children(dag, treatments)

  mediators <- outcome_parents[outcome_parents %in% treatment_children &
                                 !outcome_parents %in% treatments &
                                 !outcome_parents %in% confounders &
                                 !outcome_parents %in% latent_vars]

  # mediator-outcome confounders
  mediator_parents <- dagitty::parents(dag, mediators) # filter to include only parents of mediator variables
  moc <- mediator_parents[mediator_parents %in% outcome_parents] # include only nodes connected to both mediators and outcome (M <- MOC -> Y)
  moc <- moc[ !moc %in% c(treatments, confounders) & # remove treatment and confounder nodes
                !moc %in% treatment_parents & # double check by removing parents of treatment
                !moc %in% latent_vars ]

  nodes_trt_to_y <- nodes_between_treatment_and_outcome(dag, treatments, outcomes)

  # competing exposure
  competing_exposure <- outcome_parents[ !outcome_parents %in% mediators &
                                           !outcome_parents %in% treatments &
                                           !outcome_parents %in% confounders &
                                           !outcome_parents %in% latent_vars &
                                           !outcome_parents %in% moc &
                                           !outcome_parents %in% nodes_trt_to_y ]

  # proxy

  latent_children <- dagitty::children(dag, latent_vars)
  proxy_b <- treatment_parents[ treatment_parents %in% latent_children &
                                  !treatment_parents %in% latent_vars &
                                  !treatment_parents %in% treatments] # proxy_b
  proxy_c <- outcome_parents[ outcome_parents %in% latent_children & # proxy_c
                                !outcome_parents %in% treatments &
                                !outcome_parents %in% mediators &
                                !outcome_parents %in% latent_vars ]
  proxy <- c(proxy_b, proxy_c)

  # collider
  outcome_children <- dagitty::children(dag, outcomes)
  colliders <- outcome_children[ outcome_children %in% treatment_children &
                                   !outcome_children %in% latent_vars ]

  # instrumental variables
  instrumental_vars <- extract_instrumental_variables(dag,
                                                      treatments,
                                                      outcomes,
                                                      latent_vars,
                                                      colliders,
                                                      treatment_parents,
                                                      treatment_children,
                                                      mediator_parents,
                                                      latent_children,
                                                      outcome_children,
                                                      outcome_parents)


  instrumental_vars <- instrumental_vars[,1][instrumental_vars$role ==  "instrumental"]


  # filter latent variables
  latent_vars <- latent_vars[  !latent_vars %in% treatments &
                                 !latent_vars %in% outcomes ]

  # catch any remaining variables
  all_vars <- names(dag)

  observed <- all_vars[!all_vars %in% treatments &
                         !all_vars %in% colliders &
                         !all_vars %in% proxy &
                         !all_vars %in% competing_exposure &
                         !all_vars %in% mediators &
                         !all_vars %in% moc &
                         !all_vars %in% confounders &
                         !all_vars %in% outcomes &
                         !all_vars %in% instrumental_vars &
                         !all_vars %in% latent_vars]


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

  edges$ancestor_latent[edges$v %in% latent_vars] <- "latent"

  edges$ancestor_observed[edges$v %in% observed] <- "observed"


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

  edges$descendant_latent[edges$w %in% latent_vars] <- "latent"

  edges$descendant_observed[edges$w %in% observed] <- "observed"

  return(edges)

}

#' Edges to longer format
#'
#' @importFrom data.table as.data.table
#' @param edges Data table of edges, created using extract_node_roles().
#' @returns Data table of edges in longer format.
#' @noRd
edges_longer <- function(edges){

  edges$id <- 1:nrow(edges)

  edges_ancestors <- edges[,1:14]
  edges_descendants <- edges[,c(1:3,15:25)]

  edges_ancestors <- na.omit( reshape(edges_ancestors, varying = list(4:14), idvar = "id",
                                      v.names = "role_ancestor", direction = "long")[,c("v", "e", "w", "role_ancestor", "id")] )
  edges_ancestors <- edges_ancestors[order(edges_ancestors$id), 1:4]

  edges_descendants <- na.omit( reshape(edges_descendants, varying = list(4:14), idvar = "id",
                                        v.names = "role_descendant", direction = "long")[,c("v", "e", "w", "role_descendant", "id")] )
  edges_descendants <- edges_descendants[order(edges_descendants$id), 1:4]


  if(nrow(edges_ancestors) != nrow(edges_descendants)){ # finds missing latent edges

    edges_descendants <- find_missing_edges(edges_ancestors, edges_descendants)

  }


  edges <- do.call(cbind, lapply(list(edges_ancestors, edges_descendants[,4]), data.table::setDT))

  names(edges)[1:3] <- c("ancestor", "edge", "descendant")

  return(edges)

}

#' Find missing latent edges
#'
#' find_missing_edges() is a helper function for getEdges().
#'
#' @importFrom data.table as.data.table
#' @param dplyr anti_join semi_join
#' @param edges_ancestors Data table of edges containing ancestor node roles.
#' @param edges_descendants Data table of edges containing descendant node roles.
#' @returns A data table containing ancestor edges.
#' @noRd
find_missing_edges <- function(edges_ancestors, edges_descendants){
  .datatable.aware <- TRUE

  # this was written when I had an issue with latent variables being ignored, unsure if it has been fixed.
  # ideally this would be replaced with R core functions, or removed if no longer necessary

  missing_ancestors <- dplyr::anti_join(edges_descendants, edges_ancestors, by = c("v" = "v"))
  common_ancestors <- dplyr::semi_join(edges_descendants, edges_ancestors, by = c("v" = "v"))

  combined_ancestors <- rbind(edges_ancestors[,1:3], common_ancestors[,1:3])
  combined_ancestors[! duplicated(combined_ancestors, fromLast=TRUE) &
                       seq(nrow(combined_ancestors)) <= nrow(edges_ancestors), ]

  common_ancestors[,1:3][ common_ancestors[,1:3] %in% edges_ancestors[,1:3] ]
  missing_latent_edges <- dplyr::anti_join(edges_ancestors, common_ancestors, by = c("v", "e", "w"))

  dplyr::anti_join(edges_ancestors, edges_descendants, by = c("w" = "w"))

  missing_latent_rows <- NULL

  missing_latent_rows <- lapply( 1:nrow(missing_ancestors), function(x){
    missing_row <- missing_ancestors[x, ]
    missing_row[,"value"][missing_row[, 1] %in% latent_vars] <- "latent"
    missing_row <- as.data.frame(missing_row)

  })

  edges_ancestors <- na.omit( dplyr::bind_rows(edges_ancestors, missing_latent_rows) )

  return(edges_ancestors)

}

#' node names in path from treatment to outcome
#'
#' @importFrom dplyr bind_rows
#' @importFrom ggdag dag_paths
#' @param dag A dagitty object
#' @param treatments A vector of treatment node names
#' @param outcomes A vector of outcome node names
#' @returns Vector of nodes in path from treatment to outcome
#' @noRd
nodes_between_treatment_and_outcome <- function(dag, treatments, outcomes){

  paths_trt_to_y <- lapply(1:length(treatments), function(x){

    paths_trt_to_y <- lapply(1:length(outcomes), function(y){

      paths_trt_to_y <- na.omit(ggdag::dag_paths(dag,
                                                 from = treatments[x],
                                                 to = outcomes[y],
                                                 directed = TRUE,
                                                 paths_only = TRUE)[["data"]])

      #paths_trt_to_y <- paths_trt_to_y[!grepl(treatments[x], paths_trt_to_y$to),]
    })

  })

  paths_trt_to_y <- dplyr::bind_rows(paths_trt_to_y)

  nodes <- as.vector(unique(paths_trt_to_y$name))

  nodes <- unique( c(nodes, as.vector(na.omit(unique(paths_trt_to_y$to)))) )

  nodes <- nodes[ !nodes %in% treatments]

  return(nodes)

}


#' Add nodes to a dagitty object
#'
#' create_new_node_names() is a helper function for addNodes(). It uses existing nodes to create new node names for specified repeats/time points.
#'
#' @importFrom magrittr %>%
#' @param reference_nodes Vector of existing node names, used as a reference for the new graph nodes, e.g., c("Z1", "Z2", "Z3").
#' @param new_node_type A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param num_repeats Number of additional copies of nodes, such as time points. Each repeat number is included at the end of new node names (new_new_t1, new_node_t2, etc.).
#' @returns A vector of new node names.
#' @export
create_new_node_names <- function(reference_nodes, new_node_type, num_repeats){

  seq_repeats <- seq(num_repeats)

  new_node_names <- NULL

  if( length( reference_nodes ) > 0 ){

    if( length( seq_repeats ) == 1 ){

      new_node_names <- paste0(reference_nodes, "_", new_node_type)


    }else{

      new_node_names <- sapply(1:length(seq_repeats), function(x){

        new_node_names <- paste0(reference_nodes, "_", new_node_type, seq_repeats[x])

      })

    }

  }

  return(new_node_names)


}


#' Connect new nodes to their reference nodes
#'
#' connect_new_and_reference_nodes() is a helper function for draw_new_node_edges() and addNodes().
#'
#' @param new_node_names Inputted vector of node names to be added to the graph.
#' @param reference_node_names Inputted vector of node names, used as a reference for the new graph nodes.
#' @returns A data frame of edges connecting new nodes to their reference nodes, and each subsequent time point to the next.
#' @noRd
connect_new_and_reference_nodes <- function(new_node_names, reference_node_names){

  new_node_names_vec <- as.character(unlist(new_node_names, use.names = FALSE))

  num_ref_nodes <- length(reference_node_names)

  new_nodes_list <- suppressWarnings( lapply(1:num_ref_nodes, function(x){

    new_nodes_list <- list( v = reference_node_names[x], e = "->", w = new_node_names_vec[x] )

  }) )

  new_nodes_df <- as.data.frame( do.call( rbind, new_nodes_list ) )

  num_new_nodes <- length(new_node_names_vec)

  if( num_new_nodes > num_ref_nodes & length(new_node_names) > 1 ){

    new_nodes_list <- suppressWarnings( lapply(1:( nrow(new_node_names) - 1 ), function(x){

      new_nodes_list <- lapply( 1:nrow(new_node_names), function(y){

        new_nodes_list <- list( v = new_node_names[y, x], e = "->", w = new_node_names[y, (x+1)] )

      })

    }) )

    if( !all( is.null(unlist(new_nodes_list)) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, unlist(new_nodes_list, recursive = FALSE) ) )
    }

    new_nodes_df <- rbind(new_nodes_df, new_nodes_list)

  }else{

    new_nodes_list <- suppressWarnings( lapply( 1:( nrow(new_node_names) - 1 ), function(y){

        new_nodes_list <- list( v = new_node_names[y,], e = "->", w = new_node_names[y + 1,] )

        })
      )


    if( !all( is.null(unlist(new_nodes_list)) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, new_nodes_list) )

      new_nodes_list[] <- lapply(new_nodes_list, as.character)
    }

    new_nodes_df <- rbind(new_nodes_df, new_nodes_list)

  }




  return(new_nodes_df)

}


#' Connect new nodes to parents
#'
#' connect_new_nodes_to_parent_new_nodes() is a helper function for draw_new_node_edges() and addNodes().
#'
#' @param new_node_names Inputted vector of node names to be added to the graph.
#' @param new_node_parents_in_reference_nodes Inputted vector of new node parent names in the supplied reference nodes.
#' @returns A list of edges connecting new nodes to their parent nodes, at time point.
#' @noRd
connect_new_nodes_to_parent_new_nodes <- function(new_node_names, reference_node_names, new_node_parents_in_reference_nodes){ # needs checking single new_node_name input

  # check if more than one reference node
  if( length(reference_node_names) > 1 ){

    # connect new node to parent nodes (also included in reference nodes)
    new_nodes_list <- suppressWarnings(

      lapply(1:( length(new_node_names) ), function(t){ # t = each time point

        new_nodes_list <- lapply( 1:nrow(new_node_names), function(y){ # y = each reference (new) node e.g. Z1_a, Z2_a, Z3_a

          if( length( unlist(new_node_parents_in_reference_nodes[[y]]) ) > 0){

            new_nodes_list <- lapply( 1:length( unlist(new_node_parents_in_reference_nodes[[y]]) ), function(x){ # x = each parent node

              if( length( unlist(new_node_parents_in_reference_nodes[[y]][[x]]) ) > 0 ){

                if( length( unlist(new_node_parents_in_reference_nodes[[y]][[x]]) ) > 1){

                  new_nodes_list <- list( v =  new_node_parents_in_reference_nodes[[y]][[x,t]], e = "->", w = new_node_names[y, t] )

                }else{

                  new_nodes_list <- list( v =  new_node_parents_in_reference_nodes[[y]][[x]], e = "->", w = new_node_names[y, t] )

                }
              }

            })

          }
        })

      })
    )


    new_nodes_list <- do.call( rbind, new_nodes_list )

    if( !all( is.null(unlist(new_nodes_list)) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, unlist(new_nodes_list, recursive = FALSE) ) )
    }

    new_nodes_list[] <- lapply(new_nodes_list, as.character)


  }else{

    new_nodes_list <- suppressWarnings(

      lapply(1:( nrow(new_node_names) ), function(t){ # t = each time point

        if( length( unlist(new_node_parents_in_reference_nodes) ) > 0 ){

            new_nodes_list <- lapply( 1:length( unlist(new_node_parents_in_reference_nodes[[t]]) ), function(x){ # x = each parent node

              if( length( unlist(new_node_parents_in_reference_nodes[[t]][[x]]) ) > 0 ){

                  new_nodes_list <- list( v =  new_node_parents_in_reference_nodes[[t]][[x]], e = "->", w = new_node_names[t,] )

              }

            } )

        }

      } )
    )

    if( !all( is.null(unlist(new_nodes_list)) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, new_nodes_list) )

      new_nodes_list[] <- lapply(new_nodes_list, as.character)
    }

  }

  return(new_nodes_list)
}


#' Connect new nodes to parents
#'
#' connect_post_treatment_node_parents() is a helper function for draw_new_node_edges() and addNodes().
#'
#' @param new_node_names Inputted vector of node names to be added to the graph.
#' @param new_node_parents Inputted vector of new node parent names.
#' @returns A list of edges connecting new nodes to their parent nodes, at time point.
#' @noRd
connect_post_treatment_node_parents <- function(new_node_names, reference_node_names, new_node_parents){

  if( length(reference_node_names) > 1 ){
    # new nodes parents
    new_nodes_list <- suppressWarnings(

      lapply(1:( length(new_node_names) ), function(t){ # t = each time point

        new_nodes_list <- lapply( 1:nrow(new_node_names), function(y){ # y = each reference (new) node e.g. Z1_a, Z2_a, Z3_a

          if( length( unlist(new_node_parents[[y]]) ) > 0){

            new_nodes_list <- lapply( 1:length( unlist(new_node_parents[[y]]) ), function(x){ # x = each parent node


              if( length( unlist(new_node_parents[[y]][[x]]) ) > 0 ){

                new_nodes_list <- list( v =  new_node_parents[[y]][[x]], e = "->", w = new_node_names[y, t] )

              }else{

                new_nodes_list <- list( v =  new_node_parents[[y]], e = "->", w = new_node_names[y, t] )


              }


            })

          }


        })

      }) )

    new_nodes_list <- do.call( rbind, new_nodes_list )

    if( !all( is.null(unlist(new_nodes_list)) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, unlist(new_nodes_list, recursive = FALSE) ) )
    }

    new_nodes_list[] <- lapply(new_nodes_list, as.character)

  }else{

    # new nodes parents
    new_nodes_list <- suppressWarnings(

      lapply(1:( nrow(new_node_names) ), function(t){ # t = each time point

          if( length( unlist(new_node_parents) ) > 0 ){

            new_nodes_list <- lapply( 1:length( new_node_parents ), function(y){ # x = each parent node

              if( length( unlist(new_node_parents[[y]]) ) > 0 ){

                new_nodes_list <- lapply( 1:length( unlist(new_node_parents[[y]]) ), function(x){ # x = each parent node

                  if( length( unlist(new_node_parents[[y]][[x]]) ) > 0 ){

                    new_nodes_list <- list( v =  new_node_parents[[y]][[x]], e = "->", w = new_node_names[t,] )

                  }else{

                    new_nodes_list <- list( v =  new_node_parents[[y]], e = "->", w = new_node_names[t,] )


                  }

                })


              }

            })

          }

      }) )

    new_nodes_list <- do.call( rbind, new_nodes_list )

    #if( !all( is.null( unlist(new_nodes_list) ) ) & length( unlist(new_node_parents) ) > 1 ){
    #  new_nodes_list <-  as.data.table( do.call( rbind, new_nodes_list) )
    #  new_nodes_list[] <- lapply(new_nodes_list, as.character)
    #}else

    if( !all( is.null(unlist(new_nodes_list) ) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, unlist(new_nodes_list, recursive = FALSE) ) )

      new_nodes_list[] <- lapply(new_nodes_list, as.character)


    }


  }

  return(new_nodes_list)
}


#' Connect new nodes to children
#'
#' connect_post_treatment_node_children() is a helper function for draw_new_node_edges() and addNodes().
#'
#' @param new_node_names Inputted vector of node names to be added to the graph.
#' @param new_node_children Inputted vector of new nodes' children.
#' @returns A list of edges connecting new nodes to their child nodes, at time point.
#' @noRd
connect_post_treatment_node_children <- function(new_node_names, reference_node_names, new_node_children){


  if( length(reference_node_names) > 1 ){

    # new nodes children
    new_nodes_list <- suppressWarnings(

      lapply(1:( length(new_node_names) ), function(t){ # t = each time point

        new_nodes_list <- lapply( 1:nrow(new_node_names), function(y){ # y = each reference (new) node e.g. Z1_a, Z2_a, Z3_a

          if( length( unlist(new_node_children[[y]]) ) > 0){

            new_nodes_list <- lapply( 1:length( unlist(new_node_children[[y]]) ), function(x){ # x = each children node


              if( length( unlist(new_node_children[[y]][[x]]) ) > 0 ){

                new_nodes_list <- list( v = new_node_names[y, t], e = "->", w = new_node_children[[y]][[x]] )

              }else{

                new_nodes_list <- list( v =  new_node_names[y, t], e = "->", w =  new_node_children[[y]])


              }

            })

          }


        })

      }) )

    new_nodes_list <- do.call( rbind, new_nodes_list )

    if( !all( is.null(unlist(new_nodes_list)) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, unlist(new_nodes_list, recursive = FALSE) ) )
    }

    new_nodes_list[] <- lapply(new_nodes_list, as.character)

  }else{

    # new nodes children
    new_nodes_list <- suppressWarnings(

      lapply(1:( nrow(new_node_names) ), function(t){ # t = each time point

        if( length( unlist(new_node_children) ) > 0 ){

          new_nodes_list <- lapply( 1:length( new_node_children ), function(y){ # y = each child node

            if( length( unlist(new_node_children[[y]]) ) > 0 ){

              new_nodes_list <- lapply( 1:length( unlist(new_node_children[[y]]) ), function(x){ # x = each child node

                if( length( unlist(new_node_children[[y]][[x]]) ) > 0 ){

                  new_nodes_list <- list( v = new_node_names[t,], e = "->", w = new_node_children[[y]][[x]])

                }else{

                  new_nodes_list <- list( v = new_node_names[t,], e = "->", w = new_node_children[[y]])


                }

              })


            }

          })

        }

      }) )

    new_nodes_list <- do.call( rbind, new_nodes_list )

    #if( !all( is.null( unlist(new_nodes_list) ) ) & length( unlist(new_node_children) ) > 1 ){
    #  new_nodes_list <-  as.data.table( unlist( do.call( rbind, new_nodes_list), recursive = FALSE ) )
    #  new_nodes_list[] <- lapply(new_nodes_list, as.character)

   #}else
    if( !all( is.null(unlist(new_nodes_list) ) ) ){

      new_nodes_list <-  as.data.table( do.call( rbind, unlist(new_nodes_list, recursive = FALSE) ) )

      new_nodes_list[] <- lapply(new_nodes_list, as.character)

    }


  }

  return(new_nodes_list)
}

