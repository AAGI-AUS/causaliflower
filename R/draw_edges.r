
#' Draw graph edges
#'
#' draw_edges() is a helper function for buildGraph()
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @param dag A dagitty object.
#' @returns Data frame of edges.
#' @noRd
draw_edges <- function(type,
                       outcomes,
                       treatments,
                       confounder_vec,
                       m_o_confounder_vec,
                       mediator_vec,
                       instrumental_variables,
                       competing_exposure_vec,
                       latent_vec,
                       latent_variables,
                       collider_vec,
                       observed,
                       confounder_occurrance,
                       existing_dag = NA){

  ## outcome edges ##
  outcome_list <- c()
  # connect colliders
  if( all( complete.cases(collider_vec) ) ){

    outcome_list <- suppressWarnings( lapply(1:length(outcomes), function(x){

      outcome_list[x] <- lapply(1:length(collider_vec), function(y){

        list( c( ancestor = outcomes[x], edge = "->", descendant = collider_vec[y]) )

      })

    }) )

    outcome_list <- Filter( Negate(anyNA), unlist(unlist(outcome_list, recursive = FALSE), recursive = FALSE) )


  }
  outcome_df <- dplyr::bind_rows(outcome_list)

  ## treatment edges ##
  treatment_list <- c()
  # connect outcome
  treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

    treatment_list[x] <- lapply(1:length(outcomes), function(y){

      list( c( ancestor = treatments[x], edge = "->", descendant = outcomes[y]) )

    })

  }) )

  treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
  treatment_df <- dplyr::bind_rows(treatment_list)

  # connect mediators
  treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

    treatment_list[x] <- lapply(1:length(mediator_vec), function(y){

      list( c( ancestor = treatments[x], edge = "->", descendant = mediator_vec[y]) )

    })

  }) )

  treatment_list <- Filter( Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE) )
  treatment_df <- dplyr::bind_rows(treatment_df, treatment_list)

  # connect colliders
  treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

    treatment_list[x] <- lapply(1:length(collider_vec), function(y){

      list( c( ancestor = treatments[x], edge = "->", descendant = collider_vec[y]) )

    })

  }) )

  treatment_list <- Filter( Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE) )
  treatment_df <- dplyr::bind_rows(treatment_df, treatment_list)

  # connect all confounders (fully connected or saturated graph type)
  if( type == "full" ){

    treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

      treatment_list[x] <- lapply(1:length(treatments), function(y){

        list( c( ancestor = treatments[x], edge = "->", descendant = treatments[y]) )

      })

    }) )

    treatment_list <- Filter(Negate(anyNA), unlist(unlist(treatment_list, recursive = FALSE), recursive = FALSE))
    treatment_df <- dplyr::bind_rows(treatment_df, treatment_list)
  }

  ## confounder edges ##
  confounder_list <- c()

  # connect outcome
  confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

    confounder_list[x] <- lapply(1:length(outcomes), function(y){

      list( c( ancestor = confounder_vec[x], edge = "->", descendant = outcomes[y]) )

    })

  }) )


  confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
  confounder_df <- dplyr::bind_rows(confounder_list)

  # connect treatment
  confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

    confounder_list[x] <- lapply(1:length(treatments), function(y){

      list( c( ancestor = confounder_vec[x], edge = "->", descendant = treatments[y]) )

    })

  }) )

  confounder_list <- Filter(Negate(anyNA),  unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE))
  confounder_df <- dplyr::bind_rows(confounder_df, confounder_list)

  # connect all confounders (fully connected or saturated graph type)
  if( type == "full" | type == "saturated" ){

    confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

      confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

        list( c( ancestor = confounder_vec[x], edge = "->", descendant = confounder_vec[y]) )

      })

    }) )

    confounder_list <- Filter(Negate(anyNA), unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE))
    confounder_df <- dplyr::bind_rows(confounder_df, confounder_list)

    # connect confounders to mediators, and confounders to mediator-outcome confoundres if the inputted dag type is "full"
    if( type == "full" ){

      confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

        confounder_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = confounder_vec[x], edge = "->", descendant = mediator_vec[y]) )

        })

      }) )

      confounder_list <- Filter(Negate(anyNA), unlist(unlist(confounder_list, recursive = FALSE), recursive = FALSE))
      confounder_df <- dplyr::bind_rows(confounder_df, confounder_list)


      # connect mediator-outcome confounders
      confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

        confounder_list[x] <- lapply(1:length(m_o_confounder_vec), function(y){

          list( c( ancestor = confounder_vec[x], edge = "->", descendant = m_o_confounder_vec[y]) )

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

  if( all( complete.cases(mediator_vec) ) ){
    ## mediator edges ##
    mediator_list <- c()

    # connect outcome
    mediator_list <- suppressWarnings( lapply(1:length(mediator_vec), function(x){

      mediator_list[x] <- lapply(1:length(outcomes), function(y){

        list( c( ancestor = mediator_vec[x], edge = "->", descendant = outcomes[y]) )

      })

    }) )

    mediator_list <- Filter( Negate(anyNA), unlist(mediator_list, recursive = FALSE) )
    mediator_df <- dplyr::bind_rows(mediator_list)

    # connect all mediators if the inputted dag type is "full"
    if(type == "full"){

      mediator_list <- suppressWarnings( lapply(1:length(mediator_vec), function(x){

        mediator_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = mediator_vec[x], edge = "->", descendant = mediator_vec[y]) )

        })

      }) )

      mediator_list <- Filter( Negate(anyNA), unlist(unlist(mediator_list, recursive = FALSE), recursive = FALSE) )
      mediator_df <- dplyr::bind_rows(mediator_df, mediator_list)

      mediator_list <- suppressWarnings( lapply(1:length(mediator_vec), function(x){

        mediator_list[x] <- lapply(1:length(latent_vec), function(y){

          list( c( ancestor = mediator_vec[x], edge = "->", descendant = latent_vec[y]) )

        })

      }) )

      mediator_list <- Filter( Negate(anyNA), unlist(unlist(mediator_list, recursive = FALSE), recursive = FALSE) )
      mediator_df <- dplyr::bind_rows(mediator_df, mediator_list)

    }

  }else{

    mediator_df <- data.frame(NULL)

  }


  if( all( complete.cases(m_o_confounder_vec) ) ){

    ## mediator_outcome_confounder edges ##
    moc_list <- c()

    # connect outcome
    moc_list <- suppressWarnings( lapply(1:length(m_o_confounder_vec), function(x){

      moc_list[x] <- lapply(1:length(m_o_confounder_vec), function(y){

        list( c( ancestor = m_o_confounder_vec[x], edge = "->", descendant = outcomes[y]) )

      })

    }) )

    moc_list <- Filter(Negate(anyNA), unlist(unlist(moc_list, recursive = FALSE), recursive = FALSE))
    moc_df <- dplyr::bind_rows(moc_list)

    # connect treatments if the inputted dag type is "full"
    if(type == "full"){

      moc_list <- suppressWarnings( lapply(1:length(treatments), function(x){

        moc_list[x] <- lapply(1:length(m_o_confounder_vec), function(y){

          list( c( ancestor = treatments[x], edge = "->", descendant = m_o_confounder_vec[y]) )

        })

      }) )

      moc_list <- Filter(Negate(anyNA), unlist( unlist(moc_list, recursive = FALSE), recursive = FALSE))
      moc_df <- dplyr::bind_rows(moc_df, moc_list)

      # connect all latents if the inputted dag type is "full"
      moc_list <- suppressWarnings( lapply(1:length(m_o_confounder_vec), function(x){

        moc_list[x] <- lapply(1:length(latent_vec), function(y){

          list( c( ancestor = m_o_confounder_vec[x], edge = "->", descendant = latent_vec[y]) )

        })

      }) )

      moc_list <- Filter(Negate(anyNA), unlist( unlist(moc_list, recursive = FALSE), recursive = FALSE))
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

  }else{

    moc_df <- data.frame(NULL)

  }


  if( all( complete.cases(competing_exposure_vec) ) ){

    ## competing_exposure edges ##
    competing_exposure_list <- c()

    # connect outcome
    competing_exposure_list <- suppressWarnings( lapply(1:length(competing_exposure_vec), function(x){

      competing_exposure_list[x] <- lapply(1:length(competing_exposure_vec), function(y){

        list( c( ancestor = competing_exposure_vec[x], edge = "->", descendant = outcomes[y]) )

      })

    }) )

    competing_exposure_list <- Filter(Negate(anyNA), unlist(unlist(competing_exposure_list, recursive = FALSE), recursive = FALSE))
    competing_exposure_df <- dplyr::bind_rows(competing_exposure_list)

  }else{

    competing_exposure_df <- data.frame(NULL)

  }

  if( all( complete.cases(observed) ) ) {

    ## connect observed to ancestors and descendants ##
    observed_list <- list()

    # connect descendants
    observed_descendants <- dagitty::children(existing_dag, observed)

    observed_list <- suppressWarnings( lapply(1:length(observed), function(x){

      observed_list[x] <- lapply(1:length(observed_descendants), function(y){

        list( c( ancestor = observed[x], edge = "->", descendant = observed_descendants[y]) )

      })

    }) )

    observed_list <- Filter(Negate(anyNA), unlist(unlist(observed_list, recursive = FALSE), recursive = FALSE))
    observed_df <- dplyr::bind_rows(observed_list)

    # connect ancestors
    observed_ancestors <- dagitty::parents(existing_dag, observed)

    observed_list <- suppressWarnings( lapply(1:length(observed_ancestors), function(x){

      observed_list[x] <- lapply(1:length(observed), function(y){

        list( c( ancestor = observed_ancestors[x], edge = "->", descendant = observed[y]) )

      })

    }) )

    observed_list <- Filter(Negate(anyNA), unlist(unlist(observed_list, recursive = FALSE), recursive = FALSE))
    observed_df <- dplyr::bind_rows(observed_df, observed_list)


  }else{

    observed_df <- data.frame(NULL)

  }


  ## instrumental_variables edges ##
  instrumental_df <- draw_iv_edges(instrumental_variables, treatments)


  ## latent_variables edges ##
  latent_df <- draw_latent_edges(latent_variables)



  edges_df <- rbind(treatment_df, outcome_df, confounder_df, moc_df, mediator_df, instrumental_df, latent_df, competing_exposure_df, observed_df)  # row bind all edge data frames

  edges_df <- unique(edges_df) # remove duplicate edges

  edges_df <- na.omit(edges_df) # remove NAs

  edges_df <- edges_df[edges_df$ancestor != edges_df$descendant, ] # remove edges with identical ancestor and descendant node names


  return(edges_df)

}


#' Draw instrumental variable edges
#'
#' draw_iv_edges() is a helper function for buildGraph().
#'
#' @param instrumental_variables Inputted list or vector of instrumental variables.
#' @param treatments Vector of treatments.
#' @returns A data frame of instrumental variable edges.
#' @noRd
draw_iv_edges <- function(instrumental_variables, treatments){


  if( all( complete.cases( unlist(instrumental_variables) ) ) ){

    instrumental_list <- c()

    if( length(unlist(instrumental_variables)) > length(instrumental_variables) ){

      instrumental_list <- suppressWarnings( lapply( 1:length(instrumental_variables), function(x){

        instrumental_list[x] <- lapply(1:length( instrumental_variables[[x]][[2]] ), function(y){

          list( c( ancestor = instrumental_variables[[x]][[1]], edge = "->", descendant = instrumental_variables[[x]][[2]][[y]] ) )

        } )

      } ) )


      instrumental_list <- Filter(Negate(anyNA), unlist(unlist(instrumental_list, recursive = FALSE), recursive = FALSE))
      instrumental_df <- dplyr::bind_rows(instrumental_list)


    }else{


      instrumental_vec <- as.vector( unlist( lapply( instrumental_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

      ## instrumental variable edges ##
      instrumental_list <- c()

      # connect outcome
      instrumental_list <- suppressWarnings( lapply(1:length(instrumental_vec), function(x){

        instrumental_list[x] <- lapply(1:length(instrumental_vec), function(y){

          list( c( ancestor = instrumental_vec[x], edge = "->", descendant = treatments[y]) )

        })

      }) )

      instrumental_list <- Filter(Negate(anyNA), unlist(unlist(instrumental_list, recursive = FALSE), recursive = FALSE))
      instrumental_df <- dplyr::bind_rows(instrumental_list)


    }


  }else{

    instrumental_df <- NA

  }


  return(instrumental_df)

}


#' Draw latent variable edges
#'
#' draw_latent_edges() is a helper function for buildGraph().
#'
#' @param latent_variables Inputted list or vector of latent variables.
#' @returns A data frame of latent variable edges.
#' @noRd
draw_latent_edges <- function(latent_variables){

  if( all( complete.cases( unlist(latent_variables) ) ) ){

    latent_list <- c()

    if( length(unlist(latent_variables)) > length(latent_variables) ){

      latent_list <- suppressWarnings( lapply( 1:length(latent_variables), function(x){

        if( length(unlist(latent_variables[x])) > length(latent_variables[x]) ){

          latent_list[x] <- lapply(1:length( latent_variables[[x]][[2]] ), function(y){

            list( c( ancestor = latent_variables[[x]][[1]], edge = "->", descendant = latent_variables[[x]][[2]][[y]] ) )

          } )

        }

      } ) )


      latent_list <- Filter(Negate(anyNA), unlist(unlist(latent_list, recursive = FALSE), recursive = FALSE))
      latent_df <- dplyr::bind_rows(latent_list)


    }else{

      latent_df <- NA
    }


  }else{

    latent_df <- NA

  }


  return(latent_df)

}


#' Draw new node edges
#'
#' draw_new_node_edges() is a helper function for addNodes().
#'
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_node_names Inputted vector of node names to be added to the graph.
#' @param existing_node_names Inputted vector of node names, used as a reference for the new graph nodes.
#' @param existing_node_type A suffix added to each of the reference node names, e.g. "pre_treatment", or "t0".
#' @param new_node_type A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param temporal_reference_node Supplied in the main function, used when an alternative temporal point of reference is desired (default settings use treatment as the temporal point of reference)
#' @param num_repeats Number of additional copies of nodes, such as time points. Each repeat number is included at the end of new node names (new_new_t1, new_node_t2, etc.).
#' @returns A data frame of instrumental variable edges.
#' @noRd
draw_new_node_edges <- function(dag, new_node_names, existing_node_names, existing_node_type, new_node_type, temporal_reference_node, num_repeats){

  ## connect reference & new nodes ##

  # connect new nodes to reference nodes & consecutive new nodes (same node, at a different point in time)
  new_nodes_df <- connect_new_and_existing_nodes(new_node_names, existing_node_names)

  treatments <- dagitty::exposures(dag)


  ## fetch parent and children groups for connecting nodes ##

  # get each reference node's descendants/children
  existing_node_children <- lapply(1:length(existing_node_names), function(x){

    existing_node_children <- dagitty::children(dag, existing_node_names[x])

  })

  # get each reference node's ancestors/parents
  existing_node_parents <- lapply(1:length(existing_node_names), function(x){

    existing_node_parents <- dagitty::parents(dag, existing_node_names[x])

  })

  # get parents of new nodes, which also happen to be the existing reference nodes
  new_node_parents_in_existing_nodes <- lapply(1:length(existing_node_names), function(x){

    new_node_parents_in_existing_nodes <- existing_node_parents[[x]][ existing_node_parents[[x]] %in% existing_node_names ]

  })


  ## connect parent nodes ##

  # a list of edges connecting parent existing reference nodes to their descendants, at each time point
  new_nodes_list <- connect_new_nodes_to_parent_new_nodes(new_node_names, existing_node_names, new_node_parents_in_existing_nodes) # needs checking single new_node_name input

  if( !all( is.na(new_nodes_list) ) ){

    new_nodes_df <- rbind( new_nodes_df, new_nodes_list )

  }

  # we can remove the existing reference nodes from parents (to be connected later in this function)
  new_node_parents <- lapply(1:length(existing_node_names), function(x){

    new_node_parents <- existing_node_parents[[x]][ !existing_node_parents[[x]] %in% existing_node_names ]

  })

  treatment_parents <- dagitty::parents(dag, treatments)


  ## connecting children nodes ##
  if( existing_node_type == "pre_treatment" & !all( complete.cases(temporal_reference_node) ) ){  # if existing type = "pre_treatment", all new nodes are treated as post-treatment nodes


    new_node_children <- lapply(1:length(existing_node_names), function(x){

      new_node_children <- existing_node_children[[x]][ !existing_node_children[[x]] %in% treatments &
                                                 !existing_node_children[[x]] %in% existing_node_names[[x]] &
                                                 !existing_node_children[[x]] %in% treatment_parents ]

    })




  }else if( all( complete.cases(temporal_reference_node) ) ){

    node_names <- names(dag)

    if( length( node_names[ node_names %in% temporal_reference_node ] ) !=0 ){

      temporal_reference_node_parents <- dagitty::parents(dag, temporal_reference_node)

      new_node_children <- lapply(1:length(existing_node_names), function(x){

        new_node_children <- existing_node_children[[x]][ !existing_node_children[[x]] %in% temporal_reference_node &
                                                            !existing_node_children[[x]] %in% existing_node_names[[x]] &
                                                            !existing_node_children[[x]] %in% temporal_reference_node_parents ]

      })


    }


  }else{ # otherwise, we will just connect all new nodes to treatment and treatment parents (irrespective of whether this creates a dag)


    new_node_children <- lapply(1:length(existing_node_names), function(x){

      new_node_children <- existing_node_children[[x]][ !existing_node_children[[x]] %in% existing_node_names[[x]] ]

    } )



  }



  new_nodes_list <- connect_post_treatment_node_parents(new_node_names, existing_node_names, new_node_parents)

  if( !all( is.na(na.omit(new_nodes_list) )) ){

    new_nodes_df <- rbind( new_nodes_df, new_nodes_list )

  }


  new_nodes_list <- connect_post_treatment_node_children(new_node_names, existing_node_names, new_node_children)

  if( !all( is.na(new_nodes_list) ) ){

    new_nodes_df <- rbind( new_nodes_df, new_nodes_list )

  }







  return(new_nodes_df)

}


