
#' Draw graph edges
#'
#' draw_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param treatments Treatment variable name.
#' @param confounder_vec Vector of variable names, treated as confounders. A list can also be supplied. Order determines the assigned coordinates. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param instrumental_variables Inputted instrumental variable names.
#' @param competing_exposure_vec Vector of competing exposure names. An arrow is drawn connecting competing exposures to the outcome, with other arrows also connected depending on type of graph specified.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @param latent_variables List or vector of additional or already supplied latent (unobserved) variable names.
#' @param collider_vec Vector of collider variables, with both treatment and outcome parents.
#' @param observed Vector of variables without roles, also not unobserved (latent).
#' @param confounder_occurrance Ordered vector of confounders
#' @param existing_dag An existing dagitty object may be supplied.
#' @returns Data frame of edges.
#' @noRd
draw_edges <- function(observed_node_names,
                       type,
                       outcomes,
                       treatments,
                       confounders,
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

  .datatable.aware <- TRUE
  ## outcome edges ##
  outcome_df <- draw_outcome_edges(type, outcomes, collider_vec)
  ## treatment edges ##
  treatment_df <- draw_treatment_edges(type,
                                   outcomes,
                                   treatments,
                                   confounder_vec,
                                   mediator_vec,
                                   collider_vec)
  ## confounder edges ##
  confounder_df <- draw_confounder_edges(type,
                                    outcomes,
                                    treatments,
                                    confounders,
                                    confounder_vec,
                                    m_o_confounder_vec,
                                    mediator_vec,
                                    latent_vec,
                                    confounder_occurrance)
  ## mediator edges ##
  mediator_df <- draw_mediator_edges(type,
                                     outcomes,
                                     mediator_vec,
                                     latent_vec)

  ## mediator_outcome_confounder edges ##
  moc_df <- draw_moc_edges(type,
                           treatments,
                           outcomes,
                           confounder_vec,
                           m_o_confounder_vec,
                           mediator_vec,
                           latent_vec)

  ## competing_exposure edges ##
  competing_exposure_df <- draw_competing_exposure_edges(outcomes, competing_exposure_vec)

  ## connect observed to ancestors and descendants ##
  observed_df <- draw_observed_edges(observed, existing_dag)

  ## instrumental_variables edges ##
  instrumental_df <- draw_iv_edges(type,
                                   instrumental_variables,
                                   treatments = treatments)

  ## latent_variables edges ##
  latent_df <- draw_latent_edges(observed_node_names,
                                 latent_variables,
                                 type,
                                 outcomes,
                                 treatments,
                                 confounder_vec,
                                 m_o_confounder_vec,
                                 mediator_vec)

  ## row bind edges ##
  edges_df <- rbind(treatment_df,
                    outcome_df,
                    confounder_df,
                    moc_df,
                    mediator_df,
                    instrumental_df,
                    latent_df,
                    competing_exposure_df,
                    observed_df,
                    fill=TRUE)  # row bind all edge data frames

  edges_df <- unique(edges_df) # remove duplicate edges

  edges_df <- edges_df[ complete.cases(edges_df), ] # remove NAs

  edges_df <- edges_df[edges_df$ancestor != edges_df$descendant, ] # remove edges with identical ancestor and descendant node names


  return(edges_df)

}



#' Draw outcome edges
#'
#' draw_outcome_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param collider_vec Vector of collider variables, with both treatment and outcome parents.
#' @returns Data frame of edges.
#' @noRd
draw_outcome_edges <- function(type,
                               outcomes,
                               collider_vec
                               ){

  outcome_list <- c()

  # connect colliders
  if( all( complete.cases(collider_vec) ) ){
    ## outcome edges ##

    outcome_list <- suppressWarnings( lapply(1:length(outcomes), function(x){

      outcome_list[x] <- lapply(1:length(collider_vec), function(y){

        list( c( ancestor = outcomes[x], edge = "->", descendant = collider_vec[y]) )

      })

    }) )

    outcome_list <- Filter(Negate(anyNA), unlist(outcome_list, recursive = FALSE) )

    outcome_df <- as.data.table( do.call( rbind, unlist(outcome_list, recursive = FALSE) ) )

  }else{

    outcome_df <- NULL

  }


    # connect all confounders (fully connected or saturated graph type)
  if( length(outcomes) > 1 ){

    if( type == "full" ){

      outcome_list <- suppressWarnings( lapply(1:length(outcomes), function(x){

        outcome_list[x] <- lapply(1:length(outcomes), function(y){

          list( c( ancestor = outcomes[x], edge = "->", descendant = outcomes[y]) )

        })

      }) )

      outcome_list <- Filter(Negate(anyNA), unlist(outcome_list, recursive = FALSE))
      outcome_unlist <- as.data.table( do.call( rbind, unlist(outcome_list, recursive = FALSE) ) )

      outcome_df <- rbind(outcome_df, outcome_unlist)

    }else if( type == "saturated" ){

      outcome_occurrance <- as.numeric(order(match(outcomes, outcomes)))

      outcome_list <- suppressWarnings(lapply(1:length(outcomes), function(x){

        outcome_list[x] <- lapply(1:length(outcomes), function(y){

          list( c( ancestor = outcomes[x], edge = "->", descendant = outcomes[y], ancestor_order = outcome_occurrance[x], descendant_order = outcome_occurrance[y] ) )

        })

      }))

      outcome_list <- Filter(Negate(anyNA), unlist(outcome_list, recursive = FALSE))
      outcome_order_df <- as.data.table( do.call( rbind, unlist(outcome_list, recursive = FALSE) ) )

      outcome_order_df <- outcome_order_df[,c(1:3)][!outcome_order_df$ancestor_order > outcome_order_df$descendant_order, ] # remove rows where temporal logic is not followed
      outcome_df <- rbind(outcome_df, outcome_order_df)

    }


  }

    return(outcome_df)


}


#' Draw treatment edges
#'
#' draw_treatment_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param treatments Treatment variable name.
#' @param confounder_vec Vector of variable names, treated as confounders. A list can also be supplied. Order determines the assigned coordinates. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param collider_vec Vector of collider variables, with both treatment and outcome parents.
#' @returns Data frame of edges.
#' @noRd
draw_treatment_edges <- function(type,
                                 outcomes,
                                 treatments,
                                 confounder_vec,
                                 mediator_vec,
                                 collider_vec
                                 ){
  ## treatment edges ##
  treatment_list <- c()
  # connect outcome
  treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

    treatment_list[x] <- lapply(1:length(outcomes), function(y){

      list( c( ancestor = treatments[x], edge = "->", descendant = outcomes[y]) )

    })

  }) )

  treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
  treatment_df <- as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )


  # connect mediators
  treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

    treatment_list[x] <- lapply(1:length(mediator_vec), function(y){

      list( c( ancestor = treatments[x], edge = "->", descendant = mediator_vec[y]) )

    })

  }) )

  treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
  treatment_unlist <- as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

  treatment_df <- rbind(treatment_df, treatment_unlist)

  # connect colliders
  treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

    treatment_list[x] <- lapply(1:length(collider_vec), function(y){

      list( c( ancestor = treatments[x], edge = "->", descendant = collider_vec[y]) )

    })

  }) )

  treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
  treatment_unlist <- as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

  treatment_df <- rbind(treatment_df, treatment_unlist)

  # connect all confounders (fully connected or saturated graph type)
  if( length(treatments) > 1 ){

    if( type == "full" & length(treatments) > 1 ){

      treatment_list <- suppressWarnings( lapply(1:length(treatments), function(x){

        treatment_list[x] <- lapply(1:length(treatments), function(y){

          list( c( ancestor = treatments[x], edge = "->", descendant = treatments[y]) )

        })

      }) )

      treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
      treatment_unlist <- as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

      treatment_df <- rbind(treatment_df, treatment_unlist)

    }else if( type == "saturated" ){

      treatment_occurrance <- as.numeric(order(match(treatments, treatments)))

      treatment_list <- suppressWarnings(lapply(1:length(treatments), function(x){

        treatment_list[x] <- lapply(1:length(treatments), function(y){

          list( c( ancestor = treatments[x], edge = "->", descendant = treatments[y], ancestor_order = treatment_occurrance[x], descendant_order = treatment_occurrance[y] ) )

        })

      }))

      treatment_list <- Filter(Negate(anyNA), unlist(treatment_list, recursive = FALSE))
      treatment_order_df <- as.data.table( do.call( rbind, unlist(treatment_list, recursive = FALSE) ) )

      treatment_order_df <- treatment_order_df[,c(1:3)][!treatment_order_df$ancestor_order > treatment_order_df$descendant_order, ] # remove rows where temporal logic is not followed
      treatment_df <- rbind(treatment_df, treatment_order_df)

    }

  }


  return(treatment_df)
}


#' Draw confounder edges
#'
#' draw_confounder_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param treatments Treatment variable name.
#' @param confounder_vec Vector of variable names, treated as confounders. A list can also be supplied. Order determines the assigned coordinates. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @param confounder_occurrance Ordered vector of confounders
#' @returns Data frame of edges.
#' @noRd
draw_confounder_edges <- function(type,
                                  outcomes,
                                  treatments,
                                  confounders,
                                  confounder_vec,
                                  m_o_confounder_vec,
                                  mediator_vec,
                                  latent_vec,
                                  confounder_occurrance
                                  ){

  if( all( complete.cases(confounder_vec) ) ){

    ## confounder edges ##
    confounder_list <- c()

    # connect outcome
    confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

      confounder_list[x] <- lapply(1:length(outcomes), function(y){

        list( c( ancestor = confounder_vec[x], edge = "->", descendant = outcomes[y]) )

      })

    }) )

    confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
    confounder_df <- as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )


    # connect treatment
    confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

      confounder_list[x] <- lapply(1:length(treatments), function(y){

        list( c( ancestor = confounder_vec[x], edge = "->", descendant = treatments[y]) )

      })

    }) )

    confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
    confounder_unlist <- as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

    confounder_df <- rbind(confounder_df, confounder_unlist)


    # connect all confounders (fully connected or saturated graph type)
    if( type == "full" ){

      # check if list supplied to build_graph() 'variables'
      if( length( unlist( confounders ) ) > length( confounders ) ){

        confounder_list <- suppressWarnings( lapply(1:(length(confounders) - 1), function(x){

          lapply(1:length(confounders[[x]]), function(y){

            confounder_list[[x]] <- sapply(1:length(confounders[[x+1]]), function(z){

              list( c( ancestor = confounders[[x]][y], edge = "<->", descendant = confounders[[x+1]][z]) )

            })

          })

        }) )

        confounder_unlist <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_unlist <- data.table::as.data.table( do.call( rbind, unlist(confounder_unlist, recursive = FALSE) ) )

        colnames(confounder_df) <- c("v", "e", "w")
        confounder_df <- pdag_to_dag_edges_helper(edges = confounder_df)
        colnames(confounder_df) <- c("ancestor", "edge", "descendant")

        confounder_df <- rbind(confounder_df, confounder_unlist)

      }else{

        confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

          confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

            list( c( ancestor = confounder_vec[x], edge = "->", descendant = confounder_vec[y]) )

          })

        }) )

        confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_unlist <- as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

        confounder_df <- rbind(confounder_df, confounder_unlist)

      }

        confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

          confounder_list[x] <- lapply(1:length(mediator_vec), function(y){

            list( c( ancestor = confounder_vec[x], edge = "->", descendant = mediator_vec[y]) )

          })

        }) )

        confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_unlist <- as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

        confounder_df <- rbind(confounder_df, confounder_unlist)


        # connect mediator-outcome confounders
        confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

          confounder_list[x] <- lapply(1:length(m_o_confounder_vec), function(y){

            list( c( ancestor = confounder_vec[x], edge = "->", descendant = m_o_confounder_vec[y]) )

          })

        }) )

        confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_unlist <- as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

        confounder_df <- rbind(confounder_df, confounder_unlist)


        # connect latent variables
        confounder_list <- suppressWarnings( lapply(1:length(confounder_vec), function(x){

          confounder_list[x] <- lapply(1:length(latent_vec), function(y){

            list( c( ancestor = confounder_vec[x], edge = "->", descendant = latent_vec[y]) )

          })

        }) )

        confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_unlist <- as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

        confounder_df <- rbind(confounder_df, confounder_unlist)


    }

    if( type == "ordered" | type == "saturated" ){

      # check if list supplied to build_graph() 'variables'
      if( length( unlist( confounders ) ) > length( confounders ) ){

        confounder_list <- suppressWarnings( lapply(1:(length(confounders) - 1), function(x){

          lapply(1:length(confounders[[x]]), function(y){

            confounder_list[[x]] <- sapply(1:length(confounders[[x+1]]), function(z){

              list( c( ancestor = confounders[[x]][y], edge = "->", descendant = confounders[[x+1]][z]) )

            })

          })

        }) )

        confounder_unlist <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_unlist <- data.table::as.data.table( do.call( rbind, unlist(confounder_unlist, recursive = FALSE) ) )

        confounder_df <- rbind(confounder_df, confounder_unlist)

      }else{

        confounder_list <- suppressWarnings(lapply(1:length(confounder_vec), function(x){

          confounder_list[x] <- lapply(1:length(confounder_vec), function(y){

            list( c( ancestor = confounder_vec[x], edge = "->", descendant = confounder_vec[y], ancestor_order = confounder_occurrance[x], descendant_order = confounder_occurrance[y] ) )

          })

        }))

        confounder_list <- Filter(Negate(anyNA), unlist(confounder_list, recursive = FALSE))
        confounder_order_df <- as.data.table( do.call( rbind, unlist(confounder_list, recursive = FALSE) ) )

        confounder_order_df <- confounder_order_df[,c(1:3)][!confounder_order_df$ancestor_order > confounder_order_df$descendant_order, ] # remove rows where temporal logic is not followed
        confounder_df <- rbind(confounder_df, confounder_order_df)

      }



    }

    confounder_df <- unique(confounder_df) # remove duplicate edges
    confounder_df <- confounder_df[confounder_df$ancestor != confounder_df$descendant, ] # remove edges with identical ancestor and descendant node names

    return(confounder_df)

  }

  return( data.frame(NULL) )

}

#' Draw mediator edges
#'
#' draw_mediator_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @returns Data frame of edges.
#' @noRd
draw_mediator_edges <- function(type,
                                outcomes,
                                mediator_vec,
                                latent_vec
                                ){

  if( all( complete.cases(mediator_vec) ) ){
    ## mediator edges ##
    mediator_list <- c()

    # connect outcome
    mediator_list <- suppressWarnings( lapply(1:length(mediator_vec), function(x){

      mediator_list[x] <- lapply(1:length(outcomes), function(y){

        list( c( ancestor = mediator_vec[x], edge = "->", descendant = outcomes[y]) )

      })

    }) )

    mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
    mediator_df <- as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )


    # connect all mediators if the inputted dag type is "full"
    if(type == "full"){

      mediator_list <- suppressWarnings( lapply(1:length(mediator_vec), function(x){

        mediator_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = mediator_vec[x], edge = "->", descendant = mediator_vec[y]) )

        })

      }) )

      mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
      mediator_unlist <- as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

      mediator_df <- rbind(mediator_df, mediator_unlist)


      mediator_list <- suppressWarnings( lapply(1:length(mediator_vec), function(x){

        mediator_list[x] <- lapply(1:length(latent_vec), function(y){

          list( c( ancestor = mediator_vec[x], edge = "->", descendant = latent_vec[y]) )

        })

      }) )

      mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
      mediator_unlist <- as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

      mediator_df <- rbind(mediator_df, mediator_unlist)

    }


    if( type == "ordered" ){

      mediator_occurrance <- as.numeric( order( match( mediator_vec, mediator_vec ) ) )

      mediator_list <- suppressWarnings(lapply(1:length(mediator_vec), function(x){

        mediator_list[x] <- lapply(1:length(mediator_vec), function(y){

          list( c( ancestor = mediator_vec[x], edge = "->", descendant = mediator_vec[y], ancestor_order = mediator_occurrance[x], descendant_order = mediator_occurrance[y] ) )

        })

      }))

      mediator_list <- Filter(Negate(anyNA), unlist(mediator_list, recursive = FALSE))
      mediator_order_df <- as.data.table( do.call( rbind, unlist(mediator_list, recursive = FALSE) ) )

      mediator_order_df <- mediator_order_df[,c(1:3)][!mediator_order_df$ancestor_order > mediator_order_df$descendant_order, ] # remove rows where temporal logic is not followed
      mediator_df <- rbind(mediator_df, mediator_order_df)
    }


    mediator_df <- unique(mediator_df) # remove duplicate edges
    mediator_df <- mediator_df[mediator_df$ancestor != mediator_df$descendant, ] # remove edges with identical ancestor and descendant node names

    return(mediator_df)

  }

  return( data.frame(NULL) )

}

#' Draw mediator-outcome confounder edges
#'
#' draw_moc_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param mediator_vec Character or vector of mediator variable names.
#' @param latent_vec Vector of additional or already supplied latent (unobserved) variable names.
#' @returns Data frame of edges.
#' @noRd
draw_moc_edges <- function(type,
                           treatments,
                           outcomes,
                           confounder_vec,
                           m_o_confounder_vec,
                           mediator_vec,
                           latent_vec
                           ){

  if( all( complete.cases(m_o_confounder_vec) ) ){

    ## mediator_outcome_confounder edges ##
    moc_list <- c()

    # connect outcome
    moc_list <- suppressWarnings( lapply(1:length(m_o_confounder_vec), function(x){

      moc_list[x] <- lapply(1:length(m_o_confounder_vec), function(y){

        list( c( ancestor = m_o_confounder_vec[x], edge = "->", descendant = outcomes[y]) )

      })

    }) )

    moc_list <- Filter(Negate(anyNA), unlist(moc_list, recursive = FALSE))
    moc_df <- as.data.table( do.call( rbind, unlist(moc_list, recursive = FALSE) ) )


    # connect treatments if the inputted dag type is "full"
    if(type == "full"){

      moc_list <- suppressWarnings( lapply(1:length(treatments), function(x){

        moc_list[x] <- lapply(1:length(m_o_confounder_vec), function(y){

          list( c( ancestor = treatments[x], edge = "->", descendant = m_o_confounder_vec[y]) )

        })

      }) )

      moc_list <- Filter(Negate(anyNA), unlist(moc_list, recursive = FALSE))
      moc_unlist <- as.data.table( do.call( rbind, unlist(moc_list, recursive = FALSE) ) )

      moc_df <- rbind(moc_df, moc_unlist)

      # connect all latents if the inputted dag type is "full"
      moc_list <- suppressWarnings( lapply(1:length(m_o_confounder_vec), function(x){

        moc_list[x] <- lapply(1:length(latent_vec), function(y){

          list( c( ancestor = m_o_confounder_vec[x], edge = "->", descendant = latent_vec[y]) )

        })

      }) )

      moc_list <- Filter(Negate(anyNA), unlist(moc_list, recursive = FALSE))
      moc_unlist <- as.data.table( do.call( rbind, unlist(moc_list, recursive = FALSE) ) )

      moc_df <- rbind(moc_df, moc_unlist)

    }

    # connect mediators
    moc_list <- suppressWarnings( lapply(1:length(m_o_confounder_vec), function(x){

      moc_list[x] <- lapply(1:length(mediator_vec), function(y){

        list( c( ancestor = m_o_confounder_vec[x], edge = "->", descendant = mediator_vec[y]) )

      })

    }) )

    moc_list <- Filter(Negate(anyNA), unlist(moc_list, recursive = FALSE))
    moc_unlist <- as.data.table( do.call( rbind, unlist(moc_list, recursive = FALSE) ) )

    moc_df <- rbind(moc_df, moc_unlist)

    return(moc_df)

  }

  return( data.frame(NULL) )

}


#' Draw graph edges
#'
#' draw_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param outcomes Outcome variable name.
#' @param competing_exposure_vec Vector of competing exposure names. An arrow is drawn connecting competing exposures to the outcome, with other arrows also connected depending on type of graph specified.
#' @returns Data frame of edges.
#' @noRd
draw_competing_exposure_edges <- function(outcomes,
                                          competing_exposure_vec
                                          ){

  if( all( complete.cases(competing_exposure_vec) ) ){

    ## competing_exposure edges ##
    competing_exposure_list <- c()

    # connect outcome
    competing_exposure_list <- suppressWarnings( lapply(1:length(competing_exposure_vec), function(x){

      competing_exposure_list[x] <- lapply(1:length(competing_exposure_vec), function(y){

        list( c( ancestor = competing_exposure_vec[x], edge = "->", descendant = outcomes[y]) )

      })

    }) )

    competing_exposure_list <- Filter(Negate(anyNA), unlist(competing_exposure_list, recursive = FALSE))
    competing_exposure_df <- as.data.table( do.call( rbind, unlist(competing_exposure_list, recursive = FALSE) ) )

    return(competing_exposure_df)

  }

  return( data.frame(NULL) )

}


#' Draw graph edges
#'
#' draw_edges() is a helper function for buildGraph()
#'
#' @importFrom data.table as.data.table
#' @param observed Vector of variables without roles, also not unobserved (latent).
#' @param existing_dag An existing dagitty object may be supplied.
#' @returns Data frame of edges.
#' @noRd
draw_observed_edges <- function(observed,
                       existing_dag = NA){

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

    observed_list <- Filter(Negate(anyNA), unlist(observed_list, recursive = FALSE))
    observed_df <- as.data.table( do.call( rbind, unlist(observed_list, recursive = FALSE) ) )


    # connect ancestors
    observed_ancestors <- dagitty::parents(existing_dag, observed)

    observed_list <- suppressWarnings( lapply(1:length(observed_ancestors), function(x){

      observed_list[x] <- lapply(1:length(observed), function(y){

        list( c( ancestor = observed_ancestors[x], edge = "->", descendant = observed[y]) )

      })

    }) )

    observed_list <- Filter(Negate(anyNA), unlist(observed_list, recursive = FALSE))
    observed_unlist <- as.data.table( do.call( rbind, unlist(observed_list, recursive = FALSE) ) )

    observed_df <- rbind(observed_df, observed_unlist)

    return(observed_df)

  }

  return( data.frame(NULL) )

}


#' Draw instrumental variable edges
#'
#' draw_iv_edges() is a helper function for buildGraph().
#'
#' @importFrom data.table as.data.table
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param instrumental_variables Inputted list or vector of instrumental variables.
#' @param treatments Vector of treatments.
#' @returns A data frame of instrumental variable edges.
#' @noRd
draw_iv_edges <- function(type, instrumental_variables, treatments){

  if( all( complete.cases( unlist(instrumental_variables) ) ) ){

    instrumental_list <- c()

    if( length(unlist(instrumental_variables)) > length(instrumental_variables) ){

      instrumental_list <- suppressWarnings( lapply( 1:length(instrumental_variables), function(x){

        instrumental_list[x] <- lapply(1:length( instrumental_variables[[x]][[2]] ), function(y){

          list( c( ancestor = instrumental_variables[[x]][[1]], edge = "->", descendant = instrumental_variables[[x]][[2]][[y]] ) )

        } )

      } ) )

      instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
      instrumental_df <- as.data.table( do.call( rbind, unlist(instrumental_list, recursive = FALSE) ) )

    }else{

      instrumental_vec <- as.vector( unlist( lapply( instrumental_variables, function(x) if( identical( x, character(0) ) ) NA_character_ else x ) ) )

      if( any( type == "full" | type == "saturated" ) ){

      # connect treatments
      instrumental_list <- suppressWarnings( lapply(1:length(instrumental_vec), function(x){

        instrumental_list[x] <- lapply(1:length(treatments), function(y){

          list( c( ancestor = instrumental_vec[x], edge = "->", descendant = treatments[y]) )

        })

      }) )

      instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
      instrumental_df <- as.data.table( do.call( rbind, unlist(instrumental_list, recursive = FALSE) ) )

      }else{

        len_trt <- length(treatments)

        # connect treatments
        instrumental_list <- suppressWarnings( lapply(1:length(instrumental_vec), function(x){

          if( x <= len_trt ){

            instrumental_list[x] <- list( c( ancestor = instrumental_vec[x], edge = "->", descendant = treatments[x]) )

            }else{

            instrumental_list[x] <- list( c( ancestor = instrumental_vec[x], edge = "->", descendant = treatments[1]) )

            }

          } ) )

        }

        instrumental_list <- Filter(Negate(anyNA), unlist(instrumental_list, recursive = FALSE))
        instrumental_df <- as.data.table( do.call( rbind, instrumental_list ) )

        }

    return( instrumental_df )

    }

  return( data.frame(NULL) )

}


#' Draw latent variable edges
#'
#' draw_latent_edges() is a helper function for buildGraph().
#'
#' @importFrom data.table as.data.table
#' @param latent_variables Inputted list or vector of latent variables.
#' @param type Type of graph generated. Defaults to 'full' (fully connected graph) with arrows drawn between confounders (both directions) and from confounders to mediators. If type ='saturated', a similar saturated graph is produced except confounders are not connected to mediators, featuring bi-directional arrows between each of the confounders (follows the ESC-DAGs Mapping Stage in Ferguson et al. (2020)). When type = 'ordered', the order of supplied confounders and mediators determines the order that each node occurs, therefore directed arrows are to be connected in one direction from confounders and mediators to other confounders and mediators, respectively. This builds a saturated DAG with temporal, uni-directional arrows, based on Tennnant et al. (2021).
#' @param outcomes Outcome variable name.
#' @param treatments Treatment variable name.
#' @param confounder_vec Vector of variable names, treated as confounders. A list can also be supplied. Order determines the assigned coordinates. If type = "ordered", confounders located in the same list will be assigned similar coordinates.
#' @param m_o_confounder_vec Vector of mediator-outcome confounder names, that instead of being common causes of treatment and outcome (X <- Z -> Y) are a common cause of mediators and outcome (M <- Z -> Y). A list can also be supplied.
#' @param mediator_vec Character or vector of mediator variable names.
#' @returns A data frame of latent variable edges.
#' @noRd
draw_latent_edges <- function(observed_node_names,
                              latent_variables,
                              type,
                              outcomes,
                              treatments,
                              confounder_vec,
                              m_o_confounder_vec,
                              mediator_vec
                              ){


  if( all( complete.cases( unlist(latent_variables) ) ) ){

    ## latent edges ##
    latent_list <- c()

    if( length(unlist(latent_variables)) > length(latent_variables) ){

      latent_list <- suppressWarnings( lapply( 1:length(latent_variables), function(x){

        if( length(unlist(latent_variables[x])) > length(latent_variables[x]) ){

          latent_list[x] <- lapply(1:length( latent_variables[[x]][[2]] ), function(y){

            list( c( ancestor = latent_variables[[x]][[1]], edge = "->", descendant = latent_variables[[x]][[2]][[y]] ) )

          } )

        }

      } ) )

      latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
      latent_df <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

      return(latent_df)

    }else if( length(latent_variables > 1) ){

      # connect outcome
      latent_list <- suppressWarnings( lapply(1:length(latent_variables), function(x){

        latent_list[x] <- lapply(1:length(outcomes), function(y){

          list( c( ancestor = latent_variables[x], edge = "->", descendant = outcomes[y]) )

        })

      }) )

      latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
      latent_df <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )


      # connect treatment
      latent_list <- suppressWarnings( lapply(1:length(latent_variables), function(x){

        latent_list[x] <- lapply(1:length(treatments), function(y){

          list( c( ancestor = latent_variables[x], edge = "->", descendant = treatments[y]) )

        })

      }) )

      latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
      latent_unlist <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

      latent_df <- rbind(latent_df, latent_unlist)


      # connect all latents to confounders, other variables depending on fully connected or saturated graph type
      if( any( type == "full" | type == "saturated" ) ){
        latent_list <- suppressWarnings( lapply(1:length(latent_variables), function(x){

          latent_list[x] <- lapply(1:length(confounder_vec), function(y){

            list( c( ancestor = latent_variables[x], edge = "->", descendant = confounder_vec[y]) )

          })

        }) )

        latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
        latent_unlist <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

        latent_df <- rbind(latent_df, latent_unlist)


        # connect confounders to mediators, confounders to mediator-outcome confounders, and confounders to latents if the inputted dag type is "full"
        if( type == "full" ){

          latent_list <- suppressWarnings( lapply(1:length(latent_variables), function(x){

            latent_list[x] <- lapply(1:length(mediator_vec), function(y){

              list( c( ancestor = latent_variables[x], edge = "->", descendant = mediator_vec[y]) )

            })

          }) )

          latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
          latent_unlist <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

          latent_df <- rbind(latent_df, latent_unlist)


          # connect mediator-outcome confounders
          latent_list <- suppressWarnings( lapply(1:length(latent_variables), function(x){

            latent_list[x] <- lapply(1:length(m_o_confounder_vec), function(y){

              list( c( ancestor = latent_variables[x], edge = "->", descendant = m_o_confounder_vec[y]) )

            })

          }) )

          latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
          latent_unlist <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

          latent_df <- rbind(latent_df, latent_unlist)


          # connect latent variables
          latent_list <- suppressWarnings( lapply(1:length(latent_variables), function(x){

            latent_list[x] <- lapply(1:length(latent_variables), function(y){

              list( c( ancestor = latent_variables[x], edge = "->", descendant = latent_variables[y]) )

            })

          }) )

          latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
          latent_unlist <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

          latent_df <- rbind(latent_df, latent_unlist)

        }

      }else if( type == "ordered" ){

        latent_occurrance <- as.numeric(order(match(latent_variables, latent_variables)))

        latent_list <- suppressWarnings(lapply(1:length(latent_variables), function(x){

          latent_list[x] <- lapply(1:length(latent_variables), function(y){

            list( c( ancestor = latent_variables[x], edge = "->", descendant = latent_variables[y], ancestor_order = latent_occurrance[x], descendant_order = latent_occurrance[y] ) )

          })

        }))

        latent_list <- Filter(Negate(anyNA), unlist(latent_list, recursive = FALSE))
        latent_order_df <- as.data.table( do.call( rbind, unlist(latent_list, recursive = FALSE) ) )

        latent_order_df <- latent_order_df[,c(1:3)][!latent_order_df$ancestor_order > latent_order_df$descendant_order, ] # remove rows where temporal logic is not followed
        latent_df <- rbind(latent_df, latent_order_df)
      }

      latent_df <- unique(latent_df) # remove duplicate edges
      latent_df <- latent_df[latent_df$ancestor != latent_df$descendant, ] # remove edges with identical ancestor and descendant node names

      latent_df <- latent_df[ !unlist(latent_df[,"ancestor"]) %in% observed_node_names ] # remove latent edges containing ancestor nodes with other roles

      return(latent_df)

    }

  }

  return( data.frame(NULL) )

}



#' draws edges for copy_nodes_helper()
#'
#' draw_edges_for_copy_nodes() is called by copy_nodes() through copy_nodes_helper().
#'
#' It connects nodes by calling 'second helper' functions for each case.
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param new_node_names Inputted vector of node names to be added to the graph.
#' @param existing_node_names Inputted vector of node names, used as a reference for the new graph nodes.
#' @param existing_node_type A suffix added to each of the reference node names, e.g. "pre_treatment", or "t0".
#' @param new_node_type A suffix added to each of the new node names, e.g. "post_treatment", or "t" (a number is added for each repeat if num_repeats is specified)
#' @param temporal_reference_node Supplied in the main function, used when an alternative temporal point of reference is desired (default settings use treatment as the temporal point of reference)
#' @param num_repeats Number of additional copies of nodes, such as time points. Each repeat number is included at the end of new node names (new_new_t1, new_node_t2, etc.).
#' @returns A data frame of instrumental variable edges.
#' @noRd
draw_edges_for_copy_nodes_helper <- function(dag, new_node_names, existing_node_names, existing_node_type, new_node_type, temporal_reference_node, num_repeats){

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


#' Fully connect new nodes to others
#'
#' connect_all_nodes_to_new() is a helper function for saturate_nodes() that draws edges between new and existing nodes, in both directions.
#'
#' @importFrom data.table as.data.table
#' @param dag An existing dagitty object.
#' @param new_nodes A vector of new nodes.
#' @noRd
connect_new_nodes <- function(dag, new_nodes, ancestors, descendants, node_roles, type){
  .datatable.aware <- TRUE

  if( length(node_roles) > 0 ){

    output_list <- add_nodes_helper(dag = dag,
                            nodes = new_nodes,
                            node_role = node_role,
                            type = type)
    edges <- output_list$new_edges

    colnames(edges) <- c("v", "e", "w")

    return( edges )

  }

  node_list <- c()
  edges <- c()

  ## get node names
  node_names <- names( dag )

  if( ( length( descendants ) == 0 | length ( descendants[ descendants %in% node_names ] ) == 0  ) & length( ancestors ) == 0 ){

    descendants <- node_names

  }

  if( length( descendants ) > 0){
    ## connect new_nodes to descendants
    node_list <- suppressWarnings( lapply(1:length(new_nodes), function(x){

      node_list[x] <- lapply(1:length(descendants), function(y){

        list( c( v = new_nodes[x], e = "->", w = descendants[y]) )

      })

    }) )

    node_list <- Filter(Negate(anyNA), unlist(node_list, recursive = FALSE))
    edges <- as.data.table( do.call( rbind, unlist(node_list, recursive = FALSE) ) )
  }

  if( ( length( ancestors ) == 0 | length ( ancestors[ ancestors %in% node_names ] ) == 0  ) & length( descendants ) == 0 ){

    ancestors <- node_names

  }

  if( length( ancestors ) > 0){
    ## connect new_nodes to ancestors
    node_list <- suppressWarnings( lapply(1:length(ancestors), function(x){

      node_list[x] <- lapply(1:length(new_nodes), function(y){

        list( c( v = ancestors[x], e = "->", w = new_nodes[y]) )

      })

    }) )

    node_list <- Filter(Negate(anyNA), unlist(node_list, recursive = FALSE))
    node_unlist <- as.data.table( do.call( rbind, unlist(node_list, recursive = FALSE) ) )

    edges <- rbind(edges, node_unlist)
  }

  return( edges )
}
