
#' dagitty node edges and roles
#'
#' get_edges() filters a dagitty object and returns a data frame with edges for specified node roles.
#'
#' @importFrom data.table as.data.table
#' @param dag A dagitty object.
#' @param selected_nodes Nodes to return edges. Defaults to NULL, or can be a character or vector combination of any of the following: c("treatment", "outcome", "confounder", "mediator", "latent", "mediator-outcome-confounder", "instrumental")
#' @param output_structure Outputted data can be a "data.table", "data.frame", or "list".
#' @returns A data frame, data table, or list of edges for the roles specified in selected_nodes.
#' @export
get_edges <- function(dag,
                     selected_nodes = c("outcome",
                                        "treatment",
                                        "confounder",
                                        "mediator",
                                        "mediator_outcome_confounder",
                                        "instrumental", "proxy",
                                        "competing_exposure", "collider",
                                        "latent", "observed"),
                     output_structure = "data.table"){
  .datatable.aware <- TRUE

  all_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent", "observed")

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
  }else if(is.character(selected_nodes) & length(selected_nodes) < length(all_roles)){
    # this ensures the order of variable names in selected_nodes is consistent
    all_nodes <- all_roles

    selected_remove <- NULL
    selected_nodes <- all_nodes[all_nodes %in% selected_nodes]
  }else{

    stop("Please check the selected_nodes input and try again.")

  }


  edges <- extract_unique_node_roles(dag) # add ancestor and descendant nodes (calls a function from later in this file)


  edges <- edges_longer(edges)

  edges <- edges[role_ancestor %in% selected_nodes, ]


  if( output_structure == "data.table" ){


    return(edges)

  }else if( output_structure == "data.frame" ){


    edges <- as.data.frame(edges)

    return(edges)

  }


  edges_list <- lapply( seq_along(selected_nodes), function(x){

    edges_list <- edges[role_ancestor == selected_nodes[x], ]

  } )


  return(edges_list)

  #stop("Invalid 'output_format' - check input and try again.")

}


#' dagitty nodes grouped by role
#'
#' @param dag dagitty object
#' @return Nested list of nodes and node relationships
#' @export
get_roles <- function(dag, multiple_roles = FALSE){
  .datatable.aware <- TRUE

  if(multiple_roles == FALSE){

    edges_wide <- extract_unique_node_roles(dag)

  }else{

    edges_wide <- extract_node_roles(dag)

  }

  node_roles <- c("outcome", "treatment", "confounder", "mediator", "mediator_outcome_confounder", "instrumental", "proxy", "competing_exposure", "collider", "latent", "observed")
  num_roles <- length(node_roles)

  ## ancestor node edges to list ##
  ancestor_roles <- edges_wide[,1:14]
  ancestor_roles <- na.omit( reshape(ancestor_roles, varying = list(4:14), idvar = "id",
                                      v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  ancestor_roles <- ancestor_roles[order(ancestor_roles$id), c(1,4)]
  names(ancestor_roles)[1] <- "node"

  ## group by nodes
  unique_ancestors <- unique( ancestor_roles[,"node"] ) # vector of unique node names
  num_unique_ancestors <- nrow(unique_ancestors) # count of unique node names

  ## descendant node edges to list ##
  descendant_roles <- edges_wide[,c(1:3,15:25)]
  descendant_roles <- na.omit( reshape(descendant_roles, varying = list(4:14), idvar = "id",
                                        v.names = "role", direction = "long")[,c("v", "e", "w", "role", "id")] )
  descendant_roles <- descendant_roles[order(descendant_roles$id), c(3,4)]
  names(descendant_roles)[1] <- "node"

  ## find missing edges
  outcomes <- unique( descendant_roles[ descendant_roles$role == "outcome", ] )
  colliders <- unique( descendant_roles[ descendant_roles$role == "collider", ] )
  observed <- unique( descendant_roles[ descendant_roles$role == "observed", ] )
  latent <- unique( descendant_roles[ descendant_roles$role == "latent", ] )

  missing_outcomes <- outcomes[! unlist(outcomes[,1]) %in% unlist(unique_ancestors), ]
  missing_colliders <- colliders[! unlist(colliders[,1]) %in% unlist(unique_ancestors), ]
  missing_observed <- observed[! unlist(observed[,1]) %in% unlist(unique_ancestors), ]
  missing_latent <- latent[! unlist(latent[,1]) %in% unlist(unique_ancestors), ]

  all_roles <- rbind(ancestor_roles, missing_outcomes, missing_colliders, missing_observed, missing_latent)

  # edges grouped by each unique node in a list
  roles_list <- c()

  roles_list <- suppressWarnings( lapply(1:num_roles, function(x){

    roles_list <-  unique(unlist( all_roles[
      unlist(all_roles[,"role"]) %in% unlist(node_roles[x]) , ][, 1] ))

  }) )

  names(roles_list) <- node_roles # assign node roles as list element names

  roles_list <- lapply(roles_list, function(x) if(identical(x, character(0))) NA else x)


  return(roles_list)

}


#' Assess dagitty object edges
#'
#' assess_edges() provides ways to assess connected edges based on causal criteria and/or user inputs.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges
#' @param dag A dagitty object. Must include exposure and outcome nodes.
#' @param edges_to_keep A vector of directed arrows to be kept between (non-treatment and non-outcome) variables, e.g. c("Z1 -> Z2", "Z2 -> Z3"), c("y", "n", "y"), or c(TRUE, FALSE, TRUE).
#' @param assess_causal_criteria Defaults to FALSE. If TRUE, the user is guided through a sequence that assesses each pair of connected nodes using causal criteria. Based on the Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs) from Ferguson et al. (2020).
#' @returns A list or vector of edges.
#' @export
assess_edges <- function(dag, edges_to_keep = NA, assess_causal_criteria = FALSE){
  .datatable.aware <- TRUE
  #edges_to_keep <- initial_dag
  #dag <- saturated_graph
  edges_to_assess  <- data.table::as.data.table(dagitty::edges(dag))[, c("v", "e", "w")]


  if( dagitty::is.dagitty(edges_to_keep) ){

    edges_to_keep <- data.table::as.data.table(dagitty::edges(edges_to_keep))[, c("v", "e", "w")]

    edges_to_assess <-  edges_to_assess[!edges_to_keep, on = c("v", "e", "w")]

  }else if( all( complete.cases( unlist(edges_to_keep) ) ) ){

    if(is.vector( edges_to_keep ) ){

      edges_to_keep <- data.table::as.data.table( do.call( rbind, strsplit(edges_to_keep, " ") ) )

      colnames(edges_to_keep) <- c("v", "e", "w")

      }

    if( is.data.frame(edges_to_keep) | data.table::is.data.table(edges_to_keep) ){

      edges_to_assess <-  edges_to_assess[!edges_to_keep, on = c("v", "e", "w")]

      }else{

      stop("Edges_to_keep must be a data frame, data table, vector, or dagitty object. Please check inputs and try again.")

      }

    }



  if( assess_causal_criteria == TRUE ){

    check_skip_sequence <- FALSE

   #cat("\nThere are", nrow(edges_to_assess), "directed arrows to be assessed:", "\n", "\n", sep=" ")
   #print(edges_to_assess, quote=FALSE)
   #cat("\nAssess the posited causal relationships? (ESC-DAGs causal criteria and counterfactual thought experiment sequence)", "\n")

    edges_to_keep <- ESC_DAGs_sequence(edges_to_assess, check_skip_sequence, edges_to_keep)

    return(edges_to_keep)

  }

    if( nrow(edges_to_assess) != 0 ){


      if( all( complete.cases( unlist(edges_to_keep) ) ) ){

        ## collapse edges_to_keep to a vector
        num_edges <- nrow(edges_to_keep)

        edges_to_keep <- suppressWarnings( sapply(1:num_edges, function(x){

          edges_to_keep <- paste( edges_to_keep[x,], collapse=" ")

        }) )

      }


      ## collapse edges_to_assess to a vector and output grouped by nodes
      unique_ancestors <- unique( edges_to_assess[,"v"] )
      num_unique_ancestors <- nrow(unique_ancestors)


      edges_assess_list <- list()

      # edges_to_assess grouped by each unique node in a list
      edges_assess_list <- suppressWarnings( lapply(1:num_unique_ancestors, function(x){

        edges_to_assess[ unlist(edges_to_assess[,"v"]) %in% unlist(unique_ancestors[x]), ]


      }) )

      # nodes containing edges are collapsed, outputted in the console to allow easy assessing by copy and paste into .r file
      edges_assess_list <- lapply(1:num_unique_ancestors, function(x){

        edges_assess_list <- noquote(
          paste0( paste0( "'", sapply(1:nrow(edges_assess_list[[x]]), function(y){

            edges_assess_list[x] <- noquote( paste( edges_assess_list[[x]][y,], collapse=" "  ) )

            }), "'", collapse=", " ), sep = "") )

      })

      cat( paste("c(", paste( unlist(edges_assess_list), collapse=",\n\n" ), ")", sep = "\n", collapse = "") )

    }else{

      stop("There are no edges to assess. Please check the supplied dagitty object and try again.")


    }



  if( !all( complete.cases(edges_to_keep) ) & assess_causal_criteria == FALSE ){

    edges_to_assess <- suppressWarnings( sapply(1:nrow(edges_to_assess), function(x){

      edges_to_assess <- paste( edges_to_assess[x,], collapse=" ")

    }) )

    message("\nOutputted edges to assess.
            \n\nPrinted edges to assess - copy and paste in a .R file to use as a vector object.")

    return(edges_to_assess)

  }

  edges_list <- list(
    edges_to_keep = edges_to_keep,
    edges_to_assess = edges_to_assess
  )

  message("\nOutputted edges list (edges_to_keep & edges_to_assess).
          \n\nPrinted edges to assess - copy and paste in a .R file to use as a vector object.")

  return(edges_list)

}


#' Remove dagitty object edges
#'
#' keep_edges() removes edges based on user inputs.
#'
#' @importFrom data.table as.data.table is.data.table
#' @importFrom dagitty edges exposures outcomes latents coordinates dagitty isAcyclic
#' @param dag A saturated graph dagitty object. Exposure and outcome must be indicated, and optionally can include assigned coordinates.
#' @param edges_to_keep A vector of directed arrows to be kept between (non-treatment and non-outcome) variables, e.g. c("Z1 -> Z2", "Z2 -> Z3"), c("y", "n", "y"), or c(TRUE, FALSE, TRUE).
#' @returns A dagitty object, with directed arrows removed based on edges_to_keep.
#' @export
keep_edges <- function(dag, edges_to_keep = NA){
  .datatable.aware <- TRUE

  edges_all  <- as.data.table(dagitty::edges(dag))
  edges <- edges_all[, c("v", "e", "w")]


  if( dagitty::is.dagitty(edges_to_keep)){

    edges_to_keep <- data.table::as.data.table(dagitty::edges(edges_to_keep))[, c("v", "e", "w")]

  }else if( is.vector( edges_to_keep ) ){

    edges_to_keep <- data.table::as.data.table( do.call( rbind, strsplit(edges_to_keep, " ") ) )

    colnames(edges_to_keep) <- c("v", "e", "w")

  }

  if( all(is.null(edges_to_keep)) | all(is.na(edges_to_keep)) ){

    stop("No edges to keep were supplied.")

  }


  if( all( !is.na(edges_to_keep) )){

    if( ( is.data.frame(edges_to_keep) | data.table::is.data.table(edges_to_keep) ) ){

      edges <-  edges[ edges_to_keep, on = c("v", "e", "w")]


    }else{

      stop("Edges_to_keep must be a data frame, data table, vector, or dagitty object. Please check inputs and try again.")

    }
  }

  dag <- rebuild_dag(dag, edges)


  if(dagitty::isAcyclic(dag) == FALSE){
    warning("The outputted grpah contains cycles, and is therefore not a directed acyclic graph (DAG). Relationships may need to be further assessed.")
  }

  return(dag)

}


#' ESC-DAGs causal criteria for removinng edges
#' @param edges vector of edges whose relationships are to be assessed
#' @param num_edges number of edges to be assessed
#' @param check_skip_sequence TRUE or FALSE depending on prior inputs
#' @noRd
ESC_DAGs_sequence <- function(edges, check_skip_sequence, edges_to_keep){
  #edges <- edges_to_assess # used for debugging assess_edges()
  #check_skip_sequence <- FALSE # used for debugging assess_edges()
  #edges_to_keep <- NA # used for debugging assess_edges()
  if(check_skip_sequence == FALSE) {

    num_edges <- nrow(edges)

    cat("\nThere are ", num_edges, " directed arrows to be assessed.","\n", "\n", sep="")
    print(edges, quote=FALSE)
    cat("\nAssess the posited causal relationships using causal criteria? (ESC-DAGs protocol)", "\n")

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
        #cat("\nCode written by Aidan J Moller (2025)")
        #cat("\n===================================", "\n")

        check_ans <- TRUE

      }else if(choice == "n"){

        message("\nSkipped sequence.")

        #edges <- noquote(paste("c('", paste(edges, collapse="', '"), "')", sep = ""))

        # edges_to_keep <- noquote(paste("c('", paste(edges_to_keep, collapse="', '"), "')", sep = ""))

        message("\nOutputting causal criteria assessed edges.")

        edges_list <- list(
          edges_to_keep = edges_to_keep,
          edges_to_assess = edges
        )

        return(edges_list)

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

      arrow <- noquote( paste(edges[arrow_count], collapse=" ") )
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

          message("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

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

            #num_arrow_to_remove <- num_arrow_to_remove + 1
            #removed_arrows[num_arrow_to_remove] <- arrow_count

            #check_arrow_removed <- TRUE

            check_ans <- TRUE

            #cat("\n", "\nCausal relationship '", arrow, "' assessed; edge removed.", "\n", sep="")

            cat("\n", "\nCausal relationship '", arrow, "' assessed; answer recorded.", "\n", sep="")

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

    edges <- edges[-removed_arrows, ]

    edges_to_keep <- rbind(edges_to_keep, edges)

    return(edges_to_keep)

  }else{

    message("\nNo arrows were removed.", "\n")

    edges_list <- list(
      edges_to_keep = edges_to_keep,
      edges_to_assess = edges
    )

    return(edges_list)
  }
}
