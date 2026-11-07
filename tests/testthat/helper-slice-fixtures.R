slice_external_query <- function(..., prohibit_left = FALSE) { # Exercise external dotted identities without the bundled slice metadata
  ids <- c("trial.management.application", "toy.soil.quantity.left",
           "toy.soil.quantity.right", "harvest.output.yield")
  nodes <- as.data.frame(setNames(rep(list(rep("", length(ids))), ncol(apsim_dag_nodes)),
                                 names(apsim_dag_nodes)), stringsAsFactors = FALSE)
  nodes$node_id <- ids
  nodes$concept_id <- paste0("c:", seq_along(ids))
  nodes$node_name_readable <- c("Application", "Left soil quantity", "Right soil quantity", "Yield")
  nodes$function_of <- c("", "c:1", "c:1", "c:2;c:3")
  nodes$crop_type <- "barley"
  nodes$hidden_node <- "FALSE"
  nodes$path_dependent <- "FALSE"
  nodes$recovered_node <- "yes"
  if (prohibit_left) {
    nodes$recovered_node[2L] <- "no:Excluded fixture quantity"
  }
  nodes$tier <- c("management", "", "", "")
  aliases <- setNames(lapply(seq_along(ids), function(index) { # Give each toy quantity its own unambiguous names
    return(c(ids[index], nodes$node_name_readable[index]))
  }), ids)
  tables <- list(apsim_dag_nodes = nodes, node_aliases = aliases,
                 management_defaults = list(core = ids[1L], per_crop = list()),
                 crop_phases = crop_phases[crop_phases$crop == "barley", ],
                 phase_of = phase_of[FALSE, ])
  query <- query_graph
  environment(query) <- list2env(list(getExportedValue = function(package, name) { # Read only this fixture's tables
    if (identical(package, "causaliflower") && name %in% names(tables)) {
      return(tables[[name]])
    }
    return(base::getExportedValue(package, name))
  }), parent = environment(query_graph))
  return(query("barley", treatments = ids[1L], outcomes = ids[4L], quiet = TRUE, ...))
}
