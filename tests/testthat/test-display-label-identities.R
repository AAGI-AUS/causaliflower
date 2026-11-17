label_identity_fixture <- function(vertices){ # Build a declared chain with a canonical query ledger
  edges <- data.table::data.table(from = head(vertices, -1L), to = tail(vertices, -1L))
  graph <- dagitty::dagitty(causaliflower:::.emit_dagitty(
    vertices, edges, vertices[2L], vertices[1L], tail(vertices, 1L)))
  attr(graph, "recovered_vertices") <- vertices
  attr(graph, "edges") <- edges
  attr(graph, "latent") <- vertices[2L]
  attr(graph, "treatments") <- vertices[1L]
  attr(graph, "outcomes") <- tail(vertices, 1L)
  ledger <- data.frame(node = vertices, quantity = causaliflower:::.strip_slice(vertices),
    member_id = causaliflower:::.strip_slice(vertices), occurrence = vertices,
    member_occurrence = vertices)
  attr(graph, "query_record") <- list(version = 1L, node_map = ledger,
    source_node_map = ledger, requests = list(treatments = vertices[1L]),
    source_policy = data.frame(node_id = ledger$quantity, recovered_node = "yes"))
  attr(graph, "full_vertices") <- vertices
  attr(graph, "full_edges") <- data.table::copy(edges)
  return(graph)
}


expect_label_identity <- function(graph, shown, expected){ # Compare graph identities independently of label selection
  original <- attr(graph, "recovered_vertices")
  mapping <- stats::setNames(expected, original)
  expect_identical(attr(shown, "recovered_vertices"), expected)
  expect_setequal(names(shown), expected)
  expect_equal(anyDuplicated(expected), 0L)
  for(role in c("latent", "treatments", "outcomes")){
    expect_identical(attr(shown, role), unname(mapping[attr(graph, role)]))
  }
  expect_identical(dagitty::latents(shown), unname(mapping[attr(graph, "latent")]))
  expect_identical(dagitty::exposures(shown), unname(mapping[attr(graph, "treatments")]))
  expect_identical(dagitty::outcomes(shown), unname(mapping[attr(graph, "outcomes")]))
  inverse <- stats::setNames(original, expected)
  edges <- dagitty::edges(shown)
  edges$v <- unname(inverse[edges$v])
  edges$w <- unname(inverse[edges$w])
  edge_keys <- function(table){ return(sort(paste(table$v, table$e, table$w))) }
  expect_identical(edge_keys(edges), edge_keys(dagitty::edges(graph)))
  expected_edges <- data.table::copy(attr(graph, "edges"))
  expected_edges$from <- unname(mapping[expected_edges$from])
  expected_edges$to <- unname(mapping[expected_edges$to])
  expect_identical(attr(shown, "edges"), expected_edges)
  record <- attr(graph, "query_record")
  record$node_map$node <- unname(mapping[record$node_map$node])
  record$display_transform <- list(list(operation = "readable_dag", mapping = mapping))
  expect_identical(attr(shown, "query_record"), record)
  expect_identical(attr(shown, "full_vertices"), attr(graph, "full_vertices"))
  expect_identical(attr(shown, "full_edges"), attr(graph, "full_edges"))
  return(invisible(NULL))
}


test_that("readable labels resolve long fallback chains without changing graph identities", {
  graph <- label_identity_fixture(LETTERS[1:7])
  nodes <- data.frame(node_id = LETTERS[1:7],
    node_name_readable = c("X", "X", "A", "C", "D", "E", "Display G"))
  original <- serialize(graph, NULL)
  original_nodes <- serialize(nodes, NULL)
  set.seed(237)
  rng <- .Random.seed
  shown <- readable_dag(graph, nodes)
  expect_label_identity(graph, shown, c(LETTERS[1:6], "Display G"))
  expect_identical(serialize(graph, NULL), original)
  expect_identical(serialize(nodes, NULL), original_nodes)
  expect_identical(.Random.seed, rng)
})


test_that("readable labels resolve clashes introduced by temporal suffixes", {
  graph <- label_identity_fixture(c("A_t1", "B", "C", "D_t0"))
  nodes <- data.frame(node_id = LETTERS[1:4],
    node_name_readable = c("Plant", "Plant_t1", "A_t1", "Display D"))
  original <- serialize(graph, NULL)
  set.seed(239)
  rng <- .Random.seed
  shown <- readable_dag(graph, nodes)
  expect_label_identity(graph, shown, c("A_t1", "B", "C", "Display D_t0"))
  expect_identical(serialize(graph, NULL), original)
  expect_identical(.Random.seed, rng)
})


test_that("input labels preserve roles and arrows through both collision types", {
  cases <- list(
    list(vertices = LETTERS[1:7], bases = LETTERS[1:7],
      labels = c("X", "X", "A", "C", "D", "E", "Display G"),
      expected = c(LETTERS[1:6], "Chosen G"), override = c(G = "Chosen G")),
    list(vertices = c("A_t1", "B", "C", "D_t0"), bases = LETTERS[1:4],
      labels = c("Plant", "Plant_t1", "A_t1", "Display D"),
      expected = c("A_t1", "B", "C", "Display D_t0"), override = c(A = "Plant")))
  for(case in cases){
    graph <- label_identity_fixture(case$vertices)
    nodes <- data.frame(node_id = case$bases, node_name_readable = case$labels)
    edges <- attr(graph, "edges")
    original <- serialize(list(graph, nodes, case), NULL)
    set.seed(241)
    rng <- .Random.seed
    shown <- causaliflower:::.apply_input_labels(case$vertices, edges,
      attr(graph, "latent"), attr(graph, "treatments"), attr(graph, "outcomes"),
      nodes, case$override, "readable")
    expect_identical(shown$V, case$expected)
    mapping <- stats::setNames(case$expected, case$vertices)
    expect_identical(shown$E, data.table::data.table(
      from = unname(mapping[edges$from]), to = unname(mapping[edges$to])))
    expect_identical(shown$lat, unname(mapping[attr(graph, "latent")]))
    expect_identical(shown$X, unname(mapping[attr(graph, "treatments")]))
    expect_identical(shown$Y, unname(mapping[attr(graph, "outcomes")]))
    expect_identical(shown$input_labels,
      if( "G" %in% case$bases ){ c("Chosen G" = "G") }else{ character() })
    expect_identical(serialize(list(graph, nodes, case), NULL), original)
    expect_identical(.Random.seed, rng)
  }
})


test_that("query label fallback preserves canonical receipts and source structures", {
  ids <- c("G", LETTERS[1:6], "H", "I")
  nodes <- as.data.frame(stats::setNames(rep(list(rep("", length(ids))),
    ncol(apsim_dag_nodes)), names(apsim_dag_nodes)), stringsAsFactors = FALSE)
  nodes$node_id <- ids
  nodes$concept_id <- paste0("c:", seq_along(ids))
  nodes$node_name_readable <- c("Treatment", "X", "X", "A", "C", "D", "E",
                               "Outcome", "Visible I")
  nodes$function_of <- c("", "G", LETTERS[1:5], "F", "G")
  nodes$crop_type <- "barley"
  nodes$hidden_node <- "FALSE"
  nodes$path_dependent <- "FALSE"
  nodes$recovered_node <- "yes"
  nodes$tier[1L] <- "management"
  tables <- list(apsim_dag_nodes = nodes,
    node_aliases = stats::setNames(as.list(ids), ids),
    management_defaults = list(core = "G", per_crop = list()),
    crop_phases = crop_phases[crop_phases$crop == "barley", ], phase_of = phase_of[FALSE, ])
  query <- query_graph
  environment(query) <- list2env(list(getExportedValue = function(package, name){
    if( identical(package, "causaliflower") && name %in% names(tables) ){
      return(tables[[name]])
    }
    return(base::getExportedValue(package, name))
  }), parent = environment(query_graph))
  original_tables <- serialize(tables, NULL)
  source <- query("barley", treatments = "G", outcomes = "H", time_unroll = "none", quiet = TRUE)
  original_source <- serialize(source, NULL)
  set.seed(251)
  rng <- .Random.seed
  shown <- query("barley", treatments = "G", outcomes = "H", time_unroll = "none",
                 label = "readable", quiet = TRUE)
  mapping <- stats::setNames(ifelse(ids == "I", "Visible I", ids), ids)
  expect_setequal(names(shown), unname(mapping))
  expect_identical(dagitty::exposures(shown), "G")
  expect_identical(dagitty::outcomes(shown), "H")
  expected_edges <- data.table::copy(attr(source, "edges"))
  expected_edges$from <- unname(mapping[expected_edges$from])
  expected_edges$to <- unname(mapping[expected_edges$to])
  expect_identical(attr(shown, "edges"), expected_edges)
  source_record <- attr(source, "query_record")
  shown_record <- attr(shown, "query_record")
  expected_map <- source_record$node_map
  expected_map$node <- unname(mapping[expected_map$node])
  expect_identical(shown_record$node_map, expected_map)
  for(field in c("requests", "source_node_map", "recovery_folding", "policy",
                "source_policy", "policy_scope", "policy_applied", "aliases", "scope")){
    expect_identical(shown_record[[field]], source_record[[field]])
  }
  for(field in c("full_vertices", "full_edges", "full_latent", "contemporaneous_edges",
                "lagged_edges", "persistent_stocks", "unroll_depth")){
    expect_identical(attr(shown, field), attr(source, field))
  }
  expect_identical(serialize(source, NULL), original_source)
  expect_identical(serialize(tables, NULL), original_tables)
  expect_identical(.Random.seed, rng)
})


test_that("management options retain their column types when no choices apply", {
  columns <- c("node_id", "label", "aliases", "status", "crop_type")
  nodes <- data.frame(node_id = "A", node_name_readable = "Alpha",
    tier = "state", crop_type = "barley")
  for(crop in list(NULL, "wheat")){
    options <- get_management_options(crop, nodes, list(A = "A"), list(core = character()))
    expect_s3_class(options, "data.table")
    expect_identical(names(options), columns)
    expect_equal(nrow(options), 0L)
    expect_identical(vapply(options, typeof, character(1L)),
                     stats::setNames(rep("character", 5L), columns))
  }
  nodes$tier <- "management"
  options <- get_management_options("wheat", nodes, list(A = "A"), list(core = "A"))
  expect_identical(names(options), columns)
  expect_equal(nrow(options), 0L)
  expect_identical(vapply(options, typeof, character(1L)),
                   stats::setNames(rep("character", 5L), columns))
  expect_identical(get_management_options("barley", nodes, list(A = "A"), list(core = "A")),
    data.table::data.table(node_id = "A", label = "Alpha", aliases = "A",
                          status = "core", crop_type = "barley"))
})
