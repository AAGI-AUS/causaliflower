test_that("directed folding bypasses chains and forks without latent projection", {
  edges <- data.table::data.table(
    from = c("A", "H1", "H1", "H2", "A", "Fork", "Fork"),
    to = c("H1", "H2", "C", "B", "B", "B", "C"))
  vertices <- unique(c(edges$from, edges$to))
  original <- serialize(edges, NULL)
  f <- .query_fold_nodes(vertices, edges, c("H1", "H2", "Fork"))
  expect_setequal(f$vertices, c("A", "B", "C"))
  expect_setequal(paste(f$edges$from, f$edges$to), c("A B", "A C"))
  expect_identical(paste(f$added_edges$from, f$added_edges$to), "A C")
  expect_equal(nrow(f$edges), 2L)
  expect_identical(serialize(edges, NULL), original)
  expect_identical(f$source_edges, edges)
  # Persistence through removed occasions is a directed path too.
  temporal <- data.table::data.table(from = c("A_t0", "H_t0", "H_t1"),
                                      to = c("H_t0", "H_t1", "B_t1"))
  f <- .query_fold_nodes(unique(c(temporal$from, temporal$to)), temporal, c("H_t0", "H_t1"))
  expect_identical(paste(f$edges$from, f$edges$to), "A_t0 B_t1")
})

fold_query <- function(columns = NULL, show = TRUE, randomised = FALSE, ...) {
  query_graph("barley", treatments = "N_application", outcomes = c("grain_N", "grain_wt"),
    columns = columns, time_unroll = "selective", randomised_treatment = randomised,
    label = NULL, quiet = TRUE, params = query_params(show_hidden_nodes = show, ...))
}

test_that("query folds interpretation rows and preserves the full explanation", {
  g <- fold_query()
  f <- attr(g, "recovery_folding")
  r <- attr(g, "query_record")
  expect_equal(length(unique(.strip_slice(f$policy_removed))), 8L)
  expect_false(any(f$removed_nodes %in% names(g)))
  expect_true(all(f$removed_nodes %in% r$source_node_map$node))
  expect_false(any(f$removed_nodes %in% r$node_map$node))
  expect_true(all(f$removed_nodes %in% f$source_vertices))
  expect_true(all(.strip_slice(f$policy_removed) %in% r$policy$node_id))
  expect_gt(nrow(f$added_edges), 0L)
  expect_false(any(dagitty::edges(g)$e != "->"))
  expect_identical(attr(g, "rolled_structure"), "source_before_folding")
  expect_identical(attr(g, "edges"), attr(g, "full_edges"))
  # Reconstruct the expected adjacency with an independent elimination procedure.
  e <- as.data.frame(f$source_edges)
  for (node in f$removed_nodes) {
    parents <- e$from[e$to == node & e$from != node]
    children <- e$to[e$from == node & e$to != node]
    added <- expand.grid(from = parents, to = children, stringsAsFactors = FALSE)
    e <- unique(rbind(e[e$from != node & e$to != node, ], added))
    e <- e[e$from != e$to, ]
  }
  expect_setequal(paste(e$from, e$to), paste(f$edges$from, f$edges$to))
  expect_true(.is_acyclic(names(g), attr(g, "edges")))
  surgery <- fold_query(randomised = TRUE)
  expect_false(any(attr(surgery, "edges")$to %in% attr(surgery, "treatments")))
})

test_that("explicit interpretations survive hidden filtering but prohibited rows do not", {
  g <- fold_query(columns = c("cgr_actual", "n_translocation"), show = FALSE)
  expect_true(all(c("cgr_actual", "n_translocation") %in% .strip_slice(names(g))))
  expect_false(any(.strip_slice(attr(g, "latent")) %in% c("cgr_actual", "n_translocation")))
  expect_false("canopy_senescence" %in% .strip_slice(names(g)))
  endpoint <- query_graph("barley", treatments = "N_application", outcomes = "cgr_actual",
    time_unroll = "selective", label = NULL, quiet = TRUE,
    params = query_params(show_hidden_nodes = FALSE))
  expect_true("cgr_actual" %in% dagitty::outcomes(endpoint))
  expect_error(fold_query(columns = "spike"), "Cannot recover requested node.*spike")
  expect_error(fold_query(columns = "unspecified_organ_live_n"), "Cannot recover requested node")
})

test_that("coarse grouping cannot delete allowed members with an excluded row", {
  g <- slice_external_query(prohibit_left = TRUE,
    params = query_params(apsim_nodes_detail = "subnode", force_acyclic = TRUE))
  r <- attr(g, "query_record")
  expect_false("toy.soil.quantity.left" %in% r$node_map$member_id)
  expect_true("toy.soil.quantity.right" %in% r$node_map$member_id)
  expect_false("toy.soil.quantity" %in% r$node_map$quantity)
  expect_setequal(r$recovery_folding$refined_members,
    c("toy.soil.quantity.left", "toy.soil.quantity.right"))
})
