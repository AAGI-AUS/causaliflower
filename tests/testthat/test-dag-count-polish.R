count_polish_graph <- function(edges) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];", edges, "}"))
  ids <- .strip_slice(names(graph))
  attr(graph, "query_record") <- list(
    arguments = list(crop = "barley", params = list(water_model = "SoilWater", query_depth = Inf)),
    node_map = data.frame(node = names(graph), quantity = ids, member_id = ids,
                          occurrence = names(graph)),
    policy = data.frame(node_id = unique(ids), recovered_node = "yes")
  )
  graph
}


test_that("distinct endpoint occasions need explicit unrolling, and a partly declared quantity is marked", {
  graph <- count_polish_graph("T -> A -> Y; Y_t0 -> A; T -> Y")
  expect_error(dagby_construct(graph, keep_nodes = "Y_t0", quiet = TRUE, size = "largest"),
    "distinct occasions of an endpoint", class = "cd_static_unavailable")
  result <- dagby_construct(graph, keep_nodes = "Y_t0", quiet = TRUE, unroll = TRUE, size = "largest")
  expect_false(dagby_options(result)$summary$applied)
  expect_match(dagby_options(result)$level$static_reason, "distinct occasions of an endpoint")
  expect_identical(.gd_signature(result), .gd_signature(dagby_temporal(result)))
  mixed <- suppressWarnings(dagby_construct(count_polish_graph("T -> A_t0 -> Y; A_t1 -> A_t0; T -> Y"),
                                          latents = "A_t1", keep_nodes = "A_t1", quiet = TRUE, size = "largest"))
  expect_true(dagby_options(mixed)$summary$applied)
  expect_identical(dagitty::latents(mixed), "A")
  expect_identical(dagitty::latents(dagby_temporal(mixed)), "A_t1")
  expect_identical(dagby_options(mixed)$unmeasured$partial, list(A = list(unmeasured = "A_t1", observed = "A_t0")))
})


# Refinement must find a smaller set valid on the source, not just shrink a picture.
