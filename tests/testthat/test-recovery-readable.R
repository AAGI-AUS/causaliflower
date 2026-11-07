test_that("readable_dag reads recovered query attributes without altering its input", {
  graph <- dagitty::dagitty(paste(
    "dag { Treatment [exposure]; CropN [outcome]; SoilN_t0 [latent];",
    "Treatment -> SoilN_t1; SoilN_t0 -> SoilN_t1; SoilN_t1 -> CropN }"
  ))
  attr(graph, "recovered_vertices") <- names(graph)
  attr(graph, "edges") <- data.table::data.table(
    from = c("Treatment", "SoilN_t0", "SoilN_t1"),
    to = c("SoilN_t1", "SoilN_t1", "CropN"))
  attr(graph, "latent") <- "SoilN_t0"
  attr(graph, "treatments") <- "Treatment"
  attr(graph, "outcomes") <- "CropN"
  class(graph) <- c("dag_by_query", class(graph))
  original <- serialize(graph, NULL)
  nodes <- data.frame(node_id = c("Treatment", "SoilN", "CropN"),
    node_name_readable = c("fertiliser", "soil_N", "crop_N"))
  shown <- readable_dag(graph, nodes = nodes)
  expect_identical(serialize(graph, NULL), original)
  expect_setequal(names(shown), c("fertiliser", "soil_N_t0", "soil_N_t1", "crop_N"))
  expect_identical(dagitty::exposures(shown), "fertiliser")
  expect_identical(dagitty::outcomes(shown), "crop_N")
  expect_identical(dagitty::latents(shown), "soil_N_t0")

  # Inverting the display labels recovers the exact original directed edges.
  inverse <- c(fertiliser = "Treatment", soil_N_t0 = "SoilN_t0",
               soil_N_t1 = "SoilN_t1", crop_N = "CropN")
  shown_edges <- dagitty::edges(shown)
  recovered <- data.frame(v = unname(inverse[shown_edges$v]), e = shown_edges$e,
                          w = unname(inverse[shown_edges$w]))
  key <- function(edges) sort(paste(edges$v, edges$e, edges$w))
  expect_identical(key(recovered), key(dagitty::edges(graph)))
})
