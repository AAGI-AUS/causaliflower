test_that("query defaults retain configured context and finest node ids", {
  p <- query_params()
  expect_identical(p$query_depth, Inf)
  expect_identical(p$apsim_nodes_detail, "variable")
  expect_true(p$show_hidden_nodes)
  expect_false(p$node_recombination)
  ids <- c("sample.organ.measurement.first", "sample.organ.measurement.second")
  expect_identical(unname(.apsim_truncate(ids, p$apsim_nodes_detail)), ids)
  expect_identical(unname(.apsim_truncate(ids, "subnode")), rep("sample.organ.measurement", 2L))
  expect_error(query_params(query_depth = NA_real_), "query_depth")
})

test_that("default temporal query equals its explicit broad configuration", {
  g <- query_graph("barley", treatments = "nitrogen_applied",
                   outcomes = "grain_protein", quiet = TRUE, label = NULL)
  explicit <- query_graph("barley", treatments = "nitrogen_applied",
                          outcomes = "grain_protein", quiet = TRUE, label = NULL,
                          time_unroll = "all", randomised_treatment = FALSE,
                          params = query_params(query_depth = Inf, apsim_nodes_detail = "variable",
                                                show_hidden_nodes = TRUE, node_recombination = FALSE))
  expect_identical(attr(g, "dag_text"), attr(explicit, "dag_text"))
  expect_identical(attr(g, "time_unroll"), "all")
  expect_identical(attr(g, "crop"), "barley")
  expect_identical(attr(g, "water_model"), "SoilWater")
  expect_setequal(attr(g, "recovered_vertices"), names(g))
  expect_false(attr(g, "randomised_treatment"))
  expect_length(attr(g, "randomised_parents_cut"), 0L)
  expect_length(attr(g, "latent"), 0L)
  expect_true(.is_acyclic(attr(g, "recovered_vertices"), attr(g, "edges")))
})

test_that("the reduced query refuses an absent alternative water family", {
  expect_error(query_graph("barley", treatments = "nitrogen_applied",
                   outcomes = "grain_protein", quiet = TRUE,
                   params = query_params(water_model = "SWIM")), "bundled slice has no SWIM family rows")
})

test_that("historical compact settings and complete ancestry remain explicit", {
  compact <- query_graph("barley", treatments = "nitrogen_applied",
                         outcomes = "grain_protein", quiet = TRUE, label = NULL,
                         time_unroll = "pre_post",
                         params = query_params(query_depth = 1L, show_hidden_nodes = FALSE))
  ancestral <- query_graph("barley", treatments = "nitrogen_applied",
                           outcomes = "grain_protein", quiet = TRUE, label = NULL,
                           params = query_params(query_depth = 2L))
  expect_identical(attr(compact, "time_unroll"), "pre_post")
  expect_false(attr(compact, "show_hidden_nodes"))
  expect_identical(attr(compact, "query_depth"), 1L)
  expect_identical(attr(ancestral, "query_depth"), 2L)
  expect_true(.is_acyclic(attr(compact, "recovered_vertices"), attr(compact, "edges")))
  expect_true(.is_acyclic(attr(ancestral, "recovered_vertices"), attr(ancestral, "edges")))
})

test_that("SCM conversion retains all selected variables including isolated and hidden rows", {
  tab <- data.frame(node_id = c("X", "Hidden", "Y", "Isolated", "Management"),
                    function_of = c("NA", "X", "Hidden", "NA", "NA"),
                    hidden_node = c(FALSE, TRUE, FALSE, FALSE, FALSE),
                    tier = c("", "", "", "", "management"))
  g <- scm_to_dag(tab)
  expect_setequal(attr(g, "recovered_vertices"), tab$node_id)
  expect_identical(attr(g, "latent"), "Hidden")
  expect_setequal(paste(attr(g, "edges")$from, attr(g, "edges")$to), c("X Hidden", "Hidden Y"))
  expect_length(dagitty::exposures(g), 0L)
  expect_length(dagitty::outcomes(g), 0L)
  expect_identical(tab$function_of[2], "X")
})

test_that("SCM conversion selects edges for the crop and alternative water families", {
  tab <- data.frame(node_id = c("X", "Z", "Y", "SoilWater.sample_quantity", "SWIM.sample_quantity"),
                    function_of = c("NA", "NA", "X[wheat];Z[barley]", "NA", "NA"))
  expect_error(scm_to_dag(tab), "Select crop")
  expect_error(scm_to_dag(tab, crop = "wheat"), "Select water_model")
  wheat <- scm_to_dag(tab, crop = "wheat", water_model = "SoilWater")
  barley <- scm_to_dag(tab, crop = "barley", water_model = "SWIM")
  expect_identical(attr(wheat, "edges")$from, "X")
  expect_identical(attr(barley, "edges")$from, "Z")
  expect_true("SoilWater.sample_quantity" %in% attr(wheat, "recovered_vertices"))
  expect_false("SWIM.sample_quantity" %in% attr(wheat, "recovered_vertices"))
  expect_false("SoilWater.sample_quantity" %in% attr(barley, "recovered_vertices"))
})

test_that("SCM conversion preserves lag magnitude and a finite declared horizon", {
  tab <- data.frame(node_id = c("X", "Y"), function_of = c("NA", "X@t-2"))
  expect_error(scm_to_dag(tab), "Supply time_slices")
  expect_error(scm_to_dag(tab, time_slices = 2L), "largest retained lag")
  g <- scm_to_dag(tab, time_slices = 4L)
  expect_length(attr(g, "recovered_vertices"), 8L)
  expect_setequal(paste(attr(g, "edges")$from, attr(g, "edges")$to), c("X@t0 Y@t2", "X@t1 Y@t3"))
  expect_identical(attr(g, "lagged_edges")$lag, 2L)
  expect_true(.is_acyclic(attr(g, "recovered_vertices"), attr(g, "edges")))
})

test_that("SCM conversion includes origin and stock persistence", {
  tab <- data.frame(node_id = c("Initial", "Stock"), function_of = c("NA", "NA"),
                    init_of = c("", "Initial"), path_dependent = c(FALSE, TRUE))
  g <- scm_to_dag(tab, time_slices = 3L)
  expect_setequal(paste(attr(g, "edges")$from, attr(g, "edges")$to),
                  c("Initial@t0 Stock@t0", "Stock@t0 Stock@t1", "Stock@t1 Stock@t2"))
})

test_that("SCM conversion rejects unresolved, ambiguous and cyclic structure", {
  cycle <- data.frame(node_id = c("X", "Y"), function_of = c("Y", "X"))
  expect_error(scm_to_dag(cycle), "contemporaneous graph contains a cycle")
  expect_error(scm_to_dag(transform(cycle, function_of = c("NA", "Missing"))), "resolve to no node")
  ambiguous <- data.frame(node_id = c("X", "Y", "Z"), concept_id = c("c:C", "c:C", "c:Z"),
                          function_of = c("NA", "NA", "c:C"))
  expect_error(scm_to_dag(ambiguous), "multiple nodes")
  origin <- data.frame(node_id = c("X", "Y"), function_of = c("NA", "X"), init_of = c("Y", ""))
  expect_error(scm_to_dag(origin, time_slices = 2L), "including initial parents")
})

test_that("qualifiers are filtered before checking unresolved parents", {
  tab <- data.frame(node_id = c("X", "Y"), function_of = c("NA", "X[wheat];Missing[barley]"))
  expect_silent(scm_to_dag(tab, crop = "wheat"))
  expect_error(scm_to_dag(tab, crop = "barley"), "resolve to no node")
})

test_that("SCM conversion keeps literal ids with suffixes that resemble time indices", {
  tab <- data.frame(node_id = c("X", "State_t1"), function_of = c("NA", "X"),
                    hidden_node = c(FALSE, TRUE))
  expect_identical(attr(scm_to_dag(tab), "latent"), "State_t1")
  g <- scm_to_dag(tab, time_slices = 2L)
  expect_setequal(attr(g, "latent"), c("State_t1@t0", "State_t1@t1"))
})

test_that("selected water configuration resolves a shared concept to its active row", {
  tab <- data.frame(node_id = c("SoilWater.sample_quantity", "SWIM.sample_quantity", "Y"),
                    concept_id = c("c:Water", "c:Water", "c:Yield"),
                    function_of = c("NA", "NA", "c:Water"))
  soilwater <- scm_to_dag(tab, water_model = "SoilWater")
  swim <- scm_to_dag(tab, water_model = "SWIM")
  expect_identical(attr(soilwater, "edges")$from, "SoilWater.sample_quantity")
  expect_identical(attr(swim, "edges")$from, "SWIM.sample_quantity")
})

test_that("SCM conversion rejects conflicting concept and node identifiers", {
  tab <- data.frame(node_id = c("C", "A", "B", "Y"),
                    concept_id = c("c:Other", "c:C", "c:C", "c:Y"),
                    function_of = c("NA", "NA", "NA", "c:C"))
  expect_error(scm_to_dag(tab), "concept/node_id collision")
  expect_error(scm_to_dag(tab[ -3L, ]), "concept/node_id collision")
  identity <- data.frame(node_id = c("X", "Y"), concept_id = c("c:X", "c:Y"),
                         function_of = c("NA", "c:X"))
  expect_identical(attr(scm_to_dag(identity), "edges")$from, "X")
})
