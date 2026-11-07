test_that("runtime node metadata uses the recovery column and complete policy values", {
  nodes <- getExportedValue("causaliflower", "apsim_dag_nodes")
  expect_true("recovered_node" %in% names(nodes))
  expect_false("surfaced_node" %in% names(nodes))
  policy <- nodes$recovered_node
  expect_length(policy, nrow(nodes))
  expect_false(anyNA(policy))
  expect_true(all(policy %in% c("yes", "yes:user") |
    grepl("^no:[[:space:]]*[^[:space:]].*$", policy)))
})

test_that("recovery metadata retains the reviewed barley decisions and soil temperature", {
  nodes <- getExportedValue("causaliflower", "apsim_dag_nodes")
  expected <- c(
    "cgr_actual" = "yes:user",
    "unspecified_organ_live_n" = "no:Class-level live N does not identify an organ instance",
    "n_translocation" = "yes:user",
    "canopy_senescence" = "yes:user",
    "frac_radn_intercepted" = "yes:user",
    "spike" = "no:Whole organ label does not specify a numeric quantity",
    "soil_thermal_regime" = "yes",
    "canopy_conductance" = "yes:user"
  )
  expect_true(all(names(expected) %in% nodes$node_id))
  expect_identical(nodes$recovered_node[match(names(expected), nodes$node_id)],
                   unname(expected))
  expect_true(all(nodes$recovered_node[!nodes$node_id %in% names(expected)] == "yes"))
  # Recovery is independent of the availability fields. Soil surface remains
  # the scientific quantity's name, not terminology for graph recovery.
  expect_true(all(c("hidden_node", "observable", "apsim_obtainable", "measurable") %in% names(nodes)))
  expect_true("residue_cover" %in% nodes$node_id)
})

test_that("packaged node fields agree with the independently verified readable fixture", {
  nodes <- getExportedValue("causaliflower", "apsim_dag_nodes")
  csv <- testthat::test_path("fixtures", "slice-node-layer.csv")
  expect_true(file.exists(csv))
  source_nodes <- read.csv(csv, colClasses = "character", check.names = FALSE,
    na.strings = character(), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  attr(nodes, "cf_source_order") <- NULL
  expect_identical(nodes[names(source_nodes)], source_nodes)
  expect_true(all(c("path_dependent", "temporal_class") %in% names(nodes)))
})
