runtime_policy_nodes <- function() { # Build a minimal table with independent numbered concepts
  ids <- c(names(.cf_recovery_policy()$rules), "soil_thermal_regime", "ordinary")
  return(data.frame(node_id = ids, concept_id = paste0("c:", seq_along(ids)),
    source_field = seq_along(ids), stringsAsFactors = FALSE))
}

test_that("runtime policy preserves raw metadata and the reviewed interpretation", {
  raw <- runtime_policy_nodes()
  before <- serialize(raw, NULL)
  effective <- .cf_effective_nodes(raw)
  policy <- .cf_recovery_policy()
  expect_identical(policy$version, "dagby-recovery-20261019-v1")
  expect_identical(effective[names(raw)], raw)
  expect_identical(serialize(raw, NULL), before)
  expect_identical(effective$recovered_node,
    c("yes:user", "no:Class-level live N does not identify an organ instance",
      "yes:user", "yes:user", "yes:user",
      "no:Whole organ label does not specify a numeric quantity", "yes:user", "yes", "yes"))
  expect_identical(.cf_effective_nodes(effective), effective)
  effective$source_field[1L] <- -1L
  expect_identical(serialize(raw, NULL), before)
})

test_that("runtime policy preserves authoritative metadata and rejects invalid prepared inputs", {
  raw <- runtime_policy_nodes()
  raw$recovered_node <- rep("yes", nrow(raw))
  before <- serialize(raw, NULL)
  expect_identical(.cf_effective_nodes(raw), raw)
  expect_identical(serialize(raw, NULL), before)
  raw <- .cf_effective_nodes(runtime_policy_nodes())
  raw$recovered_node[1L] <- NA_character_
  expect_error(.cf_effective_nodes(raw), "Invalid recovered_node policy")
  raw$recovered_node[1L] <- "no:   "
  expect_error(.cf_effective_nodes(raw), "Invalid recovered_node policy")
  raw$recovered_node[1L] <- "unsupported"
  expect_error(.cf_effective_nodes(raw), "Invalid recovered_node policy")
  expect_error(.cf_effective_nodes(runtime_policy_nodes()[-1L, ]), "Reviewed policy quantities are missing.*cgr_actual")
  raw <- runtime_policy_nodes()
  raw$concept_id[1L] <- raw$concept_id[2L]
  expect_error(.cf_effective_nodes(raw), "Invalid node identity column: concept_id")
})

test_that("load hook prepares namespace and exported lazy data from an independent copy", {
  raw <- runtime_policy_nodes()
  before <- serialize(raw, NULL)
  namespace <- new.env(parent = emptyenv())
  lazydata <- new.env(parent = emptyenv())
  lazydata$apsim_dag_nodes <- raw
  expect_invisible(.cf_bind_recovery_policy(namespace, lazydata))
  expect_identical(namespace$apsim_dag_nodes, .cf_effective_nodes(raw))
  expect_identical(lazydata$apsim_dag_nodes, namespace$apsim_dag_nodes)
  expect_identical(serialize(raw, NULL), before)
  expect_invisible(.cf_bind_recovery_policy(namespace, lazydata))
  expect_identical(namespace$apsim_dag_nodes, .cf_effective_nodes(raw))
})

test_that("loaded package exposes the same authoritative policy without mutating shared data", {
  namespace <- asNamespace("causaliflower")
  exported <- getExportedValue("causaliflower", "apsim_dag_nodes")
  expect_identical(get("apsim_dag_nodes", envir = namespace, inherits = FALSE), exported)
  expect_identical(get("apsim_dag_nodes", envir = getNamespaceInfo(namespace, "lazydata"),
    inherits = FALSE), exported)
  expect_identical(.cf_effective_nodes(exported), exported)
  before <- serialize(exported, NULL)
  independent <- .cf_effective_nodes(exported)
  independent$node_name_readable[1L] <- "Changed only in the local copy"
  expect_identical(serialize(exported, NULL), before)
  expect_identical(get("apsim_dag_nodes", envir = namespace, inherits = FALSE), exported)
  expect_false(identical(independent, exported))
})
