variable_query <- function(edges, ids = NULL, exposures = "T") {
  graph <- dagitty::dagitty(paste("dag {", paste(exposures, "[exposure];", collapse = " "),
    "Y [outcome];", edges, "}"))
  nodes <- names(graph)
  if (is.null(ids)) ids <- stats::setNames(nodes, nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = NULL, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = nodes, quantity = unname(ids[nodes]),
      member_id = unname(ids[nodes]), occurrence = nodes),
    policy = data.frame(node_id = unique(unname(ids)), recovered_node = "yes"))
  graph
}

same_size_query <- function() {
  variable_query("T -> Y; T -> P_t0 -> P_t1",
    c(T = "T", Y = "Y", P_t0 = "P", P_t1 = "P"))
}

test_that("literal plain node names are separate variables and exact requests stay exact", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; P_t0 -> T; P_t0 -> Y; P_t1 -> T; P_t1 -> Y; T -> Y }")
  expect_identical(.cd_variable_count(graph), 4L)
  expect_identical(.cd_match_nodes(graph, "P_t0"), "P_t0")
  expect_error(.cd_match_nodes(graph, "P"), "absent")
  built <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  expect_setequal(names(built), names(graph))
  expect_identical(dag_options(built)$level$variables, 4L)
  expect_setequal(dag_options(built)$quantities$total$members, c("P_t0", "P_t1"))
})

test_that("source maps count canonical members and refuse fabricated identities", {
  graph <- same_size_query()
  expect_identical(.cd_variable_count(graph), 3L)
  record <- attr(graph, "query_record")
  record$requests$columns$resolved <- c(alias = "P")
  attr(graph, "query_record") <- record
  expect_setequal(.cd_match_nodes(graph, "alias"), c("P_t0", "P_t1"))
  record$node_map$quantity[record$node_map$member_id == "P"] <- "printed-combination"
  extra <- record$node_map[record$node_map$node == "P_t0", ]
  extra$member_id <- "other-source"
  record$node_map <- rbind(record$node_map, extra)
  attr(graph, "query_record") <- record
  expect_identical(.cd_variable_count(graph), 4L)
  record$node_map$member_id[[1L]] <- NA_character_
  attr(graph, "query_record") <- record
  expect_error(construct_dag(graph, quiet = TRUE), "node map")
  expect_error(.cd_variable_map(graph), "identity map")
  record$node_map <- record$node_map[record$node_map$node != record$node_map$node[[1L]], ]
  attr(graph, "query_record") <- record
  expect_error(.cd_variable_map(graph), "identity map")
})

test_that("same-variable alternatives retain exact temporal membership and stable numbers", {
  graph <- same_size_query()
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  full <- construct_dag(graph, quiet = TRUE)
  choices <- dag_options(full)$alternatives
  expect_identical(dag_options(full)$level$variables, 3L)
  expect_identical(choices$select, 1:2)
  expect_identical(choices$selected, c(TRUE, FALSE))
  expect_true(dag_options(full)$alternatives_complete)
  expect_true(all(choices$static_available))
  expect_setequal(names(temporal_dag(full)), c("T", "Y", "P_t0", "P_t1"))
  other <- select_dag(full, select = 2, quiet = TRUE)
  expect_setequal(names(temporal_dag(other)), c("T", "Y", "P_t0"))
  expect_identical(dag_options(other)$level$variables, 3L)
  expect_identical(dag_options(other)$level$select, 2L)
  expect_identical(.gd_signature(select_dag(other, quiet = TRUE)), .gd_signature(other))
  expect_identical(.gd_signature(temporal_dag(select_dag(other, quiet = TRUE))),
    .gd_signature(temporal_dag(other)))
  expect_identical(.gd_signature(temporal_dag(eval(dag_call(other)))), .gd_signature(temporal_dag(other)))
  expect_identical(dag_options(construct_dag(graph, Inf, quiet = TRUE))$candidate,
    dag_options(full)$candidate)
  .cf_cache$clear()
  expect_identical(dag_options(construct_dag(graph, quiet = TRUE))$alternatives, choices)
  withr::local_options(causaliflower.accelerate = FALSE)
  expect_identical(dag_options(construct_dag(graph, quiet = TRUE))$alternatives, choices)
})

test_that("size and strict resize operate in variables with deterministic snapping", {
  graph <- variable_query("T -> A -> B -> Y")
  small <- construct_dag(graph, 0, quiet = TRUE)
  large <- construct_dag(graph, quiet = TRUE)
  expect_identical(dag_options(small)$level$variables, 4L)
  expect_identical(prune_dag(small, 1, quiet = TRUE), small)
  expect_identical(grow_dag(large, 1, quiet = TRUE), large)
  graph <- same_size_query()
  large <- construct_dag(graph, quiet = TRUE)
  small <- prune_dag(large, 1, quiet = TRUE)
  expect_identical(dag_options(small)$level$variables, 2L)
  expect_identical(dag_options(grow_dag(small, 1, quiet = TRUE))$level$variables, 3L)
  expect_warning(tied <- construct_dag(graph, 2.5, quiet = TRUE), "equally close")
  expect_identical(dag_options(tied)$level$variables, 3L)
  other <- select_dag(large, select = 2, quiet = TRUE)
  expect_identical(grow_dag(other, 1, quiet = TRUE), other)
  expect_identical(dag_options(select_dag(other, 2, quiet = TRUE))$level$variables, 2L)
  expect_no_warning(construct_dag(graph, 1e300, quiet = TRUE))
  expect_no_warning(prune_dag(large, 1e300, quiet = TRUE))
  expect_no_warning(grow_dag(small, 1e300, quiet = TRUE))
})

test_that("view changes preserve candidate identity and defaults are static", {
  graph <- same_size_query()
  static <- construct_dag(graph, quiet = TRUE)
  temporal <- select_dag(static, unroll = TRUE, quiet = TRUE)
  expect_identical(dag_options(temporal)$candidate, dag_options(static)$candidate)
  expect_identical(dag_options(temporal)$level$variables, 3L)
  expect_identical(dag_options(temporal)$level$nodes, 4L)
  expect_true(dag_options(temporal)$arguments$unroll)
  expect_identical(.gd_signature(temporal), .gd_signature(temporal_dag(temporal)))
  expect_true(dag_options(construct_dag(temporal, quiet = TRUE))$arguments$unroll)
  expect_true(dag_options(grow_dag(prune_dag(temporal, quiet = TRUE), quiet = TRUE))$arguments$unroll)
  expect_identical(.gd_signature(prune_dag(temporal, 0, unroll = FALSE, quiet = TRUE)), .gd_signature(static))
  expect_identical(.gd_signature(eval(dag_call(temporal))), .gd_signature(temporal))
  expect_true(dag_options(grow_dag(static, 1, unroll = TRUE, quiet = TRUE))$arguments$unroll)
  endpoint <- variable_query("T_t0 -> T_t1 -> Y",
    c(T_t0 = "T", T_t1 = "T", Y = "Y"), exposures = c("T_t0", "T_t1"))
  expect_error(suppressWarnings(construct_dag(endpoint, quiet = TRUE)), "unroll = TRUE")
  allowed <- suppressWarnings(construct_dag(endpoint, unroll = TRUE, quiet = TRUE))
  expect_identical(dag_options(allowed)$level$variables, 2L)
  expect_false(dag_options(allowed)$level$static_available)
  attr(endpoint, "cf_node_map") <- attr(endpoint, "query_record")$node_map
  attr(endpoint, "query_record") <- NULL
  expect_error(suppressWarnings(construct_dag(endpoint, quiet = TRUE)), "unroll = TRUE")
})

test_that("selectors are validated and old objects migrate only by exact identity", {
  graph <- same_size_query()
  dag <- construct_dag(graph, quiet = TRUE)
  for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(select_dag(dag, select = bad, quiet = TRUE), "positive whole")
  }
  expect_error(select_dag(dag, select = 3, quiet = TRUE), "Available alternatives.*1, 2")
  expect_error(select_dag(graph), "construct_dag")
  plan <- construct_dag(graph, dry_run = TRUE, quiet = TRUE)
  expect_s3_class(select_dag(plan, select = 2, quiet = TRUE), "cf_dag")
  expect_s3_class(eval(dag_call(plan)), "cf_dag_plan")
  old <- .gd_state(dag)
  old$construction$metric <- NULL
  old$construction$candidate <- NULL
  old$construction$level$nodes <- 999L
  old <- .gd_attach(dag, old)
  expect_identical(.gd_signature(temporal_dag(select_dag(old, quiet = TRUE))), .gd_signature(temporal_dag(dag)))
  expect_identical(.gd_signature(eval(dag_call(old))), .gd_signature(dag))
  call <- dag_call(dag)
  call$candidate$binding$version <- "unavailable"
  expect_error(eval(call), "unsupported construction version")
})

test_that("plain latent declarations are inherited unless explicitly overridden", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; Z [latent]; Z -> T; Z -> Y; T -> Y }")
  inherited <- suppressWarnings(construct_dag(graph, quiet = TRUE))
  observed <- suppressWarnings(construct_dag(graph, latents = NULL, quiet = TRUE))
  expect_identical(dagitty::latents(inherited), "Z")
  expect_length(dagitty::latents(observed), 0L)
  expect_identical(dag_options(inherited)$level$variables, dag_options(observed)$level$variables)
  query <- variable_query("Z -> T; Z -> Y; T -> Y")
  record <- attr(query, "query_record")
  dagitty::latents(query) <- "Z"
  attr(query, "query_record") <- record
  omitted <- construct_dag(query, quiet = TRUE)
  cleared <- construct_dag(query, latents = NULL, quiet = TRUE)
  expect_length(dagitty::latents(omitted), 0L)
  expect_identical(.gd_signature(omitted), .gd_signature(cleared))
  expect_identical(dagitty::latents(original_query(omitted)), "Z")
})

test_that("replay keeps an unused ranked budget independent of later options", {
  withr::local_options(causaliflower.ranked_steps = "unused")
  dag <- construct_dag(same_size_query(), quiet = TRUE)
  saved <- dag_call(dag)
  options(causaliflower.ranked_steps = 20000L)
  expect_identical(.gd_signature(temporal_dag(eval(saved))), .gd_signature(temporal_dag(dag)))
})

test_that("refused rule families never claim complete alternatives or session-dependent maxima", {
  graph <- variable_query("R -> T; R -> Y; T -> Y")
  record <- attr(graph, "query_record")
  row <- record$node_map[record$node_map$node == "R", ]
  row$member_id <- "R-other"
  record$node_map <- rbind(record$node_map, row)
  attr(graph, "query_record") <- record
  .cf_cache$clear()
  expect_error(construct_dag(graph, quiet = TRUE), "several source members")
  dag <- construct_dag(graph, 0, quiet = TRUE)
  expect_false(dag_options(dag)$alternatives_complete)
  expect_false(dag_options(dag)$sizes_complete)
  expect_match(dag_options(dag)$alternatives_limit, "several source members")
  expect_identical(dag_options(dag)$level$variables, 4L)
  expect_error(construct_dag(graph, quiet = TRUE), "several source members")
  withr::local_options(causaliflower.accelerate = FALSE)
  expect_error(construct_dag(graph, quiet = TRUE), "several source members")
})
